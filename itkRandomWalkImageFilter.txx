/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkRandomWalkImageFilter.txx,v $
  Language:  C++
  Date:      $Date: 2011-03-31 12:00:00 $
  Version:   $Revision: 0.1 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/

//#ifndef __itkRandomWalkImageFilter_txx
//#define __itkRandomWalkImageFilter_txx

#define DEBUG

#include "itkRandomWalkImageFilter.h"

// ITK
#include "itkShapedNeighborhoodIterator.h"
#include "itkImageRegionIterator.h"
#include "itkNeighborhoodAlgorithm.h"
#include "itkConstantBoundaryCondition.h"

// STD
#include <vector>
#include <iostream>
#include <math.h>
#include <stdio.h>
#include <algorithm>
#include <numeric>

// Gmm++
#ifndef _SCL_SECURE_NO_DEPRECATE
#define _SCL_SECURE_NO_DEPRECATE
#include <gmm.h>
#endif

/*--------------------------------------------
  Round 
  --------------------------------------------*/
template <class T>
int round(T x) {
  return floor(x + 0.5);
}


/*--------------------------------------------
  exportCSV
  quick hack for dumping a GMM matrix into 
  a comma-separated value file (CSV)
  --------------------------------------------*/
template <class T>
void exportCSV( T A, char * filename)
{
  std::cout << "Export CSV, matrix size is " << gmm::mat_nrows(A) << 'x' << gmm::mat_ncols(A) << std::endl;

  std::ofstream file;
  file.open( filename );

  unsigned int itx, ity;
  for(itx = 0; itx < gmm::mat_nrows(A); itx++) {
    for(ity = 0; ity < (gmm::mat_ncols(A)-1); ity++) {
      file << A(itx,ity) << ",";
    }
    file << A(itx,ity) << std::endl;
  }
  file.close();
}


namespace itk
{


/*--------------------------------------------
  RandomWalkImageFilter : Constructor
  --------------------------------------------*/
template <class TInputImage, class TOutputImage>
RandomWalkImageFilter<TInputImage, TOutputImage>
::RandomWalkImageFilter()
{


}


/*--------------------------------------------
  SetSeed
  --------------------------------------------*/
template< class TInputImage, class TOutputImage>
void
RandomWalkImageFilter< TInputImage, TOutputImage>
::SetSeed( IndexType sub, int label)
{
  int ind = sub2ind( sub );
  switch (label )
  {
    case 0: 
      this->m_SeedBG = ind;
      break;
    case 1: 
      this->m_SeedFG = ind;
      break;
    default:
      // error
      std::cout << "SetSeed error!" << std::endl;
  }
}


/*--------------------------------------------
  subscript to index 
  --------------------------------------------*/
template< class TInputImage, class TOutputImage>
int
RandomWalkImageFilter< TInputImage, TOutputImage>
::sub2ind( IndexType sub )
{
  InputSizeType size = this->GetInput()->GetRequestedRegion().GetSize();

  int ind;
  ind = sub[0] + sub[1] * size[0];
  return ind;
}

/*--------------------------------------------
  Edge weight
  --------------------------------------------*/
template< class TInputImage, class TOutputImage>
double
RandomWalkImageFilter< TInputImage, TOutputImage>
::edge_weight( int g1, int g2 )
{
  double beta = 90;
  double dg = ((double) g1)/256 - ((double) g2)/256;
  double weight;
  weight =  exp(-beta * dg * dg);
  
  return std::max(weight,1e-5);
}


/*--------------------------------------------
  ComputeLaplacian
  --------------------------------------------*/
template< class TInputImage, class TOutputImage>
void
RandomWalkImageFilter< TInputImage, TOutputImage>
::ComputeLaplacian( )
{
#ifdef DEBUG
  std::cout << "Computing Laplacian (large sparse matrix)..." << std::endl;
#endif

  // Get the input and output pointers;
  typename TOutputImage::Pointer output = this->GetOutput();
  typename TInputImage::ConstPointer input  = this->GetInput();

  InputSizeType size = input->GetRequestedRegion().GetSize();

  // Boundary Conditions
  typedef itk::ConstantBoundaryCondition< TInputImage > BoundaryConditionType;
  BoundaryConditionType  BoundaryCondition ;
  BoundaryCondition.SetConstant( -1 );

  // Neigborhood
  typedef itk::ShapedNeighborhoodIterator< TInputImage, BoundaryConditionType > ShapedNeighborhoodIteratorType;
  ShapedNeighborhoodIteratorType::OffsetType top = {{0,-1}};
  ShapedNeighborhoodIteratorType::OffsetType bottom = {{0,1}};
  ShapedNeighborhoodIteratorType::OffsetType left = {{-1,0}};
  ShapedNeighborhoodIteratorType::OffsetType right = {{1,0}};
  ShapedNeighborhoodIteratorType::OffsetType centre = {{0,0}};
  ShapedNeighborhoodIteratorType::RadiusType radius;
  radius.Fill( 1 );

  // Split region into faces
  typedef itk::NeighborhoodAlgorithm::ImageBoundaryFacesCalculator< TInputImage > FaceCalculatorType;
  FaceCalculatorType faceCalculator;
  FaceCalculatorType::FaceListType faceList;
  FaceCalculatorType::FaceListType::iterator faceListIterator;
  faceList = faceCalculator( input, input->GetRequestedRegion(), radius );

  //
  // Process non-boundary regions normally.
  //
  faceListIterator=faceList.begin();

  //std::cout << "Face group 1" << std::endl;
  ShapedNeighborhoodIteratorType it( radius, input, *faceListIterator );

  it.ActivateOffset(top);
  it.ActivateOffset(bottom);
  it.ActivateOffset(left);
  it.ActivateOffset(right);
  
  // create sparse matrix
  gmm::col_matrix< gmm::wsvector<double> > A_inserter(size[0]*size[1],size[0]*size[1]);

  double w1, w2, w3, w4;
  for (it.GoToBegin(); !it.IsAtEnd(); ++it)
  {
    //std::cout << "index centre " << it.GetIndex(centre) << " index left " << it.GetIndex(left) << " value " << -edge_weight(it.GetPixel(centre),it.GetPixel(left)) << std::endl;
    w1 = edge_weight( it.GetPixel(centre),it.GetPixel(left) );
    A_inserter(sub2ind( it.GetIndex(centre) ), sub2ind( it.GetIndex(left)   )) = -w1;
    w2 = edge_weight( it.GetPixel(centre),it.GetPixel(right) );
    A_inserter(sub2ind( it.GetIndex(centre) ), sub2ind( it.GetIndex(right)  )) = -w2; 
    w3 = edge_weight( it.GetPixel(centre),it.GetPixel(bottom) );
    A_inserter(sub2ind( it.GetIndex(centre) ), sub2ind( it.GetIndex(bottom) )) = -w3;  
    w4 = edge_weight( it.GetPixel(centre),it.GetPixel(top) );
    A_inserter(sub2ind( it.GetIndex(centre) ), sub2ind( it.GetIndex(top)    )) = -w4;
    A_inserter(sub2ind( it.GetIndex(centre) ), sub2ind( it.GetIndex(centre) )) = w1 + w2 + w3 + w4;
  }

  //
  // Process each boundary region with boundary conditions.
  //
  faceListIterator++;
  for ( ; faceListIterator != faceList.end(); ++faceListIterator)
  {
    ShapedNeighborhoodIteratorType it( radius, input, *faceListIterator );

    it.OverrideBoundaryCondition(&BoundaryCondition);

    it.ActivateOffset(top);
    it.ActivateOffset(bottom);
    it.ActivateOffset(left);
    it.ActivateOffset(right);

    for (it.GoToBegin(); !it.IsAtEnd(); ++it)
    {
      double w1=0, d=0;
      if (it.GetPixel(left) >= 0) {
        w1 = edge_weight(it.GetPixel(centre),it.GetPixel(left));
        d += w1;
        A_inserter(sub2ind( it.GetIndex(centre)), sub2ind( it.GetIndex(left) )) = -w1;
      }

      if (it.GetPixel(right) >= 0) {
        w1 = edge_weight(it.GetPixel(centre),it.GetPixel(right));
        d += w1;
        A_inserter(sub2ind( it.GetIndex(centre)), sub2ind( it.GetIndex(right) )) = -w1; 
      }     

      if (it.GetPixel(bottom) >= 0) {
        w1 = edge_weight(it.GetPixel(centre),it.GetPixel(bottom));
        d += w1;
        A_inserter(sub2ind( it.GetIndex(centre)), sub2ind( it.GetIndex(bottom) )) = -w1;  
      }  

      if (it.GetPixel(top) >= 0) {
        w1 = edge_weight(it.GetPixel(centre),it.GetPixel(top));
        d += w1;
        A_inserter(sub2ind( it.GetIndex(centre)), sub2ind( it.GetIndex(top) )) = -w1;
      }

      A_inserter(sub2ind( it.GetIndex(centre)), sub2ind( it.GetIndex(centre) )) = d;

    }
    // end face
  }

  // copy writable/inserted matrix into readable matrix
  gmm::clean(A_inserter, 1E-12);
  gmm::copy(A_inserter, this->m_L);

}


/*--------------------------------------------
  Solve large sparse Linear System
  --------------------------------------------*/
template< class TInputImage, class TOutputImage>
void
RandomWalkImageFilter< TInputImage, TOutputImage>
::Solve( )
{
  //
  // Set up RHS and solve system
  //
#ifdef DEBUG
  std::cout << "Set up Right Hand Side ..." << std::endl;
#endif

  InputSizeType size = this->GetInput()->GetRequestedRegion().GetSize();

  // index of unmarked points : index_u=0:prod(size)
  std::vector<int> index_u(size[0]*size[1]);
  std::vector<int>::iterator it_index;
  int i=0;
  for (it_index=index_u.begin(); it_index != index_u.end(); ++it_index, ++i) {
    index_u[i] = i;
  }
  // BUG: erase bigger first, finish with smaller
  index_u.erase (index_u.begin() + this->m_SeedBG);
  index_u.erase (index_u.begin() + this->m_SeedFG);
  

  // set up RHS of unmarked points: B_u
  gmm::rsvector<double> B_u(size[0]*size[1]-2); 
  gmm::copy( gmm::sub_vector(gmm::mat_col( (this->m_L) , this->m_SeedFG ), gmm::sub_index(index_u) ), B_u );
  gmm::scale( B_u, -1.0 );

  // set up Laplacian matrix of unmarked points: A_u
  gmm::row_matrix< gmm::rsvector<double> > A_u(size[0]*size[1]-2,size[0]*size[1]-2);
  gmm::copy( gmm::sub_matrix( (this->m_L), gmm::sub_index(index_u), gmm::sub_index(index_u) ), A_u );

  // Set up unknown vector (unmarked points): X_u
  std::vector<double> X_u(size[0]*size[1]-2);

  /** Iterative or Direct? */
  if (size[0]*size[1] < 401)
  {
#ifdef DEBUG
    std::cout << "DIRECT Solving large linear system..." << std::endl;
#endif
    /** Direct solvers*/
    //gmm::dense_matrix<double> A_u_dense(size[0]*size[1]-2,size[0]*size[1]-2);
    //gmm:copy(A_u, A_u_dense);
    gmm::lu_solve(A_u, X_u, B_u);
  }
  else
  {
    
#ifdef DEBUG
    std::cout << "Computing preconditioner..." << std::endl;
#endif
    // Optional scalar product for cg
    gmm::identity_matrix PS; 
    /** Optional preconditioner */
    gmm::identity_matrix PR;   
    //gmm::ildlt_precond< gmm::row_matrix< gmm::rsvector<double> > > PR(A_u); 

    gmm::iteration iter(10E-6);// Iteration object with the max residu
    size_t restart = 10;       // restart parameter for GMRES
#ifdef DEBUG
    std::cout << "ITERATIVE Solving large linear system..." << std::endl;
#endif
    /** Iterative solvers */
    gmm::cg(A_u, X_u, B_u, PS, PR, iter); // Conjugate gradient
    //gmm::bicgstab(A_u, X_u, B_u, PR, iter); // BICGSTAB BiConjugate Gradient Stabilized
    //gmm::gmres(A_u, X_u, B_u, PR, restart, iter); // GMRES generalized minimum residual
    //gmm::qmr(A_u, X_u, B_u, PR, iter); // Quasi-Minimal Residual method.
    //gmm::least_squares_cg(A_u, X_u, B_u, iter); // unpreconditionned least square CG.

  }

  // insertion of marked points into solution
  // BUG: insert smaller first, finish with bigger
  X_u.insert(X_u.begin() + this->m_SeedFG, 1);
  X_u.insert(X_u.begin() + this->m_SeedBG, 0);

  //gmm::copy(X_u , this->m_X );
  this->m_X = X_u;
}

/*--------------------------------------------
  Generate Data
  --------------------------------------------*/
template< class TInputImage, class TOutputImage>
void
RandomWalkImageFilter< TInputImage, TOutputImage>
::GenerateData( )
{
  // Allocate output
  this->AllocateOutputs();

  // Get the input and output pointers;
  typename TOutputImage::Pointer output = this->GetOutput();
  typename TInputImage::ConstPointer input  = this->GetInput();

  InputSizeType size = input->GetRequestedRegion().GetSize();

#ifdef DEBUG
  std::cout << "Image size = " << size[0] << "x" << size[1] << std::endl;
#endif 

  // Compute Laplacian
  this->ComputeLaplacian();

  // Solve
  this->Solve();

#ifdef DEBUG
  std::cout << "Export solution to image..." << std::endl;
#endif

  // reshape solution
  gmm::dense_matrix<double> solution(1, gmm::vect_size(this->m_X));
  gmm::copy(gmm::row_vector(this->m_X), solution);
  gmm::reshape(solution, size[0], size[1]);

  // write matrix solution to image output
  IndexType pixelIndex;
  OutputPixelType pixelValue; 
  int it_x,it_y;
  for(it_y=0; it_y!=size[1]; ++it_y) {
    for(it_x=0; it_x!=size[0]; ++it_x) {
      pixelIndex[0] = it_x;
      pixelIndex[1] = it_y; 
      pixelValue = (OutputPixelType) round(255*solution(it_x,it_y));
//#ifdef DEBUG
//      std::cout << " " << (int)pixelValue << " " << it_x << " " << it_y << std::endl;
//#endif
      output->SetPixel( pixelIndex, pixelValue );
    }
  }

  //exportCSV(this->m_L,"L.csv");
  


}

/**
 * Standard "PrintSelf" method
 */
template <class TInputImage, class TOutput>
void
RandomWalkImageFilter<TInputImage, TOutput>
::PrintSelf(
  std::ostream& os, 
  Indent indent) const
{
  Superclass::PrintSelf( os, indent );
  os << indent << "Seeds ForeGround: " << m_SeedFG << std::endl;
  os << indent << "Seeds BackGround: " << m_SeedBG << std::endl;
}

} // end namespace itk


