/*=========================================================================

  RandomWalkTest

  

  Hugues Fontenelle

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/


#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#endif

//#define DEBUG_VNL_to_CSV
#define DEBUG

// ITK
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
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

// VNL
#include <vnl/vnl_vector.h>
#include <vnl/vnl_sparse_matrix.h>
#include <vnl/algo/vnl_sparse_lu.h>

//
template <class T>
int round(T x) {
  return floor(x + 0.5);
}


/*--------------------------------------------
  my_get_row
  --------------------------------------------*/
template<class T>
vnl_vector<T> my_get_row(vnl_sparse_matrix<T> S,unsigned int row_index)
{
  vnl_vector<T> v( S.cols(), 0 );

  for (unsigned int j = 0; j != S.cols(); ++j) {    // For each element in row, unefficient way to read all values of a sparse matrix
    v[j] = S(row_index,j);
  }
  return v;
}

/*--------------------------------------------
  exportCSV
  quick hack for dumping a VNL matrix into 
  a comma-separated value file (CSV)
  --------------------------------------------*/
template <class T>  
void exportCSV( T A, char * filename)
{
  std::cout << "Export CSV, matrix size is " << A.rows() << 'x' << A.columns() << std::endl;

  std::ofstream file;
  file.open( filename );

  unsigned int itx, ity;
  for(itx = 0; itx < A.rows(); itx++) {
    for(ity = 0; ity < (A.columns()-1); ity++) {
      file << A(itx,ity) << ",";
    }
    file << A(itx,ity) << std::endl;
  }
  file.close();
}


/*--------------------------------------------
  subscript to index 
  --------------------------------------------*/
//int sub2ind( std::vector<int> size, std::vector<int> sub )
int sub2ind( itk::Image<int, 2>::SizeType size, itk::Image<int, 2>::IndexType sub )
{
  int ind;
  ind = sub[0] + sub[1] * size[0];
  return ind;
}

/*--------------------------------------------
  index to subscript
  --------------------------------------------*/
//std::vector<int> ind2sub( std::vector<int> size, int sub )
itk::Image<int, 2>::IndexType ind2sub( itk::Image<int, 2>::SizeType size, int sub )
{
  // Warning: Untested!!!
  itk::Image<int, 2>::IndexType ind;
  int r ,c, m;
  m = size[0];
  r = (sub % m); 
  c = floor((double) sub / m);
  ind[0] = r; ind[1]=c;
  return ind;
}




/*--------------------------------------------
  Edge weight
  --------------------------------------------*/
double edge_weight( int g1, int g2 )
{
  double beta = 90;
  double dg = ((double) g1)/256 - ((double) g2)/256;
  double weight;
  weight =  exp(-beta * dg * dg);
  
  return std::max(weight,1e-5);
}

/*--------------------------------------------
  Reshape
  --------------------------------------------*/
vnl_matrix<double> reshape(vnl_vector<double> X, int M, int N)
{
  vnl_matrix<double> A(M,N);
  
  int i,j,k=0;
  for(i=0;i!=M;++i) {
    for(j=0;j!=N;++j,++k) {
      A(i,j) = X(k);
    }
  }

  return A;
}

/*--------------------------------------------
  elim element at position  pos
  --------------------------------------------*/
vnl_vector<double> elim_ele(vnl_vector<double> b, int pos)
{
  int sz = b.size();
  vnl_vector<double> b_elim(sz-1);
  int it_in,it_out;
  for(it_in=0, it_out=0; it_in!=sz; ++it_in, ++it_out)
  { 
    if (it_in==pos) {
      ++it_in;
    }
    //else {
      b_elim(it_out)=b(it_in);
    //}
  }
  return b_elim;
}

/*--------------------------------------------
  insert element at position  pos
  --------------------------------------------*/
vnl_vector<double> insert_ele(vnl_vector<double> b_elim, int pos, double value)
{
  int sz = b_elim.size();
  vnl_vector<double> b(sz+1);
  int it_in,it_out;
  for(it_in=0, it_out=0; it_in!=sz; ++it_in, ++it_out)
  { 
    if (it_in==pos) {
      b(it_out) = value;
      ++it_out;
    }
    b(it_out)=b_elim(it_in);
  }
  return b;

}


/*--------------------------------------------
  elim element at position  pos
  --------------------------------------------*/
vnl_sparse_matrix<double> elim_ele(vnl_sparse_matrix<double> A, int pos)
{
  int M = A.rows();
  int N = A.columns();
  vnl_sparse_matrix<double> A_elim(M-1,N-1);
  int it_in_x,it_out_x,it_in_y,it_out_y;

  for(it_in_y=0, it_out_y=0; it_in_y!=N; ++it_in_y, ++it_out_y) {
    if(it_in_y==pos) {
      ++it_in_y;
    }
    for(it_in_x=0, it_out_x=0; it_in_x!=M; ++it_in_x, ++it_out_x) {
      if (it_in_x==pos) {
        ++it_in_x;
      }
      A_elim(it_out_x,it_out_y)=A(it_in_x,it_in_y);
    }
  }
  return A_elim;
}

/*--------------------------------------------
  Main 
  --------------------------------------------*/
int main( int argc, char *argv[])
{
  // Read input image
  if( argc < 1 )
  {
    std::cerr << "Missing Parameters " << std::endl;
    std::cerr << "Usage: " << argv[0];
    std::cerr << " inputImage outputImage posX posY" << std::endl;
    return 1;
  } 
  const char * inputFilename  = argv[1];
  typedef short int InputPixelType;
  const unsigned int Dimension = 2;
  typedef itk::Image< InputPixelType, Dimension > ImageType;

  typedef itk::ImageFileReader< ImageType > ReaderType;
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( inputFilename );
  reader->Update();
  ImageType::Pointer inputImage = reader->GetOutput();

  // init VNL sparse matrix
  unsigned long x = inputImage->GetLargestPossibleRegion().GetSize()[0];
  unsigned long y = inputImage->GetLargestPossibleRegion().GetSize()[1];
#ifdef DEBUG
  std::cout << "image size = " << x << "x" << y << std::endl;
#endif
  vnl_sparse_matrix<double> vnlMatrix(x*y,x*y);

  // Boundary Conditions
  typedef itk::ConstantBoundaryCondition< ImageType > BoundaryConditionType;
  BoundaryConditionType  BoundaryCondition ;
  BoundaryCondition.SetConstant( -1 );

  // Neigborhood
  typedef itk::ShapedNeighborhoodIterator< ImageType, BoundaryConditionType > ShapedNeighborhoodIteratorType;
  ShapedNeighborhoodIteratorType::OffsetType top = {{0,-1}};
  ShapedNeighborhoodIteratorType::OffsetType bottom = {{0,1}};
  ShapedNeighborhoodIteratorType::OffsetType left = {{-1,0}};
  ShapedNeighborhoodIteratorType::OffsetType right = {{1,0}};
  ShapedNeighborhoodIteratorType::OffsetType centre = {{0,0}};
  ShapedNeighborhoodIteratorType::RadiusType radius;
  radius.Fill( 1 );

  // Split region into faces
  typedef itk::NeighborhoodAlgorithm::ImageBoundaryFacesCalculator< ImageType > FaceCalculatorType;
  FaceCalculatorType faceCalculator;
  FaceCalculatorType::FaceListType faceList;
  FaceCalculatorType::FaceListType::iterator faceListIterator;
  faceList = faceCalculator( inputImage, inputImage->GetRequestedRegion(), radius );

  //
  // Process non-boundary regions normally.
  //
  faceListIterator=faceList.begin();

  //std::cout << "Face group 1" << std::endl;
  ShapedNeighborhoodIteratorType it( radius, reader->GetOutput(), *faceListIterator );

  it.ActivateOffset(top);
  it.ActivateOffset(bottom);
  it.ActivateOffset(left);
  it.ActivateOffset(right);

  ImageType::SizeType size = inputImage->GetRequestedRegion().GetSize();

#ifdef DEBUG
  std::cout << "Building large sparse matrix...";
#endif
  for (it.GoToBegin(); !it.IsAtEnd(); ++it)
  {
    vnlMatrix(sub2ind( size, it.GetIndex(centre)), sub2ind(size,it.GetIndex(left))) = -edge_weight(it.GetPixel(centre),it.GetPixel(left));
    vnlMatrix(sub2ind( size, it.GetIndex(centre)), sub2ind(size,it.GetIndex(right))) = -edge_weight(it.GetPixel(centre),it.GetPixel(right));   
    vnlMatrix(sub2ind( size, it.GetIndex(centre)), sub2ind(size,it.GetIndex(bottom))) = -edge_weight(it.GetPixel(centre),it.GetPixel(bottom));  
    vnlMatrix(sub2ind( size, it.GetIndex(centre)), sub2ind(size,it.GetIndex(top))) = -edge_weight(it.GetPixel(centre),it.GetPixel(top));
    vnlMatrix(sub2ind( size, it.GetIndex(centre)), sub2ind(size,it.GetIndex(centre))) = 4;

    //ShapedNeighborhoodIteratorType::ConstIterator ci;
    //int i;
    //for (ci = it.Begin(), i = 0; ci != it.End(); ci++, i++)  
    //{
    //  std::cout << ci.Get();     
    //}    
    //std::cout << std::endl;
  }

  //
  // Process each boundary region with boundary conditions.
  //
  for ( ++faceListIterator; faceListIterator != faceList.end(); ++faceListIterator)
  {
    //std::cout << "Face group " << *faceListIterator << std::endl;  

    ShapedNeighborhoodIteratorType it( radius, reader->GetOutput(), *faceListIterator );

    it.OverrideBoundaryCondition(&BoundaryCondition);

    it.ActivateOffset(top);
    it.ActivateOffset(bottom);
    it.ActivateOffset(left);
    it.ActivateOffset(right);

    for (it.GoToBegin(); !it.IsAtEnd(); ++it)
    {
      int d=0;
      if (it.GetPixel(left) >= 0) {
        vnlMatrix(sub2ind( size, it.GetIndex(centre)), sub2ind(size,it.GetIndex(left))) = -edge_weight(it.GetPixel(centre),it.GetPixel(left));
        d++;
      }
      if (it.GetPixel(right) >= 0) {
        vnlMatrix(sub2ind( size, it.GetIndex(centre)), sub2ind(size,it.GetIndex(right))) = -edge_weight(it.GetPixel(centre),it.GetPixel(right)); 
        d++;
      }      
      if (it.GetPixel(bottom) >= 0) {
        vnlMatrix(sub2ind( size, it.GetIndex(centre)), sub2ind(size,it.GetIndex(bottom))) = -edge_weight(it.GetPixel(centre),it.GetPixel(bottom));  
        d++;
      }      
      if (it.GetPixel(top) >= 0) {
        vnlMatrix(sub2ind( size, it.GetIndex(centre)), sub2ind(size,it.GetIndex(top))) = -edge_weight(it.GetPixel(centre),it.GetPixel(top));
        d++;
      }
      vnlMatrix(sub2ind( size, it.GetIndex(centre)), sub2ind(size,it.GetIndex(centre))) = d;
    }
  }
#ifdef DEBUG
  std::cout << "Done." << std::endl;
#endif

#ifdef DEBUG_VNL_to_CSV
  exportCSV <vnl_sparse_matrix<double>> (vnlMatrix, "data.csv");
#endif

  //
  // Set up RHS and solve system
  //
  // see http://public.kitware.com/vxl/doc/release/books/core/book_6.html#SEC44
#ifdef DEBUG
  std::cout << "Set up RHS..." << std::endl;
#endif
  const int posX  = atoi(argv[3]);
  const int posY  = atoi(argv[4]);
  itk::Image<int, 2>::IndexType sub; sub[0] = posX; sub[1] = posY;
  int indRHS = sub2ind( size, sub);

  vnl_vector<double> solution(x*y); 
  vnl_vector<double> b(x*y,0); 

  b = my_get_row <double> (vnlMatrix, indRHS); // lots of CPU lost here !!! :-(

  vnl_vector<double> solution_elim(x*y-1);
  vnl_vector<double> b_elim(x*y-1,0); 
  vnl_sparse_matrix<double> vnlMatrix_elim(x*y-1,x*y-1);

#ifdef DEBUG
  std::cout << "Eliminate element in Vector..." << std::endl;
#endif
  b_elim = elim_ele(b,indRHS);
#ifdef DEBUG
  std::cout << "Eliminate element in Matrix..." << std::endl;
#endif
  vnlMatrix_elim = elim_ele(vnlMatrix,indRHS);

#ifdef DEBUG_VNL_to_CSV
  exportCSV <vnl_sparse_matrix<double>> (vnlMatrix_elim, "vnlMatrix_elim.csv");
#endif

#ifdef DEBUG
  std::cout << "Solving large linear system..." << std::endl;
#endif

  vnl_sparse_lu linear_solver(vnlMatrix_elim, vnl_sparse_lu::estimate_condition);
  linear_solver.solve( -b_elim , &solution_elim);

  solution = insert_ele(solution_elim,indRHS,1);

  // export
  vnl_matrix<double> solM(x,y);
  solM = reshape(solution, x, y);

#ifdef DEBUG_VNL_to_CSV
  exportCSV <vnl_matrix<double>> (solM, "sol.csv");
#endif
 
#ifdef DEBUG
  std::cout << "Writing solution to image..." << std::endl;
#endif
  // set up output
  const char * outputFilename  = argv[2];
  typedef unsigned char OutputPixelType;

  typedef itk::Image< OutputPixelType, Dimension > OutputImageType;
  typedef itk::ImageFileWriter< OutputImageType > WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( outputFilename );
  
  //
  // copy vnl Matrix to Image
  //
  OutputImageType::Pointer outputImage = OutputImageType::New();
  // The image region should be initialized
  outputImage->SetRegions( inputImage->GetRequestedRegion() );
  outputImage->CopyInformation( inputImage );
  outputImage->Allocate();

  // actual copy
  OutputImageType::IndexType pixelIndex;
  OutputImageType::PixelType pixelValue; 
  int it_x,it_y;
  for(it_y=0; it_y!=y; ++it_y) {
    for(it_x=0; it_x!=x; ++it_x) {
      pixelIndex[0] = it_x;
      pixelIndex[1] = it_y; 
      pixelValue = (OutputPixelType) round(255*solM(it_x,it_y));
//#ifdef DEBUG
//      std::cout << " " << (int)pixelValue << " " << it_x << " " << it_y << std::endl;
//#endif
      outputImage->SetPixel( pixelIndex, pixelValue );
    }
  }
  // ---
 
  writer->SetInput( outputImage );

  try
  {
    writer->Update();
  }
  catch( itk::ExceptionObject & err )
  {
    std::cerr << "ExceptionObject caught !" << std::endl;
    std::cerr << err << std::endl;
    return EXIT_FAILURE;
  }

  return 0;
}
