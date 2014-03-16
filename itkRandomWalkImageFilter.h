/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkRandomWalkImageFilter.h,v $
  Language:  C++
  Date:      $Date: 2011-03-31 12:00:00 $
  Version:   $Revision: 0.1 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkRandomWalkImageFilter_h
#define __itkRandomWalkImageFilter_h

#include "itkImageToImageFilter.h"
#include "itkImage.h"

#define _SCL_SECURE_NO_DEPRECATE
#include "gmm.h"

namespace itk
{
/** \class RandomWalkImageFilter
 * \brief Random Walk Segmentation
 *
 * Computes an image where a given pixel is the probability of belonging to
 * the same label as the seed provided by the user.
 *
 */
template <class TInputImage, class TOutputImage>
class ITK_EXPORT RandomWalkImageFilter :
    public ImageToImageFilter< TInputImage, TOutputImage >
{
public: /* Define methods available to everyone */

   /** Standard class typedefs. */
  typedef RandomWalkImageFilter                                      Self;
  typedef ImageToImageFilter< TInputImage, TOutputImage> Superclass;
  typedef SmartPointer<Self>                                   Pointer;
  typedef SmartPointer<const Self>                             ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(RandomWalkImageFilter, ImageToImageFilter);

  /** ImageDimension enumeration. */
  itkStaticConstMacro(ImageDimension, unsigned int, TInputImage::ImageDimension);
  itkStaticConstMacro(OutputImageDimension, unsigned int, TOutputImage::ImageDimension);

  /** Typedef to describe the output image region type. */
  //typedef typename TInputImage::RegionType  InputImageRegionType;
  //typedef typename TOutputImage::RegionType OutputImageRegionType;

  /** Pixel related typedefs. */
  typedef typename TInputImage::PixelType               InputPixelType;
  typedef typename TOutputImage::PixelType              OutputPixelType;


  typedef typename TInputImage::SizeType    InputSizeType;
  typedef typename TInputImage::IndexType    IndexType;

  /** Public Methods */
  void SetSeed ( IndexType, int );


protected: /* Define methods available only to related classes */
  RandomWalkImageFilter();
  virtual ~RandomWalkImageFilter() {}

  void GenerateData();

  void ComputeLaplacian();

  void Solve();

  int sub2ind ( IndexType );
  double edge_weight( int , int );

  void PrintSelf(std::ostream& os, Indent indent) const;


private: /* Define methods available only to this class */
  RandomWalkImageFilter(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

  // Seed indexes
  int m_SeedFG, m_SeedBG;

  // Laplacian matrix
  gmm::csc_matrix<double> m_L;

  //Solution vector
  std::vector<double> m_X;
};
  
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkRandomWalkImageFilter.txx"
#endif

#endif
