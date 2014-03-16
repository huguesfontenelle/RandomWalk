/*=========================================================================

  RandomWalkImageFilter

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#endif

#ifdef __BORLANDC__
#define ITK_LEAN_AND_MEAN
#endif

#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

#include "itkRandomWalkImageFilter.h"



int main( int argc, char * argv[] )
{
  //Image Read

  if( argc < 7 )
    {
    std::cerr << "Usage: " << std::endl;
    std::cerr << argv[0] << "  inputImageFile  outputImageFile FGseedX FGseedY BGseedX BGseedY" << std::endl;
    return EXIT_FAILURE;
    }

  typedef   short int  InputPixelType;
  typedef unsigned char  OutputPixelType;

  typedef itk::Image< InputPixelType,  2 >   TInputImage;
  typedef itk::Image< OutputPixelType, 2 >   TOutputImage;

  typedef itk::ImageFileReader< TInputImage  >  ReaderType;
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[1] );
  reader->Update();

  // Set up filter
  typedef itk::RandomWalkImageFilter<
               TInputImage, TOutputImage >  FilterType;

  FilterType::Pointer filter = FilterType::New();
  
  filter->SetInput( reader->GetOutput() );

  // Set seeds
  TInputImage::IndexType sub;
  int label;
  sub[0] = atoi(argv[3]);
  sub[1] = atoi(argv[4]);
  label = 1;
  filter->SetSeed( sub, label );
  sub[0] = atoi(argv[5]);
  sub[1] = atoi(argv[6]);
  label = 0;
  filter->SetSeed( sub, label );


  try
    {
    filter->Update();
    }
  catch( itk::ExceptionObject & excp )
    {
    std::cerr << "RandomWalkImageFilter exception thrown." << std::endl;
    std::cerr << excp << std::endl;
    return EXIT_FAILURE;
    exit( EXIT_FAILURE );
    }

  // Image Write
  typedef itk::ImageFileWriter< TOutputImage >  WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( argv[2] );
  writer->SetInput( filter->GetOutput() );
  writer->Update();

  return EXIT_SUCCESS;
}

