INCLUDE(${CMAKE_SOURCE_DIR}/IJMacros.txt)

cmake_minimum_required(VERSION 2.8)

#Change PROJECT_NAME to the name of your project
PROJECT(RandomWalk)

#The following lines are required to use Dart
#ENABLE_TESTING()
#INCLUDE(Dart)

#Declare any external dependencies that your project may have here.
#examples include: ITK, VTK, JPEG, PNG, OpenGL, ZLIB, Perl, Java
#If you're not sure what name to use, look in the Modules directory of your
#cmake install and check that a file named Find(Package).cmake exists
#
# The packages can be specified with a version number, for example:
#
# ITK 2.8.1
# ITK 3.2.0
#
# If no version is specified, the most recent release of the package
# will be used.
SET(Required_Packages
ITK
)

#this foreach loads all of the packages that you specified as required.
#It shouldn't need to be modified.
FOREACH(Package ${Required_Packages})
  LOADPACKAGE(${Package})
ENDFOREACH(Package)

SET (Libraries 
     ITKBasicFilters
	 ITKAlgorithms
     ITKCommon   
     ITKStatistics
     ITKIO
     ITKNumerics
	 )
	 
#the following block of code is an example of how to build an executable in
#cmake.  Unmodified, it will add an executable called "MyExe" to the project.
#MyExe will be built using the files MyClass.h and MyClass.cxx, and it will
#be linked to all the libraries you specified above. 
#You can build more than one executable per project
SET(CurrentExe "RandomWalkImageFilter")
ADD_EXECUTABLE(${CurrentExe} RandomWalkImageFilter.cxx) 
TARGET_LINK_LIBRARIES(${CurrentExe} ${Libraries})
SET(CurrentExe "RandomWalkTest")
ADD_EXECUTABLE(${CurrentExe} RandomWalkTest.cxx) 
TARGET_LINK_LIBRARIES(${CurrentExe} ${Libraries})

#the following line is an example of how to add a test to your project.
#Testname is the title for this particular test.  ExecutableToRun is the
#program which will be running this test.  It can either be a part of this
#project or an external executable.  After that list any args that are needed
#for this test.  Include as many tests as you like.  If your project doesn't have
#any tests you can comment out or delete the following line.
#ADD_TEST(VarianceImageFilter VarianceImageFilterTest${CMAKE_SOURCE_DIR}/images/sf4.png ${CMAKE_SOURCE_DIR}/images/sf4_test_out5.png 5)
#ADD_TEST(CompareImage ImageCompare ${CMAKE_SOURCE_DIR}/images/sf4_test_out5.png ${CMAKE_SOURCE_DIR}/images/sf4_out5.png)
#ADD_TEST(VarianceImageFilter VarianceImageFilterTest${CMAKE_SOURCE_DIR}/images/sf4.png ${CMAKE_SOURCE_DIR}/images/sf4_test_out1.png 1)
#ADD_TEST(CompareImage ImageCompare ${CMAKE_SOURCE_DIR}/images/sf4_test_out1.png ${CMAKE_SOURCE_DIR}/images/sf4_out1.png)
#Once you're done modifying this template, you should rename it to "CMakeLists.txt"

