#  README  file 

# for building the libraray libwavelet_transform.a



# This librarary  written in C  contain  functions  for computing
# wavelet coefficients   in dimension 1, dimensien 2 and dimension 3:
#
#	  wawelet_transform1D()
#  	  wawelet_transform2D()
#  	  wawelet_transform3D()
#
#
# and the corresponding inverse transforms
#
#     invwawelet_transform1D()
#     invwawelet_transform2D()
#     invwawelet_transform3D()


#   there declarations is given in "include/wavelet_transform.h"
#
#   in order to build the library  go to the build directory an

#   run      "make clean"     to remove any old objects file

#   run      "make libraries"     in order to compile all source
#   files in the directory  "source"  and collect and copy the
#    library libwavelet_transform.a  to the local top library  

#  in order to run a test go to a  the "test" directort
#    run first  "make clean"  and then  "make test"

#  and then you can run the test:
        
#      "./test"


#   all input vectors, coefficient vector and output vectors are
# assumed to be of the type  FLOAT

#  Here  FLOAT is defined to be a 32 bits  float
#  unless the macro  HIGH_PRECISION   is defined 
#  and in that case FLOAT is define as a 64 bits double

# The macro HIGH_PRECISION  will be defined at the complilations by
#  uncommenting  the line in build/makefile.setting : 
#                      #PRECISION = -dHIGH_PRECISION   


# Image data and Volume  data  as in DICOM  are often  written as
#  arrays of  unsigned 16 bits integers,
#  but most often only 12 of those 16 bits are used.
# In order to use the wavelet transforms those arrays has to
# be converted - by the user -  to arrays of type FLOATS 

# After the inverse wavelet transform a similar convertion
#  from FLOATs to unsigned integers may be done.
  

#  Note that if you want to run compile the test after changing the
#   macro PRECION  you have to first to rebuild the  library.




#  IMPORTANT  ABOUT MEMORY: IN-GOING ARRAYS are WRITTEN OVER:
#
#   All wavelet_transforms  use the input array to stor intermediate
#   results in the calculations. If you want to keep the input vector
#   results you have to make your own copy of the input vector.
#
#   The same holds for the in-going waveletkoeffisiont vector
#   when running the inverse wavelet transforms