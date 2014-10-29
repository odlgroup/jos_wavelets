#READE  file
#for a simple test program for the wavelet transform

# go to the directory "test"      (this directory)

# run  "make clean" to remove old objects file
# run "make test"   to compile  the test program   "test"

# Notice: the makefile includes the file "../build/makefile.settings"
#  where the  macro for defining HIGH_PRECISION  is set
#  if this file is changed you have to rebuild the wavelet transform
#  library  before compiling the test program.



#   This simple test program does the following.

#  1.The sizes of the volume is set by macros XLENGTH, YLENGTH, ZLENGTH

#  2 Random unsigned integers with 12 bits magnitude are generetade
#   and converted to a FLOAT arry of with size of the volume.
 
# 3. A backup copy of the volume is made for future comparation.

# 4. We run the 3D wavelet transform on the volume obtaining an
#     array of wavelet coeffisients of same size as the volume.

# 5 We run the 3D inverse wavelet transform on the array of wavelet
#   coefficients to reconstruct  an array  which should be very similar
#   to the volume array  that was generated before.

# 6  We run a program that collect  statistics from  the
#   backup copy of the generated , volume and the reconstructed volume 
#   and there difference vector.  

 
#  By changing some commenting  lines with "//" in the file test.c
#  it is easy to modify the 3D wavelet transform test
#   to a 2D test or a 1D   wavelet transform test