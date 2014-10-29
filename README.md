Libwavelet_transform
----

This library written in C contains functions for computing
wavelet coefficients in dimension 1, dimension 2 and dimension 3:

* wavelet_transform1D()
* wavelet_transform2D()
* wavelet_transform3D()

and the corresponding inverse transforms

* invwavelet_transform1D()
* invwavelet_transform2D()
* invwavelet_transform3D()

there declarations is given in the header [wavelet_transform.h](libwavelets/wavelet_transform.h).

Building
----

The build process is governed by [cmake](http://www.cmake.org/).

#### Unix:
Start by going the the directory where you want your binaries and run

    ccmake PATH_TO_SOURCE
    make

To test the code run

    test

#### Windows

To build on windows, open the CMAKE gui, run configure-generate. Then open the project with Visual Studio.

## Build flags

The build can be controlled by changing build flags in the [cmake file](CMakeLists.txt).

#### Precision
All input vectors, coefficient vector and output vectors are assumed to be of the type FLOAT. Here FLOAT is defined to be a 32 bit float unless the macro HIGH_PRECISION is defined and in that case FLOAT is define as a 64. bit double

The macro HIGH_PRECISION will be defined at the compilations by uncommenting the line:

    #add_definitions(-DHIGH_PRECISION) 

Image data and Volume data as in DICOM are often written as arrays of unsigned 16 bits integers, but most often only 12 of those 16 bits are used. In order to use the wavelet transforms those arrays has to be converted by the user to arrays of type FLOATS.

After the inverse wavelet transform a similar conversion from FLOATs to unsigned integers may be done.

Note that if you want to run compile the test after changing the precision you have to first to rebuild the library.

#### Optimization

When building on UNIX type system (with e.g. make), the build mode is controlled by uncommenting the correct flags in the block 

```cmake
# Optimized 
set(CMAKE_BUILD_TYPE Release)
# Debug
# set(CMAKE_BUILD_TYPE Debug)
```

Tests
----
Some test code with randomly generated data can be found in the [test](test) directory.

Usage notes
----
IN-GOING ARRAYS ARE WRITTEN OVER

All wavelet_transforms use the input array to stor intermediate results in the calculations. If you want to keep the input vector results you have to make your own copy of the input vector.

The same holds for the in-going wavelet coefficient vector when running the inverse wavelet transforms.
