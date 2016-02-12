#include <boost/python.hpp>
#include "boost/python/extract.hpp"
#include "boost/python/numeric.hpp"
#include <iostream>
#include <vector>

extern "C" {
#include <lib/wavelet_transform.h>
}

using namespace boost::python;

// 1D, to be removed use wavelet_decompose3Py
void wavelet_transform1DPy(uintptr_t invector,
                           int xlength,
                           int Filterlength,
                           int levels,
                           uintptr_t waveletcoefficients) {
    wavelet_transform1D(reinterpret_cast<FLOAT*>(invector),
                        xlength,
                        Filterlength,
                        levels,
                        reinterpret_cast<FLOAT*>(waveletcoefficients));
}
//to be removed, use wavelet_reconstruct3Py
void invwavelet_transform1DPy(uintptr_t waveletcoefficients,
                              int xlength,
                              int Filterlength,
                              int levels,
                              uintptr_t outvector) {
    invwavelet_transform1D(reinterpret_cast<FLOAT*>(waveletcoefficients),
                           xlength,
                           Filterlength,
                           levels,
                           reinterpret_cast<FLOAT*>(outvector));
}

// 2D, to be removed use wavelet_decompose3Py
void wavelet_transform2DPy(uintptr_t invector,
                           int xlength,
                           int ylength,
                           int Filterlength, 
                           int levels,
                           uintptr_t waveletcoefficients) {
    wavelet_transform2D(reinterpret_cast<FLOAT*>(invector),
                        xlength,
                        ylength,
                        Filterlength,
                        levels,
                        reinterpret_cast<FLOAT*>(waveletcoefficients));
}
//to be removed, use wavelet_reconstruct3Py
void invwavelet_transform2DPy(uintptr_t waveletcoefficients,
                              int xlength,
                              int ylength,
                              int Filterlength, 
                              int levels,
                              uintptr_t outvector) {
    invwavelet_transform2D(reinterpret_cast<FLOAT*>(waveletcoefficients),
                           xlength,
                           ylength,
                           Filterlength,
                           levels,
                           reinterpret_cast<FLOAT*>(outvector));
}

// 3D, to be removed use wavelet_decompose3Py
void wavelet_transform3DPy(uintptr_t invector,
                           int xlength,
                           int ylength,
                           int zlength,
                           int Filterlength, 
                           int levels,
                           uintptr_t waveletcoefficients) {
    wavelet_transform3D(reinterpret_cast<FLOAT*>(invector),
                        xlength,
                        ylength,
                        zlength,
                        Filterlength,
                        levels,
                        reinterpret_cast<FLOAT*>(waveletcoefficients));
}
//to be removed, use wavelet_reconstruct3Py
void invwavelet_transform3DPy(uintptr_t waveletcoefficients,
                              int xlength,
                              int ylength,
                              int zlength,
                              int Filterlength, 
                              int levels,
                              uintptr_t outvector) {
    invwavelet_transform3D(reinterpret_cast<FLOAT*>(waveletcoefficients),
                           xlength,
                           ylength,
                           zlength,
                           Filterlength,
                           levels,
                           reinterpret_cast<FLOAT*>(outvector));
}

void adjointinvwavelet_transform3DPy(uintptr_t invector,
                                     int xlength,
                                     int ylength,
                                     int zlength,
                                     int Filterlength, 
                                     int levels,
                                     uintptr_t waveletcoefficients) {
  adjointinvwavelet_transform3D(reinterpret_cast<FLOAT*>(invector),
                                xlength,
                                ylength,
                                zlength,
                                Filterlength,
                                levels,
                                reinterpret_cast<FLOAT*>(waveletcoefficients));
}

void adjointwavelet_transform3DPy(uintptr_t waveletcoefficients,
                                  int xlength,
                                  int ylength,
                                  int zlength,
                                  int Filterlength, 
                                  int levels,
                                  uintptr_t outvector) {
  adjointwavelet_transform3D(reinterpret_cast<FLOAT*>( waveletcoefficients),
                             xlength,
                             ylength,
                             zlength,
                             Filterlength,
                             levels,
                             reinterpret_cast<FLOAT*>(outvector));
}

void adjointinvwavelet_transform2DPy(uintptr_t invector,
                                     int xlength,
                                     int ylength,
                                     int Filterlength, 
                                     int levels,
                                     uintptr_t waveletcoefficients) {
  adjointinvwavelet_transform2D(reinterpret_cast<FLOAT*>(invector),
                                xlength,
                                ylength,
                                Filterlength,
                                levels,
                                reinterpret_cast<FLOAT*>(waveletcoefficients));
}

void adjointwavelet_transform2DPy(uintptr_t waveletcoefficients,
                                  int xlength,
                                  int ylength,
                                  int Filterlength, 
                                  int levels,
                                  uintptr_t outvector) {
  adjointwavelet_transform2D(reinterpret_cast<FLOAT*>(waveletcoefficients),
                             xlength,
                             ylength,
                             Filterlength,
                             levels,
                             reinterpret_cast<FLOAT*>(outvector));
}

void adjointinvwavelet_transform1DPy(uintptr_t invector,
                                     int xlength,
                                     int Filterlength, 
                                     int levels,
                                     uintptr_t waveletcoefficients) {
  adjointinvwavelet_transform1D(reinterpret_cast<FLOAT*>(invector),
                                xlength,
                                Filterlength,
                                levels,
                                reinterpret_cast<FLOAT*>(waveletcoefficients));
}

void adjointwavelet_transform1DPy(uintptr_t waveletcoefficients,
                                  int xlength,
                                  int Filterlength, 
                                  int levels,
                                  uintptr_t outvector) {
  adjointwavelet_transform1D(reinterpret_cast<FLOAT*>(waveletcoefficients),
                             xlength,
                             Filterlength,
                             levels,
                             reinterpret_cast<FLOAT*>(outvector));
}


// Main method
// See wavelet_dec3.c
int wavelet_decompose3Py(uintptr_t inspacevector,
                         int Xlength,
                         int Ylength,
                         int Zlength,
                         int Filterlength,
                         int Levels,
                         int Skip,
                         uintptr_t covector,
                         int ifnotSilent) {
    int colength;

    wavelet_decompose3(reinterpret_cast<FLOAT*>(inspacevector),
                       Xlength,
                       Ylength,
                       Zlength,
                       Filterlength,
                       Levels,
                       Skip,
                       reinterpret_cast<FLOAT*>(covector),
                       &colength,
                       ifnotSilent);
}

// See wavelet_rec3.c
void wavelet_reconstruct3Py(uintptr_t reccovector,
                            int colength,
                            uintptr_t outvector,
                            int Xlength,
                            int Ylength,
                            int Zlength,
                            int Filterlength,
                            int Levels,
                            int Skip,
                            int ifnotSilent) {
   wavelet_reconstruct3(reinterpret_cast<FLOAT*>(reccovector),
                        colength,
                        reinterpret_cast<FLOAT*>(outvector),
                        Xlength,
                        Ylength,
                        Zlength,
                        Filterlength,
                        Levels,
                        Skip,
                        ifnotSilent);
}

// see  wavelet_invadj3.c
void wavelet_invadjoint3Py(uintptr_t inspacevector,
                           int Xlength,
                           int Ylength,
                           int Zlength,
                           int Filterlength,
                           int Levels,
                           int Skip,
                           uintptr_t covector,
                           int ifnotSilent) {
    int colength;

    wavelet_invadjoint3(reinterpret_cast<FLOAT*>(inspacevector),
                        Xlength,
                        Ylength,
                        Zlength,
                        Filterlength,
                        Levels,
                        Skip,
                        reinterpret_cast<FLOAT*>(covector),
                        &colength,
                        ifnotSilent);
}

// see  wavelet_adj3.c
void wavelet_adjoint3Py(uintptr_t reccovector,
                        int colength,
                        uintptr_t outvector,
                        int Xlength,
                        int Ylength,
                        int Zlength,
                        int Filterlength,
                        int Levels,
                        int Skip,
                        int ifnotSilent) {
   wavelet_adjoint3(reinterpret_cast<FLOAT*>(reccovector),
                    colength,
                    reinterpret_cast<FLOAT*>(outvector),
                    Xlength,
                    Ylength,
                    Zlength,
                    Filterlength,
                    Levels,
                    Skip,
                    ifnotSilent);
}


// Expose classes and methods to Python
BOOST_PYTHON_MODULE(jos_waveletspy) {
    // 1D
    def("wavelet_transform1D", wavelet_transform1DPy);
    def("invwavelet_transform1D", invwavelet_transform1DPy);
    def("adjointinvwavelet_transform1D", adjointinvwavelet_transform1DPy);
    def("adjointwavelet_transform1D", adjointwavelet_transform1DPy);
    // 2D
    def("wavelet_transform2D", wavelet_transform2DPy);
    def("invwavelet_transform2D", invwavelet_transform2DPy);
    def("adjointinvwavelet_transform2D", adjointinvwavelet_transform2DPy);
    def("adjointwavelet_transform2D", adjointwavelet_transform2DPy);
    // 3D
    def("wavelet_transform3D", wavelet_transform3DPy);
    def("invwavelet_transform3D", invwavelet_transform3DPy);
    def("adjointinvwavelet_transform3D", adjointinvwavelet_transform3DPy);
    def("adjointwavelet_transform3D", adjointwavelet_transform3DPy);
    // Generic
    def("wavelet_decompose3", wavelet_decompose3Py);
    def("wavelet_reconstruct3", wavelet_reconstruct3Py);
    def("wavelet_invadjoint3", wavelet_invadjoint3Py);
    def("wavelet_adjoint3", wavelet_adjoint3Py);
}
