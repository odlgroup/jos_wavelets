# Copyright 2015-2016 Jan-Olov Stromberg <jostromb@kth.se>
#
# This file is part of jos_wavelets.
#
# jos_wavelets is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# jos_wavelets is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with jos_wavelets.  If not, see <http://www.gnu.org/licenses/>.

# Imports for common Python 2/3 codebase
from __future__ import print_function, division, absolute_import
from future import standard_library
standard_library.install_aliases()
from builtins import super

import numpy as np
import odl

# TODO: this is untested, check imports
from jos_wavelets import (
    wavelet_transform1D, wavelet_transform2D, wavelet_transform3D,
    adjointwavelet_transform1D, adjointwavelet_transform2D,
    adjointwavelet_transform3D,
    invwavelet_transform1D, invwavelet_transform2D, invwavelet_transform3D,
    adjointinvwavelet_transform1D, adjointinvwavelet_transform2D,
    adjointinvwavelet_transform3D)


class BiorthWaveletTransform(odl.Operator):
    """Discrete biorthogonal wavelet trafo between discrete L2 spaces.

    This class implements biorthogonal wavelet transforms using
    J.-O- Stromberg's biorthogonal wavelet library. J.-O- Stromberg's
    biorthogonal wavelets have odd length filters, the available filter
    lengths are 1, 3, 5, 7 and 9.

    In J.-O. Strombergs wavelet library (and in what follows)
    the following transforms have been implemented:

    - biorthogonal wavelet transforms
    - adjoint biorthogonal wavelet transforms
    - inverse biorthogonal wavelet transforms
    - adjoint inverse biorthogonal wavelet transforms
    """

    def __init__(self, dom, nscales, wbasis):
        """Initialize a new instance.

        Parameters
        ----------
        dom : `DiscreteLp`
            Domain of the biorthogonal wavelet transform (the "image domain").
            The exponent :math:`p` of the discrete :math:`L^p`
            space must be equal to 2.0.
        nscales : `int`
            Number of scaling levels in the wavelet transform.
        wbasis : `str`
            (``'josbiorthN'``) J.O. Stromberg's biorthogonal wavelets
            with ``N`` referring to the filter length.
            Possible values are:

            1, 3, 5, 7, 9

        Signal extention mode: reflexion
        """

        self.nscales = int(nscales)
        self.wbasis = wbasis

        if not isinstance(dom, odl.DiscreteLp):
            raise TypeError('domain {!r} is not a `DiscreteLp` instance.'
                            ''.format(dom))

        if dom.exponent != 2.0:
            raise ValueError('domain Lp exponent is {} instead of 2.0.'
                             ''.format(dom.exponent))

        max_level = int(np.ceil(np.log2(np.max(dom.shape))))

        if nscales > max_level:
            raise ValueError('Maximum useful number of scaling levels is {}, '
                             'got {}.'.format(max_level, self.nscales))
            # raiseETC ('Setting {} to {}' .format(self.nscales, max_level))

        self.filterlength = int(wbasis[-1])
        if self.filterlength not in (1, 3, 5, 7, 9):
            raise NotImplementedError('Filterlength {} not 1, 3, 5, 7 or 9'
                                      ''.format(self.filterlength))

        ran = dom.dspace_type(dom.size, dtype=dom.dtype)
        super().__init__(dom, ran, linear=True)

    def _call(self, x):
        """Compute the discrete biorthogonal wavelet transform.

        Parameters
        ----------
        x : `DiscreteLpVector`

        Returns
        -------
        arr : `numpy.ndarray`
            Flattened and concatenated coefficient array
            The length of the array depends on the size of input image to
            be transformed and on the chosen wavelet basis.
        """
        if x.space.ndim == 1:
            # A copy of x is unavoidable since bwt writes over the given inputs
            x_cpy = x.copy()
            x_cpy = x_cpy.asarray()
            coeff = self.range.element().asarray()
            wavelet_transform1D(x_cpy.ctypes.data, x_cpy.shape[0],
                                self.filterlength,
                                self.nscales,
                                coeff.ctypes.data)
            return coeff

        elif x.space.ndim == 2:
            x_cpy = x.copy()
            x_cpy = x_cpy.asarray()
            coeff = self.range.element().asarray()
            wavelet_transform2D(x_cpy.ctypes.data, x_cpy.shape[0],
                                x_cpy.shape[1], self.filterlength,
                                self.nscales, coeff.ctypes.data)
            return self.range.element(coeff)

        elif x.space.ndim == 3:
            x_cpy = x.copy()
            x_cpy = x_cpy.asarray()
            coeff = self.range.element().asarray()
            wavelet_transform3D(x_cpy.ctypes.data, x_cpy.shape[0],
                                x_cpy.shape[1], x_cpy.shape[2],
                                self.filterlength, self.nscales,
                                coeff.ctypes.data)
            return self.range.element(coeff)

    @property
    def adjoint(self):
        """The biorthogonal adjoint wavelet transform."""
        return AdjBiorthWaveletTransform(ran=self.domain,
                                         nscales=self.nscales,
                                         wbasis=self.wbasis)

    @property
    def inverse(self):
        """The biorthogonal inverse wavelet transform."""
        return InverseBiorthWaveletTransform(ran=self.domain,
                                             nscales=self.nscales,
                                             wbasis=self.wbasis)

    @property
    def adjointinverse(self):
        """The inverse of the adjoint of biorthogonal wavelet transform."""
        return InverseAdjBiorthWaveletTransform(dom=self.domain,
                                                nscales=self.nscales,
                                                wbasis=self.wbasis)


class AdjBiorthWaveletTransform(odl.Operator):
    """Discrete adjoint of biorthogonal wavelet transform between L2 space.

    This class implements adjoint biorthogonal wavelet transforms using
    J.-O- Stromberg's biorthogonal wavelet library. J.-O- Stromberg's
    biorthogonal wavelets have odd length filters, the available filter
    lengths are 1, 3, 5, 7 and 9.

    In J.-O. Strombergs wavelet library (and in what follows)
    the following transforms have been implemented:

    - biorthogonal wavelet transforms
    - adjoint biorthogonal wavelet transforms
    - inverse biorthogonal wavelet transforms
    - adjoint inverse biorthogonal wavelet transforms
    """

    def __init__(self, ran, nscales, wbasis):
        """Initialize a new instance.

        Parameters
        ----------
        ran : `DiscreteLpVector` ("image domain")

        nscales : `int`
            Number of scaling levels in the wavelet transform.

        wbasis : `str`
            (``'josbiorthN'``) J.-O. Stromberg's biorthogonal wavelets
            where ``N`` refers to the filter length.
            Possible values are:

            1, 3, 5, 7, 9

        Signal extension mode: reflexion
        """

        self.nscales = int(nscales)
        self.wbasis = wbasis

        if not isinstance(ran, odl.DiscreteLp):
            raise TypeError('domain {!r} is not a `DiscreteLp` instance.'
                            ''.format(ran))

        if ran.exponent != 2.0:
            raise ValueError('domain Lp exponent is {} instead of 2.0.'
                             ''.format(ran.exponent))
        max_level = int(np.ceil(np.log2(np.max(ran.shape))))
        if nscales > max_level:
            raise ValueError('Maximum useful number of scaling levels is {}, '
                             'got {}.'.format(max_level, self.nscales))
            # raiseETC ('Setting {} to {}' .format(self.nscales, max_level))

        self.filterlength = int(wbasis[-1])
        if self.filterlength not in (1, 3, 5, 7, 9):
            raise NotImplementedError('Filterlength {} not 1, 3, 5, 7 or 9'
                                      ''.format(self.filterlength))

        dom = ran.dspace_type(ran.size, dtype=ran.dtype)
        super().__init__(dom, ran, linear=True)

    def _call(self, coeff):
        """Discrete adjoint wavelet transform with biorthogonal wavelets.

        Parameters
        ----------
        coeff : `DiscreteLpVector`

        Returns
        -------
        arr : `DiscreteLpVector`
        """
        if len(self.range.grid.shape) == 1:
            nx = self.range.grid.shape[0]
            # A copy of coeff is unavoidable since bwt writes over
            # the given input
            coeff_in = coeff.copy().asarray()
            x = self.range.element().asarray()
            adjointwavelet_transform1D(coeff_in.ctypes.data, nx,
                                       self.filterlength,
                                       self.nscales, x.ctypes.data)
            return self.range.element(x)

        elif len(self.range.grid.shape) == 2:
            (nx, ny) = (self.range.grid.shape[0], self.range.grid.shape[1])
            coeff_in = coeff.copy().asarray()
            x = self.range.element().asarray()
            adjointwavelet_transform2D(coeff_in.ctypes.data, nx, ny,
                                       self.filterlength,
                                       self.nscales, x.ctypes.data)
            return self.range.element(x)

        elif len(self.range.grid.shape) == 3:
            (nx, ny, nz) = (self.range.grid.shape[0], self.range.grid.shape[1],
                            self.range.grid.shape[2])
            coeff_in = coeff.copy().asarray()
            x = self.range.element().asarray()
            adjointwavelet_transform3D(coeff_in.ctypes.data, nx, ny, nz,
                                       self.filterlength, self.nscales,
                                       x.ctypes.data)
            return self.range.element(x)

    @property
    def adjoint(self):
        """The biorthogonal wavelet transform."""
        return BiorthWaveletTransform(dom=self.range, nscales=self.nscales,
                                      wbasis=self.wbasis)

    @property
    def inverse(self):
        """The inverse of the adjoint of biorthogonal wavelet transform."""
        return InverseAdjBiorthWaveletTransform(ran=self.range,
                                                nscales=self.nscales,
                                                wbasis=self.wbasis)


class InverseBiorthWaveletTransform(odl.Operator):
    """Discrete inverse of biorthogonal wavelet transform between L2 space

    This class implements inverse biorthogonal wavelet transforms using
    J.-O- Stromberg's biorthogonal wavelet library. J.-O- Stromberg's
    biorthogonal wavelets have odd length filters, the available filter
    lengths are 1, 3, 5, 7 and 9.

    In J.-O. Strombergs wavelet library (and in what follows)
    the following transforms have been implemented:

    - biorthogonal wavelet transforms
    - adjoint biorthogonal wavelet transforms
    - inverse biorthogonal wavelet transforms
    - adjoint inverse biorthogonal wavelet transforms
    """

    def __init__(self, ran, nscales, wbasis):
        """Initialize a new instance.

        Parameters
        ----------
        ran : `DiscreteLpVector` ("image domain")

        nscales : `int`
            Number of scaling levels in the wavelet transform.
            The maximum number of usable scales with J.-O.S.
            biorthogonal wavelets can be determined using xxx
        wbasis : `str`
            (``'josbiorthN'``) J.-O. Stromberg's biorthogonal wavelets
            ``N`` is the filter length.
            Possible values are:

            1, 3, 5, 7 or 9

        Signal extention mode: reflection
        """
        self.nscales = int(nscales)
        self.wbasis = wbasis

        if not isinstance(ran, odl.DiscreteLp):
            raise TypeError('domain {!r} is not a `DiscreteLp` instance.'
                            ''.format(ran))

        if ran.exponent != 2.0:
            raise ValueError('domain Lp exponent is {} instead of 2.0.'
                             ''.format(ran.exponent))
        max_level = int(np.ceil(np.log2(np.max(ran.shape))))

        if nscales > max_level:
            raise ValueError('Maximum useful number of scaling levels is {}, '
                             'got {}.'.format(max_level, self.nscales))
            # raiseETC ('Setting {} to {}' .format(self.nscales, max_level))
            # nscales = max_level

        self.filterlength = int(wbasis[-1])
        if self.filterlength not in (1, 3, 5, 7, 9):
            raise NotImplementedError('Filterlength {} not 1, 3, 5, 7 or 9'
                                      ''.format(self.filterlength))

        dom = ran.dspace_type(ran.size, dtype=ran.dtype)
        super().__init__(dom, ran, linear=True)

    def _call(self, coeff):
        """Discrete inverse wavelet transform with biorthogonal wavelets.

        Parameters
        ----------
        coeff : `DiscreteLpVector`

        Returns
        -------
        arr : `DiscreteLpVector`
        """
        if len(self.range.grid.shape) == 1:
            nx = self.range.grid.shape[0]
            coeff_in = coeff.copy().asarray()
            x = self.range.element().asarray()
            invwavelet_transform1D(coeff_in.ctypes.data, nx,
                                   self.filterlength,
                                   self.nscales, x.ctypes.data)
            return self.range.element(x)

        elif len(self.range.grid.shape) == 2:
            (nx, ny) = (self.range.grid.shape[0], self.range.grid.shape[1])
            coeff_in = coeff.copy().asarray()
            x = self.range.element().asarray()
            invwavelet_transform2D(coeff_in.ctypes.data, nx, ny,
                                   self.filterlength,
                                   self.nscales, x.ctypes.data)
            return self.range.element(x)

        elif len(self.range.grid.shape) == 3:
            (nx, ny, nz) = (self.range.grid.shape[0], self.range.grid.shape[1],
                            self.range.grid.shape[2])
            coeff_in = coeff.copy().asarray()
            x = self.range.element().asarray()
            invwavelet_transform3D(coeff_in.ctypes.data, nx, ny, nz,
                                   self.filterlength,
                                   self.nscales, x.ctypes.data)
            return self.range.element(x)

    @property
    def adjoint(self):
        """The biorthogonal wavelet transform."""
        return InverseAdjBiorthWaveletTransform(dom=self.range,
                                                nscales=self.nscales,
                                                wbasis=self.wbasis)

    @property
    def inverse(self):
        """The inverse of the adjoint of biorthogonal wavelet transform."""
        return BiorthWaveletTransform(dom=self.range, nscales=self.nscales,
                                      wbasis=self.wbasis)


class InverseAdjBiorthWaveletTransform(odl.Operator):
    """Discrete biorthogonal wavelet transfrom, inverse of adjoint.

    This class implements inverse adjoint biorthogonal wavelet transforms
    using J.-O- Stromberg's biorthogonal wavelet library. J.-O- Stromberg's
    biorthogonal wavelets have odd length filters, the available filter
    lengths are 1, 3, 5, 7 and 9.

    In J.-O. Strombergs wavelet library (and in what follows)
    the following transforms have been implemented:

    - biorthogonal wavelet transforms
    - adjoint biorthogonal wavelet transforms
    - inverse biorthogonal wavelet transforms
    - adjoint inverse biorthogonal wavelet transforms
    """
    def __init__(self, dom, nscales, wbasis):
        """Initialize a new instance.

        Parameters
        ----------
        dom : `DiscreteLp` ("image domain")
        nscales : `int`
            Number of scaling levels in the wavelet transform.
            The maximum number of usable scales with J.-O.S.
            biorthogonal wavelets can be determined using xxx
        wbasis : `str`
            (``'josbiorthN'``) J.-O. Stromberg's biorthogonal wavelets
            ``N`` is the filter length.
            Possible values are:

            1, 3, 5, 7 or 9

        signal extension mode: reflection
        """

        self.nscales = int(nscales)
        self.wbasis = wbasis

        if not isinstance(dom, odl.DiscreteLp):
            raise TypeError('domain {!r} is not a `DiscreteLp` instance.'
                            ''.format(dom))

        if dom.exponent != 2.0:
            raise ValueError('domain Lp exponent is {} instead of 2.0.'
                             ''.format(dom.exponent))

        max_level = int(np.ceil(np.log2(np.max(dom.shape))))
        if nscales > max_level:
            raise ValueError('Maximum useful number of scaling levels is {}, '
                             'got {}.'.format(max_level, self.nscales))
            # raiseETC ('Setting {} to {}' .format(self.nscales, max_level))

        self.filterlength = int(wbasis[-1])
        if self.filterlength not in (1, 3, 5, 7, 9):
            raise NotImplementedError('Filterlength {} not 1, 3, 5, 7 or 9'
                                      ''.format(self.filterlength))

        ran = dom.dspace_type(dom.size, dtype=dom.dtype)
        super().__init__(dom, ran, linear=True)

    def _call(self, x):
        """Discrete biorthogonal wavelet transform, inverse of adjoint.

        Parameters
        ----------
        x : `DiscreteLpVector`

        Returns
        -------
        arr : `numpy.ndarray`
            Flattened and concatenated coefficient array
        """
        if x.space.ndim == 1:
            x_cpy = x.copy()
            x_cpy = x_cpy.asarray()
            coeff = self.range.element().asarray()
            adjointinvwavelet_transform1D(x_cpy.ctypes.data,
                                          x_cpy.shape[0],
                                          self.filterlength,
                                          self.nscales, coeff.ctypes.data)
            return self.range.element(coeff)

        elif x.space.ndim == 2:
            x_cpy = x.copy()
            x_cpy = x_cpy.asarray()
            coeff = self.range.element().asarray()
            adjointinvwavelet_transform2D(x_cpy.ctypes.data,
                                          x_cpy.shape[0], x_cpy.shape[1],
                                          self.filterlength, self.nscales,
                                          coeff.ctypes.data)
            return self.range.element(coeff)

        elif x.space.ndim == 3:
            x_cpy = x.copy()
            x_cpy = x_cpy.asarray()
            coeff = self.range.element().asarray()
            adjointinvwavelet_transform3D(x_cpy.ctypes.data,
                                          x_cpy.shape[0], x_cpy.shape[1],
                                          x_cpy.shape[2],
                                          self.filterlength,
                                          self.nscales, coeff.ctypes.data)
            return self.range.element(coeff)

    @property
    def adjoint(self):
        """The biorthogonal adjoint wavelet transform."""
        return InverseBiorthWaveletTransform(ran=self.domain,
                                             nscales=self.nscales,
                                             wbasis=self.wbasis)

    @property
    def inverse(self):
        """The biorthogonal inverse wavelet transform."""
        return AdjBiorthWaveletTransform(ran=self.domain,
                                         nscales=self.nscales,
                                         wbasis=self.wbasis)
