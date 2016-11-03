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

# Copyright 2014, 2015 The ODL development group
#
# This file is part of ODL.
#
# ODL is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# ODL is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with ODL.  If not, see <http://www.gnu.org/licenses/>.


# Imports for common Python 2/3 codebase
from __future__ import print_function, division, absolute_import
from future import standard_library
standard_library.install_aliases()

import pytest
import numpy as np
import odl
from odl.util.testutils import all_almost_equal

from jos_waveletspy import (BiorthWaveletTransform,
                            InverseAdjBiorthWaveletTransform)

# Simply modify exp_params to modify the fixture
wavelet_params = ['josbiorth1', 'josbiorth3', 'josbiorth5',
                  'josbiorth7', 'josbiorth9']
wavelet_ids = [' wavelet = {} '.format(p) for p in wavelet_params]


@pytest.fixture(scope="module", ids=wavelet_ids, params=wavelet_params)
def wbasis(request):
    return request.param


def test_bwt1d(wbasis):
    # Verify that the operator works as axpected
    # 1D test
    n = 16
    x = np.zeros(n)
    x[5:10] = 1
    nscales = 2

    # Define a discretized domain
    domain = odl.FunctionSpace(odl.Interval([-1], [1]))
    nPoints = np.array([n])
    disc_domain = odl.uniform_discr_fromspace(domain, nPoints)
    disc_phantom = disc_domain.element(x)

    # Create the discrete wavelet transform operator.
    # Only the domain of the operator needs to be defined
    Wop = BiorthWaveletTransform(disc_domain, nscales, wbasis)

    # Compute the discrete wavelet transform of discrete imput image
    coeffs = Wop(disc_phantom)

    # Compute the inverse wavelet transform
    reconstruction = Wop.inverse(coeffs)

    # Verify that reconstructions lie in correct discretized domain
    assert reconstruction in disc_domain
    assert all_almost_equal(reconstruction.asarray(), x)


def test_bwt2d():
    # 2D test
    n = 16
    x = np.zeros((n, n))
    x[5:10, 5:10] = 1
    wbasis = 'josbiorth5'
    nscales = 3

    # Define a discretized domain
    domain = odl.FunctionSpace(odl.Rectangle([-1, -1], [1, 1]))
    nPoints = np.array([n, n])
    disc_domain = odl.uniform_discr_fromspace(domain, nPoints)
    disc_phantom = disc_domain.element(x)

    # Create the discrete wavelet transform operator.
    # Only the domain of the operator needs to be defined
    Bop = BiorthWaveletTransform(disc_domain, nscales, wbasis)
    Bop2 = InverseAdjBiorthWaveletTransform(disc_domain, nscales, wbasis)

    # Compute the discrete wavelet transform of discrete imput image
    coeffs = Bop(disc_phantom)
    coeffs2 = Bop2(disc_phantom)

    reconstruction = Bop.inverse(coeffs)
    reconstruction2 = Bop2.inverse(coeffs2)

    assert all_almost_equal(reconstruction.asarray(), x)
    assert all_almost_equal(reconstruction2.asarray(), x)


def test_bwt3d():
    # 3D test
    n = 16
    x = np.zeros((n, n, n))
    x[5:10, 5:10, 5:10] = 1
    wbasis = 'josbiorth7'
    nscales = 1
    # Define a discretized domain
    domain = odl.FunctionSpace(odl.Cuboid([-1, -1, -1], [1, 1, 1]))
    nPoints = np.array([n, n, n])
    disc_domain = odl.uniform_discr_fromspace(domain, nPoints)
    disc_phantom = disc_domain.element(x)

    # Create the discrete wavelet transform operator.
    # Only the domain of the operator needs to be defined
    Bop = BiorthWaveletTransform(disc_domain, nscales, wbasis)
    Bop2 = InverseAdjBiorthWaveletTransform(disc_domain, nscales, wbasis)

    # Compute the discrete wavelet transform of discrete imput image
    coeffs = Bop(disc_phantom)
    coeffs2 = Bop2(disc_phantom)

    reconstruction = Bop.inverse(coeffs)
    reconstruction2 = Bop2.inverse(coeffs2)

    assert all_almost_equal(reconstruction.asarray(), x)
    assert all_almost_equal(reconstruction2.asarray(), x)


if __name__ == '__main__':
    pytest.main([str(__file__.replace('\\', '/')), ' -v'])
