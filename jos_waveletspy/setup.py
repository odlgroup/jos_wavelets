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

"""Setup file for jos_wavelets."""

from future import standard_library
standard_library.install_aliases()

from distutils.core import setup


setup(name='jos_wavelets',
      # version=__version__,
      version='0.1.1',
      author='Jan-Olov Stromberg, Jonas Adler, Julian Moosmann',
      author_email='odl@math.kth.se',
      url='//https://github.com/odlgroup/jos_wavelets.git',
      description='ODL plugin for the jos_wavelets library',
      license='GPLv3',
      packages=['jos_wavelets'],
      package_dir={'jos_wavelets': '.'},
      package_data={'jos_wavelets': ['*.*']},
      requires=['future', 'odl'])
