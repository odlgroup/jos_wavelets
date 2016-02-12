# -*- coding: utf-8 -*-
"""

"""

from future import standard_library
standard_library.install_aliases()

from distutils.core import setup


setup(name='jos_wavelets',
      # version=__version__,
      version='0.1.1',
      author='Jan-Olov Stromberg, Jonas Adler, Julian Moosmann',
      author_email='odl@math.kth.se',
      url='//https://github.com/odlgroup/jos_wavelets.git',
      description='Python bindings for the jos_wavelets library',
      license='GPLv3',
      packages=['jos_wavelets'],
      package_dir={'jos_wavelets': '.'},
      package_data={'jos_wavelets': ['*.*']}, requires=['scipy', 'numpy',
                                                        'matplotlib', 'numba',
                                                        'future'])
