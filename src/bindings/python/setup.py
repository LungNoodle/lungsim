""" Aether library - routines for a modelling the lung

The Aether library is an advanced modelling library for models of the lung.
"""

classifiers = """\
Development Status :: 4 - Beta
Intended Audience :: Developers
Intended Audience :: Education
Intended Audience :: Science/Research
License :: OSI Approved :: Apache Software License
Programming Language :: Python
Programming Language :: Python :: @PYTHONLIBS_MAJOR_VERSION@.@PYTHONLIBS_MINOR_VERSION@
@SETUP_PY_OPERATING_SYSTEM_CLASSIFIER@
Topic :: Scientific/Engineering :: Medical Science Apps.
Topic :: Software Development :: Libraries :: Python Modules
"""

import sys
from setuptools import setup
from setuptools.dist import Distribution


class BinaryDistribution(Distribution):
    def is_pure(self):
        return False

    def has_ext_modules(self):
        return True


doclines = __doc__#.split("\n")

PLATFORM_PACKAGE_DATA = ["*.so", "*.pyd", ]
if sys.platform.startswith('win32'):
    PLATFORM_PACKAGE_DATA.extend(["aether.dll", "aether_c.dll"])

setup(
    name='lungnoodle.aether',
    version='@Aether_VERSION@',
    author='Lung Group, Auckland Bioengineering Institute.',
    author_email='h.sorby@auckland.ac.nz',
    packages=['aether'],
    package_data={'aether': PLATFORM_PACKAGE_DATA},
    url='https://lung.bioeng.auckland.ac.nz/',
    license='http://www.apache.org/licenses/LICENSE-2.0',
    description='Aether library of routines for modelling the lung.',
    classifiers = filter(None, classifiers.split("\n")),
    long_description=doclines,
    distclass=BinaryDistribution,
    include_package_data=True,
)
