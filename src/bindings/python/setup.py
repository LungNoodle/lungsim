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
Programming Language :: Python :: 2.7
Programming Language :: Python :: 3.5
Operating System :: Microsoft :: Windows
Operating System :: Unix
Operating System :: MacOS :: MacOS X
Topic :: Scientific/Engineering :: Medical Science Apps.
Topic :: Software Development :: Libraries :: Python Modules
"""

from setuptools import setup

doclines = __doc__#.split("\n")

setup(
	name='aether',
	version='0.1.0',
	author='Lung Group, Auckland Bioengineering Institute.',
	author_email='h.sorby@auckland.ac.nz',
	packages=['aether'],
#	package_data={'aether': []},
	platforms=['any'],
	url='https://lung.bioeng.auckland.ac.nz/',
	license='http://www.apache.org/licenses/LICENSE-2.0',
	description='Aether library of routines for modelling the lung.',
	classifiers = filter(None, classifiers.split("\n")),
	long_description=doclines,
)
