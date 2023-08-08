#!/usr/bin/env python

from setuptools import setup

MAINTAINER_EMAIL = 'mesonhsupport@obs-mip.fr'
AUTHORS = 'Quentin Rodier & Meso-NH Team'

setup(
   name='MNHPy',
   version='0.3.0',
   description=('Python visulization tools for MesoNH atmospheric research model'),
   long_description=('Compatible with Meso-NH model version 5.6.0+'),
   author=AUTHORS,
   author_email=MAINTAINER_EMAIL,
   url='https://github.com/QuentinRodier/MNHPy',
   classifiers=[
       "Programming Language :: Python :: 3",
       "License :: CeCILL-C Free Software License Agreement (CECILL-C)"],
   packages=['MNHPy'],	#name of package (dir with the modules)
   # external dependencies packages
   install_requires=['numpy>=1.21.5', 'matplotlib>=3.7.0', 'cartopy>=0.21.1', 'scipy>=1.10.1', 'netCDF4>=1.5.3'],
   python_requires=">=3.8.10",
   package_dir={"": "src"},
   )
