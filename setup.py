#!/usr/bin/env python

from setuptools import setup

MAINTAINER_EMAIL = 'mesonhsupport@obs-mip.fr'
AUTHORS = 'Quentin Rodier & Meso-NH Team'

setup(
   name='MNHPy',
   version='0.2.0',
   description=('Python visulization tools for MesoNH atmospheric research model'),
   long_description=('Compatible with Meso-NH model version 5.5.1'),
   author=AUTHORS,
   author_email=MAINTAINER_EMAIL,
   url='https://github.com/QuentinRodier/MNHPy',
   classifiers=[
       "Programming Language :: Python :: 3",
       "License :: CeCILL-C Free Software License Agreement (CECILL-C)"],
   packages=['MNHPy'],	#name of package (dir with the modules)
   # external dependencies packages
   install_requires=['numpy>=1.17.2', 'matplotlib>=3.1.1', 'cartopy>=0.17', 'scipy>=1.3.1', 'netCDF4>=1.4.2'],
   python_requires=">=3.6",
   package_dir={"": "src"},
   )
