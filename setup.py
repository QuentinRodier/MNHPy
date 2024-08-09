#!/usr/bin/env python

from setuptools import setup

MAINTAINER_EMAIL = 'mesonhsupport@obs-mip.fr'
AUTHORS = 'Quentin Rodier & Meso-NH Team'

setup(
   name='MNHPy',
   version='0.3.2',
   description=('Python visulization tools for MesoNH atmospheric research model'),
   long_description=('Compatible with Meso-NH model version >=5.7.1'),
   author=AUTHORS,
   author_email=MAINTAINER_EMAIL,
   url='https://github.com/QuentinRodier/MNHPy',
   classifiers=[
       "Programming Language :: Python :: 3",
       "License :: CeCILL-C Free Software License Agreement (CECILL-C)"],
   packages=['MNHPy'],	#name of package (dir with the modules)
   # external dependencies packages
   install_requires=['numpy>=1.26.4', 'matplotlib==3.9.1', 'cartopy>=0.21.1', 'scipy>=1.14.0', 'netCDF4>=1.7.1'],
   python_requires=">=3.12.3",
   package_dir={"": "src"},
   )
