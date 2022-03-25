#!/usr/bin/env python

from setuptools import setup

MAINTAINER_EMAIL = 'qr'
AUTHORS = 'qr'

setup(
   name='MNHPy',
   version='0.0',
   description=('Python visulization tools for MesoNH'),
   license='GPL3',
   author=AUTHORS,
   author_email=MAINTAINER_EMAIL,
   url='https://github.com/QuentinRodier/MNHPy',
   packages=['MNHPy'],	#name of package (dir with the modules)
   # external dependencies packages
   install_requires=['numpy', 'matplotlib', 'cartopy'],
   )
