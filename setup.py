#!/usr/bin/env python
import os 
import sys
import setuptools  # necessary even when numpy.distutils.core.setup is used
from setuptools import dist, setup 
 
dist.Distribution().fetch_build_eggs(['numpy>=1.10'])  

 


# Include extensions only when not on readthedocs.org
if os.environ.get('READTHEDOCS', None) == 'True':
    from distutils.core import Extension
    PoroSEM = []
else:
    from numpy.distutils.core import setup, Extension 

    if sys.platform.startswith('linux'): 
        PoroSEM = [
            Extension(name='geondt.PoroSEM',
                    sources=['src/fortrangeondt/PoroSEM.f90'], 
                    extra_link_args=["-llapack"])] 
                    
    elif sys.platform.startswith('win32'):  
        PoroSEM = [
        Extension(name='geondt.PoroSEM',
                  sources=['src/fortrangeondt/PoroSEM.f90'], 
                  extra_link_args=["-lgfortran"], 
                  libraries = ["lapack"], 
                  include_dirs = ["C:\mingw64\include"],
                  library_dirs = ["C:\mingw64\lib"])]

requirements = ['numpy >=1.16.2', 
                'scipy >=1.2.1', 
                'joblib>=0.17.0', 
                'matplotlib >= 3.0.3'
] 

setup_requirements = [ ]

test_requirements = [ ] 

setup(ext_modules=PoroSEM, package_data={'geondt': ['iltcme.json']}, install_requires=requirements, include_package_data=True)
