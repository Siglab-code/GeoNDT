========
Overview
========
 
GeoNDT is a fast general-purpose computational tool for geotechnical non-destructive testing applications.  
GeoNDT is flexible, general-purpose, and can be used seamlessly for advanced signal interpretation in geophysical 
laboratory testing including the bender element (BE) and ultrasonic pulse velocity (UPV) tests, characterization of 
complex multiphase geomaterials, in-situ shallow seismic geophysics including the falling weight deflectometer (FWD) 
and multichannel analysis of surface waves (MASW) tests. The advanced physics-based signal interpretation feature of 
GeoNDT allows the quantitative characterization of geophysical and geomechanical properties of geomaterials and multilayered 
geosystems independently without making any simplified assumptions as common in the current practice.


Quick start
===========

Install (only for Linux)::

    pip install geondt

To install development version, clone this repo and install in Linux::

    git clone https://github.com/siglab/geondt

    cd geondt

    pip install -e .


To install development version, clone this repo and install in Windows::


    git clone https://github.com/siglab/geondt

    cd geondt

    python setup.py build --compiler=mingw32 

    python setup.py install  

Usage
-----

    The GeoNDT can efficiently study the three-dimensional wave propagation within soil specimens in the BE test. Sample code is given as follows: 

    >>> import numpy as np 
    >>> from geondt import one_phase_dynamic  
    >>> import  json 
    >>> with open('BE_dry.json', "r") as f:
            data = json.load(f)   
    >>> BE = one_phase_dynamic(**data["input"])   
    >>> signal = BE.run_f()  

    
 
Troubleshooting
===============

The installation procedure assumes that the Fortran compiler such as Gfortran and Lapack library are installed on your system.
To install Gfortran and Lapack in Linux::

    sudo apt install gfortran
    sudo apt-get install liblapacke-dev checkinstall 
    export gfortran="/home/kay/gcc-4.8.5/bin/gfortran"

To install Gfortran and Lapack in Windows::

* Use MinGW <https://sourceforge.net/projects/mingw-w64/> to get Gfortran. Make sure the Mingw is added to the system path. 
* Then add the liblapack.a file (can be found under lib folder in this respiratory ) in the MinGW folder (C:\mingw64\x86_64-w64-mingw32\lib). 
 

References
==========

.. [1] Liu H, Maghoul P, Mantelet, G, Shalaby A
       GeoNDT: a fast general-purpose computational tool for geotechnical non-destructive testing applications. Computers and Geotechnics.

.. [2] Liu H, Maghoul P, Shalaby A, Bahari A, Moradi F. 
       Integrated approach for the MASW dispersion analysis using the spectral element technique and trust region reflective method. 
       Computers and Geotechnics. 2020 Sep 1;125:103689.
