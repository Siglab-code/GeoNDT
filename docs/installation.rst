.. _installation:

Installation
============

This package can be installed using either PyPI or GitHub: 

PyPI
----
This is the most straightforward option!  From the command line in Linux system use
``pip`` [1]_::

    pip install geondt

The PyPI method currently supports for the installation in Linux system. For Windows users, 
GeoNDT can be installed from the source code. The package has been tested with the following setups 
(others might work, too):

* Windows (64 bit Python) and Linux (64 bit) 
* Python 3.7, 3.8, 3.9


GitHub
------
If you are intending on modifying geondt modules, it's easier to
install geondt by forking the repository and installing it locally or within
a virtual environment. After clonining the repo, you may install in linux by:: 

    git clone https://github.com/siglab/geondt

    cd geondt
    
    pip install -e .

To install development version, clone this repo and install in Windows::


    git clone https://github.com/siglab/geondt

    cd geondt

    python setup.py build --compiler=mingw32 

    python setup.py install   

The install procedure assumes that the Fortran compiler such as Gfortran and Lapack library are installed on your system.
To install Gfortran and Lapack in Linux::

    sudo apt install gfortran
    sudo apt-get install liblapacke-dev checkinstall 
    export gfortran="/home/kay/gcc-4.8.5/bin/gfortran"

To install Gfortran and Lapack in Windows::

* Use MinGW <https://sourceforge.net/projects/mingw-w64/> to get Gfortran.
* Then add the liblapack.a file (can be found under lib folder in this respiratory ) in the MinGW folder (C:\mingw64\x86_64-w64-mingw32\lib). 
 
