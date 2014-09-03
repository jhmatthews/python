# Python

Python is a (confusingly named) Monte Carlo radiative transfer code which uses the Sobolev approximation.
It has been developed by Knox Long, Christian Knigge, Stuart Sim, Nick Higginbottom and James Matthews. 

The code is not yet ready for use and should not be installed yet. If you are interested in using
Python please contact long@stsci.edu.

# Travis Build Status

Simple build checks are carried out on Travis. The latest status shows up below:

[![Build Status](https://travis-ci.org/agnwinds/python.png?branch=dev)](https://travis-ci.org/agnwinds/python)


# Getting the radiative transfer code 'Python'

You can download the required structure under the structure branch. e.g.
git clone https://github.com/agnwinds/python.git -b structure
or simply click on the 'zip' button!


Releases of progs can be found under [tags](https://github.com/agnwinds/python/tags "Wiki").

Consult the [wiki](https://github.com/agnwinds/python/wiki/_pages "Wiki") for how to install Python.


# Installation

Python and the various routines associated are set up in a self-contained directory structure. The basic directory structure and the data files that one needs to run Python need to be retrieved (and most likely recompiled).  


**If you have git installed:** To obtain the directory structure, simply retrieve it using git as follows to clone the directory structure:

    $ git clone https://github.com/agnwinds/python.git -b structure

You then need to cd to the new directory and set your environment variables
    
    $ export PYTHON = /path/to/python/
    $ cd $PYTHON 
    $ make install
    $ make clean

note that export syntax is for bash- for csh use 
  
    $ setenv PYTHON /path/to/python/


**Without git:** Use the ZIP function under the [structures](https://github.com/agnwinds/python/tree/structure "Structure") branch, and then download .tar.gz versions of the python source under [releases](https://github.com/agnwinds/python/releases).

Once you have a directory /path/to/python/ which contains the structure, place the unpacked tar.gz python source folder under /path/to/python/progs/

    $ export PYTHON = /path/to/python/
    $ cd $PYTHON 
    $ make GIT=False install
    $ make clean
    $ cd progs/python_xx #replace xx  version you download
    $ make clean
    $ make CC=gcc python       # if you want to use mpicc, ignore the CC=gcc
    $ make clean

Again, for csh use 
  
    $ setenv PYTHON /path/to/python/

As you can tell, the git install is simpler!

Please see the (wiki)[https://github.com/agnwinds/python/wiki/Installing-and-Running-Python] for how to use the code.

Any comments, email jm8g08@soton.ac.uk or long@stsci.edu.