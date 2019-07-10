# Python

Python is a (confusingly named) Monte Carlo radiative transfer code which uses the Sobolev approximation.
It has been developed by Knox Long, Christian Knigge, Stuart Sim, Nick Higginbottom, James Matthews, Sam Mangham and Edward Parkinson.
The code has been used for a variety of research projects invovling the winds of cataclysmic variables, of young stellar 
objects, and of AGN.

The code is under active development, but we are looking for beta users to test the code, and potentially use it 
for their own research. If you are interested in using Python please email the collaboration at agn-winds[at]googlegroups[dot]com or contact Knox Long via long[at]stsci[dot]edu. 

Documentation of input parameters is hosted on [ReadTheDocs](http://agnwinds.readthedocs.io/en/dev/).

## Travis Build Status

[![Build Status](https://travis-ci.org/agnwinds/python.png?branch=dev)](https://travis-ci.org/agnwinds/python)

## Installation

Python and the various routines associated are set up in a self-contained directory structure. The basic directory structure and the data files that one needs to run Python need to be retrieved and compiled. 

If you want to obtain a stable (!) release, go to the [Releases](https://github.com/agnwinds/python/releases) page.

If you want to download the latest dev version, you can zip up the git repository by clicking on the zip icon to the right of the GitHub page. Aternatively, you can clone the repository using 

    $ git clone https://github.com/agnwinds/python.git 

If you anticipate contributing to development we suggest Forking the repository and submitting pull requests with any proposed changes.

Once you have the files, you need to cd to the new directory and set your environment variables
    
    $ export PYTHON = /path/to/python/
    $ cd $PYTHON 
    $ ./configure
    $ make install  (or better make install 2>&1 | tee today.txt)
    $ make clean

If you have any difficulties with the installation, please submit an issue, along with the file today.txt

Note that the export syntax is for bash- for csh use 
  
    $ setenv PYTHON /path/to/python/

Atomic data is stored in our [data repository](https://github.com/agnwinds/data) with its own releases page. one should unzip these files and place them in a $PYTHON/data folder.

A development user may want to work on atomic data as part of their work, and pull in changes as they are made, in which case we recommend cloning the data repository:

    $ cd $PYTHON; git clone https://github.com/agnwinds/data data

## Running python

To run python you need to add the following to your $PATH variable:

    $PYTHON/bin

You can then setup your symbolic links by running 

    $ Setup_Py_Dir

and run the code by typing, e.g.

    $ py root.pf


Please see the [wiki](https://github.com/agnwinds/python/wiki/Installing-and-Running-Python) and docs folder for how to use the code.

Any comments, email one of the addresses above.
