#!/usr/bin/env python 

'''
                    Space Telescope Science Institute

Synopsis:  

A simple program to plot and smooth models made with Python.


Command line usage (if any):

    usage: plot_spec.py [-wmin 850 -wmax 1850 -smooth 21] rootname

Description:  

    This is a routine to plot the detailed spectra produced with Python.  The rootine reads
    the .spec file and produces a smoothed plot for all of the angles that have been extracted.

    Normally, the routine would be invoked from the command line with the rootnmae of the Python run.

    A plot rootname.png is produced. With none of the command line options the spectrum will comver the full
    spectral range of the detailed spectra, and will be box-car smmothed by 11 pixels.

    There are various command line options

    =h  prints a help based on __docs__
    -wmin sets the minimum wavelength of the plot
    -wmax sets the maximimy wavelength for the plot
    -smooth changes the binning in pixels

Primary routines:

    steer  Processes the command line options and then calls do_all_angles
    do_all_angles Produces the plot  (do_all_angles(root) is the routine to call from a script

Notes:
                                       
History:

140222 ksl Coding begun
1802    ksl Modified so that it worked better from the command line, and added steering routine.

'''

import sys
import numpy
import pylab
from scipy.signal import boxcar
from scipy.signal import convolve
from scipy.optimize import leastsq
# from ksl import io
from astropy.io import ascii

def read_file(filename,char=''):
    '''
    Read a file and split it into words, eliminating comments
    
    char is an optional parameter used as the delimiter for
    splitting lines into words.  Otherwise white space is
    assumed.

    History:
    
    110729    ksl    Added optional delimiters
    '''

    try:
        f=open(filename,'r')
        xlines=f.readlines()
        f.close()
    except IOError :
        print("The file %s does not exist" % filename)
        return []   
    
    lines=[]
    
    i=0
    while i<len(xlines):
        z=xlines[i].strip()
        if char=='':
            z=z.split()
        else:
            z=z.split(char)
        if len(z)>0:
            if z[0][0]!='#':
                j=0
                # print z
                while j<len(z):
                    z[j]=eval(z[j])
                    j=j+1
                lines.append(z)
        i=i+1
    return lines

def get_column_names(filename='sv.spec'):
    '''
    Silly routine to read a spec file and return the column names.
    It's silly because we need to figure out an easier way to return
    this for all files
    '''
    try:
        f=open(filename,'r')
        xlines=f.readlines()
        f.close()
    except IOError :
        print("The file %s does not exist" % filename)
        return []   
    
    words=[]
    for line in xlines:
        if line.count('# Freq.'):
            line=line.strip()
            line=line.split()
            words=line[1:]
    return words



def do_all_angles_ev(rootname='sv',smooth=21,emin=1000,emax=9000,fmax=0,fig_no=1):
    '''
    Plot each of the spectra where

    smooth is the number of bins to boxcar smooth
    fmax sets the ylim of the plot, assuming it is not zero

    and other values should be obvious.

    160405 ksl This version of the routine is intended to plot the X-ray 
        regime
    '''


    if rootname.count('.')>0:
        filename=rootname
        rootname=rootname[0:rootname.rindex('.')]
        # print('test',rootname)
    else:
        filename=rootname+'.spec'


    try:
        data=ascii.read(filename)
    except IOError:
        print('Error: Could not find %s' % filename)
        return


    print(data.colnames)

    HEV=4.136e-15

    # Now determine what the coluns containing real spectra are
    # while makeing the manes simpler
    cols=[]
    i=9
    while i<len(data.colnames):
        one=data.colnames[i]
        # print 'gotcha',one 
        new_name=one.replace('P0.50','')
        new_name=new_name.replace('A','')
        data.rename_column(one,new_name)
        cols.append(new_name)
        i=i+1



    root=rootname


    pylab.figure(fig_no,(9,6))
    pylab.clf()


    for col in cols:
        flux=data[col]
        # print(flux)
        xlabel=col+'$^{\circ}$'
        q=convolve(flux,boxcar(smooth)/float(smooth),mode='same')
        pylab.semilogx(data['Freq.']*HEV,q*data['Freq.'],'-',label=xlabel)
    pylab.xlabel(r'Energy (eV)')
    pylab.ylabel(r'$\nu F_{\nu}$')
    z=pylab.axis()
    if fmax==0:
        pylab.axis((emin,emax,0,z[3]))
    else:
        pylab.axis((emin,emax,0,fmax))

    pylab.title(root)
    pylab.legend(loc='best')
    pylab.draw()
    pylab.savefig(root+'.png')
    return



def do_all_angles(rootname='sv',smooth=21,wmin=850,wmax=1850,fmax=0,fig_no=1):
    '''
    Plot each of the spectra where

    smooth is the number of bins to boxcar smooth
    fmax sets the ylim of the plot, assuming it is not zero

    and other values should be obvious.

    140907    ksl    Updated to include fmax
    140912    ksl    Updated so uses rootname 
    141124    ksl    Updated for new formats which use astropy
    '''

    filename=rootname+'.spec'

    try:
        data=ascii.read(filename)
    except IOError:
        print('Errro: Could not find %s' % filename)
        return


    print(data.colnames)

    # Now determine what the coluns containing real spectra are
    # while makeing the manes simpler
    cols=[]
    i=9
    while i<len(data.colnames):
        one=data.colnames[i]
        # print 'gotcha',one 
        new_name=one.replace('P0.50','')
        new_name=new_name.replace('A','')
        data.rename_column(one,new_name)
        cols.append(new_name)
        i=i+1



    root=rootname

        # if wmin and wmax are not specified, use wmin and wmax from input file

    if wmin==0 or wmax==0:
        xwave=numpy.array(data['Lambda'])
        xmin=numpy.min(xwave)
        xmax=numpy.max(xwave)
        if wmin==0:
            wmin=xmin
        if wmax==0:
            wmax=xmax




    pylab.figure(fig_no,(9,6))
    pylab.clf()

    print(cols)

    for col in cols:
        flux=data[col]
        # print(flux)
        xlabel=col+'$^{\circ}$'
        q=convolve(flux,boxcar(smooth)/float(smooth),mode='same')
        pylab.plot(data['Lambda'],q,'-',label=xlabel)
    pylab.xlabel(r'Wavelength ($\AA$)')
    pylab.ylabel('Flux')
    z=pylab.axis()
    if fmax==0:
        pylab.axis((wmin,wmax,0,z[3]))
    else:
        pylab.axis((wmin,wmax,0,fmax))

    pylab.title(root)
    pylab.legend(loc='best')
    pylab.draw()
    pylab.savefig(root+'.png')
    return


def steer(argv):
    '''
    Parse the command line
    '''

    wmin=0
    wmax=0
    smooth=11

    rootname=argv[1]

    i=1
    while i<len(argv):
        if argv[i]=='-h':
            print(__doc__)
            return
        elif argv[i]=='-wmin':
            i+=1
            wmin=eval(argv[i])
        elif argv[i]=='-wmax':
            i+=1
            wmax=eval(argv[i])
        elif argv[i]=='-smooth':
            i+=1
            smooth=eval(argv[i])
        elif argv[i][0]=='-':
            print('Error: Unknown switch ---  %s' % argv[i])
            return
        else:
            rootname=argv[i]
        i+=1

    do_all_angles(rootname,smooth,wmin,wmax)




# Next lines permit one to run the routine from the command line
if __name__ == "__main__":
    import sys
    if len(sys.argv)>1:
        # doit(int(sys.argv[1]))
        steer(sys.argv)
        # doit(sys.argv[1])
    else:
        print(__doc__)
