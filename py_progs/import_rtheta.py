#!/usr/bin/env python 

'''
                    Space Telescope Science Institute

Synopsis:  

Read the master file produced by windwave2table for a
rtheta (polar-coordinate mode model and produce 
a file to be read by import in python


Command line usage (if any):

    usage: import_1d.py rootname

    e.g rtheta to read a the master.txt table

    rtheta.0.master.txt

Description:  

Primary routines:

    doit

Notes:
                                       
History:

171106 ksl Coding begun

'''

import sys
from astropy.io import ascii
import numpy


def read_file(filename,char=''):
    '''
    Read a file and split it into words, eliminating comments
    
    char is an optional parameter used as the delimiter for
    splitting lines into words.  Otherwise white space is
    assumed.

    History:
    
    110729    ksl    Added optional delimiters
    141209    ksl    Reinstalled in my standard startup
            script so there was flexibility to
            read any ascii file
    '''

    try:
        f=open(filename,'r')
        xlines=f.readlines()
        f.close()
    except IOError :
        print ("The file %s does not exist" % filename)
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
                lines=lines+[z]
        i=i+1
    return lines




def read_table(filename='foo.txt',format=''):
    '''
    Read a file using astropy.io.ascii and 
    return this 

    Description:

    Notes:

    History:


    '''
    try:
        if format=='':
            data=ascii.read(filename)
        else:
            data=ascii.read(filename,format=format)
        for col in data.itercols():
            if col.dtype.kind in 'SU':
                data.replace_column(col.name,col.astype('object'))
    except IOError:
        print ('Error: file %s does not appear to exist' % filename)
        return

    print ('Here are the column names:')
    
    print (data.colnames)

    return data


def doit(root='rtheta',outputfile=''):
    '''
    Read a master.txt file and produce a file
    which can be read in to python
    

    Description:

    Notes:

    History:

    '''

    filename=root+'.0.master.txt'
    if outputfile=='':
        outputfile=root+'.import.txt'
    


    data=read_table(filename)


    xdata=data['i','j','inwind','r','theta','v_x','v_y','v_z','rho','t_r']

    
    print (xdata)


    # This format is the easy to read back automatically
    ascii.write(xdata,outputfile,format='fixed_width_two_line')

    return


    





# Next lines permit one to run the routine from the command line
if __name__ == "__main__":
    import sys
    if len(sys.argv)>1:
        # doit(int(sys.argv[1]))
        doit(sys.argv[1])
    else:
        print ('usage: import_1d.py filename')
