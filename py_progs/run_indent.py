#!/usr/bin/env python 

'''
                    Space Telescope Science Institute

Synopsis:  

Indent in a controlled manner the .c files used by Python


Command line usage (if any):

    usage: run_indent.py filename to indent a single file
           run_indent.py *.c to indent all of the .c files in a directory
           run_indent.py *.h to indent all of the .h files
           run_indent -all to indent all of the c and .h files

Description:  

Primary routines:

    doit  processes a single file
    steer procesest the command calling either doi or do_all
    do_all processes all the .c and .h files in a directory

Notes:

    Files are not updated unless indent produces a new result
    
    If gnuindent is not found the program will do nothing
                                       
History:

180913 ksl Coding begun

'''

import sys
from astropy.io import ascii
import numpy
from glob import glob
import subprocess
from shutil import copyfile

def get_gnu():

    # Find gnu indent, if possible

    gnu=False
    indent=''

    # Brew installs gnu indent as gindent
    proc = subprocess.Popen('gindent --version', shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout,stderr=proc.communicate()
    stdout = stdout.decode().split('\n')
    if stdout[0].count('GNU'):
        gnu=True
        indent='gindent'


    # On linux machines indent is nearly always gnu
    proc = subprocess.Popen('indent --version', shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout,stderr=proc.communicate()
    stdout = stdout.decode().split('\n')
    if stdout[0].count('GNU'):
        gnu=True
        indent='indent'


    # Our original make file looked for this
    proc = subprocess.Popen('gnuindent --version', shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout,stderr=proc.communicate()
    stdout = stdout.decode().split('\n')
    if stdout[0].count('GNU'):
        gnu=True
        indent='gnuindent'

    if gnu==False:
        print('Error: could not find a gnu version of indent.')
        print('If using a Mac, one can use "brew install gnu-indent" to install gindent')
        indent =''

    return indent





def doit(filename='lines.c'):
    '''
    Properly indent a single file

    Description:

    Notes:

    If gnuindent is not found, then this function
    is a NOP.

    The indented version is first written to another
    file.  Then we check to see if the indented file
    is different from the original.  

    If the newly indented version is differnt from the
    original, then it is copied back to the original.

    This is done so that a file that is unchanged is not
    "touched', which would cause a recompilation when
    mske is uesed.

    History:


    '''
    indent=get_gnu()
    if indent=='':
        return

    indent_command='%s -gnu -l140 -bli0 -nut -o foo.txt %s' % (indent,filename)



    print(indent_command)
    proc = subprocess.Popen(indent_command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout,stderr=proc.communicate()
    stdout = stdout.decode().split('\n')
    # print(stdout)
    stderr = stderr.decode().split('\n')
    # print(stderr)


    # Now check whether anything has changed using diff

    proc = subprocess.Popen('diff %s foo.txt' % filename, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout,stderr=proc.communicate()
    stdout = stdout.decode().split('\n')
    # print(stdout)
    stderr = stderr.decode().split('\n')
    # print(stderr)

    if len(stdout)>1:
        print('There were differences after indent for file %s' % filename)
        copyfile('foo.txt',filename)
    # else:
    #     print('Indent made no changes for file %s' % filename)

    return



def do_all():
    '''
    Indent all of the .c and .h files in a directory in a standard way
    '''
    if get_gnu()=='':
        return

    files=glob('*.h')+glob('*.c')

    for one in files:
        doit(one)

    return





def steer(argv):
    '''
    Process the command line
    '''


    if get_gnu()=='':
        return

    files=[]
    i=1
    while i<len(argv):
        if argv[i]=='-h':
            print(__doc__)
            return
        if argv[i]=='-all':
            do_all()
            return
        else:
            files.append(argv[i])
        i+=1

    todo=[]
    for one in files:
        if one[len(one)-2:]=='.h' or one[len(one)-2:]=='.c':
            todo.append(one)

    if len(todo)==0:
        print('Error: Read command line but could find nothing to do')
        return

    for one in todo:
        doit(one)

    return







# Next lines permit one to run the routine from the command line
if __name__ == "__main__":
    import sys
    if len(sys.argv)>1:
        # doit(int(sys.argv[1]))
        # doit(sys.argv[1])
        steer(sys.argv)
    else:
        print (__doc__)
