#!/usr/bin/env python 

'''
                    Space Telescope Science Institute

Synopsis:  

Sumarize a model run with python, ultimately generating
an html file with various plots, etc.


Command line usage (if any):

    usage: run_check.py root

Description:  

Primary routines:

    doit

Notes:
                                       
History:

190312 ksl Coding begun

'''

import sys
from astropy.io import ascii
import numpy
import subprocess
import pylab
import xhtml
import plot_wind
import plot_spec
import plot_tot


def read_diag(root):
    '''
    Get convergence and possibly other information from the diag file
    '''

    filename='diag_%s/%s_0.diag' % (root,root)

    command="grep 'Summary  convergence' %s" % filename

    proc=subprocess.Popen(command,shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    stdout,stderr=proc.communicate()

    if len(stderr):
        print("Error - grep")
        print(stderr.decode())
        return [],[]
    else:
        print(stdout.decode())
        x=stdout.decode()
        # print(len(x))
        lines=x.split('\n')
        # print(len(lines))
        converging=[]
        converged=[]
        for line in lines:
            words=line.split()
            # print(words)
            if len(words)>5:
                converged.append(eval(words[3]))
                converging.append(eval(words[5]))
        return converged,converging


def windsave2table(root):
    '''
    Run windsave2table
    '''

    command='windsave2table %s' % root

    proc=subprocess.Popen(command,shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    stdout,stderr=proc.communicate()
    # print(stdout.decode())
    # print('Error')
    # print(stderr.decode())
    # print('return code')
    # print (proc.returncode)
    if proc.returncode:
        print('Error: running windsave2table',proc.returncode)
        return True
    elif len(stderr):
        print('Error: running windsave2table')
        print(stderr.decode())
        return True 
    else:
        return False



def plot_converged(root,converged,converging):
    '''
    Make a plot of the convergence statistics in the diag directroy
    '''

    pylab.figure(1,(6,6))
    pylab.clf()
    pylab.plot(converged,label='Converged')
    pylab.plot(converging,label='Converging')
    pylab.legend(loc='best')
    pylab.xlabel('Ionization Cycle')
    pylab.ylabel('Fraction of Cells in Wind')
    pylab.savefig('diag_%s/convergence.png' % root)

    return


    
def make_html(root,converge_plot,te_plot,tr_plot,spec_tot_plot,spec_plot):
    '''
    Make an html file that collates all the results
    '''

    string=xhtml.begin('Evaluation of how well the python run of %s succeeded' % root)

    string+=xhtml.paragraph('Provide an overview of whether the run of %s has succeeded' % root)

    string+=xhtml.hline()
    string+=xhtml.h2('Did the run converge?')
    # print(xhtml.image('file:./diag_%s/convergence.png' % root))
    string+=xhtml.image('file:./diag_%s/convergence.png' % root)

    string+=xhtml.paragraph('Which cells converged? (0 indicates success)')
    # print('Test ',converge_plot)
    string+=xhtml.image('file:%s' % (converge_plot))
    string+=xhtml.hline()

    string+=xhtml.h2('What do the temperatures look like?')
    string+=xhtml.image('file:%s' % (te_plot))
    string+=xhtml.image('file:%s' % (tr_plot))
    # string+=xhtml.image('file:./diag_%s/%s_t_e.png' % (root,root))
    # string+=xhtml.image('file:./diag_%s/%s_t_r.png' % (root,root))
    string+=xhtml.hline()
    string+=xhtml.h2('What do the total spectra look like (somewhat smoothed)?')
    string+=xhtml.image('file:%s' % (spec_tot_plot))
    string+=xhtml.hline()
    string+=xhtml.h2('What do the final spectra look like (somewhat smoothed)?')
    string+=xhtml.image('file:%s' % (spec_plot))
    string+=xhtml.hline()

    # Now add the parameter file so we can see what we have

    string+=xhtml.h2('The parameter file used:')

    x=open(root+'.pf')
    lines=x.readlines()

    string+=xhtml.add_list(lines)
    string+=xhtml.hline()



    string+=xhtml.end()

    # print(string)

    g=open(root+'.html','w')
    g.write(string)



def doit(root='ixvel',outputfile='out.txt'):
    '''
    Create a summary of a Python run, which will enough information that one can assess
    whether the run was successful

    Description:

    Notes:

    History:


    '''
    if windsave2table(root):
        print('Exiting becase windsave2table failed')
        return

    converge_plot=plot_wind.doit('%s.0.master.txt' % root,'converge',plot_dir='./diag_%s' % root)
    te_plot=plot_wind.doit('%s.0.master.txt' % root,'t_e',plot_dir='./diag_%s' % root)
    tr_plot=plot_wind.doit('%s.0.master.txt' % root,'t_r',plot_dir='./diag_%s' % root)


    converged,converging=read_diag(root)

    if len(converged):
        plot_converged(root,converged,converging)

    plot_tot.doit(root)
    spec_tot_plot=root+'.spec_tot.png'

    plot_spec.do_all_angles(root,wmin=0,wmax=0)
    spec_plot=root+'.png'

    make_html(root,converge_plot,te_plot,tr_plot,spec_tot_plot,spec_plot)


    
    



    return





# Next lines permit one to run the routine from the command line
if __name__ == "__main__":
    import sys
    if len(sys.argv)>1:
        doit(sys.argv[1])
    else:
        print (__doc__)
