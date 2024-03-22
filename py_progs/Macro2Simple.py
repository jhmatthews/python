#!/usr/bin/env python 

'''
                    Space Telescope Science Institute

Synopsis:  

Create versions of a MacroAtom imputs that can be used in 
Simple atom calculations.

This must be run in the same directory as MakeMacro.py


Command line usage (if any):

    usage: Macro2Simple.py filename

Description:  

Primary routines:

    doit

Notes:
                                       
History:

240127 ksl Coding begun

'''

from MakeMacro import get_levels
from astropy.io import ascii
import RedoPhot
import ChiantiPy.core as ch



def write_simple_levels(ion="h_1", nlevels=10):
    xtab=get_levels(ion,nlevels)
    xtab['Dtype']='LevTop'
    xtab['eqn']=-99.0
    
    ztab=xtab['Dtype','Element','Ion','islp','ilv','ion_pot','ex','g','eqn','rad_rate','config']
    outfile='%s_levels_simple.dat' % ion
    ztab.write(outfile,format='ascii.fixed_width_two_line',overwrite=True)
    return outfile



def write_simple_lines(ion="h_1"):
    infile='%s_lines.dat' % ion
    outfile='%s_lines_simple.dat' % ion
    try:
        x=ascii.read(infile)
    except:
        print('Error: write_simple_lines: Could not read %s' % infile)
        return []
    x['Dtype']='Line'
    x.write(outfile,format='ascii.fixed_width_two_line',overwrite=True)
    return outfile
    



def doit(ion='h_1',nlevels=10):
    lev_file=write_simple_levels(ion,nlevels)
    print(lev_file)
    line_file=write_simple_lines(ion)
    print(line_file)
    x = ch.ion(ion, temperature=1e5)
    nelec = x.Z - x.Ion + 1
    filename = "p%02d.%02d.dat" % (x.Z, nelec)
    outfile = "%s_phot_simple.dat" % (ion)  
    phot_tab,xtab=RedoPhot.extrap(filename,1e5)
    if len(phot_tab)==0:
        print('Error: Exiting because of previous problems')
        return    
    RedoPhot.write_phot_tab(outfile,phot_tab,xtab)
    print(outfile)
    


# Next lines permit one to run the routine from the command line
if __name__ == "__main__":
    import sys
    if len(sys.argv)==3:
        doit(sys.argv[1],int(sys.argv[2]))
    elif len(sys.argv)==2:
        doit(sys.argv[1],10)
    else:
        print (__doc__)



