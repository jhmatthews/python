#!/usr/bin/env python 

'''
Combine Chianti levels to create a more succinct MacroAtom model


Command line usage (if any):

    usage: MacroCombine.py [-guess] level_file line_file phot_file

    with
        -guess implies that the program will try to guess
            what levels should be combined.  In the absence
            of this, the program will check if the 
            levels one wants to combine have beeen
            set already.  


Description:  

    This routine attempts to reconfigure the level, line
    and phot files files to allow for the fact that one would often like
    to combine levels to make MacroAtoms smaller

    There are two basically possibilities:

    with the -guess option, the program takes the original level file, in 
    the format produced by MakeMacro, and guesses what levels to combine, 
    producing an 'intermediate level file' that contains information 
    extra columns indicating what to combine.

    It then continues to produce final level fileline annd
    photoionization files that reflect this choice.

    Alternatively, without the -guess option the routine will read
    epects to read a level file which contains the extra columns
    that indcicate what needs to combine,  and uses this to produce
    final level, line and photoionization files.

Primary routines:


Notes:
                                       
History:

231013 ksl Coding begun

'''

import sys
from astropy.io import ascii
from astropy.table import join, Table
import numpy as np
import matplotlib.pyplot as plt
import os


def read_phot(file='h_1_phot.dat'):
    '''
    Read a Python photoionization file and return 
    the information in a form the is eacy to manipulate
    
    The routine returns a table that summarizes all of the
    different x-sections, and a list of tables, one table for
    each x-section.
    '''
    f=open(file)
    lines=f.readlines()
    i=0
    summary=[]
    xsections=[] # This is going to have one table for each set of xsections
    i=0
    for line in lines:
        word=line.split()
        if word[0][0]=='#':
            pass
        elif word[0].count('PhotMacS'):
            elem=int(word[1])
            ion=int(word[2])
            level=int(word[3])
            level_up=int(word[4])
            ethresh=float(word[5])
            nlines=int(word[6])
            line_no=i
            one_record=[elem,ion,level,level_up,ethresh,nlines,line_no]
            summary.append(one_record)
            if len(summary)>1:
                # The we have already accumulated xsections
                x=Table(np.array(one_xsection),names=['e','sigma'])
                xsections.append(x)
            one_xsection=[]
        else:
            one_xsection.append([eval(word[1]),eval(word[2])])
        i+=1
    x=Table(np.array(one_xsection),names=['e','sigma'])
    xsections.append(x)
    
    summary=np.array(summary)
    # summary=np.transpose(summary)
    tab_sum=Table(summary,names=['z','ion','ll','ul','ethresh','nlines','line_no'])
    print(len(tab_sum),len(xsections))
    
    return tab_sum,xsections



def write_phot(outfile,xtab,xsec):
    '''
    Write a valid Macro photometry xsection file from what amounts
    to an internal format for this routine.  The
    internal form consists of 

    * a table (xtab) which summarizes each x-section and 
    * a list (xsec) that contains tables of the xsections
    for each x-section
    '''

    if len(xtab)!=len(xsec):
        print('Error: Trying to write a photometry file witht the number of headers (%d) != number of x-sections (%d)' % (len(xtab),len(xsec)))
        return
    
    f=open(outfile,'w')
    i=0
    while i<len(xtab):
        one=xtab[i]
        string='PhotMacS  %3d %3d %3d %3d %10.6f %3d' %  (one['z'],one['ion'],one['ll'],one['ul'],one['ethresh'],one['nlines'])
        print(string)
        f.write('%s\n' %string)
        one_x=xsec[i]
        # print('foo \n',one_x)
        j=0
        while j<len(one_x):
            f.write('%10.6f %10.6e\n' % (one_x['e'][j],one_x['sigma'][j]))
            j+=1
        i+=1
    


def plot_xsec(selection,xtab,style='-',clear=True):
    '''
    Heree 

    * xtab is a list of tables containing the xsections and 
    * selection is a list of the x-sections to plot

    Note: section is 0-based and does not necessarily
    correspond to a specific level

    Aslo not that if one wants to save the plot, one 
    must currently do that separately
    '''
    isel=np.array(selection)
    plt.figure(1,(6,6))
    if clear:
        plt.clf()
    for one in isel:
        xcross=xtab[one]
        plt.semilogy(xcross['e'],xcross['sigma'],style,label='%d' % one, alpha=0.5)
    
    plt.xlabel('E (ev)')
    plt.ylabel(r'$\sigma$ (cm$^{-2}$)')
    plt.legend()
    plt.tight_layout()
    return



def reweight(levels,g,xtab):
    '''
    Reweight a group of of photionsation xsections, based on g

    xtab contains the entire list of xsections, each
    element of which corresponds to the x=sections
    for a single transition.

    level contains a list, e.g [1,2,3]  of the xsections to consider

    and g is a listi, e.g [2,2,4] containign the g factors associated
    with the various x-sections.

    Note: As written this is zero based and this may
    not be what we want.

    '''
    # first check that the enegies are the same
    
    xsig=[]
    i=0
    for one in levels:
        # print(xtab[one])
        if len(xsig)==0:
            xsig=g[i]*np.array(xtab[one]['sigma'])
        else:
            xsig+=g[i]*np.array(xtab[one]['sigma'])
        i+=1
    xsig/=np.sum(g)
    e=np.array(xtab[levels[0]]['e'])
    qtab=Table([e,xsig],names=['e','sigma'])
    qtab['e'].format='.6f'
    qtab['sigma'].format='.7e'
    # print(qtab)
    return qtab


def redo_phot(xlevel_file='h_1_lev2phot_guess.txt',phot_file_orig='h_1_phot.dat',root_out='foo'):
    '''
    Read a file containing the levels we want to compress amd  the associted phot_file
    and produce a new phot file
    
    The photometry file reads data into a summary table (tab_sum) and a list of tables with the photoionization data (xsections)
    We need to construct a new summary tabel (tab_sum_out) and a new lst of talbes (xsections)

    This routine alos writes out a new configuration file.
    '''

    if root_out=='':
        word=phot_file_orig.split('_')
        if len(word)>2:
            root_out='%s_%s' % (word[0],word[1])
        else:
            word=xlevel_file.split('_')
            if len(word)>2:
                root_out='%s_%s' % (word[0],word[1])
            else:
                root_out=foo
                print('No root out provided and no obvious split for inputs files so using %s' % root_out)
         
    
    try:
        xlevel=ascii.read(xlevel_file)
    except:
        print('Could not read %s' % level_file)
        return
   
    tab_sum,xsections=read_phot(phot_file_orig)
    
    print('The original num of levels  is ',len(xlevel))
    print('The number of xsections     is ',len(tab_sum))
    
    # Strip off the next ion if it is there
    xlevel=xlevel[xlevel['config']!='Next']
    nlev=xlevel['xlev'][-1]
    print('There are %d new levels' % nlev)
    
    #create a place to put a revised level_file
    lev_out=xlevel[0:nlev]
    tab_sum_out=tab_sum[0:nlev]
    xsections_out=[]
    j=0
    i=1
    while i<=nlev:
        cur_lev=xlevel[xlevel['xlev']==i]
        
        # Now we need to find the corresponding set of objects 
        # in the tab_sum and xsections file
        lev_out[j]=cur_lev[0]
        
        gotcha=False
        for one in tab_sum:
            if one['ll']==cur_lev['lvl'][0]:
                if gotcha==False:
                    tab_sum_out[j]=one
                    gotcha=True
        if gotcha==False:
            print('Failed')
        # Now get the xsections
        
        if len(cur_lev)==1:
            xsections_out.append(xsections[cur_lev['lvl'][0]-1])
        else:
            xlev=np.array(cur_lev['lvl'])
            glev=np.array(cur_lev['g'])
            xlev=xlev-1
            vsections=reweight(xlev,glev,xsections)
            xsections_out.append(vsections)
                          
        i+=1
        j+=1
    
    lev_out['g']=lev_out['G']
    i=0
    while i<nlev:

        lev_out['lvl'][i]=i+1
        tab_sum_out['ll'][i]=i+1
        i+=1        
        
    
    phot_file='%s_phot_final.dat' %   root_out
    lev_file='%s_lev_final.dat' % root_out
    write_phot(phot_file,tab_sum_out,xsections_out)
    lev_out.write(lev_file,format='ascii.fixed_width_two_line',overwrite=True)



def add_gtot(filename='he_2_levels_comb.dat'):
    ''' 
    Add the final G factors to the for
    the combined levels to the level
    file
    
    The table that is produced still has 
    the same number of levels as the original 
    table.  It simply has an additional
    column labelled G that will be the 
    multiplicty associated with the final level
    
    '''
    
    x=ascii.read(filename)
    # x.info()
    x['G']=-1
    levels=np.array(np.unique(x['xlev']))
    # print(levels)
    gg=np.zeros(len(levels))
    xion=x['Ion'][0]  # This is to allow us to avoid the next level up
    for one in x:
        if one['Ion']==xion:
            # print('ok ',one)
            gg[one['xlev']-1]+=one['g']
    
    #print(gg)
    
    # Now we have the sums, we just need to add them to the file
    
    for one in x:
        if one['Ion']==xion:
            # print('ok ',one)
            one['G']=gg[one['xlev']-1]

    x.write(filename,format='ascii.fixed_width_two_line',overwrite=True)
    return filename
  
    


def redo_lines(master='h_1_levels_comb2.dat',lines='h_1_lines.dat'):
    '''
    Write out a new line file, based on the intermidate
    level file
    '''
    xmaster=ascii.read(master)
    xmaster.rename_column('Element','z')
    xmaster.rename_column('Ion','ion')
    xlines=ascii.read(lines)
    # xlines['ll_final']=-99
    # xlines['ul_final']=-99
    print(len(xlines))
    
    xlow=xmaster['z','ion','lvl','xlev','G']
    xlow.rename_column('lvl','ll')
    xlow.rename_column('xlev','xll')
    xlow.rename_column('G','Gl')

    xlow.write('low.txt',format='ascii.fixed_width_two_line',overwrite=True)
    
    xup=xmaster['z','ion','lvl','xlev','G']
    xup.rename_column('lvl','ul')
    xup.rename_column('xlev','xul')
    xup.rename_column('G','Gu')
    
    xlines=join(xlines,xlow,join_type='left')
    print(len(xlines))
        
    xlines=join(xlines,xup,join_type='left')
    
    print(len(xlines))
    
    xlines['A']=xlines['f']*xlines['gl']/xlines['gu']
    xlines['A'].format='.6f'
    
    xlines['x']=xlines['xll']*1000+xlines['xul']
    
    
    xlines.sort('ul')
    xlines.info()
    xlines.write('foo_lines.dat',format='ascii.fixed_width_two_line',overwrite=True)
    
    z,xindex=np.unique(np.array(xlines['x']),return_index=True)
    
    xfinal=xlines[xindex]

    
    print(xfinal)
    
    # Now make a new table that contains only one row for each value in z

        
    for one_row in xfinal:
        all_rows=xlines[xlines['x']==one_row['x']]
        one_row['Wave']=np.average(all_rows['Wave'])
        one_row['A']=np.sum(all_rows['A']*all_rows['gu'])/one_row['Gu']
        one_row['f']=one_row['Gu']*one_row['A']/one_row['Gl']
        # one_row['f']=np.sum(all_rows['f']*all_rows['gu'])/np.sum(all_rows['gu'])
        one_row['gu']=np.sum(all_rows['gu'])
        one_row['ul']=np.sum(all_rows['ul'])
        one_row['eu']=np.average(all_rows['eu'])
        
    xfinal['f'].format=xfinal['eu'].format='.5f'
    xfinal['Wave'].format='.5f'


    words=lines.split('.')
    out_name='%s_final.%s' % (words[-2],words[-1])

    
    
    xxx=xfinal['Dtype','z','ion','Wave','f','Gl','Gu','el','eu','xll','xul']
    xxx.write(out_name,format='ascii.fixed_width_two_line',overwrite=True)
    return xxx



    
    
def xguess(tab_name='h_1_lev2phot.txt',out_name='',frac=0.01,redo=True):
    '''
    Guess how levels for an atom should be combined and write out
    the intermediate level file that contains the various suggestions

    '''
    xtab=ascii.read(tab_name)
    cols=xtab.colnames
    for one in cols:
        if one=='xlev'and redo!=True:
            print('This file %s is already has levels to combine, so returning' % tab_name )
            print('To force a repeat, set redo to True')
            return
    xtab['xlev']=-99
    
    i=1
    ilev=1
    xpot=xtab['ion_pot'][0]
    xtab['xlev'][0]=1
    while i<len(xtab):
        cpot=xtab['ion_pot'][i]
        delta=np.fabs((cpot-xpot)/cpot)
        if delta<frac:
            xtab['xlev'][i]=ilev
        else:
            ilev+=1
            xtab['xlev'][i]=ilev
            xpot=cpot
        i+=1
        
    if out_name=='':        
        words=tab_name.split('.')
        out_name='%s_guess.%s' % (words[-2],words[-1])
        
    xtab.write(out_name,format='ascii.fixed_width_two_line',overwrite=True)
    
    return out_name


def steer(argv):
    '''
    This is just steering routine to allow operations from the command line
    '''

    level_file=''
    line_file=''
    phot_file=''
    iguess=False

    i=1
    while i<len(argv):
        if argv[i]=='-h':
            print(__doc__)
            return
        if argv[i]=='-guess':
            iguess=True
        elif argv[i][0]=='-':
            print(__doc__)
            print('Unknown option %s' % argv[i])
        elif level_file=='':
            level_file=argv[i]
        elif line_file=='':
            line_file=argv[i]
        elif phot_file=='':
            phot_file=argv[i]
        else:
            print(__doc__)
            print ('Too many arguments' % argv)
            return

        i+=1

    if level_file=='':
        print('No level file provided')
        return

    if iguess:
        level_file=xguess(level_file)
        add_gtot(level_file)

    if line_file!='':
        add_gtot(level_file)
        redo_lines(level_file,line_file)
    else:
        print('No line file provided')
        return

    if phot_file!='':
        redo_phot(xlevel_file=level_file,phot_file_orig=phot_file,root_out='')
    else:
        print('No phot file provided')

                



# Next lines permit one to run the routine from the command line
if __name__ == "__main__":
    import sys
    if len(sys.argv)>1:
        steer(sys.argv)
    else:
        print (__doc__)

