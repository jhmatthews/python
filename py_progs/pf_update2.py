#!/usr/bin/env python 
'''
	University of Southampton -- JM -- October 2014

				pf_update2.py

Synopsis:
	pf_update2.py updates a python parameter file 

Usage:
    python [-o -c] pf_update2.py new_format.pf old_list_of_files
	
Arguments:
    -o      overwrite with same filenames
    -c      copy to format New_root.pf

Notes:
    Would highly recommend having the OrderedDict object as part 
    of the collections module. Should come standard with most distributions.
'''

# we need the classes and numpy modules 
import py_classes as cls
import numpy as np
import matplotlib.pyplot as plt
import subprocess
import py_plot_util as util
import py_read_output as r
import sys




# Next lines permit one to run the routine from the command line with various options -- see docstring
if __name__ == "__main__":

    # read the arguments
    if len(sys.argv) == 3:
        print "Defaulting to copy mode (-c). Will not overwrite existing pf files."
        args = 0
    elif len(sys.argv) == 4:
        args = 1
    else:
        print "Didn't understand arguments."
        print __doc__
        sys.exit()

    # find out what mode the user wants to run in 
    copy = True
    if args > 0:
        mode = sys.argv[1]
        if mode == "-c":
            print "copy mode (-c). Will not overwrite existing pf files."
        elif mode == "-o":
            print "overwrite mode (-o). Overwrites existing pf files."
            copy = False
        else: 
            print "Didn't understand arguments."
            print __doc__
            sys.exit()



    #get the template filename and read it into an OrderedDict object if you can
    template_fname = sys.argv[1 + args]     #
    template_dictionary = r.read_pf(template_fname)

    # load the list of filenames to convert
    list_fname = sys.argv[2 + args]
    fnames = np.loadtxt(list_fname, dtype="string")


    n_fnames = len(fnames)

    print "working on %i files" % n_fnames

    # iterate over the list of files
    for i in range(n_fnames):

        # create our two OrderedDict objects
        new_dict = template_dictionary
        old_dict = r.read_pf(fnames[i])

        #replace values in new_dict with old_dict where possible
        for key,val in new_dict.iteritems():

            try:
                new_dict[key] = old_dict[key]
            except KeyError:
                print "file %s has no key %s: using template value" % (fnames[i], key)

        if copy:
            r.write_pf ( "New_%s" % fnames[i], new_dict)

        else:
            r.write_pf ( "%s" % fnames[i], new_dict)

    # all done, hopefully...
    print "All done."








