#!/usr/bin/env python

'''

Synopsis:

Find all of the python files in a directory and write out html
pydoc documentation of them in the current working directory


Command line usage (if any):

        write_docs.py [directory_name]

                if the directory name is missing then assume
                the routine is being run from the pydocs directory
                and the the directory to be documented is ../../py_progs

Description:

    The routine finds all of the .py files in the named directory and
    generates a shell script which runs pydocs to create html help
    files for each .py file.  Assuming the directory name is mydir
    It then creates a another html file, doc_mydir.html which simply
    contain links to all of the individual files



Primary routines:

Notes:

    If you want to get an html help for a specfic package the python
    command is simply
        pydocs package_name
    The help html files have links to the packages such as pyfits
    but these are incorrect because pydocs assumes all the help
    files are in the current working directory, which they are likely
    not.

    The does not delete help from routines that have been deleted from
    a package, which is a problem.  One could fix this, but one would
    then need to keep some kind of database so that you knew what files
    should exist.

    Note that you should watch the written output.  If it gives
    and odd looking error which writes out some of the code it
    means that one has not put

    if __name__ == "__main__":

    in front of the main routine.  If this happens you will need
    to indent all of the lines in main, and then rerun this
    script



History:

100805    ksl    Coding begun
111111    ksl    Modified so it actually used the entire name for the file.
        pydoc will then use that file to make the doucmentation
        Otherwise it searches for it in the path and that is less
        likely to be what I want.
180126  ksl Modified for use with python

'''

import sys
import os
import glob
from MarkupPy import markup


def make_toplevel(dirname, names):
    '''
    Create an html page to point to all
    of the individual help pages that have
    been made
    '''
    # First get the name of the directory
    dirname = dirname.replace('/', ' ')
    dirname = dirname.strip()
    dirname = dirname.split()
    dirname = dirname[len(dirname)-1]
    print(dirname)

    html_name = 'doc_'+dirname+'.html'

    # Start a page
    page = markup.page( )

    page.init(title="Documentation for %s" % dirname)

    page.h1("Documentation for python scripts in the directory  %s" % dirname)

    items = []
    for name in names:
        item = markup.oneliner.a(name, href="./%s.html" % name)
        items.append(item)

    page.ul(class_='mylist' )
    page.li(items, class_='myitem' )
    page.ul.close()

    page.p('Warning: This page is rewritten every whenever write_docs.py is run on this directoryand so this page should not be edited')
    g = open(html_name, 'w')
    g.write('%s' % page)
    g.close()


def write_docs(dirname='../../py_progs'):
    '''
    Locate all of the .py files in dirname and
    write out help in the current working directory
    using pydocs
    '''

    search_name = dirname+'/*.py'
    names = glob.glob(search_name)
    print(names)
    g = open('DoDocs', 'w')

    roots = []

    for name in names:
        words = name.split('/')
        filename = words[len(words)-1]
        words = filename.split('.')
        root = words[0]
        roots.append(root)
        # g.write('pydoc -w %s\n' % root)
        g.write('pydoc -w %s\n' % name)
    g.close()
    os.system('source DoDocs')

    # Now make a page that points to all the html pages
    # that have already been made
    make_toplevel(dirname, roots)


# Next lines permit one to run the routine from the command line
if __name__ == "__main__":
    if len(sys.argv) == 1:
        write_docs()
    elif len(sys.argv) == 2 and sys.argv[1] == '-h':
        print(__doc__)
    elif len(sys.argv) == 2:
        write_docs(sys.argv[1])
    else:
        print('usage: write_docs.py  dirname or -h for info')
