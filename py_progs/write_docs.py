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

    If you want to get an html help for a specific package the python
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
import pydoc
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

    html_name = 'doc_py_progs.html'

    # Start a page
    page = markup.page()

    page.init(title="Documentation for %s" % dirname)

    page.h1("Documentation for python scripts in the directory  %s" % dirname)

    page.p('''This page lists the scripts that exist.  Some of these scripts will be useful to users, and some
    will not.
    At present, there is no obvious way to tell, except to look read the documentation associated with
    each script by following the link, or to have seen a reference to a particular script in some other
    place in the documentation set.
    ''')

    page.p('''To use the scripts, one will need to have the py_progs directory in their PYTHONPATH.  Occasionally
    one will need to be prepared to install modules using conda or pip.''')

    page.h2('The scripts')

    items = []
    for name in names:
        item = markup.oneliner.a(name, href="./%s.html" % name)
        items.append(item)

    page.ul(class_='mylist')
    page.li(items, class_='myitem')
    page.ul.close()

    page.p('Warning: This page is rewritten every whenever write_docs.py is run on this directoryand so this page should not be edited')
    with open(html_name, 'w') as file:
        file.write('%s' % page)
        file.close()


def write_docs(dirname='../../py_progs'):
    '''
    Locate all of the .py files in dirname and
    write out help in the current working directory
    using pydocs
    '''
    # First, we delete all the existing documentation in this directory
    for item in os.listdir('.'):
        if item.endswith(".html"):
            os.remove(item)

    # Now, we write new docs
    pydoc.writedocs(dirname)

    # Now make a page that points to all the html pages
    # that have already been made
    roots = [
        item.replace('.py', '')
        for item in os.listdir(dirname)
        if item.endswith('.py')
    ]
    make_toplevel(dirname, roots)

    # Now check that we have the files we expected
    got_all = True
    for root in roots:
        if not os.path.isfile(root+'.html'):
            print('Failed to create an html file for %s.py' % root)
            got_all = False

    if got_all:
        print('html files were created for all of the .py scripts')
        return 0

    print(
        'Failed to generate documentation for some files.\n'
        'Please look at the earlier output to find more details on the error.'
    )
    return 1  # Return nonzero to indicate error


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
