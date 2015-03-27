#!/usr/bin/env python 
#specify your python interpreter path above to run as executable

'''

University of Southampton, James Matthews, 130722

"py_error.py"
	
This is the post-processing code to deal with the error logs in 
PYTHON's parallel mode.

usage and arguments:

	python py_error.py [-f folder] root
	root is the parameter files
	-f folder is the folder the diagfiles are in, defaults to diag_root

returns:

	Error summary printed to screen. Can pipe to file if required.

130722	JM	Coding began -- initial tests conducted successfully

'''

	
## import modules
import sys, subprocess

#help string for user
def help(flag):
	''' help function for user'''
	if flag=='help' or flag =='-h':
		print '''usage and arguments:

		python py_error.py [-f folder] root
		root is the parameter files
		-f folder is the folder the diagfiles are in, defaults to diag_root'''
		sys.exit(0)
	else: print 'Watchdog.py: for help type use -h flag'

# read root pf file from command line
if len(sys.argv) > 1:
	root = sys.argv[-1]
else:
	print 'Error: no arguments provided. Exiting.'
	exit (0)

help(root)


# set defaults
nprocessors = 1
diagfolder_name = 'diag_' + root


# check for args
for i in range(len(sys.argv)):
	if sys.argv[i]=="-f":
		diagfolder_name = sys.argv[i+1]


# list files in folder, put in array 'files'
command = "ls " + diagfolder_name + "/" + root + "_*.diag &"
output = subprocess.check_output ( command , shell=True )
files = output.split("\n")
if files[-1] == "": files = files[0:-1]


# number of processors and therefore diag files provided as argument
nfiles = len(files)
nprocessors = nfiles



# create some arrays to store the error logs and increment the counters
# note that we also count the number of threads an error occurs in
error_log = []
error_count = []
thread_count = []



#everything now initialised, can collate errors

for i in range(nfiles):

	filepath = files [i]

	# grep for error summary in the diag file
	command = "grep -A 100 'Error summary: End of program' " + filepath +"&"
	output = subprocess.check_output ( command , shell=True )

	# split the grep output up into lines
	data = output.split("\n")

	# now cycle through the lines in output
	for j in range(len(data)):
		No_log = True	# variable which is True if this is a new error

		# check we have a non blank entry to avoid indexing error
		if len(data[j])>0:
			error = data[j].split()

			# this next line is to check we have an actual error summary and not a string line
			if error[0].isdigit():

				# now get the err log string, and its count
				error = data[j].split()
				
				count = int(error[0])
				log_string = data[j][13:]

				# search for a match in the already fpund errors
				for k in range(len(error_log)):
					if log_string == error_log[k]:

						No_log = False	# we've found a match, increment the counts
						error_count[k] += count
						thread_count[k] += 1
				
				# this is a new error, so append to the log and count arrays
				if No_log:
					error_log.append ( log_string )
					error_count.append ( count )	
					thread_count.append (1)


n_logs=len(error_log)

# print the output to screen
print "Error summary: "+root
print "Collated errors for " + str(nprocessors) + " processors"
print "Recurrences --  Description -- number of threads in which error occurred"
for i in range(n_logs):
	print "\t%d -- %s -- %d" % (error_count[i], error_log[i], thread_count[i])

# all done.
				


