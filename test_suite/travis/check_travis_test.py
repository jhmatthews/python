'''
run simple test of python- check it runs and read sin atomic data
'''

import os
import numpy as np 


f = open("test2.out")

nerror = 0

for line in f:
	data = line.split()

	if len(data) > 0:
		if data[0] == "Get_atomic_data:":
			if data[1] == "Could" and data[2] == "not":
				nerror += 1
				print "Error: %s" % line

		elif data[0] == "Error":
			nerror += 1
			print "Error: %s" % line
	


if nerror > 0:
	print "Errors when running"
	exit(-1)

else:
	print "Ran python OK"
	exit(0)		# test was success