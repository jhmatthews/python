'''
run simple test of python- check it runs and read sin atomic data
'''


f = open("test1.out")

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
	

f = open("balmer_test.out")

n3 = 1.0
n4 = 1.0

for line in f:
	data = line.split()

	if len(data) > 0:
		if data[0] == "Get_atomic_data:":
			if data[1] == "Could" and data[2] == "not":
				nerror += 1
				print "Error: %s" % line

		#elif data[0] == "Error":
		#	nerror += 1
		#	print "Error: %s" % line


		elif data[0] == "Macro":
			if data[1] == "Atom" and data[2] == "level":
				for i in range(len(data)):
					if data[i] == "n":
						if data[i+1] == "3":
							n3 = float(data[-1])
						if data[i+1] == "4":
							n4 = float(data[-1])


	ratio = n3 / n4 



if ratio > 1.7 and ratio < 2.3:
	print "H alpha / H beta is %8.4e" % ratio
elif ratio != 1.0:
	print "Error, H alpha / H beta is %8.4e" % ratio
	nerrors += 1





if nerror > 0:
	print "Errors when running"
	exit(-1)

else:
	print "Ran python OK"
	exit(0)		# test was success
