#!/usr/bin/env python   

'''
This python script updates a .pf file for a new version of
a program with new keywords.  It inserst new keywords; it
eliminates keywords that are not in the base.pf file.  It
will not handle the situtation where the order of keywords
has been changed (because python and other programs often
do not necessarily have unique keywords.  If a new keyword
is in the same relative position in the new file as the old
The old value of the keyword will be maintained.
History
07jul	ksl	Coded
07sep	ksl	Updated to make sure a line in a slave
		file is not written more than once and
		to allow a master file line to override
		a alave line, if the 3rd word is '!!!'
'''

import sys

def process_a_file(master,slave):
	"""Update a single pf file using the master as a template"""
	try: 
		s=open(slave.rstrip(),'r')
	except IOError:
		return(0)

	print "Processing ",slave.rstrip()
	newfilename='New_'+slave.rstrip()
	sout=open(newfilename,'w')
	slines=""
	slines=s.readlines()
	# print "master: ", master[0]
	# print "slave:  ",slines[0]
	i=0
	while i < len(master) :
		j=0
		new=master[i].split() 
		while len(new) > 2 and new[2] == '!!!' :
			sout.write(master[i])
			i=i+1
			new=master[i].split()
		new_keyword=new[0]
		kk=new_keyword.find('(')
		if kk != -1:
			new_keyword=new_keyword[0:kk]
		# print 'Searching for ',new_keyword
		itest=-1
		while j < len(slines):
			old=slines[j].split()
			if len(old)==0:
				j=j+1
				break
			old_keyword=old[0]
			kk=old_keyword.find('(')
			if kk != -1:
				old_keyword=old_keyword[0:kk]
			# print 'Trying ',old_keyword
			if new_keyword == old_keyword:
				# print "Matched ", j, " to ", i
				itest=j
				jhold=itest
				break
			j=j+1
		if itest==-1:
			sout.write(master[i])
		elif len(old) > 1 :
			# sout.write(slines[itest])
			sout.write('%s	%s\n'%(new[0],old[1]))
			slines[j]='Used up'
		i=i+1
	s.close()
	return(0)

if __name__ == "__main__":		# allows one to run from command line without running automatically with write_docs.py

	argc=len(sys.argv)

	if(argc==1):
		print 'Interactive inputs'
		master_pf = raw_input("Master parameter file:  ")
		list_to_update = raw_input("List.to.update:  ")
	elif argc==3:
		print 'Getting inputs from command line'
		master_pf=sys.argv[1]
		list_to_update=sys.argv[2]
	else:
		print 'Usage is either via command line or intractively'
		print 'If interactive, simply enter pj_update.py, and answer questions'
		print 'If command line, enter pf_update.py master.pf list.of.old.pf.files'
		sys.exit(0)
		

	print 'Using ', master_pf, ' as basis for updates'
	print 'Using', list_to_update, ' as list of files to update'

	try:
		master=open(master_pf,'r')
	except IOError:
		print "File ",master_pf,"does not exist"
		sys.exit(0)


	try:
		list=open(list_to_update,'r')
	except IOError:
		print "File ",list_to_update,"does not exist"
		sys.exit(0)

	mlines=master.readlines()
	lines=list.readlines()

	for line in lines:
		process_a_file(mlines,line)




