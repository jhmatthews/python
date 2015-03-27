#!/usr/bin/env python 

'''
                    Space Telescope Science Institute


Synopsis:  
	This is a simple program which parses a .sig
	file and determines whether it is reasonable
	to restart the program to complete it

Description:  

	restart.py  -n number rootname 

Arguments:  

	None

Returns:
	A new output file 

Notes:
									   
History:
07jan	ksl	Coded


'''

import sys
import os



def parse_sigfile(root):

	sigfile=root+'.sig'

	f=open(sigfile,'r')
	lines=f.readlines()
	f.close()

	i=0
	restarts=0
	while i<len(lines):
		z=lines[i].split()
		if z[0]!='#':
			if(z[6]=='OK'):
				status='OK'
			if(z[6]=='NOK'):
				status='NOK'
			if(z[6]=='COMPLETE'):
				status='COMPLETE'
		if(z[1]=='RESTART'):
			restarts=restarts+1
		i=i+1
	return status,restarts

def append_sigfile(root,command):
	sigfile=root+'.sig'

	f=open(sigfile,'a')
	f.write('# RESTART %s\n' % command)
	f.close()

# Beginning of main routine

argc=len(sys.argv)

if(argc==1):
	print "Usage: restart.py [-n max] [-c commandfile] [-C command] root"
	sys.exit()

maxtimes=10
commandfile='Doit'
command='source '+commandfile

i=1;
while (i<argc):
	if(sys.argv[i]=='-n'):
		maxtimes=int(sys.argv[i+1])
	if(sys.argv[i]=='-c'):
		commandfile=sys.argv[i+1]
		command='source '+commandfile
	if(sys.argv[i]=='-C'):
		command=sys.argv[i+1]
	i=i+1

root=sys.argv[argc-1]

# Get the status of the previous run of the routine
status,restarts=parse_sigfile(root)

if (status=='OK' and restarts <= maxtimes):
	print '!! Restarting: Status %s and restarts %d <= maxstarts %d' % (status,restarts,maxtimes)
	append_sigfile(root,command)
	os.system(command)
else:
	print '!! Not restarting: Status %s or restarts %d > maxstarts %d' % (status,restarts,maxtimes)
sys.exit()
