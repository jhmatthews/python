#!/usr/bin/env python
'''
090301	ksl	The main routine is apparently explore, which allows
one to compare models to data.

Bugs - The problem aborst when one selects and inclination that is not
there.  It's also not obvious from the what the inputs are and particularly
how columns and inclinations are numbered.

There is some evident overlap with monte.

090301	ksl	Made some cosmetic changes to make easier to use.

'''

import numpy
import pylab
from matplotlib.ticker import FuncFormatter
from kslio import *
import carlo



def get_python_spec_angles(filename='py_ixvel/ixvel_00_00_00_00_00_00_00.spec'):
	'''
	This routine simply finds the angles for which a specific python spectrum
	has been calculated and returns them as a list
	'''
	try:
		specfile=open(filename,'r')
	except  IOError:
		print filename, " does not exist"
		return w,f

	angles=[]
	line=specfile.readline()
	while line!='':
		z=line.split('Ang:')
		if len(z)>1:
			i=1
			while i<len(z):
				angles=angles+[eval(z[i])]
				i=i+1

			return angles
		line=specfile.readline()
	return angles


def read_python_spec(filename,column):
	"""
	Read the a spectrum file created with python

	columns are properly numbered in the sense that
	column 1 refers to the first real spectrum

	The routine returns an empty spectrum on any
	error

	"""
	w=[]
	f=[]

	try:
		specfile=open(filename,'r')
	except  IOError:
		print filename, " does not exist"
		return w,f
	line=specfile.readline()
	while line!='':
		z=line.split()
		word1=z[0]
		if word1[0]!='#':
			w=w+[float(z[1])]
			try:
				f=f+[float(z[column+7])]
			except IndexError:
				print 'Spec %s does not exist in %s' % (column,filename)
				return w,f
		else:
			gotangle='no'
			if word1=='#Freq.':
				try:
					angle=z[column+7]
				except IndexError:
					print 'Angle %s does not exist in %s' % (column,filename)
					return w,f
				gotangle='yes'
			elif z[1]=='Freq.':
				try:
					angle=z[column+8]
				except IndexError:
					print 'Angle %s does not exist in %s' % (column,filename)
					return w,f
				gotangle='yes'
			else:
				print z
			if gotangle=='yes':
				print 'The angle was',angle
		line=specfile.readline()
	w=numpy.array(w)
	f=numpy.array(f)
	return w, f


def read_data_spec(filename,col_wave=0,col_flux=1,col_err=2):
	"""
	Read the ascii data spectrum file with two
	or three columns.  If there are three columnts
	it is assumed that the third column is the
	error.

	090306	ksl	Generalized so one could parse most
			spectral files
	"""
	wave=[]
	flux=[]
	error=[]

	try:
		specfile=open(filename,'r')
	except  IOError:
		print filename," does not exist"
		return wave,flux,error
	line=specfile.readline()
	while line != "":
		z=line.split()
		word1=z[0]
		if word1[0] != '#':
			wave=wave+[float(z[col_wave])]
			flux=flux+[float(z[col_flux])]
			if len(z) > col_err:
				error=error+[float(z[col_err])]
			else:
				error=error([1.0])
		line=specfile.readline()
	wave=numpy.array(wave)
	flux=numpy.array(flux)
	error=numpy.array(error)
	return wave, flux, error


def reduce_xy(x,y,xmin,xmax):
	'''
	Get the portions of the array that are within xmin and xmax
	'''
	condition= (xmin < x) & (x < xmax)
	xout=numpy.compress(condition,x)
	yout=numpy.compress(condition,y)
	return xout,yout


def ave (x,y,xmin,xmax):
	"""average a restricted portion of y between xmin and xmax"""
	xout,yout=reduce_xy(x,y,xmin,xmax)
	a=numpy.average(yout)
	return a

def rescale (zin,zout,scale):
	"""rescale a list by a scalefactor"""
	lin=len(zin)
	i=0
	while i < lin:
		zout.extend([scale*zin[i]])
		i=i+1
#        print 'len zout', len(zout)
	return len(zout)


def scale_model (wdat,fdat,wmod,fmod,wmin,wmax):
	"""Produce a rescaled model"""
	data_ave=ave(wdat,fdat,wmin,wmax)
	mod_ave=ave(wmod,fmod,wmin,wmax);
	scale=data_ave/mod_ave
	return scale*fmod



def get_data_model (data,model,ang_no):
	'''
	get_data_model just gets an observed
	spectrum, data, and a model spectrum

	data is the name of the observed spectrum
	mdoel is the filename of the model
	ang_no is the spectrum number beginning
		with 1

	Note - It's not clear that this routine
	is that useful.

	'''
	w,f,e=read_data_spec(data)
	wmod,fmod=read_python_spec(model,ang_no)
	return w,f,wmod,fmod

def plot_data_model(data,model,column,wmin,wmax):
	"""
	plot a rescaled model against the data where
	data is a standard ascii spectrum file and
	model is a python output spectrum and
	wmin and wmax are rhe regions to be included int he plot
	"""


	w,f,wmod,fmod=get_data_model(data,model,column)

	# Check that get_data_model succeeded
	if len(f)==0 or len(fmod) == 0:
		print 'data (%s) or mod (%s col  %d) were not found' % (data,model,column)
		return w,f,wmod,fmod

	ww,ff=reduce_xy(w,f,wmin,wmax)
	wwmod,ffmod=reduce_xy(wmod,fmod,wmin,wmax)


	fscale=scale_model (ww,ff,wwmod,ffmod,wmin,wmax)

	pylab.plot(wwmod,fscale,'r')
	pylab.plot(ww,ff,'b')

	reset=pylab.axis()
	xreset=list(reset)
	xreset[0]=wmin
	xreset[1]=wmax
	xreset[2]=0
	y,ylabel,yname=nice_flux_labels(xreset[3])
	pylab.yticks(y,ylabel)
	pylab.axis(xreset)
	pylab.draw()


	return w,f,wmod,fmod

def add_plot(w,f,wmod,fmod,wmin,wmax):
	'''
	This is the generic routine that plots a panel of the spectrum
	'''

	# Check that get_data_model succeeded
	if len(f)==0 or len(fmod) == 0:
		print 'data (%d) or mod (%d ) had zero length' % (len(f),len(fmod))
		return

	ww,ff=reduce_xy(w,f,wmin,wmax)
	wwmod,ffmod=reduce_xy(wmod,fmod,wmin,wmax)


	fscale=scale_model (ww,ff,wwmod,ffmod,wmin,wmax)

	pylab.plot(wwmod,fscale,'r')
	pylab.plot(ww,ff,'b')

	reset=pylab.axis()
	xreset=list(reset)
	xreset[0]=wmin
	xreset[1]=wmax
	xreset[2]=0
	xmed=numpy.median(ff)
	xreset[3]=3.*xmed
	y,ylabel,yname=nice_flux_labels(xreset[3])
	pylab.yticks(y,ylabel)
	pylab.axis(xreset)
	pylab.draw()
	return


def xadd_plot(w,f,wmod,fmod,band):
	'''
	This is the generic routine that plots a panel of the spectrum
	'''

	# Check that get_data_model succeeded
	if len(f)==0 or len(fmod) == 0:
		print 'data (%d) or mod (%d ) had zero length' % (len(f),len(fmod))
		return

	wmin=band[0]
	wmax=band[1]

	ww,ff=reduce_xy(w,f,wmin,wmax)
	wwmod,ffmod=reduce_xy(wmod,fmod,wmin,wmax)


	fscale=scale_model (ww,ff,wwmod,ffmod,wmin,wmax)

	pylab.plot(wwmod,fscale,'r')
	pylab.plot(ww,ff,'b')

	reset=pylab.axis()
	xreset=list(reset)
	xreset[0]=wmin
	xreset[1]=wmax
	xreset[2]=0
	xmed=numpy.median(ff)
	xreset[3]=3.*xmed
	y,ylabel,yname=nice_flux_labels(xreset[3])
	pylab.yticks(y,ylabel)
	print "The yname was %s ",yname

	if band[2]!='Complete':
		wmid=0.5*(wmax+wmin)
		ww=(int(wmid)/5)*5
		pylab.xticks([ww-10,ww+10])
		x=0.5*(wmin+wmid)
		y=0.8*xreset[3]
		pylab.text(x,y,band[2])


	pylab.axis(xreset)
	return yname







def xplot(data,model,column,inst='hst'):
	"""
	plot a rescaled model on an instrument
	by intrument basis

	090307 	This is a generalized versons of
	the earlier routine to reflect
	what I have learned earlier
	"""
	pylab.figure(2,figsize=(14,10))
	pylab.clf()


	hst_bands=          [[1155, 1195,'C III']]
	hst_bands=hst_bands+[[1220, 1260,'N V']]
	# 1393.7546  1402.7697
	hst_bands=hst_bands+[[1380, 1420,'Si IV']]
	hst_bands=hst_bands+[[1530, 1570,'C IV']]
	hst_bands=hst_bands+[[1170, 1650,'Complete']]


	# 1031.93 )triplet) 1037.62 (singlet)
	hut_bands=          [[1010,1050,'O VI']]
	hut_bands=hut_bands+[[1220,1260,'N V']]
	# 1393.7546  1402.7697
	hut_bands=hut_bands+[[1380, 1420,'Si IV']]
	hut_bands=hut_bands+[[1530, 1570,'C IV']]
	hut_bands=hut_bands+[[ 900, 1825,'Complete']]


	# SVI 933.378 (triplet) 944.523 singlet
	fuse_bands=          [[920,960,'S VI']]
	# CIII 977.0201
	# NIII 991.577 991.511 are slighly off the ground state, but have highes progs
	# NIII 989.799 is the ground state, but has a lower gf.
	fuse_bands=fuse_bands+[[960,1000,'CIII 977']]
	fuse_bands=fuse_bands+[[1010,1050,'O VI']]
	fuse_bands=fuse_bands+[[1140,1180,'C III 1175']]
	fuse_bands=fuse_bands+[[900,1180,'Complete']]



	if inst=='hut':
		band=hut_bands
		mytitle='HUT'
	elif inst=='fuse':
		band=fuse_bands
		mytitle='FUSE'
	else:
		mytitle='HST'
		band=hst_bands


	print band

	# Get the data and the model
	w,f,e=read_data_spec(data)
	wmod,fmod=read_python_spec(model,column)

	pylab.subplot(241)
	xadd_plot(w,f,wmod,fmod,band[0])
	# pylab.xticks([1180,1200])

	pylab.subplot(242)
	xadd_plot(w,f,wmod,fmod,band[1])
	# pylab.xticks([1230,1250])

	pylab.subplot(243)
	xadd_plot(w,f,wmod,fmod,band[2])
	# pylab.xticks([1380,1400])
	pylab.subplot(244)

	xadd_plot(w,f,wmod,fmod,band[3])
	# pylab.xticks([1540,1560])

	# The entire spectrum

	pylab.subplot(212)
	yname=xadd_plot(w,f,wmod,fmod,band[4])

	# Finally add the global labels
	pylab.figtext(0.5, 0.94,mytitle,horizontalalignment='center')
	pylab.figtext(0.5, 0.0,r'Wavelength ($\AA$) ',horizontalalignment='center')

	string=r'Flux ($ 10^{%s} ergs \/ cm^{-2} s^{-1} \AA^{-1}$)' %yname

	pylab.figtext(0.03,0.5,string,verticalalignment='center',rotation='vertical')

	pylab.draw()
	pylab.savefig(inst+'.jpg')




def get_grid(filename='sscyg_kgrid0902/Models.ls',outfile='none'):
	'''
	get_grid reads a file, containing the list of models in the format produced
	by gen_grid.py and read by pyfit3.

	It returns
		columns which is the list of the column names and
		files,which  is a list of all the files, with their associated values

	If outfile is something other than none, all of the comments are written to
	that file.  This is to prepare for writting a new file with censored values,
	see sensor below.


	Note that this does not include anything having to do with the angles for which
	spectra were generated.
	'''

	f=open(filename,'r')

	print 'get_grid ',outfile
	if outfile!='none':
		g=open(outfile,'w')



	line=f.readline()
	columns=[]
	files=[]
	n=0
	while line!='':
		z=line.split()
		if z[0]=='#':
			if z[1]=='Variable':
				columns=columns+[z[2]]
			if outfile !='none':
				print line
				g.write('%s' % line)
		else:
			nvar=len(columns)
			zz=[z[0]]
			m=1
			while m<len(z):
				zz=zz+[float(z[m])]
				m=m+1
			files=files+[zz]

		line=f.readline()
	f.close()
	if outfile !='none':
		g.close()
	return columns,files



def get_unique(mylist):
	'''
	get_unique simply return the unique values in
	a list
	'''
	myset=set(mylist)
	xlist=list(myset)
	xlist.sort()
	finallist=xlist
	return finallist

def get_chosen_model(files,possible,choices):
	'''
	get the name of a model that has the specific
	values of mdot, etc. that were requested.

	files is the list of all of the models with the
		variables that change in the grid.

	xchoices is an array contining the specific
	values of the variables that have been chosen.

	090301	ksl	Added better documentation and
			slightly modified routine
	090302	ksl	Modified to locate which particular
			spectrum as well
	'''
	m=0;
	chosen='none'
	while m<len(files):
		n=0
		good=1
		#DEBUG print files[m]
		while n<len(choices)-1:
			x=possible[n][choices[n]]
			if files[m][n+1]!=x:
				#DEBUG print "compare",files[m][n+1],x
				good=0
			n=n+1
		if good==1:
			# Have found the model
			chosen=files[m][0]
			break
		m=m+1

	print 'get_chosen: ',chosen

	# Now the spectrum is given trivially, since we are using
	# the number to indicate what spectrum we want

	return chosen,choices[n]

def get_chosen_models(files,possible,choices):
	'''
	get_chosen_models returns a new list of files that have
	been censored by choices
	'''

	print 'files :', len(files),files[0]
	print 'possible ',len(possible),possible[0]
	print 'choices: ',len(choices), choices

	good_files=[]
	m=0;
	while m<len(files):
		n=0
		good='ok'
		while n<len(choices):
			this_choice=choices[n]
			# print n,this_choice
			if this_choice[0]==-1:
				# print 'caught -1'
				n=n+1
				continue
			else:
				k=0
				good='nok'
				while k<len(this_choice):
					x=possible[n][this_choice[k]]
					if files[m][n+1]==x:
						good='ok'
						break
					k=k+1
			if good=='nok':
				break
			n=n+1


		if good=='ok':
			good_files=good_files+[files[m]]
		# print m
		m=m+1
	print 'Found %d good files ' % len(good_files)
	return good_files

def look(modellist='sscyg_kgrid0902/Models.ls',outfile='none'):
	'''
	look produces a list of the models and their parameters in
	a format that one can select individual models.

	The routine returns:
		wholestring  -- a formatted version of the unique values for each variable
		possible -- a list in which each row has the unique values for a variable
		files -- the list of files with all the associated variables for that file

	090302	ksl	Added the angles to possible and wholestring for a more consistend
			interface
	090312	ksl	Added outfile.  This really is just a pass through to
			get_grid.  It's possible that get_grid should be
			pulled out of this routine

	'''
	# Read the model list file
	columns,files=get_grid(modellist,outfile)

	# The next section looks for the values that exist of the
	# paremeters that vary an new list called possible that
	# contains the allowed values of each row.
	ncols=len(columns)
	n=1
	possible=[]
	while n<=ncols:
		m=0
		mylist=[]
		while m<len(files):
			mylist=mylist+[float(files[m][n])]
			m=m+1
		newlist=get_unique(mylist)
		row=[columns[n-1]]+newlist
		possible=possible+[row]
		n=n+1



	# Next    get the angles contained in the first python spectrum
	angles=get_python_spec_angles(filename=files[0][0])
	print 'Got the angles', angles

	columns=columns+['Angles']
	row=['Angles']
	i=0
	while i<len(angles):
		row=row+[angles[i]]
		i=i+1
	possible=possible+[row]
	# possible=possible+[angles]

	# The next section just creates a formatted version
	# of the possible values.  It seems to have been
	# done here for convenience more than anything.
	n=0

	wholestring=[]
	while n<len(possible):
		row=possible[n]
		string='%15s'% row[0]
		m=1
		while m<len(row):
			string=string+ "%8.2g " % row[m]
			m=m+1
		wholestring=wholestring+[string]
		n=n+1

	return wholestring,possible,files


def censor(modellist='py_ixvel/Models.ls',outputfilename='test.ls'):
	'''
	censor take a standard modellist and create a new modellist.

	modellist is the original list
	outputfilename is the name of the new list

	This uses raw input commands, e..g lists

	090312 	ksl	Added

	'''
	wholestring,possible,files=look(modellist,outputfilename)
	# The routine returns:
	# 	wholestring  -- a formatted version of the unique values for each variable
	# 	possible -- a list in which each row has the unique values for a variable
	# 	files -- the list of files with all the associated variables for that file

	n=0
	print 'Enter the integer vallues that you want as a python list'
	choices=[]
	while n<len(wholestring)-1:
		print wholestring[n]
		jj=len(wholestring[n].split())
		string = "Choice (1 through %s) as python list, e,g [1,2]:  " % (jj-1)
		zz=input(string)
		print zz
		if len(zz)==0:
			zz=[-1]
		choices=choices+[zz]
		n=n+1
	print choices

	censored_list=get_chosen_models(files,possible,choices)

	g=open(outputfilename,'a')

	i=0
	while i<len(censored_list):
		line=censored_list[i]
		g.write('%30s' % line[0])
		j=1
		print line
		while j<len(line):
			g.write('%9.3g' % line[j])
			j=j+1
		g.write('\n')
		i=i+1
	g.close()



def explore(modellist='py_ixvel/models.ls',spectrum='ixvel_stis2',obs='hst',ions='no'):
	'''

	The imputs are a list of files produced by the .py routine that can be used
	to generate a grid of models.  This is the same file that is read by
	pyfit3 (in the new format).

	A spectrum in a format similar to that produced by pyfit.3

	and the type of observation 'fuse' or 'hst'

	Bug - 090301 - Although the program runs, it currently produces so many diagnostics
	that one cannot easily see where one was.

	The routine fails if the inclination is not in the inclination range.

	090301	ksl	Added defaults for names, in attempt to recover what these
			routines were supposed to do.
	'''

	wholestring,possible,files=look(modellist)

	# At this point we have everything we need and so can start plotting
	choices=len(possible)*[1]
	xchoices=len(possible)*[99.99]
	ispec=1

	go_on='yes'

	while go_on=='yes':
		n=0
		while n<len(wholestring):
			print wholestring[n]
			jj=len(wholestring[n].split())
			try:
				i=int(get_input(("Choice (1 to %s)" % (jj-1)),str(choices[n])))
				choices[n]=i
				try:
					xchoices[n]=possible[n][choices[n]]
					n=n+1
				except IndexError:
					print 'choice %d out of range' % (i)
			except ValueError:
				print 'Could not parse choice, try again!'
				continue




		print ' Choices: ',choices
		print 'Xchoices: ',xchoices

		modelname,ispec=get_chosen_model(files,possible,choices)

		print spectrum,modelname,ispec


		if modelname!='none':
				xplot(spectrum,modelname,ispec,obs)
		else:
			print 'Model with these values not calculated'

		if ions!='no':
			carlo.cno(modelname[0:modelname.rindex('.spec')])

		go_on=get_input("Continue",go_on)





def do_one(data,model,column,wmin=1525.,wmax=1575.,wvel=1550.):
	"""
	plot a rescale model against hst data

	Note that the wavelengths to be plotted are
	hardcoded
	"""
	pylab.figure(3)
	pylab.draw()
	pylab.clf()
	ax=pylab.subplot(111)



	w,f,wmod,fmod=plot_data_model(data,model,column,wmin,wmax)

	x,xlabel=vel(wmin,wmax,wvel,5000.)
	pylab.xticks(x,xlabel)
	pylab.draw()
#	pylab.close()


	'''
	The purpose of this routine is to generate axis labels in
	velocity space, where the plot that is being labelled
	has tickmarks at waves, and wzero is the zero of velocity
	space.
	'''

	vmin=3e5*(wmin-wzero)/wzero
	vmax=3e5*(wmax-wzero)/wzero
	imin=int(vmin/delta_v)
	imax=int(vmax/delta_v)

	delta_w=delta_v/3e5*wzero

	x=[]
	label=[]
	i=imin
	while i<=imax:
		x=x+[wzero+i*delta_w]
		label=label+['%d' % (delta_v*i)]
		i=i+1

	return x,label


def vel_lim(wave=1550.,max_v=3000,ntics=3):
	x=[]
	label=[]
	dw=wave*max_v/3e5

	wmin=wave-dw
	wmax=wave+dw

	vmin=-max_v
	vmax=max_v

	dw=(wmax-wmin)/ntics
	dv=(vmax-vmin)/ntics


	i=0
	while i<ntics:
		x=x+[wmin+(i+0.5)*dw]
		label=label+['%d' % (vmin+(i+0.5)*dv)]
		i=i+1
	return wmin,wmax,x,label

def nice_flux_labels(ymax,nticks=5):
	'''
	Get some nicer places to put
	the flux labels
	'''


	if ymax>1e-10:
		scale=1e-10
		scaleword='-10'
	elif ymax>1e-11:
		scale=1e-11
		power='-11'
	elif  ymax>1e-12:
		scale=1e-12
		power='-12'
	elif  ymax>1e-13:
		scale=1e-13
		power='-13'
	elif  ymax>1e-14:
		scale=1e-14
		power='-14'
	else:
		scale=1e-15
		power='-15'

	smax=ymax/scale  # This is a number between 1 and 10
	print 'nice_flux', smax

	delta=3
	while 3 > smax/delta:
		delta=delta-0.1

	label=[]
	y=[]
	i=0
	while i<4:
		y=y+[(delta*scale*i)]
		label=label+['%.1f' % (delta*i)]
		i=i+1
	return y,label,power



def pyfit3_out(filename,inst='hst'):
	'''
	Plot the results of pyfit3 in a standard way.  The instrument
	defines what standard wavelength ranges are involved
	'''

	wav,flux,err=read_data_spec(filename,0,1,2)
	wav,mod,good=read_data_spec(filename,0,3,4)

	hst_bands=          [[1175.64,'C III']]
	hst_bands=hst_bands+[[1238.821,'N V']]
	# 1393.7546  1402.7697
	hst_bands=hst_bands+[[1393.754,'Si IV']]
	hst_bands=hst_bands+[[1548.203,'C IV']]


	# 1031.93 )triplet) 1037.62 (singlet)
	hut_bands=          [[1031.93,'O VI']]
	hut_bands=hut_bands+[[1238.821,'N V']]
	# 1393.7546  1402.7697
	hut_bands=hut_bands+[[1393.754,'Si IV']]
	hut_bands=hut_bands+[[1548.203,'C IV']]




	# SVI 933.378 (triplet) 944.523 singlet
	fuse_bands=          [[933.378,'S VI']]
	# CIII 977.0201
	# NIII 991.577 991.511 are slighly off the ground state, but have highes progs
	# NIII 989.799 is the ground state, but has a lower gf.
	fuse_bands=fuse_bands+[[977.02,'CIII 977']]
	fuse_bands=fuse_bands+[[1031.93,'O VI']]
	fuse_bands=fuse_bands+[[1175.64,'C III 1175']]




	if inst=='hut':
		mytitle='HUT'
		bands=hut_bands
		vel_width=9000
	elif inst=='fuse':
		mytitle='FUSE'
		bands=fuse_bands
		vel_width=9000
	else:
		mytitle='HST'
		bands=hst_bands
		vel_width=3000

	pylab.figure(3)
	pylab.clf()


	i=0
	while i<4:

		subplotno=221+i

		ax=pylab.subplot(subplotno)

		# Read the specific spectrum for this file
		wav,line,good=read_data_spec(filename,0,5+2*i,6+2*i)
		# Get the wavelengths for this band
		wmin,wmax,x,xlabel=vel_lim(bands[i][0],vel_width,3)

		# Reduce all three spectra to the wavelength ranges to be plotted
		ww,ff=reduce_xy(wav,flux,wmin,wmax)
		wwmod,ffmod=reduce_xy(wav,mod,wmin,wmax)
		wwbest,ffbest=reduce_xy(wav,line,wmin,wmax)

		# plot the spectra
		pylab.plot(wwmod,ffmod,'g')
		# pylab.plot(wwmod,ffmod,':k',linewidth=3)
		pylab.plot(ww,ff,'b')
		pylab.plot(wwbest,ffbest,'r')



		# Now label it all
		reset=pylab.axis()
		xreset=list(reset)
		xreset[0]=wmin
		xreset[1]=wmax
		xreset[2]=0
		xreset[3]=2*numpy.median(ff)
		pylab.axis(xreset)
		pylab.xticks(x,xlabel)
		y,ylabel,yname=nice_flux_labels(xreset[3])
		pylab.yticks(y,ylabel)
		# Add a label.  Note that .text is in coordinates of plt
		xpos=wmin+(wmax-wmin)*0.1
		ypos=0.8*xreset[3]
		pylab.text(xpos,ypos,bands[i][1])


		i=i+1

	string=r'Flux ($ 10^{%s} ergs \/ cm^{-2} s^{-1} \AA^{-1}$)' %yname



	# Finally add the global labels
	pylab.figtext(0.5, 0.94,mytitle,horizontalalignment='center')
	pylab.figtext(0.5, 0.0,r'Velocity ($km \/ s^{-1}$) ',horizontalalignment='center')
	pylab.figtext(0.03,0.5,string,verticalalignment='center',rotation='vertical')


#	pylab.xticks(x,xlabel)
	pylab.draw()
	pylab.savefig(filename+'.jpg')



print '''
The main routines are:

compare.explore(modellist='py_ixvel/models.ls',spectrum='ixvel_stis2',obs='hst',ions='no')
	allows one to explore models based on the variables that were used to generate the models

comapare.pyfit3_out(filename,inst='hst')
	which produces a plot that shows the results of a pyfit3 run.
	This is rougly equivalent to sm.pyfit3

compare.xplot(data,model,column,inst='hst')
	compare observed spectrum to a specific model and numbered angle on an instrument
	by instrument basis

compare.censor(modellist,newlist)
	produces a new list of files in a format for pyfit3 which have been censored to for
	example be a specific mdot.
For more extensive information, use the help facility
'''
