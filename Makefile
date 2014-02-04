#
# This is the makefile for the python related programs
#
# Make sure you have defined $PYTHON as this directory first
# then type 
# 	make [GIT=False] install
#
# if this fails, consult https://github.com/agnwinds/python/wiki/Installing-Python
# for manual install.

CFITSIO = $(PYTHON)/software/cfitsio3040

GSL = $(PYTHON)/software/gsl-1.15

GIT = True

ifeq (True, $(GIT))
	CLONE_RELEASE = cd $(PYTHON)/progs; git clone https://github.com/agnwinds/python.git; cd python; make python
	PRINT CLONE = 'Cloning Git Release'
else
	CLONE_RELEASE = 
	PRINT_CLONE = 'No git installed- have to obtain release manually from releases page on github. Exiting.'
endif

install:
	@echo 'Installing Python. the radiative transfer code'
	@echo 'Installing in directory '$(PYTHON)
	# First make cfitsio library
	@echo 'Installing CFITSIO library...'
	cd $(CFITSIO); ./configure; make; mv libcfitsio.a ../../lib; make clean
	
	# Then make GSL library
	@echo 'Installing GSL library...'
	cd $(GSL); ./configure --disable-shared --prefix=$(PYTHON)/gsl; make; make check 2>&1; make install; make clean
	
	#finally, make the latest release
	@echo 'Making latest release...'
	@echo $(PRINT_CLONE)
	$(CLONE_RELEASE)

clean: 
	rm -f *.o *~
	cd $(CFITSIO); make clean
	cd $(GSL); make clean 


test:
	mkdir $(PYTHON)/progs
	cd $(PYTHON)/progs
	git clone https://github.com/agnwinds/python.git -b Travis
	cd python
	make test
	make clean
	make py_wind
