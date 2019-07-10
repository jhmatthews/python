#
# This is the makefile for the python related programs
#
# Make sure you have defined $PYTHON as this directory first
# then type 
# 	make [GIT=False] install
#
# if this fails, consult https://github.com/agnwinds/python/wiki/Installing-Python
# for manual install.

CMAKE = mpicc
GSL = $(PYTHON)/software/gsl-2.5
GIT = True
LIBS = True


ifeq (True, $(GIT))
	CLONE_DATA =  mkdir $(PYTHON)/data; cd $(PYTHON)/data; git clone https://github.com/agnwinds/python.git -b data .
	CLONE_RELEASE = mkdir $(PYTHON)/progs; cd $(PYTHON)/progs; git clone https://github.com/agnwinds/python.git; cd python; make CC=$(CMAKE) python; make CC=$(CMAKE) py_wind; make CC=$(CMAKE) windsave2table
	PRINT_CLONE = 'Cloning Git Release'
else
	CLONE_DATA = 
	CLONE_RELEASE = 
	PRINT_CLONE = 'No git installed- have to obtain release and data manually from releases page on github. Exiting.'
endif

MAKE_SOURCE = cd $(PYTHON)/source; make CC=$(CMAKE) python; make CC=$(CMAKE) py_wind; make CC=$(CMAKE) windsave2table



ifeq (True, $(LIBS))
	INSTALL_GSL = cd $(GSL); ./configure --disable-shared --prefix=$(GSL) CC=gcc CCP=ccp; make -i; make check 2>&1; make -i install; make clean; 
	MOVE_GSL = mkdir $(PYTHON)/include/gsl/; mv $(GSL)/include/gsl/* $(PYTHON)/include/gsl; mv $(GSL)/lib/lib* $(PYTHON)/lib/;
else
	INSTALL_GSL = 
	MOVE_GSL =
endif

install:
	@echo 'Installing Python. the radiative transfer code'
	@echo 'Installing in directory '$(PYTHON)
	
	# Then make GSL library
	@echo 'Installing GSL library...'
	$(INSTALL_GSL)
	$(MOVE_GSL)

	
	#finally, make the latest release
	@echo 'Making source code...'
	$(MAKE_SOURCE)


	@echo 'all done'

clean: 
	rm -f *.o *~
	cd $(GSL); make clean;  
	cd $(PYTHON)/source/; make clean

# runs a more rigorous clean in gsl
distclean:
	rm -f *.o *~
	cd $(GSL); make distclean 
	cd $(PYTHON)/source/; make clean

# actually removes the compiled libraries
rm_lib:
	rm -rf $(PYTHON)/include/gsl/
	rm -f $(PYTHON)/lib/lib*
