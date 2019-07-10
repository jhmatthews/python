# This is the makefile for the python related programs
#
# usage      make [CC=...] [D] python
# 			 make [CC=...] [D] all to remake all of the routines, and assure they
# 			 					are properly indented.
#
# Adding D causes the routine to be run in a way that profiling and ddd can be used.
# Otherwise the run will be optimized to run as fast as possible. CC is an option to choose
# a different compiler other than mpicc.
#


#MPICC is now default compiler
CC = mpicc
# CC = gcc	can use GCC either from command line or by uncommenting this
FC = g77
# FC = gfortran
# speciify any extra compiler flags here
EXTRA_FLAGS =
LDFLAGS =


# Check a load of compiler options
# This is mostly to address GitHub issues #100
ifeq (mpicc, $(CC))
	# you're using mpicc, so we'll compile with the DMPI_ON flag
	MPI_FLAG = -DMPI_ON

	# check what underlying compiler mpi is using, and the version
	# we'll use this to print to the user and also to do some checks
	MPI_COMPILER := $(shell mpicc --showme:command)
	COMPILER_VERSION := $(shell expr `'$(CC)' -dumpversion`)
	COMPILER_PRINT_STRING = Compiling with mpicc wrapper, for compiler $(MPI_COMPILER) $(COMPILER_VERSION)

	# if it's gcc we want to check if the version is 4.8 or later
	# if it is then we'll disable aggressive loop optimizations, see #100
	ifeq (gcc, $(MPI_COMPILER))
		GCCVERSIONGTEQ4 := $(shell expr `gcc -dumpversion | cut -f1-2 -d.` \>= 4.8)
		ifeq ("$(GCCVERSIONGTEQ4)", "1")
			EXTRA_FLAGS += -fno-aggressive-loop-optimizations 	# add the flag to EXTRA_FLAGS
			COMPILER_PRINT_STRING += with -fno-aggressive-loop-optimizations
		endif
	endif

else ifeq (gcc, $(CC))
	# no mpicc,
	MPI_FLAG =

	# check the version
	# we'll use this to print to the user and also to do some checks
	COMPILER_VERSION := $(shell expr `gcc -dumpversion`)

	# if it's gcc we want to check if the version is 4.8 or later
	# if it is then we'll disable aggressive loop optimizations, see #100
	GCCVERSIONGTEQ4 := $(shell expr `gcc -dumpversion | cut -f1-2 -d.` \>= 4.8)
	COMPILER_PRINT_STRING = Compiling with $(CC) $(COMPILER_VERSION)
	ifeq ("$(GCCVERSIONGTEQ4)", "1")
		EXTRA_FLAGS += -fno-aggressive-loop-optimizations
		COMPILER_PRINT_STRING += with -fno-aggressive-loop-optimizations
	endif
else	# you must be using clang or icc
	MPI_FLAG =

	# check the version we'll use this to print to the user
	COMPILER_VERSION = $(shell expr `$(CC) -dumpversion`)
	COMPILER_PRINT_STRING = Compiling with $(CC) $(COMPILER_VERSION)
endif

# this command finds out how many files with uncommitted changes there are
GIT_DIFF_STATUS := $(shell expr `git status --porcelain 2>/dev/null| grep "^ M" | wc -l`)
GIT_COMMIT_HASH := $(shell expr `git rev-parse HEAD`)



INCLUDE = ../include
LIB = ../lib
BIN = ../bin

ifeq (D,$(firstword $(MAKECMDGOALS)))
# use pg when you want to use gprof the profiler
# to use profiler make with arguments "make D python"
# this can be altered to whatever is best
	CFLAGS = -g -pg -Wall $(EXTRA_FLAGS) -I$(INCLUDE) $(MPI_FLAG)
	FFLAGS = -g -pg
	PRINT_VAR = DEBUGGING, -g -pg -Wall flags
else
# Use this for large runs
	CFLAGS = -O3 -Wall $(EXTRA_FLAGS) -I$(INCLUDE) $(MPI_FLAG)
	FFLAGS =
	PRINT_VAR = LARGE RUNS, -03 -Wall flags
endif



# next line for debugging when concerned about memory problems and duma installed in python directory
# LDFLAGS= -L$(LIB)  -lm -lkpar  -lgsl -lgslcblas ../duma_2_5_3/libduma.a -lpthread
# next line if you want to use kpar as a library, rather than as source below
# LDFLAGS= -L$(LIB)  -lm -lkpar -lcfitsio -lgsl -lgslcblas
LDFLAGS+= -L$(LIB) -lm -lgsl -lgslcblas

#Note that version should be a single string without spaces.

VERSION = 83b

startup:
	@echo $(COMPILER_PRINT_STRING)			# prints out compiler information
	@echo 'YOU ARE COMPILING FOR' $(PRINT_VAR)	# tells user if compiling for optimized or debug
	@echo 'MPI_FLAG=' $(MPI_FLAG)
	echo "#define VERSION " \"$(VERSION)\" > version.h
	echo "#define GIT_COMMIT_HASH" \"$(GIT_COMMIT_HASH)\" >> version.h
	echo "#define GIT_DIFF_STATUS" $(GIT_DIFF_STATUS)\ >> version.h

# indent:
# 	@echo $(INDENT_STRING)
# 	$(INDENT_CMD)


foo: foo.o signal.o time.o
	$(CC) ${cfllags} foo.o signal.o time.o ${LDFLAGS}  -o foo


# these are the objects required for compiltion of python
# note that the kpar_source is now separate from this
python_objects = bb.o get_atomicdata.o photon2d.o photon_gen.o parse.o setup_files.o \
		saha.o spectra.o wind2d.o wind.o  vvector.o recipes.o \
		trans_phot.o phot_util.o resonate.o radiation.o \
		wind_updates2d.o windsave.o extract.o cdf.o roche.o random.o \
		stellar_wind.o homologous.o hydro_import.o corona.o knigge.o  disk.o\
		lines.o  continuum.o get_models.o emission.o cooling.o recomb.o diag.o \
		sv.o ionization.o  levels.o gradv.o reposition.o \
		anisowind.o wind_util.o density.o  bands.o time.o \
		matom.o estimators.o wind_sum.o yso.o cylindrical.o rtheta.o spherical.o  \
		cylind_var.o bilinear.o gridwind.o partition.o signal.o  \
		agn.o shell_wind.o compton.o zeta.o dielectronic.o \
		spectral_estimators.o matom_diag.o \
		xlog.o rdpar.o direct_ion.o pi_rates.o matrix_ion.o para_update.o \
		setup_star_bh.o setup_domains.o setup_disk.o photo_gen_matom.o macro_gov.o windsave2table_sub.o \
		import.o import_spherical.o import_cylindrical.o import_rtheta.o  \
		reverb.o paths.o setup.o run.o brem.o synonyms.o \
		setup_reverb.o setup_line_transfer.o rdpar_init.o cv.o poltest.o


# For reasons that are unclear to me.  get_models.c cannot be included in the sources
# Problems ocurr due to the prototypes that are generated.  ksl 160705
python_source= bb.c get_atomicdata.c python.c photon2d.c photon_gen.c parse.c \
		saha.c spectra.c wind2d.c wind.c  vvector.c recipes.c \
		trans_phot.c phot_util.c resonate.c radiation.c setup_files.c \
		wind_updates2d.c windsave.c extract.c cdf.c roche.c random.c \
		stellar_wind.c homologous.c hydro_import.c corona.c knigge.c  disk.c\
		lines.c  continuum.c emission.c cooling.c recomb.c diag.c \
		sv.c ionization.c  levels.c gradv.c reposition.c \
		anisowind.c wind_util.c density.c  bands.c time.c \
		matom.c estimators.c wind_sum.c yso.c cylindrical.c rtheta.c spherical.c  \
		cylind_var.c bilinear.c gridwind.c partition.c signal.c  \
		agn.c shell_wind.c compton.c zeta.c dielectronic.c \
		spectral_estimators.c matom_diag.c \
		direct_ion.c pi_rates.c matrix_ion.c para_update.c setup_star_bh.c setup_domains.c \
		setup_disk.c photo_gen_matom.c macro_gov.c windsave2table_sub.c \
		import.c import_spherical.c import_cylindrical.c import_rtheta.c\
		reverb.c paths.c setup.c run.c brem.c synonyms.c \
		setup_reverb.c setup_line_transfer.c cv.c poltest.c
#
# kpar_source is now declared separately from python_source so that the file log.h
# can be made using cproto
kpar_source = rdpar.c xlog.c synonyms.c

additional_py_wind_source = py_wind_sub.c py_wind_ion.c py_wind_write.c py_wind_macro.c py_wind.c windsave2table.c windsave2table_sub.c

prototypes:
	cp templates.h templates.h.old
	cproto -I$(INCLUDE) $(python_source) ${additional_py_wind_source} > foo.h
	cp foo.h templates.h
	rm foo.h
	cproto -I$(INCLUDE) $(kpar_source) > log.h
	cproto  -I$(INCLUDE) get_atomicdata.c > atomic_proto.h

python: startup  python.o $(python_objects)
	$(CC) ${CFLAGS} python.o $(python_objects) $(kpar_objects) $(LDFLAGS) -o python
	cp $@ $(BIN)/py
	mv $@ $(BIN)/py$(VERSION)

#This line is jsut so you can use make D python for debugging
D:
	@echo 'Debugging Mode'

py_wind_objects = py_wind.o get_atomicdata.o py_wind_sub.o windsave.o py_wind_ion.o \
		emission.o recomb.o wind_util.o  \
		cdf.o random.o recipes.o saha.o \
		stellar_wind.o homologous.o sv.o hydro_import.o corona.o knigge.o  disk.o\
		lines.o vvector.o wind2d.o wind.o  ionization.o  py_wind_write.o levels.o \
		radiation.o gradv.o phot_util.o anisowind.o resonate.o density.o \
		matom.o estimators.o yso.o photon2d.o cylindrical.o rtheta.o spherical.o \
		import.o import_spherical.o import_cylindrical.o import_rtheta.o \
		cylind_var.o bilinear.o gridwind.o py_wind_macro.o partition.o \
		spectral_estimators.o shell_wind.o compton.o zeta.o dielectronic.o \
		bb.o rdpar.o xlog.o direct_ion.o diag.o matrix_ion.o \
		pi_rates.o photo_gen_matom.o macro_gov.o \
		time.o reverb.o paths.o synonyms.o cooling.o windsave2table_sub.o \
		rdpar_init.o



py_wind: startup $(py_wind_objects)
	$(CC) $(CFLAGS) $(py_wind_objects) $(LDFLAGS) -o py_wind
	cp $@ $(BIN)
	mv $@ $(BIN)/py_wind$(VERSION)



table_objects = windsave2table.o windsave2table_sub.o get_atomicdata.o py_wind_sub.o windsave.o py_wind_ion.o \
		emission.o recomb.o wind_util.o  \
		cdf.o random.o recipes.o saha.o \
		stellar_wind.o homologous.o sv.o hydro_import.o corona.o knigge.o  disk.o\
		lines.o vvector.o wind2d.o wind.o  ionization.o  py_wind_write.o levels.o \
		radiation.o gradv.o phot_util.o anisowind.o resonate.o density.o \
		matom.o estimators.o yso.o photon2d.o cylindrical.o rtheta.o spherical.o \
		import.o import_spherical.o import_cylindrical.o import_rtheta.o  \
		cylind_var.o bilinear.o gridwind.o py_wind_macro.o partition.o \
		spectral_estimators.o shell_wind.o compton.o zeta.o dielectronic.o \
		bb.o rdpar.o rdpar_init.o xlog.o direct_ion.o diag.o matrix_ion.o \
		pi_rates.o photo_gen_matom.o macro_gov.o reverb.o paths.o time.o synonyms.o \
		cooling.o



windsave2table: $(table_objects)
	$(CC) $(CFLAGS) $(table_objects) $(LDFLAGS) -o windsave2table
	cp $@ $(BIN)
	mv $@ $(BIN)/windsave2table$(VERSION)

run_indent:
	../py_progs/run_indent.py -all


# The next line runs recompiles all of the routines after first cleaning the directory
all: clean run_indent python windsave2table py_wind


FILE = get_atomicdata.o atomic.o

libatomic.a:  get_atomicdata.o atomic.o
	ar ru libatomic.a get_atomicdata.o atomic.o
	ranlib libatomic.a
	mv libatomic.a $(LIB)
	cp atomic.h  $(INCLUDE)



clean :
	rm -f *.o  *~
