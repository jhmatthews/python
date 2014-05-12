04feb13 -This is a directory which contains all of the various types
of data that one might want to use.  This includes spectral
models for use with python as well as atomic data and linelists


* atomicxx:  The directories that contain all of the atomic data that python uses directly.
			 highest number is latest version
	* atomic77: The version of atomic data which has XSections extrapolated to higher energies
				in order to remove artificial edges in spectra caused to topbase XSections running out.
				Identical to standard73 otherwise.
	* atomic73: The version of atomic data used for the benchmark model in Higginbottom et al. 2013. 
				masterfile, standard73.
	* atomic39: The version of atomic data used for Long & Knigge 2002
				masterfile, standard39.
	* atomic70: An interim version of atomic data used during AGN development
				masterfile, standard70..
	* atomic_h20: Atomic data for the 20 level hydrogen atom. Used for macro atoms, in particular for the 
				Sim 2005 paper on YSOs. masterfile, h20.
	* atomic_old: a collection of old atomic data, mainly various macro atom data which has not been used
				since python ~58 but may be used for reference and in future.

* calculate_atomic: A directory that contains data that must be reformated and/collated so
	that it can be used by atomic.  For example, the program that create linelists
        for python from the entire Kurucz line list are here.

* hutdisk93: A set of stellar spectra calculated with TLUSTY/Sysnpec that I have used
	extenstively for calculations in the UV

* kurucz: The directories that contain the kurucz models that are used by python.  Actually
	what is used in python is in the data.91 directory.  These are binary versions
	of the Kurucz models.  The fact that the binary versions are used is an anachronism. 	     
	When I began working with these models read/write access was a lot slower than it
	is today.  

* lejeune: A set of stellar models from lejeune.  Note that I have not actually used thse.
	Their advantage is that they go further into the IR than kurucz

History:

04july -- ksl -- Created separate directories ../data44 and ../data50 that just contains
	the data needed to store the atomic files for running Python.  This directory
	contains the old data, as well as data that may not be need for standard runs
	of the program.  
	atomic49 is the version of the atomic datafile used up to and including python49
	atomic50 is a version that strips out a lot of the less used portions of the 
	atomic data.  
	Ultimately one may be able to discard both of these directoreis in favor of the
	datafiles contained in ../data44 and ../data50, but I was a little unsure that
	I had kept everything I needed for them.

1301 -- JM -- The data is now all in one folder (data). atomic data has the form 
	 atomic/atomicxx/standardxx where xx is two numbers.

1307 -- JM -- data folder still the same. standardxx now in top level directory,
	addresses bug #13
