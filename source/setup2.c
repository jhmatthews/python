/***********************************************************
                                       Space Telescope Science Institute

 Synopsis:
 
Arguments:		

Returns:
 
Description:	
	
		
Notes:

History:
	15sep	ksl	Setup and other ancillary routines that were part
			of python.c
**************************************************************/



#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "atomic.h"


#include "python.h"

/***********************************************************
                                       Space Telescope Science Institute

 Synopsis:
	print out some basic help on how to run the program
Arguments:		

Returns:
 
Description:	
		
Notes:

The easiest way to create the message, at least initially, is simply to to type
out what you want to appear on the screen and then as \n\ to all of the lines, including
the ones with nothing in them

History:
	081217	ksl	67c - Added so that ksl could remember all of the options
	09feb	ksl	68b - Added info on -v switch

**************************************************************/

int
help ()
{
  char *some_help;

  some_help = "\
\n\
This program simulates radiative transfer in a (biconical) CV, YSO, quasar or (spherical) stellar wind \n\
\n\
	Usage:  py [-h] [-r] [-t time_max] xxx  or simply py \n\
\n\
	where xxx is the rootname or full name of a parameter file, e. g. test.pf \n\
\n\
	and the switches have the following meanings \n\
\n\
	-h 	to ge this help message \n\
	-r 	restart a run of the progarm reading the file xxx.windsave \n\
	-e change the maximum number of errors before quit- don't do this unless you understand\
	the consequences! \n\
\n\
	-t time_max	limit the total time to approximately time_max seconds.  Note that the program checks \n\
		for this limit somewhat infrequently, usually at the ends of cycles, because it \n\
		is attempting to save the program outputs so that the program can be restarted with \n\
		-r if that is desired. \n\
\n\
	-v n	controls the amount of print out.  The default is 4.  Larger numbers increase  \n\
		the amount printed; smaller numbers decrease it.   \n\
	if one simply types py or pyZZ where ZZ is the version number one is queried for a name \n\
	of the parameter file. \n\
\n\
\n\
";				// End of string to provide one with help

  printf ("%s\n", some_help);

  exit (0);
}

/***********************************************************
                                       Space Telescope Science Institute

 Synopsis:
	init_geo initializes the geo structure to something that is semi-reasonable
Arguments:		

Returns:
 
Description:	
	Initial values for all of the variables that are not part of the individual
	wind descriptions that are actully read(!) into the program should be 
	created here.  The derived values are not needed.

		
Notes:

	When initializing geo, be sure to initialize to cgs units, since this cgs units
	are the working units for the program.  This is necessary for consistncy when
	one tries to restart the program.

	Note that init_geo is set up for CVs and Stars and not AGN

	XXX init_geo ought to be particularized for the system type.


History:
 	98dec	ksl	Coded and debugged.  Much of code was copied from old main routine for
			python
	04dec	ksl	This is probably still not completely up to date, but have
			added some initializations 
	080518	ksl	60a - modified to set all spectypes to SPECTYPE_BB, as part of 
			effort to get restarting models to work better
	080518	ksl	60a - modified all inputs to be in cgs units.  Added a small
      			amount of code to initialize model_list	
	081112	ksl	67 - moved more of initializaiton of geo into
			this routine as part of genearl cleanup of the main
			routine
	1508	ksl	A number of changes have been made in order to accommodate
			domains

**************************************************************/

int
init_geo ()
{
  geo.ndomain = 0;		/*ndomain is a convenience variable so we do not always
				   need to write geo.ndomain but it should nearly always
				   be set to the same value as geo.ndomain */
  geo.run_type = 0;		/* Indicates this is a run from scratch, which includes
				   the case where we already have a wind model but want
				   to change some of the parameters.  init_goe should not
				   be called at all if we are simply continuing a previous
				   run */
  geo.wind_domain_number = -1;
  geo.atmos_domain_number = -1;

/*  The domains have been created but have not been initialized at all */

  zdom[0].coord_type = 1;
  zdom[0].ndim = 30;
  zdom[0].mdim = 30;
  zdom[0].log_linear = 0;	/* Set intervals to be logarithmic */

  zdom[1].coord_type = 1;
  zdom[1].ndim = 30;
  zdom[1].mdim = 10;
  zdom[1].log_linear = 0;	/* Set intervals to be logarithmic */


  geo.disk_z0 = geo.disk_z1 = 0.0;	// 080518 - ksl - moved this up
  geo.adiabatic = 1;		// Default is now set so that adiabatic cooling is included in the wind
  geo.auger_ionization = 1;	//Default is on.


  geo.run_type = 0;		// Not a restart of a previous run

  geo.star_ion_spectype = geo.star_spectype
    = geo.disk_ion_spectype = geo.disk_spectype
    = geo.bl_ion_spectype = geo.bl_spectype = SPECTYPE_BB;
  geo.agn_ion_spectype = SPECTYPE_POW;	// 130605 - nsh - moved from python.c


  geo.rmax = 1e11;
  geo.rmax_sq = geo.rmax * geo.rmax;
  geo.rstar = 7e8;
  geo.rstar_sq = geo.rstar * geo.rstar;
  geo.mstar = 0.8 * MSOL;
  geo.m_sec = 0.4 * MSOL;
  geo.period = 3.2 * 3600;
  geo.tstar = 40000;
  geo.twind = 40000;

  geo.ioniz_mode = IONMODE_ML93;	/* default is on the spot and find the best t */
  geo.line_mode = 3;		/* default is escape probabilites */

  geo.star_radiation = 1;	/* 1 implies star will radiate */
  geo.disk_radiation = 1;	/* 1 implies disk will radiate */
  geo.bl_radiation = 0;		/*1 implies boundary layer will radiate */
  geo.wind_radiation = 0;	/* 1 implies wind will radiate */

  geo.disk_type = DISK_FLAT;		/*1 implies existence of a disk for purposes of absorption */
  geo.diskrad = 2.4e10;
  geo.disk_mdot = 1.e-8 * MSOL / YR;

  geo.t_bl = 100000.;

  geo.pl_geometry = PL_GEOMETRY_SPHERE;   // default to spherical geometry
  geo.lamp_post_height = 0.0;   // should only be used if geo.pl_geometry is PL_GEOMETRY_LAMP_POST


  strcpy (geo.atomic_filename, "data/standard78");
  strcpy (geo.fixed_con_file, "none");

  // Note that geo.model_list is initialized through get_spectype 

  /* Initialize a few other variables in python.h */
  x_axis[0] = 1.0;
  x_axis[1] = x_axis[2] = 0.0;
  y_axis[1] = 1.0;
  y_axis[0] = y_axis[2] = 0.0;
  z_axis[2] = 1.0;
  z_axis[1] = z_axis[0] = 0.0;

  geo.wcycles = geo.pcycles = 1;
  geo.wcycle = geo.pcycle = 0;

  return (0);
}

/*
Perform some simple checks on the photon distribution just produced.

History:
	01	ksl	Removed from main routine
	02jul	ksl	Loosened frequency limits to reflect the
			fact that in some cases, e.g. those in
			which the photon distribution has been split
			into small energy segments, Doppler shifts
			move photons out of that region.
	08mar	ksl	Updated slightly, ane eliminated any frequency
			checks photons generated by macro atoms since
			these often get out of range.
	090124	ksl	Modified slightly to reduce output if all
			is OK and if not debugging

*/
int
photon_checks (p, freqmin, freqmax, comment)
     char *comment;
     PhotPtr p;
     double freqmin, freqmax;
{
  int nnn, nn;
//  double lum_ioniz;  //NSH 16/2/2011 These are now declared externally to allow python to see them
//  int n_ioniz;
  int nlabel;

  geo.n_ioniz = 0;
  geo.lum_ioniz = 0.0;
  nnn = 0;
  nlabel = 0;


  /* Next two lines are to allow for fact that photons generated in
   * a frequency range may be Doppler shifted out of that range, especially
   * if they are disk photons generated right up against one of the frequency
   * limits
   * 04aug--ksl-increased limit from 0.02 to 0.03, e.g from 6000 km/s to 9000 km/s
   * 11apr--NSH-decreased freqmin to 0.4, to take account of double redshifted photons.
   * shift.
   */

  Debug ("photon_checks: %s\n", comment);

  freqmax *= (1.8);
  freqmin *= (0.6);
  for (nn = 0; nn < NPHOT; nn++)
    {
      p[nn].np = nn;		/*  NSH 13/4/11 This is a line to populate the new internal photon pointer */
      if (H * p[nn].freq > ion[0].ip)
	{
	  geo.lum_ioniz += p[nn].w;
	  geo.n_ioniz += p[nn].w / (H * p[nn].freq);
	}
      if (sane_check (p[nn].freq) != 0 || sane_check (p[nn].w))
	{
	  if (nlabel == 0)
	    {
	      Error
		("photon_checks: nphot  origin  freq     freqmin    freqmax\n");
	      nlabel++;
	    }
	  Error
	    ("photon_checks:sane_check %6d %5d %10.4e %10.4e %10.4e %5d w %10.4e \n",
	     nn, p[nn].origin, p[nn].freq, freqmin, freqmax, p[nn].w);
	  p[nn].freq = freqmax;
	  nnn++;
	}
      if (p[nn].origin < 10 && (p[nn].freq < freqmin || freqmax < p[nn].freq))
	{
	  if (nlabel == 0)
	    {
	      Error
		("photon_checks: nphot  origin  freq     freqmin    freqmax\n");
	      nlabel++;
	    }
	  Error
	    ("photon_checks: %6d %5d %10.4e %10.4e %10.4e freq out of range\n",
	     nn, p[nn].origin, p[nn].freq, freqmin, freqmax);
	  p[nn].freq = freqmax;
	  nnn++;
	}
      if (nnn > 100)
	{
	  Error
	    ("photon_checks: Exiting because too many bad photons generated\n");
	  exit (0);
	}
    }
  Log ("NSH Geo.n_ioniz=%e\n", geo.n_ioniz);

  if (nnn == 0)
    Debug ("photon_checks: All photons passed checks successfully\n");

  return (0);
}


/***********************************************************
                                       Space Telescope Science Institute

 Synopsis:
	Generalized routine to get the spectrum type and read the model files 
Arguments:		

Returns:
 
Description:	

		
Notes:

	The routine is slightly dangerous in the sense that if assumes that
	rdint always wants 0 for BB, models for 1, and 2 for uniform.  If
	we were to add another internally generated spectrum type one would
	have to carefull consider how to do this.  

	For models we have to handle two cases:
		A new model.  Here we want to start with a default value and
			to keep track of what was entered since it is likely
			we will want that together
		The continuation of an old model.  Here we need to expect the
			same choices as previously
081026 - Actual routine is still a mess.  
		* yesno is a statment of whether the component exists
		* question is really tightly associted with the way the program is 
		written
		* The progam returns the component type two ways, on through the call
		* and the other through the return.  It ssems like one would do.:




History:
	080518	ksl	Coded as part of effort to make models restart more
			easily, but also simplifies some of the code
        121025  nsh	added a mode for power law

**************************************************************/


char get_spectype_oldname[LINELENGTH] = "data/kurucz91.ls";	/*This is to assure that we read model lists in the same order everytime */
int get_spectype_count = 0;
int
get_spectype (yesno, question, spectype)
     int yesno;
     char *question;
     int *spectype;
{
  char model_list[LINELENGTH];
  int stype;
  int get_models ();		// Note: Needed because get_models cannot be included in templates.h
  if (yesno)
    {
      // XXX This is rather odd. Why are these steps needed? Why don't we fix the question here.  ksl

      // First convert the spectype to the way the questionis supposed to be answered
      if (*spectype == SPECTYPE_BB || *spectype == SPECTYPE_NONE)
	stype = 0;
      else if (*spectype == SPECTYPE_UNIFORM)
	stype = 2;
      else if (*spectype == SPECTYPE_POW)
	stype = 3;
      else
	stype = 1;
      /* Now get the response */
      rdint (question, &stype);
      /* Now convert the response back to the values which python uses */
      if (stype == 0)
	*spectype = SPECTYPE_BB;	// bb
      else if (stype == 2)
	*spectype = SPECTYPE_UNIFORM;	// uniform
      else if (stype == 3)
	*spectype = SPECTYPE_POW;	// power law
      else if (stype == 4)
	*spectype = SPECTYPE_CL_TAB;   // broken power law
	  else if (stype == 5)
	*spectype = SPECTYPE_BREM; // bremstrahlung
      else
	{
	  if (geo.run_type == SYSTEM_TYPE_PREVIOUS)
	    {			// Continuing an old model
	      strcpy (model_list, geo.model_list[get_spectype_count]);
	    }
	  else
	    {			// Starting a new model
	      strcpy (model_list, get_spectype_oldname);
	    }
	  rdstr ("Model_file", model_list);
	  get_models (model_list, 2, spectype);
	  strcpy (geo.model_list[get_spectype_count], model_list);	// Copy it to geo 
	  strcpy (get_spectype_oldname, model_list);	// Also copy it back to the old name
	  get_spectype_count++;
	}
    }
  else
    {
      *spectype = SPECTYPE_NONE;	// No radiation
    }

  return (*spectype);
}


/***********************************************************
                                       Space Telescope Science Institute

 Synopsis:
	The next couple of routines are for recording information about photons/energy impinging
	on the disk, which is stored in a disk structure called qdisk.

	qdisk_init() just initializes the structure (once the disk structue has been initialized.

	qdisk_save records the results in a file

Arguments:		

Returns:
 
Description:	

		
Notes:

	The reason the qdisk is needed as well as disk structure is that whenever the
	wavelength bands are changed the radii in the disk structure are recalibratied.
	We want qdisk to have fixed boundaries when this happens.


History:
	04mar	ksl	add section to dump heating of disk.  
 			Factor of 2 in area calculation reflects
 			fact that disk has two sides
	04dec	ksl	created variable ztot so fractional heating
			is correct if multiple subcyles
	080519	ksl	60a - Added code to calculate the irradiation of the disk
			in terms of t and w.  This should help to monitor the effect
			of irradiation on the disk
	160526	ksl	Modified output file in qdisk_save to be an
			astropy table

**************************************************************/


int
qdisk_init ()
{
  int n;
  for (n = 0; n < NRINGS; n++)
    {
      qdisk.r[n] = disk.r[n];
      qdisk.t[n] = disk.t[n];
      qdisk.g[n] = disk.g[n];
      qdisk.v[n] = disk.v[n];
      qdisk.heat[n] = 0.0;
      qdisk.nphot[n] = 0;
      qdisk.nhit[n] = 0;
      qdisk.w[n] = 0;
      qdisk.ave_freq[n] = 0;
      qdisk.t_hit[0] = 0;
    }
  return (0);
}

int
qdisk_save (diskfile, ztot)
     char *diskfile;
     double ztot;
{
  FILE *qptr;
  int n;
  double area, theat,ttot;
  qptr = fopen (diskfile, "w");
  fprintf (qptr,
	   "r         zdisk     t_disk   heat       nhit nhit/nemit  t_heat    t_irrad  W_irrad  t_tot\n");
  
  for (n = 0; n < NRINGS; n++)
    {
      area =
	(2. * PI *
	 (qdisk.r[n + 1] * qdisk.r[n + 1] - qdisk.r[n] * qdisk.r[n]));
      theat = qdisk.heat[n] / area;
      theat = pow (theat / STEFAN_BOLTZMANN, 0.25);	// theat is temperature if no internal energy production
      if (qdisk.nhit[n] > 0)
	{

	  qdisk.ave_freq[n] /= qdisk.heat[n];
	  qdisk.t_hit[n] = H * qdisk.ave_freq[n] / (BOLTZMANN * 3.832);	// Basic conversion from freq to T
	  qdisk.w[n] =
	    qdisk.heat[n] / (4. * PI * STEFAN_BOLTZMANN * area *
			     qdisk.t_hit[n] * qdisk.t_hit[n] *
			     qdisk.t_hit[n] * qdisk.t_hit[n]);
	}

      ttot=pow(qdisk.t[n],4)+pow(theat,4);
      ttot=pow(ttot,0.25);
      fprintf (qptr,
	       "%8.3e %8.3e %8.3e %8.3e %5d %8.3e %8.3e %8.3e %8.3e %8.3e\n",
	       qdisk.r[n], zdisk (qdisk.r[n]), qdisk.t[n],
	       qdisk.heat[n], qdisk.nhit[n],
	       qdisk.heat[n] * NRINGS / ztot, theat, qdisk.t_hit[n],
	       qdisk.w[n],ttot);
    }

  fclose (qptr);
  return (0);
}



/***********************************************************
	Space Telescope Science Institute

Synopsis:
	Stuart's routine to read a non-standard disk profile
	for YSO effort

Arguments:		

Returns:
 
Description:	

		
Notes:
	Originally part of main routine; moved to separate routine
	by ksl sometime in the fall of 08



History:

**************************************************************/

int
read_non_standard_disk_profile (tprofile)
     char *tprofile;
{

  FILE *fopen (), *fptr;
  int n;
  float dumflt1, dumflt2;
  int dumint;

  if ((fptr = fopen (tprofile, "r")) == NULL)
    {
      Error ("Could not open filename %s\n", tprofile);
      exit (0);
    }

  fscanf (fptr, "%d\n", &dumint);
  blmod.n_blpts = dumint;
  for (n = 0; n < blmod.n_blpts; n++)
    {
      fscanf (fptr, "%g %g", &dumflt1, &dumflt2);
      blmod.r[n] = dumflt1 * 1.e11;
      blmod.t[n] = dumflt2 * 1.e3;
    }

  fclose (fptr);

  return (0);
}


/***********************************************************
				University of Southampton

Synopsis:
	init_advanced_modes simply initialises the set of 
	advanced modes stored in the modes structure to a 
	default value. For now, this is 0 (off).

Arguments:	
    none	

Returns:
    modes is a structure declared in python.h

Description:	
	
Notes:
    see #111 and #120

History:
    1410 -- JM -- Coded
**************************************************************/


int
init_advanced_modes ()
{
  modes.iadvanced = 0;		// this is controlled by the -d flag, global mode control.
  modes.save_cell_stats = 0;	// want to save photons statistics by cell
  modes.ispy = 0;		// want to use the ispy function
  modes.keep_ioncycle_windsaves = 0;	// want to save wind file each ionization cycle
  modes.track_resonant_scatters = 0;	// want to track resonant scatters
  modes.save_extract_photons = 0;	// we want to save details on extracted photons
  modes.print_windrad_summary = 0;	// we want to print the wind rad summary each cycle
  modes.adjust_grid = 0;	// the user wants to adjust the grid scale
  modes.diag_on_off = 0;	// extra diagnostics
  modes.use_debug = 0;
  modes.print_dvds_info = 0;	// print out information on velocity gradients
  modes.quit_after_inputs = 0;	// testing mode which quits after reading in inputs
  modes.fixed_temp = 0;		// do not attempt to change temperature - used for testing
  modes.zeus_connect = 0;             // connect with zeus
  
  //note write_atomicdata  is defined in atomic.h, rather than the modes structure 
  write_atomicdata = 0;		// print out summary of atomic data 


  modes.keep_photoabs = 1;	// keep photoabsorption in final spectrum

  return (0);
}

/***********************************************************
				University of Southampton

Synopsis:
	init_observers sets up the parameter for the final parameters
	to be extracted extracted as


Arguments:	
    none	

Returns:
 
Description:	
	
Notes:

History:
    1509	ksl	Code moved from main after puttting the
    			parameters into the goe structure
**************************************************************/

int
init_observers ()
{
  int n;
  char yesno[20];


  geo.nangles = 4;
  geo.angle[0] = 10;
  geo.angle[1] = 30.;
  geo.angle[2] = 60.;
  geo.angle[3] = 80.;
  for (n = 4; n < NSPEC; n++)
    geo.angle[n] = 45;
  for (n = 0; n < NSPEC; n++)
    {
      geo.phase[n] = 0.5;
      geo.scat_select[n] = 1000;
      geo.top_bot_select[n] = 0;
    }
  geo.swavemin = 850;
  geo.swavemax = 1850;

  rddoub ("spectrum_wavemin", &geo.swavemin);
  rddoub ("spectrum_wavemax", &geo.swavemax);
  if (geo.swavemin > geo.swavemax)
    {
      geo.swavemax = geo.swavemin;
      geo.swavemin = geo.swavemax;
    }

  /* SS June 04: convert these to frequencies and store for use
     in computing macro atom and k-packet emissivities. */

  em_rnge.fmin = C / (geo.swavemax * 1.e-8);
  em_rnge.fmax = C / (geo.swavemin * 1.e-8);

  geo.matom_radiation = 0;	//initialise for ionization cycles - don't use pre-computed emissivities for macro-atom levels/ k-packets.

/* Note: Below here many of the variables which are read in are not currently part of geo stucture */

  rdint ("no_observers", &geo.nangles);

  if (geo.nangles < 1 || geo.nangles > NSPEC)
    {
      Error ("no_observers %d should not be > %d or <0\n", geo.nangles,
	     NSPEC);
      exit (0);
    }


  for (n = 0; n < geo.nangles; n++)
    rddoub ("angle(0=pole)", &geo.angle[n]);

  /* 05apr-ksl-56--For diagnostic reasons I have left questions regarding phase
   * even for systems which are not binaries.  Phase 0 in this case corresponds to
   * an extraction direction which is in the xz plane
   */
  /* JM 1502 -- change this so we only ask for phase if the system is a binary -- see #137 */

  if (geo.system_type == SYSTEM_TYPE_BINARY)
    {

      for (n = 0; n < geo.nangles; n++)
	rddoub ("phase(0=inferior_conjunction)", &geo.phase[n]);
    }
  else
    Log
      ("No phase information needed as system type %i is not a binary\n",
       geo.system_type);


  rdint ("live.or.die(0).or.extract(anything_else)", &geo.select_extract);
  if (geo.select_extract != 0)
    {
      geo.select_extract = 1;
      Log ("OK, extracting from specific angles\n");
    }
  else
    Log ("OK, using live or die option\n");

/* Select spectra with certain numbers of scatterings.  See extract 1997 aug 28 ksl 
 * 141116 - ksl The following options are clealy diagnostic and have been relegated to 
 * advanced commands*/

  if (modes.iadvanced)
    {
      strcpy (yesno, "n");
      rdstr ("Select_specific_no_of_scatters_in_spectra(y/n)", yesno);
      if (yesno[0] == 'y')
	{
	  Log
	    ("OK n>MAXSCAT->all; 0<=n<MAXSCAT -> n scatters; n<0 -> >= |n| scatters\n");
	  for (n = 0; n < geo.nangles; n++)
	    {
	      rdint ("Select_scatters", &geo.scat_select[n]);
	    }
	}
      strcpy (yesno, "n");
      rdstr ("Select_photons_by_position(y/n)", yesno);
      if (yesno[0] == 'y')
	{
	  Log
	    ("OK 0->all; -1 -> below; 1 -> above the disk, 2 -> specific location in wind\n");
	  for (n = 0; n < geo.nangles; n++)
	    {
	      rdint ("Select_location", &geo.top_bot_select[n]);
	      if (geo.top_bot_select[n] == 2)
		{
		  Log
		    ("Warning: Make sure that position will be in wind, or no joy will be obtained\n");
		  rddoub ("rho(cm)", &geo.rho_select[n]);
		  rddoub ("z(cm)", &geo.z_select[n]);
		  rddoub ("azimuth(deg)", &geo.az_select[n]);
		  rddoub ("r(cm)", &geo.r_select[n]);

		}
	    }
	}
    }

  /* Select the units of the output spectra.  This is always needed */

  rdint ("spec.type(flambda(1),fnu(2),basic(other)", &geo.select_spectype);
  if (geo.select_spectype == 1)
    {
      Log ("OK, generating flambda at 100pc\n");
    }
  else if (geo.select_spectype == 2)
    {
      Log ("OK, generating fnu at 100 pc\n");
    }
  else
    Log ("OK, basic Monte Carlo spectrum\n");

  return (0);
}

/***********************************************************
				Space Telescope Science Institute

Synopsis:
	init_photons gets information about the number of 
	cycles and how many photons there should be per cycle.
        The routine then then instantiates PhotPtr

Arguments:	
    none	

Returns:
 
Description:	
	
Notes:

History:
    1509   ksl    Moved the code from python.c
**************************************************************/
PhotPtr
init_photons ()
{
  PhotPtr p;
  double x;

  /* 140907 - ksl - Although photons_per_cycle is really an integer, 
     read in as a double so it is easier for input */

  x = 100000;
  rddoub ("photons_per_cycle", &x);
  NPHOT = x;			// NPHOT is photons/cycle

#ifdef MPI_ON
  Log ("Photons per cycle per MPI task will be %d\n", NPHOT / np_mpi_global);

  NPHOT /= np_mpi_global;
#endif

  rdint ("Ionization_cycles", &geo.wcycles);

  rdint ("spectrum_cycles", &geo.pcycles);


  if (geo.wcycles == 0 && geo.pcycles == 0)
    {
      Log
	("Both ionization and spectral cycles are set to 0; There is nothing to do so exiting\n");
      exit (0);			//There is really nothing to do!
    }

  /* Allocate the memory for the photon structure now that NPHOT is established */

  photmain = p = (PhotPtr) calloc (sizeof (p_dummy), NPHOT);

  if (p == NULL)
    {
      Error
  ("There is a problem in allocating memory for the photon structure\n");
      exit (0);
    }
  else 
    {
      /* JM 1605 -- large photon numbers can cause problems / runs to crash. Report to use (see #209) */
      Log("Allocated %10d bytes for each of %5d elements of photon structure totaling %10.1f Mb \n",
         sizeof (p_dummy), NPHOT, 1.e-6 * NPHOT * sizeof (p_dummy));
      if ( (NPHOT * sizeof (p_dummy)) > 1e9)
        Error("Over 1 GIGABYTE of photon structure allocated. Could cause serious problems.\n");
    }


  return (p);
}




/***********************************************************
				Space Telescope Science Institute

Synopsis:
	init_ioinization

Arguments:	
    none	

Returns:
 
Description:	
	
Notes:

History:
    1509   ksl    Moved the code from python.c
**************************************************************/
int
init_ionization ()
{
  int thermal_opt;


  // XXX  I is unclear to me why all of this dwon to the next XXX is not moved to a single subroutine.  It all
  // pertains to how the radiatiate tranfer is carreid out

  init_maps();
  rdmap
    ("Wind_ionization",
     &geo.ioniz_mode, maps.ion_modes);

  if (geo.ioniz_mode == IONMODE_FIXED)
    {
      rdstr ("Fixed.concentrations.filename", &geo.fixed_con_file[0]);
    }
  if (geo.ioniz_mode == IONMODE_LTE_SIM || geo.ioniz_mode == 5 || geo.ioniz_mode > 9)	/*NSH CLOUDY test - remove once done */
    {
      Log ("The allowed ionization modes are 0, 1, 2, 3, 6, 7\n");
      Error ("Unknown ionization mode %d\n", geo.ioniz_mode);
      exit (0);
    }



  /*Normally, geo.partition_mode is set to -1, which means that partition functions are calculated to take
     full advantage of the data file.  This means that in calculating the partition functions, the information
     on levels and their multiplicities is taken into account.   */

  geo.partition_mode = -1;	//?? Stuart, is there a reason not to move this earlier so it does not affect restart


  /* get_line_transfer_mode reads in the Line_transfer question from the user, 
     then alters the variables geo.line_mode, geo.scatter_mode, geo.rt_mode and geo.macro_simple */

  // XXX - Note clear that get_line_transfer_mode should be a separate routine; perhaps incoroporate here
  get_line_transfer_mode ();




  thermal_opt = 0;		/* NSH 131213 Set the option to zero - the default. The lines allow allow the
				   user to turn off mechanisms that affect the thermal balance. Adiabatic is the only one implemented
				   to start off with. */

  rdint
    ("Thermal_balance_options(0=everything.on,1=no.adiabatic)", &thermal_opt);

  if (thermal_opt == 1)
    {
      geo.adiabatic = 0;
    }

  else if (thermal_opt > 1 || thermal_opt < 0)
    {
      Error ("Unknown thermal balance mode %d\n", thermal_opt);
      exit (0);
    }


  /* 57h -- Next line prevents bf calculation of macro_estimaters when no macro atoms are present.   */

  if (nlevels_macro == 0)
    geo.macro_simple = 1;	// Make everything simple if no macro atoms -- 57h

  //SS - initalise the choice of handling for macro pops.
  if (geo.run_type == SYSTEM_TYPE_PREVIOUS)
    {
      geo.macro_ioniz_mode = 1;	// Now that macro atom properties are available for restarts
    }
  else
    {
      geo.macro_ioniz_mode = 0;
    }

  return (0);

}
