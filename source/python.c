
/***********************************************************/
/** @file  python.c
 * @author ksl
 * @date   January, 2018
 *
 * @brief  This file contains main and various related routines
 * that are central to the operation of Python
 *
 ***********************************************************/


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "atomic.h"
#include <time.h>               //To allow the used of the clock command without errors!!



#include "python.h"
#include "models.h"
#define NSPEC	20




/**********************************************************/
/**
 * @brief     The main routine for Python, which supervises the ingest of data defining a model, the initializtion
 * of the model, and actual run of the model, and writing the data to the disk
 *
 *
 * @param [in] int  argc   The number of command line arguments
 * @param [in] char *  argv[]   The command line arguments
 * @return     Always returns 0
 *
 * @details
 *
 * ### Notes ###
 *
 * The main routine of Python is fairly simple conceptually.  One gather the data, one allocates and intializes
 * all of the data structures, one runs the ionization cycles and one run the detalied spectral cycles.
 * The main ruoutin of Python supervises all of these things
 *
 * The main routine is complicated (mainly) due to the multiple ways a model can be run.
 *
 *  * A model can be run from scratch
 *  * A model can be contined from an earlier run (increasing only the number of ionization or spectral cycles)
 *  * A new model can be run beginning with an old windsave file (allowing various parameters of the radiation souces to be set differently
 *  from the old model)
 *
 *  Additionally, there is logic associated with the fact that different types of models require different data.  An AGN for example does
 *  not have (the possibility of) a secondary star.  Most of the inputs from Python come from a parameter file, but Python also has a robust
 *  set of command line switches (See parse_command_line).
 *
 *  Once all the inputs are obtained, the main routine calles routines to allocate data structures needed  to hold the model and to complete
 *  the intialization of both the data structues and various other variables.
 *
 *  The it runs the ionization cycles (See calculate_ionization );
 *
 *  Finally it runs the cycles that calulate the detailed spectra (See make_spectra);
 *
 *
 **********************************************************/

int
main (argc, argv)
     int argc;
     char *argv[];
{
  WindPtr w;

  double freqmin, freqmax;
  int n;
  char answer[LINELENGTH];
  int get_models ();            // Note: Needed because get_models cannot be included in templates.h
  int dummy_spectype;

  FILE *fopen ();

  int opar_stat, restart_stat;
  double time_max;              // The maximum time the program is allowed to run before halting
  double lstar;                 // The luminosity of the star, iv it exists

  int my_rank;                  // these two variables are used regardless of parallel mode
  int np_mpi;                   // rank and number of processes, 0 and 1 in non-parallel





#ifdef MPI_ON
  MPI_Init (&argc, &argv);
  MPI_Comm_rank (MPI_COMM_WORLD, &my_rank);
  MPI_Comm_size (MPI_COMM_WORLD, &np_mpi);
#else
  my_rank = 0;
  np_mpi = 1;
#endif

  np_mpi_global = np_mpi;       // Global variable which holds the number of MPI processes
  rank_global = my_rank;        // Global variable which holds the rank of the active MPI process
  Log_set_mpi_rank (my_rank, np_mpi);   // communicates my_rank to kpar


  opar_stat = 0;                /* Initialize opar_stat to indicate that if we do not open a rdpar file,
                                   the assumption is that we are reading from the command line */
  restart_stat = 0;             /* Assume initially that these is a new run from scratch, and not
                                   a restart */
  time_max = 13.8e9 * 3.2e7;    /* The maximum time the program will run without stopping.  This
                                   is initially set to the lifetime of the universe */
  time_max = -1;
  set_max_time (files.root, time_max);


  /* Set the verbosity level for logging.  To get more info raise the verbosity level to a higher number. To
     get less set the verbosity to a lower level. The verbosity can be reset from the comamnd line */

  verbosity = 3;
  Log_set_verbosity (verbosity);

  /* initialise the advanced mode flags (all set to 0) which is a structure in python.h */

  init_advanced_modes ();

  /* Parse the command line. Get the root. create files.diagfolder + diagfiles */

  strict = 0;

  restart_stat = parse_command_line (argc, argv);

  /* If the restart flag has been set, we check to see if a windsave file exists.  If it doues we will
     we will restart from that point.  If the windsave file does not exist we will start from scratch */

  init_log_and_windsave (restart_stat);

  Log_parallel ("Thread %d starting.\n", my_rank);

  /* Start logging of errors and comments */

  Log ("!!Python Version %s \n", VERSION);      //54f -- ksl -- Now read from version.h
  Log ("!!Git commit hash %s\n", GIT_COMMIT_HASH);

  /* warn the user if there are uncommited changes */

  int git_diff_status = GIT_DIFF_STATUS;
  if (git_diff_status > 0)
    Log ("!!Git: This version was compiled with %i files with uncommitted changes.\n", git_diff_status);

  Log ("!!Python is running with %d processors\n", np_mpi_global);
  Log_parallel ("This is MPI task number %d (a total of %d tasks are running).\n", rank_global, np_mpi_global);

  Debug ("Debug statements are on. To turn off use lower verbosity (< 5).\n");


  xsignal (files.root, "%-20s Initializing variables for %s\n", "NOK", files.root);

  opar_stat = setup_created_files ();

/* Allocate the domain structure */

  zdom = (DomainPtr) calloc (sizeof (domain_dummy), MaxDom);

  /* BEGIN GATHERING INPUT DATA */

  /* Describe the basic calculation in terms of the number of iterations which will
     be used to calculate the wind parameters and the number of iterations and wavelength
     range which will be used for the final spectrom.  Also describe the observer's views
     of the system */


  if (restart_stat == 1)        /* We want to continue a run. This is generally used
                                   because we had to limit to runtime of python or we decided
                                   we needed more ionization or spectral cycles */
  {
    Log ("Continuing a previous run of %s \n", files.root);
    strcpy (files.old_windsave, files.root);
    strcat (files.old_windsave, ".wind_save");

    if (wind_read (files.old_windsave) < 0)
    {
      Error ("python: Unable to open %s\n", files.old_windsave);        //program will exit if unable to read the file
      Exit (0);
    }
    w = wmain;

    geo.run_type = RUN_TYPE_RESTART;    // We are continuing an old run

    xsignal (files.root, "%-20s Read %s\n", "COMMENT", files.old_windsave);


    if (geo.model_count > 0)    //We have previously used models - we need to read them in again
    {
      for (n = 0; n < geo.model_count; n++)
      {
        get_models (geo.model_list[n], 2, &dummy_spectype);
      }
    }
    if (geo.pcycle > 0)
    {
      spec_read (files.specsave);
      xsignal (files.root, "%-20s Read %s\n", "COMMENT", files.specsave);
    }
  }

  else if (restart_stat == 0)   /* We are starting a new run, which is the normal mode of operation */
  {

    /* First,  establish the overall system type.  System type should be a physical system,
     * to make things easier for the user.  So really want system types to be something like
     * CV, YSO, AGN so that defaults can be set.  This is now issue #420
     */

    geo.system_type = SYSTEM_TYPE_STAR;
    geo.run_type = RUN_TYPE_NEW;

    strcpy (answer, "star");
    geo.system_type = rdchoice ("System_type(star,binary,agn,previous)", "0,1,2,3", answer);


    if (geo.system_type == SYSTEM_TYPE_PREVIOUS)
    {

      /* This option is for the confusing case where we want to start with a previous wind
         model,(presumably because that run produced a wind close to the one we are looking for,
         but we are going to change some parameters that do not affect the wind geometry,
         We will write use new filenames for the results, so all of the previous work is still saved,
       */

      strcpy (files.old_windsave, "earlier.run");
      rdstr ("Wind.old_windfile(root_only)", files.old_windsave);
      strcat (files.old_windsave, ".wind_save");


      Log ("Starting a new run from scratch starting using a previous windfile\n");

      /* Note that wind_read also reads the atomic data file that was used to create the previous run of the data. */

      if (wind_read (files.old_windsave) < 0)
      {
        Error ("python: Unable to open %s\n", files.old_windsave);      //program will exit if unable to read the file
        Exit (0);
      }

      geo.run_type = RUN_TYPE_PREVIOUS; // after wind_read one will have a different wind_type otherwise

      w = wmain;
      geo.wcycle = 0;
      geo.pcycle = 0;           /* This is a new run of an old windsave file so we set the nunber of cycles already done to 0 */
    }



    if (geo.run_type == RUN_TYPE_NEW || geo.run_type == RUN_TYPE_PREVIOUS)
    {
      /* This option is the most common one, where we are starting to define a completely new system.
       */


      if (geo.run_type == RUN_TYPE_NEW)
      {
        init_geo ();            /* Set values in the geometry structure and the domain stucture to reasonable starting
                                   values */
      }

      /* get_stellar_params gets information like mstar, rstar, tstar etc.
         it returns the luminosity of the star */

      lstar = get_stellar_params ();

      /* Describe the disk */

      get_disk_params ();


      /* describe the boundary layer / agn components to the spectrum if they exist.
         So that initial condiditions for the bl and agn are initialized sensibly this has
         to come after the disk is defined.
       */

      get_bl_and_agn_params (lstar);

      /* At this point we check whether we have any sources of radiation and exit if we do not */

      if (!geo.star_radiation && !geo.disk_radiation && !geo.bl_radiation && !geo.bl_radiation && !geo.agn_radiation)
      {
        Error ("python: No radiation sources so nothing to do but quit!\n");
        Exit (0);
      }

      /* Describe the wind (or more correctly the various domains).
       */


      rdpar_comment ("Parameters descibing the various winds or coronae in the system");

      strcpy (answer, "yes");
      geo.wind_radiation = rdchoice ("Wind.radiation(yes,no)", "1,0", answer);

      /* JM 1806 -- note that wind radiation will get "turned off" in indivisible packet/macro-atom
         mode when geo.rt_mode == RT_MODE_MACRO. This is done in get_line_transfer_mode () in
         setup_domains.c, see issue #390 */

      if (geo.run_type == RUN_TYPE_NEW)
      {
        geo.ndomain = 1;
        rdint ("Wind.number_of_components", &geo.ndomain);


        for (n = 0; n < geo.ndomain; n++)
        {

          get_domain_params (n);

        }
      }




    }

  }




/* Get the remainder of the input data.  Note that the next few lines are read from the input file whether or not the windsave file was read in,
   because these are things one would like to be able to change even if we have read in an old windsave file.  init_photons reads in
   the numbers of ionization and spectral cycles and then instatiates PhotPtr.  It is possilbe that this should be moved to else where in
   the flow of reading in data.
 */

  /* All operating modes */
  rdpar_comment ("Parameters associated with photon number, cycles,ionization and radiative transfer options");

  init_photons ();

  /* Define how ionization is going to be calculated */

  /* All operating modes */
  init_ionization ();




  /* Note: ksl - At this point, SYSTEM_TYPE_PREVIOUS refers both to a restart and to a situation where
   * one is starting from an early wind file as implemented this is quite restrictive about what one
   * can change in the previous case.   */

  if (geo.run_type == RUN_TYPE_NEW)     // Start of block to define a model for the first time
  {

    /* Describe the wind, by calling get_wind_params one or more times
       and then gets params by calling e.g. get_sv_wind_params() */

    for (n = 0; n < geo.ndomain; n++)
    {
      rdpar_comment ("Parameters for Domain %d", n);
      get_wind_params (n);
    }

  }                             // End of block to define a model for the first time
  else if (modes.zeus_connect == 1)     /* We are in rad-hydro mode, we want the new density and temperature */
  {
    /* Hydro takes the wind domain number as an argument in the current domains setup */
    Log ("We are going to read in the density and temperature from a zeus file\n");
    get_hydro (geo.hydro_domain_number);        //This line just populates the hydro structures
  }


  /* Calculate additional parameters associated with the binary star system */

  if (geo.system_type == SYSTEM_TYPE_BINARY)
    binary_basics ();

  /* Check that the parameters which have been supplied for the star, disk and boundary layer will
     allow generation of photons where that is appropriate */

  if (geo.tstar <= 0.0)
    geo.star_radiation = 0;
  if (geo.disk_mdot <= 0.0 && geo.disk_tprofile == DISK_TPROFILE_STANDARD)
    geo.disk_radiation = 0;
  if (geo.t_bl <= 0.0 || geo.lum_bl <= 0.0)
    geo.bl_radiation = 0;

  /* If the disk radius is <0, assume no disk was intended. */

  if (geo.diskrad <= 0.0)
  {
    geo.disk_type = DISK_NONE;
    geo.disk_radiation = 0;
  }

  if (geo.star_radiation)
    Log ("There is a star which radiates\n");
  else
    Log ("The star in the system does not radiate\n");

  if (!geo.disk_type)
    Log ("There is no disk in the system \n");
  else if (!geo.disk_radiation)
    Log ("The disk exists, but only absorbs photons\n");
  else
    Log ("There is a disk which radiates and absorbs\n");

  if (geo.bl_radiation)
    Log ("There is a boundary layer which radiates\n");
  else
    Log ("There is no boundary layer\n");

  if (geo.agn_radiation)
    Log ("There is a BH  which radiates\n");
  else
    Log ("There is no BH \n");

  /* Describe the spectra which will be extracted and the way it will be extracted */

  /* First initialise things to semi-reasonable values */
/* These two variables have to do with what types of spectra are created n the
 * spectrum files. They are not associated with the nature of the spectra that
 * are generated by say the boundary layer
 */

  geo.select_extract = 1;
  geo.select_spectype = 1;

/* Completed initialization of this section.  Note that get_spectype uses the source of the
 * radiation and then value given to return a spectrum type. The output is not the same
 * number as one inputs. It' s not obvious that this is a good idea. */

  if (geo.pcycles > 0 && geo.pcycle == 0)
  {
    // This should only evaluate true for when no spectrum cycles have run (I hope)

    rdpar_comment ("Parameters defining the spectra seen by observers\n");

    get_spectype (geo.star_radiation,
                  //"Rad_type_for_star(0=bb,1=models,2=uniform)_in_final_spectrum",
                  //"Central_object.rad_type_in_final_spectrum(0=bb,1=models,2=uniform)", &geo.star_spectype);
                  "Central_object.rad_type_in_final_spectrum(bb,models,uniform)", &geo.star_spectype);
    get_spectype (geo.disk_radiation,
                  //"Rad_type_for_disk(0=bb,1=models,2=uniform)_in_final_spectrum",
                  //"Disk.rad_type_in_final_spectrum(0=bb,1=models,2=uniform)", &geo.disk_spectype);
                  "Disk.rad_type_in_final_spectrum(bb,models,uniform)", &geo.disk_spectype);
    get_spectype (geo.bl_radiation,
                  //"Rad_type_for_bl(0=bb,1=models,2=uniform)_in_final_spectrum",
                  //"Boundary_layer.rad_type_in_final_spectrum(0=bb,1=models,2=uniform)", &geo.bl_spectype);
                  "Boundary_layer.rad_type_in_final_spectrum(bb,models,uniform)", &geo.bl_spectype);
    geo.agn_spectype = SPECTYPE_POW;
    get_spectype (geo.agn_radiation,
                  //"Rad_type_for_agn(3=power_law,4=cloudy_table,5=bremsstrahlung)_in_final_spectrum", &geo.agn_spectype);
                  "BH.rad_type_in_final_spectrum(power,cloudy,brems)", &geo.agn_spectype);
    if (geo.agn_radiation && geo.agn_spectype >= 0 && comp[geo.agn_spectype].nmods != 1)
    {
      Error ("python: When using models with an AGN, there should be exactly 1 model, we have %i for spectrum cycles\n",
             comp[geo.agn_ion_spectype].nmods);
      exit (0);
    }
    init_observers ();
  }

  geo.matom_radiation = 0;      //initialise for ionization cycles - don't use pre-computed emissivities for macro-atom levels/ k-packets.
  get_standard_care_factors ();
  get_meta_params ();

/* Establish limits on the frequency intervals to be used by the ionization cycles and
 * the fraquency bands for stratified sampling. These bands are alos used as the spectral
 * intervals for creating crude spectra in each of the cells
 *
 * This section of inputs might logically go earlier in the code, but
 * was put here so it would add on to existing .pf files.  It would be reasonble to consider moving
 * it to a more logical location
 */


  rdpar_comment ("Other parameters");

  bands_init (-1, &xband);
  freqmin = xband.f1[0];
  freqmax = xband.f2[xband.nbands - 1];

  if (modes.iadvanced)
  {
    /* Do we require extra diagnostics or not */
    strcpy (answer, "no");
    modes.diag_on_off = rdchoice ("@Diag.extra(yes,no)", "1,0", answer);
    if (modes.diag_on_off)
    {
      get_extra_diagnostics ();
    }
  }


  /* Wrap up and save all the inputs */


  if (strncmp (files.root, "mod", 3) == 0)
  {
    cpar ("mod.pf");
  }

  else if (opar_stat == 1)
  {
    cpar (files.input);
  }
  else
  {
    cpar (files.new_pf);
  }

  /* At this point, all inputs have been obtained at this point and the inputs have been copied to "mod.pf" or "python.pf"
   * If we have used, the -i flag, we quit; otherwise we continue on to run the model */
  if (modes.quit_after_inputs)
  {
    Log ("This was was run with the -i or --dry-run flag set, so quitting now inputs have been gathered.\n");
    error_summary ("dry run.");
    exit (0);
  }


  if (strict)
  {
    Log ("Some of the input have not been updated for the current version of Python.  Please correct and rerun\n");
    exit (0);
  }

  /* INPUTS ARE FINALLY COMPLETE */

  /* Print out some diagnositic infomration about the domains */




  Log ("There are %d domains\n", geo.ndomain);
  for (n = 0; n < geo.ndomain; n++)
  {
    Log ("%20s type: %3d  ndim: %3d mdim: %3d ndim2: %4d\n", zdom[n].name, zdom[n].wind_type, zdom[n].ndim, zdom[n].mdim, zdom[n].ndim2);
  }




  /* DFUDGE is a distance that assures we can "push through" boundaries.  setup_dfudge
     sets the push through distance depending on the size of the system.
   */

  DFUDGE = setup_dfudge ();

  /* Now define the wind cones generically. modifies the global windcone structure */
  setup_windcone ();




  /* initialize the random number generator */
  /* By default, the random number generator start with fixed seeds (differnt
   * for each processor, but this can be changed using a command line
   * switch.
   *
   * An exception is when we are in zeus mode, where it would be inappropriate
   * to use the same phtons in ezch cycle.  There we initiate the seeds unsing
   * the clock
   */
  if (modes.rand_seed_usetime == 1)
  {
    n = (unsigned int) clock () * (rank_global + 1);
    init_rand (n);

  }
  else
  {
    init_rand (1084515760 + (13 * rank_global));
  }



  /* Next line finally defines the wind if this is the initial time this model is being run */

  if (geo.run_type == RUN_TYPE_NEW)     // Define the wind and allocate the arrays the first time
  {
    define_wind ();
  }



  if (modes.zeus_connect == 1)  //We have restarted, but are in zeus connect mode, so we want to update density, temp and velocities
  {

    /* Hydro takes the wind domain number as an argument in the current domains setup */
    hydro_restart (geo.hydro_domain_number);
  }

  /* this routine checks, somewhat crudely, if the grid is well enough resolved */
  check_grid ();

  w = wmain;
  if (modes.extra_diagnostics)
  {
    /* Open a diagnostic file or files (with hardwired names) */

    init_extra_diagnostics ();
  }




  /* Start with photon history off */
  phot_hist_on = 0;

  /* If required, read in a non-standard disk temperature profile */
  if (geo.disk_tprofile == 1)
  {
    read_non_standard_disk_profile (files.tprofile);
  }



/* The next section sets up a structure qdisk to record the effects
 * of illumination on the disk.  disk_init is called primarily to get
 * a defined set of annular rings which are kept throughout the
 * ionization calculation.  A second structure qdisk is needed
 * because in the process of generating photons in various bands
 * the annular rings are changed
 *
 * disk_init calculates the flux from the disk in the energy range set by
 * freqmin and freqmax, and uses is this to identify the position of the
 * rings in the disk, so that each ring contributes the same amount to
 * the flux
 *
 * */


  disk_init (geo.rstar, geo.diskrad, geo.mstar, geo.disk_mdot, freqmin, freqmax, 0, &geo.f_disk);
  qdisk_init ();                /* Initialize a disk qdisk to store the information about photons impinging on the disk */
  xsignal (files.root, "%-20s Finished initialization for %s\n", "NOK", files.root);
  check_time (files.root);

/* XXXX - THE CALCULATION OF THE IONIZATION OF THE WIND */
  geo.ioniz_or_extract = 1;     //SS July 04 - want to compute MC estimators during ionization cycles
  //1 simply implies we are in the ionization section of the code
  //and allows routines to act accordinaly.
/* 67 -ksl- geo.wycle will start at zero unless we are completing an old run */

/* XXXX -  CALCULATE THE IONIZATION OF THE WIND */
  calculate_ionization (restart_stat);

/* XXXX - END OF CYCLE TO CALCULATE THE IONIZATION OF THE WIND */
  Log (" Completed wind creation.  The elapsed TIME was %f\n", timer ());
  /* SWM - Evaluate wind paths for last iteration */
  if (geo.reverb == REV_WIND || geo.reverb == REV_MATOM)
  {                             //If this is a mode in which we keep wind arrays, update them
    wind_paths_evaluate (w, my_rank);
  }

/* XXXX - THE CALCULATION OF A DETAILED SPECTRUM IN A SPECIFIC REGION OF WAVELENGTH SPACE */

  freqmax = C / (geo.swavemin * 1.e-8);
  freqmin = C / (geo.swavemax * 1.e-8);

  /* Perform the initilizations required to handle macro-atoms during the detailed
     calculation of the spectrum.

     Next lines turns off macro atom estimators and other portions of the code that are
     unnecessary during spectrum cycles.  */

  geo.ioniz_or_extract = 0;


/* Next step speeds up extraction stage */

  if (!modes.keep_photoabs)
  {
    DENSITY_PHOT_MIN = -1.0;    // Do not calculated photoabsorption in detailed spectrum
  }

  /*Switch on k-packet/macro atom emissivities  SS June 04 */

  if (geo.rt_mode == RT_MODE_MACRO)
  {
    geo.matom_radiation = 1;
  }

  /* Finished initializations required for macro-atom approach */

  /* Calculate and store which bf processess need to be considered in each cell
   * Note that this is not macro-specific but is just to speed the program up.
   */

  kbf_need (freqmin, freqmax);

  /* XXXX - Execute  CYCLES TO CREATE THE DETAILED SPECTRUM */
  make_spectra (restart_stat);

  return (0);
}
