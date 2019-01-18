
/***********************************************************/
/** @file  python.c
 * @author ksl
 * @date   January, 2018
 *
 * @brief  This file contains main and various related routines
 * that are central to the operation of Python
 *
 ***********************************************************/
//OLD /***********************************************************
//OLD                                        Space Telescope Science Institute
//OLD
//OLD  Synopsis:
//OLD   Python is a program designed to simulate the transfer of radiation in a wind.  It uses the
//OLD   Sobolev approximation.  It models a wind as a biconical flow.
//OLD
//OLD   This is the "main" routine for Python.  It's basic purpose is to gather the input variables
//OLD   and to control the flow of the program
//OLD
//OLD Arguments:
//OLD
//OLD   Usage:  py [-h] [-r] [-d] [-f] [-t time_max] xxx  or simply py
//OLD
//OLD   where xxx is the rootname or full name of a parameter file, e. g. test.pf
//OLD
//OLD   and the switches have the following meanings
//OLD
//OLD   -h      to get this help message
//OLD   -r      restart a run of the program reading the file xxx.windsave
//OLD
//OLD   -t time_max     limit the total time to approximately time_max seconds.  Note that the program checks
//OLD           for this limit somewhat infrequently, usually at the ends of cycles, because it
//OLD           is attempting to save the program outputs so that the program can be restarted with
//OLD           -r if theat is desired.
//OLD   -v num  determines the amount of information that is printed out.  If num is small, then
//OLD           less information is printed out; if num is large more is printed out.  Setting
//OLD           v to 5 causes the routine to print out all the information which outputs have
//OLD           included previously.  The current default is set to 3 which suppresses Debug, Log_silent
//OLD           and Error_silent
//OLD   -e nnn  The program stops normally stops after 1e6 errors of any type.  This switch allows one to
//OLD           adjust this number
//OLD   -d      Enters detailed or advanced mode. Allows one to access extra diagnositics and some
//OLD           other advanced commands
//OLD     -f    Fixed temperature mode - does not attempt to change the temperature of cells.
//OLD   -i      Diagnostic mode which quits after reading in inputs. Used for Travis test suite.
//OLD     --dry-run     Create a new .pf file and stop (equivalent to -i)
//OLD   -z      Mode to connect with zeus - it either runs two cycles in this is the first call - in order
//OLD           to obtain a good starting state, else it runs just one cycle. In both cases, it does
//OLD                   not attempt to seek a new temperature, but it does output heating and cooling rates
//OLD     --version     print out python version, commit hash and if there were files with uncommitted
//OLD           changes
//OLD       --rseed     set the random number seed to be time based, rather than fixed.
//OLD
//OLD
//OLD   if one simply types py or pyZZ where ZZ is the version number one is queried for a name
//OLD   of the parameter file.
//OLD
//OLD   NOTE - If this is modified, please also modify the help message in help() below
//OLD Returns:
//OLD
//OLD Description:
//OLD   Python is far too complicated to describe.  Basically it simulates the radiative transfer
//OLD   of photons through the wind of a cataclysmic variable or a star.  The kinematic formulation for
//OLD   the CV wind is due to Schlossman and Vitello while that of the star is due to Castor & Larmors.
//OLD
//OLD    Radiation from an optically thick disk, the WD star, a boundary layer and the wind itself
//OLD    may be included
//OLD
//OLD   There are 3 basic portions to the program which are easy to see in the main program.
//OLD
//OLD   1. A data gathering stage
//OLD
//OLD   2. A calculation of the state of ionization of the wind.
//OLD
//OLD   3. A calculation of a detailed spectrum in which the ionization is held fixed.
//OLD
//OLD
//OLD Notes:
//OLD   The program has been designed and tested both on Suns and MacIntoshes (with Symentec's
//OLD   C compiler).  Some of its peculiarities are due to memory limitations and relatively small
//OLD   stack sizes on the Mac.  When compiling, check to see that the global varible MAC is
//OLD   set properly in python.h.
//OLD
//OLD History:
//OLD   97jan   ksl     Coding on python began.
//OLD   97jul   ksl     Added capability model the spectrum as a function of orbital phase.
//OLD   97aug   ksl     Added boundary layer as additional source of radiation
//OLD   97aug   ksl     Added capability to extract spectra with a specific number of scatters
//OLD   97nov   ksl     Changed error and diagnostic logging
//OLD   98july  ksl     Added radiation from the wind itself
//OLD   98nov   ksl     Added spherical winds and greatly modified I/O
//OLD   99jan   ksl     Added ability to enter parameter file (w/out) root on command line
//OLD   99jul   ksl     Begin eliminating MAC dependencies
//OLD   00sep   ksl     Began work on adding a new coronal possibility
//OLD   01mar   ksl     Added a knigge wind
//OLD         01sept    ksl     Added a thierry wind
//OLD   01dec   ksl     Incorporated improved atomic data reading codes, and
//OLD                   made a variety of other major changes to the way in
//OLD                   which many of the subroutines operate.
//OLD   01dec   ksl     Revised way in which old windfile option works.  It
//OLD                   was not working at all originally, because of a
//OLD                   change I made in April.  The new method is a bit dangerous
//OLD                   because it does change a good bit of geo, since all of the
//OLD                   input data is read in.  In principle, his can create an inconsistency
//OLD                   between some of the info in geo and some in the wind array.
//OLD                   But it keeps the wind the same, which is what was desired.
//OLD   01dec   ksl     Added a general capability to python(40) so that during the
//OLD                   the ionization cycle photons do not necessarily have identical
//OLD                   weights.
//OLD   02jan   ksl     Moved specific inputs for KWD wind to knigge.c
//OLD   02feb   ksl     Added anisotropic wind scattering
//OLD   02feb   ksl     41.1 -- Allowed user to chose between anisotropic and
//OLD                   isotropic scattering in the wind.
//OLD   02feb   ksl     41.2 -- Added aditional options for choosing what photons
//OLD                   to count
//OLD   02apr   ksl     43.1 -- Added additional radiative transfer options.  Modified main
//OLD                   so that a spectrum is printed out after each cycle
//OLD   02jul   ksl     43.7 -- Considerable changes associated with photoionization
//OLD                   and recombination were made.  The initial attemps to allow
//OLD                   for using a more accurate detailed balance approach were
//OLD                   included in the ionizatio routines.
//OLD   03dec   ksl     Added back some timing ability
//OLD         04Mar   SS      Minor changes to set switch for macro atom method.
//OLD                         geo.line_mode=6 switches geo.rt_mode=2 (macro method)
//OLD                         and then geo.ine_mode back to =3. (SS)
//OLD         04Apr   SS      Minor change to set a switch for geo.line_mode=7
//OLD                         which is the same as 6 but with anisotropic scattering
//OLD         04May   SS      geo.line_mode=8 added: this mode is the same as mode 6
//OLD                         BUT will treat all ions as simple - it is for
//OLD                         test purposes only.
//OLD   04Jun   ksl     Added new disk + wind model for ysos.  In process, moved
//OLD                   the variables describing the wind mdots to subroutines
//OLD   04Aug   ksl     52 -- Modified input variables to allow a variety of input models
//OLD                   all ascii-based.
//OLD   04Aug   ksl     52 -- Modified input variables to allow for various disk illumination
//OLD                   models
//OLD   04Aug   ksl     52 -- Modified to allow for vertically extended disks. Note that some of
//OLD                   the input variables were rearranged in the process.
//OLD         04Aug   SS      Minor modifications to cope with finite disk thickness calculations.
//OLD   04dec   ksl     54a -- Minor mod to calculation of fraction of luminosity that hits
//OLD                   disk so that correct values are printed out in situation where there
//OLD                   are multiple subcycles.
//OLD   04dec   ksl     54e -- Minor mod to record master atomic file in geo, and, for the
//OLD                   case of fixed concentration to store the file where the concentrations
//OLD                   are contaned.
//OLD   05apr   ksl     56 -- Added new variable geo.binary_system to avoid binary system
//OLD                   questions in cases where one was modelling a single object.  Note
//OLD                   that I did not remove seemingly superfluous quesions regarding
//OLD                   phase because one can still extract photons in various azimuthal
//OLD                   directions, and if we ever move to 3d situations one might want
//OLD                   to deal with non-axially symmetric systems.
//OLD   05jul   ksl     56d -- Updated to include new windcone definitions.
//OLD   05jul   ksl     56d -- Completed basic work to allow vertically extended disk
//OLD   06may   ksl     57g -- Began work on modifying structures to lower memory requirments
//OLD                   At empty wind cells have mostly been eliminated, and macroatoms
//OLD                   have been split off into a separate structure so this is not
//OLD                   used.
//OLD   06jul   ksl     57h -- Have contined work to speed program up as much as possible
//OLD                   by giving more control over the parameters that govern the speed
//OLD                   such as SMAX_FRAC.
//OLD   06oct   ksl     58 -- This version is a cleaned up version of py57ib, which differs
//OLD                   from 57h primarily in bb.c and pdf.c.   In addition, this version
//OLD                   makes it possible to use the SV velocity law in what is otherwise
//OLD                   a KWD description of the wind.
//OLD   06nov   ksl     58c -- This version requires and updated version of kpar, that
//OLD                   keeps track of all the errors.
//OLD   07aug   ksl     58f -- Modified the way in which the program handles a situation
//OLD                   in which one attempts to run with a non-existent .pf file.
//OLD   08may   ksl     60a -- Significant modifications were made to the main routine
//OLD                   to make sure we could restart models.  Unfortunaaely rhe point
//OLD                   at which the user was asked for the atomic data needed to be changed,
//OLD                   so earlier .pf files will need to be modified.
//OLD   080808  ksl     62 -- Cleaned up wind ionization options
//OLD   081110  ksl     67 -- Revised a number of portions of the main python.c routine
//OLD                   to enable restarts.  There are still inconsistencies in the
//OLD                   way one reads information from the .pf file after having read
//OLD                   in the windsave files that could be cleaned up.
//OLD   081218  ksl     67c -- Added a switch -h to provide information about usage.
//OLD   090202  ksl     68b -- Added switch -v  to control the level of verbosity.  Updated
//OLD                   Plasma to allow routine to track scatters and absorption during
//OLD                   generation of detailed spectrum
//OLD   090402  ksl     NSPEC has been moved to the main routine, as its only remaining
//OLD                   purpose is to define some arrays in the main routine.  Note that
//OLD                   MSPEC still has meaning as the number of spectra of various types
//OLD                   that are construced without going through types.
//OLD   1108    ksl/nsh Adding code to keep track of gross spectra in a cell though xbands,
//OLD                   xj and xave_freq.  Just set up the frequence limits here.
//OLD   1112    ksl     Moved everything associated with frequency bands into bands_init
//OLD   1212    nsh     changed the way DFUDGE is defined.
//OLD         1302      jm      74b5 introduced double precision variable for photon number for calling
//OLD                   define_phot. Fixes Bug JM130302 in photon weights. It has been suggested
//OLD                   that we should make NPHOT a long integer- at the moment I have not done
//OLD                   this and simply comverted to double before calling define_phot (otherwise
//OLD                   I believe rdint would have to be redone for long integers).
//OLD   1304    ksl     75 rewrote the fix above to use a long, instead of double.
//OLD                   This has plenty of range.  Notge that there is no need to make NPHOT
//OLD                   a long, as suggested above.  We do not expect it to exceed 2e9,
//OLD                   although some kind of error check might be reasonble.
//OLD   1306    ksl     Modified to allow a power law component to a stellar spectrum.  Made
//OLD                   some changes use DEF variables instead of numbers to make choices
//OLD   1307    jm      SS Parallelized Python in June 2013, for release 76a. I have now introduced
//OLD                   slightly altered reporting to allow more succinct reporting in parallel mode.
//OLD   1307    ksl     Removed the Thierry & Hubeny O-star models as an option from the code.
//OLD                   This was never tested, and never really used.  Knox no longer even has the
//OLD                   models.  Note that Stuart is replacing this with a homologous expansion
//OLD                   model
//OLD   1308    nsh     Added a call to generate rtheta wind cones - issue #41
//OLD   1309    nsh     Changed the loop around where disk parameters are read in - issue #44
//OLD   1309    nsh     Added commands to write out warning summary - relating to issue #47
//OLD           1312    nsh     Added a new parameter file command to turn off heating and cooling mechanisms
//OLD                   at the moment it only does adiabatc
//OLD                   Also some modifications to the parallel communiactions to deal with some new
//OLD                   plasma variabales, and the min and max frequency photons seen in bands.
//OLD   1409    ksl     Added new switch -d to enable use of a new Debug logging feature
//OLD   1411    JM removed photons per cycle in favour of NPHOT. subcycles are now eliminated
//OLD   1501    JM moved some parallelization stuff to subroutines in para_update.c
//OLD                   functions are communicate_estimators_para, communicate_matom_estimators_para,
//OLD                   and gather_spectra_para
//OLD   1502            Major reorganisation of input gathering and setup. See setup.c, #136 and #139
//OLD   1508    ksl     Instroduction of the concept of domains to handle the disk and wind as separate
//OLD                   domains
//OLD     1711    ksl Removed the Elvis wind option
//OLD
//OLD   Look in Readme.c for more text concerning the early history of the program.
//OLD
//OLD **************************************************************/



#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "atomic.h"
#include <time.h>               //To allow the used of the clock command without errors!!



#include "python.h"
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
      exit (0);
    }
    w = wmain;

    geo.run_type = RUN_TYPE_RESTART;    // We are continuing an old run

    xsignal (files.root, "%-20s Read %s\n", "COMMENT", files.old_windsave);

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

    rdint ("System_type(0=star,1=binary,2=agn,3=previous)", &geo.system_type);

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
        exit (0);
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
        exit (0);
      }

      /* Describe the wind (or more correctly the various domains).
       */


      rdpar_comment ("Parameters descibing the various winds or coronae in the system");

      rdint ("Wind_radiation(y=1)", &geo.wind_radiation);
      /* JM 1806 -- note that wind radiation will get "turned off" in indivisible packet/macro-atom
         mode when geo.rt_mode == RT_MODE_MACRO. This is done in get_line_transfer_mode () in
         setup_domains.c, see issue #390 */

      if (geo.run_type == RUN_TYPE_NEW)
      {
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
 * ratiation and then value given to return a spectrum type. The output is not the same
 * number as one inputs. It' s not obvious that this is a good idea. */

  if (geo.pcycles > 0)
  {

    rdpar_comment ("Parameters defining the spectra seen by observers\n");

    get_spectype (geo.star_radiation,
                  //"Rad_type_for_star(0=bb,1=models,2=uniform)_in_final_spectrum",
                  "Central_object.rad_type_in_final_spectrum(0=bb,1=models,2=uniform)", &geo.star_spectype);
    get_spectype (geo.disk_radiation,
                  //"Rad_type_for_disk(0=bb,1=models,2=uniform)_in_final_spectrum",
                  "Disk.rad_type_in_final_spectrum(0=bb,1=models,2=uniform)", &geo.disk_spectype);
    get_spectype (geo.bl_radiation,
                  //"Rad_type_for_bl(0=bb,1=models,2=uniform)_in_final_spectrum",
                  "Boundary_layer.rad_type_in_final_spectrum(0=bb,1=models,2=uniform)", &geo.bl_spectype);
    geo.agn_spectype = 3;
    get_spectype (geo.agn_radiation, "Rad_type_for_agn(3=power_law,4=cloudy_table,5=bremsstrahlung)_in_final_spectrum", &geo.agn_spectype);
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
    rdint ("@Diag.extra(0=no,1=yes) ", &modes.diag_on_off);
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
