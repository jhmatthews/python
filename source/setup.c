/***********************************************************/
/** @file  setup.c
 * @author ksl
 * @date   March, 2018
 *
 * @brief  Various intialization routines, including a number
 * that read values from the parameter files
 *
 * ### Notes ###
 *
 * Because the input and setup of Python is relatively complex,
 * we have over time moved much of this out of main into
 * subroutines, and we have generally tried to collect
 * related portion of the initialization into different setup
 * routines, such as setup_disk or setup_files.
 *
 * This file contains a collection of these setup routines,
 * that either stand-alone or have not been moved into their
 * own files.
 ***********************************************************/


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "atomic.h"


#include "python.h"



/**********************************************************/
/**
 * @brief      initializes the geo structure to something that is semi-reasonable
 *
 * @return     Always returns 0
 *
 * @details
 * Initial values for all of the variables that are not part of the individual
 * wind descriptions that are actully read(!) into the program should be
 * created here.  The derived values are not needed.
 *
 * ### Notes ###
 *
 * In general, cgs units are the working units for Python.  Thus all
 * intialization should be converted to these units.
 *
 * @bug Currently init_geo is set up for CVs and Stars and not AGN.  We now
 * read in the system type as the first variable. The intent was to
 * allow one to use the system type to particularize how geo (aqnd
 * other variables) were intialized. But this has yet to be carreid
 * out.
 *
 *
 **********************************************************/

int
init_geo ()
{
  geo.ndomain = 0;              /*ndomain is a convenience variable so we do not always
                                   need to write geo.ndomain but it should nearly always
                                   be set to the same value as geo.ndomain */
  geo.run_type = 0;             /* Indicates this is a run from scratch, which includes
                                   the case where we already have a wind model but want
                                   to change some of the parameters.  init_goe should not
                                   be called at all if we are simply continuing a previous
                                   run */
  geo.hydro_domain_number = -1;

  if (geo.system_type == SYSTEM_TYPE_BINARY)
  {
    geo.binary = TRUE;
  }

/*  The domains have been created but have not been initialized at all */

  zdom[0].coord_type = 1;
  zdom[0].ndim = 30;
  zdom[0].mdim = 30;
  zdom[0].log_linear = 0;       /* Set intervals to be logarithmic */

  zdom[1].coord_type = 1;
  zdom[1].ndim = 30;
  zdom[1].mdim = 10;
  zdom[1].log_linear = 0;       /* Set intervals to be logarithmic */


  geo.disk_z0 = geo.disk_z1 = 0.0;
  geo.adiabatic = 1;            // Default is now set so that adiabatic cooling is included in the wind
  geo.auger_ionization = 1;     //Default is on.


  geo.run_type = 0;             // Not a restart of a previous run

  geo.star_ion_spectype = geo.star_spectype
    = geo.disk_ion_spectype = geo.disk_spectype = geo.bl_ion_spectype = geo.bl_spectype = SPECTYPE_BB;
  geo.agn_ion_spectype = SPECTYPE_POW;


  geo.rstar = 7e8;
  geo.rstar_sq = geo.rstar * geo.rstar;
  geo.mstar = 0.8 * MSOL;
  geo.m_sec = 0.4 * MSOL;
  geo.period = 3.2 * 3600;
  geo.tstar = 40000;
//OLD  geo.twind_init = 40000;

  geo.ioniz_mode = IONMODE_ML93;        /* default is on the spot and find the best t */
  geo.line_mode = 3;            /* default is escape probabilites */

  geo.star_radiation = 1;       /* 1 implies star will radiate */
  geo.disk_radiation = 1;       /* 1 implies disk will radiate */
  geo.bl_radiation = 0;         /*1 implies boundary layer will radiate */
  geo.wind_radiation = 0;       /* 1 implies wind will radiate */

  geo.disk_type = DISK_FLAT;    /*1 implies existence of a disk for purposes of absorption */
  geo.diskrad = 2.4e10;
  geo.disk_mdot = 1.e-8 * MSOL / YR;

  geo.t_bl = 100000.;

  geo.pl_geometry = PL_GEOMETRY_SPHERE; // default to spherical geometry
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

/// This is to assure that we read model lists in the same order everytime
char get_spectype_oldname[LINELENGTH] = "data/kurucz91.ls";
int get_spectype_count = 0;


/**********************************************************/
/**
 * @brief      Generalized routine to get the spectrum type for a radiation source
 * and if necessary read a set of precalculated spectra
 *
 * @param [in] int  yesno  An integer used to decide whether to ask for a spectrum type
 * @param [in, out] char *  question  The query for a spectrum type for this component
 * @param [in, out] int *  spectype   The type of spectrum to assign
 * @return    the spectype, a number that says what type of spectrum (bb, power law, etc)
 * to generate for this source
 *
 * @details
 *
 * If yesno is non-zero, the user will be asked for the type of spectrum
 * for a radation source.  If yesno is 0, we presume that this radiation souce (a star
 * a disk or the wind itself) is not to radiate in the calculation, and the spectrum
 * type is set to SPECTYPE_NONE.
 *
 * If one wants to geneate spectra from a seriels of models (such as synthetic spectra
 * caluclated for stars), this routine calls routines to read the models.
 *
 * ### Notes ###
 *
 * This routine is awkwardly constructed. For reasons that are probably
 * historical the routine converts the internally stored values to different
 * values in order to ask what type of spectrum the user wants, and then
 * translates these back to the internal value.
 *
 *
 **********************************************************/

int
get_spectype (yesno, question, spectype)
     int yesno;
     char *question;
     int *spectype;
{
  char model_list[LINELENGTH];
  int stype;
  int get_models ();            // Note: Needed because get_models cannot be included in templates.h

  if (yesno)
  {
    /* XXX The next few lines convert the value of spectype to values that correspond to what is requested
     * in the rdpar statement, and then once we have the answer they are converted back to values that the
     * program uses internally.  Presumably this was done for historical reasons, having to do with running
     * old .pf files, but it is quite awkward, and we need to fix this kind of thing.  ksl
     */

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
      *spectype = SPECTYPE_BB;  // bb
    else if (stype == 2)
      *spectype = SPECTYPE_UNIFORM;     // uniform
    else if (stype == 3)
      *spectype = SPECTYPE_POW; // power law
    else if (stype == 4)
      *spectype = SPECTYPE_CL_TAB;      // broken power law
    else if (stype == 5)
      *spectype = SPECTYPE_BREM;        // bremstrahlung
    else
    {
      if (geo.run_type == RUN_TYPE_PREVIOUS)
      {                         // Continuing an old model
        strcpy (model_list, geo.model_list[get_spectype_count]);
      }
      else
      {                         // Starting a new model
        strcpy (model_list, get_spectype_oldname);
      }
      rdstr ("Input_spectra.model_file", model_list);
      get_models (model_list, 2, spectype);
      strcpy (geo.model_list[get_spectype_count], model_list);  // Copy it to geo
      strcpy (get_spectype_oldname, model_list);        // Also copy it back to the old name
      get_spectype_count++;
    }
  }
  else
  {
    *spectype = SPECTYPE_NONE;  // No radiation
  }

  return (*spectype);
}




/**********************************************************/
/**
 * @brief      simply initialises the set of
 * advanced modes stored in the modes structure to a
 * default value.
 *
 * @return   Allways returns 0
 *
 *
 * @details
 *
 * ### Notes ###
 * modes is a structure declared in python.h
 *
 * Most of the modes are initialized to 0, which means that
 * activites that could be uddertaken for this mode, will
 * not be.
 *
 * Advanced modes are turned off by default, unless the
 * -d flag is invoked at runtime.  If it is invoked,
 * one will be asked whether to turn-on an advanced mode.
 *
 * In most cases, the advanced modes cause extra diagnostic
 * information to be saved.
 *
 * see #111 and #120 for recent work
 *
 **********************************************************/

int
init_advanced_modes ()
{
  modes.iadvanced = 0;          // this is controlled by the -d flag, global mode control.
  modes.extra_diagnostics = 0;  //  when set, want to save some extra diagnostic info
  modes.save_cell_stats = 0;    // want to save photons statistics by cell
  modes.keep_ioncycle_windsaves = 0;    // want to save wind file each ionization cycle
  modes.track_resonant_scatters = 0;    // want to track resonant scatters
  modes.save_extract_photons = 0;       // we want to save details on extracted photons
  modes.adjust_grid = 0;        // the user wants to adjust the grid scale
  modes.diag_on_off = 0;        // extra diagnostics
  modes.use_debug = 0;
  modes.print_dvds_info = 0;    // print out information on velocity gradients
  modes.quit_after_inputs = 0;  // testing mode which quits after reading in inputs
  modes.fixed_temp = 0;         // do not attempt to change temperature - used for testing
  modes.zeus_connect = 0;       // connect with zeus
  modes.report_frozen = 0;

  //note write_atomicdata  is defined in atomic.h, rather than the modes structure
  write_atomicdata = 0;         // print out summary of atomic data


  modes.keep_photoabs = 1;      // keep photoabsorption in final spectrum

  return (0);
}



/**********************************************************/
/**
 * @brief      get inputs that describe the detailed spectra that
 * one intends to extract.
 *
 * @return   Always returns 0
 *
 * @details
 *
 * ### Notes ###
 *
 * There are inputs here both for a normal extaction and in
 * advanced mode for selecting spectra from photons with specific numbers
 * of scatters, or from a certain reagion.  It is also possible
 * to extact spectra in the 'live or die" mode, in which one
 * simply tracks photons that emerge in a certain inclination
 * range.
 *
 **********************************************************/

int
init_observers ()
{
  int n;
  int ichoice;
  char answer[LINELENGTH];


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

  rdpar_comment ("The minimum and maximum wavelengths in the final spectra");
  rddoub ("Spectrum.wavemin(Angstroms)", &geo.swavemin);
  rddoub ("Spectrum.wavemax(Angstroms)", &geo.swavemax);
  if (geo.swavemin > geo.swavemax)
  {
    geo.swavemax = geo.swavemin;
    geo.swavemin = geo.swavemax;
  }

  /* convert wavelengths to frequencies and store for use
     in computing macro atom and k-packet emissivities. */

  em_rnge.fmin = C / (geo.swavemax * 1.e-8);
  em_rnge.fmax = C / (geo.swavemin * 1.e-8);

  geo.matom_radiation = 0;      //initialise for ionization cycles - don't use pre-computed emissivities for macro-atom levels/ k-packets.


  rdpar_comment ("The observers and their location relative to the system");
  rdint ("Spectrum.no_observers", &geo.nangles);

  if (geo.nangles < 1 || geo.nangles > NSPEC)
  {
    Error ("no_observers %d should not be > %d or <0\n", geo.nangles, NSPEC);
    exit (0);
  }


  for (n = 0; n < geo.nangles; n++)
    rddoub ("Spectrum.angle(0=pole)", &geo.angle[n]);

  /* Phase 0 in this case corresponds to
   * an extraction direction which is in the xz plane
   */

  if (geo.system_type == SYSTEM_TYPE_BINARY)
  {

    for (n = 0; n < geo.nangles; n++)
      rddoub ("Spectrum.orbit_phase(0=inferior_conjunction)", &geo.phase[n]);
  }
  else
    Log ("No phase information needed as system type %i is not a binary\n", geo.system_type);


  strcpy (answer, "extract");
  geo.select_extract = rdchoice ("Spectrum.live_or_die(live.or.die,extract)", "0,1", answer);
  //OLD rdint ("Spectrum.live_or_die(0=live.or.die,extract=anything_else)", &geo.select_extract);
  if (geo.select_extract != 0)
  {
    geo.select_extract = 1;
    Log ("OK, extracting from specific angles\n");
  }
  else
    Log ("OK, using live or die option\n");

/* In advanced mode, select spectra with certain numbers of scatteringsi, and or other
 * characteristics
 */

  if (modes.iadvanced)
  {
    strcpy (answer, "no");
    ichoice = rdchoice ("@Spectrum.select_specific_no_of_scatters_in_spectra(y,n)", ",1,0", answer);

    //OLD strcpy (yesno, "n");
    //OLD rdstr ("@Spectrum.select_specific_no_of_scatters_in_spectra(y/n)", yesno);
    //OLD if (yesno[0] == 'y')
    if (ichoice)

    {
      Log ("OK n>MAXSCAT->all; 0<=n<MAXSCAT -> n scatters; n<0 -> >= |n| scatters\n");
      for (n = 0; n < geo.nangles; n++)
      {
        rdint ("@Spectrum.select_scatters", &geo.scat_select[n]);
      }
    }
    strcpy (answer, "no");
    ichoice = rdchoice ("@Spectrum.select_photons_by_position(y,n)", "1,0", answer);
    //OLD strcpy (yesno, "n");
    //OLD rdstr ("@Spectrum.select_photons_by_position(y/n)", yesno);
    //OLD if (yesno[0] == 'y')
    if (ichoice)
    {
      //OLD Log ("OK 0->all; -1 -> below; 1 -> above the disk, 2 -> specific location in wind\n");
      for (n = 0; n < geo.nangles; n++)
      {
        strcpy (answer, "all");
        geo.top_bot_select[n] = rdchoice ("@Spectrum.select_location(all,below_disk,above_disk,spherical_region)", "0,-1,1,2", answer);


        //OLD rdint ("@Spectrum.select_location", &geo.top_bot_select[n]);
        if (geo.top_bot_select[n] == 2)
        {
          Log ("Warning: Make sure that position will be in wind, or no joy will be obtained\n");
          rddoub ("@Spectrum.select_rho(cm)", &geo.rho_select[n]);
          rddoub ("@Spectrum.select_z(cm)", &geo.z_select[n]);
          rddoub ("@Spectrum.select_azimuth(deg)", &geo.az_select[n]);
          rddoub ("@Spectrum.select_r(cm)", &geo.r_select[n]);

        }
      }
    }
  }

  /* Select the units of the output spectra.  This is always needed.
   * There are 3 basics choices flambda, fnu, and the internal units
   * of the program.  The first two imply that output units are scaled
   * to a distance of 100 pc. The internal units are basically a luminosity
   * within a wavelength/frequency interval. */

  strcpy (answer, "flambda");
  geo.select_spectype = rdchoice ("Spectrum.type(flambda,fnu,basic)", "1,2,3", answer);
  //OLD rdint ("Spectrum.type(flambda(1),fnu(2),basic(other)", &geo.select_spectype);

  if (geo.select_spectype == 1)
  {
    Log ("OK, generating flambda at 100pc\n");
    geo.select_spectype = SPECTYPE_FLAMBDA;
  }
  else if (geo.select_spectype == 2)
  {
    Log ("OK, generating fnu at 100 pc\n");
    geo.select_spectype = SPECTYPE_FNU;
  }
  else
  {
    Log ("OK, basic Monte Carlo spectrum\n");
    geo.select_spectype = SPECTYPE_RAW;
  }

  return (0);
}


/**********************************************************/
/**
 * @brief      gets information about the number of
 * 	cycles and how many photons there should be per cycle.
 *
 * @return     Generally returns 0
 *
 * @details
 * ??? DESCRIPTION ???
 *
 * ### Notes ###
 * The routine also allocates memory for the photon structure.
 * If the routine is unable to allocate this membory, the routine
 * will exit.
 *
 **********************************************************/

PhotPtr
init_photons ()
{
  PhotPtr p;
  double x;

  /* Although Photons_per_cycle is really an integer,
     read in as a double so it is easier for input */

  x = 100000;
  rddoub ("Photons_per_cycle", &x);
  NPHOT = x;                    // NPHOT is photons/cycle

#ifdef MPI_ON
  Log ("Photons per cycle per MPI task will be %d\n", NPHOT / np_mpi_global);

  NPHOT /= np_mpi_global;
#endif

  rdint ("Ionization_cycles", &geo.wcycles);

  rdint ("Spectrum_cycles", &geo.pcycles);


  if (geo.wcycles == 0 && geo.pcycles == 0)
  {
    Log ("Both ionization and spectral cycles are set to 0; There is nothing to do so exiting\n");
    exit (0);                   //There is really nothing to do!
  }

  /* Allocate the memory for the photon structure now that NPHOT is established */

  photmain = p = (PhotPtr) calloc (sizeof (p_dummy), NPHOT);

  if (p == NULL)
  {
    Error ("init_photons: There is a problem in allocating memory for the photon structure\n");
    exit (0);
  }
  else
  {
    /* large photon numbers can cause problems / runs to crash. Report to use (see #209) */
    Log
      ("Allocated %10d bytes for each of %5d elements of photon structure totaling %10.1f Mb \n",
       sizeof (p_dummy), NPHOT, 1.e-6 * NPHOT * sizeof (p_dummy));
    if ((NPHOT * sizeof (p_dummy)) > 1e9)
      Error ("Over 1 GIGABYTE of photon structure allocated. Could cause serious problems.\n");
  }


  return (p);
}





/**********************************************************/
/**
 * @brief      Select the ionization and line transfer modes, along
 * with various other parameters about ionization
 *
 * @return   Always returns 0
 *
 * @details
 * The routine queries for a variety of parameters which describe
 * how ionization and line transfer are handled.  This includes
 * whether surfaces reflect
 *
 * ### Notes ###
 *
 * The routine is not particular well named, and it is not
 * immediately obvious that all of the parmenters that are queried
 * for here are related.
 *
 **********************************************************/

int
init_ionization ()
{
  int thermal_opt;
  char answer[LINELENGTH];



  strcpy (answer, "matrix_bb");
  geo.ioniz_mode = rdchoice ("Wind_ionization(on.the.spot,ML93,LTE_tr,LTE_te,fixed,matrix_bb,matrix_pow)", "0,3,1,4,2,8,9", answer);

  if (geo.ioniz_mode == IONMODE_FIXED)
  {
    rdstr ("wind.fixed_concentrations_file", &geo.fixed_con_file[0]);
  }



  /*Normally, geo.partition_mode is set to -1, which means that partition functions are calculated to take
     full advantage of the data file.  This means that in calculating the partition functions, the information
     on levels and their multiplicities is taken into account.   */

  geo.partition_mode = -1;      //?? Stuart, is there a reason not to move this earlier so it does not affect restart


  /* get_line_transfer_mode reads in the Line_transfer question from the user,
     then alters the variables geo.line_mode, geo.scatter_mode, geo.rt_mode and geo.macro_simple.
     This is fairly involved and so is a separate routine  */

  get_line_transfer_mode ();


  thermal_opt = 0;              /* NSH 131213 Set the option to zero - the default. The lines allow allow the
                                   user to turn off mechanisms that affect the thermal balance. Adiabatic is the only one implemented
                                   to start off with. */

  strcpy (answer, "reflect");
  geo.absorb_reflect = rdchoice ("Surface.reflection.or.absorption(reflect,absorb,thermalized.rerad)", "1,0,2", answer);
  //OLD rdint ("Surface.reflection.or.absorption(0=no.rerad,1=high.albedo,2=thermalized.rerad)", &geo.absorb_reflect);

  // XXX Next line needs to be re-written to use rdchooce when FU ORi branch is merged
  rdint ("Thermal_balance_options(0=everything.on,1=no.adiabatic)", &thermal_opt);

  if (thermal_opt == 1)
  {
    geo.adiabatic = 0;
  }

  else if (thermal_opt > 1 || thermal_opt < 0)
  {
    Error ("Unknown thermal balance mode %d\n", thermal_opt);
    exit (0);
  }


  /* Prevent bf calculation of macro_estimaters when no macro atoms are present.   */

  if (nlevels_macro == 0)
    geo.macro_simple = 1;       // Make everything simple if no macro atoms -- 57h

  /* initialise the choice of handling for macro pops. */
  if (geo.run_type == RUN_TYPE_PREVIOUS)
  {
    geo.macro_ioniz_mode = 1;   // Now that macro atom properties are available for restarts
  }
  else
  {
    geo.macro_ioniz_mode = 0;
  }

  return (0);

}





/**********************************************************/
/**
 * @brief      Sets a global "push_through_distance" designed
 * to prevent photons from being trapped at boundaries
 *
 * @return     dfudge 	the push through distance
 *
 * @details
 *
 * dfudge is intended to be a small positive number that prevents
 * photons from being trapped at the edges of cells or wind cones.
 * It is often added to the distance to an "edge" to push the photon
 * thorough whatever boundary is being calculated.
 *
 * ### Notes ###
 *
 * setup_dfudge is used to modify the variable DFUDGE which
 * can be found in python.h.
 *
 * Originally, DFUDGE was the only number used to push through
 * boundariies, but today DFUDGE is used sparingly, if at all,
 * as push through distances within wind cells are defined
 * differently for each cell.
 *
 * There are two competing factors in defining DFUDGE.  It
 * should be short enough so that the push through distane
 * goes only a small way into s cell.  It should be large
 * enough though that round-off errors do not prevent one
 * from actually getting into a cell.
 *
 **********************************************************/

double
setup_dfudge ()
{
  double dfudge;
  double delta;
  double rmin;
  int ndom;

  rmin = VERY_BIG;
  for (ndom = 0; ndom < geo.ndomain; ndom++)
  {
    if (rmin > zdom[ndom].rmin)
    {
      rmin = zdom[ndom].rmin;
    }
  }



  delta = geo.rmax - rmin;

  if (delta < 1.e8)
  {
    dfudge = (geo.rmax - rmin) / 1000.0;
  }
  else if (delta < 1e15)
  {
    dfudge = 1e5;
  }
  else
  {
    dfudge = geo.rmax / 1.e10;
  }

  Log ("DFUDGE set to %e based on geo.rmax\n", dfudge);

  return (dfudge);
}
