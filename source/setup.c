#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "atomic.h"
#include "python.h"


/***********************************************************
             University of Southampton

Synopsis: 
  get_stellar_params sets rstar, mstar, tstar as well
  as secondary parameters based on user inputs
   
Arguments:		

Returns:
 
 
Description:	

Notes:

History:
	1502  JM 	Moved here from main()

**************************************************************/

double
get_stellar_params ()
{

  /* Describe the basic binary star system */

  geo.mstar /= MSOL;		// Convert to MSOL for ease of data entry
  rddoub ("mstar(msol)", &geo.mstar);
  geo.mstar *= MSOL;

  /* If a BH we want geo.rstar to be at least as large as the last stable orbit for
   * a non-rotating BH
   */

  if (geo.system_type == SYSTEM_TYPE_AGN)
    {
      geo.rstar = 6. * G * geo.mstar / (C * C);	//correction - ISCO is 6x Rg NSH 121025
    }

  rddoub ("rstar(cm)", &geo.rstar);


  geo.r_agn = geo.rstar;	/* At present just set geo.r_agn to geo.rstar */
  geo.rstar_sq = geo.rstar * geo.rstar;
  if (geo.system_type != SYSTEM_TYPE_AGN)
    {
      rdint ("Star_radiation(y=1)", &geo.star_radiation);
      get_spectype (geo.star_radiation,
		    "Rad_type_for_star(0=bb,1=models)_to_make_wind",
		    &geo.star_ion_spectype);

      if (geo.star_radiation)
	rddoub ("tstar", &geo.tstar_init);
    }
  else
    {
      geo.star_radiation = 0;
      geo.tstar_init = 0;
    }

  /* tstar_init and lum_star_init refer to values without the effects of backscattering */

  geo.tstar = geo.tstar_init;

  geo.lum_star = geo.lum_star_init =
    4 * PI * geo.rstar * geo.rstar * STEFAN_BOLTZMANN * pow (geo.tstar, 4.);


  /* Describe the secondary if that is required */

  if (geo.system_type == SYSTEM_TYPE_BINARY)	/* It's a binary system */
    {

      geo.m_sec /= MSOL;	// Convert units for ease of data entry
      rddoub ("msec(msol)", &geo.m_sec);
      geo.m_sec *= MSOL;

      geo.period /= 3600.;	// Convert units to hours for easy of data entry
      rddoub ("period(hr)", &geo.period);
      geo.period *= 3600.;	// Put back to cgs immediately                   
    }

  return (geo.lum_star_init);
}

/***********************************************************
             University of Southampton

Synopsis: 
  get_disk_params sets up the disk parameters according to user inputs, 
  e.g. the temperature profile, accretion rate etc.
   
Arguments:		

Returns:
  disk_illum - this is used by python.c and so needs to be returned
 
Description:	

Notes:

History:
	1502  JM 	Moved here from main()
	1510	ksl	Modified to restore illumination
			options, which were brokedn

**************************************************************/


double
get_disk_params ()
{
  geo.disk_mdot /= (MSOL / YR);	// Convert to msol/yr to simplify input
  rddoub ("disk.mdot(msol/yr)", &geo.disk_mdot);
  geo.disk_mdot *= (MSOL / YR);
  rdint ("Disk.temperature.profile(0=standard;1=readin,2=analytic)",
	 &geo.disk_tprofile);
  if (geo.disk_tprofile == DISK_TPROFILE_READIN)
    {
      rdstr ("T_profile_file", files.tprofile);
    }

  /* Set a default for diskrad for an AGN */
  if (geo.system_type == SYSTEM_TYPE_AGN)
    {
      geo.diskrad = 100. * geo.r_agn;
    }

  rddoub ("disk.radmax(cm)", &geo.diskrad);
  Log ("geo.diskrad  %e\n", geo.diskrad);

  geo.diskrad_sq = geo.diskrad * geo.diskrad;

/* If diskrad <= geo.rstar set geo.disk_type = DISK_NONE to make any disk transparent anyway. */

  if (geo.diskrad < geo.rstar)
    {
      Log ("Disk radius is less than star radius, so assuming no disk)\n");
      geo.disk_type = DISK_NONE;
    }

  if (geo.disk_type == DISK_VERTICALLY_EXTENDED)
    {				/* Get the additional variables need to describe a vertically extended disk */
      rddoub ("disk.z0(fractional.height.at.diskrad)", &geo.disk_z0);
      rddoub ("disk.z1(powerlaw.index)", &geo.disk_z1);
    }
  return (0);
}



/***********************************************************
             University of Southampton

Synopsis: 
  get_bl_and_agn_params sets up the boundary layer and agn power law parameters
  based on user input and system type
   
Arguments:		
  lstar     double 
            star luminosity as calculated by get_stellar_params
Returns:
 
Description:	

Notes:

History:
	1502  JM 	Moved here from main()

**************************************************************/

int
get_bl_and_agn_params (lstar)
     double lstar;
{
  double xbl;
  double temp_const_agn;

  rdpar_comment ("Parameters for BL or AGN");

  if (geo.system_type == SYSTEM_TYPE_AGN)	/* If it is an AGN */
    {
      geo.star_radiation = 0;	// 70b - AGN do not have a star at the center */
      geo.bl_radiation = 0;
      geo.agn_radiation = 1;
      rdint ("QSO_BH_radiation(y=1)", &geo.agn_radiation);
    }
  else
    {
      rdint ("Boundary_layer_radiation(y=1)", &geo.bl_radiation);
      geo.agn_radiation = 0;	// So far at least, our star systems don't have a BH
    }


  get_spectype (geo.bl_radiation,
		"Rad_type_for_bl(0=bb,1=models,3=pow)_to_make_wind",
		&geo.bl_ion_spectype);

  get_spectype (geo.agn_radiation,
		"Rad_type_for_agn(0=bb,1=models,3=power_law,4=cloudy_table,5=bremsstrahlung)_to_make_wind",
		&geo.agn_ion_spectype);

  Log("geo.agn_ion_spectype %d\n", geo.agn_ion_spectype);

  /* 130621 - ksl - This is a kluge to add a power law to stellar systems.  What id done
     is to remove the bl emission, which we always assume to some kind of temperature
     driven source, and replace it with a power law source

     Note that the next 3 or 4 lines just tell you that there is supposed to be a power
     law source.  They don't teel you what the parameters are.
   */

  if (geo.bl_ion_spectype == SPECTYPE_POW)
    {
      geo.agn_radiation = 1;
      geo.agn_ion_spectype = SPECTYPE_POW;
      geo.bl_radiation = 0;
      Log ("Trying to make a start with a power law boundary layer\n");
    }
  else
    {
      Log ("Not Trying to make a start with a power law boundary layer %d\n",
	   geo.bl_ion_spectype);
    }


  /* Describe the boundary layer */

  if (geo.bl_radiation && geo.bl_ion_spectype != SPECTYPE_POW)
    {
      xbl = geo.lum_bl = 0.5 * G * geo.mstar * geo.disk_mdot / geo.rstar;

      rddoub ("lum_bl(ergs/s)", &geo.lum_bl);
      Log ("OK, the bl lum will be about %.2e the disk lum\n",
	          geo.lum_bl / xbl);
      rddoub ("t_bl", &geo.t_bl);
    }
  else
    {
      geo.lum_bl = 0;
      geo.t_bl = 0;
    }

  /* Describe the agn */

  if (geo.agn_radiation && geo.system_type == SYSTEM_TYPE_AGN)	/* This peculiar line is to enamble us to add a star with a power law component */
    {
      xbl = geo.lum_agn = 0.5 * G * geo.mstar * geo.disk_mdot / geo.r_agn;

      /* If there is no disk, initilize geo.lum to the luminosity of a star */
      if (geo.disk_type == DISK_NONE)
	{
	  geo.lum_agn = lstar;
	}

      // At present we have set geo.r_agn = geo.rstar, and encouraged the user
      // set the default for the radius of the BH to be 6 R_Schwartschild.
      // rddoub("R_agn(cm)",&geo.r_agn);

      /* if we have a "blackbody agn" the luminosity is set by Stefan Boltzmann law
         once the AGN blackbody temp is read in, otherwise set by user */
      if (geo.agn_ion_spectype != SPECTYPE_BB && geo.agn_ion_spectype < 0)
	rddoub ("lum_agn(ergs/s)", &geo.lum_agn);


      Log ("OK, the agn lum will be about %.2e the disk lum\n",
	   geo.lum_agn / xbl);
      if (geo.agn_ion_spectype == SPECTYPE_POW
	  || geo.agn_ion_spectype == SPECTYPE_CL_TAB)
	{
	  geo.alpha_agn = (-1.5);
	  rddoub ("agn_power_law_index", &geo.alpha_agn);

	  if (geo.alpha_agn == -1.0)	//deal with the pathological case
	    {
	      geo.const_agn = geo.lum_agn / (log (2.42e18) - log (4.84e17));
	    }
	  else
	    {
	      geo.const_agn =
		geo.lum_agn /
		(((pow (2.42e18, geo.alpha_agn + 1.)) -
		  pow (4.84e17,
		       geo.alpha_agn + 1.0)) / (geo.alpha_agn + 1.0));
	    }
	  Log ("AGN Input parameters give a power law constant of %e\n",
	       geo.const_agn);
	}
      else if (geo.agn_ion_spectype == SPECTYPE_BREM)
	{

	  geo.brem_temp = 1.16e8;	//10kev
	  geo.brem_alpha = -0.2;	//This is the cloudy form of bremstrahlung
	  geo.const_agn = 1.0;
	  rddoub ("agn_bremsstrahlung_temp(K)", &geo.brem_temp);
	  rddoub ("agn_bremsstrahlung_alpha", &geo.brem_alpha);
	  temp_const_agn =
	    geo.lum_agn / qromb (integ_brem, 4.84e17, 2.42e18, 1e-4);
	  geo.const_agn = temp_const_agn;
	  Log ("AGN Input parameters give a Bremsstrahlung constant of %e\n",
	       temp_const_agn);

	}
      else if (geo.agn_ion_spectype == SPECTYPE_BB)
	{
	  /* note that alpha_agn holds the temperature in the case of "blackbody agn" */
	  rddoub ("agn_blackbody_temp(K)", &geo.alpha_agn);
	  geo.lum_agn =
	    4 * PI * geo.r_agn * geo.r_agn * STEFAN_BOLTZMANN *
	    pow (geo.alpha_agn, 4.);
	}
      else if (geo.agn_ion_spectype >= 0) {
        geo.lum_agn = emittance_continuum (geo.agn_ion_spectype, 1e14, 1e20, 1, 1, LUMINOSITY_UNITS);
        Log("GEO LUM AGN %8.4e\n\n", geo.lum_agn);
      }

      /* JM 1502 -- lines to add a low frequency power law cutoff. accessible
         only in advanced mode and for non broken power law. 
         default is zero which is checked before we call photo_gen_agn */
      geo.pl_low_cutoff = 0.0;
      if (modes.iadvanced && (geo.agn_ion_spectype == SPECTYPE_POW))
	rddoub ("@agn_power_law_cutoff", &geo.pl_low_cutoff);

      rdint ("geometry_for_pl_source(0=sphere,1=lamp_post)",
	     &geo.pl_geometry);

      if (geo.pl_geometry == PL_GEOMETRY_LAMP_POST)
	{
	  rddoub ("lamp_post.height(r_g)", &geo.lamp_post_height);
	  geo.lamp_post_height *= G * geo.mstar / C / C;	//get it in CGS units 
	  Log ("lamp_post_height is cm is %g\n", geo.lamp_post_height);
	}
      else if (geo.pl_geometry != PL_GEOMETRY_SPHERE)	// only two options at the moment
	{
	  Error ("Did not understand power law geometry %i. Fatal.\n",
		 geo.pl_geometry);
	  exit (0);
	}



      /* Computes the constant for the power law spectrum from the input alpha and 2-10 luminosity. 
         This is only used in the sim correction factor for the first time through. 
         Afterwards, the photons are used to compute the sim parameters. */



      if (geo.agn_ion_spectype == SPECTYPE_CL_TAB)	/*NSH 0412 - option added to allow direct comparison with cloudy power law table option */
	{
	  geo.agn_cltab_low = 1.0;
	  geo.agn_cltab_hi = 10000;
	  rddoub ("low_energy_break(ev)", &geo.agn_cltab_low);	/*lo frequency break - in ev */
	  rddoub ("high_energy_break(ev)", &geo.agn_cltab_hi);
	  geo.agn_cltab_low_alpha = 2.5;	//this is the default value in cloudy
	  geo.agn_cltab_hi_alpha = -2.0;	//this is the default value in cloudy
	}
    }
  else if (geo.agn_radiation)	/* We want to add a power law to something other than an AGN */
    {
      xbl = geo.lum_agn = 0.5 * G * geo.mstar * geo.disk_mdot / geo.r_agn;

      // At present we have set geo.r_agn = geo.rstar, and encouraged the user
      // set the default for the radius of the BH to be 6 R_Schwartschild.
      // rddoub("R_agn(cm)",&geo.r_agn);

      rddoub ("lum_agn(ergs/s)", &geo.lum_agn);
      Log ("OK, the agn lum will be about %.2e the disk lum\n",
	   geo.lum_agn / xbl);
      geo.alpha_agn = (-1.5);
      rddoub ("agn_power_law_index", &geo.alpha_agn);

      /* JM 1502 -- lines to add a low frequency power law cutoff. accessible
         only in advanced mode. default is zero which is checked before we call photo_gen_agn */
      geo.pl_low_cutoff = 0.0;
      if (modes.iadvanced)
	rddoub ("@agn_power_law_cutoff", &geo.pl_low_cutoff);


      /* Computes the constant for the power law spectrum from the input alpha and 2-10 luminosity. 
         This is only used in the sim correction factor for the first time through. 
         Afterwards, the photons are used to compute the sim parameters. */


      if (geo.alpha_agn == -1.0)	//deal with the pathological case
	{
	  geo.const_agn = geo.lum_agn / (log (2.42e18) - log (4.84e17));
	}
      else
	{
	  geo.const_agn =
	    geo.lum_agn /
	    (((pow (2.42e18, geo.alpha_agn + 1.)) -
	      pow (4.84e17, geo.alpha_agn + 1.0)) / (geo.alpha_agn + 1.0));
	}


      Log ("AGN Input parameters give a power law constant of %e\n",
	   geo.const_agn);

      if (geo.agn_ion_spectype == SPECTYPE_CL_TAB)	/*NSH 0412 - option added to allow direct comparison with cloudy power law table option */
	{
	  geo.agn_cltab_low = 1.0;
	  geo.agn_cltab_hi = 10000;
	  rddoub ("low_energy_break(ev)", &geo.agn_cltab_low);	/*lo frequency break - in ev */
	  rddoub ("high_energy_break(ev)", &geo.agn_cltab_hi);
	  geo.agn_cltab_low_alpha = 2.5;	//this is the default value in cloudy
	  geo.agn_cltab_hi_alpha = -2.0;	//this is the default value in cloudy
	}
    }

  else
    {
      geo.r_agn = 0.0;
      geo.lum_agn = 0.0;
      geo.alpha_agn = 0.0;
      geo.const_agn = 0.0;
    }
  return (0);
}





/***********************************************************
             University of Southampton
Synopsis: 
  get_meta_params reads in data pertaining to simulation meta-
  properties like reverberation mapping settings and variance
  reduction techniques.
   
Arguments:    
Returns:
 
Description:  
Notes:
History:
  1504  SWM   Added
**************************************************************/
int
get_meta_params (void)
{

  int meta_param, i, j, k, z, istate, levl, levu;
  char trackline[LINELENGTH];

  meta_param = 0;		// initialize to no reverberation tracking
  rdint ("reverb.type", &meta_param);
  switch (meta_param)
    {				//Read in reverb tyoe, if any
    case 0:
      geo.reverb = REV_NONE;
      break;
    case 1:
      geo.reverb = REV_PHOTON;
      break;
    case 2:
      geo.reverb = REV_WIND;
      break;
    case 3:
      geo.reverb = REV_MATOM;
      break;
    default:
      Error ("reverb.type: Invalid reverb mode.\n \
      Valid modes are 0=None, 1=Photon, 2=Wind, 3=Macro-atom.\n");
    }

  // ========== DEAL WITH DISK SETTINGS ==========
  if (geo.disk_type > 0 && geo.reverb != REV_NONE)
    {
      rdint ("reverb.disk_type", &meta_param);
      switch (meta_param)
	{			//Read in reverb tyoe, if any
	case 0:
	  geo.reverb_disk = REV_DISK_CORRELATED;
	  break;
	case 1:
	  geo.reverb_disk = REV_DISK_UNCORRELATED;
	  break;
	case 2:
	  geo.reverb_disk = REV_DISK_IGNORE;
	  break;
	default:
	  Error ("reverb.disk_type: Invalid reverb disk mode.\n \
        Valid modes are 0=Correlated with central source, 1=Uncorrelated, 2=Ignore.\n");
	}
    }

  // ========== DEAL WITH VISUALISATION SETTINGS ==========
  if (geo.reverb == REV_WIND || geo.reverb == REV_MATOM)
    {				//If this requires further parameters, set defaults
      geo.reverb_lines = 0;
      geo.reverb_path_bins = 100;
      geo.reverb_angle_bins = 100;
      geo.reverb_dump_cells = 0;
      geo.reverb_vis = REV_VIS_NONE;
      rdint ("reverb.path_bins", &geo.reverb_path_bins);
      rdint ("reverb.visualisation", &meta_param);
      switch (meta_param)
	{			//Select whether to produce 3d visualisation file and/or dump flat csvs of spread in cells
	case 0:
	  geo.reverb_vis = REV_VIS_NONE;
	  break;
	case 1:
	  geo.reverb_vis = REV_VIS_VTK;
	  break;
	case 2:
	  geo.reverb_vis = REV_VIS_DUMP;
	  break;
	case 3:
	  geo.reverb_vis = REV_VIS_BOTH;
	  break;
	default:
	  Error ("reverb.visualisation: Invalid mode.\n \
        Valid modes are 0=None, 1=VTK, 2=Cell dump, 3=Both.\n");
	}

      if (geo.reverb_vis == REV_VIS_VTK || geo.reverb_vis == REV_VIS_BOTH)
	//If we're producing a 3d visualisation, select bins. This is just for aesthetics
	rdint ("reverb.angle_bins", &geo.reverb_angle_bins);
      if (geo.reverb_vis == REV_VIS_DUMP || geo.reverb_vis == REV_VIS_BOTH)
	{			//If we;re dumping path arrays, read in the number of cells to dump them for
	  rdint ("reverb.dump_cells", &geo.reverb_dump_cells);
	  geo.reverb_dump_cell_x =
	    (double *) calloc (geo.reverb_dump_cells, sizeof (double));
	  geo.reverb_dump_cell_z =
	    (double *) calloc (geo.reverb_dump_cells, sizeof (double));
	  geo.reverb_dump_cell =
	    (int *) calloc (geo.reverb_dump_cells, sizeof (int));
	  for (k = 0; k < geo.reverb_dump_cells; k++)
	    {			//For each we expect, read a paired cell coord as "[i]:[j]". May need to use py_wind to find indexes.
	      rdline ("reverb.dump_cell", trackline);
	      if (sscanf
		  (trackline, "%lf:%lf", &geo.reverb_dump_cell_x[k],
		   &geo.reverb_dump_cell_z[k]) == EOF)
		{		//If this line is malformed, warn the user
		  Error ("reverb.dump_cell: Invalid position line '%s'\n \
            Expected format '[x]:[z]'\n", trackline);
		  exit (0);
		}
	    }
	}
    }

  // ========== DEAL WITH MATOM LINES ==========
  if (geo.reverb == REV_MATOM)
    {				//If this is macro-atom mode
      if (geo.rt_mode != RT_MODE_MACRO)
	{			//But we're not actually working in matom mode...
	  Error ("reverb.type: Invalid reverb mode.\n \
      Macro-atom mode selected but macro-atom scattering not on.\n");
	  exit (0);
	}

      //Read in the number of lines to be tracked and allocate space for them
      rdint ("reverb.matom_lines", &geo.reverb_lines);
      geo.reverb_line = (int *) calloc (geo.reverb_lines, sizeof (int));
      if (geo.reverb_lines < 1)
	{			//If this is <1, then warn the user and quit
	  Error ("reverb.matom_lines: \
      Must specify 1 or more lines to watch in macro-atom mode.\n");
	  exit (0);
	}

      for (i = 0; i < geo.reverb_lines; i++)
	{			//Finally, for each line we expect, read it in
	  rdline ("reverb.matom_line", trackline);
	  if (sscanf (trackline, "%d:%d:%d:%d", &z, &istate, &levu, &levl) ==
	      EOF)
	    {			//If this line is malformed, warn the user
	      Error ("reverb.matom_line: Malformed line '%s'\n \
          Expected format '[z]:[istate]:[upper level]:[lower level]'\n", trackline);
	      exit (0);
	    }
	  else
	    {			//Otherwise, sift through the line list to find what this transition corresponds to
	      for (j = 0; j < nlines_macro; j++)
		{		//And record the line position in geo for comparison purposes
		  if (line[j].z == z && line[j].istate == istate
		      && line[j].levu == levu && line[j].levl == levl)
		    {		//We're matching z, ionisation state, and upper and lower level transitions
		      geo.reverb_line[i] = line[j].where_in_list;
		    }
		}
	    }
	}
    }
  else if (geo.reverb == REV_WIND)
    {				//For wind mode...
      if (geo.wind_radiation == 0)
	{			//Warn if this data is being gathered but not used (can be useful for debug)
	  Error
	    ("reverb.type: Wind radiation is off but wind-based path tracking is enabled!\n");
	}
    }

  // ========== DEAL WITH LINE CULLING ==========
  if (geo.reverb != REV_NONE)
    {
      //Should we filter any lines out?
      //If -1, blacklist continuum, if >0 specify lines as above and whitelist
      //Automatically include matom_lines
      rdint ("reverb.filter_lines", &geo.reverb_filter_lines);
      if (geo.reverb_filter_lines > 0)
	{			//If we're given a whitelist, allocate temp storage (up to 256 lines!)
	  int temp[256], bFound;
	  for (i = 0; i < geo.reverb_filter_lines; i++)
	    {			//For each provided line, read in
	      rdint ("reverb.filter_line", &temp[i]);
	    }
	  if (geo.reverb == REV_MATOM)
	    {			//If we're in matom mode, check if those lines have already been included
	      for (i = 0; i < geo.reverb_lines; i++)
		{		//For each matom line
		  bFound = 0;
		  for (j = 0; j < geo.reverb_filter_lines; j++)
		    {		//Check if it's in the filter list
		      if (geo.reverb_line[i] == temp[j])
			bFound = 1;
		    }
		  if (!bFound)
		    {		//If it's not, add it to the filter list and increment the total lines
		      temp[geo.reverb_filter_lines++] = geo.reverb_line[i];
		    }
		}
	    }
	  //Allocate enough space for the filter list
	  geo.reverb_filter_line =
	    calloc (geo.reverb_filter_lines, sizeof (int));
	  for (i = 0; i < geo.reverb_filter_lines; i++)
	    {			//Populate the filter list from the temp list
	      geo.reverb_filter_line[i] = temp[i];
	    }
	}
    }
  return (0);
}



/***********************************************************
             University of Southampton

Synopsis: 
  get_standard_care_factors provides more control over how the program is
  run
   
Arguments:    

Returns:

Description:  

Notes:

History:
  1502  JM  Moved here from main()

**************************************************************/

int
get_standard_care_factors ()
{
  int istandard;
  istandard = 1;
  SMAX_FRAC = 0.5;
  DENSITY_PHOT_MIN = 1.e-10;

  /* 141116 - ksl - Made care factors and advanced command as this is clearly somethng that is diagnostic */

  if (modes.iadvanced)
    {
      rdint ("@Use.standard.care.factors(1=yes)", &istandard);

      if (!istandard)
	{
	  rddoub ("@Fractional.distance.photon.may.travel", &SMAX_FRAC);
	  rddoub ("@Lowest.ion.density.contributing.to.photoabsorption",
		  &DENSITY_PHOT_MIN);
	  rdint ("@Keep.photoabs.during.final.spectrum(1=yes)",
		 &modes.keep_photoabs);
	}
    }
  return (0);
}
