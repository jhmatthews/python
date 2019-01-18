
/***********************************************************/
/** @file   setup_domains.c
 * @author ksl
 * @date   January, 2018
 * @brief  Get the parameters needed to describe winds
 *
 * File containing several routines that collectively
 * define the components to a wind in python. Each component
 * of the wind is said to be a domain, and the information
 * for each domain is stored in the elements of zdom
 ***********************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "atomic.h"
#include "python.h"



/**********************************************************/
/** 
 * @brief       Get inputs that describe a particular component of the wind
 *
 * @param [in] ndom  The number (begining with 0) of this particular domain
 * @return  0 
 *
 * Sets up the one domain, which includes defining wind type, e.g whether
 * it is a shell, or a biconical flow, or an imported model, as well
 * as the type of coordinate system and its simensions.
 *
 * If the wind is to be imported from a file, it is imported in this
 * routine.
 *
 * ###Notes###
 * 1801 -   Refactored into this file in 1801.  Updates in the
 *          fall of 17 were made to allow for importing models.
 *          Note that cyl_var coordinates are not currently 
 *          working
***********************************************************/


int
get_domain_params (ndom)
     int ndom;
{
  int input_int;

  if (ndom >= geo.ndomain)
  {
    Error ("Trying to get grid params for a non-existent domain!\n");
    exit (0);
  }


  rdint ("Wind_type(0=SV,1=Star,3=Hydro,4=corona,5=knigge,6=homologous,7=yso,9=shell,10=imported)", &zdom[ndom].wind_type);

  if (zdom[ndom].wind_type == 2)
  {
    Error ("Wind_type 2, which was used to read in a previous model is no longer allowed! Use System_type instead!\n");
    exit (0);
  }


  strcat (zdom[ndom].name, "Wind");


  input_int = 1;


  /* Define the coordinate system for the grid and allocate memory for the wind structure */
  rdint ("Wind.coord_system(0=spherical,1=cylindrical,2=spherical_polar,3=cyl_var)", &input_int);
  switch (input_int)
  {
  case 0:
    zdom[ndom].coord_type = SPHERICAL;
    break;
  case 1:
    zdom[ndom].coord_type = CYLIND;
    break;
  case 2:
    zdom[ndom].coord_type = RTHETA;
    break;
  case 3:
    zdom[ndom].coord_type = CYLVAR;
    break;
  default:
    Error ("Invalid parameter supplied for 'Coord_system'. Valid coordinate types are: \n\
          0 = Spherical, 1 = Cylindrical, 2 = Spherical polar, 3 = Cylindrical (varying Z)");
  }


  if (zdom[ndom].wind_type == IMPORT)
  {
    import_wind (ndom);
  }
  else
  {
    rdint ("Wind.dim.in.x_or_r.direction", &zdom[ndom].ndim);
    if (zdom[ndom].coord_type)
    {
      rdint ("Wind.dim.in.z_or_theta.direction", &zdom[ndom].mdim);
      if (zdom[ndom].mdim < 4)
      {
        Error ("python: domain mdim must be at least 4 to allow for boundaries\n");
        exit (0);
      }
    }
    else
      zdom[ndom].mdim = 1;

  }

/* Check that NDIM_MAX is greater than NDIM and MDIM.  */

  if ((zdom[ndom].ndim > NDIM_MAX) || (zdom[ndom].mdim > NDIM_MAX))
  {
    Error ("NDIM_MAX %d is less than NDIM %d or MDIM %d. Fix in python.h and recompile\n", NDIM_MAX, zdom[ndom].ndim, zdom[ndom].mdim);
    exit (0);
  }


  /* If we are in advanced then allow the user to modify scale lengths */
  if (modes.iadvanced)
  {
    rdint ("@Diag.adjust_grid(0=no,1=yes)", &modes.adjust_grid);

    if (modes.adjust_grid)
    {
      Log ("You have opted to adjust the grid scale lengths\n");
      rddoub ("@geo.xlog_scale", &zdom[ndom].xlog_scale);
      if (zdom[ndom].coord_type != SPHERICAL)
        rddoub ("@geo.zlog_scale", &zdom[ndom].zlog_scale);
    }
  }

  zdom[ndom].ndim2 = zdom[ndom].ndim * zdom[ndom].mdim;


  return (0);
}




/**********************************************************/
/**   
 * @brief       Get detailed particular component of the wind
 *
 * @param [in] ndom  The number (begining with 0) of this particular domain
 * @return  0 
 *
 * Continues the setup of a single domain begun in get_domain_params.    
 *
 * Much of this routine is a steering routine that calls other
 * subroutines depending on the type of wind, e.g sv, for this
 * particular wind domain.
 *
 *
 * ###Notes###
 * 1801 -   Refactored into this file in 1801.  

***********************************************************/

int
get_wind_params (ndom)
     int ndom;
{
  // XXX These need to be initalized sensibly and 
  // it is not obvious that is happenning

  zdom[ndom].rmax = 1e12;

  if (geo.system_type == SYSTEM_TYPE_AGN)
  {
    zdom[ndom].rmax = 50. * geo.r_agn;
  }



  /* XXX - This should be part of the individual get_wind_parameters, not here 
   * Note that one cannot simply move the next statemnts to below the specify the various wind
   * of domains, for reasons that are unclear.  If you do this sn_1d will fail.*/

  rddoub ("Wind.radmax(cm)", &zdom[ndom].rmax);

  if (zdom[ndom].rmax <= zdom[ndom].rmin)
  {
    Error ("get_wind_parameters: rmax (%10.4e) less than or equal to rmin %10.4e in domain %d\n", zdom[ndom].rmax, zdom[ndom].rmin, ndom);
    exit (0);
  }

  rddoub ("Wind.t.init", &geo.twind_init);


  /* Next lines are to assure that we have the largest possible value of the 
   * sphere surrounding the system
   * JM 1710 -- if this is the first domain, then initialise geo.rmax see #305
   */
  if ((ndom == 0) || (zdom[ndom].rmax > geo.rmax))
  {
    geo.rmax = zdom[ndom].rmax;
  }
  geo.rmax_sq = geo.rmax * geo.rmax;


  /* Now get parameters that are specific to a given wind model

     Note: When one adds a new model, the only things that should be read in and modified
     are parameters in geo.  This is in order to preserve the ability to continue a calculation
     with the same basic wind geometry, without reading in all of the input parameters.  
   */

  if (zdom[ndom].wind_type == STAR)
  {
    get_stellar_wind_params (ndom);
  }
  else if (zdom[ndom].wind_type == SV)
  {
    get_sv_wind_params (ndom);
  }
  else if (zdom[ndom].wind_type == HYDRO)
  {
    get_hydro_wind_params (ndom);
  }
  else if (zdom[ndom].wind_type == CORONA)
  {
    get_corona_params (ndom);
  }
  else if (zdom[ndom].wind_type == KNIGGE)
  {
    get_knigge_wind_params (ndom);
  }
  else if (zdom[ndom].wind_type == HOMOLOGOUS)
  {
    get_homologous_params (ndom);
  }
  else if (zdom[ndom].wind_type == YSO)
  {
    get_yso_wind_params (ndom);
  }
  else if (zdom[ndom].wind_type == SHELL)       //NSH 18/2/11 This is a new wind type to produce a thin shell.
  {
    get_shell_wind_params (ndom);
  }
  else if (zdom[ndom].wind_type == IMPORT)      //Read in the wind model.
  {
    get_import_wind_params (ndom);
  }
  else
  {
    Error ("get_wind_parameters: Unknown wind type %d\n", zdom[ndom].wind_type);
    exit (0);
  }

  /* Get the filling factor of the wind */
  // XXX  This may  not in the right place to set the filling factor.  

  zdom[ndom].fill = 1.;

  /* JM 1606 -- the filling factor is now specified on a domain by domain basis. See #212
     XXX allows any domain to be allowed a filling factor but this should be modified when
     we know what we are doing with inputs for multiple domains. Could create confusion */

  rddoub ("Wind.filling_factor(1=smooth,<1=clumped)", &zdom[ndom].fill);

  return (0);
}




/**********************************************************/
/** 
 * @brief       Get the line transfer mode for the wind
 *
 * @param [in] None
 * @return  0 
 *
 * This rontinues simply gets the line tranfer mode for
 * all componensts of the wind.  After logging this
 * information the routine also reads in the atomic 
 * data.
 *
 *
 * ###Notes###
 * 1801 -   Refactored into this file in 1801.  It is
 *          not obvioous this is the best place for this
 *          routine since it refers to all components
 *          of the wind.

***********************************************************/



int
get_line_transfer_mode ()
{
  int user_line_mode = 0;
  rdint
    ("Line_transfer(0=pure.abs,1=pure.scat,2=sing.scat,3=escape.prob,4=anisotryopic,5=thermal_trapping,6=macro_atoms,7=macro_atoms+aniso.scattering)",
     &user_line_mode);

  /* JM 1406 -- geo.rt_mode and geo.macro_simple control different things. geo.rt_mode controls the radiative
     transfer and whether or not you are going to use the indivisible packet constraint, so you can have all simple 
     ions, all macro-atoms or a mix of the two. geo.macro_simple just means one can turn off the full macro atom 
     treatment and treat everything as 2-level simple ions inside the macro atom formalism */

  /* Set the default scattering and RT mode */
  geo.scatter_mode = SCATTER_MODE_ISOTROPIC;    // isotropic
  geo.rt_mode = RT_MODE_2LEVEL; // Not macro atom (SS)

  if (user_line_mode == 0)
  {
    Log ("Line_transfer mode:  Simple, pure absorption\n");
    geo.line_mode = user_line_mode;
  }
  else if (user_line_mode == 1)
  {
    Log ("Line_transfer mode:  Simple, pure scattering\n");
    geo.line_mode = user_line_mode;
  }
  else if (user_line_mode == 2)
  {
    Log ("Line_transfer mode:  Simple, single scattering\n");
    geo.line_mode = user_line_mode;
  }
  else if (user_line_mode == 3)
  {
    Log ("Line_transfer mode:  Simple, isotropic scattering, escape probabilities\n");
    geo.line_mode = user_line_mode;
  }
  else if (user_line_mode == 4)
  {
    Error ("get_line_transfer_mode: Line transfer mode %d is deprecated\n", user_line_mode);
    line_transfer_help_message ();
    exit (0);
  }
  else if (user_line_mode == 5)
  {
    Log ("Line_transfer mode:  Simple, thermal trapping, Single scattering \n");
    geo.scatter_mode = SCATTER_MODE_THERMAL;    // Thermal trapping model
    geo.line_mode = 3;
    geo.rt_mode = RT_MODE_2LEVEL;       // Not macro atom (SS) 
  }
  else if (user_line_mode == 6)
  {
    Log ("Line_transfer mode:  macro atoms, isotropic scattering  \n");
    geo.scatter_mode = SCATTER_MODE_ISOTROPIC;  // isotropic
    geo.line_mode = 3;
    geo.rt_mode = RT_MODE_MACRO;        // Identify macro atom treatment (SS)
    geo.macro_simple = 0;       // We don't want the all simple case (SS)
  }
  else if (user_line_mode == 7)
  {
    Log ("Line_transfer mode:  macro atoms, anisotropic  scattering  \n");
    geo.scatter_mode = SCATTER_MODE_THERMAL;    // thermal trapping
    geo.line_mode = 3;
    geo.rt_mode = RT_MODE_MACRO;        // Identify macro atom treatment (SS)
    geo.macro_simple = 0;       // We don't want the all simple case (SS)
  }
  else if (user_line_mode == 8)
  {
    Log ("Line_transfer mode: simple macro atoms, isotropic  scattering  \n");
    geo.scatter_mode = SCATTER_MODE_ISOTROPIC;  // isotropic
    geo.line_mode = 3;
    geo.rt_mode = RT_MODE_MACRO;        // Identify macro atom treatment i.e. indivisible packets
    geo.macro_simple = 1;       // This is for test runs with all simple ions (SS)
  }
  else if (user_line_mode == 9)
  {
    Log ("Line_transfer mode: simple macro atoms, anisotropic  scattering  \n");
    geo.scatter_mode = SCATTER_MODE_THERMAL;    // thermal trapping
    geo.line_mode = 3;
    geo.rt_mode = RT_MODE_MACRO;        // Identify macro atom treatment i.e. indivisible packets
    geo.macro_simple = 1;       // This is for test runs with all simple ions (SS)
  }
  else
  {
    Error ("Unknown line_transfer mode %d\n");
    line_transfer_help_message ();
    exit (0);
  }

  /* With the macro atom approach we won't want to generate photon 
     bundles in the wind so switch it off here. (SS) */
  if (geo.rt_mode == RT_MODE_MACRO)
  {
    Log ("python: Using Macro Atom method so switching off wind radiation.\n");
    geo.wind_radiation = 0;
  }


  /* read in the atomic data */
  rdstr ("Atomic_data", geo.atomic_filename);

  /* read a variable which controls whether to save a summary of atomic data
     this is defined in atomic.h, rather than the modes structure */

  if (modes.iadvanced)
  {

    rdint ("@Diag.write_atomicdata(0=no,anything_else=yes)", &write_atomicdata);
    if (write_atomicdata)
      Log ("You have opted to save a summary of the atomic data\n");
  }

  get_atomic_data (geo.atomic_filename);
  return (0);
}




/**********************************************************/
/** 
 * @brief      sets up the windcones for each domain
 *
 * @return     Always returns 0
 *
 * @details
 * 
 * This routine takes input variables, the minimim and maximum
 * radius of the wind at the disk, and the flow angles bounding
 * the wind, and calculates the variables used to define the cones
 * that bound each wind domain involving biconical flows.
 *
 * ### Notes ###
 *
 * The routine cycles through all of the existing domains, and
 * uses variables which have been read in or entered previously.
 *
 * The input variables that are used are typicall, wind_rho_min, wind_rho_max
 * and wind_thetamin and max.  They are defined in routines like,
 * get_sv_parameters.
 *
 * The angles thetamin and
 * thetamax are all defined from the z axis, so that an angle of 0
 * is a flow that is perpindicular to to the disk and one that is
 * close to 90 degrees will be parallel to the plane of the disk
 *
 * The routine files a structure of type cone, that contain, an element
 * z,  the place where the windcone intercepts the z axis, and an element
 * dzdr is the slope of the cone
 *
 * Wind cones are defined azimuthally around the z axis.
 *
 * For models that are not biconical flows, the windcones are set 
 * to include the entire domain.
 * 
 *
 **********************************************************/

int
setup_windcone ()
{
  int ndom;

  for (ndom = 0; ndom < geo.ndomain; ndom++)
  {

    if (zdom[ndom].wind_thetamin > 0.0)
    {
      zdom[ndom].windcone[0].dzdr = 1. / tan (zdom[ndom].wind_thetamin);
      zdom[ndom].windcone[0].z = (-zdom[ndom].wind_rho_min / tan (zdom[ndom].wind_thetamin));
    }
    else
    {
      zdom[ndom].windcone[0].dzdr = VERY_BIG;
      zdom[ndom].windcone[0].z = -VERY_BIG;;
    }


    if (zdom[ndom].wind_thetamax > 0.0)
    {
      zdom[ndom].windcone[1].dzdr = 1. / tan (zdom[ndom].wind_thetamax);
      zdom[ndom].windcone[1].z = (-zdom[ndom].wind_rho_max / tan (zdom[ndom].wind_thetamax));
    }
    else
    {
      zdom[ndom].windcone[1].dzdr = VERY_BIG;
      zdom[ndom].windcone[1].z = -VERY_BIG;;
    }
  }
  return (0);
}


/**********************************************************/
/** 
 * @brief      print out a help message about line transfer modes
 *
 * @return     Always returns 0
 **********************************************************/

int
line_transfer_help_message ()
{
  char *some_help;

  some_help = "\
\n\
Available line transfer modes and descriptions are: \n\
\n\
  0 Pure Absorption\n\
  1 Pure Scattering\n\
  2 Single Scattering\n\
  3 Escape Probabilities, isotropic scattering\n\
  5 Escape Probabilities, anisotropic scattering\n\
  6 Indivisible energy packets / macro-atoms, isotropic scattering\n\
  7 Indivisible energy packets / macro-atoms, anisotropic scattering\n\
  8 Indivisible energy packets, force all simple-atoms, anisotropic scattering\n\
  9 Indivisible energy packets, force all simple-atoms, anisotropic scattering\n\
\n\
  Standard mode is 5 for runs involving weight reduction and no macro-atoms\n\
  Standard macro-atom mode is 7\n\
\n\
See this web address for more information: https://github.com/agnwinds/python/wiki/Line-Transfer-and-Scattering\n\
\n\
\n\
";                              // End of string to provide one with help

  Log ("%s\n", some_help);

  return (0);
}
