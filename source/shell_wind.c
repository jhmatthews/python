

/* 
   This file was created in Feb 2011.    The purpose is to have a model where we have a single shell of material. This is different from the stellar wind because we do not want the inner surface of the wind to touch the star. This requires tight control of the grid and makes for a very prescriptive model. 
We also need a special grid, which is also stored in this file.

 */


#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "atomic.h"
#include "python.h"

/***********************************************************
                                       Space Telescope Science Institute

 Synopsis:
	Routines to produce a thin shell of material at a given radius. 


Arguments:		

Returns:
 
Description:	
	

Notes:
	


History:
 	11feb	nsh	Coded as part of the effort to put power laws into python. It allows detailed testing.
	12jan	nsh	Shell wind rewritten to use existing C+L wind model.
 		
**************************************************************/


int
get_shell_wind_params ()
{
  double vtemp[3];
  double rhotemp[200];
  double rmin, dr, r[200];
  double postemp[3];
  double speedtemp;
  double cdensity;
  double shell_vmin, shell_vmax;	//Local variables to define shellwind
  double shell_rmin, shell_rmax;	//Local variables to define shellwind
  double factor;
  int i;


  Log ("Creating wind with a single shell for testing purposes\n");


/* In order to make life simple, we are first going to check that we are in spherical coordinates, if not change!! */

  if (geo.coord_type != SPHERICAL)
    {
      Error
	("For the shell type wind, we should be in spherical coordinates - changing....\n");
      geo.coord_type = 0;
    }

/* This may not be important, but we sould make sure that NDIM is 4... */

  if (geo.ndim != 4)
    {
      Error
	("For the shell type wind, we take control of the grid, and need NDIM to be the minimum - 4 - changing\n");
      geo.ndim = 4;
      geo.mdim = 1;
    }


  geo.stellar_wind_mdot = 1.e-6;
  geo.wind_rmin = geo.rstar;
  geo.cl_beta = 1.0;
  shell_rmin = geo.cl_rmin = 2.8e9;
  shell_vmin = geo.cl_v_zero = 200e5;
  shell_vmax = geo.cl_v_infinity = 3000e5;

  geo.stellar_wind_mdot /= MSOL / YR;
  rddoub ("shell_wind_mdot(msol/yr)", &geo.stellar_wind_mdot);
  geo.stellar_wind_mdot *= MSOL / YR;

  rddoub ("shell.wind.radmin(cm)", &geo.wind_rmin);	/*Radius where wind begins */
  if (geo.wind_rmin < geo.rstar)
    {
      Error
	("get_shell_wind_params: It is unreasonable to have the wind start inside the star!\n");
      Log ("Setting geo.wind_rmin to geo.rstar\n");
      geo.wind_rmin = geo.rstar;
    }
  geo.cl_rmin = shell_rmin = geo.wind_rmin;


/*120130 NSH the next two lines have been modified to mean that the wind will end up as a CL wind, but the v_0 and v_infinity will be calulated here from these two variables, which are now local */

  rddoub ("shell.wind_v_at_rmin(cm)", &shell_vmin);	/* Velocity at base of the wind */
  rddoub ("shell.wind.v_at_rmax(cm)", &shell_vmax);	/* Final speed of wind in units of escape velocity */

/*120130 NSH This line stays as it is */
  rddoub ("shell.wind.acceleration_exponent", &geo.cl_beta);	/* Accleration scale exponent for a CL wind */
  Log ("Geo rmax = %f\n", geo.rmax);
  shell_rmax = geo.wind_rmax = geo.rmax;
/*120130 NSH These next lines invert the cl velocity equation to get the cl factors from the local shell factors */
  geo.cl_v_zero = shell_vmin;
  factor = pow ((1 - (geo.wind_rmin / geo.wind_rmax)), geo.cl_beta);
  geo.cl_v_infinity =
    (shell_vmax - shell_vmin + shell_vmin * factor) / factor;






/* Assign the generic parameters for the wind the generic parameters of the wind */
//  geo.wind_rmin = geo.rstar;

  geo.wind_thetamin = 0.0;
  geo.wind_thetamax = 90. / RADIAN;

/* define the the variables that determine the gridding */
  geo.wind_rho_min = 0;
  geo.wind_rho_max = geo.rmax;

  /* if modes.adjust_grid is 1 then we have already adjusted the grid manually */
  if (modes.adjust_grid == 0)
    {
      geo.xlog_scale = 0.3 * geo.rstar;
      geo.zlog_scale = 0.3 * geo.rstar;
    }

/* Since this is a diagnostic routine, we will write out some information to check it is doing what we think) */
  Log ("shell rmin=%f shell rmax=%f\n", shell_rmin, shell_rmax);
  dr = (shell_rmax - shell_rmin) / 100.0000;
  Log ("dr= %e, root2= %10.30e\n", dr, pow (2.0, 0.5));
  rmin = shell_rmin - (dr);

  for (i = 0; i < 103; i++)
    {
      r[i] = rmin + i * dr;

      postemp[0] = postemp[2] = r[i] / pow (2.0, 0.5);
      postemp[1] = 0.0;
      speedtemp = stellar_velocity (postemp, vtemp);
      rhotemp[i] = stellar_rho (postemp) * rho2nh;
//              dvtemp[i]=shell_dv(postemp);
      Log ("ring=%i,x=%e,r=%10.30e,speed=%10.20e,density=%10.20e\n", i,
	      r[i] / pow (2.0, 0.5), r[i], speedtemp, rhotemp[i]);
    }

  cdensity = 0.0;
  for (i = 1; i < 100; i++)
    {
      cdensity += ((rhotemp[i] + rhotemp[i + 1]) / 2.) * dr;
    }
  Log ("Column density of hydrogen=%e\n", cdensity);


  return (0);
}



/*NSH Jun2012 - 72c removed double shell_velocity(x,v) and double shell_rho(x) since these functions are no longer needed */
