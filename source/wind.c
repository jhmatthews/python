


/* The routines in this file define and summarize the properties of the wind.  It is only useful in the
   2-d version of the code.  There should now be no SV dependent code in this file, all of that having
   been moved to sv.c.  Although information about the wind is conveyed through the structure geo,
   none of the sv or cl specific values should be needed here. 

   The routines here have NO explicit dependencies on a cylindrical coordinate system
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "atomic.h"
#include "python.h"


/***********************************************************
                  Space Telescope Science Institute

 Synopsis:
	where_in_wind determines whether a position x is in the wind region.
 
Arguments:		
	double x[3]	a position
Returns:
 	where_in_wind returns 
	 W_ALL_INTORUS if the photon is in the torus
	 W_ALL_INWIND if the photon is in the wind, 
 	-1 if it is inside the inner wind cone, and 
 	-2 if it is outside the outer wind cone.
 	-3 if it is inside the minimum radius for the wind
	-4 if it is outside the maximum radius for the wind
	-5 if it is inside the disk
Description:	
	
		
Notes:
	It does not tell you whether it is in the grid or not, just whether the photon is
	between the two wind cones and/or inside minimum or maximum radius of the wind.     
	It does not update pp->grid

	11Aug. - ksl - it does not look like any of the negative return values are actually utilized
	anywhere.  It could be that one should return W_NOT_INWIND in all
	of these cases
History:
 	97jan   ksl	Coding on python began.
 	98dec	ksl	Modified so the call to where_in_wind involves only
			a position.  Old call was where_in_wind(w,p)
 	98dec	ksl	Added another option to declare the the a position is
 			inside the minimum radius of the wind. This makes it possible
 			to have a wind from a star. 
	04aug	ksl	52 -- Added possibility that a position
			is in the disk, which is now possible
			if the disk has vertical extent
        04Aug   SS      Minor modification - "else" removed.
	11Aug	ksl	Allowed for torus, and adapted so that it uses
			some of the same defined variables as w->inwind
	11Nov	ksl	Made changes to attempt to fix errors in the
			implementation of the Elvis model
	13May	nsh	Attempt to cope with an issue where code runs
			very slowly on iridis if thetamax is pi/2,
			now if it is within machine precision of
			pi/2, is returns that the photon is in wind
			and does not do the check. 
			
**************************************************************/

int
where_in_wind (x)
     double x[];
{
  double rho, rho_min, rho_max, z;
  int ireturn;

  ireturn = W_ALL_INWIND;

  //110814 - Currently geo.wind_rmin is always set to geo.rstar except for models that are 
  //calculated by others.  Thus we would all the grids to start at rmin.


  /* First check to see if photon is inside star or outside wind */
  if ((z = length (x)) < geo.wind_rmin)
    {
      ireturn = (-3);		/*x is inside the wind  radially */
    }
  else if (z > geo.wind_rmax)
    {
      ireturn = (-4);		/*the position is beyond the wind radially */
    }


  z = fabs (x[2]);		/* This is necessary to get correct answer above
				   and below plane */
  rho = sqrt (x[0] * x[0] + x[1] * x[1]);	/* This is distance from z axis */

  /* Now check to see if position is inside the disk */
  if (geo.disk_type == 2)
    {
      if (rho < geo.diskrad && z < zdisk (rho))
	return (-5);
    }

  /* 70b -ksl - Now check if the position is in the torus.  This check
     precedes the check to see if the position is in the wind
   */

  if (geo.compton_torus)
    {
      if (geo.compton_torus_rmin < rho && rho < geo.compton_torus_rmax
	  && z < geo.compton_torus_zheight)
	{
	  return (W_ALL_INTORUS);
	}
    }

  /* Now for the elvis wind model, check to see if the position is
   * inside all of the wind or if it is in the vertical section of the wind
   *
   * Note the use of sv_rmax and sv_rmin here, because thise are not the
   * same as for the windcones gescribed by geo.wind_rmin and geo.wind_rmax
   * 111124 ksl
   */

  if (geo.wind_type == 8)
    {
      if (rho < geo.sv_rmin)
	{
	  return (-1);
	}
      if (rho < geo.sv_rmax && z < geo.elvis_offset)
	{
	  return (W_ALL_INWIND);
	}
    }



  /* Check if one is inside the inner windcone */
  if (rho < (rho_min = geo.wind_rho_min + z * tan (geo.wind_thetamin)))
    {
      ireturn = (-1);
    }

  /* Finally check if positon is outside the outer windcone */
  /* NSH 130401 - The check below was taking a long time if geo.wind_thetamax was very close to pi/2.
     check inserted to simply return INWIND if geo.wind_thetamax is within machine precision of pi/2. */

  else if (fabs (geo.wind_thetamax - PI / 2.0) > 1e-6)	/* Only perform the next check if thetamax is not equal to pi/2 */
    {
      if (rho > (rho_max = geo.wind_rho_max + z * tan (geo.wind_thetamax)))
	{
	  ireturn = (-2);
	}
    }				/* If thetamax is equal to pi/2, then the photon must be in the wind if it has got to this point */


  return (ireturn);
}



/***********************************************************
                                       Space Telescope Science Institute

 Synopsis:
	 wind_check checks the wind structure for reasonability 

Arguments:		
	WindPtr w;			The structure which defines the wind in Python
 	int n				n >= 0  then an element of the array will be checked
					n<0     then the entire wind structure will be checked 

Returns:
 
Description:

The routine should be called in a number of ways

wind_check(w,-1);  to check the entire structure
wind_check(w,50);   to check element 50 of the structure
wind_check(w[50],0) to check elemetn 50 of the structure

Notes:

History:
	99dec29 ksl	Added routines for proga's wind 
	04dec	ksl	Added some additional checks, and forced
			quit when a sane_check error found
**************************************************************/

int
wind_check (www, n)
     WindPtr www;
     int n;
{
  Log ("Got to wind_check\n");
  int i, j, k, istart, istop;
  if (n < 0)
    {
      istart = 0;
      istop = NDIM2;
    }
  else
    {
      istart = n;
      istop = istart + 1;
    }

  for (i = istart; i < istop; i++)
    {
      for (j = 0; j < 3; j++)
	{
	  if (sane_check (www[i].x[j]))
	    {
	      Error ("wind_check:sane_check www[%d].x[%d] %e\n", i, j,
		     www[i].x[j]);
	    }
	  if (sane_check (www[i].xcen[j]))
	    {
	      Error ("wind_check:sane_check www[%d].xcen[%d] %e\n", i, j,
		     www[i].xcen[j]);
	    }
	  if (sane_check (www[i].v[j]))
	    {
	      Error ("wind_check:sane_check www[%d].v[%d] %e\n", i, j,
		     www[i].v[j]);
	    }
	}
      for (j = 0; j < 3; j++)
	{
	  for (k = 0; k < 3; k++)
	    {
	      if (sane_check (www[i].v_grad[j][k]))
		{
		  Error ("wind_check:sane_check www[%d].v_grad[%d][%d] %e\n",
			 i, j, k, www[i].v_grad[j][k]);
		}
	    }

	}
    }

  Log
    ("Wind_check: Punchthrough distance DFUDGE %e www[1].x[2] %e\n",
     DFUDGE, www[1].x[2]);
  Log ("Finished wind check\n");
  return (0);
}


/* model_velocity(x,v)
 * Calculate the wind velocity at a specific point in space from the original 
 * usually analytic expressions
 *
 * 04aug	ksl	52 -- adapted from wind2d.c as part of effort to 
 * 			handle multiple coordinate systems
 */

double
model_velocity (x, v)
     double x[], v[];
{
  double speed;


  if (geo.wind_type == 0)
    {
      speed = sv_velocity (x, v);
    }
  else if (geo.wind_type == 1)
    {
      speed = stellar_velocity (x, v);
    }
  else if (geo.wind_type == 3)
    {
      speed = proga_velocity (x, v);
    }
  else if (geo.wind_type == 4)
    {
      speed = corona_velocity (x, v);
    }
  else if (geo.wind_type == 5)
    {
      speed = kn_velocity (x, v);
    }
  else if (geo.wind_type == 6)
    {
      speed = homologous_velocity (x, v);
    }
  else if (geo.wind_type == 7)
    {
      speed = yso_velocity (x, v);
    }
  else if (geo.wind_type == 8)
    {
      speed = elvis_velocity (x, v);
    }
  else if (geo.wind_type == 9)
    {
      speed = stellar_velocity (x, v);
    }
  else
    {
      Error ("wind: Unknown windtype %d\n", geo.wind_type);
      exit (0);
    }

  return (speed);
}



/* model_vgrad calculates the velocity gradient at postions in the flow based on
 * the analytic wind models.  This routine is normally used only during the initialization
 * of the wind
 *
 * Notes:  This routine replaces a series of routines which had been written for each type
 * of analytic model
 *
 * History
 * 04aug	ksl	52 adapted form define_wind in wind2d.c
 */

int
model_vgrad (x, v_grad)
     double x[], v_grad[][3];
{

  double v0[3], v1[3];
  double dx[3], dv[3];
  double ds;
  int i, j;
  int vsub (), stuff_v ();

  model_velocity (x, v0);


  ds = 0.001 * length (x);
  if (ds < 1.e7)
    ds = 1.e7;

  for (i = 0; i < 3; i++)
    {
      stuff_v (x, dx);
      dx[i] += ds;

      model_velocity (dx, v1);

      if (sane_check (v1[0]) || sane_check (v1[1]) || sane_check (v1[2]))
	{
	  Error ("model_vgrad:sane_check dx %f %f %f v0 %f %f %f\n", dx[0],
		 dx[1], dx[2], v1[0], v1[1], v1[2]);
	}

      vsub (v1, v0, dv);
      for (j = 0; j < 3; j++)
	dv[j] /= ds;
      stuff_v (dv, &v_grad[i][0]);
    }

  return (0);



}

/* model_rho calculates the density of the wind in from the flow based on the
analytic wind models.  This routine is normally only used during the initialization
of the wind
History
04aug	ksl	52--adapted from define_wind in wind2d.c
11aug	ksl	70b- added call to calculate density in the torus
 */

double
model_rho (x)
     double x[];
{
  double rho;
  if (geo.wind_type == 0)
    {
      rho = sv_rho (x);
    }
  else if (geo.wind_type == 1)
    {
      rho = stellar_rho (x);
    }
  else if (geo.wind_type == 3)
    {
      rho = proga_rho (x);
    }
  else if (geo.wind_type == 4)
    {
      rho = corona_rho (x);
    }
  else if (geo.wind_type == 5)
    {
      rho = kn_rho (x);
    }
  else if (geo.wind_type == 6)
    {
      rho = homologous_rho (x);
    }
  else if (geo.wind_type == 7)
    {
      rho = yso_rho (x);
    }
  else if (geo.wind_type == 8)
    {
      rho = elvis_rho (x);
    }
  else if (geo.wind_type == 9)
    {
      rho = stellar_rho (x);
    }
  else
    {
      Error ("wind2d: Unknown windtype %d\n", geo.wind_type);
      exit (0);
    }

  /* 70b - as this is written the torus simply overlays the wind */

  if ((where_in_wind (x) == W_ALL_INTORUS)
      || (where_in_wind (x) == W_PART_INTORUS))
    {
      rho = torus_rho (x);
    }

  return (rho);

}
