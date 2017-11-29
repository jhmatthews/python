


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
	where_in_wind determines whether a position x is in a wind region.
 
Arguments:		
	double x[3]	a position
Returns:
 	where_in_wind returns 

	W_ALL_INWIND if the position is in a wind,
	and the domain number

     	


	These are the nominal returns if there is only one domain.
	If there are multiple domains these returns are not in 
	general meaning full and should be rewritten

	 W_ALL_INWIND if the photon is in the wind, 
 	-1 if it is inside the inner wind cone, and 
 	-2 if it is outside the outer wind cone.
 	-3 if it is inside the minimum radius for the wind
	-4 if it is outside the maximum radius for the wind
	-5 if it is inside the disk
Description:	
	
		
Notes:
	Where_in_wind does not tell you whether a postion is in the grid or not, just whether 
	the photon is between the two wind cones and/or inside minimum or maximum radius of the wind.     

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
	15aug	jm/ksl	Modifications to determine what domain
			a positions is in if any.  With multiple 
			domains, it is unclear that what inside
			and outside the wind mean.  
	15aug	ksl	Changed where in wind so it returns whether
			you are in the wind or not, and added a variable
			so that it also returns the domain if you are in 
			the wind.
			
**************************************************************/

int
where_in_wind (x, ndomain)
     double x[];
     int *ndomain;
{
  double rho, rad, rho_min, rho_max, z;
  int ireturn;
  int ndom, n;
  DomainPtr one_dom;


  /* First check for items, like the diks that are not domain
   * related, like the disk */

  z = fabs (x[2]);              /* Necessary to get correct answer above
                                   and below plane */
  rho = sqrt (x[0] * x[0] + x[1] * x[1]);       /* This is distance from z axis */


  /* Check if position is inside the disk for a vertically extended disk */
  if (geo.disk_type == DISK_VERTICALLY_EXTENDED)
  {
    if (rho < geo.diskrad && z < zdisk (rho))
    {
      *ndomain = -1;
      return (W_IN_DISK);
    }
  }


  /* Now check whether position is a wind region of any of the domains.
   * This is done in reverse order on the assumption that our domains
   * are layered on top of one another.  */

  ireturn = W_NOT_INWIND;
  *ndomain = -1;

  rad = length (x);

  for (ndom = geo.ndomain - 1; ndom > -1; ndom--)
  {

    one_dom = &zdom[ndom];

    /* JM 1711 -- if we have an imported model, then we don't have angles defined
       and the "inwind" variable must instead be calculated by querying the
       grid covered by this domain with where_in_grid */

    if (one_dom->wind_type == IMPORT) 
    {
        one_dom = &zdom[ndom];

        n = where_in_grid (ndom, x);
        if (n >= 0) {
          ireturn = wmain[n].inwind;
          *ndomain = ndom;
          break;
        }
    }

    else 
    {
      /* First check to see if photon is inside or outside wind */

      if (rad < one_dom->rmin)
      {
        continue;                 /*x is inside the wind  radially */
      }
      if (rad > one_dom->rmax)
      {
        continue;                 /*the position is beyond the wind radially */
      }

      if (z > one_dom->zmax)
      {
        continue;                 /*the position is beyond the wind radially */
      }




      /* Check if one is inside the inner windcone */
      if (rho < (rho_min = one_dom->wind_rho_min + z * tan (one_dom->wind_thetamin)))
      {
        continue;
      }

      /* Finally check if positon is outside the outer windcone */
      /* NSH 130401 - The check below was taking a long time if geo.wind_thetamax was very close to pi/2.
         check inserted to simply return INWIND if geo.wind_thetamax is within machine precision of pi/2. */

      if (fabs (one_dom->wind_thetamax - PI / 2.0) > 1e-6)        /* Only perform the next check if thetamax is not equal to pi/2 */
      {
        if (rho > (rho_max = one_dom->wind_rho_max + z * tan (one_dom->wind_thetamax)))
        {
          continue;
        }
      }

      /* At this point we have passed all of the tests for being in the wind */

      *ndomain = ndom;
      ireturn = W_ALL_INWIND;
      break;

    }
  }
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
        Error ("wind_check:sane_check www[%d].x[%d] %e\n", i, j, www[i].x[j]);
      }
      if (sane_check (www[i].xcen[j]))
      {
        Error ("wind_check:sane_check www[%d].xcen[%d] %e\n", i, j, www[i].xcen[j]);
      }
      if (sane_check (www[i].v[j]))
      {
        Error ("wind_check:sane_check www[%d].v[%d] %e\n", i, j, www[i].v[j]);
      }
    }
    for (j = 0; j < 3; j++)
    {
      for (k = 0; k < 3; k++)
      {
        if (sane_check (www[i].v_grad[j][k]))
        {
          Error ("wind_check:sane_check www[%d].v_grad[%d][%d] %e\n", i, j, k, www[i].v_grad[j][k]);
        }
      }

    }
  }

  Log ("Wind_check: Punchthrough distance DFUDGE %e www[1].x[2] %e\n", DFUDGE, www[1].x[2]);
  Log ("Finished wind check\n");
  return (0);
}


/* model_velocity(ndom, x,v)

Calculate the wind velocity at a specific point in space from the original 
usually analytic expressions

 History

	04aug	ksl	52 -- adapted from wind2d.c as part of effort to 
 			handle multiple coordinate systems
	15Aug	ksl	Updated for domains
 */


double
model_velocity (ndom, x, v)
     double x[], v[];
     int ndom;
{
  double speed;

  if (zdom[ndom].wind_type == SV)
  {
    speed = sv_velocity (x, v, ndom);
  }
  else if (zdom[ndom].wind_type == STAR)
  {
    speed = stellar_velocity (ndom, x, v);
  }
  else if (zdom[ndom].wind_type == HYDRO)
  {
    speed = hydro_velocity (ndom,x, v);
  }
  else if (zdom[ndom].wind_type == CORONA)
  {
    speed = corona_velocity (ndom, x, v);
  }
  else if (zdom[ndom].wind_type == KNIGGE)
  {
    speed = kn_velocity (ndom, x, v);
  }
  else if (zdom[ndom].wind_type == HOMOLOGOUS)
  {
    speed = homologous_velocity (ndom, x, v);
  }
  else if (zdom[ndom].wind_type == YSO)
  {
    speed = yso_velocity (ndom, x, v);
  }
  else if (zdom[ndom].wind_type == SHELL)
  {
    speed = stellar_velocity (ndom, x, v);
  }
  else if  (zdom[ndom].wind_type == IMPORT)
  {
      speed=import_velocity(ndom,x,v);
  }
  else
  {
    Error ("wind: Unknown windtype %d for doman %d\n", zdom[ndom].wind_type,ndom);
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
model_vgrad (ndom, x, v_grad)
     double x[], v_grad[][3];
     int ndom;
{

  double v0[3], v1[3];
  double dx[3], dv[3];
  double ds;
  int i, j;
  int vsub (), stuff_v ();

  model_velocity (ndom, x, v0);


  ds = 0.001 * length (x);
  if (ds < 1.e7)
    ds = 1.e7;

  for (i = 0; i < 3; i++)
  {
    stuff_v (x, dx);
    dx[i] += ds;

    model_velocity (ndom, dx, v1);

    if (sane_check (v1[0]) || sane_check (v1[1]) || sane_check (v1[2]))
    {
      Error ("model_vgrad:sane_check dx %f %f %f v0 %f %f %f\n", dx[0], dx[1], dx[2], v1[0], v1[1], v1[2]);
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
model_rho (ndom, x)
     int ndom;
     double x[];
{
  double rho;

  if (zdom[ndom].wind_type == SV)
  {
    rho = sv_rho (ndom, x);
  }
  else if (zdom[ndom].wind_type == STAR)
  {
    rho = stellar_rho (ndom, x);
  }
  else if (zdom[ndom].wind_type == HYDRO)
  {
    rho = hydro_rho (x);
  }
  else if (zdom[ndom].wind_type == CORONA)
  {
    rho = corona_rho (ndom, x);
  }
  else if (zdom[ndom].wind_type == KNIGGE)
  {
    rho = kn_rho (ndom, x);
  }
  else if (zdom[ndom].wind_type == HOMOLOGOUS)
  {
    rho = homologous_rho (ndom, x);
  }
  else if (zdom[ndom].wind_type == YSO)
  {
    rho = yso_rho (ndom, x);
  }
  else if (zdom[ndom].wind_type == SHELL)
  {
    rho = stellar_rho (ndom, x);
  }
  else if (zdom[ndom].wind_type == IMPORT)
  {
    rho = import_rho(ndom,x);
  }
  else
  {
    Error ("wind2d: Unknown windtype %d for domain %d\n", zdom[ndom].wind_type, ndom);
    exit (0);
  }

  return (rho);

}
