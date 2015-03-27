


/**************************************************************************
                    Space Telescope Science Institute


  Synopsis:  

  These are routines coded specifically for including a scattering
  'torus' in a system.

  Description:	

  Arguments:  


  Returns:

  Notes:

 

  History:
11aug	ksl	Coding began     

 ************************************************************************/

#include <stdio.h>
#include <strings.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "atomic.h"
#include "python.h"


double
torus_rho (x)
     double x[];
{
  double z;
  double rho;

  /* !! ERROR - this is only approximate because of rho2nh */

  z = geo.compton_torus_rmax - geo.compton_torus_rmin;
  rho = geo.compton_torus_tau / (THOMPSON * rho2nh * z);
  return (rho);
}


/* These are the routines to determine whether a photon hits the torus
as defined as a region defined as a kind of pillbox centered ont he
disk.

This is not the best place to put some of these routines.  In particular
ds_to_cylinder should probably be part of phot_util, and ds_to_torus
may belong in photon2d (where ds_to_wind) is located

ds_to_cylinder is modeled on ds_to_sphere

110930	ksl	Began coding
*/

double
ds_to_cylinder (rho, p)
     double rho;
     struct photon *p;
{
  double a, b, c, root[2];
  int i;

  a = 1. - p->lmn[2] * p->lmn[2];
  b = 2. * (p->lmn[0] * p->x[0] + p->lmn[1] * p->x[1]);
  c = (p->x[0] * p->x[0] + p->x[1] * p->x[1]) - rho * rho;

  i = quadratic (a, b, c, root);

/* root[i] is the smallest positive root unless i is
negative in which case either both roots were negative or 
both roots were imaginary */

  if (i >= 0)
    return (root[i]);
  return (VERY_BIG);
}


double
ds_to_torus (pp)
     PhotPtr pp;
{
  struct photon ptest;
  double ds, ds_best, x;
  struct plane xplane;


  ds_best = VERY_BIG;

  /* Make sure we don't mess with pp */
  stuff_phot (pp, &ptest);
  ds = ds_to_cylinder (geo.compton_torus_rmin, &ptest);
  if (ds < VERY_BIG)
    {
      move_phot (&ptest, ds);
      if (fabs (ptest.x[2]) < geo.compton_torus_zheight)
	{
	  ds_best = ds;
	}
      stuff_phot (pp, &ptest);
    }

  ds = ds_to_cylinder (geo.compton_torus_rmax, &ptest);
  if (ds < ds_best)
    {
      move_phot (&ptest, ds);
      if (fabs (ptest.x[2]) < geo.compton_torus_zheight)
	{
	  ds_best = ds;
	}
      stuff_phot (pp, &ptest);
    }

  /* At this point we know whether the photon has interecepted
   * the wall of the cylinder, but we do not know if it intercepted
   * the top or bottom of the cylinder earlier
   */

  xplane.x[0] = 0;
  xplane.x[1] = 0;
  xplane.x[2] = geo.compton_torus_zheight;
  xplane.lmn[0] = 0.0;
  xplane.lmn[1] = 0.0;
  xplane.lmn[2] = 1.0;

  ds = ds_to_plane (&xplane, &ptest);
  // Note that ds to plane can return a negative number
  if (ds > 0 && ds < ds_best)
    {
      move_phot (&ptest, ds);
      x = fabs (ptest.x[0]);
      if (geo.compton_torus_rmin < x && x < geo.compton_torus_rmax)
	{
	  ds_best = ds;
	}
      stuff_phot (pp, &ptest);
    }

  xplane.x[2] = (-geo.compton_torus_zheight);

  ds = ds_to_plane (&xplane, &ptest);
  if (ds > 0 && ds < ds_best)
    {
      move_phot (&ptest, ds);
      x = fabs (ptest.x[0]);
      if (geo.compton_torus_rmin < x && x < geo.compton_torus_rmax)
	{
	  ds_best = ds;
	}
      stuff_phot (pp, &ptest);
    }

  return (ds_best);
}
