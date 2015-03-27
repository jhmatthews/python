
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "atomic.h"
#include "python.h"


/* The reference temperature for the disk in degrees.  The maximum temperature
   actuallly seen in the disk is 0.488 tdisk
 */
double
tdisk (m, mdot, r)
     double m, mdot, r;
{
  double t;
  t = 3. * G / (8. * PI * STEFAN_BOLTZMANN) * m * mdot / (r * r * r);
  t = pow (t, 0.25);
  return (t);
}


/* effective temperature of accretion disk as a function of r
   inputs:      
   	t    reference temperature of the disk in degrees.
   	x    distance from center of disk in units of r/rmin
   where rmin is the inner edge of the disk.

   Notes: Originally (up to about python_50 only a standard steady state disk
   	was supported, for whihc a  reference is  Wade, 1984 MNRAS 208, 381.              

   History: 
   04June	SS	Modified to include a correction factor to account for
   			heating of the disk by the star following Proga et al. 1999, MNRAS, 310, 476. 
   			for now only used in the YSO model.
   04Aug	ksl	Modified so that geo.disk_illum controls how modifications
   			to the temperature are treated.  This is in preparation to allowing
			for a variety of disk illumination possibilities
   06Jan	ksl	Fixed problem with absorbing disk
*/

double
teff (t, x)
     double t, x;
{
  double q, theat, r;
  double pow ();
  double disk_heating_factor;
  int kkk;
  int n;


  q = 0.0;			/* NSH 130605 to remove o3 compile error */


  if (x < 1)
    {
      Error ("teff: x %f less than 1.0\n", x);
      return (0.0);
    }


  if ((geo.disk_tprofile != 0)
      && ((x * geo.rstar) < blmod.r[blmod.n_blpts - 1]))
    {
      if ((r = (x * geo.rstar)) < blmod.r[0])
	{
	  return (blmod.t[0]);
	}
      else
	{
	  for (n = 1; n < blmod.n_blpts; n++)
	    {
	      if ((r < blmod.r[n]) && (r > blmod.r[n - 1]))
		{
		  return (blmod.t[n]);
		}
	    }
	  Error
	    ("tdisk: inside BL profile region but failed to identify temp.\n");
	}
    }
  else
    {

      q = (1.e0 - pow (x, -0.5e0)) / (x * x * x);
      q = t * pow (q, 0.25e0);

      if (geo.disk_illum == 2)	// Absorb photons and increase t so that heat is radiated
	{
	  r = x * geo.rstar;	// 04aug -- Requires fix if disk does not extend to rstar
	  kkk = 1;		// photon cannot hit the disk at r<qdisk.r[0]
	  while (r > qdisk.r[kkk] && kkk < NRINGS - 1)
	    kkk++;
	  /* Note that disk has 2 sides */
	  theat =
	    qdisk.heat[kkk -
		       1] / (2. * PI * (qdisk.r[kkk] * qdisk.r[kkk] -
					qdisk.r[kkk - 1] * qdisk.r[kkk - 1]));

	  /* T_eff is given by T_eff**4= T_disk**4+Heating/area/STEFAN_BOLTZMANN */
	  q = pow (q * q * q * q + (theat / STEFAN_BOLTZMANN), 0.25);

	}
      else if (geo.disk_illum == 3)	// Analytic approximation for disk heating by star; implemented for YSOs
	{
	  disk_heating_factor = pow (geo.tstar / t, 4.0);
	  disk_heating_factor *=
	    (asin (1. / x) - (pow ((1. - (1. / (x * x))), 0.5) / x));
	  disk_heating_factor /= PI;

	  disk_heating_factor *= x * x * x;
	  disk_heating_factor /= (1 - sqrt (1. / x));
	  disk_heating_factor += 1;

	  q *= pow (disk_heating_factor, (1. / 4.));

	}
    }
  return (q);
}

double
gdisk (mass, mdot, rmin)
     double mass, rmin, mdot;
{
  double g0;
  g0 =
    0.625 * log10 (mass / MSOL) - 1.875 * log10 (rmin / 1.e9) +
    0.125 * log10 (mdot / 1.e16);
  g0 = 5.96e5 * pow (10., g0);
  return (g0);
}

double
geff (g0, x)
/* effective gravity of standard accretion disk as a function of r
	inputs:         g 	reference gravity in cm s**-2
			x	distance from center in units of r/rmin

*/
     double g0, x;
{
  double q;
  q = (1.0e0 - pow (x, -0.5e0));
  q = pow (x, -1.875e0) * pow (q, 0.125);
  q = g0 * q;
  return (q);
}


/* vdisk calculate the speed and velocity v at which material given
   specific position in the disk.
  
Description:

   The routine determines the speed by interpolating on stored
   values in the structure disk
  
   It then takes the xproduct of a unit vector pointing north and
   the position, to determine the direction of motion at this point
   in the disk
  
   Finally it rescales this to get the actual velocity.
  
   is travelling in the
   disk by interpolating on the velocities that are contained 

Notes"

	The routine projects the input variable x on to the xy plane
	before it calculates velocities.  

History:

	08mar	ksl	Added explanation of what the routine does. The
			routine was not changed.
        13sep	nsh	Changed linterp call to reflect a new log option.
 */
double north[] = { 0.0, 0.0, 1.0 };

double
vdisk (x, v)
     double x[];
     double v[];
{
  double xhold[3];
  double r, speed;
  int linterp ();
  stuff_v (x, xhold);
  xhold[2] = 0.0;
  r = length (xhold);
  linterp (r, disk.r, disk.v, NRINGS, &speed,0); //interpolate in linear space
  cross (north, xhold, v);	/* The velocity vector direction is given by north x r */
  renorm (v, speed);
  return (speed);
}

/***********************************************************
             Space Telescope Science Institute
 
 Synopsis:  zdisk
   
 Arguments:
        r       a radial position in the disk.
 
 
 Returns:
        The veritical height of the disk at that point in the xy
	plane. The number should be positive or zero
   
Description:
 
Notes:
 
 
History:
        04aug   ksl     52--Coded as part of effort to put 
			disks with vertical extent into python
        04Aug   SS      Multiply by "geo.diskrad" added.
                                                                                         
**************************************************************/

double
zdisk (r)
     double r;
{
  double z;
  z = geo.disk_z0 * pow (r / geo.diskrad, geo.disk_z1) * geo.diskrad;
  return (z);
}


/***********************************************************
             Space Telescope Science Institute
 
 Synopsis: ds_to_disk calculates the distance that a photon
 	would need to travel from its current position to hit 
	disk. Negative distances are allowed.
   
 Arguments:
        p       a photon pointer.
 
 
 Returns:
        The distance to the disk.  The photon is not going to
	hit the disk travelling in either direction
	returns VERY_BIG.

   
Description:
 
Notes:
	There are other routines that will return a negative distance, 
	but this is not one of them, since at present we have not 
	identified a reason to separate 
 
 
History:
        04aug   ksl     52--Coded as part of effort to put 
			disks with vertical extent into python
        04Aug   SS      Several minor modifications made.
                        ds_to_disk now also takes a second
                        input "miss_return" which tells it what
                        to return for trajectories that miss the
                        disk. 0 = return infinity while anything
                        else gives a return of the distance to 
                        the horizontal plane that touched the
                        edge of the disk. Only important for
                        non-flat disks: it influences what happens when
                        the disk is initialised.
                                                                                       
**************************************************************/

int ds_to_disk_init = 0;
struct photon ds_to_disk_photon;
struct plane diskplane, disktop, diskbottom;

double
ds_to_disk (p, miss_return)
     struct photon *p;
     int miss_return;
{
  double x1, x2;
  double s_plane, s_top, s_bottom, s_sphere, s_disk;
  double r_plane, r_top, r_bottom;
  double smin, smax;
  struct photon phit;
  void disk_deriv ();


  if (geo.disk_type == 0)
    return (VERY_BIG);		/* There is no disk! */

  if (ds_to_disk_init == 0)
    {				/* Initialize 3 structures that define
				   the plane of the disk, and two other
				   planes that encompass the disk */

      diskplane.x[0] = diskplane.x[1] = diskplane.x[2] = 0.0;
      diskplane.lmn[0] = diskplane.lmn[1] = 0.0;
      diskplane.lmn[2] = 1.0;	//changed by SS August 04


      disktop.x[0] = disktop.x[1] = 0.0;
      disktop.x[2] = geo.diskrad * geo.disk_z0;
      disktop.lmn[0] = disktop.lmn[1] = 0.0;
      disktop.lmn[2] = 1.0;	//changed by SS August 04

      diskbottom.x[0] = diskbottom.x[1] = 0.0;
      diskbottom.x[2] = (-geo.diskrad * geo.disk_z0);
      diskbottom.lmn[0] = diskbottom.lmn[1] = 0.0;
      diskbottom.lmn[2] = 1.0;	//changed by SS August 04

      ds_to_disk_init++;	// Only initialize once

    }

  /* Now calculate the place where the photon hits the diskplane */

  s_plane = ds_to_plane (&diskplane, p);
  stuff_phot (p, &phit);
  move_phot (&phit, s_plane);
  r_plane = sqrt (phit.x[0] * phit.x[0] + phit.x[1] * phit.x[1]);


  if (geo.disk_type == 1)
    {
      if (r_plane > geo.diskrad)
	return (VERY_BIG);
      return (s_plane);
    }

  /* OK now we have to deal with the hard case.  We would like to
   * avoid actually having to calculate the intercept to the disk
   * if we can because this is likely time consuming. So we first
   * determine this, by checking where the ray hits a sphere that 
   * just encompasses the disk 
   */

  s_top = ds_to_plane (&disktop, p);
  stuff_phot (p, &phit);
  move_phot (&phit, s_top);
  r_top = sqrt (phit.x[0] * phit.x[0] + phit.x[1] * phit.x[1]);

  s_bottom = ds_to_plane (&diskbottom, p);
  stuff_phot (p, &phit);
  move_phot (&phit, s_bottom);
  r_bottom = sqrt (phit.x[0] * phit.x[0] + phit.x[1] * phit.x[1]);

  /* Now if rtop and r_bottom are both greater than diskrad
   * then this photon missed the disk */
  if (r_top > geo.diskrad && r_bottom > geo.diskrad)
    {
      if (miss_return == 0)
	{
	  return (VERY_BIG);
	}
      else
	{
	  if (s_top > s_bottom)
	    {
	      return (s_top);
	    }
	  else
	    {
	      return (s_bottom);
	    }
	}
    }

  /* OK, at this point we know the photon is on a path that passes 
   * through (or passed through the disk) and we must locate it. 
   * There are two possibilities.  It hit the disk edge, or it hit
   * the face of the disk.  Most of the time it will hit the disk face
   * so we will calculate this first, and then check if s is less
   * than s_sphere. 
   *
   * To setup rtsafe, we have to find distances which bracket what
   * we want.  smax should be no more than the distance to the
   * disk_plane (assuming we have a photon headed toward the plane)
   *
   * Note that the next little section is inted to select between the
   * top and bottom face of the disk
   */

  /* I've modified the next if statement to try and get it to choose
     the correct x2 (have to be careful about the r's being less
     than diskrad. SS August 04. */

  x1 = s_plane;
  if (p->x[2] > 0.0)
    {
      if (r_top < geo.diskrad)
	{
	  x2 = s_top;
	}
      else
	{
	  x2 = s_bottom;
	}
    }
  else
    {
      if (r_bottom < geo.diskrad)
	{
	  x2 = s_bottom;
	}
      else
	{
	  x2 = s_top;
	}
    }

  if (fabs (x2) > fabs (x1))	//fabs added by SS August 04
    {
      smin = x1;
      smax = x2;
    }
  else
    {
      smin = x2;		//added by SS August 04
      smax = x1;		//added by SS August 04
    }


  stuff_phot (p, &ds_to_disk_photon);
  move_phot (&ds_to_disk_photon, smin);

  s_disk = rtsafe (disk_deriv, 0.0, smax - smin, fabs (smax - smin) / 1000.);
  //fabs added by SS August 04 - want convergence test quantity to be +ve 
  s_disk += smin;

  /* Note that s_disk can still be negative */

  /* Now we have the distance to the disk and the only question which rmains
   * is whether we encounter the outside edge of the disk first.  This is pretty
   * unlikely, but ...
   */

  s_sphere = ds_to_sphere (geo.diskrad, p);
  stuff_phot (p, &phit);
  move_phot (&phit, s_sphere);
  if (fabs (phit.x[2]) < geo.disk_z0 * geo.diskrad)
    {
      /* The photon actually hits the disk rim and you must then check whether it hits
       * the disk before the rim */
      if (fabs (s_sphere) < fabs (s_disk))
	return (s_sphere);
    }

  return (s_disk);
}

/* This is the function used by the Recipes routine rtsafe to locate the intercept 
   between the photon and the disk.  The function we are trying to zero is the
   difference between the height of the disk and the z height of the photon. 

   Minor modification by SS August 04 in case z1 = 0.

*/

void
disk_deriv (s, value, derivative)
     double s, *value, *derivative;
{
  struct photon phit;
  double z1, z2, r1, r2;
  double ds;

  stuff_phot (&ds_to_disk_photon, &phit);
  move_phot (&phit, s);
  r1 = sqrt (phit.x[0] * phit.x[0] + phit.x[1] * phit.x[1]);
  z1 = zdisk (r1) - fabs (phit.x[2]);	// this is the function

  /* OK now calculate the derivative */

  ds = (z1) / 100.;
  ds += 1;			// We must move it a little bit (in case z1 = 0) SS Aug 2004
  move_phot (&phit, ds);	// Move the photon a bit more
  r2 = sqrt (phit.x[0] * phit.x[0] + phit.x[1] * phit.x[1]);
  z2 = zdisk (r2) - fabs (phit.x[2]);	// this is the function

  *value = z1;
  *derivative = (z2 - z1) / ds;

}
