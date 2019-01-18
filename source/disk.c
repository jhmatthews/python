
/***********************************************************/
/** @file  disk.c
 * @author ksl
 * @date   January, 2018
 *
 * @brief  Various routines related setting up the disk and calculating
 * when a photon encounters the disk surface
 *
 ***********************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "atomic.h"
#include "python.h"



/**********************************************************/
/** 
 * @brief      Calculate the reference temperarure for a standard Shakura-Sunyaeev disk
 *
 * @param [in] double  m   mass of the central object
 * @param [in] double  mdot   mass accretion rate of the disk
 * @param [in] double  r   radiis of the central object
 * @return     The reference temperature for the disk
 *
 * @details
 * All of the inputs are in cgs units
 *
 * ### Notes ###
 * For an SS disk,  The maximum temperature
 *    actuallly seen in the disk is 0.488 tdisk
 *
 *
 **********************************************************/

double
tdisk (m, mdot, r)
     double m, mdot, r;
{
  double t;
  t = 3. * G / (8. * PI * STEFAN_BOLTZMANN) * m * mdot / (r * r * r);
  t = pow (t, 0.25);
  return (t);
}




/**********************************************************/
/** 
 * @brief      Calculate the temperature of the disk at a normalised distance x
 *
 * @param [in] double  t   reference temperature of the disk in degrees
 * @param [in] double  x   distance from center of disk in units of r/rmin
 * @return     The temperature
 *
 * The routine returns the effective temperature for the disk as a distance
 * which is measured in terms of the ratio of the actual distance to the inner
 * edge of the disk (rmin).
 *
 * A number of modifications to a simple steady state disk are possible, depending
 * on the variable geo.disk_tprofile:
 *
 *  - DISK_TPROFILE_READIN  the returned temperature is interpolated from values that
 *    been read in
 *  - DISK_TPROFILE_YSO an analytic profile implemented for YSOs
 *
 *  In addition for a standard profile, this is where disk heading can optionally
 *  be taken into account
 *
 *
 * ###Notes###
 *
 * A reference for the standard steady state disk is Wade, 1984 MNRAS 208, 381
 *
 **********************************************************/

double
teff (t, x)
     double t, x;
{
  double q = 0;
  double theat, r;
  double pow ();
  double disk_heating_factor;
  int kkk;
  int n;



  if (x < 1)
  {
    Error ("teff: x %f less than 1.0\n", x);
    return (0.0);
  }


  if ((geo.disk_tprofile == DISK_TPROFILE_READIN) && ((x * geo.rstar) < blmod.r[blmod.n_blpts - 1]))
  {
    /* This is the case where the temperature profile is read in as an array */
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
      Error ("tdisk: inside BL profile region but failed to identify temp.\n");
    }
  }
  else
  {
    /* This is a standard accretion disk */

    q = (1.e0 - pow (x, -0.5e0)) / (x * x * x);
    q = t * pow (q, 0.25e0);

    if (geo.absorb_reflect == BACK_RAD_ABSORB_AND_HEAT && geo.wcycle > 0)       /* Absorb photons and increase t so that heat is radiated
                                                                                   but only do this if there has been at least one
                                                                                   ionization cycle */
    {
      r = x * geo.rstar;        // 04aug -- Requires fix if disk does not extend to rstar
      kkk = 1;                  // photon cannot hit the disk at r<qdisk.r[0]
      while (r > qdisk.r[kkk] && kkk < NRINGS - 1)
        kkk++;

      /* Note that disk has 2 sides */
      theat = qdisk.heat[kkk - 1] / (2. * PI * (qdisk.r[kkk] * qdisk.r[kkk] - qdisk.r[kkk - 1] * qdisk.r[kkk - 1]));

      /* T_eff is given by T_eff**4= T_disk**4+Heating/area/STEFAN_BOLTZMANN */
      q = pow (q * q * q * q + (theat / STEFAN_BOLTZMANN), 0.25);

    }
    else if (geo.disk_tprofile == DISK_TPROFILE_YSO)    // Analytic approximation for disk heating by star; implemented for YSOs
    {
      disk_heating_factor = pow (geo.tstar / t, 4.0);
      disk_heating_factor *= (asin (1. / x) - (pow ((1. - (1. / (x * x))), 0.5) / x));
      disk_heating_factor /= PI;
      disk_heating_factor *= x * x * x;
      disk_heating_factor /= (1 - sqrt (1. / x));
      disk_heating_factor += 1;

      q *= pow (disk_heating_factor, (1. / 4.));

    }
  }
  return (q);
}


/**********************************************************/
/** 
 * @brief      Calculate a reference gravity for the disk
 *
 * @param [in] double  mass   Mass of the central object
 * @param [in] double  mdot   Mass accretion rate of the disk
 * @param [in] double  rmin   Minimum radius for the disk
 * @return     A fiducial gravity for the disk
 *
 * The gravity is needed when one constructs a disk spectrum
 * from spectra from a grid of stellar atmospheres
 *
 * The input units are all cgs
 *
 * ###Notes###
 *
 * See Long & Knigge for details
 *
 *
 **********************************************************/

double
gdisk (mass, mdot, rmin)
     double mass, rmin, mdot;
{
  double g0;
  g0 = 0.625 * log10 (mass / MSOL) - 1.875 * log10 (rmin / 1.e9) + 0.125 * log10 (mdot / 1.e16);
  g0 = 5.96e5 * pow (10., g0);
  return (g0);
}


/**********************************************************/
/** 
 * @brief      Calculate the effective gravity at a specfic position in the disk
 *
 * @param [in] double  g0  reference gravity in cm s**-2
 * @param [in] double  x   distance from center in units of r/rmin
 * @return     gravity at x in cm s**-2
 *
 * The gravity is needed when one constructs a disk spectrum
 * from spectra from a grid of stellar atmospheres
 *
 *
 * ###Notes###
 *
 * See Long & Knigge for details
 *
 **********************************************************/

double
geff (g0, x)
     double g0, x;
{
  double q;
  q = (1.0e0 - pow (x, -0.5e0));
  q = pow (x, -1.875e0) * pow (q, 0.125);
  q = g0 * q;
  return (q);
}


double north[] = { 0.0, 0.0, 1.0 };


/**********************************************************/
/** 
 * @brief      Calculate the speed and velocity v at which material given
 *    specific position in the disk
 *
 * @param [in] double  x[]   A position in the disk
 * @param [out] double  v[]   The velocity at that position
 * @return     The speed
 *
 * This routine interpolates the velocities from a grid which gives
 * the run of velocity with radius. Although the we generally
 * assume keplerian velocities for the disk those velocities
 * are not calculated here.
 *
 * ###Notes###
 *
 * The routine projects the input variable x on to the xy plane
 * before it calculates velocities
 *
 * It then takes the xproduct of a unit vector pointing north and
 * the position, to determine the direction of motion at this point
 * in the disk
 *
 *
 *  Finally it rescales this to get the actual velocity.
 *
 *
 *
 *
 **********************************************************/

double
vdisk (x, v)
     double x[];
     double v[];
{
  double xhold[3];
  double r, speed;

  stuff_v (x, xhold);
  xhold[2] = 0.0;
  r = length (xhold);
  linterp (r, disk.r, disk.v, NRINGS, &speed, 0);       //interpolate in linear space
  cross (north, xhold, v);      /* The velocity vector direction is given by north x r */
  renorm (v, speed);
  return (speed);
}




/**********************************************************/
/** 
 * @brief      Calculate the height of the disk at a certain radius
 *
 * @param [in] double  r   a radial position in the disk.
 * @return     The vertical height of the disk at that point in the xy
 * 	plane.
 *
 *
 * If the disk is vertically extended then the disk surface is
 * defined in terms of the h the fraction of the disk radius
 * at the edge of the disk, and a power law exponent.  So
 * a disk thst us a sunoke wedgge wiykd gave an exponent
 * of 1. A flat but thick disk would have an exponent of 0.
 * A flared disk wouuld generally have an exponent greater than
 * 1.
 *
 * ###Notes###
 *
 * 	zdisk returns a number that will positive or zero, so one
 * 	often needs to take this account if for example one has
 * 	a photon that hits the disk below the z=0 plane.
 *
 **********************************************************/

double
zdisk (r)
     double r;
{
  double z;
  z = geo.disk_z0 * pow (r / geo.diskrad, geo.disk_z1) * geo.diskrad;
  return (z);
}



int ds_to_disk_init = 0;
struct photon ds_to_disk_photon;
struct plane diskplane, disktop, diskbottom;


/**********************************************************/
/** 
 * @brief      Calculates the distance that a photon
 *  	would need to travel from its current position to hit disk.
 *
 * @param [in] struct photon *  p   a photon pointer.
 * @param [in] int  allow_negative   permits a negative number to
 *          be returned in certain circumstances
 * @return     The distance to the disk.
 *
 *
 * Usually, ds_to_disk returns the distance along the line of
 * sight to the disk in the direction a photon is currently travelling.
 *
 * Usually, if the photon misses the disk going in the positive
 * direction, a very large number is returned.
 *
 * If allow_negative is non-zero (true), then ds_to_disk returns
 * a negative distance if it has not hit the disk going in the
 * positive direction.
 *
 * The routine allows both for a flat disk, and a vertically
 * extended disk
 *
 *
 *
 * ###Notes###
 *
 * The need to allow for negative distances arises
 * because several of the parameterization for the wind (SV, KWD) depend
 * on the distance between the current position and footpoint
 * along a streamline.  (It is not clear
 * that this was a good approach, but that is the rationale)
 *
 * Generally, there should not be an issue with this, though
 * one can imagine that one might get into trouble for a very
 * flared (bowl-shaped) disk.
 *
 * Significant portions of the routine were rewritten in 1802 by
 * ksl to better handle flared disks.
 *
 **********************************************************/

double
ds_to_disk (p, allow_negative)
     struct photon *p;
     int allow_negative;
{
  double x1, x2;
  double s, s_negative, s_test;
  double s_plane, s_top, s_bottom, s_sphere;
  double s_disk = 0;
  double r, r_top, r_bottom;
  double smin, smax;
  struct photon phit;
  void disk_deriv ();


  if (geo.disk_type == DISK_NONE)
    return (VERY_BIG);          /* There is no disk! */

  /* Initialize 3 structures that define
     the plane of the disk, and two other
     planes that encompass the disk */

  if (ds_to_disk_init == 0)
  {
    diskplane.x[0] = diskplane.x[1] = diskplane.x[2] = 0.0;
    diskplane.lmn[0] = diskplane.lmn[1] = 0.0;
    diskplane.lmn[2] = 1.0;


    disktop.x[0] = disktop.x[1] = 0.0;
    disktop.x[2] = geo.diskrad * geo.disk_z0;
    disktop.lmn[0] = disktop.lmn[1] = 0.0;
    disktop.lmn[2] = 1.0;

    diskbottom.x[0] = diskbottom.x[1] = 0.0;
    diskbottom.x[2] = (-geo.diskrad * geo.disk_z0);
    diskbottom.lmn[0] = diskbottom.lmn[1] = 0.0;
    diskbottom.lmn[2] = 1.0;

    ds_to_disk_init++;          // Only initialize once

  }

  /* Now calculate the place where the photon hits the diskplane */

  s_test = s_plane = ds_to_plane (&diskplane, p);
  stuff_phot (p, &phit);
  move_phot (&phit, s_test);
  r = sqrt (phit.x[0] * phit.x[0] + phit.x[1] * phit.x[1]);


  s = VERY_BIG;

  if (geo.disk_type == DISK_FLAT)
  {
    if (r > geo.diskrad)
      return (VERY_BIG);
    else if (allow_negative)
    {
      return (s_plane);
    }
    else if (s_plane > 0)
    {
      return (s_plane);
      return (VERY_BIG);
    }
  }

  /* At this point, we have completed the simple case.  It is simple
   * because there is only one possible intersection with the disk
   * boundary
   *
   * For the vertically extended disk we have to keep track of
   * the smallest postive value and the smallest (in absolute
   * terms negative value.
   * OK now we have to deal with the hard case.  We would like to
   * avoid actually having to calculate the intercept to the disk
   * if we can because this is likely time consuming. So we first
   * determine this, by checking where the ray hits a sphere that
   * just encompasses the disk
   */

  s = VERY_BIG;
  s_negative = -VERY_BIG;

  s_test = s_top = ds_to_plane (&disktop, p);
  stuff_phot (p, &phit);
  move_phot (&phit, s_test);
  r = r_top = sqrt (phit.x[0] * phit.x[0] + phit.x[1] * phit.x[1]);

  if (r < geo.diskrad)
  {
    if (s_test >= 0)
    {
      s = s_top;
    }
    else
    {
      s_negative = s_top;
    }
  }



  s_test = s_bottom = ds_to_plane (&diskbottom, p);
  stuff_phot (p, &phit);
  move_phot (&phit, s_test);
  r = r_bottom = sqrt (phit.x[0] * phit.x[0] + phit.x[1] * phit.x[1]);

  if (r < geo.diskrad)
  {
    if (s_test >= 0 && s_test < s)
    {
      s = s_test;
    }
    else if (s_test < 0 && s_test > s_negative)
    {
      s_negative = s_test;
    }
  }

  s_test = s_sphere = ds_to_sphere (geo.diskrad, p);
  stuff_phot (p, &phit);
  move_phot (&phit, s_test);

  /* At this point phit is on the sphere but our check here is whether its on the diks or not
   */

  if (fabs (phit.x[2]) < geo.disk_z0 * geo.diskrad)
  {
    if (s_test >= 0 && s_test < s)
    {
      s = s_test;
    }
    else if (s_test < 0 && s_test > s_negative)
    {
      s_negative = s_test;
    }

  }

  /* So at this point we have the smallest positive and smallest negative values.  If we
   * have missed the boundary, we have s as VERY_BIG and s_negative as -VERY_BIG */


  if (s == VERY_BIG)
  {
    if (s_negative == -VERY_BIG || allow_negative == 0)
    {
      return (VERY_BIG);
    }
    else
    {
      // This is the case where we allow s to be negative, and there was no positive intercept
      s_negative = s;
    }
  }


/* At this point we must find the intercept with the disk */


  /*  First we check the unlikely case that the photon hit the edge of the
   *  disk.  This is easy because if so s will be the same as s_sphere
   *  */

  if (s == s_sphere)
  {
    return (s);
  }

  /*  Now we need to find exactly where we have hit the disk.  This
   *  is not trivial.  Basically we need to bracket the possibilites
   *
   * At this point we know the photon is on a path that passes
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

  if (fabs (x2) > fabs (x1))
  {
    smin = x1;
    smax = x2;
  }
  else
  {
    smin = x2;
    smax = x1;
  }


  stuff_phot (p, &ds_to_disk_photon);
  move_phot (&ds_to_disk_photon, smin);

  s_disk = rtsafe (disk_deriv, 0.0, smax - smin, fabs (smax - smin) / 1000.);
  s_disk += smin;

  /* Note that s_disk can still be negative */


  return (s_disk);
}



/**********************************************************/
/** 
 * @brief      Calculate the change in disk height as s function of s
 *
 * @param [in] double  s  A distance along the path of a photon
 * @param [out] double *  value   the disk height at s
 * @param [out] double *  derivative   the derivative of the disk heith at s
 * @return     N/A
 *
 * Used in ds_to_disk by rtsafe as part of the procedure to locate
 * the place a photon hits a vertically extended disk
 *
 * ###Notes###
 *
 *
 *  The function we are trying to zero is the
 *  difference between the height of the disk and the z height of the photon.
 *
 **********************************************************/

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
  z1 = zdisk (r1) - fabs (phit.x[2]);   // this is the function

  /* OK now calculate the derivative */

  ds = (z1) / 100.;
  ds += 1;                      // We must move it a little bit (in case z1 = 0) SS Aug 2004
  move_phot (&phit, ds);        // Move the photon a bit more
  r2 = sqrt (phit.x[0] * phit.x[0] + phit.x[1] * phit.x[1]);
  z2 = zdisk (r2) - fabs (phit.x[2]);   // this is the function

  *value = z1;
  *derivative = (z2 - z1) / ds;

}







/**********************************************************/
/** 
 * @brief      Initialize a structure (qdisk) for recording information about photons/energy impinging
 * 	on the disk, which is stored in a disk structure called qdisk.
 *
 * @return     Always return zero
 *
 *
 * ###Notes###
 *
 * The information stored in qdisk can be used to modify the effective temperature
 * of the disk
 *
 * The reason the qdisk is needed as well as disk structure is that whenever the
 * 	wavelength bands are changed the radii in the disk structure are recalibratied.
 * 	We want qdisk to have fixed boundaries when this happens.
 *
 **********************************************************/

int
qdisk_init ()
{
  int n;
  for (n = 0; n < NRINGS; n++)
  {
    qdisk.r[n] = disk.r[n];
    qdisk.t[n] = disk.t[n];
    qdisk.g[n] = disk.g[n];
    qdisk.v[n] = disk.v[n];
    qdisk.heat[n] = 0.0;
    qdisk.nphot[n] = 0;
    qdisk.nhit[n] = 0;
    qdisk.w[n] = 0;
    qdisk.ave_freq[n] = 0;
    qdisk.t_hit[0] = 0;
  }
  return (0);
}


/**********************************************************/
/** 
 * @brief      Save the information in the qdisk structure to a file
 *
 * @param [in] char *  diskfile   Name of the file whihc is writteen
 * @param [in] double  ztot   The total luminosity of the disk as
 *                      calculate over multiple cycles
 * @return     Always returns 0
 *
 * The routine reformats the data about disk heating which has
 * been accumulated and writes it to a file
 *
 * ###Notes###
 *
 * The data concerning heating by the disk is built up during
 * a the ionization cycles
 *
 * The file that is produced should be readable as an astropy
 * table
 *
 **********************************************************/

int
qdisk_save (diskfile, ztot)
     char *diskfile;
     double ztot;
{
  FILE *qptr;
  int n;
  double area, theat, ttot;
  qptr = fopen (diskfile, "w");
  fprintf (qptr, "r         zdisk     t_disk   heat       nhit nhit/nemit  t_heat    t_irrad  W_irrad  t_tot\n");

  for (n = 0; n < NRINGS; n++)
  {
    area = (2. * PI * (qdisk.r[n + 1] * qdisk.r[n + 1] - qdisk.r[n] * qdisk.r[n]));
    theat = qdisk.heat[n] / area;
    theat = pow (theat / STEFAN_BOLTZMANN, 0.25);       // theat is temperature if no internal energy production
    if (qdisk.nhit[n] > 0)
    {

      qdisk.ave_freq[n] /= qdisk.heat[n];
      qdisk.t_hit[n] = H * qdisk.ave_freq[n] / (BOLTZMANN * 3.832);     // Basic conversion from freq to T
      qdisk.w[n] = qdisk.heat[n] / (4. * PI * STEFAN_BOLTZMANN * area * qdisk.t_hit[n] * qdisk.t_hit[n] * qdisk.t_hit[n] * qdisk.t_hit[n]);
    }

    ttot = pow (qdisk.t[n], 4) + pow (theat, 4);
    ttot = pow (ttot, 0.25);
    fprintf (qptr,
             "%8.3e %8.3e %8.3e %8.3e %5d %8.3e %8.3e %8.3e %8.3e %8.3e\n",
             qdisk.r[n], zdisk (qdisk.r[n]), qdisk.t[n],
             qdisk.heat[n], qdisk.nhit[n], qdisk.heat[n] * NRINGS / ztot, theat, qdisk.t_hit[n], qdisk.w[n], ttot);
  }

  fclose (qptr);
  return (0);
}





/**********************************************************/
/** 
 * @brief      Read the temperature profile from a file
 *
 * @param [in, out] char *  tprofile   Name of the input file
 * @return     Always returns 0
 *
 * The input format for the file to be read in is quite spefic
 * The first line should contina the number of points n
 * the profile
 *
 * Subsequent lines should contian a radius and a temperarue
 *
 * The radius values shoule be in units of 1e11 cm
 * The temperature should be in units of 1000K
 *
 * ###Notes###
 *
 * @bug This routine which was written for the YSO study
 *  needs to be made less YSO centric. It should also be
 *  retested.
 *
 **********************************************************/

int
read_non_standard_disk_profile (tprofile)
     char *tprofile;
{

  FILE *fopen (), *fptr;
  int n;
  float dumflt1, dumflt2;
  int dumint;

  if ((fptr = fopen (tprofile, "r")) == NULL)
  {
    Error ("Could not open filename %s\n", tprofile);
    exit (0);
  }

  fscanf (fptr, "%d\n", &dumint);
  blmod.n_blpts = dumint;
  for (n = 0; n < blmod.n_blpts; n++)
  {
    fscanf (fptr, "%g %g", &dumflt1, &dumflt2);
    blmod.r[n] = dumflt1 * 1.e11;
    blmod.t[n] = dumflt2 * 1.e3;
  }

  fclose (fptr);

  return (0);
}
