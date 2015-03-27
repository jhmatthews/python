
/**************************************************************************
                    Space Telescope Science Institute


  Synopsis:

   These routines are used to describe the Roche geometry and are intended to allow one
   to calculate whether a photon intersects the surface of the star 

  Description:

   The important routines are 

   binary_basics  which uses the masses of the primary and the secondary to calculate
   the position of various Lagrange points, etc, and the size of a pillbox
   which encloses the Roche lobe of the secondary.  It must be called before
   the other routines can be used!

   ds_to_roche2   calculates the distance a photon defined by p_roche has to travel 
   before it will encounter the secondary.  If it misses, ds_to_roche2
   returns VERY_BIG.

   hit_secondary   determines whether a photon defined by p_roche will hit the secondary
   or not.  If it misses, hit_secondary returns 0, if it hits, the return is
   P_SEC.


  Arguments:		

  Returns:

  Notes:

   The secondary is located along the x, i.e. 0, axis

   To utilize the routines you must stuff the photon of interest into the external photon "p_roche".
   This is necessary because the program makes use of two Numerical Recipes Routines rtsafe
   and golden to find zeroes and minima respectively.  Note that it is possible (likely) that these are
   not the very best choices of Numerical Recipes routines to use.  In particular it is possible that
   BRENT could replace both routines.

   The routines adopt the same approach as in Keith Horne's routines dealing with the Roche
   geometry, that is a pillbox is defined around the secondary.  The radius of the pillbox is 
   determined by the maximum size of the secondary perpendicular to the line of sight between
   the two stars.  One end of the pillbox is defined by the L1 point.  The other end of the pillbox
   is defined by the backside of the secondary.  One first tests to find out whether a photon
   intersects the pillbox, and then one checks more carefully to see whether the photon actually
   hits the secondary (using the intersections with the pillbox to define the limits of the careful
   search.

   Note: ds_to_roche2 is not actually used by Python.  It would be used if one wanted to calculate
   the heating of the secondary.

  History:

   97jul   ksl	Coded       
   01jul   ksl  Created a program plot_roche, which uses these subroutines
                Improved the notes on this.  Replace the pillbox routine so
		that it was simpler and more correct, especially for photons
		which were outside the L1 plane.  It is still true that photons
		which move directly along the x axis, and perhaps others may
		present problems.

 ************************************************************************/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "atomic.h"
#include "python.h"


/* Given the mass of the primary and secondary and the period calculate some other basic parameters
   for the binary system, including the L1 point. */

/* p_roche provides a way to pass a photon structure to roche() and roche_derive which
   are in turn used by the Numerical Recipes Routine rtsafe to find the zero either to dphi_ds
   or to phi along the path lenght of the photon */

struct photon p_roche;

int
binary_basics ()
{
  double x;

  void find_l1 (), roche (), roche_deriv ();
  double rtsafe ();

  double phi (), dphi_ds (), d2phi_ds2 ();
  double roche2_width_max ();
  double pow ();

  /* Calculate a the separation of the primary and the secondary */

  x = G * geo.period * geo.period * (geo.mstar + geo.m_sec) / (16. * PI * PI);
  geo.a = pow (x, 1. / 3.);

  /*Calculate the mass ratio.  */

  geo.q = geo.m_sec / geo.mstar;

  /* Define the position of the secondary with respect to the primary.  */

  plane_sec.x[0] = geo.a;
  plane_sec.x[1] = plane_sec.x[2] = 0;
  plane_sec.lmn[0] = 1;
  plane_sec.lmn[1] = plane_sec.lmn[2] = 0;

  /* Find the position of the L1 point with respect to the primary */

  geo.l1 = 0;
  geo.phi = 0;

  /* p_roche is used here to provide an external way to pass a photon structure to find_l1 */
  /* p_roche is (near) origin of primary and pointing at secondary */

  p_roche.x[0] = 1000;
  p_roche.x[1] = p_roche.x[2] = 0;
  p_roche.lmn[0] = 1;
  p_roche.lmn[1] = p_roche.lmn[2] = 0;

  geo.l1 = x =
    rtsafe (roche_deriv, 0.01 * geo.a, 0.99 * geo.a, geo.a / 1000.);

  geo.l1_from_m2 = geo.a - geo.l1;

  plane_l1.x[0] = geo.l1 + 10000.;
  plane_l1.x[1] = plane_l1.x[2] = 0;
  plane_l1.lmn[0] = 1;
  plane_l1.lmn[1] = plane_l1.lmn[2] = 0;

  /* Similarly, find l2, the Lagrangian point behind the secondary */

  geo.l2 = x = rtsafe (roche_deriv, 1.01 * geo.a, 2.0 * geo.a, geo.a / 1000.);

  /* Now find the position of the far side of the star */

  geo.phi = phi (geo.l1);	/* Set geo.phi so phi will be zero on Roche lobes */


  /* Set geo.r2_far to be the radius of the secondary on the backside of the secondary */

  geo.r2_far = rtsafe (roche, 1.01 * geo.a, geo.l2, geo.a / 1000.) - geo.a;

  /* Define a plane on the backside of the secondary */

  plane_m2_far.x[0] = geo.r2_far + geo.a;
  plane_m2_far.x[1] = plane_m2_far.x[2] = 0;
  plane_m2_far.lmn[0] = 1.0;
  plane_m2_far.lmn[1] = plane_m2_far.lmn[2] = 0;


  Log_silent
    ("binary_basics: l1=%8.2e l2=%8.2e l1_from_m2=%8.2e r2_far %8.2e\n",
     geo.l1, geo.l2, geo.l1_from_m2, geo.r2_far);

  /* Calculate the maximum half width of the secondary in the plane of the orbit */
  geo.r2_width = roche2_width_max ();
  Log_silent ("binary_basics: r2_width=%6.2e\n", geo.r2_width);

  return (0);
}



/* Calculate the distance to the roche surface of the secondary for photon p.  
   Return -VERY_BIG if the photon does not intercept the surface.  

 */

double
ds_to_roche_2 (p)
     PhotPtr p;
{


  void roche ();
  double rtsafe ();
  double phi (), dphi_ds (), d2phi_ds2 ();
  double ds_to_sphere2 ();
  struct photon pp;

  double smin, smax, s;
  double golden (), phi ();
  double potential;

  double pillbox ();

  if ((s = pillbox (p, &smin, &smax)) == VERY_BIG)
    return (0);			/* Missed secondary */
  stuff_phot (p, &p_roche);

// So, after next line,  s is the minimum of the potential, between smin and smax
// If this value is > 0, then the photon did not hit the secondary

  potential = golden (smin, 0.5 * (smin + smax), smax, phi, 0.0001, &s);

  if (potential > 0.0)
    return (VERY_BIG);		/*Missed secondary) */


  if (smin > s)
    Log ("Error Knoxy %e %e %e\n", smin, s, smax);

  if ((s = rtsafe (roche, smin, s, geo.a / 1000.)) < 0)
    {
      Log ("Problem located Roche surface s %6.2e\n", s);
      Log ("smin %6.2e  phi %6.2e %6.2e %6.2e\n", smin, phi (smin),
	      dphi_ds (smin), d2phi_ds2 (smin));
      Log ("smax %6.2e  phi %6.2e %6.2e %6.2e\n", smax, phi (smax),
	      dphi_ds (smax), d2phi_ds2 (smax));
      stuff_phot (p, &pp);
      move_phot (&pp, smin);
      Log ("pmin %6.2e %6.2e %6.2e\n", pp.x[0], pp.x[1], pp.x[2]);
      stuff_phot (p, &pp);
      move_phot (&pp, smax);
      Log ("pmax %6.2e %6.2e %6.2e\n", pp.x[0], pp.x[1], pp.x[2]);

      return (s);		/* This is an error return */
    }

  Log ("Distance to roche surface  %6.2e\n", s);


  return (s);
}

/* Find out whether a photon will hit the secondary, but don't actually calculate the
   point where it hits */

int
hit_secondary (p)
     PhotPtr p;
{
  double smin, smax, s;
  double golden (), phi ();
  double potential;

  double pillbox ();

  if (pillbox (p, &smin, &smax) == VERY_BIG)
    return (0);			/* Missed secondary */
  stuff_phot (p, &p_roche);
  potential = golden (smin, 0.5 * (smin + smax), smax, phi, 0.0001, &s);
  if (potential > 0.0)
    return (0);			/*Missed secondary) */


  return (P_SEC);
}


/* Find out whether the photon hits the pillbox which surrounds the Roche surface of the secondary.

   The pillbox is defined by a plane at L1, a plane on the backside of the star and a cylinder whose radius
   is the maximum (half) width of the star in the plane of the orbit.

   This routine distinguishes two basic types of photons, "normal" and "abnormal" photons.  A normal photon
   is on the WD side of L1 plane; an abnormal one is on the side of the L1 plane with the secondary.
 

   pillbox returns VERY_BIG if the photon p has not intersected the pillbox.  It appears to return the distance
   to the pillbox if the photon p hits the pillbox.  smin and smax contain the entrance and exit distances.  It
   is not obvious that smin < smax.

Note:   01jul22	I revised this so that it would handle photons coming in all directions for use with plot_roche
but this instantiation of pillbox still seems to me not optical.  The basic point is that unless the photon
is in the pillbox, the only roots that matter are ones which are positive.  Only inside the pillbox will
there be a positive and a negative root.  Here we try to correct this at the end by making sure that both
roots are not negative, but one would have thought that one might establish that that photons inside the
pillbox are special and dealt with them accordindly

History:

	01jul	ksl	Completely revised in the process of creating plot_roche
	01jul22	ksl	Fixed a problem which allowed photons headed completely away from the pillbox
			to be counted as hitting the pillbox.  See note above
*/

double
pillbox (p, smin, smax)
     PhotPtr p;
     double *smin, *smax;
{
  double x1, x2;
  double a, b, c;
  double root[2], ss[2];
  int i, n;
  struct photon pp;
  double ds_to_plane ();


  n = 0;			// At this point no boundaries to the pillbox have been identified

/* If the photon is along the x-axis, then it will not intersect the cylinder anywhere, and
the only possibility is that it hit the endcaps of the cylinder. But normally, it will
not be along the x axis and one must find if there are any intersections with the cylinder
and determine where they are.  This calculation is carried out in the if section below. */

  if ((a = (p->lmn[0] * p->lmn[0])) < 1.0)
    {
      a = 1. - a;
      b = 2. * (p->x[1] * p->lmn[1] + p->x[2] * p->lmn[2]);
      c = p->x[1] * p->x[1] + p->x[2] * p->x[2] - geo.r2_width * geo.r2_width;

/* After "quadratic", root contains the distances the photon pp would need to travel to 
hit the edges of the cylinder. i is an index to the smallest positive root, if one exists.  
If both of the roots are negative or imaginary then i will be negative, and we know
that the ray did not hit the cylinder while travelling in the positive direction.  In
that case return  VERY_BIG

Note that the fact that the ray does hit the cylinder going in the positive direction 
does not necessarily mean that it hits the pillbox.  
*/

      if ((i = quadratic (a, b, c, root)) < 0)
	return (VERY_BIG);

/* If the intersection is in the part of the cylinder between the caps it is probably legitimate.
So tranport pp to the intersection of the first root and check where that lies. And then
do the same thing for the second root.  The "probably refers to the fact that a negative
root is really only possible if the photon is already in the pillbox  */

      stuff_phot (p, &pp);
      move_phot (&pp, root[0]);

      if (geo.l1 <= pp.x[0] && pp.x[0] <= plane_m2_far.x[0])
	{
	  ss[n] = root[0];
	  n++;
	}


      stuff_phot (p, &pp);
      move_phot (&pp, root[1]);	/* So pp is now located at intersection of the second root */
      if (geo.l1 <= pp.x[0] && pp.x[0] <= plane_m2_far.x[0])
	{
	  ss[n] = root[1];
	  n++;
	}

    }

/* Now find out if the photon hits the caps. This is done by translating to the
planes of the two end caps and then checking the cylindrical radius. Negative
distances are valid only if the photon is in the pillbox already */

  x1 = ds_to_plane (&plane_l1, p);	/* Calculate the distance to the L1 plane */
  stuff_phot (p, &pp);
  move_phot (&pp, x1);		/* So pp is now located at the l1 plane */
  if ((pp.x[1] * pp.x[1] + pp.x[2] * pp.x[2]) <= geo.r2_width * geo.r2_width)
    {
      ss[n] = x1;
      n++;
    }

  x2 = ds_to_plane (&plane_m2_far, p);	/* Calculate the distance to the plane behind the the star.  */
  stuff_phot (p, &pp);
  move_phot (&pp, x2);		/* So pp is now located at the r2_far plane */
  if ((pp.x[1] * pp.x[1] + pp.x[2] * pp.x[2]) <= geo.r2_width * geo.r2_width)
    {
      ss[n] = x2;
      n++;
    }

/* Finally sort it all out.  If n == 0, or if both s[0] and s[2] are negative 
then the photon did not hit the pillbox ksl 02jan */

  if ((n == 0) || (ss[0] < 0 && ss[1] < 0))
    return (VERY_BIG);		// The photon did not hit the pillbox

  if (n == 2)
    {
      if (ss[0] < ss[1])
	{
	  *smin = ss[0];
	  *smax = ss[1];
	}
      else
	{
	  *smin = ss[1];
	  *smax = ss[0];
	}

      return (*smin);
    }

  Error ("pillbox %d interfaces to pillbox is impossible\n", n);
  return (VERY_BIG);

}


/* Calculate the Roche potential at position indicated by the photon p_roche and the length s */


int phi_init = 0;
double phi_gm1, phi_gm2, phi_3, phi_4;

double
phi (s)
     double s;
{
  struct photon pp;
  double x1, x2, z, z1, z2, z3;
  double length ();

  if (phi_init == 0)
    {
      phi_gm1 = G * geo.mstar;
      phi_gm2 = G * geo.m_sec;
      phi_3 = (phi_gm1 + phi_gm2) / (geo.a * geo.a * geo.a);
      phi_3 = 0.5 * phi_3;
      phi_4 = geo.a * geo.m_sec / (geo.mstar + geo.m_sec);
//              Log_silent("phi: phi_gm1 %e phi_gm2 %e phi_3 %e phi_4 %e\n",phi_gm1,phi_gm2,phi_3,phi_4);
      phi_init++;
    }

  stuff_phot (&p_roche, &pp);
  move_phot (&pp, s);		/* So now we have the actuaal position of the photon relative to the WD */

  if ((x1 = length (&pp.x[0])) == 0)
    return (-VERY_BIG);
  z1 = -phi_gm1 / x1;
  z3 =
    -phi_3 * ((pp.x[0] - phi_4) * (pp.x[0] - phi_4) + pp.x[1] * pp.x[1]) -
    geo.phi;

  pp.x[0] -= geo.a;		/* Here we make pp refer to the positions w.r.t. the secondary */
  if ((x2 = length (&pp.x[0])) == 0)
    return (-VERY_BIG);

  z2 = -phi_gm2 / x2;

  z = z1 + z2 + z3;

  if (z < -1.e20)
    z = -1.e20;

  return (z);


}


#define EPS 10000.

double
dphi_ds (s)
     double s;
{
  double phi (), x1, x2;
  double dx, z;
  if ((dx = 0.001 * geo.a) < EPS)
    dx = EPS;

  x2 = phi (s + dx);
  x1 = phi (s);
  z = (x2 - x1) / dx;
  return (z);

}


double
d2phi_ds2 (s)
     double s;
{
  double phi (), x1, x2, x3;
  double dx;
  if ((dx = 0.001 * geo.a) < EPS)
    dx = EPS;

  x1 = phi (s);
  x2 = phi (s + dx);
  x3 = phi (s + 2. * dx);
  return ((x3 - 2. * x2 + x1) / (dx * dx));

}


/* Calculate the half width of the Roche lobe at the position x along the 
   x axis in the plane of the orbit.  There is no real guarantee that this would
   work if you were outside L2 or L3. */

double
roche_width (x)
     double x;
{
  double rho, smax;
  void roche ();
  double rtsafe ();

  if (x < geo.l1)
    smax = geo.l1;
  else
    smax = geo.l1_from_m2;

  p_roche.x[0] = x;
  p_roche.x[1] = p_roche.x[2] = 0.0;
  p_roche.lmn[0] = 0;
  p_roche.lmn[1] = 1;
  p_roche.lmn[2] = 0;

  rho = rtsafe (roche, 1000., smax, geo.a / 1000.);

  if (rho < 0)
    {
      Error ("roche_with : rtsafe failure x=%6.2e\n", x);
    }
  return (-rho);		/* This is because we are going to search for a minimun not a maximum */

}

/* Find the maximum half width of the Roche lobe of the secondary.   This routine
   uses the NR routine golden to find the point in x where the Roche lobe is maximized.  */

double
roche2_width_max ()
{
  double xmin, xmax, xmid, xbest;
  double rmin;
  double roche_width (), golden ();

  xmin = geo.l1 + 1.e5;
  xmax = geo.r2_far + geo.a - 1.e5;
  xmid = geo.a;
//      xmid=(geo.l1+geo.a)*0.5;

  Log_silent ("roche2_width_max: Search from %6.2e %6.2e %6.2e\n", xmin, xmid,
	      xmax);

  /* golden returns the value of the function it is evaluating; xbest is the position of the maximum */

  rmin = golden (xmin, xmid, xmax, roche_width, 0.0001, &xbest);

  Log_silent ("roche2_width_max: Max width at %6.2e of %6.2e\n", xbest,
	      -rmin);

  return (-rmin);		/* Because "golden" is designed to find a minimum rather than a maximum! */
}


/* This is the function used by the Recipes routine rtsafe to locate specific
   surfaces of the Roche potential.  Here *value is phi and *derivative is
   the slope of phie at s along the vector defined by p_roche.  The value of
   phi has to be set in advance through geo.phi and would normally be that of
   the Roche lobe */

void
roche (s, value, derivative)
     double s, *value, *derivative;
{
  *value = phi (s);
  *derivative = dphi_ds (s);

}

/* This is the function used by the Recipes routine rtsafe to lacate extrema of
   the Roche potential.  Thus it can be used to find the L1 point, for example.  Here
   *value is the derivative of phi at s along the vector defined by p_roche, and 
   *derivative is the second derivative */

void
roche_deriv (s, value, derivative)
     double s, *value, *derivative;
{
  *value = dphi_ds (s);
  *derivative = d2phi_ds2 (s);

}
