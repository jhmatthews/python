/* 

   This file should contain utilities which can be used to work with photon structues.  It also
   contains routines to calculate the distances to intersections of spheres and cones.  The
   routines included here should not be tightly tied to a particular wind model.   The one
   exception at present is ds_to_cone which calculates the distance to a bifurcated cone.

   1409 - Added standard headers

 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include "atomic.h"

#include "python.h"

int init_stuff_phot = 0;
size_t sizeofphot;





/**************************************************************************


  Synopsis:  

  Simply copy one photon bundle into another

  Description:	

  Arguments:  

  Returns:

  Notes:

  History:

 ************************************************************************/
int
stuff_phot (pin, pout)
     PhotPtr pin, pout;
{
  pout->x[0] = pin->x[0];
  pout->x[1] = pin->x[1];
  pout->x[2] = pin->x[2];

  pout->w = pin->w;
  pout->freq = pin->freq;
  pout->tau = pin->tau;
  pout->path = pin->path;   /* SWM 31/7/14 - Added */

  pout->lmn[0] = pin->lmn[0];
  pout->lmn[1] = pin->lmn[1];
  pout->lmn[2] = pin->lmn[2];

  pout->istat = pin->istat;
  pout->nres = pin->nres;
  pout->nrscat = pin->nrscat;
  pout->nscat = pin->nscat;
  pout->grid = pin->grid;
  pout->origin = pin->origin;
  pout->nnscat = pin->nnscat;
  pout->np = pin->np;



  return (0);
}





/**************************************************************************


  Synopsis:  

  move a photon by a distance ds

  Description:	

  Arguments:  

  Returns:

  Notes:

  History:

 ************************************************************************/

int
move_phot (pp, ds)
     PhotPtr pp;
     double ds;
{

  pp->x[0] += pp->lmn[0] * ds;
  pp->x[1] += pp->lmn[1] * ds;
  pp->x[2] += pp->lmn[2] * ds;
  pp->path += fabs(ds);         /* SWM 31/7/14 - Added */
  return (0);
}


/**************************************************************************


  Synopsis:  

	Compares two photons to determine whether they have
 	the identical position and direction.

  Description:	

  Arguments:  

  Returns:

  Notes:

  History:
	02jun	ksl	Coded

 ************************************************************************/

int
comp_phot (p1, p2)
     PhotPtr p1, p2;
{
  if (p1->x[0] != p2->x[0] || p1->lmn[0] != p2->lmn[0])
    return (1);
  if (p1->x[1] != p2->x[1] || p1->lmn[1] != p2->lmn[1])
    return (1);
  if (p1->x[2] != p2->x[2] || p1->lmn[2] != p2->lmn[2])
    return (1);
  return (0);
}


/**************************************************************************


  Synopsis:  

  phot_hist is designed to record the history of a single photon 
  bundle, in terms of a series of photon structures that are
  populated as the photon goes through the grid

  Unless phot_hist_on is true, then this routine is a NOP


  Description:	

  This is a diagnostic routine only

  Arguments:  

  Returns:

  Notes:

  History:
	090125	ksl	68b -Created to better understand where photon bundles were
			being absorbed or losing energy in the grid

 ************************************************************************/
int
phot_hist (p, iswitch)
     PhotPtr p;
     int iswitch;
{
  if (phot_hist_on == 0)
    return (0);

  if (iswitch == 0)
    {
      n_phot_hist = 0;
    }

  if (n_phot_hist < MAX_PHOT_HIST)
    {
      stuff_phot (p, &xphot_hist[n_phot_hist]);
      n_phot_hist++;
    }
  else
    {
      Error
	("phot_hist: The number of steps %d in phot_hist exceeds MAX_PHOT_HIST\n",
	 MAX_PHOT_HIST);
    }

  return (n_phot_hist);
}


/**************************************************************************


  Synopsis:  
	The next routine is destigned to update a portion of the PlasmaPtr to reflect where
 	photons along the line of sight to the observer were absorbed in the wind

  Description:	

  Arguments:  

  Returns:

  Notes:
	As photon is extracted from the wind, tau changes due to scatters and w changes due
	to absorption.  We are just recording how much energy is absorbed by scatterng processes
	here, and so the energy absorbed is the current weight (exp(-tau_before) - exp (-tau_after))

  History:
	090211	ksl	Created to store the energy removed photons headed toward the observer
			by an ion in a particular cell of the wind.
	0904	ksl	68c - Fixed problem concerning where energy was being stored

 ************************************************************************/

int
phot_history_summarize ()
{
  int n;
  PlasmaPtr xplasma;
  PhotPtr p;
  double x;
  double tau, tau_old;
  int nion;

  p = &xphot_hist[0];
  tau_old = p->tau;

  for (n = 1; n < n_phot_hist; n++)
    {
      p = &xphot_hist[n];

      tau = p->tau;		// tau is tau after the scatter

      nion = lin_ptr[p->nres]->nion;	// ion that scattered

      x = p->w * (exp (-tau_old) - exp (-tau));	// energy removed by scatter

      xplasma = &plasmamain[wmain[p->grid].nplasma];	// pointer to plasma cell where scattering occured

      xplasma->xscatters[nion] += (x);

      tau_old = tau;

      if (xplasma->xscatters[nion] < 0.0)
	{
	  Error
	    ("phot_history_summarize:  n %d n_phot_hist %d phot_hist %d nplasma %d\n",
	     n, n_phot_hist, p->grid, wmain[p->grid].nplasma);
	}
    }



  return (0);
}



/*******************************************************
            Space Telescope Science Institute
                                                                                                                                 
Synopsis:
   ds_to_cone (cc, p) This is a new version of ds_to_cone introduced in order to 
   allow cylvar coordinates  Eventually it is intended to replace ds_to_windcone 
   completely, or at least the guts of ds_to_windcone. 
							                                                                                                                                  
Returns:
                                                                                                                            
                                                                                                                            
Description:
                                                                                                                            
                                                                                                                            
Notes:

As part of the conversion, the old ds_to_cone was renamed to ds_to_windcone
   


	56d -- Cones have become more prominent in python with time.  Originally
	they were used simply to define the inner and outer edges of the
	wind, that is   to determine when a photon entered the wind.  The wind was/is 
	assumed to be a biconical flow, and therefore the most naturally was defined in
	terms of a footpoint in the disk for the innermost/outermost stream line
	and a slope.  
	
	It was expanced when rtheta coordanates were introduced
	to handle find when one had hit the theta coordinate boundaries.  And in
	56d we also want it to handle the case where the grid cells have variable
	z coordinates.  Therefore one needs to be quite careful, Since ds_to_cone
	is the ultimate name of the routine I would like to use, the first step
	I took was to change the name of the old routine to ds_to_windcone.  This
	should allow me to replace the old routine incrementally.  

	56d -- With python 56d the variable defining a cone were changed to be
	the point on the z axis intercepted by the cone (extended) through 
	the disk, and the slope.

	Our approach is as follows:

	The cone is defined by  z=z_o + dz/drho *drho
	The photon is defined by X=X_o + LMN s

	Unless dz/drho = 0, this turns into a simple quadratic equation that must
	be solved.  We work in the northen hemisphere only.
       
	There are issues concerning what to do so one crosses the disk plane.  For
        finding how far a photon can travel in a cell, we are *not* interested 
	in the possibility that the photon hits the cone on the other side of
	the disk plane, but for a photon travelling out of the wind we are 
	very interested in this possibility.	

History:
	05jul	ksl	Created so that cylvar coordinates could be incorporated
			into python.
	06sep	ksl	57h -- Corrected an error discovered by SS in the setup of the
			quadratic for a cone.  Note that if this routing ever became
			a "heavy hitter" in terms of efficiency, there are a few
			changes that could speed it up a bit since the magnitude
			of lmn is always 1.
*/


double
ds_to_cone (cc, p)
     ConePtr cc;
     struct photon *p;
{
  double dz, dzdr2;
  double a, b, c, root[2];
  double s_to_zero;		/* The path length to the xy plane */
  int i;
  struct photon pp;

  /* First of all let's work only in the "northern" hemisphere */
  stuff_phot (p, &pp);

  if (pp.x[2] < 0.0)
    {				/*move the photon to the northen hemisphere */
      pp.x[2] = -pp.x[2];
      pp.lmn[2] = -pp.lmn[2];
    }

  /* Set up and solve the quadratic equation that gives the cone intercept */

  dzdr2 = cc->dzdr * cc->dzdr;
  dz = pp.x[2] - cc->z;

  a =
    dzdr2 * (pp.lmn[0] * pp.lmn[0] + pp.lmn[1] * pp.lmn[1]) -
    (pp.lmn[2] * pp.lmn[2]);
  b =
    2. * (dzdr2 * (pp.lmn[0] * pp.x[0] + pp.lmn[1] * pp.x[1]) -
	  pp.lmn[2] * dz);
  c = dzdr2 * (pp.x[0] * pp.x[0] + pp.x[1] * pp.x[1]) - dz * dz;

  i = quadratic (a, b, c, root);	/* root[i] is the smallest positive root unless i is
					   negative in which case either both roots were negative or both roots were imaginary */

  /* Calculate the positive path length to the xy plane.  If either the
     photon is travelling in the xy plane or if the intercept is in the negative
     direction set the path length to infinity */

  if (pp.x[2] * pp.lmn[2] >= 0)
    s_to_zero = VERY_BIG;
  else
    s_to_zero = (-pp.x[2] / pp.lmn[2]);


  if (i >= 0 && root[i] < s_to_zero)
    return (root[i]);		/*Because that implies
				   the ray hit the cone before hitting the xy plane */
  return (s_to_zero);


}


/**************************************************************************


  Synopsis:  

	Calculate the pathlenth along a line of sight defined by
	a photon p to a sphere or radius r centered on the origin.  If
	the photon does not hit the sphere return a large number VERY_BIG

  Description:	

  Arguments:  

  Returns:

  Notes:

  History:

 ************************************************************************/
double
ds_to_sphere (r, p)
     double r;
     struct photon *p;
{
  double a, b, c, root[2];
  int i;
  double dot ();

  a = 1.;
  b = 2. * dot (p->x, p->lmn);
  c = dot (p->x, p->x) - r * r;

  i = quadratic (a, b, c, root);
/* root[i] is the smallest positive root unless i is
negative in which case either both roots were negative or 
both roots were imaginary */


  if (i >= 0)
    return (root[i]);		/*Because that implies
				   the ray hit the sphere */

  return (VERY_BIG);
}




/**************************************************************************


  Synopsis:  
	This is more generalized routine to find the positive distance to
   	a sphere centered at x with radius r

  Description:	

  Arguments:  

  Returns:

  Notes:

  History:

 ************************************************************************/

double
ds_to_sphere2 (x, r, p)
     double x[], r;
     struct photon *p;
{
  double a, b, c, root[2], delta[3];
  int i;
  double dot ();

  vsub (p->x, x, delta);

  a = 1.;
  b = 2. * dot (delta, p->lmn);
  c = dot (delta, delta) - r * r;

  i = quadratic (a, b, c, root);	/* root[i] is the smallest positive root unless i is
					   negative in which case either both roots were negative or both roots were imaginary */


  if (i >= 0)
    return (root[i]);		/*Because that implies
				   the ray hit the sphere */

  return (VERY_BIG);
}



/**************************************************************************


  Synopsis:  
 
	This solves a simple quadratic (or if a is zero linear equation).  The return is set up
   to make it easy to identify the smallest positive root if one exists.  The routine returns
   a negative number if both roots are negative or imaginary. 
   More specifically 
	 -1 -> both roots are imaginary
         -2 -> both roots are negative or 0
          0 -> the first root is the smallest positive  root 
          1 -> the second root is the smallest positive root


  Description:	

  Arguments:  

  Returns:

  Notes:

  History:
	05jul	ksl	56d -- Heretofore roots of 0 were on no insterest, but now that 
			quadratic is called by bilin, we would like to know about values
			of zero, specifically.  Since the simplest thing to do was to
			check both roots in this case, I added a little code to make
			sure root was defined in all cases.

 ************************************************************************/
int
quadratic (a, b, c, r)
     double a, b, c, r[];
{
  double q, z;

  if (a == 0.0)
    {				/* Then it's not really a quadratic but we can solve it
				   anyway */
      if (b == 0.0)
	{
	  r[0] = r[1] = -99.;
	  return (-1);		/* The roots are extremely imaginary, since both a a b were 0 */
	}

      r[0] = r[1] = (-c / b);	// Then it was a linear equation. Setting both roots to the same thing could be a problem ksl
      if (r[0] < 0.0)
	return (-2);		/* Generally speaking we are not interested in
				   negative distances */
      else
	return (0);
    }

  if ((q = b * b - 4. * a * c) < 0.0)
    {
      r[0] = r[1] = -99.;
      return (-1);		/* both roots are imaginary */
    }
  q = sqrt (q);
  z = 0.5 / a;

  r[0] = (-b - q) * z;
  r[1] = (-b + q) * z;


  if (r[0] > 0.0 && (r[0] < r[1] || r[1] <= 0.0))
    return (0);			/* r[0] is smallest positive root */
  if (r[1] > 0.0 && (r[1] < r[0] || r[0] <= 0.0))
    return (1);			/* r[1] is smallest positive root */
  return (-2);			/* both roots are negative */

  /* x1 should be the smallest positive root for most applications */
}



/**************************************************************************


  Synopsis:  

	ds_to_plane calculates the distance of a photon must travel to hit the plane. 
 

  Description:	

	Calcululate the distance a photon must travel to hit a plane.  The plane, just for simplicity,
	is defined as a photon structure since you need both a point and a direction to define it.  
	A plane can be defined by a position x_p and a normal lmn_p.  If the photon ray is then
	defined by x=x_v+s lmn_v and the allowed values of s are determined by the equation

	(x_v+s lmn_v - x_p)  .  lmn_p=0 where . implies the dot-product.   The routine returns
	VERY_BIG if the photon ray does not intersect the plane .


  Arguments:  

  Returns:

  Notes:

   This routine should be a more general routine than z_plane_intercept, which it replaced at
   some point.


  History:
   04aug	ksl	Changed return to +VERY_BIG if the photon cannot every hit the 
   			plane.



 ************************************************************************/

double
ds_to_plane (pl, p)
     struct plane *pl;
     struct photon *p;
{
  double denom, diff[3], numer;
  double dot ();


  if ((denom = dot (p->lmn, pl->lmn)) == 0)
    return (VERY_BIG);

  vsub (pl->x, p->x, diff);

  numer = dot (diff, pl->lmn);

  return (numer / denom);



}




/**************************************************************************


  Synopsis:  
	This routine calculates the distance a photon has to be moved to
	reach the point of closest approach to a point described by x and also
	calculates the distance of closest approach (i.e. the impact parameter).


  Description:	

  Arguments:  

  Returns:

  Notes:

  History:

 ************************************************************************/

double
ds_to_closest_approach (x, p, impact_parameter)
     double x[];		/* point for which impact parameter is calculated */
     struct photon *p;		/* Photon ptr of interest */
     double *impact_parameter;	/* distance of ray to point a closest approach */
{
  double diff[3], s, result[3];
  double length (), dot ();

  vsub (p->x, x, diff);
  s = -dot (diff, p->lmn);

  vmove (diff, p->lmn, s, result);
  *impact_parameter = length (result);



  return (s);
}
