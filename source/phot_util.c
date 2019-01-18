
/***********************************************************/
/** @file  phot_util.c
 * @author ksl
 * @date   January, 2018
 *
 * @brief
 * This file contains utilities which can be used to work with photon structues.  It also
 * contains routines to calculate the distances to intersections of spheres and cones.  The
 * routines included here should not be tightly tied to a particular wind model.   The one
 * exception at present is ds_to_cone which calculates the distance to a bifurcated cone.
 *
 *
 ***********************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include "atomic.h"

#include "python.h"

int init_stuff_phot = 0;
size_t sizeofphot;


/**********************************************************/
/**
 * @brief      Simply copy one photon bundle into another
 *
 * @param [in] PhotPtr  pin   The photon bundle to be copied
 * @param [out] PhotPtr  pout   The PhotPtr that contains the output
 * @return     Always returns 0
 *
 * @details
 * This just a simple utility to copy a photon bundle
 *
 * ### Notes ###
 * All variables are (or should be) copied.
 *
 **********************************************************/

int
stuff_phot (pin, pout)
     PhotPtr pin, pout;
{
  pout->x[0] = pin->x[0];
  pout->x[1] = pin->x[1];
  pout->x[2] = pin->x[2];

  pout->lmn[0] = pin->lmn[0];
  pout->lmn[1] = pin->lmn[1];
  pout->lmn[2] = pin->lmn[2];

  pout->w = pin->w;
  pout->w_orig = pin->w_orig;
  pout->freq = pin->freq;
  pout->freq_orig = pin->freq_orig;
  pout->tau = pin->tau;

  pout->istat = pin->istat;
  pout->nres = pin->nres;
  pout->nrscat = pin->nrscat;
  pout->nscat = pin->nscat;
  pout->nnscat = pin->nnscat;
  pout->grid = pin->grid;
  pout->origin = pin->origin;
  pout->origin_orig = pin->origin_orig;
  pout->nnscat = pin->nnscat;
  pout->np = pin->np;

  pout->path = pin->path;

  return (0);
}



/**********************************************************/
/**
 * @brief      move a photon by a distance ds
 *
 * @param [in, out] PhotPtr  pp   A photon bundle
 * @param [in] double  ds   The distance to move the photon
 * @return     Always returns 0
 *
 * @details
 * The routine simply updates the position x of the
 * photon using the current direction and the distance.
 *
 * The path length for the photon is also updated
 * for use in reverberation calculations
 *
 * ### Notes ###
 *
 **********************************************************/

int
move_phot (pp, ds)
     PhotPtr pp;
     double ds;
{

  pp->x[0] += pp->lmn[0] * ds;
  pp->x[1] += pp->lmn[1] * ds;
  pp->x[2] += pp->lmn[2] * ds;
  pp->path += fabs (ds);
  return (0);
}



/**********************************************************/
/**
 * @brief      Compares two photons to determine whether they have
 *  	the identical position and direction.
 *
 * @param [in] PhotPtr  p1   The first photon bundle
 * @param [in] PhotPtr  p2   The second phton bundle
 * @return     Returns 0 if the position and direction are identical,
 * 1 otherwise
 *
 * @details
 *
 * ### Notes ###
 * The routine is used in calculate_ds
 *
 **********************************************************/

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


//OLD /**************************************************************************
//OLD
//OLD
//OLD   Synopsis:
//OLD
//OLD   phot_hist is designed to record the history of a single photon
//OLD   bundle, in terms of a series of photon structures that are
//OLD   populated as the photon goes through the grid
//OLD
//OLD   Unless phot_hist_on is true, then this routine is a NOP
//OLD
//OLD
//OLD   Description:
//OLD
//OLD   This is a diagnostic routine only
//OLD
//OLD   Arguments:
//OLD
//OLD   Returns:
//OLD
//OLD   Notes:
//OLD
//OLD   History:
//OLD   090125  ksl     68b -Created to better understand where photon bundles were
//OLD                   being absorbed or losing energy in the grid
//OLD
//OLD  ************************************************************************/

/**********************************************************/
/**
 * @brief      Record the history of a single photon
 *   bundle, in terms of a series of photon structures that are
 *   populated as the photon goes through the grid
 *
 *   Unless phot_hist_on is true, then this routine is a NOP
 *
 * @param [in, out] PhotPtr  p   A photon to track
 * @param [in, out] int  iswitch   A switch that resets
 * @return     Returns the number of photons which have been stored
 *
 * @details
 *
 * This routine stores selected photons in an array, which one
 * can use to see what is happening in a special case.
 *
 * This is a diagnostic routine only
 *
 * ### Notes ###
 * The routine is useful (at present) only in the context of a debugger
 * since the array in which photons are stored is not written out.  The
 * basic idea is that one identifies a problem say with a specific photon
 * and then stores the history of that photon in the array so once can
 * inspect in detail what has  happened to the photon.
 *
 * @bug  It is not obvious that this routine is needed as we now have
 * another mechanism to write photons out to a file.  In any event,
 * this should be moved to some other place in the code.
 *
 **********************************************************/

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
    Error ("phot_hist: The number of steps %d in phot_hist exceeds MAX_PHOT_HIST\n", MAX_PHOT_HIST);
  }

  return (n_phot_hist);
}


//OLD /**************************************************************************
//OLD
//OLD
//OLD   Synopsis:
//OLD   The next routine is destigned to update a portion of the PlasmaPtr to reflect where
//OLD   photons along the line of sight to the observer were absorbed in the wind
//OLD
//OLD   Description:
//OLD
//OLD   Arguments:
//OLD
//OLD   Returns:
//OLD
//OLD   Notes:
//OLD   As photon is extracted from the wind, tau changes due to scatters and w changes due
//OLD   to absorption.  We are just recording how much energy is absorbed by scatterng processes
//OLD   here, and so the energy absorbed is the current weight (exp(-tau_before) - exp (-tau_after))
//OLD
//OLD   History:
//OLD   090211  ksl     Created to store the energy removed photons headed toward the observer
//OLD                   by an ion in a particular cell of the wind.
//OLD   0904    ksl     68c - Fixed problem concerning where energy was being stored
//OLD
//OLD  ************************************************************************/


/**********************************************************/
/**
 * @brief      The next routine is designed to update a portion of the PlasmaPtr to reflect where
 *  	photons along the line of sight to the observer were absorbed in the wind
 *
 * <NOTE: The [in, out] tag describes if the value of a parameter is used or altered. If it is used but not altered, delete 'OUT'. If the original value is not used and it is written to, delete 'IN'.>
 * @return     ??? RETURNS ???
 *
 * @details
 * ??? DESCRIPTION ???
 *
 * ### Notes ###
 * As photon is extracted from the wind, tau changes due to scatters and w changes due
 * 	to absorption.  We are just recording how much energy is absorbed by scatterng processes
 * 	here, and so the energy absorbed is the current weight (exp(-tau_before) - exp (-tau_after))
 *
 **********************************************************/

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

    tau = p->tau;               // tau is tau after the scatter

    nion = lin_ptr[p->nres]->nion;      // ion that scattered

    x = p->w * (exp (-tau_old) - exp (-tau));   // energy removed by scatter

    xplasma = &plasmamain[wmain[p->grid].nplasma];      // pointer to plasma cell where scattering occured

    xplasma->xscatters[nion] += (x);

    tau_old = tau;

    if (xplasma->xscatters[nion] < 0.0)
    {
      Error ("phot_history_summarize:  n %d n_phot_hist %d phot_hist %d nplasma %d\n", n, n_phot_hist, p->grid, wmain[p->grid].nplasma);
    }
  }



  return (0);
}



//OLD /*******************************************************
//OLD             Space Telescope Science Institute
//OLD
//OLD Synopsis:
//OLD    ds_to_cone (cc, p) This is a new version of ds_to_cone introduced in order to
//OLD    allow cylvar coordinates  Eventually it is intended to replace ds_to_windcone
//OLD    completely, or at least the guts of ds_to_windcone.
//OLD
//OLD Returns:
//OLD
//OLD
//OLD Description:
//OLD
//OLD
//OLD Notes:
//OLD
//OLD As part of the conversion, the old ds_to_cone was renamed to ds_to_windcone
//OLD
//OLD
//OLD
//OLD   56d -- Cones have become more prominent in python with time.  Originally
//OLD   they were used simply to define the inner and outer edges of the
//OLD   wind, that is   to determine when a photon entered the wind.  The wind was/is
//OLD   assumed to be a biconical flow, and therefore the most naturally was defined in
//OLD   terms of a footpoint in the disk for the innermost/outermost stream line
//OLD   and a slope.
//OLD
//OLD   It was expanced when rtheta coordanates were introduced
//OLD   to handle find when one had hit the theta coordinate boundaries.  And in
//OLD   56d we also want it to handle the case where the grid cells have variable
//OLD   z coordinates.  Therefore one needs to be quite careful, Since ds_to_cone
//OLD   is the ultimate name of the routine I would like to use, the first step
//OLD   I took was to change the name of the old routine to ds_to_windcone.  This
//OLD   should allow me to replace the old routine incrementally.
//OLD
//OLD   56d -- With python 56d the variable defining a cone were changed to be
//OLD   the point on the z axis intercepted by the cone (extended) through
//OLD   the disk, and the slope.
//OLD
//OLD   Our approach is as follows:
//OLD
//OLD   The cone is defined by  z=z_o + dz/drho *drho
//OLD   The photon is defined by X=X_o + LMN s
//OLD
//OLD   Unless dz/drho = 0, this turns into a simple quadratic equation that must
//OLD   be solved.  We work in the northen hemisphere only.
//OLD
//OLD   There are issues concerning what to do so one crosses the disk plane.  For
//OLD         finding how far a photon can travel in a cell, we are *not* interested
//OLD   in the possibility that the photon hits the cone on the other side of
//OLD   the disk plane, but for a photon travelling out of the wind we are
//OLD   very interested in this possibility.
//OLD
//OLD History:
//OLD   05jul   ksl     Created so that cylvar coordinates could be incorporated
//OLD                   into python.
//OLD   06sep   ksl     57h -- Corrected an error discovered by SS in the setup of the
//OLD                   quadratic for a cone.  Note that if this routing ever became
//OLD                   a "heavy hitter" in terms of efficiency, there are a few
//OLD                   changes that could speed it up a bit since the magnitude
//OLD                   of lmn is always 1.
//OLD */



/**********************************************************/
/**
 * @brief      (cc, p) This is a new version of  introduced in order to
 *    allow cylvar coordinates  Eventually it is intended to replace ds_to_windcone
 *    completely, or at least the guts of ds_to_windcone.
 *
 * <NOTE: The [in, out] tag describes if the value of a parameter is used or altered. If it is used but not altered, delete 'OUT'. If the original value is not used and it is written to, delete 'IN'.>
 * @param [in, out] ConePtr  cc   ???
 * @param [in, out] struct photon *  p   ???
 * @return     ??? RETURNS ???
 *
 * @details
 * ??? DESCRIPTION ???
 *
 * ### Notes ###
 * As part of the conversion, the old ds_to_cone was renamed to ds_to_windcone
 *
 *
 *
 * 	56d -- Cones have become more prominent in python with time.  Originally
 * 	they were used simply to define the inner and outer edges of the
 * 	wind, that is   to determine when a photon entered the wind.  The wind was/is
 * 	assumed to be a biconical flow, and therefore the most naturally was defined in
 * 	terms of a footpoint in the disk for the innermost/outermost stream line
 * 	and a slope.
 *
 * 	It was expanced when rtheta coordanates were introduced
 * 	to handle find when one had hit the theta coordinate boundaries.  And in
 * 	56d we also want it to handle the case where the grid cells have variable
 * 	z coordinates.  Therefore one needs to be quite careful, Since ds_to_cone
 * 	is the ultimate name of the routine I would like to use, the first step
 * 	I took was to change the name of the old routine to ds_to_windcone.  This
 * 	should allow me to replace the old routine incrementally.
 *
 * 	56d -- With python 56d the variable defining a cone were changed to be
 * 	the point on the z axis intercepted by the cone (extended) through
 * 	the disk, and the slope.
 *
 * 	Our approach is as follows:
 *
 * 	The cone is defined by  z=z_o + dz/drho *drho
 * 	The photon is defined by X=X_o + LMN s
 *
 * 	Unless dz/drho = 0, this turns into a simple quadratic equation that must
 * 	be solved.  We work in the northen hemisphere only.
 *
 * 	There are issues concerning what to do so one crosses the disk plane.  For
 *         finding how far a photon can travel in a cell, we are *not* interested
 * 	in the possibility that the photon hits the cone on the other side of
 * 	the disk plane, but for a photon travelling out of the wind we are
 * 	very interested in this possibility.
 *
 **********************************************************/

double
ds_to_cone (cc, p)
     ConePtr cc;
     struct photon *p;
{
  double dz, dzdr2;
  double a, b, c, root[2];
  double s_to_zero;             /* The path length to the xy plane */
  int i;
  struct photon pp;

  /* First of all let's work only in the "northern" hemisphere */
  stuff_phot (p, &pp);

  if (pp.x[2] < 0.0)
  {                             /*move the photon to the northen hemisphere */
    pp.x[2] = -pp.x[2];
    pp.lmn[2] = -pp.lmn[2];
  }

  /* Set up and solve the quadratic equation that gives the cone intercept */

  dzdr2 = cc->dzdr * cc->dzdr;
  dz = pp.x[2] - cc->z;

  a = dzdr2 * (pp.lmn[0] * pp.lmn[0] + pp.lmn[1] * pp.lmn[1]) - (pp.lmn[2] * pp.lmn[2]);
  b = 2. * (dzdr2 * (pp.lmn[0] * pp.x[0] + pp.lmn[1] * pp.x[1]) - pp.lmn[2] * dz);
  c = dzdr2 * (pp.x[0] * pp.x[0] + pp.x[1] * pp.x[1]) - dz * dz;

  i = quadratic (a, b, c, root);        /* root[i] is the smallest positive root unless i is
                                           negative in which case either both roots were negative or both roots were imaginary */

  /* Calculate the positive path length to the xy plane.  If either the
     photon is travelling in the xy plane or if the intercept is in the negative
     direction set the path length to infinity */

  if (pp.x[2] * pp.lmn[2] >= 0)
    s_to_zero = VERY_BIG;
  else
    s_to_zero = (-pp.x[2] / pp.lmn[2]);


  if (i >= 0 && root[i] < s_to_zero)
    return (root[i]);           /*Because that implies
                                   the ray hit the cone before hitting the xy plane */
  return (s_to_zero);


}


//OLD /**************************************************************************
//OLD
//OLD
//OLD   Synopsis:
//OLD
//OLD   Calculate the pathlenth along a line of sight defined by
//OLD   a photon p to a sphere or radius r centered on the origin.  If
//OLD   the photon does not hit the sphere return a large number VERY_BIG
//OLD
//OLD   Description:
//OLD
//OLD   Arguments:
//OLD
//OLD   Returns:
//OLD
//OLD   Notes:
//OLD
//OLD   History:
//OLD
//OLD  ************************************************************************/

/**********************************************************/
/**
 * @brief      Calculate the pathlenth along a line of sight defined by
 * 	a photon p to a sphere centered on the origin.
 *
 * @param [in, out] double  r   The radius of the sphere
 * @param [in, out] struct photon *  p   The photon
 * @return     The distance to the sphere.
 * 	If
 * 	the photon does not hit the sphere return a large number VERY_BIG
 *
 * @details
 * The routine simply solves the quadratic equation for calculating
 * the distance to the sphere.  It returns the smalles positive root
 *
 * ### Notes ###
 *
 *
 **********************************************************/

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
    return (root[i]);           /*Because that implies
                                   the ray hit the sphere */

  return (VERY_BIG);
}




//OLD /**************************************************************************
//OLD
//OLD
//OLD   Synopsis:
//OLD   This is more generalized routine to find the positive distance to
//OLD           a sphere centered at x with radius r
//OLD
//OLD   Description:
//OLD
//OLD   Arguments:
//OLD
//OLD   Returns:
//OLD
//OLD   Notes:
//OLD
//OLD   History:
//OLD
//OLD  ************************************************************************/


/**********************************************************/
/**
 * @brief      Find the positive distance to
 *    	a sphere centered at x with radius r
 *
 * @param [in, out] double  x[]   The center of the sphere
 * @param [in, out] double  r   The radius of the sphere
 * @param [in, out] struct photon *  p   A photon
 * @return     The distance to the sphere.
 * 	If
 * 	the photon does not hit the sphere return a large number VERY_BIG
 *
 *
 * @details
 * The routine simply solves the quadratic equation for calculating
 * the distance to the sphere.  It returns the smalles positive root
 *
 *
 * ### Notes ###
 * It is not clear we need two routines like this, one which
 * assumes that the sphere is centered on 0 and one which does
 * not.
 *
 **********************************************************/

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

  i = quadratic (a, b, c, root);        /* root[i] is the smallest positive root unless i is
                                           negative in which case either both roots were negative or both roots were imaginary */


  if (i >= 0)
    return (root[i]);           /*Because that implies
                                   the ray hit the sphere */

  return (VERY_BIG);
}



//OLD /**************************************************************************
//OLD
//OLD
//OLD   Synopsis:
//OLD
//OLD   This solves a simple quadratic (or if a is zero linear equation).  The return is set up
//OLD    to make it easy to identify the smallest positive root if one exists.  The routine returns
//OLD    a negative number if both roots are negative or imaginary.
//OLD    More specifically
//OLD    -1 -> both roots are imaginary
//OLD          -2 -> both roots are negative or 0
//OLD           0 -> the first root is the smallest positive  root
//OLD           1 -> the second root is the smallest positive root
//OLD
//OLD
//OLD   Description:
//OLD
//OLD   Arguments:
//OLD
//OLD   Returns:
//OLD
//OLD   Notes:
//OLD
//OLD   History:
//OLD   05jul   ksl     56d -- Heretofore roots of 0 were on no insterest, but now that
//OLD                   quadratic is called by bilin, we would like to know about values
//OLD                   of zero, specifically.  Since the simplest thing to do was to
//OLD                   check both roots in this case, I added a little code to make
//OLD                   sure root was defined in all cases.
//OLD
//OLD  ************************************************************************/

/**********************************************************/
/**
 * @brief      This solves a simple  quadratic (or if a is zero linear) equation.
 *
 * @param [in] double  a   The constant which multiples x**2 in the quadratic equation
 * @param [in] double  b   The constant which multiples x
 * @param [in] double  c   The constant
 * @param [out] double  r[]   The roots of the quadratic equation
 * @return   The returns depend on the a status tha depends on the nature of the solutions
 * The return is set up
 *    to make it easy to identify the smallest positive root if one exists.  The routine returns
 *    a negative number if both roots are negative or imaginary.  More specifically the routine returns
 *  *  -1 -> both roots are imaginary
 *  * -2 -> both roots are negative or 0
 *  * 0 -> the first root is the smallest positive  root
 *  * 1 -> the second root is the smallest positive root
 *
 * @details
 * The equation that is solved is ax**2+bx+c=0.
 *
 * ### Notes ###
 *
 **********************************************************/

int
quadratic (a, b, c, r)
     double a, b, c, r[];
{
  double q, z;

  if (a == 0.0)
  {                             /* Then it's not really a quadratic but we can solve it
                                   anyway */
    if (b == 0.0)
    {
      r[0] = r[1] = -99.;
      return (-1);              /* The roots are extremely imaginary, since both a a b were 0 */
    }

    r[0] = r[1] = (-c / b);     // Then it was a linear equation. Setting both roots to the same thing could be a problem ksl
    if (r[0] < 0.0)
      return (-2);              /* Generally speaking we are not interested in
                                   negative distances */
    else
      return (0);
  }

  if ((q = b * b - 4. * a * c) < 0.0)
  {
    r[0] = r[1] = -99.;
    return (-1);                /* both roots are imaginary */
  }
  q = sqrt (q);
  z = 0.5 / a;

  r[0] = (-b - q) * z;
  r[1] = (-b + q) * z;


  if (r[0] > 0.0 && (r[0] < r[1] || r[1] <= 0.0))
    return (0);                 /* r[0] is smallest positive root */
  if (r[1] > 0.0 && (r[1] < r[0] || r[0] <= 0.0))
    return (1);                 /* r[1] is smallest positive root */
  return (-2);                  /* both roots are negative */

  /* x1 should be the smallest positive root for most applications */
}



//OLD /**************************************************************************
//OLD
//OLD
//OLD   Synopsis:
//OLD
//OLD   ds_to_plane calculates the distance of a photon must travel to hit the plane.
//OLD
//OLD
//OLD   Description:
//OLD
//OLD   Calcululate the distance a photon must travel to hit a plane.  The plane, just for simplicity,
//OLD   is defined as a photon structure since you need both a point and a direction to define it.
//OLD   A plane can be defined by a position x_p and a normal lmn_p.  If the photon ray is then
//OLD   defined by x=x_v+s lmn_v and the allowed values of s are determined by the equation
//OLD
//OLD   (x_v+s lmn_v - x_p)  .  lmn_p=0 where . implies the dot-product.   The routine returns
//OLD   VERY_BIG if the photon ray does not intersect the plane .
//OLD
//OLD
//OLD   Arguments:
//OLD
//OLD   Returns:
//OLD
//OLD   Notes:
//OLD
//OLD    This routine should be a more general routine than z_plane_intercept, which it replaced at
//OLD    some point.
//OLD
//OLD
//OLD   History:
//OLD    04aug  ksl     Changed return to +VERY_BIG if the photon cannot every hit the
//OLD                           plane.
//OLD
//OLD
//OLD
//OLD  ************************************************************************/


/**********************************************************/
/**
 * @brief      calculates the distance of a photon must travel to hit the plane.
 *
 * @param [in, out] struct plane *  pl   A structure which descibes a plane in terms of a position and vector normal to the surface
 * @param [in, out] struct photon *  p   A photon
 * @return     The distance to the plane, or VERY_BIG if the photon does not hit the plane (in the positive direction)
 *
 * @details
 * Calcululate the distance a photon must travel to hit a plane.  The plane, just for simplicity,
 * is defined as a photon structure since you need both a point and a direction to define it.
 * A plane can be defined by a position x_p and a normal lmn_p.  If the photon ray is then
 * defined by x=x_v+s lmn_v and the allowed values of s are determined by the equation
 *
 * (x_v+s lmn_v - x_p)  .  lmn_p=0 where . implies the dot-product.   The routine returns
 * VERY_BIG if the photon ray does not intersect the plane .
 *
 * ### Notes ###
 *
 **********************************************************/

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




//OLD /**************************************************************************
//OLD
//OLD
//OLD   Synopsis:
//OLD   This routine calculates the distance a photon has to be moved to
//OLD   reach the point of closest approach to a point described by x and also
//OLD   calculates the distance of closest approach (i.e. the impact parameter).
//OLD
//OLD
//OLD   Description:
//OLD
//OLD   Arguments:
//OLD
//OLD   Returns:
//OLD
//OLD   Notes:
//OLD
//OLD   History:
//OLD
//OLD  ************************************************************************/


/**********************************************************/
/**
 * @brief      This routine calculates the distance a photon has to be moved to
 * 	reach the point of closest approach to a point described by x and also
 * 	calculates the distance of closest approach (i.e. the impact parameter).
 *
 * @param [in] double  x[]   A position in space
 * @param [in] struct photon *  p   A photon bundle
 * @param [out] double *  impact_parameter   The distance of closest approach of the photon
 * @return     The distance the photon would have to be moved to reach the point of closest approach,
 * which can be negative
 *
 * @details
 *
 * ### Notes ###
 *
 **********************************************************/

double
ds_to_closest_approach (x, p, impact_parameter)
     double x[];                /* point for which impact parameter is calculated */
     struct photon *p;          /* Photon ptr of interest */
     double *impact_parameter;  /* distance of ray to point a closest approach */
{
  double diff[3], s, result[3];
  double length (), dot ();

  vsub (p->x, x, diff);
  s = -dot (diff, p->lmn);

  vmove (diff, p->lmn, s, result);
  *impact_parameter = length (result);

  return (s);
}






//OLD /**************************************************************************
//OLD                     Space Telescope Science Institute
//OLD
//OLD
//OLD   Synopsis:
//OLD
//OLD   calculate the distance a photon must travel to encouter a cylinder
//OLD   that is centered on the z axis.
//OLD
//OLD   Description:
//OLD
//OLD   Given a position and a direction the distance to a cylinder can
//OLD   be solved as a simple quadratic equation
//OLD
//OLD   Arguments:
//OLD
//OLD   rho is the size of the cylinder
//OLD   p is a photon ptr
//OLD
//OLD
//OLD   Returns:
//OLD
//OLD   The smallest positive distance to the cylinder. If the
//OLD   ray does not hit the cylinder a large postive number is
//OLD   returned.
//OLD
//OLD   Notes:
//OLD
//OLD
//OLD
//OLD   History:
//OLD 11aug     ksl     Coding began
//OLD
//OLD  ************************************************************************/




/**********************************************************/
/**
 * @brief      calculate the distance a photon must travel to encouter a cylinder
 *   that is centered on the z axis.
 *
 * @param [in, out] double  rho   the size of the cylinder
 * @param [in, out] struct photon *  p   a photon ptr
 * @return     The smallest positive distance to the cylinder. If the
 *   ray does not hit the cylinder a VERY_BIG is
 *   returned.
 *
 * @details
 * Given a position and a direction the distance to a cylinder,
 * solve the simple quadratic equation that gives the distance
 * to the cylinder.
 *
 * ### Notes ###
 *
 **********************************************************/

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
