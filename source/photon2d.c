/* Coordinate system dependences removed 04aug -- ksl 
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "atomic.h"
#include "python.h"
/***********************************************************
                                       Space Telescope Science Institute

 Synopsis:   

	translate either calls translate_in_space or translate_in_wind  depending upon the
	current location of the photon.  
  
 Arguments:		

	
 Returns:
  
Description:	
 	

Notes:
	In Monte, translate controls the flow of the photon through one grid cell.  In Python,
	there are additional possibilities since there are wind free regions.  This construction is
	used so that the calling routines for translate (trans_phot and extract) can be the same
	in the 1 and 2d versions.

	where_in_wind returns the correct domain for a given position, or domain 0 if the
	position is not in the wind of a domain.  The photon structure does not have the domain
	number directly incoded, but it can be obtained from the grid number, which where in_grid updates

History:
 	1997	ksl	Coded and debugged as part of Python effort.
 	98nov	ksl	Modified call to where_in_wind 
	05apr	ksl	55d -- Minor mod to get more info on problem case
			of a photon not in wind or grid
	11aug	ksl	70b - Incorporate mulitple components
	15aug	ksl	Incorporate multiple domains

 
**************************************************************/

int
translate (w, pp, tau_scat, tau, nres)
     WindPtr w;                 //w here refers to entire wind, not a single element
     PhotPtr pp;
     double tau_scat;
     double *tau;
     int *nres;
{
  int istat;
  int ndomain;

  if (where_in_wind (pp->x, &ndomain) < 0)
  {
    istat = translate_in_space (pp);
  }
  else if ((pp->grid = where_in_grid (ndomain, pp->x)) >= 0)
  {
    istat = translate_in_wind (w, pp, tau_scat, tau, nres);
  }
  else
  {
    istat = pp->istat = -1;     /* It's not in the wind and it's not in the grid.  Bummer! */
    Error ("translate: Found photon that was not in wind or grid, istat %i\n", where_in_wind (pp->x, &ndomain));
  }

  return (istat);
}

/***********************************************************
                                       Space Telescope Science Institute

Synopsis:   

	int translate_in_space(pp) translates the photon from its current position to the 
	edge of the wind. 

Arguments:		
	PhotPtr pp; 
	
Returns:
  
Description:	
 	

Notes:
	?? Why is it necessary to check whether the photon has hit the star at 
	this point??

History:
 	1997	ksl	Coded and debugged as part of Python effort. 
	110930	ksl	Added check for torus
 
**************************************************************/


int
translate_in_space (pp)
     PhotPtr pp;
{
  double ds, x;
  int move_phot ();

  ds = ds_to_wind (pp);


/* ?? The way in which a photon is identified as hitting the star seems
a bit convoluted.  Despite the fact that it is already identified here
as being a photon that hits the star, another determination of this 
fact is made inside walls, which uses the position of the photon to see
this. One could ask a slightly different question in walls...e.g. has the
photon hit the star in its passage from pold to the current position */

  if ((x = ds_to_sphere (geo.rstar, pp)) < ds)
  {
    x = ds;
    pp->istat = P_HIT_STAR;     /* Signifying that photon is hitting star */
  }
  move_phot (pp, ds + DFUDGE);


  return (pp->istat);
}

/***********************************************************
                                       Space Telescope Science Institute

Synopsis:   

	double ds_to_wind(pp)  calculates the photon pathlength to the edge of the wind.  
   
Arguments:		
	PhotPtr pp;.
	
Returns:

 	The distance to the nearest boundary of the wind.  
  
Description:

	Python defines the boundaries of the wind in term of the intersection
	of a biconical flow (defined by an inner and outer windcone) and 
	an inner and outer radius centered on the star.	 If the user is interested
	in a biconical flow, he/she will most likely set the radii so that they
	do not truncate the wind.  Similarly, by choosing parameters for the windcones
	properly, one can assure that one has a completely spherical flow.  However,
	the user may also perversely set the parameters to that both the conical
	flow boundaries and the min and maximum radius come into play.
	
	 Usually, one would
	define the two types of The two types of boundaries are
	usually defined so that one includes the other and so the cones apply when one is
	dealing with a CV type wind, while the inner and outer radius applies for a spherical
	model.  However this routine does not require this to be the case, since it just
	calculates where the edges are.
	
	 In any event, f you are inside the wind already ds_to_wind calculates the distance to the edge of the wind. 
	If you are outside, It will also be to the nearest edge.  	

	The routine distinguishes between  two basic cases.  Either the photon is already in the wind
	in which case one looks for the nearest boundary, or the photon is not currently
	in the wind and you have to translate it until it reaches the the boundary (or
	VERY_BIG) 
Notes:
	There is no guarantee that you will still be in the region defined by the grid.

History:
 	1997	ksl	Coded and debugged as part of Python effort. 
 	98dec	ksl	Updated to add an inner and outer radial boundary.
	99oct	ksl	The existing version of this code appears to contain
			an attempt to deal with what looks to be a roundoff
			problem when the intersection point is at large radii
			but at a radii less than VERY_BIG.  I have made the 
			fudge bigger in an attempt to push through the boundary
			and logged the error.  There ought to be a better fix
			involving not letting the photon get too far from
			the WD but I am worried about interactions between
			this routine and other portions of the code. ???

        I have been trying to understand the behavior of ds_to_wind which 
        calculates the distance to the edge of the wind (from outside the
        wind.  There were two problems I believe neither of which is solved
        to my satisfaction.  The first, simpler, problem is that at very
        large distances adding DFUDGE to ds does not "punch through" since 
        DFUDGE is so small relative to DS that it is lost in roundoff errors.
        The second problem is that ds_to_wind includes calculates the distance
        to the outer radius of the wind, but that where_in_wind says that
        the photon is still in the wind at the outer boundary.  where_in_wind
        does not have a separate return value to say that the photon is still
        in the wind_cone but beyond the outer boundary.  What to do needs 
        researching because a photon can be eclipsed even if it is beyond
        the outer boundary.  AT present ds_to_wind is set up so that it 
        it keeps calculating ds_to_wind until it punches thourh this outer
        sphere.  I have simplified the logic for this from the previous program
        but exactly what should be done needs to be determined! ????		

	05jul	ksl	Changed call from ds_to_cone to ds_to_cone as
			added cylvar coord system.  Otherwise routine is
			currently unchanged.  
	15aug	ksl	Modifications for domains.  The asumption we make
			is that the poton is not in any of the wind
			regions at this point, and that we are looking
			for the closest wind boundary.  
 
**************************************************************/



double
ds_to_wind (pp)
     PhotPtr pp;
{
  struct photon ptest;
  double ds, x;
  int ndom;

  stuff_phot (pp, &ptest);

  /* First calculated the distance to the edge of the of
     all of the "computatational domain */

  ds = ds_to_sphere (geo.rmax, &ptest);

  for (ndom = 0; ndom < geo.ndomain; ndom++)
  {
    /* Check if the photon hits the inner or outer radius of the wind */
    if ((x = ds_to_sphere (zdom[ndom].rmax, &ptest)) < ds)
      ds = x;

    if ((x = ds_to_sphere (zdom[ndom].rmin, &ptest)) < ds)
      ds = x;

    /* Check if the photon hits the inner or outer windcone */

    if ((x = ds_to_cone (&zdom[ndom].windcone[0], &ptest)) < ds)
      ds = x;
    if ((x = ds_to_cone (&zdom[ndom].windcone[1], &ptest)) < ds)
      ds = x;

    if (zdom[ndom].wind_type == CORONA)
    {

      /* As currently written ds_to_plane can give a negative number */
      x = ds_to_plane (&zdom[ndom].windplane[0], &ptest);
      if (x > 0 && x < ds)
      {
        ds = x;
      }
      x = ds_to_plane (&zdom[ndom].windplane[1], &ptest);
      if (x > 0 && x < ds)
      {
        ds = x;
      }
    }

  }

  return (ds);
}



/***********************************************************
                                       Space Telescope Science Institute

 Synopsis:   
	translate_in_wind translates the photon within a single cell in the wind.  
  
 Arguments:		

	
 Returns:
  
Description:	
 	It calculates and updates the final position of the photon p, the optical depth tau after 
	having translated the photon, and, if the photon is scattered, the number of the resonance
	responsible for scattering (or absorbing) the photon bundle.  

	These last two quantities are calculated in ds_calculate and simply passed back.  

	translate_in_wind returns that status directly to the calling routine.	

Notes:


History:
 	1997	ksl	Coded and debugged as part of Python effort. 
	04aug	ksl	Removed portions of routine that depend on
			cylindrical geometry to a separate routine
			and added logic to look for disk vertical 
			extent (for python_52).  
	05apr	ksl	55d -- Enabled spherical coordinates
	05apr	ksl	55d -- Modifed the order of the processing to address issue
			of some grid cells having very small volumes in the wind
			that are not picked up by the calculation of volumes.  Basically
			there was an option to kill these photons off entirely or
			to allow than to proceed to the next cell.  
	06nov	ksl	58b -- Modified to recognize case in which
			we want to ignore a cell because the portion
			of the volume of that cell is negligibly
			small
	09mar	ksl	Provided a capability to count reduce
			the numbers of errors for a photon
			going through a region with negligibe
			volume.  
	15aug	ksl	Incorporate multiple domains
 
**************************************************************/

/* 68c - The variable below was added because there were cases where the number
 * of photons passing through a cell with a neglible volume was becoming
 * too large and stopping the program.  This is a bandaide since if this
 * is occurring a lot we should be doing a better job at calculating the
 * volume
 */

int neglible_vol_count = 0;

int
translate_in_wind (w, p, tau_scat, tau, nres)
     WindPtr w;                 //w here refers to entire wind, not a single element
     PhotPtr p;
     double tau_scat, *tau;
     int *nres;


{

  int n;
  double smax, s, ds_current;
  int istat;
  int nplasma;
  int ndom;

  WindPtr one;
  PlasmaPtr xplasma;


/* First verify that the photon is in the grid, and if not
return and record an error */

  if ((p->grid = n = where_in_grid (wmain[p->grid].ndom, p->x)) < 0)
  {
    Error ("translate_in_wind: Photon not in grid when routine entered\n");
    return (n);                 /* Photon was not in grid */
  }

/* Assign the pointers for the cell containing the photon */

  one = &wmain[n];              /* one is the grid cell where the photon is */
  nplasma = one->nplasma;
  xplasma = &plasmamain[nplasma];
  ndom = one->ndom;




/* Calculate the maximum distance the photon can travel in the cell */

  if (zdom[ndom].coord_type == CYLIND)
  {
    smax = cylind_ds_in_cell (p);       // maximum distance the photon can travel in a cell
  }
  else if (zdom[ndom].coord_type == RTHETA)
  {
    smax = rtheta_ds_in_cell (p);
  }
  else if (zdom[ndom].coord_type == SPHERICAL)
  {
    smax = spherical_ds_in_cell (p);
  }
  else if (zdom[ndom].coord_type == CYLVAR)
  {
    smax = cylvar_ds_in_cell (p);
  }
  else
  {
    Error ("translate_in_wind: Don't know how to find ds_in_cell in this coord system %d\n", zdom[ndom].coord_type);
    exit (0);
  }

  if (one->inwind == W_PART_INWIND)
  {                             // The cell is partially in the wind
    s = ds_to_wind (p);         //smax is set to be the distance to edge of the wind
    if (s < smax)
      smax = s;
    s = ds_to_disk (p, 0);      // ds_to_disk can return a negative distance
    if (s > 0 && s < smax)
      smax = s;
  }
  if (one->inwind == W_IGNORE)
  {
    if ((neglible_vol_count % 100) == 0)
      Error
        ("translate_in_wind: Photon is in cell %d with negligible volume, moving photon %.2e  Occurrences %d\n",
         n, smax, neglible_vol_count + 1);

    neglible_vol_count++;
    move_phot (p, smax);
    return (p->istat);

  }
  else if (one->inwind == W_NOT_INWIND)
  {                             //The cell is not in the wind at all

    Error ("translate_in_wind: Grid cell %d of photon is not in wind, moving photon %.2e\n", n, smax);
    Error ("translate_in_wind: photon %d position: x %g y %g z %g\n", p->np, p->x[0], p->x[1], p->x[2]);
    move_phot (p, smax);
    return (p->istat);

  }


/* At this point we now know how far the photon can travel in it's current grid cell */

  smax += one->dfudge;          /* dfudge is to force the photon through the cell boundaries. */

/* Set limits the distance a photon can travel.  There are 
a good many photons which travel more than this distance without this 
limitation, at least in the standard 30 x 30 instantiation.  It does
make small differences in the structure of lines in some cases.  
The choice of SMAX_FRAC can affect execution time.*/

  if (smax > SMAX_FRAC * length (p->x))
  {
    smax = SMAX_FRAC * length (p->x);
  }

  /* We now determine whether scattering prevents the photon from reaching the far edge of
     the cell.  calculate_ds calculates whether there are scatterings and makes use of the 
     current position of the photon and the position of the photon at the far edge of the 
     shell.  It needs a "trial photon at the maximum distance however */


/* Note that ds_current does not alter p in any way at present 02jan ksl */

  ds_current = calculate_ds (w, p, tau_scat, tau, nres, smax, &istat);

  if (p->nres < 0)
    xplasma->nscat_es++;
  if (p->nres > 0)
    xplasma->nscat_res++;


/* OK now we increment the radiation field in the cell, translate the photon and wrap 
   things up If the photon is going to scatter in this cell, radiation also reduces 
   the weight of the photon due to continuum absorption, e.g. free free */

  if (geo.rt_mode == RT_MODE_MACRO)
  {                             // Macro-method
    /* In the macro-method, b-f and other continuum processes do not reduce the photon
    weight, but are treated as as scattering processes.  Therfore most of what was in 
    subroutine radiation can be avoided. 
    */


    one = &w[p->grid];         
    nplasma = one->nplasma;
    xplasma = &plasmamain[nplasma];
    xplasma->ntot++;


    if (geo.ioniz_or_extract == 1)       
    {
     /* For an ionization cycle */
      bf_estimators_increment (one, p, ds_current);

    /*photon weight times distance in the shell is proportional to the mean intensity */
      xplasma->j += p->w * ds_current;

    /* frequency weighted by the weights and distance in the shell.  See eqn 2 ML93 */
      xplasma->ave_freq += p->freq * p->w * ds_current;

    }

  }
  else
  {
    radiation (p, ds_current);
  }


  move_phot (p, ds_current);

  p->nres = (*nres);

  return (p->istat = istat);


}

/***********************************************************
                                       Space Telescope Science Institute

 Synopsis:   
	 walls determines whether the photon has encountered the star of disk or
	 reached the edges of the grid and returns the appropriate
	status.  
  
 Arguments:		
	PhotPtr p,pold		the current and previous description of the photon bundle.
	
 Returns:
         The photon has hit the star				P_HIT_STAR
         The photon has hit the disk				P_HIT_DISK
         The photon has reached the edge of grid 		P_ESCAPE
         The status is undeterminable                    5
         
 If a photon does not fall into one of these categories, walls returns the old status, which is stored in
 p->istat 
 
Description:	
 	

Notes:

This really does determine whether the photon is in the grid, as opposed to whether the photon
has reached a radial distance geo.rwind.  I am not sure whether this is what we want???.


History:
 	1997	ksl	Coded and debugged as part of Python effort. 
 	1997nov	ksl	Corrected problem which caused mistakes in the calculation of the disk
 	 		intercept.	 
	04aug	ksl	Added checks for a vertically extended disk
    17oct   ksl Modified the way things work to "reflect" photons from
                the disk and added reflection for stars.  Now both
                are redirected as if they were emitted from the photospher,
                however note that frequencies are not changed.

**************************************************************/
double xnorth[]={0.,0.,1.};
double xsouth[]={0.,0.,-1.};

int
walls (p, pold,normal)
     PhotPtr p, pold;
     double *normal;
{
  double r, rho, rho_sq;
  double xxx[3];
  double s, z;

  /* Check to see if the photon has hit the star. If so
   * put the photon at the star surface and use that position
   * to determine the normal to the surface, the assumption
   * being that the star is located at the center of the
   * coordiante grid.
   */
  
  if ((r = dot (p->x, p->x)) < geo.rstar_sq) {
      s=ds_to_sphere(geo.rstar,pold);
      stuff_phot(pold,p);
      move_phot(p,s);
      stuff_v(p->x,normal); 
      return (p->istat = P_HIT_STAR);
  }

  /* Check to see if it has hit the disk.  */

  if (geo.disk_type == DISK_VERTICALLY_EXTENDED)
  {
    rho = sqrt (p->x[0] * p->x[0] + p->x[1] * p->x[1]);
    if ((rho * rho) < geo.diskrad_sq && fabs (p->x[2]) <= (z = zdisk (rho)))
    {
      // We are inside a vertically extended disk.  So call ds_to_disk to find out where we hit it
      s = ds_to_disk (pold, 0);
      if (s <= 0){
          Error("walls: We seem to have ben inside the disk before now\n");
      }
      stuff_phot (pold, p);
      move_phot (p, s-DFUDGE);
      return (p->istat = P_HIT_DISK);
    }
  }
  else if (geo.disk_type == DISK_FLAT && p->x[2] * pold->x[2] < 0.0)
  {                             // Then the photon crossed the xy plane and probably hit the disk
    s = (-(pold->x[2])) / (pold->lmn[2]);
    if (s < 0)
    {
      Error ("walls: distance %g<0. Position %g %g %g \n", s, p->x[0], p->x[1], p->x[2]);
      return (-1);
    }
    // Check whether it hit the disk plane beyond the geo.diskrad**2
    vmove (pold->x, pold->lmn, s, xxx);

    if (dot (xxx, xxx) < geo.diskrad_sq )
    {                           /* The photon has hit the disk */
      stuff_phot (pold, p);     /* Move the photon to the point where it hits the disk */
      move_phot (p, s - DFUDGE); 
      if (pold->x[2]>0) {
          randvcos(p->lmn,xnorth);
          stuff_v(xnorth,normal);
          }
     else {
          randvcos(p->lmn,xsouth);
          stuff_v(xsouth,normal);
          }
     return (p->istat = P_HIT_DISK);
    }

  }

  /* At this point we know the photon has not hit the disk or the star, so we now
   * need to check if it has escaped the grid.  See note above regarding whether
   * we ought to be checking this differently.  This definition is clearly coordinate
   * system dependent.
   */

  rho_sq = (p->x[0] * p->x[0] + p->x[1] * p->x[1]);
  if (rho_sq > geo.rmax_sq)
    return (p->istat = P_ESCAPE);       /* The photon is coursing through the universe */
  if (fabs (p->x[2]) > geo.rmax)
    return (p->istat = P_ESCAPE);
  return (p->istat);            /* The photon is still in the wind */
}
