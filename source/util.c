/***********************************************************
                                       Space Telescope Science Institute

Synopsis:  These are very simple utilities within python that do not
	fall into a larger catergory like those in vector.c


Description:	 

int
locate (x, xarray, xdim, nelem, frac)

double 
linterp(x,xarray,yarray,xdim)


Notes:


History:
 	02feb  ksl      Coded as part of python effort
	02mar	ksl	Eliminated ndim in favor of xdim so 
			that ndim could be used as a global 
			variable

**************************************************************/



#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "log.h"
#include "atomic.h"
#include "python.h"

#define EPS 1.e-10



/* 
This routine is intended to locate the array element that is just below the
value of x and give one the appropriate "fraction" for interpolation between
array elements.  The idea behind this is that one wants to predict the velocity
at a position x and one has a structure tha contains the velocities at specific
values of x

The array "array" must be monotonically increasing and have dimension
xdim (with array elements 0 to xdim filled.  

The program returns 0 if the desired x is contained in the array, -1
if x is below the value of the first element and -2 if the element is above
the highest element.  

the array element in question and the fractional position 
for interpolation 
*/




/* fraction is a utility to speed up the search for the bracketing
topbase photionization x-sections, but in principle it should work in
other situations within python 

Returns: 0 is the value asked for was within the bounds of the array,
	-1 if below the lower bound, +1 if above the lower bound

Notes:  The routine is similar (I found aut afterward) to the numerical
	recipes routine locate.  It is possible there should be a separate
	routine locate to avoid the division at the bottom.

	This approach may not be ideal for many applications, especially
	those in which we have a good idea of where in the array
	the new value is located.

History:
	01dec	ksl	Added limits. fraction now produces an error
			if the value is out of range.  It also now
			assumes that it should return the limits when 
			a value is out of range.
	01dec	ksl	Changed calls so that fraction returns 0 if
			OK, -1 if the point is below the limit
			1 if it is above the maximum.  Since in many
			cases you don't care that you are asking for
			things out of range, I've suppressed the errors
			after the first 10.
	02feb	ksl	Modified slightly.  Note that this routine is
			a substantial fraction of time it takes to
			run python and it needs to be absolutely as
			efficient as possible.
*/



int
fraction (value, array, npts, ival, f, mode)
     double array[];		// The array in we want to search
     int npts, *ival;		// ival is the lower point
     double value;		// The value we want to index
     double *f;			// The fractional "distance" to the next point in the array
     int mode;        // 0 = compute in lin space, 1=compute in log space
{
  int imin, imax, ihalf;

  if (value < array[0])
    {
      *ival = 0;
      *f = 0.0;
      return (-1);
    }

  imax = npts - 1;
  if (value > array[imax])
    {
      *ival = npts - 2;
      *f = 1.0;
      return (1);
    }



  imin = 0;

/* In what follows, there is a specific issue as to how to treat
the situation where the value is exactly on an array element. Right
now this is set to try to identify the array element below the one
on which the value sits, and to set the fraction to 1.  This was
to reflect the behavior of the search routine in where_in_grid. */

  while (imax - imin > 1)
    {
      ihalf = (imin + imax) >> 1;	// Compute a midpoint >> is a bitwise right shift
      if (value > array[ihalf])
	{
	  imin = ihalf;
	}
      else
	imax = ihalf;
    }

// So array[imin] just <= value
  if (mode==0)
  	*f = (value - array[imin]) / (array[imax] - array[imin]);        //linear interpolation 
  else if (mode==1)
  	*f = (log(value) - log(array[imin])) / (log(array[imax]) - log(array[imin]));     //log interpolation
  else
	{
	Error("Fraction - unknown mode %i\n",mode);
	exit(0);
	}

  *ival = imin;

  return (0);
}





/* 
Given a number x, and an array of x's in xarray, and functional
values y = f(x) in yarray and the dimension of xarray and yarray,
linterp calculates y, the linearly interpolated value of y(x). It
also returns the element in the xarray so that one can consider
not doing a search to find the right array element in certain 
circumstances

	02feb	ksl	Changed call so error could be returned
	02jul	ksl	Actually included the error in the return.
	02jul	ksl	Modified again so now it returns nelem.
*/

int
linterp (x, xarray, yarray, xdim, y, mode)
     double x;			// The value that we wish to index i
     double xarray[], yarray[];
     int xdim;
     double *y;
     int mode;      //0 = linear, 1 = log
{
  int nelem;
  double frac;
  int ierr;
  ierr = fraction (x, xarray, xdim, &nelem, &frac, mode);


  if (mode==0)
  	*y = (1. - frac) * yarray[nelem] + frac * yarray[nelem + 1];
  else if (mode==1)
	*y = exp((1. - frac) * log(yarray[nelem]) + frac * log(yarray[nelem + 1]));
  else
	{
	Error("linterp - unknown mode %i\n",mode);
	exit(0);
	}

  return (nelem);

}




/***********************************************************
           Space Telescope Science Institute

Synopsis:  coord_fraction() is used to calculate the fractional 
	contributions in interpolations of values of quantities,
	such as electron density, by various cells in the 
	grid in a coordinate system independent way. 



Description:	 
	ichoice=0 --> interpolate on vertices
	ichoice=1 --> interpolate on centers

	x --> the 3-vector position for which you want
		the fractinal position
Returns:
	ii[] --> an array that contains the 1-d element
		numbers that must be summed  
	frac[] --> the array that contains the fractional
		contribution of the corresponding 
		element
	nelem --> the number of elements that must be
		summed, nominally 2 for a spherical grid
		and 4 for a two dimensional grid.  (For
		a 3 d grid it would be 6, but we have not
		implemented this.

	 1 if  inside the grid
	-2 if outside the grid
	-1 if inside the grid

Notes:
	There are numerous times when one wants the value
	of an interpoalted  variable in the wind.  There
	is no easy way to interpolate the variable easily.
	What this routine does is calculate the fractional
	contributions of elements in the array to that
	position.  Then one must sum up the actual variable
	else where

	If positions are outside the grid, coord_fraction
	attempts to give you the value at the edge of teh
	grid.

	It's possible that coord_fraction could be used 
	to interpolate beyond the edge of the grid where
	a variable is defined, although this is not done
	at present!

History:
	04aug	ksl	52a -- Coded as a generic routine
			to get the location in the grid
			and the fractions that can be
			used for interpolation
	04dec	ksl	54a -- Minor mod to exit if don't
			understand coord type
	05apr	ksl	55d -- Made a significant change
			to the way coord_fraction works
			to concentrate variations due to different
			coordinate systems here. Note that the
			variables for coord_fraction were changed
			significantly
	05jul	ksl	56d -- Modified so that for cylvar coords
			one simple calls the a special routine
			which is in cylindvar.c.  PROBABLY THIS ROUTINE
			SHOULD BE REWRITTEN SO THIS IS THE CASE
			FOR ALL COORDINATE SYSTEMS.
	13sep	nsh	76b -- Modified calls to fraction to take
			account of new mode.
		       	

**************************************************************/

int ierr_coord_fraction = 0;

int
coord_fraction (ichoice, x, ii, frac, nelem)
     int ichoice;
     double x[];
     int ii[];
     double frac[];
     int *nelem;
{
  double r, z;
  double *xx, *zz;
  int ix, iz;
  double dr, dz;
  int cylvar_coord_fraction ();
  int n;


  /* Jump to special routine if CYLVAR coords */

  if (geo.coord_type == CYLVAR)
    {

      n = cylvar_coord_fraction (ichoice, x, ii, frac, nelem);
      if (n < 0 && ierr_coord_fraction < 1000)
	{
	  Error
	    ("coord_fraction: cylvar_coord fraction returning %d not in grid\n",
	     n);
	  ierr_coord_fraction++;
	}
      return (n);

    }

  /* Assign pointers to the xx and zz depending on whether
   * one wants to interpolate on vertex points (0) or 
   * midpoints (1)
   */

  if (ichoice == 0)
    {
      xx = wind_x;
      zz = wind_z;
    }
  else
    {
      xx = wind_midx;
      zz = wind_midz;
    }

  /* Now convert x to the appropriate coordinate system */
  if (geo.coord_type == CYLIND)
    {
      r = sqrt (x[0] * x[0] + x[1] * x[1]);
      z = fabs (x[2]);
    }
  else if (geo.coord_type == RTHETA)
    {
      r = length (x);
      z = acos (fabs (x[2]) / r) * RADIAN;
    }
  else if (geo.coord_type == SPHERICAL)
    {
      r = length (x);
      z = 0;			// To avoid -O3 warning
    }
  else
    {
      Error
	("coord_fraction: Don't know how to handle this coordinate type %d\n",
	 geo.coord_type);
      exit (0);
    }

  if (geo.coord_type == SPHERICAL)
    {				/* We are dealing with a 1d system */
      fraction (r, xx, NDIM, &ix, &dr, 0); //linear space
      ii[0] = ix;
      frac[0] = (1. - dr);
      ii[1] = ix + 1;
      frac[1] = dr;
      *nelem = 2;
      if (sane_check (dr))
	{
	  Error ("coord_frac:sane_check dr=%f for spherical coords. \n", dr);
	}
    }
  else
    {				/* We are dealing with a 2d system */
      fraction (r, xx, NDIM, &ix, &dr, 0);
      fraction (z, zz, MDIM, &iz, &dz, 0);

      ii[0] = ix * MDIM + iz;
      frac[0] = (1. - dz) * (1. - dr);

      ii[1] = (ix + 1) * MDIM + iz;
      frac[1] = (1. - dz) * dr;

      ii[2] = ix * MDIM + iz + 1;
      frac[2] = (dz) * (1. - dr);

      ii[3] = (ix + 1) * MDIM + iz + 1;
      frac[3] = (dz) * (dr);
      *nelem = 4;

      if (sane_check (dr) || sane_check (dz))
	{
	  Error ("coord_frac:sane_check dr=%f dz=%f for 2d coords\n", dr, dz);
	}
    }
  /* At this point i,j are just outside the x position */
  /* Check to see if x is outside the region of the calculation */
  /* Note that this is a very incoplethe check in the sneste that 
   * the posision could be out of the grid in other directions */

  if (r > xx[NDIM - 1])
    {
      return (-2);		/* x is outside grid */
    }
  else if (r < xx[0])
    {
      return (-1);		/*x is inside grid */
    }

  return (1);

}


/***********************************************************
           Space Telescope Science Institute

Synopsis:  where_in_2dcell calculates the fractional
	location of a position in any 2d grid.  

Description:	 
	ichoice=0 --> interpolate on vertices
	ichoice=1 --> interpolate on centers

	x --> the 3-vector position for which you want
		the fractinal position
Returns:
	fracpos[] --> the array that contains the fractional
		position 

	 0 (FALSE) if the position is in the cell
	!0 (TRUE)  if the position is not in the cell

Notes:
	This routine is general.  It does not rely on the
	predefined arrays that just contain the grid vertices
	and centers. On the other hand, it does not calculate
	which grid cell x lies in. It only tells you that x
	is not in a particular cell if that is the case.

	A cell is organized as follows:

		x01	x11
		x00	x10
	
History:
	05jul	ksl	56d -- Created in the context of incorporating
			CYLVAR coordinates.
		       	

**************************************************************/

int ierr_where_in_2dcell = 0;

int
where_in_2dcell (ichoice, x, n, fx, fz)
     int ichoice;
     double x[];
     int n;
     double *fx, *fz;
{
  double *x00, *x01, *x10, *x11;
  double z[3];
  int i;
  int bilin ();

  /* Assign the corners ot the region we want to determin
   * the fractional position of
   */

  if (n < 0 || n + MDIM + 1 >= NDIM2)
    {
      if (ierr_where_in_2dcell < 100)
	{
	  Error ("where_in_2dcell: Unreasonable n %d \n", n);
	  ierr_where_in_2dcell++;
	}
      return (n);		// And hope that calling routine knows how to handle this.
    }

  if (ichoice == 0)
    {
      x00 = wmain[n].x;
      x01 = wmain[n + 1].x;
      x10 = wmain[n + MDIM].x;
      x11 = wmain[n + MDIM + 1].x;
    }
  else
    {
      x00 = wmain[n].xcen;
      x01 = wmain[n + 1].xcen;
      x10 = wmain[n + MDIM].xcen;
      x11 = wmain[n + MDIM + 1].xcen;
    }

  /* Rotate the input vector onto the xz plane
   * which is the plane where the grid is defined
   */

  z[0] = sqrt (x[0] * x[0] + x[1] * x[1]);
  z[1] = 0.0;
  z[2] = fabs (x[2]);

  i = bilin (z, x00, x01, x10, x11, fx, fz);

  return (i);


}





/***********************************************************
          Space Telescope Science Institute
Synopsis: wind_n_to_ij(n,i,j) and wind_ij_to_n(i,j,n) are two 
routines which effectively translate between the one dimensional 
wind structure and the two dimensional physical grid 

Arguments: 

the element in the wind structure i,j the position in the 
two dimensional grid

Returns:

Description:

Notes:
For 2d, the underlying 1-d grid is organized so that as n increases, 
one goes up in z, and then steps out in rho., or in theta and then
z as the case may be.

In spherical coordinates, we assume that 
 
History:
	97jan	ksl	Coding on python began.
	02feb	ksl	Allowed for different dimensions in x and z
	05apr	ksl	Added error_check to verify routine 
			is not called with spherical coordiantes
**************************************************************/


int
wind_n_to_ij (n, i, j)
     int n, *i, *j;
{
  if (geo.coord_type == SPHERICAL)
    {
      *i = n;
      *j = 0;
/*  130624 - ksl - Removed this error as it program appears to be working as intended
 *  and there are times when this routine should be called for spherical models
 *  since we basically still have a 2d grid in this case      
       Error
	("Warning: wind_n_to_ij being called for spherical coordinates %d \n",
	 n);
*/
    }
  *i = n / MDIM;
  *j = n - (*i) * MDIM;
  return (0);
}

int
wind_ij_to_n (i, j, n)
     int *n, i, j;
{
  if (geo.coord_type == SPHERICAL)
    {
      Error
	("Warning: wind_ij_to_n being called for spherical coordinates %d %d\n",
	 i, j);
      *n = i;
      return (*n);
    }
  *n = i * MDIM + j;		// MDIM because the array is in z order
  return (*n);
}
