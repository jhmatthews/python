/* 
The overall idea for a vertically extended disk is as follows:

We will create a cylindicral coordinate system in which the rho coordinates
are established just as for a cylindrical system.  Howver, the boundaries of
the cell in the z direction will vary with radius from the star.  The 
vertices of the cells will be defined so that they have fixed offsets from
the disk surface.  The "top" and "bottom" edges are not parallel to the xy
plane, but mimic the disk surface.  

Ultimately the disk surface will be defined so that it matches the edges of
cells when a photon is in the wind.  Outise id will be defined by zdisk
as previously.

*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "atomic.h"

#include "python.h"


/***********************************************************
                                       Space Telescope Science Institute

 Synopsis:
	cylin_ds_in_cell calculates the distance to the far
        boundary of the cell in which the photon bundle resides.  	
  
 Arguments:		
 	p	Photon pointer

	
 Returns:
 	Distance to the far boundary of the cell in which the photon
	currently resides.  Negative numbers (and zero) should be
	regarded as errors.
  
Description:	

Notes:


History:
	05may	ksl	56a -- began modifications starting with same 
			routine in cylindrical
	05jun	ksl	56d -- modified so that search for top and
			bottom of each cell is determined using
			the cones that define each cell.
 
**************************************************************/



double
cylvar_ds_in_cell (p)
     PhotPtr p;


{

  int n, ix, iz, iroot;
  double a, b, c, root[2];
  double s, smax;
  int where_in_grid (), wind_n_to_ij ();
  int quadratic ();
  int radiation ();


  if ((p->grid = n = where_in_grid (p->x)) < 0)
    {
      Error ("translate_in_wind: Photon not in grid when routine entered\n");
      return (n);		/* Photon was not in wind */
    }

  wind_n_to_ij (n, &ix, &iz);	/*Convert the index n to two dimensions */

  smax = VERY_BIG;		//initialize smax to a large number

  /* Set up the quadratic equations in the radial rho direction */

  a = (p->lmn[0] * p->lmn[0] + p->lmn[1] * p->lmn[1]);
  b = 2. * (p->lmn[0] * p->x[0] + p->lmn[1] * p->x[1]);
  c = p->x[0] * p->x[0] + p->x[1] * p->x[1];

  iroot = quadratic (a, b, c - wind_x[ix] * wind_x[ix], root);	/* iroot will be the smallest positive root
								   if one exists or negative otherwise */

  if (iroot >= 0 && root[iroot] < smax)
    smax = root[iroot];

  iroot = quadratic (a, b, c - wind_x[ix + 1] * wind_x[ix + 1], root);

  if (iroot >= 0 && root[iroot] < smax)
    smax = root[iroot];

  /* At this point we have found how far the photon can travel in rho in its
     current direction.  Now we must worry about motion in the z direction  */


  s = ds_to_cone (&wmain[n].wcone, p);
  if (s < smax)
    smax = s;

  s = ds_to_cone (&wmain[n + 1].wcone, p);
  if (s < smax)
    smax = s;


  if (smax <= 0)
    {
      Error ("cylvar_ds_in_cell %f\n", smax);
    }


  return (smax);
}


/***********************************************************
                                       Space Telescope Science Institute

 Synopsis:
	cylvar_make_grid defines the cells in a cylindrical grid              
	when the disk is vertically extended.

Arguments:		
	WindPtr w;	The structure which defines the wind in Python
 
Returns:
 
Description:
	The cylvar coordinate system basically adds and off determined
	by zdisk(r) to each z position in the grid.  Beyond geo.diskrad
	the coordinate system stays fixed in the z direction (on the
	assumption that the forumula for the disk is not valid there.)

	The coordinate system is adjusted so the last element in the
	veritical direction of the grid is always the same.


History:
	05may	ksl	56a -- began modifications starting with same 
			routine in cylindrical
	05jun	ksl	56d -- Conversion for this routine appears
			complete.

**************************************************************/


int
cylvar_make_grid (w)
     WindPtr w;
{
  double dr, dz, dlogr, dlogz;
  double r, z_offset;
  int i, j, n;

  /* In order to interpolate the velocity (and other) vectors out to geo.rmax, we need
     to define the wind at least one grid cell outside the region in which we want photons
     to propagate.  This is the reason we divide by NDIM-2 here, rather than NDIM-1 */


  /* First calculate parameters that are to be calculated at the edge of the grid cell.  This is
     mainly the positions and the velocity */
  for (i = 0; i < NDIM; i++)
    {
      for (j = 0; j < MDIM; j++)
	{
	  wind_ij_to_n (i, j, &n);
	  w[n].x[1] = w[n].xcen[n] = 0;	//The cells are all defined in the xz plane

	  /*Define the grid points */
	  if (geo.log_linear == 1)
	    {			// linear intervals

	      dr = geo.rmax / (NDIM - 3);
	      w[n].x[0] = r = i * dr;	/* The first zone is at the inner radius of
					   the wind */
	      if (r < geo.diskrad)
		{
		  z_offset = zdisk (r);
		}
	      else
		z_offset = zdisk (geo.diskrad);

	      dz = (geo.rmax - z_offset) / (MDIM - 4);
	      if (j == 0)
		w[n].x[2] = 0;
	      else if (j == 1)
		w[n].x[2] = 0.9 * z_offset;	// To assure there is one comple cell in the disk
	      else
		w[n].x[2] = z_offset * (j - 1) * dz;

	      if (j == 0)
		w[n].xcen[2] = 0.5 * z_offset;
	      else
		w[n].xcen[2] = w[n].x[2] + 0.5 * dz;
	    }
	  else
	    {			//logarithmic intervals

	      dlogr = (log10 (geo.rmax / geo.xlog_scale)) / (NDIM - 3);
	      if (i == 0)
		{
		  w[n].x[0] = r = 0.0;
		  w[n].xcen[0] = 0.5 * geo.xlog_scale + zdisk (r);
		}
	      else
		{
		  w[n].x[0] = r = geo.xlog_scale * pow (10., dlogr * (i - 1));
		  w[n].xcen[0] =
		    0.5 * geo.xlog_scale * (pow (10., dlogr * (i - 1)) +
					    pow (10.,
						 dlogr * (i))) + zdisk (r);
		}

	      if (r < geo.diskrad)
		{
		  z_offset = zdisk (r);
		}
	      else
		z_offset = zdisk (geo.diskrad);

	      dlogz =
		(log10 ((geo.rmax - z_offset) / geo.zlog_scale)) / (MDIM - 4);
	      if (j == 0)
		{
		  w[n].x[2] = 0;
		  w[n].xcen[2] = 0.5 * z_offset;
		}
	      else if (j == 1)
		{
		  w[n].x[2] = 0.9 * z_offset;	// 0.9 is to keep the first grid cell out of wind
		  w[n].xcen[2] = z_offset + 0.5 * geo.zlog_scale;
		}
	      else
		{
		  w[n].x[2] = z_offset +
		    geo.zlog_scale * pow (10, dlogz * (j - 2));
		  w[n].xcen[2] = z_offset +
		    0.5 * geo.zlog_scale * (pow (10., dlogz * (j - 2)) +
					    pow (10., dlogz * (j - 1)));
		}
	    }


	}
    }

  return (0);
}


/***********************************************************
                                       Space Telescope Science Institute

 Synopsis:
	cylvar_wind_complete_grid completes the definition of some of 
	the positonal variables in the wind, and creates some subidiary
	arrays.

Arguments:		
	WindPtr w;	The structure which defines the wind in Python
 
Returns:
 
Description:
	This routine defines the windcone structures for this coordinate
	system and creates some arrays that are intended to aide
	in determining what cell one is in.  The routine is specific
	to cylvar coordinate system.  


Notes: In principle, one would like to fix problem that the cones
	are only defined out to MDIM-1 and NDIM-1, as this has the
	potential of wasting space.  However, volumes have the same problem
History:
	05may	ksl	56a -- began modifications starting with same 
			routine in cylindrical
	06jun	ksl	56d -- Completed models for cylvar

**************************************************************/


int
cylvar_wind_complete (w)
     WindPtr w;
{
  int i, j, n;
  double drho, dz;

  /* First define the windcones for each cell */
  n = 0; // silence compiler warning

  for (i = 0; i < NDIM - 1; i++)
    {
      for (j = 0; j < MDIM - 1; j++)
	{
	  /* Everything is defined in the xz plane, so life is fairly easy */
	  n = i * MDIM + j;
	  dz = w[n + MDIM].x[2] - w[n].x[2];	// change in z for one step in rho
	  drho = w[n + MDIM].x[0] - w[n].x[0];


	  if (drho > 0)
	    w[n].wcone.dzdr = dz / drho;	//new definition
	  else
	    {
	      Error ("cylvar_wind_complete: drho == 0\n");
	      exit (0);
	    }
	  w[n].wcone.z = w[n].x[2] - dz / drho * w[n].x[0];	//new definition
	}
    }
  /* Finally define some one-d vectors that make it easier to locate a photon in the wind given that we
     have adoped a "rectangular" grid of points.   

     For cylvar coordinates, we can use the standard formulation in the x direction, but not in
     the z direction

   */

  for (i = 0; i < NDIM; i++)
    {
      wind_x[i] = w[i * MDIM].x[0];
      wind_midx[i] = w[i * MDIM].xcen[0];
    }


  for (i = 0; i < NDIM; i++)
    {
      for (j = 0; j < MDIM; j++)
	{
	  wind_ij_to_n (i, j, &n);
	  wind_z_var[i][j] = w[n].x[2];
	  wind_midz_var[i][j] = w[n].xcen[2];
	}
    }

  return (0);
}

/***********************************************************
                                       Space Telescope Science Institute

 Synopsis:
 	cylvar_volume(w) calculates the wind volume of a cylindrical cell
	allowing for the fact that some cells 

 Arguments:		
	WindPtr w;    the entire wind
 Returns:

 Description:

 	This is a brute_force integration of the volume.  The technique
	is completely general.  It does not depend on the shape of the
	cell, and could be used for most of the coodinate systems.
	       
 Notes:
	Where_in grid does not tell you whether the photon is in the wind or not. 
 History:
	05may	ksl	56a -- began modifications starting with same 
			routine in cylindrical
	05jul	ksl	56d -- Made the modifications needed.
	06nov	kls	58b: Minor modifications to use W_ALL_INWIND, etc.
			instead of hardcoded values
	11aug	ksl	70b - Added ability to get volumes for multiple
			components
 
**************************************************************/

#define RESOLUTION   1000


int
cylvar_volumes (w, icomp)
     WindPtr w;
     int icomp;			// The component for which we want the volume
{
  int i, j, n;
  int jj, kk;
  double r, z;
  double rmax, rmin;
  double zmin, zmax;
  double dr, dz, x[3];
  double volume;
  double f, g;
  int bilin ();


  /* Initialize all the volumes to 0 */
  for (n = 0; n < NDIM2; n++)
    {
      w[n].vol = 0;
      w[n].inwind = W_NOT_INWIND;
    }

  for (i = 0; i < NDIM - 1; i++)
    {
      for (j = 0; j < MDIM - 1; j++)
	{
	  wind_ij_to_n (i, j, &n);

	  /* Encapsulate the grid cell with a rectangle for integrating */
	  rmin = w[n].x[0];
	  if (rmin > w[n + 1].x[0])
	    rmin = w[n + 1].x[0];
	  if (rmin > w[n + MDIM].x[0])
	    rmin = w[n + MDIM].x[0];
	  if (rmin > w[n + MDIM + 1].x[0])
	    rmin = w[n + MDIM + 1].x[0];

	  rmax = w[n].x[0];
	  if (rmax < w[n + 1].x[0])
	    rmax = w[n + 1].x[0];
	  if (rmax < w[n + MDIM].x[0])
	    rmax = w[n + MDIM].x[0];
	  if (rmax < w[n + MDIM + 1].x[0])
	    rmax = w[n + MDIM + 1].x[0];

	  zmin = w[n].x[2];
	  if (zmin > w[n + 1].x[2])
	    zmin = w[n + 1].x[2];
	  if (zmin > w[n + MDIM].x[2])
	    zmin = w[n + MDIM].x[2];
	  if (zmin > w[n + MDIM + 1].x[2])
	    zmin = w[n + MDIM + 1].x[2];

	  zmax = w[n].x[2];
	  if (zmax < w[n + 1].x[2])
	    zmax = w[n + 1].x[2];
	  if (zmax < w[n + MDIM].x[2])
	    zmax = w[n + MDIM].x[2];
	  if (zmax < w[n + MDIM + 1].x[2])
	    zmax = w[n + MDIM + 1].x[2];


	  volume = 0;
	  jj = kk = 0;
	  dr = (rmax - rmin) / RESOLUTION;
	  dz = (zmax - zmin) / RESOLUTION;
	  for (r = rmin + dr / 2; r < rmax; r += dr)
	    {
	      for (z = zmin + dz / 2; z < zmax; z += dz)
		{
		  x[0] = r;
		  x[1] = 0;
		  x[2] = z;
		  if (bilin
		      (x, w[i * MDIM + j].x, w[i * MDIM + j + 1].x,
		       w[(i + 1) * MDIM + j].x, w[(i + 1) * MDIM + j + 1].x,
		       &f, &g) == 0)
		    {
		      kk++;
		      if (where_in_wind (x) == W_ALL_INWIND)
			{
			  volume += r;
			  jj++;
			}
		    }

		}

	    }
	  w[n].vol = 4. * PI * dr * dz * volume;
	  /* OK now make the final assignement of nwind and fix the volumes */
	  if (jj == 0)
	    {
	      w[n].inwind = W_NOT_INWIND;	// The cell is not in the wind
	    }
	  else if (jj == kk)
	    {
	      //OLD 70b w[n].inwind = W_ALL_INWIND;     // All of cell is inwind
	      w[n].inwind = icomp;	// All of cell is inwind
	    }
	  else
	    {
	      //OLD 70b w[n].inwind = W_PART_INWIND;    // Some of cell is inwind
	      w[n].inwind = icomp + 1;	// Some of cell is inwind
	    }
	}
    }

  return (0);
}

/***********************************************************
                     Space Telescope Science Institute

 Synopsis:
 	cylvar_where_in_grid locates the grid position of the vector,
	when one is using cylindrical coordinates, with the z
	positions varying. 

 Arguments:		
	double x[];
	int ichoice	0 locates the grid cell using the normal cell boundaries
			1 is for the cell centers

 Returns:
 	where_in_grid normally  returns the cell number associated with
 		a position.  If the position is in the grid this will be a positive
 		integer < NDIM*MDIM.
 	x is inside the grid        -1
	x is outside the grid       -2

	fx,fz are the fractional positions
 Description:	

 	cylvar_where_in_grid returns the grid cell, and (because
	it has to be calculated) the fractional position).  
	
		
 Notes:
	Where_in grid does not tell you whether the x is in the wind or not. 

	What one means by inside or outside the grid may well be different
	for different coordinate systems.

	The basic approach is to calculate the rho postion which should tell
	you which "column" of the wind you are in, and then to calculate the
	specific cell in the column.  

 History:
	05may	ksl	56a -- began modifications starting with same 
			routine in cylindrical
	05jul	ksl	56d -- Actual modifications for cylvar added.
			This routine caused a lot of trouble, especially
			in cases where there were consequences when x
			was outside the grid.
	13sep	nsh	76b -- Modified calls to fraction to take account
			of new modes
 
**************************************************************/


int cylvar_n_approx;
int ierr_cylvar_where_in_grid = 0;

int
cylvar_where_in_grid (x, ichoice, fx, fz)
     double x[];
     int ichoice;
     double *fx, *fz;
{
  int i, j, n, ii;
  double z[3];
  double rho;
  int where_in_2dcell ();

  /* copy x to a dummy vector z, so that the z[0] component is really rho */
  z[0] = rho = sqrt (x[0] * x[0] + x[1] * x[1]);
  z[1] = 0;
  z[2] = fabs (x[2]);		/* This is necessary to get correct answer above
				   and below plane */
  if (z[2] == 0)
    z[2] = 1.e4;		//Force z to be positive  02feb ksl

  /* Attempt to find the approximate locations. Recall that ichoce == 0 
   * corresponds to the cells while ichoice == 1 cooresponds to the mid-points
   * of the cells. */

  if (ichoice == 0)
    {
      fraction (rho, wind_x, NDIM, &i, fx, 0);
      fraction (z[2], wind_z_var[i], MDIM, &j, fz, 0);	// This should get one close
    }
  else
    {
      fraction (rho, wind_midx, NDIM, &i, fx, 0);
      fraction (z[2], wind_midz_var[i], MDIM, &j, fz, 0);	// This should get one close
    }

  wind_ij_to_n (i, j, &cylvar_n_approx);

  /* Check to see if x is outside the region of the calculation.  Note that this
   * caclulation cannot be done until i is determined */

  if (rho < wind_x[0])
    return (-1);

  if (rho > wind_x[NDIM - 1] || z[2] > wind_z_var[i][MDIM - 1])
    {
      return (-2);		/* x is outside grid */
    }

  /* 
     Now check if is in this cell. If the actual j is larger than the trial j, then
     bilin returns a value of 0 the grid point is in the cell, -1 otherwise
     For cylvar coordinates g tells one whether to go up or down in the grid
   */

  n = cylvar_n_approx;

  ii = where_in_2dcell (ichoice, z, n, fx, fz);
  if (*fx < 0 || 1 < *fx)
    {
      Error
	("cylvar_where_in_grid: fraction and where_in_2d_cell incompatible %d -- %e %e\n",
	 n, z[0], z[2]);
      fraction (rho, wind_midx, NDIM, &i, fx, 0);
      ii = where_in_2dcell (ichoice, z, n, fx, fz);
    }

//  while (ii != 0 && 0 <= j && j < MDIM)
  while (ii != 0)
    {
      if (*fz < 0 && j == 0)
	{			// we are out of the grid and into the disk
	  *fz = 0;
	  if (ierr_cylvar_where_in_grid < 100)
	    {
	      Error
		("cylvar_where_in_grid: Position %f %f %f gavegrid cell in disk %d (%d %d)\n",
		 z[0], z[1], z[2], n, i, j);
	      ierr_cylvar_where_in_grid++;
	    }
	  break;
	}
      if (*fz > 1. && j == MDIM)
	{
	  // we are out of the grid at the "top" edge * fz = 0;
	  n--;
	  if (ierr_cylvar_where_in_grid < 100)
	    {
	      Error
		("cylvar_where_in_grid: Position %f %f %f gave grid cell above top edge %d (%d %d)\n",
		 z[0], z[1], z[2], n, i, j);
	      ierr_cylvar_where_in_grid++;
	    }
	  break;
	}

      /* OK we are still in the grid so we can increment or decrement j */

      if (*fz >= 1.0)
	{
	  j++;
	  n++;
	}
      else
	{
	  j--;
	  n--;
	}
      ii = where_in_2dcell (ichoice, z, n, fx, fz);
      if (*fx < 0 || 1 < *fx)
	{
	  Error
	    ("cylvar_where_in_grid: fraction and where_in_2d_cell incompatible %d -- %e %e (breaking)\n",
	     n, z[0], z[2]);
	  break;
	}
    }

  // Next error checks should not be needed

  if (j < 0 || n < 0)
    {
      if (ierr_cylvar_where_in_grid < 100)
	{
	  Error
	    ("cylvar_where_in_grid: Position %f %f %f gave unreasonable grid cell %d (%d %d)\n",
	     z[0], z[1], z[2], n, i, j);
	  ierr_cylvar_where_in_grid++;
	}
      if (n < 0)
	return (0);
      else
	n++;
    }

  return (n);

}

/***********************************************************
                     Space Telescope Science Institute

 Synopsis:

 	cylvar_get_random_location

 Arguments:		
 	int n -- Cell in which random position is to be generated
	int icomp - Component in which the random position is to
		be generated
 Returns:
 	double x -- the position
 Description:	
	
		
 Notes:
 	The encapsulation steps should probably put into a separate
	routine since this is used several times

 History:
	05may	ksl	56a -- began modifications starting with same 
			routine in cylindrical
	05jul	ksl	56d -- Updated to work for cylvar coords
	11aug	ksl	70b - Updated to include more than one
			component
 
**************************************************************/

int
cylvar_get_random_location (n, icomp, x)
     int n;			// Cell in which to create position
     int icomp;			// Component in which to create position
     double x[];		// Returned position
{
  int i, j;
  int inwind, incell;
  double r, rmin, rmax, zmin, zmax;
  double fx, fz;
  double zz;
  double phi;
  int where_in_2dcell ();

  wind_n_to_ij (n, &i, &j);

  /* Encapsulate the grid cell with a rectangle for integrating */

  rmin = wmain[n].x[0];
  if (rmin > wmain[n + 1].x[0])
    rmin = wmain[n + 1].x[0];
  if (rmin > wmain[n + MDIM].x[0])
    rmin = wmain[n + MDIM].x[0];
  if (rmin > wmain[n + MDIM + 1].x[0])
    rmin = wmain[n + MDIM + 1].x[0];

  rmax = wmain[n].x[0];
  if (rmax < wmain[n + 1].x[0])
    rmax = wmain[n + 1].x[0];
  if (rmax < wmain[n + MDIM].x[0])
    rmax = wmain[n + MDIM].x[0];
  if (rmax < wmain[n + MDIM + 1].x[0])
    rmax = wmain[n + MDIM + 1].x[0];

  zmin = wmain[n].x[2];
  if (zmin > wmain[n + 1].x[2])
    zmin = wmain[n + 1].x[2];
  if (zmin > wmain[n + MDIM].x[2])
    zmin = wmain[n + MDIM].x[2];
  if (zmin > wmain[n + MDIM + 1].x[2])
    zmin = wmain[n + MDIM + 1].x[2];

  zmax = wmain[n].x[2];
  if (zmax < wmain[n + 1].x[2])
    zmax = wmain[n + 1].x[2];
  if (zmax < wmain[n + MDIM].x[2])
    zmax = wmain[n + MDIM].x[2];
  if (zmax < wmain[n + MDIM + 1].x[2])
    zmax = wmain[n + MDIM + 1].x[2];


  /* Generate a position which is both in the cell and in the wind */
  inwind = incell = -1;
  while (inwind != icomp || incell != 0)
    {
      r =
	sqrt (rmin * rmin +
	      (rand () / (MAXRAND - 0.5)) * (rmax * rmax - rmin * rmin));

// Generate the azimuthal location
      phi = 2. * PI * (rand () / MAXRAND);
      x[0] = r * cos (phi);
      x[1] = r * sin (phi);



      x[2] = zmin + (zmax - zmin) * (rand () / (MAXRAND - 0.5));
      inwind = where_in_wind (x);	/* Some photons will not be in the wind
					   because the boundaries of the wind split the grid cell */
      incell = where_in_2dcell (0, x, n, &fx, &fz);
    }

  zz = rand () / MAXRAND - 0.5;	//positions above are all at +z distances

  if (zz < 0)
    x[2] *= -1;			/* The photon is in the bottom half of the wind */

  return (inwind);
}


/***********************************************************
                     Space Telescope Science Institute

 Synopsis:
 	cylvar_extend_density  extends the density to
	regions just outside the wind regiions so that
	extrapolations of density can be made there

 Arguments:		
 Returns:
 Description:	

     We need to updated the densities immediately outside the wind so that the density interpolation in resonate will work.
     In this case all we have done is to copy the densities from the cell which is just in the wind (as one goes outward) to the
     cell that is just inside (or outside) the wind. 

     In cylindrical coordinates, the fast dimension is z; grid positions increase up in z, and then out in r.
     In spperical polar coordinates, the fast dimension is theta; the grid increases in theta (measured)
     from the z axis), and then in r.

		
 Notes:

	Instead of copying every thing is accomplished using pointers.	

 History:
	05may	ksl	56a -- began modifications starting with same 
			routine in cylindrical
	05jul	ksl	56d -- Actually nothing was done to this
			routine. It should be OK as is.
	06may	ksl	57+ -- Updated to extend by changing the mapping into the plasma structure.
 
**************************************************************/


int
cylvar_extend_density (w)
     WindPtr w;
{

  int i, j, n, m;

  for (i = 0; i < NDIM - 1; i++)
    {
      for (j = 0; j < MDIM - 1; j++)
	{
	  wind_ij_to_n (i, j, &n);
	  if (w[n].vol == 0)

	    {			//Then this grid point is not in the wind 

	      wind_ij_to_n (i + 1, j, &m);

	      if (w[m].vol > 0)
		{
		  w[n].nplasma = w[m].nplasma;

		}
	      else if (i > 0)
		{
		  wind_ij_to_n (i - 1, j, &m);
		  if (w[m].vol > 0)
		    {
		      w[n].nplasma = w[m].nplasma;

		    }
		}
	    }
	}
    }

  return (0);

}




/***********************************************************
           Space Telescope Science Institute

Synopsis:  cylvar_coord_fraction() calculates the fractional
	position of a position in the 2d grid in a coordinate
	sytem independent way.  



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
		       	

**************************************************************/



int
cylvar_coord_fraction (ichoice, x, ii, frac, nelem)
     int ichoice;
     double x[];
     int ii[];
     double frac[];
     int *nelem;
{
  double dr, dz;
  int n;


  if ((n = cylvar_where_in_grid (x, ichoice, &dr, &dz)) < 0)
    {
      n = cylvar_n_approx;
    }


  ii[0] = n;
  frac[0] = (1. - dz) * (1. - dr);

  ii[1] = n + MDIM;
  frac[1] = (1. - dz) * dr;

  ii[2] = n + 1;
  frac[2] = (dz) * (1. - dr);

  ii[3] = n + MDIM + 1;
  frac[3] = (dz) * (dr);
  *nelem = 4;

  if (sane_check (dr) || sane_check (dz))
    {
      Error ("coord_frac:sane_check dr=%f dz=%f for 2d coords\n", dr, dz);
    }

  return (1);

}
