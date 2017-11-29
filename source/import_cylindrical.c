#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "atomic.h"
#include "python.h"

/***********************************************************
    Space Telescope Science Institute

Synopsis:
    These are general purpose routines for readin in model
    grids

Arguments:		

Returns:
 
Description:	

Notes:



History:
	17nov   ksl Began coding                           
**************************************************************/


# define LINELEN 512
# define NCELLS  512

/* The structure that holds the inputs and any subsidiary variables
 *
 * Note that i is the row number nd j is the column number*/
struct
{
  int ndim, mdim, ncell;
  int i[NDIM_MAX * NDIM_MAX], j[NDIM_MAX * NDIM_MAX],
    inwind[NDIM_MAX * NDIM_MAX];
  double x[NDIM_MAX * NDIM_MAX], z[NDIM_MAX * NDIM_MAX];
  double v_x[NDIM_MAX * NDIM_MAX], v_y[NDIM_MAX * NDIM_MAX],
    v_z[NDIM_MAX * NDIM_MAX];
  double rho[NDIM_MAX * NDIM_MAX], t[NDIM_MAX * NDIM_MAX];

  double wind_x[NDIM_MAX], wind_z[NDIM_MAX], wind_midx[NDIM_MAX],
    wind_midz[NDIM_MAX];
} xx_cyl;


/***********************************************************
    Space Telescope Science Institute

Synopsis:
    Read the an arbitray wind model in cylindrical
    coordinates


Arguments:		

Returns:
 
Description:	

Notes:

    The basic data we need to read in are

    i, j, r z  v_x v_y v_z rho (and optionally T)

    where v_x,_v_y,v_z are the velocities in the x,z plane
    and where

    i is the column number  (Thus i corresponds to ndim)
    j is the row number     (and z coresponds to mdim)

    We assume that all of the variables are centered, that is
    we are not assuming that we are giving rho at the center of
    a cell, but that r and v_r are at the edges of a cell. 
    This is someghing that would presumable be easy to change

    Note that we assume that the data are being read in in the
    same order as printed out by windsave2table, that is that
    the first "column" is read in, and then the second "column".




History:
	17nov   ksl Began coding                           
**************************************************************/

int
import_cylindrical (ndom, filename)
     int ndom;
     char *filename;
{
  FILE *fopen (), *fptr;
  char line[LINELEN];
  int n, icell, jcell, ncell, inwind;
  double q1, q2, q3, q4, q5, q6, q7;
  int jz, jx;
  double delta;

  Log ("Reading a model in cylindrical coordinates %s\n", filename);

  if ((fptr = fopen (filename, "r")) == NULL)
    {
      Error ("import_cylindrical: No such file\n");
      exit (0);
    }


  ncell = 0;
  while (fgets (line, 512, fptr) != NULL)
    {
      n =
	sscanf (line, " %d %d %d %le %le %le %le %le %le %le", &icell, &jcell,
		&inwind, &q1, &q2, &q3, &q4, &q5, &q6, &q7);
      if (n < 4)
	{
	  printf ("Error. Ignore %s \n", line);
	  continue;
	}
      else
	{
	  xx_cyl.i[ncell] = icell;
	  xx_cyl.j[ncell] = jcell;
	  xx_cyl.inwind[ncell] = inwind;
	  xx_cyl.x[ncell] = q1;
	  xx_cyl.z[ncell] = q2;
	  xx_cyl.v_x[ncell] = q3;
	  xx_cyl.v_y[ncell] = q4;
	  xx_cyl.v_z[ncell] = q5;
	  xx_cyl.rho[ncell] = q6;
	  if (n > 10)
	    {
	      xx_cyl.t[ncell] = q7;
	    }
	  else
	    {
	      xx_cyl.t[ncell] = 10000.;
	    }
	  ncell++;

	}
    }


/* Having read in the data define some initial variables concerning the model. We cannot create
 * the wind grid or other things at this point, because we do not at this point know what
 * wind cells correspnd to whate elements of the grid */

  xx_cyl.ndim = icell + 1;
  xx_cyl.mdim = jcell + 1;
  xx_cyl.ncell = ncell;



  jz = jx = 0;
  for (n = 0; n < xx_cyl.ncell; n++)
    {
      if (xx_cyl.i[n] == 0)
	{
	  xx_cyl.wind_z[jz] = xx_cyl.z[n];
	  jz++;
	}
      if (xx_cyl.j[n] == 0)
	{
	  xx_cyl.wind_x[jx] = xx_cyl.x[n];
	  jx++;
	}
    }


  /* Now fill in wind_midx and midz. Given how we construct
   * the mdpts we need to add one more on the end, though it
   * is not entirely obvious that this is needed, given the
   * assumption that we do not need extra buffer cells */

  for (n=0;n<jz-1;n++){
      xx_cyl.wind_midz[n]=0.5*(xx_cyl.wind_z[n]+xx_cyl.wind_z[n+1]);
  }

//OLD  delta=(xx_cyl.wind_z[jz-1]-xx_cyl.wind_z[jz-2]);
//OLD  xx_cyl.wind_midz[jz]=xx_cyl.wind_z[jz-1]+0.5*delta;

  delta=(xx_cyl.wind_z[n-1]-xx_cyl.wind_z[n-2]);
  xx_cyl.wind_midz[n]=xx_cyl.wind_z[n-1]+0.5*delta;



  for (n=0;n<jx-1;n++){
      xx_cyl.wind_midx[n]=0.5*(xx_cyl.wind_x[n]+xx_cyl.wind_x[n+1]);
  }

//OLD  delta=(xx_cyl.wind_x[n-1]-xx_cyl.wind_x[n-2]);
//OLD  xx_cyl.wind_midx[jx]=xx_cyl.wind_x[jx-1]+0.5*delta;


  delta=(xx_cyl.wind_x[n-1]-xx_cyl.wind_x[n-2]);
  xx_cyl.wind_midx[n]=xx_cyl.wind_x[n-1]+0.5*delta;


  /* Although the initialization of most of zdom should be postponed
   * one has to give zdom the dimensions of the array; otherwise 
   * the wrong number of elements in wmains wind will be allocated
   */

  zdom[ndom].ndim = xx_cyl.ndim;
  zdom[ndom].mdim = xx_cyl.mdim;
  zdom[ndom].ndim2=zdom[ndom].ndim * zdom[ndom].mdim;


  return (0);
}


/* The next section contains routines to make the grids for imported models.

   We make some assumptions here.  We assume that every cell is in the wind.
   and we assume that r refers to the inside edge of the cell.

 * */

int
cylindrical_make_grid_import (w, ndom)
     WindPtr w;
     int ndom;
{
  int n;
  int nn;
  double r, rmin, rmax, rho_min, rho_max, zmax;
  double x[3];

  Log ("XX Dimensions of read in model: %d %d\n", zdom[ndom].ndim,
       zdom[ndom].mdim);

/*  This is an attempt to make the grid directly.  
 *
 *  Note also that none of this will work unless a complete grid is read 
 *  in
 *  */
  for (n = 0; n < xx_cyl.ncell; n++)
    {
      wind_ij_to_n (ndom, xx_cyl.i[n], xx_cyl.j[n], &nn);
      w[nn].x[0] = xx_cyl.x[n];
      w[nn].x[1] = 0;
      w[nn].x[2] = xx_cyl.z[n];
      w[nn].v[0] = xx_cyl.v_x[n];
      w[nn].v[1] = xx_cyl.v_y[n];
      w[nn].v[2] = xx_cyl.v_z[n];
      w[nn].inwind = xx_cyl.inwind[n];

      w[nn].xcen[0] = xx_cyl.wind_midx[xx_cyl.i[n]];
      w[nn].xcen[1] = 0;
      w[nn].xcen[2] = xx_cyl.wind_midz[xx_cyl.j[n]];

      /* JM 1711 -- copy across the inwind variable to the wind pointer */
      w[nn].inwind = xx_cyl.inwind[n];

      /* JM 1711 -- you're either in, or you're out. No part in wind cells allowed! 
         there is a question here about which choice (of not in wind or all in 
         wind) is most appropriate */
      if (w[nn].inwind == W_PART_INWIND)
        w[nn].inwind = W_NOT_INWIND;
    }

// ksl - Next few lines do not do anything
//OLD  /* We now need to fill in the w[],cen */

//OLD  for (n = 0; n < zdom[ndom].ndim2; n++)
//OLD    {
//OLD      wind_ij_to_n (ndom, xx_cyl.i[n], xx_cyl.j[n], &nn);
//OLD
//OLD    }
  

  /* Now add information used in zdom */

  for (n = 0; n < zdom[ndom].ndim; n++)
    {
      zdom[ndom].wind_x[n] = xx_cyl.wind_x[n];
    }



  for (n = 0; n < zdom[ndom].mdim; n++)
    {
      zdom[ndom].wind_z[n] = xx_cyl.wind_z[n];
    }


  /* Now set up wind boundaries so they are harmless.  
   * Note that given that we have already filled out
   * the WindPtr it seems like we could do this
   * without needing the internal structure
   * I am attempting this in the rtheta model in 
   * this fashion - ksl  */


  rmax = rho_max = zmax = 0;
  rmin = rho_min = VERY_BIG;
  for (n = 0; n < xx_cyl.ncell; n++)

    {
      x[0] = xx_cyl.x[n];
      x[1] = 0;
      x[2] = xx_cyl.z[n];

      r = length (x);

      if (xx_cyl.inwind[n] >= 0)
	{
	  if (xx_cyl.x[n] > rho_max)
	    {
	      rho_max = xx_cyl.x[n];
	    }
	  if (xx_cyl.z[n] > zmax)
	    {
	      zmax = xx_cyl.z[n];
	    }
	  if (r > rmax)
	    {
	      rmax = r;
	    }
	}
      else
	{
	  if (rho_min > xx_cyl.x[n])
	    {
	      rho_min = xx_cyl.x[n];
	    }
	  if (rmin > r)
	    {
	      rmin = r;
	    }
	}
    }



  zdom[ndom].wind_rho_min = zdom[ndom].rho_min = rho_min;
  zdom[ndom].wind_rho_max = zdom[ndom].rho_max = rho_max;
  zdom[ndom].zmax = zmax;

  zdom[ndom].rmax = rmax;
  zdom[ndom].rmin = rmin;
  zdom[ndom].wind_thetamin = zdom[ndom].wind_thetamax = 0.;

  return (0);
}


/* The next section calculates velocites.  We follow the hydro approach of
 * getting those velocities from the original grid.  This is really only
 * used for setting up the grid
 */


/* velocity for cylindrical coordinates only interpolates.  One has to
 * interpoalte for vgrad. 
 *
 * This routine is dangerous because of the way it works, if one tries
 * to update wmain[].v. This is because it actually uses values in 
 * wmain[v].  Ideally one would try to avoid maing such calls, but
 * ksl has not found a good way to do this.
 * */


double
velocity_cylindrical (ndom, x, v)
     int ndom;
     double *x, *v;
{
  int j;
  int nn;
  int nnn[4], nelem;
  double frac[4];
  double vv[3];
  double speed;
  coord_fraction (ndom, 0, x, nnn, frac, &nelem);
  for (j = 0; j < 3; j++)
    {
      vv[j] = 0;
      for (nn = 0; nn < nelem; nn++)
	vv[j] += wmain[zdom[ndom].nstart + nnn[nn]].v[j] * frac[nn];
    }

  speed = length (vv);

  /* Now copy the result into v, which is very necessary if refilling wmain.v */

  v[0] = vv[0];
  v[1] = vv[1];
  v[2] = vv[2];

  return (speed);
}





/* Fill in plasma ptrs with densities.   
 *
 * For this we assume that the densities read in are 
 * given at the * midpoints of the grid
 *
 * 
 * */

/*  This routine should only be called to set up the plasma cells,
 *  and we assume that rho that was imported is the center of the 
 *  plasma cell, so there is no need to interpolate.
 */

double
rho_cylindrical (ndom, x)
     int ndom;
     double *x;
{
  double rho = 0;
  double r, z;
  int i, j, n;

  r = sqrt (x[0] * x[0] + x[1] * x[1]);
  z = fabs (x[2]);

  i = 0;
  while (z > xx_cyl.wind_z[i] && i < xx_cyl.mdim - 1)
    {
      i++;
    }
  j = 0;
  while (r > xx_cyl.wind_x[j] && j < xx_cyl.ndim - 1)
    {
      j++;
    }

  n = j * xx_cyl.mdim + i;

  rho = xx_cyl.rho[n];


  return (rho);
}
