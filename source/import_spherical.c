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

/* The next variables have to be external because we need them to be available later on */

struct
{
  int ndim;
  int element[NDIM_MAX];
  double r[NDIM_MAX], v[NDIM_MAX], rho[NDIM_MAX], t[NDIM_MAX];
} xx_1d;



/***********************************************************
    Space Telescope Science Institute

Synopsis:
    Read the an arbitray wind model intended to mimic a stellar
    wind or shell.


Arguments:		

Returns:
 
Description:	

Notes:

    The basic data we need to read in are

    i r v_r rho (and optionally T)

    We assume that all of the variables are centered, that is
    we are not assuming that we are giving rho at the center of
    a cell, but that r and v_r are at the edges of a cell. 
    This is someghing that would presumable be easy to change


History:
	17nov   ksl Began coding                           
**************************************************************/

int
import_1d (ndom, filename)
     int ndom;
     char *filename;
{
  FILE *fopen (), *fptr;
  char line[LINELEN];
  int n, icell, ncell;
  double q1, q2, q3, q4;


  Log ("Reading a 1d model %s\n", filename);


  if ((fptr = fopen (filename, "r")) == NULL)
    {
      Error ("import_1d: No such file\n");
      exit (0);
    }


  ncell = 0;
  while (fgets (line, 512, fptr) != NULL)
    {
      n = sscanf (line, " %d %le %le %le %le", &icell, &q1, &q2, &q3, &q4);
      if (n < 4)
	{
	  continue;
	}
      else
	{
	  xx_1d.element[ncell] = icell;
	  xx_1d.r[ncell] = q1;
	  xx_1d.v[ncell] = q2;
	  xx_1d.rho[ncell] = q3;
	  if (n > 4)
	    {
	      xx_1d.t[ncell] = q4;
	    }
	  else
	    {
	      xx_1d.t[ncell] = 10000.;
	    }
	  ncell++;

	}
    }


  xx_1d.ndim = ncell;

  /* Although much of the initialization of zdom can be postponeed
   * one has to define mdim and ndim of zdom here, so that the correct
   * number of wind cells will be allocated */

  zdom[ndom].ndim2=zdom[ndom].ndim = xx_1d.ndim;
  zdom[ndom].mdim = 1;

  return (0);
}


/* The next section contains routines to make the grids for imported models.

   We make some assumptions here.  We assume that every cell is in the wind.
   and we assume that r refers to the inside edge of the cell.

 * */

int
spherical_make_grid_import (w, ndom)
     WindPtr w;
     int ndom;
{

  int j, n;

//  zdom[ndom].ndim2=zdom[ndom].ndim = xx_1d.ndim;	// ADD Buffer
//  zdom[ndom].mdim = 1;
  zdom[ndom].wind_rho_min = zdom[ndom].rho_min = 0;
  zdom[ndom].rmin = xx_1d.r[0];
  zdom[ndom].wind_rho_max = zdom[ndom].zmax = zdom[ndom].rho_max =
    zdom[ndom].rmax = xx_1d.r[xx_1d.ndim - 1];
  zdom[ndom].wind_thetamin = zdom[ndom].wind_thetamax = 0.;

  for (j = 0; j < xx_1d.ndim; j++)
    {
      n = j + zdom[ndom].nstart;
      w[n].r = xx_1d.r[j];
      /* Put the radial velocity in v[0] */
      w[n].v[0]=xx_1d.v[j];
    }


      /* Need to define the midpoints of the grid */

  for (j = 0; j < zdom[ndom].ndim; j++)
    {
      n = j + zdom[ndom].nstart;
      if (j < zdom[ndom].ndim - 1)
	{
	  w[n].rcen = 0.5 * (w[n].r + w[n + 1].r);
	}
      else
	{
	  w[n].rcen = w[n].r * 1.005;
	}
      w[n].x[1] = w[n].xcen[1] = 0.0;
      w[n].x[0] = w[n].x[2] = w[n].r * sin (PI / 4.);
      w[n].xcen[0] = w[n].xcen[2] = w[n].rcen * sin (PI / 4.);
    }

  spherical_wind_complete(ndom,w);
  return (0);
}


/* The next section calculates velocites.  
 *
 * One could follow the zeus_hydro approach of getting those velocities from the original grid.  
 * but for consistency with the 2d case we get it by interpolating on values in the cells
 *
 * Note that v_r is stored in v_0
 *
 */

double
velocity_1d (ndom, x, v)
     int ndom;
     double *x, *v;
{
  double speed;
  double r;
  int nelem,nn, nnn[4];
  double frac[4];
  //
//OLD  int icell;
  r = length (x);

  
  coord_fraction (ndom, 0, x, nnn, frac, &nelem);
  speed=0;
  for (nn = 0; nn < nelem; nn++)
  {
      speed+=wmain[zdom[ndom].nstart + nnn[nn]].v[0] * frac[nn];
  }

//OLD  icell = linterp (r, xx_1d.r, xx_1d.v, xx_1d.ndim, &speed, 0);

  v[0] = x[0] / r * speed;
  v[1] = x[1] / r * speed;
  v[2] = x[2] / r * speed;
  return (speed);
}

/* Fill in plasma ptrs with densities.   
 *
 * For this we assume that the densities read in are 
 * given at the * midpoints of the grid
 *
 * 
 * */


double
rho_1d (ndom, x)
     int ndom;
     double *x;
{
  double rho = 0;
  double r;
  int n;
  r = length (x);
  n = 0;
  while (r > xx_1d.r[n] && n < xx_1d.ndim)
    {
      n++;
    }

  if (n < xx_1d.ndim)
    {
      rho = xx_1d.rho[n];
    }
  else
    {
      rho = xx_1d.rho[xx_1d.ndim - 1];
    }


  Log ("rho %e \n", rho);
  return (rho);
}

