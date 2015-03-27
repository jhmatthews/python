
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "atomic.h"
#include "python.h"

/***********************************************************
                                       Space Telescope Science Institute

 Synopsis:
	These routines define a plane parallel wind.  Initially they are intended only for use
	with balance.  
 
Arguments:		
	
Returns:
 
Description:	
	The wind is defined so that the direction in which things vary is the (+) z direction.
	Normally one would hope that photons start at - infinity or 0.
Notes:
	All of the routines should begin with letters pl, eg pl_wind_define.  Since ultimately
	we want to be able to define in a standardized fashion various types of winds.

History:
 	97jan	ksl	Coding on python began.
 	98dec	ksl	Coding on these routines begain
	07jul	ksl	58f - Still only applies to balance
 
**************************************************************/

int
pl_wind_define (w)
     WindPtr w;			// This is the entire grid
{
  double x[3], vmid[3], vbot[3], diff[3];
  double xmin, xmax, zmin, zmax;
  int i, j, n, ierr;
  double pl_rho (), pl_velocity ();
  int define_wind_grid ();
  double pl_rho (), length ();
  int pl_copy_conditions ();

  geo.log_linear = 1;		// At present do not worry about a wind which is logarithmithically
  // spaced!

  xmin = 0;
  xmax = (NDIM - 2) * 1.0;
  zmin = 0;
  zmax = (MDIM - 2) * geo.pl_vol;

  define_wind_grid (w, xmin, xmax, zmin, zmax);

  for (i = 0; i < NDIM; i++)
    {
      for (j = 0; j < MDIM; j++)
	{
	  n = i * MDIM + j;
	  x[0] = wind_midx[i];
	  x[1] = 0;
	  x[2] = wind_midz[j];
	  plasmamain[n].rho = pl_rho (x);
	  plasmamain[n].t_r = geo.pl_t_r;
	  plasmamain[n].t_e = geo.pl_t_e;
	  plasmamain[n].w = geo.pl_w;
	  x[0] = wind_x[i];
	  x[1] = 0;
	  x[2] = wind_z[j];
	  pl_velocity (x, wmain[w[n].nwind].v);
	}
    }

  for (i = 0; i < NDIM - 1; i++)
    {
      for (j = 0; j < MDIM - 1; j++)
	{
	  n = i * MDIM + j;
	  wmain[w[n].nwind].vol = geo.pl_vol;
	}
    }

  /* Calculate dvds_ave */

  for (i = 0; i < NDIM; i++)
    {
      for (j = 0; j < MDIM; j++)
	{
	  n = i * MDIM + j;
	  x[0] = wind_midx[i];
	  x[1] = 0;
	  x[2] = wind_midz[j];
	  pl_velocity (x, vmid);
	  x[0] = wind_midx[i];
	  x[1] = 0;
	  x[2] = wind_z[j];
	  pl_velocity (x, vbot);
	  vsub (vmid, vbot, diff);
	  wmain[w[n].nwind].dvds_ave =
	    length (diff) / (wind_midz[j] - wind_z[j]);
	}
    }


  if (ion_abundances (&plasmamain[0], geo.ioniz_mode) != 0)
    {
      Error
	("pl_wind_define: Error code on return from ion_abundances %d\n",
	 ierr);
    }

  /* 07jul - ksl - As written all of the plasma conditions should be the same in 
   * all of the cells
   */

  for (n = 1; n < NDIM2; n++)
    {
      pl_copy_conditions (&plasmamain[0], &plasmamain[n]);

    }

  return (0);


}


double
pl_rho (x)
     double x[];
{
  if (geo.pl_nh <= 0)
    {
      Error ("pl_rho: Surprizing value of %e for nh\n", geo.pl_nh);
    }
  return (geo.pl_nh / rho2nh);
}

double
pl_velocity (x, v)
     double x[], v[];
{
  v[0] = 0.0;
  v[1] = 0.0;
  v[2] = geo.pl_vmax * x[2] / geo.pl_vol;
  return (0);

}


/* Routines below this point are generic */

/***********************************************************
                                       Space Telescope Science Institute

 Synopsis:
	int define_wind_grid(w) defines a grid of points for use with python and balance

Arguments:		
	WindPtr w 	The previously allocated wind structure
 	
Returns:
 
Description:	
	The routine creates a rectangular grid in x and z.  The grid points can be linearly
	spaced if geo.log_linear=1, or logarithmically spaced if geo.log_linear=0;
	The routine also defines the arrays wind_x,wind_z, wind_midx, wind_midz which
	are used by where_in_grid and for calculating values at the mid point of the
	grid points.
Notes:


History:
 	97jan      ksl	Coding on python began.
 	98dec	ksl	Coding on these routines begain.  It should replace portions of wind_define
 				and all of wind_complete
 
**************************************************************/
int
define_wind_grid (w, xmin, xmax, zmin, zmax)
     WindPtr w;
     double xmin, xmax, zmin, zmax;
{
  double dr, dz;
  double logr, logz, dlogr, zscale, dlogz;
  int i, j, n;

  /* Check the inputs */
  if ((xmax - xmin <= 0.0) || (zmax - zmin <= 0.0))
    {
      Error ("define_wind_grid: xmin %e xmax %e zmin %e zmax %e\n", xmin,
	     xmax, zmin, zmax);
      exit (0);
    }

  /* Set up constants for the two types of grids */

  if (geo.log_linear == 1)
    {				//linear array

      dr = (xmax - xmin) / (NDIM - 2);	/* This assumes that theta_min >= 0 */
      dz = (zmax - zmin) / (MDIM - 2);
    }
  else if (geo.log_linear == 0)
    {				/* Use logarithmic intervals */
      dlogr = (log10 (xmax) - log10 (xmin)) / (NDIM - 2);
      zscale = 1.e7;
      dlogz = (log10 (zmax / zscale + 1.0)) / (MDIM - 2);
    }
  else
    {
      Error ("define_wind_grid: Unallowable option for geo.log_linear %d\n",
	     geo.log_linear);
      exit (0);
    }

  /* First calculate parameters that are to be calculated at the edge of the grid cell.  This is
     mainly the positions and the velocity */
  for (i = 0; i < NDIM; i++)
    {
      for (j = 0; j < MDIM; j++)
	{
	  n = i * MDIM + j;

	  /*Define the grid points */
	  if (geo.log_linear == 1)
	    {			// linear intervals

	      w[n].x[0] = xmin + i * dr;	/* The first zone is at the inner radius of
						   the wind */
	      w[n].x[1] = 0;
	      w[n].x[2] = zmin + j * dz;
	    }
	  else
	    {			//logrithmic intevals

	      logr = i * dlogr;
	      w[n].x[0] = xmin * pow (10., logr);	/* The first zone is at the inner radius of
							   the wind */
	      w[n].x[1] = 0;
	      logz = j * dlogz;
	      w[n].x[2] = zscale * (pow (10., logz) - 1.0);
	    }


	}
    }

  /* Finally define some one-d vectors that make it easier to locate a photon in the wind given that we
     have adoped a "rectangular" grid of points.  Note that rectangular does not mean equally spaced. */

  for (i = 0; i < NDIM; i++)
    wind_x[i] = w[i * MDIM].x[0];
  for (j = 0; j < MDIM; j++)
    wind_z[j] = w[j].x[2];
  for (i = 0; i < NDIM - 1; i++)
    wind_midx[i] = 0.5 * (w[i * MDIM].x[0] + w[(i + 1) * MDIM].x[0]);
  for (j = 0; j < MDIM - 1; j++)
    wind_midz[j] = 0.5 * (w[j].x[2] + w[(j + 1)].x[2]);
  /* Add something plausible for the edges */
  wind_midx[NDIM - 1] = 2. * wind_x[NDIM - 1] - wind_midx[NDIM - 2];
  wind_midz[MDIM - 1] = 2. * wind_z[MDIM - 1] - wind_midz[MDIM - 2];
  return (0);

}


/* Copy physical conditions from one wind cell to another, that is copy everything but properties
 * that directly depend on the geometry
 *
 * 1108-ksl-What is copied here is the PlasmaPtr and the variables that are copied do not 
 * appear to be complete.  This is an error
 */


int
pl_copy_conditions (win, wout)
     PlasmaPtr win, wout;
{
  int n;
  wout->rho = win->rho;
  wout->t_r = win->t_r;
  wout->t_r_old = win->t_r_old;
  wout->t_e = win->t_e;
  wout->t_e_old = win->t_e_old;
  wout->ne = win->ne;
  for (n = 0; n < NIONS; n++)
    {
      wout->density[n] = win->density[n];
      wout->partition[n] = win->partition[n];
    }
  for (n = 0; n < NIONS; n++)
    {
      wout->ioniz[n] = win->ioniz[n];
      wout->recomb[n] = win->recomb[n];
      wout->lum_ion[n] = win->lum_ion[n];
      wout->heat_ion[n] = win->heat_ion[n];
    }
  wout->j = win->j;
  wout->ave_freq = wout->ave_freq;
  wout->lum = wout->lum;
  wout->lum_rad = wout->lum_rad;
  wout->lum_rad_old = wout->lum_rad_old;
  wout->lum_lines = wout->lum_lines;
  wout->lum_ff = wout->lum_ff;
  wout->lum_adiabatic = wout->lum_adiabatic;
  wout->lum_fb = wout->lum_fb;
  wout->lum_z = wout->lum_z;
  wout->heat_tot = wout->heat_tot;
  wout->heat_tot_old = wout->heat_tot_old;
  wout->heat_lines = wout->heat_lines;
  wout->heat_ff = wout->heat_ff;
  wout->heat_photo = wout->heat_photo;
  wout->heat_z = wout->heat_z;
  wout->dmo_dt[0] = wout->dmo_dt[0];
  wout->dmo_dt[1] = wout->dmo_dt[1];
  wout->dmo_dt[2] = wout->dmo_dt[0];
  wout->w = wout->w;

  return (0);

}
