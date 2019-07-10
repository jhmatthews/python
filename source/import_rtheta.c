/***********************************************************/
/** @file  import_rtheta.c
 * @author ksl
 * @date   May, 2018
 *
 * @brief
 * These are general routines to read in a model that
 * is in polar or rtheta coordinates
 * grids
 * ###Notes###
 * There a various possibilities for how the
 * velocities could be entered.  One possibility
 * which is the way the zeus_python models work
 * is for the velocity to be given in spherical
 * polar coordinates.

 * However, internally, python uses xyz coordianes
 * for velocites (as measured in the xz plane),
 * and that is the model followed, here.  This also
 * makes these routines similar to those used
 * in imported cylindrical models.

 * This means that if the user provides a model
 * where velocities are in spherical polar coordinates
 * then one must translate them to the convention
 * here before the model is read in
 ***********************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "atomic.h"
#include "python.h"



#define LINELEN 512
#define NCELLS  512

struct
{
  int ndim, mdim, ncell;
  int i[NDIM_MAX * NDIM_MAX], j[NDIM_MAX * NDIM_MAX], inwind[NDIM_MAX * NDIM_MAX];
  double r[NDIM_MAX * NDIM_MAX], theta[NDIM_MAX * NDIM_MAX];
  double v_x[NDIM_MAX * NDIM_MAX], v_y[NDIM_MAX * NDIM_MAX], v_z[NDIM_MAX * NDIM_MAX];
  double rho[NDIM_MAX * NDIM_MAX], t[NDIM_MAX * NDIM_MAX];

  /* Note that the variables below look to be in xy, coordiantes
   * but they are really r, theta
   */

  double wind_x[NDIM_MAX], wind_z[NDIM_MAX], wind_midx[NDIM_MAX], wind_midz[NDIM_MAX];

} xx_rtheta;




/**********************************************************/
/**
 * @brief      Read the an arbitrary wind model in polar coordinates
 *
 * @param [in] int  ndom   The domain for the imported model
 * @param [in] char *  filename   The file containing the model to import
 * @return   Always returns 0
 *
 * @details
 * This routine reads the data into a set of arrays.  It's
 * only purpose is to read in the data
 *
 * ### Notes ###
 * The basic data we need to read in are
 *
 * icell, jcell, r theta inwind v_x v_y v_z  rho (and optionally T)
 *
 * where
 *
 * * r is the radial coordianate
 * * theta is the angular coordinate measured from the z axis
 * * v_x, v_y, and v_z is the velocity in cartesian coordinates
 *      as measured in the x,z plane
 * * rho is the density in cgs units
 * * inwind defines whether or not a particular cell is actually
 * in the wind
 *
 * We assume that all of the variables are centered, that is
 * we are not assuming that we are giving rho at the center of
 * a cell, but that r and v_r are at the edges of a cell.
 * This is somehing that would presumable be easy to change
 *
 **********************************************************/

int
import_rtheta (ndom, filename)
     int ndom;
     char *filename;
{
  FILE *fopen (), *fptr;
  char line[LINELEN];
  int n, icell, jcell, ncell, inwind;
  int jz, jx;
  double delta;
  double q1, q2, q3, q4, q5, q6, q7;


  Log ("Reading a model %s in polar (r,theta) coordinates \n", filename);


  if ((fptr = fopen (filename, "r")) == NULL)
  {
    Error ("import_rtheta: No such file\n");
    Exit (0);
  }


  ncell = 0;
  while (fgets (line, 512, fptr) != NULL)
  {
    n = sscanf (line, " %d %d %d %le %le %le %le %le %le %le", &icell, &jcell, &inwind, &q1, &q2, &q3, &q4, &q5, &q6, &q7);
    if (n < 4)
    {
      continue;
    }
    else
    {
      xx_rtheta.i[ncell] = icell;
      xx_rtheta.j[ncell] = jcell;
      xx_rtheta.inwind[ncell] = inwind;
      xx_rtheta.r[ncell] = q1;
      xx_rtheta.theta[ncell] = q2;
      xx_rtheta.v_x[ncell] = q3;
      xx_rtheta.v_y[ncell] = q4;
      xx_rtheta.v_z[ncell] = q5;
      xx_rtheta.rho[ncell] = q6;
      if (n > 9)
      {
        xx_rtheta.t[ncell] = q7;
      }
      else
      {
        xx_rtheta.t[ncell] = 10000.;
      }
      ncell++;

    }
  }

  /* Set and check the dimensions of the grids to be set up.
   * 
   * Note that some assumptions are built into the way the grid
   * is read in, most notably that the last cell to be read in
   * defines the dimensions of the entire grid.
   */

  zdom[ndom].ndim = xx_rtheta.ndim = icell + 1;
  zdom[ndom].mdim = xx_rtheta.mdim = jcell + 1;
  xx_rtheta.ncell = ncell;
  zdom[ndom].ndim2 = zdom[ndom].ndim * zdom[ndom].mdim;

  /* Check that the number of cells read in matches the number that was expected */

  if (ncell != xx_rtheta.ndim * xx_rtheta.mdim)
  {
    Error ("The dimensions of the imported grid seem wrong % d x %d != %d\n", xx_rtheta.ndim, xx_rtheta.mdim, xx_rtheta.ncell);
    exit (1);
  }


  jz = jx = 0;
  for (n = 0; n < xx_rtheta.ncell; n++)
  {
    if (xx_rtheta.i[n] == 0)
    {
      xx_rtheta.wind_z[jz] = xx_rtheta.theta[n];
      jz++;
    }
    if (xx_rtheta.j[n] == 0)
    {
      xx_rtheta.wind_x[jx] = xx_rtheta.r[n];
      jx++;
    }
  }

  for (n = 0; n < jz - 1; n++)
  {
    xx_rtheta.wind_midz[n] = 0.5 * (xx_rtheta.wind_z[n] + xx_rtheta.wind_z[n + 1]);
  }


  delta = (xx_rtheta.wind_z[n - 1] - xx_rtheta.wind_z[n - 2]);
  xx_rtheta.wind_midz[n] = xx_rtheta.wind_z[n - 1] + 0.5 * delta;

  for (n = 0; n < jx - 1; n++)
  {
    xx_rtheta.wind_midx[n] = 0.5 * (xx_rtheta.wind_x[n] + xx_rtheta.wind_x[n + 1]);
  }


  delta = (xx_rtheta.wind_x[n - 1] - xx_rtheta.wind_x[n - 2]);
  xx_rtheta.wind_midx[n] = xx_rtheta.wind_x[n - 1] + 0.5 * delta;

  return (0);
}




/**********************************************************/
/**
 * @brief      Use the imported data to initialize various
 * portions of the Wind and Domain structures
 *
 * @param [out] WindPtr  w   The wind structure
 * @param [in] int  ndom   The domain for the imported model
 * @return     Always returns 0
 *
 * @details
 * This routine initializes the portions of the wind structure
 * using the imported model, specirically those portions having
 * to do with positions, and velocities.
 *
 * ### Notes ###
 * The routine also initials wind_x and wind_z in the domain
 * structure, ans sets up wind_cones and other boundaries
 * intended to bound the wind.
 *
 **********************************************************/

int
rtheta_make_grid_import (w, ndom)
     WindPtr w;
     int ndom;
{
  int n, nn;
  double theta;
  double rho_max, rho_min, r_inner, r_outer, rmin, rmax;
  double zmin, zmax;

  /* As in the case of other models we assume that the grid has been
   * read in correctly and so now that the WindPtrs have been generated
   * we can put a lot of information in directly
   */

  for (n = 0; n < xx_rtheta.ncell; n++)
  {
    wind_ij_to_n (ndom, xx_rtheta.i[n], xx_rtheta.j[n], &nn);
    w[nn].r = xx_rtheta.r[n];
    w[nn].theta = theta = xx_rtheta.theta[n];

    theta /= RADIAN;

    w[nn].x[0] = w[nn].r * sin (theta);
    w[nn].x[1] = 0;
    w[nn].x[2] = w[nn].r * cos (theta);
    w[nn].v[0] = xx_rtheta.v_x[n];
    w[nn].v[1] = xx_rtheta.v_y[n];
    w[nn].v[2] = xx_rtheta.v_z[n];
    w[nn].inwind = xx_rtheta.inwind[n];

    w[nn].thetacen = xx_rtheta.wind_midz[xx_rtheta.j[n]];
    theta = w[nn].thetacen / RADIAN;

    w[nn].rcen = xx_rtheta.wind_midx[xx_rtheta.i[n]];


    w[nn].xcen[0] = w[nn].rcen * sin (theta);
    w[nn].xcen[1] = 0;
    w[nn].xcen[2] = w[nn].rcen * cos (theta);

    /* JM 1711 -- copy across the inwind variable to the wind pointer */
    w[nn].inwind = xx_rtheta.inwind[n];


    /* 1812 - ksl - For imported models, one is either in the wind or not. But we need
     * to make sure the rest of the code knows that this cell is to be ignored in
     * this case. Adapted from the code in import_cylindrical */
    if (w[nn].inwind == W_NOT_INWIND || w[nn].inwind == W_PART_INWIND)
      w[nn].inwind = W_IGNORE;
  }

  /* Now add information used in zdom */

  for (n = 0; n < zdom[ndom].ndim; n++)
  {
    zdom[ndom].wind_x[n] = xx_rtheta.wind_x[n];
  }



  for (n = 0; n < zdom[ndom].mdim; n++)
  {
    zdom[ndom].wind_z[n] = xx_rtheta.wind_z[n];
  }

  /* Now set up wind boundaries so they are harmless.
   * Note that the grid goes from near the pole
   * to the equator
   */


  rmax = rho_max = zmax = 0;
  rmin = rho_min = zmin = VERY_BIG;
  for (n = 0; n < xx_rtheta.ncell; n++)
  {
    wind_ij_to_n (ndom, xx_rtheta.i[n], xx_rtheta.j[n], &nn);

    r_inner = length (w[nn].x);

    r_outer = length (w[nn + xx_rtheta.mdim].x);


    if (w[nn].inwind >= 0)
    {
      if (w[nn + xx_rtheta.mdim].x[0] > rho_max)
      {
        rho_max = w[nn + xx_rtheta.mdim].x[0];
      }
      if (w[nn - 1].x[2] > zmax)
      {
        zmax = w[nn - 1].x[2];
      }
      if (w[nn].x[2] < zmin && w[nn].x[2] > 0)
      {
        zmin = w[nn].x[2];
      }
      if (r_outer > rmax)
      {
        rmax = r_outer;
      }
      if (rho_min > w[nn].x[0])
      {
        rho_min = w[nn].x[0];
      }
      if (rmin > r_inner)
      {
        rmin = r_inner;
      }
    }
  }

  Log ("Imported:    rmin    rmax  %e %e\n", rmin, rmax);
  Log ("Imported:    zmin    zmax  %e %e\n", zmin, zmax);
  Log ("Imported: rho_min rho_max  %e %e\n", rho_min, rho_max);


  zdom[ndom].wind_rho_min = zdom[ndom].rho_min = rho_min;
  zdom[ndom].wind_rho_max = zdom[ndom].rho_max = rho_max;
  zdom[ndom].zmax = zmax;

  zdom[ndom].rmax = rmax;
  zdom[ndom].rmin = rmin;
  zdom[ndom].wind_thetamin = zdom[ndom].wind_thetamax = 0.;

  /* The next line is necessary for calculating distances in a cell in rthota coordiatnes */

  rtheta_make_cones (ndom, w);

  return (0);
}


/* The next section calculates velocites.  We follow the hydro approach of
 * getting those velocities from the original grid.  This is really only
 * used for setting up the grid
 *
 * The code here is identical to that in velocity_cylindrical, which suggests
 * that it could be used for any regular 2d grid
 */


/**********************************************************/
/**
 * @brief      The velcity at any position in an imported
 * rtheat model
 *
 * @param [in] int  ndom   The domain for the imported model
 * @param [in] double *  x   A position (3 vector)
 * @param [out] double *  v   The calcuated velocity
 * @return     The speed at x
 *
 * @details
 * This routine interpolates on the values read in for the
 * inported model to give one a velocity
 *
 *
 * ### Notes ###
 *  In practice this routine is only used to initallize v in
 *  wind structure.  This is consistent with the way velocities
 *  are treated throughout Python.
 *
 **********************************************************/

double
velocity_rtheta (ndom, x, v)
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








/**********************************************************/
/**
 * @brief      Get the density for an imported rtheta model at x
 *
 * @param [in] int  ndom   The domain for the imported model
 * @param [in] double *  x   A postion
 * @return     The density in cgs units is returned
 *
 * @details
 * This routine finds rho from the imported model
 * at a position x.  The routine does not interpolate rho, but
 * simply locates the cell associated with x
 *
 * ### Notes ###
 * This routine is really only used to intialize rho in the
 * Plasma structure.  In reality, once the Plasma structure is
 * initialized we always interpolate within the plasma structure
 * and do not access the original data
 *
 *
 **********************************************************/

double
rho_rtheta (ndom, x)
     int ndom;
     double *x;
{
  double rho = 0;
  double r, z;
  int i, j, n;
  double ctheta, angle;

  r = length (x);

  z = fabs (x[2]);
  ctheta = z / r;
  angle = acos (ctheta) * RADIAN;

  i = 0;
  while (angle > xx_rtheta.wind_z[i] && i < xx_rtheta.mdim - 1)
  {
    i++;
  }
  j = 0;
  while (r > xx_rtheta.wind_x[j] && j < xx_rtheta.ndim - 1)
  {
    j++;
  }

  n = j * xx_rtheta.mdim + i;

  rho = xx_rtheta.rho[n];

  return (rho);
}
