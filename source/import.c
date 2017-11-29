#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "atomic.h"
#include "python.h"

/***********************************************************
    Space Telescope Science Institute

Synopsis:
    These are general purpose routines reading in model
    grids

Arguments:		

Returns:
 
Description:	

Notes:

    The routines contained here are basically steering
    routines. The real works is done in import_spherical,
    etc



History:
	17nov   ksl Began coding                           
**************************************************************/


# define LINELEN 512
# define NCELLS  512

/* Read in a model of in various coordiante systems, using the coord_type
 * to specify the type of model */

int
import_wind (ndom)
     int ndom;
{
  char filename[LINELEN];

  rdstr ("File.with.model2read", filename);

  if (zdom[ndom].coord_type == SPHERICAL)
    {
      import_1d (ndom, filename);
    }
  else if (zdom[ndom].coord_type == CYLIND)
    {
      import_cylindrical (ndom, filename);
    }
  else if (zdom[ndom].coord_type == RTHETA)
    {
      import_rtheta (ndom, filename);
    }
  else
    {
      Error
	("import_wind: Do not know how to import a model of coor_type %d\n",
	 zdom[ndom].coord_type);
      exit (0);
    }
  return (0);
}


/* Create the coordinate grids depending on the coord_type 
 *
 * */

int
import_make_grid (w, ndom)
     WindPtr w;
     int ndom;
{
  if (zdom[ndom].coord_type == SPHERICAL)
    {
      spherical_make_grid_import (w, ndom);
    }
  else if (zdom[ndom].coord_type == CYLIND)
    {
      cylindrical_make_grid_import (w, ndom);
    }
  else if (zdom[ndom].coord_type == RTHETA)
    {
      rtheta_make_grid_import (w, ndom);
    }
  else
    {
      Error
	("import_wind: Do not know how to import a model of coor_type %d\n",
	 zdom[ndom].coord_type);
      exit (0);
    }
  return (0);
}


/* Determine velocities for the various coord_types.  Note that
 * depending on the coordinate system we may alrady have calculated
 * velocities for wmain, and in that case these calls are used
 * for interpolation.
 */


double
import_velocity (ndom, x, v)
     int ndom;
     double *x, *v;
{
  double speed;
  if (zdom[ndom].coord_type == SPHERICAL)
    {
      speed = velocity_1d (ndom, x, v);
    }
  else if (zdom[ndom].coord_type == CYLIND)
    {
      speed = velocity_cylindrical (ndom, x, v);
    }
  else if (zdom[ndom].coord_type == RTHETA)
    {
      speed = velocity_rtheta (ndom, x, v);
    }
  else
    {
      Error
	("import_velocity: Do not know how to create velocities from model of coor_type %d\n",
	 zdom[ndom].coord_type);
      exit (0);
    }
  return (speed);
}


int
get_import_wind_params (ndom)
     int ndom;
{
  Log ("get_import_wind_params is currently a NOP\n");
  return (0);
}




/* Fill in plasma ptrs with densities.   
 *
 * For this we generally assume that the densities read in are 
 * given at the midpoints of the grid
 * 
 * */


double
import_rho (ndom, x)
     int ndom;
     double *x;
{
  double rho;
  if (zdom[ndom].coord_type == SPHERICAL)
    {
      rho = rho_1d (ndom, x);
    }
  else if (zdom[ndom].coord_type == CYLIND)
    {
      rho = rho_cylindrical (ndom, x);
    }
  else if (zdom[ndom].coord_type == RTHETA)
    {
      rho = rho_rtheta (ndom, x);
    }
  else
    {
      Error
	("import_rho:  Do not know how to create velocities from model of coor_type %d\n",
	 zdom[ndom].coord_type);
      exit (0);
    }
  return (rho);
}
