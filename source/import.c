
/***********************************************************/
/** @file  import.c
 * @author ksl
 * @date   May, 2018
 *
 * @brief   general purpose routines reading in model
 * grids
 *
 * 
 * The routines contained here are basically steering
 * routines. The real works is done in import_spherical,
 * etc
 *
 * For importing models, we first read in the data from a file.
 * We assume all of the data, positions, velocities and importanly
 * densities are given at the grid points of the imported model.
 *
 * We then map these models into the structures that Python uses.
 * Most of the mapping is one-to-one, but Python wants the densities
 * to be a the cell centers and not at the corners.
 ***********************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "atomic.h"
#include "python.h"



#define LINELEN 512
#define NCELLS  512

/* Read in a model of in various coordiante systems, using the coord_type
 * to specify the type of model */


/**********************************************************/
/** 
 * @brief      Import a gridded model
 *
 * @param [in] int  ndom   The domain for the model
 * @return   Always returns 0  
 *
 * @details
 *
 * This is a steering routine.  It reads the name of the file
 * to import and depending on the pre-established coordiante
 * system calls one of several coordiante system specific 
 * routines to actually read in the model
 *
 * ### Notes ###
 *
 **********************************************************/

int
import_wind (ndom)
     int ndom;
{
  char filename[LINELEN];

  rdstr ("Wind.model2import", filename);

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
    Error ("import_wind: Do not know how to import a model of coor_type %d\n", zdom[ndom].coord_type);
    exit (0);
  }
  return (0);
}


/* Create the coordinate grids depending on the coord_type 
 *
 * */


/**********************************************************/
/** 
 * @brief      Make the Python grid 
 *
 * @param [in] WindPtr  w  The entire wind structure
 * @param [in] int  ndom   The domain for the imported model
 * @return     Always returns 0
 *
 * @details
 * This is merely a steering routine for calling one
 * of the coordinate-system specific routines for creating
 * a grid from one of the imported models
 *
 * ### Notes ###
 * The fact that w is called by this routine is for constency.
 *
 **********************************************************/

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
    Error ("import_wind: Do not know how to import a model of coor_type %d\n", zdom[ndom].coord_type);
    exit (0);
  }
  return (0);
}


/* Determine velocities for the various coord_types.  Note that
 * depending on the coordinate system we may alrady have calculated
 * velocities for wmain, and in that case these calls are used
 * for interpolation.
 */



/**********************************************************/
/** 
 * @brief      Get the velocity at a position for an imported model
 *
 * @param [in] int  ndom   The domain in which the velocity is to be determined
 * @param [in] double *  x   The position where the velocity is to be determined
 * @param [out] double *  v   The velocity in cartesian coordiantes
 * @return     The speed       
 *
 * @details
 * The routine simply calls one of several coordinate-system specific routines
 * for obtaining the velocity of the wind when an imported model is involved.
 *
 * ### Notes ###
 * The routine is used to set up velocities in wmain
 *
 **********************************************************/

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
    Error ("import_velocity: Do not know how to create velocities from model of coor_type %d\n", zdom[ndom].coord_type);
    exit (0);
  }
  return (speed);
}



/**********************************************************/
/** 
 * @brief      get parameters associated with an imported wind model
 *
 * @param [in] int  ndom   The domain associated with an imported model
 * @return   Always returns 0  
 *
 * @details
 * At present this routine is simply a place holder and simply returns
 * to the call in setup_domains.  
 *
 * ### Notes ###
 *
 **********************************************************/

int
get_import_wind_params (ndom)
     int ndom;
{
  Log ("get_import_wind_params is currently a NOP\n");
  return (0);
}




/**********************************************************/
/** 
 * @brief      Get the density at an arbitray position in an imported 
 * model
 *
 * @param [in] int  ndom   The domain assosciated with an imported model
 * @param [in] double *  x   The position where we desire rho
 * @return     rho              
 *
 * @details
 * The routine simply calls coordinate system specific routines to get the
 * density for imported models
 *
 * ### Notes ###
 * The routine is used to map densities from an imported model 
 * where we assume that the density is given at the grid points.
 * In Python, we want map the grid points to the edges of wind cells,
 * but we expect the densities to be given at the centers of the cells.
 *
 **********************************************************/

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
    Error ("import_rho:  Do not know how to create velocities from model of coor_type %d\n", zdom[ndom].coord_type);
    exit (0);
  }
  return (rho);
}
