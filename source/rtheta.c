
/***********************************************************/
/** @file  rtheta.c
 * @author ksl
 * @date   March, 2018
 *
 * @brief  Generic routines for use with an azimuthally-symmetric
 * rtheta coordinate system.
 *
 ***********************************************************/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "atomic.h"
#include "python.h"





/**********************************************************/
/**
 * @brief      calculates the distance to the far
 *         boundary of the cell in which the photon bundle resides.
 *
 * @param [in] int  ndom   The domain in which the photon bundle is though to exist
 * @param [in, out] PhotPtr  p   Photon pointer
 * @return     Distance to the far boundary of the cell in which the photon
 * 	currently resides.  Negative numbers (and zero) should be
 * 	regarded as errors.
 *
 * @details
 * ??? DESCRIPTION ???
 *
 * ### Notes ###
 * ??? NOTES ???
 *
 **********************************************************/

double
rtheta_ds_in_cell (ndom, p)
     int ndom;
     PhotPtr p;


{

  int n, ix, iz;
  double s, smax;



  /* Check that that the photon is in the domain it is supposed to be
   * in.  */

  if ((p->grid = n = where_in_grid (ndom, p->x)) < 0)
  {
    Error ("rtheta_ds_in_cell: Photon not in grid when routine entered\n");
    return (n);                 /* Photon was not in wind */
  }

  wind_n_to_ij (ndom, n, &ix, &iz);     /*Convert the index n to two dimensions */


  /* Set up the quadratic equations in the radial  direction */

  smax = ds_to_sphere (zdom[ndom].wind_x[ix], p);
  s = ds_to_sphere (zdom[ndom].wind_x[ix + 1], p);
  if (s < smax)
  {
    smax = s;
  }

  /* At this point we have found how far the photon can travel in r in its
     current direction.  Now we must worry about motion in the theta direction  */

  s = ds_to_cone (&zdom[ndom].cones_rtheta[iz], p);
  if (s < smax)
  {
    smax = s;
  }

  s = ds_to_cone (&zdom[ndom].cones_rtheta[iz + 1], p);
  if (s < smax)
  {
    smax = s;
  }

  if (smax <= 0)
  {
    Error ("rtheta: ds_in_cell %f\n", smax);
  }
  return (smax);
}




/**********************************************************/
/**
 * @brief      Set up an rtheta grid
 *
 * @param [in] WindPtr  w   The structure which defines entire wind
 * @param [in] int  ndom   The current domain
 * @return     Always returns 0
 *
 * @details
 *
 *
 * 	For spherical polar, components we have defined the grid so that
 * 	i corresponds to a radial distance, and j corresponds to an angle,
 * 	e.g r theta.  increasing i corresponds to increasing distance,
 * 	and increasing theta corresponds to increasing angle measured from
 * 	the z axis. (This it should be fairly easy to implement true
 * 	spherical coordinates in the future.
 *
 * 	There are two basic options:
 *
 * 	The r-sapcing can be linear, or logarithmic.
 *
 * 	The theta spacing is always linear
 *
 *
 * ### Notes ###
 *
 * The parameters for describing the grid to be created must
 * have to have been intialized and are part of the domain structure
 *
 *
 **********************************************************/

int
rtheta_make_grid (w, ndom)
     WindPtr w;
     int ndom;
{
  double dr, theta, thetacen, dtheta, dlogr;
  int i, j, n;
  int mdim, ndim;

  /* Set zmax if this has not alrady been done.
   */

  if (zdom[ndom].zmax == 0)
  {
    /* Check if zmax has been set, and if not set it to rmax */
    zdom[ndom].zmax = zdom[ndom].rmax;
  }

  ndim = zdom[ndom].ndim;
  mdim = zdom[ndom].mdim;

  /* In order to interpolate the velocity (and other) vectors out to zdom[ndom].rmax, we need
     to define the wind at least one grid cell outside the region in which we want photons
     to propagate. */

/* Next two lines for linear intervals */
  dtheta = 90. / (mdim - 3);


  /* First calculate parameters that are to be calculated at the edge of the grid cell.  This is
     mainly the positions and the velocity */

  for (i = 0; i < ndim; i++)
  {
    for (j = 0; j < mdim; j++)
    {
      wind_ij_to_n (ndom, i, j, &n);


      /*Define the grid points */
      if (zdom[ndom].log_linear == 1)
      {                         // linear intervals

        dr = (zdom[ndom].rmax - geo.rstar) / (ndim - 3);
        w[n].r = geo.rstar + i * dr;
        w[n].rcen = w[n].r + 0.5 * dr;
      }
      else
      {                         //logarithmic intervals

        dlogr = (log10 (zdom[ndom].rmax / geo.rstar)) / (mdim - 3);
        w[n].r = geo.rstar * pow (10., dlogr * (i - 1));
        w[n].rcen = 0.5 * geo.rstar * (pow (10., dlogr * (i)) + pow (10., dlogr * (i - 1)));
      }

      /* Only the radial distance can be logarithmic */

      theta = w[n].theta = dtheta * j;
      thetacen = w[n].thetacen = w[n].theta + 0.5 * dtheta;
      if (theta > 90.)
      {
        theta = 90.;
      }
      if (thetacen > 90.)
      {
        thetacen = 90.;
      }


      /* Now calculate the positions of these points in the xz plane */
      theta /= RADIAN;
      thetacen /= RADIAN;
      w[n].x[1] = w[n].xcen[1] = 0.0;

      w[n].x[0] = w[n].r * sin (theta);
      w[n].x[2] = w[n].r * cos (theta);

      w[n].xcen[0] = w[n].rcen * sin (thetacen);
      w[n].xcen[2] = w[n].rcen * cos (thetacen);

    }
  }

  /* To define the edges of cells, one needs to define a series of
   * cone structures which are used by ds_in_cell to find the distance
   * a photon can travel
   */

  rtheta_make_cones (ndom, w);

  return (0);
}


/**********************************************************/
/**
 * @brief      defines the wind cones that are needed to calculate ds in a cell
 *
 * @param [in, out] int  ndom   The domain to populate with the cone structures
 * @param [in, out] WindPtr  w   The structure which defines the wind in Python
 * @return     Always returns 0
 *
 * @details
 *
 * This subroutine defines cones that define the edges (in the theta directio)
 * of wind cells.
 *
 * ### Notes ###
 *
 * Comment:  It is not entirely clear why we don't use fixed arrays for the cones
 * in the domains. The data volume is small, and this would simplify windsave
 * and windread.
 *
 **********************************************************/

int
rtheta_make_cones (ndom, w)
     int ndom;
     WindPtr w;
{
  int n;
  int mdim;

  mdim = zdom[ndom].mdim;


  zdom[ndom].cones_rtheta = (ConePtr) calloc (sizeof (cone_dummy), mdim);
  if (zdom[ndom].cones_rtheta == NULL)
  {
    Error ("rtheta_make_gid: There is a problem in allocating memory for the cones structure\n");
    exit (0);

  }


  for (n = 0; n < mdim; n++)
  {
    zdom[ndom].cones_rtheta[n].z = 0.0;
    zdom[ndom].cones_rtheta[n].dzdr = 1. / tan (w[n].theta / RADIAN);
  }


  return (0);
}






/**********************************************************/
/**
 * @brief      Complete the creation of an rtheta wind comain by populating certain arrays used for interpolation
 *
 * @param [in] int  ndom   The domain containing this wind component
 * @param [in] WindPtr  w   The entire wind
 * @return     Always returns 0
 *
 * @details
 *
 * ### Notes ###
 *
 * This routine sets up a few arrays in domains that have rtheta
 * coordiantes.  These arrays are used to determine where in
 * a cell one is, and how far it is to the edge of the cell
 * for photon bundles.
 *
 * The routine is called when the domain is set up, and called
 * again whenever the windsave file is read in.
 *
 **********************************************************/

int
rtheta_wind_complete (ndom, w)
     int ndom;
     WindPtr w;
{
  int i, j;
  int ndim, mdim, nstart;

  ndim = zdom[ndom].ndim;
  mdim = zdom[ndom].mdim;
  nstart = zdom[ndom].nstart;



  /* Finally define some one-d vectors that make it easier to locate a photon in the wind given that we
     have adoped a "rectangular" grid of points.  Note that rectangular does not mean equally spaced. */

  for (i = 0; i < ndim; i++)
  {
    zdom[ndom].wind_x[i] = w[nstart + i * mdim].r;
  }
  for (j = 0; j < mdim; j++)
    zdom[ndom].wind_z[j] = w[nstart + j].theta;

  for (i = 0; i < ndim - 1; i++)
    zdom[ndom].wind_midx[i] = w[nstart + i * mdim].rcen;

  for (j = 0; j < mdim - 1; j++)
    zdom[ndom].wind_midz[j] = w[nstart + j].thetacen;

  /* Add something plausible for the edges */

  zdom[ndom].wind_midx[ndim - 1] = 2. * zdom[ndom].wind_x[ndim - 1] - zdom[ndom].wind_midx[ndim - 2];
  zdom[ndom].wind_midz[mdim - 1] = 2. * zdom[ndom].wind_z[mdim - 1] - zdom[ndom].wind_midz[mdim - 2];


  /* Finally, in order to complete the r-theta wind we need to make a set of wind-cones. This is
     to allow use to use the cones routines to work out if photonds leave a cell in the theta direction */

  rtheta_make_cones (ndom, w);


  return (0);
}


//OLD /***********************************************************
//OLD                                        Space Telescope Science Institute
//OLD
//OLD  Synopsis:
//OLD   rtheta_volume(w) calculates the wind volume of a cylindrical cell
//OLD   allowing for the fact that some cells
//OLD
//OLD  Arguments:
//OLD   WindPtr w;    the entire wind
//OLD  Returns:
//OLD
//OLD  Description:
//OLD   This is a brute-froce integration of the volume
//OLD
//OLD   ksl--04aug--some of the machinations regarding what to to at the
//OLD   edge of the wind seem bizarre, like a substiture for figuring out
//OLD   what one actually should be doing.  However, volume > 0 is used
//OLD   for making certain choices in the existing code, and so one does
//OLD   need to be careful.  The existing code does not accurately calcuate
//OLD   the volume in the sense that it does not weight accroding to rho!
//OLD
//OLD  Notes:
//OLD   Where_in grid does not tell you whether the photon is in the wind or not.
//OLD  History:
//OLD   04aug   ksl     52a -- Moved from wind2d
//OLD   05apr   ksl     55d -- Modified to include the determination of whether
//OLD                   a cell was completely in the wind or not.  This
//OLD                   functionality had been in define_wind.
//OLD   06nov   ksl     58b -- Minor modification to use defined variables
//OLD                   W_ALL_INWIND, etc. instead of hardcoded vlues
//OLD   03apr   ksl     68c -- Added robust check of whether a cell was
//OLD                   in the wind or not to speed up the volume
//OLD                   calculation (and allow one to use a high resolution
//OLD                   for the numerical integration when a cell is
//OLD                   partially in the wind
//OLD   11aug   ksl     70b -- Added possibility of finding the volumes
//OLD                   in multiple components. See python.h for discussion
//OLD                   of the relationship betwedn PART and LLL
//OLD
//OLD
//OLD **************************************************************/
#define RESOLUTION   1000




/**********************************************************/
/**
 * @brief      rtheta_volume(w) calculates the wind volume of a cell in rtheta coordinate sysem
 * 	allowing for the fact that some cells may not be in the active wind region
 *
 * @param [in] int  ndom   The domain of interest
 * @param [in, out] WindPtr  w   the entire wind
 * @return     Always returns 0
 *
 * @details
 *
 * 	ksl--04aug--some of the machinations regarding what to to at the
 * 	edge of the wind seem bizarre, like a substiture for figuring out
 * 	what one actually should be doing.  However, volume > 0 is used
 * 	for making certain choices in the existing code, and so one does
 * 	need to be careful.  The existing code does not accurately calcuate
 * 	the volume in the sense that it does not weight accroding to rho!
 *
 * ### Notes ###
 *
 *
 * This is a brute-force integration of the volume
 *
 * In calculating volumes we allow for the fact that a cell
 * exists above and below the plane of the disk
 *
 *
 *
 **********************************************************/

int
rtheta_volumes (ndom, w)
     int ndom;
     WindPtr w;
{
  int i, j, n;
  int jj, kk;
  double fraction;
  double num, denom;
  double r, theta;
  double dr, dtheta, x[3];
  double rmin, rmax, thetamin, thetamax;
  double cell_volume;
  int n_inwind;
  int ndim, mdim, ndomain;

  ndim = zdom[ndom].ndim;
  mdim = zdom[ndom].mdim;

  for (i = 0; i < ndim; i++)
  {
    for (j = 0; j < mdim; j++)
    {
      wind_ij_to_n (ndom, i, j, &n);

      rmin = zdom[ndom].wind_x[i];
      rmax = zdom[ndom].wind_x[i + 1];
      thetamin = zdom[ndom].wind_z[j] / RADIAN;
      thetamax = zdom[ndom].wind_z[j + 1] / RADIAN;

      /*leading factor of 2 allows for volume above and below plane */
      cell_volume = 2. * 2. / 3. * PI * (rmax * rmax * rmax - rmin * rmin * rmin) * (cos (thetamin) - cos (thetamax));


      if (w[n].inwind == W_NOT_INWIND)
      {
        w[n].vol = 0.0;
      }

      else if (w[n].inwind == W_ALL_INWIND)
      {
        w[n].vol = cell_volume;
      }


      else if (w[n].inwind == W_NOT_ASSIGNED)
      {
        if (zdom[ndom].wind_type == IMPORT)
        {
          Error ("rtheta_volumes: Shouldn't be redefining inwind in cylind_volumes with imported model.\n");
          exit (0);
        }


        rmin = zdom[ndom].wind_x[i];
        rmax = zdom[ndom].wind_x[i + 1];
        thetamin = zdom[ndom].wind_z[j] / RADIAN;
        thetamax = zdom[ndom].wind_z[j + 1] / RADIAN;

        w[n].vol = cell_volume;

        n_inwind = rtheta_is_cell_in_wind (n);
        if (n_inwind == W_NOT_INWIND)
        {
          fraction = 0.0;       /* Force outside edge volues to zero */
          jj = 0;
          kk = RESOLUTION * RESOLUTION;
        }
        else if (n_inwind == W_ALL_INWIND)
        {
          fraction = 1.0;       /* Force outside edge volues to zero */
          jj = kk = RESOLUTION * RESOLUTION;
        }


        else
        {                       /* The grid cell is PARTIALLY in the wind */
          num = denom = 0;
          jj = kk = 0;
          dr = (rmax - rmin) / RESOLUTION;
          dtheta = (thetamax - thetamin) / RESOLUTION;
          for (r = rmin + dr / 2; r < rmax; r += dr)
          {
            for (theta = thetamin + dtheta / 2; theta < thetamax; theta += dtheta)
            {
              denom += r * r * sin (theta);;
              kk++;
              x[0] = r * sin (theta);
              x[1] = 0;
              x[2] = r * cos (theta);;
              if (where_in_wind (x, &ndomain) == W_ALL_INWIND)
              {
                num += r * r * sin (theta);     /* 0 implies in wind */
                jj++;
              }
            }
          }
          fraction = num / denom;
        }

        /* OK now assign inwind value and final volumes */

        if (jj == 0)
        {
          w[n].inwind = W_NOT_INWIND;   // The cell is not in the wind
          w[n].vol = 0.0;
        }
        else if (jj == kk)
          w[n].inwind = W_ALL_INWIND;   // The cell is completely in the wind
        else
        {
          w[n].inwind = W_PART_INWIND;  //The cell is partially in the wind
          w[n].vol *= fraction;
        }


      }
//OLD     /* JM/ksl 1711 -- the following two if statements are for if the inwind values are
//OLD        already assigned, for example by an imported model */
//OLD     /* need to zero volumes for cells not in the wind */
//OLD
//OLD     /* the logic is pretty awkward.  It seems as if we already knew a cell was in or ou
//OLD      * of the wind, then we should put these at the top rather than at the bottom of the
//OLD      * loop, here and in cylindrical_volummes
//OLD      */
//OLD
//OLD     else if (w[n].inwind == W_NOT_INWIND)
//OLD       {
//OLD         w[n].vol = 0.0;
//OLD       }
//OLD
//OLD     else if (w[n].inwind == W_ALL_INWIND)
//OLD       {
//OLD         w[n].vol = cell_volume;
//OLD       }

    }
  }

  return (0);
}



//OLD /***********************************************************
//OLD                      Space Telescope Science Institute
//OLD
//OLD  Synopsis:
//OLD   rtheta_where_in_grid locates the grid position of the vector,
//OLD   when one is using rtheta coordinates.
//OLD
//OLD  Arguments:
//OLD   double x[];
//OLD  Returns:
//OLD   where_in_grid normally  returns the cell number associated with
//OLD           a position.  If the position is in the grid this will be a positive
//OLD           integer < NDIM*MDIM.
//OLD   x is inside the grid        -1
//OLD   x is outside the grid       -2
//OLD  Description:
//OLD
//OLD
//OLD  Notes:
//OLD   Where_in grid does not tell you whether the x is in the wind or not.
//OLD
//OLD   What one means by inside or outside the grid may well be different
//OLD   for different coordinate systems.
//OLD
//OLD
//OLD
//OLD  History:
//OLD   04aug   ksl     52a -- Adapted from the same routine for cylindrical
//OLD                   systems.
//OLD           13sep   nsh     76b -- Changed calls to fraction to take account of
//OLD                   new modes.
//OLD   15aug   ksl     Added domains.
//OLD
//OLD **************************************************************/




/**********************************************************/
/**
 * @brief      locates the grid position of the vector,
 * 	when one is using rtheta coordinates.
 *
 * @param [in, out] int  ndom   The domain of interest
 * @param [in, out] double  x[]   A three-vector defining a position
 * @return     Returns the cell number associated with
 *  		a position. If x is inside the grid, the routine
 *  		returns -1, if outside -2
 *
 * @details
 * ??? DESCRIPTION ???
 *
 * ### Notes ###
 * Where_in grid does not tell you whether the x is in the wind or not.
 *
 * 	What one means by inside or outside the grid may well be different
 * 	for different coordinate systems.
 *
 **********************************************************/

int
rtheta_where_in_grid (ndom, x)
     int ndom;
     double x[];
{
  int i, j, n;
  double r, theta;
  double f;
  int ndim, mdim;

  ndim = zdom[ndom].ndim;
  mdim = zdom[ndom].mdim;

  r = length (x);
  theta = acos ((fabs (x[2] / r))) * RADIAN;

  /* Check to see if x is outside the region of the calculation */

  if (r > zdom[ndom].wind_x[ndim - 1])  /* Fixed version */
  {
    return (-2);                /* x is outside grid */
  }
  else if (r < zdom[ndom].wind_x[0])
  {
    return (-1);                /*x is inside grid */
  }

  /* Locate the position in i and j */

  fraction (r, zdom[ndom].wind_x, ndim, &i, &f, 0);
  fraction (theta, zdom[ndom].wind_z, mdim, &j, &f, 0);

  /* Convert i,j back to n */

  wind_ij_to_n (ndom, i, j, &n);

  return (n);
}

//OLD /***********************************************************
//OLD                      Space Telescope Science Institute
//OLD
//OLD  Synopsis:
//OLD   rtheta_get_random_location
//OLD
//OLD  Arguments:
//OLD   int n -- Cell in which random poition is to be generated
//OLD  Returns:
//OLD   double x -- the position
//OLD  Description:
//OLD
//OLD
//OLD  Notes:
//OLD
//OLD
//OLD
//OLD  History:
//OLD   04aug   ksl     52a -- Moved from where_in_wind as incorporated
//OLD                   multiple coordinate systems
//OLD   11aug   ksl     70b - Modifications to incoporate multiple
//OLD                   components
//OLD   15aug   ksl     Allow for mulitple domains
//OLD
//OLD **************************************************************/


/**********************************************************/
/**
 * @brief
 *
 * @param [in, out] int  n   -- Cell in which random poition is to be generated
 * @param [in, out] double  x[]   ???
 * @return     double x -- the position
 *
 * @details
 * ??? DESCRIPTION ???
 *
 * ### Notes ###
 * ??? NOTES ???
 *
 **********************************************************/

int
rtheta_get_random_location (n, x)
     int n;                     // Wind cell in which to create position
     double x[];                // Returned position
{
  int i, j;
  int inwind;
  double r, rmin, rmax, sthetamin, sthetamax;
  double theta, phi;
  double zz;
  int ndom, ndomain;

  ndom = wmain[n].ndom;
  wind_n_to_ij (ndom, n, &i, &j);

  rmin = zdom[ndom].wind_x[i];
  rmax = zdom[ndom].wind_x[i + 1];
  sthetamin = sin (zdom[ndom].wind_z[j] / RADIAN);
  sthetamax = sin (zdom[ndom].wind_z[j + 1] / RADIAN);

  /* Generate a position which is both in the cell and in the wind */

  inwind = W_NOT_INWIND;
  while (inwind != W_ALL_INWIND)
  {
    r = sqrt (rmin * rmin +
//            (rand () / (MAXRAND - 0.5)) * (rmax * rmax - rmin * rmin)); DONE
              random_number (0.0, 1.0) * (rmax * rmax - rmin * rmin));

//      theta = asin (sthetamin + (rand () / MAXRAND) * (sthetamax - sthetamin)); DONE
    theta = asin (sthetamin + random_number (0.0, 1.0) * (sthetamax - sthetamin));


//      phi = 2. * PI * (rand () / MAXRAND); DONE
    phi = 2. * PI * random_number (0.0, 1.0);

/* Project from r, theta phi to x y z  */

    x[0] = r * cos (phi) * sin (theta);
    x[1] = r * sin (phi) * sin (theta);
    x[2] = r * cos (theta);
    inwind = where_in_wind (x, &ndomain);       /* Some photons will not be in the wind
                                                   because the boundaries of the wind split the grid cell */
  }

//  zz = rand () / MAXRAND - 0.5;       //positions above are all at +z distances DONE
  zz = random_number (-1.0, 1.0);       //positions above are all at +z distances


  if (zz < 0)
    x[2] *= -1;                 /* The photon is in the bottom half of the wind */

  return (inwind);

}



//OLD /***********************************************************
//OLD                      Space Telescope Science Institute
//OLD
//OLD  Synopsis:
//OLD   rtheta_extend_density  extends the density to
//OLD   regions just outside the wind regiions so that
//OLD   extrapolations of density can be made there
//OLD
//OLD  Arguments:
//OLD  Returns:
//OLD  Description:
//OLD
//OLD
//OLD  Notes:
//OLD       Now we need to updated the densities immediately outside the wind so that the density interpolation in resonate will work.
//OLD      In this case all we have done is to copy the densities from the cell which is just in the wind (as one goes outward) to the
//OLD      cell that is just inside (or outside) the wind.
//OLD
//OLD      SS asked whether we should also be extending the wind for other parameters, especially ne.  At present we do not interpolate
//OLD      on ne so this is not necessary.  If we did do that it would be required.
//OLD
//OLD      In cylindrical coordinates, the fast dimension is z; grid positions increase up in z, and then out in r.
//OLD      In spperical polar coordinates, the fast dimension is theta; the grid increases in theta (measured)
//OLD      from the z axis), and then in r.
//OLD      In spherical coordinates, the grid increases as one might expect in r..
//OLD      *
//OLD      06may      ksl     57+ -- Trying simply to use mappings to extend the grid in density space.  The
//OLD      danger is that we will do other things, e.g update some of the wrong parameters.
//OLD
//OLD
//OLD
//OLD
//OLD  History:
//OLD   05apr   ksl     56 -- Moved functionality from wind updates
//OLD
//OLD **************************************************************/



/**********************************************************/
/**
 * @brief      extends the density to
 * 	regions just outside the wind regiions so that
 * 	extrapolations of density can be made there
 *
 * @param [in, out] int  ndom   ???
 * @param [in, out] WindPtr  w   ???
 * @return     ??? RETURNS ???
 *
 * @details
 * ??? DESCRIPTION ???
 *
 * ### Notes ###
 * Now we need to updated the densities immediately outside the wind so that the density interpolation in resonate will work.
 *      In this case all we have done is to copy the densities from the cell which is just in the wind (as one goes outward) to the
 *      cell that is just inside (or outside) the wind.
 *
 *      SS asked whether we should also be extending the wind for other parameters, especially ne.  At present we do not interpolate
 *      on ne so this is not necessary.  If we did do that it would be required.
 *
 *      In cylindrical coordinates, the fast dimension is z; grid positions increase up in z, and then out in r.
 *      In spperical polar coordinates, the fast dimension is theta; the grid increases in theta (measured)
 *      from the z axis), and then in r.
 *      In spherical coordinates, the grid increases as one might expect in r..
 *      *
 *      06may      ksl     57+ -- Trying simply to use mappings to extend the grid in density space.  The
 *      danger is that we will do other things, e.g update some of the wrong parameters.
 *
 **********************************************************/

int
rtheta_extend_density (ndom, w)
     int ndom;
     WindPtr w;
{
  int i, j, n, m;
  int ndim, mdim;

  ndim = zdom[ndom].ndim;
  mdim = zdom[ndom].mdim;

  for (i = 0; i < ndim - 1; i++)
  {
    for (j = 0; j < mdim - 1; j++)
    {
      wind_ij_to_n (ndom, i, j, &n);
      if (w[n].vol == 0)

      {                         /*Then this grid point is not in the wind */

        wind_ij_to_n (ndom, i + 1, j, &m);
        if (w[m].vol > 0)
        {                       /*Then the windcell in the +x direction is in the wind and
                                   we can copy the densities to the grid cell n  */
          w[n].nplasma = w[m].nplasma;

        }
        else if (i > 0)
        {
          wind_ij_to_n (ndom, i - 1, j, &m);
          if (w[m].vol > 0)
          {                     /*Then the grid cell in the -x direction is in the wind and
                                   we can copy the densities to the grid cell n */
            w[n].nplasma = w[m].nplasma;

          }
        }
      }
    }
  }

  return (0);

}




/**********************************************************/
/**
 * @brief      Check whether a cell in an rtheta coordiante system
 * is in the wind
 *
 * @param [in] int  n   The a cell number
 * @return     An integer indicate whether the cell is in the wind, pratially in the
 * wind, or not in the wind at all.
 *
 * @details
 * This routine performes is a robust check of whether a cell is
 * in the wind or not.
 *
 * ### Notes ###
 *
 * The routine was created to speed up the evaluation of the volumes for the wind.  It
 * checks each of the four boundaries of the wind to see whether any portions
 * of these are in the wind
 *
 * Comment:  The routine is only called by rtheta_volumes, and so we
 * already know that this cell is part of an rtheta coordiante
 * system.
 *
 *
 **********************************************************/

int
rtheta_is_cell_in_wind (n)
     int n;                     /* The wind cell number */
{
  int i, j;
  double r, theta;
  double rmin, rmax, thetamin, thetamax;
  double dr, dtheta;
  double x[3];
  int ndom, mdim, ndim;
  int ndomain;



  /* First check if the cell is in the boundary */
  ndom = wmain[n].ndom;
  wind_n_to_ij (ndom, n, &i, &j);
  ndim = zdom[ndom].ndim;
  mdim = zdom[ndom].mdim;

  if (i >= (ndim - 2) && j >= (mdim - 2))
  {
    return (W_NOT_INWIND);
  }

  /* Assume that if all four corners are in the wind that the
   * entire cell is in the wind */

  if (check_corners_inwind (n) == 4)
  {
    return (W_ALL_INWIND);
  }

  /* So at this point, we have dealt with the easy cases */

  rmin = zdom[ndom].wind_x[i];
  rmax = zdom[ndom].wind_x[i + 1];
  thetamin = zdom[ndom].wind_z[j] / RADIAN;
  thetamax = zdom[ndom].wind_z[j + 1] / RADIAN;

  dr = (rmax - rmin) / RESOLUTION;
  dtheta = (thetamax - thetamin) / RESOLUTION;

  /* Check inner and outer boundary in the z direction  */

  x[1] = 0;


  for (theta = thetamin + dtheta / 2.; theta < thetamax; theta += dtheta)
  {
    x[0] = rmin * sin (theta);
    x[2] = rmin * cos (theta);;

    if (where_in_wind (x, &ndomain) == W_ALL_INWIND)
    {
      return (W_PART_INWIND);
    }

    x[0] = rmax * sin (theta);
    x[2] = rmax * cos (theta);;
    if (where_in_wind (x, &ndomain) == W_ALL_INWIND)
    {
      return (W_PART_INWIND);
    }

  }



  for (r = rmin + dr / 2.; r < rmax; r += dr)
  {
    x[0] = r * sin (thetamin);
    x[2] = r * cos (thetamin);;
    if (where_in_wind (x, &ndomain) == W_ALL_INWIND)
    {
      return (W_PART_INWIND);
    }

    x[0] = r * sin (thetamax);
    x[2] = r * cos (thetamax);;
    if (where_in_wind (x, &ndomain) == W_ALL_INWIND)
    {
      return (W_PART_INWIND);
    }

  }


  /* If one has reached this point, then this wind cell is not in the wind */
  return (W_NOT_INWIND);


}
