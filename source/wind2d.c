
/***********************************************************/
/** @file  wind2d.c
 * @author ksl
 * @date   May, 2018
 *
 * @brief  Routines that are used to define the wind initially
 *
 * The file also contains a number of utility routines, several of
 * which are steering routines to coordinate system specific
 * routines.
 ***********************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "atomic.h"
#include "python.h"




/**********************************************************/
/**
 * @brief      initializes the structures which characterize the wind, including
 * wmain and plasmamain.
 *
 * @return     Allways returns 0
 *
 * @details
 * This is the controlling routine for initializing the wind and related structures, calculating the
 * salient properties (e.g. the volume that is in the wind) in each grid cell.
 *
 * ### Notes ###
 * This is a fairly complex routine.  In addition ot allocating space, the routine initalizes
 * the wind.  In some case the initalization is done "in-line" but more often other routines,
 * many of which coordinate system specifica are clled.
 *
 *
 **********************************************************/

int
define_wind ()
{
  int i, j, n;
  int nn;
  double rrstar;
  double x[3];
  double mdotbase, mdotwind, rr;
  int ierr;
  int n_vol, n_inwind, n_part;

  int nwind, ndom;
  int nstart, ndim, mdim;
  int nplasma;

  WindPtr w;


  /* Determine the size of the structure, we need to allocate
     from the individual domains and then allocatee the space */

  /*
     Setup the locations of starting and ending positions of the different
     domains for wmain.  The information is placed in zdom for later use.
     Note that zdom[ndom].ndim, zdom[ndom].mdim, and zdom[ndom].ndim2 shold
     have been established in get_grid_params(ndom).
   */

  n = 0;
  for (ndom = 0; ndom < geo.ndomain; ndom++)
  {
    zdom[ndom].nstart = n;
    n += zdom[ndom].ndim * zdom[ndom].mdim;
    zdom[ndom].nstop = n;
    //PLACEHOLDER  we should not really need to define NDIM and MDIM here, we may need NDIM2
    /* NDIM2 here is the total dimensions of the grid, summed over all domains
       and is used to allocate the wind pointer */
    geo.ndim2 = NDIM2 += zdom[ndom].ndim * zdom[ndom].mdim;
  }

  calloc_wind (NDIM2);

  w = wmain;

  /* Assign the domains to each cell */
  for (ndom = 0; ndom < geo.ndomain; ndom++)
  {
    for (n = zdom[ndom].nstart; n < zdom[ndom].nstop; n++)
    {
      w[n].ndom = ndom;         // Assign each wind cell to a domain
    }
  }

  /* initialize inwind and dfudge to a known state for all wind cells */
  for (n = 0; n < NDIM2; n++)
  {
    w[n].inwind = W_NOT_ASSIGNED;
    w[n].dfudge = DFUDGE;
  }

  for (ndom = 0; ndom < geo.ndomain; ndom++)
  {

    Log ("Define wind coord_type %d for domain %d\n", zdom[ndom].coord_type, ndom);

    if (zdom[ndom].wind_type == IMPORT)
    {
      import_make_grid (w, ndom);
    }
    else if (zdom[ndom].wind_type == SHELL)     /* nsh: This is the mode where we want the wind and the grid carefully
                                                   controlled to allow a very thin shell. We ensure that the coordinate type is spherical.
                                                 */
    {
      Log
        ("We are making a thin shell type grid to match a thin shell wind. This is totally aphysical and should only be used for testing purposes\n");
      shell_make_grid (w, ndom);
    }
    else if (zdom[ndom].coord_type == SPHERICAL)
    {
      spherical_make_grid (w, ndom);
    }
    else if (zdom[ndom].coord_type == CYLIND)
    {
      cylind_make_grid (ndom, w);
    }
    else if (zdom[ndom].coord_type == RTHETA)
    {
      if (zdom[ndom].wind_type == HYDRO)        /* 13jun -- nsh - 76 - This is a switch to allow one to use the
                                                   actual zeus grid in the special case of a 'proga' wind in rtheta
                                                   coordinates
                                                 */
      {
        rtheta_make_hydro_grid (w, ndom);
      }
      else
      {
        rtheta_make_grid (w, ndom);
      }
    }
    else if (zdom[ndom].coord_type == CYLVAR)
    {
      cylvar_make_grid (w, ndom);
    }
    else
    {
      Error ("define_wind: Don't know how to make coordinate type %d\n", zdom[ndom].coord_type);
    }
    for (n = zdom[ndom].nstart; n < zdom[ndom].nstop; n++)
    {
      /* For imported models we we have already set the velocities
       * at the edges of cells so this should not be done again
       */
      if (zdom[ndom].wind_type != IMPORT)
      {
        model_velocity (ndom, w[n].x, w[n].v);
      }
      model_vgrad (ndom, w[n].x, w[n].v_grad);
    }

  }
  wind_complete (w);

  /* Now determine the valid volumes of each cell and also determine whether the cells are in all
     or partially in the wind.
   */

  for (ndom = 0; ndom < geo.ndomain; ndom++)

  {
    if (zdom[ndom].coord_type == SPHERICAL)
    {
      spherical_volumes (ndom, w);
    }
    else if (zdom[ndom].coord_type == CYLIND)
    {
      cylind_volumes (ndom, w);
    }
    else if (zdom[ndom].coord_type == RTHETA)
    {
      /* 13jun -- nsh - 76 - This is a switch to allow one to use
         the actual zeus grid in the special case of a 'proga' wind
         in rtheta coordinates We dont need to work out if cells are
         in the wind, they are known to be in the wind. */
      if (zdom[ndom].wind_type == HYDRO)
      {
        rtheta_hydro_volumes (ndom, w);
      }
      else
      {
        rtheta_volumes (ndom, w);
      }
    }
    else if (zdom[ndom].coord_type == CYLVAR)
    {
      cylvar_volumes (ndom, w);
    }
    else
    {
      Error ("define_wind: Don't know how to make volumes for coordinate type %d\n", zdom[ndom].coord_type);
    }
  }




/* The routines above have established the volumes of the cells that are in the wind
 * and also assigned the variables w[].inwind at least insofar as the wind is concerned.
 */

/* Perform some consistency checks of the wind */

  NPLASMA = 0;
  for (ndom = 0; ndom < geo.ndomain; ndom++)
  {
    n_vol = n_inwind = n_part = 0;
    for (n = zdom[ndom].nstart; n < zdom[ndom].nstop; n++)
    {
      if (w[n].vol > 0.0)
        n_vol++;
      if (w[n].inwind == W_ALL_INWIND)
        n_inwind++;
      if (w[n].inwind == W_PART_INWIND)
        n_part++;

    }
    NPLASMA += n_vol;
    Log
      ("wind2d: For domain %d there are %3d cells of which %d are in inwind, %d partially in_wind, & %d with pos. vol\n",
       ndom, zdom[ndom].ndim2, n_inwind, n_part, n_vol);


    if (zdom[ndom].coord_type != SPHERICAL && zdom[ndom].wind_type != IMPORT)
    {
      for (n = zdom[ndom].nstart; n < zdom[ndom].nstop; n++)
      {
        n_inwind = check_corners_inwind (n);

        if (w[n].vol == 0 && n_inwind > 0)
        {
          wind_n_to_ij (ndom, n, &i, &j);
          Error ("wind2d: Cell %3d (%2d,%2d) in domain %d has %d corners in wind, but zero volume\n", n, i, j, ndom, n_inwind);
          w[n].inwind = W_IGNORE;
        }
        if (w[n].inwind == W_PART_INWIND && n_inwind == 4)
        {
          wind_n_to_ij (ndom, n, &i, &j);
          Error ("wind2d: Cell %3d (%2d,%2d) in domain %d has 4 corners in wind, but is only partially in wind\n", n, i, j, ndom);
        }
      }

    }
    else if (zdom[ndom].coord_type == SPHERICAL)
    {
      Log ("wind2d: Not checking corners_in_wind for SPHERICAL coordinates used in domain %d\n", ndom);
    }
    else if (zdom[ndom].wind_type == IMPORT)
    {
      /* JM 1711 -- we don't want to check corners in the case of an imported model */
      Log ("wind2d: Not checking corners_in_wind for imported model, domain %d\n", ndom);
    }
  }



  /* Allocate space for the plasma array.  To save space, this structure is sized 
   * to contain only those cells which are in the wind. */

  calloc_plasma (NPLASMA);
  calloc_dyn_plasma (NPLASMA);
  create_maps ();               /* Populate the maps between plasmamain & wmain */

  /* JM 1502 -- we want the macro structure to be allocated in geo.rt_mode = RT_MODE_MACRO. see #138  */

  if (geo.rt_mode == RT_MODE_MACRO)
  {
    calloc_macro (NPLASMA);
    calloc_estimators (NPLASMA);
  }


/* 06may -- At this point we have calculated the volumes of all of the cells and it should
be optional which variables beyond here are moved to structures othere than Wind */


/* Now calculate parameters that need to be calculated at the center of the grid cell */

  for (n = 0; n < NPLASMA; n++)
  {
    nwind = plasmamain[n].nwind;
    ndom = wmain[nwind].ndom;
    stuff_v (w[nwind].xcen, x);

    /* Next two lines allow for clumping */
    plasmamain[n].rho = model_rho (ndom, x) / zdom[ndom].fill;
    plasmamain[n].vol = w[nwind].vol * zdom[ndom].fill; // Copy volumes

    /* This is where we initialise the spectral models for the wind. */

    for (nn = 0; nn < NXBANDS; nn++)
    {
      plasmamain[n].spec_mod_type[nn] = SPEC_MOD_FAIL;  /*NSH 120817 - setting this to
                                                           a negative number means that at the outset, we assume we do not have a
                                                           suitable model for the cell */
      plasmamain[n].exp_temp[nn] = geo.tmax;    /*NSH 120817 - as an initial guess,
                                                   set this number to the hottest part of the model -
                                                   this should define where any exponential dropoff becomes important */
      plasmamain[n].exp_w[nn] = 0.0;    /* 120817 Who knows what this should be! */
      plasmamain[n].pl_alpha[nn] = geo.alpha_agn;       /*As an initial guess we assume the whole wind is
                                                           optically thin and so the spectral index for a PL illumination will be the
                                                           same everywhere.  */
      plasmamain[n].pl_log_w[nn] = -1e99;       /*131114 - a tiny weight - just to fill the variable */


      plasmamain[n].fmin_mod[nn] = 1e99;        /* Set the minium model frequency to the max frequency in the band - means it will never be used which is correct at this time - there is no model */
      plasmamain[n].fmax_mod[nn] = 1e-99;       /* Set the maximum model frequency to the min frequency in the band */

    }


/* NSH 130530 Next few lines allow the use of the temperature which can be computed from Zeus models to be
 * used as an initial guess for the wind temperature */

    if (zdom[ndom].wind_type == HYDRO)
    {
      plasmamain[n].t_r = hydro_temp (x);       //NSH 151126 - slight tidy up here - we now set t_e and t_r to hydro temp, t_e change is below//
    }
    else
    {
      plasmamain[n].t_r = zdom[ndom].twind;
    }


    /* Initialize variables having to do with converence in initial stages */
    plasmamain[n].gain = 0.5;
    plasmamain[n].dt_e_old = 0.0;
    plasmamain[n].dt_e = 0.0;

/* The next lines set the electrom temperature to 0.9 times the radiation temperature (which is a bit
	  odd since the input is the wind temperature, but is taken to be the radiation temperature). If we have
	  a fixed temprature calculation,then the wind temperature is set to be the wind temperature so the
	  user gets what they are expecting */

    if (modes.fixed_temp == 0 && modes.zeus_connect == 0)       //NSH 151126 - dont multply by 0.9 in zeus connect or fixed temp modes
      plasmamain[n].t_e = plasmamain[n].t_e_old = 0.9 * plasmamain[n].t_r;      //Lucy guess
    else
      plasmamain[n].t_e = plasmamain[n].t_e_old = plasmamain[n].t_r;
    //If we want to fix the temperature, we set it to tr which has previously been set to twind.



/* Calculate an initial guess for the weight of the PL spectrum (constant / area of a sphere / 4pi) */


    rrstar = 1. - (geo.rstar * geo.rstar) / (x[0] * x[0] + x[1] * x[1] + x[2] * x[2]);

    if (rrstar > 0)
    {
      plasmamain[n].w = 0.5 * (1 - sqrt (rrstar));
    }
    else
      plasmamain[n].w = 0.5;    //Modification to allow for possibility that grid point is inside star

    /* Determine the initial ionizations, either LTE or  fixed_concentrations */
    if (geo.ioniz_mode != IONMODE_FIXED)
    {                           /* Carry out an LTE determination of the ionization (using T_r */
      ierr = ion_abundances (&plasmamain[n], IONMODE_LTE_TR);
    }
    else
    {                           /* Set the concentrations to specified values */
      ierr = ion_abundances (&plasmamain[n], IONMODE_FIXED);
    }
    if (ierr != 0)
    {
      Error
        ("wind_define after ion_abundances: cell %d rho %8.2e t_r %8.2 t_e %8.2e w %8.2e\n",
         n, plasmamain[n].rho, plasmamain[n].t_r, plasmamain[n].t_e, plasmamain[n].w);
    }

    /* Initialize arrays for scatters 
     */

    for (j = 0; j < nions; j++)
    {
      plasmamain[n].scatters[j] = 0;
      plasmamain[n].xscatters[j] = 0;
    }
  }


/* Calculate the the divergence of the wind at the center of each grid cell */
  wind_div_v ();

/* Now calculate the adiabatic cooling and shock heating */
  for (i = 0; i < NPLASMA; i++)
  {
    if (geo.adiabatic)
    {
      nwind = plasmamain[i].nwind;
      plasmamain[i].cool_adiabatic = adiabatic_cooling (&w[nwind], plasmamain[i].t_e);
    }
    else
    {
      plasmamain[i].cool_adiabatic = 0.0;
    }

    if (geo.nonthermal)
    {
      nwind = plasmamain[i].nwind;
      plasmamain[i].heat_shock = shock_heating (&w[nwind]);
    }
    else
    {
      plasmamain[i].heat_shock = 0.0;
    }
  }




  /* Calculate one over dvds */
  dvds_ave ();


  wind_check (w, -1);           // Check the wind for reasonability

  /* zero the counters which record diagnostics for mean_intensity */
  nerr_Jmodel_wrong_freq = 0;
  nerr_no_Jmodel = 0;

/*
     Check that m_dot is correct.  This calculation is very approximate.  It only calculates mdot
     flowing up through the grid, i.e. any material that leaks out the side will not be counted.  A
     much better approach is needed, e.g. to calculate the mass flow in a spherical shell at the edge
     of the disk.  In addition, it uses a very inexact way to determine which grid points are in the
     wind.  Today, I simply modified the routine so the flux was calculated mid-way in grid space
     through the wind. ??? ksl 01mar15

     04aug -- I have not written the routine to do the check that with cooridnate systems
     other than spherical, we are getting the correct mdot.  It should be possible to do this
     in a coordinate system independent manner simply by inteegrating using interpolated values
     and so this would be a good thing to fix.  But it is not used elsewhere so I have skipped
     it for now.

     05apr -- Looked at this again.  It is clearly a straightforwrd thing to do use the rho function
     which is below. One would simply integrate over various surface integrals to make it happen.

     15aug -- ksl - It's not clear this it really correct with dommains, but the following compiles
*/

  for (ndom = 0; ndom < geo.ndomain; ndom++)
  {
    if (zdom[ndom].coord_type == CYLIND)
    {
      Log ("Warning: wind2d.c: next lines may not make sense with multiple domains\n");
      ndim = zdom[ndom].ndim;
      mdim = zdom[ndom].mdim;
      nstart = zdom[ndom].nstart;

      mdotwind = 0;
      mdotbase = 0;
      for (i = 0; i < ndim - 1; i++)
      {
        n = nstart + i * mdim + (mdim / 2);
// rr is router * router - rinner *rinner for this shell
        rr = w[n + mdim].x[0] * w[n + mdim].x[0] + w[n + mdim].x[1] * w[n + mdim].x[1] - (w[n].x[0] * w[n].x[0] + w[n].x[1] * w[n].x[1]);
        if (w[nstart + i * mdim].inwind == W_ALL_INWIND)
        {
          nplasma = w[nstart + i * mdim].nplasma;
          mdotbase += plasmamain[nplasma].rho * PI * rr * w[nstart + i * mdim].v[2];
        }
        if (w[n].inwind == W_ALL_INWIND)
        {
          nplasma = w[n].nplasma;
          mdotwind += plasmamain[nplasma].rho * PI * rr * w[n].v[2];
        }
      }

      /* Factor of 2 to allow for wind on both sides of disk */
      Log ("m dot wind: Domain %d Desired %g   Base %g Calculated %g\n", ndom, zdom[ndom].wind_mdot, 2. * mdotbase, 2. * mdotwind);
      mdot_wind (w, 1.e6, zdom[ndom].rmax / 2.);
    }
    else
    {
      Error ("wind2d.c: Not currently able to calculate mdot wind for coordtype %d in domain %d\n", zdom[ndom].coord_type, ndom);
    }
  }

  return (0);
}



int wig_n;
double wig_x, wig_y, wig_z;

/**********************************************************/
/**
 * @brief      locates the element in wmain associated with a postion
 *
 * @param [in] int  ndom   The domain number for the search
 * @param [in] double  x[]   The position
 * @return     where_in_grid normally  returns the element in wmain associated
 * with a position.  If the positions is in the grid (of one of the domains)
 * this will be a positive integer.  If the position is not in the grid of one
 * of the domains,
 * a negative number -1 will be returned if the position is inside athe grid,
 * -2 if it is outside the grid for that domain
 *
 * @details
 *
 * ### Notes ###
 * where in grid is mainly a steering routine that calls other various
 * coordinate system specific routines.
 *
 * Where_in grid does not tell you whether the position is in the wind!.
 *
 * What one means by inside or outside the grid may well be different
 * for different coordinate systems.
 *
 * This is one of the routines in Python that has a long history, in that
 * it existed prior to the implementation of domains.  This is why it
 * returns the position in wmain, rather than the position in one of
 * the plasma domains.  It would make sense to revise this.
 *
 **********************************************************/

int
where_in_grid (ndom, x)
     int ndom;
     double x[];
{
  int n;
  double fx, fz;


  if (wig_x != x[0] || wig_y != x[1] || wig_z != x[2])  // Calculate if new position
  {

    if (zdom[ndom].coord_type == CYLIND)
    {
      n = cylind_where_in_grid (ndom, x);
    }
    else if (zdom[ndom].coord_type == RTHETA)
    {
      n = rtheta_where_in_grid (ndom, x);
    }
    else if (zdom[ndom].coord_type == SPHERICAL)
    {
      n = spherical_where_in_grid (ndom, x);
    }
    else if (zdom[ndom].coord_type == CYLVAR)
    {
      n = cylvar_where_in_grid (ndom, x, 0, &fx, &fz);
    }
    else
    {
      Error ("where_in_grid: Unknown coord_type %d for domain %d\n", zdom[ndom].coord_type, ndom);
      Exit (0);
      return (0);
    }


    /* Store old positions to short-circuit calculation if asked for same position more
       than once */
    wig_x = x[0];
    wig_y = x[1];
    wig_z = x[2];
    wig_n = n;
  }

  return (wig_n);
}

int ierr_vwind = 0;


/**********************************************************/
/**
 * @brief      finds the velocity vector v for the wind in cartesian
 * 	coordinates at the position of the photon p.
 *
 * @param [in] int  ndom   The domain of interest
 * @param [in] PhotPtr  p   A photon
 * @param [out] double  v[]   The velocity at the postion given by phtoon???
 * @return     0 	on successful completion, or a negative number if
 * photon position is not in the in the grid for the domain of interest.
 *
 * @details
 * This routine carries out a bilinear interpolation of the wind velocity.  The velocities
 * at the edges of the grid cells must have been defined elsewhere (e.g in wind_define).
 *
 * If the position is outside the wind grid, then the velocities will be returned as 0.
 *
 * ### Notes ###
 * For all of the 2d coordinate systems, the positions in the WindPtr array are
 * calculated in the xz plane.  As a result, the velocities in that array are
 * calcuated there as well.  Hence, to obtain the velocity in cartesion coordinates
 * one needs to do a simple rotation of the velocity vector.
 *
 * For spherical (1d) coordinates the situation is complex, the philosopy that
 * has been adopted is to let the program run as if the wind model is intrinsically
 * 2d.   The xyz positions of the grid are defined are defined so that
 * are in the xz plane at an angle of 45 degrees to the x axis.  This was done
 * so that if one used, say an sv wind, in spherical coordinates, one would select
 * the wind velocities at a plausible location.   But it implies that the velocities
 * that have been accepted are not purely radial as one might expect.  The choice
 * that has been made is to assume that v[1] is preserved as a non zero value in
 * making the projection.  An alternative might have been to preserve the angular
 * rotation rate.
 *
 * If we ever implement a real 3d coordinate system, one will need to look at
 * this routine again.
 *
 **********************************************************/

int
vwind_xyz (ndom, p, v)
     int ndom;
     PhotPtr p;
     double v[];
{
  int i;
  double rho, r;
  double vv[3];
  double ctheta, stheta;
  double x, frac[4];
  int nn, nnn[4], nelem;



  if (ndom < 0 || ndom >= geo.ndomain)
  {
    Error ("vwind_xyz: Received invalid domain  %d\n", ndom);
  }

  coord_fraction (ndom, 0, p->x, nnn, frac, &nelem);

  for (i = 0; i < 3; i++)
  {

    x = 0;
    for (nn = 0; nn < nelem; nn++)
      x += wmain[nnn[nn]].v[i] * frac[nn];

    vv[i] = x;
  }

  rho = sqrt (p->x[0] * p->x[0] + p->x[1] * p->x[1]);

  if (zdom[ndom].coord_type == SPHERICAL)
  {                             // put the velocity into cylindrical coordinates on the xz axis
    x = sqrt (vv[0] * vv[0] + vv[2] * vv[2]);   //v_r
    r = length (p->x);
    vv[0] = rho / r * x;
    vv[2] = p->x[2] / r * x;
  }
  else if (p->x[2] < 0)         // For 2d coord syatems, velocity is reversed if the photon is in the lower hemisphere.
    vv[2] *= -1;

  if (rho == 0)
  {                             // Then we will not be able to project from a cylindrical ot a cartesian system
    Error ("vwind_xyz: Cannot determine an xyz velocity on z axis. Returnin 0,0,v[2]\n");
    v[0] = v[1] = 0;
    v[2] = vv[2];
    return (0);
  }

  /* Now we have v in cylindrical coordinates, but we would like it in cartesian coordinates.
     Note that could use project_from_cyl_xyz(p->x,vv,v) ?? */

  ctheta = p->x[0] / rho;
  stheta = p->x[1] / rho;
  v[0] = vv[0] * ctheta - vv[1] * stheta;
  v[1] = vv[0] * stheta + vv[1] * ctheta;
  v[2] = vv[2];

  if (sane_check (v[0]) || sane_check (v[1]) || sane_check (v[2]))
  {
    Error ("vwind_xyz: %f %f %f\n", v[0], v[1], v[2]);
  }


  return (0);
}




/**********************************************************/
/**
 * @brief      calculates the divergence of the velocity at the center of all the grid cells.
 *
 * @param [in,out] WindPtr  w   The entire wind domain
 * @return     Always returns 0
 *
 * The results are stored in wmain[].div_v
 *
 * @details
 * This is one of the initialization routines for the wind.
 * The divergence of velocity of the wind is used to calculated PdV cooling of a grid cell
 *
 * ### Notes ###
 *
 * Because the routine calls vwind_xyz, one must have previously
 * populated w[].v.  Also as a result of this fact, the routine
 * is not tightly tied to any particular wind model.
 *
 * The divergence, like othe scalar quantities, does not
 * need to be "rotated" differently in different coordinate
 * systems
 *
 *
 **********************************************************/
int wind_div_err = (-3);

int
wind_div_v ()
{
  int icell;
  double x_zero[3], v2[3], v1[3];
  struct photon ppp;
  double div, delta;
  double xxx[3];
  int ndom;
  double scaling;

  scaling = 1e-3;               //The scaling factor applied to 'delta' the distance away from the central point that the div_v calcs are done


  for (icell = 0; icell < NDIM2; icell++)
  {
    /* Find the center of the cell */

    stuff_v (wmain[icell].xcen, x_zero);        /*Gget the centre of the current cell in the loop */
    ndom = wmain[icell].ndom;

//    delta = 0.01 * x_zero[2];   //delta is the distance across which we measure e.g. dv_x/dx

    if (x_zero[1] != 0)
    {
      delta = fabs (fmin (wmain[icell].x[0] - x_zero[0], fmin (wmain[icell].x[1] - x_zero[1], wmain[icell].x[2] - x_zero[2])));
    }
    else
    {
      delta = fabs (fmin (wmain[icell].x[0] - x_zero[0], wmain[icell].x[2] - x_zero[2]));
    }
    delta = delta * scaling;


    if (delta == 0)
    {
      Error ("wind_div_v: Cell %d has xcen[2]==0.  This is surprising\n", icell);
      delta = wmain[icell].dfudge;
    }




    /* for each of x,y,z we first create a copy of the vector at the center. We then step 0.5*delta
       in positive and negative directions and evaluate the difference in velocities. Dividing this by
       delta gives the value of dv_x/dx, and the sum of these gives the divergence. If issues arise
       see bug report #70. */


    /* Calculate dv_x/dx at this position */
    stuff_v (x_zero, ppp.x);
    ppp.x[0] += 0.5 * delta;
    vwind_xyz (ndom, &ppp, v2);
    ppp.x[0] -= delta;
    vwind_xyz (ndom, &ppp, v1);
    div = xxx[0] = (v2[0] - v1[0]) / delta;

    /* Calculate dv_y/dy */
    stuff_v (x_zero, ppp.x);
    ppp.x[1] += 0.5 * delta;
    vwind_xyz (ndom, &ppp, v2);
    ppp.x[1] -= delta;
    vwind_xyz (ndom, &ppp, v1);
    div += xxx[1] = (v2[1] - v1[1]) / delta;


    /* Calculate dv_z/dz */
    stuff_v (x_zero, ppp.x);
    ppp.x[2] += 0.5 * delta;
    vwind_xyz (ndom, &ppp, v2);
    ppp.x[2] -= delta;
    vwind_xyz (ndom, &ppp, v1);
    div += xxx[2] = (v2[2] - v1[2]) / delta;


    /* we have now evaluated the divergence, so can store in the wind pointer */
    wmain[icell].div_v = div;

    if (div < 0 && (wind_div_err < 0 || wmain[icell].inwind == W_ALL_INWIND))
    {
      Error
        ("wind_div_v: div v %e negative in cell %d Domain %d. Major problem if inwind (%d) == 0\n",
         div, icell, wmain[icell].ndom, wmain[icell].inwind);
      wind_div_err++;
    }
  }

  return (0);
}




/**********************************************************/
/**
 * @brief      find the density of the wind at x
 *
 * @param [in, out] WindPtr  w   The entire wind
 * @param [in, out] double  x[]   A position
 * @return
 * The density at x, if the postion is in the active region of the wind.
 * If the postion is not in the active region of one of the domains, then
 * 0 is returned.  No error is reported for this
 *
 * @details
 * The routine first determines what domain the position is located in,
 * and then intepolates to find the density at a specific position
 *
 * ### Notes ###
 *
 **********************************************************/

double
rho (w, x)
     WindPtr w;
     double x[];
{
  int n;
  double dd;
  double frac[4];
  int nn, nnn[4], nelem;
  int nplasma;
  int ndom, ndomain;


  if (where_in_wind (x, &ndomain) < 0)  //note that where_in_wind is independent of grid.
    return (0.0);

  if ((ndom = ndomain) < 0)
  {
    Error ("rho: Domain not fournd %d\n", ndomain);
    return (0.0);
  }

  n = coord_fraction (ndom, 1, x, nnn, frac, &nelem);

  if (n < 0)
  {
    dd = 0;
  }
  else
  {

    dd = 0;
    for (nn = 0; nn < nelem; nn++)
    {
      nplasma = w[nnn[nn]].nplasma;
      dd += plasmamain[nplasma].rho * frac[nn];
    }

  }


  return (dd);
}



/**********************************************************/
/**
 * @brief      The routine calculates and then logs the mass loss rate in two ways
 *
 * @param [in] WindPtr  w   The entire wind domain
 * @param [in] double  z   A zheight above the disk to calculate the mass loss rate
 * @param [in] double  rmax   A radius at which to calculate the mass loss rate
 * @return     Always retruns 0
 *
 * @details
 * The routine calculates the mass loss rate in two ways, one in a plane with
 * a zheight of 0 above the disk (out to a distance rmax, and one at a spherical
 * radius rmax
 *
 * ### Notes ###
 *
 * ### Programming Comment ###
 * This routine has errors.  It is only correct if there is
 * a single domain (although it does provide a warning to this effect.
 * In particular rho(w,x) gives the correct answer
 * for rho regardless of domains, but ndom is set to 0 for vind(ndom,&p,v).
 * This should be fixed.  This is now #395
 *
 **********************************************************/
#define NSTEPS 100

int
mdot_wind (w, z, rmax)
     WindPtr w;
     double z;                  // The height (usually small) above the disk at which mdot will be calculated
     double rmax;               // The radius at which mdot will be calculated
{
  struct photon p;              //needed because vwind_xyz has a call which involves a photon
  double r, dr, rmin;
  double theta, dtheta;
  double den, rho ();
  double mdot, mplane, msphere;
  double x[3], v[3], q[3], dot ();
  int ndom;

  ndom = 0;

  Log ("For simplicity, mdot wind checks only carried out for domain 0\n");

/* Calculate the mass loss rate immediately above the disk */

  /* Check that everything is defined sensibly */
  rmin = geo.rstar;

  if (rmax <= rmin)
  {
    Error ("mdot_wind: rmax %g is less than %g, so returning\n", rmax, rmin);
    return (0.0);
  }

  dr = (rmax - rmin) / NSTEPS;
  mdot = 0;
  p.x[1] = x[1] = 0.0;
  p.x[2] = x[2] = z;
  for (r = dr / 2; r < rmax; r += dr)
  {
    p.x[0] = x[0] = r;
    den = rho (w, x);
    vwind_xyz (ndom, &p, v);
    mdot += 2 * PI * r * dr * den * v[2];
  }

  mplane = 2. * mdot;

/* Calculate the mass loss rate in a sphere */

  dtheta = PI / (2 * NSTEPS);
  mdot = 0;
  p.x[1] = x[1] = 0.0;
  q[1] = 0;
  for (theta = dtheta / 2; theta < (PI / 2.); theta += dtheta)
  {
    p.x[0] = x[0] = r * sin (theta);
    p.x[2] = x[2] = r * cos (theta);
    q[0] = sin (theta);
    q[2] = cos (theta);
    den = rho (w, x);
    vwind_xyz (ndom, &p, v);
    mdot += 2 * PI * rmax * rmax * sin (theta) * dtheta * den * dot (v, q);
  }
  msphere = 2 * mdot;           // Factor of two because wind is in both hemispheres

  Log ("Wind Mdot Desired %e Plane %e Sphere(%e) %e\n", zdom[ndom].wind_mdot, mplane, rmax, msphere);

  return (0);
}



/**********************************************************/
/**
 * @brief      is simply will produce a postion at a random place in
 * 	a cell n
 *
 * @param [in] int  n   An element number in the wind structure
 * @param [out] double  x[]   The random location that is calculaed
 * @return     Always returns 0
 *
 * @details
 * This is really just a driver routine for coordinate system specific
 * routines
 *
 * ### Notes ###
 *
 **********************************************************/

int
get_random_location (n, x)
     int n;                     // Cell in which to create position
     double x[];                // Returned position
{
  int ndom;

  ndom = wmain[n].ndom;

  if (zdom[ndom].coord_type == CYLIND)
  {
    cylind_get_random_location (n, x);
  }
  else if (zdom[ndom].coord_type == RTHETA)
  {
    rtheta_get_random_location (n, x);
  }
  else if (zdom[ndom].coord_type == SPHERICAL)
  {
    spherical_get_random_location (n, x);
  }
  else if (zdom[ndom].coord_type == CYLVAR)
  {
    cylvar_get_random_location (n, x);
  }
  else
  {
    Error ("get_random_location: Don't know this coord_type %d\n", zdom[ndom].coord_type);
    Exit (0);
  }

  return (0);
}



/**********************************************************/
/**
 * @brief      zero out the portion of plasmamain that records
 * the number of scatters in each cell
 *
 * @return     Always returns 0
 *
 * @details
 *
 * ### Notes ###
 * The routine is called at the point where we begin to
 * calculate the detailed spectrum.  It's a little unclear
 * why it is simply not incorporated into spectrum_init
 *
 **********************************************************/

int
zero_scatters ()
{
  int n, j;

  for (n = 0; n < NPLASMA; n++) /*NSH 1107 - changed loop to only run over nions to avoid running over the
                                   end of the array after the arry was dynamically allocated in 78a */
  {
    for (j = 0; j < nions; j++)
    {
      plasmamain[n].scatters[j] = 0;
    }
  }

  return (0);
}




/**********************************************************/
/**
 * @brief      The routine basically just calls where_in_wind for each of 
 ' the 4 corners of a cell in one of the 2d coordinate systems.
 *
 * @param [in] int  n   the cell number in wmain
 * @return     the number of corners of the cell that
 * are in the wind
 *
 * @details
 * This is really just a driver routine which calls coordinate system specific
 * routines
 *
 * It is used to  standardize this check for calculations of 
 * the volume of a cell that is in the wind in a domain.
 *
 * ### Notes ###
 * This is one of the routines that makes explicity assumptions
 * about guard cells, and it is not entirely clear why
 *
 * This rutine is not called for a spherical coordinate system
 *
 **********************************************************/

int
check_corners_inwind (n)
     int n;
{
  int n_inwind;
  int i, j;
  DomainPtr one_dom;
  int ndom, ndomain;

  /* find the domain */
  ndom = wmain[n].ndom;
  one_dom = &zdom[ndom];

  wind_n_to_ij (ndom, n, &i, &j);

  n_inwind = 0;

  if (i < (one_dom->ndim - 2) && j < (one_dom->mdim - 2))
  {
    if (where_in_wind (wmain[n].x, &ndomain) == W_ALL_INWIND && ndom == ndomain)
      n_inwind++;
    if (where_in_wind (wmain[n + 1].x, &ndomain) == W_ALL_INWIND && ndom == ndomain)
      n_inwind++;
    if (where_in_wind (wmain[n + one_dom->mdim].x, &ndomain) == W_ALL_INWIND && ndom == ndomain)
      n_inwind++;
    if (where_in_wind (wmain[n + one_dom->mdim + 1].x, &ndomain) == W_ALL_INWIND && ndom == ndomain)
      n_inwind++;
  }

  return (n_inwind);
}




/**********************************************************/
/**
 * @brief check whether some quantities in a grid cell are
 * changing so rapidly in the grid cell that one might
 * wish to have smaller grid cells
 *
 * @return     Always returns 0
 *
 * The results are logged.
 *
 * @details
 * loops over the entire wind points and checks quantities
 * like velocity gradients and optical depths across cells for large changes.
 * If there are large changes in cells, it tells the user how many cells and
 * suggests that the grid may need to be modified.
 *
 *
 * ### Notes ###
 * Must be called after define wind so NPLASMA and densities are set up.
 *
 * The information that is checked is very basic and it is clear that
 * someone was thinking about doing somewhat better checks.
 *
 * Errors are not generated, or rather have been commented out, and
 * so one needs to consider whether this is useful or not.  If so there
 * should be errors.
 *
 **********************************************************/

int
check_grid ()
{
  int ndom, n;
  double l_sob, vth, lambda_t, nh;
  double delta_r, delta_x, delta_z, delta_vz, delta_vx;
  double v1[3], v2[3];
  WindPtr one;
  PlasmaPtr xplasma;
  int n_dv, n_tau;

  n_dv = n_tau = 0;

  for (n = 0; n < NPLASMA; n++)
  {
    /* get the relevant pointers for this cell */
    xplasma = &plasmamain[n];
    one = &wmain[xplasma->nplasma];
    ndom = one->ndom;

    /* Hydrogen density, ne should be roughly this */
    nh = xplasma->rho * rho2nh;

    /* thermal speed */
    vth = sqrt (1.5 * BOLTZMANN * xplasma->t_e / MPROT);

    /* sobolev length -- this could be used for a check but isn't yet */
    l_sob = vth / one->dvds_ave;

    /* get approximate cell extents */
    delta_z = 2.0 * (one->xcen[2] - one->x[2]);
    delta_x = 2.0 * (one->xcen[0] - one->x[0]);

    delta_r = sqrt (delta_x * delta_x + delta_z * delta_z);

    /* Thomson mean free path 1/sigma*nh */
    lambda_t = 1.0 / THOMPSON * nh;

    /* get the velocity at cell corner and cell edge in x and z */
    model_velocity (ndom, one->x, v1);
    model_velocity (ndom, one->xcen, v2);

    delta_vx = fabs (v1[0] - v2[0]);
    delta_vz = fabs (v1[1] - v2[1]);

    /* if we have errors then increment the error counters */
    if (delta_r / lambda_t > 1)
    {
      n_tau++;
      //Error("check_grid: optical depth may be high in cell %i domain %i\n",
      //       n, ndom);
    }

    if (delta_vx > 1e8 || delta_vz > 1e8)
    {
      n_dv++;
      //Error("check_grid: velocity changes by >1,000 km/s across cell %i domain %i\n",
      //       n, ndom);
    }

  }

  /* report problems to user in a summary error message */
  if (n_dv > 0)
    Error ("check_grid: velocity changes by >1,000 km/s in %i cells\n", n_dv);

  if (n_tau > 0)
    Error ("check_grid: optical depth may be high in %i\n", n_tau);

  if (n_dv > 0 || n_tau > 0)
    Error ("check_grid: some cells have large changes. Consider modifying zlog_scale or grid dims\n");

  return (0);
}
