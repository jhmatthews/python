#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "atomic.h"
#include "python.h"


/***********************************************************
                                       Space Telescope Science Institute

 Synopsis:
	define_wind  initializes the structure which characterize the wind.  

Arguments:		
	WindPtr w;			The structure which defines the wind in Python
 
Returns:
 
Description:

	This the basic routine which initiallizes the wind structure, calculating the
	salient properties of the wind in each grid cell. 
		
Notes:  This routine has been modified a lot as multiple coordinate systems were implemented
        Now it is mostly a driver for routines found elsewhere.
	

History:
 	97jan      ksl	Coding on python began.
 	97sep13	ksl	Added the possibility of fixing the concentrations of the ions to specific
			values
 	97nov19 ksl	Fixed dmdot_da to account for the fact that mass is lost from both sides
                        of the disk.
 	98mar20	ksl	Fixed calculation of wind_midx and wind_midz.  New code is actually
 			in wind_complete
 	98may24	ksl	Modified calculations of wind volumes to force the edges to have zero volume
 			because the edge points are actually only used for interpolation and photons
 			do not interact in them.
 	98dec	ksl	Eliminated wind element .rzero, Added capability to define linear or logarithmic
 			intervals on the fly  Isolated all dependences on Shlossman and Vitello 
 			prescription of the wind to subroutine calls, which are in the file sv.c
 	98dec22	ksl	Added routines to allow definition of a one d stellar wind.  Weights inside
 			star artificially set to 0.5
	99dec29 ksl	Added routines for proga's wind 
	00oct06	ksl	Added routines for a corona
	01mar13	ksl	Added routines for knigge's wind
	01dec26	ksl	Deleted wind_summary routine since I haven't looked at this output for several 
			years (largely replaced by py_wind)
	02feb	ksl	Made changes to allow unequal dimensions in two directions
        04May   SS      Factor of two increase in the cell volumes to allow for space above and below plane.
	04aug	ksl	Began removing portions of this routine which were coordinate sytem independent
			to individial routines in wind.c.  Also redefined what is meant by .inwind.
        04Aug   SS      Modified the calulation of whether a cell is in the wind because it was 
                        missing cells when a finite disk is included (the cell corners could lie in the disk
                        but there could still be genuine wind material within the cell).
	05jul	ksl	56d -- Added checks of volume calculations for two d models as part of general
			effort to properly incoporate a vertically extneded disk.
	06may	ksl	Modified order in which items were calculated to prepare for splitting
			of the Wind structure.  Moved the creation of thw arrays here as well
			Allowed for wheich variables need to go into plasma
	06nov	ksl	58b -- Modified to identify cells which have a a calculated wind volume
			of zero, but one or more corners in the wind.  The assumption is that
			photons were intended to simply fly through these cells.  Also modified to
			use #define variables, such as W_ALL_INWIND, for describing the type of gridcell.
	07jul	kls	58f -- Added code to copy volumes from wmain[].vol to plasmamain[].vol
        12aug   nsh	73d -- Added some code to zero some parameters in the plasma structure used to 
			try and speed up the pairwise ionization scheme
	13may	nsh	added code to use temprature computed from internal energy in a
			zeus file as a first guess for temperature in a proga type model
	13sep	nsh	small change to avoid the muliplicatoin by 0.9 (lucy guess) from giving us the 				wrong zeus temperature.
	13nov	nsh	changes to take account of log formulation of power law
			also initialise the min and max frequencies seen in a band.
	14jul	nsh	added call to calloc_dyn_plasma to allocate space for arrays of variable size - currently
			those with length nion.

**************************************************************/


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
  int n_comp, n_comp_part;

  int nwind;
  int nplasma;

  WindPtr w;

  /* In order to interpolate the velocity (and other) vectors out to geo.rmax, we need
     to define the wind at least one grid cell outside the region in which we want photons
     to propagate.  This is the reason we divide by NDIM-2 here, rather than NDIM-1 */

  NDIM = ndim = geo.ndim;
  MDIM = mdim = geo.mdim;
  NDIM2 = NDIM * MDIM;
  calloc_wind (NDIM2);
  w = wmain;

  /* initialize inwind to a known state */

  for (n = 0; n < NDIM2; n++)
    {
      w[n].inwind = W_NOT_INWIND;
    }



  if (geo.wind_type == 9)	//This is the mode where we want the wind and the grid carefully controlled to allow a very thin shell. We ensure that the coordinate type is spherical. 
    {
      Log
	("We are making a thin shell type grid to match a thin shell wind. This is totally aphysical and should only be used for testing purposes\n");
      shell_make_grid (w);
    }
  else if (geo.coord_type == SPHERICAL)
    {
      spherical_make_grid (w);
    }
  else if (geo.coord_type == CYLIND)
    {
      cylind_make_grid (w);
    }
  else if (geo.coord_type == RTHETA)
    {
      if (geo.wind_type == 3) //13jun -- nsh - 76 - This is a switch to allow one to use the actual zeus grid in the special case of a 'proga' wind in rtheta coordinates
      	{
	rtheta_make_zeus_grid (w);
      	}
      else
	{
        rtheta_make_grid (w);
	}
    }
  else if (geo.coord_type == CYLVAR)
    {
      cylvar_make_grid (w);
    }
  else
    {
      Error ("define_wind: Don't know how to make coordinate type %d\n",
	     geo.coord_type);
    }
  for (n = 0; n < NDIM2; n++)
    {
      /* 04aug -- ksl -52 -- The next couple of lines are part of the changes
       * made in the program to allow more that one coordinate system in python 
       */
      model_velocity (w[n].x, w[n].v);
      model_vgrad (w[n].x, w[n].v_grad);
    }

  wind_complete (w);

/* wind_complete has to precede the volume calculations because the 1d vectors are used
in the volume calculations.  wind_complete itself is just a driver routine.  It has to
be called out as a separate routine, because the 1d vectors are not saved and have to be
recreated when a windfile is read into the program
06may ksl */


  /* Now define the valid volumes of each cell and also determine whether the cells are in all
   * or partially in the wind.
   *
   * Note - 05apr --55d -- ksl.  Previously there was a separate calculation here of whether a cell
   * was in the wind and its volume.  Given the increasingly complicated geometry when thick
   * disks were first allowed, Stuart had gone back to determining whether a cell was all or
   * partly in the wind by checking a bunch of positions in the cell.  But this is almost idetnaical 
   * to calculating the volume of the cell, and therefore I have moved this functionality
   * in the the respective volumes calculation.  At least then we will do the calculation the
   * same way both times.  . 
   */


  if (geo.coord_type == SPHERICAL)
    {
      spherical_volumes (w, W_ALL_INWIND);
    }
  else if (geo.coord_type == CYLIND)
    {
      cylind_volumes (w, W_ALL_INWIND);
    }
  else if (geo.coord_type == RTHETA)
    {
      if (geo.wind_type == 3) //13jun -- nsh - 76 - This is a switch to allow one to use the actual zeus grid in the special case of a 'proga' wind in rtheta coordinates We dont need to work out if cells are in the wind, they are known to be in the wind.
      	{
	rtheta_zeus_volumes (w);
      	}
      else
	{
        rtheta_volumes (w, W_ALL_INWIND);
	}
    }
  else if (geo.coord_type == CYLVAR)
    {
      cylvar_volumes (w, W_ALL_INWIND);
    }
  else
    {
      Error
	("wind2d.c: Don't know how to make volumes for coordinate type %d\n",
	 geo.coord_type);
    }

/* Now check if there is a second component and if so get the volumes for these cells as well */

  if (geo.compton_torus)
    {

      if (geo.coord_type == SPHERICAL)
	{
	  spherical_volumes (w, W_ALL_INTORUS);
	}
      else if (geo.coord_type == CYLIND)
	{
	  cylind_volumes (w, W_ALL_INTORUS);
	}
      else if (geo.coord_type == RTHETA)
	{
	  rtheta_volumes (w, W_ALL_INTORUS);
	}
      else if (geo.coord_type == CYLVAR)
	{
	  cylvar_volumes (w, W_ALL_INTORUS);
	}
      else
	{
	  Error
	    ("wind2d.c: Don't know how to make volumes for coordinate type %d\n",
	     geo.coord_type);
	}

    }

/* The routines above have established the volumes of the cells that are in the wind
 * and also assigned the variables w[].inwind at least insofar as the wind is concerned.
 * We now need to do the same for the torus
 */

  n_vol = n_inwind = n_part = 0;
  n_comp = n_comp_part = 0;
  for (n = 0; n < NDIM2; n++)
    {
      if (w[n].vol > 0.0)
	n_vol++;
      if (w[n].inwind == W_ALL_INWIND)
	n_inwind++;
      if (w[n].inwind == W_PART_INWIND)
	n_part++;
      if (w[n].inwind == W_ALL_INTORUS)
	n_comp++;
      if (w[n].inwind == W_PART_INTORUS)
	n_comp_part++;
    }

  Log
    ("wind2d: %3d cells of which %d are in inwind, %d partially in_wind, & %d with pos. vol\n",
     NDIM2, n_inwind, n_part, n_vol);

  if (geo.compton_torus)
    {
      Log
	("wind2d: cells of which %d are in in the torus , %d partially ini the torus\n",
	 n_comp, n_comp_part);
    }

/* 56d --Now check the volume calculations for 2d wind models 
   58b --If corners are in the wind, but there is zero_volume then ignore.
*/
  if (geo.coord_type != SPHERICAL)
    {
      for (n = 0; n < NDIM2; n++)
	{
	  n_inwind = check_corners_inwind (n, 0);
	  if (w[n].vol == 0 && n_inwind > 0)
	    {
	      wind_n_to_ij (n, &i, &j);
	      Error
		("wind2d: Cell %3d (%2d,%2d) has %d corners in wind, but zero volume\n",
		 n, i, j, n_inwind);
	      w[n].inwind = W_IGNORE;
	    }
	  if (w[n].inwind == W_PART_INWIND && n_inwind == 4)
	    {
	      wind_n_to_ij (n, &i, &j);
	      Error
		("wind2d: Cell %3d (%2d,%2d) has 4 corners in wind, but is only partially in wind\n",
		 i, j, n);
	    }
	}

    }


/* Now create the second structure, the one that is sized only to contain cells in the wind */

  if (CHOICE)
    {
      NPLASMA = n_vol;
    }
  else				/* Force NPLASMA to equal NDIM2 (for diagnostic reasons) */
    {
      NPLASMA = NDIM2;
    }

  calloc_plasma (NPLASMA);
  calloc_dyn_plasma (NPLASMA); /*78a NSH 1407 - allocate space for dynamically sized arrays*/
  create_maps (CHOICE);		// Populate the maps from plasmamain & wmain

  /* JM 1502 -- we want the macro structure to be allocated in geo.rt_mode = 2. see #138  */
  
  if (geo.rt_mode == 2)
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
      stuff_v (w[nwind].xcen, x);
      /* 140905 - ksl - Next two lines allow for clumping */
      plasmamain[n].rho = model_rho (x)/geo.fill;
      plasmamain[n].vol = w[nwind].vol*geo.fill;	// Copy volumes
/* NSH 120817 This is where we initialise the spectral models for the wind. The pl stuff is old, I've put new things in here to initialise the exponential models */
      for (nn = 0; nn < NXBANDS; nn++)
	{
	  plasmamain[n].spec_mod_type[nn] = -1;	/*NSH 120817 - setting this to a negative number means that at the outset, we assume we do not have a suitable model for the cell */
	  plasmamain[n].exp_temp[nn] = geo.tmax;	/*NSH 120817 - as an initial guess, set this number to the hottest part of the model - this should define where any exponential dropoff becomes important */
	  plasmamain[n].exp_w[nn] = 0.0;	/* 120817 Who knows what this should be! */
	  plasmamain[n].pl_alpha[nn] = geo.alpha_agn;	//As an initial guess we assume the whole wind is optically thin and so the spectral index for a PL illumination will be the same everywhere.
	  /*     plasmamain[n].pl_w[nn] = geo.const_agn / (4.0*PI*(x[0] * x[0] + x[1] * x[1] + x[2] * x[2]));  // constant / area of a sphere
	     plasmamain[n].pl_w[nn] /= 4.*PI;   // take account of solid angle NSH 120817 removed - if PL not suitable, it will be set to zero anyway, so safe to keep it at zero from the outset! */
//	  plasmamain[n].pl_w[nn] = 0.0;
	  plasmamain[n].pl_log_w[nn] = -1e99; /*131114 - a tiny weight - just to fill the variable */ 


          plasmamain[n].fmin_mod[nn] = 1e99; /* Set the minium model frequency to the max frequency in the band - means it will never be used which is correct at this time - there is no model */
          plasmamain[n].fmax_mod[nn] = 1e-99; /* Set the maximum model frequency to the min frequency in the band */
	}


/* NSH 130530 Next few lines allow the use of the temperature which can be computed from Zeus models to be used as an initial guess for the wind temperature */
      if (geo.wind_type == 3)
	{
	  plasmamain[n].t_r = proga_temp (x)/0.9; //This is a kluge, it means that we get the temperature we expect, if we read in a temperature from a zeus file - otherwise it is multiplies by 0.9 //
	}
      else
	{
	  plasmamain[n].t_r = geo.twind;
	}




      /* 70b - Initialize the temperature in the torus to a different value */
      if (w[nwind].inwind == W_ALL_INTORUS
	  || w[nwind].inwind == W_PART_INTORUS)
	{
	  plasmamain[n].t_r = geo.compton_torus_te;
	}

/* Initialize variables having to do with converence in initial stages */
      plasmamain[n].gain = 0.5;
      plasmamain[n].dt_e_old = 0.0;
      plasmamain[n].dt_e = 0.0;
      plasmamain[n].t_e = plasmamain[n].t_e_old = 0.9 * plasmamain[n].t_r;	//Lucy guess


/* Calculate an initial guess for the weight of the PL spectrum (constant / area of a sphere / 4pi) */


      rrstar =
	1. - (geo.rstar * geo.rstar) / (x[0] * x[0] + x[1] * x[1] +
					x[2] * x[2]);
      if (rrstar > 0)
	{
	  plasmamain[n].w = 0.5 * (1 - sqrt (rrstar));
	}
      else
	plasmamain[n].w = 0.5;	//Modification to allow for possibility that grid point is inside star

      /* Determine the initial ionizations, either LTE or  fixed_concentrations */
      if (geo.ioniz_mode != IONMODE_FIXED)
	{			/* Carry out an LTE determination of the ionization */
	  ierr = ion_abundances (&plasmamain[n], IONMODE_LTE);
	}
      else
	{			/* Set the concentrations to specified values */
	  ierr = ion_abundances (&plasmamain[n], IONMODE_FIXED);
	}
      if (ierr != 0)
	{
	  Error
	    ("wind_define after ion_abundances: cell %d rho %8.2e t_r %8.2 t_e %8.2e w %8.2e\n",
	     n, plasmamain[n].rho, plasmamain[n].t_r, plasmamain[n].t_e,
	     plasmamain[n].w);
	}

      /* 68b - Initialize the scatters array 73d - and the pariwise ionization denominator and temperature
       */

      for (j = 0; j < nions; j++) /* NSH 1107 - changed this loop to loop over nions rather than NIONS. Dynamic
	allocation means that these arrays are no longer of length NIONS */
	{
	  plasmamain[n].PWdenom[j] = 0.0;
	  plasmamain[n].PWnumer[j] = 0.0;
	  plasmamain[n].PWdtemp[j] = 0.0;
	  plasmamain[n].PWntemp[j] = 0.0;
	  plasmamain[n].scatters[j] = 0;
	  plasmamain[n].xscatters[j] = 0;
	}
    }


/* Calculate the the divergence of the wind at the center of each grid cell */
  wind_div_v (w);

/* Now calculate the adiabatic cooling.  Note: adiabatic cooling is not used in
 * the program at present.  There are issues concerning how to incorporate it
 * into the macro atom approach, as well as questions concerning the appropriate
 * calculation.  If changes are made to this, then they must also be made in
 * the corresponding portion of wind_updates.  04nov -- ksl
 */

/*06may -- ksl -- This is awkward because liminosities are now part of plasma structure */
  for (i = 0; i < NPLASMA; i++)
    {
      if (geo.adiabatic)
	{
	  nwind = plasmamain[i].nwind;
	  plasmamain[i].lum_adiabatic =
	    adiabatic_cooling (&w[nwind], plasmamain[i].t_e);
	}
      else
	plasmamain[i].lum_adiabatic = 0.0;
    }

  /* Calculate one over dvds */
  dvds_ave ();
  wind_check (w, -1);		// Check the wind for reasonability

  /* zero the counters which record diagnositcs from mean_intensity */
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
*/

  if (geo.coord_type == CYLIND)
    {
      mdotwind = 0;
      mdotbase = 0;
      for (i = 0; i < NDIM - 1; i++)
	{
	  n = i * MDIM + (MDIM / 2);
// rr is router * router - rinner *rinner for this shell
	  rr =
	    w[n + MDIM].x[0] * w[n + MDIM].x[0] + w[n +
						    MDIM].x[1] * w[n +
								   MDIM].x[1]
	    - (w[n].x[0] * w[n].x[0] + w[n].x[1] * w[n].x[1]);
	  if (w[i * MDIM].inwind == W_ALL_INWIND)
	    {
	      nplasma = w[i * MDIM].nplasma;
	      mdotbase +=
		plasmamain[nplasma].rho * PI * rr * w[i * MDIM].v[2];
	    }
	  if (w[n].inwind == W_ALL_INWIND)
	    {
	      nplasma = w[n].nplasma;
	      mdotwind += plasmamain[nplasma].rho * PI * rr * w[n].v[2];
	    }
	}

      //Factor of 2 to allow for wind on both sides of disk
      Log ("m dot wind: Desired %g   Base %g Calculated %g\n",
	   geo.wind_mdot, 2. * mdotbase, 2. * mdotwind);
      mdot_wind (w, 1.e6, geo.wind_rmax / 2.);
    }
  else
    {
      Error
	("wind2d.c: Don't know how to check wind mdot for this coordtype %d\n",
	 geo.coord_type);
    }

  return (0);
}


/***********************************************************
                                       Space Telescope Science Institute

 Synopsis:
 	where_in_grid locates the 1-d grid position of the photon. 

 Arguments:		
	double x[];
 Returns:
 	where_in_grid normally  returns the cell number associated with
 		a postion.  If the photon is in the grid this will be a positive
 		integer < NDIM*MDIM.
 	photon is inside the grid        -1
	photon is outside the grid       -2
 Description:	
	
		
 Notes:
	Where_in grid does not tell you whether the photon is in the wind or not. 

	What one means by inside or outside the grid may well be different
	for different coordinate systems.

	ksl--Not entirely clear why I don't make wig_x a vector.  Programming
	is not elegant at all.

 History:
 	97jan	ksl	Coding on python began.
	98jul	ksl	Removed WindPtr as one of arguments since not used
	98dec	ksl	Modified so that the argument is simply a position
	99jan	ksl	Modified to check that we are not trying to short
			circuit if we are asking for where_in_grid of same
			position as previously
	04aug	ksl	52a -- Now almos a shell, as moved to allow
			multiple coordinate systems
	05apr	ksl	55d -- Added spherical as a possiblity
 
**************************************************************/

int wig_n;
double wig_x, wig_y, wig_z;
int
where_in_grid (x)
     double x[];
{
  int n;
  double fx, fz;

  if (wig_x != x[0] || wig_y != x[1] || wig_z != x[2])	// Calculate if new position
    {

      if (geo.coord_type == CYLIND)
	{
	  n = cylind_where_in_grid (x);
	}
      else if (geo.coord_type == RTHETA)
	{
	  n = rtheta_where_in_grid (x);
	}
      else if (geo.coord_type == SPHERICAL)
	{
	  n = spherical_where_in_grid (x);
	}
      else if (geo.coord_type == CYLVAR)
	{
	  n = cylvar_where_in_grid (x, 0, &fx, &fz);
	}
      else
	{
	  Error ("where_in_grid: Unknown coord_type %d\n", geo.coord_type);
	  exit (0);
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

/***********************************************************
                                       Space Telescope Science Institute

 Synopsis:
	vwind_xyz(w,p,v) finds the velocity vector v for the wind in cartesian 
	coordinates at the position of the photon p.
 
 Arguments:		
	WindPtr w;
	PhotPtr p;
	double v[];

Returns:
	0 	on successful completion
		the value of where in grid if the photon outside the wind grid.
 
Description:
	This routine carries out a bilinear interpolation of the wind velocity.  The velocities
	at the edges of the grid cells must have been defined elsewhere (e.g in wind_define).
	
	If the position is outside the wind grid, then the velocities will be returned as 0.	
	
Notes:
	For all of the 2d coordinate systems, the positions in the WindPtr array are
	caclulated in the xz plane.  As a result, the velocities in that array are
	calcuated there as well.  Hence, to obtain the velocity in cartesion coordinates
	one needs to do a simple rotation of the velocity vector.  

	For spherical (1d) coordinates the situation is complex, the philosopy that
	has been adopted is to let the program run even if the wind model is intrinsically
	2d.   The xyz positions of the grid are defined are defined so that
	are in the xz plane at an angle of 45 degrees to the x axis.  This was done
	so that if one used, say an sv wind, in spherical coordinates, one would select
	the wind velocities at a plausible location.   But it implies that the velocities
	that have been accepted are not purely radial as one might expect.  The choice
	that has been made is to assume that v[1] is preserved as a non zero value in
	making the projection.  An alternative might have been to preserve the angular
	rotation rate.
	
	Regardless, the fundamental challenge is to make sure the program works correctly 
	for the spherically symmetric case, e.g a simple stellar wind, which it did not
	before 56b.  

	If we ever implement a real 3d coordinate system, one will need to look at
	this routine again.

History:
 	97jan	ksl	Coding on python began.
	02feb	ksl	Modified to use new general purpose fraction routine.  
	04aug	ksl	52a -- Made a shell for choosing what coordinate
			system to calculate this in
	05jun	ksl	56b -- Despite the note above nothing had been done
			to fix vwind_xyz to handle spherical coordinates
			properly.  It was OK for any two-d system.  This 
			is now fixed, but see the notes above.
	06may	ksl	57+ -- Changed call to elimatnate passing the
			Wind array.  Use wmain instead.
 
**************************************************************/
int ierr_vwind = 0;

int
vwind_xyz (p, v)
     PhotPtr p;
     double v[];
{
  int i;
  double rho, r;
  double vv[3];
  double ctheta, stheta;
  double x, frac[4];
  int nn, nnn[4], nelem;


  // gives the correct result in all coord systems.


  /* 56d -- 05jul -- I had considerable problem with the
   * next routine for cylvar coords.  Coord_fraction must
   * produce a plausible result in all cases
   */

  coord_fraction (0, p->x, nnn, frac, &nelem);



  for (i = 0; i < 3; i++)
    {

      x = 0;
      for (nn = 0; nn < nelem; nn++)
	x += wmain[nnn[nn]].v[i] * frac[nn];

      vv[i] = x;
    }

  rho = sqrt (p->x[0] * p->x[0] + p->x[1] * p->x[1]);

  if (geo.coord_type == SPHERICAL)
    {				// put the velocity into cylindrical coordinates on the xz axis
      x = sqrt (vv[0] * vv[0] + vv[2] * vv[2]);	//v_r
      r = length (p->x);
      vv[0] = rho / r * x;
      vv[2] = p->x[2] / r * x;
    }
  else if (p->x[2] < 0)		// For 2d coord syatems, velocity is reversed if the photon is in the lower hemisphere.
    vv[2] *= -1;


  if (rho == 0)
    {				// Then we will not be able to project from a cylindrical ot a cartesian system
      Error
	("vwind_xyz: Cannot determine an xyz velocity on z axis. Returnin 0,0,v[2]\n");
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


/***********************************************************
                                       Space Telescope Science Institute

 Synopsis:
  wind_div_v calculates the divergence of the velocity at the center of all the grid cells.  
Arguments:    

Returns:
 
Description:  
  This is one of the initialization routines for the wind.   
  Because the routine calls vwind_xyz, one must have previously 
  populated w[].v.  Also as a result of this fact, the routine 
  is not tightly tied to the SV prescription.
Notes:
  It is used to calculated PdV cooling of a grid cell 

  The divergence, like othe scalar quantities, does not
  need to be "rotated" differently in different coordinate
  systems

History:
  98mar20 ksl First coded and debugged 
  99dec1  ksl Explicitly typed as integer and made the return 0
      Previously it was untyped, but returned div/EPSILON
      which made no sense to me.
  04mar ksl This routine traditionally has generated
      lots of errors.  One reason for could have beeen
      that the step size was too large.  So I made this
      smaller.  But the other reason is that the div v
      is often calculated outside the wind cone, as in
      kwd winds.  In those cases, I put an extrapolation
      of the wind velocity inside the wind cone, but
      this is not necessarily physical.  If the div v
      is interpolated this could be wrong.  I have 
      made some minor modifications though to suppress
      multiple errors if the div is being calculated
      outside the wind.
  05apr ksl 55d -- Switched to using center positions for
      this calculation, rather than wind_midx, and 
      wind_midz.  Basically this was because I wanted
      to make things as independent of these two 
      arrays as possible in preparation for allowing
      more arbitrarily spaced grids.  Furthermore
      the previously formulatin was incorrect for
      rtheta components since wind_midz was actually
      theta.  It was also helpful for spherical models, 
      since it avoids the necessity of going from 1d to
      2d specification.
  13mar nsh 74b6 -- Several of the calls were addressing the
      wind as if it was a single cell, however the call
      to wind_div_v supplies the whole wind structure.
      calls that looked like w-> have been replaced by
      w[icell].
  14feb JM -- Fix for Issue #70. Instead of taking the midpoint
      and stepping ds in the +ve direction, we now step
      0.5ds in both +ve and -ve. This stops negative values of
      dvdy appearing in rotationally dominated portions of wind,
      and in some cases even leading to div_v being negative.
**************************************************************/

int wind_div_err = (-3);
int
wind_div_v (w)
     WindPtr w;
{
  int icell;
  double x_zero[3], v2[3], v1[3];
  struct photon ppp;
  double div, delta;
  double length ();
  double xxx[3];
  for (icell = 0; icell < NDIM2; icell++)
    {
      /* Find the center of the cell */

      /*   stuff_v (w->xcen, x_zero); OLD NSH 130322 - this line seems to assume w is a cell, rather than the whole wind structure */
      stuff_v (w[icell].xcen, x_zero);  /*NEW NSH 130322 - now gets the centre of the current cell in the loop */

      delta = 0.01 * x_zero[2]; //new 04mar ksl -- delta is the distance across which we measure e.g. dv_x/dx

      /* JM 1302 -- This has now been changed so instead of taking the midpoint and comparing to a point 
         delta in the +v direction, we now fo delta/2 in either direction. This is in order to correctly 
         evaluate dv/dy -- see bug report #70 */

      /* for each of x,y,z we first create a copy of the vector at the center. We then step 0.5*delta
         in positive and negative directions and evaluate the difference in velocities. Dividing this by
         delta gives the value of dv_x/dx, and the sum of these gives the divergence */


      /* Calculate dv_x/dx at this position */
      stuff_v (x_zero, ppp.x);   
      ppp.x[0] += 0.5 * delta;
      vwind_xyz (&ppp, v2);
      ppp.x[0] -= delta;
      vwind_xyz (&ppp, v1);
      div = xxx[0] = (v2[0] - v1[0]) / delta;

      /* Calculate dv_y/dy */
      stuff_v (x_zero, ppp.x);    
      ppp.x[1] += 0.5 * delta;
      vwind_xyz (&ppp, v2);
      ppp.x[1] -= delta;
      vwind_xyz (&ppp, v1);
      div += xxx[1] = (v2[1] - v1[1]) / delta;


      /* Calculate dv_z/dz */
      stuff_v (x_zero, ppp.x);  
      ppp.x[2] += 0.5 * delta;
      vwind_xyz (&ppp, v2);
      ppp.x[2] -= delta;
      vwind_xyz (&ppp, v1);
      div += xxx[2] = (v2[2] - v1[2]) / delta;


      /* we have now evaluated the divergence, so can store in the wind pointer */
      w[icell].div_v = div;


      if (div < 0 && (wind_div_err < 0 || w[icell].inwind == W_ALL_INWIND)) /*NSH 130322 another fix needed here the inwind check was w->inwind and was returning the wrong value */
  {
    Error ("wind_div_v: div v %e is negative in cell %d. Major problem if inwind (%d) == 0\n", div, icell, w[icell].inwind);  /*NSH 130222 - last fix */
    wind_div_err++;
  }
    }

  return (0);
}


/* 
find the density of the wind at x

History:
	02feb	ksl Changed the calculation to use the general purpose routine fraction.
			This made a slight difference, an improvement I believe in how
			the program performs at the boundaries.  In the new version
			if you ask for rho at a position that is below say wind_midz
			one obtains the density at the center of the cell. In the old
			version it tried to extrapolate.  As far as I can determine
			the difference was relatively modest.
	04aug	ksl	52a -- modified to carry out a coordinate system independend
			determination of rho.  
	06may	ksl	57+ -- Modified for plasma structure but this is probably not 
			what one wants in the end.  There is no reason for the call
			to include w at all, since it is not used..

*/

double
rho (w, x)
     WindPtr w;
     double x[];
{
  int n, where_in_grid ();
  double dd;
  double frac[4];
  int nn, nnn[4], nelem;
  int nplasma;


  if (where_in_wind (x) != 0)	//note that where_in_wind is independent of grid.
    return (0.0);

  n = coord_fraction (1, x, nnn, frac, &nelem);

  if (n < 0)
    {
      dd = 0;
    }
  else
    {

      dd = 0;
      //59a - ksl - 070823 - fiexed obvious error that has been there since
      //I split the structures into w and plasmamain.
      for (nn = 0; nn < nelem; nn++)
	{
	  nplasma = w[nnn[nn]].nplasma;
	  dd += plasmamain[nplasma].rho * frac[nn];
	}

    }


  return (dd);
}

/* Next is a routine to check the wind mass loss rate.  The mdot is checked in a plane running
from 0 to r, and in a sphere of radius r 

04aug -- ksl -- this looks coordinate system independnet but need to check
*/

#define NSTEPS 100
int
mdot_wind (w, z, rmax)
     WindPtr w;
     double z;			// The height (usually small) above the disk at which mdot will be calculated
     double rmax;		// The radius at which mdot will be calculated
{
  struct photon p;		//needed because vwind_xyz has a call which involves a photon
  double r, dr, rmin;
  double theta, dtheta;
  double den, rho ();
  double mdot, mplane, msphere;
  double x[3], v[3], q[3], dot ();
// Calculate the mass loss rate immediately above the disk
  rmin = geo.rstar;
  dr = (rmax - rmin) / NSTEPS;
  mdot = 0;
  p.x[1] = x[1] = 0.0;
  p.x[2] = x[2] = z;
  for (r = dr / 2; r < rmax; r += dr)
    {
      p.x[0] = x[0] = r;
      den = rho (w, x);
      vwind_xyz (&p, v);
      mdot += 2 * PI * r * dr * den * v[2];
//mdot+=2*PI*r*dr;
    }
  mplane = 2 * mdot;		// Factor of two because wind is in both hemispheres
// Calculate the mass loss rate in a sphere
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
      vwind_xyz (&p, v);
      mdot += 2 * PI * rmax * rmax * sin (theta) * dtheta * den * dot (v, q);
    }
  msphere = 2 * mdot;		// Factor of two because wind is in both hemispheres
  Log ("Wind Mdot Desired %e Plane %e Sphere(%e) %e\n",
       geo.wind_mdot, mplane, rmax, msphere);
  return (0);
}


/***********************************************************
                                       Space Telescope Science Institute

 Synopsis:
	get_random_location is simply will produce a at a random place in
	a cell n from component icomp.  
Arguments:		

Returns:
 
Description:	
	This is really just a driver routine for coordinate system specific
	routines
Notes:

	The reason you need the component here is because the boundary
	of the component may cut through the cell,and you don't want
	to generate a position outside of the component.

History:
	04aug	ksl	52a -- created as part of project to allow multiple
			coordinate systems in python
	05apr	ksl	55d -- Added spherical option
	11aug	ksl	70b -- Added option of getting a random location in
			the torus, or any new component

 
**************************************************************/

int
get_random_location (n, icomp, x)
     int n;			// Cell in which to create postion
     int icomp;			// The component we want the position in
     double x[];		// Returned position
{

  if (geo.coord_type == CYLIND)
    {
      cylind_get_random_location (n, icomp, x);
    }
  else if (geo.coord_type == RTHETA)
    {
      rtheta_get_random_location (n, icomp, x);
    }
  else if (geo.coord_type == SPHERICAL)
    {
      spherical_get_random_location (n, icomp, x);
    }
  else if (geo.coord_type == CYLVAR)
    {
      cylvar_get_random_location (n, icomp, x);
    }
  else
    {
      Error ("get_random_location: Don't know this coord_type %d\n",
	     geo.coord_type);
      exit (0);
    }

  return (0);
}


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


/* The next routine checks how many corners of a wind cell
 * are in the wind
  
 This routine returns the number of corners of a wind cell
 that are in the wind.  It is intended to standardize this
 check so that one can use such a routin in the volumes
 calculations.  

Notes:

	It has not been implemented for a spherical (1d)
	coordinate system.

	The fact that none of the corners of a cell are in
	the wind is not necessarily a guarantee that the wind
	does not pass through the cell.  This is not only
	a theoretical problem, because we generally use
	a logarithmic spacing for the wind cells and thus
	the dimensions of the cells get quite large as one
	gets far from the center.
	
	The only way to verify this is to check all four 
	surfaces of the wind.
  
History
	080403	ksl	68c - Added, or more correctly simply
			moved code into a separate routine so 
			could be called from elsewhere
 */

int
check_corners_inwind (n, icomp)
     int n;
     int icomp;			// check corners for this component
{
  int n_inwind;
  int i, j;

  wind_n_to_ij (n, &i, &j);

  n_inwind = 0;
  if (i < (NDIM - 2) && j < (MDIM - 2))
    {
      if (where_in_wind (wmain[n].x) == icomp)
	n_inwind++;
      if (where_in_wind (wmain[n + 1].x) == icomp)
	n_inwind++;
      if (where_in_wind (wmain[n + MDIM].x) == icomp)
	n_inwind++;
      if (where_in_wind (wmain[n + MDIM + 1].x) == icomp)
	n_inwind++;
    }
  return (n_inwind);
}
