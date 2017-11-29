#define MAXHYDRO 1000
#define IGHOST 0

double hydro_r_cent[MAXHYDRO];
double hydro_r_edge[MAXHYDRO];
double hydro_theta_cent[MAXHYDRO];
double hydro_theta_edge[MAXHYDRO];
int ihydro_r, ihydro_theta, j_hydro_thetamax, ihydro_mod;
double hydro_thetamax;		//The angle at which we want to truncate the theta grid
double v_r_input[MAXHYDRO * MAXHYDRO];
double v_theta_input[MAXHYDRO * MAXHYDRO];
double v_phi_input[MAXHYDRO * MAXHYDRO];
double rho_input[MAXHYDRO * MAXHYDRO];
double temp_input[MAXHYDRO * MAXHYDRO];


/*
typedef struct hydro_mod
{
  double rho;
  double v[3];
  double temp;
}
hydro_dummy, *HydroPtr;

HydroPtr hydro_ptr;
*/




/**************************************************************************
                    Space Telescope Science Institute


  Synopsis:  These routines are designed to read in files produced by
  zeus for the x-ray binary project.
  models onto sv

  Description:  These routines are designed to read the velocity and densities
  of produced by zeus  into the appropriate structures and
  allow one to calculate the density at any point in the gridded space.
  

  Arguments:    

  Returns:

  Notes:
  Daniel's models data are the results of calculations with Zeus.  The
  grid is in polar (r,theta) coordinates with theta = pi/2 corresponding
  to the disk, and 0 to the pole.  Recall that python uses cylintrical
  corrdinates with the z=0 being the disk plane.

  History:
  99dec ksl Began work
  00jan ksl Modified to read larger models and to make it slightly
      more general in terms of the input files.  Also added
      some checking on grid sizes, and set rho to zero outside
      the range of the model, and velocity at the edge of
      the model.
  04jun ksl Moved get_hydro_wind_params from python.c to this file
  13mar nsh Reopened development prior to trip to LV to work with DP.
  15aug ksl Began updates to accommodate domains

  

 ************************************************************************/
#include  <stdio.h>
#include  <stdlib.h>
#include  <strings.h>
#include  <string.h>
#include  <math.h>
#include  "log.h"
#include  "atomic.h"
#include  "python.h"

#define LINE 200




/***********************************************************
                                       Space Telescope Science Institute

 Synopsis:
  get_hydro_wind_params gets input data for Daniel Proga's wind models
Arguments:    

Returns:
 
Description:  
Notes:
History:
  99dec ksl Began work

  13mar nsh Reopened work on this file to coincide with 
      NSH working with Proga for a week.
    nsh First edit - DP says that the second value
      in the two grid files relate to the position
      where the data is defined.
    nsh Second edit - put in a new variable, hydro_hydro_thetamax.
      This defines an angle above which one disregards data.
      This was needed because Daniels data includes a flared
      accretion disk, so if you simply read in all the data
      you end up with a very dense blancket over the top of our
      disk, which causes all kinds of problems. At present,
      all that happens is that all angle data above the 
      supplied angle is replaced with the data for that last 
      angle. This should be a very small wedge of data.
  13may nsh Now read in the inteernal energy - this allows the
      computation of temperature for a cell, this makes sense
      as a first guess of temperature 
  15sept  ksl Made modifications to allow for domains, but Nick will
      need to debug this
  
**************************************************************/


int
get_hydro_wind_params (ndom)
     int ndom;
{
  Log ("Creating a wind model using a Hydro calculation = domain %i\n", ndom);

  get_hydro (ndom);

  /* updates_2d has special code related to zeus operations with zeus. The
   * next line preserves information related to the domain number
   */
  geo.hydro_domain_number = ndom;

/* Assign the generic parameters for the wind the generic parameters of the wind */
  geo.rmin = zdom[ndom].rmin = hydro_r_edge[0];

  geo.rmax = zdom[ndom].rmax = hydro_r_edge[ihydro_r] + 2.0 * (hydro_r_cent[ihydro_r] - hydro_r_edge[ihydro_r]);	//Set the outer edge of the wind to the outer edge of the final defined cell
  Log ("rmax=%e\n", geo.rmax);
  geo.rmax_sq = geo.rmax * geo.rmax;
  Log ("rmax_sq=%e\n", geo.rmax);
  zdom[ndom].wind_rho_min = 0.0;	//Set wind_rmin 0.0, otherwise wind cones dont work properly 
  Log ("rho_min=%e\n", zdom[ndom].wind_rho_min);
  zdom[ndom].wind_rho_max = zdom[ndom].rmax;	//This is the outer edge of the
  Log ("rho_max=%e\n", zdom[ndom].wind_rho_max);
  zdom[ndom].zmax = zdom[ndom].rmax;	//This is the outer edge of the
  Log ("zmax=%e\n", zdom[ndom].zmax);

  zdom[ndom].wind_thetamin = hydro_theta_edge[0];
  Log ("theta_min=%e\n", zdom[ndom].wind_thetamin);
  Log ("theta_max=%e\n", zdom[ndom].wind_thetamax);
  Log ("geo.rmin=%e\n", zdom[ndom].rmin);
  Log ("geo.rmax=%e\n", zdom[ndom].rmax);
  Log ("geo.wind_rhomin=%e\n", zdom[ndom].wind_rho_min);
  Log ("geo.wind_rhomax=%e\n", zdom[ndom].wind_rho_max);



  /* if modes.adjust_grid is 1 then we have already adjusted the grid manually */
  if (modes.adjust_grid == 0)
    {
      zdom[ndom].xlog_scale = 0.3 * geo.rstar;
      zdom[ndom].zlog_scale = 0.3 * geo.rstar;
    }

  return (0);
}



int
get_hydro (ndom)
     int ndom;
{

  FILE *fopen (), *fptr;
  char datafile[LINE];
  char aline[LINE];
  char word[LINE];
  int i, j, k;
  double r, r_edge;
  double rho;
  double theta, theta_edge, temp;
  double vr, vtheta, vphi;
  int irmax, ithetamax, itest;
  int ndim, mdim;

/*Write something into the file name strings */

  strcpy (datafile, "hdf062.dat");


  for (k = 0; k < MAXHYDRO; k++)
    {
      hydro_r_cent[k] = 0;
      hydro_theta_cent[k] = 0;
    }

  rdstr ("hydro_file", datafile);
  if ((fptr = fopen (datafile, "r")) == NULL)
    {
      Error ("Could not open %s\n", datafile);
      exit (0);
    }


  hydro_thetamax = 89.9;


  rddoub ("Hydro_thetamax(degrees:negative_means_no_maximum)",
	  &hydro_thetamax);

  //If we have set the maximum angle to a negaive value, we mean that we dont want to restrict.
  if (hydro_thetamax < 0.0)
    hydro_thetamax = VERY_BIG;
  hydro_thetamax = hydro_thetamax / RADIAN;



  ihydro_r = 0;
  j_hydro_thetamax = 0;		/* NSH 130605 to remove o3 compile error */
  irmax = 0;
  ithetamax = 0;		//Zero counters to check all cells actually have data

  while (fgets (aline, LINE, fptr) != NULL)
    {

      if (aline[0] != '#')
	{
	  sscanf (aline, "%s", word);

	  if (strncmp (word, "ir", 2) == 0)
	    {
	      Log
		("We have a hydro title line, we might do something with this in the future \n");
	    }
	  else
	    {
	      itest =
		sscanf (aline, "%d %lf %lf %d %lf %lf %lf %lf %lf %lf %lf",
			&i, &r, &r_edge, &j, &theta, &theta_edge, &vr,
			&vtheta, &vphi, &rho, &temp);
	      if (itest != 11)	//We have an line which does not match what we expect, so quit
		{
		  Error ("hydro.c data file improperly formatted\n");
		  exit (0);
		}
	      // read read the r and theta coordinates into arrays
	      hydro_r_edge[i] = r_edge;
	      hydro_r_cent[i] = r;
	      hydro_theta_cent[j] = theta;
	      hydro_theta_edge[j] = theta_edge;
	      //keep track of how many r and theta cells there are            
	      if (j > ithetamax)
		ithetamax = j;
	      if (i > irmax)
		irmax = i;
	      //If the value of theta in this cell, the edge, is greater than out theta_max, we want to make a note.
	      if (hydro_theta_edge[j] > hydro_thetamax
		  && hydro_theta_edge[j - 1] <= hydro_thetamax)
		{
		  j_hydro_thetamax = j - 1;
		  Log
		    ("current theta  (%f) > theta_max  (%f) so setting j_hydro_thetamax=%i\n",
		     theta * RADIAN, hydro_thetamax * RADIAN,
		     j_hydro_thetamax);
		}
	      //If theta is greater than thetamax, then we will replace rho with the last density above the disk
	      else if (hydro_theta_edge[j] > hydro_thetamax)
		{
		  /* NSH 130327 - for the time being, if theta is in the disk, replace with the last
		     density above the disk */
		  rho = rho_input[i * MAXHYDRO + j_hydro_thetamax];
		}
	      rho_input[i * MAXHYDRO + j] = rho;
	      temp_input[i * MAXHYDRO + j] = temp;
	      v_r_input[i * MAXHYDRO + j] = vr;
	      v_theta_input[i * MAXHYDRO + j] = vtheta;
	      v_phi_input[i * MAXHYDRO + j] = vphi;
	    }
	}
    }



  if (j_hydro_thetamax == 0 || j_hydro_thetamax == i - 1)
    {
      Log ("HYDRO j_hydro_thetamax never bracketed, using all data\n");
      ihydro_theta = ithetamax;
      zdom[ndom].wind_thetamax = 90. / RADIAN;
      hydro_thetamax = 90.0 / RADIAN;
      mdim = zdom[ndom].mdim = ihydro_theta + 2;
    }
  else
    {
      Log
	("HYDRO j_hydro_thetamax=%i, bracketing cells have theta = %f and %f\n",
	 j_hydro_thetamax, hydro_theta_cent[j_hydro_thetamax] * RADIAN,
	 hydro_theta_cent[j_hydro_thetamax + 1] * RADIAN);
      ihydro_theta = j_hydro_thetamax;
      zdom[ndom].wind_thetamax = hydro_thetamax;
      mdim = zdom[ndom].mdim = ihydro_theta + 2;
    }




  if (hydro_r_edge[0] < geo.rstar)
    {
      Error
	("Major problem, innermost edge of hydro radial grid begins inside geo.rstar\n");
      exit (0);
    }

  ihydro_r = irmax;

  Log ("Read %d r values\n", ihydro_r);
  Log ("Read %d theta values\n", ihydro_theta);
  fclose (fptr);

  /* Set a couple of last tags */

  zdom[ndom].coord_type = RTHETA;	//At the moment we only deal with RTHETA - in the future we might want to do some clever stuff
  ndim = zdom[ndom].ndim = ihydro_r + 3;	//We need an inner radial cell to bridge the star and the inside of the wind, and an outer cell
  zdom[ndom].ndim2 = zdom[ndom].ndim * zdom[ndom].mdim;	// Make ndim2 consistent with the individual dimensions


  return (0);
}



/***********************************************************
                                       Space Telescope Science Institute

 Synopsis:
  hydro_velocity calculates the wind velocity at any position x
  in cartesian coordinates
Arguments:    

Returns:
 
Description:  
Notes:
History:
  99dec ksl Began work
  04dec ksl 52a -- Mod to prevent divide by zero when 
      r=0.0, and NaN when xxx>1.0;
**************************************************************/

// ksl - Added ndom here for symmetry with other modles.  
// The hydro models do not understand domains XXX

double
hydro_velocity (ndom,x, v)
    int ndom;
     double x[];
     double v[];
{
  double length ();
  int ii, jj;
  int im, jm;
  double f1, f2;
  double r, theta;
  double v_r_interp, v_theta_interp, v_phi_interp;
  double speed;
  double xxx;

  if ((r = length (x)) == 0.0)
    {
      v[0] = v[1] = v[2] = 0.0;
      return (0.);
    };



  xxx = (sqrt (x[0] * x[0] + x[1] * x[1]) / r);
  if (xxx >= 1.)
    {
      theta = PI / 2.;
    }
  else
    theta = asin (xxx);


  /* We compute the radial velocity - in zeus this is defined at the 
  cell centre in the theta dimension and the edge in the r dimension */


  im = jm = ii = jj = 0;
  f1 = f2 = 0.0;
  hydro_frac (r, hydro_r_edge, ihydro_r, &im, &ii, &f1);
  hydro_frac (theta, hydro_theta_cent, ihydro_theta, &jm, &jj, &f2);
  v_r_interp = hydro_interp_value (v_r_input, im, ii, jm, jj, f1, f2);
  
  
  /* We compute the theta velocity - in zeus this is defined at the 
  cell centre in the r dimension and the edge in the theta dimension */

  im = jm = ii = jj = 0;
  f1 = f2 = 0.0;
  hydro_frac (r, hydro_r_cent, ihydro_r, &im, &ii, &f1);
  hydro_frac (theta, hydro_theta_edge, ihydro_theta, &jm, &jj, &f2);
  v_theta_interp = hydro_interp_value (v_theta_input, im, ii, jm, jj, f1, f2);
  
  
  /* We compute the phi velocity - in zeus this is defined at the 
  cell centre in both the r and theta dimension */
  

  im = jm = ii = jj = 0;
  f1 = f2 = 0.0;
  hydro_frac (r, hydro_r_cent, ihydro_r, &im, &ii, &f1);
  hydro_frac (theta, hydro_theta_cent, ihydro_theta, &jm, &jj, &f2);
  v_phi_interp = hydro_interp_value (v_phi_input, im, ii, jm, jj, f1, f2);
  
  
  /* Convert to carteisian velocity */


  v[0] = v_r_interp * sin (theta) + v_theta_interp * cos (theta);
  v[1] = v_phi_interp;
  v[2] = v_r_interp * cos (theta) - v_theta_interp * sin (theta);
  
  
  /*Compute the speed */

  speed = sqrt (v[0] * v[0] + v[1] * v[1] * v[2] * v[2]);
  
  

  if (sane_check (speed))
    {
      Error ("hydro_velocity:sane_check v %e %e %e\n", v[0], v[1], v[2]);
    }
  return (speed);

}



double
hydro_rho (x)
     double x[];
{
  double length ();
  int ii, jj;
  int im, jm;
  double r, theta;
  double rrho;
  double f1, f2;
  r = length (x);
  theta = asin (sqrt (x[0] * x[0] + x[1] * x[1]) / r);

  if ((hydro_r_cent[ihydro_r] - r) / r < -1e-6)
    {
      Log (" r outside hydro grid in hydro_rho %e > %e %e\n", r,
	   hydro_r_cent[ihydro_r], (hydro_r_cent[ihydro_r] - r) / r);
      rrho = 1.e-23;
      return (rrho);
    }


  /* The density is a cell centred quantity in zeus */

  im = jm = ii = jj = 0;
  f1 = f2 = 0.0;
  hydro_frac (r, hydro_r_cent, ihydro_r, &im, &ii, &f1);
  hydro_frac (theta, hydro_theta_cent, ihydro_theta, &jm, &jj, &f2);

  rrho = hydro_interp_value (rho_input, im, ii, jm, jj, f1, f2);

  if (rrho < 1e-23)
    rrho = 1e-23;

  return (rrho);
}


/***********************************************************
                                       University of Nevada Las Vegas

 Synopsis:
  hydro_temp calculates the wind temperature at any position x
  in cartesian coordiantes
Arguments:    

Returns:
 
Description:  
  This code is an exact copy of hydro_rho - except it
  maps the temperature, calculated from internal energy,
  from the hydro grid onto a cartesian grid
Notes:
History:
  13apr nsh Began work
  15oct   nsh - temperature is now computed by a helper script,
        also now common code removed to subroutines
  
**************************************************************/


double
hydro_temp (x)
     double x[];
{
  double length ();
  int ii, jj;
  int im, jm;
  double r, theta, temp;
  double f1, f2;
  r = length (x);
  theta = asin (sqrt (x[0] * x[0] + x[1] * x[1]) / r);


  im = jm = ii = jj = 0;
  f1 = f2 = 0.0;


  if ((hydro_r_cent[ihydro_r] - r) / r < -1e-6)
  {
    Log (" r outside hydro grid in hydro_temp %e > %e %e\n", r, hydro_r_cent[ihydro_r], (hydro_r_cent[ihydro_r] - r) / r);
    temp = 1e4;
    return (temp);
  }
  
  
  /* The density is a cell centred quantity in zeus */
  


  hydro_frac (r, hydro_r_cent, ihydro_r, &im, &ii, &f1);
  hydro_frac (theta, hydro_theta_cent, ihydro_theta, &jm, &jj, &f2);

  temp = hydro_interp_value (temp_input, im, ii, jm, jj, f1, f2);


  /* NSH 16/2/29 - removed lower limit - set this in the hydro translation software or hydro model */

//  if (temp < 1e4)             //Set a lower limit.
//    temp = 1e4;


  return (temp);
}




/***********************************************************
                                       Southampton

 Synopsis:
  rtheta_make_zeus_grid defines the cells in a rtheta grid based upon the coordinates that can be read in from a zeus (the hydrocode used by Proga for some simulations)            

Arguments:    
  WindPtr w;  The structure which defines the wind in Python
 
Returns:
 
Description:

  
  This is an attempt to match a zeus grid directly onto an rtheta grid.


History:
  13jun nsh 76 -- Coded and debugged.


**************************************************************/


int
rtheta_make_hydro_grid (w, ndom)
     WindPtr w;
     int ndom;
{
  double theta, thetacen, dtheta;
  int i, j, n;
  int ndim, mdim;

  ndim = zdom[ndom].ndim;
  mdim = zdom[ndom].mdim;


  for (i = 0; i < ndim; i++)
    {
      for (j = 0; j < mdim; j++)
	{
	  wind_ij_to_n (ndom, i, j, &n);
	  w[n].inwind = W_ALL_INWIND;
	  if (i == 0)		// The inner edge of the grid should be geo.rstar
	    {
	      w[n].r = geo.rstar;	//So we set the inner edge to be the stellar (or QSO) radius
	      w[n].rcen = (geo.rstar + hydro_r_edge[0]) / 2.0;	//It will be a big cell
	      w[n].inwind = W_NOT_INWIND;
	    }
	  else if (i - 1 > ihydro_r)	// We are at the radial limit of the data, this last cell will be a ghost cell of our own
	    {
	      w[n].r = geo.rmax;	// so we set the outer cell to be the edge of the wind
	      w[n].rcen = geo.rmax;	// And it has zero volume
	      w[n].inwind = W_NOT_INWIND;
	    }
	  else
	    {
	      w[n].r = hydro_r_edge[i - 1];	//The -1 is to take account of the fact that i=0 cell is the inner part of the geometry inside the wind
	      w[n].rcen = hydro_r_cent[i - 1];
	    }
	  dtheta = 2.0 * (hydro_theta_cent[j] - hydro_theta_edge[j]);



	  if (hydro_theta_cent[j] + (dtheta / 2.0) > hydro_thetamax && hydro_theta_cent[j] - (dtheta / 2.0) < hydro_thetamax)	//This cell bridges the boundary - we reset it so the lower edge is at the disk boundary
	    {
	      theta = hydro_theta_edge[j];
	      thetacen = ((hydro_thetamax + theta) / 2.0);
	    }

	  else if (j > ihydro_theta)	//We are setting up a cell past where there is any data
	    {
	      thetacen = hydro_thetamax;	//Set the center and the edge to the maximum extent of the data/interest
	      theta = hydro_thetamax;
	      w[n].inwind = W_NOT_INWIND;
	    }
	  else
	    {
	      thetacen = hydro_theta_cent[j];
	      theta = hydro_theta_edge[j];
	    }
	  w[n].theta = theta * RADIAN;
	  w[n].thetacen = thetacen * RADIAN;
	  w[n].x[1] = w[n].xcen[1] = 0.0;
	  w[n].x[0] = w[n].r * sin (theta);
	  w[n].x[2] = w[n].r * cos (theta);
	  w[n].xcen[0] = w[n].rcen * sin (thetacen);
	  w[n].xcen[2] = w[n].rcen * cos (thetacen);

	}
    }

  /* Now set up the wind cones that are needed for calclating ds in a cell */

  rtheta_make_cones (ndom, w);	//NSH 130821 broken out into a seperate routine




  /* OK finished successfuly */
  return (0);

}

/***********************************************************
                                       Southampton

 Synopsis:
  rtheta_zeus_volumes replaces rtheta_volumes for a zeus model. We know wether cells are in the wimnd
  so all we need to do is work out the volumes.
Arguments:    
  WindPtr w;  The structure which defines the wind in Python
 
Returns:
 
Description:

  
  This is an attempt to match a zeus grid directly onto an rtheta grid.


History:
  13jun nsh 76 -- Coded and debugged.
  15aug ksl Updated for domains


**************************************************************/


int
rtheta_hydro_volumes (ndom, w)
     int ndom;
     WindPtr w;
{
  int i, j, n;
  double rmin, rmax, thetamin, thetamax;
  DomainPtr one_dom;

  one_dom = &zdom[ndom];

  for (i = 0; i < one_dom->ndim; i++)
    {
      for (j = 0; j < one_dom->mdim; j++)
	{

	  wind_ij_to_n (ndom, i, j, &n);
	  if (w[n].inwind == W_ALL_INWIND)
	    {

	      rmin = one_dom->wind_x[i];
	      rmax = one_dom->wind_x[i + 1];
	      thetamin = one_dom->wind_z[j] / RADIAN;
	      thetamax = one_dom->wind_z[j + 1] / RADIAN;

	      //leading factor of 2 added to allow for volume above and below plane (SSMay04)
	      w[n].vol =
		2. * 2. / 3. * PI * (rmax * rmax * rmax -
				     rmin * rmin * rmin) * (cos (thetamin) -
							    cos (thetamax));

	      if (w[n].vol == 0.0)
		{
		  Log
		    ("Found wind cell (%i) with no volume (%e) in wind, resetting\n",
		     n, w[n].vol);
		  w[n].inwind = W_NOT_INWIND;
		}
	    }
	  else
	    w[n].vol = 0.0;
	}
    }

  return (0);
}


/***********************************************************
                                       Southampton

 Synopsis:
  hydro_frac replaces many lines of identical code, all of which
  interpolate on the input grid to get value for the python grid.
  This code find out the fraction along a coord array where the centre
  or edge of the python grid cell lies on the hydro grid
Arguments:    
    double coord;
    double coord_array[];
    int imax;
    int *cell1,*cell2;
    double *frac;
Returns:
 
Description:

XXX - ksl - It is very unclear why this routine is necessary
coor_frac should be able to do this in one go.


History:
  15sept  nsh 76 -- Coded and debugged.


**************************************************************/




int
hydro_frac (coord, coord_array, imax, cell1, cell2, frac)
     double coord;
     double coord_array[];
     int imax;
     int *cell1, *cell2;
     double *frac;
{
  int ii;
  ii = 0;
  *cell1 = 0;
  *cell2 = 0;

  while (coord_array[ii] < coord && ii < imax)	//Search through array, until array value is greater than your coordinate
    ii++;


  if (ii > imax)
    {				// r is greater than anything in Proga's model
      *frac = 1;		//If we are outside the model set fractional position to 1
      *cell1 = *cell2 = imax;	//And the bin to the outermost
      return (0);
    }
  else if (ii == 0)
    {
      *frac = 1;		//Otherwise, we must be inside the innermost bin, so again, set fration to 1
      *cell1 = *cell2 = 0;	// And the bin to the innermost. Lines below to the same for theta.
      return (0);
    }
  else if (ii > 0)
    {				//r is in the normal range

      *frac = (coord - coord_array[ii - 1]) / (coord_array[ii] - coord_array[ii - 1]);	//Work out fractional position in the ii-1th radial bin where you want to be
      *cell1 = ii - 1;		//This is the radial bin below your value of r
      *cell2 = ii;
      return (0);
    }

  return (0);
}

/***********************************************************
                                       Southampton

 Synopsis:
  hydro_interp_value replaces many lines of identical code, all of which
  interpolate on the input grid to get value for the python grid.
Arguments:    
    double array[]; - the arrayof input data
    int im,ii; //the two cells surrounding the cell in the first dim (r)
    int jm,jj; //the two cells surrounding the cell in the second dim (theta)
    double f1,f2;  //the fraction between the two values in first and second dim
Returns:
 
Description:

XXX - ksl - It is not obvious why this was really needed, assuming that
the model has been proper installed in the WindStruct.  Possibly it
has to do with where values are centered.  


History:
  15sept  nsh 76 -- Coded and debugged.


**************************************************************/



double
hydro_interp_value (array, im, ii, jm, jj, f1, f2)
     double array[];
     int im, ii;		//the two cells surrounding the cell in the first dim (r)
     int jm, jj;		//the two cells surrounding the cell in the second dim (theta)
     double f1, f2;		//the fraction between the two values in first and second dim
{
  double value;
  double d1, d2;


  d1 =
    array[im * MAXHYDRO + jm] + f1 * (array[ii * MAXHYDRO + jm] -
				      array[im * MAXHYDRO + jm]);
  d2 =
    array[im * MAXHYDRO + jj] + f1 * (array[ii * MAXHYDRO + jj] -
				      array[im * MAXHYDRO + jj]);
  value = d1 + f2 * (d2 - d1);


  return (value);
}



/***********************************************************
                                       Southampton

 Synopsis:
	hydro_restart is a subrotine which permits a previous wind save 
		file to be used in a hydro simulation. The density and temperature
		for each cell are those from the hydro simulation. The ion fractions
		are taken from the windsave file
Arguments:		
	none 
 
Returns:
 
Description:

	
	This sets up a restarted hydro model


History:
	16feb	nsh	80 -- Coded and debugged.


**************************************************************/





int
hydro_restart (ndom)
     int ndom;
{
  int n, nion;
  int nwind;
  double x[3];
  double old_density;
  int nstart, nstop, ndim2;
  WindPtr w;

  w = wmain;
  zdom[ndom].wind_type = 3;	//Temporarily set the wind type to hydro, so we can use the normal routines

  /* note that we will have passed the wind domain number as default */
  nstart = zdom[ndom].nstart;
  nstop = zdom[ndom].nstop;
  ndim2 = zdom[ndom].ndim2;

  for (n = nstart; n < nstop; n++)
    {
      /* 04aug -- ksl -52 -- The next couple of lines are part of the changes
       * made in the program to allow more that one coordinate system in python 
       */

      model_velocity (ndom, w[n].x, w[n].v);
      model_vgrad (ndom, w[n].x, w[n].v_grad);
    }

  /* JM XXX PLACEHOLDER -- unsure how we loop over the plasma cells just in one domain 
   * ksl This is an error clearly XXX */
  for (n = 0; n < NPLASMA; n++)
    {
      nwind = plasmamain[n].nwind;
      stuff_v (w[nwind].xcen, x);
      old_density = plasmamain[n].rho;
      plasmamain[n].rho = model_rho (ndom, x) / zdom[ndom].fill;
      plasmamain[n].t_r = plasmamain[n].t_e = hydro_temp (x);

      for (nion = 0; nion < nions; nion++)	//Change the absolute number densities, fractions remain the same
	{
	  plasmamain[n].density[nion] =
	    plasmamain[n].density[nion] * (plasmamain[n].rho / old_density);
	}

      plasmamain[n].ne = get_ne (plasmamain[n].density);	//get the new electron density
      partition_functions (&plasmamain[n], NEBULARMODE_LTE_GROUND);	//ensure the partition functions and level densities are correct XXX ksl why this mode

    }

  /* XXX what is going on here. this looks to be outside the grid altgogether */
  plasmamain[n].ne = get_ne (plasmamain[n].density);	//we also need to update the electron density
  partition_functions (&plasmamain[n], NEBULARMODE_LTE_GROUND);	/* WARNING fudge NSH 11/5/14 - this is as a test. We really need a better implementation
								   of partition functions and levels for a power law illuminating spectrum. We found that
								   if we didnt make this call, we would end up with undefined levels - which did really
								   crazy things */


  return (0);

}
