#define MAXPROGA 155
#define IGHOST 3 //The number of ghost cells in zeus.


double proga_r_cent[MAXPROGA];
double proga_r_edge[MAXPROGA];
double proga_dr_cent[MAXPROGA];
double proga_dr_edge[MAXPROGA];
double proga_theta_cent[MAXPROGA];
double proga_theta_edge[MAXPROGA];
double proga_dtheta_cent[MAXPROGA];
double proga_dtheta_edge[MAXPROGA];
int iproga_r, iproga_theta, j_proga_thetamax, iproga_mod;
double proga_thetamax; //The angle at which we want to truncate the theta grid

typedef struct proga_mod
{
  double rho;
  double v[3];
  double temp;
}
proga_dummy, *ProgaPtr;

ProgaPtr proga_ptr;





/**************************************************************************
                    Space Telescope Science Institute


  Synopsis:  These routines are designed to match Daniel Proga's wind
	models onto sv

  Description:	These routines are designed to read the velocity and densities
	of Daniel Proga's wind models into the appropriate structures and
	allow one to calculate the density at any point in the gridded space.
	

  Arguments:		

  Returns:

  Notes:
	Daniel's models data are the results of calculations with Zeus.  The
	grid is in polar (r,theta) coordinates with theta = pi/2 corresponding
	to the disk, and 0 to the pole.  Recall that python uses cylintrical
	corrdinates with the z=0 being the disk plane.

  History:
	99dec	ksl	Began work
	00jan	ksl	Modified to read larger models and to make it slightly
			more general in terms of the input files.  Also added
			some checking on grid sizes, and set rho to zero outside
			the range of the model, and velocity at the edge of
			the model.
	04jun	ksl	Moved get_proga_wind_params from python.c to this file
	13mar	nsh	Reopened development prior to trip to LV to work with DP.

	

 ************************************************************************/
#include	<stdio.h>
#include	<stdlib.h>
#include	<strings.h>
#include	<string.h>
#include	<math.h>
#include 	"log.h"
#include	"atomic.h"
#include 	"python.h"

#define LINE 132




/***********************************************************
                                       Space Telescope Science Institute

 Synopsis:
	get_proga_wind_params gets input data for Daniel Proga's wind models
Arguments:		

Returns:
 
Description:	
Notes:
History:
 	99dec	ksl	Began work

	13mar	nsh	Reopened work on this file to coincide with 
			NSH working with Proga for a week.
		nsh	First edit - DP says that the second value
			in the two grid files relate to the position
			where the data is defined.
		nsh	Second edit - put in a new variable, proga_proga_thetamax.
			This defines an angle above which one disregards data.
			This was needed because Daniels data includes a flared
			accretion disk, so if you simply read in all the data
			you end up with a very dense blancket over the top of our
			disk, which causes all kinds of problems. At present,
			all that happens is that all angle data above the 
			supplied angle is replaced with the data for that last 
			angle. This should be a very small wedge of data.
	13may	nsh	Now read in the inteernal energy - this allows the
			computation of temperature for a cell, this makes sense
			as a first guess of temperature	
	
**************************************************************/


int
get_proga_wind_params ()
{
  int get_proga ();






  Log
    ("Creating a wind model using a Proga Zeus Calc\n");

  get_proga ();
 // geo.wind_rmin = 8.7e8;	/*Radius where wind begins */
 // if (geo.wind_rmin < geo.rstar)
 //   {
 //     Error
//	("get_stellar_wind_params: It is unreasonable to have the wind start inside the star!\n");
//      Log ("Setting geo.wind_rmin to geo.rstar\n");
 //     geo.wind_rmin = geo.rstar;
 //   }

/* Assign the generic parameters for the wind the generic parameters of the wind */

  geo.wind_rmin = proga_r_edge[0+IGHOST];
   
  geo.wind_rmax = geo.rmax = proga_r_edge[iproga_r]+proga_dr_edge[iproga_r]; //Set the outer edge of the wind to the outer edge of the final defined cell
  Log ("rmax=%e\n",geo.rmax);
  geo.wind_rho_min = 0.0;  //Set wind_rmin 0.0, otherwise wind cones dont work properly 
  Log ("rho_min=%e\n",geo.wind_rho_min);
  geo.wind_rho_max = geo.rmax;  //This is the outer edge of the
  Log ("rho_max=%e\n",geo.wind_rho_max);
  geo.wind_thetamin = proga_theta_edge[0+IGHOST];


Log ("geo.wind_rmin=%e\n",geo.wind_rmin);
Log ("geo.wind_rmax=%e\n",geo.wind_rmax);
Log ("geo.wind_rhomin=%e\n",geo.wind_rho_min);
Log ("geo.wind_rhomax=%e\n",geo.wind_rho_min);


//      geo.wind_rho_min=4*8.7e8;
  //      geo.wind_rho_max=12.*8.7e8;
  //      geo.wind_thetamin=20.0/RADIAN;
  //      geo.wind_thetamax=65./RADIAN;

  /* if modes.adjust_grid is 1 then we have already adjusted the grid manually */
  if (modes.adjust_grid == 0)
    {
      geo.xlog_scale = 0.3 * geo.rstar;
      geo.zlog_scale = 0.3 * geo.rstar;
    }
    
  return (0);
}



int
get_proga ()
{

  FILE *fopen (), *fptr;
  char rfile[LINE], thetafile[LINE], datafile[LINE];
  char aline[LINE];
  int i, j, k;
  double r,r_edge,dr,dr_edge;
  double rho;
  double theta,theta_edge,dtheta,dtheta_edge;
  double vr, vtheta, vphi, energy;
  int irmax,ithetamax,itest;

/*Write something into the file name strings */
  strcpy (rfile, "grid1_r_big.dat");
  strcpy (thetafile, "grid_theta.dat");
  strcpy (datafile, "hdf062.221.dat");


  for (k = 0; k < MAXPROGA; k++)
    {
      proga_r_cent[k] = 0;
      proga_theta_cent[k] = 0;
    }

  proga_ptr = (ProgaPtr) calloc (sizeof (proga_dummy), MAXPROGA * MAXPROGA);

  if (proga_ptr == NULL)
    {
      Error
	("There is a problem in allocating memory for the proga structure\n");
      exit (0);
    }

  rdstr ("Proga_radii_file", rfile);
  if ((fptr = fopen (rfile, "r")) == NULL)
    {
      Error ("Could not open %s\n", rfile);
      exit (0);
    }

  iproga_r = 0;
  
  while (fgets (aline, LINE, fptr) != NULL)
    {
      itest = sscanf (aline, "%d %lf %lf %lf %lf", &i, &r_edge, &r, &dr_edge, &dr);	/*NSH 130322 - minor mod here - it is actually the second value which is the centre of the cell, where the data is defined so we ignore the first */
      if (itest != 5) //We have an line which does not match what we expect, so quit
	{
	Error("proga.c radius file improperly formatted, check for headers and footers\n");
	exit (0);
	}
      proga_r_edge[i] = r_edge;
      proga_r_cent[i] = r;
      proga_dr_edge[i] = dr_edge;
      proga_dr_cent[i] = dr;
//	printf ("PROGA %i, r=%e\n",i,r);
    }
  if (proga_r_edge[IGHOST] < geo.rstar)
 	{
	Error("Major problem, innermost edge of proga radial grid begins inside geo.rstar\n");
	exit (0);
	}

  iproga_r = i;
  Log ("Read %d values from grid_r.dat\n", i);
  fclose (fptr);

  rdstr ("Proga_theta_file", thetafile);
  if ((fptr = fopen (thetafile, "r")) == NULL)
    {
      Error ("Could not open %s\n", thetafile);
      exit (0);
    }

/* NSH 130327 - Added a parameter to define the angle at which the disk is defined
This allows one to disregard theta cells which contain the disk in Daniels model */

 	proga_thetamax=89.9;

  rddoub ("Proga_thetamax(degrees)", &proga_thetamax);

  proga_thetamax=proga_thetamax/RADIAN;
//  printf ("PROGA we will be using data down to theta=%e (%f degrees)\n",proga_thetamax,proga_thetamax*RADIAN);


  j_proga_thetamax = 0;		/* NSH 130605 to remove o3 compile error */
  while (fgets (aline, LINE, fptr) != NULL)
    {
      itest = sscanf (aline, "%d %lf %lf %lf %lf", &i, &theta_edge, &theta, &dtheta_edge, &dtheta);	/*NSH 130322 - minor mod here - it is actually the second value which is the centre of the cell, where the data is defined so we ignore the first */
      if (itest != 5) //We have an line which does not match what we expect, so quit
	{
	Error("proga.c theta file improperly formatted - check for headers and footers\n");
	exit (0);
	}
      proga_theta_cent[i] = theta;
      proga_theta_edge[i] = theta_edge;
      proga_dtheta_edge[i] = dtheta_edge;
      proga_dtheta_cent[i] = dtheta;
      if (proga_theta_edge[i] > proga_thetamax && proga_theta_edge[i - 1] <= proga_thetamax)
	{
	  j_proga_thetamax = i - 1;
	  Log
	    ("current theta  (%f) > theta_max  (%f) so setting j_proga_thetamax=%i\n",
	     theta*RADIAN, proga_thetamax*RADIAN, j_proga_thetamax);
	}
    }
  if (j_proga_thetamax==0 || j_proga_thetamax==i-1)
	{
	Log ("PROGA j_proga_thetamax never bracketed, using all data\n");
	iproga_theta=i;
	geo.wind_thetamax=90. / RADIAN;
	proga_thetamax=90.0/RADIAN;
	}
  else
	{	
 // 	printf ("PROGA j_proga_thetamax=%i, bracketing cells have theta = %f and %f\n",j_proga_thetamax,proga_theta_cent[j_proga_thetamax]*RADIAN,proga_theta_cent[j_proga_thetamax+1]*RADIAN); 
	iproga_theta=j_proga_thetamax;
	geo.wind_thetamax=proga_thetamax;
	}

  Log ("Read %d values from %s\n", iproga_theta, thetafile);
  fclose (fptr);


  rdstr ("Proga_data_file", datafile);
  if ((fptr = fopen (datafile, "r")) == NULL)
    {
      Error ("Could not open %s\n", datafile);

    }


  k = 0;
  irmax=0;
  ithetamax=0; //Zero counters to check all cells actually have data

  while (fgets (aline, LINE, fptr) != NULL)
    {
      	itest=sscanf (aline, "%d %d %lf %lf %lf %lf %lf", &i, &j, &rho, &vr, &vtheta, &vphi, &energy);
      	if (itest==7)
		{
		if (i>irmax) 
			{
			irmax=i;
			}
		if (j>ithetamax) ithetamax=j;
/* NSH 130327 - for the time being, if theta is in the disk, replace with the last
 density above the disk */
	  	proga_ptr[i * MAXPROGA + j].temp = ((2.0 / 3.0) * energy) / ((rho / MPROT) * BOLTZMANN);	//Work out the temperature from the internal energy
          	proga_ptr[i * MAXPROGA + j].rho = rho;
	  	proga_ptr[i * MAXPROGA + j].v[0] = vr;
	  	proga_ptr[i * MAXPROGA + j].v[1] = vtheta;
	  	proga_ptr[i * MAXPROGA + j].v[2] = vphi;
	  	k++;
		}
	else
		{
		Error("Proga data file not correctly formatted - check for headers and footers\n");
		exit(0);
		}

    }
//	Log ("PROGA Maximum r cell with data=%i (iproga_r=%i)\n",irmax,iproga_r);
//	Log ("PROGA Maximum theta cell with data=%i (iproga_theta=%i)\n",ithetamax,iproga_theta);
if (irmax<iproga_r)
	{
//	Log ("PROGA Maximum r cell with data=%i (iproga_r=%i) - resetting\n",irmax,iproga_r);
	iproga_r=irmax;
	}
if (ithetamax<iproga_theta)
	{
//	printf ("PROGA Maximum theta cell with data=%i (iproga_theta=%i)\n",ithetamax,iproga_theta);
	iproga_theta=ithetamax;
	}


  Log ("Read %d values from model %s\n", k, datafile);
  k = 82 * MAXPROGA + 19;


//printf("debug rho %e v %e %e %e\n",
//        proga_ptr[k].rho,
//        proga_ptr[k].v[0],
//        proga_ptr[k].v[1],
//        proga_ptr[k].v[2]);
  fclose (fptr);

 /* We need to reset the grid dimensions to those we have just read in, if we are going to try and match coordinates */
   if (geo.coord_type == RTHETA)
	{
   NDIM = geo.ndim = iproga_r-IGHOST+3; //We need an inner cell to bridge the star and the indisde of the wind, and an outer cell
   MDIM = geo.mdim = iproga_theta-IGHOST+2; //We need one outer cell
	}



  return (0);
}



/***********************************************************
                                       Space Telescope Science Institute

 Synopsis:
	proga_velocity calculates the wind velocity at any position x
	in cartesian coordinates
Arguments:		

Returns:
 
Description:	
Notes:
History:
 	99dec	ksl	Began work
	04dec	ksl	52a -- Mod to prevent divide by zero when 
			r=0.0, and NaN when xxx>1.0;
**************************************************************/


double
proga_velocity (x, v)
     double x[];
     double v[];
{
  double length ();
  int ii, jj;
  int im, jm;
  double f1, f2;
  double r, theta;
  double vr, vtheta, vphi;
  double speed;
  double xxx;
//printf ("Proga_velocity x=%e,%e,%e, v=%e,%e,%e ",x[0],x[1],x[2],v[0],v[1],v[2]);
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

 // printf (" r=%e, theta=%f\n",r,theta*RADIAN);





 // printf ("Proga_theta x %.2g %.2g %.2g  -> r= %.2g theta = %.5g\n", x[0], x[1], x[2], r,theta);
  im = jm = ii = jj = 0;



  while (proga_r_cent[ii] < r && ii <= iproga_r)	//Search through radius array, until array value is greater than your value of r
    ii++;
  while (proga_theta_cent[jj] < theta && jj <= iproga_theta)	//Search through theta array until value is greater than your value of theta
    jj++;



  if (ii > 0 && ii < iproga_r)
    {				//r is in the normal range

      f1 = (r - proga_r_cent[ii - 1]) / (proga_r_cent[ii] - proga_r_cent[ii - 1]);	//Work out fractional position in the ii-1th radial bin where you want to be
      im = ii - 1;		//This is the radial bin below your value of r

    }
  else if (ii == iproga_r)
    {				// r is greater than anything in Proga's model

      f1 = 1;			//If we are outside the model set fractional position to 1
      im = ii - 1;		//And the bin to the outermost
    }
  else
    {
      f1 = 1;			//Otherwise, we must be inide the innermost bin, so again, set fration to 1
      im = 0;			// And the bin to the innermost. Lines below to the same for theta.
    }
//	printf ("f1=%e im=%i ",f1,im);
  if (jj > 0 && jj < iproga_theta)
    {				// theta is inside the normal range

      f2 =
	(theta - proga_theta_cent[jj - 1]) / (proga_theta_cent[jj] -
					 proga_theta_cent[jj - 1]);
      jm = jj - 1;
    }
  else if (jj == iproga_theta)
    {				//theta is more than the maximum theta

      f2 = 1;
      jm = jj - 1;
    }
  else
    {				//theta is less than the expected theta

      f2 = 1;
      jm = 0;
    }
//printf ("f2=%e jm=%i \n",f1,im);
//printf ("Cells to use are %i %i %i %i\n",im * MAXPROGA + jm,im * MAXPROGA + jj,ii*MAXPROGA+jm,ii*MAXPROGA+jj);
//printf ("v_theta= %e %e %e %e\n",proga_ptr[im * MAXPROGA + jm].v[1],proga_ptr[im * MAXPROGA + jj].v[1],proga_ptr[ii
//									    *
//									    MAXPROGA
//									    +
//									    jm].v
//								  [1],proga_ptr[ii
//									    *
//									    MAXPROGA
//									    +
//									    jj].v
//								  [1]);
  vr =
    (1. - f1) * ((1. - f2) * proga_ptr[im * MAXPROGA + jm].v[0] +
		 f2 * proga_ptr[im * MAXPROGA + jj].v[0]) + f1 * ((1. -
								   f2) *
								  proga_ptr[ii
									    *
									    MAXPROGA
									    +
									    jm].v
								  [0] +
								  f2 *
								  proga_ptr[ii
									    *
									    MAXPROGA
									    +
									    jj].v
								  [0]);

  vtheta =
    (1. - f1) * ((1. - f2) * proga_ptr[im * MAXPROGA + jm].v[1] +
		 f2 * proga_ptr[im * MAXPROGA + jj].v[1]) + f1 * ((1. -
								   f2) *
								  proga_ptr[ii
									    *
									    MAXPROGA
									    +
									    jm].v
								  [1] +
								  f2 *
								  proga_ptr[ii
									    *
									    MAXPROGA
									    +
									    jj].v
								  [1]);

  vphi =
    (1. - f1) * ((1. - f2) * proga_ptr[im * MAXPROGA + jm].v[2] +
		 f2 * proga_ptr[im * MAXPROGA + jj].v[2]) + f1 * ((1. -
								   f2) *
								  proga_ptr[ii
									    *
									    MAXPROGA
									    +
									    jm].v
								  [2] +
								  f2 *
								  proga_ptr[ii
									    *
									    MAXPROGA
									    +
									    jj].v
								  [2]);



  v[0] = vr * sin (theta) - vtheta * cos (theta);
  v[1] = vphi;
  v[2] = vr * cos (theta) + vtheta * sin (theta);

  speed = sqrt (v[0] * v[0] + v[1] * v[1] * v[2] * v[2]);

  if (sane_check (speed))
    {
      Error ("proga_velocity:sane_check v %e %e %e\n", v[0], v[1], v[2]);
    }
  return (speed);

}


double
proga_rho (x)
     double x[];
{
  double length ();
  int ii, jj;
  int im, jm;
  double r, theta, rrho;
  double f1, f2;
  r = length (x);
  theta = asin (sqrt (x[0] * x[0] + x[1] * x[1]) / r);
//   printf ("x %.2g %.2g %.2g  -> r= %.2g theta = %.2g ", x[0], x[1], x[2], r,        theta);

  if (r > proga_r_cent[iproga_r])
    {
      Log (" r outside proga grid in proga_rho\n");
      rrho = 1.e-23;
      return (rrho);
    }

  im = jm = ii = jj = 0;
  while (proga_r_cent[ii] < r && ii < iproga_r)
    ii++;
  while (proga_theta_cent[jj] < theta && jj < iproga_theta)
    jj++;

  if (ii > 0)
    {
      f1 = (r - proga_r_cent[ii - 1]) / (proga_r_cent[ii] - proga_r_cent[ii - 1]);
      im = ii - 1;
    }
  else
    f1 = 1;

  if (jj > 0)
    {
      f2 =
	(theta - proga_theta_cent[jj - 1]) / (proga_theta_cent[jj] -
					 proga_theta_cent[jj - 1]);
      jm = jj - 1;
    }
  else
    f2 = 1;


//rrho=proga_ptr[ii*MAXPROGA+jj].rho;

  rrho =
    (1. - f1) * ((1. - f2) * proga_ptr[im * MAXPROGA + jm].rho +
		 f2 * proga_ptr[im * MAXPROGA + jj].rho) + f1 * ((1. -
								  f2) *
								 proga_ptr[ii
									   *
									   MAXPROGA
									   +
									   jm].rho
								 +
								 f2 *
								 proga_ptr[ii
									   *
									   MAXPROGA
									   +
									   jj].rho);

  if (rrho < 1e-23)
    rrho = 1e-23;

 //  printf ("Grid point %d %d rho %e f1=%f f2=%f\n", ii, jj, rrho,f1,f2);

  return (rrho);
}


/***********************************************************
                                       University of Nevada Las Vegas

 Synopsis:
	proga_temp calculates the wind temperature at any position x
	in cartesian coordiantes
Arguments:		

Returns:
 
Description:	
	This code is an exact copy of proga_rho - except it
	maps the temperature, calculated from internal energy,
	from the proga grid onto a cartesian grid
Notes:
History:
 	13apr	nsh	Began work
	
**************************************************************/


double
proga_temp (x)
     double x[];
{
  double length ();
  int ii, jj;
  int im, jm;
  double r, theta, temp;
  double f1, f2;
  r = length (x);
  theta = asin (sqrt (x[0] * x[0] + x[1] * x[1]) / r);
  // printf ("x %.2g %.2g %.2g  -> r= %.2g theta = %.2g\n", x[0], x[1], x[2], r,
//        theta);

  if (r > proga_r_cent[iproga_r])
    {
      Log (" r outside proga grid in proga_temp\n");
      temp = 1e4;
      return (temp);
    }

  im = jm = ii = jj = 0;
  while (proga_r_cent[ii] < r && ii < iproga_r)
    ii++;
  while (proga_theta_cent[jj] < theta && jj < iproga_theta)
    jj++;

  if (ii > 0)
    {
      f1 = (r - proga_r_cent[ii - 1]) / (proga_r_cent[ii] - proga_r_cent[ii - 1]);
      im = ii - 1;
    }
  else
    f1 = 1;

  if (jj > 0)
    {
      f2 =
	(theta - proga_theta_cent[jj - 1]) / (proga_theta_cent[jj] -
					 proga_theta_cent[jj - 1]);
      jm = jj - 1;
    }
  else
    f2 = 1;


//rrho=proga_ptr[ii*MAXPROGA+jj].rho;

  temp =
    (1. - f1) * ((1. - f2) * proga_ptr[im * MAXPROGA + jm].temp +
		 f2 * proga_ptr[im * MAXPROGA + jj].temp) + f1 * ((1. -
								   f2) *
								  proga_ptr[ii
									    *
									    MAXPROGA
									    +
									    jm].temp
								  +
								  f2 *
								  proga_ptr[ii
									    *
									    MAXPROGA
									    +
									    jj].temp);

  if (temp < 1e4)		//Set a lower limit.
    temp = 1e4;



  return (temp);
}



/***********************************************************
                                       Southampton

 Synopsis:
	rtheta_make_zeus_grid defines the cells in a rtheta grid based upon the coordinates that can be read in from a zeus (the hydrocode used by Proga for some simulations)            

Arguments:		
	WindPtr w;	The structure which defines the wind in Python
 
Returns:
 
Description:

	
	This is an attempt to match a zeus grid directly onto an rtheta grid.


History:
	13jun	nsh	76 -- Coded and debugged.


**************************************************************/


int
rtheta_make_zeus_grid (w)
     WindPtr w;
{
  double theta, thetacen;
  int i, j, n;

  



  for (i = 0; i < NDIM; i++)
    {
      for (j = 0; j < MDIM; j++)
	{
 	  wind_ij_to_n (i, j, &n);
		w[n].inwind = W_ALL_INWIND;	
	  if (i == 0)  // The inner edge of the grid should be geo.rstar
		{
//		printf ("PROGA setting inner radial grid cells (n=%i i=%i j=%i) up ",n,i,j);
		w[n].r = geo.rstar;   //So we set the inner edge to be the stellar (or QSO) radius
		w[n].rcen = (geo.rstar + proga_r_edge[IGHOST])/2.0;  //It will be a big cell
		w[n].inwind = W_NOT_INWIND;
//		printf ("inner radius =%e, center=%e\n",w[n].r,w[n].rcen);
		}
	  else if (i+IGHOST-1 > iproga_r) // We are at the radial limit of the data, this last cell will be a ghost cell of our own
		{
//		printf ("PROGA we are outside radial edge of the wind (%i > %i)",i+IGHOST,iproga_r);
		w[n].r = geo.rmax;  // so we set the outer cell to be the edge of the wind
		w[n].rcen = geo.rmax;  // And it has zero volume
		w[n].inwind = W_NOT_INWIND;
//		printf ("edge=%e, center=%e\n",w[n].r,w[n].rcen);
		}
	  else
		{
 	  	w[n].r = proga_r_edge[i+IGHOST-1]; //The -1 is to take account of the fact that i=0 cell is the inner part of the geometry inside the wind
	  	w[n].rcen = proga_r_cent[i+IGHOST-1];
		}
	  if (proga_theta_cent[j+IGHOST]+(proga_dtheta_cent[j+IGHOST]/2.0) > proga_thetamax && proga_theta_cent[j+IGHOST]-(proga_dtheta_cent[j+IGHOST]/2.0) < proga_thetamax ) //This cell bridges the boundary
		{
		theta=proga_theta_cent[j+IGHOST-1]+(proga_theta_cent[j+IGHOST-1]-proga_theta_cent[j+IGHOST-2])/2;
		thetacen=((proga_thetamax+theta)/2.0);
		}
		
	  else if (j+IGHOST >= iproga_theta)  //We are setting up a cell past where there is any data
		{
		thetacen = proga_thetamax; //Set the center and the edge to the maximum extent of the data/interest
		theta = proga_thetamax;
		w[n].inwind = W_NOT_INWIND;
		}
	  else 
		{
		thetacen = proga_theta_cent[j+IGHOST];
		theta = proga_theta_edge[j+IGHOST]; 
		}
	  w[n].theta= theta * RADIAN;
	  w[n].thetacen=thetacen * RADIAN;
	  w[n].x[1] = w[n].xcen[1] = 0.0;
	  w[n].x[0] = w[n].r * sin (theta);
	  w[n].x[2] = w[n].r * cos (theta);
	  w[n].xcen[0] = w[n].rcen * sin (thetacen);
	  w[n].xcen[2] = w[n].rcen * cos (thetacen);
//	  printf ("Cell %i r=%e, rcne=%e, theta=%f, thetacen=%f, x=%e, y=%e, z=%e, inwind=%i\n",n,w[n].r,w[n].rcen,w[n].theta,w[n].thetacen,w[n].x[0],w[n].x[1],w[n].x[2],w[n].inwind);

	}
    }

/*for (i = 0; i < NDIM; i++)
	{
	wind_ij_to_n (i, 0, &n);
	printf ("i=%i, n=%i, r=%e, rcen=%e\n",i,n,w[n].r,w[n].rcen);
	}

for (i = 0; i < MDIM; i++)
	{
	wind_ij_to_n (0, i, &n);
	printf ("j=%i,  iprogatheta=%i, n=%i, theta=%f, thetacen=%f\n",i,iproga_theta,n,w[n].theta,w[n].thetacen);
	}*/


  /* Now set up the wind cones that are needed for calclating ds in a cell */


/*
  cones_rtheta = (ConePtr) calloc (sizeof (cone_dummy), MDIM);
  if (cones_rtheta == NULL)
    {
      Error
	("rtheta_make_grid: There is a problem in allocating memory for the cones structure\n");
      exit (0);

    }


  for (n = 0; n < MDIM; n++)
    {
      cones_rtheta[n].z = 0.0;
      cones_rtheta[n].dzdr = 1. / tan (w[n].theta / RADIAN);	// New definition
    }
*/
  rtheta_make_cones(w); //NSH 130821 broken out into a seperate routine




  /* OK finished successfuly */
  return (0);

}

/***********************************************************
                                       Southampton

 Synopsis:
	rtheta_zeus_volumes replaces rtheta_volumes for a zeus model. We know wether cells are in the wimnd
	so all we need to do is work out the volumes.
Arguments:		
	WindPtr w;	The structure which defines the wind in Python
 
Returns:
 
Description:

	
	This is an attempt to match a zeus grid directly onto an rtheta grid.


History:
	13jun	nsh	76 -- Coded and debugged.


**************************************************************/


int
rtheta_zeus_volumes (w)
     WindPtr w;
{
  int i, j, n;





  double rmin, rmax, thetamin, thetamax;


  for (i = 0; i < NDIM; i++)
    {
      for (j = 0; j < MDIM; j++)
	{
	  wind_ij_to_n (i, j, &n);
	  if (w[n].inwind == W_ALL_INWIND)
	    {

	      rmin = wind_x[i];
	      rmax = wind_x[i + 1];
	      thetamin = wind_z[j] / RADIAN;
	      thetamax = wind_z[j + 1] / RADIAN;

	      //leading factor of 2 added to allow for volume above and below plane (SSMay04)
	      w[n].vol =
		2. * 2. / 3. * PI * (rmax * rmax * rmax -
				     rmin * rmin * rmin) * (cos (thetamin) -
							    cos (thetamax));
	      if (w[n].vol == 0.0)
		{
		Log ("Found wind cell (%i) with no volume (%e) in wind, resetting\n",n,w[n].vol);
		w[n].inwind = W_NOT_INWIND;
		}
	     }
	  else
		w[n].vol=0.0;
}
}
  return (0);
}


