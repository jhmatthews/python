#define MAXHYDRO 155
#define IGHOST 0

double hydro_r_cent[MAXHYDRO];
double hydro_r_edge[MAXHYDRO];
double hydro_dr_cent[MAXHYDRO];
double hydro_dr_edge[MAXHYDRO];
double hydro_theta_cent[MAXHYDRO];
double hydro_theta_edge[MAXHYDRO];
double hydro_dtheta_cent[MAXHYDRO];
double hydro_dtheta_edge[MAXHYDRO];
int ihydro_r, ihydro_theta, j_hydro_thetamax, ihydro_mod;
double hydro_thetamax; //The angle at which we want to truncate the theta grid

typedef struct hydro_mod
{
  double rho;
  double v[3];
  double temp;
}
hydro_dummy, *HydroPtr;

HydroPtr hydro_ptr;





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
	04jun	ksl	Moved get_hydro_wind_params from python.c to this file
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
 	99dec	ksl	Began work

	13mar	nsh	Reopened work on this file to coincide with 
			NSH working with Proga for a week.
		nsh	First edit - DP says that the second value
			in the two grid files relate to the position
			where the data is defined.
		nsh	Second edit - put in a new variable, hydro_hydro_thetamax.
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
get_hydro_wind_params ()
{
  int get_hydro ();






  Log
    ("Creating a wind model using a Hydro calculatioon\n");

  get_hydro ();
 // geo.wind_rmin = 8.7e8;	/*Radius where wind begins */
 // if (geo.wind_rmin < geo.rstar)
 //   {
 //     Error
//	("get_stellar_wind_params: It is unreasonable to have the wind start inside the star!\n");
//      Log ("Setting geo.wind_rmin to geo.rstar\n");
 //     geo.wind_rmin = geo.rstar;
 //   }

/* Assign the generic parameters for the wind the generic parameters of the wind */

  geo.wind_rmin = hydro_r_edge[0];
   
  geo.wind_rmax = geo.rmax = hydro_r_edge[ihydro_r]+2.0*(hydro_r_cent[ihydro_r]-hydro_r_edge[ihydro_r]); //Set the outer edge of the wind to the outer edge of the final defined cell
  Log ("rmax=%e\n",geo.rmax);
  geo.wind_rho_min = 0.0;  //Set wind_rmin 0.0, otherwise wind cones dont work properly 
  Log ("rho_min=%e\n",geo.wind_rho_min);
  geo.wind_rho_max = geo.rmax;  //This is the outer edge of the
  Log ("rho_max=%e\n",geo.wind_rho_max);
  geo.wind_thetamin = hydro_theta_edge[0];
  Log ("theta_min=%e\n",geo.wind_thetamin);

Log ("geo.wind_rmin=%e\n",geo.wind_rmin);
Log ("geo.wind_rmax=%e\n",geo.wind_rmax);
Log ("geo.wind_rhomin=%e\n",geo.wind_rho_min);
Log ("geo.wind_rhomax=%e\n",geo.wind_rho_max);




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
get_hydro ()
{

  FILE *fopen (), *fptr;
  char datafile[LINE];
  char aline[LINE];
  char word[LINE];
  int i, j, k;
  double r,r_edge;
  double rho;
  double theta,theta_edge,temp;
  double vr, vtheta, vphi;
  int irmax,ithetamax,itest;

/*Write something into the file name strings */

  strcpy (datafile, "hdf062.dat");


  for (k = 0; k < MAXHYDRO; k++)
    {
      hydro_r_cent[k] = 0;
      hydro_theta_cent[k] = 0;
    }

  hydro_ptr = (HydroPtr) calloc (sizeof (hydro_dummy), MAXHYDRO * MAXHYDRO);

  if (hydro_ptr == NULL)
    {
      Error
	("There is a problem in allocating memory for the hydro structure\n");
      exit (0);
    }

  rdstr ("hydro_file", datafile);
  if ((fptr = fopen (datafile, "r")) == NULL)
    {
      Error ("Could not open %s\n", datafile);
      exit (0);
    }

//	Log ("Hydro datafile = %s\n",datafile);

  	hydro_thetamax=89.9;


  rddoub ("Hydro_thetamax(degrees:negative_means_no_maximum)", &hydro_thetamax);

  //If we have set the maximum angle to a negaive value, we mean that we dont want to restrict.
  if (hydro_thetamax < 0.0) hydro_thetamax=VERY_BIG;
  hydro_thetamax=hydro_thetamax/RADIAN;
  


  ihydro_r = 0;
  j_hydro_thetamax = 0;		/* NSH 130605 to remove o3 compile error */
 irmax=0;
  ithetamax=0; //Zero counters to check all cells actually have data

  while (fgets (aline, LINE, fptr) != NULL)
    {

	 if (aline[0] != '#')
	{
	sscanf (aline, "%s", word);
//	printf ("first word=%s\n",word);

	if (strncmp (word, "ir", 2) == 0)
		{
		Log ("We have a hydro title line, we might do something with this in the future \n");
		}
	else
		{
      		itest = sscanf (aline, "%d %lf %lf %d %lf %lf %lf %lf %lf %lf %lf", &i, &r, &r_edge, &j, &theta, &theta_edge,  &vr, &vtheta, &vphi, &rho, &temp);	
//				printf ("aline=%s itest= %i\n",aline,itest);
      		if (itest != 11) //We have an line which does not match what we expect, so quit
			{
			Error("hydro.c data file improperly formatted\n");
			exit (0);
			}
		// read read the r and theta coordinates into arrays
		hydro_r_edge[i] = r_edge;
      		hydro_r_cent[i] = r;
      		hydro_theta_cent[j] = theta;
      		hydro_theta_edge[j] = theta_edge;
		//keep track of how many r and theta cells there are		
		if (j>ithetamax) ithetamax=j;
		if (i>irmax) irmax=i;
		//If the value of theta in this cell, the edge, is greater than out theta_max, we want to make a note.
    		if (hydro_theta_edge[j] > hydro_thetamax && hydro_theta_edge[j - 1] <= hydro_thetamax)
			{
	  		j_hydro_thetamax = j - 1;
	  		Log
	    		("current theta  (%f) > theta_max  (%f) so setting j_hydro_thetamax=%i\n",
	     			theta*RADIAN, hydro_thetamax*RADIAN, j_hydro_thetamax);
			}
		//If theta is greater than thetamax, then we will replace rho with the last density above the disk
		else if (hydro_theta_edge[j] > hydro_thetamax)
			{
//				printf ("edge is greater than thetamax %e > %e\n",hydro_theta_edge[j] , hydro_thetamax);
			/* NSH 130327 - for the time being, if theta is in the disk, replace with the last
 			density above the disk */
			rho=hydro_ptr[i * MAXHYDRO + j_hydro_thetamax].rho;
			}
		hydro_ptr[i * MAXHYDRO + j].temp = temp;	
          	hydro_ptr[i * MAXHYDRO + j].rho = rho;
//			printf ("input rho %e temp %e vr %e vtheta %e vphi %e\n",rho,temp,vr,vtheta,vphi);
	  	hydro_ptr[i * MAXHYDRO + j].v[0] = vr;
	  	hydro_ptr[i * MAXHYDRO + j].v[1] = vtheta;
	  	hydro_ptr[i * MAXHYDRO + j].v[2] = vphi;
    		}
	}
   }



  if (j_hydro_thetamax==0 || j_hydro_thetamax==i-1)
	{
	Log ("HYDRO j_hydro_thetamax never bracketed, using all data\n");
	ihydro_theta=ithetamax;
	geo.wind_thetamax=90. / RADIAN;
	hydro_thetamax=90.0/RADIAN;
	MDIM = geo.mdim = ihydro_theta+1;
	}
  else
	{	
 	Log ("HYDRO j_hydro_thetamax=%i, bracketing cells have theta = %f and %f\n",j_hydro_thetamax,hydro_theta_cent[j_hydro_thetamax]*RADIAN,hydro_theta_cent[j_hydro_thetamax+1]*RADIAN); 
	ihydro_theta=j_hydro_thetamax;
	geo.wind_thetamax=hydro_thetamax;
	MDIM = geo.mdim = ihydro_theta+2;
	}





  if (hydro_r_edge[0] < geo.rstar)
 	{
	Error("Major problem, innermost edge of hydro radial grid begins inside geo.rstar\n");
	exit (0);
	}

  ihydro_r = irmax;

  Log ("Read %d r values\n", ihydro_r);
  Log ("Read %d theta values\n", ihydro_theta);
  fclose (fptr);

/* Set a couple of last tags*/
 
	geo.coord_type=RTHETA; //At the moment we only deal with RTHETA - in the future we might want to do some clever stuff
	NDIM = geo.ndim = ihydro_r+3; //We need an inner radial cell to bridge the star and the inside of the wind, and an outer cell



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
 	99dec	ksl	Began work
	04dec	ksl	52a -- Mod to prevent divide by zero when 
			r=0.0, and NaN when xxx>1.0;
**************************************************************/


double
hydro_velocity (x, v)
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



  while (hydro_r_cent[ii] < r && ii <= ihydro_r)	//Search through radius array, until array value is greater than your value of r
    ii++;
  while (hydro_theta_cent[jj] < theta && jj <= ihydro_theta)	//Search through theta array until value is greater than your value of theta
    jj++;



  if (ii > 0 && ii < ihydro_r)
    {				//r is in the normal range

      f1 = (r - hydro_r_cent[ii - 1]) / (hydro_r_cent[ii] - hydro_r_cent[ii - 1]);	//Work out fractional position in the ii-1th radial bin where you want to be
      im = ii - 1;		//This is the radial bin below your value of r

    }
  else if (ii == ihydro_r)
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
  if (jj > 0 && jj < ihydro_theta)
    {				// theta is inside the normal range

      f2 =
	(theta - hydro_theta_cent[jj - 1]) / (hydro_theta_cent[jj] -
					 hydro_theta_cent[jj - 1]);
      jm = jj - 1;
    }
  else if (jj == ihydro_theta)
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
//printf ("Cells to use are %i %i %i %i\n",im * MAXHYDRO + jm,im * MAXHYDRO + jj,ii*MAXHYDRO+jm,ii*MAXHYDRO+jj);
//printf ("v_theta= %e %e %e %e\n",hydro_ptr[im * MAXHYDRO + jm].v[1],hydro_ptr[im * MAXHYDRO + jj].v[1],hydro_ptr[ii
//									    *
//									    MAXHYDRO
//									    +
//									    jm].v
//								  [1],hydro_ptr[ii
//									    *
//									    MAXHYDRO
//									    +
//									    jj].v
//								  [1]);
  vr =
    (1. - f1) * ((1. - f2) * hydro_ptr[im * MAXHYDRO + jm].v[0] +
		 f2 * hydro_ptr[im * MAXHYDRO + jj].v[0]) + f1 * ((1. -
								   f2) *
								  hydro_ptr[ii
									    *
									    MAXHYDRO
									    +
									    jm].v
								  [0] +
								  f2 *
								  hydro_ptr[ii
									    *
									    MAXHYDRO
									    +
									    jj].v
								  [0]);

  vtheta =
    (1. - f1) * ((1. - f2) * hydro_ptr[im * MAXHYDRO + jm].v[1] +
		 f2 * hydro_ptr[im * MAXHYDRO + jj].v[1]) + f1 * ((1. -
								   f2) *
								  hydro_ptr[ii
									    *
									    MAXHYDRO
									    +
									    jm].v
								  [1] +
								  f2 *
								  hydro_ptr[ii
									    *
									    MAXHYDRO
									    +
									    jj].v
								  [1]);

  vphi =
    (1. - f1) * ((1. - f2) * hydro_ptr[im * MAXHYDRO + jm].v[2] +
		 f2 * hydro_ptr[im * MAXHYDRO + jj].v[2]) + f1 * ((1. -
								   f2) *
								  hydro_ptr[ii
									    *
									    MAXHYDRO
									    +
									    jm].v
								  [2] +
								  f2 *
								  hydro_ptr[ii
									    *
									    MAXHYDRO
									    +
									    jj].v
								  [2]);



  v[0] = vr * sin (theta) - vtheta * cos (theta);
  v[1] = vphi;
  v[2] = vr * cos (theta) + vtheta * sin (theta);

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
  double r, theta, rrho;
  double f1, f2;
  r = length (x);
  theta = asin (sqrt (x[0] * x[0] + x[1] * x[1]) / r);
//     printf ("x %.2g %.2g %.2g  -> r= %.2g theta = %.2g ", x[0], x[1], x[2], r,        theta);

  if (r > hydro_r_cent[ihydro_r])
    {
      Log (" r outside hydro grid in hydro_rho\n");
      rrho = 1.e-23;
      return (rrho);
    }

  im = jm = ii = jj = 0;
  while (hydro_r_cent[ii] < r && ii < ihydro_r)
    ii++;
  while (hydro_theta_cent[jj] < theta && jj < ihydro_theta)
    jj++;

  if (ii > 0)
    {
      f1 = (r - hydro_r_cent[ii - 1]) / (hydro_r_cent[ii] - hydro_r_cent[ii - 1]);
      im = ii - 1;
    }
  else
    f1 = 1;

  if (jj > 0)
    {
      f2 =
	(theta - hydro_theta_cent[jj - 1]) / (hydro_theta_cent[jj] -
					 hydro_theta_cent[jj - 1]);
      jm = jj - 1;
    }
  else
    f2 = 1;


//rrho=hydro_ptr[ii*MAXHYDRO+jj].rho;
//  printf ("Data rho=%e\n",hydro_ptr[im * MAXHYDRO + jm].rho);
  rrho =
    (1. - f1) * ((1. - f2) * hydro_ptr[im * MAXHYDRO + jm].rho +
		 f2 * hydro_ptr[im * MAXHYDRO + jj].rho) + f1 * ((1. -
								  f2) *
								 hydro_ptr[ii
									   *
									   MAXHYDRO
									   +
									   jm].rho
								 +
								 f2 *
								 hydro_ptr[ii
									   *
									   MAXHYDRO
									   +
									   jj].rho);
//  printf ("Rho= %e\n",rrho);

  if (rrho < 1e-23)
    rrho = 1e-23;

 //  printf ("Grid point %d %d rho %e f1=%f f2=%f\n", ii, jj, rrho,f1,f2);

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
 	13apr	nsh	Began work
	
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
  // printf ("x %.2g %.2g %.2g  -> r= %.2g theta = %.2g\n", x[0], x[1], x[2], r,
//        theta);

  if (r > hydro_r_cent[ihydro_r])
    {
      Log (" r outside hydro grid in hydro_temp\n");
      temp = 1e4;
      return (temp);
    }

  im = jm = ii = jj = 0;
  while (hydro_r_cent[ii] < r && ii < ihydro_r)
    ii++;
  while (hydro_theta_cent[jj] < theta && jj < ihydro_theta)
    jj++;

  if (ii > 0)
    {
      f1 = (r - hydro_r_cent[ii - 1]) / (hydro_r_cent[ii] - hydro_r_cent[ii - 1]);
      im = ii - 1;
    }
  else
    f1 = 1;

  if (jj > 0)
    {
      f2 =
	(theta - hydro_theta_cent[jj - 1]) / (hydro_theta_cent[jj] -
					 hydro_theta_cent[jj - 1]);
      jm = jj - 1;
    }
  else
    f2 = 1;


//rrho=hydro_ptr[ii*MAXHYDRO+jj].rho;

  temp =
    (1. - f1) * ((1. - f2) * hydro_ptr[im * MAXHYDRO + jm].temp +
		 f2 * hydro_ptr[im * MAXHYDRO + jj].temp) + f1 * ((1. -
								   f2) *
								  hydro_ptr[ii
									    *
									    MAXHYDRO
									    +
									    jm].temp
								  +
								  f2 *
								  hydro_ptr[ii
									    *
									    MAXHYDRO
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
rtheta_make_hydro_grid (w)
     WindPtr w;
{
  double theta, thetacen, dtheta;
  int i, j, n;

  


  for (i = 0; i < NDIM; i++)
    {
      for (j = 0; j < MDIM; j++)
	{
 	  wind_ij_to_n (i, j, &n);
		w[n].inwind = W_ALL_INWIND;	
	  if (i == 0)  // The inner edge of the grid should be geo.rstar
		{
//		printf ("HYDRO setting inner radial grid cells (n=%i i=%i j=%i) up \n",n,i,j);
		w[n].r = geo.rstar;   //So we set the inner edge to be the stellar (or QSO) radius
		w[n].rcen = (geo.rstar + hydro_r_edge[0])/2.0;  //It will be a big cell
		w[n].inwind = W_NOT_INWIND;
//		printf ("inner radius =%e, center=%e\n",w[n].r,w[n].rcen);
		}
	  else if (i-1 > ihydro_r) // We are at the radial limit of the data, this last cell will be a ghost cell of our own
		{
//		printf ("HYDRO we are outside radial edge of the wind (%i > %i)\n",i-1,ihydro_r);
		w[n].r = geo.rmax;  // so we set the outer cell to be the edge of the wind
		w[n].rcen = geo.rmax;  // And it has zero volume
		w[n].inwind = W_NOT_INWIND;
//		printf ("edge=%e, center=%e\n",w[n].r,w[n].rcen);
		}
	  else
		{
 	  	w[n].r = hydro_r_edge[i-1]; //The -1 is to take account of the fact that i=0 cell is the inner part of the geometry inside the wind
	  	w[n].rcen = hydro_r_cent[i-1];
		}
		dtheta=2.0*(hydro_theta_cent[j]-hydro_theta_edge[j]);
	


//		printf ("dtheta=%e test = %e, test2= %e, hydro_thetamax = %e j = %i ihydro_theta=%i\n",dtheta,hydro_theta_cent[j]+(dtheta/2.0),hydro_theta_cent[j]-(dtheta/2.0), hydro_thetamax,j,ihydro_theta);
	  if (hydro_theta_cent[j]+(dtheta/2.0) > hydro_thetamax && hydro_theta_cent[j]-(dtheta/2.0) < hydro_thetamax ) //This cell bridges the boundary - we reset it so the lower edge is at the disk boundary
		{
//		printf ("We have reached the disk %e > %e\n",hydro_theta_cent[j]+(hydro_dtheta_cent[j]/2.0), hydro_thetamax);
		theta=hydro_theta_edge[j];
		thetacen=((hydro_thetamax+theta)/2.0);
//		printf ("Setting theta=%e thetacen=%e\n",theta,thetacen);
		}
		
	  else if (j >= ihydro_theta)  //We are setting up a cell past where there is any data
		{
//		printf ("we are past ihydro_theta , %i, %i \n",j,ihydro_theta);
		thetacen = hydro_thetamax; //Set the center and the edge to the maximum extent of the data/interest
		theta = hydro_thetamax;
		w[n].inwind = W_NOT_INWIND;
		}
	  else 
		{
		thetacen = hydro_theta_cent[j];
		theta = hydro_theta_edge[j]; 
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
/*
for (i = 0; i < NDIM; i++)
	{
	wind_ij_to_n (i, 0, &n);
	printf ("i=%i, n=%i, r=%e, rcen=%e\n",i,n,w[n].r,w[n].rcen);
	}

for (i = 0; i < MDIM; i++)
	{
	wind_ij_to_n (0, i, &n);
	printf ("j=%i,  ihydrotheta=%i, n=%i, theta=%f, thetacen=%f\n",i,ihydro_theta,n,w[n].theta,w[n].thetacen);
	}
*/

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
rtheta_hydro_volumes (w)
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


