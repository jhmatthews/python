
/***********************************************************
                                       Space Telescope Science Institute

Synopsis:  This is a diagnostic program intended to explore how
	well python mimics a velocity field. It prints out the 
	velocities and velocity gradients in a grid at a specific z height
	above or below the disk.

Arguments:		


Returns:
 
Description:	
		
Notes:

History:
**************************************************************/

//?? Note the behavior of some if not
//  all of the velocity programs have been
//    modified since this program was used.In particularly, kn_velocity returns
//    the velocity in cartesian coordinates, not cylindrical

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "atomic.h"
#include "python.h"
#define LINELENGTH 132
int
main ()
{
  WindPtr w;
  struct photon p;



  char root[LINELENGTH], input[LINELENGTH], wspecfile[LINELENGTH];
  char windsavefile[LINELENGTH];
  char outfileroot[LINELENGTH];
  char atomic[LINELENGTH];
  char parfile[LINELENGTH];
  double v[3], vs[3], vs_xyz[3], vlos, vlos_exact;
  double z;
  int i, j;
  int array_dim;
  double smin;


  FILE *fopen (), *optr;

  double ds;
  double theta, phi;
  double smax;

/* Define the sizes of the wind array */
  ndim = 30;
  mdim = 30;
  NDIM = ndim;
  MDIM = mdim;
  NDIM2 = ndim * mdim;
  dfudge = 1.e5;
  DFUDGE = dfudge;
/* End of definition of wind array sizes */



  printf
    ("This program reads a wind save file created by python and examine the wind structure\n");

  printf ("Root for wind file :");
  fgets (input, LINELENGTH, stdin);
  get_root (root, input);



  printf ("Reading data from file %s\n", root);


  /* Now create the names of all the files which will be written */

  strcpy (wspecfile, root);
  strcpy (windsavefile, root);
  strcpy (parfile, root);

  strcat (wspecfile, ".spec_tot");
  strcat (windsavefile, ".wind_save");
  strcat (parfile, ".pf");

  strcpy (outfileroot, "none");



  if (w == NULL)
    {
      printf
	("There is a problem in allocating memory for the wind structure\n");
      exit (0);
    }

/* Read in the atomic data */
  strcpy (atomic, "atomic/standard");

  opar (parfile);
  rdstr ("Atomic_data", atomic);


  get_atomic_data (atomic);

/* Read in the wind file */


  wind_read (windsavefile);
  w = wmain;

  printf ("Read wind_file %s\n", windsavefile);


/* Should only need to do next step for certain wind types */

  if (geo.wind_type == 1)
    {
//        get_stellar_wind_params ();
    }
  else if (geo.wind_type == 0)
    {
//        get_sv_wind_params ();
    }
  else if (geo.wind_type == 3)
    {
//        get_proga_wind_params ();
    }
  else if (geo.wind_type == 4)
    {
//        get_corona_params ();
    }
  else if (geo.wind_type == 5)
    {
//        get_knigge_wind_params ();
    }
  else if (geo.wind_type == 6)
    {
      get_thierry_params ();
    }


/* Initialize the remaining variables */

  z = 5e9;
  theta = 90;
  phi = 180.;
  smax = 1e10;


  rddoub ("zheight", &z);
  rddoub ("inclination", &theta);
  rddoub ("azimuth", &phi);
  rddoub ("smax", &smax);

  theta /= RADIAN;
  phi /= RADIAN;

  cpar ("py_grid.pf");

  optr = fopen ("py_grid.out", "w");

  p.lmn[0] = cos (phi) * sin (theta);
  p.lmn[1] = sin (phi) * sin (theta);
  p.lmn[2] = cos (theta);
  p.x[2] = z;

  array_dim = 101;
  smin = -smax / 2.;
  ds = smax / (array_dim - 1);
  for (i = 0; i < array_dim; i++)
    {
      p.x[0] = smin + i * ds;
      for (j = 0; j < array_dim; j++)
	{

	  p.x[1] = smin + j * ds;
	  vwind_xyz (w, &p, v);	// vzero is the initial velocity of the photon
	  vlos = dot (p.lmn, v);	// This is vlos as calulated from the grid


	  if (geo.wind_type == 1)
	    {
	      stellar_velocity (&p.x[0], vs);	//sv
	    }
	  else if (geo.wind_type == 0)
	    {
	      sv_velocity (&p.x[0], vs);	//sv
	    }
	  else if (geo.wind_type == 3)
	    {
	      proga_velocity (&p.x[0], vs);	//sv
	    }
	  else if (geo.wind_type == 4)
	    {
	      corona_velocity (&p.x[0], vs);	//sv
	    }
	  else if (geo.wind_type == 5)
	    {
	      kn_velocity (&p.x[0], vs);	//sv
	    }
	  else if (geo.wind_type == 6)
	    {
	      xthierry_velocity (&p.x[0], vs);	//thierry
	    }
// Note that kn_velocity etc return v in cylindrical coords at present
	  if (p.x[0] == 0 && p.x[1] == 0)
	    {
	      vlos_exact = 0.0;
	    }
	  else
	    {
	      project_from_cyl_xyz (p.x, vs, vs_xyz);
	      vlos_exact = dot (p.lmn, vs_xyz);	// This is vlos as calulated from the grid
	    }
	  fprintf (optr, "%d %d %8.3e %8.3e %8.3e %8.3e \n", i, j, p.x[0],
		   p.x[1], vlos, vlos_exact);

	}

    }
  return (0);
}
