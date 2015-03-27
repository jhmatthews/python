
/***********************************************************
                                       Space Telescope Science Institute

Synopsis:  This is a diagnostic program intended to explore how
	well python mimics a velocity field.

Arguments:		
	In interactive mode, one is asked for the number of photons
	one wants to generate and then repetively for the posiion
	of the photon, in cylindrical coordianates, and then the angle 
	theta, phi in which the photon is to travel.  A conversion is
	then made to the position and direction cosines.

           rdint ("nphottot", &nphottot);

	  rddoub ("rho_start(wd)", &r);
	  rddoub ("phi_start(deg", &phi_start);
	  rddoub ("z_start(wd)", &z_start);
	  rddoub ("theta_phot(deg)", &theta);
	  rddoub ("phi_phot(deg)", &phi);
	
	When reading from a file one enters the position and direction
	cosines for the photon, in cgs coordinates, one photon to a line.

	 &x[0], &x[1], &x[2], &lmn[0], &lmn[1], &lmn[2]);



Returns:
 
Description:	
	py_ray reads and existing wind_save file and then allows one
	to define rays which are then traced out of the system.  One
	can define the rays from the command line, or from a file.
		
Notes:
	This program is not very well documented, but it probably
	originated with python itself.


History:
	04nov	ksl	Updates so compiled with python_53b
**************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "atomic.h"
#include "python.h"
int nc4 = 0;


#define LINELENGTH 132

int
main ()
{
  WindPtr w;
  PhotPtr p;


  int nphot, nphottot;

  char root[LINELENGTH], input[LINELENGTH], wspecfile[LINELENGTH],
    specfile[LINELENGTH];
  char windradfile[LINELENGTH], windsavefile[LINELENGTH];
  char outfileroot[LINELENGTH];
  char photfile[LINELENGTH], atomic[LINELENGTH];
  char parfile[LINELENGTH];
  char dummy[LINELENGTH];
  double freq;


  FILE *fopen (), *optr, *iptr;

  double ds;
  double z_start, phi_start;
  double x[3], lmn[3], r, theta, phi, weight;
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
  strcpy (specfile, root);
  strcpy (windradfile, root);
  strcpy (windsavefile, root);
  strcpy (photfile, root);
  strcpy (parfile, root);

  strcat (wspecfile, ".spec_tot");
  strcat (specfile, ".spec");
  strcat (windradfile, ".wind_rad");
  strcat (windsavefile, ".wind_save");
  strcat (photfile, ".phot");
  strcat (parfile, ".pf");

  strcpy (outfileroot, "none");

//  w = (WindPtr) calloc (sizeof (wind_dummy), NDIM2);

  NPHOT = 10;
  p = (PhotPtr) calloc (sizeof (p_dummy), NPHOT);

  if (w == NULL)
    {
      printf
	("There is a problem in allocating memory for the wind structure\n");
      exit (0);
    }

/* Read in the atomic data */
  strcpy (atomic, "atomic/standard39");

  opar (parfile);
  rdstr ("Atomic_data", atomic);


  get_atomic_data (atomic);

  while ((ion[nc4].z != 6 || ion[nc4].istate != 4) && nc4 < nions)
    nc4++;
  if (nc4 == nions)
    {
      Error ("Could not locate c 4!!\n");
      exit (0);
    }

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


/* Define some photons */

  rdstr ("photfile(none.if.grid", photfile);
  if (strncmp (photfile, "none", 4) == 0)
    {				// Gen from data

      nphottot = 1;
      rdint ("nphottot", &nphottot);

      r = 5.;
      phi_start = 135;
      z_start = 0.0;
      theta = 45.;
      phi = 305;


      for (nphot = 0; nphot < nphottot; nphot++)
	{

	  rddoub ("rho_start(wd)", &r);
	  rddoub ("phi_start(deg", &phi_start);
	  rddoub ("z_start(wd)", &z_start);
	  rddoub ("theta_phot(deg)", &theta);
	  rddoub ("phi_phot(deg)", &phi);

	  phi_start /= RADIAN;
	  theta /= RADIAN;
	  phi /= RADIAN;

	  r *= geo.rstar;

	  x[0] = r * cos (phi_start);
	  x[1] = r * sin (phi_start);
	  x[2] = z_start;
	  freq = C / (1550. * 1.e-8);
	  weight = 1;
	  gen_one_phot (&p[nphot], weight, freq, x, theta, phi);
	}
    }

  else
    {
      if ((iptr = fopen (photfile, "r")) == NULL)
	{
	  Error ("py_ray: Could not open %s to read desired rays\n",
		 photfile);
	  exit (0);
	}

      nphot = 0;
      while (nphot < NPHOT && fgets (dummy, LINELENGTH, iptr) != NULL)
	{
	  sscanf (dummy, "%lf %lf %lf %lf %lf %lf", &x[0], &x[1], &x[2],
		  &lmn[0], &lmn[1], &lmn[2]);
	  stuff_v (x, p[nphot].x);
	  stuff_v (lmn, p[nphot].lmn);
	  nphot++;
	  nphottot = nphot;
	}



    }

  cpar ("py_ray.pf");

  optr = fopen ("py_ray.out", "w");

  smax = 1e11;
  ds = smax / 1000;
//  ds = smax / 100;

  for (nphot = 0; nphot < nphottot; nphot++)
    {
      process_one_phot (optr, w, &p[nphot], ds, smax);
    }
  return (0);
}


int
gen_one_phot (p, weight, freq, x, theta, phi)
     PhotPtr p;
     double weight, freq;
     double x[3];
     double theta, phi;		// relative to a cartesian system
{
  double lmn[3];

  lmn[0] = sin (theta) * cos (phi);
  lmn[1] = sin (theta) * sin (phi);
  lmn[2] = cos (theta);

  stuff_v (x, p->x);
  stuff_v (lmn, p->lmn);

  p->freq = freq;
  p->w = 1.;
  p->tau = 0.0;
  p->istat = 0;
  p->nscat = 0;
  p->nrscat = 0;
  p->grid = 0;
  p->origin = PTYPE_DISK;

  return (0);
}

int
process_one_phot (optr, w, p, ds, smax)
     FILE *optr;
     WindPtr w;			// These is the whole wind structure
     PhotPtr p;
     double ds, smax;
{
  struct photon tphot, zphot;
  double s, v[3], vs[3], zvs[3], lmn[3];
  double vels, zvels;
  double dvwind_ds ();
  double dvds, dvds_exact;
  double deltas;
  double r;
  double dot ();
  double ctheta, stheta, vf[3], vzero[3], dv[3];
  int vwind_xyz ();
  double thierry_velocity (), xthierry_velocity ();	//thierry
  double dlambda, vlos, vlos_exact;
  double vlosz, dvds_vwind, vz[3];
  double density, get_ion_density (), massden;

  s = 0;
  stuff_phot (p, &tphot);
  vwind_xyz (w, &tphot, vzero);	// vzero is the initial velocity of the photon

  while (s < smax)
    {
      move_phot (&tphot, ds);
      v[0] = v[1] = v[2] = 0;
      vwind_xyz (w, &tphot, v);	// v is the current velocity of the photon
      vsub (v, vzero, dv);
/* A photon whose which was the product of of resonace transition, would be emitted
at a wavelength of 
      lambda = 1550 * (1 - dot (tphot.lmn,vzero))
This self same photon would have a lambda in a local rest frame along the line of sight
of    lambda = lambda_zero * (1 + dot (tphot.lmn,v))

Therefore this is the wavelength in the local rest frame that a emitted at 1550 Angstroms 
would have in local rest frame later on */

      dlambda = 1550. * (1 + dot (tphot.lmn, dv) / C);
      vlos = dot (tphot.lmn, v);	// This is vlos as calculated from the grid




      if (geo.wind_type == 1)
	{
	  stellar_velocity (&tphot.x[0], vs);	//sv
	}
      else if (geo.wind_type == 0)
	{
	  sv_velocity (&tphot.x[0], vs);	//sv
	}
      else if (geo.wind_type == 3)
	{
	  proga_velocity (&tphot.x[0], vs);	//sv
	}
      else if (geo.wind_type == 4)
	{
	  corona_velocity (&tphot.x[0], vs);	//sv
	}
      else if (geo.wind_type == 5)
	{
	  kn_velocity (&tphot.x[0], vs);	//sv
	}
      else if (geo.wind_type == 6)
	{
	  xthierry_velocity (&tphot.x[0], vs);	//thierry
//        thierry_velocity (&tphot.x[0], vs);   //thierry
	}

/*We now have vs (exact) but it is in cylindrical coordinates.  
We need to convert it to xyz coordinates.  vf is exact in
xyz coordinates in the sense that it is derived from the
original specification of the velocities  */

      ctheta =
	tphot.x[0] / sqrt (tphot.x[0] * tphot.x[0] + tphot.x[1] * tphot.x[1]);
      stheta =
	tphot.x[1] / sqrt (tphot.x[0] * tphot.x[0] + tphot.x[1] * tphot.x[1]);
      vf[0] = ctheta * vs[0] - stheta * vs[1];
      vf[1] = stheta * vs[0] + ctheta * vs[1];
      vf[2] = vs[2];
      if (tphot.x[2] < 0)
	vf[2] *= -1;

      vlos_exact = dot (tphot.lmn, vf);



      dvds = dvwind_ds (w, &tphot);	/* This is dvwind_ds as calculated
					   from the gradients in the wind array */


// We know have what we need but lmn is in xyz coords, while vs is cylindrical
      project_from_xyz_cyl (tphot.x, tphot.lmn, lmn);

      vels = dot (lmn, vs);	/*So vels is the exact speeed of the wind along
				   the line of sight of the photon */

      r = length (&tphot.x[0]);
      deltas = 1e-3 * r;
      stuff_phot (&tphot, &zphot);
      move_phot (&zphot, deltas);

      vwind_xyz (w, &zphot, vz);	// v is the current velocity of the photon
      vlosz = dot (tphot.lmn, vz);	// This is vlos as calulated from the grid
      density = get_ion_density (zphot.x, nc4);
      massden = kn_rho (zphot.x);


      if (geo.wind_type == 1)
	{
	  stellar_velocity (&zphot.x[0], zvs);	//sv
	}
      else if (geo.wind_type == 0)
	{
	  sv_velocity (&zphot.x[0], zvs);	//sv
	}
      else if (geo.wind_type == 3)
	{
	  proga_velocity (&zphot.x[0], zvs);	//sv
	}
      else if (geo.wind_type == 4)
	{
	  corona_velocity (&zphot.x[0], zvs);	//sv
	}
      else if (geo.wind_type == 5)
	{
	  kn_velocity (&zphot.x[0], zvs);	//sv
	}
      else if (geo.wind_type == 6)
	{
	  xthierry_velocity (&zphot.x[0], zvs);	//thierry
//        thierry_velocity (&zphot.x, zvs);     //thierry
	}



      project_from_xyz_cyl (tphot.x, tphot.lmn, lmn);
      zvels = dot (lmn, zvs);

//zvels=lmn[0]*zvs[0]+lmn[1]*zvs[1]+lmn[2]*zvs[2];

      dvds_exact = (zvels - vels) / deltas;	/* This is exact in the
						   sense that we have gone all the way back 
						   to the original prescriptions
						   for velocity to get dvds */

      dvds_vwind = (vlos - vlosz) / deltas;	/* This is the dv_ds that would be
						   calulated from the v array of the grid */

      printf ("zvels %f vels %f deltas %f \n", zvels, vels, deltas);
      if (where_in_wind (tphot.x) == 0)
	printf
	  ("%.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.2f %.3e %.3e\n",
	   tphot.x[0], tphot.x[1], tphot.x[2], s, v[0], v[1], v[2], vs[0],
	   vs[1], vs[2], dvds, dvds_exact, dlambda, vlos, vlos_exact);

      if (where_in_wind (tphot.x) == 0)
	{
	  fprintf (optr,
		   "%.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.2f %.3e %.3e %.3e %.3e\n",
		   tphot.x[0], tphot.x[1], tphot.x[2], s, v[0], v[1], v[2],
		   vf[0], vf[1], vf[2], fabs (dvds), fabs (dvds_exact),
		   dlambda, vlos, vlos_exact, density, massden);

	  fprintf (optr,
		   "%.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.2f %.3e %.3e %.3e %.3e\n",
		   tphot.x[0], tphot.x[1], tphot.x[2], s, v[0], v[1], v[2],
		   vf[0], vf[1], vf[2], fabs (dvds_vwind), fabs (dvds_exact),
		   dlambda, vlos, vlos_exact, density, massden);
	}
      s += ds;
    }
  return (0);
}
