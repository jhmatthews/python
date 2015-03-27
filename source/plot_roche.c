/***********************************************************
                                       Space Telescope Science Institute

Synopsis: plot_roche calculates the roche service for a CV using a subset of
the pf file used by python


 
Arguments:		
	plot_roche python

Returns:
 
Description:	
	
		
Notes:
      The file program reads more of python.pf than is necessary.  This
      is in order that this program could be combined with some parts of
      dn_draw and thereby do the entire geometry.

History:
	01jul08 ksl	Began work
 	

**************************************************************/



#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "atomic.h"

#define LINELENGTH 132

#include "python.h"

int
main (argc, argv)
     int argc;
     char *argv[];
{

  int pcycles;
  struct photon ptest;
  int n, nangles;
  double angle[NSPEC - 3], phase[NSPEC - 3];
  int scat_select[NSPEC - 3], top_bot_select[NSPEC - 3];
  char yesno[20];
  int select_extract, select_spectype;
  char root[LINELENGTH], input[LINELENGTH], rochefile[LINELENGTH];
  char dummy[LINELENGTH];
  double xbl;
  double zzz, ds_to_roche_2 ();
  FILE *fptr, *fopen ();
  double theta, phi;

  printf
    ("This program simulates radiative transfer in a (biconical) CV or (spherical) stellar wind\n");

  /* Determine whether input data should be read from a file or from the terminal.  */
  if (argc == 2)
    {
      strcpy (dummy, argv[1]);
    }
  else
    {
      printf ("Input file (interactive=stdin):");
      fgets (dummy, LINELENGTH, stdin);
    }
  get_root (root, dummy);

  if (strncmp (root, "dummy", 5) == 0)
    {
      printf
	("Proceeding to create rdpar file in dummy.pf, but will not run prog\n");
    }
  else if (strncmp (root, "stdin", 5) == 0 || strncmp (root, "rdpar", 5) == 0
	   || root[0] == ' ' || strlen (root) == 0)
    {
      strcpy (root, "mod");
      printf
	("Proceeding in interactive mode\n Output files will have rootname mod\n");
    }
  else
    {
      strcpy (input, root);
      strcat (input, ".pf");

      opar (input);
      printf ("Reading data from file %s\n", input);

    }

  /* Now create the names of all the files which will be written.  Note that some files
     have the same root as the input file, while others have a generic name of python.
     This is intended so that files which you really want to keep have unique names, while
     those which are for short-term diagnostics are overwritten.  ksl 97aug. */

  strcpy (rochefile, root);
  strcat (rochefile, ".roche");


  /* Start logging of errors and comments */

//  Log_init (diagfile);
  fptr = fopen (rochefile, "w");





/* Set initial values for everything in geo struct which basically defines the overall geometry */
  init_geo ();

/* Gather input data */

  /* Describe the basic calculation in terms of the number of iterations which will
     be used to calculate the wind parameters and the number of iterations and wavelength
     range which will be used for the final spectrom.  Also describe the observer's views
     of the system */






  /* Describe the basic binary star system */

  rddoub ("mstar(msol)", &geo.mstar);
  rddoub ("rstar(cm)", &geo.rstar);
  geo.rstar_sq = geo.rstar * geo.rstar;
  if (geo.star_radiation)
    rddoub ("tstar", &geo.tstar);

  rddoub ("msec(msol)", &geo.m_sec);
  rddoub ("period(hr)", &geo.period);

  if (geo.disk_radiation)
    {
      rddoub ("disk.mdot(msol/yr)", &geo.disk_mdot);
      rddoub ("disk.radmax(cm)", &geo.diskrad);
    }
  else
    {
      geo.disk_mdot = 0;
      geo.diskrad = 0;
    }

  /* Describe the boundary layer */
  if (geo.bl_radiation)
    {
      xbl = geo.lum_bl =
	0.5 * G * geo.mstar * geo.disk_mdot * MSOL * MSOL / (geo.rstar * YR);
      rddoub ("lum_bl(ergs/s)", &geo.lum_bl);
      Log ("OK, the bl lum will be about %.2e the disk lum\n",
	   geo.lum_bl / xbl);
      rddoub ("t_bl", &geo.t_bl);
    }
  else
    {
      geo.lum_bl = 0;
      geo.t_bl = 0;
    }

  /* Describe the wind */
  rddoub ("wind.mdot(msol/yr)", &geo.wind_mdot);
  rddoub ("wind.radmax(cm)", &geo.rmax);
  geo.rmax_sq = geo.rmax * geo.rmax;
  rddoub ("wind.t.init", &geo.twind);

  /* Convert all inputs to cgs units */

  geo.mstar *= MSOL;
  geo.m_sec *= MSOL;
  geo.disk_mdot *= MSOL / YR;
  geo.wind_mdot *= MSOL / YR;
  geo.period *= 3600;

  geo.diskrad_sq = geo.diskrad * geo.diskrad;



  /* Calculate additonal parameters associated with the binary star system */

  binary_basics ();



  /* Check that the parameters which have been supplied for the star, disk and boundary layer will
     allow generation of photons where that is appropriate */

  if (geo.tstar <= 0.0)
    geo.star_radiation = 0;
  if (geo.disk_mdot <= 0.0)
    geo.disk_radiation = 0;
  if (geo.t_bl <= 0.0 || geo.lum_bl <= 0.0)
    geo.bl_radiation = 0;


  if (geo.diskrad > 0.0)
    geo.disk_type = 1;		//1 means there is a disk which absorbs

  else
    {
      geo.disk_type = 0;	//0 menas there is not an absorbing disk

      geo.disk_radiation = 0;
    }

  if (geo.star_radiation)
    Log ("There is a star which radiates\n");
  else
    Log ("The star in the system does not radiate\n");

  if (!geo.disk_type)
    Log ("There is no disk in the system \n");
  else if (!geo.disk_radiation)
    Log ("The disk exists, but only absorbs photons\n");
  else
    Log ("There is a disk which radiates and absorbs\n");

  if (geo.bl_radiation)
    Log ("There is a boundary layer which radiates\n");
  else
    Log ("There is no boundary layer\n");



/* Describe the spectra which will be extracted and the way it will be extracted */

  nangles = 4;
  angle[0] = 10;
  angle[1] = 30.;
  angle[2] = 60.;
  angle[3] = 80.;
  for (n = 4; n < NSPEC - 3; n++)
    angle[n] = 45;
  for (n = 0; n < NSPEC - 3; n++)
    {
      phase[n] = 0.5;
      scat_select[n] = 1000;
      top_bot_select[n] = 0;
    }


  select_extract = 1;
  select_spectype = 1;

  if (pcycles > 0)
    {
      if (geo.star_radiation)
	rdint ("Rad_type_for_star(0=bb,1=uniform,2=Kur)_in_final_spectrum",
	       &geo.star_spectype);
      if (geo.disk_radiation)
	rdint ("Rad_type_for_disk(0=bb,1=uniform,2=Kur)_in_final_spectrum",
	       &geo.disk_spectype);
      if (geo.bl_radiation)
	rdint ("Rad_type_for_bl__(0=bb,1=uniform,2=Kur)_in_final_spectrum",
	       &geo.bl_spectype);


      rdint ("no_observers", &nangles);
      if (nangles < 1 || nangles > NSPEC - MSPEC)
	{
	  Error ("no_observers %d should not be > %d or <0\n", nangles,
		 NSPEC - MSPEC);
	  exit (0);
	}

      for (n = 0; n < nangles; n++)
	rddoub ("angle(0=pole)", &angle[n]);
      for (n = 0; n < nangles; n++)
	rddoub ("phase(0=inferior_conjunction)", &phase[n]);

      rdint ("live.or.die(0).or.extract(anything_else)", &select_extract);
      if (select_extract != 0)
	{
	  select_extract = 1;
	  Log ("OK, extracting from specific angles\n");
	}
      else
	Log ("OK, using live or die option\n");

/* Select spectra with certain numbers of scatterings.  See extract 1997 aug 28 ksl */

      strcpy (yesno, "n");
      rdstr ("Select_specific_no_of_scatters_in_spectra(y/n)", yesno);
      if (yesno[0] == 'y')
	{
	  Log
	    ("OK n>MAXSCAT->all; 0<=n<MAXSCAT -> n scatters; n<0 -> >|n| scatters\n");
	  for (n = 0; n < nangles; n++)
	    {
	      rdint ("Select_scatters", &scat_select[n]);
	    }
	}
      strcpy (yesno, "n");
      rdstr ("Select_photons_above_below_the_disk(y/n)", yesno);
      if (yesno[0] == 'y')
	{
	  Log ("OK 0->all; <0 -> below; >0 -> above the disk\n");
	  for (n = 0; n < nangles; n++)
	    {
	      rdint ("Select_above_below", &top_bot_select[n]);
	    }
	}
    }

  /* Select the units of the output spectra.  This is always needed */

  rdint ("spec.type(flambda(1),fnu(2),basic(other)", &select_spectype);
  if (select_spectype == 1)
    {
      Log ("OK, generating flambda at 100pc\n");
    }
  else if (select_spectype == 2)
    {
      Log ("OK, generating fnu at 100 pc\n");
    }
  else
    Log ("OK, basic Monte Carlo spectrum\n");

/* Wrapup and save all the inputs */

  if (strncmp (root, "mod", 3) == 0)
    cpar ("mod.pf");
  else if (strncmp (root, "dummy", 5) == 0)
    {
      cpar ("dummy.pf");
      exit (0);
    }
  else
    cpar ("plot_roche.pf");

/* OK all inputs have been obtained at this point and the inuts have been copied to "mod.pf" or "python.pf" */


  for (theta = 0; theta <= 180; theta += 10)
    {
      for (phi = 0; phi <= 360; phi += 20)
	{
	  ptest.lmn[0] = cos (theta / RADIAN);
	  ptest.lmn[1] = sin (theta / RADIAN) * sin (phi / RADIAN);
	  ptest.lmn[2] = sin (theta / RADIAN) * cos (phi / RADIAN);
	  ptest.x[0] = geo.a;
	  ptest.x[1] = ptest.x[2] = 0;
	  move_phot (&ptest, geo.a);
	  ptest.lmn[0] *= (-1);
	  ptest.lmn[1] *= (-1);
	  ptest.lmn[2] *= (-1);
	  zzz = ds_to_roche_2 (&ptest);
	  move_phot (&ptest, zzz);
	  fprintf (fptr, "%.2le %.2le %.2le\n", ptest.x[0], ptest.x[1],
		   ptest.x[2]);
	}
    }
  return EXIT_SUCCESS;
}




/***********************************************************
                                       Space Telescope Science Institute

 Synopsis:
	init_geo initializes the geo structure to something that is semi reasonable
Arguments:		

Returns:
 
Description:	

		
Notes:


History:
 	98dec	ksl	Coded and debugged.  Much of code was copied from old main routine for
 				python
**************************************************************/

int
init_geo ()
{
  geo.wind_type = 0;		// Schlossman and Vitello

  geo.star_ion_spectype = geo.star_spectype
    = geo.disk_ion_spectype = geo.disk_spectype
    = geo.bl_ion_spectype = geo.bl_spectype = 0;

  geo.log_linear = 0;		/* Set intervals to be logarithmic */

  geo.rmax = 1e11;
  geo.rmax_sq = geo.rmax * geo.rmax;
  geo.rstar = 7e8;
  geo.rstar_sq = geo.rstar * geo.rstar;
  geo.mstar = 0.8;
  geo.m_sec = 0.4;
  geo.period = 3.2;
  geo.tstar = 40000;
  geo.twind = 40000;
  geo.wind_mdot = 1.e-9;

  geo.ioniz_mode = 3;		/* default is on the spot and find the best t */
  geo.line_mode = 3;		/* default is escape probabilites */

  geo.star_radiation = 0;	/* 1 implies star will radiate */
  geo.disk_radiation = 0;	/* 1 implies disk will radiate */
  geo.bl_radiation = 0;		/*1 implies boundary layer will radiate */
  geo.wind_radiation = 0;	/* 1 implies wind will radiate */

  geo.disk_type = 0;		/*1 implies existence of a disk for purposes of absorbtion */
  geo.diskrad = 2.4e10;
  geo.disk_mdot = 1.e-8;

  geo.t_bl = 100000.;

  geo.mstar = 0.8;


  return (0);
}
