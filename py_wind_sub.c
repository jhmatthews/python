
/***********************************************************
                                       Space Telescope Science Institute

 Synopsis:
	py_wind is a program which can be used to display various parameters of a wind
	as calculated by python.  This file contains many of the subroutiens
	used

Arguments:


Returns:

Description:
	
		
Notes:

History:
 	97jun	ksl	Coding on py_wind began.
	05jul	ksl	56d -- Moved all of the subroutines to py_wind_sub so that it was
			easier to make templates.

**************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "atomic.h"
#include "python.h"



/*
 * This routine controls what is displayed on the screen.  There are
 * basically two options, determined by the variable "determine"
 * 
 * If determine==1, then zoom sets up the variables py_wind_min, py_wind_max and
 * py_wind_delta so that the the entire wind is displayed, but it is
 * subsampled in the x direction.
 * 
 * If determine!=1, then a section of the wind is displayed at full resolution
 * 
 * History: 07jul	ksl	Made modifications to allow for the possibility that
 * the wind has fewer than 10 elemetns in the x direction
 */
int
zoom (direction)
     int direction;
{
  int center;
  if (direction == 1)
    {				/* then unzoom */
      Log ("Showing selected positions throughout wind\n");
      py_wind_min = 0;
      py_wind_max = NDIM;
      py_wind_delta = NDIM / 10;
      /*
       * Allow for the possibility that the wind has an xdimension
       * < 10
       */
      if (py_wind_delta < 1)
	py_wind_delta = 1;
    }
  else
    {
      Log ("Select part of wind to display\n");
      center = 5;
      rdint ("Center_x", &center);
      py_wind_min = center - 5;
      if (py_wind_min < 0)
	{
	  Log
	    ("zoom: this choice of center is needlessly close to origin, adjusting py_wind_min to 0\n");
	  py_wind_min = 0;
	}
      py_wind_max = py_wind_min + 10;
      if (py_wind_max > NDIM)
	{
	  py_wind_max = NDIM;
	  Log
	    ("zoom: this choice of py_wind_max is lager than NDIM, adusting py_wind_max to NDIM");
	}
      py_wind_delta = 1;
    }
  return (0);
}


/* overview */
int
overview (w, rootname)
     WindPtr w;
     char rootname[];
{
  //double lum, wind_luminosity (); JM130621: shouldn't really call this here
  int n;
  double heating, lines, ff, photo;

  heating = lines = ff = photo = 0.0;

  for (n = 0; n < NPLASMA; n++)
    {
      heating += plasmamain[n].heat_tot;
      lines += plasmamain[n].heat_lines;
      photo += plasmamain[n].heat_photo;
      ff += plasmamain[n].heat_ff;
    }
  /* lum = wind_luminosity (0., 1.e20); JM130621: shouldn't really call this here. windsave bug fix means
     we should trust what is in the geo structure */
  Log (" Total emission %8.2e heating %8.2e\n", geo.lum_ioniz, heating);
  Log ("    ff emission %8.2e heating %8.2e\n", geo.lum_ff_ioniz, ff);
  Log ("    fb emission %8.2e heating %8.2e\n", geo.lum_fb_ioniz, photo);
  Log ("  line emission %8.2e heating %8.2e\n", geo.lum_lines_ioniz, lines);
  return (0);
}

/* Summary of everything at a given position */

int
position_summary (w)
     WindPtr w;
{
  double x[3], v[3];
  struct photon p;
  int n;
  int nplasma;

  x[0] = geo.wind_rmax / 5;
  x[1] = 0.0;
  x[2] = x[0];

a:Log ("Input x=0,y=0,z=0 to return to main routine\n");
  rddoub ("x", &x[0]);
  rddoub ("y", &x[1]);
  rddoub ("z", &x[2]);

  if (length (x) == 0.0)
    return (0);

  n = where_in_grid (x);
  nplasma = wmain[n].nplasma;

  Log ("Position %8.2e  %8.2e %8.2e  Cell %5d\n", x[0], x[1], x[2], n);

  Log ("Vertex position %8.2e  %8.2e %8.2e  rtheta %8.2e %8.2e \n", w[n].x[0],
       w[n].x[1], w[n].x[2], w[n].r, w[n].theta);
  Log ("Center position %8.2e  %8.2e %8.2e  rtheta %8.2e %8.2e \n",
       w[n].xcen[0], w[n].xcen[1], w[n].xcen[2], w[n].rcen, w[n].thetacen);

  Log ("Electron density: %8.2e  rho %8.2e\n", plasmamain[nplasma].ne,
       plasmamain[nplasma].rho);
  Log ("Vel cell: %8.2e  %8.2e %8.2e\n", w[n].v[0], w[n].v[1], w[n].v[2]);

  p.x[0] = x[0];
  p.x[1] = x[1];
  p.x[2] = x[2];

  vwind_xyz (&p, v);
  Log ("Velocity: %8.2e  %8.2e %8.2e\n", v[0], v[1], v[2]);


  goto a;

  return (0);
}


/* A summary of the energy absorbed in a cell 

10nov	ksl	What was here previously was bizarre as the filename was
		being added to all the time and this caused ultimately
		a segmentation problem

*/

int
abs_summary (w, rootname, ochoice)
     WindPtr w;
     char rootname[];
     int ochoice;
{
  int n;
  double x, xtot;
  char c;
  char name[LINELENGTH];
  char filename[LINELENGTH];
  int nplasma;

  c = 't';
  xtot = 0.0;

  printf
    ("Absorption tot=t, lines=l,f=ff,b=fb,h=hydrogen,i=he1,j=he2,z=heavy elements\n");
  rdchar ("Choice", &c);


  strcpy (filename, rootname);

  switch (c)
    {
    case 't':			/* Total  Absorption */
      strcpy (name, "Total Absorbtion");
      strcat (filename, "heat_tot");
      break;
    case 'f':			/* ff */
      strcpy (name, "Free Free Absorbtion");
      strcat (filename, "heat_ff");
      break;
    case 'b':			/* Photoionization */
      strcpy (name, "Total photoionization Heating");
      strcat (filename, "heat_photo");
      break;
    case 'l':			/* Line  heating */
      strcpy (name, "Resonance Line Heating");
      strcat (filename, "heat_lines");
      break;
    case 'h':			/* H photoionization */
      strcpy (name, "H  Photoionization Heating");
      strcat (filename, "heat_h");
      break;
    case 'i':			/* He1 photo */
      strcpy (name, "He I Photoionization Heating");
      strcat (filename, "heat_he1");
      break;
    case 'j':			/* He2 photo */
      strcpy (name, "He 2 Photoionization Heating");
      strcat (filename, "heat_he2");
      break;
    case 'z':			/* Metal photo */
      strcpy (name, "Metal Photoionization Heating");
      strcat (filename, "heat_z");
      break;
    default:
      printf ("Not a valid choice\n");
      return (0);
    }



  for (n = 0; n < NDIM2; n++)
    {
      x = 0.0;
      if (w[n].vol > 0.0)
	{
	  nplasma = w[n].nplasma;
	  switch (c)
	    {
	    case 't':
	      {			/* Total heating */
		x = plasmamain[nplasma].heat_tot;
		break;
	    case 'f':		/* ff heating */
		x = plasmamain[nplasma].heat_ff;
		break;
	    case 'b':		/* photoionization heating */
		x = plasmamain[nplasma].heat_photo;
		break;
	    case 'l':		/* Line heating */
		x = plasmamain[nplasma].heat_lines;
		break;
	    case 'h':		/* H heating */
		x = plasmamain[nplasma].heat_ion[0];
		break;
	    case 'i':		/* He1 heating */
		x = plasmamain[nplasma].heat_ion[2];
		break;
	    case 'j':		/* He2 heating */
		x = plasmamain[nplasma].heat_ion[3];
		break;
	    case 'z':		/* Line heating of high z elements */
		x = plasmamain[nplasma].heat_z;
		break;
	    default:
		printf ("Not a valid choice\n");
		return (0);
	      }
	      xtot += x;

	    }
	  aaa[n] = x;
	}
    }

  display (name);
  /* Log ("Component heating %8.3e\n", xtot * py_wind_delta);
     JM130624 v76: Removed py_wind_delta term as gives too high heating */
  Log ("Component heating %8.3e\n", xtot);

  if (ochoice)
    {
      write_array (filename, ochoice);
    }
  return (0);

}

/* A summary of adiabatic cooling */
int
adiabatic_cooling_summary (w, rootname, ochoice)
     WindPtr w;
     char rootname[];
     int ochoice;
{
  int n;
  double tot;
  double adiabatic_cooling ();
  char filename[LINELENGTH];
  double t_e;



  tot = 0.0;
  for (n = 0; n < NDIM2; n++)
    {
      aaa[n] = 0;
      if (w[n].vol > 0.0)
	{
	  t_e = plasmamain[w[n].nplasma].t_e;
	  // ksl - I could not determine what the next line was supposed to do
	  // num_recomb (&plasmamain[w[n].nplasma], t_e);
	  tot += aaa[n] = adiabatic_cooling (&w[n], t_e);
	}
    }

if (geo.adiabatic==1)
	{
  display ("Adiabatic cooling");
  Log ("The total adiabatic cooling is %8.2g\n", tot);
	}
else
	{
  display ("This is only potential adiabatic cooling - it was switched off in the model");
  Log ("The total adiabatic cooling is %8.2g\n", tot);
	}

  if (ochoice)
    {
      strcpy (filename, rootname);
      strcat (filename, ".adiabatic");
      write_array (filename, ochoice);
    }
  return (0);

}


/* A summary of the lum of a cell.  The routine
 * simply reads variables that are already contained
 * in the plasma structure */

int
lum_summary (w, rootname, ochoice)
     WindPtr w;
     char rootname[];
     int ochoice;
{
  int n;
  double x, xtot;
  char c;
  char name[LINELENGTH];
  char filename[LINELENGTH];
  int nplasma;

  c = 't';
  strcpy (filename, rootname);

  printf
    ("lum tot=t,r=rad, lines=l,f=ff,b=fb,h=hydrogen,i=he1,j=he2,z=heavy elements\n");
  rdchar ("Choice", &c);
  switch (c)
    {
    case 't':			/* Total luminosity */
      strcpy (name, "Total Luminosity");
      strcat (filename, ".tot_lum");
      break;
    case 'r':			/* Radiative energo loss total */
      strcpy (name, "Total Radiative Luminosity");
      strcat (filename, ".rad_lum");
      break;
    case 'f':			/* Radiative energo loss total */
      strcpy (name, "Free Free Luminosity");
      strcat (filename, ".ff_lum");
      break;
    case 'b':			/* Radiative energo loss total */
      strcpy (name, "Free Bound (Total Recombination) Luminosity");
      strcat (filename, ".fb_lum");
      break;
    case 'l':			/* Line lumosity */
      strcpy (name, "Line Luminosity");
      strcat (filename, ".line_lum");
      break;
    case 'h':			/* H lumosity */
      strcpy (name, "H  Recombination Luminosity");
      strcat (filename, ".h1_recomb_lum");
      break;
    case 'i':			/* Line lumosity */
      strcpy (name, "He I Recombination Luminosity");
      strcat (filename, ".he1_recom_lum");
      break;
    case 'j':			/* Line lumosity */
      strcpy (name, "He 2 Recombination Luminosoity");
      strcat (filename, ".he2_recomb_lum");
      break;
    case 'z':			/* Line lumosity */
      strcpy (name, "Metal Recombination Luminosity");
      strcat (filename, ".z_recomb_lum");
      break;
    default:
      printf ("Not a valid choice\n");
      return (0);
    }



  xtot = 0.0;


  for (n = 0; n < NDIM2; n++)
    {
      x = 0.0;
      if (w[n].vol > 0.0)
	{
	  nplasma = w[n].nplasma;
	  switch (c)
	    {
	    case 't':		/* Total luminosity */
	      x = plasmamain[nplasma].lum_rad_ioniz;
	      break;
	    case 'r':		/* Radiative energo loss total */
	      x = plasmamain[nplasma].lum_rad_ioniz;
	      break;
	    case 'f':		/* Radiative energo loss total */
	      x = plasmamain[nplasma].lum_ff_ioniz;
	      break;
	    case 'b':		/* Radiative energo loss total */
	      x = plasmamain[nplasma].lum_fb_ioniz;
	      break;
	    case 'l':		/* Line luminosity */
	      x = plasmamain[nplasma].lum_lines_ioniz;
	      break;
	    case 'h':		/* H luminosity */
	      x = plasmamain[nplasma].lum_ion[0];
	      break;
	    case 'i':		/* Line luminosity */
	      x = plasmamain[nplasma].lum_ion[2];
	      break;
	    case 'j':		/* Line luminosity */
	      x = plasmamain[nplasma].lum_ion[3];
	      break;
	    case 'z':		/* Line luminosity */
	      x = plasmamain[nplasma].lum_z_ioniz;
	      break;
	    default:
	      printf ("Not a valid choice\n");
	      return (0);
	    }
	  xtot += x;

	}
      aaa[n] = x;
    }

  display (name);

  /* Log ("Luminosity total %8.3e\n", xtot * py_wind_delta);
     JM130624 v76: Removed py_wind_delta term as gives too high heating */
  Log ("Luminosity total %8.3e\n", xtot);

  if (ochoice)
    {
      write_array (filename, ochoice);
    }
  return (0);

}

/* A summary of the number of photoionizations in  a cell */

int
photo_summary (w, rootname, ochoice)
     WindPtr w;
     char rootname[];
     int ochoice;
{
  int n, ion;
  char filename[LINELENGTH];
  int nplasma;

  ion = 0;

  rdint ("ion(<5)", &ion);
  if (!(0 <= ion && ion < NIONIZ))
    {
      Log ("ion out of range\n");
      return (0);
    }
  for (n = 0; n < NDIM2; n++)
    {
      aaa[n] = 0;
      if (w[n].vol > 0.0)
	{
	  nplasma = w[n].nplasma;
	  aaa[n] =
	    plasmamain[nplasma].ioniz[ion] * plasmamain[nplasma].density[ion];
	}
    }
  display ("No of ionizations per second in cell");

  if (ochoice)
    {
      strcpy (filename, rootname);
      strcat (filename, ".photo_sum");
    }
  return (0);

}

/* A summary of the number of recombinations in  a cell */

int
recomb_summary (w, rootname, ochoice)
     WindPtr w;
     char rootname[];
     int ochoice;
{
  int n;
  int ion;
  char filename[LINELENGTH];
  int nplasma;

  ion = 0;
  rdint ("ion(<5)", &ion);
  if (!(0 <= ion && ion < NIONIZ))
    {
      Log ("recomb_summary: ioniz out of range\n");
      return (0);
    }
  for (n = 0; n < NDIM2; n++)
    {
      aaa[n] = 0;
      if (w[n].vol > 0.0)
	{
	  nplasma = w[n].nplasma;
	  num_recomb (&plasmamain[nplasma], plasmamain[nplasma].t_e);
	  aaa[n] = plasmamain[nplasma].recomb[ion];
	}
    }

  display (" Number of recombinations per second in cell");

  if (ochoice)
    {
      strcpy (filename, rootname);
      strcat (filename, ".recomb");
      write_array (filename, ochoice);
    }
  return (0);

}


/* A summary of electron densities */

int
electron_summary (w, rootname, ochoice)
     WindPtr w;
     char rootname[];
     int ochoice;
{
  int n;
  char filename[LINELENGTH];
  int nplasma;


  for (n = 0; n < NDIM2; n++)
    {
      aaa[n] = 0.0;
      if (w[n].vol > 0.0)
	{
	  nplasma = w[n].nplasma;
	  aaa[n] = plasmamain[nplasma].ne;
	}
    }
  display ("Electron densities");
  if (ochoice)
    {
      strcpy (filename, rootname);
      strcat (filename, ".ne");
      write_array (filename, ochoice);
    }
  return (0);

}


/* A summary of rho 
111002	ksl	Added to try to diagnose what was going 
		on with the torus
*/

int
rho_summary (w, rootname, ochoice)
     WindPtr w;
     char rootname[];
     int ochoice;
{
  int n;
  char filename[LINELENGTH];
  int nplasma;


  for (n = 0; n < NDIM2; n++)
    {
      aaa[n] = 0.0;
      if (w[n].vol > 0.0)
	{
	  nplasma = w[n].nplasma;
	  aaa[n] = plasmamain[nplasma].rho;
	}
    }
  display ("Rho (gm/cm**3)");
  if (ochoice)
    {
      strcpy (filename, rootname);
      strcat (filename, ".rho");
      write_array (filename, ochoice);
    }
  return (0);

}


/* A summary of rho 

Note that because of limitations in the way that display 
works cell numbers greater than 99 are not displayed as
integers unfortunately

111002	ksl	Added to try to diagnose what was going 
		on with the torus
*/

int
plasma_cell (w, rootname, ochoice)
     WindPtr w;
     char rootname[];
     int ochoice;
{
  int n;
  char filename[LINELENGTH];
  int nplasma;


  for (n = 0; n < NDIM2; n++)
    {
      aaa[n] = 0.0;
      if (w[n].vol > 0.0)
	{
	  nplasma = w[n].nplasma;
	  aaa[n] = nplasma;
	}
    }
  display ("Plasma cell number");
  if (ochoice)
    {
      strcpy (filename, rootname);
      strcat (filename, ".pnum");
      write_array (filename, ochoice);
    }
  return (0);

}




/* A summary of the average frequency */

int
freq_summary (w, rootname, ochoice)
     WindPtr w;
     char rootname[];
     int ochoice;
{
  int n;
  char filename[LINELENGTH];
  int nplasma;

  for (n = 0; n < NDIM2; n++)
    {
      aaa[n] = 0;
      if (w[n].vol > 0.0)
	{
	  nplasma = w[n].nplasma;
	  aaa[n] = plasmamain[nplasma].ave_freq;
	}
    }
  display ("Average freqency");
      if (ochoice)
	{
	  strcpy (filename, rootname);
	  strcat (filename, ".ave_freq");
	  write_array (filename, ochoice);

	}

  return (0);

}


/* A summary of the number of photons which passed through a cell.
 *
 * 111002	ksl	Modified to be able to display photons of
 * 			various types
 * */

int
nphot_summary (w, rootname, ochoice)
     WindPtr w;
     char rootname[];
     int ochoice;
{
  int n;
  char filename[LINELENGTH];
  int nplasma;
  int ichoice;
  char string[LINELENGTH];



  ichoice = 0;
  while (rdint
	 ("nphot(all=0,star=1,bl=2,disk=3,wind=4,agn=5,other=return)",
	  &ichoice) != EOF)
    {
      if (ichoice < 0 || ichoice > 5)
	{
	  return (0);
	}

      for (n = 0; n < NDIM2; n++)
	{
	  aaa[n] = 0;
	  if (w[n].vol > 0.0)
	    {
	      nplasma = w[n].nplasma;
	      if (ichoice == 0)
		{
		  aaa[n] = plasmamain[nplasma].ntot;
		  strcpy (string, "Nphot tot per cell");
		}
	      else if (ichoice == 1)
		{
		  aaa[n] = plasmamain[nplasma].ntot_star;
		  strcpy (string, "Nphot star per cell");
		}
	      else if (ichoice == 2)
		{
		  aaa[n] = plasmamain[nplasma].ntot_bl;
		  strcpy (string, "Nphot bl per cell");
		}
	      else if (ichoice == 3)
		{
		  aaa[n] = plasmamain[nplasma].ntot_disk;
		  strcpy (string, "Nphot disk per cell");
		}
	      else if (ichoice == 4)
		{
		  aaa[n] = plasmamain[nplasma].ntot_wind;
		  strcpy (string, "Nphot wind per cell");
		}
	      else if (ichoice == 5)
		{
		  aaa[n] = plasmamain[nplasma].ntot_agn;
		  strcpy (string, "Nphot agn per cell");
		}
	      else
		{
		  Error ("Unknown choice, try again\n");
		}
	    }
	}
      display (string);

      if (ochoice)
	{
	  strcpy (filename, rootname);
	  strcat (filename, ".nphot");
	  write_array (filename, ochoice);

	}
    }
  return (0);




}

/* A summary of the temperatures */

int
temp_summary (w, rootname, ochoice)
     WindPtr w;
     char rootname[];
     int ochoice;
{
  int n;
  char filename[LINELENGTH];
  int nplasma;

  for (n = 0; n < NDIM2; n++)
    {
      aaa[n] = 0;
      if (w[n].vol > 0.0)
	{
	  nplasma = w[n].nplasma;
	  aaa[n] = plasmamain[nplasma].t_e;
	}
    }
  display ("T_e");

  if (ochoice)
    {
      strcpy (filename, rootname);
      strcat (filename, ".te");
      write_array (filename, ochoice);
    }
  return (0);

}


int
temp_rad (w, rootname, ochoice)
     WindPtr w;
     char rootname[];
     int ochoice;
{
  int n;
  char filename[LINELENGTH];
  int nplasma;


  for (n = 0; n < NDIM2; n++)
    {
      aaa[n] = 0;
      if (w[n].vol > 0.0)
	{
	  nplasma = w[n].nplasma;
	  aaa[n] = plasmamain[nplasma].t_r;
	}
    }
  display ("T_rad");

  if (ochoice)
    {
      strcpy (filename, rootname);
      strcat (filename, ".tr");
      write_array (filename, ochoice);

    }
  return (0);

}

/* A summary of the radiative weights */

int
weight_summary (w, rootname, ochoice)
     WindPtr w;
     char rootname[];
     int ochoice;
{
  int n;
  int nplasma;
  char filename[LINELENGTH];

  for (n = 0; n < NDIM2; n++)
    {
      aaa[n] = 0;
      if (w[n].vol > 0.0)
	{
	  nplasma = w[n].nplasma;
	  aaa[n] = plasmamain[nplasma].w;
	}
    }
  display ("Radiative weights");

  if (ochoice)
    {
      strcpy (filename, rootname);
      strcat (filename, ".w");
      write_array (filename, ochoice);

    }

  return (0);

}

/* A summary of the velocities */

int
velocity_summary (w, rootname, ochoice)
     WindPtr w;
     char rootname[];
     int ochoice;
{
  int n;
  double x;
  char name[LINELENGTH];
  char filename[LINELENGTH];
  int ichoice;

  ichoice = 0;

  /* Note that EOF actually causes an exit via rdpar */
  while (rdint ("|v|=0;v_x=1,v_y=2,v_z=3,return=other", &ichoice) != EOF)
    {
      if (ichoice < 0 || ichoice > 3)
	return (0);

      if (ichoice == 0)
	strcpy (name, "Velocity");
      else if (ichoice == 1)
	strcpy (name, "V_x");
      else if (ichoice == 2)
	strcpy (name, "V_y");
      else
	strcpy (name, "V_z");

      for (n = 0; n < NDIM2; n++)
	{
	  x = 0;
	  if (w[n].vol > 0.0)
	    {
	      if (ichoice == 0)
		x =
		  sqrt (w[n].v[0] * w[n].v[0] + w[n].v[1] * w[n].v[1] +
			w[n].v[2] * w[n].v[2]);
	      else if (ichoice == 1)
		x = w[n].v[0];
	      else if (ichoice == 2)
		x = w[n].v[1];
	      else
		x = w[n].v[2];
	    }
	  aaa[n] = x;
	}
      display (name);


      if (ochoice)
	{
	  strcpy (filename, rootname);

	  if (ichoice == 0)
	    {
	      strcat (filename, ".vel");
	      write_array (filename, ochoice);
	    }
	  if (ichoice == 1)
	    {
	      strcat (filename, ".vx");
	      write_array (filename, ochoice);
	    }
	  if (ichoice == 2)
	    {
	      strcat (filename, ".vy");
	      write_array (filename, ochoice);
	    }
	  if (ichoice == 3)
	    {
	      strcat (filename, ".vz");
	      write_array (filename, ochoice);
	    }
	}
    }

  return (0);

}


/* A summary of the velocities */

int
mo_summary (w, rootname, ochoice)
     WindPtr w;
     char rootname[];
     int ochoice;
{
  int n;
  int ichoice;
  char name[LINELENGTH];
  char filename[LINELENGTH];
  double x;
  PlasmaPtr xplasma;


  ichoice = 0;

  /* Note that EOF actually causes an exit via rdpar */
  while (rdint
	 ("|F_rad|=0;Frac_x=1,Frad_y=2,Frad_z=3,return=other",
	  &ichoice) != EOF)
    {
      if (ichoice < 0 || ichoice > 3)
	return (0);

      if (ichoice == 0)
	strcpy (name, "|F_rad|");
      else if (ichoice == 1)
	strcpy (name, "F_rad_x");
      else if (ichoice == 2)
	strcpy (name, "F_rad_y");
      else
	strcpy (name, "F_rad_z");

      for (n = 0; n < NDIM2; n++)
	{
	  xplasma = &plasmamain[w[n].nplasma];
	  x = 0;
	  if (w[n].vol > 0.0)
	    {
	      if (ichoice == 0)
		x =
		  sqrt (xplasma->dmo_dt[0] * xplasma->dmo_dt[0] +
			xplasma->dmo_dt[1] * xplasma->dmo_dt[1] +
			xplasma->dmo_dt[2] * xplasma->dmo_dt[2]);
	      else if (ichoice == 1)
		x = xplasma->dmo_dt[0];
	      else if (ichoice == 2)
		x = xplasma->dmo_dt[1];
	      else
		x = xplasma->dmo_dt[2];
	    }
	  aaa[n] = x;
	}
      display (name);

  if (ochoice && ichoice == 0)
    {
      strcpy (filename, rootname);
      strcat (filename, ".f_rad");
      write_array (filename, ochoice);
    }

  }
  return (0);

}

/* A summary of the volumes of each cell 

	080811	ksl	Add lines from Stuart's version of this
			routine to bring this version into 
			compliance with it.
*/

int
vol_summary (w, rootname, ochoice)
     WindPtr w;
     char rootname[];
     int ochoice;
{
  int n;

  char filename[LINELENGTH];

  for (n = 0; n < NDIM2; n++)
    {
      aaa[n] = 0;
      if (w[n].vol > 0.0)
	{
	  aaa[n] = w[n].vol;
	}
    }
  display ("Volumes");

  if (ochoice)
    {
      strcpy (filename, rootname);
      strcat (filename, ".vol");
      write_array (filename, ochoice);
    }

  return (0);

}


/***********************************************************
                                       Space Telescope Science Institute

Synopsis:
	wind_element provides a detailed look at a single cell

Arguments:


Returns:

Description:

	Note that this routine does not calculate anything.  It
	simple reports on variables in the wind array.
	
		
Notes:

History:
	080804	ksl	60b -- Added reporting of partition function
	101106	ksl	69 -- Changed variable name in rdint to make
			it more obvious what was being discused here.

**************************************************************/

int
wind_element (w)
     WindPtr w;
{
  PlasmaPtr xplasma;
  int m, n, i, j, nn, mm;
  int first, last;
  n = 50;
a: printf("There are %i wind elements in this model\n",NDIM2);
rdint ("Wind.array.element", &n);

  if (n < 0)
    goto b;
  else if (n > NDIM2)
	{
	printf("No, there are %i wind elements, not %i\n",NDIM2,n);
	goto a;
	}

  wind_n_to_ij (n, &i, &j);
  xplasma = &plasmamain[w[n].nplasma];

  Log
    ("Element %d (%d,%d)  inwind %d plasma cell %d ntot %d nioniz %d nrad %d\n",
     n, i, j, w[n].inwind, xplasma->nplasma, xplasma->ntot, xplasma->nioniz,
     xplasma->nrad);
  Log ("xyz %8.2e %8.2e %8.2e vel %8.2e %8.2e %8.2e\n", w[n].x[0], w[n].x[1],
       w[n].x[2], w[n].v[0], w[n].v[1], w[n].v[2]);
  Log ("nh %8.2e ne %8.2e t_r %8.2e t_e %8.2e w %8.2e vol %8.2e\n",
       xplasma->rho * rho2nh, xplasma->ne, xplasma->t_r, xplasma->t_e,
       xplasma->w, w[n].vol);

  if (w[n].inwind < 0)
    Log ("\n# Cell is not inwind, expect all zeros to follow\n\n");
  /*70d - added compton - ksl */
  /*70g compton removed from luminosity reporting, it is now a cooling mechanism but does not produce photons
     DR cooling also added in to report */
  Log
    ("t_e %8.2e cool_tot %8.2e lum_lines  %8.2e lum_ff  %8.2e lum_fb     %8.2e cool_comp %8.2e cool_adiab %8.2e cool_DR %8.2e \n",
     xplasma->t_e, xplasma->lum_rad_ioniz+xplasma->lum_comp_ioniz+xplasma->lum_adiabatic_ioniz+xplasma->lum_dr_ioniz, xplasma->lum_lines_ioniz, xplasma->lum_ff_ioniz, xplasma->lum_fb_ioniz, xplasma->lum_comp_ioniz, xplasma->lum_adiabatic_ioniz, xplasma->lum_dr_ioniz);
  Log
    ("t_r %8.2e heat_tot %8.2e heat_lines %8.2e heat_ff %8.2e heat_photo %8.2e heat_comp %8.2e heat_icomp %8.2e\n",
     xplasma->t_r, xplasma->heat_tot, xplasma->heat_lines, xplasma->heat_ff,
     xplasma->heat_photo, xplasma->heat_comp,xplasma->heat_ind_comp);



  Log ("Recombination cooling   HII>HI %8.2e HeII>HeI %8.2e HeIII>HeII %8.2e Metals %8.2e\n",xplasma->lum_ion[0], xplasma->lum_ion[2],
     xplasma->lum_ion[3], xplasma->lum_z);
  Log ("Photoionization heating HI>HII %8.2e HeI>HeII %8.2e HeII>HeIII %8.2e Metals %8.2e\n",xplasma->heat_ion[0], xplasma->heat_ion[2],
     xplasma->heat_ion[3], xplasma->heat_z);



  Log ("The ratio of rad (total) cooling to heating is %8.2f (%8.2f) \n",
       xplasma->lum_rad_ioniz / xplasma->heat_tot,
       (xplasma->lum_rad_ioniz + xplasma->lum_adiabatic_ioniz + xplasma->lum_comp_ioniz +
	xplasma->lum_dr_ioniz) / xplasma->heat_tot);
  Log ("Adiabatic cooling %8.2e is %8.2g of total cooling\n",
       xplasma->lum_adiabatic_ioniz,
       xplasma->lum_adiabatic_ioniz / (xplasma->lum_rad + xplasma->lum_adiabatic + xplasma->lum_comp_ioniz + xplasma->lum_dr_ioniz));
  /*70g NSH compton and DR cooling are now reported seperately. */
  Log ("Compton cooling   %8.2e is %8.2g of total cooling\n",
       xplasma->lum_comp_ioniz,
       xplasma->lum_comp_ioniz / (xplasma->lum_rad + xplasma->lum_adiabatic + xplasma->lum_comp_ioniz + xplasma->lum_dr_ioniz));
  Log ("DR cooling        %8.2e is %8.2g of total cooling\n", xplasma->lum_dr_ioniz,
       xplasma->lum_dr_ioniz / (xplasma->lum_rad + xplasma->lum_adiabatic + xplasma->lum_comp_ioniz + xplasma->lum_dr_ioniz));
  Log ("Number of ionizing photons in cell nioniz %d\n", xplasma->nioniz);
  Log ("Log Ionization parameter in this cell cell based %4.2f ferland %4.2f\n", log10 (xplasma->ip), log10 (xplasma->ferland_ip));	//70h NSH computed ionizaion parameter
  Log ("ioniz %8.2e %8.2e %8.2e %8.2e %8.2e\n",
       xplasma->ioniz[0], xplasma->ioniz[1], xplasma->ioniz[2],
       xplasma->ioniz[3], xplasma->ioniz[4]);
  Log
    ("Convergence status: whole %d converging %d t_r %8.2e t_e %8.2e hc %8.2e \n",
     xplasma->converge_whole, xplasma->converging, xplasma->converge_t_r,
     xplasma->converge_t_e, xplasma->converge_hc);

  Log ("Densities:\n");
  for (nn = 0; nn < 5; nn++)
    {
      first = ele[nn].firstion;
      last = first + ele[nn].nions;
      Log ("%-5s ", ele[nn].name);
      for (m = first; m < last; m++)
	Log (" %8.2e", xplasma->density[m]);
      Log ("\n");
    }


  Log ("Partition function:\n");
  for (nn = 0; nn < 5; nn++)
    {
      first = ele[nn].firstion;
      last = first + ele[nn].nions;
      Log ("%-5s ", ele[nn].name);
      for (m = first; m < last; m++)
	Log (" %8.2e", xplasma->partition[m]);
      Log ("\n");
    }


  Log ("Non LTE Level densities:\n");


  mm = 0;
  for (nn = 0; nn < 10; nn++)
    {
      while (ion[mm].nlte <= 0 && mm < nions)
	mm++;
      if (mm == nions)
	break;

      first = ion[mm].first_levden;
      last = first + ion[mm].nlte;
      Log ("ion %3d %3d", ion[mm].z, ion[mm].istate);
      for (m = first; m < last; m++)
	Log (" %8.2e", xplasma->levden[m]);
      Log ("\n");
      mm++;
    }

  Log ("Spectral model details:\n");
  for (nn = 0; nn < geo.nxfreq; nn++)
    {
      Log ("numin= %8.2e (%8.2e) numax= %8.2e (%8.2e) Model= %2d PL_log_w= %9.2e PL_alpha= %9.2e Exp_w= %9.2e EXP_temp= %9.2e\n",xplasma-> fmin_mod[nn],geo.xfreq[nn],xplasma->fmax_mod[nn],geo.xfreq[nn+1],xplasma->spec_mod_type[nn],xplasma->pl_log_w[nn],xplasma->pl_alpha[nn],xplasma->exp_w[nn],xplasma->exp_temp[nn]);
}    


  goto a;

b:return (0);

}

int
tau_h_summary (w, rootname, ochoice)
     WindPtr w;
     char rootname[];
     int ochoice;
{
  int n;
  int nplasma;

  for (n = 0; n < NDIM2; n++)
    {
      aaa[n] = 0;
      if (w[n].vol > 0.0)
	{
	  nplasma = w[n].nplasma;
	  aaa[n] =
	    6.e-18 * plasmamain[nplasma].density[0] * pow (w[n].vol, 0.333);

	}
    }
  display ("tau_h summary");

  return (0);

}

/* A summary of the volumes of each cell */

int
coolheat_summary (w, rootname, ochoice)
     WindPtr w;
     char rootname[];
     int ochoice;
{
  int n;
  int nplasma;

  for (n = 0; n < NDIM2; n++)
    {
      aaa[n] = 0;
      if (w[n].vol > 0.0)
	{
	  nplasma = w[n].nplasma;
	  aaa[n] = plasmamain[nplasma].lum_rad_ioniz / plasmamain[nplasma].heat_tot;

	}
    }
  display ("Cooling over heating");

  return (0);

}

int
complete_file_summary (w, root, ochoice)
     WindPtr w;
     char root[];
     int ochoice;
{
  temp_summary (w, root, ochoice);
  temp_rad (w, root, ochoice);
  electron_summary (w, root, ochoice);
  convergence_summary (w, root, ochoice);
  ion_summary (w, 6, 3, 0, root, ochoice);
  ion_summary (w, 6, 4, 0, root, ochoice);
  ion_summary (w, 6, 5, 0, root, ochoice);
  ion_summary (w, 7, 4, 0, root, ochoice);
  ion_summary (w, 7, 5, 0, root, ochoice);
  ion_summary (w, 7, 6, 0, root, ochoice);
  ion_summary (w, 8, 4, 0, root, ochoice);
  ion_summary (w, 8, 5, 0, root, ochoice);
  ion_summary (w, 8, 6, 0, root, ochoice);
  ion_summary (w, 8, 7, 0, root, ochoice);
  ion_summary (w, 14, 3, 0, root, ochoice);
  ion_summary (w, 14, 4, 0, root, ochoice);
  ion_summary (w, 14, 5, 0, root, ochoice);

  /* Before 68c, this is what we printed out */
  ion_summary (w, 6, 3, 1, root, ochoice);
  ion_summary (w, 6, 4, 1, root, ochoice);
  ion_summary (w, 6, 5, 1, root, ochoice);
  ion_summary (w, 7, 4, 1, root, ochoice);
  ion_summary (w, 7, 5, 1, root, ochoice);
  ion_summary (w, 7, 6, 1, root, ochoice);
  ion_summary (w, 8, 4, 1, root, ochoice);
  ion_summary (w, 8, 5, 1, root, ochoice);
  ion_summary (w, 8, 6, 1, root, ochoice);
  ion_summary (w, 8, 7, 1, root, ochoice);
  ion_summary (w, 14, 3, 1, root, ochoice);
  ion_summary (w, 14, 4, 1, root, ochoice);
  ion_summary (w, 14, 5, 1, root, ochoice);


  ion_summary (w, 6, 3, 2, root, ochoice);
  ion_summary (w, 6, 4, 2, root, ochoice);
  ion_summary (w, 6, 5, 2, root, ochoice);
  ion_summary (w, 7, 4, 2, root, ochoice);
  ion_summary (w, 7, 5, 2, root, ochoice);
  ion_summary (w, 7, 6, 2, root, ochoice);
  ion_summary (w, 8, 4, 2, root, ochoice);
  ion_summary (w, 8, 5, 2, root, ochoice);
  ion_summary (w, 8, 6, 2, root, ochoice);
  ion_summary (w, 8, 7, 2, root, ochoice);
  ion_summary (w, 14, 3, 2, root, ochoice);
  ion_summary (w, 14, 4, 2, root, ochoice);
  ion_summary (w, 14, 5, 2, root, ochoice);


  ion_summary (w, 6, 3, 3, root, ochoice);
  ion_summary (w, 6, 4, 3, root, ochoice);
  ion_summary (w, 6, 5, 3, root, ochoice);
  ion_summary (w, 7, 4, 3, root, ochoice);
  ion_summary (w, 7, 5, 3, root, ochoice);
  ion_summary (w, 7, 6, 3, root, ochoice);
  ion_summary (w, 8, 4, 3, root, ochoice);
  ion_summary (w, 8, 5, 3, root, ochoice);
  ion_summary (w, 8, 6, 3, root, ochoice);
  ion_summary (w, 8, 7, 3, root, ochoice);
  ion_summary (w, 14, 3, 3, root, ochoice);
  ion_summary (w, 14, 4, 3, root, ochoice);
  ion_summary (w, 14, 5, 3, root, ochoice);

  return (0);
}

/* A summary of the regions in the wind */

int
wind_reg_summary (w, rootname, ochoice)
     WindPtr w;
     char rootname[];
     int ochoice;
{
  int n;


  for (n = 0; n < NDIM2; n++)
    {
      aaa[n] = -9;
      if (w[n].vol > 0.0)
	{
	  aaa[n] = w[n].inwind;

	}
    }
  display ("Regions of the wind");

  return (0);

}


/* A summary of the dvds_ave */

int
dvds_summary (w, rootname, ochoice)
     WindPtr w;
     char rootname[];
     int ochoice;
{
  char filename[LINELENGTH];
  int n;

  for (n = 0; n < NDIM2; n++)
    {
      aaa[n] = 0;
      if (w[n].vol > 0.0)
	{
	  aaa[n] = w[n].dvds_ave;

	}
    }
  display ("Average dvds");
  if (ochoice)
    {
      strcpy (filename, rootname);
      strcat (filename, ".dvds");
      write_array (filename, ochoice);
    }
  return (0);
}

/* A summary of inner shell ionization */

int
inner_shell_summary (w, rootname, ochoice)
     WindPtr w;
     char rootname[];
     int ochoice;
{
  int n;
  char filename[LINELENGTH];
  int nplasma;


  for (n = 0; n < NDIM2; n++)
    {
      aaa[n] = 0.0;
      if (w[n].vol > 0.0)
	{
	  nplasma = w[n].nplasma;
	  aaa[n] = plasmamain[nplasma].gamma_inshl[0];
	}
    }
  display ("Inner_shell");
  if (ochoice)
    {
      strcpy (filename, rootname);
      strcat (filename, ".is");
      write_array (filename, ochoice);
    }
  return (0);

}


/* A summary of the Ionization parameter - might not always be present */

int
IP_summary (w, rootname, ochoice)
     WindPtr w;
     char rootname[];
     int ochoice;
{
  int n;
  char filename[LINELENGTH];
  int nplasma;

  for (n = 0; n < NDIM2; n++)
    {
      aaa[n] = 0;
      if (w[n].vol > 0.0)
	{
	  nplasma = w[n].nplasma;
	  aaa[n] = ((plasmamain[nplasma].ferland_ip));
	}
    }
  display ("Ionization parameter (Ferland)");

  if (ochoice)
    {
      strcpy (filename, rootname);
      strcat (filename, ".f_IP");
      write_array (filename, ochoice);

    }

  for (n = 0; n < NDIM2; n++)
    {
      aaa[n] = 0;
      if (w[n].vol > 0.0)
	{
	  nplasma = w[n].nplasma;
	  aaa[n] = ((plasmamain[nplasma].ip));
	}
    }
  display ("Ionization parameter");

  if (ochoice)
    {
      strcpy (filename, rootname);
      strcat (filename, ".IP");
      write_array (filename, ochoice);

    }

  for (n = 0; n < NDIM2; n++)
    {
      aaa[n] = 0;
      if (w[n].vol > 0.0)
	{
	  nplasma = w[n].nplasma;
	  aaa[n] = ((plasmamain[nplasma].ip_direct));
	}
    }
  display ("Log Ionization parameter (direct)");

  if (ochoice)
    {
      strcpy (filename, rootname);
      strcat (filename, ".IP_direct");
      write_array (filename, ochoice);

    }

 for (n = 0; n < NDIM2; n++)
    {
      aaa[n] = 0;
      if (w[n].vol > 0.0)
	{
	  nplasma = w[n].nplasma;
	  aaa[n] = ((plasmamain[nplasma].ip_scatt));
	}
    }
  display ("Log Ionization parameter (scattered)");

  if (ochoice)
    {
      strcpy (filename, rootname);
      strcat (filename, ".IP_scatt");
      write_array (filename, ochoice);

    }




  return (0);

}



/* A summary of the Sim alpha parameter - might not always be present.

   1108	ksl	Adapted for new version of banded alpha
   1208 nsh	Changed names - reference to sim removed to make way for exponential estimators
   1401 nsh	Added more information - this now writes out all the spectral model parameters
 */

int
alpha_summary (w, rootname, ochoice)
     WindPtr w;
     char rootname[];
     int ochoice;
{
  int n, m;
  char filename[LINELENGTH];
  int nplasma;
  char word[LINELENGTH];

 

  for (m = 0; m < geo.nxfreq; m++)
    {
      for (n = 0; n < NDIM2; n++)
	{
	  aaa[n] = 0;
	  if (w[n].vol > 0.0)
	    {
	      nplasma = w[n].nplasma;
	      aaa[n] = plasmamain[nplasma].pl_alpha[m];
	    }
	}

      strcpy (word, "");
      sprintf (word, ".pl_alpha%03d", m);
      display (word);

      if (ochoice)
	{
	  strcpy (filename, rootname);
	  strcat (filename, word);
	  write_array (filename, ochoice);
	}
    }

 for (m = 0; m < geo.nxfreq; m++)
    {
      for (n = 0; n < NDIM2; n++)
	{
	  aaa[n] = 0;
	  if (w[n].vol > 0.0)
	    {
	      nplasma = w[n].nplasma;
	      aaa[n] = plasmamain[nplasma].pl_log_w[m];
	    }
	}

      strcpy (word, "");
      sprintf (word, ".pl_log_w%03d", m);
      display (word);

      if (ochoice)
	{
	  strcpy (filename, rootname);
	  strcat (filename, word);
	  write_array (filename, ochoice);
	}
    }


 for (m = 0; m < geo.nxfreq; m++)
    {
      for (n = 0; n < NDIM2; n++)
	{
	  aaa[n] = 0;
	  if (w[n].vol > 0.0)
	    {
	      nplasma = w[n].nplasma;
	      aaa[n] = plasmamain[nplasma].exp_temp[m];
	    }
	}

      strcpy (word, "");
      sprintf (word, ".exp_temp%03d", m);
      display (word);

      if (ochoice)
	{
	  strcpy (filename, rootname);
	  strcat (filename, word);
	  write_array (filename, ochoice);
	}
    }

 for (m = 0; m < geo.nxfreq; m++)
    {
      for (n = 0; n < NDIM2; n++)
	{
	  aaa[n] = 0;
	  if (w[n].vol > 0.0)
	    {
	      nplasma = w[n].nplasma;
	      aaa[n] = plasmamain[nplasma].exp_w[m];
	    }
	}

      strcpy (word, "");
      sprintf (word, ".exp_w%03d", m);
      display (word);

      if (ochoice)
	{
	  strcpy (filename, rootname);
	  strcat (filename, word);
	  write_array (filename, ochoice);
	}
    }

 for (m = 0; m < geo.nxfreq; m++)
    {
      for (n = 0; n < NDIM2; n++)
	{
	  aaa[n] = 0;
	  if (w[n].vol > 0.0)
	    {
	      nplasma = w[n].nplasma;
	      aaa[n] = plasmamain[nplasma].spec_mod_type[m];
	    }
	}

      strcpy (word, "");
      sprintf (word, ".spec_mod_type%03d", m);
      display (word);

      if (ochoice)
	{
	  strcpy (filename, rootname);
	  strcat (filename, word);
	  write_array (filename, ochoice);
	}
    }



 for (m = 0; m < geo.nxfreq; m++)
    {
      for (n = 0; n < NDIM2; n++)
	{
	  aaa[n] = 0;
	  if (w[n].vol > 0.0)
	    {
	      nplasma = w[n].nplasma;
	      aaa[n] = plasmamain[nplasma].fmin_mod[m];
	    }
	}

      strcpy (word, "");
      sprintf (word, ".fmin_mod%03d", m);
      display (word);

      if (ochoice)
	{
	  strcpy (filename, rootname);
	  strcat (filename, word);
	  write_array (filename, ochoice);
	}
    }

for (m = 0; m < geo.nxfreq; m++)
    {
      for (n = 0; n < NDIM2; n++)
	{
	  aaa[n] = 0;
	  if (w[n].vol > 0.0)
	    {
	      nplasma = w[n].nplasma;
	      aaa[n] = plasmamain[nplasma].fmax_mod[m];
	    }
	}

      strcpy (word, "");
      sprintf (word, ".fmax_mod%03d", m);
      display (word);

      if (ochoice)
	{
	  strcpy (filename, rootname);
	  strcat (filename, word);
	  write_array (filename, ochoice);
	}
    }

  return (0);

}



/* A summary of frequency banded radiation density - a crude spectrum for a cell.

 1301 nsh	Written
 */

int
J_summary (w, rootname, ochoice)
     WindPtr w;
     char rootname[];
     int ochoice;
{
  int i, n;
  char filename[LINELENGTH];
  char number[2];
  int nplasma;


  i = 1;
  rdint ("Band number for J", &i);

  for (n = 0; n < NDIM2; n++)
    {
      aaa[n] = 0;
      if (w[n].vol > 0.0)
	{
	  nplasma = w[n].nplasma;
	  aaa[n] = (plasmamain[nplasma].xj[i]);
	}
    }
  display ("J in cell");
  printf ("i=%i", i);
  sprintf (number, "%i", i);
  if (ochoice)
    {
      strcpy (filename, rootname);
      strcat (filename, ".J_band");
      strcat (filename, number);
      write_array (filename, ochoice);

    }
  return (0);

}


int
J_scat_summary (w, rootname, ochoice)
    WindPtr w;
     char rootname[];
     int ochoice;
{
  int  n;
  char filename[LINELENGTH];
  int nplasma;



  for (n = 0; n < NDIM2; n++)
    {
      aaa[n] = 0;
      if (w[n].vol > 0.0)
	{
	  nplasma = w[n].nplasma;
	  aaa[n] = (plasmamain[nplasma].j);
	}
    }
  display ("J in cell");
 
  if (ochoice)
    {
      strcpy (filename, rootname);
      strcat (filename, ".J_tot");
      write_array (filename, ochoice);
    }

  for (n = 0; n < NDIM2; n++)
    {
      aaa[n] = 0;
      if (w[n].vol > 0.0)
	{
	  nplasma = w[n].nplasma;
	  aaa[n] = (plasmamain[nplasma].j_direct);
	}
    }
  display ("J in cell from direct photons");
  if (ochoice)
    {
      strcpy (filename, rootname);
      strcat (filename, ".J_direct");
      write_array (filename, ochoice);
    }

  for (n = 0; n < NDIM2; n++)
    {
      aaa[n] = 0;
      if (w[n].vol > 0.0)
	{
	  nplasma = w[n].nplasma;
	  aaa[n] = (plasmamain[nplasma].j_scatt);
	}
    }
  display ("J in cell from scattered photons");

  if (ochoice)
    {
      strcpy (filename, rootname);
      strcat (filename, ".J_scatt");
      write_array (filename, ochoice);
    }


  return (0);

}






//Split of photons from different sources in the cell.

int
phot_split (w, rootname, ochoice)
     WindPtr w;
     char rootname[];
     int ochoice;
{
  int n;
  char filename[LINELENGTH];
  int nplasma;

  for (n = 0; n < NDIM2; n++)
    {
      aaa[n] = 0;
      if (w[n].vol > 0.0)
	{
	  nplasma = w[n].nplasma;
	  aaa[n] = (plasmamain[nplasma].ntot_wind);
	}
    }
  display ("Wind photons in cell");

  if (ochoice)
    {
      strcpy (filename, rootname);
      strcat (filename, ".ntot_wind");
      write_array (filename, ochoice);

    }
  for (n = 0; n < NDIM2; n++)
    {
      aaa[n] = 0;
      if (w[n].vol > 0.0)
	{
	  nplasma = w[n].nplasma;
	  aaa[n] = (plasmamain[nplasma].ntot_agn);
	}
    }
  display ("AGN photons in cell");

  if (ochoice)
    {
      strcpy (filename, rootname);
      strcat (filename, ".ntot_agn");
      write_array (filename, ochoice);

    }
  for (n = 0; n < NDIM2; n++)
    {
      aaa[n] = 0;
      if (w[n].vol > 0.0)
	{
	  nplasma = w[n].nplasma;
	  aaa[n] = (plasmamain[nplasma].ntot_disk);
	}
    }
  display ("Disk photons in cell");

  if (ochoice)
    {
      strcpy (filename, rootname);
      strcat (filename, ".ntot_disk");
      write_array (filename, ochoice);

    }
  for (n = 0; n < NDIM2; n++)
    {
      aaa[n] = 0;
      if (w[n].vol > 0.0)
	{
	  nplasma = w[n].nplasma;
	  aaa[n] = (plasmamain[nplasma].ntot_star);
	}
    }
  display ("Stellar photons in cell");

  if (ochoice)
    {
      strcpy (filename, rootname);
      strcat (filename, ".ntot_star");
      write_array (filename, ochoice);

    }
  return (0);

}

int
thompson (w, rootname, ochoice)
     WindPtr w;
     char rootname[];
     int ochoice;
{
  int n;
  int nplasma;
  double ne;
  char filename[LINELENGTH];

  for (n = 0; n < NDIM2; n++)
    {
      aaa[n] = 0;
      if (w[n].vol > 0.0)
	{
	  nplasma = w[n].nplasma;
	  ne = plasmamain[nplasma].ne;
	  aaa[n] = (THOMPSON * ne) * pow (w[n].vol, 0.333);
	}
    }
  display ("Thompson optical depths");

  if (ochoice)
    {
      strcpy (filename, rootname);
      strcat (filename, ".thomp");
      write_array (filename, ochoice);
    }

  return (0);

}



int
nscat_split (w, rootname, ochoice)
     WindPtr w;
     char rootname[];
     int ochoice;
{
  int n;
  int nplasma;
  char filename[LINELENGTH];

  for (n = 0; n < NDIM2; n++)
    {
      aaa[n] = 0;
      if (w[n].vol > 0.0)
	{
	  nplasma = w[n].nplasma;
	  aaa[n] = plasmamain[nplasma].nscat_es;
	}
    }
  display ("Thompson scatters");

  for (n = 0; n < NDIM2; n++)
    {
      aaa[n] = 0;
      if (w[n].vol > 0.0)
	{
	  nplasma = w[n].nplasma;
	  aaa[n] = plasmamain[nplasma].nscat_res;
	}
    }
  display ("Resonant scatters");






  if (ochoice)
    {
      strcpy (filename, rootname);
      strcat (filename, ".thomp");
      write_array (filename, ochoice);
    }

  return (0);

}

int
convergence_summary (w, rootname, ochoice)
     WindPtr w;
     char rootname[];
     int ochoice;
{
  int n;
  int nplasma;
  char filename[LINELENGTH];

  for (n = 0; n < NDIM2; n++)
    {
      aaa[n] = 0;
      if (w[n].vol > 0.0)
	{
	  nplasma = w[n].nplasma;
	  aaa[n] = plasmamain[nplasma].converge_whole;
	}
    }
  display
    ("Convergence (0=converged.  Higher numbers indicate one or more convergence tests failed)");

  if (ochoice)
    {
      strcpy (filename, rootname);
      strcat (filename, ".conv");
      write_array (filename, ochoice);
    }

  return (0);

}


/* 

   1112	ksl	Write out arrays of use in evaluating what is going on in convergence of
   		various models
*/

int
convergence_all (w, rootname, ochoice)
     WindPtr w;
     char rootname[];
     int ochoice;
{
  int n;
  int nplasma;
  char filename[LINELENGTH];


 

  for (n = 0; n < NDIM2; n++)
    {
      aaa[n] = 0;
      if (w[n].vol > 0.0)
	{
	  nplasma = w[n].nplasma;
	  aaa[n] = plasmamain[nplasma].converge_t_e;
	}
    }
  display ("t_e Convergence (0=converged)");

  if (ochoice)
    {
      strcpy (filename, rootname);
      strcat (filename, ".conv_te");
      write_array (filename, ochoice);
    }

  for (n = 0; n < NDIM2; n++)
    {
      aaa[n] = 0;
      if (w[n].vol > 0.0)
	{
	  nplasma = w[n].nplasma;
	  aaa[n] = plasmamain[nplasma].converge_hc;
	}
    }
  display ("hc Convergence (0=converged)");

  if (ochoice)
    {
      strcpy (filename, rootname);
      strcat (filename, ".conv_hc");
      write_array (filename, ochoice);
    }

  for (n = 0; n < NDIM2; n++)
    {
      aaa[n] = 0;
      if (w[n].vol > 0.0)
	{
	  nplasma = w[n].nplasma;
	  aaa[n] = plasmamain[nplasma].converge_t_r;
	}
    }
  display ("t_r Convergence (0=converged)");

  if (ochoice)
    {
      strcpy (filename, rootname);
      strcat (filename, ".conv_tr");
      write_array (filename, ochoice);
    }

 for (n = 0; n < NDIM2; n++)
    {
      aaa[n] = 0;
      if (w[n].vol > 0.0)
	{
	  nplasma = w[n].nplasma;
	  aaa[n] = plasmamain[nplasma].converge_whole;
	}
    }
  display ("Convergence (0=converged)");

  if (ochoice)
    {
      strcpy (filename, rootname);
      strcat (filename, ".conv_whole");
      write_array (filename, ochoice);
    }


  return (0);
}

/* 

   1306	nsh	Write out information pertaining to the models used in each cell/frequency band
*/
  
int
model_bands (w, rootname, ochoice)
     WindPtr w;
     char rootname[];
     int ochoice;
{
  int n, m;
  int nplasma;
  char filename[LINELENGTH];
  char word[LINELENGTH];



  for (m = 0; m < geo.nxfreq; m++)
    {
      for (n = 0; n < NDIM2; n++)
	{
	  aaa[n] = 0;
	  if (w[n].vol > 0.0)
	    {
	      nplasma = w[n].nplasma;
	      aaa[n] = plasmamain[nplasma].nxtot[m];
	    }
	}

      strcpy (word, "");
      sprintf (word, ".nxtot%03d", m);
      display (word);

      if (ochoice)
	{
	  strcpy (filename, rootname);
	  strcat (filename, word);
	  write_array (filename, ochoice);
	}
    }

  if (geo.nxfreq == 0)
    {
      printf ("This model does not use bands calculating the ionization\n");
      return (0);
    }


  /* Now write out the value of xj in each band */

  for (m = 0; m < geo.nxfreq; m++)
    {
      for (n = 0; n < NDIM2; n++)
	{
	  aaa[n] = 0;
	  if (w[n].vol > 0.0)
	    {
	      nplasma = w[n].nplasma;
	      aaa[n] = plasmamain[nplasma].xj[m];
	    }
	}

      strcpy (word, "");
      sprintf (word, ".xj%03d", m);
      display (word);

      if (ochoice)
	{
	  strcpy (filename, rootname);
	  strcat (filename, word);
	  write_array (filename, ochoice);
	}
    }


  /* Now write out the value of xave_freq in each band */

  for (m = 0; m < geo.nxfreq; m++)
    {
      for (n = 0; n < NDIM2; n++)
	{
	  aaa[n] = 0;
	  if (w[n].vol > 0.0)
	    {
	      nplasma = w[n].nplasma;
	      aaa[n] = plasmamain[nplasma].xave_freq[m];
	    }
	}

      strcpy (word, "");
      sprintf (word, ".xave_freq%03d", m);
      display (word);

      if (ochoice)
	{
	  strcpy (filename, rootname);
	  strcat (filename, word);
	  write_array (filename, ochoice);
	}
    }


  /* Now write out the value of nxtot in each band */

  for (m = 0; m < geo.nxfreq; m++)
    {
      for (n = 0; n < NDIM2; n++)
	{
	  aaa[n] = 0;
	  if (w[n].vol > 0.0)
	    {
	      nplasma = w[n].nplasma;
	      aaa[n] = plasmamain[nplasma].nxtot[m];
	    }
	}

      strcpy (word, "");
      sprintf (word, ".nxtot%03d", m);
      display (word);

      if (ochoice)
	{
	  strcpy (filename, rootname);
	  strcat (filename, word);
	  write_array (filename, ochoice);
	}
    }

  return (0);
}



/* A summary of adiabatic cooling */
int
heatcool_summary (w, rootname, ochoice)
     WindPtr w;
     char rootname[];
     int ochoice;
{
  int n;
  int nplasma;
  char filename[LINELENGTH];
	float x;

x=wind_luminosity(0.0,1e20);

  for (n = 0; n < NDIM2; n++)
    {
      aaa[n] = 0;
      if (w[n].vol > 0.0)
	{
	  nplasma = w[n].nplasma;
	  aaa[n] = plasmamain[nplasma].heat_tot;
	 if (w[n].div_v < 0.0) // add in if it is negative and hence a heating term
		{
		aaa[n] += -1.0*(plasmamain[nplasma].lum_adiabatic_ioniz);
		}
	}
    }
  display ("Total heating");

  if (ochoice)
    {
      strcpy (filename, rootname);
      strcat (filename, ".heat_tot");
      write_array (filename, ochoice);
    }

  for (n = 0; n < NDIM2; n++)
    {
      aaa[n] = 0;
      if (w[n].vol > 0.0)
	{
	  nplasma = w[n].nplasma;
	  aaa[n] = plasmamain[nplasma].heat_lines;
	}
    }
  display ("Line heating");

  if (ochoice)
    {
      strcpy (filename, rootname);
      strcat (filename, ".heat_lines");
      write_array (filename, ochoice);
    }

  for (n = 0; n < NDIM2; n++)
    {
      aaa[n] = 0;
      if (w[n].vol > 0.0)
	{
	  nplasma = w[n].nplasma;
	  aaa[n] = plasmamain[nplasma].heat_ff;
	}
    }
  display ("FF heating");

  if (ochoice)
    {
      strcpy (filename, rootname);
      strcat (filename, ".heat_ff");
      write_array (filename, ochoice);
    }

  for (n = 0; n < NDIM2; n++)
    {
      aaa[n] = 0;
      if (w[n].vol > 0.0)
	{
	  nplasma = w[n].nplasma;
	  aaa[n] = plasmamain[nplasma].heat_comp;
	}
    }
  display ("Compton heating");

  if (ochoice)
    {
      strcpy (filename, rootname);
      strcat (filename, ".heat_comp");
      write_array (filename, ochoice);
    }

  for (n = 0; n < NDIM2; n++)
    {
      aaa[n] = 0;
      if (w[n].vol > 0.0)
	{
	  nplasma = w[n].nplasma;
	  aaa[n] = plasmamain[nplasma].heat_ind_comp;
	}
    }
  display ("Induced Compton heating");

  if (ochoice)
    {
      strcpy (filename, rootname);
      strcat (filename, ".heat_ind_comp");
      write_array (filename, ochoice);
    }

  for (n = 0; n < NDIM2; n++)
    {
      aaa[n] = 0;
      if (w[n].vol > 0.0)
	{
	  nplasma = w[n].nplasma;
	  aaa[n] = plasmamain[nplasma].heat_photo;
	}
    }
  display ("Photo heating");

  if (ochoice)
    {
      strcpy (filename, rootname);
      strcat (filename, ".heat_photo");
      write_array (filename, ochoice);
    }

  for (n = 0; n < NDIM2; n++)
    {
      aaa[n] = 0;
      if (w[n].vol > 0.0)
	{
	  nplasma = w[n].nplasma;
	  aaa[n] = plasmamain[nplasma].lum_lines_ioniz;
	}
    }
  display ("Line Luminosity");

  if (ochoice)
    {
      strcpy (filename, rootname);
      strcat (filename, ".lum_lines");
      write_array (filename, ochoice);
    }



  for (n = 0; n < NDIM2; n++)
    {
      aaa[n] = 0;
      if (w[n].vol > 0.0)
	{
	  nplasma = w[n].nplasma;
	  aaa[n] = plasmamain[nplasma].lum_adiabatic_ioniz;
	}
    }
  display ("Adiabatic Luminosity");

  if (ochoice)
    {
      strcpy (filename, rootname);
      strcat (filename, ".lum_adiabatic");
      write_array (filename, ochoice);
    }

  for (n = 0; n < NDIM2; n++)
    {
      aaa[n] = 0;
      if (w[n].vol > 0.0)
	{
	  nplasma = w[n].nplasma;
	  aaa[n] = plasmamain[nplasma].lum_ff_ioniz;
	}
    }
  display ("Free Free Luminosity");

  if (ochoice)
    {
      strcpy (filename, rootname);
      strcat (filename, ".lum_ff");
      write_array (filename, ochoice);
    }

  for (n = 0; n < NDIM2; n++)
    {
      aaa[n] = 0;
      if (w[n].vol > 0.0)
	{
	  nplasma = w[n].nplasma;
	  aaa[n] = plasmamain[nplasma].lum_comp_ioniz;
	}
    }
  display ("Compton Luminosity");

  if (ochoice)
    {
      strcpy (filename, rootname);
      strcat (filename, ".lum_comp");
      write_array (filename, ochoice);
    }

  for (n = 0; n < NDIM2; n++)
    {
      aaa[n] = 0;
      if (w[n].vol > 0.0)
	{
	  nplasma = w[n].nplasma;
	  aaa[n] = plasmamain[nplasma].lum_dr_ioniz;
	}
    }
  display ("DR Luminosity");

  if (ochoice)
    {
      strcpy (filename, rootname);
      strcat (filename, ".lum_dr");
      write_array (filename, ochoice);
    }

  for (n = 0; n < NDIM2; n++)
    {
      aaa[n] = 0;
      if (w[n].vol > 0.0)
	{
	  nplasma = w[n].nplasma;
	  aaa[n] = plasmamain[nplasma].lum_fb_ioniz;
	}
    }
  display ("FB Luminosity");

  if (ochoice)
    {
      strcpy (filename, rootname);
      strcat (filename, ".lum_fb");
      write_array (filename, ochoice);
    }





 for (n = 0; n < NDIM2; n++)
    {
      aaa[n] = 0;
      if (w[n].vol > 0.0)
	{
	  nplasma = w[n].nplasma;
	  aaa[n] =
	    plasmamain[nplasma].lum_fb_ioniz + plasmamain[nplasma].lum_dr_ioniz +
	    plasmamain[nplasma].lum_comp_ioniz + plasmamain[nplasma].lum_ff_ioniz +
	    plasmamain[nplasma].lum_lines_ioniz;
	 if (w[n].div_v >= 0.0) //only add in if it is treated as a cooling term
		{
		aaa[n] += plasmamain[nplasma].lum_adiabatic_ioniz;
		}
	}
    }
  display ("Total Cooling");

  if (ochoice)
    {
      strcpy (filename, rootname);
      strcat (filename, ".lum_total");
      write_array (filename, ochoice);
    }


for (n = 0; n < NDIM2; n++)
    {
      aaa[n] = 0;
      if (w[n].vol > 0.0)
	{
	  nplasma = w[n].nplasma;
	  aaa[n] =plasmamain[nplasma].lum_rad_ioniz;
	}
    }
  display ("Total Radiating Luminosity");

  if (ochoice)
    {
      strcpy (filename, rootname);
      strcat (filename, ".lum_rad");
      write_array (filename, ochoice);
    }

  return (0);
}
