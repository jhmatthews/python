
/***********************************************************
                                       Space Telescope Science Institute

 Synopsis:
 
 	This file contains various subroutines of py_wind.  It is not part of python!
 	
	ion_summary (w, element, istate, iswitch, rootname, ochoice) calculates and 
	displays information for a specific element and ionization state.  The information
	displayed depends upon the value of iswitch  
	

Arguments:		
	WindPtr w;
	int element,istate;
	char filename[]			Output file name, if "none" then nothing is written;
	int iswitch                     0 = ion fraction
					1 = ion density
					2 = number of scatters 
     	char rootname[];                rootname of the output file
     	int ochoice;			The screen display is always the actual value. If ochoice
					is non-zere an output file is written.

Returns:
 
Description:	
	
		
Notes:

History:
 	97jun	ksl	Coding on py_wind began.
	01sep	ksl	Added switch to ion summary so that ion density rather than ion fraction could
			be printed out.  iswitch==0 --> fraction, otherwise density
	01dec	ksl	Updated for new calling structure for two_level_atom & scattering_fraction.
	04nov	ksl	Updated for multiple coordinate systems, i.e to use the routine display. 
			Note that many of these summaries seem quite dated, and I have really not 
			checked what they are supposed to do.
	05apr	ksl	Eliminated MDIM references
	090125	ksl	Added capability to display the number of scatters of an ion taking
			place in a cell.
**************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>


#include "atomic.h"
#include "python.h"


int 
ion_summary (WindPtr w, int element, int istate, int iswitch, char rootname[], int ochoice)
{
  int nion, nelem;
  int n;
  char choice[LINELENGTH], iname[LINELENGTH];
  char name[LINELENGTH];
  char filename[LINELENGTH];
  double x;
  int nplasma;



/* Find the CIV ion */
  nion = 0;
  while (nion < nions
	 && !(ion[nion].z == element && ion[nion].istate == istate))
    nion++;
  if (nion == nions)
    {
      Log ("Error--element %d ion %d not found in define_wind\n", element,
	   istate);
      return (-1);
    }
  nelem = 0;
  while (nelem < nelements && ele[nelem].z != element)
    nelem++;

  strcpy (name, "");


  for (n = 0; n < NDIM2; n++)
    {
      aaa[n] = 0;
      nplasma = w[n].nplasma;
      if (w[n].vol > 0.0 && plasmamain[nplasma].ne > 1.0)
	{
	  if (iswitch == 0)
	    {
  sprintf (name, "Element %d (%s) ion %d fractions\n", element,
	   ele[nelem].name, istate);
	      aaa[n] = plasmamain[nplasma].density[nion];
	      aaa[n] /=
		((plasmamain[nplasma].density[0] +
		  plasmamain[nplasma].density[1]) * ele[nelem].abun);
	    }
	  else if (iswitch == 1)
	    {
  sprintf (name, "Element %d (%s) ion %d density\n", element,
	   ele[nelem].name, istate);
	      aaa[n] = plasmamain[nplasma].density[nion];
	    }
	  else if (iswitch == 2)
	    {
  sprintf (name, "Element %d (%s) ion %d  #scatters\n", element,
	   ele[nelem].name, istate);
	      aaa[n] = plasmamain[nplasma].scatters[nion];
	    }
	  else if (iswitch == 3)
	    {
  sprintf (name, "Element %d (%s) ion %d scattered flux\n", element,
	   ele[nelem].name, istate);
	      aaa[n] = plasmamain[nplasma].xscatters[nion];
	    }
	  else
	    {
	      Error ("ion_summary : Unknown switch %d \n", iswitch);
	      exit (0);
	    }
	}
    }

  display (name);




  /* Store the appropriate values in a place where it does not matter */
  if (ochoice)
    {
      for (n = 0; n < NDIM2; n++)
	{
	  nplasma = w[n].nplasma;
	  if (w[n].vol > 0.0 && plasmamain[nplasma].ne > 1.0)
	    {
	      if (iswitch == 0)
		x /=
		  ((plasmamain[nplasma].density[0] +
		    plasmamain[nplasma].density[1]) * ele[nelem].abun);
	      else if (iswitch == 1)
		{
		  x = plasmamain[nplasma].density[nion];
		  x = log10 (x);
		}
	      else if (iswitch == 2)
		{
		  x = plasmamain[nplasma].scatters[nion];
		}
	      else if (iswitch == 3)
		{
		  x = plasmamain[nplasma].xscatters[nion];
		}
	      else
		{
		  Error ("ion_summary : Unknown switch %d \n", iswitch);
		  exit (0);
		}


	    }
	  else
	    x = 0.0;
	  w[n].x[1] = x;

	}

      strcpy (filename, rootname);
      if (iswitch == 0)
	strcpy (choice, ".ion");
      else if (iswitch == 1)
	strcpy (choice, ".ionc");
      else if (iswitch == 2)
	strcpy (choice, ".ions");
      else if (iswitch == 3)
	strcpy (choice, ".iona");
      else
	{
	  Error ("ion_summary : Unknown switch %d \n", iswitch);
	  exit (0);
	}

      strcat (choice, ele[nelem].name);
      sprintf (iname, "%d", istate);
      strcat (choice, iname);

      strcat (filename, choice);
      write_array (filename, ochoice);
    }

  return (0);
}

int 
tau_ave_summary (WindPtr w, int element, int istate, double freq, char rootname[], int ochoice)
{
  int nion, nelem;
  int n;
  char choice[LINELENGTH], iname[LINELENGTH];
  char name[LINELENGTH];
  char filename[LINELENGTH];
  int nplasma;



/* Find the CIV ion */
  nion = 0;
  while (nion < nions
	 && !(ion[nion].z == element && ion[nion].istate == istate))
    nion++;
  if (nion == nions)
    {
      Log ("Error--element %d ion %d not found in define_wind\n", element,
	   istate);
      return (-1);
    }
  nelem = 0;
  while (nelem < nelements && ele[nelem].z != element)
    nelem++;

  strcpy (name, "");
  sprintf (name, "Element %d (%s) ion %d fractions\n", element,
	   ele[nelem].name, istate);

  for (n = 0; n < NDIM2; n++)
    {
      aaa[n] = 0;
      if (w[n].vol > 0.0)
	{
	  nplasma = w[n].nplasma;
	  aaa[n] =
	    PI_E2_OVER_M * plasmamain[nplasma].density[nion] / freq /
	    w[n].dvds_ave;
	}
    }

  display (name);

  /* Store the appropriate values in a place where it does not matter */
  if (ochoice)
    {
      for (n = 0; n < NDIM2; n++)
	{
	  if (w[n].vol > 0.0)
	    {
	      nplasma = w[n].nplasma;
	      w[n].x[1] =
		plasmamain[nplasma].density[nion] / (0.87 *
						     plasmamain[nplasma].ne *
						     ele[nelem].abun);
	    }
	  else
	    w[n].x[1] = 0;
	}

      strcpy (filename, rootname);

      strcpy (choice, ".ion");
      strcat (choice, ele[nelem].name);
      sprintf (iname, "%d", istate);
      strcat (choice, iname);

      strcat (filename, choice);

      write_array (filename, ochoice);
    }

  return (0);
}

int 
line_summary (WindPtr w, int element, int istate, char rootname[], int ochoice)
{
  int nion, nelem;
  int n;
  double x;
  char choice[LINELENGTH], iname[LINELENGTH];
  char name[LINELENGTH];
  char filename[LINELENGTH];
  int nline;
  double freq_search, dd;

  double d1, d2, z, energy_c4, rb, tot, omega;
  int nplasma;


  element = 6;
  istate = 4;
  energy_c4 = HC / (1550e-8);

/* Find the CIV ion */
  nion = 0;
  while (nion < nions
	 && !(ion[nion].z == element && ion[nion].istate == istate))
    nion++;
  if (nion == nions)
    {
      Log ("Error--element %d ion %d not found in define_wind\n", element,
	   istate);
      return (-1);
    }
  nelem = 0;
  while (nelem < nelements && ele[nelem].z != element)
    nelem++;

/* Find the CIV line in the data */
  nline = 0;
  freq_search = C / 1548.1949e-8;

  while (fabs (1. - lin_ptr[nline]->freq / freq_search) > 0.0001
	 && nline < nlines)
    nline++;
  if (nline == nlines)
    {
      Error ("line_summary: Could not find line in linelist\n");
      exit (0);
    }

  rdint ("line_transfer(0=pure.abs,1=pure.scat,2=sing.scat,3=escape.prob)",
	 &geo.line_mode);
  if (geo.line_mode == 0)
    Log ("Pure_abs in line heating/cooling\n");
  else if (geo.line_mode == 1)
    Log ("Pure_scat_in line heating/cooling\n");
  else if (geo.line_mode == 2)
    Log ("Single scat for line heating/cooling\n");
  else if (geo.line_mode == 3)
    Log ("Escape probabilities for line heating/cooling\n");
  else
    {
      Log ("Unknown line mode\n");
      return (0);
    }

  strcpy (name, "");
  sprintf (name, "Luminosity %d (%s) ion %d fractions\n", element,
	   ele[nelem].name, istate);

  tot = 0;
  for (n = 0; n < NDIM2; n++)
    {
      aaa[n] = 0;
      if (w[n].vol > 0.0)
	{
	  nplasma = w[n].nplasma;
	  dd = plasmamain[nplasma].density[lin_ptr[nline]->nion];
	  two_level_atom (lin_ptr[nline], &plasmamain[nplasma], &d1, &d2);
	  x =
	    (d2) * a21 (lin_ptr[nline]) * H * lin_ptr[nline]->freq * w[n].vol;
	  x *= z = scattering_fraction (lin_ptr[nline], &plasmamain[nplasma]);

	  tot += x;
	  aaa[n] = x;
	}
    }

  display (name);

  tot = 2. * tot;		// Why is there a factor of 2 here??? ksl

  Log ("The total CIV luminosity (flux) is %8.2g (%8.2g)\n",
       tot, tot / (4 * PI * 1e4 * PC * PC));


  /* Store the appropriate values in a place where it does not matter */
  if (ochoice)
    {
      for (n = 0; n < NDIM2; n++)
	{
	  // Here is the calculation of the effective collisions strength
	  if (w[n].vol > 0.0)
	    {
	      nplasma = w[n].nplasma;
	      omega = 5.13 * pow (plasmamain[nplasma].t_e / 1.e5, 0.18);
	      rb =
		8.629e-6 * exp (-energy_c4 /
				(BOLTZMANN * plasmamain[nplasma].t_e)) /
		sqrt (plasmamain[nplasma].t_e) * omega;
	      w[n].x[1] =
		plasmamain[nplasma].density[nion] * plasmamain[nplasma].ne *
		rb * energy_c4 * w[n].vol;
	    }
	  else
	    w[n].x[1] = 0;
	}

      strcpy (filename, rootname);
      strcpy (choice, ".line");
      strcat (choice, ele[nelem].name);
      sprintf (iname, "%d", istate);
      strcat (choice, iname);

      strcat (filename, choice);
      write_array (filename, ochoice);
    }

  return (0);
}

int 
total_emission_summary (WindPtr w, char rootname[], int ochoice)
{
  double tot;
  int n;
  double total_emission ();
  char filename[LINELENGTH];


  tot = 0.0;
  for (n = 0; n < NDIM2; n++)
    {
      aaa[n] = 0;
      if (w[n].vol > 0.0)
	{
	  tot += aaa[n] = total_emission (&w[n], 0.0, 1e50);
	}
    }

  display ("Calculated thermal luminosities of cells");

  Log ("The total luminosity due to thermal process is %8.2g\n", tot);



  if (ochoice)
    {
      strcpy (filename, rootname);
      strcat (filename, ".abs");
      write_array (filename, ochoice);
    }


  return (0);
}

int 
modify_te (WindPtr w, char rootname[], int ochoice)
{
  int n;
  double x;
  double t_e, calc_te ();
  char filename[LINELENGTH];
  int nplasma;


  for (n = 0; n < NDIM2; n++)
    {
      nplasma = w[n].nplasma;
      aaa[n] = 0;
      if (w[n].vol > 0.0 && (x = plasmamain[nplasma].heat_tot) > 1.0)
	{
	  aaa[n] = t_e =
	    calc_te (&plasmamain[nplasma], TMIN,
		     1.2 * plasmamain[nplasma].t_r);
	}
    }

  display ("Calculate t_e");

  if (ochoice)
    {
      strcpy (filename, rootname);
      strcat (filename, ".te_mod");
      write_array (filename, ochoice);
    }

  return (0);
}



int 
partial_measure_summary (WindPtr w, int element, int istate, char rootname[], int ochoice)
{
  int nion, nelem;
  int n;
  double total;
  char choice[LINELENGTH], iname[LINELENGTH];
  char name[LINELENGTH];
  char filename[LINELENGTH];
  int nplasma;



/* Find the CIV ion */
  total = 0;
  nion = 0;
  while (nion < nions
	 && !(ion[nion].z == element && ion[nion].istate == istate))
    nion++;
  if (nion == nions)
    {
      Log ("Error--element %d ion %d not found in define_wind\n", element,
	   istate);
      return (-1);
    }
  nelem = 0;
  while (nelem < nelements && ele[nelem].z != element)
    nelem++;

  strcpy (name, "");
  sprintf (name, "Element %d (%s) ion %d fractions", element, ele[nelem].name,
	   istate);

  for (n = 0; n < NDIM2; n++)
    {
      aaa[n] = 0;
      nplasma = w[n].nplasma;
      if (plasmamain[nplasma].ne > 1.0 && w[n].vol > 0.0)
	{
	  total += aaa[n] =
	    plasmamain[nplasma].density[nion] * plasmamain[nplasma].ne *
	    w[n].vol;
	}
    }

  display ("Calculate t_e");

  Log ("The partial emission measure ne*nion*V is %8.2g\n", total);

  /* Store the appropriate values in a place where it does not matter */
  if (ochoice)
    {
      for (n = 0; n < NDIM2; n++)
	{
	  nplasma = w[n].nplasma;
	  if (plasmamain[nplasma].ne > 1.0 && w[n].vol > 0.0)
	    w[n].x[1] =
	      plasmamain[nplasma].density[nion] / (0.87 *
						   plasmamain[nplasma].ne *
						   ele[nelem].abun);
	  else
	    w[n].x[1] = 0;
	}

      strcpy (filename, rootname);
      strcpy (choice, ".part");
      strcat (choice, ele[nelem].name);
      sprintf (iname, "%d", istate);
      strcat (choice, iname);

      strcat (filename, choice);
      write_array (filename, ochoice);
    }

  return (0);
}
