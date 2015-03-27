
/* The subroutines in this file handle allocating, incrementing, and writing out the
   spectrum arrays for Python */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "atomic.h"

#include "python.h"

/***********************************************************
                                       Space Telescope Science Institute

 Synopsis:

	int spectrum_init(f1,f2,nangle,angle,phase,scat_select,top_bot_select,select_extract) 
	allocates space for and initializes the spectrum arrays

Arguments:		
	double f1,f2;			the minimimum and maximum frequency in the spectrum
	int nangle;			the number of different angles and phases for which 
					spectra must be created
	double angle[],phase[];		the corresponding inclination angles and orbital phases
	int scat_select[]		a code which allows limit the photons which will be summed
					to specific numbers of scatters
					>999 -> select everything
					0 < scat_select < MAXSCAT -> select only those photons 
					which have scattered nscat times, 
					scat_select<0 -> select those phtons with more than |nscat| 
					scatters
					This parallels the angle array
	int top_bot_select		a code which allows one to limit whether all photons or just those
					photons above or below the disk are selected
						0 -> all
						>0 -> only those above the disk
						<0-> only those below the disk 
 	int select_extract		0 for Live or Die option, extract option otherwise
 	
Returns:
 
Description:	
	The first time spectrum_init  is called (i.e. if ispec_start=0), it allocates memory 
	for the various spectrum arrays.  (This is done one time, so one needs to allocate 
	space for all the arrays even though they are not all used in the ionization step).  
	The total number of spectra created is nangle+MSPEC.
	
	Each time spectrum_init is called it rezeros all the spectrum information and 
	calculates the other information associated with each spectrum, such as the
	angle cosines and the names of each spectrum.
		
Notes:
	angle[],phase[] and scat_select[] only apply to the spectra extracted at 
	specific angles.  angle[0], phase[0], and scat_select[0] all affect spec[3]

?? I have a suspicion that these wavelengths are off by half a bin in one direction or the other ???

	Warning - Do not put anything in this routine that does anything but initialize
	or reinitialize the spectrum structure s. This is important because this routine
	is not accessed if one is continuing an old calculation of the detailed spectrum.
	It is still use on a restart where the detailed spectral cycles have not begun
	because in that case the spectra are not saved.

History:
 	97jan   ksl	Coded and debugged as part of Python effort. 
 	97jul	ksl	Updated to allow extraction of photons at specific phases
 	97aug	ksl	Updated to allow spectra to be created which sum only a
 			specific number of scatters.  (Basically this just involved
 			populating s[].nscat.  Currently it has no effect on live or
 			die option. 
 	97sep21	ksl	Provided for renormalization of spectrum through s[n].renorm.
 			This should only affect which affect Live or die option.  
	99dec29	ksl	Began modifications intended to allow the creation of additional
			spectra from the varioous sources of photons, the star, the disk
			and the wind.  Replaced the old fixed number 3 with MSPEC.  Also
			elimininated some checks which should no longer be needed.
	05apr	ksl	55d -- Eliminated initialization of shell, since no longer used
	13feb	nsh	74b5 -- Included lines to initialize the log spectrum
	1409	ksl	Added another spectrum, the spectrum of generated photons. This
			is the first spectrum in the structure

**************************************************************/

int i_spec_start = 0;

int
spectrum_init (f1, f2, nangle, angle, phase, scat_select, top_bot_select,
	       select_extract, rho_select, z_select, az_select, r_select)
     double f1, f2;
     int nangle;
     double angle[], phase[];
     int scat_select[], top_bot_select[];
     int select_extract;
     double rho_select[], z_select[], az_select[], r_select[];
{
  int i, n;
  int nspec;
  double freqmin, freqmax, dfreq;
  double lfreqmin, lfreqmax, ldfreq;	/* NSH 1302 Min, max and delta for the log spectrum */
  double x1, x2;
  char dummy[20];

  freqmin = f1;
  freqmax = f2;
  dfreq = (freqmax - freqmin) / NWAVE;

  nspec = nangle + MSPEC;

/* NSH 1302 Lines to set up a logarithmic spectrum */

  lfreqmin = log10 (freqmin);
  lfreqmax = log10 (freqmax);
  ldfreq = (lfreqmax - lfreqmin) / NWAVE;

  /* Create the spectrum arrays the first time routine is called */
  if (i_spec_start == 0)
    {
      xxspec = calloc (sizeof (spectrum_dummy), nspec);
      if (xxspec == NULL)
	{
	  Error
	    ("spectrum_init: Could not allocate memory for %d spectra with %d wavelengths\n",
	     nspec, NWAVE);
	  exit (0);
	}

      nspectra = nspec;		/* Note that nspectra is a global variable */

      i_spec_start = 1;		/* This is to prevent reallocation of the same arrays on multiple calls to spectrum_init */
    }


  Log_silent ("Zeroing or rezeroing all %d spectral matrices\n", nspec);


  for (i = 0; i <= MAXSCAT; i++)
    nscat[i] = nres[i] = 0;

  for (i = 0; i < NSTAT; i++)
    nstat[i] = 0;

  for (n = 0; n < nspec; n++)
    {
      xxspec[n].freqmin = freqmin;
      xxspec[n].freqmax = freqmax;
      xxspec[n].dfreq = dfreq;
      xxspec[n].lfreqmin = lfreqmin;
      xxspec[n].lfreqmax = lfreqmax;
      xxspec[n].ldfreq = ldfreq;
      for (i = 0; i < NSTAT; i++)
	xxspec[n].nphot[i] = 0;
      for (i = 0; i < NWAVE; i++)
	{
	  xxspec[n].f[i] = 0;
	  xxspec[n].lf[i] = 0;	/* NSH 1302 zero the logarithmic spectra */
	}
    }

  strcpy (xxspec[0].name, "Created");
  strcpy (xxspec[1].name, "Emitted");
  strcpy (xxspec[2].name, "CenSrc");
  strcpy (xxspec[3].name, "Disk   ");
  strcpy (xxspec[4].name, "Wind   ");
  strcpy (xxspec[5].name, "HitSurf");
  strcpy (xxspec[6].name, "Scattered");
  for (n = 0; n < MSPEC; n++)
    {
      xxspec[n].lmn[0] = xxspec[n].lmn[1] = xxspec[n].lmn[2] = 0.;
      xxspec[n].renorm = 1.0;
    }

  for (n = MSPEC; n < nspec; n++)
    {
/* 
We want to set up the direction cosines for extractions.  We have to be careful
about the sense of the orbital phase within the program.  Viewed from the "north"
pole in a system with the secondary along the x axis the observer moves clockwise
and phases just before 0 should be in the +x + y quadrant (since we have the disk and
wind rotating counter clockwize as viewed from the top.  Another way of saying this
is at phases just before 0, e.g. 0.98, the observer sees the receeding side of the 
disk. The minus sign in the terms associated with phase are to make this happen.
02feb ksl
*/

      sprintf (xxspec[n].name, "A%02.0f", angle[n - MSPEC]);
      xxspec[n].lmn[0] =
	sin (angle[n - MSPEC] / RADIAN) * cos (-phase[n - MSPEC] * 360. /
					       RADIAN);
      xxspec[n].lmn[1] =
	sin (angle[n - MSPEC] / RADIAN) * sin (-phase[n - MSPEC] * 360. /
					       RADIAN);
      xxspec[n].lmn[2] = cos (angle[n - MSPEC] / RADIAN);
      Log_silent ("Angle %e Angle cosines:%e %e %e\n", angle[n - MSPEC],
		  xxspec[n].lmn[0], xxspec[n].lmn[1], xxspec[n].lmn[2]);

      /* Initialize variables needed for live or die option */
      x1 = angle[n - MSPEC] - DANG_LIVE_OR_DIE;
      x2 = angle[n - MSPEC] + DANG_LIVE_OR_DIE;
      if (x1 < 0.)
	x1 = 0;
      if (x2 > 180.)
	x2 = 180.;
      x1 = fabs (cos (x1 / RADIAN));
      x2 = fabs (cos (x2 / RADIAN));
      if (x1 > x2)
	{
	  xxspec[n].mmax = x1;
	  xxspec[n].mmin = x2;
	}
      else
	{
	  xxspec[n].mmax = x2;
	  xxspec[n].mmin = x1;
	}
      if (select_extract == 0)
	{
	  xxspec[n].renorm = 1. / (xxspec[n].mmax - xxspec[n].mmin);
	}
      else
	xxspec[n].renorm = 1.;
      /* Completed initialization of variables for live or die */

/* 68g - Changed the format of the string describing the spectrum so that it was shorted
 * and so that the phase was always given.  ksl 091125
 */
      strcpy (dummy, "");
      sprintf (dummy, "P%04.2f", phase[n - MSPEC]);
      strcat (xxspec[n].name, dummy);
      xxspec[n].nscat = scat_select[n - MSPEC];
      if (xxspec[n].nscat < MAXSCAT)
	{			/* Then conditions have been place on the
				   number of scatters to be included so update the names */
	  strcpy (dummy, "");
	  if (xxspec[n].nscat > MAXSCAT)
	    sprintf (dummy, "_sc:all");
	  if (xxspec[n].nscat >= 0)
	    sprintf (dummy, "_sc:%d", xxspec[n].nscat);
	  else
	    sprintf (dummy, "_sc:>%d", -xxspec[n].nscat);
	  strcat (xxspec[n].name, dummy);
	}
      xxspec[n].top_bot = top_bot_select[n - MSPEC];
      if (xxspec[n].top_bot != 0)
	{			/* Then conditions have been placed on the last
				   location of the photon so update the names */
	  strcpy (dummy, "");
	  if (xxspec[n].top_bot == 1)
	    sprintf (dummy, "_top");
	  else if (xxspec[n].top_bot == -1)
	    sprintf (dummy, "_bot");
	  else if (xxspec[n].top_bot == 2)
	    {
	      sprintf (dummy, "_pos");
	      xxspec[n].x[0] =
		rho_select[n - MSPEC] * cos (az_select[n - MSPEC] / RADIAN);
	      xxspec[n].x[1] =
		rho_select[n - MSPEC] * sin (az_select[n - MSPEC] / RADIAN);
	      xxspec[n].x[2] = z_select[n - MSPEC];
	      xxspec[n].r = r_select[n - MSPEC];

	    }
	  else
	    {
	      Error ("spectrum_init: Unknown option %d\n", xxspec[n].top_bot);
	      exit (0);

	    }
	  strcat (xxspec[n].name, dummy);
	}
    }

  return (0);
}

/***********************************************************
                                       Space Telescope Science Institute

 Synopsis:

 int spectrum_create(p,f1,f2,nangle,select_extract) increments the spectrum arrays
 	after each flight of photons is processed.
  
Arguments:		
	PhotPtr p;
	double f1,f2;			the minimimum and maximum frequency in the spectrum
	int nangle;				the number of different angles and phases for which 
							spectra must be created
	int select_extract;		The integer stating whether the Live or Die option has 
						been chosen. (0==Live or Die)
Returns:
  
Description:	

	This routine increments the total spectrum arrays based on what has happened to each 
	photon.  In the Live or Die option, the spectra at specific angles are also created here.
 
	The routine is called after each batch of photons has been transported through the wind and 
	prints some intermediate results to assure the user that the program is still running.  

Notes:

	Summing up of the spectra in the "extract" option is done in extract.c
	
	!! To create spectra in the live or die option which are selected 
	on the number of scatters individual photons undergo, there are 
	some additional changes required to this routine.  These changes 
	should parallel those now in extract under the normal option.  97aug29
	
History:
 	97jan	ksl	Coded and debugged as part of Python effort.
 	97nov23 ksl	Modified to use new error and logging routines  
	02jan	ksl	Added live or die capability to extract specific 
			scatters and above and below plane
	02jul	ksl	Reduced printing of error messages when frequency
			seemed out of bounds since due to Doppler shifts
			at the time of photon creation, the frequencies
			can exceed the range of the formal limits somewhat.
	05apr	ksl	Removed code which summed up the points where
			the last scattering occurred, as not sufficiently
			important to retain given the desire to deal with
			both 1-d and 2-d grids simulataneouly.  
	08mar	ksl	Fixed up spectrum types to account for tracking
			of photons which had been scattered by the wind
	1212	ksl	Changed the way dealt with photons which had
			frequencies which were too low or too high
			to record the numbers and to give an error
			only if the numbers seemed large
	1409	ksl	Added code that includes a new spectrum.  The
			spectrum of generated photons which escape the
			system.  This spectrum is constructed with 
			the original weights of the photons.  Photons
			which hit a hard boudary and destroyed are
			not recorded.
**************************************************************/

int
spectrum_create (p, f1, f2, nangle, select_extract)
     PhotPtr p;
     double f1, f2;
     int nangle;
     int select_extract;

{
  int nphot, i, j, k, k1, n;
  int nspec, spectype;
  double freqmin, freqmax, dfreq;
  double lfreqmin, lfreqmax, ldfreq;
  double x1;
  int wind_n_to_ij ();
  int mscat, mtopbot;
  double delta;
  double nlow, nhigh;
  int k_orig, k1_orig;

  freqmin = f1;
  freqmax = f2;
  dfreq = (freqmax - freqmin) / NWAVE;
  nspec = nangle + MSPEC;
  nlow = 0.0;			// variable to storte the number of photons that have frequencies which are too low
  nhigh = 0.0;			// variable to storte the number of photons that have frequencies which are too high
  delta = 0.0;			// fractional frequency error allowod

/* Lines to set up a logarithmic spectrum */

  lfreqmin = log10 (freqmin);
  lfreqmax = log10 (freqmax);
  ldfreq = (lfreqmax - lfreqmin) / NWAVE;


  for (nphot = 0; nphot < NPHOT; nphot++)
    {
      if ((j = p[nphot].nscat) < 0 || j > MAXSCAT)
	nscat[MAXSCAT]++;
      else
	nscat[j]++;

      if ((j = p[nphot].nrscat) < 0 || j > MAXSCAT)
	nres[MAXSCAT]++;
      else
	nres[j]++;

    /* find out where we are in log space */
      k1 = (log10 (p[nphot].freq) - log10 (freqmin)) / ldfreq;
      if (k1 < 0)
	{
	  k1 = 0;
	}
      if (k1 > NWAVE)
	{
	  k1 = NWAVE;
	}

	/* also need to work out where we are for photon's original wavelength */
      k1_orig = (log10 (p[nphot].freq_orig) - log10 (freqmin)) / ldfreq;
      if (k1_orig < 0)
	{
	  k1_orig = 0;
	}
      if (k1_orig > NWAVE)
	{
	  k1_orig = NWAVE;
	}


    /* lines to work out where we are in a normal spectrum with linear spacing */
      k = (p[nphot].freq - freqmin) / dfreq;
      if (k < 0)
	{
	  if (((1. - p[nphot].freq / freqmin) > delta) && (geo.rt_mode != 2))
	    nlow = nlow + 1;
	  k = 0;
	}
      else if (k > NWAVE - 1)
	{
	  if (((1. - freqmax / p[nphot].freq) > delta) && (geo.rt_mode != 2))
	    nhigh = nhigh + 1;
	  k = NWAVE - 1;
	}

	/* also need to work out where we are for photon's original wavelength */
	  k_orig = (p[nphot].freq_orig - freqmin) / dfreq;
	  if (k_orig < 0)
	{
	  if (((1. - p[nphot].freq_orig / freqmin) > delta) && (geo.rt_mode != 2))
	    nlow = nlow + 1;
	  k_orig = 0;
	}
      else if (k_orig > NWAVE - 1)
	{
	  if (((1. - freqmax / p[nphot].freq_orig) > delta) && (geo.rt_mode != 2))
	    nhigh = nhigh + 1;
	  k_orig = NWAVE - 1;
	}


    xxspec[0].f[k_orig] += p[nphot].w_orig;	/* created spectrum with original weights and wavelengths */
	xxspec[0].lf[k1_orig] += p[nphot].w_orig;	/* logarithmic created spectrum */
	

      if ((i = p[nphot].istat) == P_ESCAPE)
	{
      xxspec[0].nphot[i]++;
	  xxspec[1].f[k] += p[nphot].w;	/* emitted spectrum */
	  xxspec[1].lf[k1] += p[nphot].w;	/* logarithmic emitted spectrum */
	  xxspec[1].nphot[i]++;
	  spectype = p[nphot].origin;
	  if (spectype >= 10)
	    spectype -= 10;
	  if (spectype == PTYPE_STAR || spectype == PTYPE_BL || spectype == PTYPE_AGN)	// Then it came from the bl or the star
	    {
	      xxspec[2].f[k] += p[nphot].w;	/* emitted star (+bl) spectrum */
	      xxspec[2].lf[k1] += p[nphot].w;	/* logarithmic emitted star (+bl) spectrum */
	      xxspec[2].nphot[i]++;
	    }
	  else if (spectype == PTYPE_DISK)	// Then it was a disk photon 
	    {
	      xxspec[3].f[k] += p[nphot].w;	/* transmitted disk spectrum */
	      xxspec[3].lf[k1] += p[nphot].w;	/* logarithmic transmitted disk spectrum */
	      xxspec[3].nphot[i]++;
	    }
	  else if (spectype == PTYPE_WIND)
	    {
	      xxspec[4].f[k] += p[nphot].w;	/* wind spectrum */
	      xxspec[4].lf[k1] += p[nphot].w;	/* logarithmic wind spectrum */
	      xxspec[4].nphot[i]++;
	    }
	  else
	    {
	      Error ("spectrum_create: Unknown photon type %d\n", spectype);
	    }

	  /* For Live or Die option, increment the spectra here */
	  if (select_extract == 0)
	    {
	      x1 = fabs (p[nphot].lmn[2]);
	      for (n = MSPEC; n < nspec; n++)
		{
		  /* Complicated if statement to allow one to choose whether to construct the spectrum
		     from all photons or just from photons which have scattered a specific number
		     of times.  01apr13--ksl-Modified if statement to change behavior on negative numbers 
		     to say that a negative number for mscat implies that you accept any photon with 
		     |mscat| or more scatters */
		  if (((mscat = xxspec[n].nscat) > 999 ||
		       p[nphot].nscat == mscat ||
		       (mscat < 0 && p[nphot].nscat >= (-mscat)))
		      && ((mtopbot = xxspec[n].top_bot) == 0
			  || (mtopbot * p[nphot].x[2]) > 0))

		    {
		      if (xxspec[n].mmin < x1 && x1 < xxspec[n].mmax)
			xxspec[n].f[k] += p[nphot].w;
		      xxspec[n].lf[k1] += p[nphot].w;	/* logarithmic spectrum */
		    }

		}
	    }
	}
      else if (i == P_HIT_STAR || i == P_HIT_DISK)
	{
	  xxspec[5].f[k] += p[nphot].w;	/*absorbed spectrum */
	  xxspec[5].lf[k1] += p[nphot].w;	/*logarithmic absorbed spectrum */
	  xxspec[5].nphot[i]++;
	}

      if (j > 0)
	{
	  xxspec[6].f[k] += p[nphot].w;	/* j is the number of scatters so this constructs */
	  xxspec[6].lf[k1] += p[nphot].w;	/* logarithmic j is the number of scatters so this constructs */
	  if (i < 0 || i > NSTAT - 1)
	    xxspec[6].nphot[NSTAT - 1]++;
	  else
	    xxspec[6].nphot[i]++;	/* scattering spectrum */
	}

      if (i < 0 || i > NSTAT - 1)
	nstat[NSTAT - 1]++;
      else
	nstat[i]++;


    }


  if ((nlow / nphot > 0.05) || (nhigh / nphot > 0.05))
    {
      Error
	("spectrum_create: Fraction of photons lost: %4.2f wi/ freq. low, %4.2f w/freq hi\n",
	 nlow / nphot, nhigh / nphot);
    }
  else
    {
      Log
	("spectrum_create: Fraction of photons lost:  %4.2f wi/ freq. low, %4.2f w/freq hi\n",
	 nlow / nphot, nhigh / nphot);
    }


  Log ("\nNo. of photons which have scattered n times\n");
  for (i = 0; i <= MAXSCAT; i++)
    {
      Log ("%6d", nscat[i]);
      if ((i % 10) == 9)
	Log ("\n");
    }
  Log ("\nNumber of photons resonantly scattering n times\n");
  for (i = 0; i <= MAXSCAT; i++)
    {
      Log ("%6d", nres[i]);
      if ((i % 10) == 9)
	Log ("\n");
    }
  Log ("\nNo of photons and their fates\n");
  for (i = 0; i < NSTAT; i++)
    {
      Log ("%6d", nstat[i]);
      if ((i % 10) == 9)
	Log ("\n");
    }
  Log ("\n");


  /* Now subsample the matrix which contains the positions where photons
     last scattered and print the results to the screen.  Photons which are never
     in the grid end up in cell 0,0 

     In this instance subsampling means to actually do sums, such that all of
     the photons that scattered in the wind are counted in a courser grid.

   */




  Log ("Photons contributing to the various spectra\n");
  Log ("Inwin  Scat  Esc   Star  >nscat  err  Absor  Disk  sec  Adiab(matom)\n");
  for (n = 0; n < nspectra; n++)
    {
      for (i = 0; i < NSTAT; i++)
	Log (" %5d", xxspec[n].nphot[i]);
      Log ("\n");

    }


  return (0);

}

/***********************************************************
                                       Space Telescope Science Institute

 Synopsis:
   
	spectrum_summary(filename,mode,nspecmin,nspecmax,select_spectype,renorm,loglin)  
		writes out the spectrum to a file 

Arguments:		

	char filename[]		The name of the file to write
	char mode[];		The mode in which the file should be opened, usually 'w', but
					it could be 'a'
	int nspecmin,nspecmax	These two numbers define the spectra you want to write.    
	int select_spectype	The type of spectral file you want to create, 
					0 = raw, 1= flambda,2=fnu 
        double	renorm		This is renormalization which incrementally decreases to
				one as the detailed spectral calculation goes forward.  It
				was added to allow one to print out the spectrum at the
				end of each cycle, rather than the end of the entire
				calculation.
	char loglin[]		Are we outputting a log or a linear spectrum
				

Returns:
  
Description:	
	This simple routine simply writes the spectra to a file in an easily interpretable
	ascii format. Normally one would write all of the spectra in one go, but  one can use 
	spectrum summary to write various spectra to various files by using the variables
	nspecmin and nspecmax.. 
	
	Normally s[0],s[1],and s[2] will be the escaping, scattered, and absorbed spectrum.  
	The rest will be those which have been "extracted".   
	
	It can be called multiple times. In Python, it is currently called twice, once at the
	end of the ionization stage and once when computation of the detailed spectra are
	completed.
			
Notes:

History:
 	97jan      ksl	Coded and debugged as part of Python effort.
 	97sep13	ksl	Removed some extra code which had to do with opening the files.
 	97sep21	ksl	Modified normalization of spectrum s[n].renorm which affects 
 				Live or die option.  
	02apr	ksl	Added renorm option so that the spectrum will have the
			same overall "flux" when each incremental spectrum is printed
			out.
	10nov   nsh	Added another switch if we are outputting a log or a lin spectrum

**************************************************************/



int
spectrum_summary (filename, mode, nspecmin, nspecmax, select_spectype, renorm,
		  loglin)
     char filename[], mode[];
     int loglin;		// switch to tell the code if we are outputting a log or a lin
     int nspecmin, nspecmax;
     int select_spectype;
     double renorm;		// parameter used to rescale spectrum as it is building up 

{
  FILE *fopen (), *fptr;
  int i, n;
  char string[LINELENGTH];
  double freq, freqmin, dfreq, freq1;
  double lfreqmin, lfreqmax, ldfreq;
  double x, dd;



  /* Open or reopen a file for writing the spectrum */
  if ((fptr = fopen (filename, "w")) == NULL)
    {
      Error ("spectrum_summary: Unable to open %s for writing\n", filename);
      exit (0);
    }

  /* Check that nspecmin and nspecmax are reasonable */
  if (nspecmin < 0 || nspecmax < 0 || nspecmin > nspecmax)
    {
      Error
	("spectrum_summary: nspecmin %d or nspecmax %d not reasonable \n",
	 nspecmin, nspecmax);
      exit (0);
    }

  /* Construct and write a header string  for the output file */
  fprintf (fptr, "# Python Version %s\n", VERSION);
  fprintf (fptr, "# Git commit hash %s\n", GIT_COMMIT_HASH);

  get_time (string);
  fprintf (fptr, "# Date	%s\n#  \n", string);

  /* Save all of the parameter file information to the spectrum file */

  rdpar_save (fptr);


  /* Write the rest of the header for the spectrum file */
  /* JM 1411 -- Removed comment line for column headers, see #122 */
  fprintf (fptr, "# \nFreq.        Lambda");

  for (n = nspecmin; n <= nspecmax; n++)
    {
      fprintf (fptr, " %8s", xxspec[n].name);
    }


  fprintf (fptr, "\n");


  /* Ignore the end bins because they include all photons outside the frequency range and there may be some
     as a result of the fact that the bb function generate some IR photons */
  dd = 4. * PI * (100. * PC) * (100. * PC);

  if (loglin == 0)
    {
      freqmin = xxspec[nspecmin].freqmin;
      dfreq = (xxspec[nspecmin].freqmax - freqmin) / NWAVE;
      for (i = 1; i < NWAVE - 1; i++)
	{
	  freq = freqmin + i * dfreq;
	  fprintf (fptr, "%-8e %.3f ", freq, C * 1e8 / freq);
	  for (n = nspecmin; n <= nspecmax; n++)
	    {
	      x = xxspec[n].f[i] * xxspec[n].renorm;
	      if (select_spectype == 1)
		{		/* flambda */
		  x *= (freq * freq * 1e-8) / (dfreq * dd * C);
		}
	      else if (select_spectype == 2)
		{		/*fnu */
		  x /= (dfreq * dd);
		}
	      fprintf (fptr, " %8.3g", x * renorm);
	    }


	  fprintf (fptr, "\n");
	}
    }
  else if (loglin == 1)
    {
      lfreqmin = log10 (xxspec[nspecmin].freqmin);
      freq1 = lfreqmin;
      lfreqmax = log10 (xxspec[nspecmin].freqmax);
      ldfreq = (lfreqmax - lfreqmin) / NWAVE;

      for (i = 1; i < NWAVE - 1; i++)
	{
	  freq = pow (10., (lfreqmin + i * ldfreq));
	  dfreq = freq - freq1;
	  fprintf (fptr, "%-8e %.3f ", freq, C * 1e8 / freq);
	  for (n = nspecmin; n <= nspecmax; n++)
	    {
	      x = xxspec[n].lf[i] * xxspec[n].renorm;
	      if (select_spectype == 1)
		{		/* flambda */
		  x *= (freq * freq * 1e-8) / (dfreq * dd * C);
		}
	      else if (select_spectype == 2)
		{		/*fnu */
		  x /= (dfreq * dd);
		}
	      fprintf (fptr, " %8.3g", x * renorm);	/* this really shouldn't get called if we are outputting log data */
	    }


	  fprintf (fptr, "\n");
	  freq1 = freq;
	}
    }
  fclose (fptr);

  return (0);

}






/***********************************************************
                                       University of Southampton

 Synopsis:
   
	spectrum_spectrum_restart_renormalise() renormalises the spectra

Arguments:		
				
Returns:
  
Description:	
	This simple routine renormalises the spectra when one restarts
	a run which already has some spectral cycles run.

	This is done by multiplying all the spectra we have read in by a 
	factor

	renorm_factor = (cycles in old pf file) / (cycles in new pf file)
			
Notes:
	see Issue #134

History:
 	15jan JM	Coded

**************************************************************/

int spectrum_restart_renormalise(nangle)
    int nangle;
{
  double renorm_factor;
  int n, m, nspec;

  nspec = nangle + MSPEC;
  renorm_factor = ( (double) geo.pcycle + 1.0) / ( (double) geo.pcycles + 1.0); 

  /* loop over each spectrum column and each wavelength bin */
  for (n = MSPEC; n < nspec; n++)
  {
  	for (m = 0; m < NWAVE; m++)
  	  {
        xxspec[n].f[m] *= renorm_factor;
        xxspec[n].lf[m] *= renorm_factor;
  	  } 
  }

  return (0);
}
