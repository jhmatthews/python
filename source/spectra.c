

/***********************************************************/
/** @file  spectra.c
 * @author ksl
 * @date   May, 2018
 *
 * @brief
 *
 * The subroutines in this file handle allocating, incrementing, and writing out the
 * spectrum arrays for Python
 *
***********************************************************/


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "atomic.h"

#include "python.h"


int i_spec_start = 0;

/**********************************************************/
/**
 * @brief      allocates space for and initializes the spectrum arrays, storing
 * the criteria that are later used to create each spectrum.
 *
 * @param [in] double  f1   The minimum frequency for the spectra to be created
 * @param [in] double  f2   The maximum frequency in the spectra
 * @param [in] int  nangle   The number of different inclination angles (or more properly the number of spectra) to be created
 * @param [in] double  angle[]   The inclination angles associated with each spectrum
 * @param [in] double  phase[]   The orbital phase associated with each spectrum
 * @param [in] int  scat_select[]   A parameter for each spectrum which allows one to construct spectra with specifid numbers of scatters
 * @param [in] int  top_bot_select[]   A code which allows one to select photons only from below or above the disk
 * @param [in] int  select_extract   0 for Live or Die option, non-zero fo a normal extraction
 * @param [in] double  rho_select[]   Rho coordinates for extracting only photons in a particualr region
 * @param [in] double  z_select[]   Z cooordiante for selecting only photons that scattreed or were created in a praticular region
 * @param [in] double  az_select[]  Aximuthal angle for selecting only photons in a particular region
 * @param [in] double  r_select[]   Radius of the region from which to select photons
 * @return     Always returns 0
 *
 * @details
 * The first time spectrum_init  is called (i.e. if ispec_start=0), it allocates memory
 * for the various spectrum arrays.  (This is done one time, so one needs to allocate
 * space for all the arrays even though they are not all used in the ionization step).
 * The total number of spectra created is nangle+MSPEC.
 *
 * On subsequent calls to  spectrum_init, it rezeros all the spectrum information and
 * calculates the other information associated with each spectrum, such as the
 * angle cosines and the names of each spectrum.
 *
 * ### Notes ###
 * angle[],phase[] and scat_select[] only apply to the spectra extracted at
 * specific angles.  angle[0], phase[0], and scat_select[0] all affect spec[3]
 *
 * scat_select allows one to select spectra with a specific number of scattere or range og
 * scatters.  If nscat > 999 select all.  This is the normal case. The rest are used
 * largely for diagnostic purposes.  If 0 <nscat < MAXScat selected only photons which
 * have scattered nscat times.  If nscat is negattive, then all photons with more than |nscat| re
 * included.
 *
 * If top_bot_select is 0, then all photons are counted in the spectrum; if > 0, then only
 * photons above the disk plane are selected, if <0 only those from below the disk plane are selected.
 *
 * rho_select,z_select,az_select,r_select define a region of space and can be used
 * to create spectra only from photons that scattered or were created in these particular
 * regions.  It is normally only used for diagnostic reasons
 *
 *
 * Warning - Do not put anything in this routine that does anything but initialize
 * or reinitialize the spectrum structure s. This is important because this routine
 * is not accessed if one is continuing an old calculation of the detailed spectrum.
 * It is still used on a restart where the detailed spectral cycles have not begun
 * because in that case the spectra are not saved.
 *
 **********************************************************/
int
spectrum_init (f1, f2, nangle, angle, phase, scat_select, top_bot_select, select_extract, rho_select, z_select, az_select, r_select)
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
  double lfreqmin, lfreqmax, ldfreq;    /* NSH 1302 Min, max and delta for the log spectrum */
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
      Error ("spectrum_init: Could not allocate memory for %d spectra with %d wavelengths\n", nspec, NWAVE);
      Exit (0);
    }

    nspectra = nspec;           /* Note that nspectra is a global variable */

    i_spec_start = 1;           /* This is to prevent reallocation of the same arrays on multiple calls to spectrum_init */
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
      xxspec[n].lf[i] = 0;      /* NSH 1302 zero the logarithmic spectra */
    }
  }

  strcpy (xxspec[0].name, "Created");
  strcpy (xxspec[1].name, "Emitted");
  strcpy (xxspec[2].name, "CenSrc");
  strcpy (xxspec[3].name, "Disk");
  strcpy (xxspec[4].name, "Wind");
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
    xxspec[n].lmn[0] = sin (angle[n - MSPEC] / RADIAN) * cos (-phase[n - MSPEC] * 360. / RADIAN);
    xxspec[n].lmn[1] = sin (angle[n - MSPEC] / RADIAN) * sin (-phase[n - MSPEC] * 360. / RADIAN);
    xxspec[n].lmn[2] = cos (angle[n - MSPEC] / RADIAN);
    Log_silent ("Angle %e Angle cosines:%e %e %e\n", angle[n - MSPEC], xxspec[n].lmn[0], xxspec[n].lmn[1], xxspec[n].lmn[2]);

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

    strcpy (dummy, "");
    sprintf (dummy, "P%04.2f", phase[n - MSPEC]);
    strcat (xxspec[n].name, dummy);
    xxspec[n].nscat = scat_select[n - MSPEC];
    if (xxspec[n].nscat < MAXSCAT)
    {                           /* Then conditions have been place on the
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
    {                           /* Then conditions have been placed on the last
                                   location of the photon so update the names */
      strcpy (dummy, "");
      if (xxspec[n].top_bot == 1)
        sprintf (dummy, "_top");
      else if (xxspec[n].top_bot == -1)
        sprintf (dummy, "_bot");
      else if (xxspec[n].top_bot == 2)
      {
        sprintf (dummy, "_pos");
        xxspec[n].x[0] = rho_select[n - MSPEC] * cos (az_select[n - MSPEC] / RADIAN);
        xxspec[n].x[1] = rho_select[n - MSPEC] * sin (az_select[n - MSPEC] / RADIAN);
        xxspec[n].x[2] = z_select[n - MSPEC];
        xxspec[n].r = r_select[n - MSPEC];

      }
      else
      {
        Error ("spectrum_init: Unknown option %d\n", xxspec[n].top_bot);
        Exit (0);

      }
      strcat (xxspec[n].name, dummy);
    }
  }

  return (0);
}



/**********************************************************/
/**
 * @brief      Increments the spectrum arrays
 *  	after each flight of photons is processed (during ionization
 *  	cycles and for detailed spectra in the Live or Die option).
 *
 * @param [in] PhotPtr  p   A flight of photons
 * @param [in] double  f1   The minimum frequncy in the spectrum
 * @param [in] double  f2   The maximum frequency in the spectrum
 * @param [in] int  nangle  The number of different angles and phases for which to create detailed spectra
 * @param [in] int  select_extract   The integer stating whether the Live or Die option has
 * @return     Always returns 0
 *
 * @details
 * This routine increments the total spectrum arrays based on what has happened to each
 * photon during ionization cycles.  In the Live or Die option, the spectra at specific angles
 * are also created here when detailed spectra are created.
 *
 * The routine is called after each batch of photons has been transported through the wind and
 * prints some intermediate results to assure the user that the program is still running.
 *
 * ### Notes ###
 * Summing up of the spectra in the "extract" option is done in extract.c
 *
 *
 **********************************************************/

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
  int mscat, mtopbot;
  double delta;
  double nlow, nhigh;
  int k_orig, k1_orig;
  int iwind;                    // Variable defining whether this is a wind photon
  int max_scat, max_res;

  freqmin = f1;
  freqmax = f2;
  dfreq = (freqmax - freqmin) / NWAVE;
  nspec = nangle + MSPEC;
  nlow = 0.0;                   // variable to store the number of photons that have frequencies which are too low
  nhigh = 0.0;                  // variable to store the number of photons that have frequencies which are too high
  delta = 0.0;                  // fractional frequency error allowod


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

    /* Determine whether this is a wind photon, that is was it created in the
     * wind or scattered by the wind
     */

    iwind = 0;
    if (p[nphot].origin == PTYPE_WIND || p[nphot].origin == PTYPE_WIND_MATOM || p[nphot].nscat > 0)
    {
      iwind = 1;
    }

    /* find out where we are in log space */
    k1 = (log10 (p[nphot].freq) - log10 (freqmin)) / ldfreq;
    if (k1 < 0)
    {
      k1 = 0;
    }
    if (k1 > NWAVE - 1)
    {
      k1 = NWAVE - 1;
    }

    /* also need to work out where we are for photon's original wavelength */
    k1_orig = (log10 (p[nphot].freq_orig) - log10 (freqmin)) / ldfreq;
    if (k1_orig < 0)
    {
      k1_orig = 0;
    }
    if (k1_orig > NWAVE - 1)
    {
      k1_orig = NWAVE - 1;
    }


    /* lines to work out where we are in a normal spectrum with linear spacing */
    k = (p[nphot].freq - freqmin) / dfreq;
    if (k < 0)
    {
      if (((1. - p[nphot].freq / freqmin) > delta) && (geo.rt_mode != RT_MODE_MACRO))
        nlow = nlow + 1;
      k = 0;
    }
    else if (k > NWAVE - 1)
    {
      if (((1. - freqmax / p[nphot].freq) > delta) && (geo.rt_mode != RT_MODE_MACRO))
        nhigh = nhigh + 1;
      k = NWAVE - 1;
    }

    /* also need to work out where we are for photon's original wavelength */
    k_orig = (p[nphot].freq_orig - freqmin) / dfreq;
    if (k_orig < 0)
    {
      if (((1. - p[nphot].freq_orig / freqmin) > delta) && (geo.rt_mode != RT_MODE_MACRO))
        nlow = nlow + 1;
      k_orig = 0;
    }
    else if (k_orig > NWAVE - 1)
    {
      if (((1. - freqmax / p[nphot].freq_orig) > delta) && (geo.rt_mode != RT_MODE_MACRO))
        nhigh = nhigh + 1;
      k_orig = NWAVE - 1;
    }


    xxspec[0].f[k_orig] += p[nphot].w_orig;     /* created spectrum with original weights and wavelengths */
    xxspec[0].lf[k1_orig] += p[nphot].w_orig;   /* logarithmic created spectrum */
    if (iwind)
    {
      xxspec[0].f_wind[k_orig] += p[nphot].w_orig;
      xxspec[0].lf_wind[k1_orig] += p[nphot].w_orig;
    }


    if ((i = p[nphot].istat) == P_ESCAPE)
    {
      xxspec[0].nphot[i]++;
      xxspec[1].f[k] += p[nphot].w;     /* emitted spectrum */
      xxspec[1].lf[k1] += p[nphot].w;   /* logarithmic emitted spectrum */
      if (iwind)
      {
        xxspec[1].f_wind[k] += p[nphot].w;      /* emitted spectrum */
        xxspec[1].lf_wind[k1] += p[nphot].w;    /* logarithmic emitted spectrum */
      }
      xxspec[1].nphot[i]++;
      spectype = p[nphot].origin;

      /* When a photon that originated for example in the BL which has a type of PTYPE_BL is scattered in the wind by 
       * a macro atom it's type is increased by 10.  When we want to construct a spectrum for photons originating
       * from the boundary layer we need to subtract 10 from the type.    See python.h 
       */
      if (spectype >= 10)
        spectype -= 10;

      if (spectype == PTYPE_STAR || spectype == PTYPE_BL || spectype == PTYPE_AGN)      // Then it came from the bl or the star
      {
        xxspec[2].f[k] += p[nphot].w;   /* emitted star (+bl) spectrum */
        xxspec[2].lf[k1] += p[nphot].w; /* logarithmic emitted star (+bl) spectrum */
        if (iwind)
        {
          xxspec[2].f_wind[k] += p[nphot].w;    /* emitted spectrum */
          xxspec[2].lf_wind[k1] += p[nphot].w;  /* logarithmic emitted spectrum */
        }
        xxspec[2].nphot[i]++;
      }
      else if (spectype == PTYPE_DISK)  // Then it was a disk photon
      {
        xxspec[3].f[k] += p[nphot].w;   /* transmitted disk spectrum */
        xxspec[3].lf[k1] += p[nphot].w; /* logarithmic transmitted disk spectrum */
        if (iwind)
        {
          xxspec[3].f_wind[k] += p[nphot].w;    /* emitted spectrum */
          xxspec[3].lf_wind[k1] += p[nphot].w;  /* logarithmic emitted spectrum */
        }
        xxspec[3].nphot[i]++;
      }
      else if (spectype == PTYPE_WIND)
      {
        xxspec[4].f[k] += p[nphot].w;   /* wind spectrum */
        xxspec[4].lf[k1] += p[nphot].w; /* logarithmic wind spectrum */
        if (iwind)
        {
          xxspec[4].f_wind[k] += p[nphot].w;    /* emitted spectrum */
          xxspec[4].lf_wind[k1] += p[nphot].w;  /* logarithmic emitted spectrum */
        }
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
               (mscat < 0 && p[nphot].nscat >= (-mscat))) && ((mtopbot = xxspec[n].top_bot) == 0 || (mtopbot * p[nphot].x[2]) > 0))

          {
            if (xxspec[n].mmin < x1 && x1 < xxspec[n].mmax)
            {
              xxspec[n].f[k] += p[nphot].w;
              xxspec[n].lf[k1] += p[nphot].w;   /* logarithmic spectrum */
              if (iwind)
              {
                xxspec[n].f_wind[k] += p[nphot].w;      /* emitted spectrum */
                xxspec[n].lf_wind[k1] += p[nphot].w;    /* logarithmic emitted spectrum */
              }
            }
          }

        }
      }
    }
    else if (i == P_HIT_STAR || i == P_HIT_DISK)
    {
      xxspec[5].f[k] += p[nphot].w;     /*absorbed spectrum */
      xxspec[5].lf[k1] += p[nphot].w;   /*logarithmic absorbed spectrum */
      if (iwind)
      {
        xxspec[5].f_wind[k] += p[nphot].w;      /* emitted spectrum */
        xxspec[5].lf_wind[k1] += p[nphot].w;    /* logarithmic emitted spectrum */
      }
      xxspec[5].nphot[i]++;
    }

    if (p[nphot].nscat > 0 || p[nphot].nrscat > 0)

    {
      xxspec[6].f[k] += p[nphot].w;     /* j is the number of scatters so this constructs */
      xxspec[6].lf[k1] += p[nphot].w;   /* logarithmic j is the number of scatters so this constructs */
      if (iwind)
      {
        xxspec[6].f_wind[k] += p[nphot].w;      /* emitted spectrum */
        xxspec[6].lf_wind[k1] += p[nphot].w;    /* logarithmic emitted spectrum */
      }
      if (i < 0 || i > NSTAT - 1)
        xxspec[6].nphot[NSTAT - 1]++;
      else
        xxspec[6].nphot[i]++;   /* scattering spectrum */
    }

    if (i < 0 || i > NSTAT - 1)
      nstat[NSTAT - 1]++;
    else
      nstat[i]++;


  }


  if ((nlow / nphot > 0.05) || (nhigh / nphot > 0.05))
  {
    Error ("spectrum_create: Fraction of photons lost: %4.2f wi/ freq. low, %4.2f w/freq hi\n", nlow / nphot, nhigh / nphot);
  }
  else
  {
    Log ("spectrum_create: Fraction of photons lost:  %4.2f wi/ freq. low, %4.2f w/freq hi\n", nlow / nphot, nhigh / nphot);
  }



  max_scat = max_res = 0;

  for (i = 1; i < MAXSCAT; i++)
  {
    if (nscat[i] > 0)
    {
      max_scat = i;
    }
    if (nres[i] > max_res)
    {
      max_res = i;
    }
  }

  Log ("\nNo. of photons which have scattered n times.     The max number of scatters seen was %d\n", max_scat);

  for (i = 0; i <= max_scat; i++)
  {
    Log ("%6d", nscat[i]);
    if ((i % 10) == 9)
      Log ("\n");
  }

  Log ("\nNumber of photons resonantly scattering n times.  The max number of scatters seen was %d\n", max_res);
  for (i = 0; i <= max_res; i++)
  {
    Log ("%6d", nres[i]);
    if ((i % 10) == 9)
      Log ("\n");
  }
  Log ("\nNo of photons and their fates\n!!PhotFate: ");
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
  Log ("Inwind   Scat    Esc     Star    >nscat    err    Absorb   Disk    sec    Adiab(matom)\n");
  for (n = 0; n < nspectra; n++)
  {
    for (i = 0; i < NSTAT; i++)
      Log (" %7d", xxspec[n].nphot[i]);
    Log ("\n");

  }


  return (0);

}




/**********************************************************/
/**
 * @brief      writes out the spectra to a file
 *
 * @param [in] char  filename[]   The name of the file to write
 * @param [in] int  nspecmin   The number of the first spectrum to write
 * @param [in] int  nspecmax   The number of the last spectrum to write              .
 * @param [in] int  select_spectype   The type of spectral file you want to create
 * @param [in, out] double  renorm   This is renormalization which incrementally decreases to
 * one as the detailed spectral calculation goes forward.  It
 * was added to allow one to print out the spectrum at the
 * end of each cycle, rather than the end of the entire
 * calculation.
 * @param [in] int loglin 0 to print the spectrum out in linear units, 1 in log units
 * @param [in] int iwind If false (0), print out the normal spectrum; if true (1), print
 * out only photons that were scattered or created in the wind.
 *
 * @return     Always returns 0, unless there is a major problem in which case the program
 * exits
 *
 * @details
 * This simple routine simply writes the spectra to a file in an easily interpretable
 * ascii format.  The spectra will have already been created (using spectrum_create). 
 *
 * Normally one would write all of the spectra in one go, but  one can use
 * spectrum summary to write various spectra to various files by using the variables
 * nspecmin and nspecmax.
 *
 * Normally s[0],s[1],and s[2] will be the escaping, scattered, and absorbed spectrum.
 * The rest will be those which have been "extracted".
 *
 * It is called multiple times. In Python, it is currently called at two different
 * locations in the code, once at the
 * end of each ionization cycle  and at the end of each spectrum cycle.               
 * 
 * For the the spectrum cycles, the spectrum is 
 * and renormalized so that the overall normalization of the source
 * does not change.  This enables one to plot  the spectrum as the routine
 * is continuing to calculate the detatiled spectra, improving the statistics
 *
 * The rooutine is called separately to write the spectra out with a linear and
 * logarithmic binning.
 *
 * The file includes a copy of the .pf file inputs for the program so that one can
 * track exactly what inputs wre used to create the spectrum.
 *
 * ### Notes ###
 *
 **********************************************************/

int
spectrum_summary (filename, nspecmin, nspecmax, select_spectype, renorm, loglin, iwind)
     char filename[];
     int loglin;                // switch to tell the code if we are outputting a log or a lin
     int nspecmin, nspecmax;
     int select_spectype;
     double renorm;             // parameter used to rescale spectrum as it is building up
     int iwind;

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
    Exit (0);
  }

  /* Check that nspecmin and nspecmax are reasonable */
  if (nspecmin < 0 || nspecmax < 0 || nspecmin > nspecmax)
  {
    Error ("spectrum_summary: nspecmin %d or nspecmax %d not reasonable \n", nspecmin, nspecmax);
    Exit (0);
  }

  /* Construct and write a header string  for the output file */
  fprintf (fptr, "# Python Version %s\n", VERSION);
  fprintf (fptr, "# Git commit hash %s\n", GIT_COMMIT_HASH);

  get_time (string);
  fprintf (fptr, "# Date	%s\n#  \n", string);

  if (select_spectype == SPECTYPE_RAW)
  {
    fprintf (fptr, "\n# Units: L_nu spectrum (erg/s/Hz)\n\n");
  }
  else if (select_spectype == SPECTYPE_FLAMBDA)
  {
    fprintf (fptr, "\n# Units: flambda spectrum (erg/s/cm^-2/A) at %.1f parsecs\n\n", D_SOURCE);
  }
  else if (select_spectype == SPECTYPE_FNU)
  {
    fprintf (fptr, "\n# Units: Lnu spectrum (erg/s/Hz) at %.1f parsecs\n\n", D_SOURCE);
  }
  else
  {
    Error ("spectrum_summary: Unknown select_spectype %d\n", select_spectype);
    Exit (0);
  }


  /* Save all of the parameter file information to the spectrum file */

  rdpar_save (fptr);


  /* Write the rest of the header for the spectrum file */
  /* JM 1411 -- Removed comment line for column headers, see #122 */
  fprintf (fptr, "# \nFreq.        Lambda  ");

  for (n = nspecmin; n <= nspecmax; n++)
  {
    fprintf (fptr, " %8s", xxspec[n].name);
  }


  fprintf (fptr, "\n");


  /* Ignore the end bins because they include all photons outside the frequency range and there may be some
     as a result of the fact that the bb function generate some IR photons */
  dd = 4. * PI * (D_SOURCE * PC) * (D_SOURCE * PC);

  if (loglin == 0)              /* Then were are writing out the linear version of the spectra */
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
        if (iwind)
        {
          x = xxspec[n].f_wind[i] * xxspec[n].renorm;
        }


        if (select_spectype == SPECTYPE_FLAMBDA)
        {                       /* flambda */
          x *= (freq * freq * 1e-8) / (dfreq * dd * C);
        }
        else if (select_spectype == SPECTYPE_FNU)
        {                       /*fnu */
          x /= (dfreq * dd);
        }
        else if (select_spectype == SPECTYPE_RAW)
        {                       /*generated spectrum */
          x /= (dfreq);         //With log spectra implemented, we should divide by nu, so log and lin spectra agree
        }
        fprintf (fptr, " %10.5g", x * renorm);
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
        if (iwind)
        {
          x = xxspec[n].lf_wind[i] * xxspec[n].renorm;
        }

        if (select_spectype == SPECTYPE_FLAMBDA)
        {                       /* flambda */
          x *= (freq * freq * 1e-8) / (dfreq * dd * C);
        }
        else if (select_spectype == SPECTYPE_FNU)
        {                       /*fnu */
          x /= (dfreq * dd);
        }
        else if (select_spectype == SPECTYPE_RAW)
        {                       /*generated spectrum */
          x /= (dfreq);
        }
        fprintf (fptr, " %10.5g", x * renorm);
      }


      fprintf (fptr, "\n");
      freq1 = freq;
    }
  }
  fclose (fptr);

  return (0);

}








/**********************************************************/
/**
 * @brief      renormalizes the detailed spectra in case
 * of a restart
 *
 * @param [in] int  nangle   The number of discrete angles
 * @return     Always returns 0
 *
 * @details
 *
 * This routine deals with a very special case when one is
 * restarting a previous run, and adding spectral cycles to
 * obtain a detailed spectrum with higher statistics.
 *
 * In the previous run, the spectral will have been normalized
 * so that one gets the correct flux for the simulation.  We
 * are now computing additional spectral cycles, so we have
 * to reduce the fluxes for the spectral cycles that are read
 * in, so we can continue.
 * 	renorm_factor = (cycles in old pf file) / (cycles in new pf file)
 * 
 * Restarts can also occur when one is just completing a previous
 * run, without increasing the number of cycles, and in this case
 * no renormalization is required.
 * ### Notes ###
 * see Issue #134 and more importantly #503
 *
 **********************************************************/

int
spectrum_restart_renormalise (nangle)
     int nangle;
{
  double renorm_factor;
  int n, m, nspec;

  if (geo.pcycles == geo.pcycles_renorm)
  {
    return (0);
  }

  nspec = nangle + MSPEC;

  /* If we have gooten to this point, then the number of pcycles to which the
   * spectrum has been renormalized has changed, and so we must renormalize
   * the spectrum.  The original spectrum current is normalized to g
   * geo.pcycle/geo.pcycles_renorm of the final value.  We want the psectrum
   * to be geo.pcycles/geo.pcycles of the final value, so the renormlization
   * factor is geo.pcycles_renorm/geo.pcycles.
   */


  renorm_factor = ((double) geo.pcycles_renorm) / ((double) geo.pcycles);


  Log ("spectrum_restart_renormalise:  %d %d %d %f\n", geo.pcycle, geo.pcycles_renorm, geo.pcycles, renorm_factor);

  /* loop over each spectrum column and each wavelength bin */
  for (n = 0; n < nspec; n++)
  {
    for (m = 0; m < NWAVE; m++)
    {
      xxspec[n].f[m] *= renorm_factor;
      xxspec[n].lf[m] *= renorm_factor;
    }
  }

  return (0);
}
