


/**************************************************************************
                    Space Telescope Science Institute


  Synopsis:  These are routines needed to implement AGN into python

  Description:	

  Arguments:  


  Returns:

  Notes:


  History:
10oct	nsh	oded as part of initial effort to include a power law component
		to AGN

 ************************************************************************/

#include <stdio.h>
#include <strings.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "atomic.h"
#include "python.h"

#include "log.h"

/***********************************************************
Space Telescope Science Institute

Synopsis: star_init (r, tstar, freqmin, freqmax, ioniz_or_final, f)
 Arguments:              
 Returns:
 Description:    
 
This routine calculates the luminosity of the star and the luminosity within the frequency boundaries. 
BB functions are assumed 
Notes:
History:
**************************************************************/



/* This is essentially a parallel routine that is set up for other types of sources.  It actually does not do very
 * much
 */

double
agn_init (r, lum, alpha, freqmin, freqmax, ioniz_or_final, f)
     double r, lum, alpha, freqmin, freqmax, *f;
     int ioniz_or_final;
{

  double t;
  double emit, emittance_bb (), emittance_continuum ();
  int spectype;


  if (ioniz_or_final == 1)
    spectype = geo.agn_spectype;	/* type for final spectrum */
  else
    spectype = geo.agn_ion_spectype;	/*type for ionization calculation */
  if (spectype >= 0)
    {
      /* Assume that if we simulate the continuum for an AGN that emit is the luminosity in a specific range */
      emit = emittance_continuum (spectype, freqmin, freqmax, lum, alpha);
      *f = emit;
    }
  else if (spectype == SPECTYPE_POW)
    {
      /* Emittance_pow actucally returns the specific luminosity directly */
      emit = emittance_pow (freqmin, freqmax, lum, alpha);
      *f = emit;
    }
  else if (spectype == SPECTYPE_BB)
    {

      t = alpha;		/* For this case the second variable is t */
      /* Note that emittane_bb is really the energy radiated per unit area  so
       * this must be multiplied by the area of the object */
      emit = emittance_bb (freqmin, freqmax, t);
      *f = emit;
      *f *= (4. * PI * r * r);
    }
  else if (spectype == SPECTYPE_CL_TAB)
    {
      /* Emittance_pow actucally returns the specific luminosity directly */
      emit = emittance_bpow (freqmin, freqmax, lum, alpha);
      *f = emit;
    }

  return (*f);			/* Return the luminosity    */
}


/* This routine returns the specific luminosity, that is the luminosity within the frequency interval */
double
emittance_pow (freqmin, freqmax, lum, alpha)
     double freqmin, freqmax, lum, alpha;
{
  double constant, emit;
  /* these are the frequencies over which the power law is defined - currently set to the
     equivalent of 2 to 10 keV */


#define   XFREQMIN  4.84e17
#define   XFREQMAX  2.42e18

  /* first we need to calculate the constant for the power law function */
  constant = lum / (((pow (XFREQMAX, alpha + 1.)) - pow (XFREQMIN, alpha + 1.0)) / (alpha + 1.0));	/*NSH 1205 - this seems a bit unnecessary now. The constant is calculaed elsewhere, and with the broken power law we could probably get rid of this - is it even working properly now */

  /* now we need to work out the luminosity between our limited frequency range */
  /* we may need some checking routines to make sure that the requested frequency range is within the defined range,
     or it could default to zero outside the defined range */

  emit =
    constant * ((pow (freqmax, alpha + 1.0) - pow (freqmin, alpha + 1.0)) /
		(alpha + 1.0));

  return (emit);
}




/**************************************************************************
                    Southampton University


  Synopsis:  emittance_bpow computes the emittance of a broken power law coded to try to replicate the cloudy table type power law. There are three sections - a central bit, which should cover the 2-10kev portion of the spectrum, a low frequency tail and a high frequency tail . 

  Description:	

  Arguments:  freqmin and freqmax are the min and max frequencies that we are generating photons over, lum is the 2-10kev luminosity and alpha is the 2-10kev slope.


  Returns:   The emittance of the whole power law.

  Notes:


  History:
10oct	nsh	oded as part of initial effort to include a power law component
		to AGN

 ************************************************************************/

double
emittance_bpow (freqmin, freqmax, lum, alpha)
     double freqmin, freqmax, lum, alpha;
{
  double constant, constant_low, constant_hi, emit;
  double e1, e2, e3;
  double atemp, ctemp, f1, f2;
  double pl_low, pl_hi;
  /* these are the frequencies over which the power law is defined - currently set to the
     equivalent of 2 to 10 keV */


#define   XFREQMIN  4.84e17
#define   XFREQMAX  2.42e18
  f1 = freqmin;			/* NSH - 130506 added to reomve 03 compile errors */
  f2 = freqmax;			/* NSH - 130506 added to reomve 03 compile errors */
  e1 = e2 = e3 = 0.0;		/* NSH - 130506 added to reomve 03 compile errors */
  /* first we need to calculate the constant for the 2-10 kev power law function */
  constant =
    lum / (((pow (XFREQMAX, alpha + 1.)) - pow (XFREQMIN, alpha + 1.0)) /
	   (alpha + 1.0));


/* convert broken power law bands to freq */

  pl_low = geo.agn_cltab_low / HEV;
  pl_hi = geo.agn_cltab_hi / HEV;


  constant_low = xband.pl_const[0];	//we have already worked out the constants
  constant_hi = xband.pl_const[xband.nbands - 1];



  /* now we need to work out the luminosity between our limited frequency range */
  emit = 0.0;


  /* energy in low frequency tail */
  if (freqmin >= pl_low)
    {
      e1 = 0.0;			//Our frequencies start above the cutoff
      Log
	("Broken power law: Lowest frequency is above low frequency break\n");
    }
  else if (freqmax < pl_low)	//wierd situation, where all power is in the low frequency tail
    {
      atemp = geo.agn_cltab_low_alpha;
      f1 = freqmin;
      f2 = freqmax;
      ctemp = constant_low;
      e1 =
	ctemp * ((pow (f2, atemp + 1.0) - pow (f1, atemp + 1.0)) /
		 (atemp + 1.0));
      Log ("Broken power law: Emittance below low frequency break is %e\n",
	   e1);
    }
  else
    {
      atemp = geo.agn_cltab_low_alpha;
      f1 = freqmin;
      f2 = pl_low;
      ctemp = constant_low;
      e1 =
	ctemp * ((pow (f2, atemp + 1.0) - pow (f1, atemp + 1.0)) /
		 (atemp + 1.0));
      Log ("Broken power law: Emittance below low frequency break is %e\n",
	   e1);
    }

  /* energy in main frequency band */
  if (freqmax > pl_low && freqmin < pl_hi)	//If there is any part of the required range between our central range
    {
      if (freqmin <= pl_low && freqmax >= pl_hi)
	{
	  f1 = pl_low;
	  f2 = pl_hi;
	}
      else if (freqmin <= pl_low && freqmax < pl_hi)
	{
	  f1 = pl_low;
	  f2 = freqmax;
	}
      else if (freqmin > pl_low && freqmax < pl_hi)
	{
	  f1 = freqmin;
	  f2 = freqmax;
	}
      else if (freqmin > pl_low && freqmax >= pl_hi)
	{
	  f1 = freqmin;
	  f2 = pl_hi;
	}
      atemp = alpha;
      ctemp = constant;
      e2 =
	ctemp * ((pow (f2, atemp + 1.0) - pow (f1, atemp + 1.0)) /
		 (atemp + 1.0));
      Log ("Broken power law: Emittance in centre is %e\n", e2);
    }
  else
    {
      Log ("Broken power law: No luminosity in central range\n");
    }

  /* energy in high frequency tail */
  if (freqmax <= pl_hi)
    {
      e3 = 0.0;			//Our frequencies stop below the cutoff
      Log
	("Broken power law: Highest frequency is below high frequency break\n");
    }
  else if (freqmin > pl_hi)	//Odd situation where all power is in high frequency tail
    {
      atemp = geo.agn_cltab_hi_alpha;
      f1 = freqmin;
      f2 = freqmax;
      ctemp = constant_hi;
      e3 =
	ctemp * ((pow (f2, atemp + 1.0) - pow (f1, atemp + 1.0)) /
		 (atemp + 1.0));
      Log ("Broken power law: Emittance above high frequency break is %e\n",
	   e3);
    }
  else				//normal situation where we run from boundary to the max frequency
    {
      atemp = geo.agn_cltab_hi_alpha;
      f1 = pl_hi;
      f2 = freqmax;
      ctemp = constant_hi;
      e3 =
	ctemp * ((pow (f2, atemp + 1.0) - pow (f1, atemp + 1.0)) /
		 (atemp + 1.0));
      Log ("Broken power law: Emittance above high frequency break is %e\n",
	   e3);
    }

  emit = e1 + e2 + e3;

  return (emit);
}










int
photo_gen_agn (p, r, alpha, weight, f1, f2, spectype, istart, nphot)
     PhotPtr p;
     double r, alpha, weight;
     double f1, f2;		/* The freqency mininimum and maximum if a uniform distribution is selected */
     int spectype;		/*The spectrum type to generate: 0 is bb, 1 (or in fact anything but 0)
				   is uniform in frequency space */
     int istart, nphot;		/* Respecitively the starting point in p and the number of photons to generate */
{
  double freqmin, freqmax, dfreq, t;
  int i, iend;
  int n;
  double ftest;
  double dot ();
  double planck ();
  double plaw ();
  int randvec (), randvcos ();
  double zdisk ();

  t = alpha;			/* NSH 130605, this is a slightly odd statmenent, put in to avoid 03 compilation errors, but stemming from the odd behaviour that one can give the AGN a themal spectral type, and the alpha parameter is then used to transmit the temperature. */


  if ((iend = istart + nphot) > NPHOT)
    {
      Error ("photo_gen_agn: iend %d > NPHOT %d\n", iend, NPHOT);
      exit (0);
    }
  if (f2 < f1)
    {
      Error
	("photo_gen_agn: Cannot generate photons if freqmax %g < freqmin %g\n",
	 f2, f1);
    }
  Log_silent ("photo_gen_agn creates nphot %5d photons from %5d to %5d \n",
	      nphot, istart, iend);
  freqmin = f1;
  freqmax = f2;
  dfreq = (freqmax - freqmin) / MAXRAND;
  r = (1. + EPSILON) * r;	/* Generate photons just outside the photosphere unnecessary for the AGN perhaps? */

  if (spectype == SPECTYPE_CL_TAB)	/*If we have a broken PL, we need to work out which alpha to use */
    {
      ftest = (freqmin + freqmax) / 2.0;
      for (n = 0; n < xband.nbands; n++)	/* work out what alpha and constant to use */
	{
	  if (xband.f1[n] < ftest && xband.f2[n] > ftest)
	    {
	      alpha = xband.alpha[n];
	    }
	}
    }


  for (i = istart; i < iend; i++)
    {
      p[i].origin = PTYPE_AGN;	// For BL photons this is corrected in photon_gen 
      p[i].w = weight;
      p[i].istat = p[i].nscat = p[i].nrscat = 0;
      p[i].grid = 0;
      p[i].tau = 0.0;
      p[i].nres = -1;		// It's a continuum photon
      p[i].nnscat = 1;

      if (spectype == SPECTYPE_BB)
	{
	  p[i].freq = planck (t, freqmin, freqmax);
	}
      else if (spectype == SPECTYPE_UNIFORM)
	{			/* Kurucz spectrum */
	  /*Produce a uniform distribution of frequencies */
	  p[i].freq = freqmin + rand () * dfreq;
	}
      else if (spectype == SPECTYPE_POW)	/* this is the call to the powerlaw routine we are most interested in */
	{
	  p[i].freq = get_rand_pow (freqmin, freqmax, alpha);
	}
      else if (spectype == SPECTYPE_CL_TAB)
	{
	  p[i].freq = get_rand_pow (freqmin, freqmax, alpha);
	}
      else
	{
	  p[i].freq =
	    one_continuum (spectype, t, geo.gstar, freqmin, freqmax);
	}

      if (p[i].freq < freqmin || freqmax < p[i].freq)
	{
	  Error_silent
	    ("photo_gen_agn: phot no. %d freq %g out of range %g %g\n",
	     i, p[i].freq, freqmin, freqmax);
	}

      randvec (p[i].x, r);

      /* Added by SS August 2004 for finite disk. */
      if (geo.disk_type == 2)
	{
	  while (fabs (p[i].x[2]) < zdisk (r))
	    {
	      randvec (p[i].x, r);
	    }


	  if (fabs (p[i].x[2]) < zdisk (r))
	    {
	      Error ("Photon_agn: agn photon %d in disk %g %g %g %g %g\n",
		     i, p[i].x[0], p[i].x[1], p[i].x[2], zdisk (r), r);
	      exit (0);
	    }
	}


      /* this last bit is the direction, might need a change */
      randvcos (p[i].lmn, p[i].x);
   }
  return (0);
}
