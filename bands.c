
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "atomic.h"
#include "python.h"





/***********************************************************
                Space Telescope Science Institute

Synopsis:

Python uses a form of stratified sampling in an attempt to assure
that there are photon bundles at (high, generally) frequencies

This is the routine that initializes the bands.    There are
a number of possiblilities for setting up the bands


   
Arguments:		
     double t;			A temperature which can be used to set absolute limits on the bands
     double f1, f2;		frequency limits that can overide any other limits
     int imode;			A switch used for determining how the bands are to be populated
     struct xbands *band;	A poointer to the structure that holds the band information

     The currently allow modes are 

	imode	0	Use temperature to define a single band
		1	Use f1 and f2 to define a single band
		2	Use t,f1 and f2, and hardwired bands
			to define multiple bands which have been
			tuned to be relevent to CVs and other
			systems where a hot disk (T_eff of about
			100,000 K is assumed
		3	Bands tuned for yso
		4	Query the user to specify an arbitrary
			set of bands
		5      	Hardwired very wide bands for testing 
		6 	Bands set up by nsh to test cloudy   
		7	Bands hardwired for AGN paper1 by nsh
		8       Define bands in logarithmic intervals
Returns:

	The outputs are passed to other routines through the pointer
	to xbands.  The routine itself simply returns 0 on success
 
 
Description:	

		
Notes:

	10nov - ksl - Some of the choices have checks to see whether the bands
	that are being set up lie within f1 and f2.  Others do
	not.  This seems like an error.


History:
	01dec	ksl	python40
	02jul	ksl	Adapted from an similar routine used to set
			up bands in photon_gen.  This may eventually
			replace that routine
	0810	ksl	67 - Cleaned up slightly during relook
			at balance
	0812	ksl	67c - Changed names to bands_init to make
			it clear other routines created for diagnostic
			purposes had to be updated, as modified 
			to handle yso and cv type models in same
			code
	1112	ksl	71 - Moved material that was in python having
			to do with banding to this point.  Note
			This breaks the calls that exist in balance
	1208    nsh     73 - Several new bands
	1306	ksl	Minor modifications as tried to understand how
			we should handle banding more generally
**************************************************************/


/* Actual structures are in python.h.  Here for reference only.

#define NBANDS 10
struct xbands
{
	double f1[NBANDS],f2[NBANDS];
	double min_fraction[NBANDS];
	double nat_fraction[NBANDS];          // The fraction of the accepted luminosity in this band
	double used_fraction[NBANDS];
	double flux[NBANDS];                  //The "luminosity" within a band
	double weight[NBANDS];
	int nphot[NBANDS];
	int nbands;           // Actual number of bands in use
}
xband;


*/



int
bands_init (imode, band)
     int imode;			// A switch used for determining how the bands are to be populated
     struct xbands *band;

{
  int mode;
  int nband;
  double xx;
  double tmax, freqmin, freqmax;
  double t;			// A temperature which can be used to set absolute limits on the bands
  double f1, f2;		// frequency limits that can overide any other limits
  double fmax;
  double f1_log, f2_log, df;
  int ii;

  // 59 - Increased to 20,000 A so could go further into NIR 
  freqmin = C / 12000e-8;	/*20000 A */
  tmax = TSTAR;
  if (geo.twind > tmax)
    tmax = geo.twind;
  if (geo.tstar > tmax)
    tmax = geo.tstar;
  if (geo.t_bl > tmax && geo.lum_bl > 0.0)
    tmax = geo.t_bl;
  if ((0.488 * tdisk (geo.mstar, geo.disk_mdot, geo.rstar)) > tmax)
    tmax = 0.488 * tdisk (geo.mstar, geo.disk_mdot, geo.rstar);
  freqmax = BOLTZMANN * tmax / H * 10.;
  if (freqmax < 2.0 * 54.418 / HEV)
    {
      Log ("Increasing maximum frequency to twice the Helium edge\n");
      freqmax = 2.0 * 54.418 / HEV;
    }
  else
    Log ("Maximum frequency %8.2e determined by T %8.2e\n", freqmax, tmax);
  geo.tmax = tmax;		/*NSH 120817 NSH made this a global varaible so it is available to the code to make informed guesses as to the possible location of any BB driven exponential dropoff in the spectrum */
  t = tmax;
  f1 = freqmin;
  f2 = freqmax;

  /* end of import */

  if (imode == -1)
    {
      mode = 2;
      rdint
	("Photon.sampling.approach(0=T,1=(f1,f2),2=cv,3=yso,4=user_defined,5=cloudy_test,6=wide,7=AGN,8=logarithmic)",
	 &mode);
    }
  else
    {
      mode = imode;
    }

  if (mode == 0)
    {
  /* Mode 0 is sets a single band based on the temperature given */
      band->nbands = 1;
      band->f1[0] = BOLTZMANN * t / H * 0.05;
      band->f2[0] = BOLTZMANN * t / H * 20.;
      band->min_fraction[0] = 1.0;
    }
  else if (mode == 1)
    {
	    /* Mode 1 sets a single wide band defined by f1 and f2 */
      band->nbands = 1;
      band->f1[0] = f1;
      band->f2[0] = f2;
      band->min_fraction[0] = 1.0;
    }
  else if (mode == 2)		/* Traditional cv setup */
    {
      band->nbands = 4;
      band->f1[0] = f1;
      band->f2[0] = band->f1[1] = 13.599 / HEV;
      band->f2[1] = band->f1[2] = 24.588 / HEV;
      band->f2[2] = band->f1[3] = 54.418 / HEV;
      band->f2[3] = f2;
      band->min_fraction[0] = 0;
      band->min_fraction[1] = 0.1;
      band->min_fraction[2] = 0.1;
      band->min_fraction[3] = 0.1;
      if (f1 > band->f2[0])
	{
	  Error ("bands_init: f1 (%e) > 13.599/HEV)\n", f1);
	  exit (0);
	}
      if (f2 < band->f2[2])
	{
	  Error ("bands_init: f2 (%e) < 54.418/HEV)\n", f2);
	  exit (0);
	}

    }
  else if (mode == 3)		/* YSO setup */
    {
      band->nbands = 4;
      band->f1[0] = f1;
      band->f2[0] = band->f1[1] = 1.511 / HEV;
      band->f2[1] = band->f1[2] = 3.3998 / HEV;
      band->f2[2] = band->f1[3] = 6.0000 / HEV;
      band->f2[3] = f2;
      band->min_fraction[0] = 0;
      band->min_fraction[1] = 0.1;
      band->min_fraction[2] = 0.2;
      band->min_fraction[3] = 0.2;
      if (f1 > band->f2[0])
	{
	  Error ("bands_init: f1 (%e) > 13.599/HEV)\n", f1);
	  exit (0);
	}
      if (f2 < band->f2[2])
	{
	  Error ("bands_init: f2 (%e) < 54.418/HEV)\n", f2);
	  exit (0);
	}

    }
  else if (mode == 4)
    {
      rdint ("Num.of.frequency.bands", &band->nbands);
      Log ("Lowest photon energy is ev (freq) is %f (%.2e)\n", f1 * HEV,
	      f1);
      Log ("Highest photon energy is ev (freq) is %f (%.2e)\n", f2 * HEV,
	      f2);
      Log
	("Enter band boundaries in increasing eV, and assure they are between lowest and highest energy\n");


      rddoub ("Lowest_energy_to_be_considered(eV)", &xx);
      f1 = xx / HEV;

      rddoub ("Highest_energy_to_be_considered(eV)", &xx);
      f2 = xx / HEV;

      Log ("Lowest photon energy is ev (freq) is %f (%.2e)\n", f1 * HEV, f1);
      Log ("Highest photon energy is ev (freq) is %f (%.2e)\n", f2 * HEV, f2);


      band->f1[0] = f1;

      for (nband = 0; nband < band->nbands - 1; nband++)
	{
	  rddoub ("Band.boundary(eV)", &xx);
	  band->f2[nband] = band->f1[nband + 1] = xx / HEV;

	}
      band->f2[nband] = f2;

      Log
	("Enter mimimum fraction of photons in each band.  The total must be < or = to 1\n");

      for (nband = 0; nband < band->nbands; nband++)
	{
	  rddoub ("Band.minimum_fraction)", &band->min_fraction[nband]);
	}
      for (nband = 0; nband < band->nbands; nband++)
	{
	  Log ("For band %i, f1=%10.3e, f2=%10.3e, frac=%.2f\n", nband,
	       band->f1[nband], band->f2[nband], band->min_fraction[nband]);
	}



    }
  else if (mode == 5)		/* Set up to compare with cloudy power law table command note that this also sets up the weight and photon index for the PL, to ensure a continuous distribution */
    {
      rddoub ("Lowest_energy_to_be_considered(eV)", &xx);

      if (xx > geo.agn_cltab_low)
	{
	  xx = geo.agn_cltab_low / 10.0;
	  Log ("Lowest  frequency reset to 1/10 of low frequency break\n");
	}
      f1 = xx / HEV;
      rddoub ("Highest_energy_to_be_considered(eV)", &xx);

      if (xx < geo.agn_cltab_hi)
	{
	  xx = geo.agn_cltab_hi * 10.0;
	  Log ("highest  frequency reset to 10x high frequency break\n");
	}
      f2 = xx / HEV;
      Log ("Lowest photon energy is ev (freq) is %f (%.2e)\n", f1 * HEV, f1);
      Log ("Highest photon energy is ev (freq) is %f (%.2e)\n", f2 * HEV, f2);


      band->nbands = 12;






      band->f1[0] = (geo.agn_cltab_low / HEV) / 1000.0;
      band->f2[0] = band->f1[1] = (geo.agn_cltab_low / HEV) / 100.0;
      band->f2[1] = band->f1[2] = (geo.agn_cltab_low / HEV) / 10.0;
      band->f2[2] 	        = (geo.agn_cltab_low / HEV);


/* Now set up a set of log spaced bands in the range over the central range */
      f1_log = log10 (geo.agn_cltab_low/HEV);
      f2_log = log10 (geo.agn_cltab_hi/HEV);
      df = (f2_log - f1_log) / (6);

      ii = 3;
      while (ii < 9)
	{
	  band->f1[ii] = pow (10., f1_log + (ii-3) * df);
	  band->f2[ii] = pow (10., f1_log + ((ii-3) + 1) * df);
	  ii++;
	}

      f1_log = log10 (geo.agn_cltab_hi/HEV);
      f2_log = log10 (f2);
      df = (f2_log - f1_log) / (3);

      ii = 9;
      while (ii < 12)
	{
	  band->f1[ii] = pow (10., f1_log + (ii-9) * df);
	  band->f2[ii] = pow (10., f1_log + ((ii-9) + 1) * df);
	  ii++;
	}

 //     band->f1[9] = geo.agn_cltab_hi / HEV;
 //     band->f2[9] = f2;

      //Set number of photons in each band

      band->min_fraction[0] = 0.0666;
      band->min_fraction[1] = 0.0666;
      band->min_fraction[2] = 0.0666;

      band->min_fraction[3] = 0.1;
      band->min_fraction[4] = 0.1;
      band->min_fraction[5] = 0.1;
      band->min_fraction[6] = 0.1;
      band->min_fraction[7] = 0.1;
      band->min_fraction[8] = 0.1;

      band->min_fraction[9] =  0.0666;
      band->min_fraction[10] = 0.0666;
      band->min_fraction[11] = 0.0666;

      //Set alpha for each band

      band->alpha[0] = geo.agn_cltab_low_alpha;
      band->alpha[1] = geo.agn_cltab_low_alpha;
      band->alpha[2] = geo.agn_cltab_low_alpha;
      band->alpha[3] = geo.alpha_agn;
      band->alpha[4] = geo.alpha_agn;
      band->alpha[5] = geo.alpha_agn;
      band->alpha[6] = geo.alpha_agn;
      band->alpha[7] = geo.alpha_agn;
      band->alpha[8] = geo.alpha_agn;
      band->alpha[9] = geo.agn_cltab_hi_alpha;
      band->alpha[10] = geo.agn_cltab_hi_alpha;
      band->alpha[11] = geo.agn_cltab_hi_alpha;




      //Set the constant for each band to ensure continuous distribution

      band->pl_const[0] =
	geo.const_agn * pow ((band->f2[2]),
			     geo.alpha_agn) / pow ((band->f2[2]),
						   band->alpha[0]);
      band->pl_const[1] =
	geo.const_agn * pow ((band->f2[2]),
			     geo.alpha_agn) / pow ((band->f2[2]),
						   band->alpha[0]);
      band->pl_const[2] =
	geo.const_agn * pow ((band->f2[2]),
			     geo.alpha_agn) / pow ((band->f2[2]),
						   band->alpha[0]);
      band->pl_const[3] = geo.const_agn;
      band->pl_const[4] = geo.const_agn;
      band->pl_const[5] = geo.const_agn;
      band->pl_const[6] = geo.const_agn;
      band->pl_const[7] = geo.const_agn;
      band->pl_const[8] = geo.const_agn;

      band->pl_const[9] =
	geo.const_agn * pow ((band->f2[8]),
			     geo.alpha_agn) / pow ((band->f2[8]),
						   band->alpha[9]);

      band->pl_const[10] =
	geo.const_agn * pow ((band->f2[8]),
			     geo.alpha_agn) / pow ((band->f2[8]),
						   band->alpha[9]);

      band->pl_const[11] =
	geo.const_agn * pow ((band->f2[8]),
			     geo.alpha_agn) / pow ((band->f2[8]),
						   band->alpha[9]);


      for (nband = 0; nband < band->nbands; nband++)
	Log ("f1=%e,f2=%e,alpha=%e,const=%e,lum1=%e,lum2=%e\n",
		band->f1[nband], band->f2[nband], band->alpha[nband],
		band->pl_const[nband],
		band->pl_const[nband] * pow (band->f1[nband],
					     band->alpha[nband]),
		band->pl_const[nband] * pow (band->f2[nband],
					     band->alpha[nband]));




    }

  else if (mode == 6)		//Test for balance to have a really wide frequency range
    {
      tmax = geo.tstar;
      fmax = tmax * WIEN;	//Use wiens law to get peak frequency
      band->nbands = 17;


      band->f1[0] = 1e10;
      band->f2[0] = band->f1[1] = fmax * 0.01;
      band->f2[1] = band->f1[2] = fmax * 0.1;
      band->f2[2] = band->f1[3] = fmax;
      band->f2[3] = band->f1[4] = fmax * 1.5;
      band->f2[4] = band->f1[5] = fmax * 2;
      band->f2[5] = band->f1[6] = fmax * 2.5;
      band->f2[6] = band->f1[7] = fmax * 3;
      band->f2[7] = band->f1[8] = fmax * 4;
      band->f2[8] = band->f1[9] = fmax * 6;
      band->f2[9] = band->f1[10] = fmax * 8;
      band->f2[10] = band->f1[11] = fmax * 10;
      band->f2[11] = band->f1[12] = fmax * 12;
      band->f2[12] = band->f1[13] = fmax * 14;
      band->f2[13] = band->f1[14] = fmax * 16;
      band->f2[14] = band->f1[15] = fmax * 18;
      band->f2[15] = band->f1[16] = fmax * 20;
      band->f2[16] = 1e20;	

      band->min_fraction[0] = 0.1;
      band->min_fraction[1] = 0.1;
      band->min_fraction[2] = 0.1;
      band->min_fraction[3] = 0.05;
      band->min_fraction[4] = 0.05;
      band->min_fraction[5] = 0.05;
      band->min_fraction[6] = 0.05;
      band->min_fraction[7] = 0.05;
      band->min_fraction[8] = 0.05;
      band->min_fraction[9] = 0.05;
      band->min_fraction[10] = 0.05;
      band->min_fraction[11] = 0.05;
      band->min_fraction[12] = 0.05;
      band->min_fraction[13] = 0.05;
      band->min_fraction[14] = 0.05;
      band->min_fraction[15] = 0.05;
      band->min_fraction[16] = 0.05;
    }

  else if (mode == 7)		//Test for balance matching the bands we have been using for AGN runs
    {

      band->nbands = 10;
      band->f1[0] = 1e14;
      band->f2[0] = band->f1[1] = 1e15;
      band->f2[1] = band->f1[2] = 3.162e15;
      band->f2[2] = band->f1[3] = 1e16;
      band->f2[3] = band->f1[4] = 3.162e16;
      band->f2[4] = band->f1[5] = 1e17;
      band->f2[5] = band->f1[6] = 3.162e17;
      band->f2[6] = band->f1[7] = 1e18;
      band->f2[7] = band->f1[8] = 3.162e18;
      band->f2[8] = band->f1[9] = 1e19;
      band->f2[9] = 1e20;
      band->min_fraction[0] = 0.1;
      band->min_fraction[1] = 0.1;
      band->min_fraction[2] = 0.1;
      band->min_fraction[3] = 0.1;
      band->min_fraction[4] = 0.1;
      band->min_fraction[5] = 0.1;
      band->min_fraction[6] = 0.1;
      band->min_fraction[7] = 0.1;
      band->min_fraction[8] = 0.1;
      band->min_fraction[9] = 0.1;




    }
  else if (mode == 8)		/* 1306 - ksl - Generaglized method to set up logarithmic bands */
    {
      Log("Lowest photon energy is ev (freq) is %f (%.2e)\n", f1 * HEV,
	      f1);
      Log ("Highest photon energy is ev (freq) is %f (%.2e)\n", f2 * HEV,
	      f2);

      band->nbands = 5;
      rdint ("Num.of.frequency.bands", &band->nbands);

      if (band->nbands>NBANDS){
	      Error ("bands: Asking for more bands than allowed (%d) by size of bands.  Reducing to maximum value\n",NBANDS);
	      band->nbands=NBANDS;
      }

      xx = f1 * HEV;
      rddoub ("Lowest_energy_to_be_considered(eV)", &xx);
      f1 = xx / HEV;

      xx = f1 * HEV;
      rddoub ("Highest_energy_to_be_considered(eV)", &xx);
      f2 = xx / HEV;

      f1_log = log10 (f1);
      f2_log = log10 (f2);
      df = (f2_log - f1_log) / (band->nbands + 1);
      ii = 0;
      while (ii < band->nbands)
	{
	  band->f1[ii] = pow (10., f1_log + ii * df);
	  band->f2[ii] = pow (10., f1_log + (ii + 1) * df);
	  band->min_fraction[ii] = 1. / band->nbands;
	  ii++;
	}




    }
  else
    {
      Error ("bands_init: Unknown mode %d\n", mode);
      mytrap ();
    }


  Log ("bands_init: There are %d bands\n", band->nbands);
  for (nband = 0; nband < band->nbands; nband++)
    {
      Log ("bands_init: band %i,  f1=%10.3e,  f2=%10.3e, frac=%.2f\n", nband,
	   band->f1[nband], band->f2[nband], band->min_fraction[nband]);
      Log ("bands_init: band %i, eV1=%10.3e, eV2=%10.3e, frac=%.2f\n", nband,
	   band->f1[nband] * HEV, band->f2[nband] * HEV,
	   band->min_fraction[nband]);
      Log ("bands_init: band %i, alpha1=%f, alpha2=%f, frac=%.2f\n", nband,
	   band->f1[nband] * H / (BOLTZMANN * tmax),
	   band->f2[nband] * H / (BOLTZMANN * tmax),
	   band->min_fraction[nband]);
    }


  return (0);
}


/***********************************************************
                Space Telescope Science Institute

Synopsis:

	This is the routine where the frequency 
	boundaries for course spectra are established



   
Arguments:		

Returns:

 
 
Description:	

		
Notes:
	1112 - At present everything is hardwired



History:
	1112	ksl	Moved from main routine here
	111227	ksl	Smalle modifications to reflect my moving the main
			variables into the geo.structure so that they 
			could be read by py_oind
	111227	ksl	First attempt to limit the frequency intervals to
			regions where photons are being generated
**************************************************************/

int
freqs_init (freqmin, freqmax)
     double freqmin, freqmax;
{
  int i, n, ngood, good[NXBANDS];
  double xfreq[NXBANDS];
  int nxfreq;
//  double nupeak; //Weins law preak frequency from tstar


  /* At present set up a single energy band for 2 - 10 keV */
  /*NSH 70g - bands set up to match the bands we are currently using in the.pf files. This should probably end up tied together in the long run! */
/* nxfreq = 7;
  xfreq[0] = 1.0 / HEV;
 xfreq[1] = 13.6 / HEV;
xfreq[2] = 54.42 / HEV;
 xfreq[3] = 392. / HEV;
 xfreq[4] = 739. / HEV;
 xfreq[5] = 2000 / HEV;
 xfreq[6] = 10000 / HEV;
 xfreq[7] = 50000 / HEV;*/

/* bands to match the cloudy table spectrum - needed to cover all frequencies to let induced compton work OK */

if (geo.agn_ion_spectype == SPECTYPE_CL_TAB)
	{
	nxfreq = 3;
	xfreq[0]=0.0001/HEV;
	xfreq[1]=geo.agn_cltab_low/HEV;
	xfreq[2]=geo.agn_cltab_hi/HEV;
	xfreq[3]=100000000/HEV;
	}
	
/* bands try to deal with a blackbody spectrum */
/*nupeak=WIEN*geo.tstar;
printf ("Tstar=%e, Nupeak=%e\n",geo.tstar,nupeak);*/
else
	{
  nxfreq = 10;
  xfreq[0] = freqmin;		//We need the whole range to be modelled for induced compton heating to work
  xfreq[1] = 1e15;		//This should be below the lowest threshold frequency of any element in our model
  xfreq[2] = 3.162e15;
  xfreq[3] = 1e16;
  xfreq[4] = 3.162e16;
  xfreq[5] = 1e17;
  xfreq[6] = 3.162e17;
  xfreq[7] = 1e18;
  xfreq[8] = 3.162e18;
  xfreq[9] = 1.2e19;		//This is the highest frequency defined in our ionization data
  xfreq[10] = freqmax;
	}



  Log
    ("freqs_init: Photons will be generated between %8.2f (%8.2e) and %8.2f (%8.2e)\n",
     freqmin * HEV, freqmin, freqmax * HEV, freqmax);

  ngood = 0;
  for (i = 0; i < nxfreq; i++)
    {
//      Log ("freqs_init: %10.2e %10.2e %10.2e\n", freqmin, freqmax, xfreq[i]);
      if (freqmin < xfreq[i] && xfreq[i] < freqmax)
	{
	  good[i] = 1;
	  ngood++;
	}
      else if (freqmin < xfreq[i + 1] && xfreq[i + 1] < freqmax)
	{
	  good[i] = 1;
	  ngood++;
	}
      else
	{
	  good[i] = 0;
	}
    }

  Log ("freqs_init: Of %d starting intervals, %d will have photons\n", nxfreq,
       ngood);

  n = 0;
  for (i = 0; i < nxfreq; i++)
    {
      if (good[i] == 1)
	{
	  geo.xfreq[n] = xfreq[i];
	  geo.xfreq[n + 1] = xfreq[i + 1];
	  n++;
	}
    }
  geo.nxfreq = n;

  /* OK at this point we know at least some photons will be generated in each interval, but we still don't know
   * that the we are going to have a possibilty of photons throughout the first and last intervals.
   */

  if (freqmin > geo.xfreq[0])
    {
      geo.xfreq[0] = freqmin;
    }

  Log ("freqs_init: test %e %e\n", freqmax, geo.xfreq[geo.nxfreq]);
  if (freqmax < geo.xfreq[geo.nxfreq])
    {
      geo.xfreq[geo.nxfreq] = freqmax;
    }


  Log ("freqs_init: There were %d final intervals\n", geo.nxfreq);
  for (n = 0; n < geo.nxfreq; n++)
    {
      Log ("freqs_init: %8.2f (%8.2e)    %8.2f (%8.2e)  \n",
	   geo.xfreq[n] * HEV, geo.xfreq[n], geo.xfreq[n + 1] * HEV,
	   geo.xfreq[n + 1]);
    }


  return (0);

}
