/**************************************************************************
                    Space Telescope Science Institute
                                                                                                   
                                                                                                   
  Synopsis:

The routines in this file all have to do with photoionization and 
recombination rates and emissivities.  

(Note -- These descriptions need to be cleaned up and checked at some point. 
Today 06jul27 -- I just reorganized things more sensibly, but this section
should describe all of the routines generically)

Calculate the free-bound emissivity of a cell.   This is the emissivity,
eg ne * ni *enu.  It includes the electron density

  Description:

This is the version of the program that used detailed balance to
calculate the emissivity.  The specific formulation implemented
is that described most clearly in Hazy, but that comes from
Brown and Matthews 1970.

ion_choice      0- (nions-1) causes fb to return the emissivity of
		a specific ion number, i.e. 0 will be h1, 2 will be
		he1, etc.   
		>=nions returns the total emissivity
		<0  returns the emissivity of all metals, under the
		assumption that H and He comprise the first five
		ions in the input array.
fb_choice	When an ion recombines there are three possible quantities
		of interest.
			0=the emissivity/per unit frequency of the plasma
			1=the heat lost from the electrons (that is the
				fraction of the emissivity that does not
				go into the binding energy of the level
				to which the electron+ion are recombining)
			2=the emissivity in photons/unit frequency
		0 should be used in calculations of the emission spectrum but
		1 and should be used in energy loss and gain  calculations), while
	        2 should be used in ion densities and levels
                                                                                                   
  Arguments:  
                                                                                                   
                                                                                                   
  Returns:
                                                                                                   
  Notes:
                                                                                                   

For reference here is  the freebound structures that can be used for both the
specific emissivity of a free-bound transition.  But the real structue is
in python.h

	#define NTEMPS        100             // The number of temperatures which are stored in each fbstruct
	#define NFB   10              // The maximum number of frequency intervals for which the fb emission is calculated

	struct fbstruc
	{
	  double f1, f2;
	  double emiss[NIONS][NTEMPS];
	}
	freebound[NFB];

	double xnrecomb[NIONS][NTEMPS];       // There is only one set of recombination coefficients
	double fb_t[NTEMPS];

                                                                                                   
  History:
	01oct	ksl	Began work
	01nov	ksl	Adapt to include Verland cross-sections as well.
	02jun	ksl	Modified so could be used to calculate emissivities in
			photons as well as energy.
	06jul	ksl	57h-Standardized the headers and notes.  I have not
			really checked that everything is up to date, just
			put it into a more readable format.
	12jul	nsh	73-Subroutine bad_t_rr coded to generate a total
			recombination rate from badnell type parameters
	14jan	nsh	77a - Added some checks into the integrals, to ensure
			we do not attempt to integrate over a range of 
			frequencies so large that the integrand is zero over
			an excessive range - hence causing QROMB to return
			an answer of zero.
                                                                                                   
 ************************************************************************/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "atomic.h"
#include "python.h"
#include "recipes.h"



// Next line is left in to be sure that nfb is initialized properly at
// the beginning of the program

/* Next line is required for proper initialization */
int nfb = 0;			// Actual number of freqency intervals calculated

/* FBEMISS was calculated as follows:
x= 2. * PI * MELEC * BOLTZMANN / (H*H);
x=pow(x,-1.5);
x*=8. * PI / (C*C);
x*= H;
*/
#define FBEMISS   7.67413e-62	// Calculated with constants.c



/* These are external structures used primarily because we need to call 
Numerical Recipes routines from fb_verner and fb_topbase */

struct photoionization *fb_xver;	//Verner & Ferland description of a photoionization x-section
struct topbase_phot *fb_xtop;	//Topbase description of a photoionization x-section
double fbt;			// Temperature at which thee emissivity is calculated
int fbfr;			// fb_choice (see above)



/**************************************************************************
                    Space Telescope Science Institute
                                                                                                   
                                                                                                   
  Synopsis: fb_verner returns the partial (for a specific ion) emissivity for ions 
described in terms of Verner & Ferland photoionization x-sections.
                                                                                                   
  Description:
                                                                                                   
  Arguments:  
                                                                                                   
                                                                                                   
  Returns:
                                                                                                   
  Notes:
	The routines fb_verner_partial and fb_topbase_partial return the emissivity at
	a specific freqency.  Because they are called by some NR recipes 
	routines that integrates over frequency, most of the information 
	for these routines has to be and is passed by the external structures  
                                                                                                   
                                                                                                   
                                                                                                   
  History:
	02jul	ksl	Removed all references to the wind cell.
                                                                                                   
 ************************************************************************/

double
fb_verner_partial (freq)
     double freq;
{

  int nion;
  double partial;
  double x;
  double gn, gion;
  double fthresh;
  double sigma_phot ();

  fthresh = fb_xver->freq_t;
  if (freq < fthresh)
    return (0.0);		// No recombination at frequencies lower than the threshold freq occur

  nion = fb_xver->nion;
//?? Seems like gn should actually be the multiplicity of the excited and not the ground state ???
  gn = ion[nion].g;		// This is g factor of the ion to which you are recombining
  gion = ion[nion + 1].g;	// Want the g factor of the next ion up
  x = sigma_phot (fb_xver, freq);
// Next expression from Ferland
  partial =
    FBEMISS * gn / (2. * gion) * pow (freq * freq / fbt,
				      1.5) * exp (H_OVER_K *
						  (fthresh - freq) / fbt) * x;

// 0=emissivity, 1=heat loss from electrons, 2=photons emissivity
  if (fbfr == 1)
    partial *= (freq - fthresh) / freq;
  else if (fbfr == 2)
    partial /= (H * freq);

  return (partial);
}




/**************************************************************************
                    Space Telescope Science Institute
                                                                                                   
                                                                                                   
  Synopsis: fb_topbase_partial returns the partial (for a specific ion) emissivity for ions 
described in terms of Topbase photoionization x-sections.
                                                                                                   
  Description:
                                                                                                   
  Arguments:  
                                                                                                   
                                                                                                   
  Returns:
                                                                                                   
  Notes:
	The routines fb_verner_partial and fb_topbase_partial return the emissivity 
	at a specific freqency.  Because they are called by some NR recipes routines 
	that integrates over frequency, most of the information for these routines
	 has to be and is passed by the external structures. 
                                                                                                   
                                                                                                   
                                                                                                   
  History:
	02jul	ksl	Removed all references to the wind cell.
                                                                                                   
 ************************************************************************/


double
fb_topbase_partial (freq)
     double freq;
{
  int nion;
  double partial;
  double x;
  double gn, gion;
  double fthresh;
  double sigma_phot_topbase ();

  fthresh = fb_xtop->freq[0];
  if (freq < fthresh)
    return (0.0);		// No recombination at frequencies lower than the threshold freq occur

  nion = fb_xtop->nion;
  gn = config[fb_xtop->nlev].g;
  gion = ion[nion + 1].g;	// Want the g factor of the next ion up
  x = sigma_phot_topbase (fb_xtop, freq);
// Now calculate emission using Ferland's expression


  partial =
    FBEMISS * gn / (2. * gion) * pow (freq * freq / fbt,
				      1.5) * exp (H_OVER_K *
						  (fthresh - freq) / fbt) * x;



// 0=emissivity, 1=heat loss from electrons, 2=photons emissivity
  if (fbfr == 1)
    partial *= (freq - fthresh) / freq;
  else if (fbfr == 2)
    partial /= (H * freq);



  return (partial);
}

/**************************************************************************
                    Space Telescope Science Institute
                                                                                                   
                                                                                                   
  Synopsis: integ_fb calculates the integrated emissivity of the plasma, or the number of 
recombinations per second of a particular ion.  
                                                                                                   
  Description:
                                                                                                   
  Arguments:  
	t		The temperature at which the emissivity 
			or recombination rate is calculated
	f1,f2		The frequency limits on the calculation 
			of the emissivity (ignored if the number 
			of recombinations is desired.
	nion		The ion for which the emissivity is returned
	fb_choice	A switch which determines exactly what is to
			be returned: 
			0- the full emissivity including
			the energy associated associated with the
			threshold
			1- the (reduced) emissivity, e.g. excluding
			the threshold energy.  This is the energy
			associated with kinetic energy loss
			2- the specific recombination rate.
                                                                                                   
                                                                                                   
  Returns:
	The routine returns the specific emissivity, e.g. the emissivity and
	or recombination rate per electron and per ion.
                                                                                                   
  Notes:
	As written, in July02, the idea is that if the recombination coefficients
	have been calculated the program is going to return an interpolated
	coefficient, if not, it will calculate it from scratch.  (The later is 
	much slower if one has to do it a lot of times.

	???? Error -- There is definitely an error because the program does
	not support option 0, and this is needed for calculation of the
	relative numbers of photons by fb vs free free photons.  ksl 02jul
 
                                                                                                   
                                                                                                   
  History:
	02jul	ksl	Modified to elimate need to include information about the cell.
			At this point integ_fb yields answers per electron and per ion.
			Moved DENSITY_PHOT_MIN test out of integ_fb
	02jul	ksl	Original routines has been pushed down to xinteg_fb so
			that integ_fb can be modified to use stored values when
			desired.
                                                                                                   
 ************************************************************************/


double
integ_fb (t, f1, f2, nion, fb_choice)
     double t;			// The temperature at which to calculate the emissivity
     double f1, f2;		// The frequencies overwhich to integrate the emissivity
     int nion;			// The ion for which the "specific emissivity is calculateed
     int fb_choice;		// 0=full, otherwise reduced
{
  double xinteg_fb ();
  double fnu;
  double get_fb (), get_nrecomb ();
  int n;




  if (fb_choice == 1)
    {
      for (n = 0; n < nfb; n++)
	{
	  /* See if the frequencies correspond to one previously calculated */
	  if (f1 == freebound[n].f1 && f2 == freebound[n].f2)
	    {
	      fnu = get_fb (t, nion, n);
	      return (fnu);
	    }
	}
      /* If not calculate it here */
      fnu = xinteg_fb (t, f1, f2, nion, fb_choice);
      return (fnu);
    }
  else if (fb_choice == 2)
    {
      /* See if the frequencies correspond to one previously calculated */
      if (nfb > 0)
	{
	  fnu = get_nrecomb (t, nion);
	  return (fnu);
	}
      /* If not calculate it here */
      fnu = xinteg_fb (t, f1, f2, nion, fb_choice);
      return (fnu);
    }

  Error ("integ_fb: Unknown fb_choice(%d)\n", fb_choice);
  mytrap ();
  exit (0);
}




/**************************************************************************
                    Space Telescope Science Institute
                                                                                                   
                                                                                                   
  Synopsis: total_fb returns the energy lost from the plasma due to fb emission in a
	single wind cell at a temperture t between the frequncy limits f1 and f2.  
	The energy lost is just the kinetic energy lost from the plasma 
	because it does not include the ionization potential
	associated with each recombination.  Python tracks effectively the kinetic energy
	of the plasma (not the potential energy available if everything recombined. 
                                                                                                   
  Description:
                                                                                                   
  Arguments:  
                                                                                                   
                                                                                                   
  Returns:
                                                                                                   
  Notes:
                                                                                                   
                                                                                                   
  History:
  02jul	ksl	Modified so computes fb contributions of individual ions, instead
 		of just h, he, and z, and so that one does not need to pass
 		the entire wind cell to integ_fb 
  02jul	ksl	Added call to init_freebound so that precalculated values would
  		be used where possible.  This seemed the most conservative place
		to put these calls since total_fb is always called whenever the
		band-limited luminosities are needed.
                                                                                                   
 ************************************************************************/


double
total_fb (one, t, f1, f2)
     WindPtr one;
     double t, f1, f2;
{
  double total;
  int nion;
  int nplasma;
  PlasmaPtr xplasma;

  nplasma = one->nplasma;
  xplasma = &plasmamain[nplasma];


  if (t < 1000. || f2 < f1)
    return (0);			/* It's too cold to emit */

// Initialize the free_bound structures if that is necessary
  init_freebound (1.e3, 1.e9, f1, f2); //NSH 140121 increased limit to take account of hot plasmas


// Calculate the number of recombinations whenever calculating the fb_luminosities
  num_recomb (xplasma, t);

  total = 0;
  xplasma->lum_z = 0.0;
  for (nion = 0; nion < nions; nion++)
    {
      if (xplasma->density[nion] > DENSITY_PHOT_MIN)
	{

	  total += xplasma->lum_ion[nion] =
	    one->vol * xplasma->ne * xplasma->density[nion + 1] * integ_fb (t,
									    f1,
									    f2,
									    nion,
									    1);
	  if (nion > 3)
	    xplasma->lum_z += xplasma->lum_ion[nion];
	}

    }
  return (total);
}



/**************************************************************************
                    Space Telescope Science Institute
                                                                                                   
                                                                                                   
  Synopsis: one_fb generates one free bound photon with specific frequency limits
                                                                                                   
  Description:
                                                                                                   
  Arguments:  
	one 	The wind cell in which the photon is being 
		generated
	f1,f2	The frequency limits
                                                                                                   
                                                                                                   
  Returns:
	The frequency of the fb photon that was generated.
                                                                                                   
  Notes:
	57h -- This routine was a major time sync in python57g.  Most of the problem
	was ascociated with the generation of pdfs.  So the program now has a 
	rather multilayered approach to reducing the number of pdf generation
	steps that has to be done.  For a long time, we simply stored a pdf in
 	the pdf_fb array.  For python_57h, I created a new structure, photstoremain
	which parallels plasmamain and with allows one to store photons for future
	use.  
	
	It's possible the time for this routine to take could be reduced by storing
	more than one set of extra photons for different frequncy intervals, since
	our "banding" approach causes the frequencies to shift during the photon
	generation cycle.  Alternatively, it is possible that changing the order
	of photon generation could help.  But the approach adopted eliminate
	most of the problems with the routine, and so I have not pursued that.

	060802 -- ksl
                                                                                                   
                                                                                                   
  History:
	98	ksl	Coded as part of python effort
	98oct	ksl	Removed upper limits on freqency to attempt to resolve problems with
			different frequency limits in total_fb and one_fb
	01oct	ksl	Completely rewritten for topbase x-sections
	06jul	ksl	57h -- Modified to speed the program up by checking whether
			there are free_bound frequencies that have been previsulsy
			calculated for this cell and this condition.  Also, added
			a section which creates and stores multiple photons for
			the same conditions.  This reduces very significantly
			the number of times one has to construct a pdf, which is
			the main time sink for the program
                                                                                                   
 ************************************************************************/


double fb_x[200], fb_y[200];
double fb_jumps[NLEVELS];	// There is at most one jump per level
int fb_njumps = (-1);

WindPtr ww_fb;
struct Pdf pdf_fb;
double one_fb_f1, one_fb_f2, one_fb_te;	/* Old values */

double
one_fb (one, f1, f2)
     WindPtr one;		/* a single cell */
     double f1, f2;		/* freqmin and freqmax */
{
  double freq, tt, delta;
  int n;
  double fthresh, dfreq;
  int nplasma;
  PlasmaPtr xplasma;
  PhotStorePtr xphot;

  nplasma = one->nplasma;
  xplasma = &plasmamain[nplasma];
  xphot = &photstoremain[nplasma];

  if (f2 < f1)
    {
      Error ("one_fb: f2 %g < f1 %g Something is rotten  t %g\n",
	     f2, f1, xplasma->t_e);
      exit (0);
    }

/* Check if an apprpriate photon frequency has already been generated, and 
use that instead if possible --  57h */
  tt = xplasma->t_e;
  if (xphot->n < NSTORE && xphot->f1 == f1 && xphot->f2 == f2
      && xphot->t == tt)
    {
      freq = xphot->freq[xphot->n];
//TEST      Log("one_fb:  Using precalculated fb  %d \n",xphot->n);
      (xphot->n)++;
      return (freq);
    }
//TEST  else {
//TEST  Log("one_fb %3d %3d f %8.2e %8.2e %8.2e %8.2e t %6.1f %6.1f\n",nplasma,xphot->n,f1,f2,xphot->f1,xphot->f2,tt,xphot->t);
//TEST
//TEST}


  delta = 500;			// Fudge factor to prevent generation a photon if t has changed only slightly
  /* Check to see if we have already generated a pdf */
  if (tt > (one_fb_te + delta) || tt < (one_fb_te - delta) ||
      f1 != one_fb_f1 || f2 != one_fb_f2)
    {

/* Then need to generate a new pdf */

      ww_fb = one;

      /* Create the fb_array */

      /* Determine how many intervals are between f1 and f2.  These need to be
         put in increasing frequency order */

      if (f1 != one_fb_f1 || f2 != one_fb_f2)
	{			// Regenerate the jumps 
	  fb_njumps = 0;
	  for (n = 0; n < ntop_phot; n++)
	    {			//IS THIS ADDED BRACKET CORRECT? (SS, MAY04)
	      fthresh = phot_top_ptr[n]->freq[0];
	      if (f1 < fthresh && fthresh < f2)
		{
		  fb_jumps[fb_njumps] = fthresh;
		  fb_njumps++;
		}
	    }			//IS THIS CORRECT? (SS, MAY04)
	}

      //!BUG SSMay04
      //It doesn't seem to work unless this is zero? (SS May04)
      fb_njumps = 0;		// FUDGE (SS, May04)

      /* Note -- Need to fix this to get jumps properly, that is the
         frequencies need to allow for the jumps !! ??? */

      dfreq = (f2 - f1) / 199;
      for (n = 0; n < 200; n++)
	{
	  fb_x[n] = f1 + dfreq * n;
	  fb_y[n] = fb (xplasma, xplasma->t_e, fb_x[n], nions, 0);
	}

      if (pdf_gen_from_array
	  (&pdf_fb, fb_x, fb_y, 200, f1, f2, fb_njumps, fb_jumps) != 0)
	{
	  Error ("one_fb after error: f1 %g f2 %g te %g ne %g nh %g vol %g\n",
		 f1, f2, xplasma->t_e, xplasma->ne, xplasma->density[1],
		 one->vol);
	  Error ("Giving up");
	  exit (0);
	}
      one_fb_te = xplasma->t_e;
      one_fb_f1 = f1;
      one_fb_f2 = f2;		/* Note that this may not be the best way to check for a previous pdf */
    }

/* OK, we have not created a new pdf, cdf actually.  We are in a position to
generate photons */

/* First generate the phton we need */
  freq = pdf_get_rand (&pdf_fb);

/* Now create and store for future use a set of additonal photons */

  for (n = 0; n < NSTORE; n++)
    {
      xphot->freq[n] = pdf_get_rand (&pdf_fb);

    }
  xphot->n = 0;
  xphot->t = tt;
  xphot->f1 = f1;
  xphot->f2 = f2;
  return (freq);
}





/**************************************************************************
                    Space Telescope Science Institute
                                                                                                   
                                                                                                   
  Synopsis: num_recomb calculates the total number of recombinations in (units of #/cm**2/s) 
   in the cell per second for the all ions.    
                                                                                                   
  Description:
                                                                                                   
  Arguments:  
                                                                                                   
                                                                                                   
  Returns:
                                                                                                   
  Notes:
   The calculation is made purely direct recombination using the photoionization
   x-sections.  
   
   ?? It is possible that there should be a stimulated emission correction.
   ?? It is possible that a check to see that nfb should have been made
                                                                                                   
                                                                                                   
  History:
	98jun	ck      Error in one rate corrected
	02jun	ksl	This is a completely new version of the routine and
   			uses detailed balance. (python_43.5)
	02jul	ksl	Modified so that integ_fb is the number of 
			recombinations per ne and per ion
	06may	ksl	57+ -- Modified to use plasma structure since on volume
                                                                                                   
 ************************************************************************/

int
num_recomb (xplasma, t_e)
     PlasmaPtr xplasma;
     double t_e;
{
  int nelem;
  int i, imin, imax;
  for (nelem = 0; nelem < nelements; nelem++)
    {
      imin = ele[nelem].firstion;
      imax = imin + ele[nelem].nions;
      for (i = imin; i < imax; i++)
	{
	  if (xplasma->density[i] > DENSITY_PHOT_MIN)
	    {
	      xplasma->recomb[i] =
		xplasma->ne * xplasma->density[i + 1] * integ_fb (t_e, 3e14,
								  3e17, i, 2);
	    }
	}
      xplasma->recomb[imax] = 0.0;	// Can't recombine to highest i-state

    }

  return (0);
}


/**************************************************************************
                    Space Telescope Science Institute
                                                                                                   
                                                                                                   
  Synopsis: fb calculates the free_bound emissivity of the plasma at a specific frequency
                                                                                                   
  Description:
                                                                                                   
  Arguments:  

	ion_choice	Either the total emissivity or the emissivity for a specific
			ion is caculated depending on whether ion_choice=nions, or 
			a value less than the total number of ions
                                                                                                   
                                                                                                   
  Returns:
                                                                                                   
  Notes:
                                                                                                   
                                                                                                   
  History:
	02jul	ksl	Modified to reflect desire to make fb_xtopbase and fb_verner
			independent of the wind cell.
	06may	ksl	57+ -- Switched to plasma structure since no volume
	06jul	ksl	57h -- Cleaned this routine up a bit, in part to avoid
			calling fb_verner_partial when it should not be called.
                                                                                                   
 ************************************************************************/


double
fb (xplasma, t, freq, ion_choice, fb_choice)
     PlasmaPtr xplasma;		// A cell with all its associated density data
     double t;			// The temperature at which to calculate the emissivity
     double freq;		// The frequency at which to calculate the emissivity
     int ion_choice;		// Selects which ions the emissivity is to be calculated for (see above)
     int fb_choice;		// 0=full, otherwise reduced
{
  int n;
  double fnu, x;
  int ntmin, ntmax;		// These are the Topbase photo-ionization levels that are used
  int nv;			// The Verner photo-ionization level
  int nion, nion_min, nion_max;



  if (ion_choice < nions)	//Get emissivity for this specific ion_number
    {
      nion_min = ion_choice;
      nion_max = ion_choice + 1;
    }
  else if (ion_choice == nions)	// Get the total emissivity
    {
      nion_min = 0;
      nion_max = nions;
    }
  else
    {
      Error ("fb: This choice %d for ion_choice is not supported\n",
	     ion_choice);
      exit (0);
    }


  fbt = t;			/* Externally transmitted variable */
  fbfr = fb_choice;		/* Externally transmitted variable */

  fnu = 0.0;			/* Initially set the emissivity to zero */

  for (nion = nion_min; nion < nion_max; nion++)
    {
      ntmin = ion[nion].ntop_first;
      ntmax = ntmin + ion[nion].ntop;
      nv = ion[nion].nxphot;
      x = 0.0;

/* Loop over relevent Topbase photoionization x-sections.  If 
an ion does not have Topbase photoionization x-sections then
ntmin and ntmax are the same and the loop will be skipped. */
      for (n = ntmin; n < ntmax; n++)
	{
	  fb_xtop = &phot_top[n];	/*Externally transmited to fb_topbase_partial */
	  /* We don't want to include fb transitions associated with macro atoms here
	     - they are separated out for now. (SS, Apr 04). "If" statement added. */
	  if (fb_xtop->macro_info == 0 || geo.macro_simple == 1
	      || geo.rt_mode == 1)
	    {
	      x += fb_topbase_partial (freq);
	    }
	}
// Loop over the relevant Verner x-sections   xxxx

/* 57h -- Fixed this so that only call fb_verner_partial when necessary */
      if (nv > -1)
	{
	  fb_xver = &xphot[nv];	// This assigns a specific verner x-section to fb_xver 
	  fnu += fb_verner_partial (freq);
	}
      fnu += xplasma->density[nion] * x;
    }
  fnu *= xplasma->ne;		// Correct from specific emissivity to the total fb emissivity

  return (fnu);

}


/**************************************************************************
                    Space Telescope Science Institute
                                                                                                   
                                                                                                   
  Synopsis: init_freebound initializes the structure fb_struc as well as some
associated arrays and variables (found in python.h) that describe
recombination rates and band-limited luminosities.
                                                                                                   
  Description:
                                                                                                   
  Arguments:  
                                                                                                   
	t1, t2	The lower and upper limits for which the fb 
		information is calculated.
	f1, f2	The frequency interval in which the band-limited
		fb information is calculated.
                                                                                                   
  Returns:
                                                                                                   
  Notes:
	The first time the routine is called, both recombination
	rates and band-limited luminosities are calculated.  On
	subsequent calls the routine checks to see whether it has
	already calculated the band-limited freebound emissivities, 
	and if so returns without redoing the calculation.  However, 
	if a new frequency interval is provided, the new luminosities
	are added to the free-bound structure.  To force a 
	re-initialization nfb must be set to 0.
                                                                                                   
                                                                                                   
  History:
	02jul	ksl	Coding began
	0810	ksl	67 - Modified routine so that instead of exiting
			when there are more than NFB sets of data, it 
			creates a new set of data and assumes the oldest
			set can be discarded.  This was done primarily
			to accommodate some runs of balance.
                                                                                                   
 ************************************************************************/

int init_freebound_nfb;		/*Indicates the total number of freebound sets that
				   could be used */

int
init_freebound (t1, t2, f1, f2)
     double t1, t2, f1, f2;
{
  double t;
  int i, j, nion;
  double ltmin, ltmax, dlt;
  double xinteg_fb ();
  int nput;

/* If init-freebound has never been called before initialize
   the temperatures and calculate the
    recombination rates.  Otherwise skip this section */


/* First see if there is a recombination file that already exists that can be used*/

  if (nfb == 0)
    {
      fb_read ("recomb.save");
    }

  if (nfb == 0)
    {
      if (t2 < t1)
	{
	  Error ("init_freebound: t2(%g)<t1(%g)\n", t2, t1);
	  exit (0);
	}

      ltmin = log10 (t1);
      ltmax = log10 (t2);
      dlt = (ltmax - ltmin) / (NTEMPS - 1);

      for (j = 0; j < NTEMPS; j++)
	{
	  fb_t[j] = pow (10., ltmin + dlt * j);
	}

      Log ("init_freebound: Creating recombination coefficients\n");
      for (nion = 0; nion < nions; nion++)
	{
	  for (j = 0; j < NTEMPS; j++)
	    {
	      t = fb_t[j];

	      xnrecomb[nion][j] = xinteg_fb (t, 0.0, 1.e50, nion, 2);
	    }
	}
    }
  else if (fabs (fb_t[0] - t1) > 10. || fabs (fb_t[NTEMPS - 1] - t2) > 1000.)
    {
      Error
	("init_freebound: Cannot initialize to new temps without resetting nfb");
      exit (0);

    }

/* Now check to see whether the freebound information has already
been calculated for these conditions, and if so simply return.
*/
  i = 0;
  while ((freebound[i].f1 != f1 || freebound[i].f2 != f2) && i < nfb)
    i++;

  if (i < nfb)
    {
      return (0);
    }

/* We have to calculate a new set of freebound data */
  if (i == NFB)
    {
      /* We've filled all the available space in freebound so we start recycling elements, assuming that the latest
       * ones are still likelyt to be needed
       */
      nput = init_freebound_nfb % NFB;
      init_freebound_nfb++;

      Error
	("init_freebound: Recycling freebound, storage for NFB (%d), need %d to avoid \n",
	 NFB, init_freebound_nfb);

    }
  else
    {
      nput = init_freebound_nfb = nfb;
      nfb++;
    }



/* Having reach this point, a new set of fb emissivities
must be calculated.  Note that old information is not destroyed
unless nfb had been set to 0.  The new set is added to the old
on the assumption that the fb information will be reused.
*/

  freebound[nput].f1 = f1;
  freebound[nput].f2 = f2;

/* So at this point, we know exactly what needs to be calculated */

  Log
    ("init_freebound: Creating recombination emissivites between %e and %e\n",
     f1, f2);

  for (nion = 0; nion < nions; nion++)
    {
//      Log ("init_freebound:  ion %d\n", nion);
      for (j = 0; j < NTEMPS; j++)
	{			//j covers the temps
	  t = fb_t[j];
	  freebound[nput].emiss[nion][j] = xinteg_fb (t, f1, f2, nion, 1);
	}
    }


  // OK we are done
  return (0);
}




/**************************************************************************
                    Space Telescope Science Institute
                                                                                                   
                                                                                                   
  Synopsis: Return the recombination coefficient 
                                                                                                   
  Description:
                                                                                                   
  Arguments:  
                                                                                                   
                                                                                                   
  Returns:
                                                                                                   
  Notes:
                                                                                                   
                                                                                                   
  History:
		13sep	nsh	changed call to linterp to reflect new option
                                                                                                   
 ************************************************************************/


double
get_nrecomb (t, nion)
     double t;
     int nion;
{
  int linterp ();
  double x;

  linterp (t, fb_t, xnrecomb[nion], NTEMPS, &x,0); //Interpolate in linear space
  return (x);
}


/* Return the specific emissivity due to recombination emission in an interval */

double
get_fb (t, nion, narray)
     double t;
     int nion;
     int narray;
{
  int linterp ();
  double x;

  linterp (t, fb_t, &freebound[narray].emiss[nion][0], NTEMPS, &x,0); //Interpolate in linear space
  return (x);
}



/**************************************************************************
                    Space Telescope Science Institute
                                                                                                   
                                                                                                   
  Synopsis: xinteg_fb calculates the integrated emissivity of the plasma.  
                                                                                                   
  Description:
                                                                                                   
  Arguments:  
                                                                                                   
                                                                                                   
  Returns:
                                                                                                   
  Notes:
It's unusual nature is determined
by the need to use a modififed Numerical Recipes routine for integration of fb over
a frequency range 
                                                                                                   
                                                                                                   
  History:
	02jul	ksl	Modified to elimate need to include information about the cell.
			At this point integ_fb yields answers per electron and per ion.
			Moved DENSITY_PHOT_MIN test out of integ_fb
	02jul	ksl	In attempt to store the fb coefficiencts this function
			renamed.  It actually calculates the emissivities, and integ_fb
			has become a kind of steering routine that normally, e.g. often
			reads the freebound arrays.
        04May   SS      Minor changes to exclude fb from macro atoms (which are treated elsewhere)
                                                                                                   
 ************************************************************************/


double
xinteg_fb (t, f1, f2, nion, fb_choice)
     double t;			// The temperature at which to calculate the emissivity
     double f1, f2;		// The frequencies overwhich to integrate the emissivity
     int nion;			// The ion for which the "specific emissivity is calculateed
     int fb_choice;		// 0=full, otherwise reduced
{
  int n;
  double fnu;
  double dnu;  //NSH 140120 - a parameter to allow one to restrict the integration limits.
  double fthresh, fmax;
  double den_config ();
  double sigma_phot (), sigma_phot_topbase ();
  int ntmin, ntmax;		// These are the Topbase photo-ionization levels that are used
  int nvmin, nvmax;		// These are the limits on the Verland x-sections
  double qromb ();
  double fb_topbase_partial (), fb_verner_partial ();

  dnu=0.0; //Avoid compilation errors.


  if (-1 < nion && nion < nions)	//Get emissivity for this specific ion_number
    {
      ntmin = ion[nion].ntop_first;
      ntmax = ntmin + ion[nion].ntop;
      nvmin = nion;
      nvmax = nvmin + 1;
    }
  else				// Get the total emissivity
    {
      Error ("integ_fb: %d is unacceptable value of nion\n", nion);
      mytrap ();
      return (0);
    }

// Put information where it can be used by the integrating function
  fbt = t;
  fbfr = fb_choice;

/* Limit the frequency range to one that is reasonable before integrating */

  if (f1 < 3e12)
    f1 = 3e12;			// 10000 Angstroms
  if (f2 > 3e18)		// 110819 nsh increase upper limits to include  highly ionised ions that we are now seeing in x-ray illuminated nebulas.
    f2 = 3e18;			// This is 1 Angstroms  - ksl
  if (f2 < f1)
    return (0);			/* Because there is nothing to integrate */

  fnu = 0.0;

  for (n = ntmin; n < ntmax; n++)
    {				// loop over relevent Topbase photoionzation x-sections
      fb_xtop = &phot_top[n];
      /* Adding an if statement here so that photoionization that's part of a macro atom is 
         not included here (these will be dealt with elsewhere). (SS, Apr04) */
      if (fb_xtop->macro_info == 0 || geo.macro_simple == 1 || geo.rt_mode == 1)	//Macro atom check. (SS)
	{
	  fthresh = fb_xtop->freq[0];
	  fmax = fb_xtop->freq[fb_xtop->np - 1];	// Argues that this should be part of structure
	  if (f1 > fthresh)
	    fthresh = f1;
	  if (f2 < fmax)
	    fmax = f2;

	  // Now calculate the emissivity as long as fmax exceeds xthreshold and there are ions to recombine
	  if (fmax > fthresh)
                {
//NSH 140120 - this is a test to ensure that the exponential will not go to zero in the integrations 
     		dnu = 100.0 * (fbt / H_OVER_K);
      		if (fthresh + dnu < fmax)
			{
	  		fmax = fthresh + dnu;
			}
	    	fnu += qromb (fb_topbase_partial, fthresh, fmax, 1.e-5);
	    	}
	}
    }

/* This completes the calculation of those levels for which we have Topbase x-sections, now do Verner */

  for (n = nvmin; n < nvmax; n++)
    {
      if (ion[n].phot_info == 0)
	{			// Only work on ions without Topbase and with Verner
	  fb_xver = &xphot[ion[n].nxphot];
	  fthresh = fb_xver->freq_t;
	  fmax = fb_xver->freq_max;	//So at this point these are the maximal allowable
	  if (f1 > fthresh)
	    fthresh = f1;	//So move fthresh up, if f1 was greater than this
	  if (f2 < fmax)	//Move fmax down if we wanted a narrower range
	    fmax = f2;
	  // Now integrate only if its in allowable range  && there are ions to recombine
	  if (fmax > fthresh)
//NSH 140120 - this is a test to ensure that the exponential will not go to zero in the integrations 
		{     		
		dnu = 100.0 * (fbt / H_OVER_K);
      		if (fthresh + dnu < fmax)
			{
	  		fmax = fthresh + dnu;
			}
	        fnu += qromb (fb_verner_partial, fthresh, fmax, 1.e-5);
		}
	}
    }

  return (fnu);
}


/***********************************************************
                                       Space Telescope Science Institute
                                                                                                                                      
 Synopsis:
        fb_save(w,filename)
        fb_read(filename)
                                                                                                                                      
Arguments:
                                                                                                                                      
                                                                                                                                      
Returns:
                                                                                                                                      
Description:
                                                                                                                                      
        The two routines in this file write and read the structure.
	associated with fb emission in order to save a bit of time
	calculating them.  This is semidangerous so we need to save
	some additional information to see, whether this is a good
        idea
                                                                                                                                      
Notes:
	ksl 0810 - The freebound structure could be done in a way
	that saves space on disk and in the program.  To save
	space on disk one just needs to write and read the array
	one elment at a time.  
                                                                                                                                      
History:
        06jul   ksl     57h -- Began coding
	08may	ksl	60a -- Added code to make sure that one
			was only checking the first word in VERSION.
			Previously, there was a problem with VERSION
			having an extra space compared to the version
			that was read back.  May be just an OS X
			problem.
                                                                                                                                      
**************************************************************/

int
fb_save (filename)
     char filename[];
{
  FILE *fptr, *fopen ();

  char line[LINELENGTH];
  int n;

  if ((fptr = fopen (filename, "w")) == NULL)
    {
      Error ("fb_save: Unable to open %s\n", filename);
      exit (0);
    }

  sprintf (line, "Version %s\n", VERSION);
  n = fwrite (line, sizeof (line), 1, fptr);
  sprintf (line, "%s %d %d %d %d %d %d\n", geo.atomic_filename, nelements,
	   nions, nlevels, nxphot, ntop_phot, nfb);
  n += fwrite (line, sizeof (line), 1, fptr);

  n += fwrite (fb_t, sizeof (fb_t), 1, fptr);
  n += fwrite (xnrecomb, sizeof (xnrecomb), 1, fptr);
  n += fwrite (freebound, sizeof (freebound), 1, fptr);

  fclose (fptr);

  Log ("Wrote fb information for %s to %s for %d frequency intervals \n",
       geo.atomic_filename, filename, nfb);
  Log
    ("There are  %d elem, %d ions, %d levels, %d Verner, %d topbase xsections\n",
     nelements, nions, nlevels, nxphot, ntop_phot);

  return (0);

}


int
fb_read (filename)
     char filename[];
{
  FILE *fptr, *fopen ();
  int n;
  char line[LINELENGTH];
  char version[LINELENGTH], oldversion[LINELENGTH];
  char atomic_filename[LINELENGTH];

  int xnelements, xnions, xnlevels, xnxphot, xntopphot, xnfb;

/* Initialize nfb to 0 so a return means that python will
have to calcuate the coefficients */

  nfb = 0;

  if ((fptr = fopen (filename, "r")) == NULL)
    {
      Log ("fb_read: Unable to open %s\n", filename);
      return (0);
    }

  n = fread (line, sizeof (line), 1, fptr);
  sscanf (line, "%*s %s", version);
  Log
    ("Reading Windfile %s created with python version %s with python version %s\n",
     filename, version, VERSION);
  n += fread (line, sizeof (line), 1, fptr);
  sscanf (line, "%s %d %d %d %d %d %d\n", atomic_filename, &xnelements,
	  &xnions, &xnlevels, &xnxphot, &xntopphot, &xnfb);

  Log ("Reading fb information for %s \n", atomic_filename);
  Log
    ("There are  %d elem, %d ions, %d levels, %d Verner, %d topbase, and %d freq intervals\n",
     xnelements, xnions, xnlevels, xnxphot, xntopphot, xnfb);

/* Check insofar as possible that this reconbination file was created with the same 
inputs.  The next statement is intended to handle problems with extra spaces on the
string that constitutes VERSION*/

  sscanf (VERSION, "%s", oldversion);
  if (strcmp (version, oldversion))
    {
      Log ("fb_read: Different versions of python  %s != %s or %d !=d %d \n",
	   version, oldversion, strlen (version), strlen (oldversion));
      return (0);
    }
  if (strcmp (atomic_filename, geo.atomic_filename) != 0)
    {
      Log ("fb_read: Different atomic_filename\n");
      return (0);
    }
  if (xnelements != nelements)
    {
      Log ("fb_read: old %d  new %d  nelements\n", xnelements, nelements);
      return (0);
    }
  if (xnlevels != nlevels)
    {
      Log ("fb_read: old %d  new %d  nions\n", xnions, nions);
      return (0);
    }
  if (xnlevels != nlevels)
    {
      Log ("fb_read: old %d  new %d  nlevels\n", xnlevels, nlevels);
      return (0);
    }
  if (xnxphot != nxphot)
    {
      Log ("fb_read: old %d  new %d  nnxphot\n", xnxphot, nxphot);
      return (0);
    }
  if (xnlevels != nlevels)
    {
      Log ("fb_read: old %d  new %d  ntopphot\n", xntopphot, ntop_phot);

      return (0);
    }





  nfb = xnfb;

  n += fread (fb_t, sizeof (fb_t), 1, fptr);
  n += fread (xnrecomb, sizeof (xnrecomb), 1, fptr);
  n += fread (freebound, sizeof (freebound), 1, fptr);

  fclose (fptr);


  Log ("Read geometry and wind structures from windsavefile %s\n", filename);

  return (n);


}


/***********************************************************
                                       Southampton University
                                                                                                                                      
 Synopsis:
        bad_t_rr(nion, T)
                                                                                                                                      
Arguments:
        ion - ion for which we want a recombination rate - 
		this is the upper state, so the ion which is 
		doing the recombining, there is no rate for
		H1(ion0) but there is one for H(ion1)
	temperature - the temperature we want a rate for
Returns:
	rate - the total recombination rate for this ion 
		for this temperature
                                                                                                                                      
Description:
                                                                                                                                      
        This routine generates a total recombination rate for 
		a given ion at a given temperature using 
		badnell type parameters.
                                                                                                                                      
Notes:
	
                                                                                                                                      
History:
        12jul   nsh     73 -- Began coding
	24jul	nsh	73 -- Included the shull coefficients in the chianti database
	
                                                                                                                                      
**************************************************************/

double
total_rrate (nion, T)
     int nion;
     double T;
{


  double rate;			//The returned rate
  double rrA, rrB, rrT0, rrT1, rrC, rrT2;	//The parameters
  double term1, term2, term3;	//Some temporary parameters to make calculation simpler


  rate = 0.0;			/* NSH 130605 to remove o3 compile error */

  if (ion[nion].total_rrflag != 1)
    {
      Error ("total_rrate: No T_RR parameters for ion %i\n", nion);
      return (0);
    }
  if (total_rr[ion[nion].nxtotalrr].type == RRTYPE_BADNELL)
    {
      rrA = total_rr[ion[nion].nxtotalrr].params[0];
      rrB = total_rr[ion[nion].nxtotalrr].params[1];
      rrT0 = total_rr[ion[nion].nxtotalrr].params[2];
      rrT1 = total_rr[ion[nion].nxtotalrr].params[3];
      rrC = total_rr[ion[nion].nxtotalrr].params[4];
      rrT2 = total_rr[ion[nion].nxtotalrr].params[5];


      rrB = rrB + rrC * exp ((-1.0 * rrT2) / T);	//If C=0, this does nothing


      term1 = sqrt (T / rrT0);
      term2 = 1.0 + sqrt (T / rrT0);
      term2 = pow (term2, (1 - rrB));
      term3 = 1.0 + sqrt (T / rrT1);
      term3 = pow (term3, (1 + rrB));


      rate = pow ((term1 * term2 * term3), -1.0);
      rate *= rrA;
    }
  else if (total_rr[ion[nion].nxtotalrr].type == RRTYPE_SHULL)
    {
      rate =
	total_rr[ion[nion].nxtotalrr].params[0] * pow ((T / 1.0e4),
						       total_rr[ion
								[nion].
								nxtotalrr].
						       params[1]);
    }
  else
    {
      Error ("total_rrate: unknown parameter type for ion %i\n", nion);
    }


  return (rate);


}


/***********************************************************
                                       Southampton University
                                                                                                                                      
 Synopsis:
        bad_gs_rr(nion,T)
                                                                                                                                      
Arguments:
        ion - ion for which we want a recombination rate - 
		this is the upper state, so the ion which is 
		doing the recombining, there is no rate for
		H1(ion0) but there is one for H(ion1)
	temperature - the temperature we want a rate for
Returns:
	rate - the resolved recombination rate for the ground
		state of this ion recombining into the GS of
		the lower ion.
		for this temperature
                                                                                                                                      
Description:
                                                                                                                                      
        This routine generates a total recombination rate for 
		a given ion at a given temperature using 
		badnell type parameters. We simply perform an
		interpolation 
                                                                                                                                      
Notes:
	
                                                                                                                                      
History:
        12jul   nsh     73 -- Began coding
  	14mar	nsh	77a-- Interpolaion now carried out in log space

	
                                                                                                                                      
**************************************************************/

double
badnell_gs_rr (nion, T)
     int nion;
     double T;
{
  double rate, drdt, dt;
  int i, imin, imax;
  double rates[BAD_GS_RR_PARAMS], temps[BAD_GS_RR_PARAMS];


  imin = imax = 0;		/* NSH 130605 to remove o3 compile error */


  if (ion[nion].bad_gs_rr_t_flag != 1 && ion[nion].bad_gs_rr_r_flag != 1)
    {
      Error ("bad_gs_rr: Insufficient GS_RR parameters for ion %i\n", nion);
      return (0);
    }


  for (i = 0; i < BAD_GS_RR_PARAMS; i++)
    {
      rates[i] = bad_gs_rr[ion[nion].nxbadgsrr].rates[i];
      temps[i] = bad_gs_rr[ion[nion].nxbadgsrr].temps[i];
    }

  if (T < temps[0])		//we are below the range of GS data
    {
      Log_silent
	("bad_gs_rr: Requested temp %e is below limit of data for ion %i(Tmin= %e)\n",
	 T, nion, temps[0]);
//      rate = rates[0];
     	imax=1;
	imin=0;  
  }

  else if (T >= temps[BAD_GS_RR_PARAMS - 1])	//we are above the range of GS data
    {
      Log_silent
	("bad_gs_rr: Requested temp %e is above limit (%e) of data for ion %i\n",
	 T, nion, bad_gs_rr[ion[nion].nxbadgsrr].temps[BAD_GS_RR_PARAMS - 1]);
 //     rate = rates[BAD_GS_RR_PARAMS - 1];
	imax=BAD_GS_RR_PARAMS - 1;
	imin=BAD_GS_RR_PARAMS - 2;
//We will try to extrapolate.



    }
  else				//We must be within the range of tabulated data
    {
      for (i = 0; i < BAD_GS_RR_PARAMS - 1; i++)
	{
        if (temps[i] <= T && T < temps[i + 1])	//We have bracketed the correct temperature
	    {
	      imin = i;
	      imax = i + 1;
	    }
	}
/* NSH 140313 - changed the following lines to interpolate in log space */
    }
      drdt = (log10(rates[imax]) - log10(rates[imin])) / (log10(temps[imax]) - log10(temps[imin]));
      dt = (log10(T) - log10(temps[imin]));
      rate = pow(10,(log10(rates[imin]) + drdt * dt));



    



  return (rate);


}


/***********************************************************
                                       Southampton University
                                                                                                                                      
 Synopsis:
        milne_gs_rr(nion,T)
                                                                                                                                      
Arguments:
        ion - ion for which we want a recombination rate - 
		this is the upper state, so the ion which is 
		doing the recombining, there is no rate for
		H1(ion0) but there is one for H(ion1)
	temperature - the temperature we want a rate for
Returns:
	rate - the resolved recombination rate for the ground
		state of this ion recombining into the GS of
		the lower ion.
		for this temperature
                                                                                                                                      
Description:
                                                                                                                                      
        This routine generates a GS to GS recombination rate 
		for those ions we do not have badnell data for - it
		computed the rate using the milne relation.
                                                                                                                                      
Notes:
	
                                                                                                                                      
History:
        12jul   nsh     73 -- Began coding
	
                                                                                                                                      
**************************************************************/

double
milne_gs_rr (nion, T)
     int nion;
     double T;
{
  double rate;
  int ntmin, nvmin, n;
  double fthresh, fmax, dnu;

  rate = 0.0;			/* NSH 130605 to remove o3 compile error */


  fbt = T;
  fbfr = 2;

  if (ion[nion].phot_info == 1)	//topbase
    {
      ntmin = ion[nion].ntop_ground;
      fb_xtop = &phot_top[ntmin];
      fthresh = fb_xtop->freq[0];
      fmax = fb_xtop->freq[fb_xtop->np - 1];
      dnu = 100.0 * (fbt / H_OVER_K);
      if (fthresh + dnu < fmax)
	{
	  fmax = fthresh + dnu;
	}
      rate = qromb (fb_topbase_partial, fthresh, fmax, 1e-5);
    }
  else if (ion[nion].phot_info == 0)	// verner
    {
      nvmin = nion;
      n = nvmin;
      fb_xver = &xphot[ion[n].nxphot];
      fthresh = fb_xver->freq_t;
      fmax = fb_xver->freq_max;
      dnu = 100.0 * (fbt / H_OVER_K);
      if (fthresh + dnu < fmax)
	{
	  fmax = fthresh + dnu;
	}
      rate = qromb (fb_verner_partial, fthresh, fmax, 1e-5);
    }




  return (rate);


}
