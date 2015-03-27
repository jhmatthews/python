/* These are subroutines of balance which are not related directly to ionization
calculations 
	01dec	ksl	Updated to reflect new calls to various routines, notably  
			scattering_fraction.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "atomic.h"
#include "python.h"




// Make a bunch of BB photons 
/* The band structure is intended to be used as follows:
        freq is lower limit frequency for band[n]
        min_fraction is minimum number of photos for band[n]
        nat_fraction is the fractionage of the total that would naturally occur
                if all photons were equally weighted
        used_fraction is the result of reallocation
        weight is the weight of photons in each band
        nphot is number of photons in each band
History
	02jul	ksl	Adapted to new approach to banding
*/
//02jul #define NBANDS 4
//02jul struct bands
//02jul {
//02jul   double freq;
//02jul   double min_fraction;
//02jul   double nat_fraction;          // The fraction of the accepted luminosity in this band
//02jul   double used_fraction;
//02jul   double f;                     //The "luminosity" within a band
//02jul   double weight;
//02jul   int nphot;
//02jul }
//02jul band[NBANDS + 1], uband[NBANDS + 1];

//02jul int nbands, mbands;             // Actual number of bands


int
xbb (p, t, weight, f1, f2, freq_sampling)
     PhotPtr p;
     double t;
     double weight;
     double f1, f2;
     int freq_sampling;
{

  int n, most;
  double ftot, z;
  double frac_used;
  double emittance_bb ();
  int iphot_start;
  int nphottot;
  int xmake_bb ();
  if (freq_sampling == 0)	// Old uniform approach
    {
      /* The weights are normalized to the energy density of a BB spectrum * the speed
         of light */
      weight = 4. * STEFAN_BOLTZMANN * t * t * t * t * weight / NPHOT;
      xmake_bb (p, t, f1, f2, weight, 0, NPHOT);
    }
  else				// Use band limited frequency sampling
    {
/* Extra stuff that will be moved out of here eventually */
//Initialize bands
//02jul      nbands = 4;
//02jul      band[0].freq = 0;
//02jul      band[1].freq = 13.599 / HEV;
//02jul      band[2].freq = 24.588 / HEV;
//02jul      band[3].freq = 54.418 / HEV;
//02jul      band[4].freq = 1.e50;
//02jul      band[0].min_fraction = 0;
//02jul      band[1].min_fraction = 0.1;
//02jul      band[2].min_fraction = 0.1;
//02jul      band[3].min_fraction = 0.1;
//02jul      band[4].min_fraction = 0.1;
// This ends the extra stuff that would be done on overall initialization

//02jul      memcpy (uband, band, sizeof (band));       // Copy original to a working array
//02jul      for (n = 0; n < nbands; n++)
//02jul {
//02jul   if (uband[n].freq < f1)
//02jul     uband[n].freq = f1;
//02jul   if (uband[n + 1].freq > f2)
//02jul     uband[n + 1].freq = f2;
//02jul }
// So now the frequencies are chopped off and we can calculate the band-limite luminosities 

      ftot = 0.0;
      for (n = 0; n < xband.nbands; n++)	// Now get the band limited luminosities
	{
//02jul   if (uband[n].freq < uband[n + 1].freq)
	  if (xband.f1[n] < xband.f2[n])
	    {

//02jul       ftot += uband[n].f =
//02jul         emittance_bb (band[n].freq, band[n + 1].freq, t);
	      ftot += xband.flux[n] =
		emittance_bb (xband.f1[n], xband.f2[n], t);
	      printf ("BAND %i from %e to %e at temp %f has flux %e\n", n,
		      xband.f1[n], xband.f2[n], t, xband.flux[n]);
	    }
	  else
//02jul     uband[n].f = 0.0;
//02jul   if (uband[n].f == 0.0)
//02jul     uband[n].min_fraction = 0;  //Because you will not be able to generate photons
	    xband.flux[n] = 0.0;
	  if (xband.flux[n] == 0.0)
	    xband.min_fraction[n] = 0;	//Because you will not be able to generate photons
	}
/* So now we can distribute the photons */
      frac_used = 0;
      for (n = 0; n < xband.nbands; n++)
	{
//02jul   uband[n].nat_fraction = uband[n].f / ftot;
//02jul   frac_used += uband[n].min_fraction;
	  xband.nat_fraction[n] = xband.flux[n] / ftot;
	  frac_used += xband.min_fraction[n];
	}
      nphottot = 0;
      z = 0;
      for (n = 0; n < xband.nbands; n++)
	{
//02jul   uband[n].used_fraction =
//02jul     uband[n].min_fraction + (1 - frac_used) * uband[n].nat_fraction;
//02jul   nphottot += uband[n].nphot = NPHOT * uband[n].used_fraction;
	  xband.used_fraction[n] =
	    xband.min_fraction[n] + (1 - frac_used) * xband.nat_fraction[n];

	  nphottot += xband.nphot[n] = NPHOT * xband.used_fraction[n];

//02jul   if (uband[n].used_fraction > z)
	  if (xband.used_fraction[n] > z)
	    {
//02jul       z = uband[n].used_fraction;
	      z = xband.used_fraction[n];
	      most = n;
	    }
	}

/* Because of roundoff errors nphottot may not sum to the desired value, namely NPHOT.  So
add a few more photons to the band with most photons already. It should only be a few, at most
one photon for each band.*/
      if (nphottot < NPHOT)
	{
//02jul   uband[most].nphot += (NPHOT - nphottot);
	  xband.nphot[most] += (NPHOT - nphottot);
	}

//Calculate the weight for each photon assuming no banding 
      weight = 4. * STEFAN_BOLTZMANN * t * t * t * t * weight / NPHOT;




// Now generate the photons

      iphot_start = 0;
      for (n = 0; n < xband.nbands; n++)
	{
//02jul   if (uband[n].nphot > 0)
	  if (xband.nphot[n] > 0)
	    {
//02jul       uband[n].weight =
//02jul         uband[n].nat_fraction / uband[n].used_fraction * weight;
//02jul       xmake_bb (p, t, uband[n].freq, uband[n + 1].freq,
//02jul                 uband[n].weight, iphot_start, uband[n].nphot);
//02jul       iphot_start += uband[n].nphot;
	      xband.weight[n] =
		xband.nat_fraction[n] / xband.used_fraction[n] * weight;
	      xmake_bb (p, t, xband.f1[n], xband.f2[n],
			xband.weight[n], iphot_start, xband.nphot[n]);
	      iphot_start += xband.nphot[n];
	    }
	}
    }
  return (0);
}

// This is a stripped down routine. It doesn't think. It just makes photons
int
xmake_bb (p, t_r, freqmin, freqmax, weight, iphot_start, nphot)
     PhotPtr p;
     double t_r;
     double freqmin, freqmax;
     double weight;
     int iphot_start, nphot;	// Parallels photon gen
{
  int n;
  int iphot_stop;
  double planck ();



  iphot_stop = iphot_start + nphot;
  for (n = iphot_start; n < iphot_stop; n++)
    {

      p[n].freq = planck (t_r, freqmin, freqmax);
      p[n].w = weight;
      p[n].x[0] = wind_midx[0];
      p[n].x[1] = wind_midx[0];
      p[n].x[2] = EPSILON;
      p[n].lmn[0] = 0;
      p[n].lmn[1] = 0;
      p[n].lmn[2] = 1;
    }


  return (0);

}
