
#include <stdio.h>
#include <strings.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "atomic.h"
#include "python.h"

#include "log.h"


double constant;
double T_b;



double
emittance_brem (freqmin, freqmax, lum, t)
     double freqmin, freqmax, lum, t;
{
  double emit;
  
  
  emit =  qromb(integ_brem,freqmin,freqmax,1e-4);

  return (emit);
}




double
integ_brem (freq)
     double freq;
{
  double answer;
  answer = geo.const_agn * pow(freq,geo.brem_alpha) * exp ((-1.0 * H * freq) / (BOLTZMANN * geo.brem_temp));
  return (answer);
}


double
	brem_d(alpha)
		double alpha;
{
	double answer;
	answer=pow(alpha,geo.brem_alpha)*exp(-1.0*alpha);
	return(answer);
}





#define BREM_ALPHAMIN 0.01		// Region below which we will use a low frequency approximation
#define BREM_ALPHAMAX 2.		// Region above which we will use a high frequency approximation
#define BREM_ALPHABIG 100.		//  Region over which can maximmally integrate the Planck function
#define NMAX 		1000

int ninit_brem = 0;
struct Pdf pdf_brem;

double old_brem_t = 0;
double old_brem_freqmin = 0;
double old_brem_freqmax = 0;
double brem_alphamin, brem_alphamax;
double cdf_brem_lo, cdf_brem_hi, cdf_brem_tot;	// The precise boundaries in the the bb cdf 
double cdf_brem_ylo, cdf_brem_yhi;	// The places in the CDF defined by freqmin & freqmax
double brem_lo_freq_alphamin, brem_lo_freq_alphamax, brem_hi_freq_alphamin, brem_hi_freq_alphamax;	//  the limits to use for the low and high frequency values

// bb_set is the array that pdf_gen_from_func uses to esablish the 
// specific points in the cdf of the dimensionless bb function.
double brem_set[] = {
	  0.05, 0.1, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45,
	  0.9,1.9
	};



int error_brem_hi = 0;
int error_brem_lo = 0;

double
get_rand_brem (freqmin, freqmax)
     double freqmin, freqmax;
{
  double freq, alpha, y;
  int echeck;
  
  

  /*First time through create the array containing the proper boundaries for the integral of the brem function,
     Note calling pdf_gen_from func also defines ylo and yhi */

  if (ninit_brem == 0)
    {				/* First time through p_alpha must be initialized */
      if ((echeck =
	   pdf_gen_from_func (&pdf_brem, &brem_d, BREM_ALPHAMIN, BREM_ALPHAMAX, 10,
			      brem_set)) != 0)
	{
	  Error ("get_rand_brem: on return from pdf_gen_from_func %d\n", echeck);
	}
	
      /* We need the integral of the brem function outside of the regions of interest as well */
	
      cdf_brem_tot = (pow(BREM_ALPHAMIN/100.0,0.8))/0.8+qromb (brem_d, BREM_ALPHAMIN/100.0, BREM_ALPHABIG, 1e-8);
      cdf_brem_lo = (pow(BREM_ALPHAMIN/100.0,0.8))/0.8+qromb (brem_d, BREM_ALPHAMIN/100.0, BREM_ALPHAMIN, 1e-8) / cdf_brem_tot;	//position in the full cdf of low frequcny boundary
      cdf_brem_hi = 1. - qromb (brem_d, BREM_ALPHAMAX, BREM_ALPHABIG, 1e-8) / cdf_brem_tot;	//postion in fhe full hi frequcny boundary

//      pdf_to_file (&pdf_bb, "pdf.out");
      ninit_brem++;

    }

/* If temperatures or frequencies have changed since the last call to planck
redefine various limitsi, including the region of the pdf  to be used

Note - ksl - 1211 - It is not obvious why all of these parameters need to be
reset.  A careful review of them is warranted.
*/

  if (geo.brem_temp != old_brem_t || freqmin != old_brem_freqmin || freqmax != old_brem_freqmax)
    {
      brem_alphamin = H * freqmin / (BOLTZMANN * geo.brem_temp);
      brem_alphamax = H * freqmax / (BOLTZMANN * geo.brem_temp);

      old_brem_t = geo.brem_temp;
      old_brem_freqmin = freqmin;
      old_brem_freqmax = freqmax;

      cdf_brem_ylo = cdf_brem_yhi = 1.0;
      if (brem_alphamin < BREM_ALPHABIG) //There is *some* emission - 
	{
		if (brem_alphamin>BREM_ALPHAMIN/100.) //the requested range spills into the part where we must do some qromb
		{
			cdf_brem_ylo = ((pow(BREM_ALPHAMIN/100.0,0.8))/0.8+qromb (brem_d,BREM_ALPHAMIN/100.0, brem_alphamin, 1e-8)) / cdf_brem_tot;
		}
		else
		{
			cdf_brem_ylo = ((pow(BREM_ALPHAMIN/100.0,0.8))/0.8) / cdf_brem_tot;//position in the full cdf of current low frequcny boundary		
		}
		if (cdf_brem_ylo > 1.0)
			cdf_brem_ylo = 1.0;
	}
      if (brem_alphamax < BREM_ALPHABIG)
	{
		if (brem_alphamax>BREM_ALPHAMIN/100.)
		{
			cdf_brem_yhi = ((pow(BREM_ALPHAMIN/100.0,0.8))/0.8+qromb (brem_d, BREM_ALPHAMIN/100.0, brem_alphamax, 1e-8)) / cdf_brem_tot;	//position in the full cdf of currnt hi frequcny boundary
		}
		else
		{
			cdf_brem_yhi = ((pow(BREM_ALPHAMIN/100.0,0.8))/0.8)/cdf_brem_tot;
		}
	  if (cdf_brem_yhi > 1.0)
	    cdf_brem_yhi = 1.0;
	}


/* These variables are not always used */
      brem_lo_freq_alphamin = brem_alphamin;	// Never used if 
      brem_lo_freq_alphamax = brem_alphamax;
      if (brem_lo_freq_alphamax > BREM_ALPHAMIN)
	brem_lo_freq_alphamax = BREM_ALPHAMIN;

      brem_hi_freq_alphamax = brem_alphamax;
      brem_hi_freq_alphamin = brem_alphamin;
      if (brem_hi_freq_alphamin < BREM_ALPHAMAX)
	brem_hi_freq_alphamin = BREM_ALPHAMAX;


      if (brem_alphamin < BREM_ALPHAMAX && brem_alphamax > BREM_ALPHAMIN)
	{
	  pdf_limit (&pdf_brem, brem_alphamin, brem_alphamax);
	}

    }
  /* End of section redefining limits */

  y = rand () / (MAXRAND);

  y = cdf_brem_ylo * (1. - y) + cdf_brem_yhi * y;	// y is now in an allowd place in the cdf

/* There are 3 cases to worry about
	The case where everything is in the low frequency limit
	The case where everything is in the normal limit
	The case where some photons are in the low regime and some are
	in the normal regime
*/

  if (y <= cdf_brem_lo || brem_alphamax < BREM_ALPHAMIN)
    {
      alpha = get_rand_pow (brem_lo_freq_alphamin, brem_lo_freq_alphamax,geo.brem_alpha);
    }
  else if (y >= cdf_brem_hi || brem_alphamin > BREM_ALPHAMAX)
    {
      alpha = get_rand_exp (brem_hi_freq_alphamin, brem_hi_freq_alphamax);
    }
  else
    {
      alpha = pdf_get_rand_limit (&pdf_brem);
    }

  freq = BOLTZMANN * geo.brem_temp / H * alpha;
  if (freq < freqmin || freqmax < freq)
    {
      Error ("get_rand_brem: freq %g out of range %g %g\n", freq, freqmin, freqmax);
    }
  return (freq);
}


