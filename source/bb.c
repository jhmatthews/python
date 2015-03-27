
/***********************************************************
                                       Space Telescope Science Institute

 Synopsis:
	
double planck(t,freqmin,freqmax)	returns the frequency of a "pseudo-randomly generated photon. 
double t,freqmin,freqmax;		temperature of BB, and frequency limits
  
 					The frequency is pseudo random in the following limited sense.  
					The photons are returned are weighted so that the energy distribution 
 					of a planck function is approximately reproduced. Within an energy 
 					bin the photon frequency is uniform.  

double
emittance_bb (freqmin, freqmax, t)

double freqmin, freqmax, t;
					Calculate the emittance of a bb between freqmin and freqmax

double integ_planck_d(alphamin,alphamax)
double alphamin,alphamax;

					returns the integral of the dimensionless blackbody function
					between alphamin and alphamax.  

					Because we cannot afford to keep reintegrating this integral, 
					the following procedure is used.  The first time the function is
					called the integ_planck[] is filled.  integ_plank[n]  contains the 
					integral of the dimensionless planck function from  ALPHAMIN to 
					ALPHAMAX.  Therefore if one wants to obtain the integral of the 
					dimensionless bb function, one simply interpolates this array


init_integ_planck_d() 			performs the integration of the dimensionless bb function
					in order to fill the arry integ_planck.  This routine needs only to be
					 called onetime. The actuall integration is done with the numerical 
					 recipes routine "qromb" .

check_fmax(fmin,fmax,temp)		This is a little helper routine written by NSH in August 2012. We
					were having lots of problems with trying to integrate cross sections
					x a plack function at temperatures where there is no flux at fmax. This
					little routine just checks if fmax is going to be more than 100xt/(h/k)
					which will return tiny numbers for planck, and make qromb return
					nonsense.

Arguments:		

Returns:
 
Description:	
		The way this works currently is a little obscure.  The primary routines are emittance_bb or planck.  
		They call the other routines

		The first call to either of these routines (surely emittance_bb) results in a call to integ_planck_d, 
		which in turn calls integ_planck_init.  This initialization routine in turn populates the 
		array integ_planck , which contains the integral of the bb function in an array.

		bb_emittance continues to access the array integ_plank through integ_planck_d every 
		future time is is called.

		planck does the same thing albeit more indirectly. It sets up a pdf each time new frequency 
		limits are placed on it.  planck therefore really uses the
		pdf.

		
Notes:
		It seems possible that some of these routines could be combined in a clever way
		

History:
 	97	ksl	The various portions of this file were written in 1997 as part of the python effort
 	98mar	ksl	Combined the routines here from bb_gen.c and planck.c to reduce the total number of
 			files.
 	98mar23	ksl	planck(t,freqmin,freqmax) modified to use new pdf routines
 	98apr25	ksl	integ_planck_d modified to fix an error when the ranges of alphamin and alphamax
 			were out of bounds
	98jul27	ksl	planck moddified to force denser selection of gridpoints for small and high freq relative
			to kt
	01dec	ksl	Increased number of points at high end as part of exercise to force selection
			of some high frequency photons (python_40)
	04mar	ksl	Modified to provide a better treatment when the frequency range is
			above or below the frequency limits of the cdf that was calculated 
			initially.
	05jul	ksl	56d -- Relegated writing the pdf file to a case where you were debugging
	06aug	kls	57h -- Recommented and reorganized the file, so that could complete a partial
			fix to a problem Sturat had identified.  The partial fix was generating 
			photons outsde the min and maximum frequency
	12aug	nsh	73e -- check_fmax written to check qromb calls will work properly.
	12nov	ksl	74a_ksl -- Revised  rand_exp, which is used to allow us to reasonalbly 
			approximate photons in the high alpha limit.  Eliminated some old notes about
			this routine that seemed no longer very relevant.\
	13mar	jm	74b5_JM -- fixed bug JM130303 in init_integ_planck. This bug was causing 
			the emittance in the range wavemin-wavemax in the spectral cycles to be calculated 
			wrongly. Was noticed when testing YSO models.

**************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "atomic.h"
#include "python.h"

/***********************************************************
                                       Space Telescope Science Institute

 Synopsis: planck returns the frequency of a "pseudo-randomly generated photon.  The
   frequency is pseudo random in the following limited sense.  The photons are returned
   are weighted so that the energy distribution of a planck function is approximately
   reproduced. Within an energy bin the photon frequency is uniform.  
	

Arguments:		

     temperature, 
     freqmin, freqmax  - frequency interval 

Returns:

     The frequency of a photon  drawn randomly from a BB function with 
     temperature T
 
Description:	

	The first time one enters the program, a cdf for a diminensionaless
	BB function is created.  On subseqent entries, when the temperature
	or frequency limits are changed, we use standard routines to limit
	what portion of the dimensionless cdf to use.

	Because of of the rapid fall off of the proabability density function
	at very low and high frequency, special routines exist to deal with
	photons in these regions of the spectrum.

Notes:

	12nov - ksl - I have deleted some of the old notes associated with
	this routine.  When this routine was written, computers are a lot
	slower than they are today, and I was quite concerned about the
	possibility that this routine would take a lot of computer time.
	Computers are much faster today, and this routine does not 
	consitute a significant fraction of the time of the entire 
	program.  It is possible, that a more accurate routine could
	be created more simply today.

History:
	06aug	ksl	57h --Recommented, as moved to fix a problem
			that was returning photons that were occassionally outside
			of the frequency limits.  The problem had been introduced
			to handle the YSO case, where Stuart had found the uniform
			distribution used to interpolate between frequency points
			annoyoing
	06aug	ksl	Eliminated call to Stuart's attempt to better handle
			the low frequency limit, in favor of a new routine
			that I hope fixes the problem of the frequency
			limits.  
	06sep	ksl	Hopefully correct the remaining problems with this,
			but there is a significant issue with all of our pdfs after 
			we have created pdf arrays, we simply create a uniform 
			distribution within an interval.
	06sep	ksl	57i -- Have simplfied this to elimiante for now any lo or hi
			frequency special cases
	12nov	ksl	Removed some of the old Notes assocaited with this routine.
			See version earlier than 74 for these old notes.
**************************************************************/

#define ALPHAMIN 0.4		// Region below which we will use a low frequency approximation
#define ALPHAMAX 30.		// Region above which we will use a high frequency approximation
#define ALPHABIG 100.		//  Region over which can maximmally integrate the Planck function
#define NMAX 		1000

int ninit_planck = 0;
struct Pdf pdf_bb;

double old_t = 0;
double old_freqmin = 0;
double old_freqmax = 0;
double alphamin, alphamax;
double cdf_bb_lo, cdf_bb_hi, cdf_bb_tot;	// The precise boundaries in the the bb cdf 
double cdf_bb_ylo, cdf_bb_yhi;	// The places in the CDF defined by freqmin & freqmax
double lo_freq_alphamin, lo_freq_alphamax, hi_freq_alphamin, hi_freq_alphamax;	//  the limits to use for the low and high frequency values

// bb_set is thae array that pdf_gen_from_func uses to esablish the 
// specific points in the cdf of the dimensionless bb function.
double bb_set[] = {
  0.05, 0.1, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45,
  10.3, 11.3, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0,
  19., 20., 21., 22., 23., 24., 25., 26., 27., 28., 29.
};


int error_bb_hi = 0;
int error_bb_lo = 0;

double
planck (t, freqmin, freqmax)
     double t, freqmin, freqmax;
{
  FILE *fopen ();
  double freq, alpha, y;
  double planck_d (), pdf_get_rand_limit ();
  double get_rand_pow ();
  int pdf_gen_from_func (), pdf_to_file (), echeck;
  int pdf_limit ();


  /*First time through create the array containing the proper boundaries for the integral of the BB function,
     Note calling pdf_gen_from func also defines ylo and yhi */

  if (ninit_planck == 0)
    {				/* First time through p_alpha must be initialized */
      if ((echeck =
	   pdf_gen_from_func (&pdf_bb, &planck_d, ALPHAMIN, ALPHAMAX, 29,
			      bb_set)) != 0)
	{
	  Error ("Planck: on return from pdf_gen_from_func %d\n", echeck);
	}
      /* We need the integral of the bb function outside of the regions of interest as well */


      cdf_bb_tot = qromb (planck_d, 0, ALPHABIG, 1e-8);
      cdf_bb_lo = qromb (planck_d, 0, ALPHAMIN, 1e-8) / cdf_bb_tot;	//position in the full cdf of low frequcny boundary
      cdf_bb_hi = 1. - qromb (planck_d, ALPHAMAX, ALPHABIG, 1e-8) / cdf_bb_tot;	//postion in fhe full hi frequcny boundary

//      pdf_to_file (&pdf_bb, "pdf.out");
      ninit_planck++;

    }


/* If temperatures or frequencies have changed since the last call to planck
redefine various limitsi, including the region of the pdf  to be used

Note - ksl - 1211 - It is not obvious why all of these parameters need to be
reset.  A careful review of them is warranted.
*/

  if (t != old_t || freqmin != old_freqmin || freqmax != old_freqmax)
    {

      alphamin = H * freqmin / (BOLTZMANN * t);
      alphamax = H * freqmax / (BOLTZMANN * t);

      old_t = t;
      old_freqmin = freqmin;
      old_freqmax = freqmax;

      cdf_bb_ylo = cdf_bb_yhi = 1.0;
      if (alphamin < ALPHABIG)
	{
	  cdf_bb_ylo = qromb (planck_d, 0, alphamin, 1e-8) / cdf_bb_tot;	//position in the full cdf of current low frequcny boundary
	  if (cdf_bb_ylo > 1.0)
	    cdf_bb_ylo = 1.0;
	}
      if (alphamax < ALPHABIG)
	{
	  cdf_bb_yhi = qromb (planck_d, 0, alphamax, 1e-8) / cdf_bb_tot;	//position in the full cdf of currnt hi frequcny boundary
	  if (cdf_bb_yhi > 1.0)
	    cdf_bb_yhi = 1.0;
	}


/* These variables are not always used */
      lo_freq_alphamin = alphamin;	// Never used if 
      lo_freq_alphamax = alphamax;
      if (lo_freq_alphamax > ALPHAMIN)
	lo_freq_alphamax = ALPHAMIN;

      hi_freq_alphamax = alphamax;
      hi_freq_alphamin = alphamin;
      if (hi_freq_alphamin < ALPHAMAX)
	hi_freq_alphamin = ALPHAMAX;


      if (alphamin < ALPHAMAX && alphamax > ALPHAMIN)
	{
	  pdf_limit (&pdf_bb, alphamin, alphamax);
	}

    }
  /* End of section redefining limits */


  y = rand () / (MAXRAND);

  y = cdf_bb_ylo * (1. - y) + cdf_bb_yhi * y;	// y is now in an allowd place in the cdf

/* There are 3 cases to worry about
	The case where everything is in the low frequency limit
	The case where everything is in the normal limit
	The case where some photons are in the low regime and some are
	in the normal regime
*/

  if (y <= cdf_bb_lo || alphamax < ALPHAMIN)
    {
      alpha = get_rand_pow (lo_freq_alphamin, lo_freq_alphamax, 2.);
    }
  else if (y >= cdf_bb_hi || alphamin > ALPHAMAX)
    {
      alpha = get_rand_exp (hi_freq_alphamin, hi_freq_alphamax);
    }
  else
    {
      alpha = pdf_get_rand_limit (&pdf_bb);
    }

  freq = BOLTZMANN * t / H * alpha;
  if (freq < freqmin || freqmax < freq)
    {
      Error ("planck: freq %g out of range %g %g\n", freq, freqmin, freqmax);
    }
  return (freq);
}



/***********************************************************
	Space Telescope Science Institute

 Synopsis: get_rand_pow obtains a random number between x1 and x2 
 	for a power law densiity distribution with index alpha
	

Arguments:		

Returns:
 
Description:	

Notes:
		
History:
	06sep	ksl	Coded       
**************************************************************/
double
get_rand_pow (x1, x2, alpha)
     double x1, x2, alpha;
{
  double r;
  double a;

  r = rand () / MAXRAND;

  x1 = pow (x1, alpha + 1.);
  x2 = pow (x2, alpha + 1.);

  a = (1. - r) * x1 + r * x2;

  a = pow (a, 1. / (alpha + 1.));
  return (a);
}


/***********************************************************
	Space Telescope Science Institute

 Synopsis: get_rand_exp obtains a random number between x1 and x2 
 	for an exp law density distribution with index alpha. 
	

Arguments:		

	alpha_min and alpha_max, the range of the exponential
	function over which one wants the random number 
	to be returned

Returns:

	A double precision random number
 
Description:	

	
Notes:

	The cdf for an exponential distribution can be easily
	shown to be given by solving this equation for alpha

	r*(exp(-alpha_min)-exp(-alpha_max))=(exp(-alpha_min)-exp(alpha))

	but you can recast this and solve for delta_alpha

	exp(-delta_alpha)=(1-R)+R*exp(-(alpha_max-alpha_min))

	This has the advantage that it depends only on the 
	difference of alpha_max and alpha_min and not their
	actual values, and as long as the exp of a very 
	large number turns out to be zero and not 
	not a number, it shuld not generate sane checks
		
History:
	06oct	ksl	Coded       
	12nov	ksl	Added sane check to find a problem with
			returning NaN.  The problem was due to
			the fact that alpha_min and almax was so
			large that x1 and x2 were both identically
			zero and as a result the logarithm was
			negative.  
	12nov	ksl	Changed algorithm to minimimize the effects
			of very large values of alpha_min and 
			alpha_max.  Instead of calculating 
			alpha directly, we calculate delta_alpha
			the difference between alpha_min and the
			value we want
**************************************************************/
double
get_rand_exp (alpha_min, alpha_max)
     double alpha_min, alpha_max;
{
  double r;
  //Old ksl 12nov double x1, x2;
  double x;
  double a, aa;
  double delta_alpha;

  r = rand () / MAXRAND;

  x = exp (alpha_min - alpha_max);

  //OLD 12nov ksl x1 = exp (-alpha_min);
  //OLD 12nov ksl x2 = exp (-alpha_max);

  //OLD 12nov ksl aa = (1. - r) * x1 + r * x2;

  aa = (1. - r) + r * x;
  delta_alpha = -(log (aa));

  a = alpha_min + delta_alpha;

  if (sane_check (a))
    {
      Error ("get_rand_exp:sane_check %e %e %e %e %e\n", a, aa, delta_alpha,
	     x, r);
    }
  //OLD 12 nov ksl return (-a);
  return (a);
}

/***********************************************************
	Space Telescope Science Institute

 Synopsis: integ_plank_d(alphamin,alphamax) returns the integral of the dimensionless blackbody function
   between alphamin and alphamax.  
	

Arguments:		

Returns:
 
Description:	
	To save computing time, the routine actually accesses an array that contains the
	the integral of the bb function from 0 to x.

   	The first time the function is
	called the integ_planck[] is filled.  integ_plank[n]  contains the integral of the
	dimensionless planck function from  ALPHAMIN to ALPHAMAX.  Therefore if one
	wants to obtain the integral of the dimensionless bb function, one simply interpolates
	this array. 
Notes:
		
History:
	06aug	ksl	Recommented
**************************************************************/


double integ_planck[NMAX + 1];
int i_integ_planck_d = 0;
double
integ_planck_d (alphamin, alphamax)
     double alphamin, alphamax;
{
  double x, z1, z2;
  int n;
  int init_integ_planck_d ();
  if (i_integ_planck_d == 0)
    {				/*First time through integ_planck must be defined */
      init_integ_planck_d ();
      i_integ_planck_d++;
    }

  x = (alphamin - ALPHAMIN) / (ALPHAMAX - ALPHAMIN) * NMAX;
  if (x <= 0.0)
    z1 = 0.0;
  else if (x >= (NMAX))
    {
      return (0.0);		/* Because the minimum frequency is too high */
    }
  else
    {
      n = x;
      x -= n;
      z1 = integ_planck[n] * (1. - x) + integ_planck[n + 1] * x;
    }

  x = (alphamax - ALPHAMIN) / (ALPHAMAX - ALPHAMIN) * NMAX;
  if (x < 0.0)
    {
      return (0.0);		/* Because the maximum frequency is too low */
    }
  else if (x >= (NMAX))
    {
      z2 = integ_planck[NMAX];
    }
  else
    {
      n = x;
      x -= n;
      z2 = integ_planck[n] * (1. - x) + integ_planck[n + 1] * x;
    }

  return (z2 - z1);
}


/***********************************************************
       Space Telescope Science Institute

 Synopsis: init_integ_planck_d calulates integrals of the dimensionaless bb function and
 	stores it in the array integ_planck.  Each value of integ_plank contains
	the integration from 0 to x
	

Arguments:		

Returns:
 
Description:	
   This routine needs to be called once, since ALPHAMIN and ALPHAMAX are
   hardcoded. 

   The actual integration is done with the numerical recipes routine "qromb" 

Notes:
		
History:
	97+	ksl	Originally coded
	06aug	ksl	Revised comments

**************************************************************/

int
init_integ_planck_d ()
{
  double x;
  double planck_d (), qromb ();
  int n;
  //OLD74b5 integ_planck[0] = 0;   JM130319: this should be set to ALPHAMIN- done in the for loop for simplicity (n=0).
  for (n = 0; n <= NMAX + 1; n++)
    {
      x = ALPHAMIN + n * (ALPHAMAX - ALPHAMIN) / NMAX;
// 1e-7 is the fractional accuracy in my modified version of qromb -- ksl
      integ_planck[n] = qromb (planck_d, 0.0, x, 1e-7);
    }

  return (0);
}


/***********************************************************
	Space Telescope Science Institute

 Synopsis:
"planck_d" is the dimensionless bb function.  The total emittance
   is related to planck_d as follows:
	

Arguments:		

Returns:
 
Description:	
   F_nu= 2*PI* (kT/h**3)*planck_d(h*freq/ k T)
Notes:
		
History:
	06aug	ksl	Recommented
**************************************************************/

#define EPSILON	1.e-6

double
planck_d (alpha)
     double alpha;
{
  double x;
  if (alpha < EPSILON || alpha > ALPHABIG)
    return (0);
  x = (alpha * alpha * alpha) / (exp (alpha) - 1);
  return (x);
}

// Calculate the emittance of a bb between freqmin and freqmax
// Should integrate to sigma 
double
emittance_bb (freqmin, freqmax, t)
     double freqmin, freqmax, t;
{
  double alphamin, alphamax, q1;
  double integ_planck_d ();
  q1 =
    2. * PI * (BOLTZMANN *
	       BOLTZMANN * BOLTZMANN * BOLTZMANN) / (H * H * H * C * C);


  alphamin = H * freqmin / (BOLTZMANN * t);
  alphamax = H * freqmax / (BOLTZMANN * t);
  if (alphamin > ALPHAMIN && alphamax < ALPHAMAX)
    {
      return (q1 * t * t * t * t * integ_planck_d (alphamin, alphamax));
    }
  else if (alphamax < ALPHAMIN)
    {
      return (q1 * t * t * t * t *
	      qromb (planck_d, alphamin, alphamax, 1e-7));
    }
  else
    {

      return (q1 * t * t * t * t * integ_planck_d (alphamin, alphamax));
    }
}


/***********************************************************
	Southampton University

 Synopsis:
check_fmax decides wether a maximum frequency requested for an integral is sensible.
If it is too far off the end of the planck function, qromb will malfunction. We
just have to set it to a frequency where the BB function is tiny, say where hnu/kT =100.
At this point the bb function is


	

Arguments:		

Returns:
 
Description:	
We use alphabig to define the place in the BB spectrum where we want to give up
Notes:
		
History:
	12aug	nsh	written
**************************************************************/


double
check_fmax (fmin, fmax, temp)
     double fmin, fmax, temp;
{
  double bblim;

  bblim = ALPHABIG * (temp / H_OVER_K);	/*This is the frequency at which the exponent in the 
					   planck law will be -100. This will give a *very* small b(nu). */
  if (bblim < fmax)
    {
      fmax = bblim;
    }

  return (fmax);

}




#undef NMAX
#undef ALPHAMIN
#undef ALPHAMAX
