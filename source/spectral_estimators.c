

/***********************************************************
                                       Space Telescope Science Institute

 Synopsis:
 	The routines in this file are all associated with the process of
	modelling the spectrum in a cell. It uses the mean frequency to
	find a power law model and an exponential model. It then uses
	the standard deviation to decide which is best.

	The main routine is spectral_estimators and it currently attempts
	to model the spectrum in each band in the cell an exponetial or
	a power law.  The routine is called from ionization.c

	The results are stored in parameters in the PlasmaPtr structure 
	in variables, such as spec_mod_type, pl_alpah, pl_log,exp_temp, etc. 
	
 Arguments:		
	WindPtr w;

Returns:
 
Description:	
	
Notes:

History:
	111126	ksl	This was originally coded by nsh, but it
			was inserted directly into what was intended
			to be the steering routine for all ionization
			calculations.  I have moved it to its own
			routine in python_71
	111227	ksl	Eliminated n as an external variable
			because it really is not one
	111229	ksl	Small modifications made to reflect moving
			nxfreq and xfreq into the geo structure
        120802  nsh	Recoded so all this does is compute the estimators
			this is needed for other PL modes.
        120817  nsh     Changed many things to allow for an exponential model
			to be produced. Changed name of code to spec_estimators
			and changed file it is stored in to spectral_estimators.c
	130729	nsh	Changed some sane checks that were generating hugh numbers 
			errors, when they were just checking for reasonbleness.
			Also, changed the way the code deals with no photons in	
			band. Now only an error if we actually generate photons there!
	130907	nsh	Changed two errors into warnings - no photons in band is
			not really an error - it is perfectly reasonable so we
			really dont want the code to stop if it happens too often.
			However, we do want to know about it!
	131107	nsh	Some significant changes to the spectral models. The PL
			model was changed to work in log space, so it can now cope
			with very large exponenets, which means we can model much
			more extreme edges. Also, observed upper and lower limits
			have been implemented - so the spectral models are only
			applied where photon have actually been seen. Finally,
			the EXP model has been changed to allow positive exponents, 
			and also a check has been included to avoid T becomming 
			too small for a given (high) frequency, and hence returning 0, which
			gave inf errors. All these changes should mean that we no 
l			longer discard models on the basis of numerical limitations,
			so some of the log_silents have been replaced by 
			warnings so I can see if we do still get numberical problems.


**************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "atomic.h"
#include "python.h"

/* external variables set up so zbrent can solve for alpha using sim_appha as integrand. NSH120817 Changed the names to remove reference to sim */
double spec_numin, spec_numax, spec_numean;
double lspec_numin, lspec_numax; //Log versions of numin and numax


int
spectral_estimators (xplasma)
     PlasmaPtr xplasma;
{
  double pl_alpha_min, pl_alpha_max, pl_alpha_temp, pl_w_temp, j;
  double exp_temp_min, exp_temp_max;	/*120817 the 'temperature' range we are going to search for an effective temperature for the exponential model */
  double exp_temp_temp, exp_w_temp;	/*120817 the temporary values for temperature and weight of the exponential model */
  int n,n1;
  double pl_sd, exp_sd;		/*120817 Computed standard deviations for two models for comparison with true value */
//  double ALPHAMAX = 200.0;	/*120817 Something to make it a bit more obvious as to what values of the powerlaw exponent we consider reasonable */
  int plflag, expflag;		/*120817 Two flags to say if we have a reasonable PL or EXP model, set to 1 initially, -1 means there has been some failure that means we must not use this model, +1 means it is OK */
  double genmin,genmax;
  double dfreq; /* NSH 130711 - a number to help work out if we have fully filled a band */

/* This call is after a photon flight, so we *should* have access to j and ave freq, and so we can calculate proper values for W and alpha
To avoid problems with solving, we need to find a reasonable range of values within which to search for a solution to eq18. A reasonable guess is that it is around the current value....
*/

genmin=xband.f1[0];
genmax=xband.f2[xband.nbands-1];

/* NSH 131108 The next few lines just work out which bands actually should have photons in - i.e. if we dont generate any photons in a band, we will not worry if there are no photons in that band */

for (n1 =0; n1<xband.nbands; n1++)
	{
	if (xband.nphot[n1] > 0)
	 	{
		genmin = xband.f1[n1]; //the first band with any photons will get set to 
		break;
		}
	}

for (n1 =xband.nbands-1; n1>-1; n1--)
	{
	if (xband.nphot[n1] > 0) 
		{
		genmax = xband.f2[n1]; //the first band going down will define the highest freq
		break;
		}
	}

/* We loop over all of the bands, the first band is band number 0, and the last is band nxfreq-1 */
/* 71 - 111229 - ksl - Small modification to reflect moving nxfreq, etc in the geo structure */
  for (n = 0; n < geo.nxfreq; n++)

    {
	Log_silent ("Starting out band %i in cell %i. mean=%e, sd=%e, minfreq=%e, maxfreq=%e, nphot=%i\n",n,xplasma->nplasma,xplasma->xave_freq[n],xplasma->xsd_freq[n],xplasma->fmin[n],xplasma->fmax[n],xplasma->nxtot[n]);

      plflag = expflag = 1;	//Both potential models are in the running
      if (xplasma->nxtot[n] <= 1) /*Catch the situation where there are only 1 or 0 photons in a band - we cannt reasonably try to model this situation*/ 
		{
	  	if (geo.xfreq[n] >= genmax || geo.xfreq[n+1] <= genmin) /*The band is outside where photons were generated, so not very surprising that there are no photons - just generate a log */
		  	{
			Log_silent("spectral_estimators: too few photons (1 or 0) in cell %d band %d but we weren't expecting any \n",xplasma->nplasma,n); //This is just a log, clearly it is not a problem//
			}
		else
			{
			Error("spectral_estimators: too few photons (1 or 0) in cell %d band %d to produce a model\n",xplasma->nplasma, n);
/* NSH 130709 - changed this to be a warning, there are photons produced here, so the fact that there are none getting into a cell tells us something - it may be perfectly reasonable but nice to know  */
			}
//	  	xplasma->pl_w[n] = 0;	//We also want to make sure that the weight will be zero, this way we make sure there is no contribution to the ionization balance from this frequency.
                xplasma->pl_log_w[n] = -999; //A very tiny weight - changed to take account of the new Log formulation
	  	xplasma->pl_alpha[n] = 999.9;	//Give alpha a value that will show up as an error
	  	xplasma->exp_w[n] = 0.0;	//Make sure that w is zero, s no chance of mucking up ionization balance
	  	xplasma->exp_temp[n] = 1e99;	//Give temp a value that will show up as an error
	  	xplasma->spec_mod_type[n] = -1;	//This tells the code that we have failed to model the spectrum in this band/cell/
		}
      else if (xplasma->fmax[n]==xplasma->fmin[n]) /*If all the photons in the cell are concentrated in a tiny range then we will also not expect to make a sensible model - this check could be reviewed later if lots of warning are produced */
		{
		Error("spectral_estimators: multiple photons but only one frequency seen in %d band %d\n",xplasma->nplasma,n); /*Flag as a warning, so one can see if it is an issue */
                xplasma->pl_log_w[n] = -999; //A very tiny weight
	  	xplasma->pl_alpha[n] = 999.9;	//Give alpha a value that will show up as an error
	  	xplasma->exp_w[n] = 0.0;	//Make sure that w is zero, s no chance of mucking up ionization balance
	  	xplasma->exp_temp[n] = 1e99;	//Give temp a value that will show up as an error
	  	xplasma->spec_mod_type[n] = -1;	//This tells the code that we have failed to model the spectrum in this band/cell/
		}
      else
	{

/* NSH 131108 - these line no longer needed, since we are logging the maximum and minimum frequencies in each baand, rather than globally */

//	  spec_numin = geo.xfreq[n];	/*1108 NSH n is defined in python.c, and says which band of radiation estimators we are interested in using the for power law ionisation calculation */
//	  if (xplasma->max_freq < geo.xfreq[n + 1])
//	    {
//	      Log_silent
//		("NSH resetting max frequency of band %i from %e to %e due to lack of photons\n",
//		 n, geo.xfreq[n + 1], xplasma->max_freq);
//	      spec_numax = xplasma->max_freq;
//	    }
//	  else
//	    {
//	      spec_numax = geo.xfreq[n + 1];
//	    }


/* NSH 131107 The next lines check and assign band limits. If the min/max frequency in a cell is within (fmax-fmin)/(sqrt(nphot)) of the end of
 a band, we say that the photons fill the band to that end - i.e. the fact we didnt see the minimum frequency is just because of photon numbers. */

	  dfreq = (geo.xfreq[n+1] - geo.xfreq[n])  / sqrt(xplasma->nxtot[n]); //This is a measure of the spacing between photons on average
	  if ((xplasma->fmin[n] - geo.xfreq[n]) < dfreq)
		{
		spec_numin = geo.xfreq[n];
		}
	  else
		{		
		spec_numin = xplasma->fmin[n];
		}
	  if ((geo.xfreq[n+1] - xplasma->fmax[n]) < dfreq)
		{
		spec_numax = geo.xfreq[n+1];
		}
	  else
		{		
		spec_numax = xplasma->fmax[n];
		}
	  xplasma->fmin_mod[n]=spec_numin; //This is the low frequency limit of any model we might make
	  xplasma->fmax_mod[n]=spec_numax; //This is the high frequency limit of any model we might make
	  lspec_numax = log10(spec_numax);	  
	  lspec_numin = log10(spec_numin);
	  spec_numean = xplasma->xave_freq[n];
	  j = xplasma->xj[n];


//	 Log
//	    ("NSH We are about to calculate w and alpha, band %i cell %i j=%10.2e, mean_freq=%10.2e, numin=%10.2e(%8.2fev), numax=%10.2e(%8.2fev), //number of photons in band=%i\n",
//	     n, xplasma->nplasma, j, spec_numean, spec_numin, spec_numin * HEV, spec_numax,
//	     spec_numax * HEV, xplasma->nxtot[n]);

//Log_flush();



	  pl_alpha_min = -0.1;	/*Lets just start the search around zero */
	  pl_alpha_max = +0.1;

//printf ("initial guess (lin) alpha=%f, pl_alpha_func=%e\n",pl_alpha_min,pl_alpha_func (pl_alpha_min));
//printf ("initial guess (log) alpha=%f, pl_alpha_func=%e\n",pl_alpha_min,pl_alpha_func_log (pl_alpha_min));

	  while (pl_alpha_func_log (pl_alpha_min) * pl_alpha_func_log (pl_alpha_max) >
		 0.0)
	    {
	      pl_alpha_min = pl_alpha_min - 1.0;
	      pl_alpha_max = pl_alpha_max + 1.0;
	    }


	  if (isfinite(pl_alpha_func_log(pl_alpha_min))==0
	      || isfinite(pl_alpha_func_log (pl_alpha_max))==0)
	    {
	      Error
		("spectral_estimators: Alpha cannot be bracketed (%e %e)in band %i cell %i- setting w to zero\n",
		 pl_alpha_min,pl_alpha_max,n, xplasma->nplasma); //NSH 131108 - now a warning, this should no longer happen
//	      xplasma->pl_w[n] = 0.0;
	      xplasma->pl_log_w[n] = -999.0;
	      xplasma->pl_alpha[n] = -999.0;	//Set this to a value that might let us diagnose the problem
	      plflag = -1;
	    }
	  else
	    {
/* We compute temporary values for sim alpha and sim weight. This will allow us to 
 * check that they are sensible before reassigning them */

	      pl_alpha_temp =
		zbrent (pl_alpha_func_log, pl_alpha_min, pl_alpha_max, 0.00001);
//	      if (pl_alpha_temp > ALPHAMAX)
//		pl_alpha_temp = ALPHAMAX;	//110818 nsh check to stop crazy values for alpha causing problems
//	      if (pl_alpha_temp < -1. * ALPHAMAX)
//		pl_alpha_temp = -1. * ALPHAMAX;

/*This next line computes the PL weight using an external function. Note that xplasma->j already 
 * contains the volume of the cell and a factor of 4pi, so the volume sent to sim_w is set to 1 
 * and j has a factor of 4PI reapplied to it. This means that the equation still works in balance. 
 * It may be better to just implement the factor here, rather than bother with an external call.... */


//	      pl_w_temp = pl_w (j, pl_alpha_temp, spec_numin, spec_numax);
              pl_w_temp = pl_log_w (j, pl_alpha_temp, lspec_numin, lspec_numax);

	      if ((isfinite(pl_w_temp))==0)
		{
		  Error
		    ("spectral_estimators: New PL parameters (%e) unreasonable, using existing parameters. Check number of photons in this cell\n",pl_w_temp); //NSH 131108 - now a warning, this should no longer happen
		  plflag = -1;	// Dont use this model
//		  xplasma->pl_w[n] = 0.0;
		  xplasma->pl_log_w[n] = -999.0;
		  xplasma->pl_alpha[n] = -999.0;
		}
	      else
		{
		  xplasma->pl_alpha[n] = pl_alpha_temp;
		  xplasma->pl_log_w[n] = pl_w_temp;
		}
	    }




//	  exp_temp_min = xplasma->t_r * 0.9;	/*Lets just start the search around the radiation temperature in the cell */
	  exp_temp_min = ((H * spec_numax) / (BOLTZMANN) ) * 0.9; /*NSH 131107 -  change here - we will start the search around the temperature that we know will yield a sensible answer  */
//	  exp_temp_max = xplasma->t_r * 1.1;
	  exp_temp_max = ((H * spec_numax) / (BOLTZMANN) ) / 0.9; /*NSH 131107 - and the same for the maximum temp */
	  while ((exp_temp_func (exp_temp_min) * exp_temp_func (exp_temp_max) > 0.0)  &&  ((exp_temp_func (-1.0*exp_temp_min) * exp_temp_func (-1.0*exp_temp_max) > 0.0))) /* NSH 131107 - changed to permit a negative temperautre, which will give a positive exponential */
	    {
	      if ((H * spec_numax) < (100.0 * BOLTZMANN * exp_temp_min * 0.9)) /*In this case we are going to get errors since the temperature is too to give a result in the exponential, and we will divide by zero */
		{
		exp_temp_min = exp_temp_min * 0.9; //Reduce the mininmum temperature, only if we will not end up with problems 
		}
	      	exp_temp_max = exp_temp_max * 1.1; //The maximum temperature can go up forever with no fear of numerical problems
	    }
	  if (isfinite(exp_temp_func (exp_temp_min))==0
	      || isfinite(exp_temp_func (exp_temp_max))==0)
	    {
	      Error
		("spectral_estimators: Exponential temperature cannot be bracketed (%e %e) in band %i - setting w to zero\n",
		 exp_temp_min,exp_temp_max,n);//NSH 131108 - now a warning, this should no longer happen
	      xplasma->exp_w[n] = 0.0;
	      xplasma->exp_temp[n] = -1e99;
	      expflag = -1;	//Discount an exponential model
	    }
	  else
	    {
/* We compute temporary values for sim alpha and sim weight. This will allow us to 
 * check that they are sensible before reassigning them */

/* But first see if we have a positive or negative solution. The temperatures are positive at the moment, if it was the negatives that worked, change the sign of the temperatures. */

	      if (exp_temp_func (-1.0*exp_temp_min) * exp_temp_func (-1.0*exp_temp_max) < 0.0)
		{
		exp_temp_min = -1.0* exp_temp_min;
		exp_temp_max = -1.0* exp_temp_max;
		}
//printf("trying to find correct temp between %e (%e) and %e (%e)\n",exp_temp_min,exp_temp_func(exp_temp_min),exp_temp_max,exp_temp_func(exp_temp_max));

	      exp_temp_temp = zbrent (exp_temp_func, exp_temp_min, exp_temp_max, 0.00001);	/* Solve for the effective temperature */
	      exp_w_temp = exp_w (j, exp_temp_temp, spec_numin, spec_numax);	/* Calculate the weight */


	      if ((isfinite (exp_w_temp))==0)
		{
		  Error
		    ("spectral_estimators: New exponential parameters (%e) unreasonable, using existing parameters. Check number of photons in this cell\n",exp_w_temp); //NSH 131108 - now a warning, this should no longer happen
		  expflag = -1;	//discount an exponential model
		  xplasma->exp_w[n] = 0.0;
		  xplasma->exp_temp[n] = -1e99;
		}
	      else
		{
		  xplasma->exp_temp[n] = exp_temp_temp;
		  xplasma->exp_w[n] = exp_w_temp;
		}
	    }

	  exp_sd = exp_stddev (xplasma->exp_temp[n], spec_numin, spec_numax);
	  pl_sd = pl_log_stddev (xplasma->pl_alpha[n], lspec_numin, lspec_numax);
	  Log_silent
	    ("NSH in this cell %i band %i PL estimators are log(w)=%10.2e, alpha=%5.3f giving sd=%e compared to %e\n",
	     xplasma->nplasma, n, xplasma->pl_log_w[n], xplasma->pl_alpha[n], pl_sd,
	     xplasma->xsd_freq[n]);
	  Log_silent
	    ("NSH in this cell %i band %i exp estimators are w=%10.2e, temp=%10.2e giving sd=%e compared to %e\n",
	     xplasma->nplasma, n, xplasma->exp_w[n], xplasma->exp_temp[n], exp_sd,
	     xplasma->xsd_freq[n]);
	  exp_sd =
	    fabs ((exp_sd - xplasma->xsd_freq[n]) / xplasma->xsd_freq[n]);
	  pl_sd =
	    fabs ((pl_sd - xplasma->xsd_freq[n]) / xplasma->xsd_freq[n]);
	  /* NSH 120817 These commands decide upon the best model, based upon how well the models predict the standard deviation */
	  if (expflag > 0 && plflag > 0)
	    {
	      if (exp_sd < pl_sd)
		xplasma->spec_mod_type[n] = SPEC_MOD_EXP;
	      else
		xplasma->spec_mod_type[n] = SPEC_MOD_PL;
	    }
	  else if (plflag > 0)
	    {
	      xplasma->spec_mod_type[n] = SPEC_MOD_PL;
	    }
	  else if (expflag > 0)
	    {
	      xplasma->spec_mod_type[n] = SPEC_MOD_EXP;
	    }
	  else
	    {
	      xplasma->spec_mod_type[n] = -1;	//Oh dear, there is no suitable model - this should be an error
	      Error ("No suitable model in band %i cell %i (nphot=%i fmin=%e fmax=%e)\n", n,
		     xplasma->nplasma,xplasma->nxtot[n],xplasma->fmin[n],xplasma->fmax[n]);
	  xplasma->fmin_mod[n]=spec_numax; //We will set the applicable frequency bands for the model to values that will cause errors if the model is used
	  xplasma->fmax_mod[n]=spec_numin; 
	    }
	  Log_silent ("NSH In cell %i, band %i, the best model is %i\n",
		      xplasma->nplasma, n, xplasma->spec_mod_type[n]);
	}			//End of loop that does things if there are more than zero photons in the band.
    }				//End of loop over bands 



  return (0);

}



/***********************************************************
                                       Space Telescope Science Institute

 Synopsis:
 	integrand used with zbrent to calulate alpha from
	a mean frequency and a minimum and maximum frequecy
	for a band
        The function is the integral of nu j(nu) divided by the integral of j(nu).
	The prefector w is lost in the division, so this function should equal the
	mean frequency. Zbrent searches for the value of alpha which makes this
	function equal to the mean frequency computed for the cell from photon 
	packets. So we subtract the mean frequency at the end to make the value of the
	function equal to zero at the correct alpha.
	
 Arguments:		

Returns:
 
Description:	
	
Notes:

History:  	NSH 120817 - 	Tidied up to make it a bit more clear what is going on here. There are
				now two functions, one is pl_mean, which is used elswhere, and ome is pl_alpha_func which
				is very simple, just works out the difference between the computed mean, and the measured
				mean
		NSH 131108 - 	changed into a log formulation - and changed the name.


**************************************************************/

double
pl_alpha_func_log (alpha)
     double alpha;
{
  double answer;

//  answer = pl_mean (alpha, spec_numin, spec_numax) - spec_numean;
  answer = pl_logmean (alpha, lspec_numin, lspec_numax) - spec_numean; //NSH 131106 change to deal with large number issues
  return (answer);
}

/* Function below is the old non-log formulation 
double
pl_alpha_func (alpha)
     double alpha;
{
  double answer;

  answer = pl_mean (alpha, spec_numin, spec_numax) - spec_numean;
  return (answer);
}
*/

/* This is the old, non log function to compute mean 

double
pl_mean (alpha, numin, numax)
     double alpha;
     double numin, numax;
{
  double answer, numerator, denominator;

  numerator = (pow (numax, (alpha + 2.)) - (pow (numin, (alpha + 2.)))) / (alpha + 2.);	
  denominator = (pow (numax, (alpha + 1.)) - (pow (numin, (alpha + 1.)))) / (alpha + 1.);	
  answer = numerator / denominator;
  return (answer);
}
*/

/*NSH 131108 - changed the formulation to work in log space - also changed name */   

double
pl_logmean (alpha, lnumin, lnumax)
	double alpha;
	double lnumin,lnumax;
{
  double k, answer, numerator, denominator, a,b,c,d;

  k = -1.0 * (lnumax * (alpha + 2.0)); //This is a prefector that scales everything to sensible numbers.

//printf ("lnumax=%f, lnumin=%f, k=%f\n",lnumax,lnumin,k);

a=pow (10.0,(k+(lnumax*(alpha+2.0))));
b=pow (10.0,(k+(lnumin*(alpha+2.0))));
c=pow (10.0,(k+(lnumax*(alpha+1.0))));
d=pow (10.0,(k+(lnumin*(alpha+1.0))));

//printf ("a=%e,b=%e,c=%e,d=%e\n",a,b,c,d);

  numerator = (a - b) / (alpha + 2.0);

  denominator = (c - d) / (alpha + 1.0);

  answer = numerator / denominator;


  return (answer);
}





/**************************************************************************
                    Space Telescope Science Institute


  Synopsis:  

	This program works out the scaling factor for a power law representation of a spectrum in 
	a cell. Assuming alpha is known, it computes w for E=w x integral from numin to numax of nu^alpha.
	en1 is the band limited mean radiation density. 
        

  Description:	

  Arguments:  (Input via .pf file)		


  Returns:

  Notes:


 

  History:
111227	ksl	This routine was written by Nick in 2011, but has very
		little documentations or explanation
120817  nsh     routine moved into its natural home where it is used to 
		compute the prefactor for a power law estimation of a
		spectrum in a cell. 
131107	nsh	changed to work in the new log space. Also changed the 
		name

 ************************************************************************/

/* The old - non log function */
//double
//pl_w (j, alpha, numin, numax)
//     double j;			//the band limited spectral density
//     double alpha;		//Computed spectral index for the cell
//     double numin, numax;	//Range of frequencies we are considering
//{
//  double w;			//the answer
//
//  double integral;		/*This will hold the unscaled integral of jnu from numin to numax */


//  integral = (pow (numax, (alpha + 1.0)) - pow (numin, (alpha + 1.0))) / (alpha + 1.0);	/* This is the integral */
//  w = j / integral;		/*THe scaling factor is the actual band limited J / the integral. */

//  return (w);
//}

double
pl_log_w (j, alpha, lnumin, lnumax)
     double j;			//the band limited spectral density
     double alpha;		//Computed spectral index for the cell
     double lnumin, lnumax;	//Range of frequencies we are considering
{
  double logw;			//the answer
  double logk; //scaling prefactor to permit huge numbers to be dealt with
  double log_integral;		/*This will hold the unscaled integral of jnu from numin to numax */

  logk=-1.0 * lnumax * (alpha + 1.0);

  log_integral = log10 (( pow(10,(logk+(alpha+1.0)*lnumax)) - pow(10,(logk+(alpha+1.0)*lnumin)) )/(alpha + 1.0));

  logw = log10(j) + logk - log_integral; 

  return (logw);
}


/**************************************************************************
                    Southampton University


  Synopsis:  

	This program works out the standard deviation of a photon distribution modelled
	by a power law. First we work out the integral of nu^2 jnu / int jnu.
	Then we subtract the (mean/j)^2 and then take the square root.
        

  Description:	

  Arguments:  (Input via .pf file)		


  Returns:

  Notes:


 

  History:
120817  nsh     coded
131107	nsh	Changed to work in log space

 ************************************************************************/

//This is the original non -log foormulation 
//double
//pl_stddev (alpha, numin, numax)
//     double alpha;		//Computed spectral index for the cell
//     double numin, numax;	//Range of frequencies we are considering
//{
//  double answer;		//the answer

//  double numerator;		/*This will hold the unscaled integral of nu^2 jnu from numin to numax */
//  double denominator;		/*This is the unscaled integral of jnu from numin to numax */
//  double mean;

//  numerator = (pow (numax, (alpha + 3.0)) - pow (numin, (alpha + 3.0))) / (alpha + 3.0);	/* This is the integral for nu^2 bnu */
//  denominator = (pow (numax, (alpha + 1.0)) - pow (numin, (alpha + 1.0))) / (alpha + 1.0);	/* This is the integral for nbu */
//  mean = pl_mean (alpha, numin, numax);
//  answer = sqrt ((numerator / denominator) - (mean * mean));


//  return (answer);
//}


double
pl_log_stddev (alpha, lnumin, lnumax)
     double alpha;		//Computed spectral index for the cell
     double lnumin, lnumax;	//Range of frequencies we are considering
{
  double answer;		//the answer

  double numerator;		/*This will hold the unscaled integral of nu^2 jnu from numin to numax */
  double denominator;		/*This is the unscaled integral of jnu from numin to numax */
  double mean;
  double k, a,b,c,d;

  k = -1.0 * (lnumax * (alpha + 3.0)); //This is a prefector that scales everything to sensible numbers.

//printf ("lnumax=%f, lnumin=%f, k=%f\n",lnumax,lnumin,k);

a=pow (10.0,(k+(lnumax*(alpha+3.0))));
b=pow (10.0,(k+(lnumin*(alpha+3.0))));
c=pow (10.0,(k+(lnumax*(alpha+1.0))));
d=pow (10.0,(k+(lnumin*(alpha+1.0))));

//printf ("a=%e,b=%e,c=%e,d=%e\n",a,b,c,d);

  numerator = (a - b) / (alpha + 3.0);

  denominator = (c - d) / (alpha + 1.0);

  answer = numerator / denominator;

  mean = pl_logmean (alpha, lnumin, lnumax);
  answer = sqrt ((numerator / denominator) - (mean * mean));


  return (answer);
}



/***********************************************************
                                       Southampton University

 Synopsis:
 	integrand used with zbrent to calulate effective temperature for an
	exponential model from a mean frequency and a minimum and maximum frequecy
	for a band
        The mean is the integral of nu j(nu) divided by the integral of j(nu).
	The prefector w is lost in the division, so this function should equal the
	mean frequency. Zbrent searches for the value of 'temperature' which makes this
	function equal to the mean frequency computed for the cell from photon 
	packets. So we subtract the mean frequency at the end to make the value of the
	function equal to zero at the correct temperature.
	
 Arguments:		

Returns:
 
Description:	
	
Notes:

History:  NSH 120817 - Tidied up to make it a bit more clear what is going on here. There are
	now two functions, one is pl_mean, which is used elswhere, and one is pl_alpha_func which
	is very simple, just works out the difference between the computed mean, and the measured
	mean


**************************************************************/

double
exp_temp_func (exp_temp)
     double exp_temp;
{
  double answer;

  answer = exp_mean (exp_temp, spec_numin, spec_numax) - spec_numean;
  return (answer);
}

double
exp_mean (exp_temp, numin, numax)
     double exp_temp;
     double numin, numax;
{
  double answer, numerator, denominator;
  double exp1;			/* We supply a temperature, but actually we expect the correct function to be of the form e^-hnu/kt, so this will hold -1*h/kt */
  double emin, emax;		/*The exponential evaluated at numin and numax */
  double t1, t2;

  exp1 = (-1.0 * H) / (BOLTZMANN * exp_temp);

  t1 = exp1 * numin;
  t2 = exp1 * numax;


  emin = exp (t1);
  emax = exp (t2);


  numerator = ((emax / exp1) * ((exp1 * numax) - 1.0)) - ((emin / exp1) * ((exp1 * numin) - 1.0));	/*THe integral of nu e^(-hnu/kt) from numin to numax */


  denominator = emax - emin;	/*NSH 120817 This is the integral of (J(nu)=e^(-hnu/jkt)) */
  answer = numerator / denominator;	/*NB there is a factor of exp1 top and bottom which I have cancelled out */
  return (answer);
}


/**************************************************************************
                    Southampton University


  Synopsis:  

	This program works out the scaling factor for an exponential representation of a spectrum in 
	a cell. Assuming the effective tremperature is known, it computes w for E=w x integral from numin to numax of e^(-nhu/kt).

        

  Description:	

  Arguments:  (Input via .pf file)		


  Returns:

  Notes:


 

  History:

120817  nsh     coded 

 ************************************************************************/


double
exp_w (j, exp_temp, numin, numax)
     double j;			//the band limited spectral density
     double exp_temp;		//Computed effective temperature for the cell
     double numin, numax;	//Range of frequencies we are considering
{
  double w;			//the answer

  double integral;		/*This will hold the unscaled integral of jnu from numin to numax */
  double exp1;			/* We supply a temperature, but actually we expect the correct function to be of the form e^-hnu/kt, so this will hold -1*h/kt */
  double emin, emax;		/*The exponential evaluated at numin and numax */

  exp1 = (-1.0 * H) / (BOLTZMANN * exp_temp);


  emin = exp (exp1 * numin);
  emax = exp (exp1 * numax);
  integral = (emax - emin) / exp1;	/* This is the integral */
  w = j / integral;		/*THe scaling factor is the actual band limited J / the integral. */

  return (w);
}




/**************************************************************************
                    Southampton University


  Synopsis:  

	This program works out the standard deviation of a photon distribution modelled
	by an exponential. First we work out the integral of nu^2 jnu / int jnu.
	Then we subtract the (mean/j)^2 and then take the square root.
        

  Description:	

  Arguments:  	


  Returns:

  Notes:


 

  History:
120817  nsh     coded

 ************************************************************************/


double
exp_stddev (exp_temp, numin, numax)
     double exp_temp;		//Computed spectral index for the cell
     double numin, numax;	//Range of frequencies we are considering
{
  double answer;		//the answer
  double exp1;			/* We supply a temperature, but actually we expect the correct function to be of the form e^-hnu/kt, so this will hold -1*h/kt */
  double emin, emax;		/*The exponential evaluated at numin and numax */
  double nmax, nmin;		/*The integral is a bit complicated, so we have split it up into value for numin and numax */
  double numerator, denominator;
  double mean;

  exp1 = (-1.0 * H) / (BOLTZMANN * exp_temp);
  emin = exp (exp1 * numin);
  emax = exp (exp1 * numax);

  nmax = emax * (((numax * numax) / (exp1)) - ((2 * numax) / (exp1 * exp1)) + ((2.0) / (exp1 * exp1 * exp1)));	/* Is the integral for nu^2 bnu exaluated for numax */
  nmin = emin * (((numin * numin) / (exp1)) - ((2 * numin) / (exp1 * exp1)) + ((2.0) / (exp1 * exp1 * exp1)));	/* Is the integral for nu^2 bnu exaluated for numin */
  numerator = nmax - nmin;	/*The definite integral */
  denominator = (emax - emin) / exp1;	/* This is the integral for nbu */
  mean = exp_mean (exp_temp, numin, numax);
  answer = sqrt ((numerator / denominator) - (mean * mean));


  return (answer);
}



