

/* This file was created to contain all of the routines associated
	with the attempt to calculate ionic abundances using a
	temperature fixed for each pair of ions to make the uncorrected
	abundance ration roughly equal to 1. There are then correction
	factors applied, either based upon a power law model of
	the true radiation field, or as a dilute blackbody.

	12Feb	nsh	Began coding 
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>


#include "atomic.h"
#include "python.h"



//struct photoionization *xver; //Verner & Ferland description of a photoionization x-section NSH 13Sep - should not be used any more
struct topbase_phot *xtop;	//Topbase description of a photoionization x-section 
PlasmaPtr xxxplasma;
double qromb_temp;		// This is a storage variable for the current electron temperature so it is available for qromb calls
int niterate;			// Make this a variable that all the subroutines cn see, 
				// so we can decide if we need to recompute the numerators




/***********************************************************
                                       Southampton University

  Synopsis:   

int
variable_temperature (xplasama, mode)  modifies the densities of ions, levels, and
	partition functions of ions within a cell of the wind based upon the mode, 
	and other data contained within the WindPtr itself, such as t_r, t_e, w, 
	based on the "mode".  

        Unlike nebular_concentrations, this code works on pairs of ions at a time, 
	and uses a temperature calculated to be suitable for that particular pair
	to have a ratio of abundances about equal to 1.
  
  Arguments:		
     PlasmaPtr ww;
     int mode;			// 6=correct using a dilute blackbody
     				   7=correct using a broken power law spectrum


  Returns:
 	variable temperature alters portions of the wind ptr.  Exactly how things 
	are changed depends on the mode.
 
 	variable_temperature returns 0 if it converged to an answer, -1 otherwise.  On an
 	abnormal return the density array and ne are left unchanged.
 	
  Description:	


 
  Notes:

  Section 5.3 of Nick's thesis explains the undelying ratinale for the so-called 
  variable temperature models.  The name of the routine is misleading since 
  there correctios are made for both a dilute power law and for a dilute bb. 



  History:
	2012Feb	nsh	Coded and debugged as part of QSO effort. 
        1212Dec nsh	Recoded so that the densities are computed in a 
			temporary array, and only 
			committted to the real density structure once 
			we are sure the code converges.
	2013Sep nsh	ground state fudge computed at the start, and 
			stored in an array rather
			than contiunually recomputing it - sometimes it is expensive.
	

**************************************************************/



double xxxne, xip;

int
variable_temperature (xplasma, mode)
     PlasmaPtr xplasma;
     int mode;			//   6=correct using dilute blackbody, 7=power law
{
  int nion;			
  double xnew, xsaha;
  double theta, x;
  double get_ne ();
  double t_e, t_r, xtemp, nh, xne, xxne, www;
  double a, b;
  double newden[NIONS];		//NSH 121217 - recoded so density is computed in a temperary array
  // double xip; //ionzation potential of lower ion.
  int nelem, first, last;
//  double t_r_part_correct, NSH 130401 - this variable no longer used in 75. It may come back tho!
  double t_e_part_correct;
  double sum, big;
  double pi_fudge, recomb_fudge, tot_fudge;	/*Two of the correction factors for photoionization rate, and recombination rate */
  double gs_fudge[NIONS];	/*It can be expensive to calculate this, and it only depends on t_e - which is fixed for a run. So 
				   //                 calculate it once, and store it in a temporary array */

  /* Copy xplasma to local plasma varaible, used to communicate w and alpha to the power law correction routine. 
     NSH 120703 - also used to set the  denominator calculated for a given ion for this cell last time round 
     ksl 1407 - Not obvious why this was done; the next line merely sets xxxplasma to the pointer xplasma.  It
     does not copy anythin*/

  xxxplasma = xplasma;		


  nh = xplasma->rho * rho2nh;	//LTE
  t_e = xplasma->t_e;
  t_r = xplasma->t_r;
  www = xplasma->w;


  /* Copy the current densities into the temporary array */

  for (nion = 0; nion < nions; nion++)
    {
      newden[nion] = xplasma->density[nion]; 

      /* Here we populate the recombination to ground state correction factor used in the LM and Sim ionization
       * equations.  Mode 2 imples we are include the dielectronic correction. There is no recombination fudge for 
       * the neutral ion */

      if (ion[nion].istate != 1)	
	{
	  gs_fudge[nion] = compute_zeta (t_e, nion - 1, 2);	//We call with nion-1, so when we address the array, we will want to ask for nion
	}
    }



  /* make an initial estimate of ne based on H alone,  Our guess
     assumes ion[0] is H1.  Note that x below is the fractional
     ionization of H and should vary from 0 to 1

     Note that in a pure H plasma, the left hand side of the Saha 
     equation can be written as (x*x*nh*nh)/((1-x)*nh)  and can
     then be converted into a quadratic.  That is what is done
     below.  

     Since this is an initial estimate, g factors are ignored
   *
   */

  if (t_e < MIN_TEMP)
    t_e = MIN_TEMP;		/* fudge to prevent divide by zeros */

  xsaha = SAHA * pow (t_e, 1.5);

  theta = xsaha * exp (-ion[0].ip / (BOLTZMANN * t_e)) / nh;

  if (theta < THETAMAX)
    {
      x = (-theta + sqrt (theta * theta + 4 * theta)) / 2.;
      xne = xxne = xxxne = x * nh;
    }
  else
    xne = xxne = xxxne = nh;	/*xxne is just a store so the error can report the starting value of ne. 
                                  xxxne is the shared variable so the temperature solver routine can access it */

  if (xne < 1.e-6)
    xne = xxxne = 1.e-6;	/* Set a minimum ne to assure we can calculate
				   xne the first time through the loop */


  /* At this point we have an initial estimate of ne. */

  niterate = 0;
  while (niterate < MAXITERATIONS)
    {
      /* We loop over all elements */
      for (nelem = 0; nelem < nelements; nelem++)
	{

	  first = ele[nelem].firstion;	/*first and last identify the position in the array */
	  last = first + ele[nelem].nions;	/*  So for H which has 2 ions, H1 and H2, first will 
             					    generally be 0 and last will be 2 so the for loop 
                                                    below will just be done once for nion = 1 */

	  if (first < 0 || first >= nions || last < 0 || last > nions)
	    {
	      Error
		("variable_temperature: Confusion for element %d with first ion %d and last ion %d\n",
		 nelem, first, last);
	      exit (0);
	    }


	  /* and now we loop over all the ions of this element */
	  sum = newden[first] = 1.0;	/* set the density of the first ion of this element to 1.0 
                                           -this is (always??) the neutral ion */

	  big = pow (10., 250. / (last - first));	/*make sure that we will not overflow the last ion */


	  for (nion = first + 1; nion < last; nion++)	/*nion is the upper ion of the pair, 
                                                         last is one above the actual last ion, 
                                                         so the last ion to be considered is last-1 */
	    {

	      tot_fudge = 0.0;	/* NSH 130605 to remove o3 compile error */

	      /* now we need to work out the correct temperature to use */
	      xip = ion[nion - 1].ip;	//the IP is that from the lower to the upper of the pair
	      xtemp = zbrent (temp_func, MIN_TEMP, 1e8, 10);	//work out correct temperature


              /* given this temperature, we need the pair of partition functions for these ions */
	      partition_functions_2 (xplasma, nion, xtemp, 1);	//weight of 1 give us the LTE populations.


             /* and now we need the saha equation linking the two states at our chosen temp 
                NB the partition function of the electron is included in the SAHA factor */
	      xsaha = SAHA * pow (xtemp, 1.5);
	      b = xsaha * xplasma->partition[nion]
		* exp (-ion[nion - 1].ip / (BOLTZMANN * xtemp)) / (xne *
								   xplasma->partition
								   [nion -
								    1]);


	      t_e_part_correct =
		xplasma->partition[nion - 1] / xplasma->partition[nion];
	      partition_functions_2 (xplasma, nion, t_e, 0);	/* Our only real guess here is that the 
                                                                   electron temperature might give a good 
                                                                   estimate of the partition function */

	      t_e_part_correct *=
		(xplasma->partition[nion] / xplasma->partition[nion - 1]);


	      /* we now correct b to take account of the temperature and photon field 
		 t_r and www give the actual radiation field in the cell, xtemp is the temp we used
  		 t_e is the actual electron temperature of the cell*/

	      recomb_fudge = sqrt (t_e / xtemp);

              /* Correct the SAHA equation abundance pair using an actual radiation field modelled as a dilute black body */
	      if (mode == 6)	
		{

		  if (t_r / xtemp < 0.2)	
                 /* If the true radiation temperature is much lower than the temperature at which 
                    the ion is expected to exist, we wont see it, so save time and dont bother 
                    calculating correction factors */
		    {
		      tot_fudge = 0.0;
		    }

		  else
		    {

		      pi_fudge = bb_correct_2 (xtemp, t_r, www, nion);

		      tot_fudge =
			pi_fudge * recomb_fudge * (gs_fudge[nion] +
						   www * (1 -
							  gs_fudge[nion]));
		    }
		}

              /* correct the SAHA equation abundance pair using an actual radiation field modelled as a broken power law */

	      else if (mode == 7)	
		{
		  pi_fudge = pl_correct_2 (xtemp, nion);


		  tot_fudge =
		    pi_fudge * recomb_fudge * gs_fudge[nion] *
		    t_e_part_correct;
		}

	      else
		{
		  Error ("variable_temperature: unknown mode %d\n", mode);
		  exit (0);
		}


	      /* apply correction factors */
	      b *= tot_fudge;
	      if (b > big)
		b = big;	//limit step so there is no chance of overflow
	      a = newden[nion - 1] * b;
	      sum += newden[nion] = a;

	      if (newden[nion] < 0.0)
		{
		  Error
		    ("variable_temperature: ion %i has a negative density %e\n",
		     nion, newden[nion]);
		  mytrap ();
		}
	      if (sane_check (sum))
		Error
		  ("variable_temperature:sane check failed for running density sum=%e, last ion=%i\n",
		   sum, nion);
	    }			//end of loop over ions - we now have a full set of ions for this element



	  a = nh * ele[nelem].abun / sum;	//the scaling factor to get the overall abundance right
	  for (nion = first; nion < last; nion++)
	    {
	      newden[nion] *= a;	//apply scaling
	      if (sane_check (newden[nion]))	//check nothing has gone crazy
		Error
		  ("variable_temperature:sane check failed for density newden=%e, for ion=%i\n",
		   newden[nion], nion);
	    }



         /* Re solve for the macro atom populations with the current guess for ne */
         /* JM1308 -- note that unlike theabove, here we actually modify the xplasma
	    structure for those ions which are being treated as macro ions. This means that the
	    the newden array will contain wrong values for these particular macro ions, but due
            to the if loop at the end of this subroutine they are never passed to xplasma */

	  if (geo.macro_ioniz_mode == 1)
	    {
	      macro_pops (xplasma, xne);
	    }

	  for (nion = 0; nion < nions; nion++)
	    {
             
               /* if the ion is being treated by macro_pops then use the populations just computed 
                  JM1309 -- this was missing prior to python 76c */
	      if ((ion[nion].macro_info == 1) && (geo.macro_simple == 0)
		  && (geo.macro_ioniz_mode == 1))
		{
		  newden[nion] = xplasma->density[nion];
		}

              /*Set some floor so future divisions are sensible */
	      if (newden[nion] < DENSITY_MIN)
		newden[nion] = DENSITY_MIN;
	    }

	  /* Now determine the new value of ne from the ion abundances */
	}			//end of loop over elements

      xnew = get_ne (newden);	/* determine the electron density for this density distribution */


      if (xnew < DENSITY_MIN)
	xnew = DENSITY_MIN;	/* fudge to keep a floor on ne */
      if (fabs ((xne - xnew) / (xnew)) < FRACTIONAL_ERROR || xnew < 1.e-6)
	{
	  break;
	}
      xne = xxxne = (xnew + xne) / 2.;	/*New value of ne */
      niterate++;


      if (niterate == MAXITERATIONS)
	{
	  Error
	    ("variable_temperature: failed to converge t %.2g nh %.2g xnew %.2g\n",
	     t_e, nh, xnew);
	  Error ("variable_temperature: xxne %e theta %e\n", xxne, theta);
	  return (-1);
	}
    }


/* Finally transfer the calculated densities to the real density array. This is only called if the code
 has iterated correctly, if not then the real densities will stay the same as last time*/

  xplasma->ne = xnew;
  for (nion = 0; nion < nions; nion++)
    {
      /* If statement added here to suppress interference with macro populations (SS Apr 04) */
      if (ion[nion].macro_info == 0 || geo.macro_ioniz_mode == 0
	  || geo.macro_simple == 1)
	{
	  xplasma->density[nion] = newden[nion];
	}
    }





  partition_functions (xplasma, 4);	/*WARNING fudge NSH 11/5/14 - this is as a test. 
                                         We really need a better implementation of partition functions and levels
                                        for a power law illuminating spectrum. We found that if we didnt make this call, 
                                        we would end up with undefined levels - which did really crazy things */

  return (0);
}





/***********************************************************
                                       Southampton University

  Synopsis:   
	bb_correct_2(xtemp, t_r,www,nion) calculates the 
 	correction factor of the number density of two adjacent levels calculated
	at any temperature. LM is a special case of this where the temperatures
	are the same. 
 
  Arguments:		

	double xtemp		the temperature at which the saha abundances were calculated
	double t_r		temperature of the radiation field 
	double www		dilution factor for the radiation field
	int nion;		number of the upper ion in the pair for which the correction factor is
				being calculated.

  Returns:
 	j			the correction factor
 
 
 	
 Description:	


 
  Notes:



  History:
 
	12feb		NSH	Coded as part of the quasar project
	13sep	NSH	Changed to always use the topbase style of data - verner data is now tabulated.

**************************************************************/
int bb_correct_err = 0;

double
bb_correct_2 (xtemp, t_r, www, nion)
     double t_r, xtemp, www;
     int nion;
{
  double q;
  int ion_lower, n;
  double numerator, denominator;
  int ntmin, nvmin;
  double fthresh, fmax, fmaxtemp;



  /* Initialization of fudges complete */

  ion_lower = nion - 1;		/*the call uses nion, which is the upper ion in the pair */


  if (-1 < ion_lower && ion_lower < nions)	//Get cross section for this specific ion_number
    {
      ntmin = ion[ion_lower].ntop_ground;
      nvmin = ion_lower;
    }
  else
    {
      Error ("bb_correct_2: %d is unacceptable value of nion\n", ion_lower);
      mytrap ();
      exit (0);
      return (1.0);
    }

/*The integration is to correct for photoionization rates from the lower to the upper, so
 all integrals have to be over the lower ion cross section. The integrals use the code in
 power_sub, so the temperature has to be set to an external variable. Annoyingly this is 
 called sim_te....*/

  if (ion[ion_lower].phot_info == 1)	//topbase
    {
      n = ntmin;
      xtop = &phot_top[n];	//use the topbase data
    }
  else if (ion[ion_lower].phot_info == 0)	// verner
    {
      n = nvmin;		//just the ground state ionization fraction.
      xtop = &xphot_tab[ion[n].nxphot];	//Use the tabuled verner data
    }
  else				//no data available
    {
      bb_correct_err++;
      /* 71 - 111229 - ksl - Suppress this error after 100 instances so program does not bomb */
      if (bb_correct_err < 100)
	{
	  Error
	    ("bb_correct: No photoionization xsections for ion %d (element %d, ion state %d),setting sim_frac to 1\n",
	     nion, ion[nion].z, ion[nion].istate);
	}
      else if (bb_correct_err == 100)
	{
	  Error
	    ("xinteg_sim: Suppressing photoionization xsection error, but more photo xsections are needed in atomic data\n");
	}


      q = 1.;
      return (q);
    }
  fthresh = xtop->freq[0];
  fmaxtemp = xtop->freq[xtop->np - 1];
  fmax = check_fmax (fthresh, fmaxtemp, xtemp);
  if (fthresh > fmax)
    {
      Error
	("bb_correct: After checking, fthresh has been set below fmin - we cannot compute denominator\n");
      q = 1.0;
      return (q);
    }

  qromb_temp = t_r;		//The numerator is for the actual radiation temperature

  if (pow
      (((t_r -
	 xxxplasma->PWntemp[ion_lower]) / xxxplasma->PWntemp[ion_lower]),
       2) > 1e-6)
    {
      numerator = www * qromb (tb_planck1, fthresh, fmax, 1.e-4);	//and is corrected for W
      xxxplasma->PWnumer[ion_lower] = numerator;
      xxxplasma->PWntemp[ion_lower] = t_r;
    }
  else
    {
      numerator = xxxplasma->PWnumer[ion_lower];
    }

  if (numerator == 0.0)
    return (0.0);		//There is no need to waste time evaluating the denominator


  qromb_temp = xtemp;		//The denominator is calculated for the LTE rate at our ideal temp

//      Log ("topbase n=%i,nion=%i,temp=%e,fthresh=%e,fmax=%e,hnu=%e,hnu/kt=%e,fmax=%e\n",n,ion_lower,temp,fthresh,fmax,fthresh*H,(fthresh*H)/(temp*BOLTZMANN),5.879e10*temp);



  if (pow (((xtemp - xxxplasma->PWdtemp[ion_lower]) / xxxplasma->PWdtemp[ion_lower]), 2) > 1e-6)	//Only do this integral if there is a significant change in xtemp
    {
      denominator = qromb (tb_planck1, fthresh, fmax, 1.e-4);
      xxxplasma->PWdenom[ion_lower] = denominator;
      xxxplasma->PWdtemp[ion_lower] = xtemp;	//Store the result
    }
  else
    {
      denominator = xxxplasma->PWdenom[ion_lower];
    }


  q = numerator / denominator;


  return (q);
}

/* This next little bit of code was used to calculate an ideal temperature. This gets a temperature which gives a ratio of n to n+1 of close to 1, but at this temperature, the blackbody function used for the denominator of the correction factor is almost always zero. I've left the routine in, but it doesnt seem to work...*/



double
temp_func (solv_temp)
     double solv_temp;
{
  double answer;
//      Log ("xxxne=%e, xip=%e\n",xxxne,xip);
//      Log ("1=%e,2=%e,3=%e\n",log(4.83e15/xxxne),1.5*log(temp),(xip/(BOLTZMANN*temp)));
  answer =
    log (4.83e15 / xxxne) + 1.5 * log (solv_temp) -
    (xip / (BOLTZMANN * solv_temp));

  return (answer);
}


/***********************************************************
                                       Southampton University

  Synopsis:   
	pl_correct_2(xtemp,nion) calculates the 
 	correction factor of the number density of two adjacent levels calculated
	at any temperature based upon a true radiation field modelled by a broken power law
 
  Arguments:		

	double xtemp		the temperature at which the saha abundances were calculated
	int nion;		the upper ion in an ion pair for which we want to correction factyor

  Returns:
 	q		        the correction factor
 	
 Description:	

	 

 
  Notes:
	The weights and alpha values needed for the integrals are communicated to this routine vias xxxplasma, an external varaible which is popullated in the variable_temperature routine.
	The integrals over power law require individual alphas and w, these are communicated via external variables xpl_alpha and xpl_w.	


  History:
 
	12feb	NSH	Coded as part of the quasar project
        12aug   NSH     Added code to allow an exponential model of j to be used.
	13sep	NSH	Reversed commenting out of lines that check for maximum frequency and 
			reset integration limits - unsure why this was ever commented out!
	13sep	NSH	Changed to always use the topbase style of data - verner data is now tabulated.
	13nov	NSH	Changed to use the new log version of PL, also the min and max frequencies in 
			a band are used to decide the applicability of a model.

**************************************************************/
int pl_correct_err = 0;
double xpl_alpha, xpl_w, xpl_logw;
double xexp_temp, xexp_w;


double
pl_correct_2 (xtemp, nion)
     double xtemp;
     int nion;
{
  double q;
  int ion_lower, n, j;
  double numerator, denominator;
  int ntmin, nvmin;
  double fthresh, fmax, fmaxtemp;
  double f1, f2;
  double exp_qromb, pl_qromb;


  exp_qromb = 1e-4;
  pl_qromb = 1e-4;

  ion_lower = nion - 1;		/*the call uses nion, which is the upper ion in the pair, we are interested in the PI from the lower ion */

  if (-1 < ion_lower && ion_lower < nions)	//Get cross section for this specific ion_number
    {
      ntmin = ion[ion_lower].ntop_ground;	/*We only ever use the ground state cross sections. This is for topbase */
      nvmin = ion_lower;	/*and this is for verner cross sections */
    }
  else
    {
      Error ("bb_correct_2: %d is unacceptable value of nion\n", ion_lower);
      mytrap ();
      exit (0);
      return (1.0);
    }


/*The integration is to correct for photoionization rates from the lower to the upper, so
 all integrals have to be over the lower ion cross section. The integrals use external functions, so the temperature has to be set to an external variable temp. */

  if (ion[ion_lower].phot_info == 1)	//topbase
    {
      n = ntmin;
      xtop = &phot_top[n];
    }
  else if (ion[ion_lower].phot_info == 0)	// verner
    {
      n = nvmin;		//just the ground state ionization fraction.
      xtop = &xphot_tab[ion[n].nxphot];
    }
  else
    {
      pl_correct_err++;		/* If we get here, there are no cross sections available */
      if (pl_correct_err < 100)
	{
	  Error
	    ("pl_correct: No photoionization xsections for ion %d (element %d, ion state %d),setting sim_frac to 1\n",
	     nion, ion[nion].z, ion[nion].istate);
	}
      else if (pl_correct_err == 100)
	{
	  Error
	    ("xinteg_sim: Suppressing photoionization xsection error, but more photo xsections are needed in atomic data\n");
	}
      q = 1.;			/*This is really bad actually, this will leave the abundances all wrong.... */
      return (q);
    }
  fthresh = xtop->freq[0];
  fmax = xtop->freq[xtop->np - 1];
  numerator = 0;
  if (niterate == 0)		//first time of iterating this cycle, so calculate the numerator
    {
      for (j = 0; j < geo.nxfreq; j++)	//We loop over all the bands
	{
	  xpl_alpha = xxxplasma->pl_alpha[j];
//	  xpl_w = xxxplasma->pl_w[j];
	  xpl_logw = xxxplasma->pl_log_w[j];
	  xexp_temp = xxxplasma->exp_temp[j];
	  xexp_w = xxxplasma->exp_w[j];
	  if (xxxplasma->spec_mod_type[j] > 0)	//Only bother doing the integrals if we have a model in this band
	    {
//	      f1 = geo.xfreq[j];
	      f1 = xxxplasma->fmin_mod[j]; //NSH 131114 - Set the low frequency limit to the lowest frequency that the model applies to
//	      if (geo.xfreq[j + 1] > xxxplasma->max_freq && geo.xfreq[j] < xxxplasma->max_freq)	//The maximum frequency seen in this cell is in this band, so we cannot safely use the power law estimators right up to the top of the band. Note, we hope that all the weights in bands above this will be zero!
/* NSH 130909 - this loop was commented out at some point in the past - I cannt recall why, but the lines are now reinstated - if they cause a problem - refer back to issue #50! */
//		{
//		  f2 = xxxplasma->max_freq;
//		}
//	      else
//		{
//		  f2 = geo.xfreq[j + 1];	//We can safely integrate over the whole band using the estimators for the cell/band 
		  f2= xxxplasma->fmax_mod[j]; //NSH 131114 - Set the high frequency limit to the highest frequency that the model applies to
//		}
	      if (f1 < fthresh && fthresh < f2 && f1 < fmax && fmax < f2)	//Case 1- 
		{
		  if (xxxplasma->spec_mod_type[j] == SPEC_MOD_PL)
		    {
//		      numerator += qromb (tb_pow1, fthresh, fmax, pl_qromb);
		      numerator += qromb (tb_logpow1, fthresh, fmax, pl_qromb);
		    }
		  else
		    {
		      numerator += qromb (tb_exp1, fthresh, fmax, exp_qromb);
		    }
		}
	      else if (f1 < fthresh && fthresh < f2 && f2 < fmax)	//case 2 
		{
		  if (xxxplasma->spec_mod_type[j] == SPEC_MOD_PL)
		    {
//		      numerator += qromb (tb_pow1, fthresh, f2, pl_qromb);		      
		      numerator += qromb (tb_logpow1, fthresh, f2, pl_qromb);
		    }
		  else
		    {
		      numerator += qromb (tb_exp1, fthresh, f2, exp_qromb);
		    }
		}
	      else if (f1 > fthresh && f1 < fmax && fmax < f2)	//case 3
		{
		  if (xxxplasma->spec_mod_type[j] == SPEC_MOD_PL)
		    {
//		      numerator += qromb (tb_pow1, f1, fmax, pl_qromb);
		      numerator += qromb (tb_logpow1, f1, fmax, pl_qromb);
		    }
		  else
		    {
		      numerator += qromb (tb_exp1, f1, fmax, exp_qromb);
		    }
		}
	      else if (f1 > fthresh && f2 < fmax)	// case 4
		{
		  if (xxxplasma->spec_mod_type[j] == SPEC_MOD_PL)
		    {
//		      numerator += qromb (tb_pow1, f1, f2, pl_qromb);
		      numerator += qromb (tb_logpow1, f1, f2, pl_qromb);
		    }
		  else
		    {
		      numerator += qromb (tb_exp1, f1, f2, exp_qromb);
		    }
		}
	      else		//case 5 - should only be the case where the band is outside the range for the integral.
		{
		  numerator += 0;	// Add nothing - bit of a null statement, but makes the code look nice.
		}
	    }			//End of loop to only integrate in this band if there is power
	}			//End of loop over all bands, at this point we have the numerator
      xxxplasma->PWnumer[ion_lower] = numerator;	// Store the calculated numerator for this cell - it wont change during one ionization cycle
    }				//End of if statement to decide if this is the first iteration
  else
    {
      numerator = xxxplasma->PWnumer[ion_lower];	// We are on the second iteration, so use the stored value
    }
  if (numerator == 0.0)
    {
      return (0.0);		//There is no need to waste time evaluating the denominator
    }
  fmaxtemp = xtop->freq[xtop->np - 1];
  fmax = check_fmax (fthresh, fmaxtemp, xtemp);
  if (fthresh > fmax)
    {
      Error
	("pl_correct: After checking, fthresh has been set below fmin - we cannot compute denominator\n");
      q = 1.0;
      return (q);
    }
  qromb_temp = xtemp;		//The denominator is calculated for the LTE rate at our ideal temp. If we get this wrong, then divide by zeros abound!

  if (pow (((xtemp - xxxplasma->PWdtemp[ion_lower]) / xxxplasma->PWdtemp[ion_lower]), 2) > 1e-6)	//If our guess temperature hasnt changed much, use denominator from last time
    {
      denominator = qromb (tb_planck1, fthresh, fmax, 1.e-4);
      xxxplasma->PWdenom[ion_lower] = denominator;
      xxxplasma->PWdtemp[ion_lower] = xtemp;
    }
  else
    {
      denominator = xxxplasma->PWdenom[ion_lower];
    }



  q = numerator / denominator;	/*Just work out the correction factor - hope it isnt infinite! */





  return (q);
}



/**************************************************************************
                    Space Telescope Science Institute


  Synopsis:  

  Description:	

 tb_planck is the function a\nuB\nu and is called by Qromb in order to integrate over the frequency range where the
   ionisation cross section is significant. This is the function for ions with a topbase cross section NSH 16/12/10 


  Arguments:  (Input via .pf file)		


  Returns:

  Notes:

This is almost identical to code written to compute the sim power law correction. It is copied here to make the coding easier, plus it is likely that it will supplant the earlier code so that can all be deleted.


 

  History:

12Feb NSH - written as part of the varaible temperature effort.

 ************************************************************************/


double
tb_planck1 (freq)
     double freq;
{
  double answer, bbe;
  bbe = exp ((H * freq) / (BOLTZMANN * qromb_temp));
  answer = (2. * H * pow (freq, 3.)) / (pow (C, 2));
  answer *= (1 / (bbe - 1));
//      answer*=weight;
  answer *= sigma_phot_topbase (xtop, freq);
  answer /= freq;

  return (answer);
}




/**************************************************************************
                    Space Telescope Science Institute


  Synopsis:  

  Description:	
	verner_planck is the function a\nuB\nu and is called by Qromb in order to integrate over the frequency range where the
	ionisation cross section is significant. This is the function for ions with a verner cross section NSH 16/12/10 

  Arguments:  


  Returns:

  Notes:

This is almost identical to code written to compute the sim power law correction. It is copied here to make the coding easier, plus it is likely that it will supplant the earlier code so that can all be deleted.


 

  History:

12Feb NSH - written as part of the varaible temperature effort.
13Sep	NSH - no longer called as of python76c - all data is in tabulated topbase type format - this can be deleted in time

 ************************************************************************/

/*
double
verner_planck1 (freq)
     double freq;
{
  double answer, bbe;
  bbe = exp ((H * freq) / (BOLTZMANN * qromb_temp));
  answer = (2. * H * pow (freq, 3.)) / (pow (C, 2));
  answer *= (1 / (bbe - 1));
//      answer*=weight;
  answer *= sigma_phot (xver, freq);
  answer /= freq;

  return (answer);
}
*/

/**************************************************************************
                    Space Telescope Science Institute


  Synopsis:  

  Description:	
	tb_pow is the function to allow integration of a power law photon distribution multiplied by the inoisation 
	cross section
  	to allow the numerator of the sim correction factor to be calulated for ions with topbase cross sections.

  Arguments:  


  Returns:

 Notes:

This is almost identical to code written to compute the sim power law correction. It is copied here to make the coding easier, plus it is likely that it will supplant the earlier code so that can all be deleted.


 

  History:

12Feb NSH - written as part of the varaible temperature effort.
13Nov NSH - changed to work in logspace - and name changed.
 ************************************************************************/


/* This is the old non-log version */
//double
//tb_pow1 (freq)
//     double freq;
//{
//  double answer;
//
//  answer = xpl_w * (pow (freq, (xpl_alpha - 1.0)));
//  answer *= sigma_phot_topbase (xtop, freq);	// and finally multiply by the cross section.
//
//  return (answer);
//}



double
tb_logpow1 (freq)
     double freq;
{
  double answer;

//  answer = xpl_w * (pow (freq, (xpl_alpha - 1.0)));

  answer = pow(10,xpl_logw+(xpl_alpha-1.0)*log10(freq));


  answer *= sigma_phot_topbase (xtop, freq);	// and finally multiply by the cross section.

  return (answer);
}



/**************************************************************************
                    Space Telescope Science Institute


  Synopsis:  
	verner_pow is the function to allow integration of a power law photon distribution multiplied by the inoisation cross section
   	to allow the numerator of the sim correction factor to be calulated for ions with verner cross sections.

  Description:	

  Arguments:  (Input via .pf file)		


  Returns:

  Notes:

This is almost identical to code written to compute the sim power law correction. It is copied here to make the coding easier, plus it is likely that it will supplant the earlier code so that can all be deleted.


 

  History:

12Feb NSH - written as part of the varaible temperature effort.
13Sep	NSH - no longer called as of python76c - all data is in tabulated topbase type format - this can be deleted in time
 ************************************************************************/

/*
double
verner_pow1 (freq)
     double freq;
{
  double answer;

  answer = xpl_w * (pow (freq, (xpl_alpha - 1.0)));
  answer *= sigma_phot (xver, freq);	// and finally multiply by the cross section.

  return (answer);
}
*/

/**************************************************************************
                    Space Telescope Science Institute


  Synopsis:  
	verner_exp is the function to allow integration of an exponential photon distribution multiplied by the inoisation cross section
   	to allow the numerator of the sim correction factor to be calulated for ions with verner cross sections.

  Description:	

  Arguments:  (Input via .pf file)		


  Returns:

  Notes:

This is almost identical to code written to compute the sim power law correction. It is copied here to make the coding easier, plus it is likely that it will supplant the earlier code so that can all be deleted.


 

  History:

12Aug NSH - Written as part of the effort to improve the modelling
13Sep	NSH - no longer called as of python76c - all data is in tabulated topbase type format - this can be deleted in time
 ************************************************************************/

/*
double
verner_exp1 (freq)
     double freq;
{
  double answer;

  answer = xexp_w * exp ((-1.0 * H * freq) / (BOLTZMANN * xexp_temp));
  answer *= sigma_phot (xver, freq);	// and finally multiply by the cross section.
  answer /= freq;
  return (answer);
}
*/



/**************************************************************************
                    Southampton University


  Synopsis:  

  Description:	
	tb_exp is the function to allow integration of an exponential photon distribution multiplied by the inoisation 
	cross section
  	to allow the numerator of the correction factor to be calulated for ions with topbase cross sections.

  Arguments:  


  Returns:

 Notes:

This is almost identical to code written to compute the sim power law correction. It is copied here to make the coding easier, plus it is likely that it will supplant the earlier code so that can all be deleted.


 

  History:

12Aug Written by NSH as part of the effort to improve spectral modelling
 ************************************************************************/



double
tb_exp1 (freq)
     double freq;
{
  double answer;

  answer = xexp_w * exp ((-1.0 * H * freq) / (BOLTZMANN * xexp_temp));
  answer *= sigma_phot_topbase (xtop, freq);	// and finally multiply by the cross section.
  answer /= freq;
  return (answer);
}
