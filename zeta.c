/**************************************************************************
                    Southampton University
                                                                                                   
                                                                                                   
  Synopsis:

The routines in this file are to do with computing the correction to the saha equation.



  Description:

Zeta is the term in lucy and mazzali which corrects for the proportion of recombinations
going directly to the ground state. 

nsh: I'm experimenting with moving it out of saha and stuart_sim
since the same factor will be used in both, and I'm going to add in some more options to
incorporate a correction factor for dielectronic recombination.


                                                                                                   
  Arguments:  
                                                                                                   
                                                                                                   
  Returns:
                                                                                                   
  Notes:
                                                                                                   
                                                                                                   
                                                                                                   
  History:
	11aug	nsh	Began work
	111211	ksl	Began to comment on it

                                                                                                   
 ************************************************************************/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "atomic.h"
#include "python.h"
#include "recipes.h"




/**************************************************************************
                    Southampton University
                                                                                                   
                                                                                                   
  Synopsis: compute_zeta calculates the correction term in the L+M and sim equations
        for corrected saha.
                                                                                                   
  Description:
	mode 1 is the original version, it takes ilow, ihi and interpfrac and 
        interpolates in the ground state fraction table.
        mode 2 attempts to modify zeta by application of the dielectronic recombination
 	rate
                                                                                                   
  Arguments:  
	temperature - the temperature of the cell
	nion - the ion we are recombining from
  	ilow,ihi - the bracketing elements in the interpfrac table (can be the same in at extremes)
        interpfrac - the distnce between the two bracketing temps
	f1,f2 - the frequency band over which we are interested - NB this doesnt really do that much, 
	since integ_fb ends up resetting the limits.
	mode - this will allow control here, so we can have the 
		1: interpfrac type, 
		2: interpfrac but then incorporating DR as an additional correction.
                                                                                                   
                                                                                                   
  Returns:
	zeta, the correction factor to the recombination rates. 
                                                                                                   
  Notes: When compute_zeta is called, nion is the lower ion in the pair whose abundances are being calculated.

                                                                                                   
                                                                                                   
                                                                                                   
 History:	
		Aug 2011 NSH - coding started to try to incorportate DR into python
		Feb 2012 NSH - modified to include code to find the interpfrac parameter here
				this is to improve the new code for using multiple temperature
				saha equaions - each equaion needs its own correction factor, 
				so it makes better code to call it.
		Jul 2012 NSH - modified to include use of Badnell recombination coefficients
				to allow dielectronic recombination to be included in the 
				zeta term

                                                                                                   
 ************************************************************************/

double
compute_zeta (temp, nion, mode)
     double temp;
     int mode, nion;
{
  double zeta, interpfrac, dummy;
  int ihi, ilow;

#define MIN_FUDGE  1.e-10
#define MAX_FUDGE  10.
//#define MIN_TEMP         100. //Put into python.h

  zeta = 0.0;			/* NSH 130605 to remove o3 compile error */


  /* now get the right place in the ground_frac tables  CK */
  dummy = temp / TMIN - 1.;
  ilow = dummy;			/* have now truncated to integer below */
  ihi = ilow + 1;		/*these are the indeces bracketing the true value */
  interpfrac = (dummy - ilow);	/*this is the interpolation fraction */
  if (ilow < 0)
    {
      ilow = 0;
      ihi = 0;
      interpfrac = 1.;
    }
  if (ihi > 19)
    {
      ilow = 19;
      ihi = 19;
      interpfrac = 1.;
    }


  if (mode == 1)
    {
      /* This is the old method of getting zeta - just from the ground state tables computed by Christian */

      zeta =
	ground_frac[nion].frac[ilow] +
	interpfrac * (ground_frac[nion].frac[ihi] -
		      ground_frac[nion].frac[ilow]);
    //Log ("for t_e=%f, ilow=%i, ihi=%i, interpfrac=%f, zeta=%f\n",temp,ilow,ihi,interpfrac,zeta);
    }



  /* Best try at full blown zeta including DR, if we have the data, else default to the old way of doing things */
  else if (mode == 2)		
    {
      /* We call this with the lower ion in a pair, so if that ion has only one electron, 
         (ie. Carbon 6) then we cannot have DR into this ion, so there will be no DR rate 
         associated with it. NSH 140317 We also do this if we dont have DR data at all. */	
      if (ion[nion].istate == ion[nion].z || ion[nion].drflag == 0)	
	{
	  if (ion[nion].total_rrflag == 1)	//We have total RR data
	    {

	      if (ion[nion].bad_gs_rr_t_flag == 1 && ion[nion].bad_gs_rr_r_flag == 1)	//We have tabulated gs data
		{
		  zeta =
		    badnell_gs_rr (nion, temp) / total_rrate (nion, temp);
            //Log ("We have all the data, and zeta=%f gs=%e total=%e\n",zeta,badnell_gs_rr (nion, temp),total_rrate (nion, temp));
		}

	      else		//We are going to have to integrate
		{
		  //Log
		  //  ("We do not have tabulated GS data for state %i of element %i\n",
		  //   ion[nion].istate, ion[nion].z);
		  zeta = milne_gs_rr (nion, temp) / total_rrate (nion, temp);
		}
	    }

	  else
	    {
	      zeta =
		ground_frac[nion].frac[ilow] +
		interpfrac * (ground_frac[nion].frac[ihi] -
			      ground_frac[nion].frac[ilow]);
        //Log ("We dont have the data and zeta=%f\n",zeta);
	    }

	}


      else
	{
	  if (ion[nion].drflag > 0 && ion[nion].total_rrflag == 1)	//We have the two tabulated rates
	    {
	    //Log ("We have dr and total_rr/n");
	      compute_dr_coeffs (temp);
	      if (ion[nion].bad_gs_rr_t_flag == 1 && ion[nion].bad_gs_rr_r_flag == 1)	//We also have tabulated GS data
		{
		//Log("We have gs_rr/n");
		  zeta =
		    badnell_gs_rr (nion,
				   temp) / (total_rrate (nion,
							 temp) +
					    dr_coeffs[nion]);
		}
	      else		//We are going to have to integrate
		{
		//Log ("We are ging to have to integrate");
		  zeta =
		    milne_gs_rr (nion,
				 temp) / (total_rrate (nion,
						       temp) +
					  dr_coeffs[nion]);
		}
	    }
	  else
            {
	      zeta =
		ground_frac[nion].frac[ilow] +
		interpfrac * (ground_frac[nion].frac[ihi] -
			      ground_frac[nion].frac[ilow]);
                      Log ("We dont have the data and zeta=%f\n",zeta);
	    }
	}
    }


  /* if we got here then we don't understand the mode and there must be a problem */  
  else
    {
      Error ("Compute zeta: Unkown mode %i \n", mode);
    }

  //printf("element %d ion %d zeta %g\n", ion[nion].z, ion[nion].istate, zeta);

  return (zeta);
}
