#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>


#include "atomic.h"
#include "python.h"

struct topbase_phot *xtop;	//Topbase description of a photoionization x-section 

double qromb_temp;  //Temerature used in integrations - has to be an external variable so

double xpl_alpha, xpl_w, xpl_logw;
double xexp_temp, xexp_w;



/***********************************************************
                                       West Lulworth

  Synopsis:   

int
calc_pi_rate (nion,xplasma,mode)  calculates the photoionization rarte coefficient for ion	
				nion, based upon the mean intensity stored in cell xplasma
				The mode tells the subroutine wether we are modelling the
				mean intensity as a dilute blackbody (mode2) or as a series
				of power laws and exponentials (mode1). 
  
  Arguments:		
     ion	nion;		The ion we are interested in - this is the ion which is being
				photoionized
     PlasmaPtr xplasma;		The cell in question - note that the details of the model
				is stored in this structure.
     int mode;			1 - use a power law and or exponential model for J
     				2 - use a dilute BB model for J - t_r and w are those 
				parameters in xplasma


  Returns:
 	The photioinization rate coefficient for the ion.
 	
  Description:	


 
  Notes:

  This was created in Summner 2014 in preparation for matrix ionization solver. Previously, this code
	was contained in two subroutines bb_correct_2 and pl_correct_2.The functoinsality of thses two
	have ben combined into one - hence the requirement for the mode parameter.



  History:
	2014Aug NSH - coded

**************************************************************/




double 
calc_pi_rate (nion,xplasma,mode)
	PlasmaPtr xplasma;
	int nion;
	int mode;
{
  double q;
  int  n, j;
  double pi_rate;
  int ntmin, nvmin;
  double fthresh, fmax, fmaxtemp;
  double f1, f2;
  double exp_qromb, pl_qromb;


  exp_qromb = 1e-4;
  pl_qromb = 1e-4;


  if (-1 < nion && nion < nions)	//Get cross section for this specific ion_number
    {
      ntmin = ion[nion].ntop_ground;	/*We only ever use the ground state cross sections. This is for topbase */
      nvmin = nion;	/*and this is for verner cross sections */
    }
  else
    {
      Error ("calc_pi_rate: %d is unacceptable value of nion\n", nion);
      mytrap ();
      exit (0);
      return (1.0);
    }



  if (ion[nion].phot_info == 1)	//topbase
    {
      n = ntmin;
      xtop = &phot_top[n];
    }
  else if (ion[nion].phot_info == 0)	// verner
    {
      n = nvmin;		//just the ground state ionization fraction.
      xtop = &xphot_tab[ion[n].nxphot];
    }
  else
    {
	  Error
	    ("calc_pi_rate: No photoionization xsection for ion %d (element %d, ion state %d)\n",
	     nion, ion[nion].z, ion[nion].istate);
	  exit(0); /* NSH 1408 I have decided that this is actually a really serous problem - we have no business including an ion for which we have no photoionization data.... */
    }

  fthresh = xtop->freq[0];
  fmax = xtop->freq[xtop->np - 1];
  pi_rate = 0;





if (mode==1) //Modelled version of J
	{
      for (j = 0; j < geo.nxfreq; j++)	//We loop over all the bands
	{
	  xpl_alpha = xplasma->pl_alpha[j];
	  xpl_logw = xplasma->pl_log_w[j];
	  xexp_temp = xplasma->exp_temp[j];
	  xexp_w = xplasma->exp_w[j];
	  if (xplasma->spec_mod_type[j] > 0)	//Only bother doing the integrals if we have a model in this band
	    {
	      f1 = xplasma->fmin_mod[j]; //NSH 131114 - Set the low frequency limit to the lowest frequency that the model applies to
	      f2= xplasma->fmax_mod[j]; //NSH 131114 - Set the high frequency limit to the highest frequency that the model applies to
	      if (f1 < fthresh && fthresh < f2 && f1 < fmax && fmax < f2)	//Case 1- 
		{
		  if (xplasma->spec_mod_type[j] == SPEC_MOD_PL)
		    {
		      pi_rate += qromb (tb_logpow1, fthresh, fmax, pl_qromb);
		    }
		  else
		    {
		      pi_rate += qromb (tb_exp1, fthresh, fmax, exp_qromb);
		    }
		}
	      else if (f1 < fthresh && fthresh < f2 && f2 < fmax)	//case 2 
		{
		  if (xplasma->spec_mod_type[j] == SPEC_MOD_PL)
		    {
		      pi_rate += qromb (tb_logpow1, fthresh, f2, pl_qromb);
		    }
		  else
		    {
		      pi_rate += qromb (tb_exp1, fthresh, f2, exp_qromb);
		    }
		}
	      else if (f1 > fthresh && f1 < fmax && fmax < f2)	//case 3
		{
		  if (xplasma->spec_mod_type[j] == SPEC_MOD_PL)
		    {
		      pi_rate += qromb (tb_logpow1, f1, fmax, pl_qromb);
		    }
		  else
		    {
		      pi_rate += qromb (tb_exp1, f1, fmax, exp_qromb);
		    }
		}
	      else if (f1 > fthresh && f2 < fmax)	// case 4
		{
		  if (xplasma->spec_mod_type[j] == SPEC_MOD_PL)
		    {
		      pi_rate += qromb (tb_logpow1, f1, f2, pl_qromb);
		    }
		  else
		    {
		      pi_rate += qromb (tb_exp1, f1, f2, exp_qromb);
		    }
		}
	      else		//case 5 - should only be the case where the band is outside the range for the integral.
		{
		  pi_rate += 0;	// Add nothing - bit of a null statement, but makes the code look nice.
		}
	    }			//End of loop to only integrate in this band if there is power
	}

}
else if (mode==2)  //blackbody mode
	{
  fmaxtemp = xtop->freq[xtop->np - 1];
  fmax = check_fmax (fthresh, fmaxtemp, xplasma->t_r);
  if (fthresh > fmax)
    {
      Error
	("pi_rates: temperature too low - ion %i has no PI rate\n",nion);
      pi_rate = 0.0;
    }
  else
   {
  qromb_temp = xplasma->t_r;		
  pi_rate = xplasma->w * qromb (tb_planck1, fthresh, fmax, 1.e-4);
   }
    }

pi_rate=(4*PI*pi_rate)/H;

return(pi_rate);
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






