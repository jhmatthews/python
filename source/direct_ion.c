#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "atomic.h"
#include "python.h"
#include "recipes.h"

#include <gsl/gsl_sf_expint.h>  //We need this gsl library to evaluate the first exponential integral

/* Ratio of hnu / kT beyond which we don't bother calculating 
   see #197 */
#define ALPHABIG_DIRECT_ION 100.   
/* direct_ion contains the routines relating to direct (collisional) ionization and 
threebody recombination */

/************************************************************
                                    University of Southampton
Synopsis:
  compute_di_coeffs computes the rates due to direct (collisional)
  ionization from Dere data.

Arguments:

       double t_e       electron temperature        
      
Returns:
  0 on success
       
Description:

Notes: 

History:
  
************************************************************/

int
compute_di_coeffs (T)
     double T;
{
  int n;

  for (n = 0; n < nions; n++)
    {
      if (ion[n].dere_di_flag == 0)
	{
	  di_coeffs[n] = 0.0;
	}

      else
  {
    di_coeffs[n] = q_ioniz_dere(n, T);
  }

    }   //End of loop over ions

  return (0);
}

/************************************************************
                                    University of Southampton
Synopsis:
  compute_qrecomb_coeffs computes the rates for three body recombination
  from Dere data.

Arguments:

       double t_e       electron temperature        
      
Returns:
  0 on success
       

Description:


Notes: 


History:
  1508 JM Coded

************************************************************/

int
compute_qrecomb_coeffs(T)
     double T;
{
  int n, nvmin, ntmin;
  struct topbase_phot *xtop;

  for (n = 0; n < nions; n++) //We need to generate data for the ions doing the recombining.
    {
      /* There is only any point doing this is we are not a neutral ion */

      if (ion[n].istate > 1)  //There is only any point doing this is we are not a neutral ion
    {

      if (ion[n-1].dere_di_flag == 0)
  {
    //printf ("NO COEFFS\n");
    qrecomb_coeffs[n] = 0.0;
  }

      else
  {
    /* we need to know about the bound-free jump, so we need the details
       from the ground state cross-section for for the ion below this one */

    ntmin = ion[n-1].ntop_ground;  /* We only ever use the ground state cont_ptr. 
                                      This is for topbase */
    nvmin = ion[n-1].nxphot;

    if (ion[n-1].phot_info > 0)  //topbase or hybrid VFKY (GS)+TB excited
      { 
        xtop = &phot_top[ntmin];
      }
    else if (ion[n-1].phot_info == 0)  // verner
      {   //just the ground state ionization fraction.
        xtop = &phot_top[nvmin];
      }
    else
      {
        Error ("compute_qrecomb_coeffs: no coll ionization data for ion %i\n",n-1);
      }
    
    /* this will return 0 if there aren't any coeffs */
    qrecomb_coeffs[n] = q_recomb_dere(xtop, T);
  }

    }   //End of if statement for neutral ions
	else //We are a neutral ion - so there can be no recombination
	{
		qrecomb_coeffs[n] = 0.0;
	}  
    }   //End of loop over ions

  return (0);
}



/************************************************************
                                    University of Southampton
Synopsis:
  total_di computed the total cooling due to direct (collisional)
  ionization from Dere data.

Arguments:

       WindPtr one      pointer to cell
       double t_e       electron temperature        
      
Returns:
  The total di cooling
       
Description:

Notes: 

History:
************************************************************/

double
total_di (one, t_e)
     WindPtr one;   // Pointer to the current wind cell - we need the cell volume, this is not in the plasma structure
     double t_e;    //Current electron temperature of the cell

{
  double x;            //The returned variable
  int nplasma;         //The cell number in the plasma array
  PlasmaPtr xplasma;   //pointer to the relevant cell in the plasma structure
  int n;               //loop pointers


  nplasma = one->nplasma;         //Get the correct plasma cell related to this wind cell
  xplasma = &plasmamain[nplasma]; //copy the plasma structure for that cell to local variable
  x = 0;                          //zero the luminosity


  compute_di_coeffs (t_e);        //Calculate the DR coefficients for this cell


  for (n = 0; n < nions; n++)
    {
      //We have no DI data for this ion
      if (ion[n].dere_di_flag == 0) 
  {
    x += 0.0;   //Add nothing to the sum of coefficients
  }
      else
  {
    
    x += xplasma->vol * xplasma->ne * xplasma->density[n] * di_coeffs[n] *
      dere_di_rate[ion[n].nxderedi].xi * EV2ERGS;

    //printf ("n=%i V=%e ne=%e rho=%e coeff=%e xi=%e cooling=%e\n",n, V , 
    //xplasma->ne , xplasma->density[n] , di_coeffs[n] ,
    //dere_di_rate[ion[n].nxderedi].xi*EV2ERGS,x);
  }
    }
  return (x);
}




/******************************************************************************/

/* q_ioniz_dere. This returns the collisional ionization co-efficient from Dere
data.
*/
/*could store global variables to avoid repeated interpolations */
//double t_e_last;
//double nion_last;

double
q_ioniz_dere (nion, t_e)
     int nion;
     double t_e;
{
  double coeff, t, scaled_t;
  double exp_int, dt, drdt, rate;
  int nrec, imax, imin, i;


  nrec = ion[nion].nxderedi;

  /* find the correct coefficient */


  t = (BOLTZMANN * t_e) / (dere_di_rate[nrec].xi * EV2ERGS);
  scaled_t = 1.0 - ((log (2.0)) / (log (2.0 + t)));

  /* 1/t is (hnu) / (k t_e). when this ratio is >100 we return 0, see #197 */
  if ( (1.0 / t) > ALPHABIG_DIRECT_ION)
    return 0.0;

  if (scaled_t < dere_di_rate[nrec].temps[0])  //we are below the range of DI data data
    {

      Log_silent("compute_di_coeffs: Requested temp %e is below limit of data for ion %i(Tmin= %e)\n",
                  scaled_t, nion, dere_di_rate[nrec].temps[0]);
      //rate = rates[0];
      imax = 1;
      imin = 0;
    }

  else if (scaled_t >= dere_di_rate[nrec].temps[dere_di_rate[nrec].nspline - 1]) //we are above the range of GS data
    {
      Log_silent("compute_di_coeffs: Requested temp %e is above limit (%e) of data for ion %i\n",
                  scaled_t, dere_di_rate[nrec].temps[dere_di_rate[nrec].nspline - 1], nion);

      //rate = rates[BAD_GS_RR_PARAMS - 1];
      imax = dere_di_rate[nrec].nspline - 1;
      imin = dere_di_rate[nrec].nspline - 2;
      //We will try to extrapolate.
    }

  else      //We must be within the range of tabulated data
      {
        for (i = 0; i < dere_di_rate[nrec].nspline - 1; i++)
    {
      if (dere_di_rate[nrec].temps[i] <= scaled_t && 
          scaled_t < dere_di_rate[nrec].temps[i + 1])  //We have bracketed the correct temperature
        {
          imin = i;
          imax = i + 1;
        }
    }
      /* NSH 140313 - changed the following lines to interpolate in log space */
      }


  drdt = ((dere_di_rate[nrec].rates[imax]) - (dere_di_rate[nrec].rates[imin])) / ((dere_di_rate[nrec].temps[imax]) - (dere_di_rate[nrec].temps[imin]));
  dt = ((scaled_t) - (dere_di_rate[nrec].temps[imin]));
  rate = dere_di_rate[nrec].rates[imin] + drdt * dt;

  coeff = pow (t, -0.5) * pow (dere_di_rate[nrec].xi, -1.5) * rate;

  if (exp (-1.0 / t) < (1.0 / VERY_BIG))
    {
      exp_int = 0.0;
    }
  else
    {
      exp_int = gsl_sf_expint_E1 (1.0 / t); //Evaluate the first exponential integral using gsl library.
    }

  coeff *= exp_int;
  
  return (coeff);
}


/******************************************************************************/

/* q_ioniz. This returns the collisional ionization co-efficient
Calculated following equation 5-79 of Mihalas or from data from
Dere 2007.
*/

double
q_ioniz (cont_ptr, electron_temperature)
     struct topbase_phot *cont_ptr;
     double electron_temperature;
{
  double coeff;
  double gaunt;
  double u0;
  int nion;

  /* these next two quantities only used in Hydrogen, no Dere data case */
  u0 = cont_ptr->freq[0] * H_OVER_K / electron_temperature;
  nion = cont_ptr->nion;
  gaunt = 0.1 * ion[nion].z;      //for now - from Mihalas for hydrogen and Helium

  /* when hnu / kT ratio is >100 we return 0, see #197 */
  if (u0 > ALPHABIG_DIRECT_ION)
    return (0.0);


  /* if ion[n].dere_di_flag == 1 then we have direct ionization data for this ion
     only do this if it is the ground state */
  if (ion[nion].dere_di_flag == 1 && config[cont_ptr->nlev].ilv == 1 && ion[nion].macro_info == 0)
    {
      coeff = q_ioniz_dere(nion, electron_temperature);
    }

  /* let's only apply the approximation from Mihalas for Hydogen and Helium */
  else if (ion[nion].z < 3 && ion[nion].macro_info == 1)
    {
      coeff = 1.55e13 / sqrt (electron_temperature) * gaunt * cont_ptr->x[0] *
        exp (-1. * u0) / u0;
    }

  /* otherwise return 0 */
  else
    coeff = 0.0;


  return (coeff);
}

/******************************************************************************/

/* q_recomb_dere. This returns the collisional recombination co-efficient
Calculated from inverse of q_ioniz_dere.

This equation comes from considering TE and getting the expression
q_recomb = 2.07e-16 * gl/gu * exp(E/kT) * q_ioniz
then substituting the above expression for q_ioniz.
*/

double
q_recomb_dere (cont_ptr, electron_temperature)
     struct topbase_phot *cont_ptr;
     double electron_temperature;
{
  int nion;
  double u0;
  double gaunt, coeff;
  double root_etemp;

  gaunt = 0.1;      //for now - from Mihalas for hydrogen
  u0 = cont_ptr->freq[0] * H_OVER_K / electron_temperature;
  nion = cont_ptr->nion;

  /* when hnu / kT ratio is >100 we return 0, see #197 */
  if (u0 > ALPHABIG_DIRECT_ION)
    return (0.0);

  /* if ion[n].dere_di_flag == 1 then we have direct ionization data for this ion
     only do this if it is the ground state */
  if (ion[nion].dere_di_flag == 1 && config[cont_ptr->nlev].ilv == 1)
    {
      root_etemp = sqrt(electron_temperature);     
      coeff = 2.07e-16 / (root_etemp * root_etemp * root_etemp);

      /* the original way of getting mutilplicity doesn't work for non-matoms, 
         because uplevel isn't identified */
      //coeff *= config[cont_ptr->nlev].g / config[cont_ptr->uplev].g;

      /* JM/NSH XXX -- This is the multiplicity of the ground states. 
         Should be of order 1 so it may be better to just leave it out, 
         since the collisional ionization cross section is doubtless 
         averaged over upper states... */   
      coeff *= ion[nion].g / ion[nion+1].g;   
  
      coeff *= exp(u0);
      coeff *= q_ioniz_dere(nion, electron_temperature);


      /* do a sane check here, as there's an exponential which could blow up */
      if (sane_check(coeff))
      {
        Error("q_recomb is %8.4e for ion %i at temperature %8.4e, setting to zero\n",
               coeff, nion, electron_temperature);
        coeff = 0.0;
      }
    }
  else 
    coeff = 0.0;

  return (coeff);
}


/******************************************************************************/

/* q_recomb. This returns the collisional recombination co-efficient
Calculated from inverse of q_ioniz.

JM 1301 -- Edited this to avoid need to call q_ioniz and exponential.
Should improve speed and stability

This equation comes from considering TE and getting the expression
q_recomb = 2.07e-16 * gl/gu * exp(E/kT) * (T_e**-3/2) * q_ioniz
then substituting the above expression for q_ioniz.
*/

double
q_recomb (cont_ptr, electron_temperature)
     struct topbase_phot *cont_ptr;
     double electron_temperature;
{
  double coeff;
  double gaunt, u0;
  int nion;

  nion = cont_ptr->nion;
  u0 = cont_ptr->freq[0] * H_OVER_K / electron_temperature;
  gaunt = 0.1 * ion[nion].z;      //for now - from Mihalas for hydrogen and Helium

  /* when hnu / kT ratio is >100 we return 0, see #197 */
  if (u0 > ALPHABIG_DIRECT_ION)
    return (0.0);

  /* if ion[n].dere_di_flag == 1 then we have direct ionization data for this ion
     only do this if it is the ground state */
  /* Still use the Mihalas approximation for macro-atoms */
  if (ion[nion].dere_di_flag == 1 && config[cont_ptr->nlev].ilv == 1 && ion[nion].macro_info == 0)
    {
      coeff = q_recomb_dere(cont_ptr, electron_temperature);
    }

  /* let's only apply the approximation from Mihalas for Hydogen and Helium macro-atoms */
  else if (ion[nion].z < 3 && ion[nion].macro_info == 1)
    {
      coeff = 3.2085e-3  / electron_temperature * gaunt * cont_ptr->x[0]; // normal constants * 1/T times gaunt * cross section

      coeff /= cont_ptr->freq[0] * H_OVER_K;      // divide by h nu / k
      coeff *= config[cont_ptr->nlev].g / config[cont_ptr->uplev].g;
    }

  /* otherwise return 0 */
  else
    coeff = 0.0;


  return (coeff);
}
