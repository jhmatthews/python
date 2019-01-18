
/***********************************************************/
/** @file  variable_temperature.c
 * @author nsh,ksl
 * @date   May, 2018
 *
 * @brief
 * This file was created to contain all of the routines associated
 * with calculating ionization using the so-called pair-wise
 * dilute bb and pairwise power law models.  Collectively these are
 * sometimes referred to as variable temperature models.  Both
 * involve
 * attempt to calculate ionic abundances using a
 * temperature fixed for each pair of ions to make the uncorrected
 * abundance ratio roughly equal to 1. There are then correction
 * factors applied, either based upon a piece_wise power law model of
 * the true radiation field, or as a dilute blackbody.
 *
 ***********************************************************/



#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>


#include "atomic.h"
#include "python.h"


/** Make this a variable that all the subroutines cn see,
 * so we can decide if we need to recompute the numerators */
int niterate;



//OLD /***********************************************************
//OLD                                        Southampton University
//OLD
//OLD   Synopsis:
//OLD
//OLD int
//OLD variable_temperature (xplasama, mode)  modifies the densities of ions, levels, and
//OLD   partition functions of ions within a cell of the wind based upon the mode,
//OLD   and other data contained within the WindPtr itself, such as t_r, t_e, w,
//OLD   based on the "mode".
//OLD
//OLD         Unlike nebular_concentrations, this code works on pairs of ions at a time,
//OLD   and uses a temperature calculated to be suitable for that particular pair
//OLD   to have a ratio of abundances about equal to 1.
//OLD
//OLD   Arguments:
//OLD      PlasmaPtr ww;
//OLD      int mode;                    // 6=correct using a dilute blackbody
//OLD                                      7=correct using a broken power law spectrum
//OLD
//OLD
//OLD   Returns:
//OLD   variable temperature alters portions of the wind ptr.  Exactly how things
//OLD   are changed depends on the mode.
//OLD
//OLD   variable_temperature returns 0 if it converged to an answer, -1 otherwise.  On an
//OLD   abnormal return the density array and ne are left unchanged.
//OLD
//OLD   Description:
//OLD
//OLD
//OLD
//OLD   Notes:
//OLD
//OLD   Section 5.3 of Nick's thesis explains the undelying ratinale for the so-called
//OLD   variable temperature models.  The name of the routine is misleading since
//OLD   there correctios are made for both a dilute power law and for a dilute bb.
//OLD
//OLD
//OLD
//OLD   History:
//OLD   2012Feb nsh     Coded and debugged as part of QSO effort.
//OLD         1212Dec nsh       Recoded so that the densities are computed in a
//OLD                   temporary array, and only
//OLD                   committted to the real density structure once
//OLD                   we are sure the code converges.
//OLD   2013Sep nsh     ground state fudge computed at the start, and
//OLD                   stored in an array rather
//OLD                   than contiunually recomputing it - sometimes it is expensive.
//OLD
//OLD
//OLD **************************************************************/



double xxxne, xip;


/**********************************************************/
/**
 * @brief      modifies the densities of ions, levels, and
 * partition functions of ions within a cell of the wind based upon the mode,
 * and other data contained within the WindPtr itself, such as t_r, t_e, w,
 * based on the "mode".
 *
 * @param [in,out] PlasmaPtr  xplasma   A Plasma cell
 * @param [in] int  mode   A choice of NEBULARMODE_PAIRWISE_ML93 (6) which uses
 * a dilute blackbody to represent the spectrum or
 * NEBULARMODE_PAIRWISE_SPECTRALMODEL (7) which uses a crude spectrum constructed
 * of power law or exponetial segments.
 * @return     returns 0 if the routine converged to an answer, -1 otherwise.
 * On an
 * abnormal return the density array and ne are left unchanged.
 *
 * variable temperature alters portions of the variables associated with the Plasma cellr.  
 * Exactly how things
 * are changed depends on the mode.
 *
 *
 * @details
 *
 * ### Notes ###
 * Section 5.3 of Nick's thesis explains the undelying ratinale for the so-called
 * variable temperature models.  The name of the routine is misleading since
 * there correctios are made for both a dilute power law and for a dilute bb.
 *
 * Unlike nebular_concentrations, this code works on pairs of ions at a time,
 * and uses a temperature calculated to be suitable for that particular pair
 * to have a ratio of abundances about equal to 1.
 *
 * @bug There are various questions associated with this routine about whether
 * what is being done associated with partition funcitons is correct.
 *
 **********************************************************/

int
variable_temperature (xplasma, mode)
     PlasmaPtr xplasma;
     int mode;
{
  int nion;
  double xnew, xsaha;
  double theta, x;
  double get_ne ();
  double t_e, t_r, xtemp, nh, xne, xxne;
  double a, b;
  double newden[NIONS];         //A temporaray array for recording the denisites                    
  int nelem, first, last;
  double t_e_part_correct;
  double sum, big;
  double pi_fudge, recomb_fudge, tot_fudge;     /*Two of the correction factors for photoionization rate, and recombination rate */
  double gs_fudge[NIONS];       /*It can be expensive to calculate this, and it only depends on t_e - which is fixed for a run. So
                                   //                 calculate it once, and store it in a temporary array */

  nh = xplasma->rho * rho2nh;
  t_e = xplasma->t_e;
  t_r = xplasma->t_r;


  /* Copy the current densities into the temporary array */

  for (nion = 0; nion < nions; nion++)
  {
    newden[nion] = xplasma->density[nion];

    /* Here we populate the recombination to ground state correction factor used in the LM and Sim ionization
     * equations.  Mode 2 imples we are include the dielectronic correction. There is no recombination fudge for
     * the neutral ion */

    if (ion[nion].istate != 1)
    {
      gs_fudge[nion] = compute_zeta (t_e, nion - 1, 2); //We call with nion-1, so when we address the array, we will want to ask for nion
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
    t_e = MIN_TEMP;             /* fudge to prevent divide by zeros */

  xsaha = SAHA * pow (t_e, 1.5);

  theta = xsaha * exp (-ion[0].ip / (BOLTZMANN * t_e)) / nh;

  if (theta < THETAMAX)
  {
    x = (-theta + sqrt (theta * theta + 4 * theta)) / 2.;
    xne = xxne = xxxne = x * nh;
  }
  else
    xne = xxne = xxxne = nh;    /*xxne is just a store so the error can report the starting value of ne.
                                   xxxne is the shared variable so the temperature solver routine can access it */

  if (xne < 1.e-6)
    xne = xxxne = 1.e-6;        /* Set a minimum ne to assure we can calculate
                                   xne the first time through the loop */


  /* At this point we have an initial estimate of ne. */

  niterate = 0;
  while (niterate < MAXITERATIONS)
  {
    /* We loop over all elements */
    for (nelem = 0; nelem < nelements; nelem++)
    {

      first = ele[nelem].firstion;      /*first and last identify the position in the array */
      last = first + ele[nelem].nions;  /*  So for H which has 2 ions, H1 and H2, first will
                                           generally be 0 and last will be 2 so the for loop
                                           below will just be done once for nion = 1 */

      if (first < 0 || first >= nions || last < 0 || last > nions)
      {
        Error ("variable_temperature: Confusion for element %d with first ion %d and last ion %d\n", nelem, first, last);
        exit (0);
      }


      /* and now we loop over all the ions of this element */
      sum = newden[first] = 1.0;        /* set the density of the first ion of this element to 1.0
                                           -this is (always??) the neutral ion */

      big = pow (10., 250. / (last - first));   /*make sure that we will not overflow the last ion */


      for (nion = first + 1; nion < last; nion++)       /*nion is the upper ion of the pair,
                                                           last is one above the actual last ion,
                                                           so the last ion to be considered is last-1 */
      {

        tot_fudge = 0.0;        /* NSH 130605 to remove o3 compile error */

        /* now we need to work out the correct temperature to use */
        xip = ion[nion - 1].ip; //the IP is that from the lower to the upper of the pair
        xtemp = zbrent (temp_func, MIN_TEMP, 1e8, 10);  //work out correct temperature


        /* given this temperature, we need the pair of partition functions for these ions */
        partition_functions_2 (xplasma, nion, xtemp, 1);        //weight of 1 give us the LTE populations.


        /* and now we need the saha equation linking the two states at our chosen temp
           NB the partition function of the electron is included in the SAHA factor */

        xsaha = SAHA * pow (xtemp, 1.5);
        b = xsaha * xplasma->partition[nion] * exp (-ion[nion - 1].ip / (BOLTZMANN * xtemp)) / (xne * xplasma->partition[nion - 1]);


        t_e_part_correct = xplasma->partition[nion - 1] / xplasma->partition[nion];
        partition_functions_2 (xplasma, nion, t_e, 0);  /* Get the partition functions assuming all ions
                                                           are in the ground state (hanve the zero) */

        t_e_part_correct *= (xplasma->partition[nion] / xplasma->partition[nion - 1]);


        /* we now correct b to take account of the temperature and photon field
           t_r and www give the actual radiation field in the cell, xtemp is the temp we used
           t_e is the actual electron temperature of the cell */

        recomb_fudge = sqrt (t_e / xtemp);

        /* Correct the SAHA equation abundance pair using an actual radiation field modelled as a dilute black body */
        if (mode == NEBULARMODE_PAIRWISE_ML93)
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

            pi_fudge = pi_correct (xtemp, nion, xplasma, mode);

            tot_fudge = pi_fudge * recomb_fudge * (gs_fudge[nion] + xplasma->w * (1 - gs_fudge[nion]));
          }
        }

        /* correct the SAHA equation abundance pair using an actual radiation field modelled as a broken power law */

        else if (mode == NEBULARMODE_PAIRWISE_SPECTRALMODEL)
        {
          pi_fudge = pi_correct (xtemp, nion, xplasma, mode);


          tot_fudge = pi_fudge * recomb_fudge * gs_fudge[nion] * t_e_part_correct;
        }

        else
        {
          Error ("variable_temperature: unknown mode %d\n", mode);
          exit (0);
        }


        /* apply correction factors */
        b *= tot_fudge;
        if (b > big)
          b = big;              //limit step so there is no chance of overflow
        a = newden[nion - 1] * b;
        sum += newden[nion] = a;

        if (newden[nion] < 0.0)
        {
          Error ("variable_temperature: ion %i has a negative density %e\n", nion, newden[nion]);
        }
        if (sane_check (sum))
          Error ("variable_temperature:sane check failed for running density sum=%e, last ion=%i\n", sum, nion);
      }                         //end of loop over ions - we now have a full set of ions for this element



      a = nh * ele[nelem].abun / sum;   //the scaling factor to get the overall abundance right
      for (nion = first; nion < last; nion++)
      {
        newden[nion] *= a;      //apply scaling
        if (sane_check (newden[nion]))  //check nothing has gone crazy
          Error ("variable_temperature:sane check failed for density newden=%e, for ion=%i\n", newden[nion], nion);
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
        if ((ion[nion].macro_info == 1) && (geo.macro_simple == 0) && (geo.macro_ioniz_mode == 1))
        {
          newden[nion] = xplasma->density[nion];
        }

        /*Set some floor so future divisions are sensible */
        if (newden[nion] < DENSITY_MIN)
          newden[nion] = DENSITY_MIN;
      }

      /* Now determine the new value of ne from the ion abundances */
    }                           //end of loop over elements

    xnew = get_ne (newden);     /* determine the electron density for this density distribution */


    if (xnew < DENSITY_MIN)
      xnew = DENSITY_MIN;       /* fudge to keep a floor on ne */
    if (fabs ((xne - xnew) / (xnew)) < FRACTIONAL_ERROR || xnew < 1.e-6)
    {
      break;
    }
    xne = xxxne = (xnew + xne) / 2.;    /*New value of ne */
    niterate++;


    if (niterate == MAXITERATIONS)
    {
      Error ("variable_temperature: failed to converge t %.2g nh %.2g xnew %.2g\n", t_e, nh, xnew);
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
    if (ion[nion].macro_info == 0 || geo.macro_ioniz_mode == 0 || geo.macro_simple == 1)
    {
      xplasma->density[nion] = newden[nion];
    }
  }


  partition_functions (xplasma, NEBULARMODE_LTE_GROUND);

  /* This sets all ions to the ground state.
     We really need a better implementation of partition functions and levels
     for a power law illuminating spectrum. We found that if we didnt make this call,
     we would end up with undefined levels - which did really crazy things */

  return (0);
}



int bb_correct_err = 0;

/**********************************************************/
/**
 * @brief      calculates the photoionization rate
 *  correction factor for the number density of two adjacent ions
 *
 * @param [in] double  xtemp   the temperature at which the saha abundances were calculated
 * @param [in] int  nion   number of the upper ion in the pair for which the correction factor is
 * @param [in] PlasmaPtr  xplasma   The cell currnently under test - this supplied the code with
 * @param [in] int  mode   This is the ionization mode, which should be either IONMODE_PAIRWISE_ML93 (6) or IONMODE_PAIRWISE_SPECTRALMODEL (7)
 * @return   the photoionization rate correction factor
 *
 * @details
 *
 * In IONMODE_PAIRWISE_ML93, the correction factor is calculated assuming a dilute BB function.  In
 * IONMODE_PAIRWISE_SPECTRALMODEL the factor is calculated using a piecewise power law or exponentila function.
 *
 * ### Notes ###
 * LM is a special case of this where the temperatures
 * are the same.
 *
 *
 **********************************************************/

double
pi_correct (xtemp, nion, xplasma, mode)
     double xtemp;
     int nion;
     PlasmaPtr xplasma;
     int mode;
{
  double q;
  int ion_lower;
  double numerator, denominator;
  double w_store, t_r_store;


  numerator = denominator = 0.0;
  /*nion is the upper (recombining) ion in the pair, but calc_pi_rate uses the lower (ionizing) ion  */
  ion_lower = nion - 1;



  if (niterate == 0)            /*first time of iterating this cycle, so calculate the numerator. This rate only depends on the model of J, which doesn't change during the iteration cycle */
  {
    if (mode == IONMODE_PAIRWISE_ML93)
    {
      numerator = calc_pi_rate (ion_lower, xplasma, 2, 1);      /*Call calc_pi_rate with mode 2 (dilute BB), which returns the PI rate coeficient */
    }
    else if (mode == IONMODE_PAIRWISE_SPECTRALMODEL)    /*Call calc_pi using a seriews a series of power laws and/or exponentials to represent the spectrum */
    {
      numerator = calc_pi_rate (ion_lower, xplasma, 1, 1);
    }
    xplasma->PWnumer[ion_lower] = numerator;    /* Store the calculated numerator for this cell - it wont change during one ionization cycle */
  }                             //End of if statement to decide if this is the first iteration
  else
  {
    numerator = xplasma->PWnumer[ion_lower];    // We are on the second iteration, so use the stored value
  }

  /*If the numeratior is 0, there is no PI for this ion, and so there
   * is no need to waste time evaluationg the denominatory */
  if (numerator == 0.0)
  {
    return (0.0);
  }


  if (pow (((xtemp - xplasma->PWdtemp[ion_lower]) / xplasma->PWdtemp[ion_lower]), 2) > 1e-6)    //If our guess temperature hasnt changed much, use denominator from last time
  {
    /*calc_pi_rate uses the value of t_r and w in the plasma structure to compute the rate coefficient,
     * so we need to store the actual value of t_r and w so we can substititute our guess temperature
     * and w=1 (to get LTE) */

    t_r_store = xplasma->t_r;
    w_store = xplasma->w;
    xplasma->t_r = xtemp;
    xplasma->w = 1.0;

    /*Now we can call calc_pi_rate, which will now calulate an LTE rate coefficient at our guess temperature */
    denominator = calc_pi_rate (ion_lower, xplasma, 2, 1);

    /*Put things back the way they were */
    xplasma->t_r = t_r_store;
    xplasma->w = w_store;
    xplasma->PWdenom[ion_lower] = denominator;  /*Store the LTE rate coefficient for next time */
    xplasma->PWdtemp[ion_lower] = xtemp;        /*And also the temperature at which it was calculated */
  }
  else
  {
    denominator = xplasma->PWdenom[ion_lower];  /*The temperature hadnt changed much, so use the stored value */
  }


  q = numerator / denominator;  /*This is the PI rate correction factor */


  return (q);

}




/**********************************************************/
/**
 * @brief      the function minimised by zbrent to find a temperature
 * when the ion ratios are one (so the logarithm will be 0),
 * to avoid numerical problems.
 *
 *
 * @param [in] double  solv_temp   The temperature where the ???
 * @return    The result of the temp_func calculation
 *
 * @details
 * It is the natural log of the saha equation
 * with the ne taken to the RHS. The correction
 * factors are applied after and depend on the temperature we find.
 *
 *
 * ### Notes ###
 * xxxne and xip are global variables which is declared above and
 * assigned in the main variable_temperature routine.
 *
 **********************************************************/

double
temp_func (solv_temp)
     double solv_temp;
{
  double answer;
  answer = log (SAHA / xxxne) + 1.5 * log (solv_temp) - (xip / (BOLTZMANN * solv_temp));

  return (answer);
}
