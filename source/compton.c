
/***********************************************************/
/** @file  new_compton.c
 * @author ksl
 * @date   January, 2018
 *
 * @brief  These are the routines to handle compton scattering.
 *
 * ???
 ***********************************************************/


#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "atomic.h"
#include "python.h"

PlasmaPtr xplasma;              /// Pointer to current plasma cell

/**********************************************************/
/** 
 * @brief      computes the effective compton opacity in the cell.
 *
 * @param [in] PlasmaPtr  xplasma   - pointer to the current plasma cell
 * @param [in] double  freq   - the frequency of the current photon packet being tracked
 * @return     kappa - the compton opacity for the cell.
 *
 * @details
 * This routine works out the opacity of a cell to compton (Thompson) scattering.
 * It is designed to give the correct *heating* effect when high frequency photons
 * undergo Compton scattering.
 * We use a fitted cross section given in Hazy 3 to calculate the frequency depedant
 * cross section. This cross section takes account of the energy exchange, i.e.
 * not only the cross section is freqneucy dependant, but so is the mean energy
 * transfer.- just using the Klein-Nishina cross section will give the wrong results in terms
 * of actual heating effect of these interactions.
 *
 * ### Notes ###
 * This implements the equation kappa=(sigmaT_h*ne*freq)/(me*c^2)
 *             Initally it is called by radiation
 *
 **********************************************************/

double
kappa_comp (xplasma, freq)
     PlasmaPtr xplasma;         // Pointer to current plasma cell
     double freq;               // Frequency of the current photon being tracked
{
  double x;                     // The opacity of the cell by the time we return it.
  double sigma;                 /*The cross section, thompson, or KN if hnu/mec2 > 0.01 */
  int ndom;

  ndom = wmain[xplasma->nwind].ndom;    //work out the domain we are looking at

  sigma = alpha (freq) * THOMPSON;      //obtain the energy exchange cross section

  x = (sigma * H) / (MELEC * C * C);    //Calculate the constant
  x *= xplasma->ne * freq;      //Multiply by cell electron density and frequency of the packet.

  x *= zdom[ndom].fill;         // multiply by the filling factor- should cancel with density enhancement
  return (x);
}




/**********************************************************/
/** 
 * @brief      Computes the effective opacity in the cell due to *induced* compton heating.
 *
 * @param [in] PlasmaPtr  xplasma   - pointer to the current plasma cell
 * @param [in] double  freq   - the frequency of the current photon packet being tracked
 * @return     kappa - the inducd compton opacity for the cell.
 *
 * @details
 * This routine implements a rather unusual situation, where the mean intensity in
 * a cell exceeds the blackbody value and we get induced compton heating.
 *
 * ### Notes ###
 * This implements the induced compton term in equation 6.5 in Hazy vol 3 - this is
 * the last term in the equation. It will
 * rarely if ever occur in a real model, but was included to give better agreement
 * with cloudy over a wide range of ionization parameters.
 *
 **********************************************************/

double
kappa_ind_comp (xplasma, freq)
     PlasmaPtr xplasma;         // Pointer to current plasma cell
     double freq;               // Frequency of the current photon being tracked
{
  double x;                     // The opacity of the cell by the time we return it.
  double sigma;                 /*The cross section, thompson, or KN if hnu/mec2 > 0.01 */
  double J;                     //The estimated intensity in the cell
  int ndom;

  ndom = wmain[xplasma->nplasma].ndom;


  J = 0.0;                      /* NSH 130605 to remove o3 compile error */

  /* Obtain a model for the mean intensity - we call this with mode=2, which means
     that if we have not yet completed a cycle, dont return a dilute blackbody
     estimate if we are in PL mode. */
  J = mean_intensity (xplasma, freq, 2);


  sigma = THOMPSON * alpha (freq);      //obtain the energy exchange cross section

  x = (xplasma->ne) / (MELEC);
  x *= sigma * J;
  x *= 1 / (2 * freq * freq);

  x *= zdom[ndom].fill;         // multiply by the filling factor- should cancel with density enhancement

  if (sane_check (x))           //For some reason we have a problem - we will not crash out - induced compton is unlikely to ever be dominant...
  {
    Error ("kappa_ind_comp:sane_check - undefined value for Kappa_ind_comp - setting to zero\n");
    return (0.0);
  }

  return (x);
}


/**********************************************************/
/** 
 * @brief      computes the cooling in the cell due to inverse compton scattering.
 *
 * @param [in] WindPtr  one   - pointer to the current wind cell
 * @param [in] double  t_e   - electron temperature that we need the cooling rate for
 * @return     the compton luminosity for the cell.
 *
 * @details
 * This routine computes the temperature dependant compton cooling rate in ergs/second.
 * It is mainly called as part of the search for a new temperature balancing heating
 * and cooling rates after an ionization cycle.
 *
 *
 * ### Notes ##
 * Compton cooling or 'inverse compton scattering' is the process of energy transfer
 * from hot electrons to 'cool' phortons, that is when 4kT > hnu.
 * For each interaction, the fractional energy loss from the electrons is given
 * by \frac{}\delta e}{e} = 4k_BT/mc^2.
 * Where we have a model for the mean intensity of radiation in a cell, we are
 * able to compute the cooling rate by doing an integration over the mean intensity
 * multiplied by the frequency dependant cooling cross section.
 * which looks like C=4k_BT/mc^2 * 4/pi \int{\sigma_T * beta * J_{\nu} d\nu}.
 * Beta is the effective cooling cross section from fits in Hazy 3.
 * If we are in the low frequency limit, we dont need to integrate, because the
 * cooling cross section is not frequency depedant in the Thompson limit so we
 * can just multiply by the band limited mean intensity.
 * Simiularly, if we dont have a model - we just do the best we can which is
 * to multiply by the total integrated mean indensity times the Thormpson
 * cross section.
 *
 **********************************************************/

double
total_comp (one, t_e)
     WindPtr one;               // Pointer to the current wind cell - we need the cell volume, this is not in the plasma structure
     double t_e;                //Current electron temperature of the cell
{
  double x, f1, f2;             //The returned variable
  int nplasma, j;               //The cell number in the plasma array
//  PlasmaPtr xplasma;            //pointer to the relevant cell in the plasma structure - moved to local variable to allow integration


  nplasma = one->nplasma;       //Get the correct plasma cell related to this wind cell
  xplasma = &plasmamain[nplasma];       //copy the plasma structure for that cell to local variable

  x = 0.0;

  if (xplasma->comp_nujnu < 0.0)        //Since J_nu is constant for a given cycle - we only need to compute the integral once when searching for a thermal balance
  {
    if (geo.spec_mod == 1)      //Check to see if we have generated a spectral model
    {
      for (j = 0; j < geo.nxfreq; j++)  //We loop over all the bands
      {
        if (xplasma->spec_mod_type[j] != SPEC_MOD_FAIL) //Only bother doing the integrals if we have a model in this band
        {
          f1 = xplasma->fmin_mod[j];    //NSH 131114 - Set the low frequency limit to the lowest frequency that the model applies to
          f2 = xplasma->fmax_mod[j];    //NSH 131114 - Set the high frequency limit to the highest frequency that the model applies to
          if (f1 > 1e18)        //If all the frequencies are lower than 1e18, then the cross section is constant at sigmaT
            x += qromb (comp_cool_integrand, f1, f2, 1e-6);
          else
            x += THOMPSON * xplasma->xj[j];     //In the case where we are in the thompson limit, we just multiply the band limited frequency integrated mean intensity by the Thompson cross section
        }
      }
    }

    else                        //If no spectral model - we do the best we can, and multply the mean intensity by the thompson cross section.
    {
      x = THOMPSON * xplasma->j;
    }
    xplasma->comp_nujnu = x;    //Store the result of the integral
  }
  else
    x = xplasma->comp_nujnu;    //We already have an integral, retrieve it

  x *= (16. * PI * BOLTZMANN * t_e * xplasma->ne) / (MELEC * C * C) * xplasma->vol;     //Multply by the other terms - including temperature - this gives the temperature dependance of this cooling term.


  return (x);
}


/**********************************************************/
/** 
 * @brief      computes the Klein Nishina cross section for a photon of frequency nu.
 *
 * @param [in] double  nu   - photon frequency
 * @return     the KN cross section.
 *
 * @details
 * This routine is used to compute the frequency dependant cross section
 * of a photon to electron scattering. For a frequency of less than
 * 1e18Hz this is just the Thompson scattering cross section. It is currently
 * only used to compute the chance of a photon interacting in this way. The
 * actual eneergy loss when an interaction occurs requires a cross section
 * that takes account of energy balance (provided by the alpha function.)
 *
 * ### Notes ###
 * This implements equation 7.5 in Rybicki and Lightman
 *
 **********************************************************/

double
klein_nishina (nu)
     double nu;                 //The frequency of the photon packet
{
  double x;                     //h nu / kt
  double x1, x2, x3, x4;        //variables to store intermediate results.
  double kn;                    // the final cross section

  kn = THOMPSON;                /* NSH 130605 to remove o3 compile error */
  x1 = x2 = x3 = x4 = 0.0;      /* NSH 130605 to remove o3 compile error */
  x = (H * nu) / (MELEC * C * C);       //The photon energy relative to the rest mass energy of an electron
  if (x > 0.0001)               //If the photon energy is high enough - then we need to compute the cross section - otherwise it is just the Thmopson cross section
  {
    x1 = 1. + x;
    x2 = 1. + (2. * x);
    x3 = log (x2);
    x4 = ((2. * x * x1) / x2) - x3;
    x4 *= x1 / (x * x * x);
    x4 = x4 + x3 / (2. * x);
    x4 = x4 - (1 + 3. * x) / (x2 * x2);
    kn *= 0.75 * x4;
  }

  return (kn);
}

//External variables to allow zfunc to search for the correct fractional energy change

double sigma_rand;              //The randomised cross section that our photon will see
double sigma_max;               //The cross section for the maxmimum energy loss
double x1;                      //The ratio of photon eneergy to electron energy

/**********************************************************/
/** 
 * @brief      computes a random direction for a photon undergoing compton scattering.
 *
 * @param [in] PhotPtr  p   - the photon currently being scattered
 * @param [in] PlasmaPtr  xplasma   - the cell in which the scattering is taking place
 * @return     0 if successful - the new direction is returned as part of the photon structure
 *
 * @details
 * This routine computes a new direction (and energy change) for a photon undergoing compton
 * scattering. If the freqiuency is low with respect to the electron rest mass then it is
 * totaslly elastic, isotropic Thompson scattering.
 *
 * ### Notes ###
 * EVerything is calculated in the rest frame of the electron.
 * In this frame, the photon energy change is E/E'=1+(hnu/mc^2)(1+cos\theta)
 * where \theta is the angle thruogh which the photon is deflected.
 * This is a maximum for \theta=180 - E/E'=1+2hnumc^2
 * and minimum for \theta=0 - E/E'=1
 * We compute everything by first drawing a random cross section that our photon packet will see.
 * This cross section represents an energy change and hence a direction. The
 * distribnution of angles is taken care of by using a differential cross section
 * vs energy change function. So we work out what energy change and hence direction
 * is implied by our random cross section.
 *
 **********************************************************/

int
compton_dir (p, xplasma)
     PhotPtr p;                 // Pointer to the current photon
     PlasmaPtr xplasma;         // Pointer to current plasma cell

{
  double f_min, f_max, f;       //Fractional energy changes - E_old/E_new - minimum possible, maximum possible, actual as implied by random cross section
  double n, l, m, phi, len;     //The direction cosines of the new photon direction in the frame of reference with q along the photon path
  struct basis nbasis;          //The basis function which transforms between the photon frame and the observer frame
  double lmn[3];                /* the individual direction cosines in the rotated frame */
  double x[3];                  /*photon direction in the frame of reference of the original photon */
  double dummy[3], c[3];

  x1 = H * p->freq / MELEC / C / C;     //compute the ratio of photon energy to electron energy. In the electron rest frame this is just the electron rest mass energy

  n = l = m = 0.0;              //initialise some variables to avoid warnings


  if (x1 < 0.0001)              //If the photon energy is much less than electron rest mass energy, we just have Thompson scattering - scattering is isotropic
  {
    randvec (lmn, 1.0);         //Generate a normal isotropic scatter
    f = 1.0;                    //There is no energy loss
    stuff_v (lmn, p->lmn);
  }
  else
  {
    sigma_rand = random_number (0.0, 1.0);      //Generate a random number between 0 and 1 - this represents a randomised cross section (normalised to the maximum which out photon packet sees
    f_min = 1.;                 //The minimum energy change - i.e. no energy loss - the scattering angle is zero - the photon does not chage direction
    f_max = 1. + (2. * x1);     //The maximum energy change - this occurs if the scattering angle is 180 degrees (i.e. the photon bounces straight back.) f=e_old/e_new
    sigma_max = sigma_compton_partial (f_max, x1);      //Communicated externally to the integrand function in the zbrent call below, this is the maximum cross section, used to scale the K_N function to lie between 0 and 1. This is essentually the chance of a photon scattering through 180 degrees - or the angle giving the maximum energy loss
    f = zbrent (compton_func, f_min, f_max, 1e-8);      //Find the zero point of the function compton_func - this finds the point in the KN function that represents our randomised fractional energy loss z_rand.

/*We now have the fractional energy change f - we use the 'normal' equation for compton scattering to obtain the angle cosine n=cos(\theta)	for the scattering direction*/

    n = (1. - ((f - 1.) / x1)); //This is the angle cosine of the new direction in the frame of reference of the photon - this gives a 2D scattering angle

    if (isfinite (len = sqrt (1. - (n * n))) == 0)      //Compute the length of the other angle cosines - the isfinite is to take care of the very rare occasion where n=1!
      len = 0.0;                // If n=1, then the photon has either been undeflected or bounced straight back. Both are vanishingly unlikely but need to be treated.
    phi = 0.0;                  //no need to randomise phi, the random rotation of the vector generating the basis function takes care of this

    l = len * cos (phi);        //compute the angle cosines of the other two dimensions.
    m = len * sin (phi);

    randvec (dummy, 1.0);       //Get a random vector of length 1 - this is returned in the array dummy
    cross (dummy, p->lmn, c);   //c will be perpendicular to p->lmn (the original photon direction) *and* the random vector just computed
    create_basis (p->lmn, c, &nbasis);  //create a basis with the first axis in the direction of the original photon direction, c will be perpendicular to the photon direction. Orientaion of the y/z axes are will give  arandomization of the phi axis

    x[0] = n;                   //This is the cosine direction of the new photon direction, in the frame of reference where the first axis is in the original photon direction
    x[1] = l;
    x[2] = m;

    project_from (&nbasis, x, lmn);     /* Project the vector from the FOR of the original photon into the observer frame */
    renorm (lmn, 1.0);          //Make sure the length of the direction vector is equal to 1
    stuff_v (lmn, p->lmn);      //Put the new photon direction into the photon structure
  }

  p->freq = p->freq / f;        //reduce the photon frequency by the fractional energy change
  p->w = p->w / f;              //reduce the photon weight by the same ammount to conserve photon numbers
  return (0);
}


/**********************************************************/
/** 
 * @brief      Function used to find the fractional energy change for a given cross-section in the compton scattering limit.
 *
 * @param [in out] double  f   The fractional energy change of a photon undergoing scattering
 * @return     the difference between sigma(f) and the randomised cross section we are searching for
 *
 * @details
 * This function calculates the difference between the KN cross section for a given fractional
 * energy loss f/x1 and the randomly generated cross section for our current interaction. A zero
 * finding function is used to work out where this function is equal to zero - we need to
 * do this because we cannot invert the function sigma_compton_partial.
 *
 * ### Notes ###
 * Since this function is called by zbrent, many of the variables have to be communicated
 * externally. These are
 * x1, the ratio of incoming photon energy to electron rest mass,
 * sigma_max, the KN cross section for the maximum possible energy loss
 * sigma_rand - a randomised number between 0 and 1 representing the normalised cross section we want to find
 *
 **********************************************************/

double
compton_func (f)
     double f;
{
  double ans;
  ans = (sigma_compton_partial (f, x1) / sigma_max) - sigma_rand;
  return (ans);
}

/**********************************************************/
/** 
 * @brief      is the Klein Nishina cross section as a function of the fractional energy change
 *
 * @param [in] double  f   - the fractional energy change
 * @param [in] double  x   the energy of the incoming photon divided by the rest mass energy of an electron
 * @return     the cross section for the scattering angle that gives this fractional energy change.
 *
 * @details
 * This is a functiuon which returns the Klein-Nishina cross section implied by a supplied
 * fractional energy change for an incoming photon.
 *
 * ### Notes ###
 * Stuart Sim supplied this formula, but coundn't recall where he had found it, although
 * 			he had rederived it and thinks it is correct!
 *
 **********************************************************/

double
sigma_compton_partial (f, x)
     double f;                  //This is the fractional energy change, nu/nu'
     double x;                  //h nu/mec**2 - the energy of the photon divided by the rest energy of an eectron
{
  double term1, term2, term3, tot;

  term1 = ((x * x) - (2 * x) - 2) * log (f) / x / x;
  term2 = (((f * f) - 1) / (f * f)) / 2;
  term3 = ((f - 1) / x) * ((1 / x) + (2 / f) + (1 / (x * f)));

  tot = 3 * THOMPSON * (term1 + term2 + term3) / (8 * x);

  return (tot);

}

/**********************************************************/
/** 
 * @brief      The 'heating' cross section correction to Thompson
 *
 * @param [in] double  nu   - frequency
 * @return     alpha- the effective heating cross section relative to sigma_T
 *
 * @details
 * We calculate compton heating due to photons passing through
 * a parcel of gas using an opacity. In order to work out the opacity, we
 * need the cross section that the electrons present to the photons, and
 * in order to take proper account of the energy transfer during a scatter
 * we use a 'heating crosss ection'. This is equal to 'alpha' times
 * the Thompson cross section - we calculate alpha here.
 *
 * ### Notes ###
 * The fit used here is given in Hazy 3, eq 6.6
 *
 **********************************************************/

double
alpha (nu)
     double nu;
{
  double alpha;
  if (nu < 1e17)
    alpha = 1.0;
  else
    alpha = 1. / (1. + nu * HRYD * (1.1792e-4 + 7.084e-10 * nu * HRYD));
  return (alpha);
}


/**********************************************************/
/** 
 * @brief      The 'cooling' cross section correction to Thompson
 *
 * @param [in] double  nu   - frequency
 * @return     beta - the effective cooling cross section relative to sigma_T
 *
 * @details
 * The inverse compton cooling rate in a cell is computed by ingegrating
 * the mean intensity in that cell by the frequency dependant cross section.
 * In order to take account of the energy transfer in an interaction, the
 * cross section must be the effective cooling cross section (not just the
 * Klein Nishina cross section)
 *
 * ### Notes ###
 * The fit used here is given in Hazy 3, eq 6.6
 *
 **********************************************************/

double
beta (nu)
     double nu;
{
  double alp, beta;
  if (nu < 1e17)
    beta = 1.0;
  else
  {
    alp = alpha (nu);
    beta = (1. - alp * nu * HRYD * (1.1792e-4 + (2 * 7.084e-10 * nu * HRYD)) / 4.);
  }
  return (beta);
}



/**********************************************************/
/** 
 * @brief      The integrand in the integral to obtain the total cooling rate due to inverse compton scattering.
 *
 * @param [in] double  nu   - frequency
 * @return     frequwncy dependant cross section multiplied by the mean intensity at frequency nu
 *
 * @details
 * This is is the integrand sigma x J_nu that is integrated
 * 	to obtain the compton cooling rate in a cell. The sigma in question is the
 * effective energy exchange cross section, which for coolnig is /alpha/beta/sigma_T.
 * These alpha and beta terms are computed in ther subroutines, but here we only
 * need to use \beta, because the \alpha term is taken account of when we compute
 * J_nu during the photon transport phase of the code.
 * If we are calling this function, we know we have a model so there is no need to
 * protect against the problem of not having a model
 *
 * ### Notes ###
 * See Hazy 3 eq 6.6 for more details of this approach
 *
 **********************************************************/

double
comp_cool_integrand (nu)
     double nu;
{
  double value;
  value = THOMPSON * beta (nu) * mean_intensity (xplasma, nu, 2);
  return (value);
}
