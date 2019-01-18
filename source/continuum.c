
/***********************************************************/
/** @file  continuum.c
 * @author ksl
 * @date May, 2018
 *
 * @brief  Rotines related speciically to generating MC spectra for
 * stellar atmospheres like grids of spectra
 *
 ***********************************************************/


#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "atomic.h"
#include "python.h"

#include "models.h"


//OLD /***********************************************************
//OLD                                        Space Telescope Science Institute
//OLD
//OLD  Synopsis:
//OLD
//OLD double one_continuum(spectype,t,g,freqmin,freqmax) gets a photon frequency
//OLD from a continuum grid
//OLD
//OLD Arguments:
//OLD
//OLD   spectype                                An index to the continum grid which one wants
//OLD                                           to use for obtaining a photon
//OLD   double t,g,                             Two parameters defining the model for which
//OLD                                           a photon is requried, often temperature and
//OLD                                           gravity, though this is not required
//OLD   freqmin,freqmax;                        minimum and maximum frequency of interest
//OLD
//OLD
//OLD Returns:
//OLD
//OLD
//OLD Description:
//OLD
//OLD
//OLD Notes:
//OLD   In versions of python prior to python_52, different routines
//OLD   were used to read kurucz models and hubeny models.  This was
//OLD   historical, since the routines ksl had to read the kurucz models
//OLD   read binary files (from a time when compputers were very slow).
//OLD
//OLD   The new routines are much more general, and are based on routines
//OLD   ksl had written for his fitting program kslfit.  The underlying
//OLD   routines can interpolate between grids of more than two dimensions
//OLD   although here we assume that all of the continuum grids are two
//OLD   dimensional.
//OLD
//OLD History:
//OLD   04aug   ksl     created from hub.c usied in all versions
//OLD                   of python prior to python_52.
//OLD
//OLD **************************************************************/

double old_t, old_g, old_freqmin, old_freqmax;
double jump[] = { 913.8 };


/**********************************************************/
/**
 * @brief      get a photon frequency from a grid of stored spectra
 * from a continuum grid
 *
 * @param [in] int  spectype   An index to the continum grid which one want got obtaining a photon
 * @param [in] double  t   The associated temperature (or first variable) defining the grid
 * @param [in] double  g   The associted gravity (or second variable)
 * @param [in] double  freqmin   The minimum frequency desiered for a photon
 * @param [in] double  freqmax   The maximum frequecny for a photon
 * @return   A frequency for a photon
 *
 * @details
 * This is a routine that is used to select a function from a grid of spectra
 * stored in the model structure.  For the purpose of this routine, the
 * model grid is assumed to be based on two variables t (temperature) and
 * gravity (g).  This reflects the fact that most of the grids used in Python
 * to date (e.g Kurucz models) are these kinds of grids.  The routine linearly
 * interpolates  on these to variables to produce a spectrum at the desired
 * t and g.  It then creates a coumulative distribution of the portion of
 * the spectrum bewtween fmin and fmax and then uses a random number to
 * obtain a frequncy for a single program.
 *
 * The routine should already work for any grid.
 *
 * ### Notes ###
 *
 * @bug  The model structure is general in the sense that it was intended
 * for situation with any number of variables.
 * However what is written here
 * is specific to the case of two variables. This is a problem, as we have
 * found in trying to use this in other situations.  A more general routine is needed.
 * The first step in doing this is to replace much of the code here with
 * a call to the routine model (spectype, par) in models.c
 *
 **********************************************************/


double
one_continuum (spectype, t, g, freqmin, freqmax)
     int spectype;
     double t, g, freqmin, freqmax;
{
  double lambdamin, lambdamax;
  double w_local[NCDF], f_local[NCDF];
  double f, y;
  int n, nwave;
  double par[2];
  int model ();

  /* Check if the parameters are the same as the stored ones, otherwise initialise */
  if (old_t != t || old_g != g || old_freqmin != freqmin || old_freqmax != freqmax)
  {                             /* Then we must initialize */
    if (t != comp[spectype].xmod.par[0] || g != comp[spectype].xmod.par[1])
    {
      par[0] = t;
      par[1] = g;
      model (spectype, par);
      old_t = t;
      old_g = g;
    }

    lambdamin = C * 1e8 / freqmax;
    lambdamax = C * 1e8 / freqmin;
    nwave = 0;

    /* if the first wavelength in the model is below the wavelength range in the simulation,
       interpolate on the model flux to get the flux at lambdamin. copy relevant wavelengths and
       fluxes to w_local and f_local  */
    if (comp[spectype].xmod.w[0] < lambdamin && lambdamin < comp[spectype].xmod.w[comp[spectype].nwaves - 1])
    {
      w_local[nwave] = lambdamin;
      linterp (lambdamin, comp[spectype].xmod.w, comp[spectype].xmod.f, comp[spectype].nwaves, &y, 0);
      f_local[nwave] = y;
      nwave++;
    }

    /* loop over rest of model wavelengths and fluxes and copy to w_local and f_local */
    for (n = 0; n < comp[spectype].nwaves; n++)
    {
      if (comp[spectype].xmod.w[n] > lambdamin && comp[spectype].xmod.w[n] <= lambdamax)
      {
        w_local[nwave] = comp[spectype].xmod.w[n];
        f_local[nwave] = comp[spectype].xmod.f[n];
        nwave++;
      }
    }

    /* now check if upper bound is beyond lambdamax, and if so, interpolate to get appropriate flux
       at lambda max. copy to w_local and f_local */
    if (comp[spectype].xmod.w[0] < lambdamax && lambdamax < comp[spectype].xmod.w[comp[spectype].nwaves - 1])
    {
      w_local[nwave] = lambdamax;
      linterp (lambdamax, comp[spectype].xmod.w, comp[spectype].xmod.f, comp[spectype].nwaves, &y, 0);
      f_local[nwave] = y;
      nwave++;
    }

    /* There are two pathological cases to deal with, when we only have one non zero point,
       we need to make an extra point just up/down from the penultimate/second point so we
       can make a sensible CDF. */

    if (f_local[nwave - 2] == 0.0)      //We have a zero just inside the end
    {
      nwave++;
      w_local[nwave - 1] = w_local[nwave - 2];
      f_local[nwave - 1] = f_local[nwave - 2];
      w_local[nwave - 2] = w_local[nwave - 3] / (1. - DELTA_V / (2. * C));
      linterp (w_local[nwave - 2], comp[spectype].xmod.w, comp[spectype].xmod.f, comp[spectype].nwaves, &y, 0);
      f_local[nwave - 2] = y;
    }

    /* we should now have our arrays w_local and f_local which can be used to generate a cdf */

    /*  Get_model returns wavelengths in Ang and flux in ergs/cm**2/Ang */

    if (cdf_gen_from_array (&comp[spectype].xcdf, w_local, f_local, nwave, lambdamin, lambdamax) != 0)
    {
      Error ("In one_continuum after return from cdf_gen_from_array\n");
    }
    old_freqmin = freqmin;
    old_freqmax = freqmax;
  }

  /* generate the frequency from the CDF that has been built up from the model fluxes */

  f = (C * 1.e8 / cdf_get_rand (&comp[spectype].xcdf));

  /* check if the frequency is too small or too large, and default to simulation limits */
  if (f > freqmax)
  {
    Error ("one_continuum: f too large %e\n");
    f = freqmax;
  }
  if (f < freqmin)
  {
    Error ("one_continuum: f too small %e\n");
    f = freqmin;
  }
  return (f);
}


/**********************************************************/
/** 
 * @brief      get the surface flux for Hubeny/Kurucz like
 * stellar models
 *
 * @param [in] int  spectype   an integer identifying the set of models to use
 * @param [in] double  freqmin   The mimimum frequency
 * @param [int] double  freqmax   The maximum frequency
 * @param [in] double  t   A temperature
 * @param [in] double  g   A gravity
 * @return     The surface flux for a star-like spectrum read from a grid
 *
 * @details
 * The routine gets the band-limited flux per unit area for Hubeny or Kurucz models
 * which are both in units of the Eddingting flux (H)
 *
 * ### Notes ###
 *
 * As indicated this is not a general purpose routine.  Not all spectral
 * models for stars are calculated in terms of the Eddington flux.  The
 * factor of 4 pi is the conversion from Eddington flux to phhysical flux
 *
 * For an explanation see Hubeny & Mihalas - Theory of Stellar Atmospheres,
 * equation 3.70 (or Mihalas, Stellar Atmospheres 2nd ed, eq (1-27).
 *
 * Note that the observed flux at a distance is given by
 *
 * f/H = 4pi * R**2/d**2
 *
 * and
 *
 * L=f*4pi*R**2 = 16 pi**2 R**2 H
 *
 *
 **********************************************************/


double
emittance_continuum (spectype, freqmin, freqmax, t, g)
     int spectype;
     double freqmin, freqmax, t, g;
{
  int nwav, n;
  double w, x, lambdamin, lambdamax;
  double dlambda;
  double par[2];
  int model ();

  lambdamin = C / (freqmax * ANGSTROM);
  lambdamax = C / (freqmin * ANGSTROM);
  par[0] = t;
  par[1] = g;
  model (spectype, par);
  nwav = comp[spectype].nwaves;

  if (lambdamax > comp[spectype].xmod.w[nwav - 1] || lambdamin < comp[spectype].xmod.w[0])
  {

    Error ("emittance_continum: Requested wavelengths extend beyond models wavelengths for list %s\n", comp[spectype].name);
    Error ("lambda %f %f  model %f %f\n", lambdamin, lambdamax, comp[spectype].xmod.w[0], comp[spectype].xmod.w[nwav - 1]);

  }
  x = 0;
  for (n = 0; n < nwav; n++)
  {
    w = comp[spectype].xmod.w[n];
    if (n == 0)
    {
      dlambda = comp[spectype].xmod.w[1] - comp[spectype].xmod.w[0];
    }
    else if (n == nwav - 1)
    {
      dlambda = comp[spectype].xmod.w[n] - comp[spectype].xmod.w[n - 1];
    }
    else
    {
      dlambda = 0.5 * (comp[spectype].xmod.w[n + 1] - comp[spectype].xmod.w[n - 1]);
    }
    if (lambdamin < w && w < lambdamax)
    {
      x += comp[spectype].xmod.f[n] * dlambda;
    }
  }
  x *= 4. * PI;
  return (x);
}
