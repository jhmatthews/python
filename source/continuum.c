#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "atomic.h"
#include "python.h"

#include "models.h"
//




/***********************************************************
                                       Space Telescope Science Institute

 Synopsis:

double one_continuum(spectype,t,g,freqmin,freqmax) gets a photon frequency 
from a continuum grid 

Arguments:

	spectype				An index to the continum grid which one wants
						to use for obtaining a photon
	double t,g,				Two parameters defining the model for which 
						a photon is requried, often temperature and 
						gravity, though this is not required
	freqmin,freqmax;			minimum and maximum frequency of interest


Returns:

 
Description:	
	

Notes:
	In versions of python prior to python_52, different routines
	were used to read kurucz models and hubeny models.  This was
	historical, since the routines ksl had to read the kurucz models
	read binary files (from a time when compputers were very slow).

	The new routines are much more general, and are based on routines
	ksl had written for his fitting program kslfit.  The underlying
	routines can interpolate between grids of more than two dimensions
	although here we assume that all of the continuum grids are two
	dimensional.  

History:
	04aug	ksl	created from hub.c usied in all versions
			of python prior to python_52.  
 
**************************************************************/



double old_t, old_g, old_freqmin, old_freqmax;
double jump[] = { 913.8 };

double
one_continuum (spectype, t, g, freqmin, freqmax)
     int spectype;
     double t, g, freqmin, freqmax;
{
  double par[2];		// For python we assume only two parameter models
  double lambdamin, lambdamax;
  double f;
  double pdf_get_rand ();
  int model (), nwav;

  if (old_t != t || old_g != g || old_freqmin != freqmin
      || old_freqmax != freqmax)
    {				/* Then we must initialize */
      par[0] = t;
      par[1] = g;
      nwav = model (spectype, par);
      /*  Get_model returns wavelengths in Ang and flux in ergs/cm**2/Ang */
      lambdamin = C * 1e8 / freqmax;
      lambdamax = C * 1e8 / freqmin;
      if (pdf_gen_from_array
	  (&comp[spectype].xpdf, comp[spectype].xmod.w, comp[spectype].xmod.f,
	   comp[spectype].nwaves, lambdamin, lambdamax, 1, jump) != 0)
	{
	  Error ("In one_continuum after return from pdf_gen_from_array\n");
	}
      old_t = t;
      old_g = g;
      old_freqmin = freqmin;
      old_freqmax = freqmax;
    }

  f = (C * 1.e8 / pdf_get_rand (&comp[spectype].xpdf));
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



double
emittance_continuum (spectype, freqmin, freqmax, t, g)
     int spectype;
     double freqmin, freqmax, t, g;
{
  int model (), nwav, n;
  double w, x, lambdamin, lambdamax;
  double dlambda;
  double par[2];

  lambdamin = C / (freqmax * ANGSTROM);
  lambdamax = C / (freqmin * ANGSTROM);
  par[0] = t;
  par[1] = g;
  nwav = model (spectype, par);

  if (lambdamax > comp[spectype].xmod.w[nwav - 1]
      || lambdamin < comp[spectype].xmod.w[0])
    {

      Error
	("emittance_continum: Requested wavelengths extend beyond models wavelengths for list %s\n",
	 comp[spectype].name);
      Error ("lambda %f %f  model %f %f\n", lambdamin, lambdamax,
	     comp[spectype].xmod.w[0], comp[spectype].xmod.w[nwav - 1]);

//      exit (0);
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
	  dlambda =
	    0.5 * (comp[spectype].xmod.w[n + 1] -
		   comp[spectype].xmod.w[n - 1]);
	}
      if (lambdamin < w && w < lambdamax)
	{
	  x += comp[spectype].xmod.f[n] * dlambda;
	}
    }
  x *= 4. * PI;
  return (x);
}
