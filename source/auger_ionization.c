#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "atomic.h"
#include "python.h"

int
auger_ionization (xplasma)
     PlasmaPtr xplasma;
{

  int n;
  double ionization_rate_coeff, rad_rec_rate_coeff;
  double diel_rec_rate_coeff, factor;

  for (n = 0; n < nauger; n++)
    {
      /* for now assumption for auger effect is very crude:
         1) the upper ion is made only via the Auger effect
         (e.g. OVI is made only by inner shell ionization 
         followed by autoionization of OIV)
         2) The higher ion recombination rate to the ion below 
         is its only important loss mechanism - thus OVI -> OV
         3) The high ion population is small such that it doesn't
         affect the populations of the lower ions in an important
         manner

         subject to these (highly questionable assumptions) the upper
         ion population is easy to compute from the rate at which
         inner shell ionizations occur (the MC estimator), the
         recombination rate coefficient, the electron density and the
         fixed population of the lower ion */

      ionization_rate_coeff = xplasma->gamma_inshl[n] * augerion[n].yield;
      rad_rec_rate_coeff =
	augerion[n].arad * pow (xplasma->t_e / 10000.,
				(-1. * augerion[n].etarad));
      diel_rec_rate_coeff =
	augerion[n].adi * pow ((xplasma->t_e),
			       -1.5) * exp (-1. * augerion[n].t0di /
					    xplasma->t_e) * (1. +
							     (augerion[n].bdi
							      * exp (-1. *
								     augerion
								     [n].t1di
								     /
								     xplasma->
								     t_e)));


      if ((augerion[n].nion_target != -1)
	  && (xplasma->density[augerion[n].nion] > DENSITY_MIN))
	{
	  if ((factor =
	       (ionization_rate_coeff /
		(rad_rec_rate_coeff + diel_rec_rate_coeff) / xplasma->ne)) <
	      1.0)
	    {
	      //Log("Auger_ionization: overwriting ion %d population %g\n", augerion[n].nion_target, xplasma->density[augerion[n].nion_target]);
	      //Log("Changing by factor of: %g\n",  xplasma->density[augerion[n].nion_target] / (xplasma->density[augerion[n].nion] * ionization_rate_coeff / (rad_rec_rate_coeff + diel_rec_rate_coeff) / xplasma->ne));
	      xplasma->density[augerion[n].nion_target] =
		xplasma->density[augerion[n].nion] * factor;
	      //Log("New value: %g %g\n", xplasma->density[augerion[n].nion_target], xplasma->density[augerion[n].nion]);
	    }
	  else
	    {
	      Error
		("Auger_ionization trying to set large upper population. Setting equal.\n");
	      xplasma->density[augerion[n].nion_target] =
		xplasma->density[augerion[n].nion];
	    }

	}
    }

  return (0);
}
