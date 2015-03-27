#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "atomic.h"
#include "python.h"

/***********************************************************
                                       University of Southampton

Synopsis:
	matom_emiss_report is a routine which reports the matom level and kpkt emissivities summed
	over all cells. It is called in define_phot() in the spectral cycles, after get_matom_f()

Arguments: 

Returns: 
	simply logs the level emissivites and kpkt emissivites in macro atom mode.
 
Description:	 
	prints out the level emissivities in macro atoms, summed over each cell. 
    	Should only be used in macro atom mode in the spectral cycles 
    	(when geo.matom_radiation = 1). 
	
History:
	130609 Initial coding

**************************************************************/

int
matom_emiss_report ()
{

  int n, m;
  double emiss_sum, abs_sum;

  /* Cycle over macro atom levels and log emissivities */

  for (m = 0; m < nlevels_macro; m++)
    {
      emiss_sum = 0.0;
      abs_sum = 0.0;
      for (n = 0; n < NPLASMA; n++)
	{
	  emiss_sum += macromain[n].matom_emiss[m];
	  abs_sum += macromain[n].matom_abs[m];
	}

      Log
	("Macro Atom level emissivities (summed over cells): n %d matom_abs %8.4e matom_emiss %8.4e\n",
	 m, abs_sum, emiss_sum);
    }



  /* Log kpkt emissivities as well */
  emiss_sum = 0.0;
  abs_sum = 0.0;

  for (n = 0; n < NPLASMA; n++)
    {
      emiss_sum += plasmamain[n].kpkt_emiss;
      abs_sum += plasmamain[n].kpkt_abs;

    }

  Log
    ("Kpkt emissivities (summed over cells): kpkt_abs %8.4e kpkt_emiss %8.4e\n",
     abs_sum, emiss_sum);

  /* Log totals */
  Log ("Totals: f_matom %le f_kpkt %le\n", geo.f_matom, geo.f_kpkt);

  return (0);
}
