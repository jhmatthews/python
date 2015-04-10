#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "atomic.h"
#include "python.h"
/*****************************************************************************

                        University of Southampton
Synopsis:
  choose_superlevel_deactivation decides which of the levels
  in the superlevel to jump down from

Arguments:
  xplasma 		Plasma Pointer to the cell/macro-atom

  nion 			number of level, index to config structure
  
Returns:
  n 			the new level, index to config structure

Description:
  choose_superlevel_deactivation decides which of the levels
  in the superlevel to jump down from. It does this by standard 
  RNG processes, using the quantities computed in setup_superlevels.

  It is designed for a speed improvement when upper levels in a macro-atom are approaching LTE.

Notes: 

History:

************************************************************/

int choose_superlevel_deactivation (xplasma, uplvl)
	PlasmaPtr xplasma;
	int uplvl;
{
  double threshold, run_tot, z;
  int n, nion, ground;
  MacroPtr mplasma;

  /* get some cell and ion properties */
  mplasma = &macromain[xplasma->nplasma];
  nion = config[uplvl].nion;

  /* generate random numbers and loop over all superlevel levels */
  run_tot = 0.0;
  z = ((rand () + 0.5) / MAXRAND);
  threshold = z * mplasma->superlevel_norm[nion];

  n = mplasma->superlevel_threshold[nion];
  ground = ion[nion].first_nlte_level;

  while (run_tot < threshold && n < ground + ion[nion].nlte)
	{ 
	  /* keep looping til we hit the random number. note division by g */
	  run_tot += mplasma->superlevel_lte_pops[n] / config[n].g;
	  n++;
	}
  
  /* we added on more, so take one off unless we were already at the threshold */
  if (n > mplasma->superlevel_threshold[nion])
    n--;
  
  /* throw some error messages if something has gone wrong */
  if (n == ground + ion[nion].nlte + 1)
  	Error("choose_superlevel_deactivation: failed to choose deactivation! %8.4e %8.4e %8.4e\n", 
  		   run_tot, threshold, mplasma->superlevel_threshold[nion], z);

  if (n < mplasma->superlevel_threshold[nion])
  {
  	Error("choose_superlevel_deactivation: level %i is not in superlevel!\n", n);
  }

  return n;
}

/*****************************************************************************

                        University of Southampton
Synopsis:
  setup_superlevels initialises the superlevel properties in each macro-atom

Arguments:
  
Returns:
  0 	Success

Description:
  setup_superlevels loops over all ions- if the ion has superlevels, as
  decided by the flag in the atomic data then we loop over all cells 
  and calculate lte populations and calculate the superlevel threshold.

Notes: 
  note that 

History:

************************************************************/

int setup_superlevels()
{
  double te, kt;
  int nplasma, n, ground, nion;
  MacroPtr mplasma;
  PlasmaPtr xplasma;


  for (nion = 0; nion < NIONS; nion++)
  {
    if (ion[nion].has_superlevel == 1)
    {
      for (nplasma = 0; nplasma < NPLASMA; nplasma++)
      {
      	/* get the macro and plasma pointers for this cell */
        xplasma = &plasmamain[nplasma];
  	    mplasma = &macromain[nplasma];

  	    /* temps for LTE calculation */
        te = xplasma->t_e;				// LTE at TE is what we want
  	    kt = BOLTZMANN * te;
        
        /* get the indexes for the ground. We define relative to ground,
           so the first entry is 1 */
        ground = ion[nion].first_nlte_level;
        mplasma->superlevel_lte_pops[ground] = 1.0;

        /* now calculate the LTE ratios to ground for the other elements */
        for (n = ground + 1; n < ground + ion[nion].nlte; n++)
  	      {
            mplasma->superlevel_lte_pops[n] = (config[n].g / config[ground].g) * exp ((-config[n].ex + config[ground].ex) / kt);
  	      }
        
        /* get the threshold */
        mplasma->superlevel_threshold[nion] = get_superlevel_threshold (xplasma, nion);	// this should be improved to do an on the fly calculation 

        /* the normalisation is the sum of the populations over the levels above the threshold.
           divided by the statistical weight */
        for (n = mplasma->superlevel_threshold[nion]; n < ground + ion[nion].nlte; n++)
  	      {
  	      	mplasma->superlevel_norm[nion] += mplasma->superlevel_lte_pops[n] / config[n].g;
  	      }

      }

    }
  }

  return 0;
}

/*****************************************************************************

                        University of Southampton
Synopsis:
  get_superlevel_threshold calculates the level above which we will treat
  levels as an LTE superlevel. 

Arguments:
  xplasma 		Plasma Pointer to the cell/macro-atom

  nion 			number of the ion in the ion structure
  
Returns:
  threshold 	the superlevel threshold

Description:
   The approach here is to compare the level populations from the last cycle to
  the LTE solution. If they are within a factor LTE_DEP_FRAC then we treat
  them as part of a superlevel.

Notes: 

History:
  1503 JM Coded

************************************************************/
#define LTE_DEP_FRAC 2 

int get_superlevel_threshold(xplasma, nion)
	PlasmaPtr xplasma;
	int nion;
{
  MacroPtr mplasma;
  int ground, threshold, lastden, n;
  double dep_coef, ground_den;

  ground = ion[nion].first_nlte_level;
  lastden = ion[nion].first_levden + ion[nion].nlte - 1;
  threshold = ion[nion].nlte + ground - 1;

  /* if this is the first cycle then the threshold should be the last level */
  if (geo.wcycle == 0)
   return (threshold);

  mplasma = &macromain[xplasma->nplasma];

  ground_den = xplasma->levden[ion[nion].first_levden];

  dep_coef = mplasma->superlevel_lte_pops[threshold];
  dep_coef /= (xplasma->levden[lastden] / ground_den);
  //Log ("JM Cell %i ne %8.4e pops0 %8.4e %8.4e", xplasma->nplasma, xplasma->ne, mplasma->superlevel_lte_pops[threshold], (xplasma->levden[lastden] / ground_den));
  //Log (" deps %8.4e ", dep_coef);

  n = 0;

  while (dep_coef < LTE_DEP_FRAC && dep_coef > (1./LTE_DEP_FRAC) && (threshold - ground) > LOWEST_SUPERLEVEL_THRESHOLD)
  {
  	threshold--;;
  	n++;
  	dep_coef = mplasma->superlevel_lte_pops[threshold];
    dep_coef /= xplasma->levden[lastden - n] / ground_den;
  }


  Log ("JM threshold %i ground %i upper ion ground %i\n", threshold, ground, ion[nion+1].first_nlte_level);
  
  threshold++;

  return threshold;
}