
/**************************************************************************
                    Space Telescope Science Institute


  Synopsis:
	levels (xplasma, mode) calculates the fractional occupation numbers of
	the various levels of atomic configurations as designated in
	the atomic data files

  Description:

	mode	0	LTE with t_r
		1	LTE with t_e
		2	Non-LTE (reduced by weighted BB)

  Arguments:

  Returns:

  Notes:

	0808 - ksl - levels populates the levden array in the Plasma pointer.  It
		is called from ion_abundances in python, and is called directly
		from my diagnostic routine balance.  It's closely related to
		another but separate routine which calculates the partition
		functions, callled partition

  History:
	01sep23	ksl	Began work
	01oct10	ksl	Modified so modes matched python ionization modes
			more exactly.
	01dec03	ksl	Modified to simplify so modes match those of
			nebular concentrations
	01dec12	ksl	Modified to react to changes which split "nlte"
			and "lte" levels.  Levels is explicitly for so-
			called "nlte" levels, which are tracked in the
			Wind structure
	04Apr   SS      If statement added to avoid this routine changing
                        macro atom level populations.
        04May   SS      The if statment added above is modified for the case
                        of all "simple" ions.
	06may	ksl	57+ -- Modified to make use of plasma structue
	080810	ksl	62 - Fixed problem with how the levden array was
			indexed.  The problem was that the index into the
			levden array is not the same as the index into
			the so called nlte configurations.
			Also made the actual use of variables
			like nion,n,m resemble that in partitions so the
			routine was easier to compae

 ************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "atomic.h"
#include "python.h"

int
levels (xplasma, mode)
     PlasmaPtr xplasma;
     int mode;
{
  double t, weight;
  int n, m;
  int nion;
  int m_ground;
  int nlevden;
  double kt;
  double z;

  if (mode == NEBULARMODE_TR)		// LTE with t_r
    {
      t = xplasma->t_r;
      weight = 1;
    }
  else if (mode == NEBULARMODE_TE)		// LTE with t_e
    {
      t = xplasma->t_e;
      weight = 1;
    }
  else if (mode == NEBULARMODE_ML93)		// non_LTE with t_r and weights
    {
      t = xplasma->t_r;
      weight = xplasma->w;
    }
  else if (mode == 3)		/* non_LTE with SS modification NSH 120912 - This mode is more or less defunct. It can be romoved once all the viestiges of the original PL ioinzation scheme are removed */
    {
      t = xplasma->t_e;
      weight = 1;
    }
  else if (mode == 4)		/* A test mode - this is to allow all levels to be set to GS, in the event we dont have a good idea of what the radiation field shoulb be. */
    {
      t = xplasma->t_e;
      weight = 0;
    }
  else
    {
      Error ("levels: Could not calculate levels for mode %d\n", mode);
      exit (0);
    }

  /* Next calculation should be almost identical to that contained in
     the routine partition.  An error in one place implies one in the
     other.  It seems like that this routine should be subsumed into
     that one ????  */

  kt = BOLTZMANN * t;

  for (nion = 0; nion < nions; nion++)
    {
      if (ion[nion].nlte > 0)
	/* Extra if statement added to prevent changing levden of macro atom populations (SS, Apr04) */
	{
    if (ion[nion].macro_info == 0 || geo.macro_ioniz_mode == 0)
      {     //Then calculate levels for this ion 

        z = xplasma->partition[nion];
        nlevden = ion[nion].first_levden;

        /* N.B. partition functions will most likely have been 
           calculated from "lte" levels, at least for * now ??  */

        get_boltzmann_populations(xplasma->levden, nion, weight, t, z, nlevden);
    
      }
  }
    }

  return (0);

}

/* get_boltzmann_populations calculates populations for ion nion and
   places them in an array levden_array,
   at temperature t, dilution factor w. w=1 is LTE.

   This is used by levels, above, and also by macro-atom level populations
   which calculate departure coefficients and so need LTE populations.

   nlevden will index the levden array directly in the case of levels()
   but may be set to zero if we want to store the populations in a temporary array.
*/


int
get_boltzmann_populations(levden_array, nion, w, t, z, nlevden)
    double levden_array[NLTE_LEVELS];
    int nion, nlevden;
    double w, t, z;
{
  double ground;
  int m, m_ground, n;
  double kt;
  m = ion[nion].first_nlte_level;
  m_ground = m; 
  levden_array[nlevden] = ground = config[m].g / z;  //Assumes first level is ground state

  kt = BOLTZMANN * t;

  for (n = 1; n < ion[nion].nlte; n++)
    {
      m++;
      nlevden++;
      levden_array[nlevden] = ground * w * config[m].g *
        exp ((-config[m].ex + config[m_ground].ex) / kt) / z;
    }

  return (0);
}


int get_lte_matom_populations(levden_array, nelem, xplasma)
        double levden_array[NLTE_LEVELS];
        int nelem;
        PlasmaPtr xplasma;
{
  PlasmaPtr xdummy_lte;
  plasma_dummy pdum;
  double nh, ion_fraction, z;
  int nion, n, nlevden_first;

  nh = xplasma->rho * rho2nh;

  /* create a copy of the plasma pointer and compute LTE populations- 
     we'll use this to get departure coefficients */
  xdummy_lte = &pdum;
  copy_plasma (xplasma, xdummy_lte);

  /* now compute saha ion abundances and LTE partition functions 
     and store in our copy of the plasma pointer */
  partition_functions (xdummy_lte, NEBULARMODE_TR);
  saha (xdummy_lte, xdummy_lte->ne, xdummy_lte->t_r);

  for (nion = ele[nelem].firstion;
      nion <
      (ele[nelem].firstion + ele[nelem].nions);
      nion++)
      {
        nlevden_first = ion[nion].first_levden;

        ion_fraction = xdummy_lte->density[nion] / (nh * ele[nelem].abun);
        z = xdummy_lte->partition[nion];

        /* this gives us population as an ion fraction */
        get_boltzmann_populations(xdummy_lte->levden, nion, 1.0, xdummy_lte->t_r, z, nlevden_first);
        
        /* convert the populations to fraction of entire element by multiplying by ion fraction */
        for (n = nlevden_first; n < nlevden_first+ion[nion].nlte; n++)
          {
            levden_array[n] *= ion_fraction;
          }
      }
  return (0);
}





/**************************************************************************

  Synopsis:  

  Routine to copy the necessary parts of a plasma structure for computing 
  a set of level populations. x1 points to the cell from which data is copied 
  and x2 points to the cell to which data is copied.
  
  Description:  

  Arguments:  

  Returns:

  Notes:

  History:
 ************************************************************************/

int
copy_plasma (x1, x2)
     PlasmaPtr x1, x2;
{
  x2->nwind = x1->nwind;
  x2->nplasma = x1->nplasma;
  x2->ne = x1->ne;
  x2->rho = x1->rho;
  x2->vol = x1->vol;
  x2->t_r = x1->t_r;
  x2->t_e = x1->t_e;
  x2->w = x1->w;

  /* JM 1409 -- added this for depcoef_overview_specific */
  x2->partition = x1->partition;
  x2->density = x1->density;

  /* Note this isn't everything in the cell! 
     Only the things needed for these routines */

  return (0);
}
