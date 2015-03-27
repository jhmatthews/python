#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "atomic.h"
#include "python.h"
#include "recipes.h"



/***********************************************************
                                       Space Telescope Science Institute
 
 Synopsis:
	detailed_balance (xplasma, nelem, newden)
 
 Arguments:
 
Returns:
  
Description:
 
Notes:
 
History:
        06may   ksl     57+--Modified to use new plasma structuref level populations.  
**************************************************************/

int
detailed_balance (xplasma, nelem, newden)
     PlasmaPtr xplasma;
     int nelem;
     double newden[];
{
  double rates_up[NIONS], rates_down[NIONS], fraction[NIONS];
  int n;
  int first, last;
  double tot;


  first = ele[nelem].firstion;	/*first and last identify the postion in the array */
  last = first + ele[nelem].nions;	/*  So for H which has 2 ions, H1 and H2, first will generally
					   be 0 and last will be 2 so the for loop below will just be done once for nion = 1 */


  tot = 0;
  for (n = first; n < last; n++)
    {
      tot += xplasma->density[n];
      if (xplasma->density[n] > 0.0)
	{
	  rates_up[n] = xplasma->ioniz[n] / xplasma->density[n];
	}
      else
	rates_up[n] = 0.0;
      if (xplasma->density[n + 1] > 0.0)
	{
	  rates_down[n] = xplasma->recomb[n] / xplasma->density[n + 1];
	}
      else
	rates_down[n] = 0.0;

      if (rates_up[n] < 1.e-50 || rates_down[n] < 1.e-50)
	{
	  rates_up[n] = 1e-10 * xplasma->density[n + 1];
	  rates_down[n] = 1e-10 * xplasma->density[n];
	}
    }

  rebalance (&rates_up[first], &rates_down[first], &fraction[first],
	     ele[nelem].nions);
  for (n = first; n < last; n++)
    {
      newden[n] = fraction[n] * tot;
      //Prevent the density from falling to a number likely to produce zero divides
      if (newden[n] < DENSITY_MIN)
	newden[n] = DENSITY_MIN;
    }

  return (0);
}

/* This subroutine attempts to calculate the ionization fractions
 * for a single ion in a situation in which it's assumed thatk
 * we know the rate at which ions are ionized and recombine and in
 * which all ionization occurs from the state below the current one
 * and all recombinations are from the state above
*/

int
rebalance (rates_up, rates_down, fraction, ntot)
     double rates_up[];
     double rates_down[];
     double fraction[];
     int ntot;
{
  int n;
  double sum;
  fraction[0] = sum = 1.;
  for (n = 1; n < ntot; n++)
    {

      fraction[n] = rates_up[n - 1] / rates_down[n - 1] * fraction[n - 1];
      sum += fraction[n];
    }

  for (n = 0; n < ntot; n++)
    fraction[n] /= sum;
  return (0);
}

/* After the density of ions have been changed by detailed_balance, there are
 * other parameters, e.g. the heating and cooling, the number of recombinations
 * etc. that would be expected to change with the densities assuming the same
 * photon spectrum as last passed through the wind cell.  This subroutine 
 * attempts to modify other portions of the wind structure to reflect the
 * change in ionization
 *
Notes: 
	This subroutine has to be called before the densities are updated. 

	The routine has to be called on an element by element basis to reflect
	its usage in dlucy, and also to ultimately allow one to choose which
	elements to treat in detailed balance.

	??? Initially only worried about the ionization and recombination values
	since this control the densities, but should also update the heating
	and cooling, or at least that is what I expect. ksl
	
History:
	02jul	ksl	Began coding
	06may	ksl	57+ -- Converted to use plasma structure
 */

int
wind_update_after_detailed_balance (xplasma, nelem, newden)
     PlasmaPtr xplasma;
     int nelem;
     double newden[];
{
  int nion;
  int first, last;
  first = ele[nelem].firstion;	/*first and last identify the postion in the array */
  last = first + ele[nelem].nions;	/*  So for H which has 2 ions, H1 and H2, first will generally
					   be 0 and last will be 2 so the for loop below will just be done once for nion = 1 */
  for (nion = first; nion < last; nion++)
    {
      if (xplasma->density[nion] > 0.0)
	{
	  xplasma->ioniz[nion] *= newden[nion] / xplasma->density[nion];
	  xplasma->heat_ion[nion] *= newden[nion] / xplasma->density[nion];
	}
      else
	{
	  Error
	    ("wind_update_after_detailed_balance: ioniz  Old den of nion %d (%s %d) not positive %e \n",
	     nion, ele[nelem].name, ion[nion].istate, xplasma->density[nion]);
	}
      if (xplasma->density[nion + 1] > 0.0)
	{
	  xplasma->recomb[nion] *=
	    newden[nion + 1] / xplasma->density[nion + 1];
	  xplasma->heat_ion[nion + 1] *=
	    newden[nion + 1] / xplasma->density[nion + 1];
	}

      else
	{
	  Error
	    ("wind_update_after_detailed_balance: recomb Old den of nion %d (%s %d) not positive %e \n",
	     nion + 1, ele[nelem].name, ion[nion + 1].istate,
	     xplasma->density[nion + 1]);
	}
    }

  // Next is a kluge for now to preserve heating and cooling in original notation
  // It should be deleted eventually, but one has to make sure none of the functionality
  // is lost

  if (nelem == 0)
    {

      xplasma->heat_photo =
	xplasma->heat_ion[0] + xplasma->heat_ion[2] + xplasma->heat_ion[3] +
	xplasma->heat_z;
      xplasma->heat_tot =
	xplasma->heat_lines + xplasma->heat_ff + xplasma->heat_photo;
      xplasma->lum_fb =
	xplasma->lum_ion[0] + xplasma->lum_ion[2] + xplasma->lum_ion[3] +
	xplasma->lum_z;
      xplasma->lum_rad =
	xplasma->lum_lines + xplasma->lum_ff + xplasma->lum_fb;
    }
  else if (nelem == 1)
    {
      xplasma->lum_ion[2] *= newden[3] / xplasma->density[3];
      xplasma->lum_ion[3] *= newden[4] / xplasma->density[4];
      xplasma->heat_photo =
	xplasma->heat_ion[1] + xplasma->heat_ion[2] + xplasma->heat_ion[3] +
	xplasma->heat_z;
      xplasma->heat_tot =
	xplasma->heat_lines + xplasma->heat_ff + xplasma->heat_photo;
      xplasma->lum_fb =
	xplasma->lum_ion[0] + xplasma->lum_ion[2] + xplasma->lum_ion[3] +
	xplasma->lum_z;
      xplasma->lum_rad =
	xplasma->lum_lines + xplasma->lum_ff + xplasma->lum_fb;
    }
  return (0);
}
