
/***********************************************************/
/** @file  reposition.c
 * @author ksl
 * @date   May, 2018
 *
 * @brief  Reposition a photon so it goes across a cell boundary
 *
 * ### Programming Notes ###
 *
 * This is a very short function and may not justify being in its own file. The
 * (only) function reposition could be moved into extract.c or trans_phot.c.
 *
 ***********************************************************/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "atomic.h"
#include "python.h"


/**********************************************************/
/**
 * @brief      p) attempts to assure that a photon is not scattered
 * 	a second time inappropriately by the same transition
 *
 * @param [in,out] PhotPtr  p   A photons
 * @return    Normally returns 0, but returns a negative number
 * if p is not in the wind in the domain it is supposed to be in
 *
 * @details
 * For resonant scatters, the routine moves the photon by
 * a distance dfudge. For non-resonant scattering the routine
 * simply returns.
 *
 * ### Notes ###
 *
 **********************************************************/

int
reposition (p)
     PhotPtr p;
{
  int n;

  if (p->nres < 0)
    return (0);                 /* Do nothing for non-resonant scatters */

  if ((p->grid = n = where_in_grid (wmain[p->grid].ndom, p->x)) < 0)
  {
    Error ("reposition: Photon not in grid when routine entered %d \n", n);
    return (n);                 /* Photon was not in wind */
  }

  /*
   * EP 0818: the cell specific dfudge is now used as previously the global
   * value of DFUDGE was being used
   */
  move_phot (p, wmain[p->grid].dfudge);

  return (0);
}
