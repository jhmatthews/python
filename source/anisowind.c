
/***********************************************************/
/** @file  anisowind.c
 * @author ksl
 * @date   March, 2018
 *
 * @brief  Routines to implement anisotropic scattering in
 * the wind
 *
 *
 * Python supports two ways to determine a new photon
 * direction when a photon scatters in thie wind.  These
 * include 
 *
 * * SCATTER_MODE_ISOTROPIC -isotropic scattering (randvec),
 * * SCATTER_MODE_THERMAL - thermally-broadened anisotropic scattering
 * (randwind_thermal_trapping).
 *
 ***********************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "atomic.h"
#include "python.h"

/**********************************************************/
/**
 * @brief the routine which chooses
 *   a new anisotropic direction in geo.scatter_mode = SCATTER_MODE_THERMAL
 *
 * @param [in,out] PhotPtr  p   The photon being scattered
 * @param [out] int *  nnscat   The number of times the phton
 * scattered internally before escaping the local scattering region
 *
 * @return     0 for success. Also modifies the photon ptr p
 *   to reflect new direction (p->lmn), and nnscat, which
 *   should be copied to the phoiton structure after calling
 *   this routine.
 *
 * @details
 * This routine uses a rejection method to choose a direction
 * so that the probability distribution of directions generated
 * reflects the probability of escape along each direction in accordance
 * with the sobolev optical depth.
 *
 * ### Notes ###
 * The resonance that caused the scatter must be stored in the
 * photon bundle.
 *
 * The name of the routine is something of a misnomer. Pure
 * sobolev optical depths are used.  The temperature in the
 * cell does not come into the calculation.
 *
 * 1904 - ksl - In certain situations, see issue #505, the photon can
 * get trapped in the while loop.  To prevent his from hanging
 * the routine bails after NSCAT_MAX trips through the while loop
 * and chooses the direction with the lowest tau it has found. 
 * From a physics perspective, this is not an ideal treatment of
 * the problem ans so if this occurs a lot of times, something
 * better needs to be done.
 **********************************************************/

#define NSCAT_MAX 10000

int
randwind_thermal_trapping (p, nnscat)
     PhotPtr p;
     int *nnscat;
{
  double tau_norm, p_norm;
  double tau, dvds, z, ztest;
  double z_prime[3], z_min[3];
  int nscat;
  double tau_min = VERY_BIG;
  WindPtr one;

  /* find the wind pointer for the photon */
  one = &wmain[p->grid];

  /* we want to normalise our rejection method by the escape
     probability along the vector of maximum velocity gradient.
     First find the sobolev optical depth along that vector */
  tau_norm = sobolev (one, p->x, -1.0, lin_ptr[p->nres], one->dvds_max);

  /* then turn into a probability. Note that we take account of
     this in trans_phot before calling extract */
  p_norm = p_escape_from_tau (tau_norm);

  /* Throw error if p_norm is 0 */
  if (p_norm <= 0)
    Error ("randwind_thermal_trapping: p_norm is %8.4e in cell %i", p_norm, one->nplasma);

  ztest = 1.0;
  z = 0.0;

  /* JM 1406 -- we increment nnscat here, and it is recorded in the photon
     structure. This is done because we actuall have to multiply the photon weight
     by 1/mean escape probability- which is nnscat. this is done in trans_phot.c
     before extract is called.
   */
  *nnscat = *nnscat - 1;
  nscat = 0;

  /* rejection method loop, which chooses direction and also calculated nnscat */
  while (ztest > z && nscat < NSCAT_MAX)
  {
    *nnscat = *nnscat + 1;
    randvec (z_prime, 1.0);     /* Get a new direction for the photon (isotropic */
    stuff_v (z_prime, p->lmn);  // copy to photon pointer


    /* generate random number, normalised by p_norm with a 1.2 for 20%
       safety net (as dvds_max is worked out with a sample of directions) */
    ztest = random_number (0.0, 1.0) * p_norm;

    dvds = dvwind_ds (p);
    tau = sobolev (one, p->x, -1.0, lin_ptr[p->nres], dvds);

    if (tau < tau_min)
    {
      stuff_v (z_prime, z_min);
      tau_min = tau;
    }

    z = p_escape_from_tau (tau);        /* probability to see if it escapes in that direction */
    nscat++;
  }

  if (nscat == NSCAT_MAX)
  {
    stuff_v (z_min, p->lmn);    // copy to photon pointer
    Error ("rand_wind_thermal_trapping: photon had more than  %d internal scatters\n", nscat);
  }

  return (0);
}
