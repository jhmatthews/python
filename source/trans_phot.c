
/***********************************************************/
/** @file  trans_phot.c
 * @author ksl
 * @date   April, 2018
 *
 * @brief  This file contains high level routines which carry
 * out the propagation of photons through the wind/wind domains
 * both for ionization and spectrum (extract) cycles
 *
 * ### Notes ###
 *
 * The routines contained here are central to Python, and anyone
 * who wants to understand Python in general should spend time
 * understanding how they work
 *
 * There are two basic options associated with the routines.
 *
 * * live or die is the case in which one is only concerned
 * with how a photon goes through the wind and random numbers
 * are used to generate the scattering directions.  The photon is
 * followed until it leaves the system.  This is the only option
 * used in ionization cycles (currently).
 *
 * * extract is the mode where whenever a photon scatters we
 * also calculate what the photon weight would be if it scattered
 * in a set of certain directions.  When this happens the weight
 * of the photon is changed to reflect the fact that the direction
 * is not random.
 *
 * The extract option is used
 * normaly during the spectral extraction cycles.
 * However, as an advanced option one can use the live or die
 * to construct the detailed spectrum.  One would not normally
 * want to do this, as many photons are "wasted" since they
 * don't scatter at the desired angle.  With sufficient numbers
 * of photons however the results of the two methods should
 * (by construction) be identical (or at least very very
 * similar).
 *
 *
 ***********************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include "atomic.h"
#include "python.h"


long n_lost_to_dfudge = 0;



/**********************************************************/
/**
 * @brief      int (w,p,iextract) oversees the propagation of a "flight" of photons
 *
 * @param [in] WindPtr  w   The entire wind domain
 * @param [in, out] PhotPtr  p   A pointer to a "fligh" of photons
 * @param [in] int  iextract   An integer controlling whether we are to process the
 * flight in the live or die option (0) or whether we also need to extract photons in
 * specific directions (which is usually the case in constructing spectra
 * @return   Normally returns 0
 *
 * @details
 * This routine oversees the propagation of  individual photons.  The main loop
 * covers an entire "flight" of photons.   The routine generates the random
 * optical depth a photon can travel before scattering and monitors the
 * progress of the photon through the grid.
 *
 * The real physics is done elsewhere, in lower level routines.
 *
 * ### Notes ###
 *
 * At the end of the routine the position for each of the photons in p is the
 * last point where the photon was in the wind, * not the outer boundary of
 * the radiative transfer
 *
 **********************************************************/

int
trans_phot (WindPtr w, PhotPtr p, int iextract)
{
  int nphot;
  struct photon pp, pextract;
  int absorb_reflect;           /* this is a variable used to store geo.absorb_reflect during exxtract */
  double p_norm, tau_norm;
  int nreport;
  struct timeval timer_t0;

  nreport = 100000;
  if (nreport < NPHOT / 100)
  {
    nreport = NPHOT / 100;
  }


  Log ("\n");

  timer_t0 = init_timer_t0 ();

  /* Beginning of loop over photons */

  for (nphot = 0; nphot < NPHOT; nphot++)
  {
    /* This is just a watchdog method to tell the user the program is still running */

    if (nphot % nreport == 0)
    {
      Log ("Cycle %d/%d of %s : Photon %10d of %10d or %6.1f per cent \n", geo.wcycle, geo.pcycle, basename, nphot, NPHOT,
           nphot * 100. / NPHOT);
    }

    Log_flush ();

    /* Verify that the weights are real, a check that is proably unnecessary */

    if (sane_check (p[nphot].w))
    {
      Error ("trans_phot:sane_check photon %d has weight %e\n", nphot, p[nphot].w);
    }

    stuff_phot (&p[nphot], &pp);
    absorb_reflect = geo.absorb_reflect;

    /* The next if statement is executed if we are calculating the detailed spectrum and makes sure we always run extract on
       the original photon no matter where it was generated */

    if (iextract)
    {
      // SS - for reflecting disk have to make sure disk photons are only extracted once.  Note we restore the
      // correct geo.absorb_reflect value as soon as the photons are extracted!

      if (absorb_reflect == BACK_RAD_SCATTER && p[nphot].origin == PTYPE_DISK)
      {
        geo.absorb_reflect = BACK_RAD_ABSORB_AND_DESTROY;
      }


      stuff_phot (&p[nphot], &pextract);


      /* We increase weight to account for number of scatters. This is done because in extract we multiply by the escape
         probability along a given direction, but we also need to divide the weight by the mean escape probability, which is
         equal to 1/nnscat */
      if (geo.scatter_mode == SCATTER_MODE_THERMAL && pextract.nres <= NLINES && pextract.nres > -1)
      {
        /* we normalised our rejection method by the escape probability along the vector of maximum velocity gradient.
           First find the sobolev optical depth along that vector. The -1 enforces calculation of the ion density */

        tau_norm = sobolev (&wmain[pextract.grid], pextract.x, -1.0, lin_ptr[pextract.nres], wmain[pextract.grid].dvds_max);

        /* then turn into a probability */
        p_norm = p_escape_from_tau (tau_norm);

      }
      else
      {
        p_norm = 1.0;

        /* throw an error if nnscat does not equal 1 */
        if (pextract.nnscat != 1)
          Error
            ("trans_phot: nnscat is %i for photon %i in scatter mode %i! nres %i NLINES %i\n",
             pextract.nnscat, nphot, geo.scatter_mode, pextract.nres, NLINES);
      }



      /* We then increase weight to account for number of scatters. This is done because in extract we multiply by the escape
         probability along a given direction, but we also need to divide the weight by the mean escape probability, which is
         equal to 1/nnscat */

      pextract.w *= p[nphot].nnscat / p_norm;

      if (sane_check (pextract.w))
      {
        Error ("trans_phot: sane_check photon %d has weight %e before extract\n", nphot, pextract.w);
      }
      extract (w, &pextract, pextract.origin);


      /* Restore the correct disk illumination */
      if (absorb_reflect == BACK_RAD_SCATTER && p[nphot].origin == PTYPE_DISK)
      {
        geo.absorb_reflect = BACK_RAD_SCATTER;
      }
    }                           /* End of extract loop */

    p[nphot].np = nphot;

    /* Transport a single photon */
    trans_phot_single (w, &p[nphot], iextract);

  }

  /* This is the end of the loop over all of the photons; after this the routine returns */

  /* Line to complete watchdog timer */
  Log ("\n");

  print_timer_duration ("!!python: photon transport completed in", timer_t0);

  /* sometimes photons scatter near the edge of the wind and get pushed out by DFUDGE. We record these */
  if (n_lost_to_dfudge > 0)
    Error
      ("trans_phot: %ld photons were lost due to DFUDGE (=%8.4e) pushing them outside of the wind after scatter\n",
       n_lost_to_dfudge, DFUDGE);

  n_lost_to_dfudge = 0;         // reset the counter

  return (0);
}






/**********************************************************/
/**
 * @brief      Transport a single photon photon through the wind.
 *
 * @param [in] WindPtr  w   The entire wind
 * @param [in, out] PhotPtr  p   A single photon
 * @param [in] int  iextract   If 0, then process this photon in the live or die option, without
 * calling extract
 *
 * @return     Always returns 0
 *
 * @details
 * This routine oversees the propagation of an individual photon through
 * the wind.  As the photon moves through the wind, its position, direction
 * and weight are updated.  In reverberation mode, the flight time is also
 * tracked.
 *
 * Basically what the routine does is generate a random number whihc is used to
 * determine the optical depth to a scatter, and then it calles translate
 * multiple times.   translate involves moving the photon only a single cell
 * (or alternatively a single tranfer in the windless region), and returns
 * a status.  Depending on what this status is, trans_phot_single calls
 * trans_phot again doing nothing, but if the scattering depth has been
 * reached, then trans_phot_single causes the photon to scatter,
 * which changes its direction.  This process continues until the photon
 * exits the system or hits a barrier.  If the photon hits a radiating
 * surface, the disk or star, then the photon may either be absorbed,
 * or scattered depending on the reflection/absorption mode.
 *
 *
 *
 * ### Notes ###
 *
 * This routine is called by trans_phot once for each photon in a flight of photons
 *
 *
 **********************************************************/

int
trans_phot_single (WindPtr w, PhotPtr p, int iextract)
{
  double tau_scat, tau;
  int istat;
  double rrr;
  int icell;
  int nres, *ptr_nres;
  int kkk, n;
  double weight_min;
  struct photon pp, pextract;
  int nnscat;
  int nerr;
  double p_norm, tau_norm;
  double x_dfudge_check[3];
  int ndom;
  double normal[3];

  /* Initialize parameters that are needed for the flight of the photon through the wind */

  stuff_phot (p, &pp);
  tau_scat = -log (1. - random_number (0.0, 1.0));

  weight_min = EPSILON * pp.w;
  istat = P_INWIND;
  tau = 0;
  icell = 0;

  n = 0;                        /* Avoid 03 warning */



  /* This is the beginning of the loop for a single photon and executes until the photon leaves the wind */

  while (istat == P_INWIND)
  {

    /* translate involves only a single shell (or alternatively a single tranfer in the windless region). istat as returned by
       should either be 0 in which case the photon hit the other side of the shell without scattering or 1 in which case there
       was a scattering event in the shell, 2 in which case the photon reached the outside edge of the grid and escaped, 3 in
       which case it reach the inner edge and was reabsorbed. If the photon escapes then we leave the photon at the position
       of it's last scatter.  In most other cases though we store the final position of the photon. */


    istat = translate (w, &pp, tau_scat, &tau, &nres);
    /* nres is the resonance at which the photon was stopped.  At present the same value is also stored in pp->nres, but I have
       not yet eliminated it from translate. ?? 02jan ksl */

    icell++;
    istat = walls (&pp, p, normal);
    /* pp is where the photon is going, p is where it was  */


    if (istat == -1)
    {
      Error_silent ("trans_phot: Abnormal return from translate on photon %d\n", p->np);
      break;
    }

    if (pp.w < weight_min)
    {
      istat = pp.istat = P_ABSORB;      /* This photon was absorbed by continuum opacity within the wind */
      pp.tau = VERY_BIG;
      stuff_phot (&pp, p);
      break;
    }

    if (istat == P_HIT_STAR)
    {                           /* It hit the star */
      geo.lum_star_back += pp.w;
      if (geo.absorb_reflect == BACK_RAD_SCATTER)
      {
        /* If we got here, the a new photon direction needs to be defined that will cause the photon
         * to continue in the wind.  Since this is effectively a scattering event we also have to
         * extract a photon to consturct the detailed spectrum
         */
        randvcos (pp.lmn, normal);
        move_phot (&pp, DFUDGE);
        stuff_phot (&pp, p);
        tau_scat = -log (1. - random_number (0.0, 1.0));
        istat = pp.istat = P_INWIND;    /* Set the status back to P_INWIND so the photon will continue */
        tau = 0;
        if (iextract)
        {
          stuff_phot (&pp, &pextract);
          extract (w, &pextract, PTYPE_STAR);   // Treat as stellar photon for purpose of extraction
        }
      }
      else
      {                         /*Photons that hit the star are simply absorbed so this is the end of the line for this photon */
        stuff_phot (&pp, p);
        break;
      }
    }

    if (istat == P_HIT_DISK)
    {
      /* It hit the disk */

      /* Store the energy of the photon bundle into a disk structure so that one can determine later how much and where the
         disk was heated by photons.
         Note that the disk is defined from 0 to NRINGS-2. NRINGS-1 contains the position of the outer radius of the disk. */

      rrr = sqrt (dot (pp.x, pp.x));
      kkk = 0;
      while (rrr > qdisk.r[kkk] && kkk < NRINGS - 1)
        kkk++;
      kkk--;                    /* So that the heating refers to the heating between kkk and kkk+1 */
      qdisk.nhit[kkk]++;
      geo.lum_disk_back = qdisk.heat[kkk] += pp.w;
      qdisk.ave_freq[kkk] += pp.w * pp.freq;

      if (geo.absorb_reflect == BACK_RAD_SCATTER)
      {
        /* If we got here, the a new photon direction needs to be defined that will cause the photon
         * to continue in the wind.  Since this is effectively a scattering event we also have to
         * extract a photon to consturct the detailed spectrum
         */
        randvcos (pp.lmn, normal);
        stuff_phot (&pp, p);
        tau_scat = -log (1. - random_number (0.0, 1.0));
        istat = pp.istat = P_INWIND;
        tau = 0;
        if (iextract)
        {
          stuff_phot (&pp, &pextract);
          extract (w, &pextract, PTYPE_DISK);
        }
      }
      else
      {                         /* Photons that hit the disk are to be absorbed so this is the end of the line for this photon */
        stuff_phot (&pp, p);
        break;
      }
    }

    if (istat == P_SCAT)
    {                           /* Cause the photon to scatter and reinitilize */


      pp.grid = n = where_in_grid (wmain[pp.grid].ndom, pp.x);

      if (n < 0)
      {
        Error ("trans_phot: Trying to scatter a photon which is not in the wind\n");
        Error ("trans_phot: %d grid %3d x %8.2e %8.2e %8.2e\n", pp.np, pp.grid, pp.x[0], pp.x[1], pp.x[2]);
        Error ("trans_phot: This photon is effectively lost!\n");
        istat = pp.istat = p->istat = P_ERROR;
        stuff_phot (&pp, p);
        break;
      }

      /* 1506 JM -- If the next errors reoccur, see Issue #154 for discussion */

      if (wmain[n].nplasma == NPLASMA)
      {
        Error ("trans_phot: Trying to scatter a photon which is not in a cell in the plasma structure\n");
        Error ("trans_phot: %d grid %3d x %8.2e %8.2e %8.2e\n", pp.np, pp.grid, pp.x[0], pp.x[1], pp.x[2]);
        Error ("trans_phot: This photon is effectively lost!\n");
        istat = pp.istat = p->istat = P_ERROR;
        stuff_phot (&pp, p);
        break;
      }


      if (wmain[n].vol <= 0)
      {
        Error ("trans_phot: Trying to scatter a photon in a cell with no wind volume\n");
        Error ("trans_phot: %d grid %3d x %8.2e %8.2e %8.2e\n", pp.np, pp.grid, pp.x[0], pp.x[1], pp.x[2]);
        Error ("trans_phot: istat %d\n", pp.istat);
        Error ("trans_phot: This photon is effectively lost!\n");
        istat = pp.istat = p->istat = P_ERROR;
        stuff_phot (&pp, p);
        break;

      }

      /* Add path lengths for reverberation mapping */
      if ((geo.reverb == REV_WIND || geo.reverb == REV_MATOM) && geo.ioniz_or_extract && geo.wcycle == geo.wcycles - 1)
      {
        wind_paths_add_phot (&wmain[n], &pp);
      }


      /* SS July 04 - next lines modified so that the "thermal trapping" model of anisotropic scattering is included in the
         macro atom method. What happens now is all in scatter - within that routine the "thermal trapping" model is used to
         decide what the direction of emission is before returning here.
       */


      nnscat = 0;
      nnscat++;
      ptr_nres = &nres;

      if (sane_check (pp.w))
      {
        Error ("trans_phot:sane_check photon %d has weight %e before scatter\n", p->np, pp.w);
      }
      if ((nerr = scatter (&pp, ptr_nres, &nnscat)) != 0)
      {
        Error ("trans_phot: Bad return from scatter %d at point 2", nerr);
      }
      pp.nscat++;


      if (sane_check (pp.w))
      {
        Error ("trans_phot:sane_check photon %d has weight %e after scatter\n", p->np, pp.w);
      }

      /* SS June 04: During the spectrum calculation cycles, photons are thrown away when they interact with macro atoms or
         become k-packets. This is done by setting their weight to zero (effectively they have been absorbed into either
         excitation energy or thermal energy). Since they now have no weight there is no need to follow them further. */
      /* 54b-ksl ??? Stuart do you really mean the comment above; it's not obvious to me since if true why does one need to
         calculate the progression of photons through the wind at all??? Also how is this enforced; where is pp.w set to a
         low value. */
      /* JM 1504 -- This is correct. It's one of the odd things about combining the macro-atom approach with our way of doing
         'spectral cycles'. If photons activate macro-atoms they are destroyed, but we counter this by generating photons
         from deactivating macro-atoms with the already calculated emissivities. */

      if (geo.matom_radiation == 1 && geo.rt_mode == RT_MODE_MACRO && pp.w < weight_min)
        /* Flag for the spectrum calculations in a macro atom calculation SS */
      {
        istat = pp.istat = P_ABSORB;
        pp.tau = VERY_BIG;
        stuff_phot (&pp, p);
        break;
      }

      // Calculate the line heating and if the photon was absorbed break finish up
      // XXXX ??? Need to modify line_heat for multiple scattering but not yet
      // Condition that nres < nlines added (SS)

      if (nres > -1 && nres < nlines)
      {
        pp.nrscat++;

        /* This next statement writes out the position of every resonant scattering event to a file */
        if (modes.track_resonant_scatters)
          track_scatters (&pp, wmain[n].nplasma, "Resonant");


        plasmamain[wmain[n].nplasma].scatters[line[nres].nion] += 1;

        if (geo.rt_mode == RT_MODE_2LEVEL)      // only do next line for non-macro atom case
        {
          line_heat (&plasmamain[wmain[n].nplasma], &pp, nres);
        }

        if (pp.w < weight_min)
        {
          istat = pp.istat = P_ABSORB;  /* This photon was absorbed by continuum opacity within the wind */
          pp.tau = VERY_BIG;
          stuff_phot (&pp, p);
          break;
        }
      }


      /* The next if statement causes photons to be extracted during the creation of the detailed spectrum portion of the
         program */

      /* N.B. To use the anisotropic scattering option, extract needs to follow scatter.  This is because the reweighting
         which occurs in extract needs the pdf for scattering to have been initialized. 02may ksl.  This seems to be OK at
         present. */

      if (iextract)
      {
        stuff_phot (&pp, &pextract);


        /* JM 1407 -- This next loop is required because in anisotropic scattering mode 2 we have normalised our rejection
           method. This means that we have to adjust nnscat by this factor, since nnscat will be lower by a factor of
           1/p_norm */
        if (geo.scatter_mode == SCATTER_MODE_THERMAL && pextract.nres <= NLINES && pextract.nres > -1)
        {
          /* we normalised our rejection method by the escape probability along the vector of maximum velocity gradient.
             First find the sobolev optical depth along that vector. The -1 enforces calculation of the ion density */
          tau_norm = sobolev (&wmain[pextract.grid], pextract.x, -1.0, lin_ptr[pextract.nres], wmain[pextract.grid].dvds_max);

          /* then turn into a probability */
          p_norm = p_escape_from_tau (tau_norm);

        }
        else
        {
          p_norm = 1.0;

          /* nnscat is the quantity associated with this photon being extracted */
          if (nnscat != 1)
            Error
              ("nnscat is %i for photon %i in scatter mode %i! nres %i NLINES %i\n",
               nnscat, p->np, geo.scatter_mode, pextract.nres, NLINES);
        }

        /* We then increase weight to account for number of scatters. This is done because in extract we multiply by the
           escape probability along a given direction, but we also need to divide the weight by the mean escape
           probability, which is equal to 1/nnscat */

        pextract.w *= nnscat / p_norm;

        if (sane_check (pextract.w))
        {
          Error ("trans_phot: sane_check photon %d has weight %e before extract\n", p->np, pextract.w);
        }
        extract (w, &pextract, PTYPE_WIND);     // Treat as wind photon for purpose of extraction
      }


      /* OK we are ready to continue the processing of a photon which has scattered.
       * The next steps reinitialize parameters
       so that the photon can continue throug the wind */

      tau_scat = -log (1. - random_number (0.0, 1.0));
      istat = pp.istat = P_INWIND;
      tau = 0;

      stuff_v (pp.x, x_dfudge_check);   // this is a vector we use to see if dfudge moved the photon outside the wind cone
      reposition (&pp);

      /* JM 1506 -- call walls again to account for instance where DFUDGE
         can take photon outside of the wind and into the disk or star
         after scattering. Note that walls updates the istat in pp as well.
         This may not be necessary but I think to account for every eventuality
         it should be done */
      istat = walls (&pp, p, normal);

      /* This *does not* update istat if the photon scatters outside of the wind-
         I guess P_IN_WIND is really in wind or empty space but not escaped.
         translate_in_space will take care of this next time round. All a bit
         convoluted but should work. */

      /* JM 1506 -- we don't throw errors here now, but we do keep a track
         of how many 4 photons were lost due to DFUDGE pushing them
         outside of the wind after scatter */

      if (where_in_wind (pp.x, &ndom) != W_ALL_INWIND && where_in_wind (x_dfudge_check, &ndom) == W_ALL_INWIND)
      {
        n_lost_to_dfudge++;     // increment the counter (checked at end of trans_phot)
      }

      stuff_phot (&pp, p);
      icell = 0;
    }

    /* This completes the portion of the code that handles the scattering of a photon.
     * What follows is a simple check to see if
     * this particular photon has gotten stuck in the wind */

    if (pp.nscat == MAXSCAT)
    {
      istat = pp.istat = P_TOO_MANY_SCATTERS;
      stuff_phot (&pp, p);
      break;
    }

    if (pp.istat == P_ERROR_MATOM || pp.istat == P_LOFREQ_FF || pp.istat == P_ADIABATIC)
    {
      istat = p->istat = pp.istat;
      stuff_phot (&pp, p);
      break;
    }

    /* This appears partly to be an insurance policy. It is not obvious that for example nscat
     * and nrscat need to be updated */

    p->istat = istat;
    p->nscat = pp.nscat;
    p->nrscat = pp.nrscat;
    p->w = pp.w;                // Assure that final weight of photon is returned.

  }
  /* This is the end of the loop over a photon */

  /* The next section is for diagnostic purposes.  There are two possibilities.  If you wish to know where
   * the photon was last while in the wind, you want to track p; if you wish to know where it hits the
   * outer boundary of the calculation you would want pp.  So one should keep both lines below, and comment
   * out the one you do not want. */

//OLD  if (modes.save_photons)
//OLD    {
//OLD      // save_photons (p, "Final");  // Where the last position of the photon in the wind
//OLD      save_photons (&pp, "Final"); //The position of the photon where it exits the calculation
//OLD    }
  return (0);
}
