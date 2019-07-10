/***********************************************************/
/** @file  photo_gen_matom.c
 * @author ksl, ss, jm
 * @date   January, 2018
 *
 * @brief functions for calculating emissivities and generating photons from macro-atoms and k-packets.
 *   during the spectral cycles. The actual functions which do the jumps inside an activated 
 *  macro-atom are in matom.c. This is partly done to prevent overly long files (JM1504)
 *
 ***********************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "atomic.h"
#include "python.h"
#include <gsl/gsl_block.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include "my_linalg.h"

/**********************************************************/
/** 
 * @brief      returns the specific luminosity in the band needed for the computation of the
 *        spectrum. It gets the total energy radiated by the process k-packet -> r-packet in the
 *        required wavelength range.
 *
 * @return double lum  The energy radiated by the process k-packet -> r-packet in the wind in the 
 *            wavelength range required for the specrum calculation.
 *
 * ### Notes ###
 * Literally just loops through the kpkt_emiss entries in plasmamain.
 *
 **********************************************************/

double
get_kpkt_f ()
{
  int n;
  double lum;

  lum = 0.0;

  for (n = 0; n < NPLASMA; n++)
  {
    lum += plasmamain[n].kpkt_emiss;
  }


  return (lum);
}

/**********************************************************/
/** 
 * @brief returns the specific luminosity in kpkts from nonthermal ("shock")
 *        heating. This is used to generate kpkts in the ionization cycles.
 *        This also populates the cell-by-cell kpkt luminosities in the 
 *        variable plasmamain[n].kpkt_emiss.
 *
 * @return double lum  
 *         The energy created by non-radiative heating throughout the 
 *         computational domain. 
 *
 **********************************************************/

double
get_kpkt_heating_f ()
{
  int n, nwind;
  double lum, shock_kpkt_luminosity;
  WindPtr one;

  lum = 0.0;

  for (n = 0; n < NPLASMA; n++)
  {
    nwind = plasmamain[n].nwind;
    one = &wmain[nwind];

    /* what we do depends on how the "net heating mode" is defined */
    if (KPKT_NET_HEAT_MODE)
      shock_kpkt_luminosity = (shock_heating (one) - plasmamain[n].cool_adiabatic);
    else
      shock_kpkt_luminosity = shock_heating (one);

    if (shock_kpkt_luminosity > 0)
    {
      if (geo.ioniz_or_extract)
        plasmamain[n].kpkt_emiss = shock_kpkt_luminosity;
      else
        plasmamain[n].kpkt_abs += shock_kpkt_luminosity;

      lum += shock_kpkt_luminosity;
    }
    else
      plasmamain[n].kpkt_emiss = 0.0;
  }

  return (lum);
}



/**********************************************************/
/** 
 * @brief      returns the specific band-limited luminosity in macro-atoms
 *
 * @param [in] int  mode   vvariable which controls whether or not we need to compute the
 *            emissivities (CALCULATE_MATOM_EMISSIVITIES) or use stored ones
 *            because we are restarting a spectral cycle (USE_STORED_MATOM_EMISSIVITIES)
 *            see #define statements in python.h and code in xdefine_phot().
 * @return double lum  The energy radiated by the deactivation of macro atoms in the wind in the 
 *            wavelength range required for the specrum calculation.
 *
 * @details
 * this routine calculates the luminosity in the band needed for the computation of the
 * spectrum. It gets the total energy radiated by the deactivation of macro atoms in the
 * required wavelength range. This can be a slow process, as there is no priori way 
 * (at least with our method) to work out where a photon is going to come out when a 
 * macro-atom is generated
 *
 * ### Notes ###
 * Consult Matthews thesis. 
 *
 **********************************************************/

double
get_matom_f (mode)
     int mode;
{
  int n, m;
  int mm, ss;
  double lum;
  int level_emit[NLEVELS_MACRO], kpkt_emit;
  int n_tries, n_tries_local;
  struct photon ppp;
  double contribution, norm;
  int nres, which_out;
  int my_nmin, my_nmax;         //These variables are used even if not in parallel mode


  if (mode == USE_STORED_MATOM_EMISSIVITIES)
  {
    Log ("geo.pcycle (%i) > 0, so using saved level emissivities for macro atoms\n", geo.pcycle);
  }

  else                          // we need to compute the emissivities
  {
#ifdef MPI_ON
    int num_mpi_cells, num_mpi_extra, position, ndo, n_mpi, num_comm, n_mpi2;
    int size_of_commbuffer;
    char *commbuffer;

    /* the commbuffer needs to communicate 2 variables and the number of macor levels, 
       plus the variable for how many cells each thread is doing */
    size_of_commbuffer = 8 * (3 + nlevels_macro) * (floor (NPLASMA / np_mpi_global) + 1);

    commbuffer = (char *) malloc (size_of_commbuffer * sizeof (char));
#endif

    /* add the non-radiative k-packet heating to the kpkt_abs quantity */
    get_kpkt_heating_f ();

    which_out = 0;
    n_tries = 5000000;
    geo.matom_radiation = 0;
    n_tries_local = 0;



    norm = 0;
    for (n = 0; n < NPLASMA; n++)
    {
      for (m = 0; m < nlevels_macro; m++)
      {
        norm += macromain[n].matom_abs[m];
        macromain[n].matom_emiss[m] = 0.0;
        if (sane_check (macromain[n].matom_abs[m]))
          Error ("matom_abs is %8.4e in matom %i level %i\n", macromain[n].matom_abs[m], n, m);
      }
      norm += plasmamain[n].kpkt_abs;
      plasmamain[n].kpkt_emiss = 0.0;
      if (sane_check (plasmamain[n].kpkt_abs))
        Error ("kpkt_abs is %8.4e in matom %i\n", plasmamain[n].kpkt_abs, n);
    }



    Log ("Calculating macro-atom and k-packet emissivities- this might take a while...\n");
    Log ("Number of macro-atom levels: %d\n", nlevels_macro);

    /* For MPI parallelisation, the following loop will be distributed over multiple tasks. 
       Note that the mynmim and mynmax variables are still used even without MPI on */
    my_nmin = 0;
    my_nmax = NPLASMA;

#ifdef MPI_ON

    num_mpi_cells = floor (NPLASMA / np_mpi_global);    // divide the cells between the threads
    num_mpi_extra = NPLASMA - (np_mpi_global * num_mpi_cells);  // the remainder from the above division

    /* this next loop splits the cells up between the threads. All threads with 
       rank_global<num_mpi_extra deal with one extra cell to account for the remainder */
    if (rank_global < num_mpi_extra)
    {
      my_nmin = rank_global * (num_mpi_cells + 1);
      my_nmax = (rank_global + 1) * (num_mpi_cells + 1);
    }
    else
    {
      my_nmin = num_mpi_extra * (num_mpi_cells + 1) + (rank_global - num_mpi_extra) * (num_mpi_cells);
      my_nmax = num_mpi_extra * (num_mpi_cells + 1) + (rank_global - num_mpi_extra + 1) * (num_mpi_cells);
    }
    ndo = my_nmax - my_nmin;

    Log_parallel ("Thread %d is calculating macro atom emissivities for macro atoms %d to %d\n", rank_global, my_nmin, my_nmax);

#endif




    for (n = my_nmin; n < my_nmax; n++)
    {

      /* JM 1309 -- these lines are just log statements which track progress, as this section
         can take a long time */
#ifdef MPI_ON
      if (n % 50 == 0)
        Log
          ("Thread %d is calculating  macro atom emissivity for macro atom %7d of %7d or %6.3f per cent\n",
           rank_global, n, my_nmax, n * 100. / my_nmax);
#else
      if (n % 50 == 0)
        Log ("Calculating macro atom emissivity for macro atom %7d of %7d or %6.3f per cent\n", n, my_nmax, n * 100. / my_nmax);
#endif

      for (m = 0; m < nlevels_macro + 1; m++)
      {
        if ((m == nlevels_macro && plasmamain[n].kpkt_abs > 0) || (m < nlevels_macro && macromain[n].matom_abs[m] > 0))
        {
          if (m < nlevels_macro)
          {
            if (macromain[n].matom_abs[m] > 0)
            {
              n_tries_local = (n_tries * macromain[n].matom_abs[m] / norm) + 10;
            }
            else
            {
              n_tries_local = 0;
            }
          }
          else if (m == nlevels_macro)
          {
            if (plasmamain[n].kpkt_abs > 0)
            {
              n_tries_local = (n_tries * plasmamain[n].kpkt_abs / norm) + 10;
            }
            else
            {
              n_tries_local = 0;
            }
          }
          /* We know "matom_abs" is the amount of energy absobed by
             each macro atom level in each cell. We now want to determine what fraction of
             that energy re-appears in the frequency range we want and from which macro atom
             level it is re-emitted (or if it appears via a k-packet). */

          for (mm = 0; mm < NLEVELS_MACRO; mm++)
          {
            level_emit[mm] = 0;
          }
          kpkt_emit = 0;
          if (n_tries_local > 0)
          {
            for (ss = 0; ss < n_tries_local; ss++)
            {
              if (m < nlevels_macro)
              {
                /* Dealing with excitation of a macro atom level. */
                /* First find a suitable transition in which to excite the macro atom (to get it
                   in the correct starting level. */

                /* As Knox pointed out in an earlier comment
                   here, the next loop isn't the best way to
                   do this - however I'm not sure we can
                   just replace nlines->nlines_macro since
                   we don't know for sure that the
                   macro_lines are the first lines in the
                   array. If this is really a bottle neck we
                   should set up a quicker way to identify
                   the level - the purpose is just to set
                   the macro atom state to the appropriate
                   upper level so all we need is some way to
                   tell it which state it should be. */

                nres = 0;
                while ((nres < nlines))
                {
                  if (lin_ptr[nres]->nconfigu != m)
                  {
                    nres += 1;
                  }
                  else
                  {
                    break;
                  }
                }


                if (nres > (nlines - 1))
                {
                  nres = NLINES + 1;
                  while (nres < (NLINES + 1 + nphot_total))
                  {
                    if (phot_top[nres - NLINES - 1].uplev != m)
                    {
                      nres += 1;
                    }
                    else
                    {
                      break;
                    }
                  }
                }

                if (nres > NLINES + nphot_total)
                {
                  Error ("Problem in get_matom_f (1). Abort. nres is %d, NLINES %d, nphot_total %d m %d %8.4e\n",
                         nres, NLINES, nphot_total, m, macromain[n].matom_abs[m]);
                  Exit (0);
                }

                ppp.nres = nres;
                ppp.grid = plasmamain[n].nwind;
                ppp.w = 0;
                /* this needs to be initialised because we set to istat to P_ADIABATIC
                   for adiabatic destruction */
                ppp.istat = P_INWIND;

                macro_gov (&ppp, &nres, 1, &which_out);


                /* Now a macro atom has been excited and followed until an r-packet is made. Now, if that 
                   r-packet is in the correct frequency range we record it. If not, we throw it away. */
              }
              else if (m == nlevels_macro)
              {
                /* kpkt case. */

                nres = -2;      //will do - just need
                //something that will tigger kpkt

                ppp.nres = nres;
                ppp.grid = plasmamain[n].nwind;
                ppp.w = 0;
                /* this needs to be initialised because we set to istat to P_ADIABATIC
                   for adiabatic destruction */
                ppp.istat = P_INWIND;

                macro_gov (&ppp, &nres, 2, &which_out);

                /* We have an r-packet back again. */
              }


              if (ppp.freq > geo.sfmin && ppp.freq < geo.sfmax && ppp.istat != P_ADIABATIC)
              {
                if (which_out == 1)
                {
                  if (nres < 0)
                  {
                    Error ("Negative out from matom?? Abort.\n");
                    Exit (0);
                  }

                  /* It was a macro atom de-activation. */
                  if (nres < nlines)
                  {             /* It was a line. */
                    level_emit[lin_ptr[nres]->nconfigu] += 1;
                  }
                  else
                  {
                    level_emit[phot_top[nres - NLINES - 1].uplev] += 1;
                  }
                }
                else if (which_out == 2)
                {
                  /* It was a k-packet de-activation. */
                  kpkt_emit += 1;
                }
                else
                {
                  Error ("Packet didn't emerge from matom or kpkt??? Abort. \n");
                  Exit (0);
                }
              }
            }


            /* Now we've done all the runs for this level de-activating so we can add the contributions
               to the level emissivities, the k-packet emissivity and the total luminosity. */

            for (mm = 0; mm < nlevels_macro; mm++)
            {
              contribution = 0;
              if (m < nlevels_macro)
              {
                macromain[n].matom_emiss[mm] += contribution = level_emit[mm] * macromain[n].matom_abs[m] / n_tries_local;
              }
              else if (m == nlevels_macro)
              {
                macromain[n].matom_emiss[mm] += contribution = level_emit[mm] * plasmamain[n].kpkt_abs / n_tries_local;
              }
            }

            if (m < nlevels_macro)
            {
              plasmamain[n].kpkt_emiss += kpkt_emit * macromain[n].matom_abs[m] / n_tries_local;
            }
            else if (m == nlevels_macro)
            {
              plasmamain[n].kpkt_emiss += kpkt_emit * plasmamain[n].kpkt_abs / n_tries_local;
            }
          }
        }
      }
    }


    /*This is the end of the update loop that is parallelised. We now need to exchange data between the tasks.
       This is done much the same way as in wind_update */
#ifdef MPI_ON

    /* JM Add an MPI Barrier here */
    MPI_Barrier (MPI_COMM_WORLD);

    for (n_mpi = 0; n_mpi < np_mpi_global; n_mpi++)
    {
      /* here we loop over the number of threads. If the thread is this thread then we pack the macromain information
         into memory and then broadcast it to ther other threads */
      position = 0;

      if (rank_global == n_mpi)
      {

        Log ("MPI task %d is working on matoms %d to max %d (total size %d).\n", rank_global, my_nmin, my_nmax, NPLASMA);

        MPI_Pack (&ndo, 1, MPI_INT, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        for (n = my_nmin; n < my_nmax; n++)
        {

          /* pack the number of the cell, and the kpkt and macro atom emissivites for that cell // */
          MPI_Pack (&n, 1, MPI_INT, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
          MPI_Pack (&plasmamain[n].kpkt_emiss, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
          MPI_Pack (macromain[n].matom_emiss, nlevels_macro, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);

        }

        Log_parallel ("MPI task %d broadcasting matom emissivity information.\n", rank_global);
      }


      /* Set MPI_Barriers and broadcast information to other threads */
      MPI_Barrier (MPI_COMM_WORLD);
      MPI_Bcast (commbuffer, size_of_commbuffer, MPI_PACKED, n_mpi, MPI_COMM_WORLD);
      MPI_Barrier (MPI_COMM_WORLD);
      Log_parallel ("MPI task %d survived broadcasting matom emissivity information.\n", rank_global);



      position = 0;

      /* If not this thread then we unpack the macromain information from the other threads */

      if (rank_global != n_mpi)
      {
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &num_comm, 1, MPI_INT, MPI_COMM_WORLD);
        for (n_mpi2 = 0; n_mpi2 < num_comm; n_mpi2++)
        {

          /* unpack the number of the cell, and the kpkt and macro atom emissivites for that cell */
          MPI_Unpack (commbuffer, size_of_commbuffer, &position, &n, 1, MPI_INT, MPI_COMM_WORLD);
          MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].kpkt_emiss, 1, MPI_DOUBLE, MPI_COMM_WORLD);
          MPI_Unpack (commbuffer, size_of_commbuffer, &position, macromain[n].matom_emiss, nlevels_macro, MPI_DOUBLE, MPI_COMM_WORLD);
        }
      }
    }                           // end of parallelised section

    /* add an MPI Barrier after unpacking stage */
    MPI_Barrier (MPI_COMM_WORLD);
#endif

  }                             // end of if loop which controls whether to compute the emissivities or not 

  /* this next loop just calculates lum to be the correct summed value in parallel mode */
  /* if mode == USE_STORED_MATOM_EMISSIVITIES this is all this routine does */
  lum = 0.0;                    // need to rezero, fixes segfault bug #59 

  for (n = 0; n < NPLASMA; n++)
  {

    for (mm = 0; mm < nlevels_macro; mm++)
    {
      lum += macromain[n].matom_emiss[mm];
    }
  }

  geo.matom_radiation = 1;

  return (lum);
}

/* All done. */



/**********************************************************/
/** 
 * @brief produces photon packets to account for creating of r-packets by k-packets. 

 *
 * @param [in, out] PhotPtr  p   the ptr to the entire structure for the photons
 * @param [in] double  weight   the photon weight
 * @param [in] int  photstart   The first element of the photon stucure to be populated
 * @param [in] int  nphot   the number of photons to be generated
 * @return int nphot  When it finishes it should have generated nphot photons from k-packet elliminations.
 *
 * @details produces photon packets to account for creating of r-packets by k-packets in the spectrum calculation. 
 * It should only be used once the total energy emitted in this way in the wavelength range in question is well known
 * (calculated in the ionization cycles). This routine is closely related to photo_gen_wind from which much of the code 
 * has been copied.
 *
 **********************************************************/

int
photo_gen_kpkt (p, weight, photstart, nphot)
     PhotPtr p;
     double weight;
     int photstart, nphot;
{
  int photstop;
  int icell;
  double xlum, xlumsum;
  struct photon pp;
  int nres, esc_ptr, which_out;
  int n;
  double v[3];
  double dot ();
  double test;
  int nnscat;
  double dvwind_ds (), sobolev ();
  int nplasma, ndom;
  int kpkt_mode;
  double fmin, fmax;

  photstop = photstart + nphot;
  Log ("photo_gen_kpkt  creates nphot %5d photons from %5d to %5d, weight %8.4e \n", nphot, photstart, photstop, weight);

  if (geo.ioniz_or_extract)
  {
    /* we are in the ionization cycles, so use all frequencies. kpkt_mode should allow all processes */
    fmin = xband.f1[0];
    fmax = xband.f2[xband.nbands - 1];
    kpkt_mode = KPKT_MODE_ALL;
  }
  else
  {
    /* we are in the spectral cycles, so use all the required frequency range */
    fmin = geo.sfmin;
    fmax = geo.sfmax;
    /* we only want k->r processes */
    kpkt_mode = KPKT_MODE_CONTINUUM;
  }

  for (n = photstart; n < photstop; n++)
  {
    /* locate the wind_cell in which the photon bundle originates. */

    xlum = random_number (0.0, 1.0) * geo.f_kpkt;

    xlumsum = 0;
    icell = 0;
    while (xlumsum < xlum)
    {
      if (wmain[icell].vol > 0.0)
      {
        nplasma = wmain[icell].nplasma;
        xlumsum += plasmamain[nplasma].kpkt_emiss;
      }
      icell++;
    }
    icell--;                    /* This is the cell in which the photon must be generated */

    /* Now generate a single photon in this cell */
    p[n].w = weight;

    pp.freq = 0.0;
    pp.grid = icell;

    /* This following block is a bad way of doing it - kpkt could be modified to 
       do what we want in a more elegant way. */

    test = pp.freq;

    while (test > fmax || test < fmin)
    {
      pp.w = p[n].w;
      kpkt (&pp, &nres, &esc_ptr, kpkt_mode);

      if (esc_ptr == 0 && kpkt_mode == KPKT_MODE_CONTINUUM)
      {
        test = 0.0;
      }
      else
      {
        if (esc_ptr == 0)
        {
          macro_gov (&pp, &nres, 1, &which_out);
        }
        test = pp.freq;
      }
    }

    p[n].freq = pp.freq;
    p[n].nres = nres;
    p[n].w = pp.w;


    /* The photon frequency is now known. */

    /* Determine the position of the photon in the moving frame */

    get_random_location (icell, p[n].x);

    p[n].grid = icell;

    nnscat = 1;
    // Determine the direction of the photon
    // Need to allow for anisotropic emission here
    if (p[n].nres < 0 || p[n].nres > NLINES || geo.scatter_mode == SCATTER_MODE_ISOTROPIC)
    {
      /*  It was either an electron scatter so the  distribution is isotropic, or it
         was a resonant scatter but we want isotropic scattering anyway.  */
      randvec (p[n].lmn, 1.0);  /* The photon is emitted isotropically */
    }
    else if (geo.scatter_mode == SCATTER_MODE_THERMAL)
    {                           //It was a line photon and we want the thermal trapping anisotropic model

      randwind_thermal_trapping (&p[n], &nnscat);
    }

    p[n].nnscat = nnscat;

    /* The next two lines correct the frequency to first order, but do not result in
       forward scattering of the distribution */

    ndom = wmain[icell].ndom;
    vwind_xyz (ndom, &p[n], v);
    p[n].freq /= (1. - dot (v, p[n].lmn) / C);

    p[n].istat = 0;
    p[n].tau = p[n].nscat = p[n].nrscat = 0;
    p[n].origin = PTYPE_WIND;   // Call it a wind photon

    switch (geo.reverb)
    {                           //0715 SWM - Added path generation
    case REV_MATOM:
      line_paths_gen_phot (&wmain[icell], &p[n], nres);
      break;
    case REV_WIND:
      wind_paths_gen_phot (&wmain[icell], &p[n]);
      break;
    case REV_PHOTON:
      simple_paths_gen_phot (&p[n]);
      break;
    case REV_NONE:
    default:
      break;
    }

  }



  return (nphot);               /* Return the number of photons generated */


  /* All done. */
}

/**********************************************************/
/** 
 * @brief      produces photon packets to account for creation of r-packets
 *      by deactivation of macro atoms in the spectrum calculation. It should only be used 
 *      once the total energy emitted in this way in the wavelength range in question is well known
 *      (calculated in the ionization cycles).
 *
 * @param [in, out] PhotPtr  p   the ptr to the structire for the photons
 * @param [in] double  weight   the photon weight
 * @param [in] int  photstart   The position in the photon structure of the first photon to generate
 * @param [in] int  nphot   the number of the photons to be generated 
 * @return int nphot When it finishes it should have generated nphot photons from macro atom deactivations.
 *
 * @details
 * This routine is closely related to photo_gen_kpkt from which much of the code has been copied.
 *
 * ### Notes ###
 * Consult Matthews thesis. 
 *
 **********************************************************/

int
photo_gen_matom (p, weight, photstart, nphot)
     PhotPtr p;
     double weight;
     int photstart, nphot;
{
  int photstop;
  int icell;
  double xlum, xlumsum;
  struct photon pp;
  int nres;
  int n;
  double v[3];
  double dot ();
  int emit_matom ();
  double test;
  int upper;
  int nnscat;
  double dvwind_ds (), sobolev ();
  int nplasma;
  int ndom;



  photstop = photstart + nphot;
  Log ("photo_gen_matom creates nphot %5d photons from %5d to %5d \n", nphot, photstart, photstop);

  for (n = photstart; n < photstop; n++)
  {
    /* locate the wind_cell in which the photon bundle originates. And also decide which of the macro
       atom levels will be sampled (identify that level as "upper"). */
    xlum = random_number (0.0, 1.0) * geo.f_matom;


    xlumsum = 0;
    icell = 0;
    upper = 0;
    while (xlumsum < xlum)
    {
      if (wmain[icell].vol > 0.0)
      {
        nplasma = wmain[icell].nplasma;
        xlumsum += macromain[nplasma].matom_emiss[upper];
        upper++;
        if (upper == nlevels_macro)
        {
          upper = 0;
          icell++;
        }
      }
      else
      {
        icell++;
        upper = 0;
      }
    }

    if (upper == 0)
    {
      /* If upper = 0 at this point it means the loop above finished with an increment of
         icell and setting upper = 0. In such a case the process we want was the LAST process
         in the previous cell. */
      icell--;
      upper = nlevels_macro;
    }

    upper--;                    /* This leaves the macro atom level that deactivaties. */

    /* Now generate a single photon in this cell */
    p[n].w = weight;

    pp.freq = 0.0;
    pp.grid = icell;

    /* This following block is a bad way of doing it but it'll do as
       a quick and dirty test for now. Really kpkt should be modified to 
       do what we want in a more elegant way. */

    test = pp.freq;

    while (test > geo.sfmax || test < geo.sfmin)
    {

      /* Call routine that will select an emission process for the
         deactivating macro atom. If is deactivates outside the frequency
         range of interest then ignore it and try again. SS June 04. */

      emit_matom (wmain, &pp, &nres, upper);

      test = pp.freq;
    }


    p[n].freq = pp.freq;
    p[n].nres = nres;


    /* The photon frequency is now known. */

    /* Determine the position of the photon in the moving frame */


    get_random_location (icell, p[n].x);

    p[n].grid = icell;


    /* Determine the direction of the photon
       Need to allow for anisotropic emission here
     */

    nnscat = 1;
    if (p[n].nres < 0 || p[n].nres > NLINES || geo.scatter_mode == SCATTER_MODE_ISOTROPIC)
    {
      /*  It was either an electron scatter so the  distribution is isotropic, or it
         was a resonant scatter but we want isotropic scattering anyway.  */
      randvec (p[n].lmn, 1.0);  /* The photon is emitted isotropically */
    }
    else if (geo.scatter_mode == SCATTER_MODE_THERMAL)
    {                           //It was a line photon and we want the thermal trapping anisotropic model

      randwind_thermal_trapping (&p[n], &nnscat);
    }
    p[n].nnscat = nnscat;


    /* The next two lines correct the frequency to first order, but do not result in
       forward scattering of the distribution */

    ndom = wmain[icell].ndom;
    vwind_xyz (ndom, &p[n], v);
    p[n].freq /= (1. - dot (v, p[n].lmn) / C);

    p[n].istat = 0;
    p[n].tau = p[n].nscat = p[n].nrscat = 0;
    p[n].origin = PTYPE_WIND;   // Call it a wind photon

    switch (geo.reverb)
    {                           //0715 SWM - Added path generation
    case REV_MATOM:
      line_paths_gen_phot (&wmain[icell], &p[n], nres);
      break;
    case REV_WIND:
      wind_paths_gen_phot (&wmain[icell], &p[n]);
      break;
    case REV_PHOTON:
      simple_paths_gen_phot (&p[n]);
      break;
    case REV_NONE:
    default:
      break;
    }
  }


  return (nphot);               /* Return the number of photons generated */


  /* All done. */
}
