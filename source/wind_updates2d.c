
/***********************************************************/
/** @file  wind_updates2d.c
 * @author ksl
 * @date   May, 2018
 *
 * @brief  This file contains the main routines for updating
 * and then reinitializing the wind after an ionization cycle
 *
 * The routines in this file are generic.  There is no dependence on a particlar wind model or
 * any coordinate system dependences.
 ***********************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>


#include "atomic.h"
#include "python.h"
#define LINELEN 256



int num_updates = 0;

/**********************************************************/
/**
 * @brief      updates the parameters in the wind that are
 * 	affected by radiation, including ion densities.
 *
 * @param [in] WindPtr  The entire wind
 * @return     Always returns 0
 *
 * @details
 * This is the main routine used to update the wind at the end of
 * an ionization cycle (in preparation for a new cycle).  The routine
 * is parallelized to save time
 *
 * ### Notes ###
 * At the time wind_update is called the various quantities that are accumulated
 * during the photon transfer part of the cycle has been accumulated
 * and shared between the threads.  Here certain plasma cells are assigned to
 * each thread so that the ionization structure can be updated, and each thread
 * is responsible for calculaing the updates for a certain set of cells.  At
 * the end of the routine the updates are collected and reshared.
 *
 * The real need for prallelising the routine is the work done in ion_abundances
 *
 * Once this is done, various checks are made to determined what happened as a
 * function of the updates, various variables in geo are updated,  and for
 * a hydro model the results are written to a file
 *
 * This routine is nearly 1000 lines long and might beneifit from breaking
 * it into functionl blocks, e.g by separating out the mpi communicaiton
 * into their own routines.
 *
 *
 **********************************************************/
int
wind_update (w)
WindPtr (w);
{
  int n, i, j;
  double trad, nh;

  /*1108 NSH csum added to sum compton heating 1204 NSH icsum added to sum induced compton heating */
  double wtest, xsum, psum, fsum, lsum, csum, icsum, ausum;
  double cool_sum, lum_sum, rad_sum;    //1706 - the total cooling and luminosity of the wind
  double apsum, aausum, abstot; //Absorbed photon energy from PI and auger
  double c_rec, n_rec, o_rec, fe_rec;   //1701- NSH more outputs to show cooling from a few other elements
  double c_lum, n_lum, o_lum, fe_lum;   //1708- NSH and luminosities as well
  double cool_dr_metals;
  int nn;                       //1701 - loop variable to compute recomb cooling

  double volume;
  double vol;
  char string[LINELEN];
  double t_r_old, t_e_old, dt_r, dt_e;
  double t_r_ave_old, t_r_ave, t_e_ave_old, t_e_ave;
  int iave, nmax_r, nmax_e;
  int nplasma, nshell;
  int nwind;
  int first, last, m;
  double tot, agn_ip;
  double lum_h_line, lum_he_line, lum_c_line, lum_n_line, lum_o_line, lum_fe_line;
  double h_dr, he_dr, c_dr, n_dr, o_dr, fe_dr;
  int my_nmin, my_nmax;         //Note that these variables are still used even without MPI on
  int ndom;
  FILE *fptr, *fopen ();        /*This is the file to communicate with zeus */


#ifdef MPI_ON
  int num_mpi_cells, num_mpi_extra, position, ndo, n_mpi, num_comm, n_mpi2;
  int size_of_commbuffer;
  char *commbuffer;

  /* JM 1409 -- Added for issue #110 to ensure correct reporting in parallel */
  int nmax_r_temp, nmax_e_temp;
  double dt_e_temp, dt_r_temp;


  /* the commbuffer needs to be larger enough to pack all variables in MPI_Pack and MPI_Unpack routines NSH 1407 - the
     NIONS changed to nions for the 12 arrays in plasma that are now dynamically allocated
     NSH 1703 changed NLTE_LEVELS to nlte_levels  and NTOP_PHOT to nphot_tot since they are dynamically allocated now */
  size_of_commbuffer =
//OLD    8 * (13 * nions + nlte_levels + 2 * nphot_total + 12 * NXBANDS + 2 * LPDF + NAUGER + 113) * (floor (NPLASMA / np_mpi_global) + 1);
    8 * (13 * nions + nlte_levels + 2 * nphot_total + 12 * NXBANDS + 113) * (floor (NPLASMA / np_mpi_global) + 1);
  commbuffer = (char *) malloc (size_of_commbuffer * sizeof (char));

  /* JM 1409 -- Initialise parallel only variables */
  nmax_r_temp = nmax_e_temp = -1;
  dt_e_temp = dt_r_temp = 0.0;

#endif
  dt_r = dt_e = 0.0;
  iave = 0;
  nmax_r = nmax_e = -1;
  t_r_ave_old = t_r_ave = t_e_ave_old = t_e_ave = 0.0;


  /* For MPI parallelisation, the following loop will be distributed over mutiple tasks.
     Note that the mynmim and mynmax variables are still used even without MPI on */
  my_nmin = 0;
  my_nmax = NPLASMA;
#ifdef MPI_ON
  num_mpi_cells = floor (NPLASMA / np_mpi_global);
  num_mpi_extra = NPLASMA - (np_mpi_global * num_mpi_cells);

  /* this section distributes the remainder over the threads if the cells
     do not divide evenly by thread */
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
#endif

  /* Before we do anything let's record the average tr and te from the last cycle */
  /* JM 1409 -- Added for issue #110 to ensure correct reporting in parallel */
  for (n = 0; n < NPLASMA; n++)
  {
    t_r_ave_old += plasmamain[n].t_r;
    t_e_ave_old += plasmamain[n].t_e;
  }

  /* we now know how many cells this thread has to process - note this will be
     0-NPLASMA in serial mode */

  for (n = my_nmin; n < my_nmax; n++)
  {


    nwind = plasmamain[n].nwind;
    volume = w[nwind].vol;



    /* Start with a call to the routine which normalises all the macro atom
       monte carlo radiation field estimators. It's best to do this first since
       some of the estimators include temperature terms (stimulated correction
       terms) which were included during the monte carlo simulation so we want
       to be sure that the SAME temperatures are used here. (SS - Mar 2004). */

    if (geo.rt_mode == RT_MODE_MACRO && geo.macro_simple == 0)  //test for macro atoms
    {
      mc_estimator_normalise (nwind);
      macromain[n].kpkt_rates_known = -1;
    }

    /* Store some information so one can determine how much the temps are changing */
    t_r_old = plasmamain[n].t_r;
    t_e_old = plasmamain[n].t_e;
    iave++;

    if (plasmamain[n].ntot < 100)
    {
      Log
        ("!!wind_update: Cell %4d Dom %d  Vol. %8.2e r %8.2e theta %8.2e has only %4d photons\n",
         n, w[nwind].ndom, volume, w[nwind].rcen, w[nwind].thetacen, plasmamain[n].ntot);
    }

    if (plasmamain[n].ntot > 0)
    {
      wtest = plasmamain[n].ave_freq;
      plasmamain[n].ave_freq /= plasmamain[n].j;        /* Normalization to frequency moment */
      if (sane_check (plasmamain[n].ave_freq))
      {
        Error ("wind_update:sane_check %d ave_freq %e j %e ntot %d\n", n, wtest, plasmamain[n].j, plasmamain[n].ntot);
      }

      plasmamain[n].j /= (4. * PI * volume);
      plasmamain[n].j_direct /= (4. * PI * volume);
      plasmamain[n].j_scatt /= (4. * PI * volume);

      trad = plasmamain[n].t_r = H * plasmamain[n].ave_freq / (BOLTZMANN * 3.832);
      plasmamain[n].w = PI * plasmamain[n].j / (STEFAN_BOLTZMANN * trad * trad * trad * trad);


      if (plasmamain[n].w > 1e10)
      {
        Error ("wind_update: Huge w %8.2e in cell %d trad %10.2e j %8.2e\n", plasmamain[n].w, n, trad, plasmamain[n].j);
      }
      if (sane_check (trad) || sane_check (plasmamain[n].w))
      {
        Error ("wind_update:sane_check %d trad %8.2e w %8.2g\n", n, trad, plasmamain[n].w);
        Error ("wind_update: ave_freq %8.2e j %8.2e\n", plasmamain[n].ave_freq, plasmamain[n].j);
        exit (0);
      }
    }
    else
    {                           /* It is not clear what to do with no photons in a cell */

      plasmamain[n].j = plasmamain[n].j_direct = plasmamain[n].j_scatt = 0;
      trad = plasmamain[n].t_r;
      plasmamain[n].t_e *= 0.7;
      if (plasmamain[n].t_e < MIN_TEMP)
        plasmamain[n].t_e = MIN_TEMP;
      plasmamain[n].w = 0;
    }


    /* Calculate the frequency banded j and ave_freq variables */

    for (i = 0; i < geo.nxfreq; i++)
    {                           /*loop over number of bands */
      if (plasmamain[n].nxtot[i] > 0)
      {                         /*Check we actually have some photons in the cell in this band */

        plasmamain[n].xave_freq[i] /= plasmamain[n].xj[i];      /*Normalise the average frequency */
        plasmamain[n].xsd_freq[i] /= plasmamain[n].xj[i];       /*Normalise the mean square frequency */
        plasmamain[n].xsd_freq[i] = sqrt (plasmamain[n].xsd_freq[i] - (plasmamain[n].xave_freq[i] * plasmamain[n].xave_freq[i]));       /*Compute standard deviation */
        plasmamain[n].xj[i] /= (4 * PI * volume);       /*Convert to radiation density */

      }
      else
      {
        plasmamain[n].xj[i] = 0;        /*If no photons, set both radiation estimators to zero */
        plasmamain[n].xave_freq[i] = 0;
        plasmamain[n].xsd_freq[i] = 0;  /*NSH 120815 and also the SD ???? */
      }
    }

/* 1108 NSH End of loop */



    nh = plasmamain[n].rho * rho2nh;

/* 1110 NSH Normalise IP, which at this point should be
 * the number of photons in a cell by dividing by volume
 * and number density of hydrogen in the cell
 * */

    plasmamain[n].ip /= (C * volume * nh);
    plasmamain[n].ip_direct /= (C * volume * nh);
    plasmamain[n].ip_scatt /= (C * volume * nh);

/* 1510 NSH Normalise xi, which at this point should be the luminosity of ionizing photons in a cell (just the sum of photon weights) */

    plasmamain[n].xi *= 4. * PI;
    plasmamain[n].xi /= (volume * nh);

    /* If geo.adiabatic is true, then alculate the adiabatic cooling using the current, i.e
     * previous value of t_e.  Note that this may not be  best way to determien the cooling.
     * Changes made here should also be reflected in wind2d.c.  At present, adiabatic colling
     * is not included in updates to the temperature, even if the adiabatic cooling is calculated
     * here. 04nov -- ksl
     * 05apr -- ksl -- The index being used was incorrect.  This has been fixed now
     * 11sep -- nsh -- The index for the wind (&w) for adiabatic cooling was incorrect -
     * was being called with the plasma cell rather than the approriate wind cell - fixed
     * old: adiabatic_cooling (&w[n], plasmamain[n].t_e);
     */

    if (geo.adiabatic)
      plasmamain[n].cool_adiabatic = adiabatic_cooling (&w[nwind], plasmamain[n].t_e);
    else
      plasmamain[n].cool_adiabatic = 0.0;

    if (geo.nonthermal)
      plasmamain[n].heat_shock = shock_heating (&w[nwind]);
    else
      plasmamain[n].heat_shock = 0;


    /* Calculate the densities in various ways depending on the ioniz_mode */

    ion_abundances (&plasmamain[n], geo.ioniz_mode);



    /* Perform checks to see how much temperatures have changed in this iteration */

    if ((fabs (t_r_old - plasmamain[n].t_r)) > fabs (dt_r))
    {
      dt_r = plasmamain[n].t_r - t_r_old;
      nmax_r = n;
    }
    if ((fabs (t_e_old - plasmamain[n].t_e)) > fabs (dt_e))
    {
      dt_e = plasmamain[n].t_e - t_e_old;
      nmax_e = n;
    }
    t_r_ave += plasmamain[n].t_r;
    t_e_ave += plasmamain[n].t_e;
  }



  /*This is the end of the update loop that is parallised. We now need to exchange data between the tasks. */
#ifdef MPI_ON
  for (n_mpi = 0; n_mpi < np_mpi_global; n_mpi++)
  {
    position = 0;

    if (rank_global == n_mpi)
    {
      Log ("MPI task %d is working on cells %d to max %d (total size %d).\n", rank_global, my_nmin, my_nmax, NPLASMA);
      MPI_Pack (&ndo, 1, MPI_INT, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
      for (n = my_nmin; n < my_nmax; n++)
      {
        MPI_Pack (&n, 1, MPI_INT, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].nwind, 1, MPI_INT, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].nplasma, 1, MPI_INT, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].ne, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].rho, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].vol, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (plasmamain[n].density, nions, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (plasmamain[n].partition, nions, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (plasmamain[n].levden, nlte_levels, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (plasmamain[n].PWdenom, nions, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (plasmamain[n].PWdtemp, nions, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (plasmamain[n].PWnumer, nions, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (plasmamain[n].PWntemp, nions, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].kappa_ff_factor, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].nscat_es, 1, MPI_INT, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (plasmamain[n].recomb_simple, nphot_total, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].kpkt_emiss, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].kpkt_abs, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (plasmamain[n].kbf_use, nphot_total, MPI_INT, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].kbf_nuse, 1, MPI_INT, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].t_r, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].t_r_old, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].t_e, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].t_e_old, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].dt_e, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].dt_e_old, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].heat_tot, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].abs_tot, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].heat_tot_old, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].heat_lines, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].heat_ff, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].heat_comp, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].heat_ind_comp, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].heat_lines_macro, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].heat_photo_macro, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].heat_photo, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].heat_auger, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].abs_photo, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].abs_auger, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].heat_z, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].w, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].ntot, 1, MPI_INT, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].ntot_star, 1, MPI_INT, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].ntot_bl, 1, MPI_INT, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].ntot_disk, 1, MPI_INT, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].ntot_wind, 1, MPI_INT, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].ntot_agn, 1, MPI_INT, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].mean_ds, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].n_ds, 1, MPI_INT, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].nrad, 1, MPI_INT, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].nioniz, 1, MPI_INT, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (plasmamain[n].ioniz, nions, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (plasmamain[n].recomb, nions, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (plasmamain[n].scatters, nions, MPI_INT, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (plasmamain[n].xscatters, nions, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (plasmamain[n].heat_ion, nions, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (plasmamain[n].cool_rr_ion, nions, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (plasmamain[n].lum_rr_ion, nions, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].j, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].j_direct, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].j_scatt, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].ave_freq, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].cool_tot, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (plasmamain[n].xj, NXBANDS, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (plasmamain[n].xave_freq, NXBANDS, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (plasmamain[n].xsd_freq, NXBANDS, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (plasmamain[n].nxtot, NXBANDS, MPI_INT, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].max_freq, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].lum_lines, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].lum_ff, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].cool_adiabatic, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].comp_nujnu, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].cool_comp, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].cool_dr, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].cool_di, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].cool_rr, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].lum_rr, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].cool_rr_metals, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].lum_rr_metals, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].lum_tot, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].lum_tot_old, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].cool_tot_ioniz, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].lum_lines_ioniz, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].lum_ff_ioniz, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].cool_adiabatic_ioniz, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].cool_comp_ioniz, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].cool_dr_ioniz, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].cool_di_ioniz, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].cool_rr_ioniz, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].lum_rr_ioniz, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].cool_rr_metals_ioniz, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].lum_tot_ioniz, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].heat_shock, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (plasmamain[n].dmo_dt, 3, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].gain, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].converge_t_r, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].converge_t_e, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].converge_hc, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].trcheck, 1, MPI_INT, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].techeck, 1, MPI_INT, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].hccheck, 1, MPI_INT, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].converge_whole, 1, MPI_INT, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].converging, 1, MPI_INT, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (plasmamain[n].spec_mod_type, NXBANDS, MPI_INT, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (plasmamain[n].pl_alpha, NXBANDS, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (plasmamain[n].pl_log_w, NXBANDS, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (plasmamain[n].exp_temp, NXBANDS, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (plasmamain[n].exp_w, NXBANDS, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (plasmamain[n].fmin_mod, NXBANDS, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (plasmamain[n].fmax_mod, NXBANDS, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].ip, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].ip_direct, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].ip_scatt, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].xi, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&dt_e, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&dt_r, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&nmax_e, 1, MPI_INT, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&nmax_r, 1, MPI_INT, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
      }

      Log ("MPI task %d broadcasting plasma update information.\n", rank_global);
    }

    MPI_Barrier (MPI_COMM_WORLD);
    MPI_Bcast (commbuffer, size_of_commbuffer, MPI_PACKED, n_mpi, MPI_COMM_WORLD);
    MPI_Barrier (MPI_COMM_WORLD);
    Log ("MPI task %d survived broadcasting plasma update information.\n", rank_global);

    position = 0;

    if (rank_global != n_mpi)
    {
      MPI_Unpack (commbuffer, size_of_commbuffer, &position, &num_comm, 1, MPI_INT, MPI_COMM_WORLD);
      for (n_mpi2 = 0; n_mpi2 < num_comm; n_mpi2++)
      {
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &n, 1, MPI_INT, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].nwind, 1, MPI_INT, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].nplasma, 1, MPI_INT, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].ne, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].rho, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].vol, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, plasmamain[n].density, nions, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, plasmamain[n].partition, nions, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, plasmamain[n].levden, nlte_levels, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, plasmamain[n].PWdenom, nions, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, plasmamain[n].PWdtemp, nions, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, plasmamain[n].PWnumer, nions, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, plasmamain[n].PWntemp, nions, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].kappa_ff_factor, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].nscat_es, 1, MPI_INT, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, plasmamain[n].recomb_simple, nphot_total, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].kpkt_emiss, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].kpkt_abs, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, plasmamain[n].kbf_use, nphot_total, MPI_INT, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].kbf_nuse, 1, MPI_INT, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].t_r, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].t_r_old, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].t_e, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].t_e_old, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].dt_e, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].dt_e_old, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].heat_tot, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].abs_tot, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].heat_tot_old, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].heat_lines, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].heat_ff, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].heat_comp, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].heat_ind_comp, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].heat_lines_macro, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].heat_photo_macro, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].heat_photo, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].heat_auger, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].abs_photo, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].abs_auger, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].heat_z, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].w, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].ntot, 1, MPI_INT, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].ntot_star, 1, MPI_INT, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].ntot_bl, 1, MPI_INT, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].ntot_disk, 1, MPI_INT, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].ntot_wind, 1, MPI_INT, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].ntot_agn, 1, MPI_INT, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].mean_ds, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].n_ds, 1, MPI_INT, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].nrad, 1, MPI_INT, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].nioniz, 1, MPI_INT, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, plasmamain[n].ioniz, nions, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, plasmamain[n].recomb, nions, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, plasmamain[n].scatters, nions, MPI_INT, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, plasmamain[n].xscatters, nions, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, plasmamain[n].heat_ion, nions, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, plasmamain[n].cool_rr_ion, nions, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, plasmamain[n].lum_rr_ion, nions, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].j, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].j_direct, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].j_scatt, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].ave_freq, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].cool_tot, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, plasmamain[n].xj, NXBANDS, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, plasmamain[n].xave_freq, NXBANDS, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, plasmamain[n].xsd_freq, NXBANDS, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, plasmamain[n].nxtot, NXBANDS, MPI_INT, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].max_freq, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].lum_lines, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].lum_ff, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].cool_adiabatic, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].comp_nujnu, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].cool_comp, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].cool_dr, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].cool_di, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].cool_rr, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].lum_rr, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].cool_rr_metals, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].lum_rr_metals, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].lum_tot, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].lum_tot_old, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].cool_tot_ioniz, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].lum_lines_ioniz, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].lum_ff_ioniz, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].cool_adiabatic_ioniz, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].cool_comp_ioniz, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].cool_dr_ioniz, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].cool_di_ioniz, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].cool_rr_ioniz, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].lum_rr_ioniz, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].cool_rr_metals_ioniz, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].lum_tot_ioniz, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].heat_shock, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, plasmamain[n].dmo_dt, 3, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].gain, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].converge_t_r, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].converge_t_e, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].converge_hc, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].trcheck, 1, MPI_INT, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].techeck, 1, MPI_INT, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].hccheck, 1, MPI_INT, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].converge_whole, 1, MPI_INT, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].converging, 1, MPI_INT, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, plasmamain[n].spec_mod_type, NXBANDS, MPI_INT, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, plasmamain[n].pl_alpha, NXBANDS, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, plasmamain[n].pl_log_w, NXBANDS, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, plasmamain[n].exp_temp, NXBANDS, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, plasmamain[n].exp_w, NXBANDS, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, plasmamain[n].fmin_mod, NXBANDS, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, plasmamain[n].fmax_mod, NXBANDS, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].ip, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].ip_direct, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].ip_scatt, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].xi, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &dt_e_temp, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &dt_r_temp, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &nmax_e_temp, 1, MPI_INT, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &nmax_r_temp, 1, MPI_INT, MPI_COMM_WORLD);

        /* JM 1409 -- Altered for issue #110 to ensure correct reporting in parallel */
        if (fabs (dt_e_temp) >= fabs (dt_e))
        {
          /* Check if any other threads found a higher maximum for te */
          dt_e = dt_e_temp;
          nmax_e = nmax_e_temp;
        }

        if (fabs (dt_r_temp) >= fabs (dt_r))
        {
          /* Check if any other threads found a higher maximum for tr */
          dt_r = dt_r_temp;
          nmax_r = nmax_r_temp;
        }

        t_r_ave += plasmamain[n].t_r;
        t_e_ave += plasmamain[n].t_e;
        iave++;

      }

    }

  }
  free (commbuffer);
#endif


  /* Now we need to updated the densities immediately outside the wind so that the density interpolation in resonate will work.
     In this case all we have done is to copy the densities from the cell which is just in the wind (as one goes outward) to the
     cell that is just inside (or outside) the wind.

     SS asked whether we should also be extending the wind for other parameters, especially ne.  At present we do not interpolate
     on ne so this is not necessary.  If we did do that it would be required.

     In cylindrical coordinates, the fast dimension is z; grid positions increase up in z, and then out in r.
     In spherical polar coordinates, the fast dimension is theta; the grid increases in theta (measured)
     from the z axis), and then in r.
     In spherical coordinates, the grid increases as one might expect in r..
     *
   */

  for (ndom = 0; ndom < geo.ndomain; ndom++)
  {
    if (zdom[ndom].coord_type == CYLIND)
      cylind_extend_density (ndom, w);
    else if (zdom[ndom].coord_type == RTHETA)
      rtheta_extend_density (ndom, w);
    else if (zdom[ndom].coord_type == SPHERICAL)
      spherical_extend_density (ndom, w);
    else if (zdom[ndom].coord_type == CYLVAR)
      cylvar_extend_density (ndom, w);
    else
    {
      Error ("Wind_update2d: Unknown coordinate type %d for domain %d \n", zdom[ndom].coord_type, ndom);
      exit (0);
    }
  }

  /* Finished updating region outside of wind */

  num_updates++;
  strcpy (string, "");
  sprintf (string, "# Wind update: Number %d", num_updates);

  if (modes.zeus_connect == 1 && geo.hydro_domain_number > -1)  //If we are running in zeus connect mode - we open a file for heatcool rates
  {
    Log ("Outputting heatcool file for connecting to zeus\n");
    fptr = fopen ("py_heatcool.dat", "w");
    fprintf (fptr, "i j rcen thetacen vol temp xi ne heat_xray heat_comp heat_lines heat_ff cool_comp cool_lines cool_ff rho n_h\n");

  }

  /* Check the balance between the absorbed and the emitted flux */

  //NSH 0717 - first we need to ensure the cooling and luminosities reflect the current temperature

  cool_sum = wind_cooling ();   /*We call wind_cooling here to obtain an up to date set of cooling rates */
  lum_sum = wind_luminosity (0.0, VERY_BIG);    /*and we also call wind_luminosity to get the luminosities */


  xsum = psum = ausum = lsum = fsum = csum = icsum = apsum = aausum = abstot = 0;       //1108 NSH zero the new csum counter for compton heating

  for (nplasma = 0; nplasma < NPLASMA; nplasma++)
  {
    if (sane_check (plasmamain[nplasma].heat_tot))
      Error ("wind_update:sane_check w(%d).heat_tot is %e\n", nplasma, plasmamain[nplasma].heat_tot);
    if (sane_check (plasmamain[nplasma].heat_photo))
      Error ("wind_update:sane_check w(%d).heat_photo is %e\n", nplasma, plasmamain[nplasma].heat_photo);
    if (sane_check (plasmamain[nplasma].heat_auger))
      Error ("wind_update:sane_check w(%d).heat_auger is %e\n", nplasma, plasmamain[nplasma].heat_auger);
    if (sane_check (plasmamain[nplasma].heat_photo_macro))
      Error ("wind_update:sane_check w(%d).heat_photo_macro is %e\n", nplasma, plasmamain[nplasma].heat_photo_macro);
    if (sane_check (plasmamain[nplasma].heat_ff))
      Error ("wind_update:sane_check w(%d).heat_ff is %e\n", nplasma, plasmamain[nplasma].heat_ff);
    if (sane_check (plasmamain[nplasma].heat_lines))
      Error ("wind_update:sane_check w(%d).heat_lines is %e\n", nplasma, plasmamain[nplasma].heat_lines);
    if (sane_check (plasmamain[nplasma].heat_lines_macro))
      Error ("wind_update:sane_check w(%d).heat_lines_macro is %e\n", nplasma, plasmamain[nplasma].heat_lines_macro);
    /* 1108 NSH extra Sane check for compton heating */
    if (sane_check (plasmamain[nplasma].heat_comp))
      Error ("wind_update:sane_check w(%d).heat_comp is %e\n", nplasma, plasmamain[nplasma].heat_comp);

    abstot += plasmamain[nplasma].abs_tot;
    xsum += plasmamain[nplasma].heat_tot;
    psum += plasmamain[nplasma].heat_photo;
    ausum += plasmamain[nplasma].heat_auger;
    fsum += plasmamain[nplasma].heat_ff;
    lsum += plasmamain[nplasma].heat_lines;
    csum += plasmamain[nplasma].heat_comp;      //1108 NSH Increment the compton heating counter
    icsum += plasmamain[nplasma].heat_ind_comp; //1205 NSH Increment the induced compton heating counter
    apsum += plasmamain[nplasma].abs_photo;
    aausum += plasmamain[nplasma].abs_auger;

    /* JM130621- bugfix for windsave bug- needed so that we have the luminosities from ionization
       cycles in the windsavefile even if the spectral cycles are run */

    plasmamain[nplasma].cool_tot_ioniz = plasmamain[nplasma].cool_tot;
    plasmamain[nplasma].lum_ff_ioniz = plasmamain[nplasma].lum_ff;
    plasmamain[nplasma].cool_rr_ioniz = plasmamain[nplasma].cool_rr;
    plasmamain[nplasma].lum_rr_ioniz = plasmamain[nplasma].lum_rr;
    plasmamain[nplasma].cool_rr_metals_ioniz = plasmamain[nplasma].cool_rr_metals;
    plasmamain[nplasma].lum_lines_ioniz = plasmamain[nplasma].lum_lines;
    plasmamain[nplasma].cool_comp_ioniz = plasmamain[nplasma].cool_comp;
    plasmamain[nplasma].cool_dr_ioniz = plasmamain[nplasma].cool_dr;
    plasmamain[nplasma].cool_di_ioniz = plasmamain[nplasma].cool_di;
    plasmamain[nplasma].lum_tot_ioniz = plasmamain[nplasma].lum_tot;
    plasmamain[nplasma].cool_adiabatic_ioniz = plasmamain[nplasma].cool_adiabatic;

  }



  /* JM130621- bugfix for windsave bug- needed so that we have the luminosities from ionization
     cycles in the windsavefile even if the spectral cycles are run */
  geo.lum_ff_ioniz = geo.lum_ff;
  geo.cool_rr_ioniz = geo.cool_rr;
  geo.lum_rr_ioniz = geo.lum_rr;
  geo.lum_lines_ioniz = geo.lum_lines;
  geo.cool_comp_ioniz = geo.cool_comp;
  geo.cool_dr_ioniz = geo.cool_dr;
  geo.cool_di_ioniz = geo.cool_di;
  geo.cool_adiabatic_ioniz = geo.cool_adiabatic;
  geo.lum_disk_ioniz = geo.lum_disk;
  geo.lum_star_ioniz = geo.lum_star;
  geo.lum_bl_ioniz = geo.lum_bl;
  geo.lum_wind_ioniz = geo.lum_wind;
  geo.lum_tot_ioniz = geo.lum_tot;

  /* Added this system which counts number of times two situations occur (See #91)
     We only report these every 100,000 times (one can typically get ) */
  Log ("wind_update: note, errors from mean intensity can be high in a working model\n");
  Log
    ("wind_update: can be a problem with photon numbers if there are also errors from spectral_estimators and low photon number warnings\n");
  Log ("wind_update: mean_intensity: %8.4e occurrences, this cycle, this thread of 'no model exists in a band'\n", (1.0 * nerr_no_Jmodel));
  Log
    ("wind_update: mean intensity: %8.4e occurrences, this cycle, this thread of 'photon freq is outside frequency range of spectral model'\n",
     (1.0 * nerr_Jmodel_wrong_freq));


  /* zero the counters which record diagnositics from mean_intensity */
  nerr_Jmodel_wrong_freq = 0;
  nerr_no_Jmodel = 0;






  if (modes.zeus_connect == 1 && geo.hydro_domain_number > -1)  //If we are running in zeus connect mode, we output heating and cooling rates.
  {
    for (nwind = zdom[geo.hydro_domain_number].nstart; nwind < zdom[geo.hydro_domain_number].nstop; nwind++)
    {
      if (wmain[nwind].vol > 0.0)
      {
        nplasma = wmain[nwind].nplasma;
        wind_n_to_ij (geo.hydro_domain_number, plasmamain[nplasma].nwind, &i, &j);
        i = i - 1;              //There is a radial 'ghost zone' in python, we need to make our i,j agree with zeus
        vol = w[plasmamain[nplasma].nwind].vol;
        fprintf (fptr, "%d %d %e %e %e ", i, j, w[plasmamain[nplasma].nwind].rcen, w[plasmamain[nplasma].nwind].thetacen / RADIAN, vol);        //output geometric things
        fprintf (fptr, "%e %e %e ", plasmamain[nplasma].t_e, plasmamain[nplasma].xi, plasmamain[nplasma].ne);   //output temp, xi and ne to ease plotting of heating rates
        fprintf (fptr, "%e ", (plasmamain[nplasma].heat_photo + plasmamain[nplasma].heat_auger) / vol); //Xray heating - or photoionization
        fprintf (fptr, "%e ", (plasmamain[nplasma].heat_comp) / vol);   //Compton heating
        fprintf (fptr, "%e ", (plasmamain[nplasma].heat_lines) / vol);  //Line heating 28/10/15 - not currently used in zeus
        fprintf (fptr, "%e ", (plasmamain[nplasma].heat_ff) / vol);     //FF heating 28/10/15 - not currently used in zeus
        fprintf (fptr, "%e ", (plasmamain[nplasma].cool_comp) / vol);   //Compton cooling
        fprintf (fptr, "%e ", (plasmamain[nplasma].lum_lines + plasmamain[nplasma].cool_rr + plasmamain[nplasma].cool_dr) / vol);       //Line cooling must include all recombination cooling
        fprintf (fptr, "%e ", (plasmamain[nplasma].lum_ff) / vol);      //ff cooling
        fprintf (fptr, "%e ", plasmamain[nplasma].rho); //density
        fprintf (fptr, "%e\n", plasmamain[nplasma].rho * rho2nh);       //hydrogen number density
      }
    }
    fclose (fptr);
  }
  else if (modes.zeus_connect == 1 && geo.hydro_domain_number < 0)
  {
    Error ("wind_updates2d:  Attempting to access a hydro domain in a non hydro run - not writing out hydro file\n");
  }

  /* The lines differ only in that Wind_heating adds mechanical heating, that is adiabatic heating */

  Log
    ("!!wind_update: Absorbed flux    %8.2e  (photo %8.2e ff %8.2e compton %8.2e auger %8.2e induced_compton %8.2e lines %8.2e)\n",
     abstot, apsum, fsum, csum, aausum, icsum, lsum);

  Log
    ("!!wind_update: Wind heating     %8.2e  (photo %8.2e ff %8.2e compton %8.2e auger %8.2e induced_compton %8.2e lines %8.2e adiabatic %8.2e)\n",
     xsum + geo.heat_adiabatic, psum, fsum, csum, ausum, icsum, lsum, geo.heat_adiabatic);

  /* 1108 NSH added commands to report compton cooling 1110 removed,
   * As was the case above, there are two almost identical lines.  Wind_cooling includes processes that do not produce photons,
   * not-only adiabatic cooling, but also goe.cool_comp, geo_cool_dr and geo.cool_di */
  Log
    ("!!wind_update: Wind luminosity  %8.2e (recomb %8.2e ff %8.2e lines %8.2e) after update\n",
     lum_sum, geo.lum_rr, geo.lum_ff, geo.lum_lines);


  rad_sum = wind_luminosity (xband.f1[0], xband.f2[xband.nbands - 1]);  /*and we also call wind_luminosity to get the luminosities */

  Log
    ("!!wind_update: Rad  luminosity  %8.2e (recomb %8.2e ff %8.2e lines %8.2e) after update\n",
     rad_sum, geo.lum_rr, geo.lum_ff, geo.lum_lines);

  Log
    ("!!wind_update: Wind cooling     %8.2e (recomb %8.2e ff %8.2e compton %8.2e DR %8.2e DI %8.2e lines %8.2e adiabatic %8.2e) after update\n",
     cool_sum, geo.cool_rr, geo.lum_ff, geo.cool_comp, geo.cool_dr, geo.cool_di, geo.lum_lines, geo.cool_adiabatic);


  /* Print out some diagnostics of the changes in the wind update */

  if (modes.zeus_connect == 1 || modes.fixed_temp == 1) //There is no point in computing temperature changes, because we have fixed them!
  {
    Log ("!!wind_update: We are running in fixed temperature mode - no temperature report\n");
  }
  else
  {
    t_r_ave_old /= iave;
    t_e_ave_old /= iave;
    t_r_ave /= iave;
    t_e_ave /= iave;

    if (nmax_r != -1)
    {
      wind_n_to_ij (wmain[nmax_r].ndom, nmax_r, &i, &j);
      Log ("!!wind_update: Max change in t_r %6.0f at cell %4d (%d,%d)\n", dt_r, nmax_r, i, j);
      Log ("!!wind_update: Ave change in t_r %6.0f from %6.0f to %6.0f\n", (t_r_ave - t_r_ave_old), t_r_ave_old, t_r_ave);
    }
    else
      Log ("!!wind_update: t_r did not change in any cells this cycle\n");

    if (nmax_e != -1)
    {
      wind_n_to_ij (wmain[nmax_e].ndom, nmax_e, &i, &j);
      Log ("!!wind_update: Max change in t_e %6.0f at cell %4d (%d,%d)\n", dt_e, nmax_e, i, j);
      Log ("!!wind_update: Ave change in t_e %6.0f from %6.0f to %6.0f\n", (t_e_ave - t_e_ave_old), t_e_ave_old, t_e_ave);
    }
    else
      Log ("!!wind_update: t_e did not change in any cells this cycle\n");


    Log ("Summary  t_r  %6.0f   %6.0f  #t_r and dt_r on this update\n", t_r_ave, (t_r_ave - t_r_ave_old));
    Log ("Summary  t_e  %6.0f   %6.0f  #t_e and dt_e on this update\n", t_e_ave, (t_e_ave - t_e_ave_old));
  }

  check_convergence ();

  /* Summarize the radiative temperatures (ksl 04 mar) */

  xtemp_rad (w);

/* This next block is to allow the output of data relating to the abundances of ions when python is being tested
 * with thin shell mode.We will only want this to run if the wind mode is 9, for test or thin shell mode.
 *
 * Note that this section is very dependent on the peculiar structure of the single shell model, which has only
 * one element (namely element 2) in the wind.  nshell below correspond to that particular plasma shell. Recognizing
 * this was a key element to solving bug #412.
 */
  for (ndom = 0; ndom < geo.ndomain; ndom++)
  {

    if (zdom[ndom].wind_type == SHELL)
    {

      /* nshell is the plasma cell that correspond to the second wind cell for the shell_wind model */
      nshell = wmain[zdom[ndom].nstart + 1].nplasma;
      n = plasmamain[nshell].nwind;
      for (i = 0; i < geo.nxfreq; i++)
      {                         /*loop over number of bands */
        Log
          ("Band %i f1 %e f2 %e model %i pl_alpha %f pl_log_w %e exp_t %e exp_w %e\n",
           i, geo.xfreq[i], geo.xfreq[i + 1],
           plasmamain[nshell].spec_mod_type[i],
           plasmamain[nshell].pl_alpha[i], plasmamain[nshell].pl_log_w[i], plasmamain[nshell].exp_temp[i], plasmamain[nshell].exp_w[i]);
      }
      /* Get some line diagnostics */

      lum_h_line = 0.0;
      lum_he_line = 0.0;
      lum_c_line = 0.0;
      lum_n_line = 0.0;
      lum_o_line = 0.0;
      lum_fe_line = 0.0;

      for (i = 0; i < nlines; i++)
      {
        if (lin_ptr[i]->z == 1)
          lum_h_line = lum_h_line + lin_ptr[i]->pow;
        else if (lin_ptr[i]->z == 2)
          lum_he_line = lum_he_line + lin_ptr[i]->pow;
        else if (lin_ptr[i]->z == 6)
          lum_c_line = lum_c_line + lin_ptr[i]->pow;
        else if (lin_ptr[i]->z == 7)
          lum_n_line = lum_n_line + lin_ptr[i]->pow;
        else if (lin_ptr[i]->z == 8)
          lum_o_line = lum_o_line + lin_ptr[i]->pow;
        else if (lin_ptr[i]->z == 26)
          lum_fe_line = lum_fe_line + lin_ptr[i]->pow;
      }
      agn_ip = geo.const_agn * (((pow (50000 / HEV, geo.alpha_agn + 1.0)) - pow (100 / HEV, geo.alpha_agn + 1.0)) / (geo.alpha_agn + 1.0));
      agn_ip /= (w[n].r * w[n].r);
      agn_ip /= plasmamain[nshell].rho * rho2nh;
      /* Report luminosities, IP and other diagnositic quantities */
      Log
        ("OUTPUT Lum_agn= %e T_e= %e N_h= %e N_e= %e alpha= %f IP(sim_2010)= %e Measured_IP(cloudy)= %e Measured_Xi= %e distance= %e volume= %e mean_ds=%e\n",
         geo.lum_agn, plasmamain[nshell].t_e,
         plasmamain[nshell].rho * rho2nh, plasmamain[nshell].ne,
         geo.alpha_agn, agn_ip, plasmamain[nshell].ip,
         plasmamain[nshell].xi, w[n].r, w[n].vol, plasmamain[nshell].mean_ds / plasmamain[nshell].n_ds);

      /* 1108 NSH Added commands to report compton heating */
      Log
        ("OUTPUT Absorbed_flux(ergs-1cm-3)    %8.2e  (photo %8.2e ff %8.2e compton %8.2e induced_compton %8.2e lines %8.2e auger %8.2e )\n",
         xsum / w[n].vol, psum / w[n].vol, fsum / w[n].vol, csum / w[n].vol, icsum / w[n].vol, lsum / w[n].vol, ausum / w[n].vol);

      /* 1110 NSH Added this line to report all cooling mechanisms, including those that do not generate photons. */
      Log
        ("OUTPUT Wind_cooling(ergs-1cm-3)     %8.2e (recomb %8.2e ff %8.2e compton %8.2e DR %8.2e DI %8.2e adiabatic %8.2e lines %8.2e ) after update\n",
         cool_sum / w[n].vol, geo.cool_rr / w[n].vol,
         geo.lum_ff / w[n].vol, geo.cool_comp / w[n].vol,
         geo.cool_dr / w[n].vol, geo.cool_di / w[n].vol, geo.cool_adiabatic / w[n].vol, geo.lum_lines / w[n].vol);
      Log
        ("OUTPUT Wind_luminosity(ergs-1cm-3)     %8.2e (recomb %8.2e ff %8.2e lines %8.2e ) after update\n",
         lum_sum / w[n].vol, geo.lum_rr / w[n].vol, geo.lum_ff / w[n].vol, geo.lum_lines / w[n].vol);
      /* NSH 1701 calculate the recombination cooling for other elements */

      c_rec = n_rec = o_rec = fe_rec = 0.0;
      c_lum = n_lum = o_lum = fe_lum = 0.0;
      h_dr = he_dr = c_dr = n_dr = o_dr = fe_dr = 0.0;
      cool_dr_metals = 0.0;

      for (nn = 0; nn < nions; nn++)
      {
        if (ion[nn].z == 6)
        {
          c_dr = c_dr + plasmamain[nshell].cool_dr_ion[nn];
          c_rec = c_rec + plasmamain[nshell].cool_rr_ion[nn];
          c_lum = c_lum + plasmamain[nshell].lum_rr_ion[nn];

        }
        if (ion[nn].z == 7)
        {
          n_dr = n_dr + plasmamain[nshell].cool_dr_ion[nn];
          n_rec = n_rec + plasmamain[nshell].cool_rr_ion[nn];
          n_lum = n_lum + plasmamain[nshell].lum_rr_ion[nn];
        }
        if (ion[nn].z == 8)
        {
          o_dr = o_dr + plasmamain[nshell].cool_dr_ion[nn];
          o_rec = o_rec + plasmamain[nshell].cool_rr_ion[nn];
          o_lum = o_lum + plasmamain[nshell].lum_rr_ion[nn];
        }
        if (ion[nn].z == 26)
        {
          fe_dr = fe_dr + plasmamain[nshell].cool_dr_ion[nn];
          fe_rec = fe_rec + plasmamain[nshell].cool_rr_ion[nn];
          fe_lum = fe_lum + plasmamain[nshell].lum_rr_ion[nn];
        }
        if (ion[nn].z > 2)
          cool_dr_metals = cool_dr_metals + plasmamain[nshell].cool_dr_ion[nn];
      }

      Log ("OUTPUT Wind_line_cooling(ergs-1cm-3)  H %8.2e He %8.2e C %8.2e N %8.2e O %8.2e Fe %8.2e Metals %8.2e\n", lum_h_line / w[n].vol,
           lum_he_line / w[n].vol, lum_c_line / w[n].vol, lum_n_line / w[n].vol, lum_o_line / w[n].vol, lum_fe_line / w[n].vol);
      Log ("OUTPUT Wind_recomb_cooling(ergs-1cm-3)  H %8.2e He %8.2e C %8.2e N %8.2e O %8.2e Fe %8.2e Metals %8.2e\n",
           plasmamain[nshell].cool_rr_ion[0] / w[n].vol, (plasmamain[nshell].cool_rr_ion[2] + plasmamain[nshell].cool_rr_ion[3]) / w[n].vol,
           c_rec / w[n].vol, n_rec / w[n].vol, o_rec / w[n].vol, fe_rec / w[n].vol, plasmamain[nshell].cool_rr_metals / w[n].vol);
      Log ("OUTPUT Wind_recomb_lum(ergs-1cm-3)  H %8.2e He %8.2e C %8.2e N %8.2e O %8.2e Fe %8.2e Metals %8.2e\n",
           plasmamain[nshell].lum_rr_ion[0] / w[n].vol, (plasmamain[nshell].lum_rr_ion[2] + plasmamain[nshell].lum_rr_ion[3]) / w[n].vol,
           c_lum / w[n].vol, n_lum / w[n].vol, o_lum / w[n].vol, fe_lum / w[n].vol, plasmamain[nshell].lum_rr_metals / w[n].vol);
      Log ("OUTPUT Wind_dr_cooling(ergs-1cm-3)  H %8.2e He %8.2e C %8.2e N %8.2e O %8.2e Fe %8.2e Metals %8.2e\n",
           plasmamain[nshell].cool_dr_ion[0] / w[n].vol, (plasmamain[nshell].cool_dr_ion[2] + plasmamain[nshell].cool_dr_ion[3]) / w[n].vol,
           c_dr / w[n].vol, n_dr / w[n].vol, o_dr / w[n].vol, fe_dr / w[n].vol, cool_dr_metals / w[n].vol);
      /* 1110 NSH Added this line to report all cooling mechanisms, including those that do not generate photons. */
      Log
        ("OUTPUT Balance      Cooling=%8.2e Heating=%8.2e Lum=%8.2e T_e=%e after update\n",
         cool_sum, xsum, lum_sum, plasmamain[nshell].t_e);

      for (n = 0; n < nelements; n++)
      {
        first = ele[n].firstion;
        last = first + ele[n].nions;
        Log ("OUTPUT %-5s ", ele[n].name);
        tot = 0;
        for (m = first; m < last; m++)
          tot += plasmamain[nshell].density[m];
        for (m = first; m < last; m++)
        {
          Log (" %8.2e", plasmamain[nshell].density[m] / tot);
        }
        Log ("\n");
      }
    }
  }


  return (0);
}




/**********************************************************/
/**
 * @brief      zeros those portions of the wind which contain the radiation properties
 * 	of the wind, i.e those portions which should be set to zeroed when the structure of the
 * 	wind has been changed or when you simply want to start off a calculation in a known state
 *
 * @return    Always returns 0
 *
 * @details
 * The routine is called at the beginning of each ionization calculation
 * cycle.  It should zero all heating and radiation induced cooling in the wind array.  Since
 * cooling is recalculated in wind_update, one needs to be sure that all of the appropriate
 * cooling terms are also rezeroed there as well.
 *
 * ### Notes ###
 *
 **********************************************************/

int
wind_rad_init ()
{
  int n, i;
  int njump;


  for (n = 0; n < NPLASMA; n++)
  {
    plasmamain[n].j = plasmamain[n].ave_freq = plasmamain[n].ntot = 0;
    plasmamain[n].j_direct = plasmamain[n].j_scatt = 0.0;       //NSH 1309 zero j banded by number of scatters
    plasmamain[n].ip = 0.0;
    plasmamain[n].xi = 0.0;
    plasmamain[n].ip_direct = plasmamain[n].ip_scatt = 0.0;
    plasmamain[n].mean_ds = 0.0;
    plasmamain[n].n_ds = 0;
    plasmamain[n].ntot_disk = plasmamain[n].ntot_agn = 0;       //NSH 15/4/11 counters to see where photons come from
    plasmamain[n].ntot_star = plasmamain[n].ntot_bl = plasmamain[n].ntot_wind = 0;
    plasmamain[n].heat_tot = plasmamain[n].heat_ff = plasmamain[n].heat_photo = plasmamain[n].heat_lines = 0.0;
    plasmamain[n].abs_tot = plasmamain[n].abs_auger = plasmamain[n].abs_photo = 0.0;

    plasmamain[n].heat_z = 0.0;
    plasmamain[n].max_freq = 0.0;       //NSH 120814 Zero the counter which works out the maximum frequency seen in a cell and hence the maximum applicable frequency of the power law estimators.
    plasmamain[n].cool_tot = plasmamain[n].lum_tot = plasmamain[n].lum_lines = plasmamain[n].lum_ff = 0.0;
    plasmamain[n].cool_rr = plasmamain[n].cool_rr_metals = plasmamain[n].lum_rr = 0.0;
    plasmamain[n].nrad = plasmamain[n].nioniz = 0;
    plasmamain[n].comp_nujnu = -1e99;   //1701 NSH Zero the integrated specific intensity for the cell
    plasmamain[n].cool_comp = 0.0;      //1108 NSH Zero the compton luminosity for the cell
    plasmamain[n].heat_comp = 0.0;      //1108 NSH Zero the compton heating for the cell
    plasmamain[n].heat_ind_comp = 0.0;  //1108 NSH Zero the induced compton heating for the cell
    plasmamain[n].heat_auger = 0.0;     //1108 NSH Zero the auger heating for the cell

    if (nlevels_macro > 1 && geo.macro_simple == 0)
      macromain[n].kpkt_rates_known = -1;

/* 1108 NSH Loop to zero the frequency banded radiation estimators */
/* 71 - 111279 - ksl - Small modification to reflect the fact that nxfreq has been moved into the geo structure */
    for (i = 0; i < geo.nxfreq; i++)
    {
      plasmamain[n].xj[i] = plasmamain[n].xave_freq[i] = plasmamain[n].nxtot[i] = 0;
      plasmamain[n].xsd_freq[i] = 0.0;  /* NSH 120815 Zero the standard deviation counter */
      plasmamain[n].fmin[i] = geo.xfreq[i + 1]; /* Set the minium frequency to the max frequency in the band */
      plasmamain[n].fmax[i] = geo.xfreq[i];     /* Set the maximum frequency to the min frequency in the band */
    }




    for (i = 0; i < nions; i++)
    {
      plasmamain[n].ioniz[i] = plasmamain[n].recomb[i] = plasmamain[n].heat_ion[i] = plasmamain[n].cool_rr_ion[i] =
        plasmamain[n].lum_rr_ion[i] = 0.0;

    }

    /*Block added (Dec 08) to zero the auger rate estimators */
    /* commented out by NSH 2018 - removed code */
//    for (i = 0; i < nauger; i++)
//    {
//      plasmamain[n].gamma_inshl[i] = 0.0;
//    }

    /* Next blocks added by SS Mar 2004 to zero the Macro Atom estimators. */

    /* 57h -- 0608 -- These sections actually involve enough calculations that
       they are noticeable in term sof the overall speed.  One would if possible
       like to avoid this section, since it requires the creation of macromain,
       even though macromain is not used -- ksl */


    for (i = 0; i < nlevels_macro; i++) //57h
    {
      for (njump = 0; njump < config[i].n_bbu_jump; njump++)
      {
        macromain[n].jbar[config[i].bbu_indx_first + njump] = 0.0;      // mean intensity
      }
      for (njump = 0; njump < config[i].n_bfu_jump; njump++)
      {
        macromain[n].gamma[config[i].bfu_indx_first + njump] = 0.0;
        macromain[n].gamma_e[config[i].bfu_indx_first + njump] = 0.0;
        macromain[n].alpha_st[config[i].bfd_indx_first + njump] = 0.0;  //stimulated recombination
        macromain[n].alpha_st_e[config[i].bfd_indx_first + njump] = 0.0;
      }


      /* Next block to set spontaneous recombination rates for next iteration. (SS July 04) */
      for (njump = 0; njump < config[i].n_bfd_jump; njump++)
      {
        if (plasmamain[n].t_e > 1.0)
        {
          //04Jul--ksl-modified these calls to reflect changed alpha_sp
          macromain[n].recomb_sp[config[i].bfd_indx_first + njump] = alpha_sp (&phot_top[config[i].bfd_jump[njump]], &plasmamain[n], 0);
          macromain[n].recomb_sp_e[config[i].bfd_indx_first + njump] = alpha_sp (&phot_top[config[i].bfd_jump[njump]], &plasmamain[n], 2);
        }
        else
        {
          macromain[n].recomb_sp[config[i].bfd_indx_first + njump] = 0.0;
          macromain[n].recomb_sp_e[config[i].bfd_indx_first + njump] = 0.0;
        }

      }
    }

    for (i = 0; i < ntop_phot; i++)
    {
      /* 57h -- recomb_simple is only required for we are using a macro atom approach, and only non-zero when
         this particular phot_tob xsection is treated as a simple x-section. Stuart, is this correct?? I've added
         checks so that macro_info is only 0 (false) or true (1), and so the logic of the next section can be
         simplified.  0608-ksl */
      if (geo.macro_simple || phot_top[i].macro_info)
      {

        plasmamain[n].recomb_simple[i] = 0.0;
      }
      else
      {                         //we want a macro approach, but not for this ion so need recomb_simple
        plasmamain[n].recomb_simple[i] = alpha_sp (&phot_top[i], &plasmamain[n], 2);
      }
    }


    //zero the emissivities that are needed for the spectral synthesis step.
    plasmamain[n].kpkt_emiss = 0.0;
    plasmamain[n].kpkt_abs = 0.0;
    for (i = 0; i < nlevels_macro; i++) //57h
    {
      macromain[n].matom_abs[i] = 0.0;

      macromain[n].matom_emiss[i] = 0.0;

    }

    /* End of added material. */
  }


  return (0);
}
