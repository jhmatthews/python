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
 *
 **********************************************************/
int
wind_update (WindPtr w)
{
  int n, i, j;
  double xsum, psum, fsum, lsum, csum, icsum, ausum, chexsum;
  double cool_sum, lum_sum, rad_sum;    //1706 - the total cooling and luminosity of the wind
  double apsum, aausum, abstot; //Absorbed photon energy from PI and auger
  double flux_persist_scale;
  double volume;
  double t_r_old, t_e_old, dt_r, dt_e;
  double t_r_ave_old, t_r_ave, t_e_ave_old, t_e_ave;
  int iave, nmax_r, nmax_e;
  int nplasma;
  int nwind;
  int my_nmin, my_nmax;         //Note that these variables are still used even without MPI on
  int ndom;

  dt_r = dt_e = 0.0;
  iave = 0;
  nmax_r = nmax_e = -1;
  t_r_ave_old = t_r_ave = t_e_ave_old = t_e_ave = 0.0;

  /* For MPI parallelisation, the following loop will be distributed over mutiple tasks.
     Note that the mynmim and mynmax variables are still used even without MPI on */

#ifdef MPI_ON
  int ndo = get_parallel_nrange (rank_global, NPLASMA, np_mpi_global, &my_nmin, &my_nmax);
#else
  my_nmin = 0;
  my_nmax = NPLASMA;
#endif

  /* Start with a call to the routine which normalises all the macro atom
     monte carlo radiation field estimators. It's best to do this first since
     some of the estimators include temperature terms (stimulated correction
     terms) which were included during the monte carlo simulation so we want
     to be sure that the SAME temperatures are used here. (SS - Mar 2004). */

  for (n = 0; n < NPLASMA; n++)
  {
    /* normalise macro estimators: TODO: add to comm buffer */
    if (geo.rt_mode == RT_MODE_MACRO && geo.macro_simple == FALSE)
    {
      nwind = plasmamain[n].nwind;
      normalise_macro_estimators (nwind);       // todo: update to use nplasma instead

      /* force recalculation of kpacket rates and matrices, if applicable */
      macromain[n].kpkt_rates_known = FALSE;
      macromain[n].matrix_rates_known = FALSE;
    }
  }

  /* we now know how many cells this thread has to process - note this will be
     0-NPLASMA in serial mode */

  flux_persist_scale = 0.5;     //The ammount of the latest flux that gets added into the persistent flux

  for (n = my_nmin; n < my_nmax; n++)
  {
    nwind = plasmamain[n].nwind;
    volume = w[nwind].vol;

    /* Skip cells that are partially in the wind these are not to be included
       in the calculation */

    if (modes.partial_cells == PC_EXTEND && wmain[nwind].inwind == W_PART_INWIND)
    {
      continue;
    }

    if (plasmamain[n].ntot < 100)
    {
      Log
        ("!!wind_update: Cell %4d Dom %d  Vol. %8.2e r %8.2e theta %8.2e has only %4d photons\n",
         n, w[nwind].ndom, volume, w[nwind].rcen, w[nwind].thetacen, plasmamain[n].ntot);
    }

    /* this routine normalises the unbanded and banded estimators for simple atoms in this cell */
    normalise_simple_estimators (&plasmamain[n]);

    /* If geo.adiabatic is true, then calculate the adiabatic cooling using the current, i.e
     * previous value of t_e.  Note that this may not be  best way to determine the cooling.
     * Changes made here should also be reflected in wind2d.c.  At present, adiabatic cooling
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

    /* update the persistent fluxes */
//    update_persistent_directional_flux_estimators (n, flux_persist_scale);      // TODO: this needs communicating as well
  }

  /*This is the end of the update loop that is parallised. We now need to exchange data between the tasks. */
  communicate_plasma_cells (my_nmin, my_nmax);

  /* Each rank now has the updated plasma cells, so we can now find out what the max d_t is
   * in the wind */
  for (n = 0; n < NPLASMA; ++n)
  {
    update_persistent_directional_flux_estimators (n, flux_persist_scale);      // TODO: put into parallel loop
    if ((fabs (plasmamain[n].t_r_old - plasmamain[n].t_r)) > fabs (dt_r))
    {
      dt_r = plasmamain[n].t_r - plasmamain[n].t_r_old;
      nmax_r = n;
    }
    if ((fabs (plasmamain[n].t_e_old - plasmamain[n].t_e)) > fabs (dt_e))
    {
      dt_e = plasmamain[n].t_e - plasmamain[n].t_e_old;
      nmax_e = n;
    }
    t_r_ave += plasmamain[n].t_r;
    t_e_ave += plasmamain[n].t_e;
    t_r_ave_old += plasmamain[n].t_r_old;
    t_e_ave_old += plasmamain[n].t_e_old;
//    iave++;                     // TODO: do we need to do this? iave will be NPLASMA
  }

  t_r_ave /= NPLASMA;
  t_e_ave /= NPLASMA;
  t_r_ave_old /= NPLASMA;
  t_e_ave_old /= NPLASMA;

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
      Exit (0);
    }
  }
  /* Finished updating region outside of wind */

  /* Check the balance between the absorbed and the emitted flux */
  /* NSH 0717 - ensure the cooling and luminosities reflect the current temperature */

  cool_sum = wind_cooling ();   /* We call wind_cooling here to obtain an up to date set of cooling rates */
  lum_sum = wind_luminosity (0.0, VERY_BIG, MODE_CMF_TIME);     /* and we also call wind_luminosity to get the luminosities */

  xsum = psum = ausum = lsum = fsum = csum = icsum = apsum = aausum = abstot = chexsum = 0;     //1108 NSH zero the new csum counter for compton heating

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
    csum += plasmamain[nplasma].heat_comp;
    icsum += plasmamain[nplasma].heat_ind_comp;
    apsum += plasmamain[nplasma].abs_photo;
    aausum += plasmamain[nplasma].abs_auger;
    chexsum += plasmamain[nplasma].heat_ch_ex;


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
  Log ("wind_update: mean_intensity: %8.4e occurrences, this cycle, this thread of 'no model exists in a band'\n", nerr_no_Jmodel);
  Log
    ("wind_update: mean intensity: %8.4e occurrences, this cycle, this thread of 'photon freq is outside frequency range of spectral model'\n",
     nerr_Jmodel_wrong_freq);


  /* zero the counters which record diagnositics from mean_intensity */
  nerr_Jmodel_wrong_freq = 0;
  nerr_no_Jmodel = 0;

  if (modes.zeus_connect == 1 && geo.hydro_domain_number > -1)  //If we are running in zeus connect mode, we output heating and cooling rates.
  {
    create_hydro_output_files ();
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


  rad_sum = wind_luminosity (xband.f1[0], xband.f2[xband.nbands - 1], MODE_CMF_TIME);

  Log
    ("!!wind_update: Rad luminosity  %8.2e (recomb %8.2e ff %8.2e lines %8.2e) after update\n",
     rad_sum, geo.lum_rr, geo.lum_ff, geo.lum_lines);

  Log
    ("!!wind_update: Wind cooling     %8.2e (recomb %8.2e ff %8.2e compton %8.2e DR %8.2e DI %8.2e lines %8.2e adiabatic %8.2e) after update\n",
     cool_sum, geo.cool_rr, geo.lum_ff, geo.cool_comp, geo.cool_dr, geo.cool_di, geo.lum_lines, geo.cool_adiabatic);

  if (!modes.turn_off_upweighting_of_simple_macro_atoms)
  {
    /* If we have "indivisible packet" mode on but are using the
       new BF_SIMPLE_EMISSIVITY_APPROACH then we report the flows into and out of the ion pool */
    if (geo.rt_mode == RT_MODE_MACRO)
      report_bf_simple_ionpool ();
  }

  if (modes.zeus_connect != 1 || modes.fixed_temp != 1)
  {
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
  shell_wind_diagnostics (xsum, psum, fsum, csum, icsum, lsum, ausum, chexsum, cool_sum, lum_sum);

  return (0);
}

/**********************************************************/
/**
 * @brief This summarises the flows into and out of the ionization pool for
 *        simple ions in RT_MODE_MACRO
 *
 * @return    Always returns 0
 *
 **********************************************************/
int
report_bf_simple_ionpool ()
{
  int n, m;
  int in_tot, out_tot;
  double total_in = 0.0;
  double total_out = 0.0;

  for (n = 0; n < NPLASMA; n++)
  {
    total_in += plasmamain[n].bf_simple_ionpool_in;
    total_out += plasmamain[n].bf_simple_ionpool_out;

    if (plasmamain[n].bf_simple_ionpool_out > plasmamain[n].bf_simple_ionpool_in)
    {
      Error ("The net flow out of simple ion pool (%8.4e) > than the net flow in (%8.4e) in cell %d\n",
             plasmamain[n].bf_simple_ionpool_out, plasmamain[n].bf_simple_ionpool_in, n);
    }
  }

  Log ("!! report_bf_simple_ionpool: Total flow into: %8.4e and out of: %8.4e bf_simple ion pool\n", total_in, total_out);

  total_in = total_out = 0;
  for (m = 0; m < nphot_total; m++)
  {
    in_tot = out_tot = 0;
    for (n = 0; n < NPLASMA; n++)
    {
      in_tot += plasmamain[n].n_bf_in[m];
      out_tot += plasmamain[n].n_bf_out[m];
    }


    Log ("!! report_bf:  %3d   %3d %3d %7d  %7d\n", m, phot_top[m].z, phot_top[m].istate, in_tot, out_tot);

    total_in += in_tot;
    total_out += out_tot;
  }

  Log ("!! report_bf tots:   %10.0f  %10.0f\n", total_in, total_out);


  return (0);
}


/**********************************************************/
/**
 * @brief
 *
 * @details
 *
 * ### Notes ###
 *
 **********************************************************/

static void
init_plasma (void)
{
  int i;
  int j;

  for (i = 0; i < NPLASMA; ++i)
  {
    /* Start by initialising integer fields */
    plasmamain[i].j = 0;
    plasmamain[i].ave_freq = 0;
    plasmamain[i].ntot = 0;
    plasmamain[i].n_ds = 0;
    plasmamain[i].ntot_disk = 0;
    plasmamain[i].ntot_agn = 0;
    plasmamain[i].ntot_star = 0;
    plasmamain[i].ntot_bl = 0;
    plasmamain[i].ntot_wind = 0;
    plasmamain[i].nrad = 0;
    plasmamain[i].nioniz = 0;
    for (j = 0; j < nphot_total; j++)
    {
      plasmamain[i].n_bf_in[j] = 0;
      plasmamain[i].n_bf_out[j] = 0;
    }

    /* Next we'll initialise the rest of the fields, which are doubles */
    plasmamain[i].j_direct = 0.0;
    plasmamain[i].j_scatt = 0.0;
    plasmamain[i].ip = 0.0;
    plasmamain[i].xi = 0.0;
    plasmamain[i].ip_direct = 0.0;
    plasmamain[i].ip_scatt = 0.0;
    plasmamain[i].mean_ds = 0.0;
    plasmamain[i].heat_tot = 0.0;
    plasmamain[i].heat_ff = 0.0;
    plasmamain[i].heat_photo = 0.0;
    plasmamain[i].heat_lines = 0.0;
    plasmamain[i].abs_tot = 0.0;
    plasmamain[i].abs_auger = 0.0;
    plasmamain[i].abs_photo = 0.0;
    plasmamain[i].heat_z = 0.0;
    plasmamain[i].max_freq = 0.0;
    plasmamain[i].cool_tot = 0.0;
    plasmamain[i].lum_tot = 0.0;
    plasmamain[i].lum_lines = 0.0;
    plasmamain[i].lum_ff = 0.0;
    plasmamain[i].cool_rr = 0.0;
    plasmamain[i].cool_rr_metals = 0.0;
    plasmamain[i].lum_rr = 0.0;
    plasmamain[i].comp_nujnu = -1e99;
    plasmamain[i].cool_comp = 0.0;
    plasmamain[i].heat_comp = 0.0;
    plasmamain[i].heat_ind_comp = 0.0;
    plasmamain[i].heat_auger = 0.0;
    plasmamain[i].heat_ch_ex = 0.0;
    plasmamain[i].bf_simple_ionpool_out = 0.0;
    plasmamain[i].bf_simple_ionpool_in = 0.0;

    for (j = 0; j < NUM_RAD_FORCE_DIRECTIONS; j++)
    {
      plasmamain[i].dmo_dt[j] = 0.0;
    }
    for (j = 0; j < NUM_FORCE_EST_DIRECTIONS; j++)
    {
      plasmamain[i].rad_force_es[j] = 0.0;
      plasmamain[i].rad_force_ff[j] = 0.0;
      plasmamain[i].rad_force_bf[j] = 0.0;
      plasmamain[i].F_vis[j] = 0.0;
      plasmamain[i].F_UV[j] = 0.0;
      plasmamain[i].F_Xray[j] = 0.0;
      if (geo.wcycle == 0)      // Persistent values, so only initialise for first ionisation cycle
      {
        plasmamain[i].F_vis_persistent[j] = 0.0;
        plasmamain[i].F_UV_persistent[j] = 0.0;
        plasmamain[i].F_Xray_persistent[j] = 0.0;
        plasmamain[i].rad_force_bf_persist[j] = 0.0;
      }
    }
    for (j = 0; j < NFLUX_ANGLES; j++)
    {
      if (geo.wcycle == 0)      // Persistent values, so only initialise for first ionisation cycle
      {
        plasmamain[i].F_UV_ang_x_persist[j] = 0.0;
        plasmamain[i].F_UV_ang_y_persist[j] = 0.0;
        plasmamain[i].F_UV_ang_z_persist[j] = 0.0;
      }
      plasmamain[i].F_UV_ang_x[j] = 0.0;
      plasmamain[i].F_UV_ang_y[j] = 0.0;
      plasmamain[i].F_UV_ang_z[j] = 0.0;
    }

    /* Initialise  the frequency banded radiation estimators used for estimating the coarse spectra in each i */
    for (j = 0; j < NXBANDS; j++)
    {
      plasmamain[i].nxtot[j] = 0;
      plasmamain[i].xj[j] = 0.0;
      plasmamain[i].xave_freq[j] = 0.0;
      plasmamain[i].xsd_freq[j] = 0.0;
      plasmamain[i].fmin[j] = geo.xfreq[j + 1]; /* Set the minium frequency to the max frequency in the band */
      plasmamain[i].fmax[j] = geo.xfreq[j];     /* Set the maximum frequency to the min frequency in the band */
    }
    for (j = 0; j < NBINS_IN_CELL_SPEC; ++j)
    {
      plasmamain[i].cell_spec_flux[j] = 0.0;
    }

    for (j = 0; j < nions; j++)
    {
      plasmamain[i].ioniz[j] = 0.0;
      plasmamain[i].recomb[j] = 0.0;
      plasmamain[i].heat_ion[j] = 0.0;
      plasmamain[i].cool_rr_ion[j] = 0.0;
      plasmamain[i].lum_rr_ion[j] = 0.0;
      plasmamain[i].heat_inner_ion[j] = 0.0;

    }
    for (j = 0; j < n_inner_tot; j++)
    {
      plasmamain[i].inner_ioniz[j] = 0.0;
    }
  }
}

#ifdef MPI_ON

/**********************************************************/
/**
 * @brief
 *
 * @details
 *
 * ### Notes ###
 *
 **********************************************************/

static void
communicate_alpha_sp (const int n_start, const int n_stop, const int n_cells_rank)
{
  int i;
  int int_size;
  int double_size;
  int current_rank;

  const int n_cells_max = ceil ((double) NPLASMA / np_mpi_global);
  MPI_Pack_size (1 + n_cells_max, MPI_INT, MPI_COMM_WORLD, &int_size);
  MPI_Pack_size (2 * n_cells_max * size_alpha_est + 2 * n_cells_max * nphot_total, MPI_DOUBLE, MPI_COMM_WORLD, &double_size);
  int comm_buffer_size = double_size + int_size;
  char *comm_buffer = malloc (comm_buffer_size);        // comm_buffer_size is already in bytes

  for (current_rank = 0; current_rank < np_mpi_global; ++current_rank)
  {
    if (rank_global == current_rank)
    {
      int pack_position = 0;
      MPI_Pack (&n_cells_rank, 1, MPI_INT, comm_buffer, comm_buffer_size, &pack_position, MPI_COMM_WORLD);      // how many cells to unpack
      for (i = n_start; i < n_stop; ++i)
      {
        MPI_Pack (&i, 1, MPI_INT, comm_buffer, comm_buffer_size, &pack_position, MPI_COMM_WORLD);       // which cell we're working on
        if (nlevels_macro > 0)
        {
          MPI_Pack (macromain[i].recomb_sp, size_alpha_est, MPI_DOUBLE, comm_buffer, comm_buffer_size, &pack_position, MPI_COMM_WORLD);
          MPI_Pack (macromain[i].recomb_sp_e, size_alpha_est, MPI_DOUBLE, comm_buffer, comm_buffer_size, &pack_position, MPI_COMM_WORLD);
        }
        if (nphot_total > 0)
        {
          MPI_Pack (plasmamain[i].recomb_simple, nphot_total, MPI_DOUBLE, comm_buffer, comm_buffer_size, &pack_position, MPI_COMM_WORLD);
          MPI_Pack (plasmamain[i].recomb_simple_upweight, nphot_total, MPI_DOUBLE, comm_buffer, comm_buffer_size,
                    &pack_position, MPI_COMM_WORLD);
        }
      }
    }

    // TODO: this should be non-blocking eventually
    MPI_Bcast (comm_buffer, comm_buffer_size, MPI_PACKED, current_rank, MPI_COMM_WORLD);

    if (rank_global != current_rank)
    {
      int unpack_position = 0;
      int n_cells_to_do;
      MPI_Unpack (comm_buffer, comm_buffer_size, &unpack_position, &n_cells_to_do, 1, MPI_INT, MPI_COMM_WORLD);
      for (i = 0; i < n_cells_to_do; ++i)
      {
        int current_cell;
        MPI_Unpack (comm_buffer, comm_buffer_size, &unpack_position, &current_cell, 1, MPI_INT, MPI_COMM_WORLD);

        if (nlevels_macro > 0)
        {
          MPI_Unpack (comm_buffer, comm_buffer_size, &unpack_position, macromain[current_cell].recomb_sp, size_alpha_est,
                      MPI_DOUBLE, MPI_COMM_WORLD);
          MPI_Unpack (comm_buffer, comm_buffer_size, &unpack_position, macromain[current_cell].recomb_sp_e,
                      size_alpha_est, MPI_DOUBLE, MPI_COMM_WORLD);
        }
        if (nphot_total > 0)
        {
          MPI_Unpack (comm_buffer, comm_buffer_size, &unpack_position, plasmamain[current_cell].recomb_simple,
                      nphot_total, MPI_DOUBLE, MPI_COMM_WORLD);
          MPI_Unpack (comm_buffer, comm_buffer_size, &unpack_position, plasmamain[current_cell].recomb_simple_upweight,
                      nphot_total, MPI_DOUBLE, MPI_COMM_WORLD);
        }
      }
    }
  }

  free (comm_buffer);
}

#endif

/**********************************************************/
/**
 * @brief
 *
 * @details
 *
 * ### Notes ###
 *
 **********************************************************/

static void
init_macro (void)
{
  int i;
  int k;
  int j;

  for (i = 0; i < NPLASMA; ++i)
  {
    if (geo.rt_mode == RT_MODE_MACRO)
    {
      macromain[i].kpkt_rates_known = FALSE;
    }

    plasmamain[i].kpkt_emiss = 0.0;
    plasmamain[i].kpkt_abs = 0.0;

    for (j = 0; j < nlevels_macro; ++j)
    {
      macromain[i].matom_abs[j] = 0.0;
      macromain[i].matom_emiss[j] = 0.0;

      for (k = 0; k < xconfig[j].n_bbu_jump; ++k)
      {
        macromain[i].jbar[xconfig[j].bbu_indx_first + k] = 0.0;
      }
      for (k = 0; k < xconfig[j].n_bfu_jump; ++k)
      {
        macromain[i].gamma[xconfig[j].bfu_indx_first + k] = 0.0;
        macromain[i].gamma_e[xconfig[j].bfu_indx_first + k] = 0.0;
        macromain[i].alpha_st[xconfig[j].bfd_indx_first + k] = 0.0;
        macromain[i].alpha_st_e[xconfig[j].bfd_indx_first + k] = 0.0;
      }
    }
  }

  int n_start;
  int n_stop;
  int n_cells;

#ifdef MPI_ON
  n_cells = get_parallel_nrange (rank_global, NPLASMA, np_mpi_global, &n_start, &n_stop);
#else
  n_start = 0;
  n_stop = NPLASMA;
  n_cells = NPLASMA;
#endif

  // TODO: should be reworked so we use non-blocking communication: do nlevels_macro first and then nphot_tot

  for (i = n_start; i < n_stop; ++i)
  {
    for (j = 0; j < nlevels_macro; ++j)
    {
      for (k = 0; k < xconfig[j].n_bfd_jump; ++k)
      {
        if (plasmamain[i].t_e > 1.0)
        {
          macromain[i].recomb_sp[xconfig[j].bfd_indx_first + k] = alpha_sp (&phot_top[xconfig[j].bfd_jump[k]], &plasmamain[i], 0);
          macromain[i].recomb_sp_e[xconfig[j].bfd_indx_first + k] = alpha_sp (&phot_top[xconfig[j].bfd_jump[k]], &plasmamain[i], 2);
        }
        else
        {
          macromain[i].recomb_sp[xconfig[j].bfd_indx_first + k] = 0.0;
          macromain[i].recomb_sp_e[xconfig[j].bfd_indx_first + k] = 0.0;
        }
      }
    }
    for (j = 0; j < ntop_phot; ++j)
    {
      if ((geo.macro_simple == FALSE && phot_top[j].macro_info == TRUE) || geo.rt_mode == RT_MODE_2LEVEL)
      {
        plasmamain[i].recomb_simple[j] = 0.0;
        plasmamain[i].recomb_simple_upweight[j] = 1.0;
      }
      else                      // we want a macro approach, but not for this ion so need recomb_simple instead
      {
        const double alpha_store = plasmamain[i].recomb_simple[j] = alpha_sp (&phot_top[j], &plasmamain[i], 2);
        plasmamain[i].recomb_simple_upweight[j] = alpha_sp (&phot_top[j], &plasmamain[i], 1) / alpha_store;
      }
    }
  }

#ifdef MPI_ON
  communicate_alpha_sp (n_start, n_stop, n_cells);
#endif
}

/**********************************************************/
/**
 * @brief      zeros those portions of the wind which contain the radiation properties
 * 	of the wind, i.e those portions which should be set to zeroed when the structure of the
 * 	wind has been changed or when you simply want to start off a calculation in a known state
 *
 * @details
 * The routine is called at the beginning of each ionization calculation
 * cycle.  It should zero all heating and radiation induced cooling in the Plasma structure.  Since
 * cooling is recalculated in wind_update, one needs to be sure that all of the appropriate
 * cooling terms are also rezeroed there as well.
 *
 * ### Notes ###
 *
 **********************************************************/

void
wind_rad_init ()
{
  init_plasma ();
  init_macro ();

  /* XXX Debug code --------------------------------------------------------- */
  /* TODO: let's put this into a unit test */
//  int i, j, k;
//  const int n_start = 0;
//  const int n_stop = NPLASMA;
//
//  for (i = n_start; i < n_stop; ++i)
//  {
//    for (j = 0; j < nlevels_macro; ++j)
//    {
//      Log ("macromain[%d].recomb_sp = [", i);
//      for (k = 0; k < xconfig[j].n_bfd_jump; ++k)
//      {
//        Log (" %g ", macromain[i].recomb_sp[xconfig[j].bfd_indx_first + k]);
//      }
//      Log ("]\nmacromain[%d].recomb_sp_e = [", i);
//      for (k = 0; k < xconfig[j].n_bfd_jump; ++k)
//      {
//        Log (" %g ", macromain[i].recomb_sp_e[xconfig[j].bfd_indx_first + k]);
//      }
//      Log ("]\n");
//    }
//    Log ("plasmamain[%d].recomb_simple = [", i);
//    for (j = 0; j < ntop_phot; ++j)
//    {
//      Log (" %g ", plasmamain[i].recomb_simple[j]);
//    }
//    Log ("]\nplasmamain[%d].recomb_simple_upweight = [", i);
//    for (j = 0; j < ntop_phot; ++j)
//    {
//      Log (" %g ", plasmamain[i].recomb_simple_upweight[j]);
//    }
//    Log ("]\n");
//  }
  /* XXX Debug code --------------------------------------------------------- */
}

/**********************************************************/
/**
 * @brief
 *
 * @details
 *
 * ### Notes ###
 *
 **********************************************************/

void
shell_wind_diagnostics (double xsum, double psum, double fsum, double csum, double icsum, double lsum, double ausum, double chexsum,
                        double cool_sum, double lum_sum)
{


  int ndom;
  int i;
  int n;
  int nn;
  int m;
  int nshell;

  int first;                    //ion
  int last;                     //ion

  double tot;                   //ion?
  double lum_h_line;
  double lum_he_line;
  double lum_c_line;
  double lum_n_line;
  double lum_o_line;
  double lum_fe_line;
  double agn_ip;

  for (ndom = 0; ndom < geo.ndomain; ndom++)
  {
    if (zdom[ndom].wind_type == SHELL)
    {

      /* nshell is the plasma cell that correspond to the second wind cell for the shell_wind model */
      nshell = wmain[zdom[ndom].nstart + 1].nplasma;
      n = plasmamain[nshell].nwind;
      WindPtr w = &wmain[n];
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

      Log
        ("OUTPUT Absorbed_flux(ergs-1cm-3)    %8.2e  (photo %8.2e ff %8.2e compton %8.2e induced_compton %8.2e lines %8.2e auger %8.2e charge_ex %8.2e )\n",
         xsum / w[n].vol, psum / w[n].vol, fsum / w[n].vol, csum / w[n].vol, icsum / w[n].vol, lsum / w[n].vol, ausum / w[n].vol,
         chexsum / w[n].vol);

      /* Report all cooling mechanisms, including those that do not generate photons. */
      Log
        ("OUTPUT Wind_cooling(ergs-1cm-3)     %8.2e (recomb %8.2e ff %8.2e compton %8.2e DR %8.2e DI %8.2e adiabatic %8.2e lines %8.2e ) after update\n",
         cool_sum / w[n].vol, geo.cool_rr / w[n].vol,
         geo.lum_ff / w[n].vol, geo.cool_comp / w[n].vol,
         geo.cool_dr / w[n].vol, geo.cool_di / w[n].vol, geo.cool_adiabatic / w[n].vol, geo.lum_lines / w[n].vol);
      Log ("OUTPUT Wind_luminosity(ergs-1cm-3)     %8.2e (recomb %8.2e ff %8.2e lines %8.2e ) after update\n", lum_sum / w[n].vol,
           geo.lum_rr / w[n].vol, geo.lum_ff / w[n].vol, geo.lum_lines / w[n].vol);
      /* NSH 1701 calculate the recombination cooling for other elements */

      double c_rec = 0.0;
      double n_rec = 0.0;
      double o_rec = 0.0;
      double fe_rec = 0.0;
      double c_lum = 0.0;
      double n_lum = 0.0;
      double o_lum = 0.0;
      double fe_lum = 0.0;
      double c_dr = 0.0;
      double n_dr = 0.0;
      double o_dr = 0.0;
      double fe_dr = 0.0;
      double cool_dr_metals = 0.0;

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
      Log ("radial F_es %i %e \n", nshell, plasmamain[nshell].rad_force_es[0]);
      Log ("radial F_bf %i %e \n", nshell, plasmamain[nshell].rad_force_bf[0]);
      Log ("radial F_ff %i %e \n", nshell, plasmamain[nshell].rad_force_ff[0]);
      Log ("Radial Visible flux %e \n", plasmamain[nshell].F_vis[0]);
      Log ("Radial UV      flux %e \n", plasmamain[nshell].F_UV[0]);
      Log ("Radial Xray    flux %e \n", plasmamain[nshell].F_Xray[0]);
      Log ("Total Radial   flux %e \n", plasmamain[nshell].F_vis[0] + plasmamain[nshell].F_UV[0] + plasmamain[nshell].F_Xray[0]);

    }
  }

}
