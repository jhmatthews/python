
/***********************************************************/
/** @file  resonate.c
 * @author ksl
 * @date   May, 2018
 *
 * @brief
 * These routines have to do with calculating where resonances are and
 * how strong they are.
 *
 * ###Notes###
 * These routines should be made geometry independent and therefore
 * should work with
 * very minor (header only) modifications in other configurations.
 *
 *
 ***********************************************************/



#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "atomic.h"
#include "python.h"



struct photon cds_phot_old;
double cds_v2_old, cds_dvds2_old;

/**********************************************************/
/**
 * @brief     calculate the distance a photon can travel
 * within a single shell without scattering
 *
 * @param [in] WindPtr  w   the entire wind structure
 * @param [in] PhotPtr  p   A photon (bundle)
 * @param [in] double  tau_scat   the optical depth at which the photon
 * will scatter
 * @param [in,out] double *  tau   Initially the current optical depth for
 * the photon; finally the optical depth at the distance the photon can be
 * moved.
 * @param [out] int *  nres   the number of the resonance if there was one, -1 if there was no resonant scatter
 * @param [in] double  smax   the maximum distance the photon can
 * travel in the cell
 * @param [out] int *  istat   A flag indidicating whether the
 * photon should scatter if it traels the distance estimated, 0 if no, TAU_SCAT if yes.
 * @return     The distance the photon can travel
 *
 * calculate_ds finds the distance the photon can travel subject to a
 * number of conditions, which include reach a distance where tau is the
 * distance where a scatter can occur, or a maximum distance set to ensure
 * we do not cross into another cell, or indeed to go so far that one is
 * not sure that the velocity can be approximated as a linear function of
 * distance.
 *
 * The routine returns the distance that the photon can travel subject to
 * these conditions and information about why the photon was halted at the
 * distance it was halted.
 *
 * @details
 *
 * ### Notes ###
 * Calculate_ds does not modify the p in any way!!
 *
 * Any paramaters that depend explicitly on the
 * coordinate gridding, such as the maximum distance the photon
 * can travel before hitting the edge of the
 * shell should be calculated outside of this routine.
 *
 *
 **********************************************************/
double
calculate_ds (w, p, tau_scat, tau, nres, smax, istat)
     WindPtr w;
     PhotPtr p;
     double tau_scat, *tau;
     int *nres;
     double smax;
     int *istat;
{
  int kkk;
  double kap_es;
  double freq_inner, freq_outer, dfreq, ttau, freq_av;
  double mean_freq;             //A mean freq used in compton calculations.
  int n, nn, nstart, ndelt;
  double x;
  double ds_current, ds;
  double v_inner[3], v_outer[3], v1, v2, dvds, dd;
  double v_check[3], vch, vc;
  double dvds1, dvds2;
  struct photon phot, p_now;
  int init_dvds;
  double kap_bf_tot, kap_ff, kap_cont;
  double tau_sobolev;
  WindPtr one, two;
  int check_in_grid;
  int nplasma;
  PlasmaPtr xplasma, xplasma2;
  int ndom;
  double normal[3];

  one = &w[p->grid];            //pointer to the cell where the photon bundle is located.

  nplasma = one->nplasma;
  xplasma = &plasmamain[nplasma];
  ndom = one->ndom;

  ttau = *tau;
  ds_current = 0;
  init_dvds = 0;
  dvds1 = dvds2 = 0.0;
  *nres = -1;
  *istat = P_INWIND;

  if (ttau < 0.0)
  {
    Error ("calculate_ds:  odd tau  %8.2e at %g entering calculate_ds\n", ttau, p->freq);
  }


/* Note that comp_phot compares the position and
 * direction of two photons.  If they are the same, then
 * it just takes v1 from the old value.  */

  if (comp_phot (&cds_phot_old, p))
  {
    vwind_xyz (ndom, p, v_inner);
    v1 = dot (p->lmn, v_inner);
  }
  else
  {
    v1 = cds_v2_old;
  }

  /* Initialize two photon structures phot and p_now for internal work
   * "phot"  will be a photon vector at the far edge of the cell, while p
   * remains the photon at its current positon. p_now  is located
   * at the midpoint between these tow positions.
   */

  stuff_phot (p, &phot);
  move_phot (&phot, smax);
  vwind_xyz (ndom, &phot, v_outer);
  v2 = dot (phot.lmn, v_outer);

  /* Check to see that the velocity is monotonic across the cell
   * by calculating the velocity at the midpoint of the path
   *
   * If it is not monitocinc, then reduce smax
   */
  vc = C;
  while (vc > VCHECK && smax > DFUDGE)
  {
    stuff_phot (p, &p_now);
    move_phot (&p_now, smax / 2.);
    vwind_xyz (ndom, &p_now, v_check);
    vch = dot (p_now.lmn, v_check);

    vc = fabs (vch - 0.5 * (v1 + v2));

    if (vc > VCHECK)
    {
      stuff_phot (&p_now, &phot);
      smax *= 0.5;
      v2 = vch;
    }
  }


/* This Doppler shift shifts the photon from the global to the local
 * frame of rest. Therefore multiply. See doppler notes for a discussion
 * The sign is  correct.  If the photon is moving in the same
 * direction as v, then in the rest frame of the ions,
 * then the photon frequency will be less. */


  freq_inner = p->freq * (1. - v1 / C);
  freq_outer = phot.freq * (1. - v2 / C);
  dfreq = freq_outer - freq_inner;


/* We use the doppler shifted frequency to compute the Klein-Nishina cross
 * section, if the frequency is high enough, otherwise we just use the
 * Thompson cross section.  For the time being, use the average frequency.
 * If we want true fidelity, perhaps we could compute the cross section
 * for every little path section between resonances */

  mean_freq = (freq_inner + freq_outer) / 2.0;

  /*Compute the angle averaged cross section */

  kap_es = klein_nishina (mean_freq) * xplasma->ne * zdom[ndom].fill;


/* The next section checks to see if the frequency difference on
 * the two sides is very small and if not limits the resonances
 * one has to worry about
 */

  if (fabs (dfreq) < EPSILON)
  {
    Error ("calculate_ds: v same at both sides of cell %d\n", one->nwind);
    x = -1;
    return (smax);              // This is not really the best thing to do, but it avoids disaster below

  }
  else if (dfreq > 0)
  {
    limit_lines (freq_inner, freq_outer);
    nstart = nline_min;
    ndelt = 1;
  }
  else
  {
    limit_lines (freq_outer, freq_inner);
    nstart = nline_max;
    ndelt = (-1);
  }

/*nline_min, nline_max, and nline_delt are found in atomic.h and
 * are set by limit_lines()
 */


/* Next part deals with computation of bf opacity. In the macro atom method
 * this is needed.  For the two level approximation it is not.. This section activates
 * if geo.rt_mode==RT_MODE_MACRO (switch for macro atom method). If the
 * macro atom method is not used just get kap_bf to 0 and move on). SS
 */

  kap_bf_tot = 0;
  kap_ff = 0;


  if (geo.rt_mode == RT_MODE_MACRO)
  {
/* Potentially several continuum may contribute in a given frequency
 * range so the kap_bf is an array.
 * Also need to store the total - kap_bf_tot.
 */


    freq_av = freq_inner;       //(freq_inner + freq_outer) * 0.5;  //need to do better than this perhaps but okay for star - comoving frequency (SS)


    kap_bf_tot = kappa_bf (xplasma, freq_av, 0);
    kap_ff = kappa_ff (xplasma, freq_av);

    /* Okay the bound free contribution to the opacity is now sorted out (SS) */
  }

  if (one->vol == 0)
  {
    kap_bf_tot = kap_ff = 0.0;
    Error_silent ("ds_calculate vol = 0: cell %d position %g %g %g\n", p->grid, p->x[0], p->x[1], p->x[2]);
  }


  kap_cont = kap_es + kap_bf_tot + kap_ff;      //total continuum opacity



/* Finally begin the loop over the resonances that can interact
 * with the photon in the cell
 */

  for (n = 0; n < nline_delt; n++)
  {
    nn = nstart + n * ndelt;    /* So if the frequency of resonance increases as we travel through
                                   the grid cell, we go up in the array, otherwise down */
    x = (lin_ptr[nn]->freq - freq_inner) / dfreq;

    if (0. < x && x < 1.)
    {                           /* this particular line is in resonance */
      ds = x * smax;


/* Before checking for a resonant scatter, need to check for scattering
 * due to a continuum
 * process.
 */
      if (ttau + (kap_cont) * (ds - ds_current) > tau_scat)
      {
/* then the photon was scattered by the continuum before reaching the
 * resonance.  Need to randomly select the continumm process which caused
 * the photon to scatter.  The variable threshold is used for this. */

        *nres = select_continuum_scattering_process (kap_cont, kap_es, kap_ff, xplasma);
        *istat = P_SCAT;        //flag as scattering
        ds_current += (tau_scat - ttau) / (kap_cont);   //distance travelled
        ttau = tau_scat;
        *tau = ttau;
        return (ds_current);
      }
      else
      {

/* increment tau by the continuum optical depth to this point */
        ttau += kap_cont * (ds - ds_current);


        ds_current = ds;        /* At this point ds_current is exactly the position of the resonance */
        kkk = lin_ptr[nn]->nion;


/* The density is calculated in the wind array at the center of a cell.
 * We use that as the first estimate of the density.  */
        stuff_phot (p, &p_now);
        move_phot (&p_now, ds_current); // So p_now contains the current position of the photon



        dd = get_ion_density (ndom, p_now.x, kkk);

        if (dd > LDEN_MIN)
        {
/* If we have reached this point then we have to initalize dvds1 and dvds2.
 * Otherwise there is no need to do this, especially as dvwind_ds is an
 * expensive calculation time wise */

          if (init_dvds == 0)
          {
            dvds1 = dvwind_ds (p);
            dvds2 = dvwind_ds (&phot);
            init_dvds = 1;
          }

          dvds = (1. - x) * dvds1 + x * dvds2;



          tau_sobolev = sobolev (one, p->x, dd, lin_ptr[nn], dvds);

/* tau_sobolev now stores the optical depth. This is fed into the next statement for the bb estimator calculation. SS March 2004 */

          ttau += tau_sobolev;


          if (geo.rt_mode == RT_MODE_MACRO)     //Macro Atom case (SS)
          {

/* Because push through distance may take us out of the cell we want,
 * need to make sure that the cell is correct before incrementing the
 * heating rate/estimators. So 1st check if it's still in the wind and
 * second get a pointer to the grid cell where the resonance really happens.
 */

            check_in_grid = walls (&p_now, p, normal);

            if (check_in_grid != P_HIT_STAR && check_in_grid != P_HIT_DISK && check_in_grid != P_ESCAPE)
            {
              /* The next line may be redundant.  */
              two = &w[where_in_grid (wmain[p_now.grid].ndom, p_now.x)];

              if (lin_ptr[nn]->macro_info == 1 && geo.macro_simple == 0)
              {
/* The line is part of a macro atom so increment the estimator if desired */
                if (geo.ioniz_or_extract == 1)
                {
                  bb_estimators_increment (two, p, tau_sobolev, dvds, nn);
                }
              }
              else if (two->vol == 0)
              {
                /* See issue #389 - Sometimes DFUDGE pushes a photon into a cell with no volume.  Note that this
                 * should be very rare, so if this error occurs in significant numbers the problem should be
                 * investigated further.  ksl -180626
                 */
                Error ("calculate_ds: Macro atom problem when photon moved into cell with no volume\n");
              }
              else
              {
/* The line is from a simple ion. Record the heating contribution and move on. */
                xplasma2 = &plasmamain[two->nplasma];

                bb_simple_heat (xplasma2, p, tau_sobolev, dvds, nn);

              }
            }
          }
          /* Completed special calculations for the Macro Atom case */

          /* 68b - 0902 - The next section is to track where absorption
           * is taking place along the line of sight
           * to the observer.  It is probably possibly to simplify
           * some of what is happening here, as we
           * have various photons real and imaginary in this subroutine.
           * p, the orginal photon, phot the
           * photon at the opposide edge of the cell and p_now the photon
           * at its current position.  Some
           * of these could be used to store information needed in phot_hist.
           */

          if (phot_hist_on)
          {
            p_now.tau = ttau;
            p_now.nres = nn;
            phot_hist (&p_now, 1);
          }


        }



        /* Check to see whether the photon should scatter at this point */
        if (ttau > tau_scat)
        {
          *istat = P_SCAT;
          *nres = nn;
          *tau = ttau;




          return (ds_current);
        }

        /* End of loop to process an individual resonance */
      }
      *tau = ttau;
    }
  }



/* If the photon reaches this point it was not scattered by resonances.
 * ds_current is either 0 if there were no resonances or the postion of the
 * "last" resonance if there were resonances.  But we need to check one
 * last time to see if it was scattered by continuum process.
 * Note: ksl -- It is generally a bad policy to have a bunch of repeated code
 * like this.  We should probably rewrite some of this in terms of subroutines
 * for clarity, especially the bit that calculates where the continuum
 * scattering event occurred.  04 apr
 */

  if (ttau + kap_cont * (smax - ds_current) > tau_scat)
  {
    *nres = select_continuum_scattering_process (kap_cont, kap_es, kap_ff, xplasma);

    /* A scattering event has occurred in the shell  and we
     * remain in the same shell */
    ds_current += (tau_scat - ttau) / (kap_cont);
    *istat = P_SCAT;            /* Flag for scattering (SS) */
    ttau = tau_scat;
  }
  else
  {                             /* Then we did hit the other side of the shell
                                   (or possibly the another wall of the same shell) */
    *istat = P_INWIND;
    ttau += kap_cont * (smax - ds_current);     /* kap_es replaced with kap_cont (SS) */
    ds_current = smax;

  }

  *tau = ttau;

  stuff_phot (&phot, &cds_phot_old);    // Store the final photon position
  cds_v2_old = v2;              // and the velocity along the line of sight

  return (ds_current);

}



/**********************************************************/
/**
 * @brief      determine what continuum
 * process was that caused the scatter and returns this to
 * the main routine.
 *
 * @param [in] double  kap_cont   The continuum opacity
 * @param [in] double  kap_es   The electron scattering opacity
 * @param [in] double  kap_ff   The free free opacity
 * @param [in] PlasmaPtr  xplasma   The plasma cell where everyting is being
 * calculated
 * @return     The process that cause the photon to stop/scatter at a particular
 * point
 *
 * @details
 * This routine is called to determine which of several continuum proceseses
 * cause a photon to be scattered or absorpbed.  In addition to electron scattering
 * and free-free absorption, the routine can identify which photoionization process
 * is implicated.
 *
 * ### Notes ###
 *
 * select_continuum_scattering_proces is called when it has already been
 * determined that a photon will stop or be scattered due to some continuum
 * process.  This determination is based on the total continuum opacity.  In this
 * routine a random number is generated and this is used to determine
 * which of the processes was responsible.  Data for the opacity due to
 * photonionization is passed remotely via the PlasmaPtr.
 *
 **********************************************************/

int
select_continuum_scattering_process (kap_cont, kap_es, kap_ff, xplasma)
     double kap_cont, kap_es, kap_ff;
     PlasmaPtr xplasma;
{
  int nres;
  double threshold;
  double run_tot;
  int ncont;

  threshold = random_number (0.0, 1.0) * (kap_cont);


  /* First check for electron scattering. */
  if (kap_es > threshold)
  {                             /* electron scattering event occurred (SS) */
    nres = -1;                  // flag electron scatterin (SS) */
  }
  /* Now chech for ff. */
  else if ((kap_es + kap_ff) > threshold)
  {
    nres = -2;
  }
  /* Now check for bf. */
  else
  {
/* use a running sum to find which photoionisation process it was */
/*
If a non-macro-atom run is being done this part should never be reached.
Just do a check that all is well - this can be removed eventually (SS)
*/
    if (geo.rt_mode == RT_MODE_2LEVEL)
    {
      Error ("calculate_ds: Not using macro atoms but trying to excite one? Aboort.\n");
      exit (0);
    }
    run_tot = kap_es + kap_ff;
    ncont = 0;
    while (run_tot < threshold)
    {
      run_tot += kap_bf[ncont];
      ncont++;
    }
    /* When it gets here know that excitation is in photoionisation labelled by ncont */
    nres = NLINES + 1 + xplasma->kbf_use[ncont - 1];    //modified SS Nov 04
  }
  return (nres);
}



/**********************************************************/
/**
 * @brief      calculates the bf opacity in a specific
 * 	cell.
 *
 * @param [in] PlasmaPtr  xplasma   The plasma cell of interest
 * @param [in] double  freq   The frequency at which the opacity is calculated
 * @param [in] int  macro_all   1--> macro_atoms only, 0 all topbase ions
 * @return     The bf opacity
 *
 * @details
 *
 * ### Notes ###
 * The routine allows for clumping, reducing kappa_bf by the filling
 * factor.
 *
 **********************************************************/

double
kappa_bf (xplasma, freq, macro_all)
     PlasmaPtr xplasma;
     double freq;
     int macro_all;


{
  double kap_bf_tot;
  double ft;
  double x;
  int nconf;
  double density;
  int n;
  int nn;
  int ndom;


  kap_bf_tot = 0;

  macro_all--;                  // Subtract one from macro_all to avoid >= in for loop below.


  ndom = wmain[xplasma->nwind].ndom;

  for (nn = 0; nn < xplasma->kbf_nuse; nn++)    // Loop over photoionisation processes.
  {
    n = xplasma->kbf_use[nn];
    ft = phot_top[n].freq[0];   //This is the edge frequency (SS)

    kap_bf[nn] = 0.0;

    if (freq > ft && freq < phot_top[n].freq[phot_top[n].np - 1] && phot_top[n].macro_info > macro_all)
    {
      /* Need the appropriate density at this point. */

      nconf = phot_top[n].nlev; //Returning lower level = correct (SS)

      density = den_config (xplasma, nconf);    //Need to check what this does (SS)


      if (density > DENSITY_PHOT_MIN || phot_top[n].macro_info == 1)
      {

        /* JM1411 -- added filling factor - density enhancement cancels with zdom[ndom].fill */
        kap_bf[nn] = x = sigma_phot (&phot_top[n], freq) * density * zdom[ndom].fill;   //stimulated recombination? (SS)
        kap_bf_tot += x;
      }
    }
  }


  return (kap_bf_tot);
}


/**********************************************************/
/**
 * @brief      computes and stores the set of bound-free processes which
 *         make significiant contributions to the opacity in a grid cell.
 *
 * @param [in] double  fmin   the lowest frequency of interest to the calculation
 * @param [in] double  fmax   the highest frequency of interest to the calculations
 * @return     Always returns 0
 *
 * The purpose of this routine is to speed up calculations by idenfifying which
 * bound-free x-sections are important enough to be included when calculationg the
 * bound-fee opacity, and which an be ignored because the denisty of the particular
 * ion is so low it will not contribute.
 *
 * For each cell, the routine determines what bf transitons are important
 * and stores them in one->kbf_use[n].
 *
 * The total number of such transitions
 * is given in one->kbf_nuse.
 *
 * @details
 *
 * This routine is now called before various cylces of the Monte Carlo calculation.
 * (in run.c) * It determines which
 * bf processes are worth considering during the calculation that follows.
 *
 * To do this
 * is uses the edge position (must be inside, or close to, the spectral region of
 * interest) and the edge opacity (must be greater than a threshold value, currently
 * 10^-6.
 *
 * The need for the routine is just to prevent wasting time on unimportant bf
 * transitions.  It is called (currently from resonate).
 *
 * ### Notes ###
 * The dimensionalty of kbf_use is currently set to NTOP_PHOT, which is large enough
 * to store evey transition if necessary.
 *
 **********************************************************/

int
kbf_need (fmin, fmax)
     double fmin, fmax;


{
  int nconf;
  double density;
  double tau_test, ft;
  int n;
  int nuse;

  PlasmaPtr xplasma;
  WindPtr one;
  int nplasma, nion;


  for (nplasma = 0; nplasma < NPLASMA; nplasma++)       // Loop over all the cells in the wind
  {
    xplasma = &plasmamain[nplasma];
    one = &wmain[xplasma->nwind];
    nuse = 0;

    for (n = 0; n < nphot_total; n++)   // Loop over photoionisation processes.
    {

      ft = phot_top[n].freq[0]; //This is the edge frequency (SS)

      if ((ft > (fmin / 3.)) && (ft < fmax))
      {
        nion = phot_top[n].nion;

        if (ion[nion].phot_info == 0)   // vfky
        {
          density = xplasma->density[nion];
        }
        // else if (ion[nion].phot_info > 0)  // topbase or hybrid
        else
        {
          nconf = phot_top[n].nlev;     //Returning lower level = correct (SS)
          density = den_config (xplasma, nconf);
        }

        tau_test = phot_top[n].x[0] * density * SMAX_FRAC * length (one->xcen);


        if (tau_test > 1.e-6 || phot_top[n].macro_info == 1)
        {
          /* Store the bf transition and increment nuse */
          xplasma->kbf_use[nuse] = n;
          nuse += 1;
        }
      }

    }
    xplasma->kbf_nuse = nuse;
  }


  return (0);
}



/**********************************************************/
/**
 * @brief      double (p,w,l,dvds)  calculates tau associated with a resonance, given the
 * conditions in the wind and the direction of the photon.
 *
 * It does not modify any of the variables that are passed to it, including for example
 * the photon.
 *
 * @param [in] WindPtr  one   A single wind cell
 * @param [in] double  x[]   A position
 * @param [in] double  den_ion   The density of the ion.  If less than 0, the routine calculates
 * the density at x
 * @param [in] struct lines *  lptr   A pointer to a particular ion
 * @param [in] double  dvds   the velocity gradient in the direction of travel of the photon
 * @return     The optical depth associated with a transition
 *
 * @details
 * In the sobolev approx  tau = PI e^2 / m_e c  * NL * f * lambda * 1 / (|dv/dx|) where dv/ds
 * is the velocity gradient in the direction the photon is travelling.  This can be made slightly
 * simpler by noting that lambda = freq / c and so  tau = PI e^w / m_e
 *
 * ### Notes ###
 *
 * The routine includes a correcton for the filling factor which
 * reduces tau
 *
 **********************************************************/

double
sobolev (one, x, den_ion, lptr, dvds)
     WindPtr one;               // This is a single cell in the wind
     double x[];
     double den_ion;
     struct lines *lptr;
     double dvds;
{
  double tau, xden_ion, tau_x_dvds;
  double two_level_atom (), d1, d2;
  int nion;
  double d_hold;
  int nplasma;
  int ndom;
  PlasmaPtr xplasma;

  nplasma = one->nplasma;
  xplasma = &plasmamain[nplasma];
  ndom = wmain[plasmamain->nwind].ndom;

  if ((dvds = fabs (dvds)) == 0.0)      // This forces dvds to be positive -- a good thing!
  {
    d1 = d2 = 0.;               // Elimiante warning when complied with clang
    tau = VERY_BIG;
    Error ("Sobolev: Surprise tau = VERY_BIG\n");
  }

  else if (lptr->macro_info == 1 && geo.rt_mode == RT_MODE_MACRO && geo.macro_simple == 0)
  {
    // macro atom case SS
    d1 = den_config (xplasma, lptr->nconfigl);
    d2 = den_config (xplasma, lptr->nconfigu);
  }

  else
  {
    nion = lptr->nion;

/* Next few steps to allow used of better calculation of density of this particular
ion which was done above in calculate ds.  It was made necessary by a change in the
calls to two_level atom
*/
    d_hold = xplasma->density[nion];    // Store the density of this ion in the cell

    if (den_ion < 0)
    {
      xplasma->density[nion] = get_ion_density (ndom, x, lptr->nion);   // Forced calculation of density
    }
    else
    {
      xplasma->density[nion] = den_ion; // Put den_ion into the density array
    }
    two_level_atom (lptr, xplasma, &d1, &d2);   // Calculate d1 & d2
    xplasma->density[nion] = d_hold;    // Restore w
  }



  xden_ion = (d1 - lptr->gl / lptr->gu * d2);

  if (xden_ion < 0)
  {
    Error ("sobolev: den_ion has negative density %g %g %g %g %g %g\n", d1, d2, lptr->gl, lptr->gu, lptr->freq, lptr->f);

    /*SS July 08: With macro atoms, the population solver can default to d2 = gu/gl * d1 which should
       give exactly zero here but can be negative, numerically.
       So I'm modyfying this to set tau to zero in such cases, when the populations are vanishingly small anyway. */
    tau_x_dvds = PI_E2_OVER_M * d1 * lptr->f / (lptr->freq);
    tau = tau_x_dvds / dvds;

    tau *= zdom[ndom].fill;     // filling factor is on a domain basis

    if (tau > 1.e-3)
    {
      exit (0);
    }

    else
    {
      Error ("sobolev: continuing by setting tau to zero\n");
      return (0.0);
    }
  }


  tau_x_dvds = PI_E2_OVER_M * xden_ion * lptr->f / (lptr->freq);
  tau = tau_x_dvds / dvds;

  /* JM 1411 -- multiply the optical depth by the filling factor */
  tau *= zdom[ndom].fill;

  return (tau);
}





/**********************************************************/
/**
 * @brief      calculate the  shift given the direction of the incoming
 * 	and outgoing photons and the local  velocity of the wind
 *
 * @param [in] PhotPtr  pin   The pre-scattered  photon (used for its direction)
 * @param [in,out] PhotPtr  pout   the scattered  photon (used for its direction)
 * @param [in] double  v[]   the velocity of the wind where the scatter occurred
 * @param [in] int  nres   either the number of the scatter
 * @return    Always returns 0
 *
 * pout->freq is updated
 *
 * @details
 * Given the incoming and the outgoing photon direction,
 * doppler calculates the outgoing frequency (to first order in beta).
 * There are two basic cases, depending on whether it was a resonant
 * or a nonresonant scatter.
 *
 * ### Notes ###
 * Called from scatter
 *
 **********************************************************/

int
doppler (pin, pout, v, nres)
     PhotPtr pin, pout;
     double v[];
     int nres;

{
  double dot ();

  if (nres == -1)               //Electron scattering (SS)
  {                             /*It was a non-resonant scatter */
    pout->freq = pin->freq * (1 - dot (v, pin->lmn) / C) / (1 - dot (v, pout->lmn) / C);


  }
  else if (nres > -1 && nres < nlines)
  {                             /* It was a resonant scatter. */
    pout->freq = lin_ptr[nres]->freq / (1. - dot (v, pout->lmn) / C);
  }
  else if ((nres > NLINES && nres < NLINES + nphot_total + 1) || nres == -2)
    /* It was continuum emission - new comoving frequency has been chosen by
       the matom/kpkt routine, but now need to convert in the same way
       as for lines (SS) */
  {
    /*
       If a 2-level atom run, one should never arrive here.
       Just do a check that all is well - this can be removed eventually (SS)
     */
    if (geo.rt_mode == RT_MODE_2LEVEL)
    {
      Error ("doppler: Not using macro atoms but trying to deexcite one? Abort.\n");
      exit (0);
    }
    pout->freq = pout->freq / (1. - dot (v, pout->lmn) / C);
  }
/* Now do one final check that nothing is awry.  This is another
 * check added by SS that should probably be deleted or done before this point.
 * I have made it fatal so that we will pay attention to it if it occurs. ksl */

  else
  {
    Error ("doppler: nres %d > NLINES + nphot_total %d\n", nres, NLINES + nphot_total);
    exit (0);
  }

  return (0);

}



/**********************************************************/
/**
 * @brief      determine a new direction and frequency for a photon
 * that is in the wind
 *
 * @param [in,out] PhotPtr  p   the  photon of interest
 * @param [in] int *  nres   either the number of the scatter
 * or a nonresonant scatter if nres < 0
 * @param [out] int *  nnscat   Returned from anisotropic thermal scattering model
 * @return  Always returns 0
 *
 * The results are stored in the PhotPtr p which contains the direction and frequency
 * of the scattered photon.
 *
 * @details
 * The routine calculates a new direction and frequency for a photon in both the
 * resonant and non-resonant cases.  In the frame of the wind, scattering is assumed
 * to be isotropic.
 *
 * ### Notes ###
 * This is the routine that is called when a resonant scatter does occur.  It is
 * relevant for both simple and macro atoms
 *
 * The equations for the frequency shifts are accurate only to first order in beta
 *
 * This routine should not move the photon at all, because other routines need to
 * take this photon in differing directions, and if one moves it here they may
 * encounter this resonance again
 *
 **********************************************************/

int
scatter (p, nres, nnscat)
     PhotPtr p;
     int *nres;
     int *nnscat;
{
  double v[3];
  double z_prime[3];
  int which_out;
  struct photon pold;
  int i, n;
  double p_init[3], p_final[3], dp[3], dp_cyl[3];
  WindPtr one;
  double prob_kpkt, kpkt_choice;
  double gamma_twiddle, gamma_twiddle_e, stim_fact;
  int m, llvl, ulvl;
  double v_dop;
  PlasmaPtr xplasma;
  MacroPtr mplasma;
  int ndom;



  stuff_phot (p, &pold);
  n = where_in_grid (wmain[pold.grid].ndom, pold.x);    // Find out where we are

  if (n < 0)
  {
    Error ("scatter: Trying to scatter a photon in grid cell %d\n", n);
    return (-1);
  }

  one = &wmain[p->grid];
  xplasma = &plasmamain[one->nplasma];

  /* On entering this subroutine we know that a photon packet has been
     absorbed. nres tells us which process absorbed it. There are currently
     four "absorption" processes: electron scattering (flagged -1) line
     absorption (flagged by +ve integer < NLINES), bf absorption
     (flagged by +ve integer > NLINES) and ff absorption (flagged -2). (SS) */

  /* If the macro atom method is being used then the following section must be
     performed to select a macro atom deactivation process. If not then the
     deactivation process is always the same as the activation process and so
     nothing needs to be done. */

  if (geo.rt_mode == RT_MODE_MACRO)     //check if macro atom method in use
  {


    mplasma = &macromain[xplasma->nplasma];

    /* Electron scattering is the simplest to deal with. The co-moving
       frequency is unchanged so it's just a randomisation of the direction.
       For b-b and b-f processes it is first necessary to determine the
       process by which re-emission occurs. (SS). */

    if (*nres > (-1) && *nres < NLINES)
    {
      /* It's a bb line - we can go straight to macro_gov since we know that
         we don't want a k-packet immediately. macro_gov now makes the decision
         regarding the treament (simple or full macro). */
      macro_gov (p, nres, 1, &which_out);
    }

    else if (*nres > NLINES)
    {
      /* This means that it was a photoionisation process.
         For this case we need to decide first whether to excite
         a macro atom directly or to create a k-packet. */

      /*
         The probability if creating a k-packet is given by the
         mc estimators gamma, gamma_e, alpha_st, alpha_st_e.
         Start by identifying which estimators we want and then
         by computing gamma_twiddle (in Leon's notation -
         Lucy 2003 A&A 403 261 */

      /* Now, (Apr04) I'm adding the possibility that the continuum
         is not from a macro ion but from one that we don't have/want
         a macro atom treatment. If it's non macro atom all that happens
         is an on-the-spot decision about whether to excite a fake
         bf macro atom or create a k-packet. Since we've not recorded
         monte carlo estimators for simple ions the decision about creating
         a k-packet is based only on the frequency of the scattered packet.
         (SS) */

      if (phot_top[*nres - NLINES - 1].macro_info == 1 && geo.macro_simple == 0)
      {
        /* Macro ion case (SS) */

        /* Note:  NLINES-1 in the lines below is correct.  This is becasue
           the 1st bf is identified by nres = NLINES+1 and this is
           the zeroth element of phot_top: hence the -1.  SS
         */

        llvl = phot_top[*nres - NLINES - 1].nlev;       //lower level
        ulvl = phot_top[*nres - NLINES - 1].uplev;      //upper level

        for (m = 0; m < config[llvl].n_bfu_jump; m++)
        {
          if (config[llvl].bfu_jump[m] == *nres - NLINES - 1)
          {
            break;
          }
        }

        // m should now be the label to identify which of the bf processes from llvl
        // this is. Check that it is reasonable

        if (m > config[llvl].n_bfu_jump - 1)
        {
          Error ("scatter (resonate.c): could not identify bf transition. Abort. \n");
          exit (0);
        }

        /* Need to compute the factor needed for the stimulated term. */

        stim_fact = den_config (xplasma, ulvl) / den_config (xplasma, llvl) / xplasma->ne;

        gamma_twiddle =
          mplasma->gamma_old[config[llvl].bfu_indx_first + m] - (mplasma->alpha_st_old[config[llvl].bfu_indx_first + m] * stim_fact);
        gamma_twiddle_e =
          mplasma->gamma_e_old[config[llvl].bfu_indx_first + m] - (mplasma->alpha_st_e_old[config[llvl].bfu_indx_first + m] * stim_fact);

        /* Both gamma_twiddles must be greater that zero if this is going to work. If they
           are zero then it's probably because this is the first iteration and so the've not
           been computed yet. For that first iteration k-packets will be ignored. If the
           gamma_twiddles are negative then something has gone wrong.
         */

        if (gamma_twiddle > 0 && gamma_twiddle_e > 0)
        {
          prob_kpkt = 1. - gamma_twiddle / gamma_twiddle_e;
        }
        else if (gamma_twiddle == 0 && gamma_twiddle_e == 0)
        {
          prob_kpkt = 0.;
        }
        else
        {
          Error ("scatter (resonate.c): a gamma_twiddle is negative. Abort.\n");
          exit (0);
        }

        /* Having got here we have calculated the probability of a k-packet
           being created. Now either make a k-packet or excite a macro atom. */

        kpkt_choice = random_number (0.0, 1.0); //random number for kpkt choice

        if (prob_kpkt > kpkt_choice)
        {
          macro_gov (p, nres, 2, &which_out);   //routine to deal with kpkt
        }
        else
        {
          macro_gov (p, nres, 1, &which_out);   //routine to deal with macro atom excitation
        }
      }
      else if (phot_top[*nres - NLINES - 1].macro_info == 0 || geo.macro_simple == 1)
      {
        // Simple ion case //
        /* Need to make decision about making a k-packet. Get the fraction of the energy
           that goes into the electron rather than being stored as ionisation energy: this
           fraction gives the selection probability for the packet. It's given by the
           (photon frequency / edge frequency - 1) (SS) */


        prob_kpkt = 1. - (phot_top[*nres - NLINES - 1].freq[0] / p->freq);


        /* Now choose whether or not to make a k-packet. */

        kpkt_choice = random_number (0.0, 1.0); //random number for kpkt choice


        if (prob_kpkt > kpkt_choice)
        {
          macro_gov (p, nres, 2, &which_out);   //routine to deal with kpkt
        }
        else
        {
          macro_gov (p, nres, 1, &which_out);   //routine to deal with fake macro atom bf excitation

        }
      }
      else
      {
        /* Our best-laid schemes have gang agley. It should never get here unless the input has been
           messed up in some way. (SS) */
        Error ("scatter (resonate.c): continuum scatter - seems to be neither macro nor simple. Abort.\n");
        exit (0);
      }
    }
    else if (*nres == -2)
    {                           /* This is a ff event (SS). */
      macro_gov (p, nres, 2, &which_out);       //ff always make a k-packet
    }
  }

  /* So at this point we have completed all the bits that are specific to the macro approach. 54b--ksl */

  /* SS Apr 04: I've moved this next block up so that nres is set correctly for the call to randvec */

  /* Since the error check is commented out for good reason, we should just assign
   * *nres to p->nres,  and be done with it.  ?? Stuart, assuming you agree just
   * eliminate all the clutter here.  KSL  ??
   */
  p->nres = *nres;              // Update the resonance number on a scatter

  /* SS July 04
     Next block is modified to include the thermal trapping model for anisotropic scattering.
     The code for this has been moved from trans_phot to here so that this model can work
     with macro atoms.
     For macro atoms the code above decides that emission will occur in the line - we now just need
     to use the thermal trapping model to choose the direction. */

  /* JM 1509 -- moved this here so we have ndom for compton direction code */
  ndom = wmain[p->grid].ndom;

  if (*nres == -1)              //Its an electron scatter, so we will call compton to get a direction
  {
    vwind_xyz (ndom, p, v);     //get the local velocity at the location of the photon
    v_dop = dot (p->lmn, v);    //get the dot product of the photon direction with the wind, to get the doppler velocity
    p->freq = p->freq * (1. - v_dop / C);       //This is the photon frequency in the electron rest frame
    compton_dir (p, xplasma);   //Get a new direction using the KN formula
    v_dop = dot (p->lmn, v);    //Find the dot product of the new velocity with the wind
    p->freq = p->freq / (1. - v_dop / C);       //Transform back to the observers frame

  }

  else if (*nres == -2 || *nres > NLINES || geo.scatter_mode == SCATTER_MODE_ISOTROPIC)
  {
    /*  It was either an electron scatter, bf emission or ff emission so the  distribution is isotropic,
       or it was a line photon but we want isotropic scattering anyway.  */
    randvec (z_prime, 1.0);     /* Get a new direction for the photon */
    stuff_v (z_prime, p->lmn);
  }
  else
  {                             //It was a line photon and we want to use the thermal trapping model to choose the output direction

    /* JM 1906 -- added normalisation of the below rejection method. We normalise
       to the escape probability of along the direction of dvds_max, with a safety net of
       20% in case we missed the maximum */
    randwind_thermal_trapping (p, nnscat);
  }

  /* End of modification for thermal trapping model (SS July 04) */




  vwind_xyz (ndom, p, v);       /* Get the velocity vector for the wind */

  if (*nres != -1)              //Only do this if its not an electron scatter, otherwise we have already dealt with this
    doppler (&pold, p, v, *nres);



/* We estimate velocities by interpolating between the velocities at the edges of the cell based
on the photon direction.  We have now changed the direction of the photon, and so we may not
be at the resoance as calculated this way.  reposition moves the photon to the resonance using
the new velocity

Note that one cannot fudge the frequencies, e.g. for some kind of thermal
broadening  before this or else one will defeat the purpose of reposition.

*/



/*Now calculate the momentum transfer.  What follows appears to be
correct only if there was no energy absorbed at the scattering site.
?? The rest of this is only needed in the ionization cycle.  Need to eliminate in the
detailed spectrum calculation ??
*/

  stuff_v (pold.lmn, p_init);
  renorm (p_init, pold.w / C);
  stuff_v (p->lmn, p_final);
  renorm (p_final, p->w / C);
  vsub (p_final, p_init, dp);

  project_from_xyz_cyl (pold.x, dp, dp_cyl);



  if (pold.x[2] < 0)
    dp_cyl[2] *= (-1);
  for (i = 0; i < 3; i++)
  {
    xplasma->dmo_dt[i] += dp_cyl[i];
  }


  return (0);
}
