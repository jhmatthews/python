
/***********************************************************/
/** @file  emission.c
 * @author ksl
 * @date   May, 2018
 *
 * @brief  Subroutines having to do with radiation from the wind,
 * including those that generate photons for the wind
 *
 ***********************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "atomic.h"
#include "python.h"





/**********************************************************/
/**
 * @brief      calculate the luminosity of the entire
 * wind between freqencies f1 and f2
 *
 * @param [in out] double  f1   The minimum frequency for the calculation
 * @param [in out] double  f2   The maximum frequency
 * @return     The luminosity of the entire wind
 *
 * @details
 * The routine simply calls total_emission for each wind cell with
 * a positive volume in the wind
 *
 * The routine also populates several luminosity related variables in
 * geo, which give the luminosity for separate processes, e.g free-free
 * and free_ound emission.
 *
 * ### Notes ###
 * @bug The do loop might be simpler if made over the plasma cells
 * instead of the wind, but one should be careful of the dummy cell
 * 
 * 
 * CK20180801: 
 * 
 *           in non-macro atom mode, the only continuum process treates as scattering is 
 *           electron scattering, and this is assigned nres = -1. The only valid values 
 *           of nres in non-macro-atom mode are therefore nres = -1 and 0 <= nres <= nlines-1
 *           (with the lattter range covering the lines).
 * 
 *           in macro atom mode, nres = -1 indicates electron scattering, 
 *           nres = -2 indicates ff, and nres > NLINES indicates bound-free. 
 * 	     [nres == NLINES is never used. Note also that NLINES is the *max* number of lines, whereas nlines
 *	     is the *actual* number of lines. So, actually, it's not just nres = NLINES that's never used, but 
 *	     the entire range of nlines <= nres <= NLINES]
 * 
 **********************************************************/

double
wind_luminosity (f1, f2)
     double f1, f2;             /* freqmin and freqmax */
{
  double lum, lum_lines, lum_rr, lum_ff;
  int n;
  double x;
  int nplasma;


  lum = lum_lines = lum_rr = lum_ff = 0.0;
  for (n = 0; n < NDIM2; n++)
  {

    if (wmain[n].vol > 0.0)
    {
      nplasma = wmain[n].nplasma;
      lum += x = total_emission (&wmain[n], f1, f2);
      lum_lines += plasmamain[nplasma].lum_lines;
      lum_rr += plasmamain[nplasma].lum_rr;
      lum_ff += plasmamain[nplasma].lum_ff;
    }
  }


  geo.lum_lines = lum_lines;
  geo.lum_rr = lum_rr;
  geo.lum_ff = lum_ff;

  return (lum);
}





/**********************************************************/
/**
 * @brief      Calculate the band-limited emission of a single cell
 *
 * @param [in] WindPtr  one   The wind cell of interest
 * @param [in] double  f1   The minimum frequency for the calculation
 * @param [in] double  f2   The maximum frequency for the calculation
 * @return
 * It returns the total luminosity (within frequency limits) 
 *
 * The routine also stores 
 * the luminosity due
 * to various emission processes, e.g ff, fb, lines, compton into 
 * varius variables in thea associated Plasma cell
 *
 * @details
 *
 * ### Notes ###
 * Total emission gives the total enery loss due to photons.  It does
 * not include other cooling sources, e. g. adiabatic expansion.
 *
 * The name total emission is a misnomer.  The returns a
 * band limited luminosity.  This is because the routine can be used
 * to establish the number of photons to be emitted by the wind. 
 *
 * Comment:  Compton cooling is not included here because this
 * does not result in new photons.
 *
 *
 **********************************************************/

double
total_emission (one, f1, f2)
     WindPtr one;               /* WindPtr to a specific cell in the wind */
     double f1, f2;             /* The minimum and maximum frequency over which the emission is
                                   integrated */
{
  double t_e;
  int nplasma;
  PlasmaPtr xplasma;

  nplasma = one->nplasma;
  xplasma = &plasmamain[nplasma];

  t_e = xplasma->t_e;


  if (f2 < f1)
  {
    xplasma->lum_tot = xplasma->lum_lines = xplasma->lum_ff = xplasma->lum_rr = 0;
  }
  else
  {
    if (geo.rt_mode == RT_MODE_MACRO)   //Switch for macro atoms (SS)
    {
      xplasma->lum_rr = total_fb_matoms (xplasma, t_e, f1, f2) + total_fb (one, t_e, f1, f2, FB_FULL, OUTER_SHELL);     //outer shellrecombinations

      /*
       *The first term here is the fb cooling due to macro ions and the second gives
       *the fb cooling due to simple ions.
       *total_fb has been modified to exclude recombinations treated using macro atoms.
       */
      xplasma->lum_tot = xplasma->cool_rr;
      /* Note: This the fb_matom call makes no use of f1 or f2. They are passed for
       * now in case they should be used in the future. But they could
       * also be removed.
       * (SS)
       */
      xplasma->lum_lines = total_bb_cooling (xplasma, t_e);
      xplasma->lum_tot += xplasma->lum_lines;
      /* total_bb_cooling gives the total cooling rate due to bb transisions whether they
         are macro atoms or simple ions. */
      xplasma->lum_ff = total_free (one, t_e, f1, f2);
      xplasma->lum_tot += xplasma->lum_ff;


    }
    else                        //default (non-macro atoms) (SS)
    {
      xplasma->lum_tot = xplasma->lum_lines = total_line_emission (one, f1, f2);
      xplasma->lum_tot += xplasma->lum_ff = total_free (one, t_e, f1, f2);
      /* We compute the radiative recombination luminosirty - this is not the same as the rr cooling rate and
         so is stored in a seperate variable */
      xplasma->lum_tot += xplasma->lum_rr = total_fb (one, t_e, f1, f2, FB_FULL, OUTER_SHELL);  //outer shell recombinations



    }
  }


  return (xplasma->lum_tot);


}




/**********************************************************/
/**
 * @brief      generates
 * 	photons within the wind between two freqencies and stores them in the photon array
 *
 * @param [out] PhotPtr  p   The entire photon stucture
 * @param [in] double  weight   weight of each photon to be generated
 * @param [in] double  freqmin   The minimum frequency
 * @param [in] double  freqmax   The maximum frequency
 * @param [in] int  photstart   The place in the photon stucture where the first photon will be stored
 * @param [in] int  nphot   The number of photon bundles to generate
 * @return     Normally returns the number of photons created.
 *
 * When the routine is completed new  photons will exist from p[istart]
 * to p[istart+return value].
 *
 * @details
 *
 * The routine first generates a random number which is used to determine
 * in which wind cell  should be generated, and then  determines the
 * type of photon to generate.  Once this is doen the routine cycles thourhg
 * the PlasmaCells generatating all of the photons for each cell at once.
 *
 *
 * ### Notes ###
 * This logic was adopted for speed related reasons.
 *
 * If photo_gen_wind tries to create more photons than exist in the photon structure the
 * program will stop (rather than continue incorrectly or start blasting away memory.)
 *
 * @bug This is another example where the iteration might be made over plasma cells
 * instead of looking for wind cells with positive volume
 **********************************************************/

int
photo_gen_wind (p, weight, freqmin, freqmax, photstart, nphot)
     PhotPtr p;
     double weight;
     double freqmin, freqmax;
     int photstart, nphot;
{
  int n, nn, np;
  int photstop;
  double xlum, xlumsum, lum;
  double v[3];
  int icell, icell_old;
  int nplasma = 0;
  int nnscat;
  int ndom;
  int ptype[NPLASMA][3];        //Store for the types of photons we want, ff first, fb next, line third

  for (n = 0; n < NPLASMA; n++)
  {
    for (nn = 0; nn < 3; nn++)
      ptype[n][nn] = 0;
  }

  /* Limit the lines to consider */
  limit_lines (freqmin, freqmax);


  photstop = photstart + nphot;
  Log_silent ("photo_gen_wind creates nphot %5d photons from %5d to %5d \n", nphot, photstart, photstop);

  for (n = photstart; n < photstop; n++)
  {
    /* locate the wind_cell in which the photon bundle originates.
       Note: In photo_gen, both geo.f_wind and geo.lum_wind will have been determined.
       geo.f_wind refers to the specific flux between freqmin and freqmax.  Note that
       we make sure that xlum is not == 0 or to geo.f_wind. */

    xlum = random_number (0.0, 1.0) * geo.f_wind;


    xlumsum = 0;
    icell = 0;
    while (xlumsum < xlum)
    {


      if (wmain[icell].vol > 0.0)
      {
        nplasma = wmain[icell].nplasma;
        /*increment the xlumsum by the lum_tot (the band limited luminosity in this cell */
        xlumsum += plasmamain[nplasma].lum_tot;
      }
      icell++;
    }
    icell--;

    /* At this point we know the cell in which the photon will be generated */

    nplasma = wmain[icell].nplasma;
    ndom = wmain[icell].ndom;
    plasmamain[nplasma].nrad += 1;      /* Increment the counter for the number of photons generatd in the cell */



    /*Determine the type of photon this photon will be and increment ptype, which stores the total number of
     * each photon type to be made in each cell */

    lum = plasmamain[nplasma].lum_tot;
    xlum = lum * random_number (0.0, 1.0);
    xlumsum = 0;

    p[n].nres = -1;
    p[n].nnscat = 1;
    if ((xlumsum += plasmamain[nplasma].lum_ff) > xlum)
    {
      ptype[nplasma][0]++;      /* a ff photon  */
    }
    else if ((xlumsum += plasmamain[nplasma].lum_rr) > xlum)
    {
      ptype[nplasma][1]++;      /* a fb photon */
    }
    else
    {
      ptype[nplasma][2]++;      /* a line photon */
    }
  }


/* Now actually generate the photons looping over the Plasma cells */

  photstop = photstart;

  icell_old = (-1);
  for (n = 0; n < NPLASMA; n++)
  {

    photstart = photstop;       //initially set to photstart, afterwards we start the photon number from the end of the last cell
    photstop = photstart + ptype[n][0] + ptype[n][1] + ptype[n][2];     //This is the number of photons in this cell

    icell = plasmamain[n].nwind;
    ndom = wmain[icell].ndom;

    for (np = photstart; np < photstop; np++)
    {

      if (np < photstart + ptype[n][0])
      {
        p[np].freq = one_ff (&wmain[icell], freqmin, freqmax);  /*Get the frequency of one ff photon */
        if (p[np].freq <= 0.0)
        {
          Error_silent
            ("photo_gen_wind: On return from one_ff: icell %d vol %g t_e %g\n", icell, wmain[icell].vol, plasmamain[nplasma].t_e);
          p[np].freq = 0.0;
        }
      }
      else if (np < photstart + ptype[n][0] + ptype[n][1])
      {
        p[np].freq = one_fb (&wmain[icell], freqmin, freqmax);
      }
      else
      {

        if (icell != icell_old)
        {
          lum_lines (&wmain[icell], nline_min, nline_max);      /* fill the lin_ptr->pow array. This must be done because it is not stored
                                                                   for all cells.  The if statement is intended to prevent recalculating the power if more than
                                                                   one line photon is generated from this cell in this cycle. */
          icell_old = icell;
        }
        p[np].freq = one_line (&wmain[icell], &p[np].nres);     /*And fill all the rest of the luminosity up with line photons */
        if (p[np].freq == 0)
        {
          Error ("photo_gen_wind: one_line returned 0 for freq %g %g\n", freqmin, freqmax);
        }
      }

      p[np].w = weight;
      get_random_location (icell, p[np].x);
      p[np].grid = icell;
      /*
         Determine the direction of the photon
         Need to allow for anisotropic emission here

         CK20180801: in non-macro-atom mode, the only continuum process that is treated
         as a scatter is electron scattering; but in any case, all continuum processes in non-macro-atom
         mode correspond to nres < 0.

         in macro-atom mode, bf and ff continuum processes are also treated as scattering processes, 
         and then they, too, are treated isotropically. however, in macro atom mode, nres = -1 specifically
         indicates electron scattering, nres = -2 specifically indicates ff, and nres > NLINES indicates
         bound-free. 
         [nres == NLINES is never used. Note also that NLINES is the *max* number of lines, whereas nlines
         is the *actual* number of lines. So, actually, it's not just nres = NLINES that's never used, but 
         the entire range of nlines <= nres <= NLINES]

         in both macro-atom and non-macro-atom modes, 
         line processes are treated isotropically only if geo.scatter_mode == SCATTER_MODE_ISOTROPIC;
         note that, confusingly, geo.scatter_mode == SCATTER_MODE_ANISOTROPIC is *not*
         the only anisotropic mode -- SCATTER_MODE_THERMAL is *also* anisotropic.

       */
      nnscat = 1;

      if (p[np].nres < 0 || geo.scatter_mode == SCATTER_MODE_ISOTROPIC)
      {
/*  It was either an electron scatter so the  distribution is isotropic, or it
was a resonant scatter but we want isotropic scattering anyway.  */
        randvec (p[np].lmn, 1.0);       /* The photon is emitted isotropically */
      }
      else if (geo.scatter_mode == SCATTER_MODE_THERMAL)
      {                         // It was a line photon and we want anisotropic scattering
        randwind_thermal_trapping (&p[np], &nnscat);
      }
      p[np].nnscat = nnscat;
      /* The next two lines correct the frequency to first order, but do not result in
         forward scattering of the distribution */
      vwind_xyz (ndom, &p[np], v);
      p[np].freq *= (1. + dot (v, p[np].lmn) / C);
      p[np].istat = 0;
      p[np].tau = p[np].nscat = p[np].nrscat = 0;
      p[np].origin = PTYPE_WIND;        // A wind photon
      switch (geo.reverb)
      {                         // SWM 26-3-15: Added wind paths
      case REV_WIND:
      case REV_MATOM:
        wind_paths_gen_phot (&wmain[icell], &p[np]);
        break;
      case REV_PHOTON:
        simple_paths_gen_phot (&p[np]);
        break;
      case REV_NONE:
      default:
        break;
      }
    }

  }


  return (nphot);               /* Return the number of photons generated */
}




/**********************************************************/
/**
 * @brief      gets the frequency of a
 * single collisionally excited line photon in a particular cell
 * of the wind.
 *
 * @param [in] WindPtr  one   The wind cell
 * @param [out] int *  nres   The number associated with the transition
 * @return     The frequency of the transition
 *
 * @details
 * The routine geneerates a random number and uses this to select
 * the transition that was excited
 *
 * ### Notes ###
 *
 **********************************************************/

double
one_line (one, nres)
     WindPtr one;
     int *nres;
{
  double xlum, xlumsum;
  int m;
  int nplasma;
  PlasmaPtr xplasma;
  nplasma = one->nplasma;
  xplasma = &plasmamain[nplasma];
  /* Put in a bunch of checks */
  if (xplasma->lum_lines <= 0)
  {
    Error ("one_line: requesting a line when line lum is 0\n");
    return (0);
  }
  if (nline_min == nline_max)
  {
    Error ("one_line: no lines %d %d\n", nline_min, nline_max);
    return (0);
  }

  xlum = xplasma->lum_lines * random_number (0.0, 1.0);

  xlumsum = 0;
  m = nline_min;
  while (xlumsum < xlum && m < nline_max)
  {
    xlumsum += lin_ptr[m]->pow;
    m++;
  }
  m--;
  *nres = m;
  return (lin_ptr[m]->freq);
}





/**********************************************************/
/** Next section deals with bremsstrahlung radiation
 * 4*PI* 8/3* sqrt(2*PI/3)*e**6/m**2/c**3 sqrt(m/kT) or
 * 4*PI times the normal constant to dL_nu/dnu
 **********************************************************/
#define BREMS_CONSTANT 6.85e-38


/**********************************************************/
/**
 * @brief      calculates the band-limited ff luminosity of a cell.
 *
 * @param [in] WindPtr  one   A wind cell
 * @param [in] double  t_e   The electron temperature to use
 * @param [in] double  f1   The minimum frequency
 * @param [in] double  f2   The maximum frequency
 * @return     The free-free luminosity between f1 and f2
 *
 * @details
 *
 * ### Notes ###
 * The program uses an integral formula rather than integrating on
 * the fly which is not ideal but saves time
 *
 * The gaunt factor used depends on whether data for this is read in
 * from the atomic data files.  If it is, then the gaunt factor determined
 * using data from Sutherland (1998).  If not, the gaunt factor is set
 * to 1.
 *
 **********************************************************/

double
total_free (one, t_e, f1, f2)
     WindPtr one;
     double t_e;
     double f1, f2;
{
  double g_ff_h, g_ff_he;
  double gaunt;
  double x, sum;
  double gsqrd;                 /*The scaled inverse temperature experienced by an ion - used to compute the gaunt factor */
  int nplasma, nion;
  PlasmaPtr xplasma;
  nplasma = one->nplasma;
  xplasma = &plasmamain[nplasma];
  if (f2 < f1)
  {
    Error ("total_free: band limited ff  emissivity requested by f1 %g > f2 %g\n", f1, f2);
    return (0.0);
  }

  if (ALPHA_FF * xplasma->t_e / H_OVER_K < f1)
  {
    return (0.0);
  }

  if (ALPHA_FF * xplasma->t_e / H_OVER_K < f2)
  {
    f2 = ALPHA_FF * xplasma->t_e / H_OVER_K;
  }

  if (t_e < TMIN)
  {
    return (0.0);
  }


  if (gaunt_n_gsqrd == 0)       //Maintain old behaviour because atomic data files do not include gaunt factor.
  {
    g_ff_h = g_ff_he = 1.0;
    if (nelements > 1)
    {
      x = BREMS_CONSTANT * xplasma->ne * (xplasma->density[1] * g_ff_h + 4. * xplasma->density[4] * g_ff_he) / H_OVER_K;
    }
    else
    {
      x = BREMS_CONSTANT * xplasma->ne * (xplasma->density[1] * g_ff_h) / H_OVER_K;
    }
  }
  else
  {
    sum = 0.0;
    for (nion = 0; nion < nions; nion++)
    {
      if (ion[nion].istate != 1)        //The neutral ion does not contribute
      {
        gsqrd = ((ion[nion].istate - 1) * (ion[nion].istate - 1) * RYD2ERGS) / (BOLTZMANN * t_e);
        gaunt = gaunt_ff (gsqrd);
        sum += xplasma->density[nion] * (ion[nion].istate - 1) * (ion[nion].istate - 1) * gaunt;
      }
      else
      {
        sum += 0.0;
      }
    }
    x = BREMS_CONSTANT * xplasma->ne * (sum) / H_OVER_K;
  }

  /* JM 1604 -- The reason why this is proportional to t_e**1/2,
     rather than t_e**(-1/2) as in equation 40 of LK02 is because
     one gets an extra factor of (k*t_e/h) when one does the integral */
  x *= sqrt (t_e) * xplasma->vol;
  x *= (exp (-H_OVER_K * f1 / t_e) - exp (-H_OVER_K * f2 / t_e));
  return (x);
}





/**********************************************************/
/**
 * @brief      calculate f_nu for free free emisssion
 *
 * @param [in] WindPtr  one   A wind cell
 * @param [in] double  t_e   The temperature of the plasma
 * @param [in] double  freq   The frequency at which f_nu is to be calculated
 * @return     f_nu for the specific freqency, temperature requested and densities
 * contained in the cell indicated
 *
 * @details
 *
 * ### Notes ###
 * @bug f_nu is set to 0 for t_e less than 100K.  It's not clear
 * that this limit is applied to other functions anymore.  Was this
 * missed?  Not also that the code having to do with the gaunt factor
 * is duplicated from another routine.  Should a gaunt_ff routine
 * be created for both?
 *
 * Within python, this routine is accessed through one_ff
 *
 **********************************************************/

double
ff (one, t_e, freq)
     WindPtr one;
     double t_e, freq;
{
  double g_ff_h, g_ff_he;
  double fnu;
  double gsqrd, gaunt, sum;
  int nplasma;
  int nion;
  PlasmaPtr xplasma;
  nplasma = one->nplasma;
  xplasma = &plasmamain[nplasma];
  if (t_e < TMIN)
    return (0.0);
  if (gaunt_n_gsqrd == 0)       //Maintain old behaviour because gaunt factors have not been provided in atomic data files
  {
    g_ff_h = g_ff_he = 1.0;
    if (nelements > 1)
    {
      fnu = BREMS_CONSTANT * xplasma->ne * (xplasma->density[1] * g_ff_h + 4. * xplasma->density[4] * g_ff_he);
    }
    else
    {
      fnu = BREMS_CONSTANT * xplasma->ne * (xplasma->density[1] * g_ff_h);
    }
  }
  else
  {
    sum = 0.0;
    for (nion = 0; nion < nions; nion++)
    {
      if (ion[nion].istate != 1)        //The neutral ion does not contribute
      {
        gsqrd = ((ion[nion].istate - 1) * (ion[nion].istate - 1) * RYD2ERGS) / (BOLTZMANN * t_e);
        gaunt = gaunt_ff (gsqrd);
        sum += xplasma->density[nion] * (ion[nion].istate - 1) * (ion[nion].istate - 1) * gaunt;
      }
      else
      {
        sum += 0.0;
      }
    }
    fnu = BREMS_CONSTANT * xplasma->ne * (sum);
  }


  fnu *= exp (-H_OVER_K * freq / t_e) / sqrt (t_e) * xplasma->vol;
  return (fnu);
}



/// We initialise the arrays that will contain the unscaled PDF
double ff_x[ARRAY_PDF], ff_y[ARRAY_PDF];

/// Old values
double one_ff_f1, one_ff_f2, one_ff_te;


/**********************************************************/
/**
 * @brief      randomly generate the frequency of a
 * 	ff photon within the frequency interval f1 and f2
 *
 * @param [in] WindPtr  one   A specific wind cell
 * @param [in] double  f1   The minimum frequency
 * @param [in] double  f2   The maximum frequency
 * @return     A randomly geneted frequncy for a ff photon given
 * the conditions specified
 *
 * @details
 *
 * ### Notes ###
 * one is the windcell where the photon will be created.  It is needed
 * only for the temperature.  The code would be simplified
 * if simply the temperature were transmitted.
 *
 * The routine creeates a cdf given the conditions specified and then
 * samples this.  When the routine is re-entered it checks if the
 * conditions are the same as previously.  Given this it is desirable
 * to generate all of the photons in a single cell in one go.
 *
 **********************************************************/

double
one_ff (one, f1, f2)
     WindPtr one;               /* a single cell */
     double f1, f2;             /* freqmin and freqmax */
{
  double freq, dfreq;
  int n;
  int nplasma;
  PlasmaPtr xplasma;
  int echeck;
  nplasma = one->nplasma;
  xplasma = &plasmamain[nplasma];
  if (f2 < f1)
  {
    Error ("one_ff: Bad inputs f2 %g < f1 %g returning 0.0  t_e %g\n", f2, f1, xplasma->t_e);
    return (-1.0);
  }

  /* Check to see if we have already generated a pdf */

  if (xplasma->t_e != one_ff_te || f1 != one_ff_f1 || f2 != one_ff_f2)
  {                             /* Generate a new pdf */
    dfreq = (f2 - f1) / (ARRAY_PDF - 1);
    for (n = 0; n < ARRAY_PDF - 1; n++)
    {
      ff_x[n] = f1 + dfreq * n;
      ff_y[n] = ff (one, xplasma->t_e, ff_x[n]);
    }

    ff_x[ARRAY_PDF - 1] = f2;
    ff_y[ARRAY_PDF - 1] = ff (one, xplasma->t_e, ff_x[ARRAY_PDF - 1]);
    if ((echeck = cdf_gen_from_array (&cdf_ff, ff_x, ff_y, ARRAY_PDF, f1, f2)) != 0)
    {
      Error
        ("one_ff: cdf_gen_from_array error %d : f1 %g f2 %g te %g ne %g nh %g vol %g\n",
         echeck, f1, f2, xplasma->t_e, xplasma->ne, xplasma->density[1], one->vol);
      Exit (0);
    }
    one_ff_te = xplasma->t_e;
    one_ff_f1 = f1;
    one_ff_f2 = f2;             /* Note that this may not be the best way to check for a previous pdf */
  }
  freq = cdf_get_rand (&cdf_ff);
  return (freq);
}




/**********************************************************/
/**
 * @brief      computes the frequency averaged gaunt factor for ff emissionat
 * 		scaled temperature from Sutherland (1988).
 *
 * @param [in] double  gsquared   The variable on which
 * Sutherland's ff gaunt factors are interpolated.
 * @return     An estimate of the Gaunt factor (for a particular ion)
 *
 * @details
 * It interpolates simply interpolates between the tabulated
 * factors from Table 3 of Sutherland (which is a spline fit
 * to the calculated gaunt factors at various values of gsquared)
 *
 * ### Notes ###
 * The reference for this is Sutherland, R.~S. 1998, MNRAS, 300, 321
 *
 * gsquared is the scaled inverse temperature experienced by an ion,
 * Z**2/kT(Ry).
 *
 **********************************************************/

double
gaunt_ff (gsquared)
     double gsquared;           /* the gamma squared variable */
{
  int i, index;
  double gaunt;
  double log_g2;
  double delta;                 //The log difference between our G2 and the one in the table

  delta = 0.0;                  /* NSH 130605 to remove o3 compile error */
  index = 0;                    /* NSH 130605 to remove o3 compile error */
  log_g2 = log10 (gsquared);    //The data is in log format
  if (log_g2 < gaunt_total[0].log_gsqrd || log_g2 > gaunt_total[gaunt_n_gsqrd - 1].log_gsqrd)
  {
    return (1.0);
  }

//OLD  for (i = 0; i < gaunt_n_gsqrd; i++)   /*first find the pair of parameter arrays that bracket our temperature */
//OLD  {
//OLD    if (gaunt_total[i].log_gsqrd <= log_g2 && gaunt_total[i + 1].log_gsqrd > log_g2)
//OLD    {
//OLD      index = i;                /* the array to use */
//OLD      delta = log_g2 - gaunt_total[index].log_gsqrd;
//OLD    }
//OLD  }

  i = 0;
  while (gaunt_total[i].log_gsqrd < log_g2)
  {
    i++;
  }

  index = i - 1;
  delta = log_g2 - gaunt_total[index].log_gsqrd;

  /* The outherland interpolation data is a spline fit to the gaunt function. */

  gaunt = gaunt_total[index].gff + delta * (gaunt_total[index].s1 + delta * (gaunt_total[index].s2 + gaunt_total[index].s3));
  return (gaunt);
}
