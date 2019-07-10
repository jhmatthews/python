
/***********************************************************/
/** @file  extract.c
 * @author ksl
 * @date   April, 2018
 *
 * @brief  Maing routines for extracting photons during the 
 * spectral generation phase.
 *
 ***********************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "atomic.h"
#include "python.h"



/**********************************************************/
/** 
 * @brief      A supervisory routine called to 
 * 	builds detailed spectra in the normal (extract) mode.
 *
 * @param [in, out] WindPtr  w   The entire wind
 * @param [in, out] PhotPtr  p   The photon to extract
 * @param [in, out] int  itype   An integer representing the type of photon 
 * for the purpose of being extracted.
 * @return     Always returns 0
 *
 * @details
 * extract is called when a photon begins it's flight and every time that photon
 * scatters, unless the user has exercised the "live or die" option, in
 * which case it is not called.  
 *
 * extract carries out several preperatory steps for the extraction
 * and for each spectrum one wants to build, it calls extract_one,
 * where the actual incrementing of the spectrum is done.
 *
 * itype takes on the following values:
 * * PTYPE_STAR->the photon came for the star 
 * * PTYPE_BL->the photon came from the boundary layer
 * * PTYPE_DISK->the photon being redirected arose in the disk,
 * * PTYPE_WIND->the photon being redirected arose in the wind,
 * 
 * extract uses the types to prepare the photon for extraction, including
 * doppler shifting the photon if is of PTYPE_WIND or PTYPE_STAR 
 * 
 * Usually, Python constructs a spectrum of all photons, but there are
 * advanced options which allone to restict the spectrum created to
 * those produced with a certain number of scatters or from photons 
 * that arise from the above or below the disk.  extract enforces
 * those choices before calling extract_one.
 * The parameters for this option all come in through
 * python.h, and are contained in the spectrum structure. 
 *
 * The basic choices, some of which can be used in tandom are as follows:
 *
 * *	If s[n].nscat>999; then one obtains the entire spectrum
 * *	if s[n].nscat is a positive number, then the spectrum is composed just 
 * 		of photons with that number of scatters
 * *	if s[n].nscat is a negative number, then the spectrum contains all photons
 * 		with >=  the absolute value of s[n].nscat
 * * 	s[n].top_bot=0 -> accept all photons
 * * 	s[n].top_bot    >0  -> accept photons last above the disk
 * * 	s[n].top_bot    <0  -> accept photons below the disk
 * *    One can also select photons that are to be extracted from a particular
 * spherical region.
 * ### Notes ###
 * 
 * @bug This routine as well as extract_one have options for tracking the photon
 * history.  The routines are in diag.c It is not clear that they have been used 
 * in a long time and so it may be worthwhile to remove them. Furthermore,
 * we have established a new mechanism save_phot for essentially this same
 * task.  This really should be consolidated. 
 *
 * @bug This is also commented in the text, but there is a rather bizarre separation
 * for where the photon frequency is updated and where the weight is update. The former is
 * done in extract, while the latter is done in extract_one.  This is not an error precisely
 * but makes the code more confusing than it needs to be.
 *
 **********************************************************/

int
extract (w, p, itype)
     WindPtr w;
     PhotPtr p;
     int itype;
{
  int n, mscat, mtopbot;
  struct photon pp;
  double v[3];
  double length ();
  int vsub ();
  int yep;
  double xdiff[3];
  int ndom;


  /* The next line selects the middle inclination angle for recording the absorbed energy */
  phot_history_spectrum = 0.5 * (MSPEC + nspectra);

  for (n = MSPEC; n < nspectra; n++)
  {
    /* If statement allows one to choose whether to construct the spectrum
       from all photons or just from photons that have scattered a specific number
       of times or in specific regions of the wind. A region is specified by a position
       and a radius. */

    yep = 1;                    // Start by assuming it is a good photon for extraction

    if ((mscat = xxspec[n].nscat) > 999 || p->nscat == mscat || (mscat < 0 && p->nscat >= (-mscat)))
      yep = 1;
    else
      yep = 0;

    if (yep)
    {
      if ((mtopbot = xxspec[n].top_bot) == 0)
        yep = 1;                // Then there are no positional parameters and we are done
      else if (mtopbot == -1 && p->x[2] < 0)
        yep = 1;
      else if (mtopbot == 1 && p->x[2] > 0)
        yep = 1;
      else if (mtopbot == 2)    // Then to count, the photom must originate within sn.r of sn.x
      {
        vsub (p->x, xxspec[n].x, xdiff);
        if (length (xdiff) > xxspec[n].r)
          yep = 0;

      }
      else
        yep = 0;
    }

    if (yep)                    //Then we want to extract this photon
    {


/* Create a photon pp to use here and in extract_one.  This assures we
 * have not modified p as part of extract
 */

      stuff_phot (p, &pp);
      stuff_v (xxspec[n].lmn, pp.lmn);  /* Stuff new photon direction into pp */

/* 

Need to frequency shift the disk photons as well as the wind 
photons.    

Note that split of functionality between this and extract 
one is odd. We do frequency here but weighting is carried out in  extract */

      if (itype == PTYPE_DISK)
      {
        vdisk (pp.x, v);
        doppler (p, &pp, v, -1);

      }
      if (itype == PTYPE_WIND)
      {                         /* If the photon was scattered in the wind, 
                                   the frequency also must be shifted */
        ndom = wmain[p->grid].ndom;
        vwind_xyz (ndom, &pp, v);       /*  Get the velocity at the position of pp */
        doppler (p, &pp, v, pp.nres);   /*  Doppler shift the photon -- test! */

/*  Doppler shift the photon (as nonresonant scatter) to new direction */

      }

      if (modes.save_extract_photons && 1545.0 < 2.997925e18 / pp.freq && 2.997925e18 / pp.freq < 1565.0)
      {
        save_extract_photons (n, p, &pp, v);
      }

/* 68b - 0902 - ksl - turn phot_history on for the middle spectrum.  Note that we have to wait
 * to actually initialize phot_hist because the photon bundle is reweighted in extract_one */

      if (phot_history_spectrum == n)
      {
        phot_hist_on = 1;       // Start recording the history of the photon
      }


      /* Now collect polarization information for electron scattered photons only */
//      printf("extract: nres for this photon is %d \n", p->nres);

      if (p->nres == -1 && p->nscat > 0)
      {
//        Log("extract: entering poltest now... type is %d \n p nres is %d and pp nres is %d \n p ncat is %d and pp nscat is %d \n",
//            itype, p->nres, pp.nres, p->nscat, pp.nscat);
        poltest(p, &pp);
      } else {
        /* Zero the stokes q and u for this photon as its not electron scattered */
        pp.q = 0;
        pp.u = 0;
//        Log("extract: not electron scattered. stokes q is %f and u is %f \n",pp.q, pp.u);
      }

      /* Now extract the photon */
      extract_one (w, &pp, itype, n);

      /* Make sure phot_hist is on, for just one extraction */

      phot_hist_on = 0;

    }

  }
  return (0);
}





/**********************************************************/
/** 
 * @brief      (w,pp,itype,nspec)
 *
 * @param [in] WindPtr  w   The entire wind
 * @param [in] PhotPtr  pp  The photon to be extracted
 * @param [in] int  itype   The type of photon (star, disk, wind, etc)
 * @param [in] int  nspec   the spectrum which will be incremented
 * @return     The photon status after translation
 *
 * @details
 * extract_one is analogous to the detailed portion of transphot except here the
 * basic point is to calculate the optical depth through the plasma in a certain
 * direction, and to increment the appropriate spectrum.  
 *
 * Unlike trans_phot, this routine also checks to see if 
 * whether the photon hits the secondary star, if one exists.
 *
 * ### Notes ###
 * The logic behind the weighting of the photons is described in Christian Knigge's thesis in
 * section 2.3.4.  According to equation 2.19
 * 	Pc/Pw=12 cos(theta)*(1+b cos(theta)/(3+2b) where b=1.5 corresponds to the
 * Eddington approximation.
 * 
 * In Python, and in extract and transphot in particular, tau generally refers to the tau associated
 * with scattering processes, and the weight contains the effect of dimunition of the energy of
 * the photon bundle due to pure absorption processes.  So, in extract, we add pp->w * exp(-tau)
 * to the spectrum.
 *
 * Note that both linearly and logarithmically spaced spectra are produced.
 *
 **********************************************************/

int
extract_one (w, pp, itype, nspec)
     WindPtr w;
     PhotPtr pp;
     int itype, nspec;

{
  int istat, nres;
  struct photon pstart;
  double weight_min;
  int icell;
  int k, k1;
  double x[3];
  double tau;
  double zz;
  double dvds;
  double lfreqmin, lfreqmax, ldfreq;
  int ishell;
  double normal[3];


  weight_min = EPSILON * pp->w;
  istat = P_INWIND;
  tau = 0;
  icell = 0;

/* Preserve the starting position of the photon so one can use this to determine whether the
 * photon encountered the disk or star as it tried to exist the wind.
 */

  stuff_phot (pp, &pstart);

/* Reweight the photons. Note that photons have already been frequency shifted prior 
to entering extract */

  if (itype == PTYPE_STAR || itype == PTYPE_BL)
  {                             /* It was an unscattered photon from the star */
    stuff_v (pp->x, x);
    renorm (x, 1.);
    zz = fabs (dot (x, xxspec[nspec].lmn));
    pp->w *= zz * (2.0 + 3.0 * zz);     /* Eqn 2.19 Knigge's thesis */
  }
  else if (itype == PTYPE_DISK)
  {                             /* It was an unscattered photon from the disk */
    zz = fabs (xxspec[nspec].lmn[2]);
    pp->w *= zz * (2.0 + 3.0 * zz);     /* Eqn 2.19 Knigge's thesis */
  }
  else if (pp->nres > -1 && pp->nres < NLINES)  // added < NLINES condition for macro atoms (SS)
  {

/* It was a wind photon.  In this case, what we do depends
on whether it is a photon which arose via line radiation or 
some other process.

If geo.scatter_mode==SCATTER_MODE_ISOTROPIC then there is no need 
to reweight.  This is the isotropic assumption.  Otherwise, one
needs to reweight
*/

    if (geo.scatter_mode == SCATTER_MODE_THERMAL)
    {

      dvds = dvwind_ds (pp);
      ishell = pp->grid;
      tau = sobolev (&w[ishell], pp->x, -1.0, lin_ptr[pp->nres], dvds);
      if (tau > 0.0)
        pp->w *= (1. - exp (-tau)) / tau;
      tau = 0.0;
    }

/* But in any event we have to reposition wind photons so that they don't go through
the same resonance again */

    reposition (pp);            // Only reposition the photon if it was a wind photon
  }

  if (tau > TAU_MAX)
    istat = P_ABSORB;           /* Check to see if tau already too large */
  else if (geo.binary == TRUE)
    istat = hit_secondary (pp); /* Check to see if it hit secondary */


/* 68b - 0902 - ksl If we are trying to track the history of this photon, we need to initialize the
 * phot_hist.  We had to do this here, because we have just reweighted the photon
 */

  if (phot_hist_on)
  {
    phot_hist (pp, 0);          // Initialize the photon history
  }

/* Now we can actually extract the reweighted photon */

  while (istat == P_INWIND)
  {
    istat = translate (w, pp, 20., &tau, &nres);
    icell++;

    istat = walls (pp, &pstart, normal);
    if (istat == -1)
    {
      Error ("Extract_one: Abnormal return from translate\n");
      break;
    }

    if (pp->w < weight_min)
    {
      istat = P_ABSORB;         /*This photon was absorbed within the wind */
      break;
    }

    if (istat == P_HIT_STAR)
    {                           /* It was absorbed in the photosphere */
      break;
    }
    if (istat == P_HIT_DISK)
    {                           /* It was absorbed in the disk */
      break;
    }
    if (istat == P_SCAT)
    {                           /* Cause the photon to scatter and reinitilize */
      break;
    }
  }

  if (istat == P_ESCAPE)
  {

    if (!(0 <= tau && tau < 1.e4))
      Error_silent ("Warning: extract_one: ignoring very high tau  %8.2e at %g\n", tau, pp->freq);
    else
    {
      k = (pp->freq - xxspec[nspec].freqmin) / xxspec[nspec].dfreq;

      /* Force the frequency to be in range of that recorded in the spectrum */

      if (k < 0)
        k = 0;
      else if (k > NWAVE - 1)
        k = NWAVE - 1;


      lfreqmin = log10 (xxspec[nspec].freqmin);
      lfreqmax = log10 (xxspec[nspec].freqmax);
      ldfreq = (lfreqmax - lfreqmin) / NWAVE;

      /* find out where we are in log space */
      k1 = (log10 (pp->freq) - log10 (xxspec[nspec].freqmin)) / ldfreq;
      if (k1 < 0)
      {
        k1 = 0;
      }
      if (k1 > NWAVE - 1)
      {
        k1 = NWAVE - 1;
      }

      /* Increment the spectrum.  Note that the photon weight has not been diminished
       * by its passage through the wind, even though it may have encountered a number
       * of resonance, and so the weight must be reduced by tau
       */

      xxspec[nspec].f[k] += pp->w * exp (-(tau));       //OK increment the spectrum in question
      xxspec[nspec].lf[k1] += pp->w * exp (-(tau));     //And increment the log spectrum

      /* For polarization we need to update the stokes parameters Q and U.
       * Q[k] represents the total amount of polarization. It is calculated as a degree of polarization q
       * (a fraction between -1 and +1) multiplied by the intensity.
       * Here the intensity is represented by total flux, which is represented here as
       * the photon weight (w) multiplied by exp (-tau).
       * The same calculations are applied to log Q (lQ) and U and lU */

      xxspec[nspec].Q[k] += pp->q * pp->w * exp (-(tau));        //Increment the stokes Q parameter
      xxspec[nspec].lQ[k1] += pp->q * pp->w * exp (-(tau));      //and increment the log spectrum
      xxspec[nspec].U[k] += pp->u * pp->w * exp (-(tau));        //Increment the stokes U spectrum in question
      xxspec[nspec].lU[k1] += pp->u * pp->w * exp (-(tau));      //and increment the log spectrum



      /* If this photon was a wind photon, then also increment the "reflected" spectrum */
      if (pp->origin == PTYPE_WIND || pp->origin == PTYPE_WIND_MATOM || pp->nscat > 0)
      {

        xxspec[nspec].f_wind[k] += pp->w * exp (-(tau));        //OK increment the spectrum in question
        xxspec[nspec].lf_wind[k1] += pp->w * exp (-(tau));      //OK increment the spectrum in question

      }


      /* Records the total distance travelled by extracted photon if in reverberation mode */
      if (geo.reverb != REV_NONE)
      {
        if (pstart.nscat > 0 || pstart.origin > 9 || (pstart.nres > -1 && pstart.nres < nlines))
        {                       //If this photon has scattered, been reprocessed, or originated in the wind it's important
          pstart.w = pp->w * exp (-(tau));      //Adjust weight to weight reduced by extraction
          stuff_v (xxspec[nspec].lmn, pstart.lmn);
          delay_dump_single (&pstart, nspec);   //Dump photon now weight has been modified by extraction
        }
      }



/* 68b -0902 - ksl - turn phot_history off and store the information in the appropriate locations in the PlasmaPtrs
 * The reason this is here is that we only summarizes the history if the photon actually got to the observer
 */

      if (phot_hist_on)
      {
        phot_history_summarize ();
        phot_hist_on = 0;
      }


    }

  }


  if (istat > -1 && istat < 9)
    xxspec[nspec].nphot[istat]++;
  else
    Error
      ("Extract: Abnormal photon %d %8.2e %8.2e %8.2e %8.2e %8.2e %8.2e\n",
       istat, pp->x[0], pp->x[1], pp->x[2], pp->lmn[0], pp->lmn[1], pp->lmn[2]);

  return (istat);
}






