
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "atomic.h"
#include "python.h"

/***********************************************************
                                       Space Telescope Science Institute

 Synopsis:

int extract(w,p,itype) is the supervisory routine which helps normally
	builds detailed spectra as the photons are transit the wind.

Arguments:		
	PhotPtr p;	The initial photon
	WindPtr w;
	int itype	PTYPE_STAR->the photon came for the star 
			PTYPE_BL->the photon came from the boundary layer
			PTYPE_DISK->the photon being redirected arose in the disk,
			PTYPE_WIND->the photon being redirected arose in the wind,
Returns:
 
 
Description:	

extract is called when a photon begins it's flight and every time that photon
scatters, unless the user has exercised the "live or die" option, in
which case it is not called.  extract calls extract_one for each spectrum
it wants to build, where the actual incrementing of the spectrum is done.

An option allows one to construct spectra which have undergone a certain
number of scatters.  The parameters for this option all come in through
python.h.  
	If s[n].nscat>999; then one obtains the entire spectrum
	if s[n].nscat is a positive number, then the spectrum is composed just 
		of photons with that number of scatters
	if s[n].nscat is a negative number, then the spectrum contains all photons
		with >=  the absolute value of s[n].nscat
Another option allows one to construct spectra from photons which have
	scattered above or below the plane of the disk.  As for the number of
	scatters this option comes through python.h, in this case from s[n].top_bot
	s[n].top_bot=0 -> accept all photons
	                  >0  -> accept photons last above the disk
	                  <0  -> accept photons below the disk
Notes:

History:
 	97march ksl	Coded and debugged as part of Python effort.  
 	97aug28	ksl	Added the option which allows one to construct spectra with a
			specific number of scatters.
	97sep1	ksl	Added the option which allows one to construct spectra which
			originate from a specific position above or below the disk
	02jan	ksl	Fixed problems with itypes. (Photons which were wind photons
			for the purpose of extract were not being properly labelled
			as such, and this meant they were not doppler shifted properly.
			This produced photons that were further from resonance in the
			local frame and made the wind more transmissive, especially
			in the blue wing of the line.  It suggests
			that a lot of the absorption in extract is fairly local and
			that the blue wing structure could be fairly sensitive to how
			the doppler shift is carried out.)
	09feb	ksl	68b - Added hooks to track energy deposition of extracted photons
			in the wind 

**************************************************************/


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


  /* 68b -09021 - ksl - The next line selects the middle inclination angle for recording the absorbed enery */
  phot_history_spectrum = 0.5 * (MSPEC + nspectra);

  for (n = MSPEC; n < nspectra; n++)
    {
      /* If statement allows one to choose whether to construct the spectrum
         from all photons or just from photons that have scattered a specific number
         of times or in specific regions of the wind. */

      yep = 1;			// Start by assuming it is a good photon for extraction

      if ((mscat = xxspec[n].nscat) > 999 || p->nscat == mscat
	  || (mscat < 0 && p->nscat >= (-mscat)))
	yep = 1;
      else
	yep = 0;

      if (yep)
	{
	  if ((mtopbot = xxspec[n].top_bot) == 0)
	    yep = 1;		// Then there are no positional parameters and we are done
	  else if (mtopbot == -1 && p->x[2] < 0)
	    yep = 1;
	  else if (mtopbot == 1 && p->x[2] > 0)
	    yep = 1;
	  else if (mtopbot == 2)	// Then to count, the photom must originate within sn.r of sn.x
	    {
	      vsub (p->x, xxspec[n].x, xdiff);
	      if (length (xdiff) > xxspec[n].r)
		yep = 0;

	    }
	  else
	    yep = 0;
	}



      if (yep)			//Then we want to extract this photon
	{


/* Create a photon pp to use here and in extract_one.  This assures we
 * have not modified p as part of extract
 */

	  stuff_phot (p, &pp);
	  stuff_v (xxspec[n].lmn, pp.lmn);	/* Stuff new photon direction into pp */

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
	    {			/* If the photon was scattered in the wind, 
				   the frequency also must be shifted */
	      vwind_xyz (&pp, v);	/*  Get the velocity at the position of pp */
	      doppler (p, &pp, v, pp.nres);	/*  Doppler shift the photon -- test! */

/*  Doppler shift the photon (as nonresonant scatter) to new direction */

	    }

	  if (modes.save_extract_photons && 1545.0 < 2.997925e18 / pp.freq
	      && 2.997925e18 / pp.freq < 1565.0)
	    {
	      fprintf (epltptr,
		       "%3d %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %7.2f %7.2f \n",
		       n, p->x[0], p->x[1], p->x[2], v[0], v[1], v[2],
		       p->lmn[0], p->lmn[1], p->lmn[2], pp.lmn[0],
		       pp.lmn[1], pp.lmn[2], 2.997925e18 / p->freq,
		       2.997925e18 / pp.freq);
	    }

/* 68b - 0902 - ksl - turn phot_history on for the middle spectrum.  Note that we have to wait
 * to actually initialize phot_hist because the photon bundle is reweighted in extract_one */

	  if (phot_history_spectrum == n)
	    {
	      phot_hist_on = 1;	// Start recording the history of the photon
	    }

	  /* Now extract the photon */

	  extract_one (w, &pp, itype, n);

	  /* Make sure phot_hist is on, for just one extraction */

	  phot_hist_on = 0;

	}

    }
  return (0);
}

/***********************************************************
                                       Space Telescope Science Institute

 Synopsis:

extract_one(w,pp,itype,nspec)

Arguments:		
	PhotPtr p;
	WindPtr w;
	int itype		0->the photon came from the star 
				1->the photon came from the boundary layer
				2->the photon being redirected arose in the disk,
				3->the photon being redirected arose in the wind,
	int nspec		the spectrum which will be incremented
Returns:
 
 	The photon status after translation
 
Description:	

extract_one is analogous to the detailed portion of transphot except here the
basic point is to calculate the optical depth through the plasma in a certain
direction, and to increment the appropriate spectrum.  There are a few differences
though because here we also check whether the photon hits the secondary star.
		
Notes:

The logic behind the weighting of the photons is described in Christian Knigge's thesis in
section 2.3.4.  According to equation 2.19
	Pc/Pw=12 cos(theta)*(1+b cos(theta)/(3+2b) where b=1.5 corresponds to the
Eddington approximation.

In python, and in extract and transphot in particular, tau generally refers to the tau associated
with scattering processes, and the weight contains the effect of dimunition of the energy of
the photon bundle due to pure absorption processes.  So, in extract, we add pp->w * exp(-tau)
to the spectrum.

History:
 	97march ksl	Coded and debugged as part of Python effort.  
 	97july	ksl	Included the possibility that the photon was absorbed by the secondary.
 	97sep28	ksl	Modified weightings of extracted photons to be those in Knigge's thesis for
			Eddington approximation
	02jan2	ksl	Adapted extract to use photon types

**************************************************************/



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
  int k;
  double x[3];
  double tau;
  double zz;
  double dvds;
  int ishell;


  weight_min = EPSILON * pp->w;
  istat = P_INWIND;
  tau = 0;
  icell = 0;

/* Preserve the original position of the photon so one can use this to determine whether the
 * photon encountered the disk or star as it tried to exist the wind.
 */

  stuff_phot (pp, &pstart);

/* Reweight the photons. Note that photons have already been frequency shifted prior 
to entering extract */

  if (itype == PTYPE_STAR || itype == PTYPE_BL)
    {				/* It was an unscattered photon from the star */
      stuff_v (pp->x, x);
      renorm (x, 1.);
      zz = fabs (dot (x, xxspec[nspec].lmn));
      pp->w *= zz * (2.0 + 3.0 * zz);	/* Eqn 2.19 Knigge's thesis */
    }
  else if (itype == PTYPE_DISK)
    {				/* It was an unscattered photon from the disk */
      zz = fabs (xxspec[nspec].lmn[2]);
      pp->w *= zz * (2.0 + 3.0 * zz);	/* Eqn 2.19 Knigge's thesis */
    }
  else if (pp->nres > -1 && pp->nres < NLINES)	// added < NLINES condition for macro atoms (SS)
    {

/* It was a wind photon.  In this case, what we do depends
on whether it is a photon which arose via line radiation or some other process.

If geo.scatter_mode==0 then there is no need to reweight.  This is the
isotropic assumption.

NB--It is important that reweightwind be called after scatter, as there
are variables which are set in scatter and in aniosowind that are
used by reweightwind.  02may ksl
*/

      if (geo.scatter_mode == 1)
	{			// Then we have anisotropic scattering
/* In new call it is important to realize that pp->lmn must be
the new photon direction, and that the weight of the photon will
have been changed */
	  reweightwind (pp);
	}

      else if (geo.scatter_mode == 2)	/* Then we have anisotropic
					   scattering based on a random number of scatters at the scattering
					   site */
	{

	  dvds = dvwind_ds (pp);
	  ishell = pp->grid;
	  tau = sobolev (&w[ishell], pp->x, -1.0, lin_ptr[pp->nres], dvds);
	  if (tau > 0.0)
	    pp->w *= (1. - exp (-tau)) / tau;
	  tau = 0.0;
	}
	
/* But in any event we have to reposition wind photons so thath they don't go through
the same resonance again */

      reposition (w, pp);	// Only reposition the photon if it was a wind photon
    }

  if (tau > TAU_MAX)
    istat = P_ABSORB;		/* Check to see if tau already too large */
  else if (geo.system_type == 1)	/* Changed 69 to allow for additional system types */
    istat = hit_secondary (pp);	/* Check to see if it hit secondary */


/* 68b - 0902 - ksl If we are trying to track the history of this photon, we need to initialize the
 * phot_hist.  We had to do this here, because we have just reweighted the photon
 */

  if (phot_hist_on)
    {
      phot_hist (pp, 0);	// Initialize the photon history
    }

/* Now we can actually extract the reweighted photon */

  while (istat == P_INWIND)
    {
      istat = translate (w, pp, 20., &tau, &nres);
      icell++;

      istat = walls (pp, &pstart);
      if (istat == -1)
	{
	  Error ("Extract_one: Abnormal return from translate\n");
	  break;
	}

      if (pp->w < weight_min)
	{
	  istat = P_ABSORB;	/*This photon was absorbed within the wind */
	  break;
	}

      if (istat == P_HIT_STAR)
	{			/* It was absorbed in the photosphere */
	  break;
	}
      if (istat == P_HIT_DISK)
	{			/* It was absorbed in the disk */
	  break;
	}
      if (istat == P_SCAT)
	{			/* Cause the photon to scatter and reinitilize */
	  break;
	}
    }

  if (istat == P_ESCAPE)
    {

      if (!(0 <= tau && tau < 1.e4))
	Error_silent
	  ("Warning: extract_one: ignoring very high tau  %8.2e at %g\n",
	   tau, pp->freq);
      else
	{
	  k = (pp->freq - xxspec[nspec].freqmin) / xxspec[nspec].dfreq;

	  /* Force the frequency to be in range of that recorded in the spectrum */

	  if (k < 0)
	    k = 0;
	  else if (k > NWAVE - 1)
	    k = NWAVE - 1;

	  /* Increment the spectrum.  Note that the photon weight has not been diminished
	   * by its passage through th wind, even though it may have encounterd a number
	   * of resonance, and so the weight must be reduced by tau
	   */

	  xxspec[nspec].f[k] += pp->w * exp (-(tau));	//OK increment the spectrum in question

    	/* This seems very defensive.  Is tau ever less than 0? */
	
		if (pp->nrscat > 0)		// SWM - Records total distance travelled by extract photon
		{
			stuff_v(pstart.x, pp->x);	// Restore photon to initial position (necessary for reweighting schemes)
			pp->path = pstart.path;

			if (geo.reverb > REV_NONE)		// only want to dump photon if the reverb structures are set up
			  delay_dump_single(pp, 1);	// Dump photon now weight has been modified
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
       istat, pp->x[0], pp->x[1], pp->x[2], pp->lmn[0], pp->lmn[1],
       pp->lmn[2]);

  return (istat);
}
