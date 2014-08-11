

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "atomic.h"
#include "python.h"

/***********************************************************
                                       Space Telescope Science Institute

 Synopsis:

 int trans_phot(w,p,iextract) oversees the propagation of a "flight" of photons

Arguments:		

	PhotPtr p;
	WindPtr w;
	int iextract	0  -> the live or die option and therefore no need to call extract 
                       !0  -> the normal option for python and hence the need to call "extract"
 
Returns:
  
Description:	

This routine oversees the propagation of  individual photons.  The main loop 
covers an entire "flight" of photons.   The routine generates the random
optical depth a photon can travel before scattering and monitors the
progress of the photon through the grid.  The real physics is done else
where, in "translate" or "scatter".
		
Notes:

History:
 	97jan	ksl	Coded and debugged as part of Python effort.  
 	98mar	ksl	Modified to allow for photons which are created in the wind
	00nov	ksl	Cleaned up some of the logic when a photon was scattered
			and assured that photon did not continue when pure absorption
			was used.
	01aug	ksl	Modified Error messages (and made silent) when a photon was
			detected out of the cell in which it was originally.  This
			is a possible outcome (because a photon moves as far as the
			boundary + DFUDGE, but if this happens a lot one should worry.
	01dec	ksl	Modified calls to extract to use .origin where possible
			and to change meaning of itype in extract to conform	
			to that of origin.  This mainly just separates star + bl
			into individual sources.
	02jan	ksl	Added line to assure that trans_phot returns the final weight of
			all photons, even though the position may be the position of
			last scatter, or in the absence of scatters the originating 
			position of the photon.  This was to assure that photoabsorption
			and other pure absorption processes, are properly incoporated
			into the global spectrum and in live or die option.
	04mar	ksl	Made small modifications (qdisk) to reflect fact that heating
			of disk is now stored in a separate structure.
        04Jun   SS      Small modification to kill photons when necessary during the 
                        spectral cycles of macro atom runs.
	1112	ksl	Made some changes in the logic to try to trap photons that
			had somehow escaped the wind to correct a segmenation fault
			that cropped up in spherical wind models

**************************************************************/

FILE *pltptr;
int plinit = 0;

int
trans_phot (w, p, iextract)
     WindPtr w;
     PhotPtr p;
     int iextract;		/*  0 means do not extract along specific angles; 
				   nonzero implies to extract */

{
  double tau_scat, tau;
  int nphot, istat;
  double dot (), rrr;
  int translate ();
  int icell;
  int nres, *ptr_nres;
  int kkk, n;
  double weight_min;
  struct photon pp, pextract;
  double get_ion_density ();
  int nnscat;
  int disk_illum;
  int nerr;
  double p_norm, tau_norm;


  n = 0;			// To avoid -O3 warning

  /* 05jul -- not clear whether this is needed and why it is different from DEBUG */

  if (diag_on_off && plinit == 0)
    {
      pltptr = fopen ("python.xyz", "w");
      plinit = 1;
    }


  /* 130624 -- ksl - Chaanged the way the watchdog timeer works, so that it does not continually 
   * add to what is printed out, but overwrites it.  The printf statemen below is just so th
   * we do not overwrite the last line before the routine is entered.  Note tha fprintf(stdeerr
   * is required because stderr is flused immediately, whereas printf is normally only flushed by
   * \n.  NOte that I also increased the interval where items are bing printed out to 50,000, so 
   * this would be a bit less active
   */
  /* 130718 -- jm -for the moment I've changed this back, as we agreed that printing to 
   * stderr was bad practice.
   */


  Log ("\n");

  for (nphot = 0; nphot < NPHOT; nphot++)
    {

      //This is just a watchdog method to tell the user the program is still running
      //130306 - ksl since we don't really care what the frequencies are any more
      if (nphot % 50000 == 0)
	//OLD 130718 fprintf (stderr, "\rPhoton %7d of %7d or %6.3f per cent ", nphot, NPHOT,
	Log ("Photon %7d of %7d or %6.3f per cent \n", nphot, NPHOT,
	     nphot * 100. / NPHOT);

      Log_flush ();		/*NSH June 13 Added call to flush logfile */

      /* 74a_ksl Check that the weights are real */

      if (sane_check (p[nphot].w))
	{
	  Error ("trans_phot:sane_check photon %d has weight %e\n", nphot,
		 p[nphot].w);
	}
      /* Next block added by SS Jan 05 - for anisotropic scattering 
         with extract we want to be sure that everything is initialised
         (by scatter?) before calling extract for macro atom photons.
         Insert this call to scatter which should do this. */


      if (geo.rt_mode == 2 && geo.scatter_mode == 1)
	{
	  if (p[nphot].origin == PTYPE_WIND)
	    {
	      if (p[nphot].nres > -1 && p[nphot].nres < NLINES)
		{
		  geo.rt_mode = 1;
		  /*74a_ksl Check to see when a photon weight is becoming unreal */
		  if (sane_check (p[nphot].w))
		    {
		      Error
			("trans_phot:sane_check photon %d has weight %e before scatter\n",
			 nphot, p[nphot].w);
		    }
		  if ((nerr =
		       scatter (&p[nphot], &p[nphot].nres, &nnscat)) != 0)
		    {
		      Error
			("trans_phot: Bad return from scatter %d at point 1",
			 nerr);
		    }
		  /*74a_ksl Check to see when a photon weight is becoming unreal */
		  if (sane_check (p[nphot].w))
		    {
		      Error
			("trans_phot:sane_check photon %d has weight %e aftger scatter\n",
			 nphot, p[nphot].w);
		    }
		  geo.rt_mode = 2;
		}
	    }
	}





      stuff_phot (&p[nphot], &pp);
      disk_illum = geo.disk_illum;

      /* The next if statement is executed if we are calculating the detailed spectrum and
       * makes sure we always run extract on the original photon no matter where it was
       * generated 
       */

      if (iextract)
	{
	  //SS - for reflecting disk have to make sure disk photons are only extracted once!
	  if (disk_illum == 1 && p[nphot].origin == PTYPE_DISK)
	    {
	      geo.disk_illum = 0;
	    }


	  stuff_phot (&p[nphot], &pextract);


	  /* We then increase weight to account for number of scatters.
	     This is done because in extract we multiply by the escape
	     probability along a given direction, but we also need to divide the weight
	     by the mean escape probability, which is equal to 1/nnscat */
	  if (geo.scatter_mode == 2 && pextract.nres <= NLINES
	      && pextract.nres > 0)
	    {
	      /* we normalised our rejection method by the escape 
	         probability along the vector of maximum velocity gradient.
	         First find the sobolev optical depth along that vector */
	      tau_norm = sobolev (&wmain[pextract.grid], &pextract, -1.0,
				  lin_ptr[pextract.nres],
				  wmain[pextract.grid].dvds_max);

	      /* then turn into a probability */
          p_norm = p_escape_from_tau(tau_norm);

	    }
	  else
	    {
	      p_norm = 1.0;

	      /* throw an error if nnscat does not equal 1 */
	      if (p[nphot].nnscat != 1)
		Error ("nnscat is %i for photon %i in scatter mode %i!\n",
		       p[nphot].nnscat, nphot, geo.scatter_mode);
	    }



	  /* We then increase weight to account for number of scatters.
	     This is done because in extract we multiply by the escape
	     probability along a given direction, but we also need to divide the weight
	     by the mean escape probability, which is equal to 1/nnscat */
	  pextract.w *= p[nphot].nnscat / p_norm;

	  if (sane_check (pextract.w))
		{
		  Error
			("trans_phot: sane_check photon %d has weight %e before extract\n",
			 nphot, pextract.w);
		}
	  extract (w, &pextract, pextract.origin);


	  //SS if necessary put back the correcy disk illumination
	  if (disk_illum == 1 && p[nphot].origin == PTYPE_DISK)
	    {
	      geo.disk_illum = 1;
	    }
	}

      /* Initialize parameters that are needed for the flight of the photon through the
       * wind
       */

      tau_scat = -log (1. - (rand () + 0.5) / MAXRAND);
      weight_min = EPSILON * pp.w;
      istat = P_INWIND;
      tau = 0;
      icell = 0;

      /* This is the beginning of the loop for each photon and executes until the photon leaves the wind */

      while (istat == P_INWIND)
	{

	  /* translate involves only a single shell (or alternatively a single tranfer in the windless region).
	     istat as returned by should either be 
	     0 in which case the photon hit the other side of the shell without scattering or
	     1 in which case there was a scattering event in the shell, 
	     2 in which case the photon reached the outside edge of the grid and escaped, 
	     3 in which case it reach the inner edge and was reabsorbed.  
	     If the photon escapes then we leave the photon at the position of it's last scatter.  In most other cases 
	     though we store the final position of the photon. */


	  istat = translate (w, &pp, tau_scat, &tau, &nres);
//        Log("Photon=%i,weight=%e,tauscat=%f,nres=%i,istat=%i\n",nphot,p[nphot].w,tau_scat,nres,istat);
/* nres is the resonance at which the photon was stopped.  At present the
same value is also stored in pp->nres, but I have not yet eliminated 
it from translate. ?? 02jan ksl */



	  icell++;
	  istat = walls (&pp, &p[nphot]);
//pp is where the photon is going, p is where it was



#if DEBUG
	  ispy (&pp, n);
#endif

	  if (istat == -1)
	    {
	      Error_silent
		("trans_phot: Abnormal return from translate on photon %d\n",
		 nphot);
	      break;
	    }

	  if (pp.w < weight_min)
	    {
	      istat = pp.istat = P_ABSORB;	/*This photon was absorbed by continuum opacity within the wind */
	      pp.tau = VERY_BIG;
	      stuff_phot (&pp, &p[nphot]);
	      break;
	    }

	  if (istat == P_HIT_STAR)
	    {			/* It was absorbed by the star */
	      stuff_phot (&pp, &p[nphot]);
	      break;
	    }


	  if (istat == P_HIT_DISK)
	    {
	      /* It was absorbed by the disk */

	      /* Store the energy of the photon bundle into 
	       * a disk structure so that one can determine
	       * later how much and where the disk was heated
	       * by photons
	       */
/* Note that the disk is defined from 0 to NRINGS-2. NRINGS-1 contains the 
 * position of the outer radius of the disk.
 */
	      stuff_phot (&pp, &p[nphot]);
	      rrr = sqrt (dot (pp.x, pp.x));
	      kkk = 0;
	      while (rrr > qdisk.r[kkk] && kkk < NRINGS - 1)
		kkk++;
	      kkk--;		// So that the heating refers to the heating between kkk and kkk+1
	      qdisk.nhit[kkk]++;
	      qdisk.heat[kkk] += pp.w;	// 60a - ksl - Added to be able to calculate illum of disk
	      qdisk.ave_freq[kkk] += pp.w * pp.freq;
	      break;
	    }

	  if (istat == P_SCAT)
	    {			/* Cause the photon to scatter and reinitilize */


	      /* 71 - 1112 - ksl - placed this line here to try to avoid an error
	       * I was seeing in scatter.  I believe the first if statement has
	       * a loophole that needs to be plugged, when it comes back with avalue of
	       * n = -1
	       */
	      pp.grid = n = where_in_grid (pp.x);

//OLD71       if (pp.grid > -1)
//OLD71         {
//OLD71           if ((n = where_in_grid (pp.x)) != pp.grid)
//OLD71             {
//OLD71/* This error condition happens occassionally.  The reason is that we have added a "push through 
//OLD71 * distance" to get a photon to move on to the next cell, and have not updated the grid cell 
//OLD71 * before returning from translate. If one is concerned about this restore some of the lines
//OLD71 * that can befound in trans_phot in versions up through 54a.  04dec -- ksl
//OLD71 */
//OLD71               pp.grid = n;
//OLD71
//OLD71             }
//OLD71         }
//OLD71       else

	      if (n < 0)
		{
		  Error
		    ("trans_phot: Trying to scatter a photon which is not in the wind\n");
		  Error ("trans_phot: grid %3d x %8.2e %8.2e %8.2e\n",
			 pp.grid, pp.x[0], pp.x[1], pp.x[2]);
		  Error ("trans_phot: This photon is effectively lost!\n");
		  istat = pp.istat = p[nphot].istat = P_ERROR;
		  stuff_phot (&pp, &p[nphot]);
		  break;
		}

/*57+ -- ksl -- Add check to see if there is a cell in the plasma structure for this.  This is
a problem that needs fixing */

	      if (wmain[n].nplasma == NPLASMA)
		{
		  Error
		    ("trans_phot: Trying to scatter a photon which is not in a cell in the plasma structure\n");
		  Error ("trans_phot: grid %3d x %8.2e %8.2e %8.2e\n",
			 pp.grid, pp.x[0], pp.x[1], pp.x[2]);
		  Error ("trans_phot: This photon is effectively lost!\n");
		  istat = pp.istat = p[nphot].istat = P_ERROR;
		  stuff_phot (&pp, &p[nphot]);
		  break;
		}

/*57h -- ksl -- Add check to verify the cell has some volume */

	      if (wmain[n].vol <= 0)
		{
		  Error
		    ("trans_phot: Trying to scatter a photon in a cell with no wind volume\n");
		  Error ("trans_phot: grid %3d x %8.2e %8.2e %8.2e\n",
			 pp.grid, pp.x[0], pp.x[1], pp.x[2]);
		  Log ("istat %d\n", pp.istat);
		  Error ("trans_phot: This photon is effectively lost!\n");
		  istat = pp.istat = p[nphot].istat = P_ERROR;
		  stuff_phot (&pp, &p[nphot]);
		  break;

		}


/* SS July 04 - next lines modified so that the "thermal trapping" model of anisotropic
scattering is included in the macro atom method. What happens now is all in scatter - within
that routine the "thermal trapping" model is used to decide what the direction of emission is before
returning here.  54b-ksl -- To see what the code did previously see py46.  I've confirmed that
the current version of scattering really does what the old code did for two-level lines*/


	      nnscat = 0;
	      nnscat++;
	      ptr_nres = &nres;

	      /*74a_ksl - Check added to search for error in weights */
	      if (sane_check (pp.w))
		{
		  Error
		    ("trans_phot:sane_checl photon %d has weight %e before scatter\n",
		     nphot, pp.w);
		}
	      if ((nerr = scatter (&pp, ptr_nres, &nnscat)) != 0)
		{
		  Error ("trans_phot: Bad return from scatter %d at point 2",
			 nerr);
		}
	      pp.nscat++;
	      /* 74a_ksl - Check added to search for error in weights */

	      if (sane_check (pp.w))
		{
		  Error
		    ("trans_phot:sane_check photon %d has weight %e after scatter\n",
		     nphot, pp.w);
		}

	      /* SS June 04: During the spectrum calculation cycles, photons are thrown away
	         when they interact with macro atoms or become k-packets. This is done by setting
	         their weight to zero (effectively they have been absorbed into either excitation
	         energy or thermal energy). Since they now have no weight there is no need to follow
	         them further. */
	      /* 54b-ksl ??? Stuart do you really mean the comment above; it's not obvious to me
	       * since if true why does one need to calculate the progression of photons through 
	       * the wind at all???
	       * Also how is this enforced; where is pp.w set to a low value.
	       */

	      if (geo.matom_radiation == 1 && geo.rt_mode == 2
		  && pp.w < weight_min)
		/* Flag for the spectrum calculations in a macro atom calculation SS */
		{
		  istat = pp.istat = P_ABSORB;
		  pp.tau = VERY_BIG;
		  stuff_phot (&pp, &p[nphot]);
		  break;
		}

// Calculate the line heating and if the photon was absorbed break finish up
// ??? Need to modify line_heat for multiple scattering but not yet
// Condition that nres < nlines added (SS) 

	      if (nres > -1 && nres < nlines)
		{
		  pp.nrscat++;
/* This next statement writes out the position of every resonant scattering event to a file */
		  if (diag_on_off)
		    fprintf (pltptr,
			     "Photon %i has resonant scatter at %.2e %.2e %.2e in wind cell %i (grid cell=%i). Freq=%e Weight=%e\n",
			     nphot, pp.x[0], pp.x[1], pp.x[2],
			     wmain[n].nplasma, pp.grid, pp.freq, pp.w);

		  /* 68a - 090124 - ksl - Increment the number of scatters by this ion in this cell */
		  /* 68c - 090408 - ksl - Changed this to the weight of the photon at the time of the scatter */

		  plasmamain[wmain[n].nplasma].scatters[line[nres].nion] +=
		    pp.w;

		  if (geo.rt_mode == 1)	//only do next line for non-macro atom case
		    {
		      line_heat (&plasmamain[wmain[n].nplasma], &pp, nres);
		    }

		  if (pp.w < weight_min)
		    {
		      istat = pp.istat = P_ABSORB;	/*This photon was absorbed by continuum opacity within the wind */
		      pp.tau = VERY_BIG;
		      stuff_phot (&pp, &p[nphot]);
		      break;
		    }
		}


/* The next if statement causes photons to be extracted during the creation of the detailed spectrum
 * portion of the program
 */

/* N.B. To use the anisotropic scattering option, extract needs to follow scatter.  This
is because the reweighting which occurs in extract needs the pdf for scattering to have
been initialized. 02may ksl.  This seems to be OK at present.*/

	      if (iextract)
		{
		  stuff_phot (&pp, &pextract);


		  /* JM 1407 -- This next loop is required because in anisotropic 
		     scattering mode 2 we have normalised our rejection method. 
		     This means that we have to adjust nnscat by this factor,
		     since nnscat will be lower by a factor of 1/p_norm */
		  if (geo.scatter_mode == 2 && pextract.nres <= NLINES
		      && pextract.nres > 0)
		    {
		      /* we normalised our rejection method by the escape 
		         probability along the vector of maximum velocity gradient.
		         First find the sobolev optical depth along that vector */
		      tau_norm =
			sobolev (&wmain[pextract.grid], &pextract, -1.0,
				 lin_ptr[pextract.nres],
				 wmain[pextract.grid].dvds_max);

		      /* then turn into a probability */
              p_norm = p_escape_from_tau(tau_norm);

		    }
		  else
		    {
		      p_norm = 1.0;

		      /* throw an error if nnscat does not equal 1 */
		      if (p[nphot].nnscat != 1)
			Error
			  ("nnscat is %i for photon %i in scatter mode %i!\n",
			   p[nphot].nnscat, nphot, geo.scatter_mode);
		    }

		  /* We then increase weight to account for number of scatters.
		     This is done because in extract we multiply by the escape
		     probability along a given direction, but we also need to divide the weight
		     by the mean escape probability, which is equal to 1/nnscat */
		  pextract.w *= nnscat / p_norm;

		  if (sane_check (pextract.w))
		{
		  Error
			("trans_phot: sane_check photon %d has weight %e before extract\n",
			 nphot, pextract.w);
		}
		  extract (w, &pextract, PTYPE_WIND);	// Treat as wind photon for purpose of extraction
		}




/* OK we are ready to continue the processing of a photon which has scattered. The next steps
   reinitialize parameters so that the photon can continue throug the wind */

	      tau_scat = -log (1. - (rand () + 0.5) / MAXRAND);
	      istat = P_INWIND;
	      tau = 0;
	      reposition (w, &pp);
	      stuff_phot (&pp, &p[nphot]);
	      icell = 0;
	    }



/* This completes the portion of the code that handles the scattering of a photon
 * What follows is a simple check to see if this particular photon has gotten stuck
 * in the wind  54b-ksl */

	  if (pp.nscat == MAXSCAT)
	    {
	      istat = pp.istat = P_TOO_MANY_SCATTERS;	/* Scattered too many times */
	      stuff_phot (&pp, &p[nphot]);
	      break;
	    }

	  if (pp.istat == P_ADIABATIC)
	    {
	      istat = pp.istat = p[nphot].istat = P_ADIABATIC;
	      stuff_phot (&pp, &p[nphot]);
	      break;
	    }

/* This appears partly to be an insurance policy. It is not obvious that for example nscat and nrscat need to be updated */
	  p[nphot].istat = istat;
	  p[nphot].nscat = pp.nscat;
	  p[nphot].nrscat = pp.nrscat;
	  p[nphot].w = pp.w;	// Assure that final weight of photon is returned.

	}
      /* This is the end of the loop over individual photons */

    }
  /* This is the end of the loop over all of the photons; after this the routine returns */
  //130624 ksl Line added to complete watchdog timeer,
  Log ("\n\n");

  return (0);
}
