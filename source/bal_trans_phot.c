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
				!0 -> the normal option for python and hence the need to call "extract"
 
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

**************************************************************/

int
trans_phot (w, p, iextract)
     WindPtr w;
     PhotPtr p;
     int iextract;		/*  0 means do not extract along specific angles; nonzero implies to extract */

{
  double tau_scat, tau;
  int nphot, istat;
  double dot (), rrr;
  int translate ();
  int icell;
  int nres;
  int kkk, n, nscat;
  double weight_min;
  struct photon pp;
  double tau_scat_min, tau_scat_max, tau_actual_max;

  tau_scat_min = 100.;
  tau_actual_max = tau_scat_max = 0.0;
  nscat = 0;

  for (nphot = 0; nphot < NPHOT; nphot++)
    {

      stuff_phot (&p[nphot], &pp);
      /*Limit the resonance line search to +-VMAX cm/s */
      limit_lines (pp.freq * (1. - VMAX / C), pp.freq * (1. + VMAX / C));


      tau_scat = -log (1. - (rand () + 0.5) / MAXRAND);
      if (tau_scat_min > tau_scat)
	tau_scat_min = tau_scat;
      if (tau_scat_max < tau_scat)
	tau_scat_max = tau_scat;

      weight_min = EPSILON * pp.w;
      istat = P_INWIND;
      tau = 0;
      icell = 0;

      /* translate involves only a single shell (or alternatively a single tranfer in the windless region).
         istat as returned by should either be 0 in which case the photon hit the other side of the shell without scattering or
         1 in which case there was a scattering event in the shell, 2 in which case the photon reached the outside
         edge of the grid and escaped, 3 in which case it reach the inner edge and was reabsorbed.  If the photon escapes
         then we leave the photon at the postion of it's last scatter.  In most other cases though we store the final
         position of the photon. */
      while (istat == P_INWIND)
	{
	  istat = translate (w, &pp, tau_scat, &tau, &nres);
	  if (tau > tau_actual_max)
	    tau_actual_max = tau;
	  istat = walls (&pp, &p[nphot]);	//pp is where the photon is going, p is where it was

//                      ispy(&pp,n);

	  if (istat == -1)
	    {
	      Error_silent ("Abnormal return from translate on photon %d\n",
			    nphot);
	      break;
	    }

	  if (pp.w < weight_min)
	    {
	      istat = pp.istat = P_ABSORB;	/*This photon was absorbed within the wind */
	      pp.tau = VERY_BIG;
	      stuff_phot (&pp, &p[nphot]);
	      break;
	    }

	  if (istat == P_HIT_STAR)
	    {			/* It was absorbed in the photosphere */
	      stuff_phot (&pp, &p[nphot]);
	      break;
	    }

	  if (istat == P_HIT_DISK)
	    {			/* It was absorbed in the disk */
	      stuff_phot (&pp, &p[nphot]);
	      rrr = sqrt (dot (pp.x, pp.x));
	      kkk = 0;
	      while (rrr > disk.r[kkk] && kkk < NRINGS - 1)
		kkk++;
	      disk.nhit[kkk]++;
	      break;
	    }
	  if (icell > 2 * NDIM)
	    {
	      Error
		("Went too far without a scatter or hitting the boundary! %.1e\n",
		 sqrt (dot (pp.x, pp.x)));
	      istat = pp.istat = P_ERROR;
	      stuff_phot (&pp, &p[nphot]);
	      break;
	    }
	  if (istat == P_SCAT)
	    {			/* Cause the photon to scatter and reinitilize */
	      pp.nscat++;
	      if (nres > -1)
		{
		  pp.nrscat++;
		}
//                              scatter(w,&pp,nres);
	      if ((n = pp.grid) > -1)
		{
		  if (nres > -1)
		    {
		      line_heat (&plasmamain[w[n].nplasma], &pp, nres);
//OLD                 line_heat (&w[n], &pp, nres);
		      Log
			("Scattered photon: %8.2f resonance %3d %8.2f element %2d ion %2d density %8.2e\n",
			 C / pp.freq * 1e8, nres,
			 C / lin_ptr[nres]->freq * 1e8,
			 ion[lin_ptr[n]->nion].z,
			 ion[lin_ptr[n]->nion].istate,
			 plasmamain[w[n].nplasma].density[lin_ptr[n]->nion]);
//                       w[n].density[lin_ptr[n]->nion]);
		    }
		  if (where_in_grid (pp.x) != n)
		    Error
		      ("trans_phot: reweight has moved photon somehow %d %d istat %d or something else is odd\n",
		       pp.grid, p[nphot].grid, pp.istat);
		  stuff_phot (&pp, &p[nphot]);

		  nscat++;
		  tau_scat = -log (1. - (rand () + 0.5) / MAXRAND);
		  if (tau_scat_min > tau_scat)
		    tau_scat_min = tau_scat;
		  if (tau_scat_max < tau_scat)
		    tau_scat_max = tau_scat;
		  istat = P_INWIND;
		  tau = 0;
		  icell = 0;
		}
	      else
		{
		  Error
		    ("trans_phot: Trying to scatter a photon which is not in the wind\n");
		  Error ("trans_phot: grid %d x %8.2e %8.2e %8.2e\n", pp.grid,
			 pp.x[0], pp.x[1], pp.x[2]);
		  istat = pp.istat = p[nphot].istat = P_ERROR;
		}
	    }
	  if (pp.nscat == MAXSCAT)
	    {
	      istat = pp.istat = P_TOO_MANY_SCATTERS;	/* Scattered too many times */
	      Error ("Photon scattered too many times! %d r %e z %e\n",
		     nphot, sqrt (pp.x[0] * pp.x[0] + pp.x[1] * pp.x[1]),
		     pp.x[2]);
	      stuff_phot (&pp, &p[nphot]);
	      break;
	    }
	  p[nphot].istat = istat;
	}
    }

  Log ("trans_phot: %8.2g <tau_scat <%8.2g  tau max %8.2g nscat %d\n",
       tau_scat_min, tau_scat_max, tau_actual_max, nscat);
  return (0);
}
