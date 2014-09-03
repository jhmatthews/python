
/***********************************************************
                                       Space Telescope Science Institute

 Synopsis:
	 radiation(p,ds) updates the radiation field parameters in the wind and reduces 
	 the weight of the photon as a result of the effect free free and photoabsorption.
	 radiation also keeps track of the number of photoionizations of h and he in the
	 cell. 
Arguments:

	PhotPtr p;	the photon
	double ds	the distance the photon has travelled in the cell

Returns:
	Always returns 0.  The pieces of the wind structure which are updated are
	j,ave_freq,ntot,  heat_photo,heat_ff,heat_h,heat_he1,heat_he2, heat_z,
	nioniz, and ioniz[].
 
Description:	 
Notes:
	
	The # of ionizations of a specific ion = (w(0)-w(s))*n_i sigma_i/ (h nu * kappa_tot).  (The # of ionizations
	is just given by the integral of n_i sigma_i w(s) / h nu ds, but if the density is assumed to
	be constant and sigma is also constant [therefy ignoring velocity shifts in a grid cell], then
	n_i sigma and h nu can be brought outside the integral.  Hence the integral is just over w(s),
	but that is just given by (w(0)-w(smax))/kappa_tot.)  The routine calculates the number of ionizations per
	unit volume.
	
	This routine is very time counsuming for our normal file which has a large number of x-sections.  The problem
	is going through all the do loops, not insofar as I can determine anything else.  Trying to reduce the calucatiions
	by using a density criterion is a amall effect.  One really needs to avoid the calculations, either by avoiding
	the do loops altogether or by reducing the number of input levels.  It's possible that if one were clever about
	the thresholds (as we are on the lines that one could figure out a wininning strategy as it is all brute force
        do loops..  

	57h -- This routine was very time consuming.  The time spent is/was well out of proportion to
	the affect free-bound processes have on the overall spectrum, and so signifincant changes have been made
	to the earlier versions of the routine.  The fundamental problem before 57h was that every time one
	entered the routine (which was every for every movement of a photon) in the code.  Basically there were
	3  changes made:
		1. During the detailed spectrum cycle, the code avoids the portions of the routine that are
	only needed during the ionization cycles.
		2. Switches have been installed tha skip the free-bound section altogether if PHOT_DEN_MIN is
	less than zero.  It is very reasonable to do this for the detailed spectrum calculation if one is
	doing a standard 
	were to skip portions of when they were not needed 

	?? Would it be more natural to include electron scattering here in Radiation as well ??
	?? Radiation needs a major overhaul.  A substantial portion of this routine is not needed in the extract 
	?? portion of the calculation.  In addition the do loops go through all of the ions checking one at a time
	?? whether they are above the frequencey threshold.  
	?? The solution I believe is to include some kind of switch that tells the routine when one is doing
	?? the ionization calculation and to skip the unnecessary sections when extract is being carried out.
	?? In addition, if there were a set of ptrs to the photionization structures that were orded by frequency,
	?? similar to the line list, one could then change to loops so that one only had to check up to the
	?? first x-section that had a threshold up to the photon frequency, and not all of the rest.
	?? At present, I have simply chopped the photoionizations being considered to include only thresholds
        ?? shortward of the Lyman limit...e.g. 1 Rydberg, but this makes it impossible to discuss the Balmer jump
History:
 	97jan	ksl	Coded as part of python effort
	98jul	ksl	Almost entirely recoded to eliminate arbitrary split between
			the several preexisting routines.
	98nov	ksl	Added back tracking of number of h and he ionizations
	99jan	ksl	Increased COLMIN from 1.e6 to 1.e10 and added checks so
			that one does not attempt to calculate photoionization
			cross-sections below threshold.  Both of these changes
			are intended to speed this routine.
        01oct   ksl     Modifications begun to incorporate Topbase photoionization
                        x-sections.
        01oct   ksl     Moved fb_cooling to balance_abso.  It's only used by
                        balance and probably should not be there either.
	02jun	ksl	Improved/fixed treatment of calculation of number of ionizations.
	04apr	ksl	SS had modified this routine to allow for macro-atoms in python_47, but his modifications
			left very little for radiation to accomplish.  I have returned to the old version of 
			routine, and moved the little bit that needed to be done in this routine for
			the macro approach to the calling routine.  Once we abandon the old approach this
			routine can be deleted.
	04aug	ksl	Fixed an error which had crept into the program between 51 and 51a that caused the
			disk models to be wrong.  The problem was that there are two places where the
			same frequency limit should have been used, but the limits had been set differently.
	04dec	ksl	54a -- Miniscule change to eliminate -03 warnings.  Also eliminate some variables
			that were not really being used.
	06may	ksl	57+ -- Began modifications to allow for splitting the wind and plasma structures
	06jul	ksl	57h -- Made various changes intended to speed up this routine. These include
			skipping sections of the routine in the spectrum generation
			phase that are not used, allowing control over whether fb opacities
			are calculated at all, and using frequency ordered pointer arrays
			to shorten the loops.
	11aug	nsh	70 changes made to radiation to allow compton cooling to be computed
	11aug	nsh	70 Changed printout of spectra in selected regions so it is always
			the midpoint of the wind
	12may	nsh	72 Added induced compton
	12jun 	nsh	72 Added lines to write out photon stats to a file dirung diagnostics. This is
			to allow us to see how well spectra are being modelled by power law W and alpha
	1405    JM corrected error (#73) so that photon frequency is shifted to the rest frame of the 
	        cell in question. Also added check which checks if a photoionization edge is crossed
	        along ds.
**************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "atomic.h"
#include "python.h"

#define COLMIN	0.01

int iicount = 0;

int
radiation (p, ds)
     PhotPtr p;
     double ds;
{
  PhotoionizationPtr x_ptr;
  TopPhotPtr x_top_ptr;

  WindPtr one;
  PlasmaPtr xplasma;

  double freq;
  double kappa_tot, frac_tot, frac_ff;
  double frac_z, frac_comp;	/* nsh 1108 added frac_comp - the heating in the cell due to compton heating */
  double frac_ind_comp;		/* nsh 1205 added frac_ind_comp - the heating due to induced compton heating */
  double kappa_ion[NIONS];
  double frac_ion[NIONS];
  double density, ft, tau, tau2;
  double energy_abs;
  int n, nion;
  double q, x, z;
  double w_ave, w_in, w_out;
  double den_config ();
  int nconf;
  double weight_of_packet, y;
  double v_inner[3], v_outer[3], v1, v2;
  double freq_inner, freq_outer;
  double freq_min, freq_max;
  double frac_path, freq_xs;
  struct photon phot;

  one = &wmain[p->grid];	/* So one is the grid cell of interest */
  xplasma = &plasmamain[one->nplasma];
  check_plasma (xplasma, "radiation");

  /* JM 140321 -- #73 Bugfix
     In previous version we were comparing these frequencies in the global rest fram
     The photon frequency should be shifted to the rest frame of the cell in question
     We currently take the average of this frequency along ds. In principle
     this could be improved, so we throw an error if the difference between v1 and v2 is large */

  /* calculate velocity at original position */
  vwind_xyz (p, v_inner);	// get velocity vector at new pos
  v1 = dot (p->lmn, v_inner);	// get direction cosine

  /* Create phot, a photon at the position we are moving to 
     note that the actual movement of the photon gets done after the call to radiation */
  stuff_phot (p, &phot);	// copy photon ptr
  move_phot (&phot, ds);	// move it by ds
  vwind_xyz (&phot, v_outer);	// get velocity vector at new pos
  v2 = dot (phot.lmn, v_outer);	// get direction cosine

  /* calculate photon frequencies in rest frame of cell */
  freq_inner = p->freq * (1. - v1 / C);
  freq_outer = phot.freq * (1. - v2 / C);

  /* take the average of the frequencies at original position and original+ds */
  freq = 0.5 * (freq_inner + freq_outer);



  /* calculate free-free, compton and ind-compton opacities 
     note that we also call these with the average frequency along ds */

  kappa_tot = frac_ff = kappa_ff (xplasma, freq);	/* Add ff opacity */
  kappa_tot += frac_comp = kappa_comp (xplasma, freq);	/* 70 NSH 1108 calculate compton opacity, store it in kappa_comp and also add it to kappa_tot, the total opacity for the photon path */
  kappa_tot += frac_ind_comp = kappa_ind_comp (xplasma, freq);
  frac_tot = frac_z = 0;	/* 59a - ksl - Moved this line out of loop to avoid warning, but notes 
				   indicate this is all diagnostic and might be removed */

  /* JM 1405 -- Check which of the frequencies is larger.
     if freq_max is always larger this can be removed. My checks
     indicate that it isn't */
  if (freq_outer > freq_inner)
    {
      freq_max = freq_outer;
      freq_min = freq_inner;
    }
  else
    {
      freq_max = freq_inner;
      freq_min = freq_outer;
    }


  if (freq > phot_freq_min)

    {

      if (geo.ioniz_or_extract)	// 57h -- ksl -- 060715
	{			// Initialize during ionization cycles only
	  for (nion = 0; nion < nions; nion++)
	    {
	      kappa_ion[nion] = 0;
	      frac_ion[nion] = 0;
	    }
	}
      /* Next section is for photoionization with Topbase.  There may be more
         than one x-section associated with an ion, and so one has to keep track
         of the energy that goes into heating electrons carefully.  */

      /* JM 1405 -- I've added a check here that checks if a photoionization edge has been crossed.
         If it has, then we multiply sigma*density by a factor frac_path, which is equal to the how far along 
         ds the edge occurs in frequency space  [(ft - freq_min) / (freq_max - freq_min)] */


      /* Next steps are a way to avoid the loop over photoionization x sections when it should not matter */
      if (DENSITY_PHOT_MIN > 0)	// 57h -- ksl -- 060715
	{			// Initialize during ionization cycles only


	  /* 57h -- 06jul -- ksl -- change loop to use pointers ordered by frequency */
	  for (n = 0; n < ntop_phot; n++)
	    {

	      x_top_ptr = phot_top_ptr[n];
	      ft = x_top_ptr->freq[0];

	      if (ft > freq_min && ft < freq_max)
		{
		  /* then the shifting of the photon causes it to cross an edge. 
		     Find out where between fmin and fmax the edge would be in freq space.
		     freq_xs is freq halfway between the edge and the max freq if an edge gets crossed */    
		  frac_path = (freq_max - ft) / (freq_max - freq_min);	
		  freq_xs = 0.5 * (ft + freq_max);
		}

	      else if (ft > freq_max)
		break;		// The remaining transitions will have higher thresholds

	      else if (ft < freq_min)
		{
		  frac_path = 1.0;	// then all frequency along ds are above edge
		  freq_xs = freq;	// use the average frequency
		}

	      if (freq_xs < x_top_ptr->freq[x_top_ptr->np - 1])
		{
		  /*Need the appropriate density at this point. */

		  nconf = x_top_ptr->nlev;
		  density = den_config (xplasma, nconf);

		  if (density > DENSITY_PHOT_MIN)
		    {
		      kappa_tot += x =
			sigma_phot_topbase (x_top_ptr,
					    freq_xs) * density * frac_path;

		      /* I believe most of next steps are totally diagnsitic; it is possible if 
		         statement could be deleted entirely 060802 -- ksl */

		      if (geo.ioniz_or_extract)	// 57h -- ksl -- 060715
			{	// Calculate during ionization cycles only

			  frac_tot += z = x * (freq_xs - ft) / freq_xs;
			  nion = config[nconf].nion;

			  if (nion > 3)
			    {
			      frac_z += z;
			    }

			  frac_ion[nion] += z;
			  kappa_ion[nion] += x;
			}

		    }


		}

	    }

/* Next section is for photoionization of those ions using VFKY values */

	  for (n = 0; n < nxphot; n++)
	    {
	      x_ptr = xphot_ptr[n];
	      nion = x_ptr->nion;

	      if (ion[nion].phot_info == 0)
		{		// Avoid ions with topbase x-sections
		  ft = x_ptr->freq_t;

		  if (ft > freq_min && ft < freq_max)
		    {
		      /* then the shifting of the photon causes it to cross an edge. 
		         Find out where between fmin and fmax the edge would be in freq space.
		         freq_xs is freq halfway between the edge and the max freq if an edge gets crossed */   
		      frac_path = (freq_max - ft) / (freq_max - freq_min);	
		      freq_xs = 0.5 * (ft + freq_max);
		    }

		  else if (ft > freq_max)
		    break;	// The remaining transitions will have higher thresholds

		  else if (ft < freq_min)
		    {
		      frac_path = 1.0;	// then all frequency along ds are above edge
		      freq_xs = freq;	// use the average frequency
		    }

		  density =
		    xplasma->density[nion] * ion[nion].g /
		    xplasma->partition[nion];

		  if (density > DENSITY_PHOT_MIN)
		    {
		      kappa_tot += x =
			sigma_phot (x_ptr, freq_xs) * density * frac_path;

		      /* Next if statment down to kappa_ion appers to be totally diagnostic - 060802 -- ksl */
		      if (geo.ioniz_or_extract)	// 57h -- ksl -- 060715
			{	// Calculate during ionization cycles only

			  frac_tot += z = x * (freq_xs - ft) / freq_xs;

			  if (nion > 3)
			    {
			      frac_z += z;
			    }
			  frac_ion[nion] += z;
			  kappa_ion[nion] += x;
			}
		    }
		}
	    }
	}

    }



  /* finished looping over cross-sections to calculate bf opacity 
     we can now reduce weights and record certain estimators */



  tau = kappa_tot * ds;
  w_in = p->w;

  if (sane_check (tau))
    {
      Error ("Radiation:sane_check CHECKING ff=%e, comp=%e, ind_comp=%e\n",
	     frac_ff, frac_comp, frac_ind_comp);
    }
/* Calculate the reduction in the w of the photon bundle along with the average
   weight in the cell */

  if (tau > 0.0001)
    {				/* Need differentiate between thick and thin cases */
      x = exp (-tau);
      p->w = w_out = w_in * x;
      energy_abs = w_in - w_out;
      w_ave = (w_in - w_out) / tau;
    }
  else
    {
      tau2 = tau * tau;
      p->w = w_out = w_in * (1. - tau + 0.5 * tau2);	/*Calculate to second order */
      energy_abs = w_in * (tau - 0.5 * tau2);
      w_ave = w_in * (1. - 0.5 * tau + 0.1666667 * tau2);
    }

  /*74a_ksl: 121215 -- Added to check on a problem photon */
  if (sane_check (p->w))
    {
      Error ("Radiation:sane_check photon weight is %e for tau %e\n", p->w,
	     tau);
    }

  if (geo.ioniz_or_extract == 0)
    return (0);			// 57h -- ksl -- 060715

/* Everything after this is only needed for ionization calculations */
/* Update the radiation parameters used ultimately in calculating t_r */

  xplasma->ntot++;




/* NSH 15/4/11 Lines added to try to keep track of where the photons are coming from, 
 * and hence get an idea of how 'agny' or 'disky' the cell is. */
/* ksl - 1112 - Fixed this so it records the number of photon bundles and not the total
 * number of photons.  Also used the PTYPE designations as one should as a matter of 
 * course
 */

  if (p->origin == PTYPE_STAR)
    xplasma->ntot_star++;
  else if (p->origin == PTYPE_BL)
    xplasma->ntot_bl++;
  else if (p->origin == PTYPE_DISK)
    xplasma->ntot_disk++;
  else if (p->origin == PTYPE_WIND)
    xplasma->ntot_wind++;
  else if (p->origin == PTYPE_AGN)
    xplasma->ntot_agn++;



  if (p->freq > xplasma->max_freq)	// check if photon frequency exceeds maximum frequency
    xplasma->max_freq = p->freq;

  /* JM -- 1310 -- check if the user requires extra diagnostics and
     has provided a file diag_cells.dat to store photons stats for cells they have specified
   */
  if (diag_on_off == 1 && ncstat > 0)
    {
      save_photon_stats (one, p, ds);	// save photon statistics (extra diagnostics)
    }


  /* JM 1402 -- the banded versions of j, ave_freq etc. are now updated in update_banded_estimators,
     which also updates the ionization parameters and scattered and direct contributions */

  update_banded_estimators (xplasma, p, ds, w_ave);


  if (sane_check (xplasma->j) || sane_check (xplasma->ave_freq))
    {
      Error ("radiation:sane_check Problem with j %g or ave_freq %g\n",
	     xplasma->j, xplasma->ave_freq);
    }


  if (kappa_tot > 0)
    {
      //If statement added 01mar18 ksl to correct problem of zero divide
      //  in odd situations where no continuum opacity
      z = (energy_abs) / kappa_tot;
      xplasma->heat_ff += z * frac_ff;
      xplasma->heat_tot += z * frac_ff;
      xplasma->heat_comp += z * frac_comp;	/* NSH 1108 Calculate the heating in the cell due to compton heating */
      xplasma->heat_tot += z * frac_comp;	/* NSH 1108 Add the compton heating to the total heating for the cell */
      xplasma->heat_tot += z * frac_ind_comp;	/* NSH 1205 Calculate the heating in the celldue to induced compton heating */
      xplasma->heat_ind_comp += z * frac_ind_comp;	/* NSH 1205 Increment the induced compton heating counter for the cell */
      if (freq > phot_freq_min)
//      if (freq > (CR / 100.)) //modified CR to CR/100 - SS June 04
	/* 04aug -- ksl -- 52.  Using CR/100 can speed the program up
	 * somewhat, but the limit here needs to be the same as the
	 * one above.  Setting the two limits differently can cause
	 * unpredictable and serious errors.
	 */
	{
	  xplasma->heat_photo += z * frac_tot;
	  xplasma->heat_z += z * frac_z;
	  xplasma->heat_tot += z * frac_tot;	//All of the photoinization opacities

	  /* Calculate the number of photoionizations per unit volume for H and He 
	     JM 1405 changed this to use freq_xs */
	  xplasma->nioniz++;
	  q = (z) / (H * freq * one->vol);
	  /* So xplasma->ioniz for each species is just 
	     (energy_abs)*kappa_h/kappa_tot / H*freq / volume
	     or the number of photons absorbed in this bundle per unit volume by this ion
	   */

	  for (nion = 0; nion < nions; nion++)
	    {
	      xplasma->ioniz[nion] += kappa_ion[nion] * q;
	      xplasma->heat_ion[nion] += frac_ion[nion] * z;
	    }

	}
    }

  /* Now for contribution to inner shell ionization estimators (SS, Dec 08) */
  for (n = 0; n < nauger; n++)
    {
      ft = augerion[n].freq_t;
      //Log("Auger tests: %g %g %g\n", augerion[n].freq_t, freq, p->freq);
      if (p->freq > ft)
	{
	  //Log("Adding a packet to AUGER via radiation %g \n", freq);

	  weight_of_packet = w_ave;
	  x = sigma_phot_verner (&augerion[n], freq);	//this is the cross section
	  y = weight_of_packet * x * ds;

	  xplasma->gamma_inshl[n] += y / (freq * H * one->vol);
	}
    }



  return (0);
}

/***********************************************************
                                       Space Telescope Science Institute

 Synopsis: 
	double kappa_ff(w,freq) calculates the free-free opacity allowing for stimulated emission
 	
Arguments:

Returns:
	
Description:	 
	Formula from Allen
	Includes only singly ionized H and doubly ionized he 	

Notes:

History:
	98aug	ksl	Coded as part of python effort
        04Apr   SS      If statement added for cases with hydrogen only.
                        Note: this routine assumes that H I is the first ion
                        and that He II is the fourth ion.
	05may	ksl	57+ -- Modified to eliminate WindPtr in favor
			of PlasmaPtr
   	12oct	nsh	73 -- Modified to use a prefector summed over all ions, calculated prior
			to the photon flight

**************************************************************/


double
kappa_ff (xplasma, freq)
     PlasmaPtr xplasma;
     double freq;
{
  double x;
  double exp ();
  double x1, x2, x3;

  if (gaunt_n_gsqrd == 0)	//Maintain old behaviour
    {
      if (nelements > 1)
	{
	  x = x1 =
	    3.692e8 * xplasma->ne * (xplasma->density[1] +
				     4. * xplasma->density[4]);
	}
      else
	{
	  x = x1 = 3.692e8 * xplasma->ne * (xplasma->density[1]);
	}
    }
  else
    {
      x = x1 = xplasma->kappa_ff_factor;
    }
  x *= x2 = (1. - exp (-H_OVER_K * freq / xplasma->t_e));
  x /= x3 = (sqrt (xplasma->t_e) * freq * freq * freq);

  return (x);
}

/***********************************************************
                                       Space Telescope Science Institute

 Synopsis: 
	double sigma_phot(x_ptr,freq)	calculates the photionization crossection due to the transition 
	associated with x_ptr at frequency freq
Arguments:

Returns:
	
Description:	 
	sigma_phot uses Verner et al.'s interpolation formulae for the photoionization crossection
	to calculate the bound free (or photoionization) optical depth.  The data must
	have been into the photoionization structures xphot with get_atomic_data and
	the densities of individual ions must have been calculated previously.

Notes:

History:
	98jul	ksl	Coded (actually moved from a subroutine called kappa_ds)
	06jul	ksl	57h -- Added code to avoid recalculating crossection
			if that is possible, i.e. if we are calculating
			a the x-section at the same frequency as last
			time the routine was entired for this ion. 

**************************************************************/

double
sigma_phot (x_ptr, freq)
     struct photoionization *x_ptr;
     double freq;
{
  double ft;
  double x, y;
  double f1, f2, f3;
  double xsection;


  ft = x_ptr->freq_t;		/* threshold frequency */


  if (ft < freq && freq < x_ptr->freq_max)
    {
      if (freq == x_ptr->f_last)
	return (x_ptr->sigma_last);	// Avoid recalculating xsection

      x = freq / x_ptr->freq0 - x_ptr->y0;
      y = sqrt (x * x + x_ptr->y1 * x_ptr->y1);

      /* This was line fixed by CK in 1998 Jul */
      f1 = (x - 1.0) * (x - 1.0) + x_ptr->yw * x_ptr->yw;

      //     f2 = pow (y, 0.5 * x_ptr->p - 5.5);
      f2 = exp ((0.5 * x_ptr->p - 5.5) * log (y));

//     f3 = pow (1.0 + sqrt (y / x_ptr->ya), -x_ptr->p);
      f3 = exp ((-x_ptr->p) * log ((1.0 + sqrt (y / x_ptr->ya))));


      xsection = x_ptr->sigma * f1 * f2 * f3;	// the photoinization xsection

/* Store crossesction for future use */
      x_ptr->sigma_last = xsection;
      x_ptr->f_last = freq;

      return (xsection);
    }
  else
    return (0.0);
}


/***********************************************************
				       Space Telescope Science Institute

 Synopsis:
	double sigma_phot_topbase(x_ptr,freq)	calculates the
	photionization crossection due to the transition associated with
	x_ptr at frequency freq
Arguments:
     struct topbase_phot *x_ptr;
     double freq;

Returns:

Description:
	sigma_phot uses the Topbase x-sections to calculate the bound free
	(or photoionization) xsection.	The data must have been into the
	photoionization structures xphot with get_atomic_data and the
	densities of individual ions must have been calculated previously.

Notes:

History:
	01Oct	ksl	Coded as part of general move to use Topbase data
			(especially partial xsections, which did not exist 
			in the Verner et al prescriptions
	02jul	ksl	Fixed error in the way fraction being applied.
			Sigh! and then modified program to use linterp

**************************************************************/

double
sigma_phot_topbase (x_ptr, freq)
     struct topbase_phot *x_ptr;
     double freq;
{
  int nmax;
  double xsection;
  double frac, fbot, ftop;
  int linterp ();
  int nlast;



  if (freq < x_ptr->freq[0])
    return (0.0);		// Since this was below threshold
  if (freq == x_ptr->f)
    return (x_ptr->sigma);	// Avoid recalculating xsection

  if (x_ptr->nlast > -1)
    {
      nlast = x_ptr->nlast;
      if ((fbot = x_ptr->freq[nlast]) < freq
	  && freq < (ftop = x_ptr->freq[nlast + 1]))
	{
	  frac = (log (freq) - log (fbot)) / (log (ftop) - log (fbot));
	  xsection =
	    exp ((1. - frac) * log (x_ptr->x[nlast]) +
		 frac * log (x_ptr->x[nlast + 1]));
	  //Store the results
	  x_ptr->sigma = xsection;
	  x_ptr->f = freq;
	  return (xsection);
	}
    }

/* If got to here, have to go the whole hog in calculating the x-section */
  nmax = x_ptr->np;
  x_ptr->nlast = linterp (freq, &x_ptr->freq[0], &x_ptr->x[0], nmax, &xsection, 1);	//call linterp in log space


  //Store the results
  x_ptr->sigma = xsection;
  x_ptr->f = freq;


  return (xsection);

}

/***********************************************************

  Synopsis: 
 	double sigma_phot_verner(x_ptr,freq)	calculates the photionization crossection due to the transition 
 	associated with x_ptr at frequency freq
 Arguments:
 
 Returns:
 
 Description:	 
        Same as sigma_phot but using the older compitation from Verner that includes inner shells
 
 Notes:
 
 History:
 	08dec	SS	Coded (actually mostly copied from sigma_phot)
 
**************************************************************/

double
sigma_phot_verner (x_ptr, freq)
     struct innershell *x_ptr;
     double freq;
{
  double ft;
  double y;
  double f1, f2, f3;
  double xsection;

  ft = x_ptr->freq_t;		/* threshold frequency */

  if (ft < freq)
    {
      y = freq / x_ptr->E_0 * HEV;

      f1 = ((y - 1.0) * (y - 1.0)) + (x_ptr->yw * x_ptr->yw);
      f2 = pow (y, 0.5 * x_ptr->P - 5.5 - x_ptr->l);
      f3 = pow (1.0 + sqrt (y / x_ptr->ya), -x_ptr->P);
      xsection = x_ptr->Sigma * f1 * f2 * f3;	// the photoinization xsection

      return (xsection);
    }
  else
    return (0.0);
}


/***********************************************************
				       Space Telescope Science Institute

Synopsis:

double den_config(one,nconf) returns the precalculated density
	of a particular "nlte" level.	If levels are not defined for an ion it
	assumes that all ions of are in the ground state.

Arguments:

Returns:

Description: The first case is when the density of the particular level
is in the wind array The second caseis when the density of the levels
for an ion are not tracked, and hence only the total density for the
ion is known.  In that case, we assume the ion is in it's ground state.


Notes:

History:
	01oct	ksl	Coded as part of effort to add topbase
			xsections and detailed balance to python
	05may	ksl	57+ -- Recoded to use PlasmaPtr

**************************************************************/

double
den_config (xplasma, nconf)
     PlasmaPtr xplasma;
     int nconf;
{
  double density;
  int nnlev, nion;

  nnlev = config[nconf].nden;
  nion = config[nconf].nion;

  if (nnlev >= 0)
    {				// Then a "non-lte" level with a density
      density = xplasma->levden[nnlev] * xplasma->density[nion];
    }
  else if (nconf == ion[nion].firstlevel)
    {
/* Then we are using a Topbase photoionization x-section for this ion, 
but we are not storing any densities, and so we assume it is completely in the
in the ground state */
      density = xplasma->density[nion];
    }
  else
    {
      density = 0;
    }

  return (density);


}


/***********************************************************
                Southampton University

Synopsis: pop_kappa_ff_array populates the multiplicative 	
		factor used in the FF calculation. This depends only
		on the densities of ions in the cell, and the electron
		temperature (which feeds into the gaunt factor) so it
		saves time to calculate all this just the once. This
		is called in python.c, just before the photons are 
		sent thruogh the wind.

Arguments:		

Returns:
 
Description:	

Notes:


History:
   12oct           nsh     coded 
   1212	ksl	Added sane check; note that this routine
   		is poorly documented.  Somewhere this 
		should be discribed better.  
   1407 nsh	changed loop to only go over NPLASMA cells not NPLASMA+1
 
**************************************************************/


double
pop_kappa_ff_array ()
{

  double gsqrd, gaunt, sum;
  int i, j;


  for (i = 0; i < NPLASMA ; i++)  //Changed to loop only up to NPLASMA, not NPLASMA+1
    {
      sum = 0.0;
      for (j = 0; j < nions; j++)
	{
	  if (ion[j].istate != 1)	//The neutral ion does not contribute
	    {
	      gsqrd =
		((ion[j].istate - 1) * (ion[j].istate -
					1) * RYD2ERGS) / (BOLTZMANN *
							  plasmamain[i].t_e);
	      gaunt = gaunt_ff (gsqrd);
	      sum +=
		plasmamain[i].density[j] * (ion[j].istate -
					    1) * (ion[j].istate - 1) * gaunt;
	      /* 74a_ksl  Added to diagnose problem with kappa_ff_fact producing NaN */
	      if (sane_check (sum))
		{
		  Error
		    ("pop_kappa_ff_array:sane_check sum is %e this is a problem, possible in gaunt %3\n",
		     sum, gaunt);
		}
	    }
	  else
	    {
	      sum += 0.0;	//add nothing to the sum if we have a neutral ion
	    }

	}
      plasmamain[i].kappa_ff_factor = plasmamain[i].ne * sum * 3.692e8;
    }

  return (0);
}




/***********************************************************
                Southampton University

Synopsis: 
	update_banded_estimators updates the estimators required for
	mode 7 ionization- the Power law, pairwise, modified saha ionization solver.
	It also records the values of IP.

Arguments:	
	xplasma		PlasmaPtr for the cell
	p 			Photon pointer
	ds 			ds travelled
	w_ave 		the weight of the photon. In non macro atom mode,
	            this is an average weight (passed as w_ave), but 
	            in macro atom mode weights are never reduced (so p->w 
	            is used).

Returns:
	increments the estimators in xplasma
 
Description:	
	This routine does not contain new code on initial writing, but instead 
	moves duplicated code from increment_bf_estimators and radiation() 
	to here, as duplicating code is bad practice.

Notes:


History:
   1402 JM 		Coding began
 
**************************************************************/

int
update_banded_estimators (xplasma, p, ds, w_ave)
     PlasmaPtr xplasma;
     PhotPtr p;
     double ds;
     double w_ave;
{
  int i;

  /*photon weight times distance in the shell is proportional to the mean intensity */
  xplasma->j += w_ave * ds;

  if (p->nscat == 0)
    {
      xplasma->j_direct += w_ave * ds;
    }
  else
    {
      xplasma->j_scatt += w_ave * ds;
    }



/* frequency weighted by the weights and distance       in the shell .  See eqn 2 ML93 */
  xplasma->mean_ds += ds;
  xplasma->n_ds++;
  xplasma->ave_freq += p->freq * w_ave * ds;



  /* 1310 JM -- The next loop updates the banded versions of j and ave_freq, analogously to routine inradiation
     nxfreq refers to how many frequencies we have defining the bands. So, if we have 5 bands, we have 6 frequencies, 
     going from xfreq[0] to xfreq[5] Since we are breaking out of the loop when i>=nxfreq, this means the last loop 
     will be i=nxfreq-1 */

  /* note that here we can use the photon weight and don't need to calculate anm attenuated average weight
     as energy packets are indisivible in macro atom mode */

  /* 71 - 111229 - ksl - modified to reflect fact that I have moved nxbands and xfreq into the geo structure */

  for (i = 0; i < geo.nxfreq; i++)
    {
      if (geo.xfreq[i] < p->freq && p->freq <= geo.xfreq[i + 1])
	{

	  xplasma->xave_freq[i] += p->freq * w_ave * ds;	/* 1310 JM -- frequency weighted by weight and distance */
	  xplasma->xsd_freq[i] += p->freq * p->freq * w_ave * ds;	/* 1310 JM -- input to allow standard deviation to be calculated */
	  xplasma->xj[i] += w_ave * ds;	/* 1310 JM -- photon weight times distance travelled */
	  xplasma->nxtot[i]++;	/* 1310 JM -- increment the frequency banded photon counter */

	  /* 1311 NSH lines added below to work out the range of frequencies within a band where photons have been seen */
	  if (p->freq < xplasma->fmin[i])
	    {
	      xplasma->fmin[i] = p->freq;
	    }
	  if (p->freq > xplasma->fmax[i])
	    {
	      xplasma->fmax[i] = p->freq;
	    }

	}
    }

  /* NSH 131213 slight change to the line computing IP, we now split out direct and scattered - this was 
     mainly for the progha_13 work, but is of general interest */
  /* 70h -- nsh -- 111004 added to try to calculate the IP for the cell. Note that 
   * this may well end up not being correct, since the same photon could be counted 
   * several times if it is rattling around.... */

  /* 1401 JM -- Similarly to the above routines, this is another bit of code added to radiation
     which originally did not get called in macro atom mode. */

  /* NSH had implemented a scattered and direct contribution to the IP. This doesn't really work 
     in the same way in macro atoms, so should instead be thought of as 
     'direct from source' and 'reprocessed' radiation */

  if (HEV * p->freq > 13.6)	// only record if above H ionization edge
    {

      /* IP needs to be radiation density in the cell. We sum wcontributions from
         each photon, then it is normalised in wind_update. */
      xplasma->ip += ((w_ave * ds) / (H * p->freq));

      if (p->nscat == 0)
	{
	  xplasma->ip_direct += ((w_ave * ds) / (H * p->freq));
	}
      else
	{
	  xplasma->ip_scatt += ((w_ave * ds) / (H * p->freq));
	}
    }

  return (0);
}




/***********************************************************
                Southampton University

Synopsis: 
	save_photon_stats prints photon statistics to a file

Arguments:	
	One 		WindPtr for the cell
	p 			Photon pointer
	ds 			ds travelled

Returns:
 
Description:
   the loop below is if the user requires extra diagnostics and
   has provided a file diag_cells.dat to store photons stats for cells they have specified

Notes:
   Moved here to save duplicating code between bf_estimators_increment and radiation.

History:
   1402 JM 		Coding began
 
**************************************************************/


int
save_photon_stats (one, p, ds)
     WindPtr one;
     PhotPtr p;
     double ds;
{
  int i;

  /* JM -- 1310 -- the loop below is if the user requires extra diagnostics and
     has provided a file diag_cells.dat to store photons stats for cells they have specified
   */

  for (i = 0; i < ncstat; i++)
    {
      /* check if the cell is in the specified list - ncell_stats is global variable */
      if (one->nplasma == ncell_stats[i])
	{
	  fprintf (pstatptr,
		   "PHOTON_DETAILS %3d %8.3e %8.3e %8.3e cell%3d wind cell%3d\n",
		   geo.wcycle, p->freq, p->w, ds, one->nplasma, one->nwind);
	}
    }
  return (0);
}


/*************************************************************
Synopsis: 
	mean_intensity returns a value for the mean intensity 

Arguments:	
	xplasma 		PlasmaPtr for the cell - supplies spectral model
	freq 			the frequency at which we want to get a value of J
	mode 			mode 1=use BB if we have not yet completed a cycle
				and so dont have a spectral model, mode 2=never use BB

Returns:
 
Description:
   This subroutine returns a value for the mean intensity J at a 
   given frequency, using either a dilute blackbody model
   or a spectral model depending on the value of geo.ioniz_mode. 
   to avoid code duplication.

Notes:
   This subroutine was produced
   when we started to use a spectral model to populaste the upper state of a
   two level atom, as well as to calculate induced compton heating. This was

History:
   1407 NSH 		Coding began
 
**************************************************************/



double
mean_intensity (xplasma, freq, mode)
     PlasmaPtr xplasma;		// Pointer to current plasma cell
     double freq;		// Frequency of the current photon being tracked
     int mode;			// mode 1=use BB if no model, mode 2=never use BB

{
  int i;
  double J;
  double expo;

  J = 0.0;			// Avoid 03 error


  if (geo.ioniz_mode == 5 || geo.ioniz_mode == 7)	/*If we are using power law ionization, use PL estimators */
    {
      if (geo.wcycle > 0)	/* there is only any point in worrying if we have had at least one cycle otherwise there is no model */
	{
	  for (i = 0; i < geo.nxfreq; i++)
	    {
	      if (geo.xfreq[i] < freq && freq <= geo.xfreq[i + 1])	//We have found the correct model band
		{
		  if (xplasma->spec_mod_type[i] > 0)	//Only bother if we have a model in this band
		    {

		      if (freq > xplasma->fmin_mod[i] && freq < xplasma->fmax_mod[i])	//The spectral model is defined for the frequency in question
			{

			  if (xplasma->spec_mod_type[i] == SPEC_MOD_PL)	//Power law model
			    {
			      J = pow (10, (xplasma->pl_log_w[i] + log10 (freq) * xplasma->pl_alpha[i]));
			    }

			  else if (xplasma->spec_mod_type[i] == SPEC_MOD_EXP)	//Exponential model
			    {
			      J = xplasma->exp_w[i] * exp ((-1 * H * freq) / (BOLTZMANN * xplasma->exp_temp[i]));
			    }
			  else
			    {
			      Error("mean_intensity: unknown spectral model (%i) in band %i\n",
				         xplasma->spec_mod_type[i], i);
			      J = 0.0;	//Something has gone wrong
			    }
			}

		      else
			{
			  /* We have a spectral model, but it doesnt apply to the frequency 
			     in question. clearly this is a slightly odd situation, where last
			     time we didnt get a photon of this frequency, but this time we did. 
			     Still this should only happen in very sparse cells, so induced compton 
			     is unlikely to be important in such cells. We generate a warning, just 
			     so we can see if this is happening a lot */
			   J = 0.0;

			  /* JM140723 -- originally we threw an error here. No we count these errors and 
			     in wind_updates because you actually expect 
			     it to happen in a converging run */
			  nerr_Jmodel_wrong_freq++;
			}
		    }
		  else		/* There is no model in this band - this should not happen very often  */
		    {
		      J = 0.0;	//There is no model in this band, so the best we can do is assume zero J

		      /* JM140723 -- originally we threw an error here. No we count these errors and 
			     in wind_updates because you actually expect 
			     it to happen in a converging run */
		      nerr_no_Jmodel++;
		    }



		}
	    }
	}
      else			//We have not completed an ionization cycle, so no chance of a model
	{
	  if (mode == 1)	//We need a guess, so we use the initial guess of a dilute BB
	    {
	      expo = (H * freq) / (BOLTZMANN * xplasma->t_r);
	      J = (2 * H * freq * freq * freq) / (C * C);
	      J *= 1 / (exp (expo) - 1);
	      J *= xplasma->w;
	    }
	  else			//A guess is not a good idea (i.e. we need the model for induced compton), so we return zero.
	    {
	      J = 0.0;
	    }

	}
    }

  else				/*Else, use BB estimator of J */
    {
      expo = (H * freq) / (BOLTZMANN * xplasma->t_r);
      J = (2 * H * freq * freq * freq) / (C * C);
      J *= 1 / (exp (expo) - 1);
      J *= xplasma->w;
    }

  return J;
}


