


#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "atomic.h"
#include "python.h"

/* These routines have to do with calculating where resonances are and how strong they are.  
   These routines should be made geometry independent and therefore should work with 
   very minor (header only) modifications in other configurations. */


/***********************************************************
                                       Space Telescope Science Institute

 Synopsis:

	double calculate_ds(w,p,phot,tau_scat,tau,nres,smax,istat)
	calculates the distance a photon will travel within a single shell,

Arguments:

	WindPtr w			the ptr to the structure defining the wind
	PhotPtr p,phot;			the photon at the near and the far edge of the shell
	double tau_scat			the optical depth at which the photon will scatter
	double smax			the maximum distance the photon can travel in the cell

Returns:

Calculate returns the distance the photon travels as the main result of the routine.

It also returns   
	tau, the optical depth at the calculated distance, 
	nres the number of the resonance if there was one, or -1 if there was not a resonant
		scatter in the cell.  
	istat 0 if the photon did not scatter,  1 if it did scatter.
 
Description:	 

taking
into account the resonances that appear within that shell.  Any paramaters that depend explicitly on the
coordinate gridding, such as the maximum distance the photon can travel before hitting the edge of the
shell,  should be calculated outside of this routine.


Notes:

Calculate_ds does not modify the p or phot in any way!!  It is not
entirely obvious to me that this is a good plan, but that is the
way it is.

History:
 	97jan	ksl	Coded as part of python effort
	02mar	ksl	Made small changes intended to reduce the number
			of times dvwind_ds was called in attempt to speed
			up program. 
	02jun	ksl	Added additional changes to prevent some
			of the more time consuming calculations
			when possible
        04Mar   SS      Changes for including macro atom approach.
                        Photoionisation included as a process to remove
                        photon packets are specific points.
                        tau_es is optical depth per unit length for e- scattering
                        tau_bf is optical depth per unit length for bf processes
                        tau_cont is sum of tau_es and tau_bf
	04Apr	ksl	Moved from a formalism where bf transitions are identified
			by nlines+ncont to NLINES+ncont, since presumably the latter
			will be easier to remember when debugging different types
			of input files.
			Moved the calculation of bf_opacities and of the selection
			of which continuum process causes scattering to a separate
			subroutine.  My motivation here is to keep the logic of
			calculate_ds as similar as possible to old versions of
			the code.
         04Apr   SS     Added ff as a source of opacity. ff always leads to creating of
                        a k-packet.
                        Renamed the continuum things from tau to kap since this makes
                        more sense.
                        Changed call to kappa_bf so that it includes all the topbase
                        bfs and not just the macro atom ones.
         04May   SS     Modifications to allow the macro_simple option. (All ions treated as
                        "simple").
	06may	ksl	57+ -- To allow for plasma structure.  
**************************************************************/

struct photon cds_phot_old;
double cds_v2_old, cds_dvds2_old;



double
calculate_ds (w, p, tau_scat, tau, nres, smax, istat)
     WindPtr w;			//w here refers to entire wind, not a single element
     PhotPtr p;
     double tau_scat, *tau;
     int *nres;
     double smax;
     int *istat;
{
  int kkk;
  double kap_es;
  double freq_inner, freq_outer, dfreq, ttau, freq_av;
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
  int bb_estimators_increment ();
  WindPtr one, two;
  int check_in_grid;
  int nplasma;
  PlasmaPtr xplasma, xplasma2;

  one = &w[p->grid];		//Get a pointer to the cell where the photon bundle is located.
  nplasma = one->nplasma;
  xplasma = &plasmamain[nplasma];

  kap_es = THOMPSON * xplasma->ne;
  /* This is the electron scattering opacity per unit length. For the Macro Atom approach we need an 
     equivalent opacity per unit length due to each of the b-f continuua. Call it kap_bf. (SS) */


  ttau = *tau;
  ds_current = 0;
  init_dvds = 0;
  dvds1 = dvds2 = 0.0;		// To avoid a -03 compile warning
  *nres = -1;
  *istat = P_INWIND;

  if (ttau < 0.0)
    {
      Error
	("calculate_ds:  odd tau  %8.2e at %g entering calculate_ds\n", ttau,
	 p->freq);
    }


/* So "phot"  is a photon vector at the far edge of the cell, while p remains the photon 
   vector at the near edge of the cell, and p_now is the midpoint.  Note that comp_phot
   compares the position and direction of two photons.  If they are the same, then 
   it just takes v1 from the old value.  */

  if (comp_phot (&cds_phot_old, p))
    {
      vwind_xyz (p, v_inner);
      v1 = dot (p->lmn, v_inner);
    }
  else
    {
      v1 = cds_v2_old;
    }

  /* Create phot, a photon at the far side of the cell */
  stuff_phot (p, &phot);
  move_phot (&phot, smax);
  vwind_xyz (&phot, v_outer);
  v2 = dot (phot.lmn, v_outer);

  /* Check to see that the velocity is monotonic across the cell
   * by calculating the velocity at the midpoint of the path
   *
   * If it is, then reduce smax
   */
  vc = C;
  while (vc > VCHECK && smax > DFUDGE)
    {
      stuff_phot (p, &p_now);
      move_phot (&p_now, smax / 2.);
      vwind_xyz (&p_now, v_check);
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
frame of rest. Therefore multiply. See doppler notes for a discussion

02feb: The sign appears correct.  If the photon is moving in the same 
direction as v, then in the rest frame of the ions, 
then the photon frequency will be less. */


  freq_inner = p->freq * (1. - v1 / C);
  freq_outer = phot.freq * (1. - v2 / C);
  dfreq = freq_outer - freq_inner;

/* The next section checks to see if the frequency difference on
 * the two sides is very small and if not limits the resonances
 * one has to worry about
 */

  if (fabs (dfreq) < EPSILON)
    {
      Error
	("translate: v same at both sides of cell %d\n",one->nwind); /*NSH 130724 shortened error statement, was causing issues with formatting */

/* so dfreq is %2g,\n v_inner %.2g %.2g %.2g v_outer %.2g %.2g %.2g \n",
	 one->nwind, dfreq, v_inner[0], v_inner[1], v_inner[2], v_outer[0],
	 v_outer[1], v_outer[2]); This is the part of the above error statement cut out */
      x = -1;
      return (smax);		// This is not really the best thing to do, but it avoids disaster below

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

//nline_min, nline_max, and nline_delt are found in atomic.h and are set by limit_lines()


/* Next part deals with computation of bf opacity. In the macro atom method this is needed.
In the old method it is not. This section activates if geo.rt_mode==2 (switch for macro atom
method). If the macro atom method is not used just get kap_bf to 0 and move on). SS*/

  kap_bf_tot = 0;
  kap_ff = 0;


  if (geo.rt_mode == 2)
    {
      /* Potentially several continuum may contribute in a given frequency range so the kap_bf is an array. 
         Also need to store the total - kap_bf_tot. */


      freq_av = freq_inner;	//(freq_inner + freq_outer) * 0.5;  //need to do better than this perhaps but okay for star - comoving frequency (SS)


      kap_bf_tot = kappa_bf (xplasma, freq_av, 0);
      kap_ff = kappa_ff (xplasma, freq_av);

      /* Okay the bound free contribution to the opacity is now sorted out (SS) */
    }

  if (one->vol == 0)
    {
      kap_bf_tot = kap_ff = 0.0;
      Error_silent ("ds_calculate vol = 0: cell %d position %g %g %g\n",
		    p->grid, p->x[0], p->x[1], p->x[2]);
      /* !BUG (SS June 04) */
      /* This shouldn't happen - but very occassionally it does. Seems to be a minor problem relating to
         the grid cells at the edge of the wind. Possibly requires more accurate volume calculation. */
    }


  kap_cont = kap_es + kap_bf_tot + kap_ff;	//total continuum opacity 



/* Finally begin the loop over the resonances that can interact with the
     photon in the cell */

  for (n = 0; n < nline_delt; n++)
    {
      nn = nstart + n * ndelt;	/* So if the frequency of resonance increases as we travel through
				   the grid cell, we go up in the array, otherwise down */
      x = (lin_ptr[nn]->freq - freq_inner) / dfreq;

      if (0. < x && x < 1.)
	{			/* this particular line is in resonance */
	  ds = x * smax;


/* Before checking for a resonant scatter, need to check for scattering due to a continuum
process. */
	  if (ttau + (kap_cont) * (ds - ds_current) > tau_scat)
	    {
	      /* then the photon was scattered by the continuum before reaching the resonance 
	         Need to randomly select the continumm process which caused the photon to
	         scatter.  The variable threshold is used for this. */

	      *nres =
		select_continuum_scattering_process (kap_cont, kap_es,
						     kap_ff, xplasma);
	      *istat = P_SCAT;	//flag as scattering
	      ds_current += (tau_scat - ttau) / (kap_cont);	//distance travelled
	      ttau = tau_scat;
	      *tau = ttau;
	      return (ds_current);
	    }
	  else
	    {

	      /* increment tau by the continuum optical depth to this point */
	      ttau += kap_cont * (ds - ds_current);	/*kap_cont used here rather than kap_es */


	      ds_current = ds;	/* At this point ds_current is exactly the position of the resonance */
	      kkk = lin_ptr[nn]->nion;


	      /* The density is calculated in the wind array at the center of a cell.  We use
	         that as the first estimate of the density.  So don't delete the next line.  However, we 
	         can do better except  at the edges of the grid.   01apr07 ksl.  ??? One could still do 
	         better than this since even if we can usually interpolate in one direction if not two ???
	       */
	      stuff_phot (p, &p_now);
	      move_phot (&p_now, ds_current);	// So p_now contains the current position of the photon


	      //?? It looks like this if statement could be moved up higher, possibly above the xrho line if willing to accept cruder dd estimate
	      //If the density of the ion is very small we shouldn't have to worry about a resonance, but otherwise
	      // ?? This seems like an incredibly small number; how can anything this small affect anything ??

	      dd = get_ion_density (&p_now, kkk);

	      if (dd > LDEN_MIN)
		{
		  /* If we have reached this point then we have to initalize dvds1 and dvds2. Otherwise
		     there is no need to do this, especially as dvwind_ds is an expensive calculation time wise */

		  if (init_dvds == 0)
		    {
		      dvds1 = dvwind_ds (p);
		      dvds2 = dvwind_ds (&phot);
		      init_dvds = 1;
		    }

		  dvds = (1. - x) * dvds1 + x * dvds2;

		  /* Add the line optical depth  Note that one really should translate the photon to this point 
		     before doing this (?? What is "this"??)as p-> x is being used to calculate direction of the wind */

		  ttau += tau_sobolev =
		    sobolev (one, p, dd, lin_ptr[nn], dvds);

		  /* tau_sobolev now stores the optical depth. This is fed into the next statement for the bb estimator
		     calculation. SS March 2004 */

		  /* ksl - ??? 0902 - Stuart it looks to mee as if this is being run in the Macro case even during the 
		   * exraction cycles.  Could you check. It shouldn't be necessary.  ???
		   */

		  if (geo.rt_mode == 2)	//Macro Atom case (SS)
		    {

		      /*
		         Because push through distance may take us out of the cell we want, need to make sure
		         that the cell is correct before incrementing the heating rate/estimators. So 1st check if
		         it's still in the wind and second get a pointer to the grid cell where the resonance really happens.
		       */

		      check_in_grid = walls (&p_now, p);

		      if (check_in_grid != P_HIT_STAR
			  && check_in_grid != P_HIT_DISK
			  && check_in_grid != P_ESCAPE)
			{
			  two = &w[where_in_grid (p_now.x)];
			  xplasma2 = &plasmamain[two->nplasma];

			  if (lin_ptr[nn]->macro_info == 1
			      && geo.macro_simple == 0)
			    {
			      /* The line is part of a macro atom so increment the estimator if desired (SS July 04). */
			      if (geo.ioniz_or_extract == 1)
				{
				  bb_estimators_increment (two, p,
							   tau_sobolev, dvds,
							   nn);
				}
			    }
			  else
			    {
			      /* The line is from a simple ion. Record the heating contribution and move on. */

			      bb_simple_heat (xplasma2, p, tau_sobolev, dvds,
					      nn);

			    }
			}
		    }
		  /* Completed special calculations for the Macro Atom case */

		  /* 68b - 0902 - The next section is to track where absorption is taking place along the line of sight
		   * to the observer.  It is probably possibly to simplify some of what is happening here, as we
		   * have various photons real and imaginary in this subroutine.  p, the orginal photon, phot the
		   * photon at the opposide edge of the cell and p_now hte photon at its current position.  Some
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
ds_current is either 0 if there were no resonances or the postion of the 
"last" resonance if there were resonances.  But we need to check one
last time to see if it was scattered by continuum process.
Note: ksl -- It is generally a bad policy to have a bunch of repeated code
like this.  We should probably rewrite some of this in terms of subroutines
for clarity, especially the bit that calculates where the continuum scattering
event occurred.  04 apr

 */

  if (ttau + kap_cont * (smax - ds_current) > tau_scat)
    {
      *nres =
	select_continuum_scattering_process (kap_cont, kap_es, kap_ff,
					     xplasma);

      /* A scattering event has occurred in the shell  and we remain in the same shell */
      ds_current += (tau_scat - ttau) / (kap_cont);
      *istat = P_SCAT;		/* Flag for scattering (SS) */
      ttau = tau_scat;
    }
  else
    {				/* Then we did hit the other side of the shell  
				   (or possibly the another wall of the same shell) */
      *istat = P_INWIND;
      ttau += kap_cont * (smax - ds_current);	/* kap_es replaced with kap_cont (SS) */
      ds_current = smax;

    }

  *tau = ttau;

  stuff_phot (&phot, &cds_phot_old);	// Store the final photon position
  cds_v2_old = v2;		// and the velocity along the line of sight

  return (ds_current);

}


/***********************************************************
                                       Space Telescope Science Institute

Synopsis:
	select_continuum_scattering_process() what the scattering
	process was that caused the scatter and returns this to
	the main routine.

Arguments:


Returns:


Description:	 


Notes:
	ksl--Accuracy might be improved if one used a specific position to
	calculate the opacity.  It is not obvious that these are the
	variables we want to transmit in the routine but this should
	be OK for now



History:
	04apr	ksl	Moved code that had been in calculate_ds to
			separate routine.  
        04Apr   SS      Modified to include possibility of ff
                        For a ff event nres = -2
                        Changed some variable names to kap rather
                        than tau so that it makes more sense.
        04Nov   SS      Modified to take the wind pointer argument since
                        the meaning of the elements of kap_bf could now
                        vary from cell to cell.

**************************************************************/
int
select_continuum_scattering_process (kap_cont, kap_es, kap_ff, xplasma)
     double kap_cont, kap_es, kap_ff;
     PlasmaPtr xplasma;
{
  int nres;
  double threshold;
  double run_tot;
  int ncont;

  threshold = ((rand () + 0.5) / MAXRAND) * (kap_cont);

  /* First check for electron scattering. */
  if (kap_es > threshold)
    {				/* electron scattering event occurred (SS) */
      nres = -1;		// flag electron scatterin (SS) */
    }
  /* Now chech for ff. */
  else if ((kap_es + kap_ff) > threshold)
    {
      nres = -2;
    }
  /* Now check for bf. */
  else
    {				/* use a running sum to find which photoionisation process it was */
      /* 
         If a non-macro-atom run is being done this part should never be reached.
         Just do a check that all is well - this can be removed eventually (SS)
       */
      if (geo.rt_mode == 1)
	{
	  Error
	    ("calculate_ds: Not using macro atoms but trying to excite one? Aboort.\n");
	  exit (0);		//hopefully this will never happen and this check can be deleted at some
	  //point (SS)
	}
      run_tot = kap_es + kap_ff;
      ncont = 0;
      while (run_tot < threshold)
	{
	  run_tot += kap_bf[ncont];
	  ncont++;
	}
      /* When it gets here know that excitation is in photoionisation labelled by ncont */
/* ksl 04apr: I don't particularly like this approach for labelling a photoionization transition.
It implies that coding and then decoding needs to be done accurately.  It is quite likely that
we should expand the types of scattering events, e.g P_SCAT --> P_RES, P_CONT. It is a bit of
a conundrum though because we have previously only had two types of scatters, and so for example
doppler used nres to define how to doppler shift. In the end, I decided to keep Stuart's 
approach but but to use a fixed offset of NLINES rather than nlines.  This has the modest
advantage that photoionizations always start at the same place.  ksl*/
      nres = NLINES + 1 + xplasma->kbf_use[ncont - 1];	//modified SS Nov 04
    }
  return (nres);
}


/***********************************************************
                                       Space Telescope Science Institute

Synopsis:
	kappa_bf(one, freq,macro_all) calculates the bf opacity in a specific
	cell. 


Arguments:

	macro_all	1--> macro_atoms only
			0--> all topbase ions


Returns:


Description:	 


Notes:
	ksl--Accuracy might be improved if one used a specific position to
	calculate the opacity.  


History:
	04apr	ksl	Moved code that had been in calculate_ds to
			separate routine.  
        04Apr   SS      Changed some variable names to kap rather than tau to 
                        make more sense.

**************************************************************/
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


  kap_bf_tot = 0;		//initalise to 0 (SS)

  macro_all--;			// Subtract one from macro_all to avoid >= in for loop below.


  //if (freq > CR)
  //  {

  for (nn = 0; nn < xplasma->kbf_nuse; nn++)	// Loop over photoionisation processes. 
    // This is mostly copied from old radiation.c (SS)
    {
      n = xplasma->kbf_use[nn];
      ft = phot_top[n].freq[0];	//This is the edge frequency (SS)

      kap_bf[nn] = 0.0;

      if (freq > ft && freq < phot_top[n].freq[phot_top[n].np - 1]
	  && phot_top[n].macro_info > macro_all)
	{
	  /*Need the appropriate density at this point. */

	  nconf = phot_top[n].nlev;	//Returning lower level = correct (SS)

	  density = den_config (xplasma, nconf);	//Need to check what this does (SS)


	  if (density > DENSITY_PHOT_MIN || phot_top[n].macro_info == 1)
	    {

	      /*          kap_tot += x = (delete) */
	      kap_bf[nn] = x = sigma_phot_topbase (&phot_top[n], freq) * density;	//stimulated recombination? (SS)
	      kap_bf_tot += x;
	    }
	}
    }

  //}

  return (kap_bf_tot);
}

/***********************************************************
                                       Space Telescope Science Institute

Synopsis:
	kbf_need(one) computes and stores the set of bound-free processes which
        make significiant contributions to the opacity in a grid cell.


Arguments:

       w - a pointer to the wind
       fmin - the lowest frequency of interest to the calculation
       fmax - the highest frequency of interest to the calculations


Returns:
	For each cell, the routine determines what bf transitons are inportant
	and stores them in one->kbf_use[n].  The total number of such transitions
	is given in one->kbf_nuse.


Description:	 

      This routine is now called before a Monte Carlo calculation. It determines which 
      bf processes are worth considering during the calculation that follows. To do this
      is uses the edge position (must be inside, or close to, the spectral region of 
      interest) and the edge opacity (must be greater than a threshold value, currently
      10^-6.

      The need for the routine is just to prevent wasting time on unimportant bf
      transitions.  It is called (currently from resonate).


Notes:
	The dimensionalty of kbf_use is currently set to NTOP_PHOT, which is large enough
	to store evey transition if necessary.


History:
        04Nov   SS	work began
	04dec	ksl	54b-- Small modification.  There is no need to search in 2d, one d
			is fine, since vol for edge cells is zero
	06may	ksl	57 -- Began work to match to plama structure.  Uses vol but only
			to decide which cells to calculate for so I switched to summing
			over the plaama structure.
                   

**************************************************************/
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
  int nplasma;


  for (nplasma = 0; nplasma < NPLASMA; nplasma++)	// Loop over all the cells in the wind
    {
      xplasma = &plasmamain[nplasma];
      one = &wmain[xplasma->nwind];
      nuse = 0;

      for (n = 0; n < ntop_phot; n++)	// Loop over photoionisation processes. 
	{

	  ft = phot_top[n].freq[0];	//This is the edge frequency (SS)

	  if ((ft > (fmin / 3.)) && (ft < fmax))
	    {


	      nconf = phot_top[n].nlev;	//Returning lower level = correct (SS)

	      density = den_config (xplasma, nconf);

	      tau_test =
		phot_top[n].x[0] * density * SMAX_FRAC * length (one->xcen);


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

/***********************************************************
                                       Space Telescope Science Institute

 Synopsis:

double sobolev(p,w,l,dvds)  calculates tau associated with a resonance, given the 
conditions in the wind and the direction of the photon.

It does not modify any of the variables that are passed to it, including for example
the photon.  

Arguments:

	PhotPtr p		p is a single photon, whose position is that of the resonanace
	WindPtr w		w is a WindPtr to the cell in which the photon resides; it is not 
				the whole array 
	struct lines *l	ptr to the transition  which has a resonance at this point
	double dvds	the velocity gradient in the direction of travel of the photon

Returns:

	The optical depth associated with a transition

Description:	 

In the sobolev approx  tau = PI e^2 / m_e c  * NL * f * lambda * 1 / (|dv/dx|) where dv/ds 
is the velocity gradient in the direction the photon is travelling.  This can be made slightly 
simpler by noting that lambda = freq / c and so  tau = PI e^w / m_e   

Notes:

	One worries that the Sobolev length, e.g delta lambda/lambda c/dvds might lead to 
	problems if it is too long, e.g. greater than the size of the system.  I have done tests 
	to see whether limiting the minimum dvds made much difference in generating knigge9b.  At
	very high values e.g. 3, it essentially makes the wind non scattering.  At values of
	10**-3 or less, it made no difference in the models I looked at.   02apr-ksl

History:
 	97jan	ksl	Coded as part of python effort
	00dec	ksl	Modified to account for population
			of lines, i.e. two_level_atom
	01dec	ksl	Modified to account for changed call to two_level_atom.  Note
			that the calling routine resonate interpolates the density
			of the ion in question, which accounts for the shuffling
			of densities.  Made it really work correctly 02jan
	02may	ksl	Added possibility of having sobolev calculate an accurate
			density
	06may	ksl	57+ -- Began mods for plasma structure

**************************************************************/

double
sobolev (one, p, den_ion, lptr, dvds)
     WindPtr one;		// This is a single cell in the wind
     PhotPtr p;
     double den_ion;
     struct lines *lptr;
     double dvds;
{
  double tau, xden_ion;
  double two_level_atom (), d1, d2;
  int nion;
  double d_hold;
  int nplasma;
  PlasmaPtr xplasma;

  nplasma = one->nplasma;
  xplasma = &plasmamain[nplasma];

  if ((dvds = fabs (dvds)) == 0.0)	// This forces dvds to be positive -- a good thing!
    {
      tau = VERY_BIG;
      Error ("Sobolev: Surprize tau = VERY_BIG\n");
    }

  else if (lptr->macro_info == 1 && geo.rt_mode == 2 && geo.macro_simple == 0)
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
      d_hold = xplasma->density[nion];	// Store the density of this ion in the cell
      if (den_ion < 0)
	{
	  xplasma->density[nion] = get_ion_density (p, lptr->nion);	// Forced calculation of density 
	}
      else
	{
	  xplasma->density[nion] = den_ion;	// Put den_ion into the density array
	}
      two_level_atom (lptr, xplasma, &d1, &d2);	// Calculate d1 & d2
      xplasma->density[nion] = d_hold;	// Restore w
    }



  xden_ion = (d1 - lptr->gl / lptr->gu * d2);

  if (xden_ion < 0)
    {
      Error ("sobolev: den_ion has negative density %g %g %g %g %g %g\n",
	     d1, d2, lptr->gl, lptr->gu, lptr->freq, lptr->f);

      /*SS July 08: With macro atoms, the population solver can default to d2 = gu/gl * d1 which should
         give exactly zero here but can be negative, numerically. 
         So I'm modyfying this to set tau to zero in such cases, when the populations are vanishingly small anyway. */
      tau_x_dvds = PI_E2_OVER_M * d1 * lptr->f / (lptr->freq);
      tau = tau_x_dvds / dvds;
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

//  tau = PI_E2_OVER_M * xden_ion * lptr->f ;
//  tau /= (lptr->freq*dvds);

// N.B. tau_x_dvds is an external variable and is used in anisotropic
// scattering  0902 - ksl - I believe the question here, which has
// been here for a long time  is whether this is 
// good programming practice ????
  tau_x_dvds = PI_E2_OVER_M * xden_ion * lptr->f / (lptr->freq);
  tau = tau_x_dvds / dvds;

  return (tau);
}


/***********************************************************
                                       Space Telescope Science Institute

 Synopsis:

	doppler(pin,pout,v,nres) calculate the doppler shift given the direction of the incoming 
	and outgoing photons and the local  velocity of the wind 

Arguments:

	PhotPtr pin,pout	the incoming and outgoing photon
	double v[];		the velocity of the wind where the scatter occurred
	int nres;			either the number of the scatter, or a nonresonant scatter
						if nres < 0

Returns:

	pout->freq is updated

Description:	 

	Given the incoming and the outgoing photon direction, 
	doppler calculates the outgoing frequency (to first order in beta).  
	There are two basic cases, depending on whether it was a resonant 
	or a nonresonant scatter.  

Notes:
	This is more obscure than I first thought because second order
	terms can matter.  So for example, if as originally (up through
	python_43.0) one tries to convert from the global frame to the
	local rest frame and back using the formulas
		freq2= freq1 * (1-v/c)
		freq3= freq2 * (1-v/c)
	One finds that freq3=freq1 (1-v**2/c**2)...which is not accurate
	enough.  Therfore if you want to use a first order calculation
	you have to say.
		freq2= freq1 * (1-v/c)
		freq3= fre12 / (1-v/c)
	Alternatively of course, you could use the special relativistic
	formulat of
		freq2=freq1 * (1-v/c)/sqrt(1-v/c**2)
	but maybe this is OK for now.

History:
 	97jan	ksl	Coded as part of python effort
	02apr	ksl	Made changes to attempt to assure that
			doppler shifts are exact in first order.
			See. notes.
        04May   SS      Changes made to include macro atom approach
                        (doppler shift of bf and ff emission).
***********************************************************/


int
doppler (pin, pout, v, nres)
     PhotPtr pin, pout;
     double v[];
     int nres;

{
  double dot ();
//  double ftemp;
// double beta;
//  double q[3];

  if (nres == -1)		//Electron scattering (SS)
    {				/*It was a non-resonant scatter */
      pout->freq =
	pin->freq * (1 - dot (v, pin->lmn) / C) / (1 -
						   dot (v, pout->lmn) / C);
//    beta=(dot (v, pin->lmn) / C);
      //   ftemp=pin->freq*sqrt((1-beta)/(1+beta));
      //  beta=(dot (v, pout->lmn) / C);
      // pout->freq=ftemp/sqrt((1-beta)/(1+beta));


    }
  else if (nres > -1 && nres < nlines)
    {				/* It was a resonant scatter. */
      pout->freq = lin_ptr[nres]->freq / (1. - dot (v, pout->lmn) / C);
    }
  else if ((nres > NLINES && nres < NLINES + ntop_phot + 1) || nres == -2)
    /* It was continuum emission - new comoving frequency has been chosen by
       the matom/kpkt routine, but now need to convert in the same way 
       as for lines (SS) */
    {
      /* 
         If a 2-level atom run, one should never arrive here.
         Just do a check that all is well - this can be removed eventually (SS)
       */
      if (geo.rt_mode == 1)
	{
	  Error
	    ("doppler: Not using macro atoms but trying to deexcite one? Abort.\n");
	  exit (0);
	}
      pout->freq = pout->freq / (1. - dot (v, pout->lmn) / C);
    }
/* Now do one final check that nothing is awry.  This is another
 * check added by SS that should probably be deleted or done before this point.  
 * I have made it fatal so that we will pay attention to it if it occurs. ksl */

  else
    {
      Error ("doppler: nres %d > NLINES+ntop_phot %d\n", nres,
	     NLINES + ntop_phot);
      exit (0);
    }

  return (0);

}

/***********************************************************
                                       Space Telescope Science Institute

 Synopsis:

	int scatter(w,p,nres) determines a new direction and frequency for a photon 
	that scatters in the wind 

Arguments:

	WindPtr w
	PhotPtr p			the  photon of interest
	int nres;			either the number of the scatter, or a nonresonant scatter
						if nres < 0

Returns:

	The results are stored in the PhotPtr p which contains the direction and frequency
	of the scattered photon.	
	
Description:	 

	The routine calculates a new direction and frequency for a photon in both the
	resonant and non-resonant cases.  In the frame of the wind, scattering is assumed
	to be isotropic.  

Notes:
	This is the routine that is called when a resonant scatter does occur.
	
	The equations for the frequency shifts are accurate only to first order in beta

	This routine should not move the photon at all, because other routines need to
	take this photon in differing directions, and if one moves it here they may
	encounter this resonance again
	
		
History:
	96dec	ksl	Coded as part of monte effort
	97jan	ksl 	updated for python
	02jan	ksl	Completely rewritten.  Current version (and previous version)
			is first order in v/c and assumes that the distribution is isotropic
			in the Eulerian frame.
	02jan	ksl	Added update of resonance number
	02feb	ksl	Included an approach for anisotropic scattering
        04feb   SS      Substantially modified to incorporate Macro Atom
                        treatment
        04Mar   SS      Creating of k-packets following bf absn is now allowed.
        04Apr   SS      Modified to branch when simple ions are involved 
                        (i.e. those for which a full macro atom treatment is
                        not required)
                        Also ff is now included: nres = -2 flags ff which leads
                        to creating a k-packet always.
        04May   SS      Modifications to allow the macro_simple option (i.e. all
                        ions treated as "simple").
        04Jun   SS      Modifications so simplify the interation of macro atoms and
                        k-packets. This routine does not now access kpkt and matom 
                        directly but via macro_gov which supervises what the k-packets
                        and r-packets do. Hopefully this make the logic easier to follow.

        04Dec  SS     Added nnscat pointer to call so that the anisotropic thermal
                        scattering model would work as before for "simple" calculations.
                        Previously nnscat was always just = 1.

        1406 	JM 		Added normalisation of rejection method for anisotropic scattering
        				'thermal trapping' model.
        				See Issue #82.




***********************************************************/

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
  PlasmaPtr xplasma;
  MacroPtr mplasma;



  stuff_phot (p, &pold);
  n = where_in_grid (pold.x);	// Find out where we are

  //71 - 1112 Check added to test out spherical wind models 
  if (n < 0)
    {
      Error ("scatter: Trying to scatter a photon in grid cell %d\n", n);
      return (-1);
    }

  one = &wmain[p->grid];
  xplasma = &plasmamain[one->nplasma];
  //OLD - did not trap a problem if (xplasma==NULL){
  //OLD - did not trap a problem          Error("Houston, we have a null pointer at %d %d",p->grid,one->nplasma);
  //OLD - did not trap a problem }

  /* On entering this subroutine we know that a photon packet has been 
     absorbed. nres tells us which process absorbed it. There are currently
     four "absorption" processes: electron scattering (flagged -1) line
     absorption (flagged by +ve integer < NLINES), bf absorption 
     (flagged by +ve integer > NLINES) and ff absorption (flagged -2). (SS) */

  /* If the macro atom method is being used then the following section must be
     performed to select a macro atom deactivation process. If not then the
     deactivation process is always the same as the activation process and so
     nothing needs to be done. */

  if (geo.rt_mode == 2)		//check if macro atom method in use
    {

      /* 1112 - 71 - ksl - Moved to avoid trying to reference mplasma if there are no 
         macro atoms.  This was to fix a segmentation fault that appeared
         when compiling with a new version of gcc.   It's possible that the error below
         could occur if we were in a macro atom approach but had no macro atoms.  Need
         to fix all this up with a thorough review of macro atoms. !!!
       */

      if (geo.nmacro > 0)
	{
	  mplasma = &macromain[xplasma->nplasma];
	}
      else
	{
	  mplasma = NULL;
	  Error
	    ("Resonate: In macro atom section, but no macro atoms.  Seems very odd\n");
	}

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

	  if (phot_top[*nres - NLINES - 1].macro_info == 1
	      && geo.macro_simple == 0)
	    {
	      /* Macro ion case (SS) */

	      // Note:  NLINES-1 in the lines below is correct.  This is becasue
	      // the 1st bf is identified by nres = NLINES+1 and this is 
	      // the zeroth element of phot_top: hence the -1.  SS

	      llvl = phot_top[*nres - NLINES - 1].nlev;	//lower level
	      ulvl = phot_top[*nres - NLINES - 1].uplev;	//upper level

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
		  Error
		    ("scatter (resonate.c): could not identify bf transition. Abort. \n");
		  exit (0);
		}

	      // Need to compute the factor needed for the stimulated term.

	      stim_fact =
		den_config (xplasma, ulvl) / den_config (xplasma,
							 llvl) / xplasma->ne;

	      gamma_twiddle =
		mplasma->gamma_old[config[llvl].bfu_indx_first + m] -
		(mplasma->alpha_st_old[config[llvl].bfu_indx_first + m] *
		 stim_fact);
	      gamma_twiddle_e =
		mplasma->gamma_e_old[config[llvl].bfu_indx_first + m] -
		(mplasma->alpha_st_e_old[config[llvl].bfu_indx_first + m] *
		 stim_fact);

	      // Both gamma_twiddles must be greater that zero if this is going to work. If they 
	      // are zero then it's probably because this is the first iteration and so the've not
	      //been computed yet. For that first iteration k-packets will be ignored. If the
	      // gamma_twiddles are negative then something has gone wrong.

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
		  Error
		    ("scatter (resonate.c): a gamma_twiddle is negative. Abort.\n");
		  exit (0);
		}

	      /* Having got here we have calculated the probability of a k-packet
	         being created. Now either make a k-packet or excite a macro atom. */

	      kpkt_choice = ((rand () + 0.5) / MAXRAND);	//random number for kpkt choice

	      if (prob_kpkt > kpkt_choice)
		{
		  macro_gov (p, nres, 2, &which_out);	//routine to deal with kpkt
		}
	      else
		{
		  macro_gov (p, nres, 1, &which_out);	//routine to deal with macro atom excitation
		}
	    }
	  else if (phot_top[*nres - NLINES - 1].macro_info == 0
		   || geo.macro_simple == 1)
	    {
	      // Simple ion case //
	      /* Need to make decision about making a k-packet. Get the fraction of the energy
	         that goes into the electron rather than being stored as ionisation energy: this
	         fraction gives the selection probability for the packet. It's given by the
	         (photon frequency / edge frequency - 1) (SS) */


	      prob_kpkt =
		1. - (phot_top[*nres - NLINES - 1].freq[0] / p->freq);


	      /* Now choose whether or not to make a k-packet. */

	      kpkt_choice = ((rand () + 0.5) / MAXRAND);	//random number for kpkt choice

	      if (prob_kpkt > kpkt_choice)
		{
		  macro_gov (p, nres, 2, &which_out);	//routine to deal with kpkt
		}
	      else
		{
		  macro_gov (p, nres, 1, &which_out);	//routine to deal with fake macro atom bf excitation

		}
	    }
	  else
	    {
	      /* Our best-laid schemes have gang agley. It should never get here unless the input has been
	         messed up in some way. (SS) */
	      Error
		("scatter (resonate.c): continuum scatter - seems to be neither macro nor simple. Abort.\n");
	      exit (0);
	    }
	}
      else if (*nres == -2)
	{			/* This is a ff event (SS). */
	  macro_gov (p, nres, 2, &which_out);	//ff always make a k-packet
	}
    }

  /* So at this point we have completed all the bits that are specific to the macro approach. 54b--ksl */

  /* SS Apr 04: I've moved this next block up so that nres is set correctly for the call to randvec */

  /* Since the error check is commented out for good reason, we should just assign 
   * *nres to p->nres,  and be done with it.  ?? Stuart, assuming you agree just
   * eliminate all the clutter here.  KSL  ??
   */
  p->nres = *nres;		// Update the resonance number on a scatter

  /* SS July 04
     Next block is modified to include the thermal trapping model for anisotropic scattering.
     The code for this has been moved from trans_phot to here so that this model can work 
     with macro atoms.
     For macro atoms the code above decides that emission will occur in the line - we now just need
     to use the thermal trapping model to choose the direction. */


  if (*nres == -1 || *nres == -2 || *nres > NLINES || geo.scatter_mode == 0
      || geo.scatter_mode > 2)
    //geo.scatter_mode > 2 should never happen but to keep consistency with what was here before I've
    //added it as a possibility  SS.  ?? I'm not sure about the geo.scatter mode > 2 stuff.  Seems
    // as if we should put an error check at the end and bail.  ksl 04dec.
    {
      /*  It was either an electron scatter, bf emission or ff emission so the  distribution is isotropic, 
         or it was a line photon but we want isotropic scattering anyway.  */
      randvec (z_prime, 1.0);	/* Get a new direction for the photon */
      stuff_v (z_prime, p->lmn);
    }

  else if (geo.scatter_mode == 1)
    {				// It was a line photon and we want anisotropic scattering model 1

      randwind (p, z_prime, wmain[n].lmn);
      stuff_v (z_prime, p->lmn);

    }
  else
    {				//It was a line photon and we want to use the thermal trapping model to choose the output direction

      /* JM 1906 -- added normalisation of the below rejection method. We normalise
         to the escape probability of along the direction of dvds_max, with a safety net of 
         20% in case we missed the maximum */
      randwind_thermal_trapping(p, nnscat);
    }

  /* End of modification for thermal trapping model (SS July 04) */



  //stuff_v (z_prime, p->lmn);

  vwind_xyz (p, v);		/* Get the velocity vector for the wind */


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
