#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "atomic.h"
#include "python.h"


/* A couple of external things for use in the routines for computing gamma's below. */
struct topbase_phot *cont_ext_ptr2;	//continuum pointer passed externally
double temp_ext2;		//temperature passed externally
double temp_ext_rad;		//radiation temperature passed externally 

/************************************************************
                                    Imperial College London
Synopsis:
	bf_estimators_increment is the routine to increment the bound-free
	estimators as needed for the macro atom calculations. Whenever a packet
	travels a path length ds through a shell the estimator for all bound-free
	processes where the packet frequency is above the threshold frequency
	is incremented. The estimator is not normalised here. That is done later.

Arguments:

       WindPtr one                 pointer to cell
       PhotPtr p                   the packet
       double ds                   the path length

Returns:
       

Description:


Notes: The estimators are not normalised by this routine.
       This routine also now computes the total heating rate due
       to simple ions.  Therefore the routine needs to be called
       even for the simple non-macro atom case.

History:
          Mar 04  SS   Coding began.
          Apr 04  SS   Computation of heating rate due to simple ions added.
          Apr 04  SS   Computation of heating rate due to ff added.
          May 04  SS   Added lines for macro_simple option (i.e. treat all ions
                       as simple)
          Sep 04  SS   Modified to record energy absorped by macro atom levels and k-packets.
	06may	ksl	57+ -- Replaced wind with plasma structure, mainly.  Note that
			one is already assigned here, and so I did not switch everything
			but it may be that this should be done
************************************************************/

int
bf_estimators_increment (one, p, ds)
     WindPtr one;
     PhotPtr p;
     double ds;

{
  double freq_av;
  double weight_of_packet;
  double x, ft;
  double y, yy;
  double exponential, heat_contribution;
  int n, m, llvl, nn;
  double sigma_phot_topbase ();
  double density;
  double abs_cont;
  int nplasma;
  PlasmaPtr xplasma;
  MacroPtr mplasma;


  nplasma = one->nplasma;
  xplasma = &plasmamain[nplasma];
  mplasma = &macromain[xplasma->nplasma];


  freq_av = p->freq;
  // the continuum neglect variation of frequency along path and
  // take as a single "average" value.  

  if (p->freq > xplasma->max_freq)  // check if photon frequency exceeds maximum frequency
    xplasma->max_freq = p->freq;


  /* JM -- 1310 -- check if the user requires extra diagnostics and
     has provided a file diag_cells.dat to store photons stats for cells they have specified
   */
  if (diag_on_off == 1 && ncstat > 0)
    {
      save_photon_stats(one, p, ds);  // save photon statistics (extra diagnostics)
    }


  
  /* JM 1402 -- the banded versions of j, ave_freq etc. are now updated in update_banded_estimators,
     which also updates the ionization parameters and scattered and direct contributions */

  update_banded_estimators(xplasma, p, ds, p->w);

  /* check that j and ave freq give sensible numbers */
  if (sane_check (xplasma->j) || sane_check (xplasma->ave_freq))
    {
      Error ("radiation:sane_check Problem with j %g or ave_freq %g\n",
	     xplasma->j, xplasma->ave_freq);
    }




  for (nn = 0; nn < xplasma->kbf_nuse; nn++)
    {
      n = xplasma->kbf_use[nn];
      ft = phot_top[n].freq[0];	//This is the edge frequency (SS)

      llvl = phot_top[n].nlev;	//Returning lower level = correct (SS)

      density = den_config (xplasma, llvl);

      /* JM130729 Bugfix 31: This if loop causes the else statement for simple ions to be 
       * entered in macro atom mode- it appeared to be introduced sometime between 58 and 68.
       *
       * if (kap_bf[nn] > 0.0 && (freq_av > ft) && phot_top[n].macro_info == 1
       *          && geo.macro_simple == 0)
       */
      if (kap_bf[nn] > 0.0 && (freq_av > ft))	// does the photon cause bf heating?
	{

	  if (phot_top[n].macro_info == 1 && geo.macro_simple == 0)	// it is a macro atom
	    {
	      x = kap_bf[nn] / density;	//this is the cross section

	      /* Now identify which of the BF processes from this level this is. */

	      m = 0;
	      while (m < config[llvl].n_bfu_jump
		     && config[llvl].bfu_jump[m] != n)
		m++;

	      // m should now be the label to identify which of the bf processes from llvl
	      // this is. Check that it is reasonable

	      if (m > config[llvl].n_bfu_jump - 1)
		{
		  Error
		    ("bf_estimators_increment: could not identify bf transition. Abort. \n");
		  exit (0);
		}

	      // Now calculate the contributions and add them on.
	      weight_of_packet = p->w;
	      y = weight_of_packet * x * ds;
	      exponential =
		y * exp (-(freq_av - ft) / BOLTZMANN / xplasma->t_e);
	      mplasma->gamma[config[llvl].bfu_indx_first + m] += y / freq_av;
	      mplasma->alpha_st[config[llvl].bfu_indx_first + m] +=
		exponential / freq_av;
	      mplasma->gamma_e[config[llvl].bfu_indx_first + m] += y / ft;
	      mplasma->alpha_st_e[config[llvl].bfu_indx_first + m] +=
		exponential / ft;

	      /* Now record the contribution to the energy absorbed by macro atoms. */
	      yy = y * den_config (xplasma, llvl);
	      mplasma->matom_abs[phot_top[n].uplev] += abs_cont =
		yy * ft / freq_av;
	      xplasma->kpkt_abs += yy - abs_cont;
	      /* the following is just a check that flags packets that appear to travel a 
	         suspiciously large optical depth in the continuum */
	      if ((yy / weight_of_packet) > 50)
		{
		  Log
		    ("bf_estimator_increment: A packet survived an optical depth of %g\n",
		     yy / weight_of_packet);
		  Log ("bf_estimator_increment: freq_av %g, ft %g\n", freq_av,
		       ft);
		}
	    }

	  else			// it is a simple ion
	    {
	      /* Now we are dealing with the heating due to the bf continua of simple ions. No stimulated
	         recombination is included here. (SS, Apr 04) */
	      if (density > DENSITY_PHOT_MIN)
		{
		  x = sigma_phot_topbase (&phot_top[n], freq_av);	//this is the cross section
		  weight_of_packet = p->w;
		  y = weight_of_packet * x * ds;

		  xplasma->heat_photo += heat_contribution =
		    y * density * (1.0 - (ft / freq_av));
		  xplasma->heat_tot += heat_contribution;
		  /* This heat contribution is also the contibution to making k-packets in this volume. So we record it. */

		  xplasma->kpkt_abs += heat_contribution;
		}
	    }
	}
    }

  /* Now for contribution to heating due to ff processes. (SS, Apr 04) */

  weight_of_packet = p->w;

  y = weight_of_packet * kappa_ff (xplasma, freq_av) * ds;

  xplasma->heat_ff += heat_contribution = y;	// record ff hea        


  /* Now for contribution to heating due to compton processes. (JM, Sep 013) */

  y = weight_of_packet * kappa_comp (xplasma, freq_av) * ds;

  xplasma->heat_comp += y;	// record the compton heating
  heat_contribution += y;	// add compton to the heat contribution


  /* Now for contribution to heating due to induced compton processes. (JM, Sep 013) */

  y = weight_of_packet * kappa_ind_comp (xplasma, freq_av) * ds;

  xplasma->heat_ind_comp += y;	// record the induced compton heating
  heat_contribution += y;	// add induced compton to the heat contribution



  xplasma->heat_tot += heat_contribution;	// heat contribution is the contribution from compton, ind comp and ff processes



  /* This heat contribution is also the contibution to making k-packets in this volume. So we record it. */

  xplasma->kpkt_abs += heat_contribution;


  /* Now for contribution to inner shell ionization estimators (SS, Dec 08) */

  for (n = 0; n < nauger; n++)
    {
      ft = augerion[n].freq_t;
      if (freq_av > ft)
	{
	  Log ("estimators: Adding a packet to AUGER %g \n", freq_av);

	  weight_of_packet = p->w;
	  x = sigma_phot_verner (&augerion[n], freq_av);	//this is the cross section
	  y = weight_of_packet * x * ds;

	  xplasma->gamma_inshl[n] += y / freq_av / H / one->vol;
	}
    }



  return (0);
  /* All done. (SS) */
}



/************************************************************
                                    Imperial College London
Synopsis:
	bb_estimators_increment is the routine to increment the bound-bound
	estimators as needed for the macro atom calculations. Whenever a packet
	comes into resonance with a line this routine is called WHETHER OR NOT
	the packet interacts with the line.
	The estimator is not normalised here. That is done later.

Arguments:

       WindPtr one                 pointer to cell
       PhotPtr p                   the packet
       double tau_sobolev          optical depth of line
       double dvds                 velocity gradient
       int nn                      the label for the line in question

Returns:
	0 on success
        1 if the line was a "simple" line.  The routine should not have
	  been called if that was the case.
       

Description:


Notes: The estimators are not normalised by this routine
       For the moment, the increment to the estimator is taken to be
       weight * (1 - exp(-tau))/ tau / dvds
       This is slightly different from what I do in my code (not
       accurate to v/c corrections). If calculations where accuracy to
       order v/c is needed then this routine should be improved.

History:
          Mar 04  SS   Coding began.
	04apr	ksl	Modified to check if line was a macro-line
			Fixed error calls to refer to routine and not
			filename.  
          May 04  SS  modifications to allow macro_simple option
          May 04  SS  modifications to cope with input line list that is not ordered
                      (use of where_in_list)
          Sep 04  SS  modified to record rate of energy absorbed by macro atom levels
                      ans k-packets
	06may	ksl	57+ -- Modifications to accommodate plasma structue
          
************************************************************/



int
bb_estimators_increment (one, p, tau_sobolev, dvds, nn)
     WindPtr one;
     PhotPtr p;
     double tau_sobolev;
     double dvds;
     int nn;

{
  int llvl;
  int n;
  struct lines *line_ptr;
  double weight_of_packet;
  double y;
  int nmax;

  PlasmaPtr xplasma;
  MacroPtr mplasma;

  xplasma = &plasmamain[one->nplasma];
  mplasma = &macromain[xplasma->nplasma];


  /* 04apr ksl: Start by checking that this was a macro-line */
  /* There should now be a diversion in place before this routine is called
     to identify macro/simple ions. But I'll leave the check here too - backup (SS, Apr04). */

  line_ptr = lin_ptr[nn];

  if (line_ptr->macro_info == 0 || geo.macro_simple == 1)
    return (-1);

  /* Start by identifying which estimator we want. Get lower level of
     transition . */

  llvl = line_ptr->nconfigl;

  /* Now find which of the bb transition from level llvl this is. */


  nmax = config[llvl].n_bbu_jump;

  n = 0;
  while (n < nmax && line[config[llvl].bbu_jump[n]].where_in_list != nn)
    n++;

  if (n == nmax)
    {
      Error
	("bb_estimators_increment: could not identify bb transition. Abort. \n");
      exit (0);
    }




  /* Okay now know which estimator we wish to increment so do it. */

  weight_of_packet = p->w;
  dvds = fabs (dvds);		//make sure that it is positive


  if (tau_sobolev > 0.00001)
    {
      y = weight_of_packet * (1. - exp (-tau_sobolev)) / tau_sobolev / dvds;
    }
  else				//To avoid tau_sobolev = 0
    {
      y = weight_of_packet / dvds;
    }

  if (y >= 0)
    {
      mplasma->jbar[config[llvl].bbu_indx_first + n] += y;
    }
  else
    {
      Error
	("bb_estimators_increment: trying to add negative contribution to jbar. Abort. \n");
      exit (0);
    }

  /* Record contribution to energy absorbed by macro atoms. */

  mplasma->matom_abs[line_ptr->nconfigu] +=
    weight_of_packet * (1. - exp (-tau_sobolev));

  return (0);
  /* All done. (SS) */
}

/************************************************************
                                    Imperial College London
Synopsis:
	mc_estimator_normalise (n) is the routine to normalise the 
	bb and bf estimators. It is performed as part of the wind update stage
	of Python. The estimators should only be changed during the ionisation
	cycles - after that they should be fixed. 

Arguments:

       w            WindPtr
       n            the cell for which the normalisation is to be done

Returns:
       

Description:


Notes: This routine normalises the bb and bf mc estimators needed
       for the macro atom jumping probabilities. During the main mc
       simulation these were stored unnormalised. They are now
       normalised and moved into the "old" slots for use in the next
       iteration. The current slots are reset to zero for use in the
       next iteration.

       It also now computes the bf heating rate using the normalised
       estimators. heat_tot and heat_photo are incremented but nothing else
       is for now.

       ksl -- This routine loops over nlte_levels, which in principle
       could include non-macro ions, but that should not matter 
       since nbfu_jumps will be zero for these.

History:
          Mar 04  SS   Coding began.
          Mar 04  SS   Added steps for computation of bf heating rate.
          Jun 04  SS   Added compution of collisional bf heating rate (3 body recombination)
	06may	ksl	57+ -- Changes to allow for plaama structure. Eliminate passage
			of entire wind array
          
************************************************************/

int
mc_estimator_normalise (n)
     int n;

{
  double volume;
  int i, j;
  double stimfac, line_freq, stat_weight_ratio;
  double heat_contribution;
  WindPtr one;
  PlasmaPtr xplasma;
  MacroPtr mplasma;

  one = &wmain[n];
  xplasma = &plasmamain[one->nplasma];
  mplasma = &macromain[xplasma->nplasma];

  /* All the estimators need the volume so get that first. */
  volume = one->vol;

  /* bf estimators. During the mc calculation the quantity stored
     was weight * cross-section * path-length / frequency.
     To normalise this we need to put in:
     1 / h  (Planck constant)
     1 / Volume
     1 / Time 
     I think that the weights are chosen such that Time = 1. 
     So the normalisation of gamma and gamma_e is easy. 
   */


  /* For the alpha_st estimators (stimulated recombination) the
     estimator also needs to be multiplied by a temperature
     dependent factor which is almost, but not quite, given by
     the LTE population ratio. The muliplicative factor
     is given by: */
  stimfac =
    0.5 * pow (H * H / 2. / PI / MELEC / BOLTZMANN / xplasma->t_e, 3. / 2.);

  for (i = 0; i < nlte_levels; i++)
    {
      for (j = 0; j < config[i].n_bfu_jump; j++)
	{

	  mplasma->gamma_old[config[i].bfu_indx_first + j] = mplasma->gamma[config[i].bfu_indx_first + j] / H / volume;	//normalise
	  mplasma->gamma[config[i].bfu_indx_first + j] = 0.0;	//re-initialise for next iteration
	  mplasma->gamma_e_old[config[i].bfu_indx_first + j] = mplasma->gamma_e[config[i].bfu_indx_first + j] / H / volume;	//normalise
	  mplasma->gamma_e[config[i].bfu_indx_first + j] = 0.0;	//re-initialise for next iteration

	  /* For the stimulated recombination parts we need the the
	     ratio of statistical weights too. 
	     For free electron statistical weight = 2 is included in
	     stimfac above. */

	  stat_weight_ratio =
	    config[phot_top[config[i].bfu_jump[j]].uplev].g / config[i].g;

	  mplasma->alpha_st_old[config[i].bfu_indx_first + j] =
	    mplasma->alpha_st[config[i].bfu_indx_first +
			      j] * stimfac * stat_weight_ratio / H / volume;
	  mplasma->alpha_st[config[i].bfu_indx_first + j] = 0.0;

	  mplasma->alpha_st_e_old[config[i].bfu_indx_first + j] =
	    mplasma->alpha_st_e[config[i].bfu_indx_first +
				j] * stimfac * stat_weight_ratio / H / volume;
	  mplasma->alpha_st_e[config[i].bfu_indx_first + j] = 0.0;

	  /* For continuua whose edges lie beyond freqmin assume that gamma
	     is given by a black body. */

	  /* For now place the limit at 7.5e12 which is 400000AA */
	  /* Try also doing it for very high energy ones - greater than 50eV: 1.2e16 since up there the statistics of the estimators are very poor at the moment. Ideally we don't want to have this so should probably switch this back sometime (SS August 05) !!!BUG */

	  if (phot_top[config[i].bfu_jump[j]].freq[0] < 7.5e12
	      || phot_top[config[i].bfu_jump[j]].freq[0] > 1.2e16)
	    {
	      mplasma->gamma_old[config[i].bfu_indx_first + j] =
		get_gamma (&phot_top[config[i].bfu_jump[j]], xplasma);
	      mplasma->gamma_e_old[config[i].bfu_indx_first + j] =
		get_gamma_e (&phot_top[config[i].bfu_jump[j]], xplasma);
	      mplasma->alpha_st_e_old[config[i].bfu_indx_first + j] =
		get_alpha_st_e (&phot_top[config[i].bfu_jump[j]], xplasma);
	      mplasma->alpha_st_old[config[i].bfu_indx_first + j] =
		get_alpha_st (&phot_top[config[i].bfu_jump[j]], xplasma);
	    }

	}

      /* That deals with the bf jumps. Now need to sort out the bb jumps. */

      /* For bb jumps the normalisation requires:
         1/ 4 PI
         1/ Volume
         1/ Time
         c / nu  - to convert dvds^-1 to dnuds^-1

         Also, I'm putting the correction factor for stimulated emission in here - so that it doesn't have to be
         computed in the macro atom jumping probabilities. (SS)


       */

      for (j = 0; j < config[i].n_bbu_jump; j++)
	{

	  /* The correction for stimulated emission is (1 - n_lower * g_upper / n_upper / g_lower) */

	  stimfac =
	    den_config (xplasma,
			line[config[i].bbu_jump[j]].nconfigu) /
	    den_config (xplasma, i);
	  stimfac =
	    stimfac * config[i].g /
	    config[line[config[i].bbu_jump[j]].nconfigu].g;
	  if (stimfac < 1.0)
	    {
	      stimfac = 1. - stimfac;	//all's well
	    }
	  else
	    {
	      Error
		("mc_estimator_normalise: bb stimulated correction factor is out of bound. Abort.\n");
	      Error
		("stimfac %g, i %d, line[config[i].bbu_jump[j]].nconfigu %d\n",
		 stimfac, i, line[config[i].bbu_jump[j]].nconfigu);
	      Log
		("estimators: den_config (xplasma, i) %g  den_config (xplasma, line[config[i].bbu_jump[j]].nconfigu) %g \n",
		 den_config (xplasma, i), den_config (xplasma,
						      line[config[i].bbu_jump
							   [j]].nconfigu));
	      stimfac = 0.0;
	      //exit (0);
	    }

	  //get the line frequency
	  line_freq = line[config[i].bbu_jump[j]].freq;
	  mplasma->jbar_old[config[i].bbu_indx_first + j] =
	    mplasma->jbar[config[i].bbu_indx_first +
			  j] * C * stimfac / 4. / PI / volume / line_freq;

	  mplasma->jbar[config[i].bbu_indx_first + j] = 0.0;
	}
    }

  /* bb and bf now normalised. Done. */
  /* Get the heating contribution from macro atom bb transitions (the
     line heating). */

  xplasma->heat_lines += heat_contribution =
    macro_bb_heating (xplasma, xplasma->t_e);
  xplasma->heat_lines_macro = heat_contribution;
  xplasma->heat_tot += heat_contribution;

  /* Get the bf heating contributions here too. (SS June 04) */

  xplasma->heat_photo += heat_contribution =
    macro_bf_heating (xplasma, xplasma->t_e);
  xplasma->heat_photo_macro = heat_contribution;
  xplasma->heat_tot += heat_contribution;


  /* Now that we have estimators, set the plag to use them for the level populations */


  geo.macro_ioniz_mode = 1;


  return (0);
}


/************************************************************
                                    Imperial College London
Synopsis: total_fb_matoms computes the cooling rate due to free-bound
	recombinations using mc estimators. It is modelled on total_fb
	but makes use of the mc estimators for stimulated recombination.

Arguments:

       ww            WindPtr
       t_e           electron temperature
       f1            lower frequency
       f2            upper frequency

Returns:
       

Description:


Notes: This returns the cooling rate for bf recombination.

       ksl -- This routine loops over nlte_levels, which in principle
       could include non-macro ions, but that should not matter since
       since nbfu_jumps will be zero for these.

History:
          Mar 04  SS   Coding began.
          Jun 04  SS   Modified to include the cooling due to
                       collisional ionization.
	04Jul	ksl	Modified calls to alpha_sp to reflect
			fact that alpha_sp has been generalized
			to cover the energy averaged and 
			spontaneous case.
	06may	ksl	57+ -- Modified to use plasma structure
			Uses volume so did not delete Wind structure
			but if we put volume in both places could do this

          
************************************************************/

double
total_fb_matoms (xplasma, t_e, f1, f2)
     PlasmaPtr xplasma;
     double t_e;
     double f1, f2;
{
  double cool_contribution;
  double t_e_store;
  struct topbase_phot *cont_ptr;
  double total, density;
  int i, j;
  double q_ioniz ();
  MacroPtr mplasma;

  mplasma = &macromain[xplasma->nplasma];

  t_e_store = xplasma->t_e;	//store the temperature - will put it back at the end
  xplasma->t_e = t_e;		//for use in calls to alpha_sp below

  total = 0;			// initialise

  if (geo.macro_simple == 0)	//allow for "only-simple" calculations (SS May04)
    {
      for (i = 0; i < nlte_levels; i++)
	{
	  for (j = 0; j < config[i].n_bfu_jump; j++)
	    {
	      /* Need the density for the upper level in the recombination
	         process. */
	      cont_ptr = &phot_top[config[i].bfu_jump[j]];
	      density = den_config (xplasma, cont_ptr->uplev);
	      cool_contribution =
		(mplasma->alpha_st_e_old[config[i].bfu_indx_first + j] +
		 alpha_sp (cont_ptr, xplasma, 1)
		 - mplasma->alpha_st_old[config[i].bfu_indx_first + j]
		 - alpha_sp (cont_ptr, xplasma, 0))
		* H * phot_top[config[i].bfu_jump[j]].freq[0] * density *
		xplasma->ne * wmain[xplasma->nwind].vol;

	      /* Now add the collisional ionization term. */
	      density = den_config (xplasma, cont_ptr->nlev);
	      cool_contribution +=
		q_ioniz (cont_ptr,
			 t_e) * density * xplasma->ne * H *
		phot_top[config[i].bfu_jump[j]].freq[0] *
		wmain[xplasma->nwind].vol;

	      /* That's the bf cooling contribution. */
	      total += cool_contribution;
	    }
	}

      xplasma->t_e = t_e_store;	//restore the original value

    }

  return (total);

  //All done (SS).

}

/************************************************************
                                    Imperial College London
Synopsis:    
       total_bb_cooling computes the total cooling in bb transitions
       (due to collisions) for both macro atoms and simple ions. 
       It is used by the heating/cooling calculation to get the temperature.


Arguments:

       ww            WindPtr
       t_e           electron temperature


Returns:
       

Description:


Notes: This returns the total cooling rate for bb collisions.
       

History:
          Apr 04  SS   Coding began.
          May 04  SS   Corrections made to get energy radiated rather than number of cooling events.
          May 04  SS   Macro_simple option added (for all ions to be "simple")
	06may	ksl	57+ -- Adapted to plasma structure.  Calls probably need to modified
			Uses volume
 
          
************************************************************/

double
total_bb_cooling (xplasma, t_e)
     PlasmaPtr xplasma;
     double t_e;
{
  double cool_contribution;
  struct lines *line_ptr;
  double total, upper_density, lower_density;
  int i;
  double coll_rate, rad_rate;

  total = 0;			// initialise
  for (i = 0; i < nlines; i++)
    {
      line_ptr = &line[i];
      if (line_ptr->macro_info == 1 && geo.macro_simple == 0)
	{			//This is a line from a macro atom for which we know
	  //the upper and lower level populations
	  lower_density = den_config (xplasma, line_ptr->nconfigl);
	  cool_contribution =
	    (lower_density * q12 (line_ptr, t_e)) * xplasma->ne *
	    wmain[xplasma->nwind].vol * line_ptr->freq * H;
	}
      else
	{			//It's a simple line - don't know the level populations
	  // - just use a two-level-atom approach

	  //The cooling rate is computed using the scattering probability formalism in KSL's notes on Python.

	  two_level_atom (line_ptr, xplasma, &lower_density, &upper_density);
	  coll_rate =
	    q21 (line_ptr,
		 t_e) * xplasma->ne * (1. -
				       exp (-H_OVER_K * line_ptr->freq /
					    t_e));

	  cool_contribution =
	    (lower_density * line_ptr->gu / line_ptr->gl -
	     upper_density) * coll_rate / (exp (H_OVER_K * line_ptr->freq /
						t_e) -
					   1.) * wmain[xplasma->nwind].vol *
	    line_ptr->freq * H;


	  rad_rate = a21 (line_ptr) * p_escape (line_ptr, xplasma);
	  cool_contribution *= rad_rate / (rad_rate + coll_rate);
	}

      /* That's the bb cooling contribution. */

      total += cool_contribution;

    }

  return (total);

  //All done (SS).

}

/************************************************************
                                    Imperial College London
Synopsis:    
       macro_bb_heating computes the total heating due to bb transitions
       (due to collisions) for macro atoms. The heating in simple ions
       is taken care of elsewhere. 
       It is used by the heating/cooling calculation to get the temperature.


Arguments:

       ww            WindPtr
       t_e           electron temperature


Returns:
       

Description:


Notes: This returns the heating rate for bb collisions in macro atoms.
       

History:
          Apr 04  SS   Coding began.
          May 04  SS   Corrections made to get energy rather than number of heating events.
          May 04  SS   Macro_simple option added (for all ions to be "simple")
	06may	ksl	57+ -- Modified for plasma structue.  Note that volume is used
          
************************************************************/

double
macro_bb_heating (xplasma, t_e)
     PlasmaPtr xplasma;
     double t_e;
{
  double heat_contribution;
  struct lines *line_ptr;
  double total, upper_density;
  int i;


  total = 0;			// initialise

  for (i = 0; i < nlines; i++)
    {
      line_ptr = &line[i];
      if (line_ptr->macro_info == 1 && geo.macro_simple == 0)
	{			//This is a line from a macro atom for which we know
	  //the upper and lower level populations
	  upper_density = den_config (xplasma, line_ptr->nconfigu);
	  heat_contribution =
	    upper_density * q21 (line_ptr,
				 t_e) * xplasma->ne *
	    wmain[xplasma->nwind].vol * line_ptr->freq * H;
	  total += heat_contribution;
	}
    }

  return (total);

  //All done (SS).

}

/************************************************************
                                    Imperial College London
Synopsis:    
       macro_bf_heating computes the total heating due to bf transitions
       for macro atoms. The heating in simple ions
       is taken care of elsewhere. 
       It is used by the heating/cooling calculation to get the temperature.


Arguments:

       ww            WindPtr
       t_e           electron temperature


Returns:
       

Description:


Notes: This returns the heating rate for bf transitions in macro atoms.
       

History:
         Jun 04 - SS coding began: previously the computation of the 
                     bf heating was done in the normalisation of the estimators
                     but with the inclusion of three body recombination it makes
                     more sense to put it all in one subroutine.
	06may	ksl	57+ -- Modified for plama structue
          
************************************************************/

double
macro_bf_heating (xplasma, t_e)
     PlasmaPtr xplasma;
     double t_e;
{
  double heat_contribution;
  double total, upper_density, lower_density;
  int i, j;
  double q_recomb ();
  MacroPtr mplasma;

  mplasma = &macromain[xplasma->nplasma];

  total = 0;			// initialise

  for (i = 0; i < nlte_levels; i++)
    {
      for (j = 0; j < config[i].n_bfu_jump; j++)
	{
	  heat_contribution = 0.0;
	  /* Photoionization part. */
	  lower_density =
	    den_config (xplasma, phot_top[config[i].bfu_jump[j]].nlev);
	  heat_contribution =
	    (mplasma->gamma_e_old[config[i].bfu_indx_first + j] -
	     mplasma->gamma_old[config[i].bfu_indx_first + j]) * H *
	    phot_top[config[i].bfu_jump[j]].freq[0] * lower_density *
	    wmain[xplasma->nwind].vol;

	  /* Three body recombination part. */
	  upper_density =
	    den_config (xplasma, phot_top[config[i].bfu_jump[j]].uplev);
	  heat_contribution +=
	    q_recomb (&phot_top[config[i].bfu_jump[j]],
		      t_e) * xplasma->ne * xplasma->ne * H * upper_density *
	    wmain[xplasma->nwind].vol *
	    phot_top[config[i].bfu_jump[j]].freq[0];

	  total += heat_contribution;

	}
    }

  return (total);

  //All done (SS).

}

/************************************************************
                                    Imperial College London
Synopsis:
	bb_simple_heat records the heating contribution from lines of
        simple elements. It is called whenever a packet comes into 
        resonance with the line and records the heating contribution 
        from that line and packet

Arguments:

       WindPtr one                 pointer to cell
       PhotPtr p                   the packet
       double tau_sobolev          optical depth of line
       double dvds                 velocity gradient
       int nn                      the label for the line in question
     
Returns:
	0 on success
       

Description:


Notes: 


History:
       04 Apr  SS: coding began
       04 May  SS: major re-write (it wasn't correct before)
       04 Nov  SS: modified to record the heating contribution as a means
                   of making k-pkts
	06may	ksl	57+ -- Modified for plasma structure.  There
			is no volume here, so have changed the entire
			routine to use the plasma structure

************************************************************/



int
bb_simple_heat (xplasma, p, tau_sobolev, dvds, nn)
     PlasmaPtr xplasma;
     PhotPtr p;
     double tau_sobolev;
     double dvds;
     int nn;

{
  double heat_contribution;
  double weight_of_packet;
  struct lines *line_ptr;
  double electron_temperature;
  double rad_rate, coll_rate, normalisation;
  double d1, d2;		//densities of lower and upper level
  double b12 ();

  /* The heating contribution is modelled on the macro atom bb estimator 
     calculations for the radiative excitation rate. This (energy) excitation
     rate is multiplied by the destruction probability to get the heating. 
     The destruction probability is obtained following the discussion in KSL's notes
     on Python. */


  weight_of_packet = p->w;
  line_ptr = lin_ptr[nn];
  electron_temperature = xplasma->t_e;
  two_level_atom (line_ptr, xplasma, &d1, &d2);	//get level densities

  rad_rate = a21 (line_ptr) * p_escape (line_ptr, xplasma);

  coll_rate = q21 (line_ptr, electron_temperature) * xplasma->ne
    * (1. - exp (-H_OVER_K * line_ptr->freq / electron_temperature));

  normalisation = rad_rate + coll_rate;


  /* Now add the heating contribution. */

  xplasma->heat_lines += heat_contribution =
    weight_of_packet * (coll_rate / normalisation) * (1. -
						      exp (-1. *
							   tau_sobolev));

  xplasma->heat_tot += heat_contribution;
  xplasma->kpkt_abs += heat_contribution;

  return (0);

}


/**************************************************
  get_gamma - to get the energy weighted photoionization rate
  for a black body radiation field with known dilution factor
  and temperature


	06may	ksl	57+ -- Changed to reflect plasma
			structure.  Changed call directly
			since volume is not involved.
*****************************************************/

double
get_gamma (cont_ptr, xplasma)
     struct topbase_phot *cont_ptr;
     PlasmaPtr xplasma;
{
  double gamma_value;
  double fthresh, flast;
  double qromb ();
  double gamma_integrand ();

  temp_ext2 = xplasma->t_r;	//external temperature
  cont_ext_ptr2 = cont_ptr;	//external cont pointer
  fthresh = cont_ptr->freq[0];	//first frequency in list
  flast = cont_ptr->freq[cont_ptr->np - 1];	//last frequency in list

  gamma_value = qromb (gamma_integrand, fthresh, flast, 1e-4);

  gamma_value *= 8 * PI / C / C * xplasma->w;

  return (gamma_value);

}

/********************************************************
 Function to give the integrand for gamma at frequency freq
**************************************************/

double
gamma_integrand (freq)
     double freq;
{
  double fthresh;
  double x;
  double integrand;
  double tt;

  fthresh = cont_ext_ptr2->freq[0];
  tt = temp_ext2;

  if (freq < fthresh)
    return (0.0);		// No photoionization at frequencies lower than the threshold freq occur

  x = sigma_phot_topbase (cont_ext_ptr2, freq);	//this is the cross-section
  integrand = x * freq * freq / (exp (H_OVER_K * freq / tt) - 1);

  return (integrand);
}

/*****************************************************************************/
/**************************************************
  get_gamma_e - to get the energy weighted photoionization rate
  for a black body radiation field with known dilution factor
  and temperature

	06may	ksl	57+ -- Changed to reflect plasma structure
			No need for volume so eleminated Wind altogther
*****************************************************/

double
get_gamma_e (cont_ptr, xplasma)
     struct topbase_phot *cont_ptr;
     PlasmaPtr xplasma;
{
  double gamma_e_value;
  double fthresh, flast;
  double qromb ();
  double gamma_e_integrand ();

  temp_ext2 = xplasma->t_r;	//external temperature
  cont_ext_ptr2 = cont_ptr;	//external cont pointer
  fthresh = cont_ptr->freq[0];	//first frequency in list
  flast = cont_ptr->freq[cont_ptr->np - 1];	//last frequency in list

  gamma_e_value = qromb (gamma_e_integrand, fthresh, flast, 1e-4);

  gamma_e_value *= 8 * PI / C / C * xplasma->w;

  return (gamma_e_value);

}

/********************************************************
 Function to give the integrand for gamma_e at frequency freq
**************************************************/

double
gamma_e_integrand (freq)
     double freq;
{
  double fthresh;
  double x;
  double integrand;
  double tt;

  fthresh = cont_ext_ptr2->freq[0];
  tt = temp_ext2;

  if (freq < fthresh)
    return (0.0);		// No photoionization at frequencies lower than the threshold freq occur

  x = sigma_phot_topbase (cont_ext_ptr2, freq);	//this is the cross-section
  integrand =
    x * freq * freq * freq / (exp (H_OVER_K * freq / tt) - 1) / fthresh;

  return (integrand);
}

/*****************************************************************************/

/******************************************* 
get_alpha_st - to get the stimulated recombination estimator 

	06may	ksl	57+ -- Modified for new structure
			Changed call to eliminate WindPtr altogether
*********************************************/
#define ALPHA_SP_CONSTANT 5.79618e-36

double
get_alpha_st (cont_ptr, xplasma)
     struct topbase_phot *cont_ptr;
     PlasmaPtr xplasma;
{
  double alpha_st_value;
  double fthresh, flast;
  double qromb ();
  double alpha_st_integrand ();

  temp_ext2 = xplasma->t_e;	//external for use in integrand
  temp_ext_rad = xplasma->t_r;
  cont_ext_ptr2 = cont_ptr;	//"
  fthresh = cont_ptr->freq[0];	//first frequency in list
  flast = cont_ptr->freq[cont_ptr->np - 1];	//last frequency in list
  alpha_st_value = qromb (alpha_st_integrand, fthresh, flast, 1e-4);

  /* The lines above evaluate the integral in alpha_sp. Now we just want to multiply 
     through by the appropriate constant. */
  if (cont_ptr->macro_info == 1 && geo.macro_simple == 0)
    {
      alpha_st_value = alpha_st_value * config[cont_ptr->nlev].g
	/ config[cont_ptr->uplev].g * pow (xplasma->t_e, -1.5);
    }
  else				//case for simple element
    {
      alpha_st_value = alpha_st_value * config[cont_ptr->nlev].g / ion[cont_ptr->nion + 1].g * pow (xplasma->t_e, -1.5);	//g for next ion up used
    }

  alpha_st_value = alpha_st_value * ALPHA_SP_CONSTANT * xplasma->w;

  return (alpha_st_value);
}



/******************************************************************************/

/* alpha_st_integrand. This returns the integrand for alpha_st at a chosen
   frequency*/

double
alpha_st_integrand (freq)
     double freq;		//frequency 
{
  double fthresh;
  double x;
  double integrand;
  double tt;
  double ttrr;

  fthresh = cont_ext_ptr2->freq[0];
  tt = temp_ext2;		//this is the electron temperature
  /* Also need the radiation temperature here */
  ttrr = temp_ext_rad;		//will do for now

  if (freq < fthresh)
    return (0.0);		// No recombination at frequencies lower than the threshold freq occur

  x = sigma_phot_topbase (cont_ext_ptr2, freq);	//this is the cross-section
  integrand =
    x * freq * freq * exp (H_OVER_K * (fthresh - freq) / tt) /
    (exp (H_OVER_K * freq / ttrr) - 1);

  return (integrand);
}

/*****************************************************************************/
/******************************************* 
get_alpha_st_e - to get the stimulated recombination estimator 

	06may	ksl	57+ -- Modified for new plasma structure.  In this
			case I eleimanted Wind altogether
*********************************************/
#define ALPHA_SP_CONSTANT 5.79618e-36

double
get_alpha_st_e (cont_ptr, xplasma)
     struct topbase_phot *cont_ptr;
     PlasmaPtr xplasma;
{
  double alpha_st_e_value;
  double fthresh, flast;
  double qromb ();
  double alpha_st_e_integrand ();

  temp_ext2 = xplasma->t_e;	//external for use in integrand
  temp_ext_rad = xplasma->t_r;	//"
  cont_ext_ptr2 = cont_ptr;	//"
  fthresh = cont_ptr->freq[0];	//first frequency in list
  flast = cont_ptr->freq[cont_ptr->np - 1];	//last frequency in list
  alpha_st_e_value = qromb (alpha_st_e_integrand, fthresh, flast, 1e-4);

  /* The lines above evaluate the integral in alpha_sp. Now we just want to multiply 
     through by the appropriate constant. */
  if (cont_ptr->macro_info == 1 && geo.macro_simple == 0)
    {
      alpha_st_e_value = alpha_st_e_value * config[cont_ptr->nlev].g
	/ config[cont_ptr->uplev].g * pow (xplasma->t_e, -1.5);
    }
  else				//case for simple element
    {
      alpha_st_e_value = alpha_st_e_value * config[cont_ptr->nlev].g / ion[cont_ptr->nion + 1].g * pow (xplasma->t_e, -1.5);	//g for next ion up used
    }

  alpha_st_e_value = alpha_st_e_value * ALPHA_SP_CONSTANT * xplasma->w;

  return (alpha_st_e_value);
}



/******************************************************************************/

/* alpha_st_e_integrand. This returns the integrand for alpha_st at a chosen
   frequency*/

double
alpha_st_e_integrand (freq)
     double freq;		//frequency 
{
  double fthresh;
  double x;
  double integrand;
  double tt;
  double ttrr;

  fthresh = cont_ext_ptr2->freq[0];
  tt = temp_ext2;		//this is the electron temperature
  /* Also need the radiation temperature here */
  ttrr = temp_ext_rad;		//will do for now

  if (freq < fthresh)
    return (0.0);		// No recombination at frequencies lower than the threshold freq occur

  x = sigma_phot_topbase (cont_ext_ptr2, freq);	//this is the cross-section
  integrand =
    x * freq * freq * exp (H_OVER_K * (fthresh - freq) / tt) /
    (exp (H_OVER_K * freq / ttrr) - 1) * freq / fthresh;

  return (integrand);
}

/*****************************************************************************/
