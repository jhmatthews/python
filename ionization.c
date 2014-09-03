

/***********************************************************
                                       Space Telescope Science Institute

 Synopsis:
 	ionization routines for the wind one cell at a time
	
 Arguments:		
	WindPtr w;

Returns:
 
Description:	
	The intent is that the routine ion_abundances is the steering routine for
	all calculations of the abundances
	
Notes:

History:
	98jul	ksl	Incorporated Christians mods with my comments 
	98nov	ksl	Created the steering routine ion_abundances
	00jan	ksl	Homogenized calls so that ionization_mod_on_the_spot
			and ionization_mod_on_the_spot_exact do not call the
			entire wind structure.  NOTE THAT THESE LAST ROUTINES
			HAVE NOT BEEN SERIOUSLY RECHECKED, SINCE THEY ARE NOT
			MUCH USED.
	00jan	ksl	Moved ionization_mod_on_the_spot and ionization_mod_on_the_spot_exact
			to "Legacy code" Legacy1.c
	01oct	ksl	Add calls to levels for calculation of level populations.  The
			calls are currently hardwired.
	08aug	ksl	60b	Evidently ksl modified the calls to ion_abundances
				previously, but leaving w as the variable was
				confusing.  Fixed this
	080808	ksl	62  - Removed option 4 for a partial detailed balance
			which seems totally obsolete at this point as
			we have not followed this up.  This was part of
			the cleanup of the ionizaton balance routines
			Also removed option 5, as this was not supported
			later in the program
	11nov	ksl	71 - Moved the so-called Sim power law approximation
			to its own routine to restore this routine to its intent 
			to be a steering routine and not something with a lot
			of code particular to one calculation
        12feb   nsh	71c - added options for mode 6 and 7, the two new ioinzation
			codes that compute saha abundances for pairs of ions, based on a
			good temperature for that pair, then corrects. 
**************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "atomic.h"
#include "python.h"




int
ion_abundances (xplasma, mode)
     PlasmaPtr xplasma;
     int mode;
{
  int ireturn;


  if (mode == IONMODE_ML93_FIXTE)
    {
/* on-the-spot approximation using existing t_e.   This routine does not attempt 
to match heating and cooling in the wind element! */

      if ((ireturn = nebular_concentrations (xplasma, NEBULARMODE_ML93)))
	{
	  Error
	    ("ionization_abundances: nebular_concentrations failed to converge\n");
	  Error
	    ("ionization_abundances: j %8.2e t_e %8.2e t_r %8.2e w %8.2e\n",
	     xplasma->j, xplasma->t_e, xplasma->w);
	}
    }
  else if (mode == IONMODE_LTE)
    {
/* LTE using t_r  (ksl - checked - 080808 */

      ireturn = nebular_concentrations (xplasma, NEBULARMODE_TR);
    }
  else if (mode == IONMODE_FIXED)
    {				//  Hardwired concentrations

      ireturn = fix_concentrations (xplasma, 0);
    }
  else if (mode == IONMODE_ML93)
    {
/* On the spot, with one_shot at updating t_e before calculating densities */

/* Shift values to old */
      xplasma->dt_e_old = xplasma->dt_e;
      xplasma->dt_e = xplasma->t_e - xplasma->t_e_old;	//Must store this before others
      xplasma->t_e_old = xplasma->t_e;
      xplasma->t_r_old = xplasma->t_r;
      xplasma->lum_rad_old = xplasma->lum_rad;

      ireturn = one_shot (xplasma, mode);

/* Convergence check */
      convergence (xplasma);
    }
  else if (mode == 4)
    {				//LTE with SIM correction this is called from define_wind, sim_alpha and sim_w are set to geo values in define_wind. Not sure this is ever called now, we thought it best to set the values to LTE when in define_wind.
      ireturn = nebular_concentrations (xplasma, 5);
    }
 // else if (mode == 5)     NSH 130626 removed mode
 //   {				// One shot at updating t_e before calculating densities using Stuart's power law correction
 //     ireturn = spectral_estimators (xplasma);	//Slight recoding - now this just computes the estimators - this was the code clogging up this mode. Now it looks like the others.

 //     xplasma->dt_e_old = xplasma->dt_e;
 //     xplasma->dt_e = xplasma->t_e - xplasma->t_e_old;	//Must store this before others
 //     xplasma->t_e_old = xplasma->t_e;
 //     xplasma->t_r_old = xplasma->t_r;
 //     xplasma->lum_rad_old = xplasma->lum_rad;


      /* Log ("NSH in this cell, we have %e AGN photons and %e disk photons\n",
         xplasma->ntot_agn, xplasma->ntot_disk); Removed, no longer gives reasonable answers due to banding */

   //   ireturn = one_shot (xplasma, mode);


/* Convergence check */
    //  convergence (xplasma);



  //  }
  else if (mode == IONMODE_PAIRWISE_ML93)
    {
      /* Feb 2012 new for mode 6. New abundances have been computed using pairwise Saha equation
         approach. We can now attempt to balance heating and cooling with the new abundance in the
         same way as mode 3. */

/* Shift values to old */
      xplasma->dt_e_old = xplasma->dt_e;
      xplasma->dt_e = xplasma->t_e - xplasma->t_e_old;	//Must store this before others
      xplasma->t_e_old = xplasma->t_e;
      xplasma->t_r_old = xplasma->t_r;
      xplasma->lum_rad_old = xplasma->lum_rad;

      ireturn = one_shot (xplasma, mode);

/* Convergence check */
      convergence (xplasma);
    }
  else if (mode == IONMODE_PAIRWISE_SPECTRALMODEL)
    {
/* Feb 2012 NSH - new for mode 7. KSL has moved a lot of the mechanics that used to be here into
 power_abundances. This, once called, calculates the weight and alpha for each band in this cell. There is a lot of code that was clogging up this routine. Once this is done, one_shot gets called from within that routine. */
      ireturn = spectral_estimators (xplasma);	/*Aug 2012 NSH - slight change to help integrate this into balance, power_estimators does the work of getting banded W and alpha. Then oneshot gets called. */
      xplasma->dt_e_old = xplasma->dt_e;
      xplasma->dt_e = xplasma->t_e - xplasma->t_e_old;	//Must store this before others
      xplasma->t_e_old = xplasma->t_e;
      xplasma->t_r_old = xplasma->t_r;
      xplasma->lum_rad_old = xplasma->lum_rad;


      ireturn = one_shot (xplasma, mode);


/* Convergence check */
      convergence (xplasma);
    }


  else
    {
      Error
	("ion_abundances: Could not calculate abundances for mode %d\n",
	 mode);
      exit (0);
    }

  /* If we want the Auger effect deal with it now. Initially, this is
     put in here, right at the end of the ionization calculation -
     the assumption is that the Auger effect is only for making minor
     ions so that the ionization balance of the other ions is not
     affected in an important way. */

  if (geo.auger_ionization == 1)
    {
      auger_ionization (xplasma);
    }


  return (ireturn);

}


/***********************************************************
              Space Telescope Science Institute

 Synopsis: convergence checks to see whether a single cell
	is or is not converging
	
 Arguments:		
	WindPtr w;

Returns:
 
Description:	

	The routine attempts to determine whether a cell is
	on the track to a final solution by checking whether 
	the electron and radiation temperatures are getting 
	smaller with cycle and whetehr the difference between
	heating and cooling is drroping.

	The routine returns a number between 0 and 3, 
	depending on the number of convergence checks that
	are passed.  If all convergence tests are pssed
	then the number returned will be 0

	The routine also adjust the gain which controls how
	far the electron temperture can change in a cycle.

	
Notes:

History:
	06may	ksl	57+: Began modifications to reflect new
			structure definition.  Note that it
			is likelely that this entire routine
			will ultimatetely be replaced because
			everything here should only be in the wind
        11sep  nsh     70f: Added lines to track which of the convergence criteria
			in each cell was being met
**************************************************************/
int
convergence (xplasma)
     PlasmaPtr xplasma;
{
  int trcheck, techeck, hccheck, whole_check, converging;
  double epsilon;

  trcheck = techeck = hccheck = converging = 0;
  xplasma->trcheck = xplasma->techeck = xplasma->hccheck = 0;	//NSH 70g - zero the global variables
  epsilon = 0.05;

  /* Check the fractional change in tempperatature and if is less than 
   * epsiolong increment trcheck and techeck 
   */

  if ((xplasma->converge_t_r =
       fabs (xplasma->t_r_old - xplasma->t_r) / (xplasma->t_r_old +
						 xplasma->t_r)) > epsilon)
    xplasma->trcheck = trcheck = 1;
  if (xplasma->t_e < TMAX)
	{
  	if ((xplasma->converge_t_e =
       		fabs (xplasma->t_e_old - xplasma->t_e) / (xplasma->t_e_old +
						 xplasma->t_e)) > epsilon)
    		xplasma->techeck = techeck = 1;
	if ((xplasma->converge_hc =
       		fabs (xplasma->heat_tot -
	     	(xplasma->lum_adiabatic + xplasma->lum_rad + xplasma->lum_dr +
	      	xplasma->lum_comp)) / fabs(xplasma->heat_tot + xplasma->lum_comp +
				     xplasma->lum_adiabatic +
				     xplasma->lum_dr + xplasma->lum_rad )) > epsilon)
    		xplasma->hccheck = hccheck = 1;
	}
  else //If the cell has reached the maximum temperature
	{
	xplasma->techeck = techeck = xplasma->hccheck = hccheck = 2; //we mark it as overlimit
	}

//110919 nsh modified line below to include the adiabatic cooling in the check that heating equals cooling
//111004 nsh further modification to include DR and compton cooling, now moved out of lum_rad

  /* Check whether the heating and colling balance to within epsilon and if so set hccheck to 1 */
 /* 130722 added a fabs to the bottom, since it is now conceivable that this could be negative if 
   lum_adiabatic is large and negative - and hence heating */

/* NSH 130711 - also changed to have fabs on top and bottom, since heating can now be negative!) */

/* NSH 130725 - moved the hc check to be within the if statement about overtemp - we cannot expect hc to converge if we are hitting the maximum temperature */
  /* whole_check is the sum of the temperature checks and the heating check */

  xplasma->converge_whole = whole_check = trcheck + techeck + hccheck;

  /* Converging is a situation where the change in electron
   * temperature is dropping with time and the cell is oscillating
   * around a temperature.  If that is the case, we drop the 
   * amount by which the temperature can change in this cycle
   */

  if (xplasma->dt_e_old * xplasma->dt_e < 0
      && fabs (xplasma->dt_e) > fabs (xplasma->dt_e_old))
    converging = 1;
  xplasma->converging = converging;

  if (converging == 1)
    {				// Not converging
      xplasma->gain *= 0.7;
      if (xplasma->gain < 0.1)
	xplasma->gain = 0.1;
    }
  else
    {
      xplasma->gain *= 1.1;
      if (xplasma->gain > 0.8)
	xplasma->gain = 0.8;
    }

  return (whole_check);
}





/***********************************************************
              Space Telescope Science Institute

 Synopsis:
   check_convergence -- Do a global check on how well the wind is converging
	
 Arguments:		
	WindPtr w;

Returns:
 
Description:	
	
Notes:
Notes--Eventually should absorb some of the calculations in wind_updates

History:
	05apr	ksl	55d -- Trivially eliminated MDIM from this
			in effort to make coordinate system independent.
			Whether something is in the wind or not is
			now totally determeined by the volume.
        06may	ksl	57+ -- Now using plasma structure
        11sep   nsh     70e added in lines to write out the split of temperature and luminosity convergence
**************************************************************/
int
check_convergence ()
{
  int n;
  int nconverge, nconverging, ntot;
  int nte, ntr, nhc;		//NSH 70g - three new counters for the different convergence criteria
  int nmax;			//NSH 130725 - counter for cells which are marked as converged, but over temp
  double xconverge, xconverging;

  nconverge = nconverging = ntot = 0;
  ntr = nte = nhc = nmax = 0;		//NSH 70i zero the counters

  for (n = 0; n < NPLASMA; n++)
    {
      ntot++;
      if (plasmamain[n].converge_whole == 0)
	nconverge++;
      if (plasmamain[n].trcheck == 0)	//NSH 70g - count up the three individual convergence criteria
	ntr++;
      if (plasmamain[n].techeck == 0)
	nte++;
      if (plasmamain[n].hccheck == 0)
	nhc++;
      if (plasmamain[n].techeck == 2)
	nmax++;
      if (plasmamain[n].converging == 0)
	nconverging++;

    }

  xconverge = ((double) nconverge) / ntot;
  xconverging = ((double) nconverging) / ntot;
  Log
    ("!!Check_converging: %4d (%.3f) converged and %4d (%.3f) converging of %d cells\n",
     nconverge, xconverge, nconverging, xconverging, ntot);
  Log ("!!Check_convergence_breakdown: t_r %4d t_e(real) %4d t_e(maxed) %4d hc(real) %4d\n", ntr, nte, nmax, nhc);	//NSH 70g split of what is converging
  Log
    ("Summary  convergence %4d %.3f  %4d  %.3f  %d  #  n_converged fraction_converged  converging fraction_converging total cells\n",
     nconverge, xconverge, nconverging, xconverging, ntot);
  Log_flush ();			/*NSH June 13 Added call to flush logfile */
  return (0);
}

/***********************************************************
                                       Space Telescope Science Institute

 Synopsis:
	one_shot calculates new densities of ions in a single element of the wind
	according to equation 11 of Lucy & Mazzali, after having found the
	temperature which matches heating and cooling for the previous
	densities
	
 Arguments:		
	PlasmaPtr xplasma;

Returns:
 
Description:
	
Notes:
	This routine attempts to match heating and cooling in the wind element!
	To do this it calls calc_te.  Based on the returned value of te, the
	routine then calculates the concentrations in the on-the-spot approximation.

	IT SEEMS LIKELY some code could be eliminated by simply having this routine
	call the on-the-spot routine directly.

History:
	98	ksl	Coded as part of python effort
	02jul	ksl	Added mode variable so could try detailed balance
	06may	ksl	57+ -- Switched to use plasma structue

**************************************************************/



int
one_shot (xplasma, mode)
     PlasmaPtr xplasma;
     int mode;

{
  double te_old, te_new;
  double gain;



  gain = xplasma->gain;


  te_old = xplasma->t_e;
  te_new = calc_te (xplasma, 0.7 * te_old, 1.3 * te_old);

  xplasma->t_e = (1 - gain) * te_old + gain * te_new;

/* NSH 130722 - NOTE - at this stage, the cooling terms are still those computed from the 'ideal' t_e, not the new t_e - this may be worth investigatiing. */
  if (xplasma->t_e > TMAX)
	{	
	xplasma->t_e = TMAX;
	}

//  Log ("One_shot: %10.2f %10.2f %10.2f\n", te_old, te_new, w->t_e);

/* Modes in the driving routines are not identical to those in nebular concentrations.
The next lines are an attempt to mediate this problem.  It might be better internally
at least to define a flag for using one shot, and have the modes take on the
meaning in nebular concentrations.
*/

  if (mode == IONMODE_ML93)
    mode = NEBULARMODE_ML93;
  else if (mode <= 1 || mode ==5 || mode >= 8)	/* modification to cope with mode 5 - SIM + two new modes in Feb 2012  + mode 5 now removed*/
    {

      Error ("one_shot: Sorry, Charlie, don't know how to process mode %d\n",
	     mode);
      exit (0);
    }

  if (xplasma->t_r > 10.)
    {				/* Then modify to an on the spot approx */
      if (nebular_concentrations (xplasma, mode))
	{
	  Error
	    ("ionization_on_the_spot: nebular_concentrations failed to converge\n");
	  Error
	    ("ionization_on_the_spot: j %8.2e t_e %8.2e t_r %8.2e w %8.2e\n",
	     xplasma->j, xplasma->t_e, xplasma->w);
	}
      if (xplasma->ne < 0 || VERY_BIG < xplasma->ne)
	{
	  Error ("ionization_on_the_spot: ne = %8.2e out of range\n",
		 xplasma->ne);
	}
    }
  else
    {
      Error ("ionization_on_the_spot: t_r exceptionally small %g\n",
	     xplasma->t_r);
      mytrap ();
      exit (0);
    }


  return (0);
}


/* calc_te determines and returns the electron temperature in the wind such that the energy emitted
   by the wind is equal to energy emitted.

   Description:
   calc_te does not modify any abundances.  It simply takes the current value of the heating in the
   cell and attempts to find the value of the electron temperature which will result in cooling which
   matches the heating.

   This approach is not entirely self consistent because if te is actually different then the
   abundances will be different and the heating will change as well.

   This routine is a kluge because it does not really deal with what happens if the cooling curve 
   has maxima and minima.

   xxxplasma is just a way to tranmit information to zero_emit

   History:

   98dec        ksl     Updated calls so that tmin and tmax were communicated externally,
			rather than hardwired
   01dec	ksl	Reversed 98dec decision as a result of massive changes for python38
   01dec	ksl	Assured that ww->t_e is updated here
   01dec	ksl	Added capability to modify the desired goal of calc_te from the full
			heating to something intermediate between the current value and the
			ultimate goal
   01dec	ksl	Rewrote to assure that boundaries will be bracketed properly and if
			not calc_te will handle
   04June       SS      Modified so that changes in the heating rate due to changes in the
                        temperature are included for macro atoms.
	06may	ksl	Modified for plasma structue
 */

PlasmaPtr xxxplasma;

double
calc_te (xplasma, tmin, tmax)
     PlasmaPtr xplasma;
     double tmin, tmax;
{
  double z1, z2;
  int macro_pops ();
		

  /* 110916 - ksl - Note that we assign a plasma pointer here to a fixed structure because
   * we need to call zbrent and we cannot pass the xplasma ptr directly
   */

  xxxplasma = xplasma;

  xplasma->t_e = tmin;
  z1 = zero_emit (tmin);
  xplasma->t_e = tmax;
  z2 = zero_emit (tmax);

  /* The way this works is that if we have a situation where the cooling
   * at tmax and tmin brackets the heating, then we use zbrent to improve
   * the estimated temperature, but if not we chose the best direction
   */

  if ((z1 * z2 < 0.0))
    {				// Then the interval is bracketed 
      xplasma->t_e = zbrent (zero_emit, tmin, tmax, 50.);
    }
  else if (fabs (z1) < fabs (z2))
    {
      xplasma->t_e = tmin;
    }
  else
    {
      xplasma->t_e = tmax;
    }
  /* With the new temperature in place for the cell, get the correct value of heat_tot.
     SS June  04 */

  /* ksl - I basically don't undestand what is going on here.  If we start using
   * macro atoms a lot we need to understand them better ??? - 
   * Look at zero emit as well 091611 */

  xplasma->heat_tot -= xplasma->heat_lines_macro;
  xplasma->heat_lines -= xplasma->heat_lines_macro;
  xplasma->heat_lines_macro = macro_bb_heating (xplasma, xplasma->t_e);
  xplasma->heat_tot += xplasma->heat_lines_macro;
  xplasma->heat_lines += xplasma->heat_lines_macro;

  xplasma->heat_tot -= xplasma->heat_photo_macro;
  xplasma->heat_photo -= xplasma->heat_photo_macro;
  xplasma->heat_photo_macro = macro_bf_heating (xplasma, xplasma->t_e);
  xplasma->heat_tot += xplasma->heat_photo_macro;
  xplasma->heat_photo += xplasma->heat_photo_macro;


  return (xplasma->t_e);

}



/* This is just a function which has a zero when total energy loss is equal to total energy gain */

double
zero_emit (t)
     double t;
{
  double difference;
  double total_emission ();
  int macro_pops ();
  double macro_bb_heating (), macro_bf_heating ();


  /* NOTE - IMPORTANT

     SS May 04: I'm removing adiabatic cooling for now. I'll need to be put
     back in once the heating/cooling tests are all sorted out. 

     SS June 04: Adiabatic cooling is still not switched on. But it is now 
     switched off in emission.c rather than here.

     SS July 04: Adiabatic cooling is now switched back on. These comments can
     all be deleted once this is tested (soon).
   */

  /* This block is not needed now - SS July 04
     if (xxxplasma->lum_adiabatic > 0)
     {
     Error("zero_emit: adiabatic cooling is switched on.\n");
     }
   */

  /*Original method */
  xxxplasma->t_e = t;


  /* Correct heat_tot for the change in temperature. SS June 04. */
  //macro_pops (xxxplasma, xxxplasma->ne);
  xxxplasma->heat_tot -= xxxplasma->heat_lines_macro;
  xxxplasma->heat_lines -= xxxplasma->heat_lines_macro;
  xxxplasma->heat_lines_macro = macro_bb_heating (xxxplasma, t);
  xxxplasma->heat_tot += xxxplasma->heat_lines_macro;
  xxxplasma->heat_lines += xxxplasma->heat_lines_macro;

  xxxplasma->heat_tot -= xxxplasma->heat_photo_macro;
  xxxplasma->heat_photo -= xxxplasma->heat_photo_macro;
  xxxplasma->heat_photo_macro = macro_bf_heating (xxxplasma, t);
  xxxplasma->heat_tot += xxxplasma->heat_photo_macro;
  xxxplasma->heat_photo += xxxplasma->heat_photo_macro;

  //  difference = (xxxplasma->heat_tot - total_emission (xxxplasma, 0., VERY_BIG));


   /* 70d - ksl - Added next line so that adiabatic cooling reflects the temperature we
   * are testing.  Adiabatic cooling is proportional to temperature
   */

  if (geo.adiabatic)
     {
     if (wmain[xxxplasma->nwind].div_v >= 0.0) //This is the case where we have adiabatic cooling - we want to retain the old behaviour, so we use the 'test' temperature to compute it. If div_v is less than zero, we don't do anything here, and so the existing value of adiabatic cooling is used - this was computed in wind_updates2d before the call to ion_abundances.
	{	     
	xxxplasma->lum_adiabatic = adiabatic_cooling (&wmain[xxxplasma->nwind], t);
	}     
     }
   else 
     {
     xxxplasma->lum_adiabatic=0.0; 
     }


  /* difference =
     xxxplasma->heat_tot - xxxplasma->lum_adiabatic -
     total_emission (&wmain[xxxplasma->nwind], 0., VERY_BIG); */


  /* 70g - nsh adding this line in next to calculate dielectronic recombination cooling without generating photons */
  compute_dr_coeffs (t);
  xxxplasma->lum_dr = total_dr (&wmain[xxxplasma->nwind], t);

/* 70g compton cooling calculated here to avoid generating photons */
  xxxplasma->lum_comp = total_comp (&wmain[xxxplasma->nwind], t);

  difference = xxxplasma->heat_tot - xxxplasma->lum_adiabatic - xxxplasma->lum_dr - xxxplasma->lum_comp - total_emission (&wmain[xxxplasma->nwind], 0., VERY_BIG);	//NSH 1110 - total emission no longer computes compton.*/




  return (difference);
}
