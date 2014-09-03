

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "atomic.h"
#include "python.h"
#include <gsl/gsl_block.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
//#include <gsl/gsl_blas.h>
#include "my_linalg.h"

/*****************************************************************************

                                    Imperial College London
Synopsis:
	Matom is the core of the implementation of Leon Lucy's Macro Atom
	approach to dealing with radiation-matter interactions. It is called whenever
	a photon packet activates a macro atom. As input it takes the process by
	which activation occurred (identified by the label "nres") which allow it 
	to deduce the level that has been excited. It then calculates all the 
	Macro Atom jumping/deactivaion probabilities following Lucy and so deduces
	the process by which deactivation occurs. At output "nres" identifies this 
	process and the packet information has been updates. 

Arguments:

       WindPtr w                   the ptr to the structure defining the wind
       PhotPtr p                   the packet at the point of activation
       int nres                    the process which activates the Macro Atom
  
Returns:
       int nres                    the process by which deactivation occurs
       PhotPtr p                   the packet following deactivation
       int escape                  to tell us whether the de-activation of the
                                   macro atom is via an r-packet (escape = 1) or
                                   a k-packet (escape = 0)

Description:


Notes: 

ksl-- There is a lot of arithmetic required to keep track of indices.  I would
be inclined to figure out a way to avoid this.  One possibility would be to record the
go_to level in an array when one is calculating the emission probabilities. This would
make it easier to add new processes, I suspect. 

It would be possible to convert the for loop to a while statement.  This would avoid the
break statement in the middle.  

History:
          Feb 04  SS   Coding began.
          Mar 04  SS   Minor changes made (based on comments from ksl)
	04apr	ksl	A. Eliminated some unnecessary if statements.  These are indicated
			by the words extra.  Stuart should remove the lines assuming he
			agrees.
			B. Changed the logic of the major do loop so that the break is moved
			higher in the do loop.  Stuart, the intent is to make the
			code more readable, and to reduce the level of nesting. 
			C. Corrected an errror that on upward bb transitions seemed to leave
			the macro-atom in the same state as previously.
			D. Shamelessly modified some of the comments to make it more
			straightforward for me to understand.
			E. m = m + 1 --> m++ in several places
			F. Modified selections to allow simple lines to be read in but
			not to affect matom.  nlines_macro is the number on macro-lines
			read in and these are the first elements in the line structure
			G. Modified the way matom determines that a transition is a photo
			ionization transition bo be consitent with change made in
			resonate.
          Apr 04   SS   Modifications to include collisional jumps/deactivation included.
          May 04   SS   Bug corrected for identification of line (nres is the place of the
                        line in the ORDERED list, not the input list).
          Jun 04   SS   Modified to carry the "escape" variable to identify de-activation
                        via r-packets and k-packets - avoids the need to call kpkt from
                        within this routine
          Jun 04   SS   Adding collisional ionization as an upwards jump possibility. No collisional
                        recombination for the moment.
          Jun 04   SS   Now putting in collisional recombination too.
          July04   SS   Modifying so that this routine does not call alpha_sp but rather uses 
                        pre-computed (stored values) for the spontaneous recombination coefficient.
                        This is an effort to speed up the code.
	06may	ksl	57+ -- Adapted for use with plsama structure.  Changed call to
			eliminate passing entire w array
	06jun	ksl	57g -- Split macro variables into a separate structure. The structue
			which is createdin gridwind, is only created if there are macroatoms.
	06jul	ksl	57h -- In the process of speeding up the program in the simple 
			non-macro atom case, I tried to make sure that matom and the 
			derivative routiens were not called at all.  Note that at present
			the macroatom case is quite slow, due to what is happening in
			matom.  I am suspicious that it could be speeded up a lot.
        07jul     SS    Experimenting with retaining jumping/emission probabilities to save time.

************************************************************/

int
matom (p, nres, escape)
     PhotPtr p;
     int *nres;
     int *escape;
{
  struct lines *line_ptr;
  struct topbase_phot *cont_ptr;
  int uplvl, uplvl_old;
  double jprbs[2 * (NBBJUMPS + NBFJUMPS)];
  double eprbs[NBBJUMPS + NBFJUMPS];
  double pjnorm, penorm;
  double threshold, run_tot;
  double sp_rec_rate;
  int n, m;
  int njumps;
  int nbbd, nbbu, nbfd, nbfu;
  double t_e, ne;
  double bb_cont, choice, bf_cont;
  WindPtr one;
  double rad_rate, coll_rate;
  PlasmaPtr xplasma;
  MacroPtr mplasma;
  double jprbs_known[NLEVELS_MACRO][2 * (NBBJUMPS + NBFJUMPS)],
    eprbs_known[NLEVELS_MACRO][2 * (NBBJUMPS + NBFJUMPS)];
  double pjnorm_known[NLEVELS_MACRO], penorm_known[NLEVELS_MACRO];
  int prbs_known[NLEVELS_MACRO];


  for (n = 0; n < NLEVELS_MACRO; n++)
    {
      prbs_known[n] = -1;	//flag all as unknown
    }


  one = &wmain[p->grid];	//This is to identify the grid cell in which we are
  xplasma = &plasmamain[one->nplasma];
  check_plasma (xplasma, "matom");

  mplasma = &macromain[xplasma->nplasma];


  t_e = xplasma->t_e;		//electron temperature 
  ne = xplasma->ne;		//electron number density


  /* The first step is to identify the configuration that has been excited. */


  if (*nres < NLINES)		//this means that it was a line excitation CHECK
    {
      uplvl = lin_ptr[*nres]->nconfigu;

    }
  else if (*nres > NLINES)	//this means it was photoionisation
    {
      uplvl = phot_top[*nres - NLINES - 1].uplev;
    }
  else
    {
      Error ("matom: upper level not identified. nres = %d\n", *nres);
      exit (0);
    }

  /* Now follows the main loop to govern the macro atom jumps. Keeps jumping until
     an emission occurs (at which point it breaks). */


  for (njumps = 0; njumps < MAXJUMPS; njumps++)
    {
      /*  The excited configuration is now known. Now compute all the probabilities of deactivation
         /jumping from this configuration. Then choose one. */

      nbbd = config[uplvl].n_bbd_jump;	//store these for easy access -- number of bb downward jumps
      nbbu = config[uplvl].n_bbu_jump;	// number of bb upward jump from this configuration
      nbfd = config[uplvl].n_bfd_jump;	// number of bf downward jumps from this transition
      nbfu = config[uplvl].n_bfu_jump;	// number of bf upward jumps from this transiion

      if (prbs_known[uplvl] != 1)
	{
	  /* Start by setting everything to 0  */

	  m = 0;		//m counts the total number of possible ways to leave the level
	  pjnorm = 0.0;		//stores the total jump probability
	  penorm = 0.0;		//stores the total emission probability

	  for (n = 0; n < nbbd + nbfd; n++)
	    {
	      eprbs_known[uplvl][n] = eprbs[n] = 0;	//stores the individual emission probabilities SS
	      jprbs_known[uplvl][n] = jprbs[n] = 0;	//stores the individual jump probabilities SS
	    }
	  for (n = nbbd + nbfd; n < nbbd + nbfd + nbbu + nbfu; n++)
	    {
	      jprbs_known[uplvl][n] = jprbs[n] = 0;	/*slots for upward jumps */
	    }
	  /* Finished zeroing. */


	  /* bb */

	  /* First downward jumps. (I.e. those that have emission probabilities. */

	  /* For bound-bound decays the jump probability is A-coeff * escape-probability * energy */
	  /* At present the escape probability is only approximated (p_escape). This should be improved. */
	  /* The collisional contribution to both the jumping and deactivation probabilities are now added (SS, Apr04) */

	  for (n = 0; n < nbbd; n++)
	    {
	      line_ptr = &line[config[uplvl].bbd_jump[n]];

	      rad_rate = (a21 (line_ptr) * p_escape (line_ptr, xplasma));
	      coll_rate = q21 (line_ptr, t_e);	// this is multiplied by ne below

	      if (coll_rate < 0)
		{
		  coll_rate = 0;
		}

	      bb_cont = rad_rate + (coll_rate * ne);
	      jprbs_known[uplvl][m] = jprbs[m] = bb_cont * config[line_ptr->nconfigl].ex;	//energy of lower state

	      eprbs_known[uplvl][m] = eprbs[m] = bb_cont * (config[uplvl].ex - config[line[config[uplvl].bbd_jump[n]].nconfigl].ex);	//energy difference

	      if (jprbs[m] < 0.)	//test (can be deleted eventually SS)
		{
		  Error ("Negative probability (matom, 1). Abort.");
		  exit (0);
		}
	      if (eprbs[m] < 0.)	//test (can be deleted eventually SS)
		{
		  Error ("Negative probability (matom, 2). Abort.");
		  exit (0);
		}

	      pjnorm += jprbs[m];
	      penorm += eprbs[m];
	      m++;
	    }

	  /* bf */
	  for (n = 0; n < nbfd; n++)
	    {

	      cont_ptr = &phot_top[config[uplvl].bfd_jump[n]];	//pointer to continuum
	      if (n < 25)
		{
		  sp_rec_rate = mplasma->recomb_sp[config[uplvl].bfd_indx_first + n];	//need this twice so store it
		  bf_cont =
		    (sp_rec_rate + q_recomb (cont_ptr, t_e) * ne) * ne;
		}
	      else
		{
		  bf_cont = 0.0;
		}

	      jprbs_known[uplvl][m] = jprbs[m] = bf_cont * config[phot_top[config[uplvl].bfd_jump[n]].nlev].ex;	//energy of lower state
	      eprbs_known[uplvl][m] = eprbs[m] = bf_cont * (config[uplvl].ex - config[phot_top[config[uplvl].bfd_jump[n]].nlev].ex);	//energy difference
	      if (jprbs[m] < 0.)	//test (can be deleted eventually SS)
		{
		  Error ("Negative probability (matom, 3). Abort.");
		  exit (0);
		}
	      if (eprbs[m] < 0.)	//test (can be deleted eventually SS)
		{
		  Error ("Negative probability (matom, 4). Abort.");
		  exit (0);
		}
	      pjnorm += jprbs[m];
	      penorm += eprbs[m];
	      m++;
	    }

	  /* Now upwards jumps. */

	  /* bb */
	  /* For bound-bound excitation the jump probability is B-coeff times Jbar with a correction 
	     for stimulated emission. To avoid the need for recalculation all the time, the code will
	     be designed to include the stimulated correction in Jbar - i.e. the stimulated correction
	     factor will NOT be included here. (SS) */
	  /* There is no emission probability for upwards transitions. */
	  /* Collisional contribution to jumping probability added. (SS,Apr04) */

	  for (n = 0; n < nbbu; n++)
	    {
	      line_ptr = &line[config[uplvl].bbu_jump[n]];
	      rad_rate =
		(b12 (line_ptr) *
		 mplasma->jbar_old[config[uplvl].bbu_indx_first + n]);

	      coll_rate = q12 (line_ptr, t_e);	// this is multiplied by ne below

	      if (coll_rate < 0)
		{
		  coll_rate = 0;
		}

	      jprbs_known[uplvl][m] = jprbs[m] = ((rad_rate) + (coll_rate * ne)) * config[uplvl].ex;	//energy of lower state



	      if (jprbs[m] < 0.)	//test (can be deleted eventually SS)
		{
		  Error ("Negative probability (matom, 5). Abort.");
		  exit (0);
		}
	      pjnorm += jprbs[m];
	      m++;
	    }

	  /* bf */
	  for (n = 0; n < nbfu; n++)
	    {
	      /* For bf ionization the jump probability is just gamma * energy
	         gamma is the photoionisation rate. Stimulated recombination also included. */
	      cont_ptr = &phot_top[config[uplvl].bfu_jump[n]];	//pointer to continuum

	      jprbs_known[uplvl][m] = jprbs[m] = (mplasma->gamma_old[config[uplvl].bfu_indx_first + n] - (mplasma->alpha_st_old[config[uplvl].bfu_indx_first + n] * xplasma->ne * den_config (xplasma, cont_ptr->uplev) / den_config (xplasma, cont_ptr->nlev)) + (q_ioniz (cont_ptr, t_e) * ne)) * config[uplvl].ex;	//energy of lower state
	      if (jprbs[m] < 0.)	//test (can be deleted eventually SS)
		{
		  Error ("Negative probability (matom, 6). Abort?\n");
		  Log
		    ("mplasma->gamma_old[uplvl][n] %g, (mplasma->alpha_st_old[uplvl][n] * xplasma->ne * den_config (xplasma, cont_ptr->uplev) / den_config (xplasma, cont_ptr->nlev)) %g\n",
		     mplasma->gamma_old[config[uplvl].bfu_indx_first + n],
		     (mplasma->alpha_st_old[config[uplvl].bfu_indx_first +
					    n] * xplasma->ne *
		      den_config (xplasma,
				  cont_ptr->uplev) / den_config (xplasma,
								 cont_ptr->
								 nlev)));
		  jprbs_known[uplvl][m] = jprbs[m] = 0.0;

		}
	      pjnorm += jprbs[m];
	      m++;
	    }

	  pjnorm_known[uplvl] = pjnorm;
	  penorm_known[uplvl] = penorm;
	  prbs_known[uplvl] = 1;

	}

      /* Probabilities of jumping (j) and emission (e) are now known. The normalisation factor pnorm has
         also been recorded. the integer m now gives the total number of possibilities too. 
         now select what happens next. Start by choosing the random threshold value at which the
         event will occur. */

      threshold = ((rand () + 0.5) / MAXRAND);

      if ((pjnorm_known[uplvl] + penorm_known[uplvl]) <= 0.0)
	{
	  Error ("matom: macro atom level has no way out %d %g %g\n", uplvl,
		 pjnorm_known[uplvl], penorm_known[uplvl]);
	  exit (0);
	}

      if (((pjnorm_known[uplvl] /
	    (pjnorm_known[uplvl] + penorm_known[uplvl])) < threshold)
	  || (pjnorm_known[uplvl] == 0))
	break;			// An emission occurs and so we leave the for loop.

      uplvl_old = uplvl;

// Continue on if a jump has occured 

      /* Use a running total to decide where event occurs. */
      run_tot = 0;

      n = 0;
      threshold = ((rand () + 0.5) / MAXRAND);
      threshold = threshold * pjnorm_known[uplvl_old];
      while (run_tot < threshold)
	{
	  run_tot += jprbs_known[uplvl_old][n];
	  n++;
	}
      // This added to prevent case where theshold is essentially 0. It is already checked that the jumping probability is not zero
      if (n > 0)
	{
	  n = n - 1;
	}

      /* n now identifies the jump that occurs - now set the new level. */
      if (n < nbbd)
	{			/* bb downwards jump */
	  uplvl = line[config[uplvl].bbd_jump[n]].nconfigl;
	}
      else if (n < (nbbd + nbfd))
	{			/* bf downwards jump */
	  uplvl = phot_top[config[uplvl].bfd_jump[n - nbbd]].nlev;
	}
      else if (n < (nbbd + nbfd + nbbu))
	{			/* bb upwards jump */
	  uplvl = line[config[uplvl].bbu_jump[n - nbbd - nbfd]].nconfigu;
	}
      else if (n < (nbbd + nbfd + nbbu + nbfu))
	{			/* bf upwards jump */
	  uplvl =
	    phot_top[config[uplvl].bfu_jump[n - nbbd - nbfd - nbbu]].uplev;
	}
      else
	{
	  Error ("Trying to jump but nowhere to go! Matom. Abort");
	  exit (0);
	}

/* ksl: Check added to verify that the level actually changed */
      if (uplvl_old == uplvl)
	{
	  Error ("matom: uplvl did not change with jump: %d %d\n", uplvl, n);

	}

    }

  /* When is gets here either the sum has reached maxjumps: didn't find an emission: this is
     an error and stops the run OR an emission mechanism has been chosen in which case all is well. SS */

  if (njumps == MAXJUMPS)
    {
      Error ("Matom: jumped %d times with no emission. Abort.\n", MAXJUMPS);
      exit (0);
    }


  /* If it gets here then an emission has occurred. SS */
  /* Use a running total to decide where event occurs. SS */

  run_tot = 0;
  n = 0;
  threshold = ((rand () + 0.5) / MAXRAND);
  threshold = threshold * penorm_known[uplvl];	//normalise to total emission prob.
  while (run_tot < threshold)
    {
      run_tot += eprbs_known[uplvl][n];
      n++;
    }
  n = n - 1;
  /* n now identifies the jump that occurs - now set nres for the return value. */
  if (n < nbbd)
    {				/* bb downwards jump */
      /* With collisions included we need to decide whether the deactivation is radiative (r-packet)
         or collisional (k-packet). Get a random number and then use the ratio of the collisional
         emission probability to the (already known) collisional+radiative probability to decide whether
         collisional or radiative deactivation occurs. */
      choice = ((rand () + 0.5) / MAXRAND);	// the random number

      line_ptr = &line[config[uplvl].bbd_jump[n]];	//pointer for the bb transition

      rad_rate = a21 (line_ptr) * p_escape (line_ptr, xplasma);

      /* JM130716 in some old versions the coll_rate was set incorrectly here- 
         it needs to be multiplied by electron density */
      coll_rate = q21 (line_ptr, t_e) * ne;

      if (coll_rate < 0)
	{
	  coll_rate = 0.;
	}

      if (choice > (coll_rate / (rad_rate + coll_rate)))
	{
	  *escape = 1;		//don't need to follow this routine with a k-packet
	  /* This is radiative deactivation. */
	  *nres = line[config[uplvl].bbd_jump[n]].where_in_list;
	  p->freq = line[config[uplvl].bbd_jump[n]].freq;
	  /* This is the comoving frequency - changed to rest frequency by doppler */
	}
      else
	{			/* This is collisional deactivation. In this case a k-packet is created and so it must be processed by
				   the appropriate routine. (SS, Apr04) */
	  /* This is done by returning escape = 0 */
	  *escape = 0;

	}
    }
  else if (n < (nbbd + nbfd))
    {				/* bf downwards jump */
      /* With collisional recombination included we need to decide whether to deactivate
         radiatively or make a k-packet. */

      cont_ptr = &phot_top[config[uplvl].bfd_jump[n - nbbd]];

      rad_rate = mplasma->recomb_sp[config[uplvl].bfd_indx_first + n - nbbd];	//again using recomb_sp rather than alpha_sp (SS July 04)
      coll_rate = q_recomb (cont_ptr, t_e) * ne;

      choice = ((rand () + 0.5) / MAXRAND);	// the random number

      if (choice > (coll_rate / (rad_rate + coll_rate)))
	{			//radiative deactivation
	  *escape = 1;
	  *nres = config[uplvl].bfd_jump[n - nbbd] + NLINES + 1;
	  /* continuua are indicated by nres > NLINES */
	  p->freq = phot_top[config[uplvl].bfd_jump[n - nbbd]].freq[0]
	    -
	    (log (1. - (rand () + 0.5) / MAXRAND) * xplasma->t_e / H_OVER_K);
	  /* Co-moving frequency - changed to rest frequency by doppler */
	  /*Currently this assumed hydrogenic shape cross-section - Improve */
	}
      else
	{			//collisional deactivation
	  *escape = 0;
	}

    }
  else
    {
      Error
	("Trying to emitt from Macro Atom but no available route (matom). Abort.");
      exit (0);
    }

  return (0);
}

/********************************************************************************/

/*
  b12  very similar to a21 - returns the b12 Einstein coefficient.
  History:
  2004feb       coded by S Sim
*/

// !! ksl OK, Probably should be moved to atomic.c for consistency, but that can be done
// later.  
#define B12_CONSTANT 5.01983e25

struct lines *b12_line_ptr;
double b12_a;

double
b12 (line_ptr)
     struct lines *line_ptr;
{
  double freq;

  if (b12_line_ptr != line_ptr)
    {
      freq = line_ptr->freq;
      b12_a = B12_CONSTANT * line_ptr->f / freq;
      b12_line_ptr = line_ptr;
    }

  return (b12_a);
}

/************************************************************/


/* As for similar routines in recomb.c, in order to use the integrator the 
   following external structures are used (SS)*/

/* This relates to the alpha_sp routines at the end of this file */

struct topbase_phot *cont_ext_ptr;	//continuum pointer passed externally
double temp_ext;		//temperature passed externally
int temp_choice;		//choice of type of calcualation for alpha_sp

/*****************************************************************************/

/*
alpha_sp - to govern the calculation of the spontaneous recombination rate.

The rate is given by (4 pi /c2) (gu/gl) (h2/2 pi m k T)^(3/2) 
times the integral of   a(nu) nu2 exp [(chi- h nu)/kT].

04jul30	ksl	Modified so that one does not need to have multiple versions
		of the code depending on how whether the integrand is 
		multiplied by 1, f/fthresh, or f/fthresh-1.  This was
		done eliminate alpha_sp_e as a separate set of routines
		and to assure that bf rates are positive 
			ichoice = 0   --> spontanous recombination
			ichoice = 1   --> energy weighted recombination 
			ichoice = 2   --> the difference between energy_weighted
					and spontaneous

	06may	ksl	57+ -- Modified to use plasma structure
*/
#define ALPHA_SP_CONSTANT 5.79618e-36

double
alpha_sp (cont_ptr, xplasma, ichoice)
     struct topbase_phot *cont_ptr;
     PlasmaPtr xplasma;
     int ichoice;
{
  double alpha_sp_value;
  double fthresh, flast;
  double qromb ();
  double alpha_sp_integrand ();

  temp_choice = ichoice;
  temp_ext = xplasma->t_e;	//external for use in alph_sp_integrand
  cont_ext_ptr = cont_ptr;	//"
  fthresh = cont_ptr->freq[0];	//first frequency in list
  flast = cont_ptr->freq[cont_ptr->np - 1];	//last frequency in list
  alpha_sp_value = qromb (alpha_sp_integrand, fthresh, flast, 1e-4);

  /* The lines above evaluate the integral in alpha_sp. Now we just want to multiply 
     through by the appropriate constant. */
  if (cont_ptr->macro_info == 1 && geo.macro_simple == 0)
    {
      alpha_sp_value = alpha_sp_value * config[cont_ptr->nlev].g
	/ config[cont_ptr->uplev].g * pow (xplasma->t_e, -1.5);
    }
  else				//case for simple element
    {
      alpha_sp_value = alpha_sp_value * config[cont_ptr->nlev].g / ion[cont_ptr->nion + 1].g * pow (xplasma->t_e, -1.5);	//g for next ion up used
    }

  alpha_sp_value = alpha_sp_value * ALPHA_SP_CONSTANT;

  return (alpha_sp_value);
}



/******************************************************************************/

/* alpha_sp_integrand. This returns the integrand for alpha_sp at a chosen
   frequency*/

double
alpha_sp_integrand (freq)
     double freq;		//frequency 
{
  double fthresh;
  double x;
  double integrand;
  double tt;

  fthresh = cont_ext_ptr->freq[0];
  tt = temp_ext;

  if (freq < fthresh)
    return (0.0);		// No recombination at frequencies lower than the threshold freq occur

  x = sigma_phot_topbase (cont_ext_ptr, freq);	//this is the cross-section
  integrand = x * freq * freq * exp (H_OVER_K * (fthresh - freq) / tt);

  if (temp_choice == 1)
    return (integrand * freq / fthresh);	//energy weighed case
  if (temp_choice == 2)
    return (integrand * (freq - fthresh) / fthresh);	// difference case
  return (integrand);		//spontanoues case
}

/*****************************************************************************/
/****************************************************************************/


/* kpkt
This deals with the elimination of k-packets. Whenever a k-packet is created it is
immediately destroyed insitu by this routine. At output "nres" identified the process
that destroys the k-packet and the packet information has been updated in the same
way as at the end of matom. */

/************************************************************
                                    Imperial College London
Synopsis:

Arguments:

       WindPtr w                   the ptr to the structure defining the wind
       PhotPtr p                   the packet at the point of activation
       int nres                    the process which activates the Macro Atom

Returns:
       int nres                    the process by which deactivation occurs
       PhotPtr p                   the packet following deactivation
       int escape                  indicated whether removal of k-packet was by creating
                                   an r-packet (escape = 1) or by exciting a macro atom
                                   (escape = 0)

Description:


Notes: 


History:
          Mar 04  SS   Coding began.
          Apr 04  SS   Various improvements including the addition of ff and collisions.
          May 04  SS   Minor changes made to collisional cooling rates (bug fixed)
                       and fb cooling for simple continua.
          May 04  SS   Modified to work for case with all "simple" ions.
          May 04  SS   Modified to use "scattering probability" formalism for 
                       simple ion cooling rates.
          Jun 04  SS   Modified to include the "escape" variable to identify excitation of
                       macro atoms. This removes the need to call matom from within this routine.
          Jun 04  SS   Modified to include collisional ionization as a cooling term.
          July04  SS   Modified to use recomb_sp(_e) rather than alpha_sp(_e) to speed up.
	06may	ksl	57+ -- Modified to use new plasma array.  Eliminated passing
			entire w array
	131030	JM 		-- Added adiabatic cooling as possible kpkt destruction choice
          
************************************************************/

int
kpkt (p, nres, escape)
     PhotPtr p;
     int *nres;
     int *escape;
{

  int i;
  int ulvl;
  double cooling_bf[NTOP_PHOT];
  double cooling_bf_col[NTOP_PHOT];	//collisional cooling in bf transitions
  double cooling_bb[NLINES];
  double cooling_adiabatic;
  struct topbase_phot *cont_ptr;
  struct lines *line_ptr;
  double cooling_normalisation;
  double destruction_choice;
  double electron_temperature;
  double cooling_bbtot, cooling_bftot, cooling_bf_coltot;
  double lower_density, upper_density;
  double cooling_ff;
  WindPtr one;
  PlasmaPtr xplasma;
  MacroPtr mplasma;

  double coll_rate, rad_rate;
  double q_ioniz ();


  /* Idea is to calculated the cooling
     terms for all the processes and then choose one at random to destroy the k-packet and
     turn it back into a photon bundle. 
     The routine considers bound-free, collision excitation and ff
     emission. */

  /* The loop here runs over all the bf processes, regardless of whether they are macro or simple
     ions. The reason that we can do this is that in the macro atom method stimulated recombination
     has already been considered before creating the k-packet (by the use of gamma-twiddle) in "scatter"
     and so we only have spontaneous recombination to worry about here for ALL cases. */


  one = &wmain[p->grid];
  xplasma = &plasmamain[one->nplasma];
  check_plasma (xplasma, "kpkt");
  mplasma = &macromain[xplasma->nplasma];

  electron_temperature = xplasma->t_e;


  /* ksl 091108 - If the kpkt destruction rates for this cell are not known they are calculated here.  This happens
   * every time the wind is updated */

  if (mplasma->kpkt_rates_known != 1)
    {
      cooling_normalisation = 0.0;
      cooling_bftot = 0.0;
      cooling_bbtot = 0.0;
      cooling_ff = 0.0;
      cooling_bf_coltot = 0.0;

      for (i = 0; i < ntop_phot; i++)
	{
	  cont_ptr = &phot_top[i];
	  ulvl = cont_ptr->uplev;

	  if (cont_ptr->macro_info == 1 && geo.macro_simple == 0)
	    {
	      upper_density = den_config (xplasma, ulvl);
	      /* SS July 04 - for macro atoms the recombination coefficients are stored so use the
	         stored values rather than recompue them. */
	      cooling_bf[i] = mplasma->cooling_bf[i] =
		upper_density * H * cont_ptr->freq[0] *
		(mplasma->recomb_sp_e[config[ulvl].bfd_indx_first +
				      cont_ptr->down_index]);
	      // _sp_e is defined as the difference 
	    }
	  else
	    {
	      upper_density = xplasma->density[cont_ptr->nion + 1];

	      cooling_bf[i] = mplasma->cooling_bf[i] =
		upper_density * H * cont_ptr->freq[0] *
		(xplasma->recomb_simple[i]);
	    }


	  /* Note that the electron density is not included here -- all cooling rates scale
	     with the electron density so I've factored it out. */
	  if (cooling_bf[i] < 0)
	    {
	      Error ("kpkt: bf cooling rate negative. Density was %g\n",
		     upper_density);
	      Error ("alpha_sp(cont_ptr, xplasma,2) %g \n",
		     alpha_sp (cont_ptr, xplasma, 2));
	      Error ("i, ulvl, ntop_phot, nion %d %d %d %d\n", i, ulvl,
		     ntop_phot, cont_ptr->nion);
	      Error ("nlev, z, istate %d %d %d \n", cont_ptr->nlev,
		     cont_ptr->z, cont_ptr->istate);
	      Error ("freq[0] %g\n", cont_ptr->freq[0]);
	      cooling_bf[i] = mplasma->cooling_bf[i] = 0.0;
	    }
	  else
	    {
	      cooling_bftot += cooling_bf[i];
	    }
	  cooling_normalisation += cooling_bf[i];
	  if (cont_ptr->macro_info == 1 && geo.macro_simple == 0)
	    {
	      /* Include collisional ionization as a cooling term in macro atoms. Don't include
	         for simple ions for now.  SS */

	      lower_density = den_config (xplasma, cont_ptr->nlev);
	      cooling_bf_col[i] = mplasma->cooling_bf_col[i] =
		lower_density * H * cont_ptr->freq[0] * q_ioniz (cont_ptr,
								 electron_temperature);

	      cooling_bf_coltot += cooling_bf_col[i];

	      cooling_normalisation += cooling_bf_col[i];

	    }



	}

      /* end of loop over ntop_phot */

      for (i = 0; i < nlines; i++)
	{
	  line_ptr = &line[i];
	  if (line_ptr->macro_info == 1 && geo.macro_simple == 0)
	    {			//It's a macro atom line and so the density of the upper level is stored
	      cooling_bb[i] = mplasma->cooling_bb[i] =
		den_config (xplasma, line_ptr->nconfigl) * q12 (line_ptr,
								electron_temperature)
		* line_ptr->freq * H;

	      /* Note that the electron density is not included here -- all cooling rates scale
	         with the electron density so I've factored it out. */
	    }
	  else
	    {			//It's a simple line. Get the upper level density using two_level_atom

	      two_level_atom (line_ptr, xplasma, &lower_density,
			      &upper_density);

	      /* the collisional rate is multiplied by ne later */
	      coll_rate =
		q21 (line_ptr,
		     electron_temperature) * (1. -
					      exp (-H_OVER_K *
						   line_ptr->freq /
						   electron_temperature));

	      cooling_bb[i] =
		(lower_density * line_ptr->gu / line_ptr->gl -
		 upper_density) * coll_rate / (exp (H_OVER_K *
						    line_ptr->freq /
						    electron_temperature) -
					       1.) * line_ptr->freq * H;

	      rad_rate = a21 (line_ptr) * p_escape (line_ptr, xplasma);

	      /* Now multiply by the scattering probability - i.e. we are only going to consider bb cooling when
	         the photon actually escapes - we don't to waste time by exciting a two-level macro atom only so that
	         it makes another k-packet for us! (SS May 04) */


	      cooling_bb[i] *=
		rad_rate / (rad_rate + (coll_rate * xplasma->ne));
	      mplasma->cooling_bb[i] = cooling_bb[i];
	    }

	  if (cooling_bb[i] < 0)
	    {
	      cooling_bb[i] = mplasma->cooling_bb[i] = 0.0;
	    }
	  else
	    {
	      cooling_bbtot += cooling_bb[i];
	    }
	  cooling_normalisation += cooling_bb[i];
	}

      /* end of loop over nlines  */


      /* 57+ -- This might be modified later since we "know" that xplasma cannot be for a grid with zero
         volume.  Recall however that vol is part of the windPtr */
      if (one->vol > 0)
	{
	  cooling_ff = mplasma->cooling_ff =
	    total_free (one, xplasma->t_e, 0.0,
			VERY_BIG) / one->vol / xplasma->ne;
	}
      else
	{
	  /* SS June 04 - This should never happen, but sometimes it does. I think it is because of 
	     photons leaking from one cell to another due to the push-through-distance. It is sufficiently
	     rare (~1 photon in a complete run of the code) that I'm not worrying about it for now but it does
	     indicate a real problem somewhere. */

	  /* SS Nov 09: actually I've not seen this problem for a long
	     time. Don't recall that we ever actually fixed it,
	     however. Perhaps the improved volume calculations
	     removed it? We delete this whole "else" if we're sure
	     volumes are never zero. */

	  cooling_ff = mplasma->cooling_ff = 0.0;
	  Error ("kpkt: A scattering event in cell %d with vol = 0???\n",
		 one->nwind);
	  //Diagnostic      return(-1);  //57g -- Cannot diagnose with an exit
	  exit (0);
	}


      if (cooling_ff < 0)
	{
	  Error ("kpkt: ff cooling rate negative. Abort.");
	  exit (0);
	}
      else
	{
	  cooling_normalisation += cooling_ff;
	}


      /* JM -- 1310 -- we now want to add adiabatic cooling as another way of destroying kpkts
         this should have already been calculated and stored in the plasma structure. Note that 
         adiabatic cooling does not depend on type of macro atom excited */

      /* note the units here- we divide the total luminosity of the cell by volume and ne to give cooling rate */

      cooling_adiabatic = xplasma->lum_adiabatic / one->vol / xplasma->ne;

      if (geo.adiabatic == 0 && cooling_adiabatic > 0.0)
        {
      	  Error("Adiabatic cooling turned off, but non zero in cell %d", xplasma->nplasma);
        }


      /* JM 1302 -- Negative adiabatic coooling- this used to happen due to issue #70, where we incorrectly calculated dvdy, 
         but this is now resolved. Now it should only happen for cellspartly in wind, because we don't treat these very well.
         Now, if cooling_adiabatic < 0 then set it to zero to avoid runs exiting for part in wind cells. */
      if (cooling_adiabatic < 0)
        {
          Error("kpkt: Adiabatic cooling negative! Major problem if inwind (%d) == 0\n",
          	     one->inwind);
	      Log ("kpkt: Setting adiabatic kpkt destruction probability to zero for this matom.\n");
	      cooling_adiabatic = 0.0;
        }


      cooling_normalisation += cooling_adiabatic;


 
      mplasma->cooling_bbtot = cooling_bbtot;
      mplasma->cooling_bftot = cooling_bftot;
      mplasma->cooling_bf_coltot = cooling_bf_coltot;
      mplasma->cooling_adiabatic = cooling_adiabatic;
      mplasma->cooling_normalisation = cooling_normalisation;
      mplasma->kpkt_rates_known = 1;

    }





  /* The cooling rates for the recombination and collisional processes are now known. 
     Choose which process destroys the k-packet with a random number. */

  destruction_choice =
    ((rand () + 0.5) / MAXRAND) * mplasma->cooling_normalisation;


  if (destruction_choice < mplasma->cooling_bftot)
    {				//destruction by bf
      for (i = 0; i < ntop_phot; i++)
	{
	  if (destruction_choice < mplasma->cooling_bf[i])
	    {
	      /* Having got here we know that desturction of the k-packet was via the process labelled
	         by i. Let's just check that i is a sensible number. */

	      if (i > ntop_phot - 1)
		{
		  Error
		    ("kpkt (matom.c): trying to destroy k-packet in unknown process. Abort.\n");
		  exit (0);
		}

	      /* If it gets here, all seems fine. Now set nres for the destruction process. */

	      *nres = i + NLINES + 1;
	      *escape = 1;	//record that an r-packet is made - no need to excite a macro atom again


	      /* Now (as in matom) choose a frequency for the new packet. */

	      p->freq = phot_top[i].freq[0]
		-
		(log (1. - (rand () + 0.5) / MAXRAND) * xplasma->t_e /
		 H_OVER_K);
	      /* Co-moving frequency - changed to rest frequency by doppler */
	      /*Currently this assumed hydrogenic shape cross-section - Improve */

	      /* k-packet is now eliminated. All done. */
	      return (0);
	    }
	  else
	    {
	      destruction_choice =
		destruction_choice - mplasma->cooling_bf[i];
	    }
	}
    }
  else if (destruction_choice <
	   (mplasma->cooling_bftot + mplasma->cooling_bbtot))
    {				//this means that a collisional destruction has occurred - this results in 
      //a macro atom being excited. Choose which macro atom and level to excite
      destruction_choice = destruction_choice - mplasma->cooling_bftot;
      for (i = 0; i < nlines; i++)
	{
	  if (destruction_choice < mplasma->cooling_bb[i])
	    {			//This is the bb collision which removes the k-packet
	      *nres = line[i].where_in_list;	//label for bb process 
	      if (line[i].macro_info == 1 && geo.macro_simple == 0)	//line is for a macro atom
		{

		  /* escape = 0 flag returned to tell macro_gov that
		     a macro atom should be excited, rather than making a call to matom here. */

		  *escape = 0;

		}
	      else		//line is not for a macro atom - use simple method
		{
		  /* Since the cooling rate accounts for the scattering fraction we know that if we
		     get here we want a line emission, not just an excited macro atom. (SS May 04) */
		  *escape = 1;	//No need for re-exciting a macro atom.
		  p->freq = line[i].freq;
		}
	      /* When it gets here the packet is back to an
	         r-packet and the emission mechanism is identified by nres
	         i.e. that's it finished. (SS, Apr 04). */
	      return (0);
	    }
	  else
	    {
	      destruction_choice =
		destruction_choice - mplasma->cooling_bb[i];
	    }
	}
    }
  else if (destruction_choice <
	   (mplasma->cooling_bftot +
	    mplasma->cooling_bbtot + mplasma->cooling_ff))
    {				//this is a ff destruction

      /* The limits for ff emission are hard-wired: 40 microns ->
         twice energy of He II edge. Shouldn't be a problem unless
         we're in plasma with temperatures that gives significant ff
         emission outside this range. */

      *escape = 1;		//we are making an r-packet not exciting a macro atom

      *nres = -2;

      p->freq = one_ff (one, 7.5e12, 2.626e16);	//get frequency of resulting energy packet

      return (0);
    }


  /* JM 1310 -- added loop to check if destruction occurs via adiabatic cooling */
  else if (destruction_choice <
	   (mplasma->cooling_bftot +
	    mplasma->cooling_bbtot + mplasma->cooling_ff + mplasma->cooling_adiabatic))
    {

      if (geo.adiabatic == 0)
      {
      	Error("Destroying kpkt by adiabatic cooling even though it is turned off.");
      }
      *escape = 1;		// we want to escape but set photon weight to zero
      *nres = -2;		
      //p->w = 0.0;		// JM131030 set photon weight to zero as energy is taken up in adiabatic expansion

      p->istat = P_ADIABATIC; 	// record that this photon went into a kpkt destruction from adiabatic cooling

      return (0);
    }


  else
    {
      /* We want destruction by collisional ionization in a macro atom. */
      destruction_choice =
	destruction_choice - mplasma->cooling_bftot -
	mplasma->cooling_bbtot - mplasma->cooling_ff - mplasma->cooling_adiabatic;

      for (i = 0; i < ntop_phot; i++)
	{
	  if (destruction_choice < mplasma->cooling_bf_col[i])
	    {
	      /* Having got here we know that destruction of the k-packet was via the process labelled
	         by i. Let's just check that i is a sensible number. */

	      if (i > ntop_phot - 1)
		{
		  Error
		    ("kpkt (matom.c): trying to destroy k-packet in unknown process. Abort.\n");
		  exit (0);
		}

	      /* Now set nres for the destruction process. */

	      *nres = i + NLINES + 1;
	      *escape = 0;	//collisional ionization makes an excited macro atom


	      /* k-packet is now eliminated. All done. */
	      return (0);
	    }
	  else
	    {
	      destruction_choice =
		destruction_choice - mplasma->cooling_bf_col[i];
	    }
	}
    }



  Error ("matom.c: Failed to select a destruction process in kpkt. Abort.\n");
  Error
    ("matom.c: cooling_bftot %g, cooling_bbtot %g, cooling_ff %g, cooling_bf_coltot %g cooling_adiabatic %g\n",
     mplasma->cooling_bftot, mplasma->cooling_bbtot,
     mplasma->cooling_ff, mplasma->cooling_bf_coltot, mplasma->cooling_adiabatic);

  exit (0);

  return (0);
}



/************************************************************
                                    Imperial College London
Synopsis:
       fake_matom_bb is the macro atom routine that deals with line events involving
       simple ions (i.e. ions for which a full macro atom treatment is not employed.
       When this routine is called a simple line has absorbed a packet. This routine 
       creates a fake two-level macro atom and determines whether the packet energy
       is simply re-emitted in the line or is thermalised. If it is thermalised it
       turns into a k-packet and the appropriate routine is called. 


Arguments:

       WindPtr w                   the ptr to the structure defining the wind
       PhotPtr p                   the packet at the point of activation
       int nres                    the process which activates the Macro Atom

Returns:
       int nres                    the process by which deactivation occurs
       PhotPtr p                   the packet following deactivation
       int escape                  identifies whether the macro atom deactivated via an
                                   r-packet (escape = 1) or a k-packet (escape = 0). If
                                   a k-packet then the call to this routine should be
                                   followed by a call to kpkt.

Description:


Notes: 

History:
          Apr 04  SS   Coding began.
          Jun 04  SS   Modified to return escape = 1 for r-packet and 2 for k-packet
                       to avoid call to k-packet within this routine.
	06may	ksl	57+ -- Eliminationg passing entire w structure

************************************************************/

int
fake_matom_bb (p, nres, escape)
     PhotPtr p;
     int *nres;
     int *escape;
{
  double kprb, rprb;
  WindPtr one;
  PlasmaPtr xplasma;
  struct lines *line_ptr;
  double electron_temperature;
  double normalisation;
  double choice;

  one = &wmain[p->grid];	//record where we are
  xplasma = &plasmamain[one->nplasma];

  line_ptr = lin_ptr[*nres];	//record the line pointer
  electron_temperature = xplasma->t_e;

  /* Upon calling we know that the upper level of our fake two level macro
     atom is excited. Since it's only two-levels there are no jumping probabilities
     and the only decision that needs to be made is whether it de-excites by an r-packet
     in the line or a collisionally produced k-packet. To decide this we compute the
     two emission probabilities. rprb is the probability of the r-packet and kprb
     for the k-packet. */

  /* Since both k- and r-packet emission have the same delta energy associated with them then
     we can ignore the energy factor that is otherwise needed in the de-activation probabilities.
     The upper level population is also factored out. (SS) */

  rprb = a21 (line_ptr) * p_escape (line_ptr, xplasma);
  // Einstein-A * escape probability

  kprb =
    q21 (line_ptr,
	 electron_temperature) * xplasma->ne * (1. -
						exp (-H_OVER_K *
						     line_ptr->freq /
						     electron_temperature));
  // Collision de-excitation rate coeff. * electron density 
  /* The extra factor of (1 - exp (-h nu / k T)) is explained in KSL's notes on Python. It appears because we are
     using the "scattering fraction" formalism for simple ions. (SS May 04). */

  normalisation = kprb + rprb;

  kprb = kprb / normalisation;

  rprb = rprb / normalisation;

  /* Now just use a random number to decide what happens. */

  choice = ((rand () + 0.5) / MAXRAND);

  /* If "choice" is less than rprb then we have chosen a radiative decay - for this fake macro atom there 
     is only one line so there's nothing to do - the energy is re-radiated in the line and that's it. We
     don't need to change nres. Otherwise we've chosen a collisional destruction and so we need to get
     kpkt to deal with the resulting k-packet and give us a new value of nres following thermalisation and
     re-emission of the energy. */

  if (choice < rprb)
    {
      *escape = 1;		//it's an r-packet de-activation
      p->freq = line_ptr->freq;	// As in matom, set this to the comoving frequency
    }
  else
    {
      /* Return the instruction to macro_gov to
         make a k-packet. */

      *escape = 0;

    }

  /* That's it - we've dealt with the simple line. */

  return (0);

}


/************************************************************
                                    Imperial College London
Synopsis:
       fake_matom_bf is the macro atom routine that deals with photoionisation 
       events involving simple ions (i.e. ions for which a full macro atom treatment 
       is not employed).
       When this routine is called a photoionisation has absorbed a packet. 
       The idea of this routine is to deal with the subsequenct
       by creating a fake two-level atom. 
       However, in the absense of collisional recombination (or something
       similar) there's only one thing that can happen - radiative 
       recombination. Therefore there's no need to do anything here unless
       collisional recombination is to be introduced at some point in the
       future.
       All this routine does for now is choose a new frequency for the 
       emitted r-packet.


Arguments:

       WindPtr w                   the ptr to the structure defining the wind
       PhotPtr p                   the packet at the point of activation
       int nres                    the process which activates the Macro Atom

Returns:
       int nres                    the process by which deactivation occurs
       PhotPtr p                   the packet following deactivation
       int escape                  in principle this tells us whether de-activation is
                                   via an r-packet or a k-packet. For this routine at the 
                                   moment only r-packets are possible so it always returns
                                   escape = 1

Description:


Notes:    

History:
          Apr 04  SS   Coding began.
          Jun 04  SS   Modified to include "escape" being set to 1
	06may	ksl	57+ -- Modified for new structure.  Have not
			fixed the call to fake_atom

************************************************************/

int
fake_matom_bf (p, nres, escape)
     PhotPtr p;
     int *nres;
     int *escape;
{
  WindPtr one;
  PlasmaPtr xplasma;

  one = &wmain[p->grid];	//record where we are
  xplasma = &plasmamain[one->nplasma];

  /* All that can happen is radiative recombination (no collisional recombination
     for now. Just choose frequency for re-emission). */

  *escape = 1;			//always an r-packet here

  p->freq = phot_top[*nres - NLINES - 1].freq[0]
    - (log (1. - (rand () + 0.5) / MAXRAND) * xplasma->t_e / H_OVER_K);

  /*Currently this assumes hydrogenic shape cross-section - Improve */

  return (0);

}



  /************************************************************
                                    Imperial College London
Synopsis:

       macro_pops uses the Monte Carlo estimators to compute a set
       of level populations for levels of macro atoms.


Arguments: one -> pointer to the wind 
           xne -> current value for electron density in this shell

Returns:   Should compute the fractional level populations for 
          macro atoms and store them in "levden" array. The ion fractions
          are also computed and stored in w[n].density[nion]


Description:


Notes:  
         This routine uses a matrix inversion method to get the level populations.
         For now the matrix solver used is the LU decomposition method provided
         by the Gnu Scientific Library (GSL, which is free). This requires the
         include files that I've added to the top of this file and access to the
         GSL "library" file (I've added the library into the Makefile too). I found
         GSL to be very easy to install but if there are problems in the future we
         may need to switch to another matrix solver. (SS, Apr 04)
         

History:
          Apr 04  SS   Coding began.
          May 04  SS   Checking/testing began. Corrected some problems with the
                       array indexing in the last part of the calculation.
          June04  SS   Adding collisional ionization/recombination rate.
	  Aug 05  SS   Modified to include fast collisional transition between
	               ground states and levels that should be metastable
                       - some treatment of metastables is needed for 
                       doing CNO macro atoms.
	06may	ksl	57+ -- Modified to reflect plasma structue
    1404 	jm  77a -- Added test to check that populations are solution to rate matrix equation in macro_pops
                       + more robust erro reporting in general. See pull request #75.
    1404  	jm  77a -- Now only clean for population inversions if levels have a radiative transition between them #76
************************************************************/

int
macro_pops (xplasma, xne)
     PlasmaPtr xplasma;
     double xne;
{

  int index_element, index_ion, index_lvl;
  int n_macro_lvl;
  double rate;
  double rate_matrix[NLEVELS_MACRO][NLEVELS_MACRO];
  int radiative_flag[NLEVELS_MACRO][NLEVELS_MACRO];	// 140423 JM flag if two levels are radiatively linked
  int conf_to_matrix[NLEVELS_MACRO];
  struct lines *line_ptr;
  struct topbase_phot *cont_ptr;
  int nn, mm;
  int index_bbu, index_bbd, index_bfu, index_bfd;
  int lower, upper;
  double this_ion_density, level_population;
  double inversion_test;
  double q_ioniz (), q_recomb ();
  int s;			//NEWKSL
  double *a_data, *b_data;
  gsl_permutation *p;		//NEWKSL
  gsl_matrix_view m;		//NEWKSL
  gsl_vector_view b;		//NEWKSL
  gsl_vector *populations, *test_vector;	//NEWKSL
  gsl_matrix *test_matrix;
  int index_fast_col, ierr;
  double test_val;

  MacroPtr mplasma;
  mplasma = &macromain[xplasma->nplasma];

  /* Start with an outer loop over elements: there are no rates that couple
     levels of different elements so we can always separate them out. */

  for (index_element = 0; index_element < nelements; index_element++)
    {

      /* Zero all elements of the matrix before doing anything else. */

      for (nn = 0; nn < NLEVELS_MACRO; nn++)
	{
	  for (mm = 0; mm < NLEVELS_MACRO; mm++)
	    {
	      rate_matrix[mm][nn] = 0.0;

	      /* 140423 JM -- new int array to flag if  two levels are radiatively linked
	         initialize to 0 */
	      radiative_flag[mm][nn] = 0; 
	    }
	}

      /* See if this element uses a macro atom treatment or is a simple element.
         For now I'm assuming that either all ions of a given element are
         treated using the macro atom method, or else none are (mixing and 
         matching is probably a bad idea because of the way in which bf 
         processes couple different ionisation stages). */

      /* The check is against the first ion of the element. */

      if (ion[ele[index_element].firstion].macro_info == 1
	  && geo.macro_simple == 0)
	{

	  /* Having established that the ion requires a macro atom treatment we
	     are going to construct a matrix of rates between the levels and 
	     invert that matrix to get the level populations. The first thing we need 
	     to do is work out how many levels we are dealing with in total. This is
	     easily done by summing up the number of levels of each ion. */

	  n_macro_lvl = 0;

	  for (index_ion = ele[index_element].firstion;
	       index_ion <
	       (ele[index_element].firstion + ele[index_element].nions);
	       index_ion++)
	    {
	      for (index_lvl = ion[index_ion].first_nlte_level;
		   index_lvl <
		   ion[index_ion].first_nlte_level + ion[index_ion].nlte;
		   index_lvl++)
		{
		  /* I want to be able to easily go from knowing the index of a level in the 
		     configurations structure to its position in the rates matrix. So I'm making
		     two arrays here that allow the mapping between these indices to be done easily.
		   */
		  conf_to_matrix[index_lvl] = n_macro_lvl;
		  n_macro_lvl++;
		}
	    }

	  /* We now know how many levels there are and therefore how big the matrix we
	     need to invert will be. */

	  /* Now we want to populate the matrix with all the rates between the levels. */

	  for (index_ion = ele[index_element].firstion;
	       index_ion <
	       (ele[index_element].firstion + ele[index_element].nions);
	       index_ion++)
	    {
	      index_lvl = ion[index_ion].first_nlte_level;

	      /* The next loop judges whether or not a level is to be fixed in population relative to ground
	         star. The input radiative lifetime is used to judge this at the moment. If the lifetime was set
	         to be long (essentially infite) then a very fast collisional transition is put in to dominate
	         all other rates into and out of this level. 

	         Whether this is really the best thing to do I don't know, but it's an improvement over ignoring 
	         this issue altogether! SS Aug 05 */

	      for (index_fast_col = index_lvl;
		   index_fast_col <
		   ion[index_ion].first_nlte_level + ion[index_ion].nlte - 1;
		   index_fast_col++)
		{
		  {
		    if (config[index_fast_col + 1].rad_rate > 1.e15)
		      {
			fast_line.gl = config[index_lvl].g;
			fast_line.gu = config[index_fast_col + 1].g;
			fast_line.freq =
			  (config[index_fast_col + 1].ex -
			   config[index_lvl].ex) / H;
			fast_line.f = 1e4;
			rate = q12 (&fast_line, xplasma->t_e) * xne;
			lower = conf_to_matrix[index_lvl];
			upper = conf_to_matrix[index_fast_col + 1];
			rate_matrix[lower][lower] += -1. * rate;
			rate_matrix[upper][lower] += rate;
			rate = q21 (&fast_line, xplasma->t_e) * xne;
			rate_matrix[upper][upper] += -1. * rate;
			rate_matrix[lower][upper] += rate;
		      }
		  }
		}

	      for (index_lvl = ion[index_ion].first_nlte_level;
		   index_lvl <
		   ion[index_ion].first_nlte_level + ion[index_ion].nlte;
		   index_lvl++)
		{
		  /* Now add contribution for all the bb and bf processes between the levels. This is
		     done by looping over the numbers "bbu, bbd, bfu, bfd" which tell us how many
		     processes there are. */


		  for (index_bbu = 0;
		       index_bbu < config[index_lvl].n_bbu_jump; index_bbu++)
		    {
		      /* These are bb upwards jumps. The rate in such a jump is given by
		         Jbar which has been computed as a Monte Carlo estimator. I'm also
		         including a collisional term (which depends on ne). */

		      line_ptr = &line[config[index_lvl].bbu_jump[index_bbu]];
		      rate =
			b12 (line_ptr) *
			mplasma->jbar_old[config[index_lvl].bbu_indx_first +
					  index_bbu];
		      rate += q12 (line_ptr, xplasma->t_e) * xne;

		      /* This is the rate out of the level in question. We need to add it
		         to the matrix in two places: firstly as a -ve contribution to the
		         diagonal and secondly as a +ve contribution for the off-diagonal
		         corresponding to the level populated by this process. */

		      /* Get the matix indices for the upper and lower states of the jump. */

		      lower = conf_to_matrix[index_lvl];
		      upper = conf_to_matrix[line_ptr->nconfigu];

		      rate_matrix[lower][lower] += -1. * rate;
		      rate_matrix[upper][lower] += rate;

		      /* There's a radiative jump between these levels, so we want to clean
		         for popualtion inversions. Flag this jump */
		      radiative_flag[index_lvl][line_ptr->nconfigu] = 1;

		      if (rate < 0.0 || sane_check(rate))
		      {
		      	Error("macro_pops: bbu rate is %8.4e in cell/matom %i\n", rate, xplasma->nplasma);
		      }
		    }

		  for (index_bbd = 0;
		       index_bbd < config[index_lvl].n_bbd_jump; index_bbd++)
		    {
		      /* These are bb downwards jumps. The rate in such a jump is given by
		         the A-value. I'm also
		         including a collisional term (which depends on ne). */

		      line_ptr = &line[config[index_lvl].bbd_jump[index_bbd]];
		      rate = (a21 (line_ptr) * p_escape (line_ptr, xplasma));
		      rate += q21 (line_ptr, xplasma->t_e) * xne;

		      /* This is the rate out of the level in question. We need to add it
		         to the matrix in two places: firstly as a -ve contribution to the
		         diagonal and secondly as a +ve contribution for the off-diagonal
		         corresponding to the level populated by this process. */

		      /* Get the matix indices for the upper and lower states of the jump. */

		      upper = conf_to_matrix[index_lvl];
		      lower = conf_to_matrix[line_ptr->nconfigl];

		      rate_matrix[upper][upper] += -1. * rate;
		      rate_matrix[lower][upper] += rate;


		      /* There's a radiative jump between these levels, so we want to clean
		         for popualtion inversions. Flag this jump */
		      radiative_flag[line_ptr->nconfigl][index_lvl] = 1;

		      if (rate < 0.0 || sane_check(rate))
		      {
		      	Error("macro_pops: bbd rate is %8.4e in cell/matom %i\n", rate, xplasma->nplasma);
		      }
		    }



		  for (index_bfu = 0;
		       index_bfu < config[index_lvl].n_bfu_jump; index_bfu++)
		    {
		      /* These are bf upwards jumps. The rate in such a jump is given by
		         gamma which has been computed as a Monte Carlo estimator. */

		      cont_ptr =
			&phot_top[config[index_lvl].bfu_jump[index_bfu]];
		      rate =
			mplasma->gamma_old[config[index_lvl].bfu_indx_first +
					   index_bfu];
		      rate += q_ioniz (cont_ptr, xplasma->t_e) * xne;

		      /* This is the rate out of the level in question. We need to add it
		         to the matrix in two places: firstly as a -ve contribution to the
		         diagonal and secondly as a +ve contribution for the off-diagonal
		         corresponding to the level populated by this process. */

		      /* Get the matix indices for the upper and lower states of the jump. */

		      lower = conf_to_matrix[index_lvl];
		      upper = conf_to_matrix[cont_ptr->uplev];

		      rate_matrix[lower][lower] += -1. * rate;
		      rate_matrix[upper][lower] += rate;

		      if (rate < 0.0 || sane_check(rate))
		      {
		      	Error("macro_pops: bfu rate is %8.4e in cell/matom %i\n", rate, xplasma->nplasma);
		      }

		      /* Now deal with the stimulated emission. */
		      /* Lower and upper are the same, but now it contributes in the
		         other direction. */

		      rate =
			mplasma->
			alpha_st_old[config[index_lvl].bfu_indx_first +
				     index_bfu] * xne;

		      rate_matrix[upper][upper] += -1. * rate;
		      rate_matrix[lower][upper] += rate;

		      if (rate < 0.0 || sane_check(rate))
		      {
		      	Error("macro_pops: st. recomb rate is %8.4e in cell/matom %i\n", rate, xplasma->nplasma);
		      }

		    }



		  for (index_bfd = 0;
		       index_bfd < config[index_lvl].n_bfd_jump; index_bfd++)
		    {
		      /* These are bf downwards jumps. The rate in such a jump is given by
		         the alpha value. */
		      cont_ptr =
			&phot_top[config[index_lvl].bfd_jump[index_bfd]];
		      /* Get new values of the recombination rates and store them. */
		      mplasma->recomb_sp[config[index_lvl].bfd_indx_first +
					 index_bfd] =
			alpha_sp (cont_ptr, xplasma, 0);
		      mplasma->recomb_sp_e[config[index_lvl].bfd_indx_first +
					   index_bfd] =
			alpha_sp (cont_ptr, xplasma, 2);
		      rate =
			mplasma->recomb_sp[config[index_lvl].bfd_indx_first +
					   index_bfd] * xne;
		      rate += q_recomb (cont_ptr, xplasma->t_e) * xne * xne;


		      /* This is the rate out of the level in question. We need to add it
		         to the matrix in two places: firstly as a -ve contribution to the
		         diagonal and secondly as a +ve contribution for the off-diagonal
		         corresponding to the level populated by this process. */

		      /* Get the matix indices for the upper and lower states of the jump. */

		      upper = conf_to_matrix[index_lvl];
		      lower = conf_to_matrix[cont_ptr->nlev];

		      rate_matrix[upper][upper] += -1. * rate;
		      rate_matrix[lower][upper] += rate;

		      if (rate < 0.0 || sane_check(rate))
		      {
		      	Error("macro_pops: bfd rate is %8.4e in cell/matom %i\n", rate, xplasma->nplasma);
		      }
		    }
		}
	    }

	  /* The rate matrix is now filled up. Since the problem is not closed as it stands, the next
	     thing is to replace one of the rows of the matrix (say the first row) with the constraint
	     that the sum of all the populations is 1.0 (this will let us get population fractions). */

	  for (index_lvl = 0; index_lvl < n_macro_lvl; index_lvl++)
	    {
	      rate_matrix[0][index_lvl] = 1.0;
	    }

	  /* Now we can just invert the matrix to get the fractional level populations. */


	  /********************************************************************************/
	  /* The block that follows (down to next line of ***s) is to do the
	     matrix inversion. It uses LU decomposition - the code for doing this is 
	     taken from the GSL manual with very few modifications. */
	  /* here we solve the matrix equation M x = b, where x is our vector containing
	     level populations as a fraction w.r.t the whole element */

	  /* Replaced inline array allocaation with calloc, which will work with older version of c compilers */

	  a_data =
	    (double *) calloc (sizeof (rate), n_macro_lvl * n_macro_lvl);

	  for (nn = 0; nn < n_macro_lvl; nn++)
	    {
	      for (mm = 0; mm < n_macro_lvl; mm++)
		{
		  a_data[nn * n_macro_lvl + mm] = rate_matrix[nn][mm];
		}
	    }


	  /* Replaced inline array allocaation with calloc, which will work with older version of c compilers 
	     calloc also sets the elements to zero, which is required */

	  b_data = (double *) calloc (sizeof (rate), n_macro_lvl);

	  /* replace the first entry with 1.0- this is part of the normalisation constraint */
	  b_data[0] = 1.0;

      /* create gsl matrix/vector views of the arrays of rates */
	  m = gsl_matrix_view_array (a_data, n_macro_lvl, n_macro_lvl);	//KSLNEW

	  /* these are used for testing the solution below */
	  test_matrix = gsl_matrix_alloc(n_macro_lvl, n_macro_lvl);
	  test_vector = gsl_vector_alloc (n_macro_lvl);

	  gsl_matrix_memcpy(test_matrix, &m.matrix); 	// create copy for testing 

	  b = gsl_vector_view_array (b_data, n_macro_lvl);	//KSLNEW

	  /* the populations vector will be a gsl vector which stores populations */
	  populations = gsl_vector_alloc (n_macro_lvl);	//KSLNEW


	  p = gsl_permutation_alloc (n_macro_lvl);	//NEWKSL

	  gsl_linalg_LU_decomp (&m.matrix, p, &s);

	  gsl_linalg_LU_solve (&m.matrix, p, &b.vector, populations);

	  gsl_permutation_free (p);


	  /**********************************************************/

      /* JM 140414 -- before we clean, we should check that the populations vector we 
         have just created really is a solution to the matrix equation */
      
      /* declaration contained in my_linalg.h, taken from gsl library */
      ierr = gsl_blas_dgemv(CblasNoTrans, 1.0, test_matrix, populations, 0.0, test_vector);

      if (ierr != 0)
        {
      	  Error("macro_pops: bad return when testing matrix solution to rate equations.\n");
        }

      /* now cycle through and check the solution to y = m * populations really is
         (1, 0, 0 ... 0) */

      for (mm = 0; mm < n_macro_lvl; mm++)
        {

          /* get the element of the vector we want to check */
          test_val = gsl_vector_get (test_vector, mm);
          
          if (mm == 0)		// 0th element should be 1
          {
      	    if ( (fabs(test_val - 1.0)) > EPSILON)
      	      Error("macro_pops: test vector has first element = %8.4e, should be 1.0\n",
      	      	     test_val);    
          }     

          else 				// all other elements should be zero
          {
          	if ( fabs(test_val)  > EPSILON)
      	      Error("macro_pops: test vector has element %i = %8.4e, should be 0.0\n",
      	      	     mm, test_val);    
          }
        }

        /* free memory */
        free (a_data);	
	    free (b_data);	
	    free (test_vector);
	    free (test_matrix);



	  /* MC noise can cause population inversions (particularly amongst highly excited states)
	     which are never a good thing and most likely unphysical.
	     Therefor let's follow Leon's procedure (Lucy 2003) and remove inversions. */



	  for (index_ion = ele[index_element].firstion;
	       index_ion <
	       (ele[index_element].firstion + ele[index_element].nions);
	       index_ion++)
	    {
	      for (index_lvl = ion[index_ion].first_nlte_level;
		   index_lvl <
		   ion[index_ion].first_nlte_level + ion[index_ion].nlte;
		   index_lvl++)
		{		/* Start loop with lowest level of the ion. For each level in turn check to see if there's a population 
				   inversion i.e. is  upper_pop > lower_pop * g_upper / g_lower. If it is then replace upper_pop with
				   lower_pop * g_upper / g_lower. We loop over all levels higher than the currently chosen lower level. */
		  for (nn = index_lvl + 1;
		       nn <
		       (ion[index_ion].first_nlte_level +
			ion[index_ion].nlte); nn++)
		    {

              /* this if statement means we only clean if there's a radiative jump between the levels */
              if (radiative_flag[index_lvl][nn])
              {
		        inversion_test = gsl_vector_get (populations, conf_to_matrix[index_lvl]) * config[nn].g / config[index_lvl].g * 0.999999;	//include a correction factor 
		        if (gsl_vector_get (populations, conf_to_matrix[nn]) >
			    inversion_test)
			  {
			    gsl_vector_set (populations, conf_to_matrix[nn],
					  inversion_test);
			  }
		      }
		    }
		}
	    }


	  /* The populations are now known. The populations need to be stored
	     firstly as ion populations and secondly as fractional
	     level populations within an ion. Get the ion
	     populations and write them to one->density[nion]. The level populations
	     are to be put in "levden". */

	  nn = 0;
	  mm = 0;
	  for (index_ion = ele[index_element].firstion;
	       index_ion <
	       (ele[index_element].firstion + ele[index_element].nions);
	       index_ion++)
	    {
	      this_ion_density = 0.0;
	      for (index_lvl = ion[index_ion].first_nlte_level;
		   index_lvl <
		   ion[index_ion].first_nlte_level + ion[index_ion].nlte;
		   index_lvl++)
		{
		  level_population = gsl_vector_get (populations, conf_to_matrix[index_lvl]);	
		  this_ion_density += level_population;

		  /* JM 140409 -- add error check for negative populations */
		  if (level_population < 0.0 || sane_check(level_population))
		    Error("macro_pops: level %i has frac. pop. %8.4e in cell %i\n", 
		    	   index_lvl, level_population, xplasma->nplasma);

		  nn++;
		}

	      /* Write the level populations to the 
	         levden array. These are fractional level populations within an ion. */

	      for (index_lvl = ion[index_ion].first_nlte_level;
		   index_lvl <
		   ion[index_ion].first_nlte_level + ion[index_ion].nlte;
		   index_lvl++)
		{

		  xplasma->levden[config[index_lvl].nden] =
		    gsl_vector_get (populations,
				    conf_to_matrix[index_lvl]) /
		    this_ion_density;

		  mm++;

		}

	      xplasma->density[index_ion] =
		this_ion_density * ele[index_element].abun * xplasma->rho *
		rho2nh;




	    }

	  gsl_vector_free (populations);
	}
    }

  return (0);
  /* All done. (SS, Apr 04) */
}

  /************************************************************
                                    Imperial College London
Synopsis:

       macro_gov is a new routine that will sit at a higher level in the code than either matom or kpkt
       and will govern the passage of packets between these routines. At the moment, since matom and kpkt 
       call each other it call all get rather confusing and loads of nested subroutine calls ensue. 
       macro_gov removes this by calling matom and kpkt which on return tell macro_gov to either return
       an r-packet to resonate or make another call to either kpkt or matom as appropriate.



Arguments:   
       WindPtr w                   the ptr to the structure defining the wind
       PhotPtr p                   the packet at the point of activation
       int nres                    the process which activates the Macro Atom
       int matom_or_kpkt           to tell us if we should initially excite a
                                   Macro Atom or create a k-packet. For Macro Atom
                                   it is set to 1 for a k-packet it is set to 2.
       int which_out               set to 1 if return is via macro atom and 2 if via kpkt

           

Returns:   Will return an r-packet after (possibly) several calls to matom and kpkt

       int nres                    the process by which re-emission occurs
       PhotPtr p                   the packet following re-emission



Description:


Notes:   I've written this to be as general as possible so that if we want to improve the treatment of simple
         ions it should not need to be re-written.

         During the spectrum calculation the emission of r-packets within the spectral region of interest
         is done using emissivities obtained during the ionization cycles. Therefore whenever an r-packet
         is converted into a k-packet or an excited macro atom that ends the need to follow the packet any
         further. To deal with this, this routine sets the weights of such packets to zero. Upon return to 
         trans_phot these packets will then be thrown away.

         
History:
          SS   June 04 - coding began
          SS   June 04 - modified to kill photons that interact with macro atoms or
                         become k-packets during the specrum calculation part of the code.
          SS   June 04 - returns 1 if last interaction was with full macro atom. 0 otherwise.
          SS   Sep  04 - added return of "which_out"
	  ksl	0803	Added lines which modify p->origin to reflect the fact
	  		that a photn has been processed by an interaction with a
			macro atom.  Basically, what I have done is simply to 
			add 10 to the origin.  This was done to account for the fact
			that its hard to constrain the frequence range of a photon
			processed by a macro atom
************************************************************/

int
macro_gov (p, nres, matom_or_kpkt, which_out)
     PhotPtr p;
     int *nres;
     int matom_or_kpkt;
     int *which_out;

{
  int escape;			//this tells us when the r-packet is escaping 

  escape = 0;			//start with it not being ready to escape as an r-packet


  while (escape == 0)
    {
      if (matom_or_kpkt == 1)	//excite a macro atom - depending on simple/macro choice call different routines
	{

	  if (*nres > (-1) && *nres < NLINES && geo.macro_simple == 0
	      && lin_ptr[*nres]->macro_info == 1)
	    /* This is a bb line of a macro atoms and we want the full macro atom treatment. */
	    {

	      if (geo.matom_radiation == 1)
		{
		  /* During the spectrum cycles we want to throw these photons away. */
		  p->w = 0.0;
		  escape = 1;	//This doesn't matter but it breaks us out of this loop
		}
	      else
		{
		  matom (p, nres, &escape);

		  if (escape == 1)
		    {
		      /* It is going to escape as a r-packet that was created by de-activation of a macro atom. 
		         Therefore, if the frequency is suitable it should be recorded as a macro atom emission event for use in the 
		         computation of the k-packet emissivity needed for the final spectrum calculation. */
		      *which_out = 1;
		      /* 0803 - ksl - 60 - Added code to modify the photon origin to indicate the packet has been processed
		         by a macro atom */
		      if (p->origin < 10)
			p->origin += 10;
		      return (1);
		    }
		}
	    }
	  else if (*nres > (-1) && *nres < NLINES
		   && (geo.macro_simple == 1
		       || lin_ptr[*nres]->macro_info == 0))
	    /* This is a bb line for which we don't want a full macro atom treatment. */
	    {
	      fake_matom_bb (p, nres, &escape);
	    }
	  else if (*nres > NLINES
		   && phot_top[*nres - NLINES - 1].macro_info == 1
		   && geo.macro_simple == 0)
	    /* This is a bf continuum of a macro atom and we want the full treatment. */
	    {

	      if (geo.matom_radiation == 1)
		{
		  /* During the spectrum cycles we want to throw these photons away. */
		  p->w = 0.0;
		  escape = 1;	//This doesn't matter but it breaks us out of this loop
		}
	      else
		{
		  matom (p, nres, &escape);

		  if (escape == 1)
		    {
		      /* It is going to escape as a r-packet that was created by de-activation of a macro atom. 
		         Therefore, if the frequency is suitable it should be recorded as a macro atom emission event for use in the 
		         computation of the k-packet emissivity needed for the final spectrum calculation. */
		      *which_out = 1;
		      /* 0803 - ksl - 60 - Added code to modify the photon origin to indicate the packet has been processed
		         by a macro atom */
		      if (p->origin < 10)
			p->origin += 10;
		      return (1);
		    }
		}
	    }
	  else if (*nres > NLINES
		   && (phot_top[*nres - NLINES - 1].macro_info == 0
		       || geo.macro_simple == 1))
	    /* This is a bf continuum but we don't want the full macro atom treatment. */
	    {
	      fake_matom_bf (p, nres, &escape);
	    }

	  matom_or_kpkt = 2;	//if it did not escape then it must have had a
	  //de-activation by collision processes -> need a k-packet 
	}
      else if (matom_or_kpkt == 2)	//deal with a k-packet
	{
	  if (geo.matom_radiation == 1)
	    {
	      /* During the spectrum cycles we want to throw these photons away. */
	      p->w = 0.0;
	      escape = 1;	//This doesn't matter but it breaks us out of the loop.
	    }
	  else
	    {
	      kpkt (p, nres, &escape);

	    }

	  matom_or_kpkt = 1;	//if it did not escape then the k-packet must have been 
	  //destroyed by collisionally exciting a macro atom - 
	  //excite a macro atom next
	}
      else
	{
	  Error ("macro_gov: Unknown choice for next action. Abort.\n");
	  exit (0);
	}
    }

  //When it gets here an escpae has taken place. That's it done.

  *which_out = 2;

  /* 0803 - ksl - 60 - Added code to modify the photon origin to indicate the packet has been processed
     by a macro atom */

  if (p->origin < 10)
    p->origin += 10;

  return (0);

  //All done. 
}





/************************************************************
                                    Imperial College London
Synopsis:

       get_kpkt_f returns the specific luminosity in the band needed for the computation of the
       spectrum. It gets the total energy radiated by the process k-packet -> r-packet in the
       required wavelength range.


Arguments:   
       WindPtr w                   the ptr to the structure defining the wind

                                  
           

Returns:   The energy radiated by the process k-packet -> r-packet in the wind in the 
           wavelength range required for the specrum calculation.



Description:


Notes:  
        
         
History:
          June 04 SS - coding began
	06may	ksl	57+ -- Recoded to use plasma array

************************************************************/

double
get_kpkt_f ()
{
  int n;
  double lum;

  lum = 0.0;

  for (n = 0; n < NPLASMA; n++)
    {
      lum += plasmamain[n].kpkt_emiss;
    }


  return (lum);
}

/* All done. */

/************************************************************
                                    Imperial College London
Synopsis:

       get_matom_f returns the specific luminosity in the band needed for the computation of the
       spectrum. It gets the total energy radiated by the deactivation of macro atoms in the
       required wavelength range.


Arguments:   
       WindPtr w                   the ptr to the structure defining the wind

                                              

Returns:   The energy radiated by the deactivation of macro atoms in the wind in the 
           wavelength range required for the specrum calculation.



Description:


Notes:  
        
         
History:
          June 04 SS - coding began
          Sep  04 SS - significant modification to improve the treatment of macro
                       atoms in spectral synthesis steps. 
	06may	ksl	57+ -- Recoded to use plasma structure

************************************************************/

double
get_matom_f ()
{
  int n, m;
  int mm, ss;
  double lum;
  int level_emit[NLEVELS_MACRO], kpkt_emit;
  int n_tries, n_tries_local;
  struct photon ppp;
  double contribution, norm;
  int nres, which_out;
  int my_nmin, my_nmax;		//These variables are used even if not in parallel mode

#ifdef MPI_ON
  int num_mpi_cells, num_mpi_extra, position, ndo, n_mpi, num_comm, n_mpi2;
  int size_of_commbuffer;
  char *commbuffer;

  /* the commbuffer needs to communicate 2 variables and the number of macor levels, 
     plus the variable for how many cells each thread is doing */
  size_of_commbuffer =
    8 * (3 + nlevels_macro) * (floor (NPLASMA / np_mpi_global) + 1);

  commbuffer = (char *) malloc (size_of_commbuffer * sizeof (char));
#endif



  which_out = 0;
  lum = 0.0;
  n_tries = 5000000;
  geo.matom_radiation = 0;
  n_tries_local = 0;



  norm = 0;
  for (n = 0; n < NPLASMA; n++)
    {
      for (m = 0; m < nlevels_macro; m++)
	{
	  norm += macromain[n].matom_abs[m];
	}
      norm += plasmamain[n].kpkt_abs;
    }



  Log ("Calculating macro atom emissivities- this might take a while...\n");

  /* For MPI parallelisation, the following loop will be distributed over multiple tasks. 
     Note that the mynmim and mynmax variables are still used even without MPI on */
  my_nmin = 0;
  my_nmax = NPLASMA;

#ifdef MPI_ON

  num_mpi_cells = floor(NPLASMA/np_mpi_global);  // divide the cells between the threads
  num_mpi_extra = NPLASMA - (np_mpi_global*num_mpi_cells);  // the remainder from the above division
  
  /* this next loop splits the cells up between the threads. All threads with 
     rank_global<num_mpi_extra deal with one extra cell to account for the remainder */
  if (rank_global < num_mpi_extra)
    {
      my_nmin = rank_global*(num_mpi_cells+1);
      my_nmax = (rank_global+1)*(num_mpi_cells+1);     
    }
  else
    {
      my_nmin = num_mpi_extra*(num_mpi_cells+1) + (rank_global-num_mpi_extra)*(num_mpi_cells);
      my_nmax = num_mpi_extra*(num_mpi_cells+1) + (rank_global-num_mpi_extra+1)*(num_mpi_cells);
    }
  ndo = my_nmax-my_nmin;

  Log_parallel
    ("Thread %d is calculating macro atom emissivities for macro atoms %d to %d\n",
     rank_global, my_nmin, my_nmax);

#endif




  for (n = my_nmin; n < my_nmax; n++)
    {

      /* JM 1309 -- these lines are just log statements which track progress, as this section
         can take a long time */
#ifdef MPI_ON
      if (n % 50 == 0)
	Log
	  ("Thread %d is calculating  macro atom emissivity for macro atom %7d of %7d or %6.3f per cent\n",
	   rank_global, n, my_nmax, n * 100. / my_nmax);
#else
      if (n % 50 == 0)
	Log
	  ("Calculating macro atom emissivity for macro atom %7d of %7d or %6.3f per cent\n",
	   n, my_nmax, n * 100. / my_nmax);
#endif

      for (m = 0; m < nlevels_macro + 1; m++)
	{
	  if ((m == nlevels_macro && plasmamain[n].kpkt_abs > 0)
	      || (m < nlevels_macro && macromain[n].matom_abs[m] > 0))
	    {
	      if (m < nlevels_macro)
		{
		  if (macromain[n].matom_abs[m] > 0)
		    {
		      n_tries_local =
			(n_tries * macromain[n].matom_abs[m] / norm) + 10;
		    }
		  else
		    {
		      n_tries_local = 0;
		    }
		}
	      else if (m == nlevels_macro)
		{
		  if (plasmamain[n].kpkt_abs > 0)
		    {
		      n_tries_local =
			(n_tries * plasmamain[n].kpkt_abs / norm) + 10;
		    }
		  else
		    {
		      n_tries_local = 0;
		    }
		}
	      /* We know "matom_abs" is the amount of energy absobed by
	         each macro atom level in each cell. We now want to determine what fraction of
	         that energy re-appears in the frequency range we want and from which macro atom
	         level it is re-emitted (or if it appears via a k-packet). */

	      for (mm = 0; mm < NLEVELS_MACRO; mm++)
		{
		  level_emit[mm] = 0;
		}
	      kpkt_emit = 0;
	      if (n_tries_local > 0)
		{
		  for (ss = 0; ss < n_tries_local; ss++)
		    {
		      if (m < nlevels_macro)
			{
			  /* Dealing with excitation of a macro atom level. */
			  /* First find a suitable transition in which to excite the macro atom (to get it
			     in the correct starting level. */

			  /* As Knox pointed out in an earlier comment
			     here, the next loop isn't the best way to
			     do this - however I'm not sure we can
			     just replace nlines->nlines_macro since
			     we don't know for sure that the
			     macro_lines are the first lines in the
			     array. If this is really a bottle neck we
			     should set up a quicker way to identify
			     the level - the purpose is just to set
			     the macro atom state to the appropriate
			     upper level so all we need is some way to
			     tell it which state it should be. */

			  nres = 0;
			  while ((nres < nlines))
			    {
			      if (lin_ptr[nres]->nconfigu != m)
				{
				  nres += 1;
				}
			      else
				{
				  break;
				}
			    }


			  if (nres > (nlines - 1))
			    {
			      nres = NLINES + 1;
			      while (nres < (NLINES + 1 + ntop_phot))
				{
				  if (phot_top[nres - NLINES - 1].uplev != m)
				    {
				      nres += 1;
				    }
				  else
				    {
				      break;
				    }
				}
			    }

			  if (nres > NLINES + ntop_phot)
			    {
			      Error ("Problem in get_matom_f (1). Abort. \n");
			      exit (0);
			    }

			  ppp.nres = nres;
			  ppp.grid = plasmamain[n].nwind;
			  ppp.w = 0;

			  macro_gov (&ppp, &nres, 1, &which_out);


			  /* Now a macro atom has been excited and followed until an r-packet is made. Now, if that 
			     r-packet is in the correct frequency range we record it. If not, we throw it away. */
			}
		      else if (m == nlevels_macro)
			{
			  /* kpkt case. */

			  nres = -2;	//will do - just need
			  //something that will tigger kpkt

			  ppp.nres = nres;
			  ppp.grid = plasmamain[n].nwind;
			  ppp.w = 0;

			  macro_gov (&ppp, &nres, 2, &which_out);

			  /* We have an r-packet back again. */
			}


		      if (ppp.freq > em_rnge.fmin && ppp.freq < em_rnge.fmax)
			{
			  if (which_out == 1)
			    {
			      if (nres < 0)
				{
				  Error
				    ("Negative out from matom?? Abort.\n");
				  exit (0);
				}

			      /* It was a macro atom de-activation. */
			      if (nres < nlines)
				{	/* It was a line. */
				  level_emit[lin_ptr[nres]->nconfigu] += 1;
				}
			      else
				{
				  level_emit[phot_top
					     [nres - NLINES - 1].uplev] += 1;
				}
			    }
			  else if (which_out == 2)
			    {
			      /* It was a k-packet de-activation. */
			      kpkt_emit += 1;
			    }
			  else
			    {
			      Error
				("Packet didn't emerge from matom or kpkt??? Abort. \n");
			      exit (0);
			    }
			}
		    }


		  /* Now we've done all the runs for this level de-activating so we can add the contributions
		     to the level emissivities, the k-packet emissivity and the total luminosity. */

		  for (mm = 0; mm < nlevels_macro; mm++)
		    {
		      contribution = 0;
		      if (m < nlevels_macro)
			{
			  macromain[n].matom_emiss[mm] += contribution =
			    level_emit[mm] * macromain[n].matom_abs[m] /
			    n_tries_local;
			}
		      else if (m == nlevels_macro)
			{
			  macromain[n].matom_emiss[mm] += contribution =
			    level_emit[mm] * plasmamain[n].kpkt_abs /
			    n_tries_local;
			}
		      lum += contribution;
		    }

		  if (m < nlevels_macro)
		    {
		      plasmamain[n].kpkt_emiss +=
			kpkt_emit * macromain[n].matom_abs[m] / n_tries_local;
		    }
		  else if (m == nlevels_macro)
		    {
		      plasmamain[n].kpkt_emiss +=
			kpkt_emit * plasmamain[n].kpkt_abs / n_tries_local;
		    }
		}
	    }
	}
    }


  /*This is the end of the update loop that is parallelised. We now need to exchange data between the tasks.
    This is done much the same way as in wind_update */
#ifdef MPI_ON
  for (n_mpi = 0; n_mpi < np_mpi_global; n_mpi++)
    {
      /* here we loop over the number of threads. If the thread is this thread then we pack the macromain information
         into memory and then broadcast it to ther other threads */
      position = 0;

      if (rank_global == n_mpi)
	{

	  Log
	    ("MPI task %d is working on matoms %d to max %d (total size %d).\n",
	     rank_global, my_nmin, my_nmax, NPLASMA);

	  MPI_Pack (&ndo, 1, MPI_INT, commbuffer,
		    size_of_commbuffer, &position, MPI_COMM_WORLD);
	  for (n = my_nmin; n < my_nmax; n++)
	    {

	      /* pack the number of the cell, and the kpkt and macro atom emissivites for that cell // */
	      MPI_Pack (&n, 1, MPI_INT, commbuffer,
			size_of_commbuffer, &position, MPI_COMM_WORLD);
	      MPI_Pack (&plasmamain[n].kpkt_emiss, 1, MPI_DOUBLE,
			commbuffer, size_of_commbuffer, &position,
			MPI_COMM_WORLD);
	      MPI_Pack (macromain[n].matom_emiss, nlevels_macro, MPI_DOUBLE,
			commbuffer, size_of_commbuffer, &position,
			MPI_COMM_WORLD);

	    }

	  Log_parallel
	    ("MPI task %d broadcasting matom emissivity information.\n",
	     rank_global);
	}


      /* Set MPI_Barriers and broadcast information to other threads */
      MPI_Barrier (MPI_COMM_WORLD);
      MPI_Bcast (commbuffer, size_of_commbuffer, MPI_PACKED,
		 n_mpi, MPI_COMM_WORLD);
      MPI_Barrier (MPI_COMM_WORLD);
      Log_parallel
	("MPI task %d survived broadcasting matom emissivity information.\n",
	 rank_global);



      position = 0;

      /* If not this thread then we unpack the macromain information from the other threads */

      if (rank_global != n_mpi)
	{
	  MPI_Unpack (commbuffer, size_of_commbuffer, &position,
		      &num_comm, 1, MPI_INT, MPI_COMM_WORLD);
	  for (n_mpi2 = 0; n_mpi2 < num_comm; n_mpi2++)
	    {

	      /* unpack the number of the cell, and the kpkt and macro atom emissivites for that cell */
	      MPI_Unpack (commbuffer, size_of_commbuffer,
			  &position, &n, 1, MPI_INT, MPI_COMM_WORLD);
	      MPI_Unpack (commbuffer, size_of_commbuffer, &position,
			  &plasmamain[n].kpkt_emiss, 1, MPI_DOUBLE,
			  MPI_COMM_WORLD);
	      MPI_Unpack (commbuffer, size_of_commbuffer, &position,
			  macromain[n].matom_emiss, nlevels_macro, MPI_DOUBLE,
			  MPI_COMM_WORLD);
	    }
	}
    }


  /* this next loop just corrects lum to be the correct summed value in parallel mode */
  lum = 0.0;   // need to rezero, fixes segfault bug #59 

  for (n = 0; n < NPLASMA; n++)
    {

      for (mm = 0; mm < nlevels_macro; mm++)
	{
	  lum += macromain[n].matom_emiss[mm];
	}
    }
#endif



  geo.matom_radiation = 1;

  return (lum);
}

/* All done. */



/************************************************************
                                    Imperial College London
Synopsis:

     photo_gen_kpkt produces photon packets to account for creating of r-packets
     by k-packets in the spectrum calculation. It should only be used once the total 
     energy emitted in this way in the wavelength range in question is well known
     (calculated in the ionization cycles).


Arguments:   
       WindPtr w                   the ptr to the structure defining the wind
       PhotPtr p;                  the ptr to the structire for the photons
       double weight;              the photon weight
       int photstart, nphot;       the number of the first photon to be generated and
                                   the total number of photons to be generated by this
                                   routine

                                           

Returns:   
      When it finishes it should have generated nphot photons from k-packet elliminations.


Description:
      This routine is closely related to photo_gen_wind from which much of the code has been copied.

Notes:  
        
         
History:
          June 04 SS - coding began
	  04aug	ksl	Modified so that uses coordinate system
	  		independent routine get_random_location
			to set location of generated kpkt
          Sep  04 SS - minor bug fixed
	06may	ksl	57+ -- Updated to use plasma structure. 
			and to elimate passing entrie wind structure

************************************************************/
int
photo_gen_kpkt (p, weight, photstart, nphot)
     PhotPtr p;
     double weight;
     int photstart, nphot;
{
  int photstop;
  int icell;
  double xlum, xlumsum;
  struct photon pp;
  int nres, esc_ptr;
  int n;
  double v[3];
  double dot ();
  int wind_n_to_ij (), stuff_v (), randvec ();
  int get_random_location ();
  double test;
  int nnscat;
  double dvwind_ds (), sobolev ();
  int nplasma;



  photstop = photstart + nphot;
  Log ("photo_gen_kpkt creates nphot %5d photons from %5d to %5d \n", nphot,
       photstart, photstop);

  for (n = photstart; n < photstop; n++)
    {
      /* locate the wind_cell in which the photon bundle originates. */

      xlum = (rand () + 0.5) / (MAXRAND) * geo.f_kpkt;

      xlumsum = 0;
      icell = 0;
      while (xlumsum < xlum)
	{
	  if (wmain[icell].vol > 0.0)
	    {
	      nplasma = wmain[icell].nplasma;
	      xlumsum += plasmamain[nplasma].kpkt_emiss;
	    }
	  icell++;
	}
      icell--;			/* This is the cell in which the photon must be generated */

      /* Now generate a single photon in this cell */
      p[n].w = weight;

      pp.freq = 0.0;
      pp.grid = icell;

      /* This following block is a bad way of doing it - kpkt could be modified to 
         do what we want in a more elegant way. */

      test = pp.freq;

      while (test > em_rnge.fmax || test < em_rnge.fmin)
	{
	  kpkt (&pp, &nres, &esc_ptr);
	  if (esc_ptr == 0)
	    {
	      test = 0.0;
	    }
	  else
	    {
	      test = pp.freq;
	    }
	}


      p[n].freq = pp.freq;
      p[n].nres = nres;

      /* The photon frequency is now known. */

      /* Determine the position of the photon in the moving frame */

      /* ! Need to account for possibility that photon is in the torus */
      get_random_location (icell, 0, p[n].x);

      p[n].grid = icell;

      nnscat = 1;
      // Determine the direction of the photon
      // Need to allow for anisotropic emission here
      if (p[n].nres < 0 || p[n].nres > NLINES || geo.scatter_mode == 0)
	{
	  /*  It was either an electron scatter so the  distribution is isotropic, or it
	     was a resonant scatter but we want isotropic scattering anyway.  */
	  randvec (p[n].lmn, 1.0);	/* The photon is emitted isotropically */
	}
      else if (geo.scatter_mode == 1)
	{			// It was a line photon and we want anisotropic scattering

	  // -1. forces a full reinitialization of the pdf for anisotropic scattering

	  randwind (&p[n], p[n].lmn, wmain[icell].lmn);

	}
      else if (geo.scatter_mode == 2)
	{			//It was a line photon and we want the thermal trapping anisotropic model

	  randwind_thermal_trapping(&p[n], &nnscat);
	}

      p[n].nnscat = nnscat;





      /* The next two lines correct the frequency to first order, but do not result in
         forward scattering of the distribution */

      vwind_xyz (&p[n], v);
      p[n].freq /= (1. - dot (v, p[n].lmn) / C);

      p[n].istat = 0;
      p[n].tau = p[n].nscat = p[n].nrscat = 0;
      p[n].origin = PTYPE_WIND;	// Call it a wind photon

    }


  return (nphot);		/* Return the number of photons generated */


  /* All done. */
}

/************************************************************
                                    Imperial College London
Synopsis:

     photo_gen_matom produces photon packets to account for creation of r-packets
     by deactivation of macro atoms in the spectrum calculation. It should only be used 
     once the total energy emitted in this way in the wavelength range in question is well known
     (calculated in the ionization cycles).


Arguments:   
       WindPtr w                   the ptr to the structure defining the wind
       PhotPtr p;                  the ptr to the structire for the photons
       double weight;              the photon weight
       int photstart, nphot;       the number of the first photon to be generated and
                                   the total number of photons to be generated by this
                                   routine

                                           

Returns:   
      When it finishes it should have generated nphot photons from macro atom deactivations.


Description:
      This routine is closely related to photo_gen_kpkt from which much of the code has been copied.

Notes:  
        
         
History:
          June 04 SS - coding began
          Aug  04 SS - modified as for gen_kpkt above by KSL
          Sep  04 SS - minor bug fixed
	06may	ksl	57+ -- Initial adaptation to plasma structure
	0812	ksl	67c -- Changed call to photo_gen_matom to
			eliminate WindPtr w, since we have access
			to this through wmain.

************************************************************/
int
photo_gen_matom (p, weight, photstart, nphot)
     PhotPtr p;
     double weight;
     int photstart, nphot;
{
  int photstop;
  int icell;
  double xlum, xlumsum;
  struct photon pp;
  int nres;
  int n;
  double v[3];
  double dot ();
  int wind_n_to_ij (), stuff_v (), randvec ();
  int emit_matom ();
  double test;
  int upper;
  int nnscat;
  double dvwind_ds (), sobolev ();
  int nplasma;



  photstop = photstart + nphot;
  Log ("photo_gen_matom creates nphot %5d photons from %5d to %5d \n", nphot,
       photstart, photstop);

  for (n = photstart; n < photstop; n++)
    {
      /* locate the wind_cell in which the photon bundle originates. And also decide which of the macro
         atom levels will be sampled (identify that level as "upper"). */

      xlum = (rand () + 0.5) / (MAXRAND) * geo.f_matom;

      xlumsum = 0;
      icell = 0;
      upper = 0;
      while (xlumsum < xlum)
	{
	  if (wmain[icell].vol > 0.0)
	    {
	      nplasma = wmain[icell].nplasma;
	      xlumsum += macromain[nplasma].matom_emiss[upper];
	      upper++;
	      if (upper == nlevels_macro)
		{
		  upper = 0;
		  icell++;
		}
	    }
	  else
	    {
	      icell++;
	      upper = 0;
	    }
	}

      if (upper == 0)
	{
	  /* If upper = 0 at this point it means the loop above finished with an increment of
	     icell and setting upper = 0. In such a case the process we want was the LAST process
	     in the previous cell. */
	  icell--;
	  upper = nlevels_macro;
	}

      upper--;			/* This leaves the macro atom level that deactivaties. */

      /* Now generate a single photon in this cell */
      p[n].w = weight;

      pp.freq = 0.0;
      pp.grid = icell;

      /* This following block is a bad way of doing it but it'll do as
         a quick and dirty test for now. Really kpkt should be modified to 
         do what we want in a more elegant way. */

      test = pp.freq;

      while (test > em_rnge.fmax || test < em_rnge.fmin)
	{

	  /* Call routine that will select an emission process for the
	     deactivating macro atom. If is deactivates outside the frequency
	     range of interest then ignore it and try again. SS June 04. */

	  emit_matom (wmain, &pp, &nres, upper);

	  test = pp.freq;
	}


      p[n].freq = pp.freq;
      p[n].nres = nres;


      /* The photon frequency is now known. */

      /* Determine the position of the photon in the moving frame */

      /* !! ERROR - need to account for posibilyt that photon is in torus */

      get_random_location (icell, 0, p[n].x);

      p[n].grid = icell;


      // Determine the direction of the photon
      // Need to allow for anisotropic emission here
      nnscat = 1;
      if (p[n].nres < 0 || p[n].nres > NLINES || geo.scatter_mode == 0)
	{
	  /*  It was either an electron scatter so the  distribution is isotropic, or it
	     was a resonant scatter but we want isotropic scattering anyway.  */
	  randvec (p[n].lmn, 1.0);	/* The photon is emitted isotropically */
	}
      else if (geo.scatter_mode == 1)
	{			// It was a line photon and we want anisotropic scattering

	  // -1. forces a full reinitialization of the pdf for anisotropic scattering

	  randwind (&p[n], p[n].lmn, wmain[icell].lmn);

	}
      else if (geo.scatter_mode == 2)
	{			//It was a line photon and we want the thermal trapping anisotropic model

      randwind_thermal_trapping(&p[n], &nnscat);
	}
      p[n].nnscat = nnscat;


      /* The next two lines correct the frequency to first order, but do not result in
         forward scattering of the distribution */

      vwind_xyz (&p[n], v);
      p[n].freq /= (1. - dot (v, p[n].lmn) / C);

      p[n].istat = 0;
      p[n].tau = p[n].nscat = p[n].nrscat = 0;
      p[n].origin = PTYPE_WIND;	// Call it a wind photon

    }


  return (nphot);		/* Return the number of photons generated */


  /* All done. */
}

/************************************************************
                                    Imperial College London
Synopsis:
	emit_matom is a scaled down version of matom which deals with the emission due
        to deactivating macro atoms in the detailed spectrum part of the calculation.


Arguments:

       WindPtr w                   the ptr to the structure defining the wind
       PhotPtr p                   the packet at the point of activation
       int upper                   the upper level that we deactivate from

Returns:
       int nres                    the process by which deactivation occurs
       PhotPtr p                   the packet following deactivation
    
Description:


Notes: 

History:
        June 04  -   coding began
	06may	ksl	57+ -- Initial adaptation to plasma structure

************************************************************/

int
emit_matom (w, p, nres, upper)
     WindPtr w;
     PhotPtr p;
     int *nres;
     int upper;
{
  struct lines *line_ptr;
  struct topbase_phot *cont_ptr;
  int uplvl;
  double alpha_sp ();
  double eprbs[NBBJUMPS + NBFJUMPS];
  double penorm;
  double threshold, run_tot;
  double sp_rec_rate;
  int n, m;
  int nbbd, nbfd;
  double t_e, ne;
  double bb_cont;
  WindPtr one;
  PlasmaPtr xplasma;


  one = &w[p->grid];		//This is to identify the grip cell in which we are
  xplasma = &plasmamain[one->nplasma];

  t_e = xplasma->t_e;		//electron temperature 
  ne = xplasma->ne;		//electron number density


  /* The first step is to identify the configuration that has been excited. */

  uplvl = upper;

  /* Unlike the proper matom routine, we don't want any jumps, only deactivations here. */


  /*  The excited configuration is now known. Now compute all the probabilities of deactivation
     /jumping from this configuration. Then choose one. */


  nbbd = config[uplvl].n_bbd_jump;	//store these for easy access -- number of bb downward jumps
  nbfd = config[uplvl].n_bfd_jump;	// number of bf downared jumps from this transition


  // Start by setting everything to 0

  m = 0;			//m counts the total number of possible ways to leave the level
  penorm = 0.0;			//stores the total emission probability

  for (n = 0; n < nbbd + nbfd; n++)
    {
      eprbs[n] = 0;		//stores the individual emission probabilities SS
    }
  /* Finished zeroing. */


  /* bb */

  /* First downward jumps. */

  for (n = 0; n < nbbd; n++)
    {
      line_ptr = &line[config[uplvl].bbd_jump[n]];
      /* Since we are only interested in making an r-packet here we can (a) ignore collisional
         deactivation and (b) ignore lines outside the frequency range of interest. */
      if ((line_ptr->freq > em_rnge.fmin) && (line_ptr->freq < em_rnge.fmax))	// correct range
	{
	  bb_cont = (a21 (line_ptr) * p_escape (line_ptr, xplasma));

	  eprbs[m] = bb_cont * (config[uplvl].ex - config[line[config[uplvl].bbd_jump[n]].nconfigl].ex);	//energy difference

	  if (eprbs[m] < 0.)	//test (can be deleted eventually SS)
	    {
	      Error ("Negative probability (matom, 2). Abort.");
	      exit (0);
	    }

	  penorm += eprbs[m];
	}
      m++;
    }

  /* bf */
  for (n = 0; n < nbfd; n++)
    {
      cont_ptr = &phot_top[config[uplvl].bfd_jump[n]];	//pointer to continuum

      /* If the edge is above the frequency range we are interested in then we need not consider this
         bf process. */
      if (cont_ptr->freq[0] < em_rnge.fmax)	//means that it may contribute
	{
	  sp_rec_rate = alpha_sp (cont_ptr, xplasma, 0);
	  eprbs[m] = sp_rec_rate * ne * (config[uplvl].ex - config[phot_top[config[uplvl].bfd_jump[n]].nlev].ex);	//energy difference
	  if (eprbs[m] < 0.)	//test (can be deleted eventually SS)
	    {
	      Error ("Negative probability (matom, 4). Abort.");
	      exit (0);
	    }
	  penorm += eprbs[m];
	}
      m++;
    }



  /* No need for upwards jumps - they don't provide emission mechanisms. */

  /* Probabilities of emission (e) are now known. The normalisation factor pnorm has
     also been recorded. the integer m now gives the total number of possibilities too. 
     now select what happens next. Start by choosing the random threshold value at which the
     event will occur. */

  threshold = ((rand () + 0.5) / MAXRAND);

  run_tot = 0;
  n = 0;
  threshold = threshold * penorm;	//normalise to total emission prob.
  while (run_tot < threshold)
    {
      run_tot += eprbs[n];
      n++;
    }
  n = n - 1;
  /* n now identifies the jump that occurs - now set nres for the return value. */
  if (n < nbbd)
    {				/* bb downwards */
      line_ptr = &line[config[uplvl].bbd_jump[n]];	//pointer for the bb transition
      *nres = line[config[uplvl].bbd_jump[n]].where_in_list;
      p->freq = line[config[uplvl].bbd_jump[n]].freq;
      /* This is the comoving frequency - changed to rest frequency by doppler */
    }
  else if (n < (nbbd + nbfd))
    {				/* bf downwards jump */
      *nres = config[uplvl].bfd_jump[n - nbbd] + NLINES + 1;
      /* continuua are indicated by nres > NLINES */
      p->freq = phot_top[config[uplvl].bfd_jump[n - nbbd]].freq[0]
	- (log (1. - (rand () + 0.5) / MAXRAND) * t_e / H_OVER_K);
      /* Co-moving frequency - changed to rest frequency by doppler */
      /*Currently this assumed hydrogenic shape cross-section - Improve */
    }
  else
    {
      Error
	("Trying to emit from Macro Atom but no available route (emit_matom). Abort.");
      exit (0);
    }
  return (0);
}

/* The frequency and the value of nres have been set correctly. All done. */

/********************************************************************************/


/******************************************************************************/

/* q_ioniz. This returns the collisional ionization co-efficient
Calculated following equation 5-79 of Mihalas.
*/

double
q_ioniz (cont_ptr, electron_temperature)
     struct topbase_phot *cont_ptr;
     double electron_temperature;
{
  double coeff;
  double gaunt;
  double u0;

  u0 = cont_ptr->freq[0] * H_OVER_K / electron_temperature;
  gaunt = 0.1;			//for now - from Mihalas for hydrogen

  coeff = 1.55e13 / sqrt (electron_temperature) * gaunt * cont_ptr->x[0] *
    exp (-1. * u0) / u0;

  return (coeff);
}

/******************************************************************************/

/* q_recomb. This returns the collisional recombination co-efficient
Calculated from inverse of q_ioniz.

JM 1301 -- Edited this to avoid need to call q_ioniz and exponential.
Should improve speed and stability

This equation comes from considering TE and getting the expression
q_recomb = 2.07e-16 * gl/gu * exp(E/kT) * q_ioniz
then substituting the above expression for q_ioniz.
*/

double
q_recomb (cont_ptr, electron_temperature)
     struct topbase_phot *cont_ptr;
     double electron_temperature;
{
  double coeff;
  double gaunt;


  gaunt = 0.1;			//for now - from Mihalas for hydrogen

  coeff = 3.2085e-3  / electron_temperature * gaunt * cont_ptr->x[0];	// normal constants * 1/T times gaunt * cross section

  coeff /= cont_ptr->freq[0] * H_OVER_K;			// divide by h nu / k
  coeff *= config[cont_ptr->nlev].g / config[cont_ptr->uplev].g;

  return (coeff);
}


