

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
	     
          /* this error condition can happen in unconverged hot cells where T_R >> T_E.
             for the moment we set to 0 and hope spontaneous recombiantion takes care of things */
	      /* note that we check and report this in check_stimulated_recomb() in estimators.c once a cycle */
	      if (jprbs[m] < 0.)	//test (can be deleted eventually SS)
		{
		  //Error ("Negative probability (matom, 6). Abort?\n");
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

  x = sigma_phot (cont_ext_ptr, freq);	//this is the cross-section
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
  double cooling_bf[nphot_total];
  double cooling_bf_col[nphot_total];	//collisional cooling in bf transitions
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
  double freqmin, freqmax;


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

  /* JM 1511 -- Fix for issue 187. We need band limits for free free packet
     generation (see call to one_ff below) */
  if (geo.ioniz_or_extract)
    {
      /* in spectral cycles, so use the boundaries of the photon generation bands */
      freqmin = xband.f1[0];
      freqmax = xband.f2[xband.nbands-1];
    }
  else
    {
      /* in spectral cycles, use the frequency range of the final spectrum */
      freqmin = em_rnge.fmin;
      freqmax = em_rnge.fmax;
    }


  /* ksl 091108 - If the kpkt destruction rates for this cell are not known they are calculated here.  This happens
   * every time the wind is updated */

  if (mplasma->kpkt_rates_known != 1)
    {
      cooling_normalisation = 0.0;
      cooling_bftot = 0.0;
      cooling_bbtot = 0.0;
      cooling_ff = 0.0;
      cooling_bf_coltot = 0.0;

      /* JM 1503 -- we used to loop over ntop_phot here, 
         but we should really loop over the tabulated Verner Xsections too
         see #86, #141 */
      for (i = 0; i < nphot_total; i++)
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
	      Error ("i, ulvl, nphot_total, nion %d %d %d %d\n", i, ulvl,
		     nphot_total, cont_ptr->nion);
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

      /* end of loop over nphot_total */

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
			VERY_BIG) / xplasma->vol / xplasma->ne;		// JM 1411 - changed to use filled volume
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

      cooling_adiabatic = xplasma->lum_adiabatic / xplasma->vol / xplasma->ne; // JM 1411 - changed to use filled volume

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

      /* JM 1503 -- we used to loop over ntop_phot here, 
         but we should really loop over the tabulated Verner Xsections too
         see #86, #141 */
      for (i = 0; i < nphot_total; i++)
	{
	  if (destruction_choice < mplasma->cooling_bf[i])
	    {
	      /* Having got here we know that desturction of the k-packet was via the process labelled
	         by i. Let's just check that i is a sensible number. */

	      if (i > nphot_total - 1)
		{
		  Error
		    ("kpkt (matom.c): trying to destroy k-packet in unknown process. Abort.\n");
		  exit (0);
		}

	      /* If it gets here, all seems fine. Now set nres for the destruction process. */

	      *nres = i + NLINES + 1;
	      *escape = 1;	//record that an r-packet is made - no need to excite a macro atom again


	      /* Now (as in matom) choose a frequency for the new packet. */

	      p->freq = phot_top[i].freq[0] -
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

      /* used to do one_ff (one, 7.5e12, 2.626e16) here,
         but now use the band boundaries, see #187.*/
      p->freq = one_ff (one, freqmin, freqmax);	//get frequency of resulting energy packet

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

      for (i = 0; i < nphot_total; i++)
	{
	  if (destruction_choice < mplasma->cooling_bf_col[i])
	    {
	      /* Having got here we know that destruction of the k-packet was via the process labelled
	         by i. Let's just check that i is a sensible number. */

	      if (i > nphot_total - 1)
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


