/* macro_gov.c is a file which contains the functions which govern macro-atoms and obtain
   their level populations. The actual functions which do the jumps inside an activated 
   macro-atom are in matom.c. This is partly done to prevent overly long files (JM1504)
*/


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
  double *a_data, *b_data;
  double *populations;
  int index_fast_col, ierr, insane, sane_populations;

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

	  sane_populations = 0;
      while (sane_populations == 0)
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
		      rate = b12 (line_ptr) *
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

		      if (rate < 0.0 || sane_check(rate))
              {
		      	Error("macro_pops: bbu rate is %8.4e in cell/matom %i\n", rate, xplasma->nplasma);
		      }

		      /* There's a radiative jump between these levels, so we want to clean
		         for popualtion inversions. Flag this jump */
		      radiative_flag[index_lvl][line_ptr->nconfigu] = 1;
		    }

		  for (index_bbd = 0;
		       index_bbd < config[index_lvl].n_bbd_jump; index_bbd++)
		    {
		      /* These are bb downwards jumps. The rate in such a jump is given by
		         the A-value. I'm also
		         including a collisional term (which depends on ne). */

		      line_ptr = &line[config[index_lvl].bbd_jump[index_bbd]];
		      rate = (a21 (line_ptr) * p_escape (line_ptr, xplasma));
		      //rate =0.0;
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
	  populations = (double *) calloc (sizeof (rate), n_macro_lvl);

	  /* replace the first entry with 1.0- this is part of the normalisation constraint */
	  b_data[0] = 1.0;

      /* this next routine is a general routine which solves the matrix equation
         via LU decomposition */
	  ierr = solve_matrix(a_data, b_data, n_macro_lvl, populations, xplasma->nplasma);

	  if (ierr != 0)
	  	Error("macro_pops: bad return from solve_matrix\n");

      /* free memory */
      free (a_data);	
	  free (b_data);


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
		        inversion_test = populations[conf_to_matrix[index_lvl]] * config[nn].g / config[index_lvl].g * 0.999999;	//include a correction factor 

		        if (populations[conf_to_matrix[nn]] > inversion_test)
			  {
			    populations[conf_to_matrix[nn]] = inversion_test;
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
	  insane = 0;


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
		  level_population = populations[conf_to_matrix[index_lvl]];	
		  this_ion_density += level_population;

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
		    populations[conf_to_matrix[index_lvl]] /
		    this_ion_density;

		  if (xplasma->levden[config[index_lvl].nden] < 0.0 || sane_check(xplasma->levden[config[index_lvl].nden]))
		  {
		  	Error("macro_pops: level %i has frac. pop. %8.4e in cell %i\n", 
		    	   index_lvl, xplasma->levden[config[index_lvl].nden], xplasma->nplasma);
		  	insane = 1;
		  }

          
		  mm++;

		}

	      xplasma->density[index_ion] =
		this_ion_density * ele[index_element].abun * xplasma->rho *
		rho2nh;

		  if (sane_check(xplasma->density[index_ion]) ||  xplasma->density[index_ion] < 0.0 )
		  {
		  	Error("macro_pops: ion %i has frac. pop. %8.4e in cell %i\n", 
		    	   index_ion, xplasma->density[index_ion], xplasma->nplasma);
		  	insane = 1;
		  }




	    }

      /* if the variable insane has been set to 1 then that means we had either a negative or 
         non-finite level population somewhere. If that is the case, then set all the estimators
         to dilute blackbodies instead and go through the solution again */
      if (insane)
        {
          Error("macro_pops: found unreasonable populations in cell %i, so adopting dilute BBody excitation\n", xplasma->nplasma);
          get_dilute_estimators(xplasma);
        }
      /* if we didn't set insane to 1 then we have a realistic set of populations, so set sane_populations to 1 to break
         the while loop */
      else
      	sane_populations = 1;	

	  }  // end of while sane loop


	}  // end of if statement for macro-atoms
    }  // end of elements loop


  return (0);
  /* All done. (SS, Apr 04) */
}
