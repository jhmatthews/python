

/* This file contains various ways in which to calculate the
   ionization of a plasma.  The guiding difference between
   whether to put routines here or in ionization.c appears to
   be that here one takes the conditions, radiative weight as
   given and calculate the densities of each ionic species, while
   ionization.c is involved in trying to match heating and cooling.

	00Jan02	ksl	Removed some of the routines which allowed
   			special calculation of H and He 
	01dec	ksl	Major modifications to all routines.  (a) update
			calls to concentrations and nebular_concentrations
			so WindPtrs are passed.
	02jun	ksl	Modified calls yet again.  This time the goal
			is to have all the modes etc interpreted in
			nebular concentrations, so one can actually tell what 
			is happening in nebular concentrations, and to
			avoid doing any real work in nebular concentrations.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>


#include "atomic.h"
#include "python.h"

/***********************************************************
                                       Space Telescope Science Institute

  Synopsis:   

int
nebular_concentrations (xplasama, mode)  modifies the densities of ions, levels, and
	partition functions of ions within a cell of the wind based upon the mode, 
	and other data contained within the WindPtr itself, such as t_r, t_e, w, 
	based on the "mode".  

        This version of the code is a first, quick and dirty, attempt to include the
        fraction of recombinations going into the ground state as a correction factor
        following Mazzali and Lucy (CK).  In other words, ionization equilibrium here
        is calculated from equation 11 of Mazzali & Lucy.
  
  Arguments:		
     PlasmaPtr ww;
     int mode;			// 0=saha using tr, 1=saha using te, 2= Lucy & Mazzali


  Returns:
 	nebular concentrations alters portions of the wind ptr.  Exactly how things 
	are changed depends on the mode.
 
 	nebular_concentrations returns 0 if it converged to an answer, -1 otherwise.  On an
 	abnormal return the density array and ne are left unchanged.
 	
  Description:	

	As of 01dec (python39), nebular concentrations serves as the steering routine
	for all ionization calculations.  This change was made to push down all of the
	various choices about calculation of ionization to a single routine, so that
	no longer needs to include case statements or logic for how to handle the various
	choices elsewhre.  The Lucy & Mazzali calculation is still carried out within
	nebular concentrations although that is historical at present.

	It remains true that for the Lucy and Mazalli calculation, concentrations is
	called initially, a choice which today seems bizarre...since the LM calculation
	can be carried equally well straight away.
	
	My pre-01dec notes said:

	This is a first attempt to fudge the ionization equation in the wind to account 
	for the fact that one is not sitting in a BB radiation field.  The routine assumes one
	has already put the concentrations of of the ions under LTE conditions for a temperature
	t_r into density.  It uses the same dumb procedure as for "concentrations."  

 
  Notes:
	It's probably OK but it would be worthwhile worrying whether I have introduced 
	any problems by setting the floor to any species and ne of DENSITY_MIN	

	Mode is not consistent across all of the parts of python, so don't assume that
	mode here is geo.ionization_mode...for example.


  History:
	1997	ksl	Coded and debugged as part of Python effort. 
	97aug27	ksl	Eliminated the goto statement in the main iteration loop and
 			added numerous comments.
	97oct2	ck	Implemented Mazzali & Lucy eqn 11.
	98may	ksl	Moved ne calculation to separate subroutine
	98may	ksl	Moved actual calculation for Mazzali & Lucy to a separate routine for
 			each element
	01dec	ksl	Changed the way in which nebular concentrations was called, so that all
			the variables are taken from the windptr except for the mode. In addition
			nebular_concentrations was made the sole call from other parts of python, and
			so for example saha should never be called elsewhere.  This last step wasw
			to encapsulate the determination of abundances (given known conditions.)
	02jun	ksl	Moved Lucy-Mazzali calculation to a separate routine.  Now nebular
			concentrations is simply a driving routine
	06may	ksl	57+ -- Modified to use plasma structure.  There are no volumes involved
			so use plasma stucture in call to routines.
	080808	ksl	60b -- Removed mode 4 which attempted to carry out detailed balance
			for H and He and ussd LM for other elements.  (This was carried out
			in the routine dlucy) This is wholly replaced by the macro atom approach.  
	11	ksl	Main changes here made by Nick to incorporate the power law approximation
			to ionization for AGN
        12Feb   nsh	More changes made to allow for two new modes - for multiple saha equaions.
			mode 6 corrects using a dilute blackbody and should be almost the same as the 				LM method
			mode 7 corrects using a power law spectrum, and show end up supplanting mode 5.


**************************************************************/

int
nebular_concentrations (xplasma, mode)
     PlasmaPtr xplasma;
     int mode;			// 0=saha using tr, 1=saha using te, 2= Lucy & Mazzali
{
  double get_ne ();
  int lucy_mazzali1 ();
  int m;

  if (mode == NEBULARMODE_TR)
    {				// LTE all the way -- uses tr

      partition_functions (xplasma, mode);

      m = concentrations (xplasma, mode);

    }
  else if (mode == NEBULARMODE_TE)
    {				//LTE but uses temperature te

      partition_functions (xplasma, mode);

      m = concentrations (xplasma, mode);

    }
  else if (mode == NEBULARMODE_ML93)		// This is the standard LM method
    {


      partition_functions (xplasma, mode);	// calculate partition functions, using t_r with weights

      /* JM 1308 -- concentrations then populates xplasma with saha abundances. 
	 in macro atom mode it also call macro_pops, which is done incorrectly, 
	 the escape probabilities are calculated with saha ion densities each time,
	 because saha() repopulates xplasma in each cycle. */

      m = concentrations (xplasma, NEBULARMODE_TR);	// Saha equation using t_r

      /* JM 1308 -- lucy then applies the lucy mazzali correction factors to the saha abundances. 
	 in macro atom mode it also call macro_pops which is done correctly in this case, as lucy_mazzali, 
	 which is called in lucy, does repopulate xplasma like saha does. This means that in this mode,
	 it doesn't actually matter that concentrations does macro level populations wrong, as that is 
         corrected here. It should be sorted veyr soon, however. */

      m = lucy (xplasma);	// Main routine for running LucyMazzali


    }
 // else if (mode == 5)		NSH 130626 - mode removed		
				/* This replicates Sim's (2008) power
				   law method for ionization in a non-BB radiation
				   field.  */
 //   {

 //     partition_functions (xplasma, 1);	//lte partition function using t_e and no weights

 //     m = concentrations (xplasma, 1);	// Saha equation using t_e 

 //     m = sim_driver (xplasma);
 //   }

  /* Two new modes, they could probably be combined into one if statement, but having 
     two adds little complexity and allows for other   modifications if required. No 
     call to partition functions is required, since this is done on a pairwise basis 
     in the routine. Similarly there is no call to concentrations, since this is also 
     done pairwise inside the routine. */

  else if (mode == NEBULARMODE_PAIRWISE_ML93)		/* Pairwise calculation of abundances, using a 
				   temperature computed to ensure a reasonable
				   ratio between the two, and then corrected for
				   a dilute blackbody radiation field. */
    {
      m = variable_temperature (xplasma, mode);
    }

  else if (mode == NEBULARMODE_PAIRWISE_SPECTRALMODEL)		/* Pairwise calculation of abundances, using a 
				   temperature computed to ensure a reasonable
				   ratio between the two, and then corrected for
				   a radiation field modelled by a power law */
    {
      m = variable_temperature (xplasma, mode);
    }

  else
    {
// If reached this point the program does not understand what is desired.
      Error ("nebular_concentrations: Unknown mode %d\n", mode);
      exit (0);

    }


  return (m);
}


/***********************************************************
                                       Space Telescope Science Institute

 Synopsis:   
int
concentrations (xplasma, mode)
     PlasmaPtr xplasma;
     int mode;			//   0=saha using tr, 1=saha using te
  	determines the concentrations of individual ions assuming a form of the 
	Saha equantion.  The atomic data required must have been read in with the 
	routine "get_atomic_data".  Most of the actual information to calculate
	the concentrations is already contained in ww
  
 Arguments:		

     PlasmaPtr xplasma;
     int mode;			//   0=saha using tr, 1=saha using te

 Returns:
	concentrations returns 0 if it converged; -1 if it did not converge.  If it does
	not converge, both density and *ne are likely to be garbage.

	Other results...e.g. the densities.. are stored in ww.
 
 Description:	
 	
	The Saha eqn. is  (e.g. Motz, p 110, eqn 4.14

  	 N_r+1 N_e / N_r =  (2*pi*m*kT)**1.5/h**3 2 Z_r+1(T)/Z_r * exp (-X/kT)
   
	where N_r+1 is the number density of the r+1 ionization state
	and   Z_r+1 is the partition function of that state.  

	The program works by perfoming a crude iteration on ne.  It could undoubtedly 
	be improved, if one wanted better accuracy or higher speed.


 
 Notes:
	
	So that one doesn't get zero-divides in other programs, e.g. nebular_concentrations, 
	a floor is set on the density of any ion or the electron density.



 History:
 	1997      ksl	Coded and debugged as part of Python effort.  All of the partition functions
 				have been set to 1.
 	97aug27 ksl	Eliminated the goto statement in the iteration of the Saha equation and
 				added numberous comments.
 	97aug27 ksl      Updated code to include an estimate of the partition function based on
 				the ground state multiplities.  Add an error return if the function failed
 				to return.
 	98may	ksl	Moved ne calculation to separate subroutine
 	98may	ksl	Moved saha calculation for an individual element to a separate routine
	00mar	ksl	Fixed error in initial estimate of ne so that it is correct at high temperatures
			Basically the problem was due to the fact that at high enough temperatures
			a subtraction was being performed which could lead to an initial value of
			0 for ne from which the code never recovered.  Note--I believe there still
			may be a small problem with this section of the code because it does not
			take account of the partition function in the initial estimate.  On the
			other hand the effect of this error should be minor.
	01jul	ksl	Added a real partition function to routine
	01dec	ksl	Changed variables used to call saha, and made a subroutine of nebular
			concentrations (python_39)
	02jul	ksl	Included the electron multiplicity in the SAHA constant.
      
        04Apr   SS      Added a call to macro_pops so that macro atom ionization fractions are used
                        in place of others.
	06may	ksl	57+ -- Switched to plasma structue wind no volumes
	07mar	ksl	Convert warnings to errors to stop many writes 
**************************************************************/



//#define SAHA 4.82907e15               /* 2* (2.*PI*MELEC*k)**1.5 / h**3  (Calculated in constants) */
//#define MAXITERATIONS 200
//#define FRACTIONAL_ERROR 0.03
//#define THETAMAX       1e4
//#define MIN_TEMP         100. NSH 0712 - moved into python.h because this is used in severalpaces


int
concentrations (xplasma, mode)
     PlasmaPtr xplasma;
     int mode;			//   0=saha using tr, 1=saha using te
{
  int nion, niterate;
  double xne, xxne, xnew, xsaha;
  double theta, x;
  double get_ne ();
  double t, nh;
  int saha ();


  // This needs to be moved up into nebular_concentrations given that we
  // do not want these routines to do mode shifting.
  //
  if (mode == NEBULARMODE_TR)
    {
      t = xplasma->t_r;
    }
  else if (mode == NEBULARMODE_TE)
    {
      t = xplasma->t_e;
    }
  else if (mode == NEBULARMODE_ML93)
    {
      t = sqrt (xplasma->t_e * xplasma->t_r);
    }
  else
    {
      Error ("Concentrations: Unknown mode %d\n", mode);
      mytrap ();
      exit (0);
    }

  nh = xplasma->rho * rho2nh;	//LTE

  /* make an initial estimate of ne based on H alone,  Our guess
     assumes ion[0] is H1.  Note that x below is the fractional
     ionization of H and should vary from 0 to 1

     Note that in a pure H plasma, the left hand side of the Saha 
     equation can be written as (x*x*nh*nh)/((1-x)*nh)  and can
     then be converted into a quadratic.  That is what is done
     below.  

     Since this is an initial estimate g factors have been ignored
     in what follows.
   */

  if (t < MIN_TEMP)
    t = MIN_TEMP;		/* fudge to prevent divide by zeros */

  xsaha = SAHA * pow (t, 1.5);

  theta = xsaha * exp (-ion[0].ip / (BOLTZMANN * t)) / nh;

  if (theta < THETAMAX)
    {
      x = (-theta + sqrt (theta * theta + 4 * theta)) / 2.;
      xne = xxne = x * nh;
    }
  else
    xne = xxne = nh;

  if (xne < 1.e-6)
    xne = 1.e-6;		/* fudge to assure we can actually calculate
				   xne the first time through the loop */

  /* At this point we have an initial estimate of ne. */


  niterate = 0;
  while (niterate < MAXITERATIONS)
    {

      /* Assuming a value of ne calculate the relative densities of each ion for each element */
      /* JM 1308 -- Here we  actually populate the plasma structure with saha abundances. In Lucy Mazzali mode, 
	 we apply a correction to these abundances. Note that in macro atom mode this call should 
	 really only happen before the first ionization cycle */	

      saha (xplasma, xne, t);


      /* New (SS Apr 04) call the routine that will replace the macro atoms ionization fractions and
         level populations with those computed with the Monte Carlo estimators. geo.macro_iniz_mode
         controls the use of macro_pops. */
      /* JM 1308 -- At the moment, the escape probabilities in macro_pops are calculated using saha ion 
	 densities which is wrong. Fortunately, macro_pops is called again in lucy() and converges on the 
         correct value, but this should be fixed. I am not sure if it even needs to be called at all. */

      if (geo.macro_ioniz_mode == 1)
	{
	  macro_pops (xplasma, xne);
	}

      /*Set some floor so future divisions are sensible */
      for (nion = 0; nion < nions; nion++)
	{
	  if (xplasma->density[nion] < DENSITY_MIN)
	    xplasma->density[nion] = DENSITY_MIN;
	}

      /* Now determine the new value of ne from the ion abundances */
      xnew = get_ne (xplasma->density);	/* determine the electron density for this density distribution */

      if (xnew < DENSITY_MIN)
	xnew = DENSITY_MIN;	/* fudge to keep a floor on ne */

      if (fabs ((xne - xnew) / (xnew)) < FRACTIONAL_ERROR || xnew < 1.e-6)
	break;

      xne = (xnew + xne) / 2.;	/* Make a new estimate of xne */
      niterate++;
    }

  if (niterate == MAXITERATIONS)
    {
      Error ("concentrations: failed to converge t %.2g nh %.2g xnew %.2g\n",
	     t, nh, xnew);
      Error ("concentrations: xxne %e theta %e\n", xxne, theta);
      return (-1);
    }

  xplasma->ne = xnew;
  return (0);
}



/***********************************************************
                                       Space Telescope Science Institute

  Synopsis:   

   	Calculate the saha densities for all of the ions in a single
   	cell (given ne). 
  
  Arguments:
	xplasma, ne, t
	the plasma cell to populate, the electron density and the temperature	


  Returns:
 	
  Description:	

 
  Notes:

  	080808 - In the new version of this routine, all of the partition
  	functions are assumed to have been calculated before entering
  	the routine.
	
	Note that this routine populates the actual xplasma structure
   	with saha abundances.

	This routine assumes ne is known.  Iteration to calculate ne
	takes place in concentrations.


  History:
   98may        ksl     Abstracted from concentrations
   01jul	ksl	Added partition functions more properly
   080808	ksl	Modified so that calls parallel those of other 
   			functions, that is so that this uses
			the plasma array


**************************************************************/

int
saha (xplasma, ne, t)
     PlasmaPtr xplasma;
     double ne, t;

{
  double nh;
  int nelem;
  double *density, *partition;

  int first, last, nion;
  double xsaha;
  double sum, a, b;
  double big;

  density = xplasma->density;
  partition = xplasma->partition;

  nh = xplasma->rho * rho2nh;	//LTE
  xsaha = SAHA * pow (t, 1.5);

  for (nelem = 0; nelem < nelements; nelem++)
    {
      first = ele[nelem].firstion;	/*first and last identify the position in the array */

      last = first + ele[nelem].nions;	/*  So for H which has 2 ions, H1 and H2, first will generally
					   be 0 and last will be 2 so the for loop below will just be done once for nion = 1 */

      if (first < 0 || first >= nions || last < 0 || last > nions)
	{
	  Error
	    ("saha: Confusion for element %d with first ion %d and last ion %d\n",
	     nelem, first, last);
	  exit (0);
	}

/*    These lines were put in to make sim work properly, ideally there should be a switch so 
      if we are doing things the old way, we keep the old numbers. But, saha doesnt know the mode....
      sum = density[first] = 1e-250;
      big = pow (10., 250. / (last - first));
      big=big*1e6;   */

      sum = density[first] = 1.0;
      big = pow (10., 250. / (last - first));

      for (nion = first + 1; nion < last; nion++)
	{ 

	  /* JM 1309 -- this next if statement is to ensure that saha densities are only calculated if the 
             ion in question is not a macro ion. Otherwise, this will affect the escape 
             probabilities that are calculated in macro_pops. The exception to this is prior
             to the first ionization cycle when we need to populate saha densities as a first guess */

          if ((ion[nion].macro_info == 0) || (geo.macro_ioniz_mode == 0))
            {
	      b = xsaha * partition[nion]
	          * exp (-ion[nion - 1].ip / (BOLTZMANN * t)) / (ne *
							   partition[nion -
								     1]);
//        if (b > big && nh < 1e5) this is a line to only modify things if the density is high enough to matter
	      if (b > big)
	        b = big;		//limit step so there is no chance of overflow

	      a = density[nion - 1] * b;

	      sum += density[nion] = a;
  
	      if (density[nion] < 0.0)
	        mytrap ();

	      if (sane_check (sum))
	        Error ("saha:sane_check failed for density summation\n");
  
            }
	}

      a = nh * ele[nelem].abun / sum;

      for (nion = first; nion < last; nion++)
	{

          /* JM 1309 -- this next if statement is to ensure that saha densities are only calculated if the 
             ion in question is not a macro ion. Otherwise, this will affect the escape 
             probabilities that are calculated in macro_pops. The exception to this is prior
             to the first ionization cycle when we need to populate saha densities as a first guess */

          if ((ion[nion].macro_info == 0) || (geo.macro_ioniz_mode == 0))
            {

	      density[nion] *= a;
	      if (sane_check (density[nion]))
	        Error ("saha:sane_check failed for density of ion %i\n", nion);
            }    
	}
    }


  return (0);
}



/***********************************************************
                                       Space Telescope Science Institute

  Synopsis:   
	lucy() applied the Lucy Mazzali correction factors to the saha abundances
	contained in xplasma. It thus assumes that saha() has been called before it, 
	in the concentrations routine.

	Note that lucy also applied the macro_pops routine, and required that the 
	ne convergence criteria is rigorous enough that the macro atom level populations
	converge on the correct values.

  Arguments:
	xplasma, the plasma pointer for the cell in question		

  Returns:
	iterates on ne to return the corrected lucy mazzali abundances and an improved
	ne value.
 	
  Description:	

 
  Notes:

	Concentrations should have been called before this routine is executed.

	Procedurally, the routine is analogous to concentrations()

  History:

	02jun	ksl	Made separate routine, removing it from nebular concentrations
        04Apr   SS      If statement added to stop this routine from changing macro atom
	                populations.
        04May   SS      "If" statement modified for compatibility with the "macro_simple" option
	                (i.e. all ions treated as simple).
	07mar	ksl	Convert warnings to errors to stop many repeads



**************************************************************/


#define MIN_FUDGE  1.e-10
#define MAX_FUDGE  10.

int
lucy (xplasma)
     PlasmaPtr xplasma;
{
  int nelem, nion, niterate;
  double xne, xnew;
  double newden[NIONS];
  double t_r, nh;
  double t_e, www;

  t_r = xplasma->t_r;
  t_e = xplasma->t_e;
  www = xplasma->w;


  /* Initally assume electron density from the LTE densities */

  xne = xplasma->ne;
  if (xne < DENSITY_MIN)
    {
      Error
	("nebular_concentrations: Very low ionization: ne initially %8.2e\n",
	 xne);
      xne = DENSITY_MIN;
    }

  nh = xplasma->rho * rho2nh;	//LTE -- Not clear needed at this level

  /* Begin iteration loop to find ne */
  niterate = 0;
  while (niterate < MAXITERATIONS)
    {
      for (nelem = 0; nelem < nelements; nelem++)
	{

	  /* JM1308 -- Here we apply the lucy/mazzali correction factors to the saha abundances
	     which are contained in xplasma. These corrected abundances are copied to newden, which
	     is transferred over to xplasma when we converge on ne */

	  lucy_mazzali1 (nh, t_r, t_e, www, nelem, xplasma->ne,
			 xplasma->density, xne, newden);
  }

    /* Re solve for the macro atom populations with the current guess for ne */
    /* JM1308 -- note that unlike lucy mazzali above, here we actually modify the xplasma
       structure for those ions which are being treated as macro ions. This means that the
       the newden array will contain wrong values for these particular macro ions, but due
       to the if loop at the end of this subroutine they are never passed to xplasma */
	  if (geo.macro_ioniz_mode == 1)
	    {
	      macro_pops (xplasma, xne);
	    }


      for (nion = 0; nion < nions; nion++)
	{

	  /* if the ion is being treated by macro_pops then use the populations just computed */
	  if ((ion[nion].macro_info == 1) && (geo.macro_simple == 0)
	      && (geo.macro_ioniz_mode == 1))
	    {
	      newden[nion] = xplasma->density[nion];
	    }

	  /*Set some floor so future divisions are sensible */
	  if (newden[nion] < DENSITY_MIN)
	    newden[nion] = DENSITY_MIN;
	}
      xnew = get_ne (newden);
      if (xnew < DENSITY_MIN)
	xnew = DENSITY_MIN;

      /* Check to see whether the search for xne has converged and if so exit loop */
      if (fabs ((xne - xnew) / (xnew)) < FRACTIONAL_ERROR || xnew < 1.e-6)
	break;

      /* else start another iteration of the main loop */
      xne = (xnew + xne) / 2.;	/* Make a new estimate of xne */
      niterate++;
    }
  /* End of main iteration loop */

  if (niterate == MAXITERATIONS)
    {
      Error
	("nebular_concentrations: failed to converge:nh %8.2e www %8.2e t_e %8.2e  t_r %8.2e \n",
	 nh, www, t_e, t_r);
      return (-1);
    }

/* Finally transfer the calculated densities to the real density array */

  xplasma->ne = xnew;
  for (nion = 0; nion < nions; nion++)
    {
      /* If statement added here to suppress interference with macro populations (SS Apr 04) */
      if (ion[nion].macro_info == 0 || geo.macro_ioniz_mode == 0
	  || geo.macro_simple == 1)
	{
	  xplasma->density[nion] = newden[nion];
	}
    }
  return (0);
}





/***********************************************************
                                       Space Telescope Science Institute

  Synopsis:   
	lucy_mazzali1(nh,t_r,t_e,www,nelem,ne,density,xne,newden) calculates ion
	densities of a single element (defined by nelem) according to a modified on
	the spot aprroximation, in this case equation 11 of Mazzali and Lucy which
	attempts to consider ionzations from excieted and ground state levels.
 
  Arguments:		

	double nh,ne		number density of H atoms and the electron density
	double t_r,t_e		temperature of the radiation field and the electrons in degrees
	double density[]	densities of individual ions under pure Saha conditions
	double www		dilution factor for the radiation field
	int nelem;		the element number for which the calculation is carried out

  Returns:
 	double newden[]		the modified concentrations of the individual ions. Array 
 				elements of density correspond to those ions in the structure ions.
 
 	nebular_concentrations returns 0 if it converged to an answer, -1 otherwise.  On an
 	abnormal return the density array and ne are left unchanged.
 	
 Description:	

	This is a first attempt to fudge the ionization equation in the wind to account 
	for the fact that one is not sitting in a BB radiation field.  The routine assumes one
	has already put the concentrations of of the ions under LTE conditions for a temperature
	t_r into density.  It uses the same dumb procedure as for "concentrations."  

 
  Notes:
	It's probably OK but it would be worthwhile worrying whether I have introduced 
	any problems by setting the floor to any species and ne of DENSITY_MIN	


  History:
 
	98may	ksl	Coded as a a separate routine for each element
	07mar	ksl	Convert warnings to errors to stop many repeads
	12july 	nsh	Modified to use external code to compute zeta.

**************************************************************/


int
lucy_mazzali1 (nh, t_r, t_e, www, nelem, ne, density, xne, newden)
     double nh, t_r, t_e, www;
     int nelem;
     double ne, density[], xne, newden[];
{
  double fudge;			// dummy, interpfrac; 0712 zeta calculation moved into zeta.c
  double fudge2, q;
  double sum, a;
//  int ilow, ihi; 0712 nsh zeta moved into subroutine
  int first, last, nion;
  double numerator, denominator;
  if (t_r > MIN_TEMP)
    {
      fudge = www * sqrt (t_e / t_r);

//NSH 0712These lines now moved into zeta.c we retain the caculation of sqrt t_e/t_r
      /* now get the right place in the ground_frac tables  CK */
//      dummy = t_e / TMIN - 1.;
//      ilow = dummy;           /* have now truncated to integer below */
//      ihi = ilow + 1;         /*these are the indeces bracketing the true value */
//      interpfrac = (dummy - ilow);    /*this is the interpolation fraction */
//      if (ilow < 0)
//      {
//        ilow = 0;
//        ihi = 0;
//        interpfrac = 1.;
//      }
//    if (ihi > 19)
//      {
//        ilow = 19;
//        ihi = 19;
//        interpfrac = 1.;
//      }

    }
  else
    {
      Error_silent
	("lucy_mazzali1: t_r too low www %8.2e t_e %8.2e  t_r %8.2e \n", www,
	 t_e, t_r);
      fudge = 0.0;
//      interpfrac = 0.0;
//      ihi = ilow = 0;
    }

  if (fudge < MIN_FUDGE || MAX_FUDGE < fudge)
    {
      Error_silent
	("lucy_mazzali1: fudge>10 www %8.2e t_e %8.2e  t_r %8.2e \n",
	 www, t_e, t_r);
    }

  /* Initialization of fudges complete */

  first = ele[nelem].firstion;	/*identify the position of the first and last ion in the array */
  last = first + ele[nelem].nions;	/*  So for H which has 2 ions, H1 and H2, first will generally
					   be 0 and last will be 2 so the for loop below will just be done once for nion = 1 */


  while (density[first] < 1.1 * DENSITY_MIN)
    {
      newden[first] = DENSITY_MIN;
      first++;
    }
  while (density[last - 1] < 1.1 * DENSITY_MIN)
    {
      newden[last - 1] = DENSITY_MIN;
      last--;
    }

  sum = newden[first] = 1.;
  for (nion = first + 1; nion < last; nion++)
    {
      numerator = newden[nion - 1] * fudge * (ne) * density[nion];
      denominator = density[nion - 1] * xne;
      q = numerator / denominator;

      /* now apply the Mazzali and  Lucy fudge, i.e. zeta + w (1-zeta)
         first find fraction of recombinations going directly to
         the ground state for this temperature and ion */

      fudge2 = compute_zeta (t_e, nion - 1, 2);	/* NSH 1207 - call external function, mode 2 uses chianti and badnell data to try to take account of DR in zeta - if atomic data is not read in, then the old interpolated zeta will be returned */


/* NSH 0712 - Old way of doing it, not moved into routines in zeta.c 
	ground_frac[nion - 1].frac[ilow] +
	interpfrac * (ground_frac[nion - 1].frac[ihi] -
		      ground_frac[nion - 1].frac[ilow]); */
      /*Since nion-1 points at state i-1 (saha: n_i/n_i-1) we want ground_frac[nion-1].
         Note that we NEVER access the last ion for any element in that way
         which is just as well since you can't have recombinations INTO
         a fully ionized state -- The corresponding lines in recombination.out 
         are just dummies to make reading the data easier */
      fudge2 = fudge2 + www * (1.0 - fudge2);
      newden[nion] = fudge2 * q;
      sum += newden[nion];
// This is the equation being calculated
//              sum+=newden[nion]=newden[nion-1]*fudge*(*ne)*density[nion]/density[nion-1]/xne;
    }
  a = nh * ele[nelem].abun / sum;
  for (nion = first; nion < last; nion++)
    newden[nion] *= a;



  return (0);
}

/***********************************************************
                       Space Telescope Science Institute

 Synopsis:   
  	fix_concentrations(xplasma,mode) sets the concentrations of individual ions
  	assuming hard coded ionization fractions that are read from a file
  
 Arguments:		


 Returns:
	
Description:	

	The first time the routine is called it reads the file specified by
	geo.con_file which populates the force_con array.

	On subsequent calls the values of one of an xplasma element are 
	modified.

	The format of the file is fairly obvious

	Element.z   Ion   ion_fraction

	6            4       1

	would set the CIV fraction to 1.  
	
	Note there is nothing forcing the totals to be anything sensible.


 	

Notes:
	To change the ions that you want to fix, simply update con_force,
	AND MAKE SURE NFORCE IS SET TO THE NUMBER OF IONS

History:
 	97sep13	ksl	Adapted from concentrations
 	98may	ksl	Moved ne calculation to separate subroutine
	01mar	ksl	Moved specification of fixed concentrations to a file
	04dec	ksl	54e -- previously a specific file fixed.dat was
			read.  With this change, the concentration file
			name is read as part of the inputs
	080802	ksl	60b -- Made fix_concentrations calls resemble
			those of other routines of the same type

**************************************************************/

int fix_con_start = 0;
struct force_con
{
  int z, istate;
  float frac;
}
con_force[10];
int nforce;

int
fix_concentrations (xplasma, mode)
     PlasmaPtr xplasma;
     int mode;			// 0=saha using tr, 1=saha using te, 2= Lucy & Mazzali

{
  int nelem, nion;
  int n;
  double get_ne ();
  FILE *fopen (), *cptr;
  char line[LINELENGTH];

  double nh;

  nh = xplasma->rho * rho2nh;


  /* Define the element and ion which will be present in the wind */


  if (fix_con_start == 0)
    {
      if ((cptr = fopen (geo.fixed_con_file, "r")) == NULL)
	{
	  Error
	    ("fix_concentrations: Could not open %s to read concentrations\n",
	     geo.fixed_con_file);
	  exit (0);
	}

      nforce = 0;
      while (fgets (line, LINELENGTH, cptr) != NULL && nforce < 10)
	{
	  if ((n =
	       sscanf (line, "%d %d %f", &con_force[nforce].z,
		       &con_force[nforce].istate,
		       &con_force[nforce].frac)) != 3)
	    {
	      Error ("fix_concentrations: Read %d items from %s\n", n, line);
	    }
	  Log ("Fixing element %d istate %d to fraction %f\n",
	       con_force[nforce].z, con_force[nforce].istate,
	       con_force[nforce].frac);
	  nforce++;
	}

      fix_con_start = 1;
      fclose (cptr);
    }

  /* Set all the ion abundances to 0 and *ne to 0 */
  for (nion = 0; nion < nions; nion++)
    xplasma->density[nion] = 0;

  /* Search for matches and if one is found set the density of that ion to be the density of that
     atom */
  for (n = 0; n < nforce; n++)
    {
      nion = 0;
      while (nion < nions &&
	     !(ion[nion].z == con_force[n].z
	       && ion[nion].istate == con_force[n].istate))
	nion++;
      if (nion < nions)
	{			/* There was a match */
	  /* Find the element associated with the ion */
	  nelem = 0;
	  while (nelem < nelements && ele[nelem].z != con_force[n].z)
	    nelem++;
	  /* Increment the ion density and the electron density */
	  xplasma->density[nion] = nh * ele[nelem].abun * con_force[n].frac;
	}
    }

  xplasma->ne = get_ne (xplasma->density);

  partition_functions (xplasma, 0);

  return (0);
}




/***********************************************************
                                       Space Telescope Science Institute

  Synopsis:   
	get_ne simple calculates the electron density given the 
	densities of each ion.
  
  Arguments:		


  Returns:
 	
  Description:	

 
  Notes:

  	This makes use of the fact that the densities are stored
	in ion order.

  History:


**************************************************************/

double
get_ne (density)
     double density[];
{
  int n;
  double ne;
  ne = 0;
  for (n = 0; n < nions; n++)
    {
      ne += density[n] * (ion[n].istate - 1);
    }
  return (ne);
}
