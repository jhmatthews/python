/***********************************************************
                                       West Lulworth

  Synopsis:   

int   matrix_ion_populations (xplasma,mode)
  
  Arguments:		
     PlasmaPtr xplasma;		The cell in question -
     int mode;			1 - use a power law and or exponential model for J
     				2 - use a dilute BB model for J - t_r and w are those 
				parameters in xplasma


  Returns:
	

 	
  Description:	

	The general structure of matrix_ion_populations is as follows:
	First we compute all of the rates that we have data for linking all ionization stages 
	We make an initial guess at the electron density
	We attempt to solve the ionization rates matrix, and produce a new guess at the electron density
	We now calculate new rates, resolve and proceed until the electron density converges.


 
  Notes:

  This routine is a first try at implementing a matrix solver for the ionization state in a cell. 
	This was motivated by the desire to incorporate more ionization and recombination processes
	into python. 

	



  History:
	2014Aug NSH - coded
	2014 Nov NSH - tidied up

**************************************************************/



#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "atomic.h"
// struct topbase_phot *xtop;
// double g1, g2, xne, xip;
#define SAHA 4.82907e15
// double xip, xxxne, qromb_temp;
#include "python.h"


#include <gsl/gsl_block.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
// #include <gsl/gsl_blas.h>
#include "my_linalg.h"




int
matrix_ion_populations (xplasma, mode)
     PlasmaPtr xplasma;
     int mode;

{
  int nn, mm, nrows;
  double rate_matrix[nions][nions];
  double newden[NIONS];
  double nh, t_e, t_r, www;
  double xne, xxne, xxxne;
  //double xsaha, x, theta;
  //int s;                                                                                                
  double b_temp[nions];
  double *b_data, *a_data;
  double *populations;
  int ierr, niterate;
  double xnew;
  //double xden[nelements];
  int xion[nions];		// This array keeps track of what ion is in each line
  int xelem[nions];		// This array keeps track of the element for each ion
  double pi_rates[nions];
  double rr_rates[nions];
  double inner_rates[n_inner_tot]; //This array contains the rates for each of the inner shells. Where they go to requires the electron yield array

  /* Copy some quantities from the cell into local variables */

  nh = xplasma->rho * rho2nh;	// The number density of hydrogen ions - computed from density
  t_e = xplasma->t_e;		// The electron temperature in the cell - used for collisional processes 
  t_r = xplasma->t_r;		// The radiation temperature - used for PI if we have a BB approximation
  www = xplasma->w;		// The radiative weight in the cell - again for BB approximation for PI

  /* Dielectronic recombination and direct ionization coefficients depend only on electron temperature, calculate them now -
     they will not change */

  compute_dr_coeffs (t_e);
  compute_di_coeffs (t_e);

  /* JM 1508 -- also compute direct recombination coefficients */
  compute_qrecomb_coeffs(t_e);


  /* In the following loop, over all ions in the simulation, we compute the radiative recombination rates, and photionization
     rates OUT OF each ionization stage. The PI rates are calculated either using the modelled mean intensity in a cell, or
     using the dilute blackbody approximation, depending on which mode we are in. At the same time, we copy the ion densities
     from the plasma strcuture into a local array. We will only overwrite the numbders in the structure if we believe the
     results are an improvement on what is there. */

  for (mm = 0; mm < nions; mm++)
    {
      newden[mm] = xplasma->density[mm];	// newden is our local density array
      xion[mm] = mm;		// xion is an array we use to track which ion is in which row of the matrix
      if (ion[mm].istate != 1)	// We can recombine since we are not in the first ionization stage
	{
	  rr_rates[mm] = total_rrate (mm, xplasma->t_e);	// radiative recombination rates
	}
      if (ion[mm].istate != ion[mm].z + 1)	// we can photoionize, since we are not in the z+1th ionization state (bare)
	{
	  if (mode == NEBULARMODE_MATRIX_BB)
	    {
	      pi_rates[mm] = calc_pi_rate (mm, xplasma, 2, 1);	// PI rate for the BB model
	    }
	  else if (mode == NEBULARMODE_MATRIX_SPECTRALMODEL)
	    {
	      pi_rates[mm] = calc_pi_rate (mm, xplasma, 1, 1);	// PI rate for an explicit spectral model
	    }
	  else
	    {
	      // If reached this point the program does not understand what type of spectral model to apply 
	      Error ("matrix_ion_populations: Unknown mode %d\n", mode);
	      exit (0);
	    }		
	}

    for (nn = 0; nn < nelements; nn++)
{
  if (ion[mm].z == ele[nn].z)
    {
      xelem[mm] = nn;	/* xelem logs which element each ow in the arrays refers to. This is important because we need
			   to know what the total density will be for a group of rows all representing the same
			   element. */
    }
	}
  }
	
	
	/* The next loop generates the inner shell ionization rates, if they are present in the atomic data and
	we wish to compute auger ionizaion rates. This only computes the rates out of each ion, we also need to 
    consult the electron yield array if we are to compute the change in electron number*/
  
	if (geo.auger_ionization==1)
	{
		for (mm=0; mm<n_inner_tot; mm++)
		{
			if (mode==NEBULARMODE_MATRIX_BB)
			{
				inner_rates[mm]=calc_pi_rate(mm, xplasma, 2, 2);
			}
			else if (mode==NEBULARMODE_MATRIX_SPECTRALMODEL)
			{
				inner_rates[mm]=calc_pi_rate(mm,xplasma,1,2);
			}
		}
	}

	
		


  /* This next line sets the partition function for each ion. This has always been the place here python calculates the
     partition funcions and sets the level densities for each ion. It needs to be done, or other parts of the code which rely
     on sensible level populations don't work properly. In the case of the dilute blackbody, the code works well, however we do 
     not currently (v78) have a procedure to calucate the levels for a spectral model case. We therefore call partition
     functions with mode 4 - this is a special mode which forces the ions to be in the ground state. This is reasonable for a
     radiation dominated plasma, since any excitations into higher states will quickly be followed by radative de-excitation.
     We should really do better though... */

  if (mode == NEBULARMODE_MATRIX_BB)
    {
      partition_functions (xplasma, NEBULARMODE_ML93);	// We use t_r and the radiative weight
    }
  else if (mode == NEBULARMODE_MATRIX_SPECTRALMODEL)
    {
      partition_functions (xplasma, 4);	// Set to ground state 
    }

  /* Next we need to obtain an initian guess for the electron density. In the past this has been done by calculating the
     hydrogen density directly from the Saha equation at the current electron temperature. In initial testing of this mode -
     this seemed to be a little unstable. At the moment, we use the last ionization state to compute n_e. For the first
     iteraion, this will just be calculated in LTE from the initial temperature structure of the wind, which will give almost
     the same result as the original procedure, or for successive calculations, it should be a better guess. I've leftin the
     original code, commented out...  */

  xne = xxne = xxxne = get_ne (newden);	// Set n_e to the current value. 

  /* xne is the current working number xxne */


  /* These following commented out lines calculate the electron density from Saha given the current electron temperature */
  // if (t_e < MIN_TEMP)
  // t_e = MIN_TEMP;  fudge to prevent divide by zeros 
  // xsaha = SAHA * pow (t_e, 1.5);
  // theta = xsaha * exp (-ion[0].ip / (BOLTZMANN * t_e)) / nh;
  // if (theta < THETAMAX)
  // {
  // x = (-theta + sqrt (theta * theta + 4 * theta)) / 2.;
  // xne = xxne = xxxne = x * nh;
  // }
  // else
  // xne = xxne = xxxne = nh; /*xxne is just a store so the error can report the starting value of ne. */
  /* xxxne is the shared variable so the temperature solver routine can access it */
  // if (xne < 1.e-6)
  // xne = xxxne = 1.e-6; /* Set a minimum ne to assure we can calculate
  // xne the first time through the loop 



  /* We are now going to iterate on the electron density - MAXITERATIONS is set in python.h and is currently (78) set to 200.
     We would normally expect to converge much fater than this */

  niterate = 0;
  while (niterate < MAXITERATIONS)
    {


      populate_ion_rate_matrix (xplasma, rate_matrix, pi_rates, inner_rates, rr_rates,
				b_temp, xne, xelem);


      /* The array is now fully populated, and we can begin the process of solving it */

      nrows = nions;		/* This is a placeholder, we may end up removing rows and columns that have no rates (or very
				   low rates) connecting them to other rows. This may improve stability but will need to be
				   done carefully */

      /* The solver routine was taken largely wholesale from the matom routine. I have left in most of the original comments, and 
         added a few of my own for clarification */



		/********************************************************************************/
      /* The block that follows (down to next line of ***s) is to do the matrix inversion. It uses LU decomposition - the code
         for doing this is taken from the GSL manual with very few modifications. */
      /* here we solve the matrix equation M x = b, where x is our vector containing level populations as a fraction w.r.t the
         whole element */
      /* JM 1411 -- the actual LU decomposition- the process of obtaining a solution - is done by the routine solve_matrix() */

      /* Replaced inline array allocaation with calloc, which will work with older version of c compilers */

      /* This next line produces an array of the correct size to hold the rate matrix */
	
      a_data = (double *) calloc (sizeof (double), nrows * nrows);

      /* We now copy our rate matrix into the prepared matrix */
      for (mm = 0; mm < nrows; mm++)
	{
	  for (nn = 0; nn < nrows; nn++)
	    {
	      a_data[mm * nrows + nn] = rate_matrix[mm][nn];
	    }
	}



      /* Replaced inline array allocaation with calloc, which will work with older version of c compilers calloc also sets the
         elements to zero, which is required */

      b_data = (double *) calloc (sizeof (double), nrows);
      populations = (double *) calloc (sizeof (double), nrows);

      /* This b_data column matrix is the total number density for each element, placed into the row which relates to the neutral 
         ion. This matches the row in the rate matrix which is just 1 1 1 1 for all stages. NB, we could have chosen any line for 
         this. */

      for (nn = 0; nn < nrows; nn++)
	{
	  b_data[nn] = b_temp[nn];
	}

	     ierr = solve_matrix (a_data, b_data, nrows, populations, xplasma->nplasma);

      if (ierr != 0)
	      Error ("matrix_ion_populations: bad return from solve_matrix\n",ierr);
	    if (ierr ==2)
		    Error ("matrix_ion_populations: some matrix rows failing relative error check\n");
	    else if (ierr ==3)
			  Error ("matrix_ion_populations: some matrix rows failing absolute error check\n");

      /* free memory */
      free (a_data);
      free (b_data);

      /* Calculate level populations for macro-atoms */
      if (geo.macro_ioniz_mode == 1)
	{
	  macro_pops (xplasma, xne);
	}

      /* We now have the populations of all the ions stored in the matrix populations. We copy this data into the newden array
         which will temperarily store all the populations. We wont copy this to the plasma structure until we are sure thatwe
         have made things better. We just loop over all ions, and copy. The complexity in the loop is to future proof us against
         the possibility that there are some ions that are not included in the matrix scheme becuse there is no route into or
         out of it. */

      for (nn = 0; nn < nions; nn++)
	{
	  newden[nn] = 0.0;	// initialise the arrayelement

	  /* if the ion is being treated by macro_pops then use the populations just computed */
	  if ((ion[nn].macro_info == 1) && (geo.macro_simple == 0)
	      && (geo.macro_ioniz_mode == 1))
	    {
	      newden[nn] = xplasma->density[nn];
	    }

	  for (mm = 0; mm < nrows; mm++)	// inner loop over the elements of the population array
	    {
	      if (xion[mm] == nn)	// if this element contains the population of the ion is question
		{
		  newden[nn] = populations[mm];	// get the population
		}
	    }


	  if (newden[nn] < DENSITY_MIN)	// this wil also capture the case where population doesnt have a value for this ion
	    newden[nn] = DENSITY_MIN;
	}
	free (populations);
      xnew = get_ne (newden);	/* determine the electron density for this density distribution */


      if (xnew < DENSITY_MIN)
	xnew = DENSITY_MIN;	/* fudge to keep a floor on ne */


      if (fabs ((xne - xnew) / (xnew)) < FRACTIONAL_ERROR || xnew < 1.e-6)	/* We have converged, or have the situation where we
										   have a neutral plasma */
	{
	  break;		/* Break out of the while loop - we have finished our iterations */
	}
      xne = xxxne = (xnew + xne) / 2.;	/* New value of ne */

      niterate++;


      if (niterate == MAXITERATIONS)
	{
	  Error
	    ("matrix_ion_populations: failed to converge for cell %i t %e nh %e xnew %e\n",
	     xplasma->nplasma, t_e, nh, xnew);


	  for (nn = 0; nn < geo.nxfreq; nn++)
	    {
	      Log
		("numin= %e (%e) numax= %e (%e) Model= %2d PL_log_w= %e PL_alpha= %e Exp_w= %e EXP_temp= %e\n",
		 xplasma->fmin_mod[nn], geo.xfreq[nn], xplasma->fmax_mod[nn],
		 geo.xfreq[nn + 1], xplasma->spec_mod_type[nn],
		 xplasma->pl_log_w[nn], xplasma->pl_alpha[nn],
		 xplasma->exp_w[nn], xplasma->exp_temp[nn]);
	    }

	  Error ("matrix_ion_populations: xxne %e theta %e\n", xxne);

	  return (-1);		/* If we get to MAXITERATIONS, we return without copying the new populations into plasma */
	}
    }				/* This is the end of the iteration loop */


  xplasma->ne = xnew;
  for (nn = 0; nn < nions; nn++)
    {
      /* If statement added here to suppress interference with macro populations (SS Apr 04) */
      if (ion[nn].macro_info == 0 || geo.macro_ioniz_mode == 0
	  || geo.macro_simple == 1)
	{
	  xplasma->density[nn] = newden[nn];
	}
      if ( (sane_check(xplasma->density[nn])) || (xplasma->density[nn] < 0.0) )
        Error("matrix_ion_populations: ion %i has population %8.4e in cell %i\n", 
               nn, xplasma->density[nn], xplasma->nplasma);
    }


  partition_functions (xplasma, 4);	/* WARNING fudge NSH 11/5/14 - this is as a test. We really need a better implementation
					   of partition functions and levels for a power law illuminating spectrum. We found that
					   if we didnt make this call, we would end up with undefined levels - which did really
					   crazy things */



  return (0);
}


/***********************************************************
                                       West Lulworth

  Synopsis:   
    populate_ion_rate_matrix populates a rate_matrix of shape nions x nions
    with the pi_rates and rr_rates supplied at the density xne in question.

  
  Arguments:	
    xplasma 
      PlasmaPtr to the plasma cell for physical properties

    pi_rates
      PI rates for each ion

    rr_rates 
      radiative recombination rates for each ion

    b_temp
      array which is mostly zeros, but first element is set to abundances
      so the equation is soluble

    xne 
      our guess of ne, which has not been copied to xplasma yet

    xelem
      array which stores the element corresponding to each ion in the matrix

  Returns:
	

 	
  Description:
 
  Notes:
    This routine includes the process of replacing the first row of the matrix with
    1s in order to make the problem soluble. 

  History:
	2014Aug JM - moved code here from main routine

**************************************************************/

int
populate_ion_rate_matrix (xplasma, rate_matrix, pi_rates, inner_rates, rr_rates, b_temp,
			  xne, xelem)
     PlasmaPtr xplasma;
     double rate_matrix[nions][nions];
     double pi_rates[nions];
	 double inner_rates[n_inner_tot];
     double rr_rates[nions];
     double xne;
     double b_temp[nions];
     int xelem[nions];

{
  int nn, mm;
  double nh;
 int n_elec,d_elec,ion_out; //The number of electrons left in a current ion


  nh = xplasma->rho * rho2nh;	// The number density of hydrogen ions - computed from density

  /* First we initialise the matrix */
  for (nn = 0; nn < nions; nn++)
    {
      for (mm = 0; mm < nions; mm++)
	{
	  rate_matrix[nn][mm] = 0.0;
	}
    }


  /* The next block of loops populate the matrix. For simplicity of reading the code each process has its own loop. Some rates
     actually dont change during each iteration, but those that depend on n_e will. All are dealt with together at the moment,
     but this could be streamlined if it turns out that there is a bottleneck. */

  /* Now we populate the elements relating to PI depopulating a state */


  for (mm = 0; mm < nions; mm++)
    {
      if (ion[mm].istate != ion[mm].z + 1)	// we have electrons
	{
	  rate_matrix[mm][mm] -= pi_rates[mm];
	}
    }


  /* Now we populate the elements relating to PI populating a state */

  for (mm = 0; mm < nions; mm++)
    {
      for (nn = 0; nn < nions; nn++)
	{

	  if (mm == nn + 1 && ion[nn].istate != ion[nn].z + 1
	      && ion[mm].z == ion[nn].z)
	    {
	      rate_matrix[mm][nn] += pi_rates[nn];
	    }
	}
    }

  /* Now we populate the elements relating to direct ionization depopulating a state */

  for (mm = 0; mm < nions; mm++)
    {
      if (ion[mm].istate != ion[mm].z + 1 && ion[mm].dere_di_flag > 0)	// we have electrons and a DI rate
	{
	  rate_matrix[mm][mm] -= (xne * di_coeffs[mm]);
	}
    }

  /* Now we populate the elements relating to direct ionization populating a state - this does depend on the electron density */

  for (mm = 0; mm < nions; mm++)
    {
      for (nn = 0; nn < nions; nn++)
	{
	  if (mm == nn + 1 && ion[nn].istate != ion[nn].z + 1
	      && ion[mm].z == ion[nn].z && ion[nn].dere_di_flag > 0)
	    {
	      rate_matrix[mm][nn] += (xne * di_coeffs[nn]);
	    }
	}
    }


  /* Now we populate the elements relating to radiative recomb depopulating a state */

  for (mm = 0; mm < nions; mm++)
    {
      if (ion[mm].istate != 1)	// we have space for electrons
	{
	  rate_matrix[mm][mm] -= xne * (rr_rates[mm] + xne * qrecomb_coeffs[mm]);
	}
    }


  /* Now we populate the elements relating to radiative recomb populating a state */

  for (mm = 0; mm < nions; mm++)
    {
      for (nn = 0; nn < nions; nn++)
	{
	  if (mm == nn - 1 && ion[nn].istate != 1 && ion[mm].z == ion[nn].z)
	    {
	      rate_matrix[mm][nn] += xne * (rr_rates[nn] + xne * qrecomb_coeffs[nn]);
	    }
	}
    }

  /* Now we populate the elements relating to dielectronic recombination depopulating a state */

  for (mm = 0; mm < nions; mm++)
    {
      if (ion[mm].istate != 1 && ion[mm].drflag > 0)	// we have space for electrons
	{
	  rate_matrix[mm][mm] -= (xne * dr_coeffs[mm]);
	}
    }


  /* Now we populate the elements relating to dielectronic recombination populating a state */



  for (mm = 0; mm < nions; mm++)
    {
      for (nn = 0; nn < nions; nn++)
	{
	  if (mm == nn - 1 && ion[nn].istate != 1 && ion[mm].z == ion[nn].z
	      && ion[mm].drflag > 0)
	    {
	      rate_matrix[mm][nn] += (xne * dr_coeffs[nn]);
	    }
	}
    }
	
	
	
	for (mm=0; mm<n_inner_tot; mm++)  //There mare be several rates for each ion, so we loop over all the rates
	{
		if (inner_cross[mm].n_elec_yield != -1)  //we only want to treat ionization where we have info about the yield
		{
			ion_out=inner_cross[mm].nion;  //this is the ion which is being depopulated
			rate_matrix[ion_out][ion_out]-=inner_rates[mm]; //This is the depopulation
			n_elec=ion[ion_out].z-ion[ion_out].istate+1;
			if (n_elec>11)
				n_elec=11;
			for (d_elec=1;d_elec<n_elec;d_elec++) //We do a loop over the number of remaining electrons
			{
				nn=ion_out+d_elec; //We will be populating a state d_elec stages higher
				rate_matrix[nn][ion_out] += inner_rates[mm]*inner_elec_yield[inner_cross[mm].n_elec_yield].prob[d_elec-1];
			}
		}
	}
	



  /* Now, we replace the first line for each element with 1's and 0's. This is done because we actually have more equations
     than unknowns. We replace the array elements relating to each ion stage in this element with a 1, and all the other array
     elements (which relate to other elements (in the chemicalsense) with zeros. This is equivalent to the equation
     1*n1+1*n2+1*n3 = n_total - i.e. the sum of all the partial number densities adds up to the total number density for that
     element. This loop also produces the 'b matrix'. This is the right hand side of the matrix equation, and represents the
     total number density for each element. This can be a series of 1's for each row in the matrix relating to the ground
     state of the repsective element, however one can just use the total number density for that elements and then the
     densities which the solver computes are just the actual number densities for each ion. This does mean that different rows
     are orders of magnitude different, so one can imagine numerical issues. However each element is essentially solved for
     seperately.. Something to watch */

  for (nn = 0; nn < nions; nn++)
    {
      if (ion[nn].istate == 1)
	{
	  b_temp[nn] = nh * ele[xelem[nn]].abun;
	  for (mm = 0; mm < nions; mm++)
	    {
	      if (ion[mm].z == ion[nn].z)
		{
		  rate_matrix[nn][mm] = 1.0;
		}
	      else
		{
		  rate_matrix[nn][mm] = 0.0;
		}
	    }
	}
      else
	{
	  b_temp[nn] = 0.0;
	}
    }

  return (0);
}

/***********************************************************
                                       AMNH, New York

  Synopsis: 
    solve_matrix solves the matrix equation A x = b for the vector x.  
     
  Arguments:
    a_data 
    	1d array of length nrows * nrows, containing entries
    	for matrix A	

    b_data
    	1d array of length nrows, containing entries
    	for vector b

    nrows
        dimensions of matrices/vectors. Should be consistent
        with arrays passed to routine.

    populations	 
    	1d array of length nrows, containing entries
        for vector x - modified by this routine		

  Returns:
 	
  Description:
 
  Notes:
    
  History:
	2014Aug JM - moved code here from main routine

**************************************************************/

int
solve_matrix (a_data, b_data, nrows, x, nplasma)
     double *a_data, *b_data;
     int nrows;
     double *x;
     int nplasma;
{
  int mm, ierr, s;
  /* s is the 'sign' of the permutation - is had the value -1^n where n is the number of
     permutations. We dont use it anywhere, but in principle it can be used to refine the
     solution via gsl_linalg_LU_refine */
  double test_val;
  double det;
  
  gsl_permutation *p;
  gsl_matrix_view m;
  gsl_vector_view b;
  gsl_vector *test_vector, *populations;
  gsl_matrix *test_matrix;

  ierr = 0;
  test_val = 0.0;

  /* create gsl matrix/vector views of the arrays of rates */
  m = gsl_matrix_view_array (a_data, nrows, nrows);

  /* these are used for testing the solution below */
      test_matrix = gsl_matrix_alloc (nrows, nrows);
	  test_vector = gsl_vector_alloc (nrows);

	   gsl_matrix_memcpy (test_matrix, &m.matrix);	// create copy for testing 


  b = gsl_vector_view_array (b_data, nrows);

  /* the populations vector will be a gsl vector which stores populations */
  populations = gsl_vector_alloc (nrows);


  p = gsl_permutation_alloc (nrows);	// NEWKSL

  gsl_linalg_LU_decomp (&m.matrix, p, &s);

  det = gsl_linalg_LU_det(&m.matrix, s); // get the determinant to report to user

  if (det == 0)
    Error("Rate Matrix Determinant is %8.4e for cell %i", det, nplasma);

  gsl_linalg_LU_solve (&m.matrix, p, &b.vector, populations);

  gsl_permutation_free (p);

  /* JM 140414 -- before we clean, we should check that the populations vector we have just created really is a solution to
     the matrix equation */

  /* gsl_blas_dgemv declaration contained in my_linalg.h, taken from gsl library */
  /* The following line does the matrix multiplication test_vector = 1.0 * test_matrix * populations The CblasNoTrans
     statement just says we do not do anything to test_matrix, and the 0.0 means we do not add a second matrix to the result
     If the solution has worked, then test_vector should be equal to b_temp */

    ierr =
		  gsl_blas_dgemv (CblasNoTrans, 1.0, test_matrix, populations, 0.0,
			    test_vector);

				if (ierr != 0)
				{
					Error
					("solve_matrix: bad return when testing matrix solution to rate equations.\n");
				}

  /* now cycle through and check the solution to y = m * populations really is (1, 0, 0 ... 0) */

				for (mm = 0; mm < nrows; mm++)
				{

      /* get the element of the vector we want to check */
					test_val = gsl_vector_get (test_vector, mm);

      /* b_data is (1,0,0,0..) when we do matom rates. test_val is normally something like
         1e-16 if it's supposed to be 0. We have a different error check if b_data[mm] is 0 */

					if (b_data[mm] > 0.0)
					{
						if (fabs ((test_val - b_data[mm]) / test_val) > EPSILON)
						{
							Error("solve_matrix: test solution fails relative error for row %i %e != %e\n",
								mm, test_val, b_data[mm]);
							ierr = 2;
						}
					}
					else if (fabs (test_val - b_data[mm]) > EPSILON)	// if b_data is 0, check absolute error
						
					{
						Error("solve_matrix: test solution fails absolute error for row %i %e != %e\n",
							mm, test_val, b_data[mm]);
						ierr = 3;
					}
				}

  /* copy the populations to a normal array */
  for (mm = 0; mm < nrows; mm++)
    x[mm] = gsl_vector_get (populations, mm);

  /* free memory */
  gsl_vector_free (test_vector);
  gsl_matrix_free (test_matrix);
  gsl_vector_free (populations);

  return (ierr);
}
