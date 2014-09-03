
/**************************************************************************
                    Space Telescope Science Institute


  Synopsis:
  	The routines here provide a way to read grids of models and
	then to interpolate between them.

	There are two main routines

		get_models (modellist, npars, spectype)

	which reads a list of models

	and
		model (spectype, par)

	which given a set of parameters will interpolate between modles
	given the set of models one wants to interpolate betwee, which
	is determined by spectype, and the set of parameters (par) that
	define the models
 

  History:
	13jul	ksl	Added these notes about the purpose of the routines in this file
************************************************************************/
 


/**************************************************************************
                    Space Telescope Science Institute


  Synopsis:

	get_models (modellist, npars, spectype)


  Description:	These are a generic set of routines to read a grid or grids
	of models, e.g. a set of Kurucz stellar atmospheres which can
	used to simulate a disk, given a run of temperature and gravity
   	
	
	The grid can be any reasonable number of dimension

  Arguments:		
     char modellist[];		// filename containing location and associated parameters of models
     int npars;			// Number of parameters which vary for these models

     The modellist is an ascii file listing one model on each line of the file. Comments
     lines begining with # are ignored.  For example

     kurucz91/fp00t4000g30k2c125.txt         4000          3.0

     refers to a Kurucz model with a temperature of 4000 K and a log gravity of 3.0


  Returns:
     int *spectype;		//  The returned spectrum type 

     spectype really is just an integer that is incremented each time a new set of models
     is read in.  

  Notes:
	DO NOT INCLUDE ANYTHING THAT IS NOT GENERIC HERE!  

	The models general are also ascii files with wavelengths and fluxes (flambda
	or something proportional to flambda, although since the routine is generic
	that is not something that is enforced by this routine.  

	For 1d models we assume the models are read in increasing order of
	the parameter of interest.  

	2 or more parameter models can be read in any order.  Bilinear interpolation
	is done for the 2 parmeter models (as long as the grid is fully filled).
	To deal with cases in which the grid is not fully filled, we follow 
	the same approach as for get_hub, i.e. for parameters t and g (i.e.
	the first parameter being t, and the second g., we first search the grid
	for all models with t just below (above) the desired value.  We then 
	interpolate on g for the model just below, and on g for the model just
	above, and finally interpolate between these two.  
	The algorithm for interpolating between models 


	Note that each time a new set of models is read in, the spectype is
	incremented.

	These routines were originally parrt of as set of routines, kslfit,  intended
	to fit UV models to UV spectra. Subsequently, 080t17, the kslfit routines of the
	same name were modified and so the two sets of routines are no longer
	in sync, even though they have the same name.



  History:
04jul	ksl	Adapted for use in python
04dec	ksl	54a -- miniscule change to eliminate warning with 03
13jul	ksl	Another small change to elimated a warning on some compileres
	

 ************************************************************************/
#include	<math.h>
#include	<stdio.h>
#include        <stdlib.h>
#include	<strings.h>
#include	<string.h>
#include 	"atomic.h"
#include	"python.h"	//This needs to come before modlel.h so that what is in models.h is used
#include         "models.h"
#define    	BIG 1e32


/* Get all the models of one type and regrid them onto the wavelength grid of the data */
int get_models_init = 0;	// needed so can initialize nmods_tot the first time this routine is called

int
get_models (modellist, npars, spectype)
     char modellist[];		// filename containing location and associated parameters of models
     int npars;			// Number of parameters which vary for these models
     int *spectype;		//  The returned spectrum type 


{
  FILE *mptr, *fopen ();
  char dummy[LINELEN];
  int n, m, mm, nxpar;
  double xpar[NPARS], xmin[NPARS], xmax[NPARS];
  int get_one_model ();
  int nw, nwaves;


  nwaves = 0;

  if (get_models_init == 0)
    {
      nmods_tot = 0;
      ncomps = 0;		// The number of different sets of models that have been read in
      get_models_init = 1;
    }

  /* Now check if this set of models has been read in previously.  If so return the 
   * spectype when this was read.
   */

  n = 0;
  while (n < ncomps)
    {
      if (strcmp (modellist, comp[n].name) == 0)
	{
	  *spectype = n;
	  return (0);
	}
      n++;
    }


  if ((mptr = fopen (modellist, "r")) == NULL)
    {
      Error ("get_models:Could not open file %s containing list of models \n",
	     modellist);
      exit (0);
    }

/* Now initialize the model summary structure */
  strcpy (comp[ncomps].name, modellist);
  comp[ncomps].modstart = nmods_tot;
  comp[ncomps].npars = npars;
  comp[ncomps].xpdf.limit1 = -99.;
  comp[ncomps].xpdf.limit2 = -99.;

/* Now get all the models of this type */
  n = nmods_tot;		// This is the starting point since may have read models in before
/* Establish initial limits on xmin and xmax so that they can be properly populated */
  for (m = 0; m < NPARS; m++)
    {
      xmin[m] = (BIG);
      xmax[m] = (-BIG);
      comp[ncomps].xmod.par[m] = -99;
    }

  nw = -1;			// Initiallize nw
  while (n < NMODS && (fgets (dummy, LINELENGTH, mptr)) != NULL)
    {
      if (dummy[0] == '#' || dummy[0] == '!')
	{
	}			//skip comment lines in models
      else
	{
	  nxpar =
	    sscanf (dummy, "%s %lf %lf %lf %lf %lf %lf %lf %lf %lf",
		    mods[n].name, &xpar[0], &xpar[1], &xpar[2], &xpar[3],
		    &xpar[4], &xpar[5], &xpar[6], &xpar[7], &xpar[8]);
	  if (nxpar < npars)
	    {
	      Error ("get_models: nxpar (%d) < npars (%d) in line %s\n",
		     nxpar, npars, dummy);
	      exit (0);
	    }
	  for (m = 0; m < npars; m++)
	    {
	      mods[n].par[m] = xpar[m];
	      if (xpar[m] > xmax[m])
		xmax[m] = xpar[m];
	      if (xpar[m] < xmin[m])
		xmin[m] = xpar[m];
	    }
	  //130727 - ksl - added a new variable mm to prevent a warning thrown by some compilers 
	  for (mm = m; mm < NPARS; mm++)
	    mods[n].par[mm] = -99;

	  nwaves = get_one_model (mods[n].name, &mods[n]);
	  if (nw > 0 && nwaves != nw)
	    {
	      Error
		("get_models: file %s has %d wavelengths, others have %d\n",
		 mods[n].name, nwaves, nw);
	      exit (0);
	    }

	  if ((n % 100) == 0)
	    printf ("Model n %d %s\n", n, mods[n].name);
	  n++;
	}
    }

  if (n == NMODS)
    {
      Error
	("get_models: Reached maximum number of models %d. Please increase NMODS in .h file \n",
	 n);
      exit (0);
    }
/* Now complete the initialization of the modsum structure */
  comp[ncomps].modstop = nmods_tot = n;
  comp[ncomps].nmods = comp[ncomps].modstop - comp[ncomps].modstart;
  comp[ncomps].nwaves = nwaves;
  for (n = 0; n < nwaves; n++)
    {
      comp[ncomps].xmod.w[n] = mods[comp[ncomps].modstart].w[n];
    }

  if (comp[ncomps].nmods == 0)
    {
      Error
	("get_models: No models from %s were read. Please check list of models!\n",
	 comp[ncomps].name);
      exit (0);
    }

  /* The next 3 lines set a normalization that is used by kslfit.  They are mostly
   * not relevant to python, where comp[ncomp[.min[0] refers to a normalization for 
   * a model.  I've kept this initialization for now */
  comp[ncomps].min[0] = 0;
  comp[ncomps].max[0] = 1000;

  for (m = 0; m < npars; m++)
    {
      comp[ncomps].min[m] = xmin[m];
      comp[ncomps].max[m] = xmax[m];
    }

  *spectype = ncomps;		// Set the spectype 
  ncomps++;
  return (*spectype);
}

/* Get a single model model 
   This routine simple reads a model file from disk and puts the result into the structure
   onemod.  This file need not have the same wavelengths as the data nor other models.

080915  ksl     Added error to catch the case where the model being read in has more
                wavelengths than allowed, ie. more than NWAVES
 */
int
get_one_model (filename, onemod)
     char filename[];
     struct Model *onemod;
{
  FILE *ptr;
  char dummy[LINELEN];
  int n;
  double w, f;

  if ((ptr = fopen (filename, "r")) == NULL)
    {
      Error ("Could not open filename %s\n", filename);
      exit (0);
    }
  n = 0;
  while (n < NWAVES && (fgets (dummy, LINELEN, ptr)) != NULL)
    {
      if ((dummy[0] != '#'))
	{
	  sscanf (dummy, "%le %le", &w, &f);
	  onemod->w[n] = w;
	  onemod->f[n] = f;

	  n++;
	}
    }
  onemod->nwaves = n;


  if (n >= NWAVES)
    {
      Error ("get_one_model: model %s has more than %d wavelengths\n",
	     filename, NWAVES);
    }


  fclose (ptr);
  return (n);
}

/**************************************************************************
                    Space Telescope Science Institute


  Synopsis:  model interpolates between models and places the results in 
  	comp[spectype].xmod

  Description:	

  Arguments:		

  Returns:
  	Returns nwaves in the the spectrum if a new spectrum was calculated
	0 if it was the old spectrum.

  Notes:

  This routine has always been a kluge. In 02aug, I have attempted to improve it.
  The problem has always been what to do in situations where the model grid does
  not fill out a regular, rectangular grid in parameter space.  This arises, for 
  example, in a model grid of stellar atmospheres with T & g, since the range of
  g for which a model can be calculated (i.e. for which Tlusty converges) depends 
  on T.  

  The solution which I try to implement in the new version of the code is to attempt
  to find where models should exist in a completely filled grid.  By completely filled,
  I mean that if you are fitting t and g, you will have 4 models which are just 
  above and below the desired t and g, e.g. mod(t1,g1), mod(t1,g2), mod(t2,g1), 
  mod(t2,g2).  In that case we carry out bi-linear interpolation to find the
  weights.  The problem arises when mod(t2,g2) does not exist.

  The "practical" solution this routine now implements is to ignore the missing model,
  treating each of the others as if bi-linear interpolation was being used, and
  then renormalizing to make sure the weights of the models that one does have
  sum to 1.  The solution is "practical" because it always gives a result, and
  because it gives the an appropriate answer when the grid is filled.  It is also
  easily extensible to more than two dimensions.  But it should be noted, the
  approach is also dangerous.

  History:
04aug	ksl	Adapted from similar routine used by kslfit.  In this
		case the results are placed in comp[spectype]xmod
06jul	ksl	57g -- fixed bug when trying to retrieve the best model
		when the first parameter (T) was exactly at a model value
		but this modeld would be excluded by the gravity calculation.
		Also did a little clean-up of comments

 ************************************************************************/
int nmodel_error = 0;
int nmodel_terror = 0;


int
model (spectype, par)
     int spectype;
     double par[];


{
  int j, n;
  int good_models[NMODS];	// Used to establish which models are to be included in creating output model
  double xmin[NPARS], xmax[NPARS];	// The vertices of a completely filled grid
  double weight[NMODS];		// The weights assigned to the models
  double hi, lo, delta, wtot;
  int ngood;
  double f;
  int nwaves;
  double flux[NWAVES];
  double q1, q2, lambda, tscale, xxx;	// Used for rescaleing according to a bb



/* First determine whether we already have interpolated this
exact model previously */

  n = 0;
  while (n < comp[spectype].npars && comp[spectype].xmod.par[n] == par[n])
    {
      n++;
    }
  if (n == comp[spectype].npars)
    {
      return (0);		// This was the model stored in comp already
    }


  /* First identify the models of interest */
  n = 0;
  while (n < comp[spectype].modstart)
    {
      weight[n] = good_models[n] = 0;
      n++;
    }
  while (n < comp[spectype].modstop)
    {
      weight[n] = good_models[n] = 1;
      n++;
    }
  while (n < nmods_tot)
    {
      weight[n] = good_models[n] = 0;
      n++;
    }

  for (j = 0; j < comp[spectype].npars; j++)
    {
      xmax[j] = comp[spectype].max[j];
      xmin[j] = comp[spectype].min[j];
      hi = BIG;
      lo = -BIG;
      for (n = comp[spectype].modstart; n < comp[spectype].modstop; n++)
	{
	  if (good_models[n])
	    {
	      delta = mods[n].par[j] - par[j];
	      if (delta > 0.0 && delta < hi)
		{
		  xmax[j] = mods[n].par[j];
		  hi = delta;
		}
	      if (delta <= 0.0 && delta >= lo)
		{
		  xmin[j] = mods[n].par[j];
		  lo = delta;
		}
	    }
	}
      // So at this point we know what xmin[j] and xmax[j] and we 
      // need to prune good_models
      for (n = comp[spectype].modstart; n < comp[spectype].modstop; n++)
	{
	  // Next lines excludes the models which are out of range.
	  if (mods[n].par[j] > xmax[j] || mods[n].par[j] < xmin[j])
	    good_models[n] = 0;
	  // Next line modifies the weight of this model assuming a regular grid
	  // If xmax==xmin, then par[j] was outside of the range of the models and
	  // so we need to weight the remaining models fully.
	  if (good_models[n] && xmax[j] > xmin[j])
	    {
	      f = (par[j] - xmin[j]) / (xmax[j] - xmin[j]);
	      if (mods[n].par[j] == xmax[j])
		{
		  // Then the model is at the maximum for this parameter
		  weight[n] *= f;
		}
	      else
		weight[n] *= (1. - f);

/* 57g -- If the weight given to a model is going to be zero, it needs to be
excluded from furthur consideration -- 07jul ksl */
	      if (weight[n] == 0.0)
		good_models[n] = 0;

	    }
	}
    }

  // At this point, we should have all the input models we want to include in the
  // final output weighting, as well as the relative weighting of the models.  

  wtot = 0;
  ngood = 0;
  for (n = comp[spectype].modstart; n < comp[spectype].modstop; n++)
    {
      if (good_models[n])
	{
	  wtot += weight[n];
	  ngood++;
	}
    }
  if (wtot == 0)
    {
      Error
	("model: Wtot must be greater than 0 or something is badly wrong\n");
      exit (0);
    }
  for (n = comp[spectype].modstart; n < comp[spectype].modstop; n++)
    {
      if (good_models[n])
	weight[n] /= wtot;
    }

// So now we know the absolute weighting.

  if (ngood == 0)
    {
      Error ("model: No models from %s survived pruning\n",
	     comp[spectype].name);
      exit (0);
    }
  else if (ngood == 1 && nmodel_error < 20)
    {
      Error
	("model: Only one model after pruning for parameters, consider larger model grid\n");
      for (n = comp[spectype].modstart; n < comp[spectype].modstop; n++)
	{
	  if (good_models[n])
	    {
	      Error ("model: %s %8.2f %8.2f\n", mods[n].name, par[0], par[1]);
	    }
	}
      nmodel_error++;
    }

  nwaves = comp[spectype].nwaves;

// Now create the spectrum 
  for (j = 0; j < nwaves; j++)
    {
      flux[j] = 0;
    }


  for (n = comp[spectype].modstart; n < comp[spectype].modstop; n++)
    {
      if (good_models[n])
	{
	  for (j = 0; j < nwaves; j++)
	    {
	      flux[j] += weight[n] * mods[n].f[j];
	    }
	}
    }



/* 081103 ksl -  Next section is new to deal with models with Temperatures less than in grid 
   140212 jm -- I've now included a fix for higher temperatures as well. See Pull request #67    

if our temperature is higher or lower than that in the grid then we want to scale it by 

scaling factor = B_nu ( T ) / B_nu (Tmax) = (exp(hnu / kTmax) - 1) / (exp(hnu / kT) - 1)

=> f_out =  1/(e**hnu/kT_out -1) /  1/(e**hnu/kT_in -1) * f_in = (1/q1)/(1/q2) * f_in= q2/q1 * f_in

Note that in the algorithm below we have to worry that the exp can become quite large, infinity
even, and so for those cases we want to make sure to calculate the ratio of qs directly
*/

  if (par[0] < comp[spectype].min[0] || par[0] > comp[spectype].max[0])   // is temp outside grid range
    {
      for (j = 0; j < nwaves; j++)              // cycle through wavelength bins
	{
	  lambda = comp[spectype].xmod.w[j] * 1.e-8;	// Convert lamda to cgs


    /* tscale is temperature to use in the BB function by which we need to scale the flux. 
       tscale can be larger or smaller than our actual temperature */
    if (par[0] < comp[spectype].min[0])
      tscale = comp[spectype].min[0];            // lowest temperature model

    else if (par[0] > comp[spectype].max[0])
      tscale = comp[spectype].max[0];            // highest temperature model


    /* calculate h*nu/kT for both temperatures */

	  q1 = H_OVER_K * C / (lambda * par[0]);	   //  h*nu/kT for model desired 

	  q2 = H_OVER_K * C  / (lambda * tscale);	   //  h*nu/kT for model that exists


	  /* q can be large- line below is attempt to keep exponents in range in that case */
	  if (q1 > 50. || q2 > 50.)
	    {
	      xxx = exp (q2 - q1);	// q2 - q1 should be negative since q1 is model desired
	    }
	  else
	    {
	      q2 = exp (q2) - 1.0;	// Model that exists has higher temperature
	      q1 = exp (q1) - 1.0;	// Model desired has lower T, hence q1 is larger
	      xxx = q2 / q1;
	    }

    /* multiply flux by scaling factor */
	  flux[j] *= xxx;

    /* nmodel_terror counts number of models where this is true. */
	  if ( nmodel_terror < 20 )
	    {
	      Error
		("model: Taking corrective action because parameter %f outside bound %f -> scaling factor %f \n",
		 par[0], tscale, xxx);

	      nmodel_terror++;
	    }
	}
    }



/* End of section to reweight the spectra. we can now copy fluxes to structure */

  for (j = 0; j < nwaves; j++)
    {
      comp[spectype].xmod.f[j] = flux[j];
    }
  for (j = 0; j < comp[spectype].npars; j++)
    {
      comp[spectype].xmod.par[j] = par[j];
    }

  return (nwaves);
}
