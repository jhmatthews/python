
/***********************************************************/
/** @file  get_models.c
 * @author ksl
 * @date   March, 2018
 *
 * @brief  The routines here provide a way to read grids of models and
 *   then to interpolate between them
 *
 * @details
 *    	There are two main routines
 *
 *    		get_models (modellist, npars, spectype)
 *
 *    	which reads a list of models
 *
 *    	and
 *
 *    		model (spectype, par)
 *
 *    	which given a set of parameters will interpolate between modles
 *    	given the set of models one wants to interpolate betwee, which
 *    	is determined by spectype, and the set of parameters (par) that
 *    	define the models
 *
 *  ### Notes ###
 *
 *  A detailed understanding of these routines requires an understandig
 *  of the associated data structures which are in models.h
 *
 *  The individual models are contained in a stucture called Models
 *
 *  A summary of a set of models is contained in ModSum, along with
 *  a interpolated model generated by the routine model
 *
 *  get_models is called as part of the initialization process.
 *
 *  model is called as one generates photons for specific raddiation
 *  sources, e.g the disk.  It stores a spectrum/cdf of that spectrum
 *  in ModSum, which can be used to genreate photons for that spectrum.
 *
 *
 ***********************************************************/




#include	<math.h>
#include	<stdio.h>
#include        <stdlib.h>
#include	<strings.h>
#include	<string.h>
#include 	"atomic.h"
#include	"python.h"      //This needs to come before modlel.h so that what is in models.h is used
#include         "models.h"
#define    	BIG 1e32


/* Get all the models of one type and regrid them onto the wavelength grid of the data */

/// Needed so can initialize nmods_tot the first time this routine is called
int get_models_init = 0;

/**********************************************************/
/**
 * @brief     This routine reads in a set of models for use in generation of spectra,
 *     e.g. a set of Kurucz stellar atmospheres which can
 *     used to simulate a disk, given a run of temperature and gravity
 *
 * @param [in] char  modellist[]   filename containing location and associated parameters of models
 * @param [in] int  npars   Number of parameters which vary for these models
 * @param [out] int * spectype   spectype an integer that is incremented each time a new set of models is read in
 * @return     *spectype is also returned
 *
 *
 * @details
 *
 * These are a generic set of routines to read a grid or grids
 * 	of models, e.g. a set of Kurucz stellar atmospheres which can
 * 	used to simulate a disk, given a run of temperature and gravity
 *
 * 	The grid can be any reasonable number of dimension
 *
 * ### Notes ###
 *
 * Basically, the way this works is as followsL
 *
 * modellist is a file that cantains a list of models, and the parameters that
 * are associated with each model. (There can be mulitple parameters associated
 * with each file. Typically these are temperature, and gravity, but these routines
 * are generic, and so one could image that the the spectra are organized quite differently
 *
 * Each model, that is read in contains a list of wavelengths and a flux, that should
 * be proportional to flambda.
 *
 * There are some assumptions that are made in terms of the order in the modellist file.
 * Specifically"
 *
 *
 * 	For 1d models we assume the models are read in increasing order of
 * 	the parameter of interest.
 *
 * 	2 or more parameter models can be read in any order.  Bilinear interpolation
 * 	is done for the 2 parameter models (as long as the grid is fully filled).
 *
 *
 * 	To deal with cases in which the grid is not fully filled, we follow
 * 	the same approach as for get_hub, i.e. for parameters t and g (i.e.
 * 	the first parameter being t, and the second g., we first search the grid
 * 	for all models with t just below (above) the desired value.  We then
 * 	interpolate on g for the model just below, and on g for the model just
 * 	above, and finally interpolate between these two.
 * 	The algorithm for interpolating between models
 *
 *
 * 	Note that each time a new set of models is read in, the spectype is
 * 	incremented.
 *
 *
 *
 **********************************************************/

int
get_models (modellist, npars, spectype)
     char modellist[];          // filename containing location and associated parameters of models
     int npars;                 // Number of parameters which vary for these models
     int *spectype;             //  The returned spectrum type


{
  FILE *mptr, *fopen ();
  char dummy[LINELENGTH];
  int n, m, mm, nxpar;
  double xpar[NPARS], xmin[NPARS], xmax[NPARS];
  int get_one_model ();
  int nw, nwaves;

  printf ("BLAH in get_models get_models_init=%i\n",get_models_init);
  nwaves = 0;

  if (get_models_init == 0)
  {
    nmods_tot = 0;
    ncomps = 0;                 // The number of different sets of models that have been read in
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
    Error ("get_models:Could not open file %s containing list of models \n", modellist);
    Exit (0);
  }

/* Now initialize the model summary structure */
  strcpy (comp[ncomps].name, modellist);
  comp[ncomps].modstart = nmods_tot;
  comp[ncomps].npars = npars;
  comp[ncomps].xcdf.limit1 = -99.;
  comp[ncomps].xcdf.limit2 = -99.;
  


/* Now get all the models of this type */
  n = nmods_tot;                // This is the starting point since may have read models in before
/* Establish initial limits on xmin and xmax so that they can be properly populated */
  for (m = 0; m < NPARS; m++)
  {
    xmin[m] = (BIG);
    xmax[m] = (-BIG);
    comp[ncomps].xmod.par[m] = -99;
  }
  


  nw = -1;                      // Initiallize nw
  while (n < NMODS && (fgets (dummy, LINELENGTH, mptr)) != NULL)
  {
    if (dummy[0] == '#' || dummy[0] == '!')
    {
    }                           //skip comment lines in models
    else
    {
      nxpar =
        sscanf (dummy, "%s %lf %lf %lf %lf %lf %lf %lf %lf %lf",
                mods[n].name, &xpar[0], &xpar[1], &xpar[2], &xpar[3], &xpar[4], &xpar[5], &xpar[6], &xpar[7], &xpar[8]);
      if (nxpar < npars)
      {
        Error ("get_models: nxpar (%d) < npars (%d) in line %s\n", nxpar, npars, dummy);
        Exit (0);
      }
      for (m = 0; m < npars; m++)
      {
        mods[n].par[m] = xpar[m];
        if (xpar[m] > xmax[m])
          xmax[m] = xpar[m];
        if (xpar[m] < xmin[m])
          xmin[m] = xpar[m];
      }
      for (mm = m; mm < NPARS; mm++)
        mods[n].par[mm] = -99;

      nwaves = get_one_model (mods[n].name, &mods[n]);
      if (nw > 0 && nwaves != nw)
      {
        Error ("get_models: file %s has %d wavelengths, others have %d\n", mods[n].name, nwaves, nw);
        Exit (0);
      }

      if ((n % 100) == 0)
        Log ("Model n %d %s\n", n, mods[n].name);
      n++;
    }
  }


  if (n == NMODS)
  {
    Error ("get_models: Reached maximum number of models %d. Please increase NMODS in .h file \n", n);
    Exit (0);
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
    Error ("get_models: No models from %s were read. Please check list of models!\n", comp[ncomps].name);
  }
  
  printf ("Got here ncomps=%i\n",ncomps);

  /* The next 3 lines set a normalization that is used by kslfit.  They are mostly
   * not relevant to python, where comp[ncomp[.min[0] refers to a normalization for
   * a model.  I've kept this initialization for now */
  comp[ncomps].min[0] = 0;
  comp[ncomps].max[0] = 1000;
  
  printf ("and here\n");

  for (m = 0; m < npars; m++)
  {
    comp[ncomps].min[m] = xmin[m];
    comp[ncomps].max[m] = xmax[m];
  }
  
  printf ("and lastly here ncomps=%i\n",ncomps);


  *spectype = ncomps;           // Set the spectype
  
  printf ("really????\n");

  ncomps++;

  
  return (*spectype);
}


/**********************************************************/
/**
 * @brief      reads a single model file from disk and puts the result into the structure
 *    onemod
 *
 * @param [in] char  filename[]   The name of the file containing the model
 * @param [in] struct Model *  onemod   A pointer to the structure where the model is to be stored
 * @return     The number of wavelengths and fluxes read from the file
 *
 * @details
 * This is simply a utility routine called by get_models to read in a single model into the modle stucture
 *
 *
 **********************************************************/

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
    Exit (0);
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
    Error ("get_one_model: model %s has more than %d wavelengths\n", filename, NWAVES);
  }


  fclose (ptr);
  return (n);
}

int nmodel_error = 0;
int nmodel_terror = 0;

/**********************************************************/
/**
 * @brief      interpolates between spectra in a grid  placing the results in
 *   	comp[spectype].xmod
 *
 * @param [in] int  spectype   A number which refers to a specific collecton of models, which were read in
 * @param [in] double  par[]   An array that contains the parameters for the models we want to interpoalted, e.g [T, grav]
 * @return     Returns nwaves in the the spectrum if a new spectrum was calculated
 * 	0 if it was the old spectrum.
 *
 * @details
 * This is an intepolation routine.  Having previously read in a set of models which are associated
 * with a set of parameter, e.g T and g, we want the spectrum at some t and g, which is not necessaily
 * in the grid.   This routien does the interpolation and stores the result  in comp[spectype].mod
 *
 * The routine is generic as long as the requested parameters are within the grid of models,
 * but the routine makes the assumption that the first parameter is a temperature if the
 * value of the first parameter outside the range of that parameter.
 *
 *
 * ### Notes ###
 *
 * If the parameters for the interpolated routine are within the grid of models, and the grid
 * is fully filled, then a straighfoward interpolation is carried out.
 *
 * If models are missing from the grid or if the requested model has parameters which are
 * outside the grid, an "interpolated" model is still returned but assumptions are made
 * which are intended to provide a reasonable solution for the case, where the first parameter is
 * a temperature, and the second parameter is a gravity.
 *
 * Specfically, for a missing model, the first thing that happens is that the routine
 * first ignores the missing model treating each of the others as if bi-linear interpolation was being used, and
 * then renormalizing to make sure the weights of the models that one does have
 * sum to 1.  The solution is "practical" because it always gives a result, and
 * because it gives the an appropriate answer when the grid is filled.  It is also
 * easily extensible to more than two dimensions.  But it should be noted, the
 * approach is also dangerous.
 *
 * Lastly, if the first parameter is out of range, it assumes that the first parameter was a
 * temperature, and it modifies the spectrum by the ratio of a BB speccturm for the desired
 * spectrum to the temperature of the max/min temperature in the grid.  So for example, suppose the
 * models which were read in were for Kuruce models, which run out a temperature of about 50,000k
 * and you wnat the spectrum of a star which is 100,000K. Then the routine will deliver a spectrum
 * which has line features of a 50,000K star, but an overall specturm which has the shape of a 100,000 K
 * BB. This is not a good situation for this extreme case, but may be reasonable if the temperature
 * differences are less extreme.
 *
 *
 **********************************************************/
int
model (spectype, par)
     int spectype;
     double par[];


{
  int j, n;
  int good_models[NMODS];       // Used to establish which models are to be included in creating output model
  double xmin[NPARS], xmax[NPARS];      // The vertices of a completely filled grid
  double weight[NMODS];         // The weights assigned to the models
  double hi, lo, delta, wtot;
  int ngood;
  double f;
  int nwaves;
  double flux[NWAVES];
  double q1, q2, lambda, tscale, xxx;   // Used for rescaleing according to a bb



/* First determine whether we already have interpolated this
exact model previously */

  n = 0;
  while (n < comp[spectype].npars && comp[spectype].xmod.par[n] == par[n])
  {
    n++;
  }
  if (n == comp[spectype].npars)
  {
    return (0);                 // This was the model stored in comp already
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
    /*   So at this point we know what xmin[j] and xmax[j] and we
       need to prune good_models
     */
    for (n = comp[spectype].modstart; n < comp[spectype].modstop; n++)
    {
      // Next lines excludes the models which are out of range.
      if (mods[n].par[j] > xmax[j] || mods[n].par[j] < xmin[j])
        good_models[n] = 0;
      /* Next line modifies the weight of this model assuming a regular grid
         If xmax==xmin, then par[j] was outside of the range of the models and
         so we need to weight the remaining models fully.
       */

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

/*If the weight given to a model is going to be zero, it needs to be
excluded from furthur consideration */
        if (weight[n] == 0.0)
          good_models[n] = 0;

      }
    }
  }

  /* At this point, we should have all the input models we want to include in the
     final output weighting, as well as the relative weighting of the models.
   */

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
    Error ("model: Wtot must be greater than 0 or something is badly wrong\n");
    Exit (0);
  }
  for (n = comp[spectype].modstart; n < comp[spectype].modstop; n++)
  {
    if (good_models[n])
      weight[n] /= wtot;
  }

// So now we know the absolute weighting.

  if (ngood == 0)
  {
    Error ("model: No models from %s survived pruning\n", comp[spectype].name);
    Exit (0);
  }
  else if (ngood == 1 && nmodel_error < 20)
  {
    Error ("model: Only one model after pruning for parameters, consider larger model grid\n");
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



/* Next section is new to deal with models with Temperatures less than in grid

if our temperature is higher or lower than that in the grid then we want to scale it by

scaling factor = B_nu ( T ) / B_nu (Tmax) = (exp(hnu / kTmax) - 1) / (exp(hnu / kT) - 1)

=> f_out =  1/(e**hnu/kT_out -1) /  1/(e**hnu/kT_in -1) * f_in = (1/q1)/(1/q2) * f_in= q2/q1 * f_in

Note that in the algorithm below we have to worry that the exp can become quite large, infinity
even, and so for those cases we want to make sure to calculate the ratio of qs directly
*/

  if (par[0] < comp[spectype].min[0] || par[0] > comp[spectype].max[0]) // is temp outside grid range
  {
    for (j = 0; j < nwaves; j++)        // cycle through wavelength bins
    {
      lambda = comp[spectype].xmod.w[j] * 1.e-8;        // Convert lamda to cgs


      /* tscale is temperature to use in the BB function by which we need to scale the flux.
         tscale can be larger or smaller than our actual temperature */
      if (par[0] < comp[spectype].min[0])
        tscale = comp[spectype].min[0]; // lowest temperature model

      else if (par[0] > comp[spectype].max[0])
        tscale = comp[spectype].max[0]; // highest temperature model


      /* calculate h*nu/kT for both temperatures */

      q1 = H_OVER_K * C / (lambda * par[0]);    //  h*nu/kT for model desired

      q2 = H_OVER_K * C / (lambda * tscale);    //  h*nu/kT for model that exists


      /* q can be large- line below is attempt to keep exponents in range in that case */
      if (q1 > 50. || q2 > 50.)
      {
        xxx = exp (q2 - q1);    // q2 - q1 should be negative since q1 is model desired
      }
      else
      {
        q2 = exp (q2) - 1.0;    // Model that exists has higher temperature
        q1 = exp (q1) - 1.0;    // Model desired has lower T, hence q1 is larger
        xxx = q2 / q1;
      }

      /* multiply flux by scaling factor */
      flux[j] *= xxx;
    }

    /* nmodel_terror counts number of models where this is true. */
    if (nmodel_terror < 20)
    {
      Error ("model: Rescaling spectra because parameter %f outside bound %f of spectra in grid\n", par[0], tscale);

      nmodel_terror++;
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
