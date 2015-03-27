
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "atomic.h"
#include "python.h"
#include "balance_templates.h"

//OLD #define LINELENGTH 132
#define NP    10000

/* xsaha determines the abundances using the electron temperature and the Saha equation
   and then calculates the heating and cooling 

  02jun	ksl	Moved initializations of heading etc. after ionization calculation in
 		an attempt to get xsaha to behave more like python since some of ionization
	        calculations use the current value of the the total heating.	
	06nov	ksl	58b - Began work to get this to work again
 */

int
xsaha (one, p, nh, t_r, t_e, weight, ioniz_mode, freq_sampling)
     WindPtr one;
     PhotPtr p;
     double nh;
     double t_r;
     double t_e;
     double weight;
     int ioniz_mode;
     int freq_sampling;
{
  double total_fb ();
  double total_free ();
  double total_line_emission ();
  double total_emission ();
  double s;
  double alpha_recomb ();
  double planck (), kappa_ff ();
  double luminosity;
  double fb_cooling ();
  double vol;
//  double freqmin, freqmax;
  int n;
  PlasmaPtr xplasma;


  xplasma = &plasmamain[one->nplasma];

/* Now determine the abundances, concentrations, and partition functions based
on the conditions as currently defined.

ion_abundances is a python routine.  Currently (Sept 01) legal values are

0   on-the-spot approximation using existing t_e
1   LTE (e.g. concentrations)
2   Fixed concentrations
3   On the spot, with one_shot at updating t_e before calculating densities

It is far from clear that all of these will work at present with balance.  

01sept ksl

*/

  if (ioniz_mode < 0)
    {
      Log
	("xsaha: Calculating heating and cooling without changing abundances\n");
    }
  else
    {
      if (ioniz_mode == 0)
	{
	  Log
	    ("xsaha: Calculating abundances from On_the_spot or nebular approx\n");
	}
      else if (ioniz_mode == 1)
	{
	  Log ("xsaha: Calculating abundances from LTE\n");
	}
      else if (ioniz_mode == 2)
	{
	  Log ("xsaha: Fixed concentrations\n");
	}
      else if (ioniz_mode == 3)
	{
	  Log ("xsaha: Abundances using one shot (one_shot at t update\n");
	}
      else if (ioniz_mode == 4)
	{
	  Log ("xsaha: Test of detailed balance\n");
	}

      else
	{
	  Error ("xsaha: ioniz_mode unknown. Good luck!\n");
	}

      ion_abundances (xplasma, ioniz_mode);
    }


/* initialize some parts of the wind ptr (see wind_rad_init) */

  xplasma->j = xplasma->ave_freq = xplasma->lum = xplasma->heat_tot =
    xplasma->ntot = 0;
  xplasma->heat_ff = xplasma->heat_photo = xplasma->heat_lines = 0.0;
  xplasma->heat_z = 0;
  xplasma->lum_z = 0.0;
  xplasma->ntot = xplasma->nioniz = 0;
  for (n = 0; n < nions; n++)
    {
      xplasma->ioniz[n] = 0;
      xplasma->recomb[n] = 0;
      xplasma->heat_ion[n] = 0;
      xplasma->lum_ion[n] = 0;
    }

  if (xplasma->t_e != t_e)
    Error ("xsaha xplasma->t_e %f != t_e %f. Why?\n", xplasma->t_e, t_e);
  if (xplasma->t_r != t_r)
    Error ("xsaha xplasma->t_r %f != t_r %f. Why?\n", xplasma->t_r, t_r);
  if (xplasma->w != weight)
    Error ("xsaha xplasma->w %f != weight %f. Why?\n", xplasma->w, weight);

  vol = wmain[xplasma->nwind].vol;
  s = vol;
//  s = xplasma->vol;

/* Now calculate the heating  based on these abundances, etc.

There appear to be 3 types of heating, one in which the line is assumed
to have a simple flat topped profile, one in which the heating is calculated
using Python (i.e. there is no absorption until the photon scatters in a MC
sense), and one with a simple sobolev_line assumption (in which case one calculates
tau to define how much got absorbed.)

There appearted to be an error in the routine below which I fixed up, just so that
the choices actually matched what was said.  I fixed this.  01sept ksl

*/

/* Create photons which comprise a weighted BB distribution */

/* 0810 - ksl - with next call xband->nbands, is set to 1, and the 
 * frequency limites depned on t_r 
 */

  //0LD init_bands (t_r, 0.0, 1e50, 0, &xband);
  //OLD 120626 bands_init (t_r, 0.0, 1e50, 0, &xband);
  bands_init (t_r, &xband);

//???? New  make_bb (p, t_r, weight);
  xbb (p, t_r, weight, xband.f1[0], xband.f2[0], 1);

  if (geo.rt_mode == 2)
    {
      for (n = 0; n < NPHOT; n++)
	{
	  sobolev_line_heating (xplasma, &p[n], s);
	  radiation (&p[n], s);
	}
    }
  else if (geo.rt_mode == 1)
    {
      trans_phot (one, p, 0);
      xplasma[0].heat_tot += xplasma[0].heat_lines;
    }
  else
    {
      for (n = 0; n < NPHOT; n++)
	{
	  line_heating (xplasma, &p[n], s);
	  radiation (&p[n], s);
	}
    }

  luminosity = total_emission (one, 0.0, VERY_BIG);
  num_recomb (&xplasma[0], xplasma->t_e);

  summary (xplasma);

  return (0);
}


/***********************************************************
               Space Telescope Science Institute

Synopsis: 
	The routine cycle is intended to mimic what happens in Python in
	a very simple way, updating the ion abundances, etc. in a single
	cell, calculating the heating of the cell based on the updated
	abundances and then calling one shot to try to adjust the temperature.

   
Arguments:		

Returns:
 
 
Description:	

		
Notes:
	If t_r or t_e are changed then the routine will initialize everything
	to LTE, but otherwise it will go straight to the heating, and updating
	of abundances based on the new value of the heating.

	NB: t_r, t_e, weight are only used for generating the photon field. 
	The ionization calculation is based solely on what is contained in
	xplasma.


History:
	01dec	ksl	Updated to reflect new ways of calling routines
		    	like concentrations.	
	0810	ksl	67 - Relook and slight cleanup
 
**************************************************************/

double cy_tr = 0.0;
double cy_nh = 0.0;

int
cycle (xplasma, p, nh, t_r, t_e, weight, mode, freq_sampling, radmode)
     PlasmaPtr xplasma;
     PhotPtr p;
     double nh;
     double t_r;
     double t_e;
     double weight;
     int mode;
     int freq_sampling;
     int radmode;
{
  double total_fb ();
  double total_free ();
  double total_line_emission ();
  double total_emission ();
  double ne, s;
  double alpha_recomb ();
  double planck (), kappa_ff ();
  double luminosity;
  double agn_weight;		/*The weight of an agn photon */
  double vol;
  double agn_ip;		/* the ionization parameter */
  int i;
//  double spec[1000],fmin,fmax,dfreq;
//  int ipoint;
  double en1, en2, mean_freq, energy_density, sim_w_new;	/*Sim estimators */
  int n;
//  FILE *fopen(), *fptr;
  WindPtr one;

  one = &wmain[xplasma->nwind];
  printf ("Arriving in cycle, nh=%e, t_r =%e cy_tr=%e \n", nh, t_r, cy_tr);
/* initialize some parts of the wind ptr (see wind_rad_init) */



  xplasma->j = xplasma->ave_freq = xplasma->lum = xplasma->heat_tot =
    xplasma->ntot = xplasma->max_freq = 0;
  xplasma->heat_ff = xplasma->heat_photo = xplasma->heat_lines = 0.0;
  xplasma->heat_comp = xplasma->heat_ind_comp = 0.0;	//NSH 1208 73 added to zero compton heating
  xplasma->heat_z = 0;
  xplasma->lum_z = 0.0;
  xplasma->ntot = xplasma->nioniz = 0;
  for (n = 0; n < nions; n++)
    {
      xplasma->PWdenom[n] = 0.0;	//NSH 1208 73 added to zero pairwise denominator store
      xplasma->PWdtemp[n] = 0.0;	//NSH 1208 73 added to zero pariwise denom temp store
      xplasma->ioniz[n] = 0;
      xplasma->recomb[n] = 0;
      xplasma->heat_ion[n] = 0;
      xplasma->lum_ion[n] = 0;
    }

  if (xplasma->t_e != t_e)
    Error ("cycle: xsaha xplasma->t_e %f != t_e %f. Why?\n", xplasma->t_e,
	   t_e);
  if (xplasma->t_r != t_r)
    Error ("cycle: xsaha xplasma->t_r %f != t_r %f. Why?\n", xplasma->t_r,
	   t_r);
  if (xplasma->w != weight)
    Error ("cycle: xsaha xplasma->w %f != weight %f. Why?\n", xplasma->w,
	   weight);

  vol = wmain[xplasma->nwind].vol;
//  s = pow(vol,1./3.); //NSH 21/1/11 - changed from s=vol to s=pow(vol,1/3) - assumes plane illuminated cube
  s = vol;			// NSH 120802 Changed back, note, this means that the cell is 1cm2 area, with depth equal to volume/area = numerically equal to volume.

/* 
Now set up the initial concentrations.  This is to parallel what is
done in python where we initially begin with an LTE like assumption.
On successive calls to cycle this section will be skipped as a result
of the if statement
 */

  if (t_r != cy_tr || nh != cy_nh)
    {
      Log
	("Cycle: Fixing abundances to Saha values since new t_r %g or nh %g\n",
	 t_r, nh);

      if (mode != 5)		/*try to maintain original behiaviour */
	{
	  xplasma->t_e = 0.9 * t_e;	// Lucy guess
	  nebular_concentrations (xplasma, 2);
	}
      else
	{
	  printf
	    ("Were are going to try to use ioniz_mode 5 - going to NC\n");
	  nebular_concentrations (xplasma, 5);
	}
      ne = xplasma->ne;		//NSH 14/1/2011 - seemed necessary in order for the following line to produce a result.
      Log ("Cycle: On the spot estimate  t_r of %g gives ne %g\n", t_r, ne);
      cy_tr = t_r;
      cy_nh = nh;
      xplasma->dt_e = 0.0;	// Allow restart of oneshot limits
    }


/* Now calculate the heating */

/* First make some photons which comprise a BB distribution */

/* Set the frequency limits and then generate the photons */

  //OLD init_bands (t_r, 0.0, 1e50, 0, &xband);


  Error
    ("Cycle: Routines have been fixed to compile, but band_init is not properly sorted out\n");
  if (radmode == 1)
    {
//old1206  bands_init (t_r, &xband);
      geo.tstar = t_r;		//nsh 1208 put in to ensure we get a reasonable band even if we have a very high t_r.
      bands_init (7, &xband);
      xbb (p, t_r, weight, xband.f1[0], xband.f2[0], freq_sampling);
    }



  else if (radmode == 2)
    {
      agn_ip =
	emittance_pow (100 / HEV, 50000 / HEV, geo.lum_agn, geo.alpha_agn);
      printf
	("Calculating heating and cooling - I think the luminosity from 2-10 keV is %e Giving an ionisation parameter of %e\n",
	 geo.lum_agn, agn_ip / (geo.d_agn * geo.d_agn * nh));
//old1206  bands_init (t_r, 1.23e15,1.21e19,1,&xband);     /* at the moment we will have one big band */
      bands_init (1, &xband);	/* at the moment we will have one big band */
      agn_weight =
	emittance_pow (xband.f1[0], xband.f2[0], geo.lum_agn, geo.alpha_agn);

      agn_weight =
	(agn_weight) / (NPHOT * (4.0 * PI * geo.d_agn * geo.d_agn));

      /* This next line will generate the photons */
      photo_gen_agn (p, geo.r_agn, geo.alpha_agn, agn_weight, xband.f1[0],
		     xband.f2[0], -4, 0, NPHOT);
      printf ("Each pl photon has a weight equal to %e\n", agn_weight);
    }


  if (mode == 7)
    {
      freqs_init (xband.f1[0], xband.f2[xband.nbands - 1]);	//Set up estimators in case we want to do power law things.
      for (i = 0; i < geo.nxfreq; i++)
	{
	  xplasma->xj[i] = xplasma->xave_freq[i] = xplasma->xsd_freq[i] = xplasma->nxtot[i] = 0.0;	//zero estimators
	  //       xplasma->spec_mod_type[i]=-1; //Tell the code we dont have a model at the moment
	}
    }


  /* The lines below are used to generate an output file to show what photons are produced currently commented out because they cause a segmentaion fault very occasionally on my mac. The problem is caused by very high frequency photons generating an ipoint of 1000, quick fix would be to increase size of array, but it probably needs a more clever fix */
  /*
     fmin=xband.f1[0];
     fmax=xband.f2[0];
     dfreq=(fmax-fmin)/1001;

     printf("fmin=%e,fmax=%e,dfreq=%e\n",fmin,fmax,dfreq);

     for (n=0; n<999; n++)
     {
     spec[n]=0.0; 
     }
     for (n = 0; n < NPHOT; n++)
     {
     ipoint=(p[n].freq-fmin)/dfreq;
     spec[ipoint]=spec[ipoint]+p[n].w;  
     }



     fptr=fopen("bal.spec","w");
     for (n = 0; n < 998; n++)  
     {    
     fprintf(fptr,"freq= %e nphot= %e\n",(fmin+(n+1)*dfreq),spec[n]);  
     }  
     fclose(fptr);
   */
/*for (n=0; n<NPHOT; n++)
	{
	printf ("PHOT_DAT %e %e \n",p[n].freq,p[n].w);
	}*/

/* Shift values to old */
  xplasma->heat_tot_old = xplasma->heat_tot;
  xplasma->dt_e_old = xplasma->dt_e;

/* Set up sim estimators */
  en1 = 0.0;
  en2 = 0.0;
/* Next calculate the heating for this distribution */
  for (n = 0; n < NPHOT; n++)
    {
//      printf ("FLYING PHOTONS0 f=%e w=%e\n",p[n].freq,p[n].w);
/* sum sim estimators */
      en1 = en1 + p[n].w * s;
      en2 = en2 + p[n].w * s * p[n].freq;
      line_heating (xplasma, &p[n], s);
//      printf ("FLYING PHOTONS1 w=%e\n",p[n].w);
      radiation (&p[n], s);
//      printf ("FLYING PHOTONS2 w=%e\n",p[n].w);

    }
  xplasma->j /= (4. * PI * vol);	//Factor of 2 has been removed
  if (mode == 7)		//Generate PL estimators
    {
      for (i = 0; i < geo.nxfreq; i++)	/*loop over number of bands */
	{
	  if (xplasma->nxtot[i] > 0)	/*Check we actually have some photons in the cell in this band */
	    {
	      xplasma->xave_freq[i] /= xplasma->xj[i];	/*Normalise the average frequency */
	      xplasma->xsd_freq[i] /= xplasma->xj[i];	/*Normalise the mean square frequency */
	      xplasma->xsd_freq[i] = sqrt (xplasma->xsd_freq[i] - (xplasma->xave_freq[i] * xplasma->xave_freq[i]));	/*Compute standard deviation */
	      xplasma->xj[i] /= (4 * PI * vol);	/*Convert to radiation density */
	      printf
		("NSH Band %i - We have calculated xave_freq[i]=%e, xj[i]=%e xsd[i]=%e \n",
		 i, xplasma->xave_freq[i], xplasma->xj[i],
		 xplasma->xsd_freq[i]);
	    }
	  else
	    {
	      xplasma->xj[i] = 0;	/*If no photons, set both radiation estimators to zero */
	      xplasma->xave_freq[i] = 0;
	      xplasma->xsd_freq[i] = 0;
	    }
	}
    }





/*	if (radmode==2)
	{

	printf ("E1=%e,E2=%e",en1,en2);
	mean_freq=en2/en1;
	energy_density=en1/(4*PI*vol);
	printf ("Mean freq=%e\n",mean_freq);
	Error("Cycle: This may be an error\n");
	for (i=0;i<NXBANDS;i++){
		xplasma->pl_alpha[i]=sim_alphasolve(en2/en1,1.23e15,1.21e19);
		printf("I think alpha=%f\n",xplasma->pl_alpha[i]);



	sim_w_new=pl_w(en1,wmain[xplasma->nwind].vol,1,xplasma->pl_alpha[i],1.23e15,1.21e19);
	printf("sim W computed as %e compared to current weight used %e \n",sim_w_new,xplasma->pl_w[i]);

        xplasma->pl_w[i]=sim_w_new;
	}

	mean_freq=en2/en1;   

//	  trad = xplasma->t_r =
//	    H * mean_freq / (BOLTZMANN * 3.832);
//	  xplasma->w =
//	    PI * energy_density / (STEFAN_BOLTZMANN * trad * trad * trad *
//				    trad);	
	} NSH 120817 Commented all this PL stuff out - needs recoding and is just confusing things at the moment */


/* Now calculate the luminosity for these conditions */
  luminosity = total_emission (one, xband.f1[0], xband.f2[xband.nbands - 1]);
  num_recomb (&xplasma[0], xplasma->t_e);
  summary (xplasma);

/* Shift values to old */
  xplasma->dt_e = xplasma->t_e - xplasma->t_e_old;	//Must store this before others
  xplasma->t_e_old = xplasma->t_e;
  xplasma->t_r_old = xplasma->t_r;
  xplasma->lum_rad_old = xplasma->lum_rad;
  if (mode == 7)
    {
      spectral_estimators (xplasma);
    }
  printf ("About to go off to oneshot - mode = %i\n", mode);
  one_shot (xplasma, mode);




  Log ("Cycle: one_shot old t_r t_e %8.2g %8.2g, new t_r t_e %8.2g %8.2g\n",
       t_r, t_e, xplasma->t_r, xplasma->t_e);

//  Error ?? -- Next statement is probably superfluous, since calc_te, called by oneshot calculated total emission
//  luminosity = total_emission (xplasma, freqmin, freqmax);
  luminosity = total_emission (one, xband.f1[0], xband.f2[xband.nbands - 1]);
  num_recomb (&xplasma[0], xplasma->t_e);

  summary (xplasma);

/* Convergence check */
  convergence (xplasma);

  return (0);

}



/* 

The routine dumb_step increments te upward or downward depending on whether the total luminosity at
this temperature is less than or greater than the current value of the heating. It does not
recalculate the heating at any point.   The resulting te is incorporated into the windptr
	01dec	ksl	Updated to reflect new calls to routines like total_emission.  The routine
			seemed to work following these changes.
*/

int
dumb_step (one, te)
     WindPtr one;
     double *te;
{
  double lum_current, lum_hotter, lum_cooler;
  double tstep, hardness;
  double total_emission ();
  PlasmaPtr xplasma;

  xplasma = &plasmamain[one->nplasma];

  if (*te > 20000.)
    tstep = 0.1 * (*te);
  else
    tstep = 0.05 * (*te);

  xplasma->t_e = *te;
  lum_current = total_emission (one, 0.0, VERY_BIG);
  num_recomb (&xplasma[0], xplasma->t_e);

  Log
    ("Current values of luminosities compared to previously calculated heating\n");
  Log
    ("t_e %8.2e lum_tot  %8.2e lum_lines  %8.2e lum_ff  %8.2e lum_fb     %8.2e %8.2e %8.2e %8.2e %8.2e\n",
     xplasma->t_e, xplasma->lum_rad, xplasma->lum_lines, xplasma->lum_ff,
     xplasma->lum_fb, xplasma->lum_ion[0], xplasma->lum_ion[2],
     xplasma->lum_ion[3], xplasma->lum_z);
  Log
    ("t_r %8.2e heat_tot %8.2e heat_lines %8.2e heat_ff %8.2e heat_photo %8.2e %8.2e %8.2e %8.2e %8.2e\n",
     xplasma->t_r, xplasma->heat_tot, xplasma->heat_lines, xplasma->heat_ff,
     xplasma->heat_photo, xplasma->heat_ion[0], xplasma->heat_ion[2],
     xplasma->heat_ion[3], xplasma->heat_z);

  xplasma->t_e = *te + tstep;
  lum_hotter = total_emission (one, 0.0, VERY_BIG);
  num_recomb (&xplasma[0], xplasma->t_e);

  Log ("Results of raising t_e\n");
  Log
    ("t_e %8.2e lum_tot  %8.2e lum_lines  %8.2e lum_ff  %8.2e lum_fb     %8.2e %8.2e %8.2e %8.2e %8.2e\n",
     xplasma->t_e, xplasma->lum_rad, xplasma->lum_lines, xplasma->lum_ff,
     xplasma->lum_fb, xplasma->lum_ion[0], xplasma->lum_ion[2],
     xplasma->lum_ion[3], xplasma->lum_z);
  Log
    ("t_r %8.2e heat_tot %8.2e heat_lines %8.2e heat_ff %8.2e heat_photo %8.2e %8.2e %8.2e %8.2e %8.2e\n",
     xplasma->t_r, xplasma->heat_tot, xplasma->heat_lines, xplasma->heat_ff,
     xplasma->heat_photo, xplasma->heat_ion[0], xplasma->heat_ion[2],
     xplasma->heat_ion[3], xplasma->heat_z);

  xplasma->t_e = *te - tstep;
  lum_cooler = total_emission (one, 0.0, VERY_BIG);
  num_recomb (&xplasma[0], xplasma->t_e);

  Log ("Results of lowering t_e\n");
  Log
    ("t_e %8.2e lum_tot  %8.2e lum_lines  %8.2e lum_ff  %8.2e lum_fb     %8.2e %8.2e %8.2e %8.2e %8.2e\n",
     xplasma->t_e, xplasma->lum_rad, xplasma->lum_lines, xplasma->lum_ff,
     xplasma->lum_fb, xplasma->lum_ion[0], xplasma->lum_ion[2],
     xplasma->lum_ion[3], xplasma->lum_z);
  Log
    ("t_r %8.2e heat_tot %8.2e heat_lines %8.2e heat_ff %8.2e heat_photo %8.2e %8.2e %8.2e %8.2e %8.2e\n",
     xplasma->t_r, xplasma->heat_tot, xplasma->heat_lines, xplasma->heat_ff,
     xplasma->heat_photo, xplasma->heat_ion[0], xplasma->heat_ion[2],
     xplasma->heat_ion[3], xplasma->heat_z);

  Log ("te %8.2g lum  %8.2e\n", *te, lum_current);
  Log ("te %8.2g lum  %8.2e\n", *te + tstep, lum_hotter);
  Log ("te %8.2g lum  %8.2e\n", *te - tstep, lum_cooler);

  Log ("tr %8.2g heat %8.2e\n", xplasma->t_r, xplasma->heat_tot);

  hardness =
    (lum_current - xplasma->heat_tot) / (lum_current + xplasma->heat_tot);

  if (fabs (hardness) < 0.1)
    {
      Log ("Convergence is achieved!\n");
    }
  else if (hardness > 0.0)
    {				// luminosity is greater the heating 

      if (lum_cooler < lum_hotter)
	{
	  *te -= tstep;
	  Log ("Dropping t_e to %8.2g\n", *te);
	}
      else
	{
	  *te += tstep;
	  Log ("Raising t_e to % 8.2g\n", *te);
	}
    }
  else
    {				//luminosity is less than heating

      if (lum_cooler > lum_hotter)
	{
	  *te -= tstep;
	  Log ("Dropping t_e to %8.2g\n", *te);
	}
      else
	{
	  *te += tstep;
	  Log ("Raising t_e to % 8.2g\n", *te);
	}
    }

  xplasma->t_e = *te;
  return (0);
}

int
find_te (one)
     WindPtr one;
{
  double calc_te ();
  double luminosity, total_emission ();
  PlasmaPtr xplasma;

  xplasma = &plasmamain[one->nplasma];

  luminosity = total_emission (one, 0.0, VERY_BIG);
  num_recomb (&xplasma[0], xplasma->t_e);

  summary (xplasma);

  xplasma->t_e = calc_te (xplasma, TMIN, 1.2 * xplasma->t_r);

  luminosity = total_emission (one, 0.0, VERY_BIG);
  num_recomb (&xplasma[0], xplasma->t_e);

  summary (xplasma);

  return (0);
}
