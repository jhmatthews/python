/* 
	01dec	ksl	Updated for changes in calls to various routines, including nebular_concentrations,
			scattering_fraction.  Note that abso_two_level_atom was not updated
			and this likely is a problem in terms of accurate comparisons to python.
	06nov	ksl	58b -- began to try to get balance to work once abain
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "atomic.h"
#include "python.h"
#include "balance_templates.h"

//OLD #define LINELENGTH 132
#define NP    10000


int
absolute_saha (www, p, nh, t_r, t_e, weight, freq_sampling)
     PlasmaPtr www;
     PhotPtr p;
     double nh;
     double t_r;
     double t_e;
     double weight;
     int freq_sampling;
{
  double total_fb ();
  double total_free ();
  double absolute_total_line_emission ();
  double absolute_total_emission ();
  double xheat_lines, absolute_line_heating ();
  double ne, s;
  double alpha_recomb ();
  double planck (), kappa_ff ();
  double luminosity;
  double fb_cooling (), fb_h, fb_he1, fb_he2;
  double vol;
  int n;

/* OK make some photons which comprise a BB distribution */



  //OLD init_bands (t_r, 0.0, 1e50, 0, &xband);
  Error
    ("absolute_saha: Next line is likely an error in for inputs to bands_init\n");
  //OLD0620 bands_init (t_r, 0.0, 1e50, 0, &xband);
  bands_init (1, &xband);
  xbb (p, t_r, weight, xband.f1[0], xband.f2[0], freq_sampling);


/* Now set up the initial concentrations */
  nebular_concentrations (www, 1);


  Log ("Saha abundances FIXEd at %g gives ne %g\n", t_e, ne);
  for (n = 0; n < 2; n++)
    Log ("%2d %8.2g ", n, www->density[n]);
  Log ("\n");
  for (n = 2; n < 5; n++)
    Log ("%2d %8.2g ", n, www->density[n]);
  Log ("\n");
  for (n = 5; n < 10; n++)
    Log ("%2d %8.2g ", n, www->density[n]);
  Log ("\n");


/* Now calculate the heating */


/* initialize some parts of the wind ptr (see wind_rad_init) */

  www->j = www->ave_freq = www->lum = www->heat_tot = www->ntot = 0;
  www->heat_ff = www->heat_photo = www->heat_lines = 0.0;
  www->heat_z = 0;
  www->lum_z = 0.0;
  www->ntot = www->nioniz = 0;
  www->t_e = t_e;
  www->t_r = t_r;
  vol = wmain[www->nwind].vol;
  s = vol;
//OLD  s = www->vol;

  for (n = 0; n < nions; n++)
    {
      www->ioniz[n] = 0;
      www->recomb[n] = 0;
      www->heat_ion[n] = 0;
      www->lum_ion[n] = 0;
    }


/* Now calculate the ionization coefficients */
  for (n = 0; n < NPHOT; n++)
    {
// this call will certainly only work if p[n].grid==0
      xheat_lines = absolute_line_heating (www, &p[n], s);
      radiation (&p[n], s);
    }



  luminosity =
    absolute_total_emission (&wmain[www->nwind], t_e, xband.f1[0],
			     xband.f2[0]);
//OLD  luminosity = absolute_total_emission (www, xband.f1[0], xband.f2[0]);
//  luminosity = absolute_total_emission (www, t_e, freqmin, freqmax);

  summary (www);


//  fb_cooling (www->t_e, freqmin, freqmax, &fb_h, &fb_he1, &fb_he2);
  fb_cooling (www->t_e, xband.f1[0], xband.f2[0], &fb_h, &fb_he1, &fb_he2);

  Log
    ("Warning: UNCLEAR WHY fb_cooling and not new functions being used here?????\n");

  Log ("Detailed balance emissivity: h  %e he1 %e he2 %e\n", fb_h, fb_he1,
       fb_he2);
  Log ("Detailed balance cooling:    h  %e he1 %e he2 %e\n",
       www->ne * www->density[1] * fb_h, www->ne * www->density[3] * fb_he1,
       www->ne * www->density[4] * fb_he2);
  return (0);
}

//??? I have changed the calls to total_emission but not absolute_total_emission ksl 
double
absolute_total_emission (one, t_e, f1, f2)
     WindPtr one;		/* WindPtr to a specific cell in the wind */
     double t_e;		/* The electron temperature of the gas, which can be different from
				   the value stored in ww */
     double f1, f2;		/* The minimum and maximum frequency over which the emission is
				   integrated */
{
  double absolute_total_line_emission (), total_free (), total_fb ();
  double ffmax;
  int nplasma;
  PlasmaPtr xplasma;


  nplasma = one->nplasma;
  xplasma = &plasmamain[nplasma];


  if (f1 < 0.05 * BOLTZMANN * t_e / H)
    f1 = 0.05 * BOLTZMANN * t_e / H;
  ffmax = 20. * BOLTZMANN * t_e / H;
//      ffmax=10.*BOLTZMANN*t_e/H;  !! Increased to get fb emission at lower T
  // Why is it limited at all here ksl
  //      if(ffmax>f2)ffmax=f2;
  ffmax = f2;

  if (ffmax < f1)
    {
      xplasma->lum_rad = xplasma->lum_lines = xplasma->lum_ff =
	xplasma->lum_fb = 0;
    }
  else
    {
      xplasma->lum_rad = xplasma->lum_lines =
	absolute_total_line_emission (xplasma, t_e, f1, ffmax);
      xplasma->lum_rad += xplasma->lum_ff = total_free (one, t_e, f1, ffmax);
      xplasma->lum_rad += xplasma->lum_fb = total_fb (one, t_e, f1, ffmax);
    }

  return (xplasma->lum_rad);


}

double
absolute_total_line_emission (ww, t_e, f1, f2)
     PlasmaPtr ww;		/* WindPtr to a specific cell in the wind */
     double t_e;		/* The electron temperature of the gas, which can be different from
				   the value stored in ww */
     double f1, f2;		/* Minimum and maximum frequency */
{

  int n;
  double lum, x, z;
  double dd, d1, d2, sqrt_t;
  double f1_t, f2_t;
  double q21 (), scattering_fraction (), sf, absolute_two_level_atom ();
  double a, a21 ();
  char dummy[20];
  double weight_hold, te_hold, tr_hold;
  double vol;


  if (t_e <= 0 || f2 < f1)
    return (0);
  f1_t = 0.05 * BOLTZMANN * t_e / H;
  f2_t = 20.0 * BOLTZMANN * t_e / H;
  if (f1_t > f1)
    f1 = f1_t;
  if (f2_t < f2)
    f2 = f2_t;
  if (f2 <= f1)
    {
      Log ("Error f2<f1 in absolute_line heating\n");
      exit (0);
    }


  limit_lines (f1, f2);

  lum = 0;
  sqrt_t = sqrt (t_e);

  for (n = nline_min; n < nline_max; n++)
    {

      absolute_two_level_atom (lin_ptr[n], ww, &d1, &d2);

      a = a21 (lin_ptr[n]);

      dd = ww->density[lin_ptr[n]->nion];


// The lines below are to accommodate a new scattering_fraction call
// It is not obvious that all of these things need to be held in this
// way, but this should leave the operation of the code unchanged
      te_hold = ww->t_e;
      tr_hold = ww->t_r;
      weight_hold = ww->w;
      ww->t_e = t_e;
      ww->t_r = t_e;
      ww->w = 1;
      sf = scattering_fraction (lin_ptr[n], ww);
      ww->w = weight_hold;
      ww->t_e = te_hold;
      ww->t_r = tr_hold;

      vol = wmain[ww->nwind].vol;
      x = d2 * a * H * lin_ptr[n]->freq * vol * (1. - sf);
//OLD x = d2 * a * H * lin_ptr[n]->freq * ww->vol * (1. - sf);

      lum += lin_ptr[n]->pow = x;

      sprintf (dummy, "%e", x);
      if (dummy[0] == 'N')
	{
	  Log ("total_line_emission %s %e\n", dummy, z);
	}

//              x=2.*H*lin_ptr[n]->freq*lin_ptr[n]->freq*lin_ptr[n]->freq/(C*C);
      //              x/=(exp(H*lin_ptr[n]->freq/(BOLTZMANN*ww->t_e))-1.);
    }

  return (lum);

}

int ialh = 0;

/* Calculate line heating in the cell for one photon */
double
absolute_line_heating (w, p, ds)
     PlasmaPtr w;
     PhotPtr p;
     double ds;
{
  double f1, f2, phi, tau;
  double dd;
  double heating;
  double absolute_line_nsigma (), nsigma;
  double sf, scattering_fraction ();
  int n;
  double f1_t, f2_t;
  double x;

/* Carefully assure that we adopt the same frequency limits as in make_bb */

  f1 = 0.95 * p->freq;
  f2 = 1.05 * p->freq;
  f1_t = 0.05 * BOLTZMANN * w->t_e / H;
  f2_t = 20.0 * BOLTZMANN * w->t_e / H;
  if (f1_t > f1)
    f1 = f1_t;
  if (f2_t < f2)
    f2 = f2_t;
  if (f2 <= f1)
    {
      Log ("Error f2<f1 in absolute_line heating\n");
      exit (0);
    }
  else
    phi = 1. / (f2 - f1);

  limit_lines (f1, f2);

  nsigma = 0;
  for (n = nline_min; n < nline_max; n++)
    {
      dd = w->density[lin_ptr[n]->nion];
      sf = scattering_fraction (lin_ptr[n], w);

      nsigma += x = absolute_line_nsigma (lin_ptr[n], w) * (1. - sf);
    }
  tau = nsigma * ds * phi;
// Don't calculate exponents if tau is very small since you won't get a good result
  //      if(tau<1e-5) z=tau;
  //      else z=(1.-exp(-tau));
  //      heating=p->w*(z);
  //Try everything assuming we are optically thin
  heating = p->w * tau;
  w->heat_lines += heating;
  w->heat_tot += heating;
  return (w->heat_lines);
}

/* Calculate the total line absorption crossection for a specific transition
   allowing for stimulated emission */

double
absolute_line_nsigma (line_ptr, w)
     struct lines *line_ptr;
     PlasmaPtr w;
{
  double d1, d2, x;
  int ion;
  double f;
  double a21 (), absolute_two_level_atom ();

  ion = line_ptr->nion;
  f = line_ptr->freq;

  absolute_two_level_atom (line_ptr, w, &d1, &d2);

  x = a21 (line_ptr);
  x *= C * C / (8. * PI * f * f);
  x *= d1 * line_ptr->gu / line_ptr->gl - d2;
  return (x);
}

/* Calculate the ratio of populations of a two level atom
   in lte */

double
absolute_two_level_atom (line_ptr, www, d1, d2)
     struct lines *line_ptr;
     PlasmaPtr www;
     double *d1, *d2;
{
  double freq, te;
  int gu, gl;
  double dd;
  double z, q;

  gu = line_ptr->gu;
  gl = line_ptr->gl;
  freq = line_ptr->freq;

  te = www->t_e;
  dd = www->density[line_ptr->nion];

  z = H_OVER_K * freq / te;

  q = gu / gl * exp (-z);

  *d1 = dd / (1. + q);
  *d2 = dd - (*d1);

  return (q);
}





/***********************************************************
                                       Space Telescope Science Institute

 Synopsis: 
	double fb_cooling(t_e,f1,f2,fb_h,fb_he1,fb_he2)
Arguments:

Returns:
	
Description:	 

Notes:

History:
	98jul	ksl	Coded as part of python effort
	01oct	ksl	Moved to balance_abso.  It is unclear
			why this is needed!!!!

**************************************************************/

/* The FB_Constant is 4/sqrt(pi) * (m_e/2k)**1.5 h**4 / (m_e**3 c**2) */
#define FB_CONSTANT 3.839436e-62
/* Varibles needed for doing the integration */
double fb_alpha;
struct photoionization *fb_xptr;

double
fb_cooling (t_e, f1, f2, fb_h, fb_he1, fb_he2)
     double t_e, f1, f2;
     double *fb_h, *fb_he1, *fb_he2;
{
  int n;
  double h, he1, he2;
  double ff1, ff2;
  double zz;
  double qromb (), fb_cooling_d ();

  //printf("fb test t_e %e f1 %e f2 %e \n",t_e,f1,f2);

  fb_alpha = -H_OVER_K / t_e;

  (h) = (he1) = (he2) = 0.0;

  for (n = 0; n < nxphot; n++)
    {
      fb_xptr = &xphot[n];
      if (fb_xptr->nion == 0 && fb_xptr->freq_t < f2)
	{
	  ff1 = fb_xptr->freq_t;
	  if (f1 > ff1)
	    ff1 = f1;
	  printf ("h ratio f2/ft %e\n", f2 / fb_xptr->freq_t);
	  ff2 = f2;
	  if (ff2 > 20. * fb_xptr->freq_t)
	    ff2 = 20. * fb_xptr->freq_t;
	  h += qromb (fb_cooling_d, ff1, ff2, 1.e-4);
	}
      if (fb_xptr->nion == 2 && fb_xptr->freq_t < f2)
	{
	  ff1 = fb_xptr->freq_t;
	  if (f1 > ff1)
	    ff1 = f1;
	  he1 += qromb (fb_cooling_d, ff1, f2, 1.e-4);
	}
      if (fb_xptr->nion == 3 && fb_xptr->freq_t < f2)
	{
	  ff1 = fb_xptr->freq_t;
	  if (f1 > ff1)
	    ff1 = f1;
	  he2 += qromb (fb_cooling_d, ff1, f2, 1.e-4);
	}
    }

  zz = FB_CONSTANT * pow (t_e, -1.5);

  *fb_h = zz * ion[0].g / ion[1].g * (h);
  *fb_he1 = zz * ion[2].g / ion[3].g * (he1);
  *fb_he2 = zz * ion[3].g / ion[4].g * (he2);

//      printf("tester %d %d should be 2 and 1 \n",ion[3].g,ion[4].g);


  return (0);
}


double
fb_cooling_d (f)
     double f;
{
  double ft, x;
  double sigma_phot ();
  ft = fb_xptr->freq_t;
  x = (f - ft) * f * f * exp (fb_alpha * (f - ft)) * sigma_phot (fb_xptr, f);
  return (x);

}
