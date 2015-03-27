/* These are subroutines of balance which are not related directly to ionization
calculations 
	01dec	ksl	Updated to reflect new calls to various routines, notably  
			scattering_fraction.
	06nov	ksl	58b -- Began to try to get this to work again
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "atomic.h"
#include "python.h"

//OLD #define LINELENGTH 132



int
summary (one)
     PlasmaPtr one;
{
  int n, first, last, nn, m;
  int nmax, n_tenth, n_hundredth;
  double pmax, cool_tot;
  FILE *fopen (), *heatcoolfile;

  heatcoolfile = fopen ("heat_cool.out", "a");

  /* Write out the element densities */
  for (nn = 0; nn < 5; nn++)
    {
      first = ele[nn].firstion;
      last = first + ele[nn].nions;
      Log ("%-5s ", ele[nn].name);
      for (m = first; m < last; m++)
	Log (" %8.2e", one->density[m]);
      Log ("\n");
    }
  cool_tot = one->lum_rad + one->lum_comp + one->lum_dr;	//dr and compton cooling do not generate photons
  Log
    ("t_e %8.2e cool_tot %8.2e lum_lines  %8.2e lum_ff  %8.2e lum_comp  %8.2e lum_dr %8.2e        lum_fb     %8.2e %8.2e %8.2e %8.2e %8.2e\n",
     one->t_e, cool_tot, one->lum_lines, one->lum_ff, one->lum_comp,
     one->lum_dr, one->lum_fb, one->lum_ion[0], one->lum_ion[2],
     one->lum_ion[3], one->lum_z);

  fprintf (heatcoolfile,
	   "t_e %8.2e cool_tot %8.2e lum_lines  %8.2e lum_ff  %8.2e lum_comp  %8.2e lum_dr %8.2e        lum_fb     %8.2e %8.2e %8.2e %8.2e %8.2e\n",
	   one->t_e, cool_tot, one->lum_lines, one->lum_ff, one->lum_comp,
	   one->lum_dr, one->lum_fb, one->lum_ion[0], one->lum_ion[2],
	   one->lum_ion[3], one->lum_z);

  Log
    ("t_r %8.2e heat_tot %8.2e heat_lines %8.2e heat_ff %8.2e heat_comp %8.2e heat_ind_comp %8.2e heat_photo %8.2e %8.2e %8.2e %8.2e %8.2e\n",
     one->t_r, one->heat_tot, one->heat_lines, one->heat_ff, one->heat_comp,
     one->heat_ind_comp, one->heat_photo, one->heat_ion[0], one->heat_ion[2],
     one->heat_ion[3], one->heat_z);

  fprintf (heatcoolfile,
	   "t_r %8.2e heat_tot %8.2e heat_lines %8.2e heat_ff %8.2e heat_comp %8.2e heat_ind_comp %8.2e heat_photo %8.2e %8.2e %8.2e %8.2e %8.2e\n",
	   one->t_r, one->heat_tot, one->heat_lines, one->heat_ff,
	   one->heat_comp, one->heat_ind_comp, one->heat_photo,
	   one->heat_ion[0], one->heat_ion[2], one->heat_ion[3], one->heat_z);
  fprintf (heatcoolfile, "\n");
  fclose (heatcoolfile);


  Log ("Heating/cooling: Tot %8.2g lines %8.2g \n",
       one->heat_tot / cool_tot, one->heat_lines / one->lum_lines);
  Log ("Number of ionizing photons in cell nioniz %d\n", one->nioniz);

  /* Write out the number of ionizations and recombinations */
  for (nn = 0; nn < 5; nn++)
    {
      first = ele[nn].firstion;
      last = first + ele[nn].nions;
      Log ("%-5s : Ioniz  ", ele[nn].name);
      for (m = first; m < last; m++)
	Log (" %8.2e", one->ioniz[m]);
      Log ("\n");
      Log ("%-5s : Recomb ", ele[nn].name);
      for (m = first; m < last; m++)
	Log (" %8.2e", one->recomb[m]);
      Log ("\n");
      Log ("%-5s : Heat   ", ele[nn].name);
      for (m = first; m < last; m++)
	Log (" %8.2e", one->heat_ion[m]);
      Log ("\n");
      Log ("%-5s : Cool   ", ele[nn].name);
      for (m = first; m < last; m++)
	Log (" %8.2e", one->lum_ion[m]);
      Log ("\n");
    }



  nmax = 0;
  pmax = 0;
  for (n = 0; n < nlines; n++)
    {
      if (lin_ptr[n]->pow > pmax)
	{
	  pmax = lin_ptr[n]->pow;
	  nmax = n;
	}
    }

  n_hundredth = n_tenth = 0;
  for (n = 0; n < nlines; n++)
    {
      if (lin_ptr[n]->pow > 0.1 * pmax)
	{
	  n_tenth++;
	  Log
	    ("Strong lines: lambda %8.2f pow %8.2g element %2d ion %2d density %8.2e\n",
	     C * 1e8 / lin_ptr[n]->freq, lin_ptr[n]->pow,
	     ion[lin_ptr[n]->nion].z, ion[lin_ptr[n]->nion].istate,
	     one->density[lin_ptr[n]->nion]);
	}
      if (lin_ptr[n]->pow > 0.01 * pmax)
	{
	  n_hundredth++;
	}
    }

  Log
    ("Line sum: pmax %8.2e lambda_max %8.2f frac_max %8.2g n_tenth %d n_hundreth %d\n",
     pmax, C * 1e8 / lin_ptr[nmax]->freq, pmax / one->lum_lines, n_tenth,
     n_hundredth);



  return (0);
}


/* 
Calculate line heating in the cell for one photon.  This routine assumes that the
line is described by a square profile

ksl 01 aug
*/

double
line_heating (w, p, ds)
     PlasmaPtr w;
     PhotPtr p;
     double ds;
{
  double f1, f2, phi, tau, dvds;
  double dd, heating;
  double line_nsigma (), scattering_fraction (), sf, nsigma;
  int n;

  phi = 1. / (0.1 * p->freq);
  f1 = 0.95 * p->freq;
  f2 = 1.05 * p->freq;

  dvds = 1.0;

  if (limit_lines (f1, f2) == 0)
    {
      return (w->heat_lines);	// There were no lines of interest

    }

  nsigma = 0;
  for (n = nline_min; n <= nline_max; n++)
    {
      dd = w->density[lin_ptr[n]->nion];
      sf = scattering_fraction (lin_ptr[n], w);

      nsigma += line_nsigma (lin_ptr[n], w) * (1. - sf);
    }
  tau = nsigma * ds * phi;
  heating = p->w * tau;
  w->heat_lines += heating;
  w->heat_tot += heating;
  return (w->heat_lines);
}

/* 
Calculate line heating in the cell for one photon assuming the Sobolev approximation
Note that this avoids the MC aspects of python in which heating only occurs if the
photon is absorbed.

ksl 01sep

*/
double
sobolev_line_heating (w, p, ds)
     PlasmaPtr w;
     PhotPtr p;
     double ds;
{
  double f1, f2, tau, dvds;
  double dd, heating;
  double line_nsigma (), scattering_fraction (), sf;
  int n;

  f1 = p->freq;
  f2 = p->freq * (1 + geo.pl_vmax / C);

  dvds = 1.0;

  if (limit_lines (f1, f2) == 0)
    {
      return (w->heat_lines);	//there were no interesting lines

    }
//      limit_lines(0,VERY_BIG);

  tau = 0;
  for (n = nline_min; n <= nline_max; n++)
    {
      if (f1 < lin_ptr[n]->freq && lin_ptr[n]->freq < f2)
	{
	  dd = w->density[lin_ptr[n]->nion];
	  sf = scattering_fraction (lin_ptr[n], w);
	  tau +=
	    line_nsigma (lin_ptr[n],
			 w) * (1. -
			       sf) * C / (lin_ptr[n]->freq * geo.pl_vmax /
					  geo.pl_vol);
	}
    }
  heating = p->w * tau;
  w->heat_lines += heating;
  w->heat_tot += heating;
  return (w->heat_lines);
}
