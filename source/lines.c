
/***********************************************************
                                       Space Telescope Science Institute

 Synopsis:
	Subroutines in this file all have to do with resonance line radiation
	
	double	total_line_emission(ww,t_e,f1,f2) -- calculates the total luminosity due to line 
						emission in a cell
	double	q12(line_ptr,t)/ q21(line_ptr,t)  return the collisional excitation/ dexcitation coefficient
	double	a21(line_ptr)  returns the Einstein A coefficient for the line designated by line_ptr
	double	two_level_atom(line_ptr,ne,te,w,tr,dd,d1,d2)  calulates the population of ions in the 
				upper and lower state of the transition referred to by line_ptr
	double	line_nsigma(line_ptr,w)  calculates kappa_tot for a line in the two level approximation
	double	scattering_fraction(line_ptr,ne,te,dd,dvds,w,tr) calulates the fraction of ions which decay
				by radiation
	int	line_heat(www,pp,nres) calculates the amount heating which occurs due to a line 
				resonance
		
Arguments:		


Returns:
 
Description:	
	

Notes:

History:
 	97jun	ksl	Coding on py_wind began.
 	98feb	ksl	Coding of these subroutines began.
 	98apr25	ksl	Added checks in several routines to fix problems when the maximum freq
 				was less than the minimum frequency
 	98nov	ksl	Corrected various bugs associated with line heating, and modified to incorporate
			pdf's in line heating.  This change was desirable to speed the program when
			large numbers of lines are in the line list, because it cuts down the times
			the power of a line or set of lines is calculated..

 
**************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "atomic.h"
#include "python.h"
#include "recipes.h"


/* Calculate the total line emission from a wind cell

   History:
	98	ksl	Coded
	98aug29	ksl	Removed code which limits calculation to strong lines.  This calculation
			had the effect of eliminating all line emission at low density (10**8) under
			certain situations.  As a result, line emission and line absorption were not
			guaranteed to balance in conditions that were supposed to approach lte 
	98nov25	ksl	Corrected formulae to properly account for stimulated emission.
	99jan	ksl	Added limits to avoid calculation of extremely weak lines in
			attempt to speed up program.  Note that this is fairly dangerous
	06may	ksl	57+ -- Mods for change to  plasma.  Downsream programs need volume
 */

double
total_line_emission (one, f1, f2)
     WindPtr one;		/* WindPtr to a specific cell in the wind */
     double f1, f2;		/* Minimum and maximum frequency */
{

  double lum;
  double t_e;

  t_e = plasmamain[one->nplasma].t_e;

  if (t_e <= 0 || f2 < f1)
    return (0);

  limit_lines (f1, f2);

//  lum = lum_lines (ww, t_e, nline_min, nline_max);
  lum = lum_lines (one, nline_min, nline_max);

//Now populate the crude pdf for this wind element


  if (xxxpdfwind == 1)
    lum_pdf (&plasmamain[one->nplasma], lum);


  return (lum);

}

double
lum_lines (one, nmin, nmax)
     WindPtr one;		/* WindPtr to a specific cell in the wind */
     int nmin, nmax;		/* The min and max index in lptr array for which the power is to be calculated */
{
  int n;
  double lum, x, z;
  double dd, d1, d2;
  double q;
  double t_e;			/* The electron temperature of the gas, which can be different from
				   the value stored in ww */
  int nplasma;
  double foo1, foo2, foo3, foo4;
  PlasmaPtr xplasma;



  nplasma = one->nplasma;
  xplasma = &plasmamain[nplasma];
  t_e = xplasma->t_e;
  lum = 0;
  for (n = nmin; n < nmax; n++)
    {
      dd = xplasma->density[lin_ptr[n]->nion];

      if (dd > LDEN_MIN)
	{			/* potentially dangerous step to avoid lines with no power */
	  two_level_atom (lin_ptr[n], xplasma, &d1, &d2);
	  x = foo1 = lin_ptr[n]->gu / lin_ptr[n]->gl * d1 - d2;

	  z = exp (-H_OVER_K * lin_ptr[n]->freq / t_e);


//Next lines required if want to use escape probabilities               

	  q = 1. - scattering_fraction (lin_ptr[n], xplasma);

	  x *= foo2 = q * a21 (lin_ptr[n]) * z / (1. - z);

    /* JM 1411 -- corrected to use filled volume, rather than cell volume */
	  x *= foo3 = H * lin_ptr[n]->freq * xplasma->vol;
	  if (geo.line_mode == 3)
	    x *= foo4 = p_escape (lin_ptr[n], xplasma);	// Include effects of line trapping 
	  else
	    {
	      foo4 = 0.0;	// Added to prevent compilation warning
	    }

	  lum += lin_ptr[n]->pow = x;
	  if (x < 0)
	    {
	      Log
		("lum_lines: foo %10.3g (%10.3g %10.3g %10.3g) %10.3g %10.3g %10.3g %10.3g %10.3g %10.3g %10.3g\n",
		 foo1, d1, d2, dd, foo2, foo3, foo4, lin_ptr[n]->el,
		 xplasma->t_r, t_e, xplasma->w);
	      //mytrap ();  JM 1410 -- mytrap is deprecated
	    }
	  if (sane_check (x) != 0)
	    {
	      Error ("total_line_emission:sane_check %e %e\n", x, z);
	    }
	}
      else
	lin_ptr[n]->pow = 0;
    }


  return (lum);
}

/* This routine creates a luminosty pdf */
int
lum_pdf (xplasma, lumlines)
     PlasmaPtr xplasma;
     double lumlines;
{
  int n, m;
  double xsum, vsum;

  xplasma->pdf_x[0] = nline_min;
  xplasma->pdf_y[0] = 0;

  n = nline_min;
  vsum = 0.0;
  for (m = 1; m < LPDF; m++)
    {
      xsum = m * lumlines / (LPDF - 1);	/* This is the target */
      while ((vsum += lin_ptr[n]->pow) < xsum && n < nline_max)
	n++;
      n++;			// otherwise one will add lin_ptr[n]->pow twice
/* Why this is done this way is tricky.  The important point is that

xplasma->pdf_y[m]= sum lin_ptr[mm]->pow  where mm runs from pdf_x[m-1] to
pdf_x[m]-1

*/
      xplasma->pdf_x[m] = n;
      xplasma->pdf_y[m] = vsum;


    }

  return (0);
}



#define ECS_CONSTANT 4.773691e16	//(8*PI)/(sqrt(3) *nu_1Rydberg

/* 

   q21 calculates and returns the collisional de-excitation coefficient q21

   Notes:  This uses a simple approximation for the effective_collision_strength omega.
   ?? I am not sure where this came from at present and am indeed no sure that
   omega is omega(2->1) as it should be.  This should be carefully checked.??

   c21=n_e * q21
   q12 = g_2/g_1 q21 exp(-h nu / kT )

   History:
   98aug        ksl     Recoded from other routines so that one always calculates q21 in
			the same way.
	01nov	ksl	Add tracking mechanism so if called with same conditions it
			returns without recalculation
	12oct	nsh	Added, then commented out approximate gaunt factor given in
			hazy 2.

 */
struct lines *q21_line_ptr;
double q21_a, q21_t_old;

double
q21 (line_ptr, t)
     struct lines *line_ptr;
     double t;
{
  double gaunt;
  double omega;


  if (q21_line_ptr != line_ptr || t != q21_t_old)
    {

/*NSH 121024 - the followinglines implement the approximate gaunt factor as described in eq 4.21 in hazy 2*/
/* NSH 121026 commented out in py74a - not certain that this approximate gaunt factor actually inmproves anything */
      /*     if (line_ptr->istate == 1) //Neutral
         {
         gaunt = ((BOLTZMANN*t)/(H*line_ptr->freq))/10.0;
         }
         else
         {
         gaunt = 0.2;
         } */



      gaunt = 1;
      omega =
	ECS_CONSTANT * line_ptr->gl * gaunt * line_ptr->f / line_ptr->freq;
      q21_a = 8.629e-6 / (sqrt (t) * line_ptr->gu) * omega;
      q21_t_old = t;
    }

  return (q21_a);
}

double
q12 (line_ptr, t)
     struct lines *line_ptr;
     double t;
{
  double x;
  double q21 ();
  double exp ();

  x =
    line_ptr->gu / line_ptr->gl * q21 (line_ptr,
				       t) * exp (-H_OVER_K * line_ptr->freq /
						 t);

  return (x);
}


/* 
   a21 alculates and returns the Einstein A coefficient 
   History:
   98aug        ksl     Coded and debugged
   99jan        ksl Modified so would shortcircuit calculation if 
   called multiple times for same a
 */
#define A21_CONSTANT 7.429297e-22	// 8 * PI * PI * E * E / (MELEC * C * C * C)

struct lines *a21_line_ptr;
double a21_a;

double
a21 (line_ptr)
     struct lines *line_ptr;
{
  double freq;

  if (a21_line_ptr != line_ptr)
    {
      freq = line_ptr->freq;
      a21_a =
	A21_CONSTANT * line_ptr->gl / line_ptr->gu * freq * freq *
	line_ptr->f;
      a21_line_ptr = line_ptr;
    }

  return (a21_a);
}


/* 
   two_level_atom(line_ptr,ne,te,w,tr,dd,d1,d2) calculates the ratio
   n2/n1 and gives the individual densities for the states of a two level
   atom.  

   Description:
   In the two level approximation,

   n2 / n1 = ( c12 + g2/g1 c*c / (2 h nu * nu * nu )  * A21 J12) / 
   (c21 + A21 + c*c / (2 h nu * nu * nu )  * A21 J21)  

   In the on-the-spot approx. we assume,. 

   J21= W  (2 h nu * nu * nu )  / (exp (h nu / k T(rad)) -1)

   and this is what is calulated here

   History:
	98Aug	ksl	Coded and debugged
	99jan	ksl	Tried to improve speed by assuring that it does not
			calculate exactly the same atom twice in a row
	00nov	ksl	Adopted a very simple radiated domininated model for 
			transitions from excited state to excited state. 
			This is essentially a modified LTE approach.  Currently
			at least there is no collisional correction for either
			level, although it would be easy to add for the excited
			level.  The difficulty here is that one does not readily
			have all of the level information at hand.  
	01dec	ksl	Modified calling scheme and added partionion functions
	06may	ksl	57+ -- Modified to use plasma structue
	07mar	ksl	58c -- Tried to address problems associated with our 
			inconsistent treatment of level populations for two
			level atoms, This is at best a bandaide.
	14jul	nsh	78 -- changed to allow the use of a computed model for
			the mean intensity in a cell to calualate influence of radiation
			on the upper state population of a two level atom.
 */

struct lines *old_line_ptr;
double old_ne, old_te, old_w, old_tr, old_dd;
double old_d1, old_d2, old_n2_over_n1;

double
two_level_atom (line_ptr, xplasma, d1, d2)
     struct lines *line_ptr;
     PlasmaPtr xplasma;
     double *d1, *d2;
{
  double a, a21 ();
  double q, q21 (), c12, c21;
  double freq;
  double g2_over_g1;
  double n2_over_n1;
  double n1_over_ng;
  double n2_over_ng;
  double z;
  double exp ();
  int gg;
  double xw;
  double ne, te, w, tr, dd;
  int nion;
  double J; //Model of the specific intensity


  //Check and exit if this routine is called for a macro atom, since this should never happen

  if (line_ptr->macro_info == 1 && geo.rt_mode == 2 && geo.macro_simple == 0)
    {
      Error ("Calling two_level_atom for macro atom line. Abort.\n");
      exit (0);
    }

/* Move variables used in the calculation from the xplasma structure into subroutine variables */
  ne = xplasma->ne;
  te = xplasma->t_e;
  tr = xplasma->t_r;
  w = xplasma->w;
  nion = line_ptr->nion;
  dd = xplasma->density[nion];

  /* Calculate the number density of the lower level for the transition using the partition function */
  ;
  if (ion[nion].nlevels > 0)
    {
      dd *= config[ion[nion].firstlevel].g / xplasma->partition[nion];
    }

  if (old_line_ptr == line_ptr
      && old_ne == ne
      && old_te == te && old_w == w && old_tr == tr && old_dd == dd)
    {				// Then there is no need to recalculate eveything
      *d1 = old_d1;
      *d2 = old_d2;
    }
  else
    {
      if (line_ptr->el == 0.0)
	{			// Then the lower level is the ground state

/* For a ground state connected transition we correct for the partition
function in calculating the density of the lower level, and then we 
make an on-the-spot approximation for the upper level.  The only reason
this is a "improvement" over the numbers available from the partition
function directly is the allowance for collisions, and also for the
possibility that not all lines have upper levels that are included
in the configuration structure. 01dec ksl */

	  a = a21 (line_ptr);
	  q = q21 (line_ptr, te);
	  freq = line_ptr->freq;
	  g2_over_g1 = line_ptr->gu / line_ptr->gl;


	  c21 = ne * q;
	  c12 = c21 * g2_over_g1 * exp (-H_OVER_K * freq / te);

    /* NSH this if block removed to simplify code - it can be reinstated 
      if runtimes are a problem, but it will need an additional if statement
      to avoid missing it out if we are using a modelled specific intensity. */
//	  if (w < 1.e-6) 
//	    {			// Radiation is unimportant
//	      n2_over_n1 = c12 / (c21 + a);
//	    }
//	  else
//	    {			//Include effects of stimulated emission
//	      z = w / (exp (H_OVER_K * freq / tr) - 1.); //original

	      z= (C * C) / (2. * H * freq * freq * freq); //This is the factor which relates the A coefficient to the b coefficient

//	      n2_over_n1 = (c12 + g2_over_g1 * a * z) / (c21 + a * (1. + z)); //original

     /* we call mean intensity with mode 1 - this means we are happy to use the 
        dilute blackbody approximation even if we havent run enough spectral cycles 
        to have a model for J */
	   J = mean_intensity (xplasma, freq, 1);

     /* this equation is equivalent to equation 4.29 in NSH's thesis with the 
        einstein b coefficients replaced by a multiplied by suitable conversion 
        factors from the einstein relations. */
     n2_over_n1 = (c12 + g2_over_g1 * a * z * J) / (c21 + a*(1. + (J * z)));  


//	    }


	  *d1 = dd;
	  *d2 = *d1 * n2_over_n1;

	}
      else
	{			// The transition has both levels above the ground state

/* 
In the event that both levels are above the ground state, we assume
that the upper level population is given by an on-the-spot approximation.
We make the same assumption for the lower level, unless the lower level
is matastable in which case we set the weight to 1 and force equlibrium 
*/

	  gg = ion[line_ptr->nion].g;
	  z = w / (exp (line_ptr->eu / (BOLTZMANN * tr)) + w - 1.);
	  n2_over_ng = line_ptr->gu / gg * z;

/* For lower level, use an on the spot approximation if the lower level has a short radiative lifetive;
Othewise, assert that the lower level is metastable and set the radiative weight to 1 
ERROR -- At present I don't believe what one should do with metastable lower levels is
ERROR -- worked out, either in the program... where are we getting radiative rates
ERROR -- or conceptually
07mar - ksl - We still need to determine whether this makes sense at all !!
*/

	  xw = w;		// Assume all lower levels are allowed at present

	  z = xw / (exp (line_ptr->el / (BOLTZMANN * tr)) + xw - 1.);
	  n1_over_ng = line_ptr->gl / gg * z;

	  *d1 = dd * n1_over_ng;
	  *d2 = dd * n2_over_ng;
	  n2_over_n1 = n2_over_ng / n1_over_ng;
	}


      old_line_ptr = line_ptr;
      old_ne = ne;
      old_te = te;
      old_w = w;
      old_tr = tr;
      old_dd = dd;
      old_d1 = (*d1);
      old_d2 = (*d2);
      old_n2_over_n1 = n2_over_n1;
    }

  return (old_n2_over_n1);

}

/* Calculate the total line absorption crossection for a specific transition
   allowing for stimulated emission */

double
line_nsigma (line_ptr, xplasma)
     struct lines *line_ptr;
     PlasmaPtr xplasma;
{
  double d1, d2, x;

  two_level_atom (line_ptr, xplasma, &d1, &d2);	//xxxx

  x = (d1 - line_ptr->gl / line_ptr->gu * d2);
  x *= PI_E2_OVER_MC * line_ptr->f;
  return (x);
}





/***********************************************************
                                       Space Telescope Science Institute

scattering fraction(line_ptr,ne,te,dd,dvds,w,tr) calculate the fraction of excited 
state atoms which correspond to scattered photons, i.e. the portion which are
excited by radiation and return to the ground state via spontaneous emission.

Description:

In the radiative transfer equation for lines, we separate the fraction of photons
which are "absorbed" and the fraction which are "scattered".  

If line_mode==0, the atomosphere is a completely absorbing and no photons
		will be scattered.  In this mode, assuming the wind is a source
		of emission, the emissivity will be the Einstein A coefficient
If line_mode==1, the atmosphere is completely scattering, there will be no
		interchange of energy between the photons and the electrons
		as a result of radiation transfer
If line_mode==2, then a simple single scattering approximation is applied in which
		case the scattered flux is just  A21/(C21+A21*(1-exp(-h nu/k T_e). 
If line_mode==3, then radiation trapping is included as well.  The basic idea
		is to calculate the average number of scatters in a single
		interaction and there is heat lost in each of these scatters.
Notes: 
	It may be more efficient to combine several of these routines in view
	of the fact that the exp is calculated several times

History:
	98	ksl	Coded
	98aug	ksl	Rewritten to allow for stimulated decays in a
			two-level atom population.  Also rewrote the
			radiative trapping section.  It was clearly not
			correct previously! But is it correct now???
	98nov	ksl	Corrected scattering fraction to agree with Python notes.
			It was not correct previously 
	06may	ksl	57+ -- Modified to use plasma structure
*/


double
scattering_fraction (line_ptr, xplasma)
     struct lines *line_ptr;
     PlasmaPtr xplasma;
{
  double q, escape;
  double a, c, z;
  double sf;
  double ne, te;
  double w;			/* the radiative weight, and radiation tempeature */

  if (geo.line_mode == 0)
    return (0.0);		//purely absorbing atmosphere

  else if (geo.line_mode == 1)
    return (1.);		//purely scattering atmosphere

  //Populate variable from previous calling structure
  ne = xplasma->ne;
  te = xplasma->t_e;
  w = xplasma->w;

  c = (-H_OVER_K * line_ptr->freq / te);
  a = exp (c);
  z = 1.0 - a;
  a = a21 (line_ptr);
  c = ne * q21 (line_ptr, te) * z;
  q = c / (a + c);		//q == epsilon in Rybicki and elsewhere

  if (geo.line_mode == 2)
    return (1 - q);		//single scattering atmosphere

  else if (geo.line_mode == 3)
    {				// atmosphere with  line trapping 

      escape = p_escape (line_ptr, xplasma);
      //The following is exact
      sf = (1. - q) * escape / (q + escape * (1. - q));
      if (sane_check (sf))
	{
	  Error
	    ("scattering fraction:sane_check sf %8.2e q %8.2e escape %8.2e w %8.2e\n",
	     sf, q, escape, w);
	  //mytrap ();  JM 1410 -- mytrap is deprecated
	}
      return (sf);
    }

  else
    {				// Unknown treatment of line radiation

      Error ("scattering_fraction: Cannot handle %d line_mode\n",
	     geo.line_mode);
      exit (0);
    }

}


/* p_escape approximates the escape probability based on dvds

   History:
   98dec        ksl     Coded
   99jan        ksl     Added code so that it shortcircuits if 
   asked to claculate the same escape probability

	06may	ksl	57+ -- Modify for plasma structue
  1411 JM -- changed to use the sobolev function to calculate tau.
 */
struct lines *pe_line_ptr;
double pe_ne, pe_te, pe_dd, pe_dvds, pe_w, pe_tr;
double pe_escape;

double
p_escape (line_ptr, xplasma)
     struct lines *line_ptr;
     PlasmaPtr xplasma;
{
  double tau, two_level_atom ();
  double escape;
  double ne, te;
  double dd;			/* density of the relevent ion */
  double dvds;
  double w, tr;			/* the radiative weight, and radiation tempeature */
  WindPtr one;

//Populate variable from previous calling structure
  ne = xplasma->ne;
  te = xplasma->t_e;
  tr = xplasma->t_r;	//JM1308 in pre 76b versions this was incorrectly set to xplasma->t_e
  w = xplasma->w;
  dd = xplasma->density[line_ptr->nion];

  one = &wmain[xplasma->nwind];
  dvds = one->dvds_ave;

// Band-aid to prevent divide by zero in calculation of tau below
  if (dvds <= 0.0)
    {
      Error ("Warning: p_escape: dvds <=0 \n");
      return (0.0);
    }
  if (pe_line_ptr != line_ptr
      || pe_ne != ne
      || pe_te != te
      || pe_dd != dd || pe_dvds != dvds || pe_w != w || pe_tr != tr)
    {

      /* JM 1411 -- we used to have duplicated code here, but 
         now we call the sobolev function itself */
      tau = sobolev (one, one->x, dd, line_ptr, dvds);

      /* JM 1408 -- moved calculation of p_escape to subroutine below */
      escape = p_escape_from_tau (tau);


      pe_line_ptr = line_ptr;
      pe_ne = ne;
      pe_te = te;
      pe_dd = dd;
      pe_dvds = dvds;
      pe_w = w;
      pe_tr = tr;

      pe_escape = escape;
    }


  return (pe_escape);
}

/* p_escape_from_tau calculates the probability of escape
   given an actual tau. It simple returns the equation

   (1. - exp (-tau)) / tau;

   Except for high and low tau where it returns 1/tau and
   1.0 respectively. This is used by p_escape above, which 
   calculates the sobolev escape probability, and 
   also by the anisotropic scattering routines. 

   History:
   1408 JM  Moved here to avoid code duplication
 */


double
p_escape_from_tau(tau)
double tau;
{
  double escape;

  /* TAU_MIN is defined in python.h
     this is to stop numerical problems when tau is low */
  if (tau < TAU_MIN)
    escape = 1.;

  else if (tau < 10.0)
    escape = (1. - exp (-tau)) / tau;

  else
    escape = 1. / tau;

  return (escape);
}


/* line_heat calculates the amount of line heating that occurs after a resonance. It is called
   by trans_phot in python 

   History
   98sept       ksl     Coded
   98dec        ksl     Updated so that both heat_lines and heat_tot are included
	06my	ksl	57+ Updated for new structure approach
 */

int
line_heat (xplasma, pp, nres)
     PlasmaPtr xplasma;
     PhotPtr pp;
     int nres;
{
  double x, sf;


  check_plasma (xplasma, "line_heat");

  sf = scattering_fraction (lin_ptr[nres], xplasma);

  if (sane_check (sf))
    {
      Error ("line_heat:sane_check scattering fraction %g\n", sf);
    }
  x = pp->w * (1. - sf);
  xplasma->heat_lines += x;
  xplasma->heat_tot += x;

  // Reduce the weight of the photon bundle


  pp->w *= sf;

  return (0);

}
