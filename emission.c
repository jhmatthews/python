/***********************************************************
                                       Space Telescope Science Institute

 Synopsis:
	Subroutines in this file all have to do with radiation from the wind
Arguments:		


Returns:
 
Description:	
	

Notes:

	Comment: Some but not all of the calls to routines were changed here
	when the Plasma structure was introduced.  A general cleanup would
	make this consistent.  ksl 110805

History:
 	97jun	ksl	Coding on py_wind began.
 	98feb	ksl	Coding of these subroutines began.
 	98apr25	ksl	Added checks in several routines to fix problems when the maximum freq
 				was less than the minimum frequency
	01oct	ksl	Recombination related routines removed to recomb.c while
			topbase modifications made.
	01nov	ksl	Implement pdf functionality for line luminosity spped up
        11Aug   nsh     Changes to total_emission and wind_luminosity to allow compton heating to be computed
	11Oct	nsh	Compton heating removed from total_emission, DR cooling added to wind_luminosity
 
**************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "atomic.h"
#include "python.h"
//#include "recipes.h"




/***********************************************************
                                       Space Telescope Science Institute

Synopsis:  wind_luminosity (f1, f2) calculate the luminosity of the entire 
wind between freqencies f1 and f2 

Arguments:		


Returns:
 
Description:	
	

Notes:

History:
	06may	ksl	57+ -- Changed to accommodate plasma structure.  
			Note that the call for wind_luminostiy has changed
			to remove the pointer call.  Ultimately rewrite
			using plasma structure alone.
	11aug	nsh	70 Modifications made to incorporate compton cooling
        11sep   nsh     70 Modifications in incorporate DR cooling (very approximate at the moment)
	12sep	nsh	73 Added a counter for adiabatic luminosity (!)
 	13jul	nsh	76 Split up adiabatic luminosity into heating and cooling.
 
**************************************************************/

double
wind_luminosity (f1, f2)
     double f1, f2;		/* freqmin and freqmax */
{
  double lum, lum_lines, lum_fb, lum_ff, lum_comp, lum_dr, lum_adiab, heat_adiab;	//1108 NSH Added a new variable for compton cooling
//1109 NSH Added a new variable for dielectronic cooling
//1307 NSH Added a new variable to split out negtive adiabatic cooling (i.e. heating).
  int n;
  double x;
  int nplasma;


  lum = lum_lines = lum_fb = lum_ff = lum_comp = lum_dr = lum_adiab = heat_adiab = 0;	//1108 NSH Zero the new counter 1109 including DR counter
  for (n = 0; n < NDIM2; n++)
    {
      if (wmain[n].vol > 0.0)
	{
	  nplasma = wmain[n].nplasma;
	  lum += x = total_emission (&wmain[n], f1, f2);
	  lum_lines += plasmamain[nplasma].lum_lines;
	  lum_fb += plasmamain[nplasma].lum_fb;
	  lum_ff += plasmamain[nplasma].lum_ff;
	  lum_comp += plasmamain[nplasma].lum_comp;	//1108 NSH Increment the new counter by the compton luminosity for that cell.
	  lum_dr += plasmamain[nplasma].lum_dr;	//1109 NSH Increment the new counter by the DR luminosity for the cell.
	  if (geo.adiabatic) //130722 NSH - slight change to allow for adiabatic heating effect - now logged in a new global variable for reporting.
		{
		if (plasmamain[nplasma].lum_adiabatic >= 0.0)
			{
			lum_adiab += plasmamain[nplasma].lum_adiabatic;
			}
		else
			{
			heat_adiab += plasmamain[nplasma].lum_adiabatic;
			}
		}	
	  else
		{
		lum_adiab = 0.0;
		}
	  if (x < 0)
	    mytrap ();
	  if (recipes_error != 0)
	    Error ("wind_luminosity: Received recipes error on cell %d\n", n);
	}
    }

/* Deciding which of these to fill is tricky because geo contains
   values for geo.lum_wind and geo.f_wind separately.  Could be cleaner.
   ksl 98mar6 */

//      geo.lum_wind=lum;
  geo.lum_lines = lum_lines;
  geo.lum_fb = lum_fb;
  geo.lum_ff = lum_ff;
  geo.lum_comp = lum_comp;	//1108 NSH The total compton luminosity of the wind is stored in the geo structure
  geo.lum_dr = lum_dr;		//1109 NSH the total DR luminosity of the wind is stored in the geo structure
  geo.lum_adiabatic = lum_adiab;
  geo.heat_adiabatic = heat_adiab;

  return (lum);
}



/***********************************************************
             Space Telescope Science Institute

Synopsis:  total_emission (one, f1, f2) Calculate the total emission of a single cell 

Arguments:		


Returns:
 
Description:	
	

Notes:
	Total emission gives the total enery loss due to photons.  It does
	not include other coooling sources, e. g. adiabatic expansion.

	It returns the total luminosity, but also stores the luminosity due
	to various types of emssion, e.g ff, fb, lines, compton into the
	Plasms cells

	Comment: The call to this routine was changed when PlasmaPtrs
	were introduced, but it appears that the various routines 
	that were called were not changed.

History:
	97	ksl	Coded
	98oct	ksl     Removed temperature limits within total_emission in
			and attempt to concentrate the frequency limits elsewhere.
	01nov	ksl	Removed superflous code
	01nov	ksl	Changed call so t_e passed through wind structure.  Other
			routines should be changed ALSO, BUT NOT NOW
	04Mar   SS      Modified to divert if Macro Atoms being used.
	04Apr   SS      Modified to include a collisional cooling term
	                - this is getting very messy and confusing - I think
                        I should probably put together a completely separte
			set of routines for doing the heating and cooling in the
			macro atom cases. (SS) since this is now not really a
			total emission but total cooling.
	05may	ksl	57+ -- To use plasma structure for most things.  If vol
			is added to plasma then one can change the call
	11aug	nsh	70 changes made to allow compton heating to be calculated.
 
 
**************************************************************/


double
total_emission (one, f1, f2)
     WindPtr one;		/* WindPtr to a specific cell in the wind */
     double f1, f2;		/* The minimum and maximum frequency over which the emission is
				   integrated */
{
  double total_line_emission (), total_free (), total_fb ();
  double total_fb_matoms (), total_bb_cooling ();
  double t_e;
  int nplasma;
  PlasmaPtr xplasma;

  nplasma = one->nplasma;
  xplasma = &plasmamain[nplasma];

  t_e = xplasma->t_e;		// Change so calls to total emission are simpler

  if (f2 < f1)
    {
      xplasma->lum_rad = xplasma->lum_lines = xplasma->lum_ff = xplasma->lum_fb = 0;	//NSH 1108 Zero the new lum_comp variable NSH 1101 - removed
    }
  else
    {
      if (geo.rt_mode == 2)	//Switch for macro atoms (SS)
	{
	  xplasma->lum_fb = total_fb_matoms (xplasma, t_e, f1, f2) +
	    total_fb (one, t_e, f1, f2);
	  //The first term here is the fb cooling due to macro ions and the second gives
	  //the fb cooling due to simple ions.
	  //total_fb has been modified to exclude recombinations treated using macro atoms.
	  xplasma->lum_rad = xplasma->lum_fb;
	  //Note: This the fb_matom call makes no use of f1 or f2. They are passed for
	  //now in case they should be used in the future. But they could
	  //also be removed.
	  // (SS)
	  xplasma->lum_lines = total_bb_cooling (xplasma, t_e);
	  xplasma->lum_rad += xplasma->lum_lines;
	  /* total_bb_cooling gives the total cooling rate due to bb transisions whether they
	     are macro atoms or simple ions. */
	  xplasma->lum_ff = total_free (one, t_e, f1, f2);
	  xplasma->lum_rad += xplasma->lum_ff;


	}
      else			//default (non-macro atoms) (SS)
	{
	  xplasma->lum_rad = xplasma->lum_lines =
	    total_line_emission (one, f1, f2);
	  xplasma->lum_rad += xplasma->lum_ff = total_free (one, t_e, f1, f2);
	  xplasma->lum_rad += xplasma->lum_fb = total_fb (one, t_e, f1, f2);


	}
      /* NSH 1108 - This line calls the function total_comp which returns the compton luminosity for the cell
         this is then added to lum_rad, the total luminosity of the cell */
      /* NSH 1110 - Ths line has now been commented out. It was adding the compton luminosity to the lum_rad 
         variable. This was then generating line photons to fill the compton luminosity. We need to do something
         better, but at the moment, simply removing this line, and putting the calculation of compton luminosity 
         into calc_te with the adiabatic cooling and the new DR cooling is the way to make things a little more stable */
      //OLD     xplasma->lum_rad += xplasma->lum_comp = total_comp (one, t_e); 
    }


  return (xplasma->lum_rad);


}



/***********************************************************
                                       Space Telescope Science Institute

Synopsis:  
	adiabatic_cooling (one, t) determines the amount of 
	adiabatic cooling in a cell, in units of luminosity.

Arguments:		
	WindPtr one;	pointer to wind cell
	double t;		electron temperature

Returns:
 
Description:	
   Adiabatic cooling is clearly the amount of PdV work done by
   a fluid element, per unit time dt. Thus it is equal to P dV/dt.  
   The only real question is whether dV/dt is given by the volume * div v.
   div v here is the divergence of the velocity field.
	
Notes:
  JM 1401 -- I've rederived the expression for dv/dT. It follows
        directly from the continuity equation and is indeed equal 
        to volume * div_v. 

        Note also that this function should only be called
        if geo.adiabatic == 1, in which case it populates
        xplasma->lum_adiabatic. This is used in heating and cooling
        balance. We also use it as a potential destruction choice for 
        kpkts in which case the kpkt is thrown away by setting its istat 
        to P_ADIABATIC.


History:
	04nov	ksl	Stuart had switched adiabatic cooling off
  			by commenting out the entire program.  I have
  			reinstated the calculation, but created a new
  			variable geo.adiabatic to determine whether this
  			routine is called in python.  At present, even
  			if the routine is called it has no effect.  The
  			or more properly, one issue is how to meld adiabatic 
  			cooling with the macro atom approach.
	06may	ksl	57+ -- Adapted to include plasma structure
	11aug	ksl	70 - Adiabatic cooling returns the cooling but
			does not store it.
 
 
**************************************************************/


double
adiabatic_cooling (one, t)
     WindPtr one;
     double t;
{
  double cooling;
  int nplasma;
  PlasmaPtr xplasma;

  nplasma = one->nplasma;
  xplasma = &plasmamain[nplasma];

  //JM 1401 -- here was an old factor of 3/2 which KSL and JM believe to be incorrect. 
  cooling = xplasma->ne * BOLTZMANN * t * one->vol * one->div_v;

  return (cooling);
}


/***********************************************************
                                       Space Telescope Science Institute

Synopsis: photo_gen_wind(p,w,weight,freqmin,freqmax,istart) generates 
	photons within the wind and stores them in the photon array

Arguments:		
	WindPtr w;
	PhotPtr p;
	double weight;				weight of each photon
	double freqmin,freqmax;          frequency interval within which photons are created
	int istart;					photon number in which first photon is placed

Returns:
	Normally returns the number of photons created.  So photons will exist from p[istart]
	to p[istart+return value].
	
	If photo_gen_wind tries to create more photons than exist in the photon structure the
	program will stop (rather than continue incorrectly or start blasting away memory.)
 
Description:	
	Calls photo_gen_wind_one for each cell.

Notes:

History:
 	98feb	ksl	Coding began
 	98may16	ksl	Added line to keep track of the number of recombinations
 	98dec21	ksl	Added the variable inwind and clarified loop which finds a position in
 				the wind in the grid cell
	99oct19	ksl	Fixed error in loop above.
	06may	ksl	57+ -- Began to modify to use plasma, but this is a
			stopgap for now since should not need to be addressing
			wind array.  Changed call to remove wind since entire
			grid was tramsmitted.
 
**************************************************************/

int
photo_gen_wind (p, weight, freqmin, freqmax, photstart, nphot)
     PhotPtr p;
     double weight;
     double freqmin, freqmax;
     int photstart, nphot;
{
  int n;
  int photstop;
  double xlum, xlumsum, lum;
  double v[3];
  int icell;
  int nplasma;
  int nnscat;


  photstop = photstart + nphot;
  Log_silent ("photo_gen_wind creates nphot %5d photons from %5d to %5d \n", nphot,
       photstart, photstop);

  for (n = photstart; n < photstop; n++)
    {
      /* locate the wind_cell in which the photon bundle originates.
         Note: In photo_gen, both geo.f_wind and geo.lum_wind will have been determined.
         geo.f_wind refers to the specific flux between freqmin and freqmax.  Note that
         we make sure that xlum is not == 0 or to geo.f_wind. */
      xlum = (rand () + 0.5) / (MAXRAND) * geo.f_wind;

      xlumsum = 0;
      icell = 0;
      while (xlumsum < xlum)
	{
//?? 57+ This is not the best way to do this.  We should be able to sum over the plasma structure
	  if (wmain[icell].vol > 0.0)	//only consider cells with volume greater than zero
	    {
	      nplasma = wmain[icell].nplasma;	//get the plasma cell in this wind cell
	      xlumsum += plasmamain[nplasma].lum_rad;	/*increment the xlumsum by the lum_rad in this plasma cell - note that due to the way wind_luminosity gets called, this is actually the band limited flux not the luminosity. */
	    }
	  icell++;		/*If we are not yet up to xlum, go to the next cell */
	}
      icell--;			/* We have got up to xlum, so this is the cell in which the photon must be generated */

      nplasma = wmain[icell].nplasma;


      /* Now generate a single photon in this cell */

      /*Get the total luminosity and MORE IMPORTANT populate xcol.pow and other parameters */
      lum = plasmamain[nplasma].lum_rad;	/* Whilst this says lum - I'm (nsh) pretty sure this is actually a flux between two frequency limits) */

      xlum = lum * (rand () + 0.5) / (MAXRAND);	/*this makes a small test luminosity */

      xlumsum = 0;

      p[n].nres = -1;
      p[n].nnscat = 1;
      if ((xlumsum += plasmamain[nplasma].lum_ff) > xlum)	/*Add the free free luminosity of the cell to the running total. If it is more than our small test luminosity, then we need to make some ff photons */
	{
	  p[n].freq = one_ff (&wmain[icell], freqmin, freqmax);	/*Get the frequency of one ff photon */
	  if (p[n].freq <= 0.0)
	    {
	      Error_silent
		("photo_gen_wind: On return from one_ff: icell %d vol %g t_e %g\n",
		 icell, wmain[icell].vol, plasmamain[nplasma].t_e);
	      p[n].freq = 0.0;
	    }
	}
      else if ((xlumsum += plasmamain[nplasma].lum_fb) > xlum)	/*Do the same for fb */
	{
	  p[n].freq = one_fb (&wmain[icell], freqmin, freqmax);
	}
      else
	{
	  p[n].freq = one_line (&wmain[icell], freqmin, freqmax, &p[n].nres);	/*And fill all the rest of the luminosity up with line photons */
	}
      p[n].w = weight;
      /* Determine the position of the photon in the moving frame */

      /* !! ERROR - Need to account for emission from torus if it exists */

      if (wmain[icell].inwind > 1)
	{
	  get_random_location (icell, 2, p[n].x);	/* NSH 1110 Added this if statement to take account of photons being generated from the torus. Hope I've done it correctly!! */
	}
      else
	{
	  get_random_location (icell, 0, p[n].x);
	}

      p[n].grid = icell;





	  // Determine the direction of the photon
	  // ?? Need to allow for anisotropic emission here
      // JM 1406 -- I think there's a mistake here. I believe this should be
      // if (p[n].nres < 0 || p[n].nres > NLINES || geo.scatter_mode == 0)
      // to allow for isotropic BF continuum emission
      nnscat = 1;
      if (p[n].nres < 0 || geo.scatter_mode != 1)
	{
/*  It was either an electron scatter so the  distribution is isotropic, or it
was a resonant scatter but we want isotropic scattering anyway.  */
	  randvec (p[n].lmn, 1.0);	/* The photon is emitted isotropically */
	}
      else if (geo.scatter_mode == 1) 
	{			// It was a line photon and we want anisotropic scattering mode 1

// -1. forces a full reinitialization of the pdf for anisotropic scattering
      
	  randwind (&p[n], p[n].lmn, wmain[icell].lmn);

	}
	  else if (geo.scatter_mode == 2) 
	{			// It was a line photon and we want anisotropic scattering mode 2
	  randwind_thermal_trapping (&p[n], &nnscat);
	}
	p[n].nnscat = nnscat;


      /* The next two lines correct the frequency to first order, but do not result in
         forward scattering of the distribution */

      vwind_xyz (&p[n], v);
      p[n].freq *= (1. + dot (v, p[n].lmn) / C);

      p[n].istat = 0;
      p[n].tau = p[n].nscat = p[n].nrscat = 0;
      p[n].origin = PTYPE_WIND;	// A wind photon

    }


  return (nphot);		/* Return the number of photons generated */

}


/***********************************************************
                                       Space Telescope Science Institute

Synopsis: one_line (one, freqmin, freqmax, nres) gets the frequency of a 
single collisionally excited line photon in a particular cell
of the wind.

Arguments:		

Returns:
 
Description:	

Notes:

History:
	01nov	ksl	Modified to use line pdfs
	02jan	ksl	Modified to return number of resonance line from linptr
	06may	ksl	57+ -- Modified to partially account for plasma structue
			but a further change will be needed when move volume to
			plasma
 
**************************************************************/


double
one_line (one, freqmin, freqmax, nres)
     WindPtr one;
     double freqmin, freqmax;
     int *nres;
{
  double xlum, xlumsum;
  int m, n;
  int lnmin, lnmax;
  int nplasma;
  PlasmaPtr xplasma;

  nplasma = one->nplasma;
  xplasma = &plasmamain[nplasma];

  xlum = xplasma->lum_lines * (rand () / (MAXRAND - 0.5));


// Using the precaclatuted line luminosity pdf's locate the portion of the lin_ptr array 
// from which the photon will will be drawn

  n = 1;
  while (xlum > xplasma->pdf_y[n] && n < LPDF)
    n++;
  lnmin = xplasma->pdf_x[n - 1];
  lnmax = xplasma->pdf_x[n];
  xlumsum = xplasma->pdf_y[n - 1];

/* At this point we know from the line pdf roughly where the line will be located in
frequency space but we need to recalculate the line luminosities within this frequency
space because the line luminosities contained in the line ptr arrays are not current */

  if ((xlumsum + lum_lines (one, lnmin, lnmax)) < xlum)
    {
      Error ("one_line: lumin from lum_lines insufficient\n");
      mytrap ();
    }


  m = lnmin;
  while (xlumsum <= xlum)
    {
      xlumsum += lin_ptr[m]->pow;
      m++;
    }

  *nres = m - 1;

  return (lin_ptr[m - 1]->freq);
}

/* Next section deals with bremsstrahlung radiation */
#define BREMS_CONSTANT 6.85e-38	/*4*PI* 8/3* sqrt(2*PI/3)*e**6/m**2/c**3 sqrt(m/kT) or
				   4*PI times the normal constant to dL_nu/dnu */


/***********************************************************
                                       Space Telescope Science Institute

Synopsis:  total_free calculates the ff luminosity of a cell.

Arguments:		

Returns:
 
Description:	

Notes:

History:
 
**************************************************************/

/* 

Note: program uses an integral formula rather than integrating on
   the fly which is not ideal but saves time 

 	04Apr	SS 	added if statement for case when running hydrogen only. 
	06may	ksl	57+ -- Modified to partially account for plasma structue
			but a further change will be needed when move volume to
			plasma
	07jul	ksl	58f - This routine requires a WindPtr because we still need the
			volume and that is still part of WindPtr
        12sep	nsh	73g - increased the number of ions we will use to all of them!!
	12sep	nsh	73g - incorporated sutherlands data for gaunt factor
	12dec	nsh	74b - put in code to cope with the case where gaunt factor data is not read in
        13jul	nsh	76  - fixed major bugs in this code relating to bug number 29!
*/

double
total_free (one, t_e, f1, f2)
     WindPtr one;
     double t_e;
     double f1, f2;
{
  double g_ff_h, g_ff_he;
  double gaunt;
  double x, sum;
  double gsqrd;			/*The scaled inverse temperature experienced by an ion - used to compute the gaunt factor */
  int nplasma, nion;
  PlasmaPtr xplasma;

  nplasma = one->nplasma;
  xplasma = &plasmamain[nplasma];
  if (t_e < 100.)
    return (0.0);
  if (f2 < f1)
    {
      return (0.0);
    }

//  gaunt=1.0; /*NSH 120920 - this is a placeholder, we need to calculate the gaunt factor at some point */

//  if (nelements > 1)
//    {
/*NSH 120924 - this summation works out the z^2 times number density term for all ions the gaunt factor is calculated for each ion */


  if (gaunt_n_gsqrd == 0)	//Maintain old behaviour
    {
      g_ff_h = g_ff_he = 1.0;
      if (nelements > 1)
	{
	  x =
	    BREMS_CONSTANT * xplasma->ne * (xplasma->density[1] * g_ff_h +
					    4. * xplasma->density[4] *
					    g_ff_he) / H_OVER_K;
	}
      else
	{
	  x =
	    BREMS_CONSTANT * xplasma->ne * (xplasma->density[1] * g_ff_h) /
	    H_OVER_K;
	}
    }
  else
    {
      sum = 0.0;		/*NSH 120920 - zero the summation over all ions */
      for (nion = 0; nion < nions; nion++)
		{
		if (ion[nion].istate !=1) //The neutral ion does not contribute
			{			
			gsqrd=((ion[nion].istate-1)*(ion[nion].istate-1)*RYD2ERGS)/(BOLTZMANN*t_e);
			gaunt=gaunt_ff(gsqrd);
			sum += xplasma->density[nion] * (ion[nion].istate-1) * (ion[nion].istate-1) * gaunt;
			}
		else
			{
			sum += 0.0;
			}
		}
      x = BREMS_CONSTANT * xplasma->ne * (sum) / H_OVER_K;
    }

  x *= sqrt (t_e) * one->vol;
  x *= (exp (-H_OVER_K * f1 / t_e) - exp (-H_OVER_K * f2 / t_e));

  return (x);
}




/***********************************************************
                                       Space Telescope Science Institute

Synopsis: 

Arguments:		

Returns:
 
Description:	

Notes:

History:
 
**************************************************************/
/* Calculate the free-free emissivity in a cell.  Just consider H2+He3 EDIT NSH - now do all ions - its fast so why not?

SS Apr 04: added an "if" statement to deal with case where there's only H. 
	06may	ksl	57+ -- Modified to partially account for plasma structue
			but a further change will be needed when move volume to
			plasma
        12sep	nsh	73g - increased the number of ions we will use to all of them!!
	12sep	nsh	73g - incorporated sutherlands data for gaunt factor
	12dec	nsh	74b - put code in to allow for the case where gaunt data is not available.
        13jul	nsh	76  - fixed major bugs in this code relating to bug number 29!
	
*/

double
ff (one, t_e, freq)
     WindPtr one;
     double t_e, freq;
{
  double g_ff_h, g_ff_he;
  double fnu;
  double gsqrd, gaunt, sum;
  int nplasma;
  int nion;
  PlasmaPtr xplasma;

  nplasma = one->nplasma;
  xplasma = &plasmamain[nplasma];


  if (t_e < 100.)
    return (0.0);



  if (gaunt_n_gsqrd == 0)	//Maintain old behaviour
    {
      g_ff_h = g_ff_he = 1.0;
      if (nelements > 1)
	{
	  fnu =
	    BREMS_CONSTANT * xplasma->ne * (xplasma->density[1] * g_ff_h +
					    4. * xplasma->density[4] *
					    g_ff_he);
	}
      else
	{
	  fnu = BREMS_CONSTANT * xplasma->ne * (xplasma->density[1] * g_ff_h);
	}
    }
  else
    {
      sum = 0.0;
      for (nion = 0; nion < nions; nion++)
		{
		if (ion[nion].istate !=1) //The neutral ion does not contribute
			{			
			gsqrd=((ion[nion].istate-1)*(ion[nion].istate-1)*RYD2ERGS)/(BOLTZMANN*t_e);
			gaunt=gaunt_ff(gsqrd);
			sum += xplasma->density[nion] * (ion[nion].istate-1) * (ion[nion].istate-1) * gaunt;
			}
		else
			{
			sum += 0.0;
			}
		}
      fnu = BREMS_CONSTANT * xplasma->ne * (sum);
    }


  fnu *= exp (-H_OVER_K * freq / t_e) / sqrt (t_e) * one->vol;

  return (fnu);

}


/***********************************************************
                Space Telescope Science Institute

Synopsis: one_ff(one, f1, f2)  determines the frequency of a 
	ff photon within the frequency interval f1 and f2

Arguments:		

Returns:
 
Description:	

Notes:
	one is the windcell where the photon will be created.  It is needed
	only for the temperature.  ?? The code would be simplified
	if simply the temperature were transmitted.

History:
   98           ksl     coded as part of python effort
   98oct        ksl     Removed the internal frequency limits to assure 
			that total ff and one_ff were using
   			the same limits
 
**************************************************************/

struct Pdf pdf_ff;
double ff_x[200], ff_y[200];
double one_ff_f1, one_ff_f2, one_ff_te;	/* Old values */

double
one_ff (one, f1, f2)
     WindPtr one;		/* a single cell */
     double f1, f2;		/* freqmin and freqmax */
{
  double dummy, freq, dfreq;
  int n;
  int nplasma;
  PlasmaPtr xplasma;
  int echeck;

  nplasma = one->nplasma;
  xplasma = &plasmamain[nplasma];

  if (f2 < f1)
    {
      Error ("one_ff: Bad inputs f2 %g < f1 %g returning 0.0  t_e %g\n",
	     f2, f1, xplasma->t_e);
      return (-1.0);
    }

  /* Check to see if we have already generated a pdf */

  if (xplasma->t_e != one_ff_te || f1 != one_ff_f1 || f2 != one_ff_f2)
    {				/* Generate a new pdf */

      dfreq = (f2 - f1) / 199;
      for (n = 0; n < 200; n++)
	{
	  ff_x[n] = f1 + dfreq * n;
	  ff_y[n] = ff (one, xplasma->t_e, ff_x[n]);
	}



      if ((echeck =
	   pdf_gen_from_array (&pdf_ff, ff_x, ff_y, 200, f1, f2, 0,
			       &dummy)) != 0)
	{
	  Error
	    ("one_ff: pdf_gen_from_array error %d : f1 %g f2 %g te %g ne %g nh %g vol %g\n",
	     echeck, f1, f2, xplasma->t_e, xplasma->ne, xplasma->density[1],
	     one->vol);
	  exit (0);
	}
      one_ff_te = xplasma->t_e;
      one_ff_f1 = f1;
      one_ff_f2 = f2;		/* Note that this may not be the best way to check for a previous pdf */
    }
  freq = pdf_get_rand (&pdf_ff);
  return (freq);
}




/***********************************************************
                Southampton university

Synopsis: gaunt_ff computes the frequency averaged gaunt factor at
		scaled temperature. It interpolates between the tabulated
		factors

Arguments:		

Returns:
 
Description:	

Notes:


History:
   12           nsh     coded 
 
**************************************************************/


double
gaunt_ff (gsquared)
     double gsquared;		/* the gamma squared varaiable */
{
  int i, index;
  double gaunt;
  double log_g2;
  double delta;			//The log difference between our G2 and the one in the table

  delta = 0.0;			/* NSH 130605 to remove o3 compile error */
  index = 0;			/* NSH 130605 to remove o3 compile error */

  log_g2 = log10 (gsquared);	//The data is in log format

  if (log_g2 < gaunt_total[0].log_gsqrd
      || log_g2 > gaunt_total[gaunt_n_gsqrd - 1].log_gsqrd)
    {
//              Error ("gaunt ff - request gsqrd outside range - returning gaunt factor =1\n"); /*Removed this as an error - it happens a lot, and is probably not important! */
      return (1.0);
    }
  for (i = 0; i < gaunt_n_gsqrd; i++)	/*first find the pair of parameter arrays that bracket our temperature */
    {
      if (gaunt_total[i].log_gsqrd <= log_g2
	  && gaunt_total[i + 1].log_gsqrd > log_g2)
	{
	  index = i;		/* the array to use */
	  delta = log_g2 - gaunt_total[index].log_gsqrd;
	}
    }

  gaunt =
    gaunt_total[index].gff + delta * (gaunt_total[index].s1 +
				      delta * (gaunt_total[index].s2 +
					       gaunt_total[index].s3));





  return (gaunt);
}
