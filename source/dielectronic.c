/**************************************************************************
                    Southampton University
                                                                                                   
                                                                                                   
  Synopsis:

The routines in this file all have to do with dielectronic recombination.


  Description:

The first routine is compute_dr_coeffs
It takes as its input a temperature, and it populates the fd_coeff array.
This is an initial test, and it will be called every time a cell needs a set
of DR coefficients. It could well be quicker in the future to generate one,
temperature split array, and then sample it. Lets see...



                                                                                                   
  Arguments:  
                                                                                                   
                                                                                                   
  Returns:
                                                                                                   
  Notes:
                                                                                                   
                                                                                                   

                                                                                                   
  History:
	11aug	nsh	Began work

                                                                                                   
 ************************************************************************/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "atomic.h"
#include "python.h"
#include "recipes.h"




/**************************************************************************
                    Space Telescope Science Institute
                                                                                                   
                                                                                                   
  Synopsis: compute_dr_coeffs returns the volumetric dielectronic rate
 	coefficients for a given temperature. 
                                                                                                   
  Description:
                                                                                                   
  Arguments:  
	temperature
                                                                                                   
                                                                                                   
  Returns:
	nothing, but populates the array dr_coeffs
                                                                                                   
  Notes:
 	the rates are associated with the ion being recombined into. 
                                                                                                   
                                                                                                   
                                                                                                   
  History:
	11sep	nsh	Written as part of python70 effort to incorporate DR.
                                                                                                   
 ************************************************************************/

int 
compute_dr_coeffs (double temp)
{
  int n, n1, n2;
  double Adi, Bdi, T0, T1;
  for (n = 1; n < nions+1; n++)
    {
      if (ion[n].drflag == 0)
	{
	  dr_coeffs[n] = 0.0;
	}
      else
	{
	  n1 = ion[n].nxdrecomb;
	  dr_coeffs[n] = 0.0;
	  if (drecomb[n1].type == DRTYPE_BADNELL)
	    {
	      for (n2 = 0; n2 < drecomb[n1].nparam; n2++)
		{
		  dr_coeffs[n] +=
		    (drecomb[n1].c[n2] *
		     exp (-1 * (drecomb[n1].e[n2] / temp)));
		}
	      dr_coeffs[n] *= pow (temp, -1.5);
	    }
	  else if (drecomb[n1].type == DRTYPE_SHULL)
	    {
	      Adi = drecomb[n1].shull[0];
	      Bdi = drecomb[n1].shull[1];
	      T0 = drecomb[n1].shull[2];
	      T1 = drecomb[n1].shull[3];
	      dr_coeffs[n] =
		Adi * pow (temp, -1.5) * exp ((-1.0 * T0) / temp);
	      dr_coeffs[n] *= (1 + Bdi * exp ((-1.0 * T1) / temp));
	    }
	  else
	    {
	      Error ("Compute_dr_coeffs: Unknown DR data rtype for ion %i\n",
		     n);
	      dr_coeffs[n] = 0.0;
	    }
	}
    }
  return (0);
}


/**************************************************************************
                    Space Telescope Science Institute
                                                                                                   
                                                                                                   
  Synopsis: total_dr calculates the total luminosity from DR. 
                                                                                                   
  Description:
                                                                                                   
  Arguments:  
	pointer to grid cell we are interested in	
	temperature
	
                                                                                                   
                                                                                                   
  Returns:
  	the total luminosity of this cell due to dielectronic recombinaions
                                                                                                   
  Notes:
 	
                                                                                         
                                                                                                   
                                                                                                   
  History:
	11sep	nsh	Written as part of python70 effort to incorporate DR. Initally we are just doing a ROM calculation by multiplying the volumetric rate be the ion density, the electron density and mean eelectron energy.
        12jul 	nsh	Changed to take account of the fact tha we are now assiciating a rate with the ion being recombined into, this means that for rate[nion] we need density [nion+1].
                                                                                                   
 ************************************************************************/

double 
total_dr (WindPtr one, double t_e)		//Current electron temperature of the cell
{
  double x;			//The returned variable
  double meanv, meanke;		//The mean velocity and kinetic energy of electrons in the cell
  int nplasma;			//The cell number in the plasma array
  PlasmaPtr xplasma;		//pointer to the relevant cell in the plasma structure
  int n;			//loop pointers


  nplasma = one->nplasma;	//Get the correct plasma cell related to this wind cell
  xplasma = &plasmamain[nplasma];	//copy the plasma structure for that cell to local variable
  x = 0;			//zero the luminosity


  compute_dr_coeffs (t_e);	//Calculate the DR coefficients for this cell
  meanv = pow ((2 * BOLTZMANN * t_e / MELEC), 0.5);
  meanke = 0.5 * MELEC * meanv * meanv;

  for (n = 1; n < nions+1; n++)
    {
      if (ion[n].drflag == 0)	//We have no DR for this ion.
	{
	  x += 0.0;		//Add nothing to the sum of coefficients
	}
      else
	{
	  x +=
	    xplasma->vol * xplasma->ne * xplasma->density[n +
						      1] * dr_coeffs[n] *
	    meanke;
	}
    }
  return (x);
}
