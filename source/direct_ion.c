#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "atomic.h"
#include "python.h"
#include "recipes.h"

#include <gsl/gsl_sf_expint.h> //We need this gsl library to evaluate the first exponential integral





int
compute_di_coeffs (T)
     double T;
{
  int i,n,nrec;
  double rate, drdt, t, dt, scaled_t;
  int imin,imax;
  double rates[DERE_DI_PARAMS],temps[DERE_DI_PARAMS];
  double exp_int;


  for (n = 0; n < nions; n++)
    {
      if (ion[n].dere_di_flag == 0)
	{
//printf ("NO COEFFS\n");
	  di_coeffs[n] = 0.0;
	}
      else
	{
nrec=ion[n].nxderedi;
/* find the correct coefficient */
	  


  for (i = 0; i < dere_di_rate[nrec].nspline; i++)
    {
      rates[i] = dere_di_rate[nrec].rates[i];
      temps[i] = dere_di_rate[nrec].temps[i];
    }


t=(BOLTZMANN*T)/(dere_di_rate[nrec].xi*EV2ERGS);
scaled_t=1.0-((log(2.0))/(log(2.0+t)));

  if (scaled_t < temps[0])		//we are below the range of DI data data
    {
      Log_silent
	("bad_gs_rr: Requested temp %e is below limit of data for ion %i(Tmin= %e)\n",
	 scaled_t, n, temps[0]);
//      rate = rates[0];
     	imax=1;
	imin=0;  
  }

  else if (scaled_t >= temps[dere_di_rate[nrec].nspline - 1])	//we are above the range of GS data
    {
      Log_silent
	("bad_gs_rr: Requested temp %e is above limit (%e) of data for ion %i\n",
	 scaled_t, temps[dere_di_rate[nrec].nspline - 1],n);
 //     rate = rates[BAD_GS_RR_PARAMS - 1];
	imax=dere_di_rate[nrec].nspline - 1;
	imin=dere_di_rate[nrec].nspline - 2;
//We will try to extrapolate.



    }
  else				//We must be within the range of tabulated data
    {
      for (i = 0; i < dere_di_rate[nrec].nspline - 1; i++)
	{
        if (temps[i] <= scaled_t && scaled_t < temps[i + 1])	//We have bracketed the correct temperature
	    {
	      imin = i;
	      imax = i + 1;
	    }
	}
/* NSH 140313 - changed the following lines to interpolate in log space */
    }
//if (ion[n].z==6 && ion[n].istate==4){
//
//printf ("T=%e Interploating between %e(%e) and %e(%e)\n",scaled_t,rates[imin],temps[imin],rates[imax],temps[imax]); 
//}
      drdt = ((rates[imax]) - (rates[imin])) / ((temps[imax]) - (temps[imin]));
      dt = ((scaled_t) - (temps[imin]));
      rate = rates[imin] + drdt * dt;


//printf ("t=%e xi=%e scaled_rate=%e\n",t,dere_di_rate[nrec].xi,rate);


di_coeffs[n]=pow(t,-0.5)*pow(dere_di_rate[nrec].xi,-1.5)*rate;
//printf ("Doing 1/t=%e\n",1.0/t);

if (exp(-1.0/t) < (1.0/VERY_BIG))
	{
	exp_int=0.0;
	}
	else
	{
exp_int=gsl_sf_expint_E1(1.0/t); //Evaluate the first exponential integral using gsl library.
	}

di_coeffs[n]*=exp_int;
}



} //End of loop over ions



return(0);
}	 

double
total_di (one, t_e)
     WindPtr one;		// Pointer to the current wind cell - we need the cell volume, this is not in the plasma structure
     double t_e;		//Current electron temperature of the cell

{
  double x;			//The returned variable
  int nplasma;			//The cell number in the plasma array
  PlasmaPtr xplasma;		//pointer to the relevant cell in the plasma structure
  int n;			//loop pointers


  nplasma = one->nplasma;	//Get the correct plasma cell related to this wind cell
  xplasma = &plasmamain[nplasma];	//copy the plasma structure for that cell to local variable
  x = 0;			//zero the luminosity


  compute_di_coeffs (t_e);	//Calculate the DR coefficients for this cell


  for (n = 0; n < nions; n++)
    {
      if (ion[n].dere_di_flag == 0)	//We have no DR for this ion.
	{
	  x += 0.0;		//Add nothing to the sum of coefficients
	}
      else
	{

	  x +=
	    xplasma->vol * xplasma->ne * xplasma->density[n] * di_coeffs[n] *
	    dere_di_rate[ion[n].nxderedi].xi*EV2ERGS;
//printf ("n=%i V=%e ne=%e rho=%e coeff=%e xi=%e cooling=%e\n",n, V , xplasma->ne , xplasma->density[n] , di_coeffs[n] ,
//	    dere_di_rate[ion[n].nxderedi].xi*EV2ERGS,x);
	}
    }
  return (x);
}

 
