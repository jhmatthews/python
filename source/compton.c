


/**************************************************************************
                    Space Telescope Science Institute


  Synopsis:  These are the routines to handle compton scattering.

  Description:	There are two routines in here:
		kappa_compton (xplasma,freq) this calculates the opacity in the cell xplasma due to compton cooling

  Arguments: 		


  Returns:

  Notes:

  History:
2011	nsh	Coded as part of the effort to include compton scattering in August 2011 at Southampton.

 ************************************************************************/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "atomic.h"
#include "python.h"
#include <gsl/gsl_rng.h>


/**************************************************************************
                    Space Telescope Science Institute


  Synopsis:  kappa_compton computes the opacity in the cell due to compton heating.

  Description:	

  Arguments:   xplasma - pointer to the current plasma cell
		freq - the frequency of the current photon packet being tracked

  Returns:   kappa - the compton opacity for the cell.

  Notes:   This implements the equation kappa=(sigmaT*ne*freq)/(me*c^2)
            Initally it is called by radiation

  History:
2011	nsh	Coded as part of the effort to include compton scattering in August 2011 at Southampton.
feb 2013 - nsh - approximate KN cross section replaced by correct value

 ************************************************************************/


double
kappa_comp (xplasma, freq)
     PlasmaPtr xplasma;		// Pointer to current plasma cell
     double freq;		// Frequency of the current photon being tracked
{
  double x;			// The opacity of the cell by the time we return it.
  double sigma;			/*The cross section, thompson, or KN if hnu/mec2 > 0.01 */

/*	alpha=1/(1+freq*HRYD*(1.1792e-4+(7.084e-10*freq*HRYD))); NSH 130214 This is the approximate way of doing it.*/

 //	sigma=THOMPSON/(1+freq*HRYD*(1.1792e-4+(7.084e-10*freq*HRYD)));


  sigma = klein_nishina (freq);	//NSH 130214 - full KN formula

  x = (sigma * H) / (MELEC * C * C);	//Calculate the constant
  x *= xplasma->ne * freq;	//Multiply by cell electron density and frequency of the packet.

  x *= geo.fill;    // multiply by the filling factor- should cancel with density enhancement
  return (x);
}



/**************************************************************************
                    Space Telescope Science Institute


  Synopsis:  kappa_ind_compton computes the opacity in the cell due to induced compton heating.

  Description:	

  Arguments:   xplasma - pointer to the current plasma cell
		freq - the frequency of the current photon packet being tracked

  Returns:   kappa - the inducd compton opacity for the cell.

  Notes:   This implements the induced compton term in equation 6.5 in cloudy 10.3

  History:
2011	nsh	Coded as part of the effort to include compton scattering in August 2011 at Southampton.
feb 2013 - nsh - approximate KN cross section replaced by correct value
1309 	JM	Removed w and ds as arguments as no longer required
1312	NSH	Recoded to use the new exponential and power law models.

 ************************************************************************/


double
kappa_ind_comp (xplasma, freq)
     PlasmaPtr xplasma;		// Pointer to current plasma cell
     double freq;		// Frequency of the current photon being tracked
     //double w;			// The weight of the photon packet
     //double ds;			//The distance the photon travels
{
  double x;			// The opacity of the cell by the time we return it.
  double sigma;			/*The cross section, thompson, or KN if hnu/mec2 > 0.01 */
  double J;		//The estimated intensity in the cell
//  int i;

  /* Previously, NSH had used the following formula whixch required ds and w to work  
     J=(4*PI*w*ds)/(C*xplasma->vol); //Calcuate the intensity NSH This works for a thin shell... Why? Dont know.
  */ 

  J = 0.0;			/* NSH 130605 to remove o3 compile error */

  /* Obtain a model for the mean intensity - we call this with mode=2, which means 
     that if we have not yet completed a cycle, dont return a dilute blackbody 
     estimate if we are in PL mode. */
  J = mean_intensity (xplasma, freq, 2); 

  /* 1407 -- JM -- There was a lot of commented out code here which I've deleted-
     NSH moved it into the mean_intensity subroutine. See Pull Request #88 */


  sigma = klein_nishina (freq);	//NSH 130214 - full KN formula


  x = (xplasma->ne) / (MELEC);
  x *= sigma * J;		// NSH 130214 factor of THOMPSON removed, since alpha is now the actual compton cross section
  x *= 1 / (2 * freq * freq);

  x *= geo.fill;    // multiply by the filling factor- should cancel with density enhancement

  if (sane_check (x)) //For some reason we have a problem
    {
      Error
	("kappa_ind_comp:sane_check - undefined value for Kappa_ind_comp - setting to zero\n");
      return (0.0);
    }

  return (x);
}

/**************************************************************************
                    Space Telescope Science Institute


  Synopsis:  total_comp computes the cooling in the cell due to compton cooling.

  Description:	

  Arguments:   one - pointer to the current wind cell
		t_e - electron temperature of the current cell

  Returns:   the compton luminosity for the cell.

  Notes:   This implements the equation C=16piVsigmatJ(kTe/me)

  History:
2011	nsh	Coded as part of the effort to include compton scattering in August 2011 at Southampton.

 ************************************************************************/

double
total_comp (one, t_e)
     WindPtr one;		// Pointer to the current wind cell - we need the cell volume, this is not in the plasma structure
     double t_e;		//Current electron temperature of the cell
{
  double x;			//The returned variable
  int nplasma;			//The cell number in the plasma array
  PlasmaPtr xplasma;		//pointer to the relevant cell in the plasma structure


  nplasma = one->nplasma;	//Get the correct plasma cell related to this wind cell
  xplasma = &plasmamain[nplasma];	//copy the plasma structure for that cell to local variable

  x = 16. * PI * THOMPSON * BOLTZMANN / (MELEC * C * C);	//Keep all the constants together
  x *= xplasma->ne * xplasma->vol * xplasma->j * t_e;	//multiply by the volume (from wind) and j (from plasma) and t_e


  return (x);
}



/**************************************************************************
                    Space Telescope Science Institute


  Synopsis:  klein_nishina computes the KN cross section for a photon of frequency nu.

  Description:	

  Arguments:   nu - photon frequency


  Returns:   the KN cross section.

  Notes:   This implements equation 7.5 in Rybicki and Lightman

  History:
2013	nsh	Coded

 ************************************************************************/

double
klein_nishina (nu)
     double nu;			//The frequency of the photon packet
{
  double x;			//h nu / kt
  double x1, x2, x3, x4;	//variables to store intermediate results.
  double kn;			// the final cross section

  kn = THOMPSON;		/* NSH 130605 to remove o3 compile error */
  x1 = x2 = x3 = x4 = 0.0;	/* NSH 130605 to remove o3 compile error */
  x = (H * nu) / (MELEC * C * C);
  if (x > 0.0001)
    {
      x1 = 1. + x;
      x2 = 1. + (2. * x);
      x3 = log (x2);
      x4 = ((2. * x * x1) / x2) - x3;
      x4 *= x1 / (x * x * x);
      x4 = x4 + x3 / (2. * x);
      x4 = x4 - (1 + 3. * x) / (x2 * x2);
      kn *= 0.75 * x4;
    }

  return (kn);
}

double z_rand,sigma_tot,x1; //External variables to allow zfunc to search for the correct fractional energy change

/**************************************************************************
                    Southampton University


  Synopsis:  compton_dir computes a random direction (and hence frequency change)
			for a photon undergoing compton scattering. At low frequencies, this
			is just thompson scattering.

  Description:	

  Arguments:  
			p - the photon currently being scattered, this gives us the current direction and the frequency
			xplasma- pointer to the current plasma cell

  Returns:   kappa - nothing, but updates the photon frequency and direction

  Notes:   

  History:
2015	NSH coded as part of teh summer 2015 code sprint

 ************************************************************************/



int
	compton_dir (p,xplasma)
		PhotPtr p;
        PlasmaPtr xplasma;		// Pointer to current plasma cell
		
{
	double f_min,f_max,f; //The theoretical maxmimum energy change
	double n,l,m,phi,len;  //The direction cosines of the new photon direction in the frame of reference with q along the photon path
	struct basis nbasis;  //The basis function which transforms between the photon frame and the pbserver frame	
	double lmn[3]; /* the individual direction cosines in the rotated frame */
	double x[3]; /*photon direction in the frame of reference of the original photon*/
	double dummy[3],c[3];
	
	x1=H*p->freq/MELEC/C/C; //compute the ratio of photon energy to electron energy. In the electron rest frame this is just the electron rest mass energ 
		
	n=l=m=0.0;  //initialise some variables to avoid warnings
	
	
	if (x1<0.0001) //If the photon energy is much less than electron mass, we just have thompson scattering.
	{
		randvec(lmn,1.0); //Generate a normal isotropic scatter
		f=1.0;   //There is no energy loss
		stuff_v (lmn, p->lmn);
	}
	else
	{
		z_rand=rand()/MAXRAND; //Generate a random number between 0 and 1 - this is the random location in the klein nishina scattering distribution - it gives the energy loss and also direction.
		f_min=1.;         //The minimum energy loss - i.e. no energy loss
		f_max=1.+(2.*x1);  //The maximum energy loss
	
		sigma_tot=sigma_compton_partial(f_max,x1);  //Communicated externally to the integrand function in the zbrent call below, this is the maximum cross section, used to scale the K_N function to lie between 0 and 1.
	
		f=zbrent(compton_func,f_min,f_max,1e-8);  //Find the zero point of the function compton_func - this finds the point in the KN function that represents our random energy loss.
		n=(1.-((f-1.)/x1));   //This is the angle cosine of the new direction in the frame of reference of the photon
//		printf ("f=%e n=%e fmin=%e fmax=%e\n",f,n,f_min,f_max);
	
		if (isfinite(len=sqrt(1.-(n*n)))==0)   //Compute the length of the other sides - the isfinite is to take care of the very rare occasion where n=1!
		len=0.0;
		phi = 0.0;           //no need to randomise phi, the random rotation of the vector generating the basis function takes care of this

		l = len*cos (phi);   //compute the angle cosines of the other two dimensions.
		m = len*sin (phi);
		
		randvec(dummy,1.0);  //Get a random vector
	   
		cross(dummy, p->lmn, c);  //c will be perpendicular to p->lmn
				
//		printf ("c= %e %e %e n=%e acos(n)=%f\n",c[0],c[1],c[2],n,acos(n)*RADIAN);
		
		create_basis (p->lmn, c, &nbasis); //create a basis with the first axis in the direction of the original photon direction, c will be perpendicular to the photon direction. Orientaion of the y/z axes are will give the randomization of the phi axis
		
		
		x[0] = n;  //This is the cosine direction of the new photon direction, in the frame of reference where the first axis is in the original photon direction
		x[1] = l;
		x[2] = m;
	  

      project_from (&nbasis, x, lmn);	/* Project the vector from the FOR of the original photon into the observer frame */
		renorm(lmn,1.0);
//		Log("f=%e freq=%e  n=%e l_old=%e m_old=%e n_old=%e l_new=%e m_new=%e n_new=%e len=%e\n",f,p->freq,n,p->lmn[0],p->lmn[1],p->lmn[2],lmn[0],lmn[1],lmn[2],length(lmn));
//	  Log_flush();
//	  renorm (lmn,1.0);

		stuff_v (lmn,p->lmn);
//		randvec(a,1.0); //Generate a normal isotropic scatter
//				f=1.0;   //There is no energy loss
//		stuff_v (a, p->lmn);
//	printf ("Original %e %e %e theta %e new %e %e %e\n",pold.lmn[0],pold.lmn[1],pold.lmn[2],n,p->lmn[0],p->lmn[1],p->lmn[2]);		
		
	}
	
	p->freq=p->freq/f;  //reduce the photon frequency
	p->w=p->w/f;  //reduce the photon weight by the same ammount to conserve photon numbers
return(0);
}

/**************************************************************************
                    Southampton University


  Synopsis:  compton_func is a simple function that is equal to zero when the 
			external variable z_rand is equal to the normalised KN cross section.
			It is used in a call to zbrent in the function compton_dir

  Description:	

  Arguments:  
			f- the fractional energy change which is supplied to sigma_compton_partial

  Returns:   the difference between sigma(f) and the randomised cross section we are searching for

  Notes:   Since this function is called by zbrent, many of the variables have to be communicated
			externally. These are 
				x1, the ratio of photon energy to electron rest mass,
				sigma_tot, the total angle integrated KN cross section
				z_rand - a randomised number between 0 and 1 representing the normalised cross section we want to find
		

  History:
2015	NSH coded as part of the summer 2015 code sprint

 ************************************************************************/


double
	compton_func(f)
		double f;
{
	double ans;
	ans=(sigma_compton_partial(f,x1)/sigma_tot)-z_rand;
		return(ans);
}

/**************************************************************************
                    Southampton University


  Synopsis:  sigma_compton_partial is the KN cross section as a function of the 
			fractional energy change

  Description:	

  Arguments:  
			f - the fractional energy change
			x	the enerrgy of the photon divided by the rest mass energy of an electron

  Returns:   the cross section for the scattering angle that gives this fractional energy change.

  Notes:   
			Stuart Sim supplied this formula, but coundn't recall where he had found it, although
			he had rederived it and thinks it is correct!

  History:
2015	NSH coded as part of the summer 2015 code sprint

 ************************************************************************/


double 
   sigma_compton_partial(f,x)
      double f; //This is the fractional energy change, nu/nu'
   double x; //h nu/mec**2 - the energy of the photon divided by the rest energy of an eectron
   {
      double term1,term2,term3,tot;
   
   term1 = ( (x*x) - (2*x) - 2 ) * log(f) / x / x;
   term2 = ( ((f*f) -1) / (f * f)) / 2;
   term3 = ( (f - 1) / x) * ( (1/x) + (2/f) + (1/(x*f)));

   tot = 3 * THOMPSON * (term1 + term2 + term3) / (8 * x);

   return(tot);
   
}
   
   
   
   
   
    


