#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "atomic.h"

#include "python.h"
/***********************************************************
                                       Space Telescope Science Institute

 Synopsis:   
	translate either calls translate_in_space or translate_in_wind  depending upon the
	current location of the photon.  
  
 Arguments:		

	
 Returns:
  
Description:	
 	

Notes:
	In Monte, translate controls the flow of the photon through one grid cell.  In Python,
	there are additional possibilities since there are wind free regions.  This construction is
	used so that the calling routines for translate (trans_phot and extract) can be the same
	in the 1 and 2d versions.

History:
 	1997	ksl	Coded and debugged as part of Python effort.
 	98nov	ksl	Modified call to where_in_wind 
 
**************************************************************/

int
translate (w, pp, tau_scat, tau, nres)
     WindPtr w;
     PhotPtr pp;
     double tau_scat;
     double *tau;
     int *nres;
{
  int istat;


//      if(where_in_wind(w,pp)!=0) {
  //      if(where_in_wind(pp->x)!=0) {
  //              istat=translate_in_space(pp);
  //      }
  //      else if ((pp->grid=where_in_grid(pp->x))>=0) {
  istat = translate_in_wind (w, pp, tau_scat, tau, nres);
//      }
  //      else istat=pp->istat=-1;  /* It's not in the wind and it's not in the grid.  Bummer! */


  return (istat);


}

/***********************************************************
                                       Space Telescope Science Institute

Synopsis:   

	int translate_in_space(pp) translates the photon from its current position to the 
	edge of the wind. 

Arguments:		
	PhotPtr pp; 
	
Returns:
  
Description:	
 	

Notes:
	?? Why is it necessary to check whether the photon has hit the star at 
	this point??

History:
 	1997	ksl	Coded and debugged as part of Python effort. 
 
**************************************************************/


int
translate_in_space (pp)
     PhotPtr pp;
{
  double ds, x;

  double ds_to_wind (), ds_to_sphere ();

  ds = ds_to_wind (pp);

  if ((x = ds_to_sphere (geo.rstar, pp)) < ds)
    {
      x = ds;
      pp->istat = P_HIT_STAR;	/* Signifying that photon is hitting star */
    }
  move_phot (pp, ds + DFUDGE);
  return (pp->istat);
}

/***********************************************************
                                       Space Telescope Science Institute

 Synopsis:   
	double ds_to_wind(pp)  calculates the photon pathlength to the edge of the wind.  

  
  
 Arguments:		
	PhotPtr pp;.
	
 Returns:
 	The distance to the boundary of the wind
  
Description:	
	If you are inside the wind already ds_to_wind calculates the distance to the edge of the wind. 
	If you are outside, It will also be to the edge.  	

Notes:
	There is no guarantee that you will still be in the region defined by the grid.

History:
 	1997	ksl	Coded and debugged as part of Python effort. 
 
**************************************************************/


double
ds_to_wind (pp)
     PhotPtr pp;
{
  double ds, ds_to_cone (), x;

  ds = ds_to_cone (&windcone[0], pp);

  if ((x = ds_to_cone (&windcone[1], pp)) < ds)
    ds = x;

  return (ds);
}

/***********************************************************
                                       Space Telescope Science Institute

 Synopsis:   
	translate_in_wind translates the photon within a single cell in the wind.  
  
  
 Arguments:		

	
 Returns:
  
Description:	
 	It
	calculates and updates the final position of the photon p, the optical depth tau after 
	having translated the photon, and the number of the resonance.  These last two quantities 
	are caculated in ds_calculate and simply passed back.  translate_in_wind stores
	the final status of the photon in p and returns that status directly to the calling routine.	

Notes:


History:
 	1997	ksl	Coded and debugged as part of Python effort. 
 
**************************************************************/



int
translate_in_wind (w, p, tau_scat, tau, nres)
     WindPtr w;
     PhotPtr p;
     double tau_scat, *tau;
     int *nres;


{

  struct photon phot;
  int n, ix, iz;
  double smax, ds_current;
  double calculate_ds ();
  double ds_to_wind ();
  int istat;



  if ((p->grid = n = where_in_grid (p->x)) < 0)
    {
      Error ("translate_in_wind: Photon not in grid when routine entered\n");
      return (n);		/* Photon was not in wind */
    }

  wind_n_to_ij (n, &ix, &iz);	/*Convert the index n to two dimensions */

  smax = geo.pl_vol - p->x[2];


  /* At this point we now know how far the photon can travel in it's current
     grid cell */


  smax += DFUDGE;		/* DFUDGE is to force the photon through the cell boundaries */

  /* We now determine whether scattering prevents the photon from reaching the far edge of
     the cell.  calculate_ds calculates whether there are scatterings and makes use of the current
     position of the photon and the position of the photon at the far edge of the shell.  It needs a
     "trial photon at the maximimum distance however */

  stuff_phot (p, &phot);
  move_phot (&phot, smax);

//OLD  ds_current = calculate_ds (w, p, &phot, tau_scat, tau, nres, smax, &istat);
  ds_current = calculate_ds (w, p, tau_scat, tau, nres, smax, &istat);

/* OK now we increment the radiation field in the cell, translate the photon and wrap things up */

  radiation (p, ds_current);

  move_phot (p, ds_current);

  return (p->istat = istat);


}

/***********************************************************
                                       Space Telescope Science Institute

 Synopsis:   
	 walls determines whether the photon has encountered the star of disk or
	 reached the edges of the grid and returns the appropriate
	status.  
  
 Arguments:		
	PhotPtr p,pold		the current and previous descripiton of the photon bundle.
	
 Returns:
         The photon has hit the star                       P_HIT_STAR
         The photon has hit the disk                       P_HIT_DISK
         The photon has reached the edge of grid    P_ESCAPE
         The status is undeterminable                    5
         
 If a photon does not fall into one of these categories, walls returns the old status, which is stored in
 p->istat 
 
Description:	
 	

Notes:


History:
 	1997	ksl	Coded and debugged as part of Python effort. 
 	1997nov	ksl	Corrected problem which caused mistakes in the calculation of the disk
 	 			intercept.	 

**************************************************************/

int
walls (p, pold)
     PhotPtr p, pold;
{
  double dot ();

//      /* Check to see if the photon has hit the star */
  //      if((r=dot(p->x,p->x))<geo.rstar_sq) return(p->istat=P_HIT_STAR);

  /* Check to see if it has hit the disk */

//      if(geo.idisk==0 && p->x[2]*pold->x[2]<0.0){
  //              s=(-(pold->x[2]))/(pold->lmn[2]);
  //              if(s<0) {
  //                      Error("walls: distance %g<0\n",s);
  //                      exit(0);
  //              }
  //              vmove(pold->x,pold->lmn,s,xxx);
  //              if (dot(xxx,xxx) <geo.diskrad_sq) { /* The photon has hit the disk*/
  //                      stuff_phot(pold,p);     /* Move the photon to the point where it hits the disk */
  //                      move_phot(p,s);
  //                      return(p->istat=P_HIT_DISK);  
  //              }
  //      }

//      rho_sq=(p->x[0]*p->x[0]+p->x[1]*p->x[1]);
  //      if(rho_sq>geo.rmax_sq)      return(p->istat=P_ESCAPE);  /* The photon is coursing through the universe */
  //      if(fabs(p->x[2])>geo.rmax)  return(p->istat=P_ESCAPE); 

  if (p->x[2] >= geo.pl_vol)
    return (p->istat = P_ESCAPE);
  return (p->istat);		/* The photon is still in the wind */

}
