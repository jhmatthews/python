/*

                                       Space Telescope Science Institute

Synopsis:
	These are all of the routines necessary to define a YSO wind, which comprises
	a KWD wind in some places and a stellar wind elsewhere.  The model is for
	a stellar wind within the inner windcone.
Arguments:		

Returns:
 
Description:	

Notes:


History:
 	04jun	ksl	Added this possibility, beginning with the same code for a kwd wind
			Where possible, I have simply chosen to get the data needed
			from the original kwd or stellar wind routines.

 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "atomic.h"
#include "python.h"


/***********************************************************
                                       Space Telescope Science Institute

Synopsis:
	get_yso_wind_params gets input data that are necessary for a YSO
	comprising a region that resembles a KWD wind and a region outside
	this that contins a speherical wind
Arguments:		

Returns:
 
Description:	
	Basically this routine just reads the inputs from the two separate
	wind models, since this is a combination of the two.
Notes:
	The model has a spherical wind on the interior.  The windcones
	are set to include the region from the polar axis to the outside
	edge of the wind.  (Note that when filling out the grid regions
	outside the wind may have non-zero values of quantities, such
	as density.  This allows for interpolation to occur. But there
	will not be (should not be) any scattering or whatever in this
	region since that wind cones define the regions where photons
	are scattered.  


History:
 	04jun	ksl	Added this possibility
        04Sep   SS      Changed finite disk case to cut off wind
                        at the top (rather than bottom) outer edge of disk.
**************************************************************/

int
get_yso_wind_params ()
{

/* The approach to get the input parameters is to call both input parameter routines
one after the other*/

  get_stellar_wind_params ();
  get_knigge_wind_params ();

/* Assign the generic parameters for the wind the generic parameters of the wind */

  geo.wind_rmin = geo.rstar;
  geo.wind_rmax = geo.rmax;
  geo.wind_thetamin = 0.0;
/* Somewhat paradoxically diskrad is in cm, while dn_ratio which is really d in KWD95 is 
in units of WD radii */
  geo.wind_thetamax = atan (geo.diskrad / (geo.kn_dratio * geo.rstar));
// Line above would be 90 degeres if we want a stellar wind outside the windcone
  /* Next lines added by SS Sep 04. Changed the wind shape so that the boundary touches the outer 
     corner of the disk rather than the intersection of the disk edge with the xy-plane. */

  if (geo.disk_type == 2)
    {
      geo.wind_thetamax =
	atan (geo.diskrad /
	      (((geo.kn_dratio * geo.rstar) + zdisk (geo.diskrad))));
    }


  geo.wind_rho_min = geo.rstar;
  geo.wind_rho_max = geo.diskrad;

  /* The change in the boundary of the wind (as corner of disk -- see above) 
     means that wind_rho_max nees to be redefined so that it is used correctly
     to compute the boundary of the wind elsewhere. */

  if (geo.disk_type == 2)
    {
      geo.wind_rho_max =
	geo.diskrad - (zdisk (geo.diskrad) * tan (geo.wind_thetamax));
    }

  /* if modes.adjust_grid is 1 then we have already adjusted the grid manually */
  if (modes.adjust_grid == 0)
    {
      geo.xlog_scale = geo.rstar;
      geo.zlog_scale = 1e7;
    }

  return (0);
}



/***********************************************************
                                       Space Telescope Science Institute

 Synopsis:
	double yso_velocity(x,v) calculates the v in cartesian coordinates
	of a yso wind from a position x in cartesian coordinates.  
Arguments:		
	double x[]		the position where for the which one desires the velocity
Returns:
	double v[]		the calculated velocity in cartesina coordinates
	
	The amplitude of the velocity is returned 
	
Description:	
		
Notes:

History:
	04jun	ksl	The routine simple chooses between the way to calculate
                        velocity gradients and uses the original routines
	04aug	ksl	52 -- Modified underlying routines so that the velocity
			return is in cartesian coordinates.
 
**************************************************************/

double
yso_velocity (x, v)
     double x[], v[];
{
  double r, rzero;
  double zzz;
  double dd;
  double fabs ();

/* First step is to determine whether we are inside the windcone of the kwd wind. 
This is done as it was done in kn_velocity by seeing if the kwd streamline hits
the disk.  If it does not, then we calculate the velocity for a stellar wind */

  dd = geo.rstar * geo.kn_dratio;	//dd is the distance of the focus point in cm 
  r = sqrt (x[0] * x[0] + x[1] * x[1]);	// rho of the point we have been given
  rzero = r / (1. + fabs (x[2] / dd));	// Nothing prevents this from being less than R_wd

  if (rzero < geo.rstar)
    zzz = stellar_velocity (x, v);
  else
    zzz = kn_velocity (x, v);


  return (zzz);

}

/***********************************************************
                                       Space Telescope Science Institute

 Synopsis:
	double yso_rho(x) calculates the density of an kn_wind at a position x
Arguments:		
	double x[]	the position where for the which one desires the density
Returns:
	The density at x is returned in gram/cm**3
	
Description:	
		
Notes:

History:
 	04jun	ksl	Began work, using kn_rho as starting point.  The routine simple
                        chooses between the way to calculate velocity gradients and 
			uses the original routines

 
**************************************************************/

double
yso_rho (x)
     double x[];
{
  double r, rzero;
  double dd;
  double rho;

  dd = geo.rstar * geo.kn_dratio;
  r = sqrt (x[0] * x[0] + x[1] * x[1]);	//rho coordinate of the point we have been given
  rzero = r / (1. + fabs (x[2] / dd));	//rho at the base for this streamline

  if (rzero < geo.rstar)
    {
      rho = stellar_rho (x);
    }
  else
    rho = kn_rho (x);

  return (rho);
}


/*
   Next routine is eliminated in py52 in favor of a generic routine 04 aug

yso_vel_grad calculates the velocity radient tensor at any point in
the flow

The velocity gradient is defined as a 3 x 3 tensor such that

velgrad[i][j]= dv_i/dx_j


        04jun   ksl     Began work on this.  The routine simple
			chooses between the way to calculate
			velocity gradients and uses the original
			routines

*/
/*
int
yso_vel_grad (x, velgrad)
     double x[], velgrad[][3];
{
  double dd, r, rzero;

  dd = geo.rstar * geo.kn_dratio;
  r = sqrt (x[0] * x[0] + x[1] * x[1]);	//rho coordinate of the point we have been given
  rzero = r / (1. + fabs (x[2] / dd));	//rho at the base for this streamline

  if (rzero < geo.rstar)
    {
      stellar_vel_grad (x, velgrad);
    }
  else
    kn_vel_grad (x, velgrad);

  return (0);
}
*/
