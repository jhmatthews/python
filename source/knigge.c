
/*

                                       Space Telescope Science Institute

Synopsis:
	These are all of the routines necessary to define a KWD wind
Arguments:		

Returns:
 
Description:	
Notes:


History:
 	01mar	ksl	Added this possibility
	02jan	ksl	Rechecked all of the routines in an atttempt to
			determine why I do not produce KWD figure 7 well.
        11may   nsh     Added the option to have a multiple of sound speed as initial velocity

 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "atomic.h"
#include "python.h"


/***********************************************************
	Space Telescope Science Institute

 Synopsis:
	get_knigge_wind_params gets input data that is necessary for a Knigge's 
	description of the wind
Arguments:		

Returns:
 
Description:	
	The parameters, geo.kn...,  obtained here are only used in the routines in stellar_winds.c
	which calculate the velocity and density of the wind during the initialization process.
	Other portions of the structure, geo defined here are more general purpose.		
Notes:


History:
 	01mar	ksl	Added this possibility
	04jun	ksl	Moved readding wind_mdot to this routine
        04Sep   SS      Changed finite disk case to cut off wind
                        at the top (rather than bottom) outer edge of disk.
	080518	ksl	60a - Fixed to account for fact that geo.wind_mdot 
			has been initialized in cgs units
**************************************************************/

int
get_knigge_wind_params ()
{
  double kn_wind_mdot_integral ();
  double qromb ();
  double dmin;
  double disktheta, test;

  Log ("Creating Knigge's wind description for Cataclysmic Variable\n");

  geo.wind_mdot /= MSOL / YR;
  rddoub ("wind.mdot(msol/yr)", &geo.wind_mdot);
  geo.wind_mdot *= MSOL / YR;


  geo.kn_r_scale = 7e10;	/*Accleration length scale for wind */
  geo.kn_alpha = 1.5;		/* Accleration scale exponent for wind */
  geo.kn_v_infinity = 3;	/* Final speed of wind in units of escape velocity */
  geo.kn_lambda = 0.0;		/* Mass loss rate exponent */
  geo.kn_dratio = dmin = 0.5 * sqrt (geo.diskrad / geo.rstar);	/* Center of collimation in units of the WD 
								   radius. The value set here is for the minimum collimation, see KWD95.  The coefficient 0.5 is 
								   approximate */
  geo.kn_v_zero = 1.0;		/* NSH 19/04/11 New parameter - multiple of sound speed to be used as initial velocity of wind */
  geo.wind_rho_min=1;
  geo.wind_rho_max=geo.diskrad / geo.rstar;


/* There is confusion in various papers concerning whether to use d or d/dmin.  In KWD95, d/dmin was
used but in later papers, e.g KD97 d in WD radii was used.  I believe d is more natural and so will use it, 
but one should remember that this differs from KWD.  To emphasize this we will calculate and log d/dmin.
Terminolgy is awful here. -- ksl 

As now represented geo.kn_dratio is the distance to the focus point in stellar radii!

*/
  rddoub ("kn.d(in_wd_radii)", &geo.kn_dratio);
  Log_silent ("dmin = %f so the ratio d/dmin here is %f  (%.2e %.2e) \n",
	      dmin, geo.kn_dratio / dmin, geo.diskrad, geo.rstar);


  rddoub ("kn.mdot_r_exponent", &geo.kn_lambda);	/* Mass loss rate exponent */
  rddoub ("kn.v_infinity(in_units_of_vescape)", &geo.kn_v_infinity);	/* Final speed of wind in units of escape velocity */
  if (geo.kn_v_infinity < 0)
    {
      Log
	("Since geo.kn_v_infinity is less than zero, will use SV prescription for velocity law.\n Velocity at base remains the soundspeed\n");
    }

  rddoub ("kn.acceleration_length(cm)", &geo.kn_r_scale);	/*Accleration length scale for wind */
  rddoub ("kn.acceleration_exponent", &geo.kn_alpha);	/* Accleration scale exponent */
  rddoub ("kn.v_zero(multiple_of_sound_speed_at_base)", &geo.kn_v_zero);

/* Assign the generic parameters for the wind the generic parameters of the wind. Note that one should not really
use these generic parameters in the rest of the routines here, as especially for the yso model they may have
to be modified -- ksl 04jun */


  /* JM 1502 -- added capability for user to adjust launch radii of KWD wind
     required to reproduce Stuart's X-ray models (Sim+ 2008,2010) 
     units are in stellar radii / WD radii / grav radii as in SV model 
     
     1509 -- ksl -- This had not be implemented quite right.  Now fixed.  Defaults
     are for standard KWD model */

  rddoub ("kn.rmin(in_wd_radii)", &geo.wind_rho_min);
  rddoub ("kn.rmax(in_wd_radii)", &geo.wind_rho_max ); 

  geo.wind_rho_min *= geo.rstar;
  geo.wind_rho_max *= geo.rstar;
  geo.wind_thetamin = atan (1. / geo.kn_dratio);

/* Somewhat paradoxically diskrad is in cm, while dn_ratio which is really d in KWD95 is 
in units of WD radii */
  geo.wind_thetamax = atan (geo.diskrad / (geo.kn_dratio * geo.rstar));

  /* Next lines added by SS Sep 04. Changed the wind shape so that the boundary touches the outer 
     corner of the disk rather than the intersection of the disk edge with the xy-plane. */

  if (geo.disk_type == 2)
    {
      geo.wind_thetamax =
	atan (geo.diskrad /
	      (((geo.kn_dratio * geo.rstar) + zdisk (geo.diskrad))));
    }

  geo.wind_rmin=geo.wind_rho_min;
  geo.wind_rmax=geo.rmax;

  /* The change in the boundary of the wind (as corner of disk -- see above) 
     means that wind_rho_max nees to be redefined so that it is used correctly
     to compute the boundary of the wind elsewhere. */

  if (geo.disk_type == 2) // If disk_type==2, then the disk is vertically extended
    {
      geo.wind_rho_max =
	geo.diskrad - (zdisk (geo.diskrad) * tan (geo.wind_thetamax));
    }

  /* Some changes were made here to make this more general by Nick */
//OLD70d  geo.xlog_scale = geo.rstar/10.0;
//geo.zlog_scale = 1e8;
//OLD70d  geo.zlog_scale = geo.rstar/10.0;
/* 70d - ksl - I changed the scaling to something that produced more cells
 * in the wind, at the cost of slightly less spatial resolution at the inner
 * edge of the wind
 */

  /* if modes.adjust_grid is 1 then we have already adjusted the grid manually */
  if (modes.adjust_grid == 0)
    {
      geo.xlog_scale = geo.rstar;
      geo.zlog_scale = geo.rstar;
    }


/*Now calculate the normalization factor for the wind*/

  test = geo.rstar;

  /* For non-flat disk some streamlines are missing (SS). */

  if (geo.disk_type == 2)
    {
      disktheta = atan (zdisk (geo.diskrad) / geo.diskrad);
      test =
	geo.kn_dratio * geo.rstar * sin (geo.wind_thetamin) *
	cos (disktheta) / sin ((PI / 2.) - geo.wind_thetamin - disktheta);
    }

  geo.mdot_norm = qromb (kn_wind_mdot_integral, test, geo.diskrad, 1e-6);

  return (0);
}


/* The routines in this file define and summarize the properties of the wind.  The routines here are
   specific to the Knigge description of a wind. 
 */

/***********************************************************
	Space Telescope Science Institute

 Synopsis:
	double kn_velocity(x,v) calculates the v in cartesion coordinates
	of a Knigge wind from a position x in cartesian coordinates.  
Arguments:		
	double x[]		the postion where for the which one desires the velocity
Returns:
	double v[]		the calculated velocity
	
	The amplitude of the velocity is returned 
	
Description:	

	This implements equations 17 & 18 of KWD95.  To get the constant
	velocity law set kn_alpha to 0.  

	58 -- If v_inf<0, then an SV law is used for the velocity along
	the streamline.	
		
Notes:

At the inner edge of the accretion disk the temperature of the disk
goes to zero, and so in principle the thermal velocity goes to zero.
The kfudge below keeps this from happening. What did Christian actually
do here?? 

04mar	ksl	The way in which the wind velocity is made to fill
		space, particular the space inside the wind cone,
		is not obvious.  The way in which this is done at
		present causes the div v calulation to be negative
		inside the windcone.  I have suppressed this after
		the first error in wind_div_v, but this could be
		a mistake.
History:
	01mar	ksl	Began coding
	02jan	ksl	Removed the implicit assumption that we always
			wanted the velocity "above" the disk plane.
	04aug	ksl	Modified to account for a disk with vertical
			extent, and to return the .
	05jul	ksl	Corrected error in the direction of ptest
	06oct	ksl	58 -- provided a very simple way to use a 
			SV velocity law (that is a law in which
			the scale height is always the place where
			the velocity is half the peak velocity) in
			what otherwise looks like s KWD model.  The
			SV law is used whenever v_infinity is less 
			than 0. Note that vzero is not changed from
			the KWD prescription.  
		
**************************************************************/

#define kfudge 1.01


double
kn_velocity (x, v)
     double x[], v[];
{
  double r, rzero, theta;
  double ldist, zzz, v_escape, vl;
  double dd;
  double vzero, kn_vzero ();
  double fabs ();
  double tdisk ();
  double teff ();
  struct photon ptest;
  double xtest[3];
  double s;
  double ds_to_disk ();

/*dd is the distance of the focus point along the z axis in cm */

  dd = geo.rstar * geo.kn_dratio;
  r = sqrt (x[0] * x[0] + x[1] * x[1]);	// rho of the point we have been given

/* rzero is the position where the stream line  hits the xy plane. This could less than R_wd. */

  rzero = r / (1. + fabs (x[2] / dd));
  theta = atan (rzero / dd);

/* ldist is the poloidal distance along the streamline */

  ldist = sqrt ((r - rzero) * (r - rzero) + x[2] * x[2]);

/* 04aug -- ksl -- 52 Take the thickness of the disk into account if that is necessary.  */

  if (geo.disk_type == 2)
    {
      xtest[0] = r;		// Define xtest in the +xz plane
      xtest[1] = 0;
      xtest[2] = fabs (x[2]);
      ptest.x[0] = rzero;	// Define ptest at the intersection of the streamline and x axis
      ptest.x[1] = 0.0;
      ptest.x[2] = EPSILON;
      ptest.lmn[0] = sin (theta);	//56d -- ksl -- lmn is along the streamline toward xtest
      ptest.lmn[1] = 0.0;
      ptest.lmn[2] = cos (theta);
      s = ds_to_disk (&ptest, 1);
      move_phot (&ptest, s);	// Now test photon is at disk surface
      vsub (ptest.x, xtest, xtest);
      ldist = length (xtest);
      rzero = length (ptest.x);
    }

/* Note, Knigge's velocity law is slightly different from SV */

/* If statement limits the poloidal component of the velocity to be no more than 
the poloidal velocity at the inner edge of the wind. It is there for continuity reasons */

  if (rzero < geo.rstar)
    v_escape = sqrt (2. * G * geo.mstar / geo.rstar);
  else
    v_escape = sqrt (2. * G * geo.mstar / rzero);

  vzero = kn_vzero (rzero);
  vzero *= geo.kn_v_zero;	/*NSH 19/4/11 Multiply sound speed by new multiplyer to allow faster velocity at base of wind */

/* 578 -- 06oct -- ksl -- The next lines are modified to allow one to create a SV style
velocity law if kn_v_infinity is less than 0 */

  if (geo.kn_v_infinity > 0.0)
    {
      zzz = ldist / (ldist + geo.kn_r_scale);
      zzz = pow (zzz, geo.kn_alpha);	// In Knigge's notation this is actually beta
      vl = vzero + (geo.kn_v_infinity * v_escape - vzero) * zzz;
    }
  else
    {
      zzz = pow (ldist / geo.kn_r_scale, geo.kn_alpha);
      vl =
	vzero + ((-geo.kn_v_infinity) * v_escape - vzero) * zzz / (1. + zzz);
    }


  v[0] = vl * sin (theta);

/* Unlike the polidal component we had been treating the WD as a point source for 
the azimuthal velocity.  One could do this more consistently, with the following
code.: 
  if (r > 0)
    v[1] = sqrt (G * geo.mstar * rzero) / r;  // Eqn 8 KWD
  else
    v[1] = 0;

As far as I could tell it made no difference  02apr-ksl  */

  if (rzero > geo.rstar)
    v[1] = sqrt (G * geo.mstar * rzero) / r;	// Eqn 8 KWD
  else if (r > 0)
    v[1] = sqrt (G * geo.mstar * geo.rstar) / r;
  else
    v[1] = 0;
  v[2] = vl * cos (theta);
  if (x[2] < 0)
    v[2] *= (-1);

/* 

04aug -- ksl --52 At this point we have calculated the velocity in the xz plane, which
is identical to the statement that we have calculated it in cylindrical coordinates.  
We now project back to xyz coordinates if necessary.  

06oct -- ksl -- In practice, the current version of python currently should never need to carry out  
this projection.  However, it is needed, if we were to move to full 3d calculation and for certain 
test programs.

*/
  if (x[1] != 0.0)
    {
      project_from_cyl_xyz (x, v, xtest);
      stuff_v (xtest, v);
    }

  return (sqrt (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]));
}

/***********************************************************
	Space Telescope Science Institute

 Synopsis:
	double kn_rho(x) calculates the density of an kn_wind at a position x
Arguments:		
	double x[]	the position where for the which one desires the density
Returns:
	The density at x is returned in gram/cm**3
	
Description:	
		
Notes:

History:
 	01mar	ksl	Began work.
 
**************************************************************/

double
kn_rho (x)
     double x[];
{
  double r, rzero;
  double dd;
  double vqr, vzero;
  double v[3], rho;
  double kn_velocity ();
  double kn_vzero ();
  double kn_rho_zero ();
  struct photon ptest;
  double s, theta;
  double ds_to_disk ();

  dd = geo.rstar * geo.kn_dratio;
  r = sqrt (x[0] * x[0] + x[1] * x[1]);	//rho coordinate of the point we have been given
  rzero = r / (1. + fabs (x[2] / dd));	//rho at the base for this streamline
  /* If the disk is thick we need to modify the position of rzero */

  if (geo.disk_type == 2)
    {
      theta = atan (rzero / dd);
      ptest.x[0] = rzero;
      ptest.x[1] = 0.0;
      ptest.x[2] = EPSILON;
      ptest.lmn[0] = cos (theta);
      ptest.lmn[1] = 0.0;
      ptest.lmn[2] = sin (theta);
      s = ds_to_disk (&ptest, 1);
      move_phot (&ptest, s);	// Now test photon is at disk surface
      rzero = sqrt (ptest.x[0] * ptest.x[0] + ptest.x[1] * ptest.x[1]);
    }
  kn_velocity (x, v);
  vqr = sqrt (v[0] * v[0] + v[2] * v[2]);	//poloidal velocity
  vzero = kn_vzero (rzero);	//polidal velocity at base
  if ((r * vqr) > 0)		// If statement to assure denominator is not zero
    {
      rho = kn_rho_zero (rzero) * (rzero * rzero * vzero) / (r * r * vqr);
    }
  else
    rho = 0.0;
  return (rho);
}



/***********************************************************
	Space Telescope Science Institute

 Synopsis:
	Calculate the sound speed (launch speed) of the wind at r
Arguments:		
 
Returns:
 
Description:	
	For KWD, the launch speed is the sound velocity
	
		
Notes:

History:
 	01mar      ksl	Coding began.
 
**************************************************************/


double
kn_vzero (r)
     double r;
{
  double tref, tdisk ();
  double t, teff ();
  double ratio, v;
  tref = tdisk (geo.mstar, geo.disk_mdot, geo.rstar);
  ratio = r / geo.rstar;
  if (ratio < kfudge)
    ratio = kfudge;
  t = teff (tref, ratio);
  v = 1.e6 * sqrt (t / 1e4);	// Frank, King & Raine 1985
  return (v);
}


/***********************************************************
	Space Telescope Science Institute

 Synopsis:
	double kn_wind_mdot_integral(r) is the integrand of SV model 
	for mdot as a function of radius
Arguments:		
	double r;
 
Returns:
 
Description:	
	
		
Notes:
       There is an extra factor of 2 because the wind emerges 
       from both sides of the disk

History:
 	01mar      ksl	Coding began.
 
**************************************************************/

double
kn_wind_mdot_integral (r)
     double r;
{
  double tref, tdisk ();
  double t, teff ();
  double x, ratio;
  tref = tdisk (geo.mstar, geo.disk_mdot, geo.rstar);
  ratio = r / geo.rstar;
  if (ratio < kfudge)
    ratio = kfudge;
  t = teff (tref, ratio);
  x = 4. * PI * pow (t, 4. * geo.kn_lambda) * r;
  return (x);
}


/***********************************************************
	Space Telescope Science Institute

 Synopsis:
	Calculate rho at the base of the wind
Arguments:		
 
Returns:
 
Description:	
	
		
Notes:

History:
 	01mar      ksl	Coding began.
 
**************************************************************/

double
kn_rho_zero (r)
     double r;
{
  double tref, tdisk ();
  double t, teff ();
  double x, ratio;
  double kn_vzero (), vzero, dd, cosd;
  tref = tdisk (geo.mstar, geo.disk_mdot, geo.rstar);
  ratio = r / geo.rstar;
  if (ratio < kfudge)
    ratio = kfudge;
  t = teff (tref, ratio);
  x = geo.wind_mdot * pow (t, 4. * geo.kn_lambda) / geo.mdot_norm;
  vzero = kn_vzero (r);
  dd = geo.rstar * geo.kn_dratio;
  cosd = dd / sqrt (r * r + dd * dd);
  x /= (vzero * cosd);
  return (x);
}
