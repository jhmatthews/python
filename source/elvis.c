
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "atomic.h"
#include "python.h"

/***********************************************************
                     Space Telescope Science Institute

 Synopsis:
	These are the routines that are specific to an offiset wind model 
	intended to mimic a proposed model for AGN winds.  It is very similar
	to an SV wind except for an offset
Arguments:		

Returns:
 
Description:	

Notes:

	These routines were all developed from the original sv routines by SS and so 
	changes in one may well effect the other.  The velocity flow assumes that one
	calculates a "poloidal distance" starting with a vertical flow until you get
	to the offset height at which point the remaining distance is calculated according
	to that prescribed by the SV model 

History:
         Oct06 SS
	 1111	ksl	Began to try to understand in detail what Stuart intended prior to any
	 		real code modifications.  
**************************************************************/


/***********************************************************
                                       Space Telescope Science Institute

 Synopsis:
	get_elvis_wind_params gets input data which is necessary for a Shlossman & Vitello 
	description of the wind with a vertical offset as requested for the Elvis AGN model
	This is all copied almost directly from the equivalent SV routines.
Arguments:		

Returns:
 
Description:	
	The parameters, geo.sv...,  needed to define the wind.
	The parameters are identical to that nneed for a SV wind with an additional offset.
Notes:


History:
         Oct06 SS
	 1111	ksl	Began to try to understand in detail what Stuart intended prior to any
	 		real code modifications
**************************************************************/

int 
get_elvis_wind_params (void)
{
  double windmin, windmax, theta_min, theta_max;
  double qromb (), elvis_wind_mdot_integral ();

  Log ("Creating an SV/Elvis wind model for an AGN\n");

  rddoub ("wind.mdot(msol/yr)", &geo.wind_mdot);
  geo.wind_mdot *= MSOL / YR;


  geo.sv_rmin = 2.8e9;
  geo.sv_rmax = 8.4e9;
  geo.sv_thetamin = 20. / RADIAN;
  geo.sv_thetamax = 65. / RADIAN;
  geo.sv_gamma = 1.;
  geo.sv_v_zero = 6e5;		/* velocity at base of wind */
  geo.sv_r_scale = 7e10;	/*Accleration length scale for wind */
  geo.sv_alpha = 1.5;		/* Accleration scale exponent */
  geo.sv_v_infinity = 3;	/* Final speed of wind in units of escape velocity */
  geo.sv_lambda = 0.0;		/* Mass loss rate exponent */
  geo.elvis_offset = 1.e14;

  windmin = geo.sv_rmin / geo.rstar;
  windmax = geo.sv_rmax / geo.rstar;
  rddoub ("sv.diskmin(wd_rad)", &windmin);
  rddoub ("sv.diskmax(wd_rad)", &windmax);


  geo.sv_rmin = windmin * geo.rstar;
  geo.sv_rmax = windmax * geo.rstar;

  theta_min = geo.sv_thetamin * RADIAN;
  theta_max = geo.sv_thetamax * RADIAN;
  rddoub ("sv.thetamin(deg)", &theta_min);
  rddoub ("sv.thetamax(deg)", &theta_max);
  geo.sv_thetamin = theta_min / RADIAN;
  geo.sv_thetamax = theta_max / RADIAN;

  rddoub ("sv.mdot_r_exponent", &geo.sv_lambda);	/* Mass loss rate exponent */
  rddoub ("sv.v_infinity(in_units_of_vescape", &geo.sv_v_infinity);	/* Final speed of wind in units of escape velocity */

  rddoub ("sv.acceleration_length(cm)", &geo.sv_r_scale);	/*Accleration length scale for wind */
  rddoub ("sv.acceleration_exponent", &geo.sv_alpha);	/* Accleration scale exponent */

  rddoub ("elvis_offset(cm)", &geo.elvis_offset);	/* This is the vertical offset - the length over which the wind rises vertically */

/* Assign the generic parameters for the wind the generic parameters of the wind */


  geo.wind_rmin = geo.rstar;
  geo.wind_rmax = geo.rmax;
//Old  geo.wind_rho_min = geo.sv_rmin;
//Old  geo.wind_rho_max = geo.sv_rmax;
  geo.wind_rho_min = geo.sv_rmin - (geo.elvis_offset * tan (geo.sv_thetamin));
  geo.wind_rho_max = geo.sv_rmax - (geo.elvis_offset * tan (geo.sv_thetamin));
  geo.wind_thetamin = geo.sv_thetamin;
  geo.wind_thetamax = geo.sv_thetamax;
//OLD  geo.xlog_scale = geo.sv_rmin + (geo.elvis_offset * tan (geo.sv_thetamin));

  /* if modes.adjust_grid is 1 then we have already adjusted the grid manually */
  if (modes.adjust_grid == 0)
    {
      geo.xlog_scale = geo.sv_rmin;
      geo.zlog_scale = 1e15;	/* Big number - for AGN */
    }

/*Now calculate the normalization factor for the wind*/

  geo.mdot_norm =
    qromb (elvis_wind_mdot_integral, geo.sv_rmin, geo.sv_rmax, 1e-6);
  return (0);
}

/***********************************************************
                                       Space Telescope Science Institute

 Synopsis:
	double elvis_velocity(x,v) calulates the v of an offset Schlossman Vitello wind from a position
	x
Arguments:		
	double x[]		the position where for the which one desires the velocity
Returns:
	double v[]		the calculated velocity
	
	The amplitude of the velocity is returned 
	
Description:	
		
Notes:

History:
 
**************************************************************/

double 
elvis_velocity (double x[], double v[])
{
  double r, rzero, theta, speed;
  double ldist, zzz, v_escape, vl;
  double elvis_find_wind_rzero ();
  double elvis_theta_wind ();
  struct photon ptest;
  double xtest[3];
  double s;
  double ds_to_disk ();

  zzz = v_escape = -99.;

  rzero = elvis_find_wind_rzero (x);
  theta = elvis_theta_wind (rzero);

  /* Calculate the poloidal distance, assuming a stream line that is vertical 
   * up to the elvis offset and then along an SV stream line after that.  ksl 111123 
   */

  r = sqrt (x[0] * x[0] + x[1] * x[1]);
  if (fabs (x[2]) > geo.elvis_offset)
    {
      ldist =
	geo.elvis_offset + sqrt ((r - rzero) * (r - rzero) +
				 (x[2] - geo.elvis_offset) * (x[2] -
							      geo.elvis_offset));
    }
  else
    {
      ldist = fabs (x[2]);
    }

  /* If the disk is vertically extended ksl 111124 
   * ERROR? - This calculation of the poloical distance replaces the one above.  It does not take the offset into account*/
  if (geo.disk_type == 2)
    {
      xtest[0] = r;		// Define xtest in the +z plane
      xtest[1] = 0;
      xtest[2] = fabs (x[2]);
      ptest.x[0] = rzero;
      ptest.x[1] = 0.0;
      ptest.x[2] = EPSILON;
      ptest.lmn[0] = sin (theta);	// 56d -- ptest direction is along stream line
      ptest.lmn[1] = 0.0;
      ptest.lmn[2] = cos (theta);
      s = ds_to_disk (&ptest, 1);
      move_phot (&ptest, s);	// Now test photon is at disk surface
      vsub (ptest.x, xtest, xtest);
      ldist = length (x);
    }


  /* Having calculated the "poloidal distance" calculate the velocity ksl 111124 */
  vl = geo.sv_v_zero;
  if (ldist > 0)
    {
      zzz = pow (ldist / geo.sv_r_scale, geo.sv_alpha);

      if (rzero < geo.rstar)
	v_escape = sqrt (2. * G * geo.mstar / geo.rstar);
      else
	v_escape = sqrt (2. * G * geo.mstar / rzero);

      vl =
	geo.sv_v_zero + (geo.sv_v_infinity * v_escape -
			 geo.sv_v_zero) * zzz / (1. + zzz);
    }

  /* Now calculate the direction of the velocity */
  if (fabs (x[2]) > geo.elvis_offset)
    {
      v[0] = vl * sin (theta);
      v[2] = vl * cos (theta);
    }
  else
    {
      v[0] = 0;
      v[2] = vl;
    }


  if (r > 0)
    v[1] = sqrt (G * geo.mstar * rzero) / r;
  else
    v[1] = 0;

  if (x[2] < 0)
    v[2] *= (-1);
  /* 04aug -- ksl --52 At this point we have calculated the velocity in the xz plane, which
   * is identical to the statement that we have calculated it in
   * cylindrical coordinates.  The next bit projects back to xyz 
   * coordinates if x was not originally in the xz plane.
   */
  if (x[1] != 0.0)
    {
      project_from_cyl_xyz (x, v, xtest);
      stuff_v (xtest, v);
    }


  speed = (sqrt (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]));
  if (sane_check (speed))
    {
      Error ("elvis_velocity:sane_check x %f %f %f v %f %f %f\n", x[0], x[1],
	     x[2], v[0], v[1], v[2]);
      Error
	("elvis_velocity: rzero %f theta %f ldist %f zzz %f v_escape %f vl %f\n",
	 rzero, theta, ldist, zzz, v_escape, vl);
    }


  return (speed);

}

/***********************************************************
                                       Space Telescope Science Institute

 Synopsis:
	double elvis_rho(x) calculates the density of an offset sv_wind at a position x
Arguments:		
	double x[]	the position where for the which one desires the denisty
Returns:
	The density at x is returned in gram/cm**3
	
Description:	
		
Notes:

History:
 
**************************************************************/

double 
elvis_rho (double x[])
{
  double r, rzero, theta;
  double ldist;
  double elvis_find_wind_rzero ();
  double elvis_theta_wind ();
  double dmdot_da;
  double dtheta_drzero, dr_drzero;

  double v[3], rho;
  double elvis_velocity ();
//  double elvis_find_wind_rzero (), elvis_theta_wind ();
  struct photon ptest;
  double xtest[3];
  double s;
  double ds_to_disk ();


  elvis_velocity (x, v);

  rzero = elvis_find_wind_rzero (x);

  /* ERROR - This is a mistake, and is a root cause of some of the problems we are seeing in the Elvis wind 111124 */
  if (rzero > geo.sv_rmax)
    {
      rho = 1.e-20;		//This is a fudge which keeps the wind boundaries conical by filling in excess space with
      // "empty" space - a waste of memory and time but to do otherwise would require a 
      // more substantial change to the way wind boundaries are defined SSOct06
      return (rho);
    }
  if (rzero < (geo.sv_rmin + (geo.elvis_offset * tan (geo.sv_thetamin))))
    {
      rho = 1.e-20;		//fudge, as above
      return (rho);
    }




  theta = elvis_theta_wind (rzero);

  r = sqrt (x[0] * x[0] + x[1] * x[1]);
  if (fabs (x[2]) > geo.elvis_offset)
    {
      ldist =
	geo.elvis_offset + sqrt ((r - rzero) * (r - rzero) +
				 (x[2] - geo.elvis_offset) * (x[2] -
							      geo.elvis_offset));
    }
  else
    {
      ldist = fabs (x[2]);
    }

  if (geo.disk_type == 2)
    {
      xtest[0] = r;		// Define xtest in the +z plane
      xtest[1] = 0;
      xtest[2] = fabs (x[2]);
      ptest.x[0] = rzero;
      ptest.x[1] = 0.0;
      ptest.x[2] = EPSILON;
      ptest.lmn[0] = cos (theta);
      ptest.lmn[1] = 0.0;
      ptest.lmn[2] = sin (theta);
      s = ds_to_disk (&ptest, 1);
      move_phot (&ptest, s);	// Now test photon is at disk surface
      vsub (ptest.x, xtest, xtest);
      ldist = length (xtest);
      rzero = length (ptest.x);
    }

  if (fabs (x[2]) < geo.elvis_offset)
    {
      dmdot_da =
	geo.wind_mdot * pow (rzero, geo.sv_lambda) / geo.mdot_norm / 2.;
      rho = dmdot_da / v[2];

      return (rho);
    }

  /* Reduced by a factor of 2 to account for radiation on both sides of the disk */
  dmdot_da =
    geo.wind_mdot * pow (rzero,
			 geo.sv_lambda) * cos (theta) / geo.mdot_norm / 2.;

  /* Although the current definition of sv_theta_wind is continuous, the derivative is not continuous accross the
     outer boundary of the wind and thus dtheta_drzero would appear to change discontinuously.   This created
     large apparent density jumps at the outside edge of the wind. We can't allow that and so we force the derivative to equal
     that at the edge of the wind if you try to calculate the density rho.  ksl 97apr23 */

  if (rzero > geo.sv_rmax)
    rzero = geo.sv_rmax;
  dtheta_drzero =
    (elvis_theta_wind (rzero) -
     elvis_theta_wind ((1. - EPSILON) * rzero)) / (EPSILON * rzero);

  dr_drzero = 1. + ldist * dtheta_drzero / cos (theta);
  /* Note VS93 eqn 8 is dr/drzero but equation  7 is drzero/dr   ksl 97 apr 19 */
  rho = rzero * dmdot_da / (dr_drzero * r * v[2]);

  return (rho);
}


/***********************************************************
                                       Space Telescope Science Institute

 Synopsis:
	elvis_find_wind_rzero(p) locates the radial position in the disk from which the 
	stream line arises for an offset SV wind (the "elvis" model).  
 
Arguments:		
	double p[]	    A 3 vector position, and not a photon structure
Returns:
	elvis_find_wind_rzero returns the radius of the footpoint of the stream line
 
Description:	
		
Notes:

History:
 
**************************************************************/


double 
elvis_find_wind_rzero (
    double p[]		/* Note that p is a 3 vector and not a photon structure */
)
{
  double x, z;
  double elvis_zero_r ();
  double zbrent ();
  double rho_min, rho_max, rho;
  int elvis_zero_init ();

  /* thetamin and theta max are defined w.r.t  z axis */


  z = fabs (p[2]);		/* This is necessary to get correct answer above
				   and below plane */

  /* If the position is underneath the elvis_offset, the stream line is assumed
   * to be vertical and so you just return rho 111124 ksl */
  if (z <= geo.elvis_offset)
    {
      x = (sqrt (p[0] * p[0] + p[1] * p[1]));	// If p is in the xy plane, there is no need to think further
      return (x);
    }


  elvis_zero_init (p);		/* This initializes the routine sv_zero_r.  It is not
				   actually needed unless zbrent is called, but it
				   does allow you to check your answer otherwize
				 */

  /* The next lines provide a value for the footpoint position when the positions
   * we want the footpoint for is outside of the wind 111124 ksl */
  rho_min = geo.sv_rmin + z * tan (geo.sv_thetamin);
  rho_max = geo.sv_rmax + z * tan (geo.sv_thetamax);
  rho = sqrt (p[0] * p[0] + p[1] * p[1]);

  /* Note that theta_min and theta_max are measured from the z axis */

  if (rho <= rho_min)
    {
      x = geo.sv_rmin * rho / rho_min;
      return (x + (geo.elvis_offset * tan (geo.sv_thetamin)));
    }
  if (rho >= rho_max)
    {
      x = geo.sv_rmax + rho - rho_max;
      return (x + (geo.elvis_offset * tan (geo.sv_thetamax)));
    }

  /* 100 here means that zbrent will end if one has a guess of rzero which is
     correct ot 100 cm */

  x =
    zbrent (elvis_zero_r,
	    geo.sv_rmin + (geo.elvis_offset * tan (geo.sv_thetamin)),
	    geo.sv_rmax + (geo.elvis_offset * tan (geo.sv_thetamax)), 100.);
  return (x);

}

/***********************************************************
                                       Space Telescope Science Institute

 Synopsis:
	The next two subroutines are used to guess the position in the disk
	from which a streamline arises.  They work in the same manner as sv_zero_r 
	and sv_zero_init.
 
Arguments:		

Returns:
 
Description:	
	
		
Notes:

History:
 
**************************************************************/

double zero_p[3];

int 
elvis_zero_init (double p[])
{
  stuff_v (p, zero_p);
  zero_p[2] = fabs (zero_p[2]);	/* Required to get correct 
				   answer below (in -z ) the disk */
  return (0);
}

/* This routine is used to test whether a guess of r_zero is correct.  If
   you have the answer exactly then elvis_zero_r will return 0 */

double 
elvis_zero_r (double r)
{
  double theta;
  double rho, rho_guess;
  double elvis_theta_wind ();

  theta = elvis_theta_wind (r);

  rho = sqrt (zero_p[0] * zero_p[0] + zero_p[1] * zero_p[1]);
  if (zero_p[2] < geo.elvis_offset)
    {
      rho_guess = r;
    }
  else
    {
      rho_guess = r + tan (theta) * (zero_p[2] - geo.elvis_offset);
    }
  return (rho_guess - rho);


}

/***********************************************************
                                       Space Telescope Science Institute

 Synopsis:
	double elvis_theta_wind(r) finds the angle at which the wind that rises from a specific point on the disk bend to
 
Arguments:
	double r	a radial distance on the surface of the disk		

Returns:
 
Description:	
		
Notes:

History:
 
**************************************************************/

double 
elvis_theta_wind (double r)
{
  double theta;
  if (r <= (geo.sv_rmin + geo.elvis_offset * tan (geo.sv_thetamin)))
    return (atan (tan (geo.sv_thetamin * r / geo.sv_rmin)));
  if (r >= (geo.sv_rmax + geo.elvis_offset * tan (geo.sv_thetamax)))
    return (geo.sv_thetamax);
  theta = geo.sv_thetamin +
    (geo.sv_thetamax -
     geo.sv_thetamin) * pow ((r - geo.sv_rmin -
			      geo.elvis_offset * tan (geo.sv_thetamin)) /
			     (geo.sv_rmax +
			      geo.elvis_offset * tan (geo.sv_thetamax) -
			      geo.sv_rmin -
			      geo.elvis_offset * tan (geo.sv_thetamin)),
			     geo.sv_gamma);
  return (theta);

}

/***********************************************************
                                       Space Telescope Science Institute

 Synopsis:
	double elvis_wind_mdot_integral(r) is the integrand of SV model for mdot as a function
	of radius - works in the same way as the SV routine but calls the elvis version of the theta routines. 
Arguments:		
	double r;
 
Returns:
 
Description:	
	
		
Notes:

History:

 
**************************************************************/

double 
elvis_wind_mdot_integral (double r)
{
  double x;
  double elvis_theta_wind ();

  if (r < (geo.sv_rmin + geo.elvis_offset * tan (geo.sv_thetamin)))
    {
      x = 0.0;
      return (x);
    }
  if (r > geo.sv_rmax)
    {
      x = 0.0;
      return (x);
    }


  x = 2 * PI * pow (r, geo.sv_lambda + 1.) * cos (elvis_theta_wind (r));
  return (x);

}



/***********************************************************
	Space Telescope Science Institute

 Synopsis:
 	ds_to_pillbox calculates the distance to a pillbox, really
	an annular region with an inner radius and outer radius
	and a height.  The pillbox is symmetric with respect to the xy
Arguments:		
	pp	photon ptr containing a position and direction
	rmin	The minimum rho for the pillbox
	rmax	The maximum rho for the pillbox
	height 	Height of the pillbox
 
Returns:
 
Description:	
	
		
Notes:
	No allowance is made for a geometrically thick disk.  This is
	a purely geometrical construction.

	This is modeled on ds_to_torus, but here all the parameters
	are passed.  It's likely that this routine should replace
	that one

	This should almost surely be moved to phot_util

History:
	11nov	ksl	Began coding to deal with issues that
			currently exist with the elvis model

 
**************************************************************/

double 
ds_to_pillbox (PhotPtr pp, double rmin, double rmax, double height)
{

  struct photon ptest;
  double ds, ds_best, x;
  struct plane xplane;


  ds_best = VERY_BIG;

  /* Make sure we don't mess with pp */
  stuff_phot (pp, &ptest);
  ds = ds_to_cylinder (rmin, &ptest);

  /* Calculate the distance to the innner cylinder */
  if (ds < VERY_BIG)
    {
      /* Check whether we encounted the
       * part of the cylinder we are interested in
       */
      move_phot (&ptest, ds);
      if (fabs (ptest.x[2]) < height)
	{
	  ds_best = ds;
	}
      /* Now reinitialize ptest */
      stuff_phot (pp, &ptest);
    }

  /* Similarly calculate the distance to the outer
   * cylinder
   */
  ds = ds_to_cylinder (rmax, &ptest);
  if (ds < ds_best)
    {
      move_phot (&ptest, ds);
      if (fabs (ptest.x[2]) < height)
	{
	  ds_best = ds;
	}
      stuff_phot (pp, &ptest);
    }

  /* At this point we know whether the photon has interecepted
   * the wall of the cylinder, but we do not know if it intercepted
   * the top or bottom of the cylinder earlier
   */

  xplane.x[0] = 0;
  xplane.x[1] = 0;
  xplane.x[2] = height;
  xplane.lmn[0] = 0.0;
  xplane.lmn[1] = 0.0;
  xplane.lmn[2] = 1.0;

  ds = ds_to_plane (&xplane, &ptest);
  // Note that ds to plane can return a negative number
  if (ds > 0 && ds < ds_best)
    {
      move_phot (&ptest, ds);
      x = fabs (ptest.x[0]);
      if (rmin < x && x < rmax)
	{
	  ds_best = ds;
	}
      stuff_phot (pp, &ptest);
    }

  xplane.x[2] = (-height);

  ds = ds_to_plane (&xplane, &ptest);
  if (ds > 0 && ds < ds_best)
    {
      move_phot (&ptest, ds);
      x = fabs (ptest.x[0]);
      if (rmin < x && x < rmax)
	{
	  ds_best = ds;
	}
      stuff_phot (pp, &ptest);
    }

  return (ds_best);
}
