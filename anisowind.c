#include <stdio.h>
#include <stdlib.h>
#include <math.h>


#include "atomic.h"
#include "python.h"


/* 
randwind (xyz, lmn, north) is a routine to calculate a new direction
for the photon for specific presciptions regarding the anisotropy
of scattering.  To use a different prescription the only thing that
should have to be modified is vrandwind.

Here 

 xyz is the position of the photon in cartesion coordinates
 north is the the direction of maxium velocity gradient in 
	cylindrical coordinates

and 
 lmn is the new direction of the photon in cartesian coordinates 
	which is returned is returned.

02feb	ksl	Modified jumps in randwind so evenly spaced.
02feb	ksl	Mad major modifications to the way in which
		the direction was projected back to the 
		cartesian frame, and also to the generation
		of the photon direction in the azimuthal
		direction.  The routine was seriously 
		screwed up.
02may	ksl	Made modifications to generalize anisowind
		and to incorporate changes associated
		with creating a new (constant source function)
		based scattering function.  The routines are
		now much more integrated in the sense that
		only one routine vrandwind needs to be
		changed to incorporate a new scattering
		model (as long as the model is axisymmetric)
		about some direction..presumable the direction
		of tau_min.  N.B. There are some dependences
		on the order in which the routines are called,
		i.e. you must call randwind before vrandwind,
		since randwind initializes some variables that
		one needs in vrandwind.  
02may	ksl	Changed calls to avoid having to pass tau_min
		to the routine.  Hopefully this will help
		eliminate questions of what has been calculated
		when, or at least isolate it to the routines
		in this file.
02may	ksl	Added a check to prevent tau from going too low
		when calculating pdf_rand.  ?? It might be more
		computationally efficient to switch to isotropic
		in this case ??
*/

#define  TAU_TOP   10.0		/* tau above which we assume the angular distribution function
				   no longer changes significantly, e.g above which the surface behaves
				   like a solid surface */
#define  TAU_BOT    0.01


double tau_randwind = -1000.;
//struct Pdf pdf_randwind;    // Moved into python.h for python_43.2

int
randwind (p, lmn, north)
     PhotPtr p;
     double lmn[3], north[3];
{

  double xyz[3];
  double n;			/* the individual direction cosines in the rotated frame */
  double q;
  double phi;
  int create_basis (), project_from ();
  double vrandwind ();
  double pdf_get_rand ();
  double xlmn[3], dummy[3];
  int project_from_cyl_xyz ();
  double north_xyz[3];
  struct basis nbasis;
  int k;
  double tau;


  if (p->nres < 0)
    {
      Error
	("randwind: Cannot calculate line scattering without nres %d of line\n",
	 p->nres);
      return (-1);
    }

// Next statement gets information needed to recalculate everything
  k = where_in_grid (p->x);
  tau = sobolev (&wmain[k], p, -1., lin_ptr[p->nres], wmain[k].dvds_max);

  stuff_v (p->x, xyz);

  if (tau > TAU_TOP)
    tau = TAU_TOP;		// 
  if (tau < TAU_BOT)
    tau = TAU_BOT;

  if (fabs (tau - tau_randwind) > 0.01)
    {				// (Re)create pdf
      make_pdf_randwind (tau);
      stuff_phot (p, &phot_randwind);
    }

  xlmn[0] = n = pdf_get_rand (pdf_randwind);
  if (sane_check (n))
    {
      Error ("anisowind:sane_check of pdf_get_rand returned %f\n", n);
      xlmn[0] = n = pdf_get_rand (pdf_randwind);
    }

  q = sqrt (1. - n * n);

  phi = 2. * PI * (rand () / MAXRAND);
  xlmn[1] = q * cos (phi);
  xlmn[2] = q * sin (phi);


/* So at this point we have the direction cosines of the photon in a 
   coordinate system in which "north" is along the FIRST axis. "north"
   as presented is a vector in the cylindrical coordinate system.  So
   we need to express "north" in cartesian coordinates.
*/


  project_from_cyl_xyz (xyz, north, north_xyz);

/* Given north in cartesion coordinates, we can now follow the
same procedure as in randvcos, namely, create a basis in which
north is along the FIRST axis.  This is the basis in which we
have calculated the photon direction.  The other two directions
don't really matter, but to create a basis we need to ensure
that the two axes given in create_basis are not parallel to
one another.  */

  cross (north_xyz, z_axis, dummy);
  if (length (dummy) > 0.0)
    create_basis (north_xyz, z_axis, &nbasis);
  else
    create_basis (north_xyz, x_axis, &nbasis);

/* Now project from the coordinate system described by nbasis to
the cartesian frame */

  project_from (&nbasis, xlmn, lmn);

  if (sane_check (lmn[0]) || sane_check (lmn[1]) || sane_check (lmn[2]))
    {
      Error ("Randwind:sane_check NAN problem lmn %f %f %f\n", lmn[0], lmn[1],
	     lmn[2]);
    }

  return (0);
}

/* This is the function that is used to generate the pdf for scattering.  It
basically is calculating a function that is proportional to the probability
density.  dP/dcos(theta) for the photon. 

double vrandwind (x)

     double x   The direction cosine with respect to the direction 
		or tau_min


	02feb14	ksl	Made changes to set a minimum for vrandwind in
			in an attempt to prevent problems with pdf
			generation
*/

#define VRANDWIND_FLOOR 1e-5
#define VRANDWIND_TAU 11.51

double
vrandwind (x)
     double x;
{
//  double a, b;
  double abscos, z;

  abscos = fabs (x);
  if (abscos == 0.0)		// Then at 90 degrees
    {
      return (VRANDWIND_FLOOR);
    }


// The next lines generate a model for scattering in which
// the source function is assumed to be constant through the
// scattering region.  It is also related to escape probs.

  z = (1. - exp (-tau_randwind / abscos)) * abscos;

// These generate the probability density for the Eddington approximation
//  a = 0.5;
//  b = 1.5;
//  z = x * (a * (1. + b * x));


// The next line would generate a completely isotropic distribution
// z=0.1;  // Completely isotropic

  if (sane_check (z))
    {
      Error ("vrandwind:sane_check tau_randwind %f x %f\n", tau_randwind, x);
      return (0.0);
    }


  return (z);


}

/* 
reweightwind calculates the weight of a photon in the wind that
is forced to scatter along a specific line of sight.  

Here 
 p must contain the position and NEW direction for the photon
 being scattered.  In general, p must also be up to date w. r. t.
 the number of the resonance that generated the photon.

The reweighting is returned, but it is also incorporated into the photon.

Notes:  
 xyz is the position of the photon in cartesion coordinates
 lmn is the direction of the photon in cartesian coordinates
 north is the the direction of maxium velocity gradient in 
	cylindrical coordinates (whose basis is referenced to x)

tau_scatter_min, the tau in the direction of the maximum
velocity gradient, for this photon will in general have already
been calculated.

and the reweighting is returned.

02feb	ksl	Coded (or revised) to allow model with the photon
		to be scattered at the center of its resonant
		survace
02may	ksl	Changed significantly to allow for easy changes
		to the scattering function.  The new routine
		is based on the fact that the weight of 
		a photon is essentially the a properly normalized
		probability density
02may	ksl	Changed calling function for reweightwind, in part
		has already been calculated for this position, i.e.
		to better encapsulate this routine.
04deck	ksl	Miniscule mod to make compile cleanly with 03
*/

#define REWEIGHTWIND_TAU_MAX 100.
int reweightwind_init = 1;	//TRUE to start
double reweightwind_zmax;

double
reweightwind (p)
     PhotPtr p;
{
  double ctheta, x, z;
  double north_xyz[3];
  int project_from_cyl_xyz ();
  double xyz[3], lmn[3], north[3], tau;
  double delta[3];
  int k;

  vsub (p->x, phot_randwind.x, delta);
// Idea here is that if photon has moved from position where it was you must recalc.
  if ((x = length (delta)) > DFUDGE)
    {
      k = where_in_grid (p->x);
      tau = sobolev (&wmain[k], p, -1., lin_ptr[p->nres], wmain[k].dvds_max);
      make_pdf_randwind (tau);	// Needed for the normalization
//      Log
//      ("Had to calculate the pdf for this photon %d %e tau_scat_min %10.2f norm %10.2f \n",
//       p->nres, x, tau, pdf_randwind->norm);
      stuff_phot (p, &phot_randwind);
    }
  else
    tau = 0.0;

  stuff_v (p->x, xyz);
  stuff_v (p->lmn, lmn);
  k = where_in_grid (p->x);
  stuff_v (wmain[k].lmn, north);

/* We need to know the cos of the angle between the direction
of maximum velocity gradient and the new direction of
the photon.  north is a vector expressed in terms of
a cylindrical coordinate system (at the point xyz) so one
way to proceed is to project north to xyz coords which is
what we do here */

  project_from_cyl_xyz (xyz, north, north_xyz);

  ctheta = fabs (dot (lmn, north_xyz));


// Factor of 2 needed because the interval is from -1 to 1
// ?? It's definitely needed to make a uniform distribution work
// but is this reason really right
  z = 2. / pdf_randwind->norm;
  x = vrandwind (ctheta) * z;

  if (sane_check (x) || x > 2.0)
    {
      Error
	("Reweightwind:sane_check x %f tau %f ctheta %f z %e pdf_randwind->norm %f\n",
	 x, tau, ctheta, z, pdf_randwind->norm);
      x = 2.0;
    }

  p->w *= x;
  return (x);

}


/*
 *
 * 02june	ksl	Modified make_pdf_randwind so that the first time
 * 			the program is entered an array of cumulative
 * 			distribution functions is created.  This is 
 * 			designed to speed the program up significantly
*/

int init_make_pdf_randwind = 1;
int make_pdf_randwind_njumps;
double make_pdf_randwind_jumps[180];
#define LOGTAUMIN -2.
#define LOGTAUMAX 1.
double pdf_randwind_dlogtau;

int
make_pdf_randwind (tau)
     double tau;
{
  int jj;
  int echeck;
  int pdf_gen_from_func ();
  double vrandwind ();
  double xtau;
  double log10 ();

// Initalize jumps the first time routine is called
  if (init_make_pdf_randwind)
    {
      make_pdf_randwind_njumps = 0;
      for (jj = -88; jj <= 88; jj += 2)
	{
	  make_pdf_randwind_jumps[make_pdf_randwind_njumps] =
	    sin (jj / 57.29578);
	  make_pdf_randwind_njumps++;
	}
      pdf_randwind_dlogtau = (LOGTAUMAX - LOGTAUMIN) / 99.;
      for (jj = 0; jj < 100; jj++)
	{
	  xtau = pow (10., LOGTAUMIN + pdf_randwind_dlogtau * jj);
	  tau_randwind = xtau;	// This is passed to vrandwind by an external variable
	  if ((echeck =
	       pdf_gen_from_func (&pdf_randwind_store[jj], &vrandwind, -1.0,
				  1.0, make_pdf_randwind_njumps,
				  make_pdf_randwind_jumps)) != 0)
//Old ksl 04mar gave warning      &make_pdf_randwind_jumps)) != 0)
	    {
	      Error ("Randwind: return from pdf_gen_from_func %d\n", echeck);
	    }
	}
      init_make_pdf_randwind = 0;
    }

  if (sane_check (tau))
    {
      Error
	("make_pdf_randwind:sane_check Need proper tau (%e) to make pdf_randwind\n",
	 tau);
      tau = 10.;		// Forces something close to isotropic
    }

  jj = (log10 (tau) - LOGTAUMIN) / pdf_randwind_dlogtau + 0.5;

  if (jj < 0)
    jj = 0;
  if (jj > 99)
    jj = 99;

  pdf_randwind = &pdf_randwind_store[jj];
  tau_randwind = tau;		// This is passed to vrandwind by an external variable

  return (0);
}


/***************************************************************
                      
                      University of Southampton

Synopsis:   
  randwind_thermal_trapping is the routine which chooses
  a new anisotropic direction in geo.scatter_mode = 2.
  
Arguments:   

  
Returns:
  0 for success. Also modifies the photon ptr p
  to reflect new direction (p->lmn), and nnscat, which
  should be copied to the phoiton structure after calling
  this routine.
  
Description:  
  This routine uses a rejection method to choose a direction
  so that the probability distribution of directions generated 
  reflects the probability of escape along each direction in accordance
  with the sobolev optical depth. 

Notes:
  
History:
  1406  Moved code here from photo_gen_matom and scatter to avoid 
        duplication


****************************************************************/


int 
randwind_thermal_trapping(p, nnscat)
  PhotPtr p;
  int *nnscat;
{
  double tau_norm, p_norm;
  double tau, dvds, z, ztest;
  double z_prime[3];
  WindPtr one;

  /* find the wind pointer for the photon */
  one = &wmain[p->grid];

  /* we want to normalise our rejection method by the escape 
    probability along the vector of maximum velocity gradient.
    First find the sobolev optical depth along that vector */
  tau_norm = sobolev (one, p, -1.0, lin_ptr[p->nres], one->dvds_max);

  /* then turn into a probability. Note that we take account of
     this in trans_phot before calling extract */
  p_norm = p_escape_from_tau(tau_norm);

  /* Throw error if p_norm is 0 */
  if (p_norm <= 0)
    Error("randwind_thermal_trapping: p_norm is %8.4e in cell %i",
           p_norm, one->nplasma);

  ztest = 1.0;
  z = 0.0;
 
  /* JM 1406 -- we increment nnscat here, and it is recorded in the photon
    structure. This is done because we actuall have to multiply the photon weight 
    by 1/mean escape probability- which is nnscat. this is done in trans_phot.c
    before extract is called. 
  */  
  *nnscat = *nnscat - 1;

   /* rejection method loop, which chooses direction and also calculated nnscat*/
   while (ztest > z)
  {
    *nnscat = *nnscat + 1; //- JM - see above 
    randvec (z_prime, 1.0);       /* Get a new direction for the photon (isotropic */
    stuff_v (z_prime, p->lmn);    // copy to photon pointer


    /* generate random number, normalised by p_norm with a 1.2 for 20% 
       safety net (as dvds_max is worked out with a sample of directions) */
    ztest = (rand () + 0.5) / MAXRAND * p_norm;   
    dvds = dvwind_ds (p);
    tau = sobolev (one, p, -1.0, lin_ptr[p->nres], dvds);

    z = p_escape_from_tau (tau);  /* probability to see if it escapes in that direction */
  }
  
  /* one could copy to the photon pointer here, 
     but for the moment this is done after calling this routine */
  //p->nnscat = nnscat;

  return (0);
}



