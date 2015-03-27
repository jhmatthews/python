#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <string.h>
#include "atomic.h"
#include "python.h"

/***********************************************************
                                       Space Telescope Science Institute

 Synopsis:
        double
        dvwind_ds(w,p) finds the gradient velocity vector v, e.g. dvds for the wind 
        coordinates at the position of the photon p.

 Arguments:
        WindPtr w;
        PhotPtr p;

Returns:
        dvds    on successful completion
                the value of where in grid if the photon outside the wind grid.

Description:

Notes:

History:
    01dec	ksl	Coding began, as a result of finding 
			features in output spectra that were 
			due to a discontinuous dv/ds.
	02jan	ksl	Fixed problem in which photons below the
			disk plane were not being dealt with 
			consistently to those above the plane
	04aug	ksl	52a -- Revised to allow for multiple
			coordinate systems and to use the same
			fractional position routines as most
			other portions of python.  Also rewrote
			so more nearly resembles those programs.
	05apr	ksl	55d -- Modified to accept new version
			of coord_fraction, which is much more
			coordinate system independent
	06may	ksl	57+ -- There is no point in transferring
			entire wind to dvwind_ds when we have
			wmain
    1411    JM -- use on the fly method for spherical coordinates
            to avoid problems. This should be improved. See notes
            below and #118 

**************************************************************/

double
dvwind_ds (p)
     PhotPtr p;
{
  double v_grad[3][3];
  double lmn[3], dvel_ds[3], dvds;
  int j, k, n, nn;
  double dot_tensor_vec ();
  struct photon pp;
  int nnn[4], nelem;		// At present the largest number of dimenssion in the grid is 2
  double frac[4];
  double x;


  /* We want the change in velocity along the line of sight, but we
  need to be careful because of the fact that we have elected to 
  combine the upper and lower hemispheres in the wind array.  Since
  we are only concerned with the the scalar dv_ds, the safest thing
  to do is to create a new photon that is only in the upper hemisphere 
  02jan ksl */

  stuff_phot (p, &pp);
  if (pp.x[2] < 0.0)
    {				/*move the photon to the northen hemisphere */
      pp.x[2] = -pp.x[2];
      pp.lmn[2] = -pp.lmn[2];
    }
  
  /* JM 1411 -- ideally, we want to do an interpolation on v_grad here. However, 
  the interpolation was incorrect in spherical coordinates (see issue #118).
  For the moment, I've adopted an on the fly method for spherical coordinates.
  This should be improved by figuring out out how the velocity gradient 
  tensor ought to be rotated in order to give the right answer for spherical 
  coordinates */

  if (geo.coord_type == SPHERICAL)
    {
      struct photon pnew;
      double v1[3], v2[3], diff[3];
      double ds;

      /* choose a small distance which is dependent on the cell size */	
      vsub (pp.x, wmain[pp.grid].x, diff);
      ds = 0.001 * length(diff);

      /* calculate the velocity at the position of the photon */
      /* note we use model velocity, which could potentially be slow,
      but avoids interpolating (see #118) */
      model_velocity(pp.x, v1);

      /* copy the photon and move it by ds, and evaluate the velocity
      at the new point */
      stuff_phot(&pp, &pnew);
      move_phot(&pnew, ds);
      model_velocity(pnew.x, v2);

      /* calculate the relevant gradient */
      dvds = fabs(dot(v1, pp.lmn) - dot(v2, pp.lmn)) / ds;
    }

  else // for non spherical coords we interpolate on v_grad
    {

      n = coord_fraction (0, pp.x, nnn, frac, &nelem);


      for (j = 0; j < 3; j++)
        {
          for (k = 0; k < 3; k++)
	  {
	      x = 0;
	      for (nn = 0; nn < nelem; nn++)
	        x += wmain[nnn[nn]].v_grad[j][k] * frac[nn];

	      v_grad[j][k] = x;

	  }
        }

      // ??? Not clear this is correct !!!  for multiple coordinate systems 
      /* v_grad is in cylindrical cordinates, or more precisely intended
      to be azimuthally symmetric.  One could either 
      (a) convert v_grad to cartesian coordinates at the position of the
      photon or (b) 
      convert the photon direction to cylindrical coordinates 
      at the position of the photon. Possibility b is more straightforward
      given the existing code*/

      project_from_xyz_cyl (pp.x, pp.lmn, lmn);

      dvds = dot_tensor_vec (v_grad, lmn, dvel_ds);

      /* Note that dvel_ds is also in cylindrical coordinates should it ever be
      needed.

      It is not clear that vwind_xyz should not be modified simply to return v
      along the line of sight of the photon as well...since vwind the results are vwind
      are immediately dotted with the photon direction in resonate which is the only place
      the routine is used I believe . ksl
      */
    }

  if (sane_check (dvds))
    {
      Error ("dvwind_ds: sane_check %f\n", dvds);
    }

  return (dvds);

}




/* Calculate the average tau assuming an oscillator strength of 1

??? Not clear to me that this and everything should not be done with gradv above 

   History:
	98	ksl	Coded
	98dec	ksl	Modified so it calulates a length ds which is
			appropriate to the cell size
	02feb	ksl	Added code to find the direction of maximum velocity
			gradient, and also to assure that in calculating dvds_ave
			the pseudo photon pp did not cross the plane of the
			disk.
	04aug	ksl	52a -- Corrected basic error in dvds_ave.  Previous values
			were 10x too large
	05apr	ksl	55d -- Modified to use w[].xcen rather than mid_x and
			mid_z.  This was necessary for rtheta and spherical
			coordinates.  
	05jul	ksl	56d -- Routine was attempting to calcuate vwind_xyz
			outside of the grid on occassion.  This can be a problem, 
			if the result is really wrong.  I have simply tried to
			address this in vwind, which really means vwind always
			needs to produce a plausible result, which in turn
			means that coord_frac must always produce a reasonable
			result.
	05jul	ksl	Added DEBUG statements to limit creation of dvds.diag
			in cases where not needed.  
	06may	ksl	57+ -- Have not changed except to fix call since it is
			claimed that one needs this everywhere.
	0903	ksl	68c -- This routine appears to have been wrotng since
			at least the last time it was changed, although it is
			not clear what the effect of the error was.  pp was not
			intialized in the statements below
 */

#define N_DVDS_AVE	10000
int
dvds_ave ()
{
  struct photon p, pp;
  double v_zero[3], delta[3], vdelta[3], diff[3];
  double sum, dvds, ds;
  double dvds_max, lmn[3];
  int n;
  int icell;
  double dvds_min, lmn_min[3];	//TEST
  char filename[LINELENGTH];

  strcpy (filename, basename);
  strcat (filename, ".dvds.diag");
  if (modes.print_dvds_info)
    {
      optr = fopen (filename, "w");	//TEST
    }

  for (icell = 0; icell < NDIM2; icell++)
    {

      dvds_max = 0.0;		// Set dvds_max to zero for the cell.
      dvds_min = 1.e30;		// TEST


      /* Find the center of the cell */

      stuff_v (wmain[icell].xcen, p.x);

      /* Define a small length */

      vsub (p.x, wmain[icell].x, diff);
      ds = 0.001 * length (diff);

      /* Find the velocity at the center of the cell */
      vwind_xyz (&p, v_zero);

      sum = 0.0;
      for (n = 0; n < N_DVDS_AVE; n++)
	{
	  randvec (delta, ds);
	  if (p.x[2] + delta[2] < 0)
	    {			// Then the new position would punch through the disk
	      delta[0] = (-delta[0]);	// So we reverse the direction of the vector
	      delta[1] = (-delta[1]);
	      delta[2] = (-delta[2]);
	    }
	  vadd (p.x, delta, pp.x);
	  vwind_xyz (&pp, vdelta);
	  vsub (vdelta, v_zero, diff);
	  dvds = length (diff);

	  /* Find the maximum and minimum values of dvds and the direction
	   * for this
	   */

	  if (dvds > dvds_max)
	    {
	      dvds_max = dvds;
	      renorm (delta, 1.0);
	      stuff_v (delta, lmn);
	    }
	  if (dvds < dvds_min)
	    {
	      dvds_min = dvds;
	      renorm (delta, 1.0);
	      stuff_v (delta, lmn_min);
	    }

	  sum += dvds;

	}

      wmain[icell].dvds_ave = sum / (N_DVDS_AVE * ds);	/* 10000 pts and ds */

//Store the maximum and the direction of the maximum
      wmain[icell].dvds_max = dvds_max / ds;
      stuff_v (lmn, wmain[icell].lmn);

      if (modes.print_dvds_info) 
        { 
          fprintf (optr,
	         "%d %8.3e %8.3e %8.3e %8.3e %8.3e %8.3e %8.3e %8.3e %8.3e %8.3e %8.3e %8.3e \n",
	         icell, p.x[0], p.x[1], p.x[2], dvds_max / ds,
	         dvds_min / ds, lmn[0], lmn[1], lmn[2], lmn_min[0],
	         lmn_min[1], lmn_min[2], dot (lmn, lmn_min));
        }

    }


  if (modes.print_dvds_info) 
    fclose (optr);

  return (0);
}
