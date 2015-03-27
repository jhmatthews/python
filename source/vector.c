/***********************************************************
                                       Space Telescope Science Institute

 Synopsis:

A series of subroutines designed to carry out common 3 vector operations
		
	double dot(a,b)		returns the dot product of the 3 vectors a and b
	double length(x)	returns the length of the 3 vector x
	int renorm(a,scalar)		renormalizes the 3-vector a to a length of
							a; retuns -1 if a has length 0.
	int cross(a,b,c)		c = a x b  (where x is the cross-product)
		
	vmove(u,lmn,s,result)
	double u[],lmn[],s,result[];
					Returns the 3-vec result = u + s * lmn

	vsub(u,v,result)
	double u[],v[],result[];
					Returns the 3 vec result = u -v

	vadd(u,v,result)
	double u[],v[],result[];
					Returns the 3 vec result = u +v

	stuff_v(vin,vout)
	double vin[],vout[];
					Returns the 3-vec  vout=vin.
					

	int project_from_xyz_cyl(a,b,result)
	double a[],b[],result[];
					Return the 3 vector result which is the vector b in cylindrical coords
					at a.   a and b are assumed to have been in cartesian coordinates.  They
					are not necessarily the same since (for example) one might want the
					velocity at some position in a flow.

					
	The next few routines are designed to allow one to handle relatively arbitrary basis manipulations
	

 	int create_basis(u,v,basis_new)
	double u[],v[];
	struct basis *basis_new;		
					creates an orthnormal basis with the first axis in the u direction, the
					y axis in the uv plane, and the the third axis in the perpendicular direction.  
					The system is right-handed 	

	project_from(basis_from, v_in,v_out)
	struct basis *basis_from;     direction cosines to go from rotated to unrotated frame 
	double v_in[],v_out[];       v_in here is in rotated frame, v_out in unrotated frame 
					projects a vector in the locally rotated frame onto the unrotated frame 

	int project_to(basis_from, v_in,v_out)
	struct basis *basis_from;     direction cosines to go from rotated to unrotated frame 
	double v_in[],v_out[];       v_in here is in unrotated frame, v_out in rotated frame 
					projects a vector in the unrotated frame onto the rotated frame 

	reorient(basis_from,basis_to,v_from,v_to)
	struct basis *basis_from,*basis_to;
	double v_from[],v_to[];
					projects a vector in one rotated frame to another rotated frame.  Note that
					a[i][j] is always the project from a rotated frame to the unrotated (Master) 
					frame.  The projection in the other direction is just the transpose of this


Description:	 



Notes:


History:
 	97jan	ksl	Coded as part of python effort
 	98feb	ksl	Added a routine to project vectors from cartesian coordinates to
 				cylintrical coordinates.

**************************************************************/



#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "log.h"
#include "atomic.h"
#include "python.h"

#define EPS 1.e-10



/* A basis is defined such that if x is a 3 vector as expressed an unprimed cartesian coordinate
   frame, and if y is the same vector in some rotated frame, then
   x[i] = a[i][j] y[j]

History:
	02jan	ksl	Encountered some weird Nan errors with py_ray and so rewrote to
			include a sane_check
 */

// This exists in python.h and so is commented out here.
//struct basis
//{
//  double a[3][3];
//
//};

double
dot (a, b)
     double a[], b[];
{
  double x;

  x = a[0] * b[0] + a[1] * b[1] + a[2] * b[2];

  if (sane_check (x))
    {
      Error ("dot:sane_check of x: a %f %f %f b %f %f %f\n",
	     a[0], a[1], a[2], b[0], b[1], b[2]);
    }
  return (x);
}

/* Return the length of a 3 vector */

double
length (a)
     double a[];
{
  double x, y;
  double sqrt ();
  y = (a[0] * a[0] + a[1] * a[1] + a[2] * a[2]);
  if (sane_check (y))
    Error ("length:sane_check of y: a %f %f %f \n", a[0], a[1], a[2]);
  x = sqrt (y);
  if (sane_check (x))
    Error ("length:sane_check of x=sqrt(y): y= %f\n", y);
  return (x);
}


/* Renormalize a vector so that it has length "scalar" */

int
renorm (a, scalar)
     double a[], scalar;
{
  double x;
  double dot ();
  x = (dot (a, a));
  if (x < EPS)
    {
      Log ("renorm: Cannot renormalize a vector of length 0\n");
      Log ("renorm: %e %e %e\n", a[0], a[1], a[2]);
      Log("renorm returning -1\n");
      return (-1);
    }
  x = scalar / sqrt (x);
  a[0] *= x;
  a[1] *= x;
  a[2] *= x;
  return (0);
}

int
cross (a, b, c)
     double a[], b[], c[];
{
  c[0] = a[1] * b[2] - a[2] * b[1];
  c[1] = a[2] * b[0] - a[0] * b[2];
  c[2] = a[0] * b[1] - a[1] * b[0];

  return (0);
}

/* translate a vector in the direction lmn (must already be normalized) by a distance
   s and report the answer in result */

int
vmove (u, lmn, s, result)
     double u[], lmn[], s, result[];
{
  result[0] = lmn[0] * s + u[0];
  result[1] = lmn[1] * s + u[1];
  result[2] = lmn[2] * s + u[2];
  return (0);
}

int
vsub (u, v, result)
     double u[], v[], result[];
{
  result[0] = u[0] - v[0];
  result[1] = u[1] - v[1];
  result[2] = u[2] - v[2];
  return (0);

}

int
vadd (u, v, result)
     double u[], v[], result[];
{
  result[0] = u[0] + v[0];
  result[1] = u[1] + v[1];
  result[2] = u[2] + v[2];
  return (0);

}

int
stuff_v (vin, vout)
     double vin[], vout[];
{
  vout[0] = vin[0];
  vout[1] = vin[1];
  vout[2] = vin[2];
  return (0);

}

/* Dot product of a 2d tensor and a vector to produce
a vector

The program actually returns the dot product of the 
resultant, so if for example the tensor is the velocity
gradient tensor, and if vin is a direction (of unit length)
then the program will return the change in speed along
the direction of vin. 

	01dec	ksl	Added for python40

*/

double
dot_tensor_vec (tensor, vin, vout)
     double tensor[3][3], vin[3], vout[3];
{
  double dot ();
  vout[0] = dot (tensor[0], vin);
  vout[1] = dot (tensor[1], vin);
  vout[2] = dot (tensor[2], vin);
  return (dot (vin, vout));
}

/* Project a vector b in xyz coords from position a in xyz coords
into cylindrical coordinates */

int
project_from_xyz_cyl (a, b, result)
     double a[], b[], result[];
{

  double n_rho[3], n_phi[3], n_z[3];

  n_rho[0] = a[0];
  n_rho[1] = a[1];
  n_rho[2] = 0;
  if (renorm (n_rho, 1.0) == -1)
    {
      Log
	("Position on z axis; conversion cylintrical coords indeterminate\n");
      return (-1);
    }
  n_z[0] = n_z[1] = 0;
  n_z[2] = 1;
  cross (n_z, n_rho, n_phi);

  result[0] = dot (b, n_rho);
  result[1] = dot (b, n_phi);
  result[2] = b[2];


  return (0);
}


/* Project a vector "b" in which is expressed using a cyl basis from position "a" in 
cartesion (xyz) coords into cartesion xyz coordinates.  Note this is clearly different
from the situation where both a and b are in cylindrical coordinates. */

int
project_from_cyl_xyz (a, b, result)
     double a[], b[], result[];
{
  double x, ctheta, stheta;

  if ((x = sqrt (a[0] * a[0] + a[1] * a[1])) == 0)
    {
      Error ("project_from_cyl_xyz: indeterminate a =%f %f %f\n", a[0], a[1],
	     a[2]);
    }
  ctheta = a[0] / x;
  stheta = a[1] / x;

  result[0] = b[0] * ctheta - b[1] * stheta;
  result[1] = b[0] * stheta + b[1] * ctheta;
  result[2] = b[2];

  return (0);
}



/*  The next few routine have to do with various basis manipulations */

/* create_basis creates an orthnormal basis with the first axis in the u direction, the
   y axis in the uv plane, and the the third axis in the perpendicular direction.  The
   system is right-handed */

int
create_basis (u, v, basis_new)
     double u[], v[];
     struct basis *basis_new;
{
  int i;
  double x[3], y[3], z[3];
  double mu_x;
  double dot ();

  for (i = 0; i < 3; i++)
    {
      x[i] = u[i];
      y[i] = v[i];
    }
  if (renorm (x, 1.) || renorm (y, 1.))
    {
      Log ("Problem creating basis: Either u or v had length 0\n");
      Log ("create_basis: u %e %e %e\n", u[0], u[1], u[2]);
      Log ("create_basis: v %e %e %e\n", v[0], v[1], v[2]);
      return (-1);
    }
  if ((mu_x = dot (x, y)) > 1. - EPS)
    {
      Log ("Problem creating basis u,v parallel\n");
      Log ("create_basis: u %e %e %e\n", u[0], u[1], u[2]);
      Log ("create_basis: v %e %e %e\n", v[0], v[1], v[2]);
      return (-1);
    }
  for (i = 0; i < 3; i++)
    y[i] -= mu_x * x[i];	/* y=y-dot(x,y)*x is orthogonal to x */

  renorm (y, 1.);


  cross (x, y, z);


  for (i = 0; i < 3; i++)
    {
      basis_new->a[i][0] = x[i];
      basis_new->a[i][1] = y[i];
      basis_new->a[i][2] = z[i];
    }

  return (0);
}

/* This routine projects a vector in the locally rotated frame onto the unrotated frame */

int
project_from (basis_from, v_in, v_out)
     struct basis *basis_from;	/* direction cosines to go from rotated to unrotated frame */
     double v_in[], v_out[];	/*v_in here is in rotated frame, v_out in unrotated frame */

{
  int i, j;
  for (i = 0; i < 3; i++)
    {
      v_out[i] = 0;
      for (j = 0; j < 3; j++)
	{
	  v_out[i] += basis_from->a[i][j] * v_in[j];
	}
    }
  return (0);
}

/* This routine projects a vector in the unrotated frame onto the rotated frame */
int
project_to (basis_from, v_in, v_out)
     struct basis *basis_from;	/* direction cosines to go from rotated to unrotated frame */
     double v_in[], v_out[];	/*v_in here is in unrotated frame, v_out in rotated frame */

{
  int i, j;
  for (i = 0; i < 3; i++)
    {
      v_out[i] = 0;
      for (j = 0; j < 3; j++)
	{
	  v_out[i] += basis_from->a[j][i] * v_in[j];
	}
    }
  return (0);
}


/* This routine projects a vector in one rotated frame to another rotated frame.  Note that
   a[i][j] is always the project from a rotated frame to the unrotated (Master) frame.  The projection
   in the other direction is just the transpose of this */

int
reorient (basis_from, basis_to, v_from, v_to)
     struct basis *basis_from, *basis_to;
     double v_from[], v_to[];
{
  double a[3][3];
  int i, j, k;
  for (i = 0; i < 3; i++)
    {
      for (j = 0; j < 3; j++)
	{
	  a[i][j] = 0;
	  for (k = 0; k < 3; k++)
	    {
	      a[i][j] += basis_to->a[k][i] * basis_from->a[k][j];
	    }
	}
    }
  for (i = 0; i < 3; i++)
    {
      v_to[i] = 0;
      for (j = 0; j < 3; j++)
	{
	  v_to[i] += a[i][j] * v_from[j];
	}
    }
  return (0);
}
