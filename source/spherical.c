

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "atomic.h"
#include "python.h"

/* Notes on spherical coordinates

   In spherical coordinates, we set NDIM = to whatever and MDIM = 1.  

   Some of the variables in the Wind_ptr array, like the divergence
   have to be calculated at a certain place.  Logically, one would
   do this along one of the principle axes.  But at present we have
   calculated this along a 45 degree angle.  The reason we did this
   was because I want the code to run for normal bipolar winds.  But
   this can be an issue for anisotropic scattering!!  It should not
   affect items like the divergence, but it certainly affects the
   gradients. 

*/
/***********************************************************
                                       Space Telescope Science Institute

 Synopsis:
	spherical_ds_in_cell calculates the distance to the far
        boundary of the cell in which the photon bundle resides.  	
  
 Arguments:		
 	p	Photon pointer

	
 Returns:
 	Distance to the far boundary of the cell in which the photon
	currently resides.  Negative numbers (and zero) should be
	regarded as errors.
  
Description:	

Notes:

History:
 	05apr	ksl	55d: Adapted from rtheta.c
	15aug	ksl	Domains incorporated
 
**************************************************************/

double
spherical_ds_in_cell (p)
     PhotPtr p;

{

  int n, ix;
  double s, smax;
  int ndom;

  ndom = wmain[p->grid].ndom;

  if ((p->grid = n = where_in_grid (ndom, p->x)) < 0)
  {
    Error ("translate_in_wind: Photon not in grid when routine entered\n");
    return (n);                 /* Photon was not in wind */
  }

  ix = n;

  /* Set up the quadratic equations in the radial  direction */

  smax = ds_to_sphere (zdom[ndom].wind_x[ix], p);
  s = ds_to_sphere (zdom[ndom].wind_x[ix + 1], p);
  if (s < smax)
    smax = s;

  if (smax <= 0)
  {
    Error ("spherical: ds_in_cell %f\n", smax);
  }
  return (smax);
}



/***********************************************************
               Space Telescope Science Institute

 Synopsis:
	spherical_make_grid defines the cells in a spherical grid              

Arguments:		
	WindPtr w;	The structure which defines the wind in Python
 
Returns:
 
Description:

	In spherical coordinates w runs from 0 to NDIM.  Note that the 
	centers of the grid cells are defined in the xz plane at a 45 
	degree angle.  This was done so that one would be in a plausible
	region of a biconical wind.


History:
 	05apr	ksl	55d: Adapted from rtheta.c
	05jun	ksl	56a: Eliminated some superfluous variables, and
			added better comments.

**************************************************************/


int
spherical_make_grid (w, ndom)
     WindPtr w;
     int ndom;
{
  double dr, dlogr;
  int j, n;
  int ndim;


  ndim = zdom[ndom].ndim;


  for (j = 0; j < ndim; j++)
  {
    {
      n = j + zdom[ndom].nstart;        // This is the element in wmain

      /*Define the grid points */
      if (zdom[ndom].log_linear == 1)
      {                         // linear intervals

        dr = (zdom[ndom].rmax - geo.rstar) / (ndim - 3);
        w[n].r = geo.rstar + j * dr;
        w[n].rcen = w[n].r + 0.5 * dr;
      }
      else
      {                         //logarithmic intervals

        dlogr = (log10 (zdom[ndom].rmax / zdom[ndom].rmin)) / (ndim - 3);
        w[n].r = zdom[ndom].rmin * pow (10., dlogr * (j - 1));
        w[n].rcen = 0.5 * zdom[ndom].rmin * (pow (10., dlogr * (j)) + pow (10., dlogr * (j - 1)));
        Log ("New W.r = %e, w.rcen = %e\n", w[n].r, w[n].rcen);
      }

      /* Now calculate the positions of these points in the xz plane.
         There is a choice about how one does this.  Here we  have elected
         to calculate this at a 45 degree angle.  in the hopes this will 
         be a reasonable portion of the wind in a biconical flow.
       */

      w[n].x[1] = w[n].xcen[1] = 0.0;
      w[n].x[0] = w[n].x[2] = w[n].r * sin (PI / 4.);
      w[n].xcen[0] = w[n].xcen[2] = w[n].rcen * sin (PI / 4.);

    }
  }

  return (0);

}



/***********************************************************
                        Space Telescope Science Institute

 Synopsis:
	spherical_wind_complete (w)

 Arguments:		
	WindPtr w;    the entire wind
 Returns:

 Description:
 	This simple little routine just populates one dimensional 
	arrays that are used for interpolation.  

 Notes:
 History:
 	05apr	ksl	55d: Adapted from rtheta.c
	15aug	ksl	Add domains
	16mar	ksl	Removed reference to mdim as this is
			a spherical grid
 
**************************************************************/


int
spherical_wind_complete (ndom, w)
     int ndom;
     WindPtr w;
{
  int i;
  int ndim, nstart;

  ndim = zdom[ndom].ndim;
  nstart = zdom[ndom].nstart;


  for (i = 0; i < ndim; i++)
    zdom[ndom].wind_x[i] = w[nstart + i].r;
  for (i = 0; i < ndim - 1; i++)
    zdom[ndom].wind_midx[i] = w[nstart + i].rcen;
  /* Add something plausible for the edges */
  zdom[ndom].wind_midx[ndim - 1] = 2. * zdom[ndom].wind_x[ndim - 1] - zdom[ndom].wind_midx[ndim - 2];
  zdom[ndom].wind_midz[ndim - 1] = 2. * zdom[ndom].wind_z[ndim - 1] - zdom[ndom].wind_midz[ndim - 2];

  return (0);
}


/***********************************************************
                                       Space Telescope Science Institute

 Synopsis:
 	spherical_volume(w) calculates the wind volume of cell
	allowing for the fact that some cells are partially in the wind

 Arguments:		
	WindPtr w;    the entire wind
 Returns:

 Description:
		
 Notes:
 	This is like not the best way to do this integration. It
	should probaly be done interns of sin theta, rather than
	theta, since this would evenly sample the volume which after
	all is what we are trying to calculate. ksl 05apr 

 History:
 	05apr	ksl	55d: Adapted from rtheta.c
	06nov	ksl	58b: Minor modification to use W_ALL_INWIND
			etc., instead of hardcoded values
	11aug	ksl	70b Add the ability to find a different 
			component.  Note that this makes explicit
			use of the way components are defined,
			See python.h
	15aug	ksl	Added support for domains
 
**************************************************************/
#define RESOLUTION   100


int
spherical_volumes (ndom, w)
     int ndom;
     WindPtr w;
{
  int i, n;
  double fraction;
  double num, denom;
  double r, theta;
  double dr, dtheta, x[3];
  double rmin, rmax;
  double thetamin, thetamax;
  int jj, kk;
  int ndim, nstart, ndomain;

  ndim = zdom[ndom].ndim;
  nstart = zdom[ndom].nstart;


  thetamin = 0.0;
  thetamax = 0.5 * PI;

  for (i = 0; i < ndim; i++)
  {
    {
      n = i + nstart;           // nstart is the offset into the wind cell
      rmin = zdom[ndom].wind_x[i];
      rmax = zdom[ndom].wind_x[i + 1];

      w[n].vol = 4. / 3. * PI * (rmax * rmax * rmax - rmin * rmin * rmin);

      if (i == ndim - 1)
      {
        fraction = 0.0;         /* Force outside edge volues to zero */
        jj = 0;
        kk = RESOLUTION;
      }
      else if (i == ndim - 2)
      {
        fraction = 0.0;         /* Force outside edge volues to zero */
        jj = 0;
        kk = RESOLUTION;
      }
      else
      {                         /* Determine whether the cell is in the wind */
        num = denom = 0;
        jj = kk = 0;
        dr = (rmax - rmin) / RESOLUTION;
        dtheta = (thetamax - thetamin) / RESOLUTION;
        for (r = rmin + dr / 2; r < rmax; r += dr)
        {
          for (theta = thetamin + dtheta / 2; theta < thetamax; theta += dtheta)
          {
            denom += r * r * sin (theta);;
            kk++;
            x[0] = r * sin (theta);
            x[1] = 0;
            x[2] = r * cos (theta);;
            if (where_in_wind (x, &ndomain) == W_ALL_INWIND)
            {
              num += r * r * sin (theta);       /* 0 implies in wind */
              jj++;
            }
          }
        }
        fraction = num / denom;
      }
      if (jj == 0)
      {
        w[n].inwind = W_NOT_INWIND;     // The cell is not in the wind
        w[n].vol = 0.0;
      }
      else if (jj == kk)
        w[n].inwind = W_ALL_INWIND;     // The cell is completely in the wind
      else
      {
        w[n].inwind = W_PART_INWIND;    //The cell is partially in the wind
        w[n].vol *= fraction;
      }


    }
  }

  return (0);
}


/***********************************************************
                     Space Telescope Science Institute

 Synopsis:
 	spherical_where_in_grid locates the element in wmain,
	when one is using spherical coordinates. 

 Arguments:		
	double x[];
 Returns:
 	where_in_grid normally  returns the wind element  associated with
 		a position.  If the position is in the grid this will be a positive
 		integer < NDIM*MDIM.
 	x is inside the grid        -1
	x is outside the grid       -2
 Description:	
	
		
 Notes:
	Where_in grid does not tell you whether the x is in the wind or not. 

	What one means by inside or outside the grid may well be different
	for different coordinate systems.

	This routine is not normally used directly.  Instead it is called bia
	where_in_grid which determines which coordinate specific where_in_grid
	routine to call.

 History:
 	05apr	ksl	55d: Adapted from rtheta.c.  
	13sep	nsh	76b Changed call to fraction to take account of new mode
	1605	ksl	Updated so it returns the element in wmain to
			make this routine consistent with that of
			the same routine for other coordiante systems
 
**************************************************************/



int
spherical_where_in_grid (ndom, x)
     int ndom;
     double x[];
{
  int n;
  double r;
  double f;
  int ndim;


  ndim = zdom[ndom].ndim;

  r = length (x);

  /* Check to see if x is outside the region of the calculation */
  if (r > zdom[ndom].wind_x[ndim - 1])
  {
    return (-2);                /* x is outside grid */
  }
  else if (r < zdom[ndom].wind_x[0])
  {
    return (-1);                /*x is inside grid */
  }

  fraction (r, zdom[ndom].wind_x, ndim, &n, &f, 0);

  // n is the position with this domain, so zdom[ndom].nstart is added 
  // get to wmain

  return (n + zdom[ndom].nstart);
}

/***********************************************************
                     Space Telescope Science Institute

 Synopsis:
 	spherical_get_random_location

 Arguments:		
 	int n -- Cell in which random position is to be generated
 Returns:
 	double x -- the position
 Description:	
	
		
 Notes:



 History:
 	05apr	ksl	55d: Adapted from rtheta.c
	11aug	ksl	70b - Modified to account for torus
 
**************************************************************/

int
spherical_get_random_location (n, x)
     int n;                     // Cell in which to create position
     double x[];                // Returned position
{
  int i, j;
  int inwind;
  double r, rmin, rmax;
  double theta, phi;
  int ndom, ndomain;

  ndom = wmain[n].ndom;
  wind_n_to_ij (ndom, n, &i, &j);
  rmin = zdom[ndom].wind_x[i];
  rmax = zdom[ndom].wind_x[i + 1];


  /* Generate a position which is both in the cell and in the wind */
  inwind = W_NOT_INWIND;
  while (inwind != W_ALL_INWIND)
  {
    r = (rmin * rmin * rmin) + (rmax * rmax * rmax - rmin * rmin * rmin) * (rand () / (MAXRAND - 0.5));
    r = pow (r, (1. / 3.));
    theta = acos (2. * (rand () / MAXRAND) - 1);

    phi = 2. * PI * (rand () / MAXRAND);
/* Project from r, theta phi to x y z  */
    x[0] = r * cos (phi) * sin (theta);
    x[1] = r * sin (phi) * sin (theta);
    x[2] = r * cos (theta);
    inwind = where_in_wind (x, &ndomain);       /* Some photons will not be in the wind */
  }

  return (inwind);
}



/***********************************************************
                     Space Telescope Science Institute

 Synopsis:
 	spherical_extend_density  extends the density to
	regions just outside the wind regiions so that
	extrapolations of density can be made there

 Arguments:		
 Returns:
 Description:	
 Notes:
 	There are several reasons for this code in a spherical
	wind, but the main one is probably to make sure models
	that really should not be modelled in python still actually
	run without producing silly computational errors.

 	SS asked whether we should also be extending the wind for other 
	parameters, especially ne.  At present we do not interpolate
	on ne so this is not necessary.  If we did do that it would be required.

	In cylindrical coordinates, the fast dimension is z; grid positions 
	increase up in z, and then out in r.  In spperical polar coordinates, 
	the fast dimension is theta; the grid increases in theta (measured) from 
	the z axis), and then in r.  In spherical coordinates, 
	the grid increases as one might expect in r..
	
		



 History:
	05apr	ksl	56 -- Moved functionality from wind updates   
	06may	ksl	57+ -- Using mappings in attempt to push everything into
			plasma structure.  OK as long as we do not update
			something incorrectly as a result.
 
**************************************************************/


int
spherical_extend_density (ndom, w)
     int ndom;
     WindPtr w;
{

  int j, n, m;
  int ndim, nstart;

  ndim = zdom[ndom].ndim;
  nstart = zdom[ndom].nstart;
  /* 
     Now we need to updated the densities immediately outside the wind so that the density interpolation in resonate will work.
     In this case all we have done is to copy the densities from the cell which is just in the wind (as one goes outward) to the
     cell that is just inside (or outside) the wind. 
   */

  for (j = 0; j < ndim - 1; j++)
  {
    n = nstart + j;
    if (w[n].vol == 0)          // Then the grid point is not in the wind

    {
      m = n + 1;
      if (w[m].vol > 0)         // Then grid point n is just inside the wind
      {                         // To extend we copy  copy the densities to the grid cell n
        w[n].nplasma = w[m].nplasma;

      }
      else if (j > 0)
      {
        m = n - 1;
        if (w[m].vol > 0)       // The grid point is just outside the wind
        {
          w[n].nplasma = w[m].nplasma;

        }
      }
    }
  }

  return (0);

}


/***********************************************************
               Space Telescope Science Institute

 Synopsis:
	shell_make_grid defines the cells in a thin shell. One shell inside the shell, one outside and one shell exactly fitting the shell.            

Arguments:		
	WindPtr w;	The structure which defines the wind in Python
 
Returns:
 
Description:

	In spherical coordinates w runs from 0 to NDIM.  Note that the 
	centers of the grid cells are defined in the xz plane at a 45 
	degree angle.  This was done so that one would be in a plausible
	region of a biconical wind.
        NSH - Ive changed it to use 1/root2 rather than sin 45. More accurate.


History:
 	11feb	nsh	Adapted from spherical_make_grid.c


**************************************************************/


int
shell_make_grid (w, ndom)
     WindPtr w;
     int ndom;
{
  int n;
  int ndim;

  ndim = zdom[ndom].ndim;


  w[0].r = zdom[ndom].rmin - (zdom[ndom].rmax - zdom[ndom].rmin);
  w[1].r = zdom[ndom].rmin;
  w[2].r = zdom[ndom].rmax;
  w[3].r = zdom[ndom].rmax + (zdom[ndom].rmax - zdom[ndom].rmin);



  w[0].rcen = (w[0].r + w[1].r) / 2;
  w[1].rcen = (w[1].r + w[2].r) / 2;
  w[2].rcen = (w[2].r + w[3].r) / 2;
  w[3].rcen = w[2].rcen + (zdom[ndom].rmax - zdom[ndom].rmin);




  /* Now calculate the positions of these points in the xz plane.
     There is a choice about how one does this.   I have elected
     to assume that we want to calculate this at a 45 degree angle.
     in the hopes this will be a reasonable portion of the wind in
     a biconical flow.
   */
  for (n = 0; n < ndim; n++)
  {
    Log ("Cell %i:  inner edge = %2.20e, centre = %2.20e\n", n, w[n].r, w[n].rcen);
    w[n].x[1] = w[n].xcen[1] = 0.0;

    //NSH Slight change here, using 1/root2 give more accurate results than sin45.


    w[n].x[0] = w[n].x[2] = w[n].r / pow (2.0, 0.5);
    w[n].xcen[0] = w[n].xcen[2] = w[n].rcen / pow (2.0, 0.5);
  }


  return (0);

}
