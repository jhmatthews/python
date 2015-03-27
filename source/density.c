/***********************************************************
                                       Space Telescope Science Institute

Synopsis:
	get_ion_density finds the density of a specific ion designated by nion, at a position held in
the photon pointer p

Arguments:		

Returns:
 
Description:	
Notes:

	Densities (rho) are defined at the centers of cells. rho is defined in wind2d.c
	an etends from 0 to NDIM-1 (inclusive), and 0 to MDIM-1 (inclusive) . However, 
	in wind_updates2d.c, the ionization fractions are only calculated from 0 to NDIM-2, 
        and 0 to MDIM-2 inclusive.


History:
 	01mar	ksl	Added this possibility
	04nov	ksl	53b: Made modifications to allow for 
			different coordinate systems
	05apr	ksl	55d: Adapted to new version of coord_fraction
			in continuing effort to enable more 
			coordinate systems.
	06may	ksl	57+ -- Modified to use plasma structure
  1411 JM Modified to use a general vector x, rather than a PhotPtr
**************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "atomic.h"
#include "python.h"



double
get_ion_density (x, nion)
     double x[];
     int nion;
{
  double dd;
  int nn, nnn[4], nelem;
  double frac[4];
  int nplasma;

/* nnn is an array which contains the elements of the Wind (coordinate)
structure that must be summed. We need to convert that to a position
in the plasma structure*/

  if ((coord_fraction (1, x, nnn, frac, &nelem)) > 0)
    {

      dd = 0;

      for (nn = 0; nn < nelem; nn++)
	{
	  nplasma = wmain[nnn[nn]].nplasma;
	  dd += plasmamain[nplasma].density[nion] * frac[nn];
	}
    }
  else
    {
      dd = 0;
    }


  return (dd);
}
