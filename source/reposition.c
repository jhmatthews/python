

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "atomic.h"

#include "python.h"

/***********************************************************
                                       Space Telescope Science Institute

 Synopsis:   
	reposition(w,p) attempts to assure that a photon is not scattered 
	a second time inappropriately by the same transition
  
 Arguments:		

	
 Returns:
  
Description:	

Notes:


History:
 	02	ksl	Created to handle problem which was discovered
			in python_40
 
**************************************************************/



int
reposition (w, p)
     WindPtr w;			//w here refers to entire wind, not a single element
     PhotPtr p;

{

  int n;
  int stuff_phot ();
  int where_in_grid ();


  if (p->nres < 0)
    {				// Do nothing for non-resonant scatters
      return (0);
    }

  if ((p->grid = n = where_in_grid (p->x)) < 0)
    {
      Error ("reposition: Photon not in grid when routine entered %d \n", n);
      return (n);		/* Photon was not in wind */
    }


  move_phot (p, DFUDGE);
  return (0);


}
