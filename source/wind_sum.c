
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>


#include "atomic.h"
#include "python.h"

/***********************************************************
                         Space Telescope Science Institute

 Synopsis:
	xtemp_rad (w) prints out t_r and t_e in the grid

 Arguments:		
	WindPtr w;

Returns:
 
Description:	
	
Notes:
	This was added in spring 2004 to try to understand how the
	temperature was evolving in the models

	This routine is adapted from a very simillar routine in py_wind 

	The routine won't fail but the output may be less useful
	when coodinates sysstems other than cylindrical are
	used.
History:
 	04mar	ksl	Coding on python began.
        04May	SS	Modified so that Te is given too.
	06may	ksl	57+ -- to reflect plasma structure
	09feb	ksl	68b - Added back subsampling when dimensions
			are very large, and a warning about
			the display when not using cylindrical 
			coords.

**************************************************************/

int
xtemp_rad (w)
     WindPtr w;
{
  int i, j, n;
  double x;
  int py_wind_min, py_wind_max, py_wind_delta;
  int nplasma;

  py_wind_min = 0;
  py_wind_max = NDIM;

  /* py_wind_delta can be used to subsample the array */
  py_wind_delta = 1;
  if (MDIM > 30)
    py_wind_delta = 1 + MDIM / 30;


  if (geo.coord_type != 1)
    {
      Log
	("Warning: Since coord type is not cylindrical, next print out may look odd\n");
    }

  Log ("\n T rad\n");

  Log ("   z\\x   ");
  for (i = py_wind_min; i < py_wind_max; i += py_wind_delta)
    Log ("%8.2e ", w[i * MDIM].x[0]);
  Log ("\n");

  for (j = 0; j < MDIM; j++)
    {
      Log ("%8.2e ", w[j].x[2]);
      for (i = py_wind_min; i < py_wind_max; i += py_wind_delta)
	{
	  n = i * MDIM + j;
	  if (w[n].vol > 0.0)
	    {
	      nplasma = w[n].nplasma;
	      x = plasmamain[nplasma].t_r;
	    }
	  else
	    x = 0.0;
	  Log ("%8.2g ", x);
	}
      Log ("\n");
    }

  Log ("\n T e\n");

  Log ("   z\\x   ");
  for (i = py_wind_min; i < py_wind_max; i += py_wind_delta)
    Log ("%8.2e ", w[i * MDIM].x[0]);
  Log ("\n");

  for (j = 0; j < MDIM; j++)
    {
      Log ("%8.2e ", w[j].x[2]);
      for (i = py_wind_min; i < py_wind_max; i += py_wind_delta)
	{
	  n = i * MDIM + j;
	  if (w[n].vol > 0.0)
	    {
	      nplasma = w[n].nplasma;
	      x = plasmamain[nplasma].t_e;
	    }
	  else
	    x = 0.0;
	  Log ("%8.2g ", x);
	}
      Log ("\n");
    }

  return (0);

}
