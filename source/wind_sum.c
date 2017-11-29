
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
	
	This is just a logging routine
	
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
	15aug	ksl	Modified for domains, on the assumption
			that we wanted to print out results for
			all of the domians

**************************************************************/

int
xtemp_rad (w)
     WindPtr w;
{
  int i, j, n;
  double x;
  int py_wind_min, py_wind_max, py_wind_delta;
  int nplasma;
  int ndom, ndim, mdim, nstart;
  int ntot;



  for (ndom = 0; ndom < geo.ndomain; ndom++)
    {

      Log ("Results for Domain %d:\n", ndom);
      ndim = zdom[ndom].ndim;
      mdim = zdom[ndom].mdim;
      nstart = zdom[ndom].nstart;



      py_wind_min = 0;
      py_wind_max = ndim;

      /* py_wind_delta can be used to subsample the array */
      py_wind_delta = 1;
      if (mdim > 30)
	py_wind_delta = 1 + mdim / 30;

      /* XXX - Nick had commented the next secion out for reasons which are totally unclear */

      if (zdom[ndom].coord_type == SPHERICAL)
	{
	  Log ("\n R\n");
	  j = 1;
	  for (i = py_wind_min; i < py_wind_max; i += 1)
	    {
	      Log ("%8.2e ", w[nstart + i * mdim].r);
	      if (j % 10 == 0)
		{
		  Log ("\n");
		}
	      j++;
	    }
	  Log ("\n T rad\n");

	  j = 1;
	  for (i = py_wind_min; i < py_wind_max; i += 1)
	    {
	      n = nstart + i;
	      Log ("%8.2e ", plasmamain[n].t_r);
	      if (j % 10 == 0)
		{
		  Log ("\n");
		}
	      j++;
	    }

	  Log ("\n T e\n");

	  j = 1;
	  for (i = py_wind_min; i < py_wind_max; i += 1)
	    {
	      n = nstart + i;
	      Log ("%8.2e ", plasmamain[n].t_e);
	      if (j % 10 == 0)
		{
		  Log ("\n");
		}
	      j++;
	    }

	  Log ("\n nphot\n");

	  j = 1;
	  for (i = py_wind_min; i < py_wind_max; i += 1)
	    {
	      n = nstart + i;
	      if (plasmamain[n].vol > 0.0)
		{
		  ntot = plasmamain[n].ntot;
		}
	      else
		ntot = 0;
	      Log ("%8d ", ntot);
	      if (j % 10 == 0)
		{
		  Log ("\n");
		}
	      j++;
	    }


	  return (0);
	}
      /* 2d coord systems */
      if (zdom[ndom].coord_type != CYLIND)
	{
	  Log
	    ("Warning: Since coord type is not cylindrical, next print out may look odd\n");
	}


      Log ("\n T rad\n");

      Log ("   z\\x   ");
      for (i = py_wind_min; i < py_wind_max; i += py_wind_delta)
	Log ("%8.2e ", w[nstart + i * mdim].x[0]);
      Log ("\n");

      for (j = 0; j < mdim; j++)
	{
	  Log ("%8.2e ", w[j].x[2]);
	  for (i = py_wind_min; i < py_wind_max; i += py_wind_delta)
	    {
	      n = nstart + i * mdim + j;
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
	Log ("%8.2e ", w[nstart + i * mdim].x[0]);
      Log ("\n");

      for (j = 0; j < mdim; j++)
	{
	  Log ("%8.2e ", w[j].x[2]);
	  for (i = py_wind_min; i < py_wind_max; i += py_wind_delta)
	    {
	      n = nstart + i * mdim + j;
	      nplasma = w[n].nplasma;
	      if (plasmamain[nplasma].vol > 0.0)
		{
		  x = plasmamain[nplasma].t_e;
		}
	      else
		x = 0.0;
	      Log ("%8.2g ", x);
	    }
	  Log ("\n");
	}

      Log ("\n ntot \n");

      Log ("   z\\x   ");
      for (i = py_wind_min; i < py_wind_max; i += py_wind_delta)
	Log ("%8.2e ", w[nstart + i * mdim].x[0]);
      Log ("\n");

      for (j = 0; j < mdim; j++)
	{
	  Log ("%8.2e ", w[j].x[2]);
	  for (i = py_wind_min; i < py_wind_max; i += py_wind_delta)
	    {
	      n = nstart + i * mdim + j;
	      nplasma = w[n].nplasma;
	      if (plasmamain[nplasma].vol > 0.0)
		{
		  ntot = plasmamain[nplasma].ntot;
		}
	      else
		ntot = 0;
	      Log ("%8d ", ntot);
	    }
	  Log ("\n");
	}
    }



  return (0);

}
