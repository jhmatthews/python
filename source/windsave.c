/* Coordinate system independent, except for wind_complete, which is fixed 03aug */

/***********************************************************
                                       Space Telescope Science Institute

 Synopsis:
	wind_save(w,filename)
	wind_read(filename)
	spec_save(filename)
	spec_read(filename)

Arguments:		



Returns:
 
Description:	
	
	The first two routines in this file write and read the wind structure.  		
	The second two routines do the same thing for the spectrum structure

Notes:


History:
 	98mar	ksl 	Replaced old ascii routine with binary read and write
 	98mar20	ksl	Fixed problem in createing wind_midx and wind_midz
	02apr	ksl	Changed windread so that it would dynamically allocate
			the WindPtr array after reading geo.  In the process
			changed the call to windread since one cannot assign
			w inside the routine and have that assigment be durable
			outside of the subroutine.  Instead we set wmain.
	06may	ksl	Added read & write statement for plasma structure
	08mar	ksl	60 - Added read & write statements for macro structure
	08may	ksl	60a - Fixed write statement for macro structure so not
			written when no macro atoms
	08dec	ksl	67c - Added routines to read and write the 
			spec structure as part of the general effort
			to allow python to restart. Modified the call to
			wind_save to eliminate superfluous passing of
			the wind ptr.
	14jul	nsh	78a - Added code to allow dynamically allocated arrays
			in the plasma structure to be read in and written out.
 
**************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "atomic.h"

#include "python.h"


int
wind_save (filename)
     char filename[];
{
  FILE *fptr, *fopen ();
  char line[LINELENGTH];
  int n, m;

  if ((fptr = fopen (filename, "w")) == NULL)
    {
      Error ("wind_save: Unable to open %s\n", filename);
      exit (0);
    }

  sprintf (line, "Version %s\n", VERSION);
  n = fwrite (line, sizeof (line), 1, fptr);
  n += fwrite (&geo, sizeof (geo), 1, fptr);
  n += fwrite (wmain, sizeof (wind_dummy), NDIM2, fptr);
  n += fwrite (plasmamain, sizeof (plasma_dummy), NPLASMA, fptr);

/* NSH 1407 - The following loop writes out the variable length arrays
in the plasma structure */
      for (m = 0; m < NPLASMA; m++)
{
  n += fwrite (plasmamain[m].density,sizeof(double),nions,fptr);
  n += fwrite (plasmamain[m].partition,sizeof(double),nions,fptr);

  n += fwrite (plasmamain[m].PWdenom,sizeof(double),nions,fptr);
  n += fwrite (plasmamain[m].PWdtemp,sizeof(double),nions,fptr);
  n += fwrite (plasmamain[m].PWnumer,sizeof(double),nions,fptr);
  n += fwrite (plasmamain[m].PWntemp,sizeof(double),nions,fptr);

  n += fwrite (plasmamain[m].ioniz,sizeof(double),nions,fptr);
  n += fwrite (plasmamain[m].recomb,sizeof(double),nions,fptr);

  n += fwrite (plasmamain[m].scatters,sizeof(int),nions,fptr);
  n += fwrite (plasmamain[m].xscatters,sizeof(double),nions,fptr);

  n += fwrite (plasmamain[m].heat_ion,sizeof(double),nions,fptr);
  n += fwrite (plasmamain[m].lum_ion,sizeof(double),nions,fptr);
}




  if (geo.nmacro)
    {
      n += fwrite (macromain, sizeof (macro_dummy), NPLASMA, fptr);
      for (m = 0; m < NPLASMA; m++)
	{
	  n +=
	    fwrite (macromain[m].jbar, sizeof (double), size_Jbar_est, fptr);
	  n +=
	    fwrite (macromain[m].jbar_old, sizeof (double), size_Jbar_est,
		    fptr);
	  n +=
	    fwrite (macromain[m].gamma, sizeof (double), size_gamma_est,
		    fptr);
	  n +=
	    fwrite (macromain[m].gamma_old, sizeof (double), size_gamma_est,
		    fptr);
	  n +=
	    fwrite (macromain[m].gamma_e, sizeof (double), size_gamma_est,
		    fptr);
	  n +=
	    fwrite (macromain[m].gamma_e_old, sizeof (double), size_gamma_est,
		    fptr);
	  n +=
	    fwrite (macromain[m].alpha_st, sizeof (double), size_gamma_est,
		    fptr);
	  n +=
	    fwrite (macromain[m].alpha_st_old, sizeof (double),
		    size_gamma_est, fptr);
	  n +=
	    fwrite (macromain[m].alpha_st_e, sizeof (double), size_gamma_est,
		    fptr);
	  n +=
	    fwrite (macromain[m].alpha_st_e_old, sizeof (double),
		    size_gamma_est, fptr);
	  n +=
	    fwrite (macromain[m].recomb_sp, sizeof (double), size_alpha_est,
		    fptr);
	  n +=
	    fwrite (macromain[m].recomb_sp_e, sizeof (double), size_alpha_est,
		    fptr);
	  n +=
	    fwrite (macromain[m].matom_emiss, sizeof (double), nlevels_macro,
		    fptr);
	  n +=
	    fwrite (macromain[m].matom_abs, sizeof (double), nlevels_macro,
		    fptr);

	}

    }

  fclose (fptr);

  Log
    ("wind_write sizes: NPLASMA %d size_Jbar_est %d size_gamma_est %d size_alpha_est %d nlevels_macro %d\n",
     NPLASMA, size_Jbar_est, size_gamma_est, size_alpha_est, nlevels_macro);

  return (n);

}

/*

11dec	ksl	Updated so returns -1 if it cannot open the windsave file.  This
		was done to enable one to handle missing files differently in
		different cases
14jul	nsh	Code added to read in variable length arrays in plasma structure
*/
int
wind_read (filename)
     char filename[];
{
  FILE *fptr, *fopen ();
  int n, m;
  char line[LINELENGTH];
  char version[LINELENGTH];

  if ((fptr = fopen (filename, "r")) == NULL)
    {
      return (-1);
    }

  n = fread (line, sizeof (line), 1, fptr);
  sscanf (line, "%*s %s", version);
  Log
    ("Reading Windfile %s created with python version %s with python version %s\n",
     filename, version, VERSION);

  n += fread (&geo, sizeof (geo), 1, fptr);

  /* 091105 - ksl - Now read the atomic data file.  This is necessary to do here in order to establish the 
   * values for the dimensionality of some of the variable length structures, associated 
   * with macro atoms, especially but likely to be a good idea ovrall
   */

  get_atomic_data (geo.atomic_filename);




/* Now allocate space for the wind array */

  ndim = geo.ndim;
  mdim = geo.mdim;
  NDIM = ndim;
  MDIM = mdim;
  NDIM2 = ndim * mdim;
  NPLASMA = geo.nplasma;

  calloc_wind (NDIM2);
  n += fread (wmain, sizeof (wind_dummy), NDIM2, fptr);

  calloc_plasma (NPLASMA);

  n += fread (plasmamain, sizeof (plasma_dummy), NPLASMA, fptr);

  calloc_dyn_plasma (NPLASMA); /*NSH 1407 allocate space for the dynamically allocated plasma arrays */
  for (m = 0; m < NPLASMA; m++)
{

  n += fread (plasmamain[m].density,sizeof(double),nions,fptr);
  n += fread (plasmamain[m].partition,sizeof(double),nions,fptr);

  n += fread (plasmamain[m].PWdenom,sizeof(double),nions,fptr);
  n += fread (plasmamain[m].PWdtemp,sizeof(double),nions,fptr);
  n += fread (plasmamain[m].PWnumer,sizeof(double),nions,fptr);
  n += fread (plasmamain[m].PWntemp,sizeof(double),nions,fptr);

  n += fread (plasmamain[m].ioniz,sizeof(double),nions,fptr);
  n += fread (plasmamain[m].recomb,sizeof(double),nions,fptr);

  n += fread (plasmamain[m].scatters,sizeof(int),nions,fptr);
  n += fread (plasmamain[m].xscatters,sizeof(double),nions,fptr);

  n += fread (plasmamain[m].heat_ion,sizeof(double),nions,fptr);
  n += fread (plasmamain[m].lum_ion,sizeof(double),nions,fptr);
}



  if (geo.nmacro > 0)
    {
      calloc_macro (NPLASMA);
      n += fread (macromain, sizeof (macro_dummy), NPLASMA, fptr);
      calloc_estimators (NPLASMA);

      for (m = 0; m < NPLASMA; m++)
	{
	  n +=
	    fread (macromain[m].jbar, sizeof (double), size_Jbar_est, fptr);
	  n +=
	    fread (macromain[m].jbar_old, sizeof (double), size_Jbar_est,
		   fptr);
	  n +=
	    fread (macromain[m].gamma, sizeof (double), size_gamma_est, fptr);
	  n +=
	    fread (macromain[m].gamma_old, sizeof (double), size_gamma_est,
		   fptr);
	  n +=
	    fread (macromain[m].gamma_e, sizeof (double), size_gamma_est,
		   fptr);
	  n +=
	    fread (macromain[m].gamma_e_old, sizeof (double), size_gamma_est,
		   fptr);
	  n +=
	    fread (macromain[m].alpha_st, sizeof (double), size_gamma_est,
		   fptr);
	  n +=
	    fread (macromain[m].alpha_st_old, sizeof (double), size_gamma_est,
		   fptr);
	  n +=
	    fread (macromain[m].alpha_st_e, sizeof (double), size_gamma_est,
		   fptr);
	  n +=
	    fread (macromain[m].alpha_st_e_old, sizeof (double),
		   size_gamma_est, fptr);
	  n +=
	    fread (macromain[m].recomb_sp, sizeof (double), size_alpha_est,
		   fptr);
	  n +=
	    fread (macromain[m].recomb_sp_e, sizeof (double), size_alpha_est,
		   fptr);
	  n +=
	    fread (macromain[m].matom_emiss, sizeof (double), nlevels_macro,
		   fptr);
	  n +=
	    fread (macromain[m].matom_abs, sizeof (double), nlevels_macro,
		   fptr);

	  /* Force recalculation of kpkt_rates */

	  macromain[m].kpkt_rates_known = 0;
	}

    }

  fclose (fptr);

  wind_complete (wmain);

  Log ("Read geometry and wind structures from windsavefile %s\n", filename);

  return (n);

}


/* This routine now just calls routines that create 1 dimensional arrays used in
 * various interpolations.  It is separate from coordinate_generation, because these
 * arrays are not part of wind_save, and hence must be regenerated in py_wind.
 *
 * It might make more sense just to save the interpolation arrays and then one might
 * be able to reduce the total number of subroutines, by absorbint the calculations
 * done by wind_complete into another subroutine.
 *
 * History:
 * 	04aug	ksl	Routine rewritten just to be a driver 
 */

int
wind_complete (w)
     WindPtr w;
{

  if (geo.coord_type == SPHERICAL)
    {
      spherical_wind_complete (w);
    }
  else if (geo.coord_type == CYLIND)
    {
      cylind_wind_complete (w);
    }
  else if (geo.coord_type == RTHETA)
    {
      rtheta_wind_complete (w);
    }
  else if (geo.coord_type == CYLVAR)
    {
      cylvar_wind_complete (w);
    }
  else
    {
      Error ("wind_complete: Don't know how to complete coord_type %d\n",
	     geo.coord_type);
      exit (0);
    }

  return (0);
}

int
spec_save (filename)
     char filename[];
{

  FILE *fptr, *fopen ();
  char line[LINELENGTH];
  int n;

  if ((fptr = fopen (filename, "w")) == NULL)
    {
      Error ("spec_save: Unable to open %s\n", filename);
      exit (0);
    }

  sprintf (line, "Version %s  nspectra %d\n", VERSION, nspectra);
  n = fwrite (line, sizeof (line), 1, fptr);
  n += fwrite (xxspec, sizeof (spectrum_dummy), nspectra, fptr);
  fclose (fptr);

  return (n);
}


int
spec_read (filename)
     char filename[];
{
  FILE *fptr, *fopen ();
  int n;

  char line[LINELENGTH];
  char version[LINELENGTH];

  if ((fptr = fopen (filename, "r")) == NULL)
    {
      Error ("spec_read: Unable to open %s\n", filename);
      exit (0);
    }
  n = fread (line, sizeof (line), 1, fptr);
  sscanf (line, "%*s %s %*s %d", version, &nspectra);
  Log
    ("Reading specfile %s with %d spectra created with python version %s with python version %s\n",
     filename, nspectra, version, VERSION);


  /* First allocate space */

  xxspec = calloc (sizeof (spectrum_dummy), nspectra);
  if (xxspec == NULL)
    {
      Error
	("spectrum_init: Could not allocate memory for %d spectra with %d wavelengths\n",
	 nspectra, NWAVE);
      exit (0);
    }

/* Now read the rest of the file */

  n += fread (xxspec, sizeof (spectrum_dummy), nspectra, fptr);

  fclose (fptr);

  Log ("Read spec structures from specfile %s\n", filename);

  return (n);

}
