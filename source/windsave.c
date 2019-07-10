
/***********************************************************/
/** @file  windsave.c
 * @author ksl
 * @date   March, 2018
 *
 * @brief  Routines to save and read in the structues the constitute
 * a model and and spectra which have been genrerated
 *
 * 
 * The first two routines in this file write and read the wind structure.  		
 * The second two routines do the same thing for the spectrum structure
 * 
 * ### Notes ###
 *
 * The files here are all written out as binary files.  They are
 * used for restars, and also by routines like py_wind and windsave2talbe
 * which inspect what is happening in the wind.
 *
 * There are separate ascii_writing 
 * routines for writing the spectra out for plotting.)
 * 
 ***********************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "atomic.h"
#include "python.h"



/**********************************************************/
/** 
 * @brief      Save all of the strutures associated with the 
 * wind to a file
 *
 * @param [in] char  filename[]   The name of the file to write to
 * @return     The number of successful writes
 *
 * @details
 *
 * ### Notes ###
 *
 * For the most part, adding a variable to the structures geo,
 * or plasma, does not require changes to this routine, unless
 * new variable length arrays are involved.
 *
 **********************************************************/

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
    Exit (0);
  }

  sprintf (line, "Version %s\n", VERSION);
  n = fwrite (line, sizeof (line), 1, fptr);
  n += fwrite (&geo, sizeof (geo), 1, fptr);
  n += fwrite (zdom, sizeof (domain_dummy), geo.ndomain, fptr);
  n += fwrite (wmain, sizeof (wind_dummy), NDIM2, fptr);
  n += fwrite (&disk, sizeof (disk), 1, fptr);
  n += fwrite (&qdisk, sizeof (disk), 1, fptr);
  n += fwrite (plasmamain, sizeof (plasma_dummy), NPLASMA, fptr);

/* Write out the variable length arrays
in the plasma structure */

  for (m = 0; m < NPLASMA; m++)
  {
    n += fwrite (plasmamain[m].density, sizeof (double), nions, fptr);
    n += fwrite (plasmamain[m].partition, sizeof (double), nions, fptr);

    n += fwrite (plasmamain[m].ioniz, sizeof (double), nions, fptr);
    n += fwrite (plasmamain[m].recomb, sizeof (double), nions, fptr);
    n += fwrite (plasmamain[m].inner_recomb, sizeof (double), nions, fptr);


    n += fwrite (plasmamain[m].scatters, sizeof (int), nions, fptr);
    n += fwrite (plasmamain[m].xscatters, sizeof (double), nions, fptr);

    n += fwrite (plasmamain[m].heat_ion, sizeof (double), nions, fptr);
    n += fwrite (plasmamain[m].cool_rr_ion, sizeof (double), nions, fptr);
    n += fwrite (plasmamain[m].cool_dr_ion, sizeof (double), nions, fptr);
    n += fwrite (plasmamain[m].lum_rr_ion, sizeof (double), nions, fptr);


    n += fwrite (plasmamain[m].levden, sizeof (double), nlte_levels, fptr);
    n += fwrite (plasmamain[m].recomb_simple, sizeof (double), nphot_total, fptr);
    n += fwrite (plasmamain[m].recomb_simple_upweight, sizeof (double), nphot_total, fptr);
    n += fwrite (plasmamain[m].kbf_use, sizeof (double), nphot_total, fptr);
  }

/* Now write out the macro atom info */

  if (geo.nmacro)
  {
    n += fwrite (macromain, sizeof (macro_dummy), NPLASMA, fptr);
    for (m = 0; m < NPLASMA; m++)
    {
      n += fwrite (macromain[m].jbar, sizeof (double), size_Jbar_est, fptr);
      n += fwrite (macromain[m].jbar_old, sizeof (double), size_Jbar_est, fptr);
      n += fwrite (macromain[m].gamma, sizeof (double), size_gamma_est, fptr);
      n += fwrite (macromain[m].gamma_old, sizeof (double), size_gamma_est, fptr);
      n += fwrite (macromain[m].gamma_e, sizeof (double), size_gamma_est, fptr);
      n += fwrite (macromain[m].gamma_e_old, sizeof (double), size_gamma_est, fptr);
      n += fwrite (macromain[m].alpha_st, sizeof (double), size_gamma_est, fptr);
      n += fwrite (macromain[m].alpha_st_old, sizeof (double), size_gamma_est, fptr);
      n += fwrite (macromain[m].alpha_st_e, sizeof (double), size_gamma_est, fptr);
      n += fwrite (macromain[m].alpha_st_e_old, sizeof (double), size_gamma_est, fptr);
      n += fwrite (macromain[m].recomb_sp, sizeof (double), size_alpha_est, fptr);
      n += fwrite (macromain[m].recomb_sp_e, sizeof (double), size_alpha_est, fptr);
      n += fwrite (macromain[m].matom_emiss, sizeof (double), nlevels_macro, fptr);
      n += fwrite (macromain[m].matom_abs, sizeof (double), nlevels_macro, fptr);

    }

  }

  fclose (fptr);

  Log_silent
    ("wind_write sizes: NPLASMA %d size_Jbar_est %d size_gamma_est %d size_alpha_est %d nlevels_macro %d\n",
     NPLASMA, size_Jbar_est, size_gamma_est, size_alpha_est, nlevels_macro);

  return (n);

}

/*

   wind_read (filename)

   History
	11dec	ksl	Updated so returns -1 if it cannot open the windsave file.  This
			was done to enable one to handle missing files differently in
			different cases
	14jul	nsh	Code added to read in variable length arrays in plasma structure
	15aug	ksl	Updated to read domain structure
	15oct	ksl	Updated to read disk and qdisk stuctures
*/


/**********************************************************/
/** 
 * @brief      Read back the windsavefile 
 *
 * @param [in] char  filename[]   The full name of the windsave file
 * @return     The number of successful reads, or -1 if the file cannot 
 * be opened
 *
 * @details
 * 
 * The routine reads in both the windsave file and the
 * associated atomic data files for a model. It also reads the
 * disk and qdisk structures.
 *
 *
 * ### Notes ###
 *
 * ### Programming Comment ### 
 * This routine calls wind_complete. This looks superfluous, since 
 * wind_complete and its subsidiary routines but it
 * also appears harmless.  ksl 
 *
 **********************************************************/

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
  Log ("Reading Windfile %s created with python version %s with python version %s\n", filename, version, VERSION);

  /* Now read in the geo structure */

  n += fread (&geo, sizeof (geo), 1, fptr);


  /* Read the atomic data file.  This is necessary to do here in order to establish the 
   * values for the dimensionality of some of the variable length structures, associated 
   * with macro atoms, especially but likely to be a good idea ovrall
   */

  get_atomic_data (geo.atomic_filename);


/* Now allocate space for the wind array */

  NDIM2 = geo.ndim2;
  NPLASMA = geo.nplasma;

  zdom = (DomainPtr) calloc (sizeof (domain_dummy), MaxDom);
  n += fread (zdom, sizeof (domain_dummy), geo.ndomain, fptr);

  calloc_wind (NDIM2);
  n += fread (wmain, sizeof (wind_dummy), NDIM2, fptr);

  /* Read the disk and qdisk structures */

  n += fread (&disk, sizeof (disk), 1, fptr);
  n += fread (&qdisk, sizeof (disk), 1, fptr);

  calloc_plasma (NPLASMA);

  n += fread (plasmamain, sizeof (plasma_dummy), NPLASMA, fptr);

  /*Allocate space for the dynamically allocated plasma arrays */

  calloc_dyn_plasma (NPLASMA);


  /* Read in the dynamically allocated plasma arrays */

  for (m = 0; m < NPLASMA; m++)
  {

    n += fread (plasmamain[m].density, sizeof (double), nions, fptr);
    n += fread (plasmamain[m].partition, sizeof (double), nions, fptr);

    n += fread (plasmamain[m].ioniz, sizeof (double), nions, fptr);
    n += fread (plasmamain[m].recomb, sizeof (double), nions, fptr);
    n += fread (plasmamain[m].inner_recomb, sizeof (double), nions, fptr);

    n += fread (plasmamain[m].scatters, sizeof (int), nions, fptr);
    n += fread (plasmamain[m].xscatters, sizeof (double), nions, fptr);

    n += fread (plasmamain[m].heat_ion, sizeof (double), nions, fptr);
    n += fread (plasmamain[m].cool_rr_ion, sizeof (double), nions, fptr);
    n += fread (plasmamain[m].cool_dr_ion, sizeof (double), nions, fptr);
    n += fread (plasmamain[m].lum_rr_ion, sizeof (double), nions, fptr);

    n += fread (plasmamain[m].levden, sizeof (double), nlte_levels, fptr);
    n += fread (plasmamain[m].recomb_simple, sizeof (double), nphot_total, fptr);
    n += fread (plasmamain[m].recomb_simple_upweight, sizeof (double), nphot_total, fptr);
    n += fread (plasmamain[m].kbf_use, sizeof (double), nphot_total, fptr);
  }


  /*Allocate space for macro-atoms and read in the data */

  if (geo.nmacro > 0)
  {
    calloc_macro (NPLASMA);
    n += fread (macromain, sizeof (macro_dummy), NPLASMA, fptr);
    calloc_estimators (NPLASMA);

    for (m = 0; m < NPLASMA; m++)
    {
      n += fread (macromain[m].jbar, sizeof (double), size_Jbar_est, fptr);
      n += fread (macromain[m].jbar_old, sizeof (double), size_Jbar_est, fptr);
      n += fread (macromain[m].gamma, sizeof (double), size_gamma_est, fptr);
      n += fread (macromain[m].gamma_old, sizeof (double), size_gamma_est, fptr);
      n += fread (macromain[m].gamma_e, sizeof (double), size_gamma_est, fptr);
      n += fread (macromain[m].gamma_e_old, sizeof (double), size_gamma_est, fptr);
      n += fread (macromain[m].alpha_st, sizeof (double), size_gamma_est, fptr);
      n += fread (macromain[m].alpha_st_old, sizeof (double), size_gamma_est, fptr);
      n += fread (macromain[m].alpha_st_e, sizeof (double), size_gamma_est, fptr);
      n += fread (macromain[m].alpha_st_e_old, sizeof (double), size_gamma_est, fptr);
      n += fread (macromain[m].recomb_sp, sizeof (double), size_alpha_est, fptr);
      n += fread (macromain[m].recomb_sp_e, sizeof (double), size_alpha_est, fptr);
      n += fread (macromain[m].matom_emiss, sizeof (double), nlevels_macro, fptr);
      n += fread (macromain[m].matom_abs, sizeof (double), nlevels_macro, fptr);

      /* Force recalculation of kpkt_rates */

      macromain[m].kpkt_rates_known = 0;
    }

  }

  fclose (fptr);

  wind_complete (wmain);

  Log ("Read geometry and wind structures from windsavefile %s\n", filename);

  return (n);

}




/**********************************************************/
/** 
 * @brief      A driver routine that calls coordinate-system specific routines
 * that complete the descirption of the wind
 *
 * @param [in] WindPtr w  The entire wind
 * @return     Always returns 0
 *
 * @details
 *
 * For the most point, the various routines that are called
 * just recalculate some of the various arrays used for 
 * finding the boundaries of an individual cell.
 *
 * ### Notes ###
 *
 * The need to call this routine whenever a windsave file is read
 * in could be obviated if all of the variable length arrays
 * in domains were actually written out.  But since these
 * variables are easy to calculate again, it is done here.
 * 
 *
 **********************************************************/

int
wind_complete (w)
     WindPtr w;
{
  int ndom;

  /* JM Loop over number of domains */


  for (ndom = 0; ndom < geo.ndomain; ndom++)
  {
    if (zdom[ndom].coord_type == SPHERICAL)
    {
      spherical_wind_complete (ndom, w);
    }
    else if (zdom[ndom].coord_type == CYLIND)
    {
      cylind_wind_complete (ndom, w);
    }
    else if (zdom[ndom].coord_type == RTHETA)
    {
      rtheta_wind_complete (ndom, w);
    }
    else if (zdom[ndom].coord_type == CYLVAR)
    {
      cylvar_wind_complete (ndom, w);
    }
    else
    {
      Error ("wind_complete: Don't know how to complete coord_type %d\n", zdom[ndom].coord_type);
      Exit (0);
    }

  }
  return (0);
}


/**********************************************************/
/** 
 * @brief      Save all the spectra to a binary file
 *
 * @param [in] char  filename[]   The file to write to 
 * @return     The number of successful writes
 *
 * @details
 *
 * ### Notes ###
 * 
 * The program exits if one cannot write the file
 *
 **********************************************************/

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
    Exit (0);
  }

  sprintf (line, "Version %s  nspectra %d\n", VERSION, nspectra);
  n = fwrite (line, sizeof (line), 1, fptr);
  n += fwrite (xxspec, sizeof (spectrum_dummy), nspectra, fptr);
  fclose (fptr);

  return (n);
}



/**********************************************************/
/** 
 * @brief      Read the binary file containing all the spectra
 *
 * @param [in] char  filename[]   The name of the file
 * @return     The number of spectra that were read in
 *
 * @details
 *
 * The routine allocates space for the spectra and reads
 * then in.
 *
 * ### Notes ###
 *
 * The first line of the file contains the Python version
 * and the number of spectra to be read in
 * 
 * The program exits if the file does not exist
 *
 **********************************************************/

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
    Exit (0);
  }

  n = fread (line, sizeof (line), 1, fptr);

  sscanf (line, "%*s %s %*s %d", version, &nspectra);
  Log ("Reading specfile %s with %d spectra created with python version %s with python version %s\n", filename, nspectra, version, VERSION);


  /* First allocate space */

  xxspec = calloc (sizeof (spectrum_dummy), nspectra);
  if (xxspec == NULL)
  {
    Error ("spectrum_init: Could not allocate memory for %d spectra with %d wavelengths\n", nspectra, NWAVE);
    Exit (0);
  }

/* Now read the rest of the file */

  n += fread (xxspec, sizeof (spectrum_dummy), nspectra, fptr);

  fclose (fptr);

  Log ("Read spec structures from specfile %s\n", filename);

  return (n);

}
