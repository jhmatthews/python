
/***********************************************************/
/** @file  diag.c
 * @author ksl
 * @date   May, 2018
 *
 * @brief
 * This routine consolidates an approach to providing extra diagnostics
 * for debugging python. It also
 * contains initialization routines
 * that control which diagnosics are tracked.
 *
 * The first few functions are associated with gather information about what
 * one wants to track
 *
 * The basic idea for the extra diagnositics routines are as follows:
 *
 * There is one file to which all of the diagnostics are written.
 * Which diagnositica are written depends on how varaious of the modes
 * variables are set.
 *
 * All of the various diagnostic routines should be line oriented and begin
 * with a unique identifier to allow easy parsing of the outputs.
 *
 *
 * To add a new diagnostic one should
 *
 *  * (0) check that one of the existing diagnostics is insuffient for your needs
 *  * (1) add a routine here that writes to the diagnostic file.
 *  * (2) add a rdchoice statement in get_extra_diagnostics that allows you to track
 *  this diagnositic
 *  * (3) Assure that your new routine has the appropriate doxygen description,
 *  * (4) Add a yaml file for the new rdpar statement), and
 *  * (5) call the routine with the approriate condidtions from wherever you wish.
 *
 * As written, you must put your routine here.  This is done on purpose
 * so that we will not proliferate extra diagnostics around the various
 * subroutines.
 *
 *
 * ### Notes ###
 *
 * The unified approach makes most sense when one wants similar information
 * spread out over several routines.  There will be places where one
 * wishes to print out detailed information within a certain routine.
 * A current example is print_dvds_info, which is totally contained within
 * gradv.c
 *
 * This approached was adopted, based onthe discussion isssue #338
 * A detailed look at the subroutines indicates that this, especially the
 * part having to do with writing to a single diagnostic file is a work
 * in progress.
 *
 * Right now there are a number of very similar routines, some of
 * which might be combined.
 *
 * Note that as most of these routines are currently writen, one should
 * check whether the that appropriate flag is set, before calling 
 * the routine.  Otherwise, one may write to a file that has not been
 * opened.
 *
 *
 ***********************************************************/





#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "atomic.h"
#include "python.h"




/**********************************************************/
/**
 * @brief      get inputs that provides more control over how the program is
 *   run
 *
 * @return     Always returns 0
 *
 * @details
 * In advanced mode, this routine provides ways to change
 * several somewhat unrelated parameters in the routines
 *
 * The parameters are
 * * the fractional distance a photons may travel
 * * the lowest ion allowed to contribute to photabosrption
 * * whether photoabsorption is considered in the final spectrum
 *
 * ### Notes ###
 * @bug It is not obvious that much recent thought has been
 * given to the choices that are here.  The fractional distance
 * that a photon travel is intended to make sure the velocity
 * along the line of sight can be approximated linearly.  If
 * a photon travels too far in an azimuthal direction the sense
 * of the velocity can change and this prevensts this
 *
 * The lowest ion density contributing to photoionization is used
 * to determine what ions one has to calculate the photoionzation
 * xsection for.  The lower this density; the more that have to be
 * calculated, and as a result the slower the program.
 *
 * Keeping photoionizaion during final spectrum allows one to
 * check the contribution of photoabsorption.
 *
 **********************************************************/

int
get_standard_care_factors ()
{
  int istandard;
  char answer[LINELENGTH];
  istandard = 1;
  SMAX_FRAC = 0.5;
  DENSITY_PHOT_MIN = 1.e-10;

  /* 141116 - ksl - Made care factors and advanced command as this is clearly somethng that is diagnostic */

  if (modes.iadvanced)
  {
    strcpy (answer, "no");
    istandard = rdchoice ("@Diag.use_standard_care_factors(yes,no)", "1,0", answer);

    if (!istandard)
    {
      rddoub ("@Diag.fractional_distance_photon_may_travel", &SMAX_FRAC);
      rddoub ("@Diag.lowest_ion_density_for_photoabs", &DENSITY_PHOT_MIN);
      strcpy (answer, "no");
      modes.keep_photoabs = rdchoice ("@Diag.keep_photoabs_in_final_spectra(yes,no)", "1,0", answer);
    }
  }
  return (0);
}



/**********************************************************/
/**
 * @brief      Allows the user to specify what extra diagnositcs
 * to write out
 *
 * @return     Always returns 0
 *
 * Based on the inputs, the rougines modifies variables in the modes structure to turn on whichever modes
 * are asked for.
 *
 * @details
 * The routine is just for setting up the options.  The routine is called
 * only if one has requested extra diagnostics.
 *
 * ### Notes ###
 * see #111 and #120
 *
 * @bug This routine should be combined with get_standard_care_factors
 *
 **********************************************************/

int
get_extra_diagnostics ()
{
  char answer[LINELENGTH];
  if (modes.iadvanced == 0)
    Error ("Getting extra_diagnostics but advanced mode is off!\n");

  Log ("get_extra_diagnostics: Getting extra diagnostics as requested...\n");

  /* read the options. */
  strcpy (answer, "no");
  modes.save_cell_stats = rdchoice ("@Diag.save_cell_statistics(yes,no)", "1,0", answer);

  strcpy (answer, "no");
  modes.keep_ioncycle_windsaves = rdchoice ("@Diag.keep_ioncycle_windsaves(yes,no)", "1,0", answer);

  strcpy (answer, "no");
  modes.make_tables = rdchoice ("@Diag.make_ioncycle_tables(yes,no)", "1,0", answer);

  strcpy (answer, "no");
  modes.save_photons = rdchoice ("@Diag.save_photons(yes,no)", "1,0", answer);

  strcpy (answer, "no");
  modes.save_extract_photons = rdchoice ("@Diag.save_extract_photons(yes,no)", "1,0", answer);

  strcpy (answer, "no");
  modes.print_dvds_info = rdchoice ("@Diag.print_dvds_info(yes,no)", "1,0", answer);

  strcpy (answer, "no");
  modes.track_resonant_scatters = rdchoice ("@Diag.track_resonant_scatters(yes,no)", "1,0", answer);

  if (modes.save_cell_stats || modes.save_photons || modes.save_extract_photons | modes.track_resonant_scatters)
  {
    modes.extra_diagnostics = 1;
  }
  else
  {
    modes.extra_diagnostics = 0;
  }

  return 0;
}





int eplinit = 0;

/// To say if we have checked to see if we need to log photons
int pstatinit = 0;

/// Extra diagnostics file
FILE *epltptr;



/**********************************************************/
/**
 * @brief      opon a file for writing out
 *  extra diagnostics.  In some cases reads a file
 *  that specifies in which cells ones wants diagnostics
 *
 * @return     Always retuns 0
 *
 * @details
 * This routine gets some residual inforamtion needed
 * to specify exactlay what one wants to track, and opens
 * files that will be used to write the diagnostics.
 *
 * ### Notes ###
 * see #111 and #120
 *
 * The diagnostic filenames are  hardwired
 *
 * @bug Ultimately we would like to write the extra diagnositcs to
 * a single file
 *
 **********************************************************/

int
init_extra_diagnostics ()
{
  FILE *cellfile;               /*File that may or may not exist, pointing to cells we want to write out photon stats for */
  int cell;                     /*Temporary storage of cell to use */

  if (eplinit == 0 && modes.extra_diagnostics)
  {
    epltptr = fopen ("python.ext.txt", "w");
    eplinit = 1;
  }

  ncstat = 0;                   /*Zero the counter for the number of cells to be tracked */
  if (pstatinit == 0 && modes.save_cell_stats)  /* Check we havent already done this */
  {
    cellfile = fopen ("diag_cells.dat", "r");   /*This is the file containing cells to track */
    if (cellfile != NULL)       /*If there actually *is* a file read it */
    {
      while (fscanf (cellfile, "%d", &cell) == 1)       /*If the line contains only one integer number read it in, otherwise quit reading */
      {
        Log ("open_diagfile: Cell diagnostics - we have a cell - %i, ncstat=%i, NCSTAT=%i\n", cell, ncstat, NCSTAT);
        if (-1 < cell && cell < geo.nplasma && ncstat < NCSTAT) /*if the cells are real */
        {
          Log ("open_diagfile: Cell numbers have been accepted as real.\n");
          ncell_stats[ncstat] = cell;
          ncstat = ncstat + 1;
        }
        else
        {
          Error ("open_diagfile: %i is an unacceptable cell number for photon tracking\n", cell);
        }
      }
      fclose (cellfile);
      pstatptr = fopen ("cell_phot_stats.dat", "w");
    }
    else
    {
      Log ("open_diagfile: We have no file of cells to track, so we wont be doing any cell tracking\n");
    }
    pstatinit = 1;              /* We have initialised this routine */
  }

  return (0);
}




/**********************************************************/
/**
 * @brief      prints photon statistics to a file
 *
 * @param [in] WindPtr  one   WindPtr for the cell
 * @param [in] PhotPtr  p   Photon pointer (for one photon)
 * @param [in] double  ds   ds travelled
 * @param [in] double  w_ave   The average weight of the photon as it travelled
 * @return     Always returns 0
 *
 * @details
 * The routine checks whether this is one of a number of cells in which one
 * wishes to record the photons statitics and if so writes out information
 * about the photon to a file.
 *
 * ### Notes ###
 * Called from radiation
 *
 **********************************************************/

int
save_photon_stats (one, p, ds, w_ave)
     WindPtr one;
     PhotPtr p;
     double ds, w_ave;
{
  int i;

  /* JM -- 1310 -- the loop below is if the user requires extra diagnostics and
     has provided a file diag_cells.dat to store photons stats for cells they have specified
   */

  for (i = 0; i < ncstat; i++)
  {
    /* check if the cell is in the specified list - ncell_stats is global variable */
    if (one->nplasma == ncell_stats[i])
    {
      fprintf (pstatptr,
               "PHOTON_DETAILS cycle %3d n_photon %d freq %8.3e  w %8.3e ave_w %8.3e ds %8.3e nscat %d plasma cell %3d wind cell %3d\n",
               geo.wcycle, p->np, p->freq, p->w, w_ave, ds, p->nscat, one->nplasma, one->nwind);
    }
  }
  return (0);
}



/**********************************************************/
/**
 * @brief      saves informations about photons in
 *     a particulare wavelength range
 *
 * @param [in] int  n   The number of the spectrum
 * @param [in] PhotPtr  p   The photon before being doppler shifted
 * @param [in] PhotPtr  pp   The phtoon after being doppler shiftd
 * @param [in] double *  v   The velocity vector at the possibon of p/pp
 * @return     Always returns 0
phase(0=inferior_conjunction)                   0.5
 *
 * @details
 * This diagnostic routine was proably written to address a concern
 * about whether the photons were being properly doppler shifted on extract
 *
 * ### Notes ###
 * Called from extract.c, where currently one restricts the wavelenght
 * range of the photons to around CIV
 *
 **********************************************************/

int
save_extract_photons (n, p, pp, v)
     int n;
     PhotPtr p, pp;
     double *v;
{
  fprintf (epltptr,
           "EXTRACT %3d %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %7.2f %7.2f \n",
           n, p->x[0], p->x[1], p->x[2], v[0], v[1], v[2],
           p->lmn[0], p->lmn[1], p->lmn[2], pp->lmn[0], pp->lmn[1], pp->lmn[2], 2.997925e18 / p->freq, 2.997925e18 / pp->freq);

  return (0);
}

int save_photon_number = 0;


/**********************************************************/
/**
 * @brief      save_photon
 *
 * @param [in] PhotPtr  p   Photon pointer
 * @param [in] char  comment[]   A comment indicating why/at what point the photon information
 * needed to be recorded.
 * @return     Always returns 0
 *
 * @details
 * This is a diagnostic that allows one to print out information about
 * a photon at any time.
 *
 * ### Notes ###
 * It was written as part of the effort to
 *    debug imports for the fu_ori project.
 *
 **********************************************************/

int
save_photons (p, comment)
     PhotPtr p;
     char comment[];
{
  save_photon_number += 1;

  fprintf (epltptr,
           "PHOTON %3d %3d %10.4e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %3d %3d %3d %3d %s \n",
           geo.wcycle, p->np, p->freq, p->x[0], p->x[1], p->x[2], p->lmn[0], p->lmn[1],
           p->lmn[2], p->grid, p->istat, p->origin, p->nres, comment);

  return (0);
}


/**********************************************************/
/**
 * @brief
 *
 * @param [in] PhotPtr  p   Photon pointer
 * @param [in] int  nplasma   the cell in which the photon resides
 * @param [in] char *  comment   An arbitrary comment
 * @return     Always returns 0
 *
 * @details
 * Another way to track photons. This one is focussed on what happens
 * scattering events.
 *
 * ### Notes ###
 *
 **********************************************************/

int
track_scatters (p, nplasma, comment)
     PhotPtr p;
     int nplasma;
     char *comment;
{

  fprintf (epltptr, "Scattter %i %.2e %.2e %.2e  %i %e %e %i %s\n", p->np,
           p->x[0], p->x[1], p->x[2], p->grid, p->freq, p->w, nplasma, comment);

  return (0);
}
