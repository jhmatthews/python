/* diag.c contains a number of routines relating to extra diagnositics
   and related file i/o. 
   JM 1410 -- moved a few routines here relating to extra diagnostics
*/

/***********************************************************
                                       Southampton

 Synopsis:
	 open_diagfile sets up diagnostic files for use in python when the extra.diagnostics flag is 

Arguments:

	PhotPtr p;	the photon
	double ds	the distance the photon has travelled in the cell

Returns:
	Always returns 0.  .
 
Description:	 
Notes:
	

History:

	12jun 	nsh	72 Added lines to set up a file for outputting photons in given cells. The cells are read in from a file called diag_cells.dat. The existance of this file defined wether the diagnostic takes place.
	13jul	jm	changed print statements to logs and made more descriptive
	14oct   jm  changed to reflect move to 'advanced mode', see #111 and #120

**************************************************************/




#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "atomic.h"
#include "python.h"

int eplinit = 0;
int pstatinit = 0;		/*To say if we have checked to see if we need to log photons */

int
open_diagfile ()
{
  FILE *cellfile;		/*File that may or may not exist, pointing to cells we want to write out photon stats for */
  int cell;			/*Temporary storage of cell to use */

  if (eplinit == 0 && modes.save_extract_photons)
    {
      epltptr = fopen ("python.ext", "w");
      eplinit = 1;
    }

  ncstat = 0;			/*Zero the counter for the number of cells to be tracked */
  if (pstatinit == 0 && modes.save_cell_stats)		/* Check we havent already done this */
    {
      cellfile = fopen ("diag_cells.dat", "r");	/*This is the file containing cells to track */
      if (cellfile != NULL)	/*If there actually *is* a file read it */
	{
	  while (fscanf (cellfile, "%d", &cell) == 1)	/*If the line contains only one integer number read it in, otherwise quit reading */
	    {
	      Log ("open_diagfile: Cell diagnostics - we have a cell - %i, ncstat=%i, NCSTAT=%i\n", cell,
		      ncstat, NCSTAT);
	      if (-1 < cell && cell < geo.nplasma && ncstat < NCSTAT)	/*if the cells are real */
		{
		  Log ("open_diagfile: Cell numbers have been accepted as real.\n");
		  ncell_stats[ncstat] = cell;
		  ncstat = ncstat + 1;
		}
	      else
		{
		  Error
		    ("open_diagfile: %i is an unacceptable cell number for photon tracking\n",
		     cell);
		}
	    }
	  fclose (cellfile);
	  pstatptr = fopen ("cell_phot_stats.dat", "w");
	}
      else
	{
	  Log
	    ("open_diagfile: We have no file of cells to track, so we wont be doing any cell tracking\n");
	}
      pstatinit = 1;		/* We have initialised this routine */
    }

  // diag_on_off = 0;            // 0=off everything else is on
  return (0);
}






/***********************************************************
				University of Southampton

Synopsis:
	get_extra_diagnostics reads in extra diagnostics if 
	the user has asked for them. It uses rd_extra() in rdpar.c
	to get the actual lines, and then checks them against 
	hardwired options via a number of if loops.

Arguments:	
    none	

Returns:
    modifies the modes structure to turn on whichever modes
    are asked for.
 
Description:	
	
Notes:
    see #111 and #120

History:
    1410 -- JM -- Coded
**************************************************************/

int get_extra_diagnostics()
{
  int wordlength;
  char firstword[LINELENGTH];
  int noptions, rdstat;
  double answer;


  noptions = 0;
  wordlength = 0;
  rdstat = 0;

  if (modes.iadvanced == 0)
  	Error("Getting extra_diagnostics but advanced mode is off!\n");

  Log("get_extra_diagnostics: Getting extra diagnostics as requested...\n");

  while (noptions < NMAX_OPTIONS && rdstat == 0)
  {

    rdstat = rd_extra(firstword, &answer, &wordlength);


    /* if rdstat is 1 we reached EOF, so exit this routine */
    if (rdstat)
      {
      	if (noptions == 0)
      	  Error("get_extra_diagnostics: EOF: %i options read, but extra diagnostics on!\n", noptions);

        return(0);
      }

    Log("%s %8.4e\n", firstword, answer);

    /* would you like to save cell photon statistics */
	if (strncmp ("save_cell_statistics", firstword, wordlength) == 0)
	  {
	    modes.save_cell_stats = answer;
  	    Log("You are tracking photon statistics by cell\n");
	  }
	  
    /* would you like to use ispy mode */
	else if (strncmp ("ispymode", firstword, wordlength) == 0)
      {
	    modes.ispy = answer;
	    Log("ISPY mode is on\n");
	  }

	/* would you like to track resonant scatters */
	else if (strncmp ("track_resonant_scatters", firstword, wordlength) == 0)
      {
	    modes.track_resonant_scatters = answer;
	    Log("You are tracking resonant scatters\n");
	  }

    /* would you like keep windsave files for each cycle */
	else if (strncmp ("keep_ioncycle_windsaves", firstword, wordlength) == 0)
      {
	    modes.keep_ioncycle_windsaves = answer;
	    Log("You are keeping windsave files for each cycle\n");
	  }

	/* would you like to save data on extract */
	else if (strncmp ("save_extract_photons", firstword, wordlength) == 0)
      {
	    modes.save_extract_photons = answer;
	    Log("You are saving data on extract\n");
	  }

	/* would you like to print wind_rad_summary*/
	else if (strncmp ("print_windrad_summary", firstword, wordlength) == 0)
      {
	    modes.print_windrad_summary = answer;
	    Log("You are printing wind_rad_summary\n");
	  }

	/* would you like to print dvds_info */
	else if (strncmp ("print_dvds_info", firstword, wordlength) == 0)
      {
	    modes.print_dvds_info = answer;
	    Log("You are printing dvds_info\n");
	  }

    else
      {
      	Error("get_extra_diagnostics: didn't understand question %s, continuing!\n", firstword); 
      	noptions--;		// this isn't a real option, so decrement it before we increment it...!
      }

    noptions++;	// increment noptions
  }

  return 0;
}




/***********************************************************
                Southampton University

Synopsis: 
	save_photon_stats prints photon statistics to a file

Arguments:	
	One 		WindPtr for the cell
	p 			Photon pointer
	ds 			ds travelled

Returns:
 
Description:
   the loop below is if the user requires extra diagnostics and
   has provided a file diag_cells.dat to store photons stats for cells they have specified

Notes:
   Moved here to save duplicating code between bf_estimators_increment and radiation.

History:
   1410 JM 		Coding began
   1410 JM      Moved here from python.c	
 
**************************************************************/



int
save_photon_stats (one, p, ds, w_ave)
     WindPtr one;
     PhotPtr p;
     double ds,w_ave;
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
   geo.wcycle, p->np, p->freq, p->w, w_ave,ds, p->nscat, one->nplasma, one->nwind);
	}
    }
  return (0);
}
