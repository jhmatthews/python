
/***********************************************************/
/** @file  windsave2table.c
 * @author ksl
 * @date   April, 2018
 *
 * @brief  A standalone routine for writing a standard set of data
 * from a windsavefile to ascii tables all readable with as astropy.io.ascii
 *
 * This routine is run from the command line, as follows
 *
 * windsave2table  windsave_root
 *
 * where windsave_root is the root name for a python run, or more precisel
 * the rootname of a windsave file.
 *
 * The routine reads the windsavefile and then writes out a selected 
 * set of variables into a variety of number of files, each of which
 * are readable as astropy tables, where each row corresponds to one
 * element of the wind.   
 *
 * All of the files begin
 * with the rootname and the domain number.
 * If Python run that created the windsave file
 * has multiple domains a sepearate output file is created for each
 * domain
 *
 * All of the files, begin with
 *
 *        x        z    i    j inwind 
 *
 * where x and z correspond to the position of the cell, i and j correspond to 
 * the cell number and inwind says whether the cell was in the wind, so that it is
 * fairly straightforward to create plotting routines for various parameters
 * that are contained in the remaining columns.  
 *
 * The files include a so-called masterfile which records basic parameters
 * like n_e, velocity,  rho, t_e, t_r, and ionization fractions
 * for a few key ions in each cell when the windsave file was written
 *
 * ### Notes ###
 *
 * Whereas py_wind is intended to be run interactively, windsave2table is
 * entirely hardwired so that it produces a standard set of output
 * files.  To change the outputs one has to modify the routine
 *
 * This file just contains the driving routine.  All of the 
 * real work is carried out in windsave2table_sub.c  Indeed so 
 * little of the real work is done here that it might be sensible
 * to move some portion of that code here.
 *
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
 * @brief      windsave2table writes key variables in a windsave file 
 * to an astropy table calculated Python.  This is the  main routine.
 *
 * @param [in] int  argc   The number of argments in the command line
 * @param [in] char *  argv[]   The command line
 * @return     Always returns 0, unless the windsave file is not
 * found in which case the routine will issue and error before
 * exiting.  
 *
 * @details
 * argc and argv[] are the standard variables provided to main 
 * in a c-program.  Only the first command line argument is 
 * parsed and this should be rootname of the windsave file
 *
 * ### Notes ###
 *
 * This routine is a supervisory routine. The real work 
 * is in do_windsave2table
 *
 * The routine does not read the .pf file.  It does read  
 * the windsave file and the associated atomic data file
 * before calling do_windsave2table.
 *
 **********************************************************/

int
main (argc, argv)
     int argc;
     char *argv[];
{



  char root[LINELENGTH], input[LINELENGTH];
  char outputfile[LINELENGTH];
  char windsavefile[LINELENGTH];
  char parameter_file[LINELENGTH];
  int create_master_table (), create_ion_table ();
  int do_windsave2table ();



  strcpy (parameter_file, "NONE");

  /* Next command stops Debug statements printing out in py_wind */
  Log_set_verbosity (3);

  if (argc == 1)
  {
    printf ("Root for wind file :");
    fgets (input, LINELENGTH, stdin);
    get_root (root, input);
  }
  else
  {
    strcpy (input, argv[argc - 1]);
    get_root (root, input);
  }


  printf ("Reading data from file %s\n", root);

  /* Now create the names of all the files which will be written */

  strcpy (windsavefile, root);
  strcpy (outputfile, root);

  strcat (windsavefile, ".wind_save");
  strcat (outputfile, ".txt");


/* Read in the wind file */

  if (wind_read (windsavefile) < 0)
  {
    Error ("py_wind: Could not open %s", windsavefile);
    exit (0);
  }


  printf ("Read wind_file %s\n", windsavefile);

  get_atomic_data (geo.atomic_filename);

  printf ("Read Atomic data from %s\n", geo.atomic_filename);

  do_windsave2table (root);

  return (0);
}
