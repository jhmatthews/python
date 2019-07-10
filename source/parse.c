
/***********************************************************/
/** @file  parse.c
 * @author ksl
 * @date   February, 2018
 *
 * @brief  Routines associated with parsing the command line
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
 * @brief      parses the command line options
 *
 * @param [in]  int  argc   the number of command line arguments
 * @param [in]  char *  argv[]   The command line arguments
 * @return      restart_stat   1 if restarting a previous model,
 * 0 in all other cases.
 *
 * Python has a fairly rich set of command line options, which
 * are parsed by this routine
 *
 * The routine also creates the diag folder, which is where most
 * of the log files are written 
 *
 * ###Notes###
 *
 * The general purpose of each of the command line options
 * should be fairly obvious from reading the code.
 *
 * If changes to the command line interface are made they should
 * be described in the routine help 
 *
 * Although this routine uses the standard Log and Error commands
 * the diag files have not been created yet and so this information
 * is really simply written to the terminal.  
 *
 **********************************************************/

int
parse_command_line (argc, argv)
     int argc;
     char *argv[];
{
  int restart_stat, verbosity, max_errors, i;
  int j = 0;
  char dummy[LINELENGTH];
  int mkdir ();
  double time_max;
  char *fgets_rc;

  restart_stat = 0;

  if (argc == 1)
  {
    printf ("Parameter file name (e.g. my_model.pf, or just my_model):");
    fgets_rc = fgets (dummy, LINELENGTH, stdin);
    if (!fgets_rc)
    {
      Error ("Input rootname is NULL or invalid\n");
      exit (1);
    }
    get_root (files.root, dummy);
    strcpy (files.diag, files.root);
    strcat (files.diag, ".diag");
  }
  else
  {

    for (i = 1; i < argc; i++)
    {
      if (strcmp (argv[i], "-h") == 0)
      {
        help ();
      }
      else if (strcmp (argv[i], "-r") == 0)
      {
        Log ("Restarting %s\n", files.root);
        restart_stat = 1;
        j = i;
      }
      else if (strcmp (argv[i], "-t") == 0)
      {
        if (sscanf (argv[i + 1], "%lf", &time_max) != 1)
        {
          Error ("python: Expected time after -t switch\n");
          exit (1);
        }
        set_max_time (files.root, time_max);
        i++;
        j = i;
        Log ("Program will stop after time %f\n", time_max);

      }
      else if (strcmp (argv[i], "-v") == 0)
      {
        if (sscanf (argv[i + 1], "%d", &verbosity) != 1)
        {
          Error ("python: Expected verbosity after -v switch\n");
          exit (0);
        }
        Log_set_verbosity (verbosity);
        i++;
        j = i;
        Log ("Verbosity level set to %d\n", verbosity);

      }
      else if (strcmp (argv[i], "-e") == 0)
      {
        if (sscanf (argv[i + 1], "%d", &max_errors) != 1)
        {
          Error ("python: Expected max errors after -e switch\n");
          exit (1);
        }
        Log_quit_after_n_errors (max_errors);
        i++;
        j = i;
        Log ("Setting the maximum number of errors of a type before quitting to %d\n", max_errors);

      }
      else if (strcmp (argv[i], "-e_write") == 0)
      {
        if (sscanf (argv[i + 1], "%d", &max_errors) != 1)
        {
          Error ("python: Expected max errors after -e switch\n");
          exit (1);
        }
        Log_print_max (max_errors);
        i++;
        j = i;
        Log ("Setting the maximum number of errors of a type to print out to  %d\n", max_errors);

      }
      else if (strcmp (argv[i], "-d") == 0)
      {
        modes.iadvanced = 1;
        Log ("Enabling advanced/diagnositic inputs (@ commands)\n");
        j = i;
      }
      else if (strcmp (argv[i], "-f") == 0)
      {
        modes.fixed_temp = 1;
        Log ("Invoking fixed temperature mode\n");
        j = i;
      }

      else if (strcmp (argv[i], "--rseed") == 0)
      {
        modes.rand_seed_usetime = 1;
        j = i;
        Log ("Using a random seed in random number generator\n");
      }
      else if (strcmp (argv[i], "-z") == 0)
      {
        modes.zeus_connect = 1;
        Log ("Setting zeus_connect to %i\n", modes.zeus_connect);
        j = i;
      }
      else if (strcmp (argv[i], "-i") == 0)
      {
        modes.quit_after_inputs = 1;
        j = i;
      }
      else if (strcmp (argv[i], "-p") == 0)
      {
        Log ("Logarithmic photon stepping enabled\n");
        modes.photon_speedup = 1;
        if (sscanf (argv[i + 1], "%lf", &PHOT_RANGE) == 1)
        {
          i++;
        }
        else
        {
          PHOT_RANGE = 1.;
        }
        j = i;
      }
      else if (strcmp (argv[i], "--dry-run") == 0)
      {
        modes.quit_after_inputs = 1;
        j = i;
      }

      else if (strcmp (argv[i], "--version") == 0)
      {
        /* give information about the pyhon version, such as commit hash */
        Log ("Python Version %s \n", VERSION);  //54f -- ksl -- Now read from version.h
        Log ("Built from git commit hash %s\n", GIT_COMMIT_HASH);
        /* warn the user if there are uncommited changes */
        int git_diff_status = GIT_DIFF_STATUS;
        if (git_diff_status > 0)
          Log ("This version was compiled with %i files with uncommitted changes.\n", git_diff_status);
        exit (1);
      }

      else if (strncmp (argv[i], "-", 1) == 0)
      {
        Error ("python: Unknown switch %s\n", argv[i]);
        help ();
      }
    }

    /* The last command line variable is always the .pf file */

    if (j + 1 == argc)
    {
      Error ("All of the command line has been consumed without specifying a parameter file name, so exiting\n");
      exit (1);
    }


    strcpy (dummy, argv[argc - 1]);
    get_root (files.root, dummy);

    /* This completes the parsing of the command line */

    /* Create a subdirectory to store diaganostic files */

    sprintf (files.diagfolder, "diag_%s/", files.root);
    mkdir (files.diagfolder, 0777);
    strcpy (files.diag, files.diagfolder);
    sprintf (dummy, "_%d.diag", rank_global);
    strcat (files.diag, files.root);
    strcat (files.diag, dummy);


  }

  return (restart_stat);
}




/**********************************************************/
/** 
 * @brief      print out some basic help concering command line options for the program
 *
 * @return     0 in all cases
 *
 * This simply prints out a hardwired multi-line string to the command line.
 *
 * ###Notes###
 *
 * An improved version of this routine would simply read and print and external 
 * file.
 *
 * If one chooses to add to what is printed out then one needs to be careful 
 * concerning the ends of lines.
 *
 **********************************************************/

void
help ()
{
  char *some_help;

  some_help = "\
\n\
This program simulates radiative transfer in a (biconical) CV, YSO, quasar or (spherical) stellar wind \n\
\n\
Usage:  py [-h] [-r] [-t time_max] [-v n] [--dry-run] [-i] [--version] [--rseed] [-p n_steps] xxx  or simply py \n\
\n\
where xxx is the rootname or full name of a parameter file, e. g. test.pf \n\
\n\
and the switches have the following meanings \n\
\n\
 -h             Print this help message and stop \n\
 -r             Restart a run of the progarm reading the file xxx.windsave \n\
 -t time_max    Limit the total time to approximately time_max seconds.  Note that the program checks \n\
                for this limit somewhat infrequently, usually at the ends of cycles, because it \n\
                is attempting to save the program outputs so that the program can be restarted with \n\
                -r if that is desired. \n\
 -v n           Increase or decrease the amount of print out.  The default is 4.  Larger numbers increase  \n\
                the amount printed; smaller numbers decrease it.   \n\
 --dry-run      Create a new .pf file and stop \n\
 -i             Same as --dry-run \n\
 --version      Print out python version, commit hash and if there were files with uncommitted \n\
                changes and stop \n\
 --rseed        Set the random number seed to be time-based, rather than fixed. \n\
\n\
Other switches exist but these are not intended for the general user.\n\
These are largely diagnostic or for special cases. These include\n\
 -d             Enable advanced/diagnostic inputs (normally for debugging purposes) \n\
                Python will then query the user for information about what to do with a series of \n\
                inputs beginning with @ \n\
 -e             Change the maximum number of errors before the program will quit\n\
 -e_write 	Change the maximum number of errors to print out before recording errors silently\n\
 -f             Invoke a fixed temperature mode, used for runs with Zeus \n\
 -z             Invoke a special mode for that causes Python to start with a run from Zeus\n\
 -p range       Invoke the photon logarithmic stepping algorithm which in some cases can result in a speed up\n\
                Range is in powers of 10, the difference beween the number of photons in the first cycle \n\
                compared to the last \n\
\n\
If one simply types py or pyZZ where ZZ is the version number, one is queried for a name \n\
of the parameter file and inputs will be requested from the command line. \n\
\n\
\n\
";                              // End of string to provide one with help

  printf ("%s\n", some_help);

  exit (0);                     // Note that here we simply do want to exit, not use Exit
}
