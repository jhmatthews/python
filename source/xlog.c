/***********************************************************
                                       Space Telescope Science Institute

 Synopsis:
	These are a simple series of routines designed to store comments and errors
	in a diagnostic file.
		
		int Log_init(filename)				Open a logfile
		int Log_append(filename)			Append to a logfile

		int Shout( char *format, ...)			Send a message to the screen and logfile that cannot be
								suppressed.

		int Log ( char *format, ...)			Send a message to the screen and logfile
		int Log_silent ( char *format, ...)		Send a message to the logfile


		int Error ( char *format, ...)			Send a message prefaced by the word "Error:" to
										the screen and to the logfile.
		int Error_silent ( char *format, ...)		Send a message prefaced by the word "Error:" to
										the logfile
		
		int Log_close()					Close the current logfile

		int error_summary(char *format)			Summarize all of the erors that have been
								logged to this point in time
		int Log_set_verbosity(vlevel)			Set the verbosity of the what is printed to
								the cren and the log file

		int Log_print_max(print_max)			Set the number of times a single error will be
								output to the scren and the log file

   		int Log_flush()					simply flushes the logfile to disk


                int Log_set_mpi_rank(rank, n_mpi)		Tells kpar the rank of the parallel process in parallel mode,
								and divides max errors by n_mpi

		int Log_parallel(message)			Log statement for parallel reporting

    JM- Debug statements are controlled by verbosity now, so no need for Log_Debug
		int Log_debug(value)				Turn on logging of Debug statements if value is non-zero

		int Debug( char *format, ...) 			Log an statement to the screen.  This is essentially a 
								intended to replace a printf statement in situations where
								one is debugging code.  The use of Debug instead of log
								means that a future developer should be free to remove the
								Debug statement from the code


Arguments:		


Returns:
 
Description:	
	
	Normally, one would begin by issuing a Log_init (or Log_append) command.  This will open a
	file which will be used for logging.  If the user does not issue a Log_init command
	but attempts to use one of the other routines in the file, then Log_init will
	be issued internally but the and the which will used for logging will be called
	"logfile".
	
	All of the Log... and Error... allow one to send what is essentially fprintf and/if desired
	printf commands.  They are designed to handle variable numbers of arguments.  So
	for example if you want to send a message to the log file only, which includes the
	variables i, and j, one would say:
		
		Log_silent("This message writes i (%d) and j (%d) to the screen\n",i,j);
		
	Log_close will close the existing log file.  If one actually uses Log_close, one should
	be sure not to try and reopen the same file for logging because this will overwrite
	the original file.

	What is actually recorded depends on the verbosity, which is currently set to 5.  Lower levels
	of verbosity will print or record less.  In particular:

	The longer term plan is to introduce something that will always be recorded.  
	
	

Notes:
	These is used more generally than in python.  

History:
  	98feb	ksl	Coding of these subroutines began.
	99dec	ksl	Rewrote sane_check for linux using routine finite
	06nov	ksl	Added error logging
	07aug	ksl	Added a total error count to cause the program to die
			I had a situation that the number of times an error
			occured seemed to get to be so many that the
			integer values wrapped.
	08nov	ksl	Added a new Log_append command in order to allow
			restarts.
	09feb	ksl	Added a new capability to control the verbosity of what 
			is printed or recorded
	09mar	ksl	Added capability to change the number of times the same
			error is logged to the screen, and to cause the 
			program to quit.
	12dec	nsh	Added Log_flush() which flushed the logfile to disk. This
			is to help in following jobs on the Soton Iridis cluster 
			which seems to have a very large write buffer meaning 	
			nothing gets written out for a long time
	12dec	nsh	Added some lines to error - since ap not preserved, and 
			error logging to diag file was giving crazy results.
	13jul	jm	Added Log_set_mpi_rank and Log_parallel statements for parallel 
			reporting. Also added new verbosity setting, LOG_PARALLEL, and 
			divided max errors by n_mpi processors which is communicated
			via log_set_mpi_rank. The verbosity and rank are also communicated
			to rdpar from this thread.
	13sep	nsh	Added a new class of reporting - warnings. Things we would like to 
			know about (no photons in band, no model in band) but do not want 
			to crash the code.
  1407 JM removed warning - we would like to throw errors
  1411 JM debug statements are controlled by verbosity now, 
          so no need for Log_Debug

 
**************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <math.h>
#include "log.h"

#define LINELENGTH 132
#define NERROR_MAX 500          // Number of different errors that are recorded

/* definitions of what is logged at what verboisty level */

#define SHOW_PARALLEL		1
#define SHOW_LOG  		2
#define SHOW_ERROR		2
#define SHOW_DEBUG	  	4
#define SHOW_LOG_SILENT  	5
#define SHOW_ERROR_SILENT	5


int my_rank = 0;                // rank of mpi process, set to zero

int log_print_max = 100;        // Maximum number of times a single error will be reported.  
                                // Note that it will still be counted.
int max_errors = 100000;        // Maximum number of times an error can occur before giving up

typedef struct error_log
{
  char description[LINELENGTH];
  int n;
} error_dummy, *ErrorPtr;

ErrorPtr errorlog;


int nerrors;

FILE *diagptr;
int init_log = 0;
int log_verbosity = 5;          // A parameter which can be used to suppress what would normally be logged or printed

//1411 JM debug statements are controlled by verbosity now, so no need for Log_Debug
//int log_debug=0;      //A parameter which is set to cause Debug commands to be logged

int
Log_init (filename)
     char *filename;
{
  FILE *fopen ();

  if ((diagptr = fopen (filename, "w")) == NULL)
  {
    printf ("Yikes: could not even open log file %s\n", filename);
    exit (0);
  }
  init_log = 1;

  nerrors = 0;
  errorlog = (ErrorPtr) calloc (sizeof (error_dummy), NERROR_MAX);

  if (errorlog == NULL)
  {
    printf ("There is a problem in allocating memory for the errorlog structure\n");
    exit (0);
  }

  return (0);
}


int
Log_append (filename)
     char *filename;
{
  FILE *fopen ();

  if ((diagptr = fopen (filename, "a")) == NULL)
  {
    printf ("Yikes: could not even open log file %s\n", filename);
    exit (0);
  }
  init_log = 1;

  nerrors = 0;
  errorlog = (ErrorPtr) calloc (sizeof (error_dummy), NERROR_MAX);

  if (errorlog == NULL)
  {
    printf ("There is a problem in allocating memory for the errorlog structure\n");
    exit (0);
  }

  return (0);
}

int
Log_close ()
{
  fclose (diagptr);
  init_log = 0;
  free (errorlog);              // Release the error summary structure
  return (0);
}

/* The next routine allows the user to change the amount of reporting that is
 * carried out through the logging subroutines
 */

int
Log_set_verbosity (vlevel)
     int vlevel;
{
  log_verbosity = vlevel;
  rdpar_set_verbose (vlevel);
  return (0);
}


/* The next routine allows the user to change the amount of reporting that is
 * carried out through the logging subroutines
 */

int
Log_print_max (print_max)
     int print_max;
{
  log_print_max = print_max;
  return (0);
}


/* The next routine allows the user to change the amount of reporting that is
 * carried out through the logging subroutines
 */

int
Log_quit_after_n_errors (n)
     int n;
{
  max_errors = n;
  return (0);
}

int
Log (char *format, ...)
{
  va_list ap, ap2;
  int result;

  if (init_log == 0)
    Log_init ("logfile");

  if (log_verbosity < SHOW_LOG)
    return (0);

  va_start (ap, format);
  va_copy (ap2, ap);            /* ap is not necessarily preserved by vprintf */

  if (my_rank == 0)
    result = vprintf (format, ap);
  result = vfprintf (diagptr, format, ap2);
  va_end (ap);
  return (result);
}

int
Log_silent (char *format, ...)
{
  va_list ap, ap2;
  int result;

  if (init_log == 0)
    Log_init ("logfile");

  if (log_verbosity < SHOW_LOG_SILENT)
    return (0);
  va_start (ap, format);
  va_copy (ap2, ap);            /* ap is not necessarily preserved by vprintf */

  result = vfprintf (diagptr, format, ap);
  va_end (ap);
  return (result);
}

int
Error (char *format, ...)
{
  va_list ap, ap2;
  int result;

  if (init_log == 0)
    Log_init ("logfile");

  if (error_count (format) > log_print_max || log_verbosity < SHOW_ERROR)
    return (0);

  va_start (ap, format);
  va_copy (ap2, ap);            /*NSH 121212 - Line added to allow error logging to work */
  if (my_rank == 0)             // only want to print errors if master thread
    result = vprintf (format, ap);

  fprintf (diagptr, "Error: ");
  result = vfprintf (diagptr, format, ap2);
  va_end (ap);
  return (result);
}



int
Error_silent (char *format, ...)
{
  va_list ap, ap2;
  int result;
  if (init_log == 0)
    Log_init ("logfile");

  if (error_count (format) > log_print_max || log_verbosity < SHOW_ERROR_SILENT)
    return (0);

  va_start (ap, format);
  va_copy (ap2, ap);            /* ap is not necessarily preserved by vprintf */

  if (my_rank == 0)             // only want to print errors if master thread
    result = vprintf (format, ap);
  fprintf (diagptr, "Error: ");
  result = vfprintf (diagptr, format, ap2);
  va_end (ap);
  return (result);
}



int
Shout (char *format, ...)
{
  va_list ap, ap2;
  int result;

  if (init_log == 0)
    Log_init ("logfile");

  if (error_count (format) > log_print_max)
    return (0);

  printf ("Error: ");
  va_start (ap, format);
  va_copy (ap2, ap);            /* ap is not necessarily preserved by vprintf */
  result = vprintf (format, ap);
  fprintf (diagptr, "Error: ");
  result = vfprintf (diagptr, format, ap);
  va_end (ap);
  return (result);
}

int
sane_check (x)
     double x;
{
  int i;
  if ((i = isfinite (x)) == 0)
  {
    Error ("sane_check: %d %e\n", i, x);
    return (-1);
  }
  return (0);
}


int
error_count (char *format)
{
  int n;
  n = 0;
  while (n < nerrors)
  {
    if (strcmp (errorlog[n].description, (format)) == 0)
      break;
    n++;
  }

  if (n == nerrors)
  {
    strcpy (errorlog[nerrors].description, format);
    errorlog[n].n = 1;
    if (nerrors < NERROR_MAX)
    {
      nerrors++;
    }
    else
    {
      printf ("Exceeded number of different errors that can be stored\n");
      error_summary ("Quitting because there are too many differnt types of errors\n");
      exit (0);
    }
  }
  else
  {
    n = errorlog[n].n++;
    if (n == log_print_max)
      Error ("error_count: This error will no longer be logged: %s\n", format);
    if (n == max_errors)
    {
      error_summary ("Something is drastically wrong for any error to occur so much!\n");
      exit (0);
    }
  }
  return (n + 1);
}



int
error_summary (message)
     char *message;
{
  int n;
  Log ("\nError summary: %s\n", message);
  Log ("Recurrences --  Description\n");
  for (n = 0; n < nerrors; n++)
  {
    Log ("%9d -- %s", errorlog[n].n, errorlog[n].description);
  }

  return (0);
}


/*NSH 121107 added a routine to flush the diagfile*/

int
Log_flush ()
{
  if (init_log == 0)
    Log_init ("logfile");

  fflush (diagptr);
  return (0);
}


/* JM130717 the next set of routines are designed to clean up the code 
 * in parallel mode
 */

/* The next routine simply sets the rank of the process 
 * if not in parallel mode then we set my_rank to zero
 */

int
Log_set_mpi_rank (rank, n_mpi)
     int rank, n_mpi;
{
  my_rank = rank;
  rdpar_set_mpi_rank (rank);    //this just communicates the rank to rdpar      

  /* if in parallel mode we divide by the number of parallel processes for max errors */
  //     max_errors /= n_mpi;                
  //log_print_max /= n_mpi; for the moment we won'y change this as only master thread prints errors anyway
  return (0);
}


/* log statement that prints parallel statements to screen */

int
Log_parallel (char *format, ...)
{
  va_list ap, ap2;
  int result;

  if (log_verbosity < SHOW_PARALLEL)
    return (0);

  va_start (ap, format);
  va_copy (ap2, ap);            /* ap is not necessarily preserved by vprintf */

  if (my_rank == 0)
    result = vprintf (format, ap);

  fprintf (diagptr, "Para: ");
  result = vfprintf (diagptr, format, ap2);

  return (result);
}


/* Set a flag to which cause Debug statements to be logged */
// int Log_debug(value)
//      int value;
// {
//      log_debug=value;
//      return(0);
// }
/* JM- Debug statements are controlled by verbosity now, so no need for Log_Debug */
/* Log a debug statement if the external varialbe log_debug
 * has been set to a non-zero value
 */

int
Debug (char *format, ...)
{
  va_list ap, ap2;
  int result;

  if (log_verbosity < SHOW_DEBUG)
    return (0);

  if (init_log == 0)
    Log_init ("logfile");

  va_start (ap, format);
  va_copy (ap2, ap);            /*NSH 121212 - Line added to allow error logging to work */
  if (my_rank == 0)             // only want to print errors if master thread
    vprintf ("Debug: ", ap);
  result = vprintf (format, ap);
  fprintf (diagptr, "Debug: ");
  result = vfprintf (diagptr, format, ap2);
  va_end (ap);
  return (result);
}
