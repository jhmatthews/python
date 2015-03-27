


/**************************************************************************
                    Space Telescope Science Institute


  Synopsis:  

  These are routines which provide a very brief summary of 
  how the program is proceeding, mainly information about
  how many cycles have completed.  They are used for restarting
  python, primarily.

  Description:	

  Arguments:  


  Returns:

  Notes:

 

  History:
08nov	ksl	Coded as part of effort to be able to restart jobs
		on the Royal cluster

 ************************************************************************/

#include <stdio.h>
#include <strings.h>
#include <stdarg.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "log.h"

#include "atomic.h"
#include "python.h"


/* 

xsignal generates a single line message to a file names root.sig

All of the messages hae the format of

Mon Nov 10 09:05:34 2008     10.0  message 

where the message is determined by the format and the extra variables 
that are passed to the program.  This portion of the routine
operates in the same way that an fprintf statement operates
and was derived from the routines used for logging.


0811	ksl	Created as part of effort to make the program work on
		on Royal where there are time limits for how long a
		single process can run
 
*/

int
xsignal (char *root, char *format, ...)
{

  va_list ap, ap2;
  int result;

  char curtime[LINELENGTH];
  char message[LINELENGTH];
  FILE *fopen (), *sptr;
  char filename[LINELENGTH];
  double elapsed_time;


  /* if we are in MPI mode we only want the master process to write to sig file */
#ifdef MPI_ON
  if (rank_global==0)
    {
#endif


  /* Make the filemne */
  strcpy (filename, "");
  strcpy (filename, root);
  strcat (filename, ".sig");

  /* Open the file so that it will append if the file exists */


  if ((sptr = fopen (filename, "a")) == NULL)
    {
      Error ("xsignal: Could not even open signal file %s\n", filename);
      exit (0);
    }

  /* Now generate the message */

  /* Get the current time */
  get_time (curtime);


  elapsed_time = timer ();

  /* Get the time since the time was initiated */


  fprintf (sptr, "%s %8.1f ", curtime, elapsed_time);


  va_start (ap, format);
  va_copy (ap2, ap);		/* Added because vfprintf can change ap */
  result = vfprintf (sptr, format, ap);
  va_end (ap);


  vsprintf (message, format, ap2);
  Log ("xxx %s %8.1f %s", curtime, elapsed_time, message);


  fclose (sptr);

#ifdef MPI_ON
  }
#endif

  return (0);
}

/* 
 * rm the old signal file so one can begin again
 */

int
xsignal_rm (char *root)
{
  char filename[LINELENGTH];
  char command[LINELENGTH];

  /* Make the filemne */
  strcpy (filename, "");
  strcpy (filename, root);
  strcat (filename, ".sig");

  strcpy (command, "rm ");
  strcat (command, filename);
  system (command);

  return (0);

}


/* 
max_time is the amount of time in seconds that one wantts the program to run 
without stopping.  It can be updated at any point.

Note

check_time is the routine that halts the process if it runs too long.

The first call to the routine timer() sets the time.  It is fine to call this
directly but that is not necessary, as timeer will be called every time xsignal
is invoked.

0811	ksl	Created as part of efYort to make the program work on
		on Royal where there are time limits for how long a
		single process can run
*/


double max_time = -1.0;

int
set_max_time (char *root, double t)
{
  Log ("Setting maximum time for program to run to be %d s\n", t);
  xsignal (root, "MAXTIME set to %.0f\n", t);
  max_time = t;
  return 0;
}


/*

   check_time checks whether the elapsed time is greater than the max_time
   and terminates the program if that is the case.

   If the current time is greater than the allowed time. Then the program is
   stopped, and a signal is sent to the .sig file that the program can be 
   restarted

   If the maximum time has not been set then this routine is a NOP

Note:
	Generally speaking one should send a message to xsignal before
	invoking check_time to record the status of the signal file

0811	ksl	Created as part of efYort to make the program work on
		on Royal where there are time limits for how long a
		single process can run
*/

int
check_time (char *root)
{
  double t;
  if (max_time > 0.0 && (t = timer () > max_time))
    {
      error_summary ("Time allowed has expired expired\n");
      xsignal (root, "COMMENT max_time %.1f exceeded\n", max_time);
      exit (0);
    };

  return (0);
}
