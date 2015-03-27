
#include <stdio.h>
#include <strings.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <sys/time.h>
#include <time.h>


/*
Return the time in seconds since the timer was initiated
*/

double t_init = 0.0;

double
timer ()
{
  struct timeval tv;
  struct timezone tz;
  double t;

  gettimeofday (&tv, &tz);
  t = tv.tv_sec + 1.e-6 * tv.tv_usec;
  if (t_init == 0.0)
    {
      t_init = t;
      return (0.0);
    }
  else
    {
      return (t - t_init);
    }
}


/* 
Get the current time.  as an ascii string

Notes:

ctime returns a string with a trailing \n, which we need
to strip off.  This accounts for the rather bizarre handling
of the string.


0811	ksl	67 - As part of effort to make python more 
		friendly to running on Royal

 */

int
get_time (curtime)
     char curtime[];
{
  time_t tloc;
  time (&tloc);
  strcpy (curtime, "");
  strncpy (curtime, ctime (&tloc), 24);
  curtime[24] = '\0';		// We need to end the string properly
  return (0);
}
