
/**************************************************************************
                    Space Telescope Science Institute


  Synopsis:

  Description:	This subroutine retrieves a line from a previously
	opened file.  Lines beginnig with # are assumed to be comments
	and are ignored.  Blank lines are also ingnored. 

  Arguments:		

  Returns:
	0   after a normal read.
	EOF at the end of file 

  Notes:
	The real reason for this routine is to provide a common way to
	incorporate comments into ascii files, since if that were not
	the case one could use fgets and then sscanf.

	sscanf would then be used to parse the file.
	
	If the file is not opened properly one currently gets a segmentation
	error (on linux/gcc).  fgets is supposed to return EOF when a file
	is not open, but his did not happen.  An alternative would be to
	move the opening of the file into this routine, but I have not
	done this.

  History:
	99aug7	ksl	Began work
	06jul12	ksl	Increased LINELEN to 256
	

 ************************************************************************/
#include	<math.h>
#include	<stdio.h>
#include	<string.h>
#include "log.h"

#define LINELEN	256

FILE *get_line_ptr;

int
get_line (fptr, line)
     FILE *fptr;
     char line[];
{

  char dummy[LINELEN];
  int j;


a:if (fgets (line, LINELEN, fptr) == NULL)
    {
      return (EOF);
    }
  j = sscanf (line, "%s", dummy);
  if (dummy[0] == '#' || j == 0)
    {
      goto a;
    }


  return (0);

}
