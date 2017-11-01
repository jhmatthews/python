


/**************************************************************************
                    Space Telescope Science Institute


  Synopsis:  

  This is a routine which is intended to help with program updates that
  change the name of a keyword in a parameter file.  It is called from
  rdpar.string_process_from_file

  Description:	

  The routine simply matches the string in new_question to a string in
  the array new_names, below.  If it finds a match, the old_name is 
  returned in old_question.



  Arguments:  	


  Returns:

  	0 if there was no match
	1 if there was a match of a name in the parameter file to
	  one of the new_names, and in old_question the name of
	  the old_value

  Notes:

  To add another variable to the list, just record the new_name and the
  old_name in the arrays below, and increase the number of names by 1

  Do not include the material that is in paren, that is for

  xxx(many_choices) -->  xxxx

  The routine is completely hardwired as wrtten though clearly this 
  could be changed.
  	


  History:
	16sept	ksl	Coded as part of an effort to make python more 
			robust to changes in parameter file names

 ************************************************************************/

#include <stdio.h>
#include <strings.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "log.h"

#define	LINELEN 132


char *old_names[] = { "mstar", "rstar","Disk.illumination.treatment" };
char *new_names[] = { "Central.object.mass", "Central.object.radius","Surface.reflection.or.absorption" };

int nunber_of_names = 3;



int
check_synonyms (new_question, old_question)
     char new_question[], old_question[];
{
  int n;
  char *line;
  char firstword[LINELEN];
  int nwords, wordlength;
  char *ccc, *index ();


// Strip off any extra bits in the new question
  line = new_question;
  strcpy (firstword, "");
  nwords = sscanf (line, "%s ", firstword);
  wordlength = strlen (firstword);


  if (nwords == 0)
  {
    return (0);
  }


// Strip off everthing in paren, or rather find out how much of the string to comapre
  if ((ccc = index (firstword, '(')) != NULL)
  {
    wordlength = (int) (ccc - firstword);
    if (wordlength == 0)
      return (0);
  }



  for (n = 0; n < nunber_of_names; n++)
  {
    if (strncmp (new_names[n], firstword, wordlength) == 0)
    {
      Log ("Matched keyword %s in .pf file to %s in current python version\n", new_question, old_names[n]);
      strcpy (old_question, old_names[n]);
      return (1);
    }

  }

  return (0);

}
