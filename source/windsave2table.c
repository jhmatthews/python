
/***********************************************************
                                       Space Telescope Science Institute

 Synopsis:
	windsave2table writes key variables in a wind save file to an astropy table    
		as calculated by python.  This is the main routine.

Arguments:		

	py_wind  windsave_root



Returns:
 
Description:	
	

	
Notes:



History:
	150428	ksl	Adpated from routines in py_wind.c

**************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "atomic.h"
#include "python.h"



int
main (argc, argv)
     int argc;
     char *argv[];
{


  int ochoice;

  char root[LINELENGTH], input[LINELENGTH];
  char outputfile[LINELENGTH];
  char windsavefile[LINELENGTH];
  char parameter_file[LINELENGTH];
  double xdoub;
  int create_master_table(), create_ion_table();


  // py_wind uses rdpar, but only in an interactive mode. As a result 
  // there is no associated .pf file

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

/* Note that wind_read allocates the space for the WindPtr array.  The
reason that this is done in wind_read is because wind_read also reads
the geo structure.  The geo struc contains the dimensions of the wind
array and so until it is read the space for the WindPtr structure cannot
be allocated.  Howver, one cannot assign w in a subroutine.  Therefore
the old call, namely wind_read(w,windsavefile, will not work.  The statement
w=wmain is a way around this since wmain is external and therefore can
be assigned.  Then w can be set to this ptr and all is restored. The
use of w is endemic in the program. and it is always called through main.
I did not change this now.  Though it could be done.  02apr ksl */

  if (wind_read (windsavefile) < 0)
    {
      Error ("py_wind: Could not open %s", windsavefile);
      exit (0);
    }

/* aaa is used to store variable for writing to files for the purpose of plotting*/
  aaa = calloc (sizeof (xdoub), NDIM2);

  printf ("Read wind_file %s\n", windsavefile);

  get_atomic_data (geo.atomic_filename);

  printf ("Read Atomic data from %s\n", geo.atomic_filename);


  ochoice = 1;


  create_master_table (root);
  create_ion_table (root, 6);
  create_ion_table (root, 7);
  create_ion_table (root, 8);
  create_ion_table (root, 14);
  create_ion_table (root, 26);
  return (0);
}




/***********************************************************
                                       Space Telescope Science Institute

Synopsis:

	create_master_table writes a selected variables of in the windsaave
	file to an astropy table

	It is intended to be easily modifible.

Arguments:		

	rootname of the file that will be written out


Returns:
 
Description:	

	The routine reads data directly from wmain, and then calls 
	get_one or get_ion multiple times to read info from the Plasma
	structure.  
	
	It then writes the data to an  astropy file
Notes:

	To add a variable one just needs to define the column_name
	and send the appropriate call to either get_one or get_ion.

	There is some duplicated code in the routine that pertains
	to whether one is dealing with a spherecial or a 2d coordinate
	system.  It should be possible to delete this



History:
	150428	ksl	Adpated from routines in py_wind.c
	150501	ksl	Cleaned this routine up, and added a few more variables

**************************************************************/

int
create_master_table (rootname)
     char rootname[];
{
  char filename[132];
  double *get_one ();
  double *get_ion ();
  double *c[50],*converge;
  char column_name[50][20];
  char one_line[1024], start[132], one_value[20];


  int i, ii, jj;
  int n, ncols;
  FILE *fopen (), *fptr;

  strcpy (filename, rootname);
  strcat (filename, "_master.txt");

  fptr = fopen (filename, "w");

  /* Get the variables that one needs */

  c[0] = get_one ("ne");
  strcpy (column_name[0], "ne");

  c[1] = get_one ("t_e");
  strcpy (column_name[1], "t_e");

  c[2] = get_one ("t_r");
  strcpy (column_name[2], "t_r");

  c[3] = get_ion (1, 1, 0);
  strcpy (column_name[3], "h1");

  c[4] = get_ion (2, 2, 0);
  strcpy (column_name[4], "he2");

  c[5] = get_ion (6, 4, 0);
  strcpy (column_name[5], "c4");

  c[6] = get_ion (7, 5, 0);
  strcpy (column_name[6], "n5");

  c[7] = get_ion (8, 6, 0);
  strcpy (column_name[7], "o6");

  c[8]=get_one("dmo_dt_x");
  strcpy (column_name[8], "dmo_dt_x");


  c[9]=get_one("dmo_dt_y");
  strcpy (column_name[9], "dmo_dt_y");

  c[10]=get_one("dmo_dt_z");
  strcpy (column_name[10], "dmo_dt_z");

  c[11] = get_one ("ntot");
  strcpy (column_name[11], "ntot");

  ncols = 12;


  converge=get_one("converge");

  /* At this point oll of the data has been collected */





  if (geo.coord_type == SPHERICAL)
    {


      /* First assemble the heder line
       */

      sprintf (start, "%8s %4s %8s %6s %8s %8s %8s ", "r", "i", "inwind", "converge","v_x",
	       "v_y", "v_z");
      strcpy (one_line, start);
      n = 0;
      while (n < ncols)
	{
	  sprintf (one_value, "%9s ", column_name[n]);
	  strcat (one_line, one_value);

	  n++;
	}
      fprintf (fptr, "%s\n", one_line);





      /* Now assemble the lines of the table */
      for (i = 0; i < NDIM2; i++)
	{
	  // This line is different from the two d case
	  sprintf (start, "%8.2e %4d %6d %8.0f %8.2e %8.2e %8.2e ", wmain[i].r, i,
		   wmain[i].inwind,converge[i],wmain[i].v[0], wmain[i].v[1],
		   wmain[i].v[2]);
	  strcpy (one_line, start);
	  n = 0;
	  while (n < ncols)
	    {
	      sprintf (one_value, "%9.2e ", c[n][i]);
	      strcat (one_line, one_value);
	      n++;
	    }
	  fprintf (fptr, "%s\n", one_line);
	}
    }
  else
    {

      Log ("Processing 2d wind\n");

      /* First assemble the header line
       */

      sprintf (start, "%8s %8s %4s %4s %6s %8s %8s %8s %8s ", "x", "z", "i", "j",
	       "inwind", "converge","v_x", "v_y", "v_z");
      strcpy (one_line, start);
      n = 0;
      while (n < ncols)
	{
	  sprintf (one_value, "%9s ", column_name[n]);
	  strcat (one_line, one_value);

	  n++;
	}
      printf ("%s\n", one_line);
      fprintf (fptr, "%s\n", one_line);


      /* Now assemble the lines of the table */
      for (i = 0; i < NDIM2; i++)
	{
	  wind_n_to_ij (i, &ii, &jj);
	  sprintf (start, "%8.2e %8.2e %4d %4d %6d %8.0f %8.2e %8.2e %8.2e ",
		   wmain[i].xcen[0], wmain[i].xcen[2], ii, jj,
		   wmain[i].inwind, converge[i],wmain[i].v[0], wmain[i].v[1],
		   wmain[i].v[2]);
	  strcpy (one_line, start);
	  n = 0;
	  while (n < ncols)
	    {
	      sprintf (one_value, "%9.2e ", c[n][i]);
	      strcat (one_line, one_value);
	      n++;
	    }
	  fprintf (fptr, "%s\n", one_line);
	}
    }
  return (0);
}




/***********************************************************
                                       Space Telescope Science Institute

 Synopsis:
	windsave2table writes key variables in a wind save file to an astropy table    
		as calculated by python.  This is the main routine.

Arguments:		

	py_wind  windsave_root



Returns:
 
Description:	
	

	
Notes:



History:
	150428	ksl	Adpated from routines in py_wind.c

**************************************************************/
int
create_ion_table (rootname, iz)
     char rootname[];
     int iz;			// Where z is the element 
{
  char filename[132];
  double *get_one ();
  double *get_ion ();
  double *c[50];
  int first_ion, number_ions;
  char element_name[20];
  int istate[50];
  char one_line[1024], start[132], one_value[20];


  int i, ii, jj, n;
  FILE *fopen (), *fptr;
// First we actually need to determine what ions exits, but we will ignore this for now

  i = 0;
  while (i < nelements)
    {
      if (ele[i].z == iz)
	{
	  break;
	}
      i++;
    }


  first_ion = ele[i].firstion;
  number_ions = ele[i].nions;
  strcpy (element_name, ele[i].name);

  Log ("element %d %d %s\n", first_ion, number_ions, element_name);

/* Open the output file */

  sprintf (filename, "%s_%s.txt", rootname, element_name);

  fptr = fopen (filename, "w");


  i = 0;
  while (i < number_ions)
    {
      istate[i] = ion[first_ion + i].istate;

      c[i] = get_ion (iz, istate[i], 0);
      i++;
    }



  if (geo.coord_type == SPHERICAL)
    {


      /* First assemble the heder line
       */

      sprintf (start, "%8s %4s %6s ", "r", "i", "inwind");
      strcpy (one_line, start);
      n = 0;
      while (n < number_ions)
	{
	  sprintf (one_value, "     i%02d ", istate[n]);
	  strcat (one_line, one_value);

	  n++;
	}
      fprintf (fptr, "%s\n", one_line);





      /* Now assemble the lines of the table */
      for (i = 0; i < NDIM2; i++)
	{
	  // This line is different from the two d case
	  sprintf (start, "%8.2e %4d %6d ", wmain[i].r, i, wmain[i].inwind);
	  strcpy (one_line, start);
	  n = 0;
	  while (n < number_ions)
	    {
	      sprintf (one_value, "%8.2e ", c[n][i]);
	      strcat (one_line, one_value);
	      n++;
	    }
	  fprintf (fptr, "%s\n", one_line);
	}
    }
  else
    {

      Log ("Processing 2d wind\n");

      /* First assemble the heder line
       */

      sprintf (start, "%8s %8s %4s %4s %6s ", "x", "z", "i", "j", "inwind");
      strcpy (one_line, start);
      n = 0;
      while (n < number_ions)
	{
	  sprintf (one_value, "     i%02d ", istate[n]);
	  strcat (one_line, one_value);

	  n++;
	}
      fprintf (fptr, "%s\n", one_line);

      /* Now assemble the lines of the table */
      for (i = 0; i < NDIM2; i++)
	{
	  wind_n_to_ij (i, &ii, &jj);
	  sprintf (start, "%8.2e %8.2e %4d %4d %6d ", wmain[i].xcen[0],
		   wmain[i].xcen[2], ii, jj, wmain[i].inwind);
	  strcpy (one_line, start);
	  n = 0;
	  while (n < number_ions)
	    {
	      sprintf (one_value, "%8.2e ", c[n][i]);
	      strcat (one_line, one_value);
	      n++;
	    }
	  fprintf (fptr, "%s\n", one_line);
	}
    }
  return (0);
}


/*


*/

/***********************************************************
                                       Space Telescope Science Institute

Synopsis:

	Get get denensity, etc for one particular ion

Arguments:		




Returns:
 
Description:	
	

   This will return an array with all zeros if there is no such ion
	
Notes:



History:
	150428	ksl	Adpated from routines in py_wind.c

**************************************************************/

double *
get_ion (element, istate, iswitch)
     int element, istate, iswitch;
{
  int nion, nelem;
  int n;
  // char choice[LINELENGTH], iname[LINELENGTH];
  char name[LINELENGTH];
  // char filename[LINELENGTH];
  int nplasma;
  double *x;



  x = (double *) calloc (sizeof (double), NDIM2);
/* Find the ion */
  nion = 0;
  while (nion < nions
	 && !(ion[nion].z == element && ion[nion].istate == istate))
    nion++;
  if (nion == nions)
    {
      Log ("Error--element %d ion %d not found in define_wind\n", element,
	   istate);
      return (x);
      // return (-1);
    }
  nelem = 0;
  while (nelem < nelements && ele[nelem].z != element)
    nelem++;

  strcpy (name, "");


  // Now populate the array
  for (n = 0; n < NDIM2; n++)
    {
      x[n] = 0;
      nplasma = wmain[n].nplasma;
      if (wmain[n].vol > 0.0 && plasmamain[nplasma].ne > 1.0)
	{
	  if (iswitch == 0)
	    {
	      sprintf (name, "Element %d (%s) ion %d fractions\n", element,
		       ele[nelem].name, istate);
	      x[n] = plasmamain[nplasma].density[nion];
	      x[n] /=
		((plasmamain[nplasma].density[0] +
		  plasmamain[nplasma].density[1]) * ele[nelem].abun);
	    }
	  else if (iswitch == 1)
	    {
	      sprintf (name, "Element %d (%s) ion %d density\n", element,
		       ele[nelem].name, istate);
	      x[n] = plasmamain[nplasma].density[nion];
	    }
	  else if (iswitch == 2)
	    {
	      sprintf (name, "Element %d (%s) ion %d  #scatters\n", element,
		       ele[nelem].name, istate);
	      x[n] = plasmamain[nplasma].scatters[nion];
	    }
	  else if (iswitch == 3)
	    {
	      sprintf (name, "Element %d (%s) ion %d scattered flux\n",
		       element, ele[nelem].name, istate);
	      x[n] = plasmamain[nplasma].xscatters[nion];
	    }
	  else
	    {
	      Error ("xion_summary : Unknown switch %d \n", iswitch);
	      exit (0);
	    }
	}
    }



  return (x);
}

/**************************************************************************


  Synopsis:  
	Get a simple variable from the PlasmaPtr array

  Description:	

  Arguments:  

  Returns:

  	The values in the plasma pointer for this variable. A double
	will be returned even if the PlasmaPtr varible is an integer

  Notes:
  	Getting any simple variable from the plama structure should
	follow this template.

  History:
  	150429 ksl Adapted from te_summary in py_wind

 ************************************************************************/
double *
get_one (variable_name)
     char variable_name[];
{
  int n;
  int nplasma;
  double *x;

  x = (double *) calloc (sizeof (double), NDIM2);
  for (n = 0; n < NDIM2; n++)
    {
      x[n] = 0;
      if (wmain[n].vol > 0.0)
	{
	  nplasma = wmain[n].nplasma;

	  if (strcmp (variable_name, "ne") == 0)
	    {
	      x[n] = plasmamain[nplasma].ne;
	    }
	  else if (strcmp (variable_name, "t_e") == 0)
	    {
	      x[n] = plasmamain[nplasma].t_e;
	    }
	  else if (strcmp (variable_name, "t_r") == 0)
	    {
	      x[n] = plasmamain[nplasma].t_r;
	    }
	  else if (strcmp (variable_name, "converge")==0)
	  {
		  x[n]=plasmamain[nplasma].converge_whole;
	  }
	  else if (strcmp (variable_name, "dmo_dt_x")==0)
	  {
		  x[n]=plasmamain[nplasma].dmo_dt[0];
	  }
	  else if (strcmp (variable_name, "dmo_dt_y")==0)
	  {
		  x[n]=plasmamain[nplasma].dmo_dt[1];
	  }
	  else if (strcmp (variable_name, "dmo_dt_z")==0)
	  {
		  x[n]=plasmamain[nplasma].dmo_dt[2];
	  }
	  else if (strcmp (variable_name, "ntot")==0)
	  {
		  x[n]=plasmamain[nplasma].ntot;
	  }





			  
	  else
	    {
	      Error ("get_one: Unknown variable %s\n", variable_name);
	    }
	}
    }

  return (x);

}
