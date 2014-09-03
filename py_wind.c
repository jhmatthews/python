
/***********************************************************
                                       Space Telescope Science Institute

 Synopsis:
	py_wind is a program which can be used to display various parameters of a wind
		as calculated by python.  This is the main routine.

Arguments:		

	py_wind [-h] [-s] [-p parameter_file] [root]

	where
		-h 	prints out a short help file and exits (see help routine below)
		-s	causes certain parameters to be printed out in individual files
			after which the program exits
	 	-d	looks for and cylcles through the diagnostic files printing
			out all of the ascii files for each diagnostic file
		-p parameter file
			Instead of reading the choices from the command line read them
			from a parameter file
		root	optional root name of wind_save file.  If this is not given,
			the user is queried for this


Returns:
 
Description:	
	
	Py_wind simply reads and then displays the wind file created by python.  It can
	select various parameters from the wind and write them to a file so that one
	can create displays of them.  

	One can enter the commands interactively, which is the usual case, or one can 
	use a parameter file for the commands.  Note however, that one needs to be 
	sure that the parameter file responds to requests for choices int he correct order.
	This can be accomplished by first going through a file interactively.  The results
	of the interactive section will be stored in py_wind.pf (if you end with a q, and
	not an EOF response to the choice question)  EOF terminates the program at that point
	before the command file is written to py_wind.pf
	
	The file can contain either the original gridding which was used by python, in which 
	case the file prefix will be "x.", or "z", in which case it will be regridded to a 
	linear array.  This option is intended so one can create a contour plot more easily.


		
Notes:

History:
 	97jun	ksl	Coding on py_wind began.
 	97sep13	ksl	Modified slightly to reflect the fact that EOF exits program through rdpar
 	97oct2	ksl	Modified so that t_e and t_r can be printed separately
 	98feb	ksl	Modified to add the amount of energy and momentum flux absorbed per cell
 				as output possibilities. Also in py_ion added the ability to print out the
 				luminosity in CIV if the wind was assumed to be a thin plasma.
	98jul	ck	Modified so the number of ionzations and recombinations can be checked
			against one another
	98oct	ksl	Modified abs_sum and lum_suim so that the various processes which go
				into heating and cooling the plasma can be displayed individually.
	99nov	ksl	Added an additional option to create a fixed set of output files. 
	00mar	ksl	Updated main routine and write_array subroutine to allow one to choose
			whether to relinearize the grid or to leave it in as it was in the windfile
			Added output files for velocity.
	01sep	ksl	Updated ion_summary so that either concentrations or ion fractions could be
			examined.
	02apr	ksl	Fixed so could read wind_save files with various sizes 
	02jun	ksl	Fixed and so could read diagnostic outputs when one had printed out 
			results at each iteration of the ionization cycle
	04nov	ksl     Made major modification to eliminate duplicative code, by adding new
			task display.  Also introduced variable py_wind_project to switch
			between raw mode (regardless of coord system), and a mode which projects
			on to a yz plane the system of coordinates.
	04nov	ksl	Added a choice to explore what is the situation at a specific position
			in the grid.  This was added to see how different coordinate systems
			compared.
	05jan	ksl	54f -- Modified			
	05apr	ksl	56 -- Eliminated MDIM from the routines in this file.  
			MDIM is isolated to routines like display.  Added swithches to
			change choices regarding what, if anything is printed to a 
			file.
	05jul	ksl	56d -- Moved all of the subroutines to a separate routine py_wind_sub.c
			This is intended to facilitate the create of templates by Makefile.
	06jul	ksl	57g -- Made modifications to reflect the change in the way structures
			had been defined.
	090117	ksl	68a -- Added command line functionality so could run more easily
			in a non interactive mode
	090125	ksl	Changed functionality of I and i switches to allow for printing
			out information about the number of scatters, etc, by an ion
			in the wind.  This reflects a desire to be able to understand
			better where in the wind ions are "active" in creating the
			spectrum.
	111002	ksl	Added rho to one of the options
	111125	ksl	Modified the way the choice options are displayed, and put the various
			choices into a subroutine.  The goal here is to make py_wind easier
			to run from the command line.  The output files also been modifed
			so that we see which cells are in the wind and so the cell numbers are
			printed out.
	111125	ksl	Incorporated a crude capability to use a parameter file that takes
			a series of commands and prints out the files associated with them.
	111227	ksl	Added the capability to read a new wind file while the program
			is running

**************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "atomic.h"
#include "python.h"


//char *choice_options;

/* 111125 - ksl - Replaced print statements giving choices with a string. The point is to be able to include
 * the choices in the help string  Note carefullly the format if you revise this
 * lines should end with \n\  to make the string continue.  Do not leave any trailing spaces after the last
 * \ to avoid warnings
 */

char *choice_options = "\n\
   n=ne,  R=rho,  v=vel,        i=ion info, j=ave_tau, f=ave_freq, p=nphot, S=sim_alpha\n\
   r=t_r, t=t_e,  w=rad_weight,  s=vol,     l=lum,     C=cooling/heating,  b=adiabatic cooling\n\
   a=abs, c=c4,   g=photo,       h=recomb,  k=tau H,   l=lum,     m=F_rad, x=total, y=mod_te,\n\
   o=overview,    e=everything, P=Partial emission meas, I=Ionisation parameter\n\
   W=wind_region, D=dvds_ave, X=position summary, M=macro atom info, G=inner shell\n\
   d=convergence status  E=convergence_all_info   B=PlasmaPtr  J=Radiation density\n\
   H=All Heating and Cooling mechanisms in one shot  O=Spectral model parameters S=Spectral models\n\
   z=Zoom,u=unZoom,Z=switch to/from raw and yz projected modes, F=Create files, A=Change file write defaults\n\
   N=new.windfile q=quit (preferred over EOF)\n";

int
main (argc, argv)
     int argc;
     char *argv[];
{

//  WindPtr w;

  int i;
  int ochoice;
  char c;

  char root[LINELENGTH], input[LINELENGTH], wspecfile[LINELENGTH],
    specfile[LINELENGTH];
  char windradfile[LINELENGTH], windsavefile[LINELENGTH];
  char parameter_file[LINELENGTH];
  char photfile[LINELENGTH];
  double freq;
  int interactive;


  // py_wind uses rdpar, but only in an interactive mode. As a result 
  // there is no associated .pf file

  interactive = 1;		/* Default to the standard operating mofe for py_wind */
  strcpy (parameter_file, "NONE");

  if (argc == 1)
    {
      printf ("Root for wind file :");
      fgets (input, LINELENGTH, stdin);
      get_root (root, input);
    }
  else
    {
      for (i = 1; i < argc; i++)
	{
	  if (strcmp (argv[i], "-h") == 0)
	    {
	      py_wind_help ();
	    }
	  else if (strcmp (argv[i], "-d") == 0)
	    {
	      interactive = -1;
	    }

	  else if (strcmp (argv[i], "-s") == 0)
	    {
	      interactive = 0;
	    }
	  else if (strcmp (argv[i], "-p") == 0)
	    {
	      interactive = 0;
	      i = i + 1;
	      strcpy (parameter_file, argv[i]);
	    }
	  else if (strncmp (argv[i], "-", 1) == 0)
	    {
	      Error ("py_wind: unknown switch %s\n", argv[i]);
	      py_wind_help ();
	    }
	}

      strcpy (input, argv[argc - 1]);
      get_root (root, input);
    }


  printf ("Reading data from file %s\n", root);

  /* Now create the names of all the files which will be written */

  strcpy (wspecfile, root);
  strcpy (specfile, root);
  strcpy (windradfile, root);
  strcpy (windsavefile, root);
  strcpy (photfile, root);

  strcat (wspecfile, ".spec_tot");
  strcat (specfile, ".spec");
  strcat (windradfile, ".wind_rad");
  strcat (windsavefile, ".wind_save");
  strcat (photfile, ".phot");
//  strcat (parfile, ".pf");


  /* Initialize other variables here */

  py_wind_project = 1;		// The default is to try to project onto a yz plane 

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
  aaa = calloc (sizeof (freq), NDIM2);

  printf ("Read wind_file %s\n", windsavefile);

  get_atomic_data (geo.atomic_filename);

  printf ("Read Atomic data from %s\n", geo.atomic_filename);


/* Produce a standard set of output files and exit*/
  if (interactive == 0 && strcmp (parameter_file, "NONE") == 0)
    {
      zoom (1);			/* This affects the logfile */
      ochoice = 1;
      complete_file_summary (wmain, root, ochoice);
      exit (0);
    }
  else if (interactive == -1)
    {
      /* Write the sumary ascii files for all of the diagnostic files 
       * as well as the original */
      zoom (1);			/* This affects the logfile */
      ochoice = 1;
      complete_file_summary (wmain, root, ochoice);
      i = 0;
      strcpy (root, "");
      sprintf (root, "python%02d", i);
      strcpy (windsavefile, "");
      sprintf (windsavefile, "python%02d.wind_save", i);
      while (wind_read (windsavefile) > 0)
	{
	  Log ("Trying %s %s\n", windsavefile, root);
	  complete_file_summary (wmain, root, ochoice);
	  strcpy (root, "");
	  sprintf (root, "python%02d", i);
	  strcpy (windsavefile, "");
	  sprintf (windsavefile, "python%02d.wind_save", i);
	  i++;
	}
      exit (0);
    }




  if (strcmp (parameter_file, "NONE") != 0)
    {
      zoom (1);			/* This affects the logfile */
      ochoice = 1;
      opar (parameter_file);
    }


/* Choices */
  ochoice = 0;
  rdint ("Make_files(0=no,1=original,2=regrid_to_linear)", &ochoice);
  c = 'i';
  zoom (1);



/* 111126 - Consolodated choice statements into a string that is an external variable 
 * so that it can be printed out whenever necessary.
 */

  printf ("%s\n", choice_options);
  printf ("Model %s   :\n", root);
  rdchar ("Choice", &c);

  while (c != EOF)
    {
      one_choice (c, root, ochoice);
      printf ("%s\n", choice_options);
      rdchar ("Choice", &c);
    }

  return (0);
}


/***********************************************************
                        Space Telescope Science Institute

Synopsis:

Arguments:		



Returns:
 
Description:	


		
Notes:

History:
	111125	ksl	In an attempt to make it easier to operate py_wind
			from the command line the huge choice loop has been
			moved into a separate routine

**************************************************************/

int
one_choice (choice, root, ochoice)
     char choice;
     char *root;
     int ochoice;
{
  double lambda, freq;
  int n, istate, iswitch;
  char windsavefile[LINELENGTH];

  iswitch = 0;

  /* JM 1312 --initialise variables to avoid compilation warnings */
  istate = 0;
  n = 0;


  switch (choice)
    {
    case 'a':			/* Energy absorbed */
      abs_summary (wmain, root, ochoice);
      break;
    case 'A':			// Change the file defaults
      rdint ("Make_files(0=no,1=original,2=regrid_to_linear)", &ochoice);
      break;
    case 'b':			/*Adiabatic cooling */
      adiabatic_cooling_summary (wmain, root, ochoice);
      break;
    case 'B':
      plasma_cell (wmain, root, ochoice);
      break;
    case 'c':			/*C4 emission */
      line_summary (wmain, n, istate, root, ochoice);
      break;
    case 'C':			/*the ratio cooling to heating */
      coolheat_summary (wmain, root, ochoice);
      break;
    case 'd':
      convergence_summary (wmain, root, ochoice);
      break;
    case 'D':			/* dvds summary */
      dvds_summary (wmain, root, ochoice);
      break;
    case 'E':
      convergence_all (wmain, root, ochoice);
      break;
    case 'e':			/* print out everything about an element */
      wind_element (wmain);
      break;
    case 'f':			/* Electron summary */
      freq_summary (wmain, root, ochoice);
      break;
    case 'F':			/* Complete file summary */
      complete_file_summary (wmain, root, ochoice);
      break;
    case 'g':			/*n photo */
      photo_summary (wmain, root, ochoice);
      break;
    case 'G':			/* inner shell summary */
      inner_shell_summary (wmain, root, ochoice);
      break;
    case 'h':			/*n photo */
      Log ("Don't get discouraged.  This takes a little while!");
      recomb_summary (wmain, root, ochoice);
      break;
    case 'H':			/* heating and cooling mechanisms breakdown */
      heatcool_summary (wmain, root, ochoice);
      break;
    case 'i':			/* Allow user to display information about ions in the wind */

      rdint ("Ion_info_type(0=fraction,1=density,2=scatters,3=abs", &iswitch);

      n = 6;
      istate = 4;

      while (rdint ("element(0=return)", &n) != EOF)
	{
	  if (n <= 0)
	    break;
	  rdint ("ion", &istate);
	  ion_summary (wmain, n, istate, iswitch, root, ochoice);	// 0 implies ion fractions
	}
      break;
    case 'I':
      IP_summary (wmain, root, ochoice);
      break;

    case 'j':			/* Calculate the average tau at the center of a cell */
      n = 6;
      istate = 4;
      lambda = 1550;

      rddoub ("wavelength", &lambda);
      freq = C / (lambda * 1.e-8);

      while (rdint ("element(0=return)", &n) != EOF)
	{
	  if (n <= 0)
	    break;
	  rdint ("ion", &istate);
	  tau_ave_summary (wmain, n, istate, freq, root, ochoice);
	}
      break;
    case 'J':			/* radiation density in cell */
      J_summary (wmain, root, ochoice);
      break;
    case 'k':			/* tau at H edge */
      tau_h_summary (wmain, root, ochoice);
      break;
    case 'K':			/* cell J split by direct photons and scattered photons */
      J_scat_summary (wmain, root, ochoice);
      break;
    case 'l':			/* Lum of shell */
      lum_summary (wmain, root, ochoice);
      break;
    case 'm':			/* Radiation force */
      mo_summary (wmain, root, ochoice);
      break;
    case 'M':
      macro_summary (wmain, root, ochoice);
      break;
    case 'n':			/* Electron summary */
      electron_summary (wmain, root, ochoice);
      break;
    case 'N':			/* Read a different wind save file */
      rdstr ("New.rootname", root);
      strcpy (windsavefile, root);
      strcat (windsavefile, ".wind_save");
      if (wind_read (windsavefile) < 0)
	{
	  Error ("one_choice: Could not read %s", windsavefile);
	}


/* aaa is used to store variable for writing to files for the purpose of plotting*/
      if (aaa != NULL)
	{
	  free (aaa);
	}
      aaa = calloc (sizeof (freq), NDIM2);

      printf ("Read wind_file %s\n", windsavefile);

      get_atomic_data (geo.atomic_filename);

      printf ("Read Atomic data from %s\n", geo.atomic_filename);

      break;
    case 'o':			/* overview */
      overview (wmain, root);
      break;
    case 'O':			/* spectral model parameters */
      model_bands (wmain, root, ochoice);
      break;
    case 'p':			/* nphot summary */
      nphot_summary (wmain, root, ochoice);
      break;
    case 'P':			/* Allow user to display information about the wind */

      n = 6;
      istate = 4;

      while (rdint ("element(0=return)", &n) != EOF)
	{
	  if (n <= 0)
	    break;
	  rdint ("ion", &istate);
	  partial_measure_summary (wmain, n, istate, root, ochoice);
	}
      break;
    case 'r':			/* Temp summary */
      temp_rad (wmain, root, ochoice);
      break;
    case 'R':			/* Rho summary */
      rho_summary (wmain, root, ochoice);
      break;
    case 's':			/* Volume summary */
      vol_summary (wmain, root, ochoice);
      break;
    case 'S':
      alpha_summary (wmain, root, ochoice);
      break;
    case 't':			/* Temp summary */
      temp_summary (wmain, root, ochoice);
      break;
    case 'T':
      thompson (wmain, root, ochoice);
      break;
    case 'v':			/* Velocity summary */
      velocity_summary (wmain, root, ochoice);
      break;
    case 'V':			/* Split of scatters in the cell between electron and resonant */
      nscat_split (wmain, root, ochoice);
      break;
    case 'w':			/* inten weight summary */
      weight_summary (wmain, root, ochoice);
      break;
    case 'W':			/*Show regions in the wind */
      wind_reg_summary (wmain, root, ochoice);
      break;
    case 'x':			/*Total emission */
      total_emission_summary (wmain, root, ochoice);
      break;
    case 'X':			/* Position summary */
      position_summary (wmain);
      break;
    case 'y':			/* Recalculate temperatures */
      modify_te (wmain, root, ochoice);
      break;
    case 'Y':			/* Split of photons from different sources */
      phot_split (wmain, root, ochoice);
      break;
    case 'z':			/* inspect a specific region */
      zoom (0);
      break;
    case 'Z':			/* Switch between raw and projected mode */
      if (py_wind_project == 0)
	{
	  py_wind_project = 1;
	  Log ("Switching to raw display, whatever the coordinate system");
	}
      else if (py_wind_project == 1)
	{
	  py_wind_project = 0;
	  Log ("Switching to projected y z display");
	}
      break;
    case 'u':			/* Go back to full image */
      zoom (1);
      break;
    case 'q':			/* quit */
      /* Write out a parameterfile that gives all of the commands used in this run */
      cpar ("py_wind.pf");
      exit (0);
      break;

    }

  return (0);
}

//  goto a;

//}



/***********************************************************
                        Space Telescope Science Institute

Synopsis:
	py_wind_help simply prints out help to the screen
	and exits

Arguments:		



Returns:
 
Description:	


		
Notes:
	Unfortunately unlike python the program language
	c is not self documenting

History:
	111125	ksl	Modified so routine also prints out
			the string that contains all the choices

**************************************************************/


int
py_wind_help ()
{

  char *some_help;

/* Beginning of a string to provide some help for py_wind */
  some_help = "\
\n\
This program reads a wind save file created by python and examine the wind structure\\
\n\
	Usage: py_wind [-h] [-s] [-p parameter_file] [root] \n\
\n\
	where\n\
		-h 	prints out a short help file and exits (see help routine below) \n\
		-s	causes certain parameters to be printed out in individual files \n\
			after which the program exits \n\
		-d	searches for the diagnostic wind save files and prints information \n\
			to ascii files \n\
		-p	Gets information a parameter file \n\
		root	optional root name of wind_save file.  If this is not given, \n\
			the user is queried for this \n\
\n\
";

  printf ("%s\n", some_help);

  printf ("Choices are indicated below\n%s\n", choice_options);

  exit (0);
}
