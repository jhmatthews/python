#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "atomic.h"
#include "python.h"



/***********************************************************
             University of Southampton

Synopsis: 
  init_files 
   
Arguments:   
  argc            command line arguments 

Returns:
  restart_start   1 if restarting
 
 
Description:  

Notes:

History:
  1502  JM  Moved here from main()

**************************************************************/

int
init_log_and_windsave (restart_stat)
     int restart_stat;
{
  FILE *fopen (), *qptr;

  if (restart_stat == 0)
    {				// Then we are simply running from a new model
      xsignal_rm (files.root);	// Any old signal file
      xsignal (files.root, "%-20s %s \n", "START", files.root);
      Log_init (files.diag);
    }
  else
    {
      /* Note that alghough we chekc that we dan open the windsave file, it is not read here.   */

      strcpy (files.windsave, files.root);
      strcat (files.windsave, ".wind_save");
      qptr = fopen (files.windsave, "r");

      if (qptr != NULL)
	{
	  /* Then the file does exist and we can restart */
	  fclose (qptr);
	  xsignal (files.root, "%-20s %s\n", "RESTART", files.root);
	  Log_append (files.diag);
	}
      else
	{
	  /* It does not exist and so we start from scratch */
	  restart_stat = 0;
	  xsignal_rm (files.root);	// Any old signal file
	  xsignal (files.root, "%-20s %s \n", "START", files.root);
	  Log_init (files.diag);
	}
    }

  return (0);
}


/***********************************************************
             University of Southampton

Synopsis: 

  setup_dfudge works out dfudge and returns it to the user.
  the global variable DFUDGE is not altered here.
   
Arguments:		

Returns:
  dfudge 	the push through distance 

Description:	

  DFUDGE is the push through distance when photons are not
  travelling within wind cells.  (Inside a cell in a domain
  the push through distance is defined on a cell by cell
  basis)
  
Notes:

  There are two competing factors in defining DFUDGE.  It
  should be short enough so that the push through distane
  goes only a small way into s cell.  It should be large
  enough though that round-off errors do not prevent one
  from actually getting into a cell.  When this happens
  one can get "stuck photons".

  Prior to domains, we had effectively 3 separate ways of
  defining dfudge, one for the SHELL model, one for normal
  stellar systems and one for AGN,  The shell model was
  different because it can be thin but large scale.

  The current version of setup_dfudge preserves the onld
  value of dfudge as much as possible (except the it
  senses the old SHELL case bu the differene between
  gero.rmax and geo.rmin


History:
	1502	JM 	Moved here from main()
	1605	ksl	Revised to remove the dependence on
			a specific geometry namely SHELL

**************************************************************/

double
setup_dfudge ()
{
  double dfudge;
  double delta;

  delta = geo.rmax - geo.rmin;

  if (delta < 1.e8)
    {
      dfudge = (geo.rmax - geo.rmin) / 1000.0;
    }
  else if (delta < 1e15)
    {
      dfudge = 1e5;
    }
  else
    {
      dfudge = geo.rmax / 1.e10;
    }

  Log ("DFUDGE set to %e based on geo.rmax\n", dfudge);

  return (dfudge);
}



/***********************************************************
             University of Southampton

Synopsis: 
  setup_windcone sets up the windcone 
   
Arguments:		

Returns:

Description:	

Notes:
  The angles thetamin and
  thetamax are all defined from the z axis, so that an angle of 0
  is a flow that is perpeindicular to to the disk and one that is
  close to 90 degrees will be parallel to the plane of the disk
  geo.wind_thetamin and max are defined in the routines that initialize
  the various wind models, e. g. get_sv_wind_parameters. These
  have been called at this point.  

  z is the place where the windcone intercepts the z axis
  dzdr is the slope 

  111124 fixed notes on this - ksl
History:
	1502  JM 	Moved here from main()
	1508	ksl	Modified to construct wind cones  fo
			all domains

**************************************************************/

int
setup_windcone ()
{
  int ndom;

  for (ndom = 0; ndom < geo.ndomain; ndom++)
    {

      if (zdom[ndom].wind_thetamin > 0.0)
	{
	  zdom[ndom].windcone[0].dzdr = 1. / tan (zdom[ndom].wind_thetamin);
	  zdom[ndom].windcone[0].z =
	    (-zdom[ndom].wind_rho_min / tan (zdom[ndom].wind_thetamin));
	}
      else
	{
	  zdom[ndom].windcone[0].dzdr = VERY_BIG;
	  zdom[ndom].windcone[0].z = -VERY_BIG;;
	}


      if (zdom[ndom].wind_thetamax > 0.0)
	{
	  zdom[ndom].windcone[1].dzdr = 1. / tan (zdom[ndom].wind_thetamax);
	  zdom[ndom].windcone[1].z =
	    (-zdom[ndom].wind_rho_max / tan (zdom[ndom].wind_thetamax));
	}
      else
	{
	  zdom[ndom].windcone[1].dzdr = VERY_BIG;
	  zdom[ndom].windcone[1].z = -VERY_BIG;;
	}
    }
  return (0);
}












/***********************************************************
             University of Southampton

Synopsis: 
  setup_created_files 
   
Arguments:    

Returns:

Description:  

Notes:

History:
  1502  JM  Moved here from main()

**************************************************************/

int
setup_created_files ()
{
  int opar_stat;

  opar_stat = 0;		/* 59a - ksl - 08aug - Initialize opar_stat to indicate that if we do not open a rdpar file, 
				   the assumption is that we are reading from the command line */


  if (strncmp (files.root, "stdin", 5) == 0
      || strncmp (files.root, "rdpar", 5) == 0 || files.root[0] == ' '
      || strlen (files.root) == 0)
    {
      strcpy (files.root, "mod");
      Log
	("Proceeding in interactive mode\n Output files will have rootname mod\n");
    }

  else
    {
      strcpy (files.input, files.root);
      strcat (files.input, ".pf");

      if ((opar_stat = opar (files.input)) == 2)
	{
	  Log ("Reading data from file %s\n", files.input);
	}
      else
	{
	  Log ("Creating a new parameter file %s\n", files.input);
	}

    }


  /* Now create the names of all the files which will be written.  Note that some files
     have the same root as the input file, while others have a generic name of python.
     This is intended so that files which you really want to keep have unique names, while
     those which are for short-term diagnostics are overwritten.  ksl 97aug. */

  strcpy (basename, files.root);	//56d -- ksl --Added so filenames could be created by routines as necessary

  strcpy (files.wspec, files.root);	//generated photons
  strcpy (files.lwspec, files.root);	//generated photon in log space

  strcpy (files.wspec_wind, files.root);
  strcpy (files.lwspec_wind, files.root);

  strcpy (files.spec, files.root);
  strcpy (files.lspec, files.root);

  strcpy (files.spec_wind, files.root);
  strcpy (files.lspec_wind, files.root);

  strcpy (files.new_pf, files.root);
  strcat (files.new_pf, ".out.pf");


  strcpy (files.windrad, "python");
  strcpy (files.windsave, files.root);
  strcpy (files.specsave, files.root);

  /* 130722 JM we now save python.phot and disk.diag files under diag_root folder */
  strcpy (files.phot, files.diagfolder);
  strcpy (files.disk, files.diagfolder);
  strcat (files.phot, "python");
  strcat (files.disk, files.root);

  strcat (files.wspec, ".spec_tot");
  strcat (files.lwspec, ".log_spec_tot");

  strcat (files.wspec_wind, ".spec_tot_wind");
  strcat (files.lwspec_wind, ".log_spec_tot_wind");


  strcat (files.spec, ".spec");
  strcat (files.lspec, ".log_spec");

  strcat (files.spec_wind, ".spec_wind");
  strcat (files.lspec_wind, ".log_spec_wind");


  strcat (files.windrad, ".wind_rad");
  strcat (files.windsave, ".wind_save");
  strcat (files.specsave, ".spec_save");
  strcat (files.phot, ".phot");
  strcat (files.disk, ".disk.diag");


  return (opar_stat);
}





