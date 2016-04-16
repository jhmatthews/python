
/**************************************************************************
                    Space Telescope Science Institute


  Synopsis:

  Description:	These are subroutines used by py_wind to write ascii
	files of various variables for use with smongo or normal display
	software.

  Arguments:		

  Returns:

  Notes:

  History:
	01aug	ksl	Replaced ndim with xdim so that ndim could be used
			in python.h
	05apr	ksl	56 -- eliminated fitswriting capability as part of
			a general reorganization of py_wind to accommodate
			mulitple files and for use with tecplot.  56d is
			the last version with fits capability
	

 ************************************************************************/

#include        <math.h>
#include        <stdio.h>
#include        <strings.h>
#include	<string.h>
#include	<stdlib.h>
#include	"atomic.h"
#include	"python.h"
#include        "templates.h"




/**************************************************************************
                    Space Telescope Science Institute


  Synopsis:  write_array writes out the arrays in a way which can be read easily
   into a routine which can do contouring or to make a fits file 



  Description:

  Arguments:		
  	choice		0-> don't write array 
			1--> write without interpolation
			2-> interpolate to linear array

  Returns:

  Notes:

  History:
	02feb	ksl	Fixed to allow different values of NDIM and MDIM.  This
			is tricky, but I believe that ain is defined such that
			the second axis is in the z direction so that the proper
			dimension is ain {NDIM] [MDIM]
	05apr	ksl	56 -- Completely rewritten to enable a variety of 
			coordinate systems and to produce tecplot type 
			output files.
	1111	ksl	71 -- Modified format of files to ease ability
			to plot the data in a variety of formats
	1111	ksl	71 - Tecplot dependencies removed

	

 ************************************************************************/
#define ODIM 256
int init_write_array = 1;
float aout[ODIM][ODIM];

int
write_array (filename, choice)
     char filename[];
     int choice;
{
  //Dynamical allocation is allowed, although I generally avoid it -- 05apr ksl
  float rin[NDIM], zin[MDIM];
  float r, z;
  float rmin, rmax, zmin, zmax;
  int ii, jj;
  FILE *fopen (), *fptr;
  char outfile[LINELENGTH];
  double length ();
  double xx[3];
  int i;
  int nn, nnn[4], nelem;
  double frac[4];

  if (choice == 0)		// Then we don't really want output files
    {
      return (0);
    }

// Open the appropriate file

  strcpy (outfile, filename);
  strcat (outfile, ".dat");	// Add standard extension to filenames, i.e. the one used by tecplot
  fptr = fopen (outfile, "w");


  /* Write out the header information for the file */
  fprintf (fptr, "# TITLE= \"%s\"\n", outfile);
  if (geo.coord_type == SPHERICAL)
    {
      fprintf (fptr, "# Coord_Sys SPHERICAL\n");
    }
  else if (geo.coord_type == CYLIND)
    {
      fprintf (fptr, "# Coord_Sys CYLIND\n");
    }
  else if (geo.coord_type == RTHETA)
    {
      fprintf (fptr, "# Coord_Sys RTHETA\n");
    }
  else if (geo.coord_type == CYLVAR)
    {
      fprintf (fptr, "# Coord_Sys CYLVAR\n");
    }
  else
    {
      Error ("write_array: Unknown coordinaate system type: %d\n",
	     geo.coord_type);
    }


// Put the r and z coord. grid  into easier to understand arrays
  for (i = 0; i < NDIM; i++)
    {
      zin[i] = wmain[i].x[2];
      rin[i] = wmain[i * MDIM].x[0];
    }


/* In py_wind the filenames are set to start with z if the output coordinates
are linear, and x otherwise.  This is not particularly transparent ?? ksl */

  if (choice == 1)
    {				// print out the original array elements

      if (geo.coord_type == SPHERICAL)
	{
    
    /* JM 1411 -- added header for reading with astropy ascii module */
    fprintf (fptr, "r var inwind i\n");

	  for (i = 0; i < NDIM; i++)
	    {
	      fprintf (fptr, "%8.2e %8.2e %3d %3d \n", wmain[i].r, aaa[i],
		       wmain[i].inwind, i);
	    }
	}
      else
	{

    /* JM 1411 -- added header for reading with astropy ascii module */
    fprintf (fptr, "x z var inwind i j\n");

	  for (i = 0; i < NDIM2; i++)
	    {
	      wind_n_to_ij (i, &ii, &jj);
	      fprintf (fptr, "%8.4e %8.4e %8.5e %3d %3d %3d\n",
		       wmain[i].xcen[0], wmain[i].xcen[2], aaa[i],
		       wmain[i].inwind, ii, jj);
	    }
	}

    }
  else if (choice == 2)
    {				// Then regrid to a linear grid
/* Wind array is organized so that it increments up along z and then over */

/* Want linearly spaced output array */
      rmin = zmin = 0;
      rmax = zmax = geo.rmax;

      xx[1] = 0.0;
      for (ii = 0; ii < ODIM; ii++)
	{
	  xx[0] = r = rmin + ii * (rmax - rmin) / (ODIM - 1);
	  for (jj = 0; jj < ODIM; jj++)
	    {
	      xx[2] = z = zmin + jj * (zmax - zmin) / (ODIM - 1);
	      //OLD 70B if (where_in_wind (xx) == 0)
	      if (where_in_wind (xx) >= 0)
		{		// Then the position is in the wind region
		  coord_fraction (0, xx, nnn, frac, &nelem);
		  for (nn = 0; nn < nelem; nn++)
		    {
		      aout[ii][jj] += aaa[nnn[nn]] * frac[nn];
		    }

		}
	    }
	}

      /* Now print out the array */
      // fprintf (fptr, "VARIABLES= \"X\" \"Z\" \"Var1\" \n");
      // fprintf (fptr, "ZONE I=%d J=%d DATAPACKING=POINT\n", ODIM, ODIM);
      fprintf (fptr, "# Resampled outputs\n");
      for (jj = 0; jj < ODIM; jj++)
	{
    /* JM 1411 -- added header for reading with astropy ascii module */
    fprintf (fptr, "r z var\n");

	  z = zmin + (zmax - zmin) * jj / (ODIM - 1);
	  for (ii = 0; ii < ODIM; ii++)
	    {
	      r = rmin + (rmax - rmin) * ii / (ODIM - 1);
	      fprintf (fptr, "%8.2e %8.2e %8.2e\n", r, z, aout[ii][jj]);
	    }
	}
      /*Finished printing out the array */
    }
  else
    {
      Error ("write_array: Unknown choice %d\n", choice);
    }

  fclose (fptr);
  return (0);
}


/**************************************************************************
                    Space Telescope Science Institute


  Synopsis:

  Description:
	This is a generalized display routine, intended to allow various
 	simplifications of py_wind

  Arguments:		

  Returns:

  Notes:

	For cylindrical coordinates, one goes up in z, and over in x. To cylindrical coordinates
	we want to display so that each colum represents a constant line on the axis
	w[mdim][ndim].  So for a system with mdim=20 and ndim=25.  There are 25 elements
	in the z axis and 20 elements in the xaxis

	Therefore a row, constant z, is displayed by by incrementing by 20 or mdim

	For rtheta coordinates, one goes around in theta, up in r.  The fasted moving
	coordinate is r (when thinking of 1-d versions of the array.  Unless we
	are going to write a separate routine, the simplest thing to do
	is to make each column represent a constant r

	w[mdim][ndim]  So, in spherical polar coorcinates, a system with mdim 20 has
	20 angles, and 25 radii.  

	So for spherical polar, constant r is displayed by incrementing by 1.  As a result
	it is unclear that one can easily use the same routine for the two situations
	since you seem to be incrementing the opposite axes, since in the one case one 
	wants MDIM rows and other case one wants NDIM rows.  

	It would be possible if you plot theta lines in each row, but what this means
	is that the first row is closest to the z axis


  History:
  	111125	ksl	Put standard headers before routine prior to updating slightly
	

 ************************************************************************/



int
display (name)
     char name[];
{
  int i, j, n;
  Log ("\n %s \n", name);
  Log ("z/theta \\x/r");
  for (i = py_wind_min; i < py_wind_max; i += py_wind_delta)
    if (geo.coord_type == 1 || py_wind_project == 1)
      Log ("%8.2e ", wmain[i * MDIM].x[0]);
    else
      Log ("%8.2e ", wmain[i * MDIM].rcen);

  Log ("\n");

  for (j = 0; j < MDIM; j++)
    {
      if (geo.coord_type == 1 || py_wind_project == 1)
	Log ("%8.2e ", wmain[j].x[2]);
      else
	Log ("%8.2e ", wmain[j].thetacen);

      for (i = py_wind_min; i < py_wind_max; i += py_wind_delta)
	{
	  n = i * MDIM + j;
	  Log ("%8.2g ", aaa[n]);
	}
      Log ("\n");
    }

  return (0);
}
