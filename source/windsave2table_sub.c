
/***********************************************************/
/** @file  windsave2table_sub.c
 * @author ksl
 * @date   April, 2018
 *
 * @brief  Subroutines for windsave2table
 *
 * These are various routines used by windsave2table which reads
 * a windsave file and writes various variables of that windsave
 * file to ascii files which can be read with astropy.io ascii
 * as tables.
 *
 * Unlike py_wind, windsave2table is hardwired and to change what
 * is written one must actually modify the routines themselves.
 ***********************************************************/


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "atomic.h"
#include "python.h"




/**********************************************************/
/**
 * @brief      The main routine associated with windsave2table, this routine calls
 * various other routines which write individual files
 *
 * @param [in] char *  root   The rootname of the windsave file
 * @return     Always returns 0
 *
 * @details
 *
 *
 *
 * ### Notes ###
 *
 * The routine cycles through the various domains and calles subroutines
 * that write the individual files for each domain.
 *
 *
 *
 **********************************************************/

int
do_windsave2table (root, ion_switch)
     char *root;
     int ion_switch;
{
  int ndom;
  char rootname[LINELENGTH];


  for (ndom = 0; ndom < geo.ndomain; ndom++)
  {

    sprintf (rootname, "%s.%d", root, ndom);

    create_master_table (ndom, rootname);
    create_heat_table (ndom, rootname);
    create_ion_table (ndom, rootname, 1, ion_switch);
    create_ion_table (ndom, rootname, 2, ion_switch);
    create_ion_table (ndom, rootname, 6, ion_switch);
    create_ion_table (ndom, rootname, 7, ion_switch);
    create_ion_table (ndom, rootname, 8, ion_switch);
    create_ion_table (ndom, rootname, 11, ion_switch);
    create_ion_table (ndom, rootname, 14, ion_switch);
    create_ion_table (ndom, rootname, 26, ion_switch);
  }
  return (0);
}






/**********************************************************/
/**
 * @brief      writes specific  variables of a windsaave
 * which are intended to be of general interest to
 * file which has the format of an astropy table
 *
 * 	It is intended to be easily modifible.
 *
 * @param [in] int  ndom   A domain number
 * @param [in] char  rootname   The rootname of the master file
 * @return   Always returns 0
 *
 * @details
 *
 * The master table is contains basic information for each
 * cell in the wind, such as the electron density, the density,
 * the ionization parameter, and the radiative and electron temperature
 *
 * The routine takes data directly from wmain, and then calls
 * get_one or get_ion multiple times to fet info from the Plasma
 * structure.
 *
 * It then writes the data to an ascii file which can be read as
 * an  astropy table
 *
 * ### Notes ###
 * To add a variable one just needs to define the column_name
 * and send the appropriate call to either get_one or get_ion.
 *
 * There is some duplicated code in the routine that pertains
 * to whether one is dealing with a spherecial or a 2d coordinate
 * system.  It should be possible to delete this
 *
 **********************************************************/

int
create_master_table (ndom, rootname)
     int ndom;
     char rootname[];
{
  char filename[132];
  double *c[50], *converge;
  char column_name[50][20];
  char one_line[1024], start[132], one_value[20];


  int i, ii, jj;
  int nstart, ndim2;
  int n, ncols;
  FILE *fptr;

  strcpy (filename, rootname);
  strcat (filename, ".master.txt");


  fptr = fopen (filename, "w");

  /* Get the variables that one needs */

  c[0] = get_one (ndom, "vol");
  strcpy (column_name[0], "vol");

  c[1] = get_one (ndom, "rho");
  strcpy (column_name[1], "rho");

  c[2] = get_one (ndom, "ne");
  strcpy (column_name[2], "ne");

  c[3] = get_one (ndom, "t_e");
  strcpy (column_name[3], "t_e");

  c[4] = get_one (ndom, "t_r");
  strcpy (column_name[4], "t_r");

  c[5] = get_ion (ndom, 1, 1, 0);
  strcpy (column_name[5], "h1");

  c[6] = get_ion (ndom, 2, 2, 0);
  strcpy (column_name[6], "he2");

  c[7] = get_ion (ndom, 6, 4, 0);
  strcpy (column_name[7], "c4");

  c[8] = get_ion (ndom, 7, 5, 0);
  strcpy (column_name[8], "n5");

  c[9] = get_ion (ndom, 8, 6, 0);
  strcpy (column_name[9], "o6");

  c[10] = get_one (ndom, "dmo_dt_x");
  strcpy (column_name[10], "dmo_dt_x");


  c[11] = get_one (ndom, "dmo_dt_y");
  strcpy (column_name[11], "dmo_dt_y");

  c[12] = get_one (ndom, "dmo_dt_z");
  strcpy (column_name[12], "dmo_dt_z");

  c[13] = get_one (ndom, "ip");
  strcpy (column_name[13], "ip");

  c[14] = get_one (ndom, "ntot");
  strcpy (column_name[14], "ntot");

  c[15] = get_one (ndom, "nrad");
  strcpy (column_name[15], "nrad");

  c[16] = get_one (ndom, "nioniz");
  strcpy (column_name[16], "nioniz");


  /* This should be the maxium number above +1 */
  ncols = 17;


  converge = get_one (ndom, "converge");

  /* At this point oll of the data has been collected */


  nstart = zdom[ndom].nstart;
  ndim2 = zdom[ndom].ndim2;


  if (zdom[ndom].coord_type == SPHERICAL)
  {


    /* First assemble the header line
     */

    sprintf (start, "%9s %9s %4s %6s %6s %8s %8s %8s ", "r", "rcen", "i", "inwind", "converge", "v_x", "v_y", "v_z");
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

    for (i = 0; i < ndim2; i++)
    {
      // This line is different from the two d case
      sprintf (start, "%9.3e %9.3e %4d %6d %8.0f %8.2e %8.2e %8.2e ",
               wmain[nstart + i].r, wmain[nstart + i].rcen, i, wmain[nstart + i].inwind,
               converge[i], wmain[nstart + i].v[0], wmain[nstart + i].v[1], wmain[nstart + i].v[2]);
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
  else if (zdom[ndom].coord_type == CYLIND)
  {

    /* First assemble the header line */

    sprintf (start, "%8s %8s %8s %8s %4s %4s %6s %8s %8s %8s %8s ", "x", "z", "xcen", "zcen", "i", "j", "inwind", "converge", "v_x", "v_y",
             "v_z");
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

    for (i = 0; i < ndim2; i++)
    {
      wind_n_to_ij (ndom, nstart + i, &ii, &jj);
      sprintf (start,
               "%8.2e %8.2e %8.2e %8.2e %4d %4d %6d %8.0f %8.2e %8.2e %8.2e ",
               wmain[nstart + i].x[0], wmain[nstart + i].x[2], wmain[nstart + i].xcen[0], wmain[nstart + i].xcen[2], ii,
               jj, wmain[nstart + i].inwind, converge[i], wmain[nstart + i].v[0], wmain[nstart + i].v[1], wmain[nstart + i].v[2]);
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

  else if (zdom[ndom].coord_type == RTHETA)
  {

    /* First assemble the header line */

    sprintf (start, "%8s %8s %8s %9s %8s %8s %8s %8s %4s %4s %6s %8s %8s %8s %8s ", "r", "theta", "r_cen", "theta_cen", "x", "z", "xcen",
             "zcen", "i", "j", "inwind", "converge", "v_x", "v_y", "v_z");
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

    for (i = 0; i < ndim2; i++)
    {
      wind_n_to_ij (ndom, nstart + i, &ii, &jj);
      sprintf (start,
               "%8.2e %8.2e %8.2e %9.2e %8.2e %8.2e %8.2e %8.2e %4d %4d %6d %8.0f %8.2e %8.2e %8.2e ",
               wmain[nstart + i].r, wmain[nstart + i].theta, wmain[nstart + i].rcen, wmain[nstart + i].thetacen,
               wmain[nstart + i].x[0], wmain[nstart + i].x[2], wmain[nstart + i].xcen[0], wmain[nstart + i].xcen[2], ii,
               jj, wmain[nstart + i].inwind, converge[i], wmain[nstart + i].v[0], wmain[nstart + i].v[1], wmain[nstart + i].v[2]);
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
    printf ("Error: Cannot print out files for coordinate system type %d\n", zdom[ndom].coord_type);
  }


  return (0);
}




/**********************************************************/
/**
 * @brief      writes a selected variables related to heating and cooling
 * processes to an ascii file which can be read as an astropy table
 *
 *
 * @param [in] int  ndom   The domain of interest
 * @param [in, out] char  rootname[]   The rootname of the windsave file
 * @return     Always returns 0
 *
 * @details
 *
 * ### Notes ###
 *
 * To add a variable one just needs to define the column_name
 * and send the appropriate call to either get_one or get_ion.
 *
 *
 **********************************************************/

int
create_heat_table (ndom, rootname)
     int ndom;
     char rootname[];
{
  char filename[132];
  double *c[50], *converge;
  char column_name[50][20];
  char one_line[1024], start[132], one_value[20];


  int i, ii, jj;
  int nstart, ndim2;
  int n, ncols;
  FILE *fptr;

  strcpy (filename, rootname);
  strcat (filename, ".heat.txt");


  fptr = fopen (filename, "w");

  /* Get the variables that one needs */

  c[0] = get_one (ndom, "vol");
  strcpy (column_name[0], "vol");

  c[1] = get_one (ndom, "rho");
  strcpy (column_name[1], "rho");

  c[2] = get_one (ndom, "ne");
  strcpy (column_name[2], "ne");

  c[3] = get_one (ndom, "t_e");
  strcpy (column_name[3], "t_e");

  c[4] = get_one (ndom, "t_r");
  strcpy (column_name[4], "t_r");

  c[5] = get_one (ndom, "w");
  strcpy (column_name[5], "w");

  c[6] = get_one (ndom, "heat_tot");
  strcpy (column_name[6], "heat_tot");

  c[7] = get_one (ndom, "heat_comp");
  strcpy (column_name[7], "heat_comp");

  c[8] = get_one (ndom, "heat_lines");
  strcpy (column_name[8], "heat_lines");

  c[9] = get_one (ndom, "heat_ff");
  strcpy (column_name[9], "heat_ff");

  c[10] = get_one (ndom, "heat_photo");
  strcpy (column_name[10], "heat_photo");

  c[11] = get_one (ndom, "heat_auger");
  strcpy (column_name[11], "heat_auger");

  c[12] = get_one (ndom, "cool_tot");
  strcpy (column_name[12], "cool_tot");

  c[13] = get_one (ndom, "cool_comp");
  strcpy (column_name[13], "cool_comp");

  c[14] = get_one (ndom, "lum_lines");
  strcpy (column_name[14], "lum_lines");

  c[15] = get_one (ndom, "cool_dr");
  strcpy (column_name[15], "cool_dr");

  c[16] = get_one (ndom, "lum_ff");
  strcpy (column_name[16], "lum_ff");


  c[17] = get_one (ndom, "cool_rr");
  strcpy (column_name[17], "cool_rr");

  c[18] = get_one (ndom, "cool_adiab");
  strcpy (column_name[18], "cool_adiab");

  c[19] = get_one (ndom, "heat_shock");
  strcpy (column_name[19], "heat_shock");

  c[20] = get_one (ndom, "heat_lines_macro");
  strcpy (column_name[20], "heat_lines_macro");

  c[21] = get_one (ndom, "heat_photo_macro");
  strcpy (column_name[21], "heat_photo_macro");

  /* This should be the maxium number above +1 */
  ncols = 22;


  converge = get_one (ndom, "converge");

  /* At this point all of the data has been collected */


  nstart = zdom[ndom].nstart;
  ndim2 = zdom[ndom].ndim2;


  if (zdom[ndom].coord_type == SPHERICAL)
  {


    /* First assemble the header line
     */

    sprintf (start, "%9s %4s %6s %6s %8s %8s %8s ", "r", "i", "inwind", "converge", "v_x", "v_y", "v_z");
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

    for (i = 0; i < ndim2; i++)
    {
      // This line is different from the two d case
      sprintf (start, "%9.3e %4d %6d %8.0f %8.2e %8.2e %8.2e ",
               wmain[nstart + i].r, i, wmain[nstart + i].inwind,
               converge[i], wmain[nstart + i].v[0], wmain[nstart + i].v[1], wmain[nstart + i].v[2]);
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

    /* First assemble the header line */

    sprintf (start, "%8s %8s %4s %4s %6s %8s %8s %8s %8s ", "x", "z", "i", "j", "inwind", "converge", "v_x", "v_y", "v_z");
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

    for (i = 0; i < ndim2; i++)
    {
      wind_n_to_ij (ndom, nstart + i, &ii, &jj);
      sprintf (start,
               "%8.2e %8.2e %4d %4d %6d %8.0f %8.2e %8.2e %8.2e ",
               wmain[nstart + i].xcen[0], wmain[nstart + i].xcen[2], ii,
               jj, wmain[nstart + i].inwind, converge[i], wmain[nstart + i].v[0], wmain[nstart + i].v[1], wmain[nstart + i].v[2]);
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




/**********************************************************/
/**
 * @brief      makes an astropy table containing the relative abundances
 *  of a given element as a function of the position in the grid
 *
 * @param [in] int  ndom   The domain number
 * @param [in] char  rootname[]   The rootname of the windsave file
 * @param [in] int  iz   atomic number of the element
 * @return     0 on completion
 *
 * @details
 *
 * ### Notes ###
 *
 * The routine calls get_ion_density mulitiple times, once for
 * each ion of an atom
 *
 **********************************************************/

int
create_ion_table (ndom, rootname, iz, ion_switch)
     int ndom;
     char rootname[];
     int iz;                    // Where z is the element
     int ion_switch;            // Determines what is actually printed out
{
  char filename[132];
  double *c[50];
  int first_ion, number_ions;
  char element_name[20];
  int istate[50];
  char one_line[1024], start[132], one_value[20];
  int nstart, ndim2;


  int i, ii, jj, n;
  FILE *fptr;

/* First we actually need to determine what ions exits, but we will ignore this for now */

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

// Log ("element %d %d %s\n", first_ion, number_ions, element_name);

/* Open the output file */

  sprintf (filename, "%s.%s.txt", rootname, element_name);

  fptr = fopen (filename, "w");



  i = 0;
  while (i < number_ions)
  {
    istate[i] = ion[first_ion + i].istate;

    c[i] = get_ion (ndom, iz, istate[i], ion_switch);
    i++;
  }

  nstart = zdom[ndom].nstart;
  ndim2 = zdom[ndom].ndim2;



  if (zdom[ndom].coord_type == SPHERICAL)
  {


    /* First assemble the header line
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

    for (i = 0; i < ndim2; i++)
    {
      // This line is different from the two d case
      sprintf (start, "%8.2e %4d %6d ", wmain[nstart + i].r, i, wmain[nstart + i].inwind);
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


    /* First assemble the header line */

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

    for (i = 0; i < ndim2; i++)
    {
      wind_n_to_ij (ndom, nstart + i, &ii, &jj);
      sprintf (start, "%8.2e %8.2e %4d %4d %6d ", wmain[nstart + i].xcen[0], wmain[nstart + i].xcen[2], ii, jj, wmain[nstart + i].inwind);
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



/**********************************************************/
/**
 * @brief      Get get density, etc for one particular ion
 *
 * @param [in] int  ndom   the domain number
 * @param [in] int  element   the element number
 * @param [in] int  istate   the ionization state
 * @param [in] int  iswitch   a switch controlling exactly what is returned for that ion
 * @return     Normally returns an array with values associated with what is requested
 *    	This will return an array with all zeros if there is no such ion
 *
 * @details
 *
 * ### Notes ###
 *
 **********************************************************/

double *
get_ion (ndom, element, istate, iswitch)
     int ndom, element, istate, iswitch;
{
  int nion, nelem;
  int n;
//OLD  char name[LINELENGTH];
  int nplasma;
  double *x;
  int nstart, ndim2;
  double nh;


  nstart = zdom[ndom].nstart;
  ndim2 = zdom[ndom].ndim2;

  x = (double *) calloc (sizeof (double), ndim2);

/* Find the ion */

  nion = 0;
  while (nion < nions && !(ion[nion].z == element && ion[nion].istate == istate))
    nion++;
  if (nion == nions)
  {
    Log ("Error--element %d ion %d not found in define_wind\n", element, istate);
    return (x);
  }
  nelem = 0;
  while (nelem < nelements && ele[nelem].z != element)
    nelem++;

//OLD  strcpy (name, "");

  /* Now populate the array */

  for (n = 0; n < ndim2; n++)
  {
    x[n] = 0;
    nplasma = wmain[nstart + n].nplasma;
    if (wmain[nstart + n].vol > 0.0 && plasmamain[nplasma].rho > 0.0)
    {
      if (iswitch == 0)
      {
//OLD        sprintf (name, "Element %d (%s) ion %d fractions\n", element, ele[nelem].name, istate);
        x[n] = plasmamain[nplasma].density[nion];
        nh = rho2nh * plasmamain[nplasma].rho;
        x[n] /= (nh * ele[nelem].abun);
      }
      else if (iswitch == 1)
      {
//OLD        sprintf (name, "Element %d (%s) ion %d density\n", element, ele[nelem].name, istate);
        x[n] = plasmamain[nplasma].density[nion];
      }
      else if (iswitch == 2)
      {
//OLD        sprintf (name, "Element %d (%s) ion %d  #scatters\n", element, ele[nelem].name, istate);
        x[n] = (double) plasmamain[nplasma].scatters[nion] / plasmamain[nplasma].vol;
      }
      else if (iswitch == 3)
      {
//OLD        sprintf (name, "Element %d (%s) ion %d scattered flux\n", element, ele[nelem].name, istate);
        x[n] = plasmamain[nplasma].xscatters[nion];
      }
      else
      {
        Error ("get_ion : Unknown switch %d \n", iswitch);
        exit (0);
      }
    }
  }

  return (x);
}



/**********************************************************/
/**
 * @brief      Get a simple variable from the PlasmaPtr array
 *
 * @param [in] int  ndom   The domain in question
 * @param [in] char  variable_name[]   The name of the variable
 * @return     The values in the plasma pointer for this variable. A double
 * 	will be returned even if the PlasmaPtr varible is an integer
 *
 *
 * @details
 *
 * A simple variable is a variable that is just a number, not an array
 *
 * The routine performes a simple tranlation of the character name
 * to a variable in the PlasmaPtr.
 *
 * ### Notes ###
 * Only selected variables are returned, but new variables are easy
 * to add using the template of the other variables
 *
 **********************************************************/

double *
get_one (ndom, variable_name)
     int ndom;
     char variable_name[];
{
  int n;
  int nplasma;
  double *x;
  int ndim2;
  int nstart;

  nstart = zdom[ndom].nstart;
  ndim2 = zdom[ndom].ndim2;

  x = (double *) calloc (sizeof (double), ndim2);

  for (n = 0; n < ndim2; n++)
  {
    x[n] = 0;
    if (wmain[n + nstart].vol > 0.0)
    {
      nplasma = wmain[n + nstart].nplasma;


      if (strcmp (variable_name, "ne") == 0)
      {
        x[n] = plasmamain[nplasma].ne;
      }
      else if (strcmp (variable_name, "rho") == 0)
      {
        x[n] = plasmamain[nplasma].rho;
      }
      else if (strcmp (variable_name, "vol") == 0)
      {
        x[n] = plasmamain[nplasma].vol;
      }
      else if (strcmp (variable_name, "t_e") == 0)
      {
        x[n] = plasmamain[nplasma].t_e;
      }
      else if (strcmp (variable_name, "t_r") == 0)
      {
        x[n] = plasmamain[nplasma].t_r;
      }

      else if (strcmp (variable_name, "converge") == 0)
      {
        x[n] = plasmamain[nplasma].converge_whole;
      }
      else if (strcmp (variable_name, "dmo_dt_x") == 0)
      {
        x[n] = plasmamain[nplasma].dmo_dt[0];
      }
      else if (strcmp (variable_name, "dmo_dt_y") == 0)
      {
        x[n] = plasmamain[nplasma].dmo_dt[1];
      }
      else if (strcmp (variable_name, "dmo_dt_z") == 0)
      {
        x[n] = plasmamain[nplasma].dmo_dt[2];
      }
      else if (strcmp (variable_name, "ntot") == 0)
      {
        x[n] = plasmamain[nplasma].ntot;
      }
      else if (strcmp (variable_name, "ip") == 0)
      {
        x[n] = plasmamain[nplasma].ip;
      }
      else if (strcmp (variable_name, "heat_tot") == 0)
      {
        x[n] = plasmamain[nplasma].heat_tot;
      }
      else if (strcmp (variable_name, "heat_comp") == 0)
      {
        x[n] = plasmamain[nplasma].heat_comp;
      }
      else if (strcmp (variable_name, "heat_lines") == 0)
      {
        x[n] = plasmamain[nplasma].heat_lines;
      }
      else if (strcmp (variable_name, "heat_ff") == 0)
      {
        x[n] = plasmamain[nplasma].heat_ff;
      }
      else if (strcmp (variable_name, "heat_photo") == 0)
      {
        x[n] = plasmamain[nplasma].heat_photo;
      }
      else if (strcmp (variable_name, "heat_auger") == 0)
      {
        x[n] = plasmamain[nplasma].heat_auger;
      }
      else if (strcmp (variable_name, "cool_comp") == 0)
      {
        x[n] = plasmamain[nplasma].cool_comp;
      }
      else if (strcmp (variable_name, "lum_lines") == 0)
      {
        x[n] = plasmamain[nplasma].lum_lines;
      }
      else if (strcmp (variable_name, "lum_ff") == 0)
      {
        x[n] = plasmamain[nplasma].lum_ff;
      }
      else if (strcmp (variable_name, "cool_rr") == 0)
      {
        x[n] = plasmamain[nplasma].cool_rr;
      }
      else if (strcmp (variable_name, "cool_dr") == 0)
      {
        x[n] = plasmamain[nplasma].cool_dr;
      }
      else if (strcmp (variable_name, "cool_tot") == 0)
      {
        x[n] = plasmamain[nplasma].cool_tot;
      }
      else if (strcmp (variable_name, "w") == 0)
      {
        x[n] = plasmamain[nplasma].w;
      }
      else if (strcmp (variable_name, "nrad") == 0)
      {
        x[n] = plasmamain[nplasma].nrad;
      }
      else if (strcmp (variable_name, "nioniz") == 0)
      {
        x[n] = plasmamain[nplasma].nioniz;
      }
      else if (strcmp (variable_name, "heat_shock") == 0)
      {
        x[n] = plasmamain[nplasma].heat_shock;
      }
      else if (strcmp (variable_name, "cool_adiab") == 0)
      {
        x[n] = plasmamain[nplasma].cool_adiabatic;
      }
      else if (strcmp (variable_name, "heat_lines_macro") == 0)
      {
        x[n] = plasmamain[nplasma].heat_lines_macro;
      }
      else if (strcmp (variable_name, "heat_photo_macro") == 0)
      {
        x[n] = plasmamain[nplasma].heat_photo_macro;
      }


      else
      {
        Error ("get_one: Unknown variable %s\n", variable_name);
      }
    }
  }

  return (x);

}
