
/***********************************************************
                                       Space Telescope Science Institute

 Synopsis:
	py_wind is a program which can be used to display various parameters of a wind
	as calculated by python.  This file contains subroutines intended to
	display information about macro atoms
	used

Arguments:		


Returns:
 
Description:	
		
Notes:
	There two things that one would like to do: Display everything in a single
	cell and display a quantity in each cell.

  Note that level 3 here has index 2- we refer to levels by nconfig, so 0 
  is the first level.

History:
 	97jun	ksl	Coding on py_wind began.
	08mar	ksl	60 - Began this effort
	08jul	ss	61 - Modifications made to get more information
			about departure coeff
	08aug	ksl	62 - Brought back in line with concurrent changes
			which ksl had made to some of the
			routines, e.g do_partitions and saha
  1312 JM  Coded some new routines to report more information on level emissivities
	

**************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "atomic.h"
#include "python.h"


/**************************************************************************

  Synopsis:  

  xadiabatic_cooling_summary prints the adiabatic cooling in each cell,
  and the total adiabatic cooling. 
  
  Description:  

  Arguments:  

  Returns:

  Notes:

  History:

 ************************************************************************/
int
xadiabatic_cooling_summary (w, rootname, ochoice)
     WindPtr w;
     char rootname[];
     int ochoice;
{
  int n;
  double tot;
  double adiabatic_cooling ();
  char filename[LINELENGTH];
  double t_e;



  tot = 0.0;
  for (n = 0; n < NDIM2; n++)
    {
      aaa[n] = 0;
      if (w[n].vol > 0.0)
	{
	  t_e = plasmamain[w[n].nplasma].t_e;
	  num_recomb (&plasmamain[w[n].nplasma], t_e);
	  tot += aaa[n] = adiabatic_cooling (&w[n], t_e);
	}
    }


  display ("Adiabatic cooling");
  Log ("The total adiabatic cooling is %8.2g\n", tot);

  if (ochoice)
    {
      strcpy (filename, rootname);
      strcat (filename, ".adiabatic");
      write_array (filename, ochoice);
    }
  return (0);

}




/**************************************************************************

  Synopsis:  

  macro_summary is a routine intended to allow one to 
  display information about macro atoms. It is the top level routine
  for macro information and requires further input from the user.
  
  Description:  

  Arguments:  

  Returns:

  Notes:

  History:
 ************************************************************************/

int
macro_summary (w, rootname, ochoice)
     WindPtr w;
     char rootname[];
     int ochoice;
{
  int nmacro;
  int nconfig;
  int icell;
  int choose;
  int version;
  int nlev, nplasma;

  choose = 0;
  nlev = 0;
  nplasma=0;

  /* get input from the user */
  rdint ("Detailed cell info (0), levels (1) \
emissivities (2) P_escapes (3) Detailed Pops (4) taus (5) estimators (6)", &choose);

  if (choose == 0)
    {
      icell = 301;
      rdint ("Cell.no", &icell);



      ion_overview (icell);

      /* What configurations are macro atom configurations */
      nconfig = -1;
      rdint ("Configuration.Number (-1 for all b-values; -2 to rollup)",
	     &nconfig);

      while (nconfig >= -1)
	{
	  if (nconfig >= 0)
	    {
	      nmacro = config_overview (nconfig, icell);
	    }
	  else
	    {
	      depcoef_overview (icell);
	    }
	  rdint ("Configuration.Number(-1 for all b-values; -2 to rollup)",
		 &nconfig);
	}
    }

  /* specific level information for each cell */
  else if (choose == 1)
    {
      nconfig = -1;
      version = 0;
      rdint ("Number density (0), Fractional Pop (1) or Departure coeff (2)", &version);
      rdint ("Level.Number(-1 to rollup)", &nconfig);
      while (nconfig >= 0)
	{
	  depcoef_overview_specific (version, nconfig, w, rootname, ochoice);
	  rdint ("Level.Number(-1 to rollup)", &nconfig);
	}
    }


  /* JM 1311 -- new loop added to report level emissivities for each cell */
  else if (choose == 2)
    {
      rdint ("Level emissivity to view (0 - kpkts, -1 - back):", &nlev);
      while (nlev >= 0)
        {
          level_emissoverview (nlev, w, rootname, ochoice);
          rdint ("Level emissivity to view (0 - kpkts, -1 - back):", &nlev);
        }
    }

  /* JM 1311 -- new loop added to report escape probabilities 
     Note this currently only works for Balmer lines */
  else if (choose == 3)
    {
      rdint ("Upper level escape to view (-1 - back, Halpha = 3):", &nlev);
      while (nlev >= 0)
        {
          level_escapeoverview (nlev, w, rootname, ochoice);
          rdint ("Upper level escape to view (-1 - back, Halpha = 3):", &nlev);
        }
    }

  /* JM 1311 -- new loop added to report level populations in a cell */
  else if (choose == 4)
    {
      rdint ("Cell to view (-1 to rollup):", &nplasma);
      while (nplasma >= 0)
        {
          level_popsoverview (nplasma, w, rootname, ochoice);
          rdint ("Cell to view (-1 to rollup):", &nplasma);
        }
    }

  /* JM 1311 -- new loop added to report taus 
     Note this currently only works for Balmer lines */
  else if (choose == 5)
    {
      rdint ("Upper level tau to view (-1 - back, Halpha = 3):", &nlev);
      while (nlev >= 0)
        {
          level_tauoverview (nlev, w, rootname, ochoice);
          rdint ("Upper level tau to view (-1 - back, Halpha = 3):", &nlev);
        }
    }

  /* JM 1311 -- add stuff for estimators here */
  // else 
  //   {
  //     rdint ("Estimator to view: jbar (0) back (-1)", &version);
  //     while (version >= 0)
  //       {
  //         estimator_overview (version, w, rootname, ochoice);
  //         rdint ("Upper level tau to view (-1 - back, Halpha = 3):", &nlev);
  //       }
  //   }



  return (0);
}


/**************************************************************************

  Synopsis:  

  macro_summary is a routine intended to allow one to 
  display information about macro atoms. It is the top level routine
  for macro information and requires further input from the user.
  
  Description:  

  Arguments:  

  Returns:

  Notes:

  History:
 ************************************************************************/

int
ion_overview (icell)
     int icell;
{
  int n;
  IonPtr p;


  printf (" Z Ion I.Potent f.lte  n.lte f.nlte n.nlte  macro\n");
  for (n = 0; n < nions; n++)
    {
      p = &ion[n];
      printf ("%2d %2d %8.2e %5d %5d %6d %6d %5d\n", p->z, p->istate, p->ip,
	      p->firstlevel, p->nlevels, p->first_levden, p->nlte,
	      p->macro_info);
    }

  return (p->macro_info);
}

/**************************************************************************

  Synopsis:  

  macro_summary is a routine intended to allow one to 
  display information about macro atoms. It is the top level routine
  for macro information and requires further input from the user.
  
  Description:  

  Arguments:  

  Returns:

  Notes:

  History:
 ************************************************************************/

int
config_overview (n, icell)
     int n, icell;
{
  ConfigPtr p;
  PlasmaPtr x;
  MacroPtr m;
  double xden;
  int i, ii;

  p = &config[n];
  /* initialize the density */
  xden = -1;
  if (icell >= 0 && icell < NDIM2)
    {
      x = &plasmamain[icell];
      xden = x->levden[p->nden];
    }


  printf
    (" Z Ion nden macro n_bbu_jump  n_bbd_jump n_bfu_jump nbfd_jump  lev_den\n");

  printf ("%2d %2d %4d %5d %10d %10d %10d %10d %8.2e\n", p->z, p->istate,
	  p->nden, p->macro_info, p->n_bbu_jump, p->n_bbd_jump, p->n_bfu_jump,
	  p->n_bfd_jump, xden);

  if (icell >= 0 && icell < NDIM2 && p->macro_info == 1)
    {
      /* Print information about the jumps */
      printf ("bbu_jump: ");
      for (i = 0; i < p->n_bbu_jump; i++)
	printf (" %3d", p->bbu_jump[i]);
      printf ("\n");

      printf ("bbd_jump: ");
      for (i = 0; i < p->n_bbd_jump; i++)
	printf (" %3d", p->bbd_jump[i]);
      printf ("\n");

      printf ("bfu_jump: ");
      for (i = 0; i < p->n_bfu_jump; i++)
	printf (" %3d", p->bfu_jump[i]);
      printf ("\n");

      printf ("bfd_jump: ");
      for (i = 0; i < p->n_bfd_jump; i++)
	printf (" %3d", p->bfd_jump[i]);
      printf ("\n");

      m = &macromain[icell];
      printf ("matom_emis: %8.2e\n", m->matom_emiss[n]);
      printf ("matom_abs : %8.2e\n", m->matom_abs[n]);

      /* Detailed information on the bb transitions */
      printf ("bbu_jump:\n");
      for (i = 0; i < p->n_bbu_jump; i++)
	{
	  //ii=p->bbu_jump[i];
	  ii = i;
	  printf (" %3d %8.2e %8.2e\n", ii,
		  (m->jbar[config[n].bbu_indx_first + ii]),
		  (m->jbar_old[config[n].bbu_indx_first + ii]));
	}

      printf ("bbd_jump:\n");
      for (i = 0; i < p->n_bbd_jump; i++)
	{
	  //                          ii=p->bbd_jump[i];
	  ii = i;
	  printf (" %3d %8.2e %8.2e\n", ii,
		  (m->jbar[config[n].bbu_indx_first + ii]),
		  (m->jbar_old[config[n].bbu_indx_first + ii]));
	}

      /* Detailed information on the fb transitions */

      printf ("bfu_jump:\n");
      for (i = 0; i < p->n_bfu_jump; i++)
	{
	  //      ii=p->bfu_jump[i];
	  ii = i;
	  printf (" %3d %g %g %g %g %g %g %g %g\n", ii,
		  (m->gamma[config[n].bfu_indx_first + ii]),
		  (m->gamma_old[config[n].bfu_indx_first + ii]),
		  (m->alpha_st[config[n].bfu_indx_first + ii]),
		  (m->alpha_st_old[config[n].bfu_indx_first + ii]),
		  (m->alpha_st_e[config[n].bfu_indx_first + ii]),
		  (m->alpha_st_e_old[config[n].bfu_indx_first + ii]),
		  (m->recomb_sp[config[n].bfd_indx_first + ii]),
		  (m->recomb_sp_e[config[n].bfd_indx_first + ii]));
	}



      printf ("bfd_jump:\n");
      for (i = 0; i < p->n_bfd_jump; i++)
	{
	  //      ii=p->bfd_jump[i];
	  ii = i;
	  printf (" %3d %g %g %g %g %g %g %g %g\n", ii,
		  (m->gamma[config[n].bfu_indx_first + ii]),
		  (m->gamma_old[config[n].bfu_indx_first + ii]),
		  (m->alpha_st[config[n].bfu_indx_first + ii]),
		  (m->alpha_st_old[config[n].bfu_indx_first + ii]),
		  (m->alpha_st_e[config[n].bfu_indx_first + ii]),
		  (m->alpha_st_e_old[config[n].bfu_indx_first + ii]),
		  (m->recomb_sp[config[n].bfd_indx_first + ii]),
		  (m->recomb_sp_e[config[n].bfd_indx_first + ii]));
	}



    }
  else
    {
      printf ("This is not a macro level\n");
    }
  return (p->macro_info);
}



/**************************************************************************

  Synopsis:  

  depcoef_overview is a routine which prints information about departure 
  coefficients and level populations in a cell. There is probably a crossover
  between this and level_popsoverview
  
  Description:  

  Arguments:  

  Returns:

  Notes:

  History:
        080812  ksl 62 - Updated to reflect new versions of programs to do 
        partition functions and calculate saha densities
 ************************************************************************/

int
depcoef_overview (icell)
     int icell;
{
  ConfigPtr p;
  PlasmaPtr x, xdummy;
  double xden, lteden;
  int n;

  int copy_plasma ();
  plasma_dummy pdum;


  /* initialize the density */
  xden = -1;
  if (icell >= 0 && icell < NDIM2)
    {
      x = &plasmamain[icell];
      xdummy = &pdum;

      copy_plasma (x, xdummy);
      geo.macro_ioniz_mode = 0;

      partition_functions (xdummy, 1);
      saha (xdummy, xdummy->ne, xdummy->t_e);
      geo.macro_ioniz_mode = 1;
    }

  printf
    (" The cell is at (%g, %g, %g) and has volume %g.\n",
     wmain[x->nwind].xcen[0], wmain[x->nwind].xcen[1],
     wmain[x->nwind].xcen[2], wmain[x->nwind].vol);
  printf (" Z Ion nden macro  b       fpop    lte_fpop    t_e\n");

  for (n = 0; n < NLEVELS; n++)
    {
      p = &config[n];
      if (icell >= 0 && icell < NDIM2 && p->macro_info == 1)
	{
	  xden = den_config (x, n);
	  lteden = den_config (xdummy, n);
	  //  xden = x->levden[p->nden];
	  //lteden = xdummy->levden[p->nden];
	  printf ("%2d %2d %4d %5d %8.2e %8.2e %8.2e %8.2e\n", p->z,
		  p->istate, p->nden, p->macro_info, (xden / lteden), xden,
		  lteden, x->t_e);
	}
    }

  return (0);
}


/**************************************************************************

  Synopsis:  

  Routine to copy the necessary parts of a plasma structure for computing 
  a set of level populations. x1 points to the cell from which data is copied 
  and x2 points to the cell to which data is copied.
  
  Description:  

  Arguments:  

  Returns:

  Notes:

  History:
 ************************************************************************/


int
copy_plasma (x1, x2)
     PlasmaPtr x1, x2;
{
  x2->nwind = x1->nwind;
  x2->nplasma = x1->nplasma;
  x2->ne = x1->ne;
  x2->rho = x1->rho;
  x2->vol = x1->vol;
  x2->t_r = x1->t_r;
  x2->t_e = x1->t_e;
  x2->w = x1->w;

  /* JM 1409 -- added this for depcoef_overview_specific */
  x2->partition = x1->partition;
  x2->density = x1->density;

  /* Note this isn't everything in the cell! 
     Only the things needed for these routines */

  return (0);
}




/**************************************************************************

  Synopsis:  

  depcoef_overview_specific gives populations of nconfig in each cell
  
  Description:  

  Arguments: 
    version 0 gives number densities, 1 gives departure coefficients 

  Returns:

  Notes:

  History:
        080812  ksl 62 - Updated to reflect new versions of programs to do 
        partition functions and calculate saha densities
 ************************************************************************/

int
depcoef_overview_specific (version, nconfig, w, rootname, ochoice)
     int version;
     int nconfig;
     WindPtr w;
     char rootname[];
     int ochoice;
{
  int n;
  char filename[LINELENGTH], lname[LINELENGTH];
  ConfigPtr p;
  PlasmaPtr xplasma, xdummy;
  double xden, lteden, ion_density;
  int copy_plasma ();
  plasma_dummy pdum;


  nconfig = nconfig - 1;

  for (n = 0; n < NDIM2; n++)
    {
      aaa[n] = 0;
      if (w[n].vol > 0.0)
	{
	  xden = -1;

	  xplasma = &plasmamain[w[n].nplasma];
	  xdummy = &pdum;

	  copy_plasma (xplasma, xdummy);
	  geo.macro_ioniz_mode = 0;

	  partition_functions (xdummy, 1);
	  saha (xdummy, xdummy->ne, xdummy->t_e);

	  geo.macro_ioniz_mode = 1;

	  p = &config[nconfig];
	  xden = den_config (xplasma, nconfig);
	  lteden = den_config (xdummy, nconfig);

    ion_density = xplasma->density[p->nion];

	  if (version == 0)
	    {
	      aaa[n] = xden;
	    }
    else if (version == 1)
      {
        aaa[n] = xden / ion_density;
      }
	  else
	    {
	      aaa[n] = xden / lteden;
	    }
	}
    }

  if (version == 0)
    display ("Level Population Number Densities");
  else if (version == 1)
    display ("Fractional Level Population");
  else
    display ("Departure Coefficient");

  if (ochoice)
    {
      strcpy (filename, rootname);

      sprintf (lname, ".lev%d_", nconfig+1);
      strcat (filename, lname);

      if (version == 0)
        strcat (filename, "den");
      else if (version == 1)
        strcat (filename, "frac_pop");
      else
        strcat (filename, "dep_coef");

      write_array (filename, ochoice);
    }

  return (0);
}






/**************************************************************************

  Synopsis:  

  level_popsoverview is a routine which prints information about departure 
  coefficients and level populations in a cell. There is probably a crossover
  between this and depcoef_overview 
  
  Description:  

  Arguments:  

  Returns:

  Notes:

  History:
      1312  JM  Coded
 ************************************************************************/


int 
level_popsoverview (nplasma, w, rootname, ochoice)
  int nplasma;
  WindPtr w;
  char rootname[];
  int ochoice;
{
  int i;
  PlasmaPtr xplasma, xdummy;
  plasma_dummy pdum;
  char filename[LINELENGTH];
  char lname[LINELENGTH];
  double xden, lteden;
  int copy_plasma ();

  strcpy (filename, rootname);
  strcpy (filename, rootname);
  strcat (filename, ".lev_cell");

  sprintf (lname, "%d", nplasma);
  strcat (filename, lname);
  strcat(filename,".dat");

  FILE *f = fopen(filename, "w");

  if (f == NULL)
  {
      printf("Error opening file!\n");
      exit(1);
  }

  xplasma = &plasmamain[nplasma];
  xdummy = &pdum;

  copy_plasma (xplasma, xdummy);
  geo.macro_ioniz_mode = 0;
  
  Log("# Level, Pops, Dep coefs\n");
  if (ochoice)
    fprintf(f, "# Level, Pops, Dep coefs\n");

  Log("# Cell %i Radiation field t_r %8.4e w %8.4e\n",
        nplasma, xplasma->t_r, xplasma->w);
  Log("# Cell %i Physical t_e %8.4e ne %8.4e\n",
      nplasma, xplasma->t_e, xplasma->ne);

  if (ochoice)
  {
    fprintf(f, "# Cell %i Radiation field t_r %8.4e w %8.4e\n",
        nplasma, xplasma->t_r, xplasma->w);
    fprintf(f, "# Cell %i Physical t_e %8.4e ne %8.4e\n",
      nplasma, xplasma->t_e, xplasma->ne);
  }

  for (i=0; i < nlevels_macro; i++)
  {
    partition_functions (xdummy, 1);
    saha (xdummy, xdummy->ne, xdummy->t_e);

    geo.macro_ioniz_mode = 1;

    //p = &config[i];
    xden = den_config (xplasma, i);
    lteden = den_config (xdummy, i);
    Log("%i %8.4e %8.4e\n", i+1, xplasma->levden[i], xden/lteden);
    if (ochoice)
      fprintf(f, "%i %8.4e %8.4e\n", i+1, xplasma->levden[i], xden/lteden);
  }
  fclose(f);

  return (0);
}




/**************************************************************************

  Synopsis:  

  Routine to print level emissivities in each cell.
  The emissivity is the quantity
  A21 * n2 * beta_sobolev * nu_21 * H * volume
  And thus has units of ergs/s
  
  Description:  

  Arguments:  

  Returns:

  Notes:

  History:
      1312  JM  Coded to track 
                e.g. H-alpha emission through wind
 ************************************************************************/


int 
level_emissoverview (nlev, w, rootname, ochoice)
  int nlev;
  WindPtr w;
  char rootname[];
  int ochoice;
{
  int n, nplasma;
  char name[LINELENGTH], lname[LINELENGTH];
  char filename[LINELENGTH];

  strcpy (name, "");

  if (nlev != 0)
    {
      sprintf (name, "Level emissivities for Level %i", nlev);
    }
  else
    {
      sprintf (name, "Kpkt emissivities");
    }


  for (n = 0; n < NDIM2; n++)
    {
      aaa[n] = 0;
      nplasma = w[n].nplasma;

      if (w[n].vol > 0.0 && plasmamain[nplasma].ne > 1.0)
        {
          if (nlev != 0)
            {
              aaa[n] = macromain[nplasma].matom_emiss[nlev - 1];
            }
          else
            {
              aaa[n] = plasmamain[nplasma].kpkt_emiss;
            }
        }
    }

  display(name);

  if (ochoice)
    {
      strcpy (filename, rootname);
      strcat (filename, ".lev");

      if (nlev != 0)
        {
          sprintf (lname, "%d", nlev);
          strcat (filename, lname);
        }
      else
        {
          strcat (filename, "k");
        }

      strcat (filename, "_emiss");
      write_array (filename, ochoice);

    }

  return (0);
}




/**************************************************************************

  Synopsis:  

  Routine to print escape probabilities for a given 
  Balmer line in each cell.

  Description:  

  Arguments:  

  Returns:

  Notes:

  History:
      1312  JM  Coded
 ************************************************************************/



int level_escapeoverview (nlev, w, rootname, ochoice)
  int nlev;
  WindPtr w;
  char rootname[];
  int ochoice;
{
  PlasmaPtr xplasma;
  int n, nplasma, nline, found;
  char name[LINELENGTH], lname[LINELENGTH];
  char filename[LINELENGTH];
  double lambda;

  
  found = 0;
  nline = 0;



  /* Find the line in line list */
  while (found != 1 && nline < nlines)
  {

    if (lin_ptr[nline]->levl == 2 && lin_ptr[nline]->levu == nlev && 
       lin_ptr[nline]->z == 1 && lin_ptr[nline]->istate == 1)
      { 
        found = 1;
      }

    nline++; 
  }


  if (nline == nlines)
    {
      Error ("level_escapeoverview: Could not find line in linelist\n");
      exit (0);
    }

  nline--;
  
  lambda = ( C/ lin_ptr[nline]->freq ) * 1e8;

  strcpy (name, "");
  sprintf (name, "Balmer series P_escapes for Level %i, Lambda %.1f", nlev, lambda);


  for (n = 0; n < NDIM2; n++)
    {
      aaa[n] = 0;
      nplasma = w[n].nplasma;

      if (w[n].vol > 0.0 && plasmamain[nplasma].ne > 1.0)
        {
          xplasma = &plasmamain[nplasma];
          aaa[n] = p_escape (lin_ptr[nline], xplasma);

          
        }
    }

  display(name);

  if (ochoice)
    {
      strcpy (filename, rootname);
      strcat (filename, ".lev");

      sprintf (lname, "%d", nlev);
      strcat (filename, lname);

      strcat (filename, "_esc");
      write_array (filename, ochoice);

    }

  return (0);
}




/**************************************************************************

  Synopsis:  

  Routine to print taus for a given 
  Balmer line in each cell.

  Description:  

  Arguments:  

  Returns:

  Notes:

  History:
      1312  JM  Coded
 ************************************************************************/



int level_tauoverview (nlev, w, rootname, ochoice)
  int nlev;
  WindPtr w;
  char rootname[];
  int ochoice;
{
  PlasmaPtr xplasma;
  WindPtr one;
  int n, nplasma, nline, found;
  char name[LINELENGTH], lname[LINELENGTH];
  char filename[LINELENGTH];
  double lambda;

  
  found = 0;
  nline = 0;



  /* Find the line in line list */
  while (found != 1 && nline < nlines)
  {

    if (lin_ptr[nline]->levl == 2 && lin_ptr[nline]->levu == nlev && 
       lin_ptr[nline]->z == 1 && lin_ptr[nline]->istate == 1)
      { 
        found = 1;
      }

    nline++; 
  }


  if (nline == nlines)
    {
      Error ("level_tauoverview: Could not find line in linelist\n");
      exit (0);
    }

  nline--;
  
  lambda = ( C/ lin_ptr[nline]->freq ) * 1e8;

  strcpy (name, "");
  sprintf (name, "Balmer series taus for Level %i, Lambda %.1f", nlev, lambda);


  for (n = 0; n < NDIM2; n++)
    {
      aaa[n] = 0;
      nplasma = w[n].nplasma;

      if (w[n].vol > 0.0 && plasmamain[nplasma].ne > 1.0)
        {
          xplasma = &plasmamain[nplasma];
          one = &wmain[xplasma->nwind];
          aaa[n] = sobolev (one, one->x, xplasma->density[lin_ptr[nline]->nion], lin_ptr[nline], one->dvds_ave);

          
        }
    }

  display(name);

  if (ochoice)
    {
      strcpy (filename, rootname);
      strcat (filename, ".lev");

      sprintf (lname, "%d", nlev);
      strcat (filename, lname);

      strcat (filename, "_tau");
      write_array (filename, ochoice);

    }

  return (0);
}
