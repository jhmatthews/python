
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
	There are probably two things that one would like to do: Display everything in a single
	cell and display 

History:
 	97jun	ksl	Coding on py_wind began.
	08mar	ksl	60 - Began this effort
	08jul	ss	61 - Modifications made to get more information
			about departure coeff
	08aug	ksl	62 - Brought back in line with concurrent changes
			which ksl had made to some of the
			routines, e.g do_partitions and saha
	

**************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "atomic.h"
#include "python.h"



/* A summary of adiabatic cooling */
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


/* A routine intended to allow one to display information about macro atoms */
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

  choose = 0;
  rdint ("Detailed cell info (0) or specific level info (other)", &choose);

  if (choose == 0)
    {
      icell = 301;
      rdint ("Cell.no", &icell);



      ion_overview (icell);
      /* What configurations are macro atom configurations */

      nconfig = -1;
      rdint ("Configuration.Number(-1 for all b-values; -2 to rollup)",
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
  else
    {
      nconfig = -1;
      version = 0;
      rdint ("Configuration.Number(-1 to rollup)", &nconfig);
      rdint ("Departure coeff (1) or number density (0)", &version);
      while (nconfig >= 0)
	{
	  depcoef_overview_specific (version, nconfig, w, rootname, ochoice);
	  rdint ("Configuration.Number(-1 to rollup)", &nconfig);
	}
    }



  return (0);
}


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
	  printf (" %3d %8.2e %8.2e\n", ii, (m->jbar[n][ii]),
		  (m->jbar_old[n][ii]));
	}

      printf ("bbd_jump:\n");
      for (i = 0; i < p->n_bbd_jump; i++)
	{
	  //                          ii=p->bbd_jump[i];
	  ii = i;
	  printf (" %3d %8.2e %8.2e\n", ii, (m->jbar[n][ii]),
		  (m->jbar_old[n][ii]));
	}

      /* Detailed information on the fb transitions */

      printf ("bfu_jump:\n");
      for (i = 0; i < p->n_bfu_jump; i++)
	{
	  //      ii=p->bfu_jump[i];
	  ii = i;
	  printf (" %3d %g %g %g %g %g %g %g %g\n", ii, (m->gamma[n][ii]),
		  (m->gamma_old[n][ii]), (m->alpha_st[n][ii]),
		  (m->alpha_st_old[n][ii]), (m->alpha_st_e[n][ii]),
		  (m->alpha_st_e_old[n][ii]), (m->recomb_sp[n][ii]),
		  (m->recomb_sp_e[n][ii]));
	}



      printf ("bfd_jump:\n");
      for (i = 0; i < p->n_bfd_jump; i++)
	{
	  //      ii=p->bfd_jump[i];
	  ii = i;
	  printf (" %3d %g %g %g %g %g %g %g %g\n", ii, (m->gamma[n][ii]),
		  (m->gamma_old[n][ii]), (m->alpha_st[n][ii]),
		  (m->alpha_st_old[n][ii]), (m->alpha_st_e[n][ii]),
		  (m->alpha_st_e_old[n][ii]), (m->recomb_sp[n][ii]),
		  (m->recomb_sp_e[n][ii]));
	}



    }
  else
    {
      printf ("This is not a macro level\n");
    }
  return (p->macro_info);
}

/*********************************************************************************************/

/*

   080812	ksl	62 - Updated to reflect new versions of programs to do 
   			partition functions and calculate saha densities
*/

int
depcoef_overview (icell)
     int icell;
{
  ConfigPtr p;
  PlasmaPtr x, xdummy;
//OLD  double xden, lteden, nh;
//OLD  int n, nelem;
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
//OLD      do_partitions(xdummy, 1);
//OLD      levels(xdummy, 1);
//OLD      nh = x->rho * rho2nh;   
//OLD      for (nelem = 0; nelem < nelements; nelem++)
//OLD   {
//OLD     saha(xdummy->ne, nh, xdummy->t_e, nelem,  xdummy->density);
//OLD   }
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

/*********************************************************************************************/

int
copy_plasma (x1, x2)
     PlasmaPtr x1, x2;
{
  /*Routine to copy the necessary parts of a plasma structure for computing a set of level populations.
     x1 points to the cell from which data is copied and x2 points to the cell to which data is copied */

  x2->nwind = x1->nwind;
  x2->nplasma = x1->nplasma;
  x2->ne = x1->ne;
  x2->rho = x1->rho;
  x2->vol = x1->vol;
  x2->t_r = x1->t_r;
  x2->t_e = x1->t_e;
  x2->w = x1->w;

  /*Note this isn't everything in the cell! Only the things needed to get a set of
     LTE populations! */

  return (0);

}

/*********************************************************************************************/


/*

   080812	ksl	62 - Updated to reflect new versions of programs to do 
   			partition functions and calculate saha densities
*/

int
depcoef_overview_specific (version, nconfig, w, rootname, ochoice)
     int version;
     int nconfig;
     WindPtr w;
     char rootname[];
     int ochoice;
{
  int n;
  char filename[LINELENGTH];
//OLD  int nelem;
  ConfigPtr p;
  PlasmaPtr x, xdummy;
//OLD  double xden, lteden, nh;
  double xden, lteden;
  int copy_plasma ();
  plasma_dummy pdum;


  for (n = 0; n < NDIM2; n++)
    {
      aaa[n] = 0;
      if (w[n].vol > 0.0)
	{
	  xden = -1;

	  x = &plasmamain[w[n].nplasma];
	  xdummy = &pdum;

	  copy_plasma (x, xdummy);
	  geo.macro_ioniz_mode = 0;

	  partition_functions (xdummy, 1);
	  saha (xdummy, xdummy->ne, xdummy->t_e);

//OLD     do_partitions(xdummy, 1);
//OLD     levels(xdummy, 1);
//OLD     nh = x->rho * rho2nh;   
//OLD     for (nelem = 0; nelem < nelements; nelem++)
//OLD       {
//OLD         saha(xdummy->ne, nh, xdummy->t_e, nelem,  xdummy->density);
//OLD       }

	  geo.macro_ioniz_mode = 1;



	  p = &config[nconfig];
	  xden = den_config (x, nconfig);
	  lteden = den_config (xdummy, nconfig);

	  if (version == 0)
	    {
	      aaa[n] = xden;
	    }
	  else
	    {
	      aaa[n] = xden / lteden;
	    }
	}
    }
  display ("Dep_coef");

  if (ochoice)
    {
      strcpy (filename, rootname);
      strcat (filename, ".dep_coef");
      write_array (filename, ochoice);

    }

  return (0);

}
