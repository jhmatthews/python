

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "atomic.h"

#include "python.h"

/***********************************************************
                                       Space Telescope Science Institute

 Synopsis:
	This file contains routines for handling the interaction
	between the GridPtr and WindPtr array

 Arguments:

	
 Returns:

Description:

Notes:


History:
	1112	ksl	Began to make changes so that one could 
			read a new wind file from within one 
			of the python routines

**************************************************************/



/***********************************************************
                                       Space Telescope Science Institute

 Synopsis: create_maps (ichoice)

 Arguments:
     int ichoice	!0 (or true) then the size of the plasma structure is
			just given by the cells with non-zero volume
			0 (or false) then the size of the palma structure
			is the same as the wind structue

	
 Returns:

Description:

Notes:


History:
	06may	ksl	Created as part of the attempt to reduce the
			overall size of the structures

**************************************************************/


int
create_maps (ichoice)
     int ichoice;
{
  int i, j;
  j = 0;
  if (ichoice)
    {
      //map to reduce space
      for (i = 0; i < NDIM2; i++)
	{
	  wmain[i].nwind = i;
	  if (wmain[i].vol > 0)
	    {
	      wmain[i].nplasma = j;
	      plasmamain[j].nplasma = j;
	      plasmamain[j].nwind = i;
	      j++;
	    }
	  else
	    {
	      wmain[i].nplasma = NPLASMA;
	    }
	}
      if (j != NPLASMA)
	{
	  Error
	    ("create_maps: Problems with matching cells -- Expected %d Got %d\n",
	     NPLASMA, j);
	  exit (0);
	}
    }
  else
    //one plama cell for each wind cell
    {
      for (i = 0; i < NDIM2; i++)
	{
	  wmain[i].nwind = i;
	  wmain[i].nplasma = i;
	  plasmamain[i].nplasma = i;
	  plasmamain[i].nwind = i;

	}

    }

  plasmamain[NPLASMA].nplasma = NPLASMA;
  plasmamain[NPLASMA].nwind = -1;
  return (0);
}



/***********************************************************
                                       Space Telescope Science Institute

 Synopsis: calloc_wind (nelem)

 Arguments:

	
 Returns:

 Description:

 Notes:


 History:
	06may	ksl	57a -- Coded
	11dec	ksl	71 - Modified so that the memory would be
			reallocated if necessary

**************************************************************/


int
calloc_wind (nelem)
     int nelem;
{

  if (wmain != NULL)
    {
      free (wmain);
    }

  wmain = (WindPtr) calloc (sizeof (wind_dummy), nelem + 1);

  if (wmain == NULL)
    {
      Error
	("There is a problem in allocating memory for the wind structure\n");
      exit (0);
    }
  else
    {
      Log_silent
	("Allocated %10d bytes for each of %5d elements of             totaling %10.1f Mb\n",
	 sizeof (wind_dummy), nelem, 1.e-6 * nelem * sizeof (wind_dummy));
    }

  return (0);
}



/***********************************************************
                                       Space Telescope Science Institute

 Synopsis: calloc_plasma (nelem)

 Arguments:

	
 Returns:

Description:

Notes:

	This only allocates elements.  It does not populate them
	with any information.


History:
11dec	ksl	71 - Modified so that the memory would be
		reallocated if necessary

**************************************************************/


int
calloc_plasma (nelem)
     int nelem;
{

  if (plasmamain != NULL)
    {
      free (plasmamain);
    }

  //Allocate one extra element to store data where there is no volume

  plasmamain = (PlasmaPtr) calloc (sizeof (plasma_dummy), (nelem + 1));
  geo.nplasma = nelem;

  if (plasmamain == NULL)
    {
      Error
	("There is a problem in allocating memory for the plasma structure\n");
      exit (0);
    }
  else
    {
      Log_silent
	("Allocated %10d bytes for each of %5d elements of      plasma totaling %10.1f Mb \n",
	 sizeof (plasma_dummy), (nelem + 1),
	 1.e-6 * (nelem + 1) * sizeof (plasma_dummy));
    }

  /* Now allocate space for storing photon frequencies -- 57h */
  if (photstoremain != NULL)
    {
      free (photstoremain);
    }
  photstoremain =
    (PhotStorePtr) calloc (sizeof (photon_store_dummy), (nelem + 1));

  if (photstoremain == NULL)
    {
      Error
	("There is a problem in allocating memory for the photonstore structure\n");
      exit (0);
    }
  else
    {
      Log_silent
	("Allocated %10d bytes for each of %5d elements of photonstore totaling %10.1f Mb \n",
	 sizeof (photon_store_dummy), (nelem + 1),
	 1.e-6 * (nelem + 1) * sizeof (photon_store_dummy));
    }

  return (0);
}


int
check_plasma (xplasma, message)
     PlasmaPtr xplasma;
     char message[];
{
  if (xplasma->nplasma == NPLASMA)
    {
      Error ("check_plasma -- %s \n", message);
      return (1);
      //true
    }
  else
    return (0);
}



/***********************************************************
               Space Telescope Science Institute

 Synopsis: calloc_macro (nelem)

 Arguments:

	
 Returns:

Description:

Notes:

060805 -- ksl -- Begin notes on changing the macro structures

At present, we allocate macromain if there are any macrolevels,
but the size of macromain is determined by NLEVELS_MACRO and
NBBJUMPS or NBFJUMPS.  These are not set dynamically, and so
one should consider how to do this, since for the current values
these are set to 500 for NLEVELS_MACRO, 20 for NBBJUMPS, and
NBFJUMPS, which means they are roughly 1000 per level.  The
first thing to do about this is to reduce them to values that
are more like those in our datafiles.

If further work is needed there are several possibilties,
each of which has some disadvantages relative to the current
appraoch.  The advantage of the current approch is that each
plasma cell has a separate macro structure element and that
these are addressed rather logically, e.g
	macromain[n].jbar[macro_level][nbf]
It is pretty easy to see what is going on in arrays.  If one
wants to dynamically allocate say the number of macro_levels,
then one will create a new structure called for example
macro2 to contain jbar, jbar_old etc, and then to map to the
appropriate elements of it  from the macro_array.  There are
several ways to do this.

One would be to follow the approach used in other portions of
the atomic structures to add t add the place to start in the
macro2 structure as an integer.  So for example, if there
were 35 macrolevels read in, one could simply add a number
that told one where to start in the macro2 array, much as
is done in looking for firstion and lastion for an element.

In that case, the addressing would look somelike
	macro2[macromain[n].start_level]].jbar[nbf]

This is not crazy since we know how many levels have
been read in, and so the start_levels will increase
by the number of levels.  An advantage of this approach
is that you can dump the structues to disk, and read
them back with another program and eveything will be
properly indexed.

Another appraoch that would work would be to create
pointers instead of indexes.  This may be faster. In
this case the addresing would look something like.
	macro[n].lev_ptr[nlevel]->jbar[nbf]
where lev_ptr was defined as *Macro2Ptr within the
macro structure.

A final approach would be to dynamically resize the
way the macro array itselve is addressed within each
routine. This can also be done.  One simply creates a that
parallels the macro structue within each routine, e.g

struct macro xmacro[][nlevel]

and addrsses things as

        xmacro[n][nlev].jbar[nbf]

This also looks fairly straightforward, in c, and reads pretty
esaily.


At present, I have not done anything about this

End notes on resizing -- 060805



History:
	06jun	ksl	Coded
	0608	ksl	57h -- Recoded to eliminate the creation
			of the array when there are no macro_atoms.
			Previously the array had been created but
			not written out.
	1112	ksl	71 - Added checks to see if macromain, had
			been allocated previously and if so to
			reallocate

**************************************************************/


int
calloc_macro (nelem)
     int nelem;
{


  if (nlevels_macro == 0 && geo.nmacro == 0)
    {
      
      Log
	("calloc_macro: Allocated no space for macro since nlevels_macro==0 and geo.nmacro==0\n");
      return (0);
    }
  if (macromain != NULL)
    {
      free (macromain);
    }

  //Allocate one extra element to store data where there is no volume

  macromain = (MacroPtr) calloc (sizeof (macro_dummy), (nelem + 1));
  geo.nmacro = nelem;

  if (macromain == NULL)
    {
      Error
	("calloc_macro: There is a problem in allocating memory for the macro structure\n");
      exit (0);
    }
  else if (nlevels_macro > 0 || geo.nmacro > 0)
    {
      Log
	("Allocated %10d bytes for each of %5d elements of       macro totaling %10.1f Mb \n",
	 sizeof (macro_dummy), (nelem + 1),
	 1.e-6 * (nelem + 1) * sizeof (macro_dummy));
    }
  else
    {
      Log
	("calloc_macro: Allocated no space for macro since nlevels_macro==0\n");
    }

  return (0);
}


/**************************************************************i

This sections seems to have been added by Stuart in the summer of 2009
but it is not documented.  It was an attempt to reduce the size
of the windsave file when macro atoms are used.


Note that nelem here refers to an element of the macro array, not
the number of elements


1112	ksl	71 - Added code that is intended to allow one to realloate the memory
		if necessary, but the way this is construced makes it easy to
		cause errors and it is not obvious how to check this until
		we put a macro model back in
130625  JM      Commented out free statements due to PYWIND MALLOC MATOM BUG
 */


int
calloc_estimators (nelem)
     int nelem;
{
  int n;

  if (nlevels_macro == 0 && geo.nmacro == 0)
    {
      geo.nmacro = 0;
      Log_silent
	("Allocated no space for MA estimators since nlevels_macro==0 and geo.nmacro==0\n");
      return (0);
    }
  //Allocate one extra element to store data where there is no volume


  size_Jbar_est = 0;
  size_gamma_est = 0;
  size_alpha_est = 0;

  for (n = 0; n < nlevels_macro; n++)
    {
      Log
	("calloc_estimators: level %d has n_bbu_jump %d  n_bbd_jump %d n_bfu_jump %d n_bfd_jump %d\n",
	 n, config[n].n_bbu_jump, config[n].n_bbd_jump, config[n].n_bfu_jump,
	 config[n].n_bfd_jump);
      config[n].bbu_indx_first = size_Jbar_est;
      size_Jbar_est += config[n].n_bbu_jump;
      config[n].bfu_indx_first = size_gamma_est;
      size_gamma_est += config[n].n_bfu_jump;
      config[n].bfd_indx_first = size_alpha_est;
      size_alpha_est += config[n].n_bfd_jump;
    }




  Log
    ("calloc_estimators: size_Jbar_est %d size_gamma_est %d size_alpha_est %d\n",
     size_Jbar_est, size_gamma_est, size_alpha_est);


  for (n = 0; n < nelem; n++)
    {
      /* JM130625: Commented out free statements due to PYWIND MALLOC MATOM BUG
	 if (macromain[n].jbar != NULL)
	{
	  free (macromain[n].jbar);
	}*/
      if ((macromain[n].jbar =
	   calloc (sizeof (double), size_Jbar_est)) == NULL)
	{
	  Error
	    ("calloc_estimators: Error in allocating memory for MA estimators\n");
	  exit (0);
	}

      /*if (macromain[n].jbar_old != NULL)
	{
	  free (macromain[n].jbar_old);
	}*/
      if ((macromain[n].jbar_old =
	   calloc (sizeof (double), size_Jbar_est)) == NULL)
	{
	  Error
	    ("calloc_estimators: Error in allocating memory for MA estimators\n");
	  exit (0);
	}

      /*if (macromain[n].gamma != NULL)
	{
	  free (macromain[n].gamma);
	}*/
      if ((macromain[n].gamma =
	   calloc (sizeof (double), size_gamma_est)) == NULL)
	{
	  Error
	    ("calloc_estimators: Error in allocating memory for MA estimators\n");
	  exit (0);
	}

      /*if (macromain[n].gamma_old != NULL)
	{
	  free (macromain[n].gamma_old);
	}*/
      if ((macromain[n].gamma_old =
	   calloc (sizeof (double), size_gamma_est)) == NULL)
	{
	  Error
	    ("calloc_estimators: Error in allocating memory for MA estimators\n");
	  exit (0);
	}

      /*if (macromain[n].gamma_e != NULL)
	{
	  free (macromain[n].gamma_e);
	}*/
      if ((macromain[n].gamma_e =
	   calloc (sizeof (double), size_gamma_est)) == NULL)
	{
	  Error
	    ("calloc_estimators: Error in allocating memory for MA estimators\n");
	  exit (0);
	}

      /*if (macromain[n].gamma_e_old != NULL)
	{
	  free (macromain[n].gamma_e_old);
	}*/
      if ((macromain[n].gamma_e_old =
	   calloc (sizeof (double), size_gamma_est)) == NULL)
	{
	  Error
	    ("calloc_estimators: Error in allocating memory for MA estimators\n");
	  exit (0);
	}

      /*if (macromain[n].alpha_st != NULL)
	{
	  free (macromain[n].alpha_st);
	}*/
      if ((macromain[n].alpha_st =
	   calloc (sizeof (double), size_gamma_est)) == NULL)
	{
	  Error
	    ("calloc_estimators: Error in allocating memory for MA estimators\n");
	  exit (0);
	}

      /*if (macromain[n].alpha_st_old != NULL)
	{
	  free (macromain[n].alpha_st_old);
	}*/
      if ((macromain[n].alpha_st_old =
	   calloc (sizeof (double), size_gamma_est)) == NULL)
	{
	  Error
	    ("calloc_estimators: Error in allocating memory for MA estimators\n");
	  exit (0);
	}

      /*if (macromain[n].alpha_st_e != NULL)
	{
	  free (macromain[n].alpha_st_e);
	}*/
      if ((macromain[n].alpha_st_e =
	   calloc (sizeof (double), size_gamma_est)) == NULL)
	{
	  Error
	    ("calloc_estimators: Error in allocating memory for MA estimators\n");
	  exit (0);
	}

      /*if (macromain[n].alpha_st_e_old)
	{
	  free (macromain[n].alpha_st_e_old);
	}*/
      if ((macromain[n].alpha_st_e_old =
	   calloc (sizeof (double), size_gamma_est)) == NULL)
	{
	  Error
	    ("calloc_estimators: Error in allocating memory for MA estimators\n");
	  exit (0);
	}

      /*if (macromain[n].recomb_sp != NULL)
	{
	  free (macromain[n].recomb_sp);
	}*/
      if ((macromain[n].recomb_sp =
	   calloc (sizeof (double), size_alpha_est)) == NULL)
	{
	  Error
	    ("calloc_estimators: Error in allocating memory for MA estimators\n");
	  exit (0);
	}

      /*if (macromain[n].recomb_sp_e != NULL)
	{
	  free (macromain[n].recomb_sp_e);
	}*/
      if ((macromain[n].recomb_sp_e =
	   calloc (sizeof (double), size_alpha_est)) == NULL)
	{
	  Error
	    ("calloc_estimators: Error in allocating memory for MA estimators\n");
	  exit (0);
	}

      /*if (macromain[n].matom_emiss != NULL)
	{
	  free (macromain[n].matom_emiss);
	}*/
      if ((macromain[n].matom_emiss =
	   calloc (sizeof (double), nlevels_macro)) == NULL)
	{
	  Error
	    ("calloc_estimators: Error in allocating memory for MA estimators\n");
	  exit (0);
	}

      /*if (macromain[n].matom_abs != NULL)
	{
	  free (macromain[n].matom_abs);
	}*/
      if ((macromain[n].matom_abs =
	   calloc (sizeof (double), nlevels_macro)) == NULL)
	{
	  Error
	    ("calloc_estimators: Error in allocating memory for MA estimators\n");
	  exit (0);
	}

      /* Added ksl 091103 59e */
      /*if (macromain[n].cooling_bf != NULL)
	{
	  free (macromain[n].cooling_bf);
	}*/
      if ((macromain[n].cooling_bf =
	   calloc (sizeof (double), ntop_phot)) == NULL)
	{
	  Error
	    ("calloc_estimators: Error in allocating memory for MA estimators\n");
	  exit (0);
	}

      /*if (macromain[n].cooling_bf_col != NULL)
	{
	  free (macromain[n].cooling_bf_col);
	}*/
      if ((macromain[n].cooling_bf_col =
	   calloc (sizeof (double), ntop_phot)) == NULL)
	{
	  Error
	    ("calloc_estimators: Error in allocating memory for MA estimators\n");
	  exit (0);
	}

      /*if (macromain[n].cooling_bb != NULL)
	{
	  free (macromain[n].cooling_bb);
	}*/
      if ((macromain[n].cooling_bb =
	   calloc (sizeof (double), nlines)) == NULL)
	{
	  Error
	    ("calloc_estimators: Error in allocating memory for MA estimators\n");
	  exit (0);
	}
    }



  if (nlevels_macro > 0 || geo.nmacro > 0)
    {
      Log_silent
	("Allocated %10.1f Mb for MA estimators \n",
	 1.e-6 * (nelem + 1) * (2. * nlevels_macro + 2. * size_alpha_est +
				8. * size_gamma_est +
				2. * size_Jbar_est) * sizeof (double));
    }
  else
    {
      Log_silent ("Allocated no space for macro since nlevels_macro==0\n");
    }

  return (0);
}

/***********************************************************
                                       West Lulworth

 Synopsis:
	This subroutine allocates space for dynamic arrays in the plasma 
	structure

 Arguments:
	nelem - the number of plasma cells (actually NPLASMA+1) to allow for empty cell

	
 Returns:


Description:
	This subroutine allocates space for variable length arrays in the plasma structure.
	This is to save space over the previous versions of python, which allocated arrays
	based on guesses as to the lengths. 

Notes:
	Arrays sized to the number of ions are largest,
	and dominate the size of nplasma, so these were first to be dynamically allocated.

History:
	1407	nsh	Started out allocating arrays that have length nion

**************************************************************/


int
calloc_dyn_plasma (nelem)
     int nelem;
{
  int n;

for (n=0;n<nelem+1;n++)  //We loop over all elements in the plasma array, adding one for an empty cell used for extrapolations.
	{
  	if ((plasmamain[n].density =calloc (sizeof (double), nions)) == NULL)
	{
	  Error
	    ("calloc_dyn_plasma: Error in allocating memory for density\n");
	  exit (0);
	}
	if ((plasmamain[n].partition =calloc (sizeof (double), nions)) == NULL)
	{
	  Error
	    ("calloc_dyn_plasma: Error in allocating memory for partition\n");
	  exit (0);
	}
	if ((plasmamain[n].PWdenom =calloc (sizeof (double), nions)) == NULL)
	{
	  Error
	    ("calloc_dyn_plasma: Error in allocating memory for PWdenom\n");
	  exit (0);
	}
	if ((plasmamain[n].PWdtemp =calloc (sizeof (double), nions)) == NULL)
	{
	  Error
	    ("calloc_dyn_plasma: Error in allocating memory for PWtemp\n");
	  exit (0);
	}
	if ((plasmamain[n].PWnumer =calloc (sizeof (double), nions)) == NULL)
	{
	  Error
	    ("calloc_dyn_plasma: Error in allocating memory for PWnumer\n");
	  exit (0);
	}
	if ((plasmamain[n].PWntemp =calloc (sizeof (double), nions)) == NULL)
	{
	  Error
	    ("calloc_dyn_plasma: Error in allocating memory for PWntemp\n");
	  exit (0);
	}
	if ((plasmamain[n].ioniz =calloc (sizeof (double), nions)) == NULL)
	{
	  Error
	    ("calloc_dyn_plasma: Error in allocating memory for ioniz\n");
	  exit (0);
	}
	if ((plasmamain[n].recomb =calloc (sizeof (double), nions)) == NULL)
	{
	  Error
	    ("calloc_dyn_plasma: Error in allocating memory for recomb\n");
	  exit (0);
	}
	if ((plasmamain[n].scatters =calloc (sizeof (int), nions)) == NULL)
	{
	  Error
	    ("calloc_dyn_plasma: Error in allocating memory for scatters\n");
	  exit (0);
	}
	if ((plasmamain[n].xscatters =calloc (sizeof (double), nions)) == NULL)
	{
	  Error
	    ("calloc_dyn_plasma: Error in allocating memory for xscatters\n");
	  exit (0);
	}
	if ((plasmamain[n].heat_ion =calloc (sizeof (double), nions)) == NULL)
	{
	  Error
	    ("calloc_dyn_plasma: Error in allocating memory for heat_ion\n");
	  exit (0);
	}
	if ((plasmamain[n].lum_ion =calloc (sizeof (double), nions)) == NULL)
	{
	  Error
	    ("calloc_dyn_plasma: Error in allocating memory for lum_ion\n");
	  exit (0);	
	}
}

     Log
	("Allocated %10d bytes for each of %5d elements variable length plasma arrays totaling %10.1f Mb \n",
	 sizeof (double)*nions*12, (nelem+1),
	 1.e-6 * (nelem+1 ) * sizeof (double)*nions*12);

  return (0);
}
