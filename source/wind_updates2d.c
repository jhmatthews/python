/* The routines in this file are generic.  There is no dependence on a particlar wind model or
 * any coordinate system dependences.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>


#include "atomic.h"
#include "python.h"
#define LINELEN 132

/***********************************************************
                                       Space Telescope Science Institute

 Synopsis: wind_update(w) updates the parameters in the wind that are 
	affected by radiation.

 Arguments:		
	WindPtr w;

Returns:
 
Description:	
	
Notes:
	There is not much point in calling this until you have propagated a few photons
	The routine should not be called until at least some photons have been propagated 
	through the wind because the radiation and electron temperatures are calculated from 
	the frequency moment
History:
 	97jan      ksl	Coding on python began.
 	97sep13	ksl	Added the possibility of fixing the concentrations of the ions to specific
 			values, but this really should not be called because in that case there
 			is little point in updating the wind since the only thing which wind_update
 			really does is to change the ionization fractions.  This change was really
 			for consistency.
 	97oct2	ksl	Modified relation between t_e and t_r to use CK's prescription derived from
 			Cloudy comparisons.
 	98may3	ksl	Fixed problems with the portion of the routine which placed new densities in
 			the regions just outside the wind.  These densities are used for interpolation
 			purposes.
 	98may16	ksl	Added some checks to indicate how great the changes in t_e and t_r in each
 			wind_update
 	98may17	ksl	Changed some of the logic to (hopefully) eliminate some calculations outside
 			the wind which should not be necessary.   Added some additional tracking of 
 			changes to see how things converge, and where the biggest changes in t_e 
 			and t_r are.
 	98may24	ksl	Changed code so that w[].vol > 0 is the primary criterion for attempting to
 			recalculate the ionization abundances
 	98may31 ksl	Corrected error in calculating dilution factor because cells are actually two 
 			cells one above and one below the plane of the disk.
 	98dec	ksl	Relocated section of code which chose between various ways of calculating
 			ion abundances 	to a subroutine ion_abundances.  Removed section of code
 			in which contained estimates of what t_e ought to be based on t_r, w, etc.
	99jan	ksl	Modified how t_e and t_r change in the event there were no photons in
			a cell.  Old solution was to set t_e and t_r to zero; new solution was
			to leave t_r at old value and drop t_e by the maximum amount 0.7*t_e_old
			being used in one_shot.  
        04May   SS      Removed extra factor of two in calculation of the radiation density since 
			the cell
                        volume now includes ares both above and below the plane.
	05apr	ksl	55d: Began work to make more compatible with 1-d coordinate systems.
	05apr	ksl	56: Moved functionality that extends the ion densityies
			to cylindrical.c, rtheta.c, and spherical.c in 
			continuing effort to isolate difference between
			coordinate system types.
	06may	ksl	57+: Began splitting of wind and plasma structure.  This routine loops
			over the entire grid and ultimately eveything here is likely to be
			contained in the plasma structure.  Therefore here we want to use
			plasmamain
	11aug	nsh	70: modifications made to wind_update  and wind_rad_init to incorporate 
			compton heating and cooling.
        12apr	nh	72: modifications makde to wind_update and wind_rad_init to incoprorate
			induced copmton heating.
        130621  jm      76: added lines to store luminosities from ionization cycles in geo and plasma
			structure to solve windsave bug
   	13jul	nsh	76: added lines to deal with the case when adiacaitc cooling becomes heating in
			complex wind geometries
	13nov	nsh	77: Some changes in the parallel communications for new variables.
	14jul	nsh	78a: Changed the length of the communications array to take account of
			dynamically allocated arrays in plasma structure. Also changed some of the pack
			and unpack commands for the same reason.
	14sept	nsh	78b: Changes to deal with the inclusion of direct recombination
	14nov 	JM 78b: Changed volume to be the filled volume


**************************************************************/


int num_updates = 0;

int
wind_update (w)
WindPtr (w);
{
  int n, i, j;
  double trad, nh;
  double wtest, xsum, asum, psum, fsum, lsum, csum, icsum, ausum;	/*1108 NSH csum added to sum compton heating 1204 NSH icsum added to sum induced compton heating */
  double volume;
  double vol;
  char string[LINELEN];
  double t_r_old, t_e_old, dt_r, dt_e;
  double t_r_ave_old, t_r_ave, t_e_ave_old, t_e_ave;
  int iave, nmax_r, nmax_e;
  int nplasma;
  int nwind;
  int first, last, m;
  double tot, agn_ip;
  double nsh_lum_hhe;
  double nsh_lum_metals;
  int my_nmin, my_nmax;	//Note that these variables are still used even without MPI on
  FILE *fptr, *fopen (); /*This is the file to communicate with zeus */
  
#ifdef MPI_ON
  int num_mpi_cells, num_mpi_extra, position, ndo, n_mpi, num_comm, n_mpi2;
  int size_of_commbuffer;
  char *commbuffer;

  /* JM 1409 -- Added for issue #110 to ensure correct reporting in parallel */
  int nmax_r_temp, nmax_e_temp;
  double dt_e_temp, dt_r_temp;


  /* the commbuffer needs to be larger enough to pack all variables in MPI_Pack and MPI_Unpack routines NSH 1407 - the 
  NIONS changed to nions for the 12 arrays in plasma that are now dynamically allocated */
  size_of_commbuffer = 8 * (12*nions + NLTE_LEVELS + 2*NTOP_PHOT + 12*NXBANDS + 2*LPDF + NAUGER + 106)*(floor(NPLASMA/np_mpi_global)+1);
  commbuffer = (char *) malloc(size_of_commbuffer*sizeof(char));

  /* JM 1409 -- Initialise parallel only variables */
  nmax_r_temp = nmax_e_temp = -1;
  dt_e_temp = dt_r_temp = 0.0;
      
#endif
  dt_r = dt_e = 0.0;
  iave = 0;
  nmax_r = nmax_e = -1;
  t_r_ave_old = t_r_ave = t_e_ave_old = t_e_ave = 0.0;


  /* For MPI parallelisation, the following loop will be distributed over mutiple tasks. 
     Note that the mynmim and mynmax variables are still used even without MPI on */
  my_nmin = 0;
  my_nmax = NPLASMA;
#ifdef MPI_ON
  num_mpi_cells = floor(NPLASMA/np_mpi_global);
  num_mpi_extra = NPLASMA - (np_mpi_global*num_mpi_cells);

  /* this section distributes the remainder over the threads if the cells
     do not divide evenly by thread */
  if (rank_global < num_mpi_extra)
    {
      my_nmin = rank_global*(num_mpi_cells+1);
      my_nmax = (rank_global+1)*(num_mpi_cells+1);     
    }
  else
    {
      my_nmin = num_mpi_extra*(num_mpi_cells+1) + (rank_global-num_mpi_extra)*(num_mpi_cells);
      my_nmax = num_mpi_extra*(num_mpi_cells+1) + (rank_global-num_mpi_extra+1)*(num_mpi_cells);
    }
  ndo = my_nmax-my_nmin;
#endif

  /* Before we do anything let's record the average tr and te from the last cycel */
  /* JM 1409 -- Added for issue #110 to ensure correct reporting in parallel */
  for (n=0; n<NPLASMA; n++)
  {
  	t_r_ave_old += plasmamain[n].t_r;
    t_e_ave_old += plasmamain[n].t_e;
  }

  /* we now know how many cells this thread has to process - note this will be 
     0-NPLASMA in serial mode */

  for (n = my_nmin; n < my_nmax; n++)
    {


      nwind = plasmamain[n].nwind;
      volume = w[nwind].vol;

      //volume = plasmamain[n].vol; 	// 1411 JM -- we want volume to be the filled volume
      //1504 -- JM No we don't! The mean intensity shouldn't change if we *just* clump


      /* Start with a call to the routine which normalises all the macro atom 
         monte carlo radiation field estimators. It's best to do this first since
         some of the estimators include temperature terms (stimulated correction
         terms) which were included during the monte carlo simulation so we want 
         to be sure that the SAME temperatures are used here. (SS - Mar 2004). */

      if (geo.rt_mode == 2 && geo.macro_simple == 0)	//test for macro atoms
	{		
	  mc_estimator_normalise (nwind);  
	  macromain[n].kpkt_rates_known = -1;
	}

      /* Store some information so one can determine how much the temps are changing */
      t_r_old = plasmamain[n].t_r;
      t_e_old = plasmamain[n].t_e;
      iave++;


      if (plasmamain[n].ntot < 100)
	{
	  Log
	    ("!!wind_update: Cell %d  with volume %8.2e has only %d photons\n",
	     n, volume, plasmamain[n].ntot);
	}

      if (plasmamain[n].ntot > 0)
	{
	  wtest = plasmamain[n].ave_freq;
	  plasmamain[n].ave_freq /= plasmamain[n].j;	/* Normalization to frequency moment */
	  if (sane_check (plasmamain[n].ave_freq))
	    {
	      Error ("wind_update:sane_check %d ave_freq %e j %e ntot %d\n",
		     n, wtest, plasmamain[n].j, plasmamain[n].ntot);
	    }

	  plasmamain[n].j /= (4. * PI * volume);	//Factor of 2 has been removed from this line (SS, May04)
	  plasmamain[n].j_direct /= (4. * PI * volume);
	  plasmamain[n].j_scatt /= (4. * PI * volume);



	  trad = plasmamain[n].t_r =
	    H * plasmamain[n].ave_freq / (BOLTZMANN * 3.832);
	  plasmamain[n].w =
	    PI * plasmamain[n].j / (STEFAN_BOLTZMANN * trad * trad * trad *
				    trad);


	  if (plasmamain[n].w > 1e10)
	    {
	      Error
		("wind_update: Huge w %8.2e in cell %d trad %10.2e j %8.2e\n",
		 plasmamain[n].w, n, trad, plasmamain[n].j);
	    }
	  if (sane_check (trad) || sane_check (plasmamain[n].w))
	    {
	      Error ("wind_update:sane_check %d trad %8.2e w %8.2g\n", n,
		     trad, plasmamain[n].w);
	      Error ("wind_update: ave_freq %8.2e j %8.2e\n",
		     plasmamain[n].ave_freq, plasmamain[n].j);
	      exit (0);
	    }
	}
      else
	{			// It is not clear what to do with no photons in a cell

	  plasmamain[n].j = plasmamain[n].j_direct = plasmamain[n].j_scatt = 0;
	  trad = plasmamain[n].t_r;
	  plasmamain[n].t_e *= 0.7;
	  if (plasmamain[n].t_e < TMIN)
	    plasmamain[n].t_e = TMIN;
	  plasmamain[n].w = 0;
	}


      /* 1108 NSH/KSL  This loop is to calculate the frequency banded j and ave_freq variables */

      /* 71 - 111279 - ksl - Small modification to reflect the fact that nxfreq has been moved into the geo structure */
      for (i = 0; i < geo.nxfreq; i++)	/*loop over number of bands */
	{
	  if (plasmamain[n].nxtot[i] > 0)	/*Check we actually have some photons in the cell in this band */
	    {

	      plasmamain[n].xave_freq[i] /= plasmamain[n].xj[i];	/*Normalise the average frequency */
	      plasmamain[n].xsd_freq[i] /= plasmamain[n].xj[i];	/*Normalise the mean square frequency */
	      plasmamain[n].xsd_freq[i] = sqrt (plasmamain[n].xsd_freq[i] - (plasmamain[n].xave_freq[i] * plasmamain[n].xave_freq[i]));	/*Compute standard deviation */
	      plasmamain[n].xj[i] /= (4 * PI * volume);	/*Convert to radiation density */

	    }
	  else
	    {
	      plasmamain[n].xj[i] = 0;	/*If no photons, set both radiation estimators to zero */
	      plasmamain[n].xave_freq[i] = 0;
	      plasmamain[n].xsd_freq[i] = 0;	/*NSH 120815 and also the SD ???? */
	    }
	}

/* 1108 NSH End of loop */



      nh = plasmamain[n].rho * rho2nh;

/* 1110 NSH Normalise IP, which at this point should be the number of photons in a cell by dividing by volume and number density of hydrogen in the cell */

      plasmamain[n].ip /= (C * volume * nh);
      plasmamain[n].ip_direct /= (C * volume * nh);
      plasmamain[n].ip_scatt /= (C * volume * nh);

/* 1510 NSH Normalise xi, which at this point should be the luminosity of ionizing photons in a cell (just the sum of photon weights) */

		plasmamain[n].xi *= 4.*PI;
		plasmamain[n].xi /= ( volume * nh);

      /* If geo.adiabatic is true, then alculate the adiabatic cooling using the current, i.e 
       * previous value of t_e.  Note that this may not be  best way to determien the cooling. 
       * Changes made here should also be reflected in wind2d.c.  At present, adiabatic colling
       * is not included in updates to the temperature, even if the adiabatic cooling is calculated
       * here. 04nov -- ksl 
       * 05apr -- ksl -- The index being used was incorrect.  This has been fixed now 
       * 11sep -- nsh -- The index for the wind (&w) for adiabatic cooling was incorrect - 
       * was being called with the plasma cell rather than the approriate wind cell - fixed 
       * old: adiabatic_cooling (&w[n], plasmamain[n].t_e);
       */

      if (geo.adiabatic)
	plasmamain[n].lum_adiabatic = adiabatic_cooling (&w[nwind], plasmamain[n].t_e);
      else
	plasmamain[n].lum_adiabatic = 0.0;


      /* Calculate the densities in various ways depending on the ioniz_mode */

      ion_abundances (&plasmamain[n], geo.ioniz_mode);



      /* Perform checks to see how much temperatures have changed in this iteration */

      if ((fabs (t_r_old - plasmamain[n].t_r)) > fabs (dt_r))
	{
	  dt_r = plasmamain[n].t_r - t_r_old;
	  nmax_r = n;
	}
      if ((fabs (t_e_old - plasmamain[n].t_e)) > fabs (dt_e))
	{
	  dt_e = plasmamain[n].t_e - t_e_old;
	  nmax_e = n;
	}
      t_r_ave += plasmamain[n].t_r;
      t_e_ave += plasmamain[n].t_e;
    }



  /*This is the end of the update loop that is parallised. We now need to exchange data between the tasks. */
#ifdef MPI_ON
  for (n_mpi=0; n_mpi < np_mpi_global; n_mpi++)
    {
      position = 0;
      
      if (rank_global == n_mpi)
	{
	  Log("MPI task %d is working on cells %d to max %d (total size %d).\n", rank_global, my_nmin, my_nmax, NPLASMA);
	  MPI_Pack(&ndo, 1, MPI_INT, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
	  for (n = my_nmin; n < my_nmax; n++)
	    {
	      MPI_Pack(&n, 1, MPI_INT, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
	      MPI_Pack(&plasmamain[n].nwind, 1, MPI_INT, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
	      MPI_Pack(&plasmamain[n].nplasma, 1, MPI_INT, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
	      MPI_Pack(&plasmamain[n].ne, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
	      MPI_Pack(&plasmamain[n].rho, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
	      MPI_Pack(&plasmamain[n].vol, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
	      MPI_Pack(plasmamain[n].density, nions, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
	      MPI_Pack(plasmamain[n].partition, nions, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
	      MPI_Pack(plasmamain[n].levden, NLTE_LEVELS, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
	      MPI_Pack(plasmamain[n].PWdenom, nions, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
	      MPI_Pack(plasmamain[n].PWdtemp, nions, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
	      MPI_Pack(plasmamain[n].PWnumer, nions, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
	      MPI_Pack(plasmamain[n].PWntemp, nions, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
	      MPI_Pack(&plasmamain[n].kappa_ff_factor, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
	      MPI_Pack(&plasmamain[n].nscat_es, 1, MPI_INT, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD); 
	      MPI_Pack(plasmamain[n].recomb_simple, NTOP_PHOT, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
	      MPI_Pack(&plasmamain[n].kpkt_emiss, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
	      MPI_Pack(&plasmamain[n].kpkt_abs, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
	      MPI_Pack(plasmamain[n].kbf_use, NTOP_PHOT, MPI_INT, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
	      MPI_Pack(&plasmamain[n].kbf_nuse, 1, MPI_INT, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
	      MPI_Pack(&plasmamain[n].t_r, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
	      MPI_Pack(&plasmamain[n].t_r_old, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
	      MPI_Pack(&plasmamain[n].t_e, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
	      MPI_Pack(&plasmamain[n].t_e_old, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
	      MPI_Pack(&plasmamain[n].dt_e, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
	      MPI_Pack(&plasmamain[n].dt_e_old, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
	      MPI_Pack(&plasmamain[n].heat_tot, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
	      MPI_Pack(&plasmamain[n].heat_tot_old, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
	      MPI_Pack(&plasmamain[n].heat_lines, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
	      MPI_Pack(&plasmamain[n].heat_ff, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
	      MPI_Pack(&plasmamain[n].heat_comp, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
	      MPI_Pack(&plasmamain[n].heat_ind_comp, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
	      MPI_Pack(&plasmamain[n].heat_lines_macro, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
	      MPI_Pack(&plasmamain[n].heat_photo_macro, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
	      MPI_Pack(&plasmamain[n].heat_photo, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
	      MPI_Pack(&plasmamain[n].heat_auger, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
	      MPI_Pack(&plasmamain[n].heat_z, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
	      MPI_Pack(&plasmamain[n].w, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
	      MPI_Pack(&plasmamain[n].ntot, 1, MPI_INT, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
	      MPI_Pack(&plasmamain[n].ntot_star, 1, MPI_INT, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
	      MPI_Pack(&plasmamain[n].ntot_bl, 1, MPI_INT, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
	      MPI_Pack(&plasmamain[n].ntot_disk, 1, MPI_INT, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
	      MPI_Pack(&plasmamain[n].ntot_wind, 1, MPI_INT, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
	      MPI_Pack(&plasmamain[n].ntot_agn, 1, MPI_INT, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
	      MPI_Pack(&plasmamain[n].mean_ds, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
	      MPI_Pack(&plasmamain[n].n_ds, 1, MPI_INT, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
	      MPI_Pack(&plasmamain[n].nrad, 1, MPI_INT, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
	      MPI_Pack(&plasmamain[n].nioniz, 1, MPI_INT, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
	      MPI_Pack(plasmamain[n].ioniz, nions, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
	      MPI_Pack(plasmamain[n].recomb, nions, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
	      MPI_Pack(plasmamain[n].scatters, nions, MPI_INT, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
	      MPI_Pack(plasmamain[n].xscatters, nions, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
	      MPI_Pack(plasmamain[n].heat_ion, nions, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
	      MPI_Pack(plasmamain[n].lum_ion, nions, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
	      MPI_Pack(&plasmamain[n].j, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
	      MPI_Pack(&plasmamain[n].j_direct, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
	      MPI_Pack(&plasmamain[n].j_scatt, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
	      MPI_Pack(&plasmamain[n].ave_freq, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
	      MPI_Pack(&plasmamain[n].lum, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
	      MPI_Pack(plasmamain[n].xj, NXBANDS, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
	      MPI_Pack(plasmamain[n].xave_freq, NXBANDS, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
	      MPI_Pack(plasmamain[n].xsd_freq, NXBANDS, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
	      MPI_Pack(plasmamain[n].nxtot, NXBANDS, MPI_INT, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
	      MPI_Pack(&plasmamain[n].max_freq, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
	      MPI_Pack(&plasmamain[n].lum_lines, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
	      MPI_Pack(&plasmamain[n].lum_ff, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
	      MPI_Pack(&plasmamain[n].lum_adiabatic, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
	      MPI_Pack(&plasmamain[n].lum_comp, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
	      MPI_Pack(&plasmamain[n].lum_dr, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
	      MPI_Pack(&plasmamain[n].lum_di, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
	      MPI_Pack(&plasmamain[n].lum_fb, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
	      MPI_Pack(&plasmamain[n].lum_z, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
	      MPI_Pack(&plasmamain[n].lum_rad, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
	      MPI_Pack(&plasmamain[n].lum_rad_old, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
	      MPI_Pack(&plasmamain[n].lum_ioniz, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
	      MPI_Pack(&plasmamain[n].lum_lines_ioniz, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
	      MPI_Pack(&plasmamain[n].lum_ff_ioniz, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
	      MPI_Pack(&plasmamain[n].lum_adiabatic_ioniz, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
	      MPI_Pack(&plasmamain[n].lum_comp_ioniz, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
	      MPI_Pack(&plasmamain[n].lum_dr_ioniz, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
	      MPI_Pack(&plasmamain[n].lum_di_ioniz, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
	      MPI_Pack(&plasmamain[n].lum_fb_ioniz, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
	      MPI_Pack(&plasmamain[n].lum_z_ioniz, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
	      MPI_Pack(&plasmamain[n].lum_rad_ioniz, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
	      MPI_Pack(plasmamain[n].dmo_dt, 3, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
	      MPI_Pack(&plasmamain[n].npdf, 1, MPI_INT, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
	      MPI_Pack(plasmamain[n].pdf_x, LPDF, MPI_INT, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
	      MPI_Pack(plasmamain[n].pdf_y, LPDF, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
	      MPI_Pack(&plasmamain[n].gain, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
	      MPI_Pack(&plasmamain[n].converge_t_r, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
	      MPI_Pack(&plasmamain[n].converge_t_e, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
	      MPI_Pack(&plasmamain[n].converge_hc, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
	      MPI_Pack(&plasmamain[n].trcheck, 1, MPI_INT, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
	      MPI_Pack(&plasmamain[n].techeck, 1, MPI_INT, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
	      MPI_Pack(&plasmamain[n].hccheck, 1, MPI_INT, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
	      MPI_Pack(&plasmamain[n].converge_whole, 1, MPI_INT, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
	      MPI_Pack(&plasmamain[n].converging, 1, MPI_INT, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
	      MPI_Pack(plasmamain[n].gamma_inshl, NAUGER, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
	      MPI_Pack(plasmamain[n].spec_mod_type, NXBANDS, MPI_INT, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
	      MPI_Pack(plasmamain[n].pl_alpha, NXBANDS, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
	      MPI_Pack(plasmamain[n].pl_log_w, NXBANDS, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
	      MPI_Pack(plasmamain[n].exp_temp, NXBANDS, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
	      MPI_Pack(plasmamain[n].exp_w, NXBANDS, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
	      MPI_Pack(plasmamain[n].fmin_mod, NXBANDS, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
	      MPI_Pack(plasmamain[n].fmax_mod, NXBANDS, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
	      MPI_Pack(&plasmamain[n].sim_ip, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
	      MPI_Pack(&plasmamain[n].ferland_ip, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
	      MPI_Pack(&plasmamain[n].ip, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
	      MPI_Pack(&plasmamain[n].ip_direct, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
	      MPI_Pack(&plasmamain[n].ip_scatt, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
	      MPI_Pack(&plasmamain[n].xi, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
	      MPI_Pack(&dt_e, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
	      MPI_Pack(&dt_r, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
	      MPI_Pack(&nmax_e, 1, MPI_INT, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
	      MPI_Pack(&nmax_r, 1, MPI_INT, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
	    }
	  
	  Log("MPI task %d broadcasting plasma update information.\n", rank_global);
	}
      MPI_Barrier(MPI_COMM_WORLD);
      MPI_Bcast(commbuffer, size_of_commbuffer, MPI_PACKED, n_mpi, MPI_COMM_WORLD);
      MPI_Barrier(MPI_COMM_WORLD);
      Log("MPI task %d survived broadcasting plasma update information.\n", rank_global);

      position = 0;
      
      if (rank_global != n_mpi)
	{
	  MPI_Unpack(commbuffer, size_of_commbuffer, &position, &num_comm, 1, MPI_INT, MPI_COMM_WORLD);
	  for (n_mpi2 = 0; n_mpi2 < num_comm; n_mpi2++)
	    {
	      MPI_Unpack(commbuffer, size_of_commbuffer, &position, &n, 1, MPI_INT, MPI_COMM_WORLD);
	      MPI_Unpack(commbuffer, size_of_commbuffer, &position, &plasmamain[n].nwind, 1, MPI_INT, MPI_COMM_WORLD);
	      MPI_Unpack(commbuffer, size_of_commbuffer, &position, &plasmamain[n].nplasma, 1, MPI_INT, MPI_COMM_WORLD);
	      MPI_Unpack(commbuffer, size_of_commbuffer, &position, &plasmamain[n].ne, 1, MPI_DOUBLE, MPI_COMM_WORLD);
	      MPI_Unpack(commbuffer, size_of_commbuffer, &position, &plasmamain[n].rho, 1, MPI_DOUBLE, MPI_COMM_WORLD);
	      MPI_Unpack(commbuffer, size_of_commbuffer, &position, &plasmamain[n].vol, 1, MPI_DOUBLE, MPI_COMM_WORLD);
	      MPI_Unpack(commbuffer, size_of_commbuffer, &position, plasmamain[n].density, nions, MPI_DOUBLE, MPI_COMM_WORLD);
	      MPI_Unpack(commbuffer, size_of_commbuffer, &position, plasmamain[n].partition, nions, MPI_DOUBLE, MPI_COMM_WORLD);
	      MPI_Unpack(commbuffer, size_of_commbuffer, &position, plasmamain[n].levden, NLTE_LEVELS, MPI_DOUBLE, MPI_COMM_WORLD);
	      MPI_Unpack(commbuffer, size_of_commbuffer, &position, plasmamain[n].PWdenom, nions, MPI_DOUBLE, MPI_COMM_WORLD);
	      MPI_Unpack(commbuffer, size_of_commbuffer, &position, plasmamain[n].PWdtemp, nions, MPI_DOUBLE, MPI_COMM_WORLD);
	      MPI_Unpack(commbuffer, size_of_commbuffer, &position, plasmamain[n].PWnumer, nions, MPI_DOUBLE, MPI_COMM_WORLD);
	      MPI_Unpack(commbuffer, size_of_commbuffer, &position, plasmamain[n].PWntemp, nions, MPI_DOUBLE, MPI_COMM_WORLD);	  
	      MPI_Unpack(commbuffer, size_of_commbuffer, &position, &plasmamain[n].kappa_ff_factor, 1, MPI_DOUBLE, MPI_COMM_WORLD);
	      MPI_Unpack(commbuffer, size_of_commbuffer, &position, &plasmamain[n].nscat_es, 1, MPI_INT, MPI_COMM_WORLD);
	      MPI_Unpack(commbuffer, size_of_commbuffer, &position, plasmamain[n].recomb_simple, NTOP_PHOT, MPI_DOUBLE, MPI_COMM_WORLD);	  
	      MPI_Unpack(commbuffer, size_of_commbuffer, &position, &plasmamain[n].kpkt_emiss, 1, MPI_DOUBLE, MPI_COMM_WORLD);
	      MPI_Unpack(commbuffer, size_of_commbuffer, &position, &plasmamain[n].kpkt_abs, 1, MPI_DOUBLE, MPI_COMM_WORLD);
	      MPI_Unpack(commbuffer, size_of_commbuffer, &position, plasmamain[n].kbf_use, NTOP_PHOT, MPI_INT, MPI_COMM_WORLD);	  
	      MPI_Unpack(commbuffer, size_of_commbuffer, &position, &plasmamain[n].kbf_nuse, 1, MPI_INT, MPI_COMM_WORLD);
	      MPI_Unpack(commbuffer, size_of_commbuffer, &position, &plasmamain[n].t_r, 1, MPI_DOUBLE, MPI_COMM_WORLD);
	      MPI_Unpack(commbuffer, size_of_commbuffer, &position, &plasmamain[n].t_r_old, 1, MPI_DOUBLE, MPI_COMM_WORLD);
	      MPI_Unpack(commbuffer, size_of_commbuffer, &position, &plasmamain[n].t_e, 1, MPI_DOUBLE, MPI_COMM_WORLD);
	      MPI_Unpack(commbuffer, size_of_commbuffer, &position, &plasmamain[n].t_e_old, 1, MPI_DOUBLE, MPI_COMM_WORLD);
	      MPI_Unpack(commbuffer, size_of_commbuffer, &position, &plasmamain[n].dt_e, 1, MPI_DOUBLE, MPI_COMM_WORLD);
	      MPI_Unpack(commbuffer, size_of_commbuffer, &position, &plasmamain[n].dt_e_old, 1, MPI_DOUBLE, MPI_COMM_WORLD);
	      MPI_Unpack(commbuffer, size_of_commbuffer, &position, &plasmamain[n].heat_tot, 1, MPI_DOUBLE, MPI_COMM_WORLD);
	      MPI_Unpack(commbuffer, size_of_commbuffer, &position, &plasmamain[n].heat_tot_old, 1, MPI_DOUBLE, MPI_COMM_WORLD);
	      MPI_Unpack(commbuffer, size_of_commbuffer, &position, &plasmamain[n].heat_lines, 1, MPI_DOUBLE, MPI_COMM_WORLD);
	      MPI_Unpack(commbuffer, size_of_commbuffer, &position, &plasmamain[n].heat_ff, 1, MPI_DOUBLE, MPI_COMM_WORLD);
	      MPI_Unpack(commbuffer, size_of_commbuffer, &position, &plasmamain[n].heat_comp, 1, MPI_DOUBLE, MPI_COMM_WORLD);
	      MPI_Unpack(commbuffer, size_of_commbuffer, &position, &plasmamain[n].heat_ind_comp, 1, MPI_DOUBLE, MPI_COMM_WORLD);
	      MPI_Unpack(commbuffer, size_of_commbuffer, &position, &plasmamain[n].heat_lines_macro, 1, MPI_DOUBLE, MPI_COMM_WORLD);
	      MPI_Unpack(commbuffer, size_of_commbuffer, &position, &plasmamain[n].heat_photo_macro, 1, MPI_DOUBLE, MPI_COMM_WORLD);
	      MPI_Unpack(commbuffer, size_of_commbuffer, &position, &plasmamain[n].heat_photo, 1, MPI_DOUBLE, MPI_COMM_WORLD);
	      MPI_Unpack(commbuffer, size_of_commbuffer, &position, &plasmamain[n].heat_auger, 1, MPI_DOUBLE, MPI_COMM_WORLD);
	      MPI_Unpack(commbuffer, size_of_commbuffer, &position, &plasmamain[n].heat_z, 1, MPI_DOUBLE, MPI_COMM_WORLD);
	      MPI_Unpack(commbuffer, size_of_commbuffer, &position, &plasmamain[n].w, 1, MPI_DOUBLE, MPI_COMM_WORLD);
	      MPI_Unpack(commbuffer, size_of_commbuffer, &position, &plasmamain[n].ntot, 1, MPI_INT, MPI_COMM_WORLD);
	      MPI_Unpack(commbuffer, size_of_commbuffer, &position, &plasmamain[n].ntot_star, 1, MPI_INT, MPI_COMM_WORLD);
	      MPI_Unpack(commbuffer, size_of_commbuffer, &position, &plasmamain[n].ntot_bl, 1, MPI_INT, MPI_COMM_WORLD);
	      MPI_Unpack(commbuffer, size_of_commbuffer, &position, &plasmamain[n].ntot_disk, 1, MPI_INT, MPI_COMM_WORLD);
	      MPI_Unpack(commbuffer, size_of_commbuffer, &position, &plasmamain[n].ntot_wind, 1, MPI_INT, MPI_COMM_WORLD);
	      MPI_Unpack(commbuffer, size_of_commbuffer, &position, &plasmamain[n].ntot_agn, 1, MPI_INT, MPI_COMM_WORLD);
	      MPI_Unpack(commbuffer, size_of_commbuffer, &position, &plasmamain[n].mean_ds, 1, MPI_DOUBLE, MPI_COMM_WORLD);
	      MPI_Unpack(commbuffer, size_of_commbuffer, &position, &plasmamain[n].n_ds, 1, MPI_INT, MPI_COMM_WORLD);
	      MPI_Unpack(commbuffer, size_of_commbuffer, &position, &plasmamain[n].nrad, 1, MPI_INT, MPI_COMM_WORLD);
	      MPI_Unpack(commbuffer, size_of_commbuffer, &position, &plasmamain[n].nioniz, 1, MPI_INT, MPI_COMM_WORLD);
	      MPI_Unpack(commbuffer, size_of_commbuffer, &position, plasmamain[n].ioniz, nions, MPI_DOUBLE, MPI_COMM_WORLD);	  
	      MPI_Unpack(commbuffer, size_of_commbuffer, &position, plasmamain[n].recomb, nions, MPI_DOUBLE, MPI_COMM_WORLD);	  
	      MPI_Unpack(commbuffer, size_of_commbuffer, &position, plasmamain[n].scatters, nions, MPI_INT, MPI_COMM_WORLD);	  
	      MPI_Unpack(commbuffer, size_of_commbuffer, &position, plasmamain[n].xscatters, nions, MPI_DOUBLE, MPI_COMM_WORLD);	  
	      MPI_Unpack(commbuffer, size_of_commbuffer, &position, plasmamain[n].heat_ion, nions, MPI_DOUBLE, MPI_COMM_WORLD);	  
	      MPI_Unpack(commbuffer, size_of_commbuffer, &position, plasmamain[n].lum_ion, nions, MPI_DOUBLE, MPI_COMM_WORLD);	  
	      MPI_Unpack(commbuffer, size_of_commbuffer, &position, &plasmamain[n].j, 1, MPI_DOUBLE, MPI_COMM_WORLD);
	      MPI_Unpack(commbuffer, size_of_commbuffer, &position, &plasmamain[n].j_direct, 1, MPI_DOUBLE, MPI_COMM_WORLD);
	      MPI_Unpack(commbuffer, size_of_commbuffer, &position, &plasmamain[n].j_scatt, 1, MPI_DOUBLE, MPI_COMM_WORLD);
	      MPI_Unpack(commbuffer, size_of_commbuffer, &position, &plasmamain[n].ave_freq, 1, MPI_DOUBLE, MPI_COMM_WORLD);
	      MPI_Unpack(commbuffer, size_of_commbuffer, &position, &plasmamain[n].lum, 1, MPI_DOUBLE, MPI_COMM_WORLD);
	      MPI_Unpack(commbuffer, size_of_commbuffer, &position, plasmamain[n].xj, NXBANDS, MPI_DOUBLE, MPI_COMM_WORLD);	  
	      MPI_Unpack(commbuffer, size_of_commbuffer, &position, plasmamain[n].xave_freq, NXBANDS, MPI_DOUBLE, MPI_COMM_WORLD);	  
	      MPI_Unpack(commbuffer, size_of_commbuffer, &position, plasmamain[n].xsd_freq, NXBANDS, MPI_DOUBLE, MPI_COMM_WORLD);	  
	      MPI_Unpack(commbuffer, size_of_commbuffer, &position, plasmamain[n].nxtot, NXBANDS, MPI_INT, MPI_COMM_WORLD);	  
	      MPI_Unpack(commbuffer, size_of_commbuffer, &position, &plasmamain[n].max_freq, 1, MPI_DOUBLE, MPI_COMM_WORLD);
	      MPI_Unpack(commbuffer, size_of_commbuffer, &position, &plasmamain[n].lum_lines, 1, MPI_DOUBLE, MPI_COMM_WORLD);
	      MPI_Unpack(commbuffer, size_of_commbuffer, &position, &plasmamain[n].lum_ff, 1, MPI_DOUBLE, MPI_COMM_WORLD);
	      MPI_Unpack(commbuffer, size_of_commbuffer, &position, &plasmamain[n].lum_adiabatic, 1, MPI_DOUBLE, MPI_COMM_WORLD);
	      MPI_Unpack(commbuffer, size_of_commbuffer, &position, &plasmamain[n].lum_comp, 1, MPI_DOUBLE, MPI_COMM_WORLD);
	      MPI_Unpack(commbuffer, size_of_commbuffer, &position, &plasmamain[n].lum_dr, 1, MPI_DOUBLE, MPI_COMM_WORLD);
	      MPI_Unpack(commbuffer, size_of_commbuffer, &position, &plasmamain[n].lum_di, 1, MPI_DOUBLE, MPI_COMM_WORLD);
	      MPI_Unpack(commbuffer, size_of_commbuffer, &position, &plasmamain[n].lum_fb, 1, MPI_DOUBLE, MPI_COMM_WORLD);
	      MPI_Unpack(commbuffer, size_of_commbuffer, &position, &plasmamain[n].lum_z, 1, MPI_DOUBLE, MPI_COMM_WORLD);
	      MPI_Unpack(commbuffer, size_of_commbuffer, &position, &plasmamain[n].lum_rad, 1, MPI_DOUBLE, MPI_COMM_WORLD);
	      MPI_Unpack(commbuffer, size_of_commbuffer, &position, &plasmamain[n].lum_rad_old, 1, MPI_DOUBLE, MPI_COMM_WORLD);
	      MPI_Unpack(commbuffer, size_of_commbuffer, &position, &plasmamain[n].lum_ioniz, 1, MPI_DOUBLE, MPI_COMM_WORLD);
	      MPI_Unpack(commbuffer, size_of_commbuffer, &position, &plasmamain[n].lum_lines_ioniz, 1, MPI_DOUBLE, MPI_COMM_WORLD);
	      MPI_Unpack(commbuffer, size_of_commbuffer, &position, &plasmamain[n].lum_ff_ioniz, 1, MPI_DOUBLE, MPI_COMM_WORLD);
	      MPI_Unpack(commbuffer, size_of_commbuffer, &position, &plasmamain[n].lum_adiabatic_ioniz, 1, MPI_DOUBLE, MPI_COMM_WORLD);
	      MPI_Unpack(commbuffer, size_of_commbuffer, &position, &plasmamain[n].lum_comp_ioniz, 1, MPI_DOUBLE, MPI_COMM_WORLD);
	      MPI_Unpack(commbuffer, size_of_commbuffer, &position, &plasmamain[n].lum_dr_ioniz, 1, MPI_DOUBLE, MPI_COMM_WORLD);
	      MPI_Unpack(commbuffer, size_of_commbuffer, &position, &plasmamain[n].lum_di_ioniz, 1, MPI_DOUBLE, MPI_COMM_WORLD);
	      MPI_Unpack(commbuffer, size_of_commbuffer, &position, &plasmamain[n].lum_fb_ioniz, 1, MPI_DOUBLE, MPI_COMM_WORLD);
	      MPI_Unpack(commbuffer, size_of_commbuffer, &position, &plasmamain[n].lum_z_ioniz, 1, MPI_DOUBLE, MPI_COMM_WORLD);
	      MPI_Unpack(commbuffer, size_of_commbuffer, &position, &plasmamain[n].lum_rad_ioniz, 1, MPI_DOUBLE, MPI_COMM_WORLD);
	      MPI_Unpack(commbuffer, size_of_commbuffer, &position, plasmamain[n].dmo_dt, 3, MPI_DOUBLE, MPI_COMM_WORLD);
	      MPI_Unpack(commbuffer, size_of_commbuffer, &position, &plasmamain[n].npdf, 1, MPI_INT, MPI_COMM_WORLD);
	      MPI_Unpack(commbuffer, size_of_commbuffer, &position, plasmamain[n].pdf_x, LPDF, MPI_INT, MPI_COMM_WORLD);
	      MPI_Unpack(commbuffer, size_of_commbuffer, &position, plasmamain[n].pdf_y, LPDF, MPI_DOUBLE, MPI_COMM_WORLD);
	      MPI_Unpack(commbuffer, size_of_commbuffer, &position, &plasmamain[n].gain, 1, MPI_DOUBLE, MPI_COMM_WORLD);
	      MPI_Unpack(commbuffer, size_of_commbuffer, &position, &plasmamain[n].converge_t_r, 1, MPI_DOUBLE, MPI_COMM_WORLD);
	      MPI_Unpack(commbuffer, size_of_commbuffer, &position, &plasmamain[n].converge_t_e, 1, MPI_DOUBLE, MPI_COMM_WORLD);
	      MPI_Unpack(commbuffer, size_of_commbuffer, &position, &plasmamain[n].converge_hc, 1, MPI_DOUBLE, MPI_COMM_WORLD);
	      MPI_Unpack(commbuffer, size_of_commbuffer, &position, &plasmamain[n].trcheck, 1, MPI_INT, MPI_COMM_WORLD);
	      MPI_Unpack(commbuffer, size_of_commbuffer, &position, &plasmamain[n].techeck, 1, MPI_INT, MPI_COMM_WORLD);
	      MPI_Unpack(commbuffer, size_of_commbuffer, &position, &plasmamain[n].hccheck, 1, MPI_INT, MPI_COMM_WORLD);
	      MPI_Unpack(commbuffer, size_of_commbuffer, &position, &plasmamain[n].converge_whole, 1, MPI_INT, MPI_COMM_WORLD);
	      MPI_Unpack(commbuffer, size_of_commbuffer, &position, &plasmamain[n].converging, 1, MPI_INT, MPI_COMM_WORLD);
	      MPI_Unpack(commbuffer, size_of_commbuffer, &position, plasmamain[n].gamma_inshl, NAUGER, MPI_DOUBLE, MPI_COMM_WORLD);
	      MPI_Unpack(commbuffer, size_of_commbuffer, &position, plasmamain[n].spec_mod_type, NXBANDS, MPI_INT, MPI_COMM_WORLD);
	      MPI_Unpack(commbuffer, size_of_commbuffer, &position, plasmamain[n].pl_alpha, NXBANDS, MPI_DOUBLE, MPI_COMM_WORLD);
	      MPI_Unpack(commbuffer, size_of_commbuffer, &position, plasmamain[n].pl_log_w, NXBANDS, MPI_DOUBLE, MPI_COMM_WORLD);
	      MPI_Unpack(commbuffer, size_of_commbuffer, &position, plasmamain[n].exp_temp, NXBANDS, MPI_DOUBLE, MPI_COMM_WORLD);
	      MPI_Unpack(commbuffer, size_of_commbuffer, &position, plasmamain[n].exp_w, NXBANDS, MPI_DOUBLE, MPI_COMM_WORLD);
	      MPI_Unpack(commbuffer, size_of_commbuffer, &position, plasmamain[n].fmin_mod, NXBANDS, MPI_DOUBLE, MPI_COMM_WORLD);
	      MPI_Unpack(commbuffer, size_of_commbuffer, &position, plasmamain[n].fmax_mod, NXBANDS, MPI_DOUBLE, MPI_COMM_WORLD);
	      MPI_Unpack(commbuffer, size_of_commbuffer, &position, &plasmamain[n].sim_ip, 1, MPI_DOUBLE, MPI_COMM_WORLD);
	      MPI_Unpack(commbuffer, size_of_commbuffer, &position, &plasmamain[n].ferland_ip, 1, MPI_DOUBLE, MPI_COMM_WORLD);
	      MPI_Unpack(commbuffer, size_of_commbuffer, &position, &plasmamain[n].ip, 1, MPI_DOUBLE, MPI_COMM_WORLD);
	      MPI_Unpack(commbuffer, size_of_commbuffer, &position, &plasmamain[n].ip_direct, 1, MPI_DOUBLE, MPI_COMM_WORLD);
	      MPI_Unpack(commbuffer, size_of_commbuffer, &position, &plasmamain[n].ip_scatt, 1, MPI_DOUBLE, MPI_COMM_WORLD);
	      MPI_Unpack(commbuffer, size_of_commbuffer, &position, &plasmamain[n].xi, 1, MPI_DOUBLE, MPI_COMM_WORLD);
	      MPI_Unpack(commbuffer, size_of_commbuffer, &position, &dt_e_temp, 1, MPI_DOUBLE, MPI_COMM_WORLD);
	      MPI_Unpack(commbuffer, size_of_commbuffer, &position, &dt_r_temp, 1, MPI_DOUBLE, MPI_COMM_WORLD);
	      MPI_Unpack(commbuffer, size_of_commbuffer, &position, &nmax_e_temp, 1, MPI_INT, MPI_COMM_WORLD);
	      MPI_Unpack(commbuffer, size_of_commbuffer, &position, &nmax_r_temp, 1, MPI_INT, MPI_COMM_WORLD);

          /* JM 1409 -- Altered for issue #110 to ensure correct reporting in parallel */
	      if (fabs(dt_e_temp) >= fabs(dt_e))
	      {
	      	/* Check if any other threads found a higher maximum for te */
	      	dt_e = dt_e_temp;
	      	nmax_e = nmax_e_temp;
	      }

	      if (fabs(dt_r_temp) >= fabs(dt_r))
	      {
	      	/* Check if any other threads found a higher maximum for tr */
	      	dt_r = dt_r_temp;
	      	nmax_r = nmax_r_temp;
	      }

	      t_r_ave += plasmamain[n].t_r;
	      t_e_ave += plasmamain[n].t_e;
	      iave++;
	      
	    }

	}

    }
free (commbuffer);
#endif


  /* Now we need to updated the densities immediately outside the wind so that the density interpolation in resonate will work.
     In this case all we have done is to copy the densities from the cell which is just in the wind (as one goes outward) to the
     cell that is just inside (or outside) the wind. 

     SS asked whether we should also be extending the wind for other parameters, especially ne.  At present we do not interpolate
     on ne so this is not necessary.  If we did do that it would be required.

     In cylindrical coordinates, the fast dimension is z; grid positions increase up in z, and then out in r.
     In spperical polar coordinates, the fast dimension is theta; the grid increases in theta (measured)
     from the z axis), and then in r.
     In spherical coordinates, the grid increases as one might expect in r..
     *
   */

  if (geo.coord_type == CYLIND)
    cylind_extend_density (w);
  else if (geo.coord_type == RTHETA)
    rtheta_extend_density (w);
  else if (geo.coord_type == SPHERICAL)
    spherical_extend_density (w);
  else if (geo.coord_type == CYLVAR)
    cylvar_extend_density (w);
  else
    {
      Error ("Wind_update2d: Unknown coordinate type %d\n", geo.coord_type);
      exit (0);
    }

  /* Finished updatating region outside of wind */

  num_updates++;
  strcpy (string, "");
  sprintf (string, "# Wind update: Number %d", num_updates);

  if (modes.zeus_connect==1) //If we are running in zeus connect mode - we open a file for heatcool rates
  {
	  Log("Outputting heatcool file for connecting to zeus\n");
      fptr = fopen ("py_heatcool.dat", "w");
  	 fprintf(fptr,"i j rcen thetacen vol temp xi ne heat_xray heat_comp heat_lines heat_ff cool_comp cool_lines cool_ff rho n_h\n");
	  
  }

  /* Check the balance between the absorbed and the emitted flux */

  xsum = psum = ausum = lsum = fsum = csum = icsum = 0;	//1108 NSH zero the new csum counter for compton heating

  for (nplasma = 0; nplasma < NPLASMA; nplasma++)
    {
      if (sane_check (plasmamain[nplasma].heat_tot))
	Error ("wind_update:sane_check w(%d).heat_tot is %e\n", nplasma,
	       plasmamain[nplasma].heat_tot);
      if (sane_check (plasmamain[nplasma].heat_photo))
	Error ("wind_update:sane_check w(%d).heat_photo is %e\n", nplasma,
	       plasmamain[nplasma].heat_photo);
	  if (sane_check (plasmamain[nplasma].heat_auger))
	 	Error ("wind_update:sane_check w(%d).heat_auger is %e\n", nplasma,
	 	       plasmamain[nplasma].heat_auger);
      if (sane_check (plasmamain[nplasma].heat_photo_macro))
	Error ("wind_update:sane_check w(%d).heat_photo_macro is %e\n",
	       nplasma, plasmamain[nplasma].heat_photo_macro);
      if (sane_check (plasmamain[nplasma].heat_ff))
	Error ("wind_update:sane_check w(%d).heat_ff is %e\n", nplasma,
	       plasmamain[nplasma].heat_ff);
      if (sane_check (plasmamain[nplasma].heat_lines))
	Error ("wind_update:sane_check w(%d).heat_lines is %e\n", nplasma,
	       plasmamain[nplasma].heat_lines);
      if (sane_check (plasmamain[nplasma].heat_lines_macro))
	Error ("wind_update:sane_check w(%d).heat_lines_macro is %e\n",
	       nplasma, plasmamain[nplasma].heat_lines_macro);
      /* 1108 NSH extra Sane check for compton heating */
      if (sane_check (plasmamain[nplasma].heat_comp))
	Error ("wind_update:sane_check w(%d).heat_comp is %e\n", nplasma,
	       plasmamain[nplasma].heat_comp);
      xsum += plasmamain[nplasma].heat_tot;
      psum += plasmamain[nplasma].heat_photo;
      ausum += plasmamain[nplasma].heat_auger;
      fsum += plasmamain[nplasma].heat_ff;
      lsum += plasmamain[nplasma].heat_lines;
      csum += plasmamain[nplasma].heat_comp;	//1108 NSH Increment the compton heating counter
      icsum += plasmamain[nplasma].heat_ind_comp;	//1205 NSH Increment the induced compton heating counter

// Comment - ksl - Next line generated an inappopriate amount of output in the .diag file
// OLD           Log ("OUTPUT logIP(cloudy_thoeretical)= %e logIP(cloudy_actual)=%e\n",
// OLD           log10(plasmamain[nplasma].ferland_ip),log10(plasmamain[nplasma].ip)); 
      


      /* JM130621- bugfix for windsave bug- needed so that we have the luminosities from ionization
         cycles in the windsavefile even if the spectral cycles are run */
      plasmamain[nplasma].lum_ioniz = plasmamain[nplasma].lum;
      plasmamain[nplasma].lum_ff_ioniz = plasmamain[nplasma].lum_ff;
      plasmamain[nplasma].lum_fb_ioniz = plasmamain[nplasma].lum_fb;
      plasmamain[nplasma].lum_z_ioniz = plasmamain[nplasma].lum_z;
      plasmamain[nplasma].lum_lines_ioniz = plasmamain[nplasma].lum_lines;
      plasmamain[nplasma].lum_comp_ioniz = plasmamain[nplasma].lum_comp;
      plasmamain[nplasma].lum_dr_ioniz = plasmamain[nplasma].lum_dr;
      plasmamain[nplasma].lum_di_ioniz = plasmamain[nplasma].lum_di;
      plasmamain[nplasma].lum_rad_ioniz = plasmamain[nplasma].lum_rad;
      plasmamain[nplasma].lum_adiabatic_ioniz = plasmamain[nplasma].lum_adiabatic;




    }
	


  /* JM130621- bugfix for windsave bug- needed so that we have the luminosities from ionization
     cycles in the windsavefile even if the spectral cycles are run */
  geo.lum_ff_ioniz = geo.lum_ff;
  geo.lum_fb_ioniz = geo.lum_fb;
  geo.lum_lines_ioniz = geo.lum_lines;
  geo.lum_comp_ioniz = geo.lum_comp;
  geo.lum_dr_ioniz = geo.lum_dr;
  geo.lum_di_ioniz = geo.lum_di;
  geo.lum_adiabatic_ioniz = geo.lum_adiabatic;
  geo.lum_disk_ioniz = geo.lum_disk;
  geo.lum_star_ioniz = geo.lum_star;
  geo.lum_bl_ioniz = geo.lum_bl;
  geo.lum_wind_ioniz = geo.lum_wind;
  geo.lum_tot_ioniz = geo.lum_tot;

  /* Added this system which counts number of times two situations occur (See #91)
     We only report these every 100,000 times (one can typically get ) */
  Log ("wind_update: note, errors from mean intensity can be high in a working model\n");
  Log ("wind_update: can be a problem with photon numbers if there are also errors from spectral_estimators and low photon number warnings\n");
  Log ("wind_update: mean_intensity: %8.4e occurrences, this cycle, this thread of 'no model exists in a band'\n",
  	  	      (1.0*nerr_no_Jmodel));
  Log ("wind_update: mean intensity: %8.4e occurrences, this cycle, this thread of 'photon freq is outside frequency range of spectral model'\n",
  	  	     (1.0*nerr_Jmodel_wrong_freq));


  /* zero the counters which record diagnositics from mean_intensity */
  nerr_Jmodel_wrong_freq = 0;
  nerr_no_Jmodel = 0;



  asum = wind_luminosity (0.0, VERY_BIG);  /*We call wind_luminosity here to obtain an up to date set of cooling rates*/


  if (modes.zeus_connect==1) //If we are running in zeus connect mode, we output heating and cooling rates.
  {
	  for (nplasma = 0; nplasma < NPLASMA; nplasma++)
	  {
	  
 		 wind_n_to_ij (plasmamain[nplasma].nwind, &i, &j);
		 vol=w[plasmamain[nplasma].nwind].vol;
	  fprintf(fptr,"%d %d %e %e %e ",i,j,w[plasmamain[nplasma].nwind].rcen,w[plasmamain[nplasma].nwind].thetacen/RADIAN,vol); //output geometric things
	  fprintf(fptr,"%e %e %e ",plasmamain[nplasma].t_e,plasmamain[nplasma].xi,plasmamain[nplasma].ne); //output temp, xi and ne to ease plotting of heating rates
	  fprintf(fptr,"%e ",(plasmamain[nplasma].heat_photo+plasmamain[nplasma].heat_auger)/vol); //Xray heating - or photoionization
	  fprintf(fptr,"%e ",(plasmamain[nplasma].heat_comp)/vol); //Compton heating
	  fprintf(fptr,"%e ",(plasmamain[nplasma].heat_lines)/vol); //Line heating 28/10/15 - not currently used in zeus
	  fprintf(fptr,"%e ",(plasmamain[nplasma].heat_ff)/vol); //FF heating 28/10/15 - not currently used in zeus
	  fprintf(fptr,"%e ",(plasmamain[nplasma].lum_comp)/vol); //Compton cooling
	  fprintf(fptr,"%e ",(plasmamain[nplasma].lum_lines+plasmamain[nplasma].lum_fb+plasmamain[nplasma].lum_dr)/vol); //Line cooling must include all recombinatiobs cooling
	  fprintf(fptr,"%e ",(plasmamain[nplasma].lum_ff)/vol); //ff cooling
	  fprintf(fptr,"%e ",plasmamain[nplasma].rho); //density
	  fprintf(fptr,"%e\n",plasmamain[nplasma].rho*rho2nh); //hydrogen number density
   }
    fclose(fptr);
}

  /* 1108 NSH Added commands to report compton heating */
  Log ("!!wind_update: Absorbed flux    %8.2e  (photo %8.2e ff %8.2e compton %8.2e auger %8.2e induced_compton %8.2e lines %8.2e)\n", 
        xsum, psum, fsum, csum, ausum, icsum, lsum);
	
  /* 1306 Added line to split out absorbed flux from wind heating */
  Log ("!!wind_update: Wind heating     %8.2e  (photo %8.2e ff %8.2e compton %8.2e auger %8.2e induced_compton %8.2e lines %8.2e adiabatic %8.2e)\n", 
        xsum+geo.heat_adiabatic, psum, fsum, csum, ausum, icsum, lsum, geo.heat_adiabatic);
	
  /* 1108 NSH added commands to report compton cooling 1110 removed, 
   * this line now just reports cooling mechanisms that will generate photons */
  Log ("!!wind_update: Wind luminosity  %8.2e (recomb %8.2e ff %8.2e lines %8.2e) after update\n", 
        asum, geo.lum_fb, geo.lum_ff, geo.lum_lines);	

  /* 1110 NSH Added this line to report all cooling mechanisms, including those that do not generate photons. */
  Log ("!!wind_update: Wind cooling     %8.2e (recomb %8.2e ff %8.2e compton %8.2e DR %8.2e DI %8.2e lines %8.2e adiabatic %8.2e) after update\n", 
        asum + geo.lum_comp + geo.lum_dr + geo.lum_di + geo.lum_adiabatic, geo.lum_fb, geo.lum_ff, geo.lum_comp, 
        geo.lum_dr, geo.lum_di, geo.lum_lines, geo.lum_adiabatic);	



  /* Print out some diagnositics of the changes in the wind update */
		if (modes.zeus_connect==1 || modes.fixed_temp==1)	     //There is no point in computing temperature changes, because we have fixed them!
   	{
   		Log ("!!wind_update: We are running in fixed temperature mode - no temperature report\n");
		
   	} 
		else
		{
 	t_r_ave_old /= iave;
  t_e_ave_old /= iave;
  t_r_ave /= iave;
  t_e_ave /= iave;
  wind_n_to_ij (nmax_r, &i, &j);
  Log ("!!wind_update: Max change in t_r %6.0f at cell %4d (%d,%d)\n", dt_r,
       nmax_r, i, j);
  Log ("!!wind_update: Ave change in t_r %6.0f from %6.0f to %6.0f\n",
       (t_r_ave - t_r_ave_old), t_r_ave_old, t_r_ave);
  wind_n_to_ij (nmax_e, &i, &j);
  Log ("!!wind_update: Max change in t_e %6.0f at cell %4d (%d,%d)\n", dt_e,
       nmax_e, i, j);
  Log ("!!wind_update: Ave change in t_e %6.0f from %6.0f to %6.0f\n",
       (t_e_ave - t_e_ave_old), t_e_ave_old, t_e_ave);

  Log ("Summary  t_e  %6.0f   %6.0f  #t_e and dt_e on this update\n", t_e_ave,
       (t_e_ave - t_e_ave_old));
  Log ("Summary  t_r  %6.0f   %6.0f  #t_r and dt_r on this update\n", t_r_ave,
       (t_r_ave - t_r_ave_old));
   } 
 

  check_convergence ();
  /* Summarize the radiative temperatures (ksl 04 mar)*/
  xtemp_rad (w);




/* This next block is to allow the output of data relating to the abundances of ions when python is being tested
 * with thin shell mode.We will only want this to run if the wind mode is 9, for test or thin shell mode.
 */

  if (geo.wind_type == 9)
    {
      n = plasmamain[0].nwind;
      for (i = 0; i < geo.nxfreq; i++)	/*loop over number of bands */
	{
	  Log
	    ("Band %i f1 %e f2 %e model %i pl_alpha %f pl_log_w %e exp_t %e exp_w %e\n",
	     i, geo.xfreq[i], geo.xfreq[i + 1],
	     plasmamain[0].spec_mod_type[i], plasmamain[0].pl_alpha[i],
	     plasmamain[0].pl_log_w[i], plasmamain[0].exp_temp[i],
	     plasmamain[0].exp_w[i]);
	}
      /* Get some line diagnostics */
      nsh_lum_hhe = 0.0;
      nsh_lum_metals = 0.0;

      for (i = 0; i < nlines; i++)
	{
	  if (lin_ptr[i]->z < 3)
	    nsh_lum_hhe = nsh_lum_hhe + lin_ptr[i]->pow;
	  else
	    nsh_lum_metals = nsh_lum_metals + lin_ptr[i]->pow;
	}
      agn_ip =
	geo.const_agn *
	(((pow (50000 / HEV, geo.alpha_agn + 1.0)) -
	  pow (100 / HEV, geo.alpha_agn + 1.0)) / (geo.alpha_agn + 1.0));
      agn_ip /= (w[n].r * w[n].r);
      agn_ip /= plasmamain[0].rho * rho2nh;
      /* Report luminosities, IP and other diagnositic quantities */
      Log
	("OUTPUT Lum_agn= %e T_e= %e N_h= %e N_e= %e alpha= %f IP(sim_2010)= %e Meaured_IP(cloudy)= %e Measured_Xi= %e distance= %e volume= %e mean_ds=%e\n",
	 geo.lum_agn, plasmamain[0].t_e, plasmamain[0].rho * rho2nh,
	 plasmamain[0].ne, geo.alpha_agn, agn_ip, plasmamain[0].ip, plasmamain[0].xi,w[n].r,
	 w[n].vol, plasmamain[0].mean_ds / plasmamain[0].n_ds);

      /* 1108 NSH Added commands to report compton heating */
      Log ("OUTPUT Absorbed_flux(ergs-1cm-3)    %8.2e  (photo %8.2e ff %8.2e compton %8.2e induced_compton %8.2e lines %8.2e auger %8.2e )\n", 
            xsum / w[n].vol, psum / w[n].vol,  fsum / w[n].vol, csum / w[n].vol, icsum / w[n].vol, lsum / w[n].vol, ausum / w[n].vol);	

      /* 1110 NSH Added this line to report all cooling mechanisms, including those that do not generate photons.*/
      Log ("OUTPUT Wind_cooling(ergs-1cm-3)     %8.2e (recomb %8.2e ff %8.2e compton %8.2e DR %8.2e DI %8.2e adiabatic %8.2e lines %8.2e ) after update\n",
           (asum + geo.lum_comp + geo.lum_dr + geo.lum_di +geo.lum_adiabatic) / w[n].vol, geo.lum_fb / w[n].vol, geo.lum_ff / w[n].vol, 
           geo.lum_comp / w[n].vol, geo.lum_dr / w[n].vol,  geo.lum_di / w[n].vol,geo.lum_adiabatic / w[n].vol, geo.lum_lines / w[n].vol);	

      Log ("OUTPUT Wind_line_cooling(ergs-1cm-3)  HHe %8.2e Metals %8.2e\n",
	   nsh_lum_hhe / w[n].vol, nsh_lum_metals / w[n].vol);

    /* 1110 NSH Added this line to report all cooling mechanisms, including those that do not generate photons. */
      Log ("OUTPUT Balance      Cooling=%8.2e Heating=%8.2e Lum=%8.2e T_e=%e after update\n", asum + geo.lum_comp + geo.lum_dr + geo.lum_di + geo.lum_adiabatic, xsum, asum, plasmamain[0].t_e);

      for (n = 0; n < nelements; n++)
	{
	  first = ele[n].firstion;
	  last = first + ele[n].nions;
	  Log ("OUTPUT %-5s ", ele[n].name);
	  tot = 0;
	  for (m = first; m < last; m++)
	    tot += plasmamain[0].density[m];
	  for (m = first; m < last; m++)
	    {
	      Log (" %8.2e", plasmamain[0].density[m] / tot);
	    }
	  Log ("\n");
	}
    }







  return (0);
}


/***********************************************************
                                       Space Telescope Science Institute

 Synopsis: wind_rad_init(w) zeros those portions of the wind which contain the radiation properties 
	of the wind, i.e those portions which should be set to zeroed when the structure of the 
	wind has been changed or when you simply want to start off a calculation in a known state

 Arguments:		
	WindPtr w;

Returns:
 
Description:	
	
Notes:
	The routine is called by the main program at the beginning of each ionization calculation
	cycle.  It should zero all heating and radiation induced cooling in the wind array.  Since
	cooling is recalculated in wind_update, one needs to be sure that all of the appropriate
	cooling terms are also rezeroed there as well.
	

History:
 	97jan	ksl	Coding on python began.
 	98apr	ksl	Fixed major problem with this routine, previously I had only been 
 			reinitializing the first NDIM portions of the wind.
 	98may	ksl	Added initialization for .ioniz and .recomb
 	98sept	ksl	Added initialization for heating of h, he1, he2,z. and also for the remaining
 			luminosities.
	03jul	ksl	Updated to reflect modified wind structure with heat and cooling
			for all ions
	04mar	ss	Added additional quantities to zero for macro-atoms
        04July  SS      Modified to include setting the spontaneous recombination rates.
        04 Nov  SS      Modified to include setting the recomb_simple rates (needed for cooling in
                        simple continuua).
	06may	ksl	Modified to use wmain and plasmamain, in process of redefining all of
			the structures.  Most if not all of this will end up in plasmamain
	06jul	ksl	57+ -- Added to changes to allow for a separate macro structure
	06aug	ksl	57h -- Additional changes to allow for the fact that marcomain
			is not created at all if no macro atoms.
	13dec	nsh	77 zero various new plasma variables

**************************************************************/


int
wind_rad_init ()
{
  int n, i;
  int njump;


  for (n = 0; n < NPLASMA; n++)
    {
      plasmamain[n].j = plasmamain[n].ave_freq = plasmamain[n].ntot = 0;
      plasmamain[n].j_direct = plasmamain[n].j_scatt = 0,0;  //NSH 1309 zero j banded by number of scatters
      plasmamain[n].ip = 0.0;
      plasmamain[n].xi = 0.0;

      plasmamain[n].ip_direct = plasmamain[n].ip_scatt = 0.0;
      plasmamain[n].mean_ds = 0.0;
      plasmamain[n].n_ds = 0;
      plasmamain[n].ntot_disk = plasmamain[n].ntot_agn = 0;	//NSH 15/4/11 counters to see where photons come from
      plasmamain[n].ntot_star = plasmamain[n].ntot_bl =
	plasmamain[n].ntot_wind = 0;
      plasmamain[n].heat_tot = plasmamain[n].heat_ff =
	plasmamain[n].heat_photo = plasmamain[n].heat_lines = 0.0;
      plasmamain[n].heat_z = 0.0;
      plasmamain[n].max_freq = 0.0;	//NSH 120814 Zero the counter which works out the maximum frequency seen in a cell and hence the maximum applicable frequency of the power law estimators.
      plasmamain[n].lum = plasmamain[n].lum_rad = plasmamain[n].lum_lines =
	plasmamain[n].lum_ff = 0.0;
      plasmamain[n].lum_fb = plasmamain[n].lum_z = 0.0;
      plasmamain[n].nrad = plasmamain[n].nioniz = 0;
      plasmamain[n].lum_comp = 0.0;	//1108 NSH Zero the compton luminosity for the cell
      plasmamain[n].heat_comp = 0.0;	//1108 NSH Zero the compton heating for the cell
      plasmamain[n].heat_ind_comp = 0.0;	//1108 NSH Zero the induced compton heating for the cell
      plasmamain[n].heat_auger = 0.0;	//1108 NSH Zero the auger heating for the cell
	  
      if (nlevels_macro > 1)
	macromain[n].kpkt_rates_known = -1;

/* 1108 NSH Loop to zero the frequency banded radiation estimators */
/* 71 - 111279 - ksl - Small modification to reflect the fact that nxfreq has been moved into the geo structure */
      for (i = 0; i < geo.nxfreq; i++)
	{
	  plasmamain[n].xj[i] = plasmamain[n].xave_freq[i] =
	    plasmamain[n].nxtot[i] = 0;
	  plasmamain[n].xsd_freq[i] = 0.0;	/* NSH 120815 Zero the standard deviation counter */
          plasmamain[n].fmin[i] = geo.xfreq[i+1]; /* Set the minium frequency to the max frequency in the band */
          plasmamain[n].fmax[i] = geo.xfreq[i]; /* Set the maximum frequency to the min frequency in the band */
	}




      for (i = 0; i < nions; i++)
	{
	  plasmamain[n].ioniz[i] = plasmamain[n].recomb[i] =
	    plasmamain[n].heat_ion[i] = plasmamain[n].lum_ion[i] = 0.0;
	  
	}

      /*Block added (Dec 08) to zero the auger rate estimators */
      for (i = 0; i < nauger; i++)
	{
	  plasmamain[n].gamma_inshl[i] = 0.0;
	}

      /* Next blocks added by SS Mar 2004 to zero the Macro Atom estimators. */

      /* 57h -- 0608 -- These sections actually involve enough calculations that
         they are noticeable in term sof the overall speed.  One would if possible
         like to avoid this section, since it requires the creation of macromain, 
         even though macromain is not used -- ksl */


      for (i = 0; i < nlevels_macro; i++)	//57h
	{
	  for (njump = 0; njump < config[i].n_bbu_jump; njump++)
	    {
	      macromain[n].jbar[config[i].bbu_indx_first + njump] = 0.0;	// mean intensity
	    }
	  for (njump = 0; njump < config[i].n_bfu_jump; njump++)
	    {
	      macromain[n].gamma[config[i].bfu_indx_first + njump] = 0.0;
	      macromain[n].gamma_e[config[i].bfu_indx_first + njump] = 0.0;
	      macromain[n].alpha_st[config[i].bfd_indx_first + njump] = 0.0;	//stimulated recombination
	      macromain[n].alpha_st_e[config[i].bfd_indx_first + njump] = 0.0;
	    }


	  /* Next block to set spontaneous recombination rates for next iteration. (SS July 04) */
	  for (njump = 0; njump < config[i].n_bfd_jump; njump++)
	    {
	      if (plasmamain[n].t_e > 1.0)
		{
		  //04Jul--ksl-modified these calls to reflect changed alpha_sp
		  macromain[n].recomb_sp[config[i].bfd_indx_first + njump] =
		    alpha_sp (&phot_top[config[i].bfd_jump[njump]],
			      &plasmamain[n], 0);
		  macromain[n].recomb_sp_e[config[i].bfd_indx_first + njump] =
		    alpha_sp (&phot_top[config[i].bfd_jump[njump]],
			      &plasmamain[n], 2);
		}
	      else
		{
		  macromain[n].recomb_sp[config[i].bfd_indx_first + njump] =
		    0.0;
		  macromain[n].recomb_sp_e[config[i].bfd_indx_first + njump] =
		    0.0;
		}

	    }
	}

      for (i = 0; i < ntop_phot; i++)
	{
	  /* 57h -- recomb_simple is only required for we are using a macro atom approach, and only non-zero when
	     this particular phot_tob xsection is treated as a simple x-section. Stuart, is this correct?? I've added
	     checks so that macro_info is only 0 (false) or true (1), and so the logic of the next section can be 
	     simplified.  0608-ksl */
	  if (geo.macro_simple || phot_top[i].macro_info)
	    {

	      plasmamain[n].recomb_simple[i] = 0.0;
	    }
	  else
	    {			//we want a macro approach, but not for this ion so need recomb_simple
	      plasmamain[n].recomb_simple[i] =
		alpha_sp (&phot_top[i], &plasmamain[n], 2);
	    }
	}


      //zero the emissivities that are needed for the spectral synthesis step.
      plasmamain[n].kpkt_emiss = 0.0;
      plasmamain[n].kpkt_abs = 0.0;
      for (i = 0; i < nlevels_macro; i++)	//57h
	{
	  macromain[n].matom_abs[i] = 0.0;

	  macromain[n].matom_emiss[i] = 0.0;

	}

      /* End of added material. */
    }


  return (0);
}

/***********************************************************
                                       Space Telescope Science Institute

 Synopsis: wind_rad_summary(w,filename,mode) writes out a file which sumarizes the
	the properties of the radiation field in each cell in the wind. 
Arguments:		
	WindPtr w;
	char filename[]		Name of the output file
	char mode[]		Either "w" to write a new file or "a" to append to
						and existing file
Returns:
 
Description:	
	
Notes:
	There is not much point in calling this until you have propagated a few photons
History:
 	97jan      ksl	Coding on python began.
	06may	ksl	57+: Began modificatioins to split structures.  Much of this
			will end up in the plasma structure.

**************************************************************/

int istat_wind_rad_summary = 0;

int
wind_rad_summary (w, filename, mode)
     WindPtr w;
     char filename[], mode[];
{
  FILE *fopen (), *fptr;
  int n;
  int nplasma;


  /* Open or reopen the file where you want to put the summary of the radiation properties of the wind */
  if (istat_wind_rad_summary == 0)
    {				/*First time through recreate the file regardless */
      if ((fptr = fopen (filename, "w")) == NULL)
	{
	  Error ("wind_rad_summary: Unable to open %s\n", filename);
	  exit (0);
	}

    }
  else
    {
      if ((fptr = fopen (filename, mode)) == NULL)
	{
	  Error ("wind_rad_summary: Unable to reopen %s with mode %s\n",
		 filename, mode);
	  exit (0);
	}
    }

  istat_wind_rad_summary++;

  fprintf (fptr,
	   "n       x         z         j       ntot    ave_freq    T_rad      W       Lum     Lum_abs\n");
  for (n = 0; n < NDIM2 - 1; n++)
    {
      nplasma = w[n].nplasma;
      fprintf (fptr,
	       "%-3d %8.3e %8.3e %8.3e %8d %8.3e %8.3e %8.3e %8.3e %8.3e\n",
	       n, w[n].x[0], w[n].x[2], plasmamain[n].j, plasmamain[n].ntot,
	       plasmamain[nplasma].ave_freq, plasmamain[nplasma].t_r,
	       plasmamain[nplasma].w, plasmamain[nplasma].lum,
	       plasmamain[nplasma].heat_tot);

    }

  fclose (fptr);

  return (0);
}








/***********************************************************
                                       Space Telescope Science Institute

 Synopsis: wind_ip() populates the plasma object ferland_ip which is intended to be an 
      estimate of the ionization parameter for that cell. It assumes all ionizaing photons
      are produces from the origin.

Arguments:		

Returns:
 
Description:	
	
Notes:
	There is not much point in calling this until you have propagated a few photons
History:
 	11Oct - NSH Coded to try and provide a 'correct' ionisation parameter for the wind,
               calculated exactly as per the ionization parameter in hazy1 (eq 5.4)

**************************************************************/


int
wind_ip ()
{
  int n;
  float r;
  for (n = 0; n < NPLASMA; n++)
    {
      r =
	sqrt ((wmain[plasmamain[n].nwind].xcen[0] *
	       wmain[plasmamain[n].nwind].xcen[0] +
	       wmain[plasmamain[n].nwind].xcen[1] *
	       wmain[plasmamain[n].nwind].xcen[1] +
	       wmain[plasmamain[n].nwind].xcen[2] *
	       wmain[plasmamain[n].nwind].xcen[2]));

      plasmamain[n].ferland_ip =
	geo.n_ioniz / (4 * PI * C * plasmamain[n].rho * rho2nh * (r * r));

//OLD     Log ("NSH log(ferland_ip) for cell %i = %e (r=%e nh=%e nioniz=%e)\n",n,log10(plasmamain[n].ferland_ip),r,plasmamain[n].rho*rho2nh,geo.n_ioniz);
      r =
	sqrt ((wmain[plasmamain[n].nwind].x[0] *
	       wmain[plasmamain[n].nwind].x[0] +
	       wmain[plasmamain[n].nwind].x[1] *
	       wmain[plasmamain[n].nwind].x[1] +
	       wmain[plasmamain[n].nwind].x[2] *
	       wmain[plasmamain[n].nwind].x[2]));

      plasmamain[n].ferland_ip =
	geo.n_ioniz / (4 * PI * C * plasmamain[n].rho * rho2nh * (r * r));

//OLD     Log  ("NSH log(ferland_ip) for cell %i = %e (r=%e nh=%e nioniz=%e)\n",n,log10(plasmamain[n].ferland_ip),r,plasmamain[n].rho*rho2nh,geo.n_ioniz);
    }
  return (0);
}
