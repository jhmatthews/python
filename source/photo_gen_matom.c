/* photo_gen_matom.c contains functions for calculating emissivities and generating photons
    during the spectral cycles. The actual functions which do the jumps inside an activated 
   macro-atom are in matom.c. This is partly done to prevent overly long files (JM1504)
*/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "atomic.h"
#include "python.h"
#include <gsl/gsl_block.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
//#include <gsl/gsl_blas.h>
#include "my_linalg.h"

/************************************************************
                                    Imperial College London
Synopsis:

       get_kpkt_f returns the specific luminosity in the band needed for the computation of the
       spectrum. It gets the total energy radiated by the process k-packet -> r-packet in the
       required wavelength range.


Arguments:   
       WindPtr w                   the ptr to the structure defining the wind

                                  
           

Returns:   The energy radiated by the process k-packet -> r-packet in the wind in the 
           wavelength range required for the specrum calculation.



Description:


Notes:  
        
         
History:
          June 04 SS - coding began
	06may	ksl	57+ -- Recoded to use plasma array

************************************************************/

double
get_kpkt_f ()
{
  int n;
  double lum;

  lum = 0.0;

  for (n = 0; n < NPLASMA; n++)
    {
      lum += plasmamain[n].kpkt_emiss;
    }


  return (lum);
}

/* All done. */

/************************************************************
                                    Imperial College London
Synopsis:

       get_matom_f returns the specific luminosity in the band needed for the computation of the
       spectrum. It gets the total energy radiated by the deactivation of macro atoms in the
       required wavelength range.


Arguments:   
       mode  variable which controls whether or not we need to compute the
             emissivities (CALCULATE_MATOM_EMISSIVITIES) or use stored ones
             because we are restarting a spectral cycle (USE_STORED_MATOM_EMISSIVITIES)
             see #define statements in python.h and code in xdefine_phot().

                                              

Returns:   The energy radiated by the deactivation of macro atoms in the wind in the 
           wavelength range required for the specrum calculation.



Description:


Notes:  
        
         
History:
          June 04 SS - coding began
          Sep  04 SS - significant modification to improve the treatment of macro
                       atoms in spectral synthesis steps. 
	06may	ksl	57+ -- Recoded to use plasma structure

************************************************************/

double
get_matom_f (mode)
    int mode;
{
  int n, m;
  int mm, ss;
  double lum;
  int level_emit[NLEVELS_MACRO], kpkt_emit;
  int n_tries, n_tries_local;
  struct photon ppp;
  double contribution, norm;
  int nres, which_out;
  int my_nmin, my_nmax;		//These variables are used even if not in parallel mode


  if (mode == USE_STORED_MATOM_EMISSIVITIES) 
    {
  	  Log("geo.pcycle (%i) > 0, so using saved level emissivities for macro atoms\n", 
  		   geo.pcycle);
    }

  else   // we need to compute the emissivities
    {
#ifdef MPI_ON
  int num_mpi_cells, num_mpi_extra, position, ndo, n_mpi, num_comm, n_mpi2;
  int size_of_commbuffer;
  char *commbuffer;

  /* the commbuffer needs to communicate 2 variables and the number of macor levels, 
     plus the variable for how many cells each thread is doing */
  size_of_commbuffer =
    8 * (3 + nlevels_macro) * (floor (NPLASMA / np_mpi_global) + 1);

  commbuffer = (char *) malloc (size_of_commbuffer * sizeof (char));
#endif



  which_out = 0;
  n_tries = 5000000;
  geo.matom_radiation = 0;
  n_tries_local = 0;



  norm = 0;
  for (n = 0; n < NPLASMA; n++)
    {
      for (m = 0; m < nlevels_macro; m++)
	{
	  norm += macromain[n].matom_abs[m];
	  if (sane_check(macromain[n].matom_abs[m]))
	  	Error("matom_abs is %8.4e in matom %i level %i\n",macromain[n].matom_abs[m], n, m);
	}
      norm += plasmamain[n].kpkt_abs;
      if (sane_check(plasmamain[n].kpkt_abs))
        Error("kpkt_abs is %8.4e in matom %i\n",plasmamain[n].kpkt_abs, n);
    }



  Log ("Calculating macro atom emissivities- this might take a while...\n");

  /* For MPI parallelisation, the following loop will be distributed over multiple tasks. 
     Note that the mynmim and mynmax variables are still used even without MPI on */
  my_nmin = 0;
  my_nmax = NPLASMA;

#ifdef MPI_ON

  num_mpi_cells = floor(NPLASMA/np_mpi_global);  // divide the cells between the threads
  num_mpi_extra = NPLASMA - (np_mpi_global*num_mpi_cells);  // the remainder from the above division
  
  /* this next loop splits the cells up between the threads. All threads with 
     rank_global<num_mpi_extra deal with one extra cell to account for the remainder */
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

  Log_parallel
    ("Thread %d is calculating macro atom emissivities for macro atoms %d to %d\n",
     rank_global, my_nmin, my_nmax);

#endif




  for (n = my_nmin; n < my_nmax; n++)
    {

      /* JM 1309 -- these lines are just log statements which track progress, as this section
         can take a long time */
#ifdef MPI_ON
      if (n % 50 == 0)
	Log
	  ("Thread %d is calculating  macro atom emissivity for macro atom %7d of %7d or %6.3f per cent\n",
	   rank_global, n, my_nmax, n * 100. / my_nmax);
#else
      if (n % 50 == 0)
	Log
	  ("Calculating macro atom emissivity for macro atom %7d of %7d or %6.3f per cent\n",
	   n, my_nmax, n * 100. / my_nmax);
#endif

      for (m = 0; m < nlevels_macro + 1; m++)
	{
	  if ((m == nlevels_macro && plasmamain[n].kpkt_abs > 0)
	      || (m < nlevels_macro && macromain[n].matom_abs[m] > 0))
	    {
	      if (m < nlevels_macro)
		{
		  if (macromain[n].matom_abs[m] > 0)
		    {
		      n_tries_local =
			(n_tries * macromain[n].matom_abs[m] / norm) + 10;
		    }
		  else
		    {
		      n_tries_local = 0;
		    }
		}
	      else if (m == nlevels_macro)
		{
		  if (plasmamain[n].kpkt_abs > 0)
		    {
		      n_tries_local =
			(n_tries * plasmamain[n].kpkt_abs / norm) + 10;
		    }
		  else
		    {
		      n_tries_local = 0;
		    }
		}
	      /* We know "matom_abs" is the amount of energy absobed by
	         each macro atom level in each cell. We now want to determine what fraction of
	         that energy re-appears in the frequency range we want and from which macro atom
	         level it is re-emitted (or if it appears via a k-packet). */

	      for (mm = 0; mm < NLEVELS_MACRO; mm++)
		{
		  level_emit[mm] = 0;
		}
	      kpkt_emit = 0;
	      if (n_tries_local > 0)
		{
		  for (ss = 0; ss < n_tries_local; ss++)
		    {
		      if (m < nlevels_macro)
			{
			  /* Dealing with excitation of a macro atom level. */
			  /* First find a suitable transition in which to excite the macro atom (to get it
			     in the correct starting level. */

			  /* As Knox pointed out in an earlier comment
			     here, the next loop isn't the best way to
			     do this - however I'm not sure we can
			     just replace nlines->nlines_macro since
			     we don't know for sure that the
			     macro_lines are the first lines in the
			     array. If this is really a bottle neck we
			     should set up a quicker way to identify
			     the level - the purpose is just to set
			     the macro atom state to the appropriate
			     upper level so all we need is some way to
			     tell it which state it should be. */

			  nres = 0;
			  while ((nres < nlines))
			    {
			      if (lin_ptr[nres]->nconfigu != m)
				{
				  nres += 1;
				}
			      else
				{
				  break;
				}
			    }


			  if (nres > (nlines - 1))
			    {
			      nres = NLINES + 1;
			      while (nres < (NLINES + 1 + nphot_total))
				{
				  if (phot_top[nres - NLINES - 1].uplev != m)
				    {
				      nres += 1;
				    }
				  else
				    {
				      break;
				    }
				}
			    }

			  if (nres > NLINES + nphot_total)
			    {
			      Error ("Problem in get_matom_f (1). Abort. \n");
			      exit (0);
			    }

			  ppp.nres = nres;
			  ppp.grid = plasmamain[n].nwind;
			  ppp.w = 0;

			  macro_gov (&ppp, &nres, 1, &which_out);


			  /* Now a macro atom has been excited and followed until an r-packet is made. Now, if that 
			     r-packet is in the correct frequency range we record it. If not, we throw it away. */
			}
		      else if (m == nlevels_macro)
			{
			  /* kpkt case. */

			  nres = -2;	//will do - just need
			  //something that will tigger kpkt

			  ppp.nres = nres;
			  ppp.grid = plasmamain[n].nwind;
			  ppp.w = 0;

			  macro_gov (&ppp, &nres, 2, &which_out);

			  /* We have an r-packet back again. */
			}


		      if (ppp.freq > em_rnge.fmin && ppp.freq < em_rnge.fmax)
			{
			  if (which_out == 1)
			    {
			      if (nres < 0)
				{
				  Error
				    ("Negative out from matom?? Abort.\n");
				  exit (0);
				}

			      /* It was a macro atom de-activation. */
			      if (nres < nlines)
				{	/* It was a line. */
				  level_emit[lin_ptr[nres]->nconfigu] += 1;
				}
			      else
				{
				  level_emit[phot_top
					     [nres - NLINES - 1].uplev] += 1;
				}
			    }
			  else if (which_out == 2)
			    {
			      /* It was a k-packet de-activation. */
			      kpkt_emit += 1;
			    }
			  else
			    {
			      Error
				("Packet didn't emerge from matom or kpkt??? Abort. \n");
			      exit (0);
			    }
			}
		    }


		  /* Now we've done all the runs for this level de-activating so we can add the contributions
		     to the level emissivities, the k-packet emissivity and the total luminosity. */

		  for (mm = 0; mm < nlevels_macro; mm++)
		    {
		      contribution = 0;
		      if (m < nlevels_macro)
			{
			  macromain[n].matom_emiss[mm] += contribution =
			    level_emit[mm] * macromain[n].matom_abs[m] /
			    n_tries_local;
			}
		      else if (m == nlevels_macro)
			{
			  macromain[n].matom_emiss[mm] += contribution =
			    level_emit[mm] * plasmamain[n].kpkt_abs /
			    n_tries_local;
			}      
		    }

		  if (m < nlevels_macro)
		    {
		      plasmamain[n].kpkt_emiss +=
			kpkt_emit * macromain[n].matom_abs[m] / n_tries_local;
		    }
		  else if (m == nlevels_macro)
		    {
		      plasmamain[n].kpkt_emiss +=
			kpkt_emit * plasmamain[n].kpkt_abs / n_tries_local;
		    }
		}
	    }
	}
    }


  /*This is the end of the update loop that is parallelised. We now need to exchange data between the tasks.
    This is done much the same way as in wind_update */
#ifdef MPI_ON
  for (n_mpi = 0; n_mpi < np_mpi_global; n_mpi++)
    {
      /* here we loop over the number of threads. If the thread is this thread then we pack the macromain information
         into memory and then broadcast it to ther other threads */
      position = 0;

      if (rank_global == n_mpi)
	{

	  Log
	    ("MPI task %d is working on matoms %d to max %d (total size %d).\n",
	     rank_global, my_nmin, my_nmax, NPLASMA);

	  MPI_Pack (&ndo, 1, MPI_INT, commbuffer,
		    size_of_commbuffer, &position, MPI_COMM_WORLD);
	  for (n = my_nmin; n < my_nmax; n++)
	    {

	      /* pack the number of the cell, and the kpkt and macro atom emissivites for that cell // */
	      MPI_Pack (&n, 1, MPI_INT, commbuffer,
			size_of_commbuffer, &position, MPI_COMM_WORLD);
	      MPI_Pack (&plasmamain[n].kpkt_emiss, 1, MPI_DOUBLE,
			commbuffer, size_of_commbuffer, &position,
			MPI_COMM_WORLD);
	      MPI_Pack (macromain[n].matom_emiss, nlevels_macro, MPI_DOUBLE,
			commbuffer, size_of_commbuffer, &position,
			MPI_COMM_WORLD);

	    }

	  Log_parallel
	    ("MPI task %d broadcasting matom emissivity information.\n",
	     rank_global);
	}


      /* Set MPI_Barriers and broadcast information to other threads */
      MPI_Barrier (MPI_COMM_WORLD);
      MPI_Bcast (commbuffer, size_of_commbuffer, MPI_PACKED,
		 n_mpi, MPI_COMM_WORLD);
      MPI_Barrier (MPI_COMM_WORLD);
      Log_parallel
	("MPI task %d survived broadcasting matom emissivity information.\n",
	 rank_global);



      position = 0;

      /* If not this thread then we unpack the macromain information from the other threads */

      if (rank_global != n_mpi)
	{
	  MPI_Unpack (commbuffer, size_of_commbuffer, &position,
		      &num_comm, 1, MPI_INT, MPI_COMM_WORLD);
	  for (n_mpi2 = 0; n_mpi2 < num_comm; n_mpi2++)
	    {

	      /* unpack the number of the cell, and the kpkt and macro atom emissivites for that cell */
	      MPI_Unpack (commbuffer, size_of_commbuffer,
			  &position, &n, 1, MPI_INT, MPI_COMM_WORLD);
	      MPI_Unpack (commbuffer, size_of_commbuffer, &position,
			  &plasmamain[n].kpkt_emiss, 1, MPI_DOUBLE,
			  MPI_COMM_WORLD);
	      MPI_Unpack (commbuffer, size_of_commbuffer, &position,
			  macromain[n].matom_emiss, nlevels_macro, MPI_DOUBLE,
			  MPI_COMM_WORLD);
	    }
	}
    }  // end of parallelised section
  #endif	

  } // end of if loop which controls whether to compute the emissivities or not 

  /* this next loop just calculates lum to be the correct summed value in parallel mode */
  /* if mode == USE_STORED_MATOM_EMISSIVITIES this is all this routine does */
  lum = 0.0;   // need to rezero, fixes segfault bug #59 

  for (n = 0; n < NPLASMA; n++)
    {

      for (mm = 0; mm < nlevels_macro; mm++)
	{
	  lum += macromain[n].matom_emiss[mm];
	}
    }

  geo.matom_radiation = 1;

  return (lum);
}

/* All done. */



/************************************************************
                                    Imperial College London
Synopsis:

     photo_gen_kpkt produces photon packets to account for creating of r-packets
     by k-packets in the spectrum calculation. It should only be used once the total 
     energy emitted in this way in the wavelength range in question is well known
     (calculated in the ionization cycles).


Arguments:   
       WindPtr w                   the ptr to the structure defining the wind
       PhotPtr p;                  the ptr to the structire for the photons
       double weight;              the photon weight
       int photstart, nphot;       the number of the first photon to be generated and
                                   the total number of photons to be generated by this
                                   routine

                                           

Returns:   
      When it finishes it should have generated nphot photons from k-packet elliminations.


Description:
      This routine is closely related to photo_gen_wind from which much of the code has been copied.

Notes:  
        
         
History:
          June 04 SS - coding began
	  04aug	ksl	Modified so that uses coordinate system
	  		independent routine get_random_location
			to set location of generated kpkt
          Sep  04 SS - minor bug fixed
	06may	ksl	57+ -- Updated to use plasma structure. 
			and to elimate passing entrie wind structure

************************************************************/
int
photo_gen_kpkt (p, weight, photstart, nphot)
     PhotPtr p;
     double weight;
     int photstart, nphot;
{
  int photstop;
  int icell;
  double xlum, xlumsum;
  struct photon pp;
  int nres, esc_ptr;
  int n;
  double v[3];
  double dot ();
  int wind_n_to_ij (), stuff_v (), randvec ();
  int get_random_location ();
  double test;
  int nnscat;
  double dvwind_ds (), sobolev ();
  int nplasma;



  photstop = photstart + nphot;
  Log ("photo_gen_kpkt creates nphot %5d photons from %5d to %5d \n", nphot,
       photstart, photstop);

  for (n = photstart; n < photstop; n++)
    {
      /* locate the wind_cell in which the photon bundle originates. */

      xlum = (rand () + 0.5) / (MAXRAND) * geo.f_kpkt;

      xlumsum = 0;
      icell = 0;
      while (xlumsum < xlum)
	{
	  if (wmain[icell].vol > 0.0)
	    {
	      nplasma = wmain[icell].nplasma;
	      xlumsum += plasmamain[nplasma].kpkt_emiss;
	    }
	  icell++;
	}
      icell--;			/* This is the cell in which the photon must be generated */

      /* Now generate a single photon in this cell */
      p[n].w = weight;

      pp.freq = 0.0;
      pp.grid = icell;

      /* This following block is a bad way of doing it - kpkt could be modified to 
         do what we want in a more elegant way. */

      test = pp.freq;

      while (test > em_rnge.fmax || test < em_rnge.fmin)
	{
	  kpkt (&pp, &nres, &esc_ptr);
	  if (esc_ptr == 0)
	    {
	      test = 0.0;
	    }
	  else
	    {
	      test = pp.freq;
	    }
	}


      p[n].freq = pp.freq;
      p[n].nres = nres;

      /* The photon frequency is now known. */

      /* Determine the position of the photon in the moving frame */

      /* ! Need to account for possibility that photon is in the torus */
      get_random_location (icell, 0, p[n].x);

      p[n].grid = icell;

      nnscat = 1;
      // Determine the direction of the photon
      // Need to allow for anisotropic emission here
      if (p[n].nres < 0 || p[n].nres > NLINES || geo.scatter_mode == 0)
	{
	  /*  It was either an electron scatter so the  distribution is isotropic, or it
	     was a resonant scatter but we want isotropic scattering anyway.  */
	  randvec (p[n].lmn, 1.0);	/* The photon is emitted isotropically */
	}
      else if (geo.scatter_mode == 1)
	{			// It was a line photon and we want anisotropic scattering

	  // -1. forces a full reinitialization of the pdf for anisotropic scattering

	  randwind (&p[n], p[n].lmn, wmain[icell].lmn);

	}
      else if (geo.scatter_mode == 2)
	{			//It was a line photon and we want the thermal trapping anisotropic model

	  randwind_thermal_trapping(&p[n], &nnscat);
	}

      p[n].nnscat = nnscat;





      /* The next two lines correct the frequency to first order, but do not result in
         forward scattering of the distribution */

      vwind_xyz (&p[n], v);
      p[n].freq /= (1. - dot (v, p[n].lmn) / C);

      p[n].istat = 0;
      p[n].tau = p[n].nscat = p[n].nrscat = 0;
      p[n].origin = PTYPE_WIND;	// Call it a wind photon

    }


  return (nphot);		/* Return the number of photons generated */


  /* All done. */
}

/************************************************************
                                    Imperial College London
Synopsis:

     photo_gen_matom produces photon packets to account for creation of r-packets
     by deactivation of macro atoms in the spectrum calculation. It should only be used 
     once the total energy emitted in this way in the wavelength range in question is well known
     (calculated in the ionization cycles).


Arguments:   
       WindPtr w                   the ptr to the structure defining the wind
       PhotPtr p;                  the ptr to the structire for the photons
       double weight;              the photon weight
       int photstart, nphot;       the number of the first photon to be generated and
                                   the total number of photons to be generated by this
                                   routine

                                           

Returns:   
      When it finishes it should have generated nphot photons from macro atom deactivations.


Description:
      This routine is closely related to photo_gen_kpkt from which much of the code has been copied.

Notes:  
        
         
History:
          June 04 SS - coding began
          Aug  04 SS - modified as for gen_kpkt above by KSL
          Sep  04 SS - minor bug fixed
	06may	ksl	57+ -- Initial adaptation to plasma structure
	0812	ksl	67c -- Changed call to photo_gen_matom to
			eliminate WindPtr w, since we have access
			to this through wmain.

************************************************************/
int
photo_gen_matom (p, weight, photstart, nphot)
     PhotPtr p;
     double weight;
     int photstart, nphot;
{
  int photstop;
  int icell;
  double xlum, xlumsum;
  struct photon pp;
  int nres;
  int n;
  double v[3];
  double dot ();
  int wind_n_to_ij (), stuff_v (), randvec ();
  int emit_matom ();
  double test;
  int upper;
  int nnscat;
  double dvwind_ds (), sobolev ();
  int nplasma;



  photstop = photstart + nphot;
  Log ("photo_gen_matom creates nphot %5d photons from %5d to %5d \n", nphot,
       photstart, photstop);

  for (n = photstart; n < photstop; n++)
    {
      /* locate the wind_cell in which the photon bundle originates. And also decide which of the macro
         atom levels will be sampled (identify that level as "upper"). */

      xlum = (rand () + 0.5) / (MAXRAND) * geo.f_matom;

      xlumsum = 0;
      icell = 0;
      upper = 0;
      while (xlumsum < xlum)
	{
	  if (wmain[icell].vol > 0.0)
	    {
	      nplasma = wmain[icell].nplasma;
	      xlumsum += macromain[nplasma].matom_emiss[upper];
	      upper++;
	      if (upper == nlevels_macro)
		{
		  upper = 0;
		  icell++;
		}
	    }
	  else
	    {
	      icell++;
	      upper = 0;
	    }
	}

      if (upper == 0)
	{
	  /* If upper = 0 at this point it means the loop above finished with an increment of
	     icell and setting upper = 0. In such a case the process we want was the LAST process
	     in the previous cell. */
	  icell--;
	  upper = nlevels_macro;
	}

      upper--;			/* This leaves the macro atom level that deactivaties. */

      /* Now generate a single photon in this cell */
      p[n].w = weight;

      pp.freq = 0.0;
      pp.grid = icell;

      /* This following block is a bad way of doing it but it'll do as
         a quick and dirty test for now. Really kpkt should be modified to 
         do what we want in a more elegant way. */

      test = pp.freq;

      while (test > em_rnge.fmax || test < em_rnge.fmin)
	{

	  /* Call routine that will select an emission process for the
	     deactivating macro atom. If is deactivates outside the frequency
	     range of interest then ignore it and try again. SS June 04. */

	  emit_matom (wmain, &pp, &nres, upper);

	  test = pp.freq;
	}


      p[n].freq = pp.freq;
      p[n].nres = nres;


      /* The photon frequency is now known. */

      /* Determine the position of the photon in the moving frame */

      /* !! ERROR - need to account for posibilyt that photon is in torus */

      get_random_location (icell, 0, p[n].x);

      p[n].grid = icell;


      // Determine the direction of the photon
      // Need to allow for anisotropic emission here
      nnscat = 1;
      if (p[n].nres < 0 || p[n].nres > NLINES || geo.scatter_mode == 0)
	{
	  /*  It was either an electron scatter so the  distribution is isotropic, or it
	     was a resonant scatter but we want isotropic scattering anyway.  */
	  randvec (p[n].lmn, 1.0);	/* The photon is emitted isotropically */
	}
      else if (geo.scatter_mode == 1)
	{			// It was a line photon and we want anisotropic scattering

	  // -1. forces a full reinitialization of the pdf for anisotropic scattering

	  randwind (&p[n], p[n].lmn, wmain[icell].lmn);

	}
      else if (geo.scatter_mode == 2)
	{			//It was a line photon and we want the thermal trapping anisotropic model

      randwind_thermal_trapping(&p[n], &nnscat);
	}
      p[n].nnscat = nnscat;


      /* The next two lines correct the frequency to first order, but do not result in
         forward scattering of the distribution */

      vwind_xyz (&p[n], v);
      p[n].freq /= (1. - dot (v, p[n].lmn) / C);

      p[n].istat = 0;
      p[n].tau = p[n].nscat = p[n].nrscat = 0;
      p[n].origin = PTYPE_WIND;	// Call it a wind photon

    }


  return (nphot);		/* Return the number of photons generated */


  /* All done. */
}