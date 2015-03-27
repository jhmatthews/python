/***********************************************************
                                       Space Telescope Science Institute

 Synopsis: balance is a test program to allow one to understand how and whether various
     radiative process approach lte
   
 Arguments:		

Returns:
 
 
Description:	

	G:   An option to calculate the abundances of the light elements in a variety of situations, for comparison
	     with other routines like Cloudy
		
Notes:
	Insofar as possible the program uses python routines directly

History:
	98jul	ksl	Coding began during a visit to Oxford
	01mar	ksl	Verified program still runs basically, but note
			that a number of the underlying routines have
			changed.
	01sept	ksl	Attempted to really sort out some of the various options
			that no longer mattered and eliminate these sections
			of balance, and added options for levels.
	02jun	ksl	Added an option to set the total heating
	02jun	ksl	Added a subroutine to calculate what the derived values
			for the radiation field would be after photon passage.
	06nov	ksl	58b -- Began work to upgrade this so it will work
			again.  Note that we have not tried to keep this program
			working for some time.
	07jul	ksl	58e -- It is unclear where I ended up with this, but
			program does compile and seem to run on OSX, except
			for VERY_BIG warning.
	10nov	ksl	69 -- Fixed several small problems and added a good bit
			of documentation.  The outputs look plausible.
	12jun	ksl	72 - Restored to "operational" condition.  There
			is sill work to do to see whather it is actually
			producing what we want
 
**************************************************************/


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "atomic.h"
#include "python.h"
#include "balance_templates.h"





int
main (argc, argv)
     int argc;
     char *argv[];

{
  FILE *fopen (), *fptr;
  WindPtr w;
  PlasmaPtr wplasma;
  PhotPtr p;
  double t_r, t_e;
  double alpha_agn, lum_agn, r_agn, d_agn;
  double nh, heat;
  double weight, vol, vmax;
  char choice;
  char atomic_filename[LINELENGTH];
  int lmode;
  int freq_sampling;
  int check_field ();
  int do_detail ();
  int cmode;
  int radmode;			/* 1=blackbody, 2=powerlaw from agn */
  int i;
  int nelement, nele, z, nn, first, last, m;
  double tot, agn_ip;
  NPHOT = 100000;
  DENSITY_PHOT_MIN = 1.0e-10;	//This turns on photoionization calculations during heating and cooling.
  geo.ioniz_or_extract = 1;	// Set the overall pointer to tell the code we are in an ionization cycle.

/* These next lines define the grid size and give the push through distance */

  ndim = 3;
  mdim = 3;
  NDIM = ndim;
  MDIM = mdim;
  NDIM2 = ndim * mdim;
  NPLASMA = NDIM2;
  dfudge = 0.0001;
  DFUDGE = dfudge;

/* End of definition of wind array size */

/* Get the atomic data */

  strcpy (atomic_filename, "atomic/standard73");
  rdstr ("Atomic_data", atomic_filename);
  get_atomic_data (atomic_filename);
  strcpy (geo.atomic_filename, atomic_filename);

/* Create the wind structure */
  wmain = w = (WindPtr) calloc (sizeof (wind_dummy), NDIM2);
  plasmamain = wplasma = (PlasmaPtr) calloc (sizeof (plasma_dummy), NPLASMA);
//      printf ("The first random number is %i",rand());
  for (i = 0; i < NDIM2; i++)
    {
      w[i].nplasma = wplasma[i].nwind = i;
    }

  if (w == NULL)
    {
      Error
	("There is a problem in allocating memory for the wind structure\n");
      exit (0);
    }

/* 
 Note above that we have created two variables that do the same thing, e.g wmain and w, 
 and plasmamain and wplasam.  The difference is that wamain and plasmamain are held in 
 common, wherease the other two are not
*/

/* Generate the photon structure */

  p = (PhotPtr) calloc (sizeof (p_dummy), NPHOT);
  if (p == NULL)
    {
      Log ("There was a problem allocating p\n");
      exit (0);
    }

/* Create a set of inital choices */

  geo.line_mode = 2;
  geo.ioniz_mode = 0;
  geo.rt_mode = 0;
  geo.pl_nh = nh = 1e10;
  geo.pl_t_r = t_r = 20000;
  geo.pl_t_e = t_e = 20000;
  geo.pl_w = weight = 1;
  geo.pl_vol = vol = 1.0;
  geo.pl_vmax = vmax = 1e6;
  geo.lum_agn = lum_agn = 1e43;
  geo.alpha_agn = alpha_agn = -1.5;
  geo.r_agn = r_agn = 1e14;
  geo.d_agn = d_agn = 1e16;
  freq_sampling = 1;
  radmode = 1;

  choice = 'l';

/* Initialize a plane parallel wind */

  pl_wind_define (w);

/* Initialize the occupation number of the levels in the wind to LTE using t_e */
  /* 080802 - ksl - It's unclear to me why call to levels is needed here, as I would have thought
   * it would have been called by pl_wind_define, which does the ionzation balance
   * 101104 - ksl - While that may be logincal, as far as I can tell pl_wind_define does
   * not actually call partion_fucntion, even though it does set up the ionization levels.
   * It might be useful to test this with ddd, but that's beyond the scope of what I am
   * doing now.  It should not be harmful to call this again here in any event.
   * Note that, in python levels is called from partition function.
   */


  levels (wplasma, 1);
  print_levels (wplasma);

/* Completed the initialization */

a:				/* This is the return point.  One could replace all this whith some kind of while statement */

  /* Print out the current values of various variables */

  Log
    ("\nnh %8.2g t_r %8.2g t_e %8.2g weight %8.2e vol %8.2g vmax %8.2g line_mode %2d ioniz_mode %d rt_mode %d freq_sampling %d\n\n",
     nh, t_r, t_e, weight, vol, vmax, geo.line_mode, geo.ioniz_mode,
     geo.rt_mode, freq_sampling);

  /* Print out the choices that the user can make to change various variables or to cause a calculation to be made */

  Log
    ("h=Nh, r=t_r, e=t_e, w=weight, v=vol, s=vmax, p=line_mode, k=rad_trans_mode\n");
  Log ("i=ionization mode, j=calculate heating/cooling and abundances\n");
  Log ("l=lte(force t_e=t_r w=1),  a=absolute_lte(ions+levels), A=levels\n");
  Log
    ("b=frequency sampling d=dumb, f=find_te, ,g=go  H=Set heating, D -- modify to detailed rates\n");
  Log ("B=blackbody spectrum(default), P=powerlaw spectrum\n");
  Log
    ("Other: a: Force t_r=t_e and LTE level populations, and calculate saha abundances\n");
  Log
    ("Other: A: Calculate levels with various options.  Just sets lmode and recaclates partition func and levels\n");
  Log ("Other  b: Set freqqence sampling frequency_sampling(0=uniform)\n");
  Log ("Other: B: Setting radiation mode to blackbody");
  Log
    ("Other:  d: Increment t_e up or down to attempt to match the current value of the heating.\n");
  Log
    ("Other:  D: Modify the current abundances to produce detail balance, and then recalculate the heating and cooling without modifying the abundances\n");
  Log
    ("Other:  f: Adjust t_e to balance heating and cooling without changing the abundances\n");
  Log
    ("Other   g: Adjust t_e and abundances using python's one_shot routine.\n");
  Log
    ("Other   G: Calculate ionization fractions as a function of T and W. Remake figure in LK02\n");
  Log ("Other   H: Set the total heating\n");
  Log ("Other   I: Run multi_concen\n");
  Log
    ("Other   j: Calculate abundances and then heating and cooling with current ionization mode, etc. Do not update afterwards\n");
  Log
    ("Other   l:Force Saha ion abundances and t_e=t_r and weight to 1 and then calculate heating and cooling. (uses python routines)\n");
  Log
    ("Other   o:  output the ion fractions of a particular element for offline investigation\n");
  Log ("Other   p: Set the line transfer mode mode\n");
  Log ("Other  P: get the parameters for a power law spectrum \n");
  Log
    ("Other  S:Carry out one calculation of the sim parameter and apply to current densities\n");
  Log ("Other  T: This is just for testing of routines\n");






  rdchar ("Choice", &choice);
  Log_silent ("Choice %c\n", choice);

/* Initialize the wind structure to match the input variables */

  wplasma[0].w = weight;
  wplasma[0].t_e = t_e;
  wplasma[0].t_r = t_r;
  w[0].vol = wplasma[0].vol = vol;
//      w[0].dvds_ave=1.0;
  wplasma[0].rho = nh / rho2nh;
  wplasma[0].gain = 0.5;

/* Now exeute the choices */

  switch (choice)
    {


    case 'a':			/* Force Saha ion densities, and LTE level populations, and set t_e=t_r and weight to 1
				   This differs from the case l (lte option) in that LTE level populations are forced explicitly */

      rddoub ("LTE temp", &t_e);
      Log_silent ("LTE temp %8.2g\n", t_e);

      geo.pl_nh = nh;
      geo.pl_t_r = t_r = t_e;
      geo.pl_t_e = t_e;
      geo.pl_w = weight = 1.;
      geo.pl_vol = vol;
      geo.pl_vmax = vmax;
      pl_wind_define (w);

      Log
	("Calculating heating and cooling with t_e=t_r, LTE ionization & LTE exited state populations\n");
      absolute_saha (&wplasma[0], p, nh, t_e, t_e, weight, freq_sampling);
      break;

    case 'A':
/* Calculate the levels with various options */

      lmode = 1;

      rdint ("Level_calc(0=on.the.spot,1=lte)", &lmode);

      partition_functions (wplasma, lmode);
      levels (wplasma, lmode);
      print_levels (wplasma);
      break;


    case 'b':			// freqmode
      rdint ("frequency_sampling(0=uniform)", &freq_sampling);
      if (freq_sampling == 0)
	Log ("OK -- using normal sampling with no required photon numbers\n");
      else
	{
	  freq_sampling = 1;
	  Log
	    ("OK-- setting energy bands to assure substantial photons above ionization edges\n");
	}
      break;

    case 'B':
      Log ("Setting radiation mode to blackbody");
      radmode = 1;

      break;


    case 'd':
/* 
Increment t_e up or down to attempt to match the current value of the heating.  
It's very dumb though because the heating is not automatically calculated.
*/

      dumb_step (&w[0], &t_e);
      break;

    case 'D':
/* Modify the current abundances to produce detail balance, and then recalculate
the heating and cooling without modifying the abundances
*/
      do_detail (&wplasma[0]);
      xsaha (&w[0], p, nh, t_r, t_e, weight, -1, freq_sampling);
      break;
    case 'e':			// Set the electron temperature

      rddoub ("t_e", &t_e);
      Log_silent ("t_e %8.2g\n", t_e);
      break;

    case 'f':

/*Adjust t_e to balance heating and cooling without changing the abundances */

      Log
	("Adjusting t_e to balance heating and cooling without changing the abundances\n");
      find_te (&w[0]);
      t_e = wplasma[0].t_e;
      break;

    case 'g':			// Adjust t_e and abundances using python's one_shot routine.  

/* The first time this routine is called after t_r or nh have been called a saha density
distribution is set up based on concentrations, and then one_shot is called.  If this
case is repeatedly executed it should eventually match the abundances and temperatures.
Note that the ionization mode is not changed by this routine.  */

      if (geo.ioniz_mode == 0)
	cmode = 3;
      if (geo.ioniz_mode == 4)
	cmode = 4;
      if (geo.ioniz_mode == 5)
	cmode = 5;
      printf ("wplasma[0].t_r=%f tr=%f\n", wplasma[0].t_r, t_r);
      cycle (&wplasma[0], p, nh, t_r, t_e, weight, cmode, freq_sampling,
	     radmode);
      t_e = wplasma[0].t_e;
      break;


    case 'G':
/* 
Calculate ionization fractions as a function of T and W.  This is the option that
was used to calculate ionization fractions in the figure which is in
Long and Knigge 2002
 */

      if (geo.ioniz_mode == 0)
	cmode = 3;
      if (geo.ioniz_mode == 3)
	cmode = 3;		/*NSH 120823 - fixed an oddity that if you *told* balance to run in mode 3, it broke! */
      if (geo.ioniz_mode == 4)
	cmode = 4;
      if (geo.ioniz_mode == 5)
	cmode = 5;
      if (geo.ioniz_mode == 6)
	cmode = 6;
      if (geo.ioniz_mode == 7)
	cmode = 7;
      multicycle (&wplasma[0], p, cmode, freq_sampling, radmode);
      break;

    case 'h':			/*Set the hydrogen density */

      rddoub ("NH", &nh);
      Log_silent ("NH %8.2g\n", nh);
      geo.pl_nh = nh;
      break;

    case 'H':			/*Set the total heating */
      rddoub ("Total_heating", &heat);
      wplasma[0].heat_tot = heat;
      break;

    case 'i':			/*Set the mode for which ionization abundances will be calculated, e.g. saha abundances */

      printf
	("0=on_the_spot,1=saha,2=fix to std vals, 3=_on_the_spot (with one t update), 4=detailed, 10=same as 3, but ground mult.\n");
      rdint ("ioniz_mode", &geo.ioniz_mode);
      geo.partition_mode = -1;

      if (geo.ioniz_mode == 10)	// Use recalc but only ground state multiplicities
	{
	  geo.ioniz_mode = 3;
	  geo.partition_mode = 0;
	}
      printf ("Ionization mode %d", geo.ioniz_mode);
      Log_silent ("Ionization mode %d", geo.ioniz_mode);
      break;

    case 'I':
      multi_concen (&wplasma[0]);
      break;

    case 'j':

/* Calculate abundances and then heating and cooling with current ionization mode, etc.
but do not update basic parameters aftewards.  Multiple calls to j should give identical
answers, except for statistical changes */

      geo.pl_nh = nh;
      geo.pl_t_r = t_r;
      geo.pl_t_e = t_e;
      geo.pl_w = weight;
      geo.pl_vol = vol;
      geo.pl_vmax = vmax;
      pl_wind_define (w);

/* Assure that www is initialized before going into xsaha */
      wplasma[0].t_e = t_e;
      wplasma[0].t_r = t_r;
      wplasma[0].w = weight;



/* Currently (Sept 01), python takes ioniz_mode to mean

0  on-the-spot approximation using existing t_e
1  LTE
2  Hardwired 
3  On the spot, with one_shot at updating t_e before calculating densities

*/
      xsaha (&w[0], p, nh, t_r, t_e, weight, geo.ioniz_mode, freq_sampling);


      check_field (&wplasma[0]);

      break;

    case 'k':			/*Set the radiative tranfer mode.  USED ONLY BY BALANCE */

/* Simple here means that the line shape is set to a simple square wave
   MC means we have run the MC code and so lines are only absorbed as closely to
	the way they are handled in python as possible, i.e. only if tau exceeds a
	defined value
   Simple sobolev calculates the heating using the Sobolev approx, but without the 
   MC aspects of deciding whether any part of a photon was absorbed
*/
      rdint ("Radiative transfer(0=simple,1=MC,2=simple_sobolev)",
	     &geo.rt_mode);
      Log_silent
	("Radiative transfer(0=simple,1=MC,2=simple_sobobolev): %d\n",
	 geo.rt_mode);
      break;

    case 'l':

/* Force Saha ion abundances and t_e=t_r and weight to 1 and then calculate heating and cooling.
In this option, the levels are calculated with the same routines as in python.  The results should be identical
to case a unless there is an error somewhere */

      rddoub ("LTE temp", &t_e);
      Log_silent ("LTE temp %8.2g\n", t_e);

      geo.pl_nh = nh;
      geo.pl_t_r = t_r = t_e;
      geo.pl_t_e = t_e;
      geo.pl_w = weight = 1.0;
      geo.pl_vol = vol;
      geo.pl_vmax = vmax;
      pl_wind_define (w);
//              w[0].dvds_ave=1.0;
      geo.ioniz_mode = 1;


      Log
	("Calculating heating and cooling with t_e=t_r and Saha=LTE ionization\n");
/* Assure that www is initialized before going into xsaha */
      wplasma[0].t_e = t_e;
      wplasma[0].t_r = t_e;
      wplasma[0].w = weight;

      xsaha (&w[0], p, nh, t_e, t_e, weight, 1, freq_sampling);
      break;

    case 'o':			/* output the ion fractions of a particular element for offline investigation */

      rdint ("Element atomic number?", &nele);	// what element are we interested in
      z = nele;
      nelement = 999;		//set nelement to crazy value to check if we have an element in the database
      for (nn = 0; nn < nelements; nn++)	//loop to see if we know about this element
	{
	  if (ele[nn].z == z)	// does the atomic number of element nn in out database match what we are looking for?
	    nelement = nn;	// if so, set nelement to nn
	}
      if (nelement == 999)	//if. after loop nelement has not been reset, we dont have data for the requested element
	{
	  printf ("Your element is not in the database");	//console user with caring words
	  break;		// give up
	}
      first = ele[nelement].firstion;	//link into array for first ion of this element
      last = first + ele[nelement].nions;	//last array element will be first plus 
      printf ("%-5s ", ele[nelement].name);
      for (m = first; m < last; m++)
	tot += wplasma[0].density[m];
      for (m = first; m < last; m++)
	{
	  printf ("%8.2e ", (wplasma[0].density[m]));
	}
      printf ("\n");


      break;

    case 'O':
      fptr = fopen ("bal_summary.out", "w");
      fprintf (fptr, "Mod   %.1f %.0f %.3f %8.0f \n", log10 (nh),
	       t_r, agn_ip, t_e);

      for (nn = 0; nn < nelements; nn++)
	{
	  first = ele[nn].firstion;
	  last = first + ele[nn].nions;
	  fprintf (fptr, "%-5s ", ele[nn].name);
	  Log ("%-5s ", ele[nn].name);
	  tot = 0;
	  for (m = first; m < last; m++)
	    tot += wplasma[0].density[m];
	  for (m = first; m < last; m++)
	    {
	      fprintf (fptr, " %8.2e", wplasma[0].density[m]);
	      Log (" %8.2e", wplasma[0].density[m]);
	    }
	  fprintf (fptr, "\n");
	  Log ("\n");
	}
      fclose (fptr);


      break;



    case 'p':			/* Set the line transfer mode mode */

      rdint ("line_mode", &geo.line_mode);
      Log_silent ("line_mode %d\n", geo.line_mode);
      break;

    case 'P':			/* get the parameters for a power law spectrum */
      rddoub ("Power law luminosity between 2 and 10 KeV (ergs/s)", &lum_agn);
      rddoub ("Power law alpha - if you want negative, put in negative",
	      &alpha_agn);
      geo.lum_agn = lum_agn;
      geo.alpha_agn = alpha_agn;
      geo.const_agn =
	geo.lum_agn /
	(((pow (2.42e18, geo.alpha_agn + 1.)) -
	  pow (4.84e17, geo.alpha_agn + 1.0)) / (geo.alpha_agn + 1.0));
      Log ("Input parameters give a power law constant of %e\n",
	   geo.const_agn);
      rddoub
	("Radius of power law region - used to calculate the geometric dilution factor",
	 &r_agn);
      rddoub ("Distance from cell to power law region", &d_agn);
      geo.lum_agn = lum_agn;
      geo.alpha_agn = alpha_agn;
      geo.r_agn = r_agn;
      geo.d_agn = d_agn;



      wplasma[0].w = weight =
	0.5 * (1 -
	       sqrt (1 -
		     ((geo.r_agn * geo.r_agn) / (geo.d_agn * geo.d_agn))));
      for (i = 0; i < NXBANDS; i++)
	{
	  wplasma[0].pl_alpha[i] = alpha_agn;	//Set the alpha for this cell, in this case it is just the global version
	  wplasma[0].pl_w[i] = geo.const_agn / ((4. * PI) * (4. * PI * geo.r_agn * geo.r_agn));	//Set the pl w factor for this cell alone. 
	  i = i + 1;
	}
      Log ("Input parameters give a BB weight of %e and a SIM weight of %e\n",
	   wplasma[0].w, wplasma[0].pl_w);
      agn_ip =
	geo.const_agn *
	(((pow (50000 / HEV, geo.alpha_agn + 1.0)) -
	  pow (100 / HEV, geo.alpha_agn + 1.0)) / (geo.alpha_agn + 1.0));
      agn_ip /= (d_agn * d_agn);
      agn_ip /= nh;
      Log ("Input parameters give an Ionisation Parameter of %e\n", agn_ip);
      radmode = 2;		/* set the mode to use power law radiation */
      geo.ioniz_mode = 5;	//if we are here we can assume that we probably are interested in power laws, lets help out the user...
      rdint
	("What ionization mode should we use? 1=LTE, 3=Modified OTS(Lucy), 5=Modified LTE(Sim)",
	 &geo.ioniz_mode);


      break;

    case 'r':			/* Set the radiation temperature */

      rddoub ("t_r", &t_r);
      Log_silent ("t_r %8.2g\n", t_r);
      break;

    case 's':			/* Set the maximum velocity across the cell */

      rddoub ("vmax", &vmax);
      geo.pl_nh = nh;
      geo.pl_t_r = t_r;
      geo.pl_t_e = t_e;
      geo.pl_w = weight;
      geo.pl_vol = vol;
      geo.pl_vmax = vmax;
      pl_wind_define (w);
      Log_silent ("vmax %8.2g\n", vmax);
      break;


    case 'S':			/* Carry out one calculation of the sim parameter and apply to current densities */
      Log ("Calling sim correction with T_e=%e, r=%e, alpha=%e, and C=%e\n",
	   wplasma[0].t_e, wplasma[0].w, geo.alpha_agn, geo.const_agn);
      sim_driver (&wplasma[0]);

      break;


    case 'T':			//This is just for testing of routines

      //OLD init_bands (t_e, 0.0, 1e50, 1, &xband);
      //Old 120626 bands_init (t_e, 0.0, 1e50, 1, &xband);
      bands_init (t_e, &xband);
      init_freebound (1e3, 1e6, xband.f1[0], xband.f2[0]);

      break;


    case 'v':			/* Set the volume of the cell */

      rddoub ("volume", &vol);
      w[0].vol = wplasma[0].vol;
      geo.pl_nh = nh;
      geo.pl_t_r = t_r;
      geo.pl_t_e = t_e;
      geo.pl_w = weight;
      geo.pl_vol = vol;
      geo.pl_vmax = vmax;
      pl_wind_define (w);
      Log_silent ("vol %8.2g\n", vol);
      break;

    case 'w':			/* Set the radiative weight of the spectrum */

      rddoub ("w", &weight);
      Log_silent ("weight %8.2g\n", weight);
      break;


    case 'q':			/* Wrap up and quit */
      goto b;
      break;

    default:			/* Oops this was not a choice that did anything */
      goto a;
      break;

    }

  goto a;

b:
  fb_save ("recomb.save");
  exit (0);

}


/***********************************************************
                                       Space Telescope Science Institute

 Synopsis:  multicycle calculates abundances of the light elements in a variety of situations for comparison
            with other routines like Cloudy
   
 Arguments:		

Returns:
 
 
Description:	

	Results of this routine always appear in bal_summary.out
		
Notes:
	Most of the options are hardcoded.  At some level we should be reluctant to change this
	since it provides a trail back to Long and Knigge 2002

	In this routine, www is a PlasmaPtr

History:
	0810	ksl	67- Relooked at routine in effort to 
			begin studying ionization balance again
 
**************************************************************/

int
multicycle (www, p, mode, freq_sampling, radmode)
     PlasmaPtr www;
     PhotPtr p;
     int mode;
     int freq_sampling;
     int radmode;

{
  FILE *fopen (), *fptr, *fudgefile, *numfile, *denomfile, *heatcoolfile;
  double nh;
  double t_r;
  double t_e;
  double weight;
  double tot;
  double agn_ip, lum_agn_50, logip;

  double t_rmin, t_rmax, t_rdelt;
  double lt_rmin, lt_rmax, lt_rdelt;
  double t_rtemp;
  double ww, w_min, w_max, w_delt;

  int n, m, first, last, nn;
  int i;

  printf ("Ionization mode in multicycle %d", mode);

  /* cycle in AGN  we vary lum_agn to replicate the ionisation parameters in Sim (2010) */



  /* Finished with inputs */

  fptr = fopen ("bal_summary.out", "w");
  fudgefile = fopen ("fudge_summary.out", "w+");
  numfile = fopen ("num_summary.out", "w+");
  denomfile = fopen ("denom_summary.out", "w+");
  heatcoolfile = fopen ("heat_cool.out", "w+");
  if (radmode == 1)
    {
      /* Set hardcoded limits */
      nh = 1e12;
      t_r = 30000;
      t_e = 20000;
      weight = 0.1;

      t_rmin = 2000;
      t_rmax = 100000;
      t_rdelt = 2000;

      lt_rmin = log10 (1000.0) * 100;
      lt_rmax = log10 (10000000.0) * 100;
      lt_rdelt = 1;



      w_max = -2;
      w_min = -2;
      w_delt = 0.25;
      for (ww = w_max; ww >= w_min; ww -= w_delt)	// Note w starts at 1
	{
	  weight = pow (10, ww);
//                      for (t_r = t_rmin; t_r <= t_rmax; t_r += t_rdelt)
	  for (t_rtemp = lt_rmin; t_rtemp <= lt_rmax; t_rtemp += lt_rdelt)
	    {
	      t_r = pow (10, t_rtemp / 100.0);
	      www->gain = 0.5;
	      www->dt_e_old = 0.0;
	      www->dt_e = 0.0;
	      www->t_e = www->t_e_old = 0.9 * t_r;	//lucy guess


	      www->w = weight;
	      www->t_r = t_r;
	      printf ("NSH We are going into cycle, t_r=%f\n", www->t_r);
	      t_e = www->t_e;


	      /* Find the ionization fractions and electron temperature */
	      for (n = 0; n < 20; n++)
		{
		  Log ("multicycle %d\n", n);
		  www->rho = nh / rho2nh;
//                              printf ("wplasma[0].tr=%f tr=%f\n",www->t_r,t_r);
//                              printf ("wplasma[0].w=%f w=%f\n",www->w,weight);
		  cycle (www, p, nh, t_r, t_e, weight, mode, freq_sampling,
			 radmode);
		  t_e = www->t_e;
		  Log ("Mo  %2d %.1f %.0f %.1f %8.0f %5.2g\n", n, log10 (nh),
		       t_r, (ww), t_e, www->heat_tot / www->lum_rad);
		}

	      /* Write out the results */
	      fprintf (fptr, "Mod   %.1f %.0f %.1f %8.0f %5.2f \n",
		       log10 (nh), t_r, (ww), t_e,
		       www->heat_tot / www->lum_rad);
	      Log ("Mod    %.1f %.0f %.1f %8.0f %5.2g\n", log10 (nh), t_r,
		   (ww), t_e, www->heat_tot / www->lum_rad);
	      for (nn = 0; nn < nelements; nn++)
		{
		  first = ele[nn].firstion;
		  last = first + ele[nn].nions;
		  fprintf (fptr, "%-5s ", ele[nn].name);
		  Log ("%-5s ", ele[nn].name);
		  tot = 0;
		  for (m = first; m < last; m++)
		    tot += www->density[m];
		  for (m = first; m < last; m++)
		    {
		      fprintf (fptr, " %8.2e", log10 (www->density[m] / tot));
		      Log (" %8.2e", log10 (www->density[m] / tot));
		    }
		  fprintf (fptr, "\n");
		  Log ("\n");
		}
	    }
	}
      fclose (fptr);
      return (0);
    }
  else if (radmode == 2)
    {
      nh = geo.pl_nh;
      printf ("Starting multicycle, nh=%e\n", nh);
      /* first we need to set the agn parameters up to get the right range of ionisation parameter 100-10,000 */
//      agn_ip=1.0;     /* this will be our first value of ionisation parameter */
//      lum_agn_50=agn_ip*nh;    /*we will need to set the luminosity over 2-10KeV to implies a 2-10keV luminosity of, and agive the right value 0.1-50KeV */
//      geo.d_agn=1e15;    /* Set some reasonable distance - we have chosen to vary L */





//      lum_agn_50=lum_agn_50*geo.d_agn*geo.d_agn; /* we now know what the required luminosity is */
//      lum_agn_std=emittance_pow (0.1/HEV, 50/HEV,1.0,geo.alpha_agn); /* this is what a luminosity of 1 gives us */
//      lum_agn_scale=lum_agn_50/lum_agn_std;
//      agn_ip=emittance_pow (0.1/HEV, 50/HEV,lum_agn_scale*100.,geo.alpha_agn);
//      printf ("agn_lum goes from %e giving IP %e ",lum_agn_scale,agn_ip/(nh*geo.d_agn*geo.d_agn));
//      agn_ip=emittance_pow (0.1/HEV, 50/HEV,lum_agn_scale*10000.,geo.alpha_agn);
//      printf ("to %e giving IP %e \n",lum_agn_scale*100,agn_ip/(nh*geo.d_agn*geo.d_agn));
      /* we now know lum_agn_scale, the luminosity required to give an ionisation parameter of 1, we will go between 100 and 10000 */
      t_r = 30000.;
      for (logip = 5; logip <= 5.01; logip += 0.1)
	{
	  fudgefile = fopen ("fudge_summary.out", "a");
	  numfile = fopen ("num_summary.out", "a");
	  denomfile = fopen ("denom_summary.out", "a");
	  heatcoolfile = fopen ("heat_cool.out", "a");

	  fprintf (fudgefile, "Logip=%e \n", logip);
	  fprintf (numfile, "Logip=%e \n", logip);
	  fprintf (denomfile, "Logip=%e \n", logip);
	  fprintf (heatcoolfile, "Logip=%e \n", logip);

	  fclose (fudgefile);
	  fclose (numfile);
	  fclose (denomfile);
	  fclose (heatcoolfile);

	  agn_ip = pow (10, logip);

/* Next line calculates the luminosity over the range 0.1 to 50 KeV, this will allow us to calculate a constant for the power law */
	  lum_agn_50 = agn_ip * nh * geo.d_agn * geo.d_agn;

/* Now we have the luminosity, calculate the constant assuming L(0.1 to 50) = const x integral of nu^alpha from 0.1 to 50 */
	  geo.const_agn =
	    lum_agn_50 /
	    (((pow (50000 / HEV, geo.alpha_agn + 1.0)) -
	      pow (100 / HEV, geo.alpha_agn + 1.0)) / (geo.alpha_agn + 1.0));

/* The next line calculates the luminosity between 2 and 10KeV. This number is needed to compute the photon weight later on. This is because we are going to use the python code to do this, and when this was written, we decided to standardise our AGN luminosity to the luminosity over that range */
	  geo.lum_agn =
	    geo.const_agn * ((pow (10000 / HEV, geo.alpha_agn + 1.0)) -
			     pow (2000 / HEV,
				  geo.alpha_agn + 1.0)) / (geo.alpha_agn +
							   1.0);
	  printf
	    ("Multicycle - an IP of %e, power law coefficient of %e, which gives a 2-10keV luminosity of %e\n",
	     agn_ip, geo.const_agn, geo.lum_agn);
	  printf ("The power law coefficient =%e\n", geo.const_agn);


//              printf("logip=%f agn_ip=%e, luminosity=%e\n",logip,agn_ip,lum_agn_scale*agn_ip);
//              geo.lum_agn=lum_agn_scale*agn_ip;
//              printf ("We are setting the AGN luminosity to %e",geo.lum_agn);
	  www->w = weight = 0.5 * (1 - sqrt (1 - ((geo.r_agn * geo.r_agn) / (geo.d_agn * geo.d_agn))));	/* set the weight for BB type calculations */

	  /* 120626 -- Added loop for new plasma structure */
	  for (i = 0; i < NXBANDS; i++)
	    {
	      www->pl_w[i] = (geo.const_agn) / ((4 * PI) * (4.0 * PI * geo.d_agn * geo.d_agn));	/* Set the sim_w parameter - added 7/2/11 as part of the effort to get the sim code into python - it needs the ability to have a different weight for each cell */
	      www->pl_alpha[i] = geo.alpha_agn;	/* Set the pl_alpha parameter - added 7/2/11 as part of the effort to get the sim code into python - needs the ability to very alpha for each cell */
	    }

	  t_r = t_r + 1.;	//makes no difference, but sparks a new whole cycle, so stops using the same old densities from last loop When multicycle was first written, either w or tr would be changed in the loop.
	  www->gain = 0.5;
	  www->dt_e_old = 0.0;
	  www->dt_e = 0.0;
	  www->t_e = www->t_e_old = 1.0e6;	//Set electron temperature to sim guess
	  www->t_r = t_r;
	  t_e = www->t_e;


	  /* Find the ionization fractions and electron temperature */
	  for (n = 0; n < 20; n++)
	    {
	      Log ("multicycle %d IP= %f current t_e= %f current t_r= %f\n",
		   n, logip, www->t_e, www->t_r);
	      www->rho = nh / rho2nh;
	      printf ("wplasma[0].tr=%f tr=%f\n", www->t_r, t_r);
	      printf ("wplasma[0].w=%f w=%f\n", www->w, weight);
	      printf ("just going into cycle, nh=%e\n", nh);
	      cycle (www, p, nh, t_r, t_e, weight, mode, freq_sampling,
		     radmode);
	      t_e = www->t_e;
	      Log ("Mo  %2d %.1f %.0f %.1f %8.0f %5.2g\n", n, log10 (nh),
		   t_r, (ww), t_e, www->heat_tot / www->lum_rad);
	    }

	  /* Write out the results */
	  fprintf (fptr, "Mod   %.1f %.0f %.3f %8.0f %5.2f \n", log10 (nh),
		   t_r, logip, t_e, www->heat_tot / www->lum_rad);
	  Log ("Mod    %.1f %.0f %.1f %8.0f %5.2g\n", log10 (nh), t_r, (ww),
	       t_e, www->heat_tot / www->lum_rad);
	  for (nn = 0; nn < nelements; nn++)
	    {
	      first = ele[nn].firstion;
	      last = first + ele[nn].nions;
	      fprintf (fptr, "%-5s ", ele[nn].name);
	      Log ("%-5s ", ele[nn].name);
	      tot = 0;
	      for (m = first; m < last; m++)
		tot += www->density[m];
	      for (m = first; m < last; m++)
		{
		  fprintf (fptr, " %8.2e", www->density[m] / tot);
		  Log (" %8.2e", www->density[m] / tot);
		}
	      fprintf (fptr, "\n");
	      Log ("\n");
	    }


	}
      return (0);
    }
  else
    {
      printf ("something has gone badly wrong - radmode is not 1 or 2");
      return (0);
    }
}

int
print_levels (w)
     PlasmaPtr w;
{
  int n, m;
  for (n = 0; n < nions; n++)
    {
      if (ion[n].nlte > 0)
	{
	  m = 0;
	  while (ion[n].z != ele[m].z)
	    m++;
	  printf ("%-5s %3d %3d  ", ele[m].name, ion[n].nlevels, ion[n].nlte);

	  for (m = 0; m < ion[n].nlte; m++)
	    printf ("%8.2e ", w->levden[ion[n].first_levden + m]);
	  printf ("\n");
	}
    }
  return (0);
}

int
multi_concen (www)
     PlasmaPtr www;

{
  FILE *fopen (), *fptr;
  double nh;
  double t_r;
  double t_e;
  double weight;
  double tot;

  double t_rmin, t_rmax, t_rdelt;
  double ww, w_min, w_max, w_delt;

  int n, m, first, last, nn;

  nh = 1e12;
  t_r = 30000;
  t_e = 20000;
  weight = 0.1;

  t_rmin = 10000;
  t_rmax = 100000;
  t_rdelt = 2000;

  w_max = 0;
  w_min = -6;
  w_delt = 2;

  fptr = fopen ("con_summary.out", "w");

  for (ww = w_max; ww >= w_min; ww -= w_delt)	// Note w starts at 1
    {
      www->w = pow (10, ww);
      for (t_r = t_rmin; t_r <= t_rmax; t_r += t_rdelt)
	{
	  www->gain = 0.5;
	  www->dt_e_old = 0.0;
	  www->dt_e = 0.0;
	  www->t_e = www->t_e_old = 0.9 * t_r;	//Lucy guess
	  for (n = 0; n < 10; n++)
	    {
	      Log ("multiCycle %d\n", n);

	      www->rho = nh / rho2nh;
	      nebular_concentrations (www, 2);
	      t_e = www->t_e;

	      Log ("Mo  %2d %.1f %.0f %.1f %8.0f %5.2f\n", n, log10 (nh),
		   t_r, (ww), t_e, www->heat_tot / www->lum_rad);

	    }

	  fprintf (fptr, "Mod   %.1f %.0f %.1f %8.0f %5.2f \n", log10 (nh),
		   t_r, (ww), t_e, www->heat_tot / www->lum_rad);
	  Log ("Mod    %.1f %.0f %.1f %8.0f %5.2f\n", log10 (nh), t_r, (ww),
	       t_e, www->heat_tot / www->lum_rad);

	  for (nn = 0; nn < 5; nn++)
	    {
	      first = ele[nn].firstion;
	      last = first + ele[nn].nions;
	      fprintf (fptr, "%-5s ", ele[nn].name);
	      Log ("%-5s ", ele[nn].name);
	      tot = 0;
	      for (m = first; m < last; m++)
		tot += www->density[m];
	      for (m = first; m < last; m++)
		{
		  fprintf (fptr, " %8.2e", log10 (www->density[m] / tot));
		  Log (" %8.2e", log10 (www->density[m] / tot));
		}
	      fprintf (fptr, "\n");
	      Log ("\n");
	    }
	}
    }
  fclose (fptr);
  return (0);
}

int
check_field (one)
     PlasmaPtr one;
{

  double trad;
  double vol;

  vol = one->vol;


  one->ave_freq /= one->j;
  one->j /= (4. * PI * vol);
  one->t_r = H * one->ave_freq / (BOLTZMANN * 3.832);
  one->w = PI * one->j / (STEFAN_BOLTZMANN * trad * trad * trad * trad);

  Log ("Check_field: trad %.2f  w %.2e\n", trad, one->w);
  return (0);
}

int
do_detail (one)
     PlasmaPtr one;
{
  double newden[NIONS];
  int n, nelem;

  for (nelem = 0; nelem < nelements; nelem++)
    {
      detailed_balance (one, nelem, newden);
      wind_update_after_detailed_balance (one, nelem, newden);

    }
  for (n = 0; n < NIONS; n++)
    {
      one->density[n] = newden[n];
    }
  return (0);
}
