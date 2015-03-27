


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>


#include "atomic.h"
#include "python.h"


/***********************************************************
                                       Space Telescope Science Institute

Synopsis:
	partition_functions calculates the partition functions
	for a single ceell in the grid

Arguments:


Returns:
	
Description:

Mode here is identical to that used by nebular_concentrations, e.g

0 = LTE using t_r
1 = LTE using t_e
2 = Lucy and Mazzali



Notes:

	0800802 - This is a rewritten version of two routines, do_partitions
	and partition.   It is still true that levels is almost identical
	to some of the code here and therefore it is unclear why it is needed.
	levels is always called from partition_functions


History:
	080802	ksl	60b -- Brought levels routine into
			do_partitions since they are always
			called one after the other.  The routines
			are almost identical and it remains
			unclear why this was coded twice
	080802	ksl	Combined two routines do_partitions
			and partition into one in order to
			make this routine very similar to
			levels
	080804	ksl	Removed hubeny calculation of g from
			this code.  The advantage of the 
			Hubeny calculation was that it had
			a correction for density that we
			no longer had, but it was not being
			accessed by the code at all.

**************************************************************/



int
partition_functions (xplasma, mode)
     PlasmaPtr xplasma;
     int mode;
{
  int nion;
  double partition ();
  double t, weight;

  int n, m;
  int m_ground;			/* added by SS Jan 05 */
  double z, kt;


  if (mode == NEBULARMODE_TR)
    {
      //LTE using t_r
      t = xplasma->t_r;
      weight = 1;
    }
  else if (mode == NEBULARMODE_TE)
    {
      //LTE using t_e
      t = xplasma->t_e;
      weight = 1;
    }
  else if (mode == NEBULARMODE_ML93)
    {
      //Non LTE calculation with radiative weights
      t = xplasma->t_r;
      weight = xplasma->w;
    }
  else if (mode == 3)		/*NSH 120912 This mode is more or less defunct. When the last vestigies of the mode 3 ionizetion scheme (the original sim PL correction) is removed, this can go too */
    {
      //Non LTE calculation with non BB radiation field. Use T_e to get partition functions, same as mode 1- 
      t = xplasma->t_e;
      weight = 1;
    }
  else if (mode == 4)		/*NSH 120912 This is a test mode, used to set partition functions to ground state only. This is achieved by setting W to 0. At this point, the temperature is a moot point, so lest go with t_e, since this is only going to be called if we are doing a power law calculation */
    {
      t = xplasma->t_e;
      weight = 0;
    }

  else
    {
      Error ("partition_functions: Unknown mode %d\n", mode);
      exit (0);
    }

  /* Calculate the partition function for each ion in turn */
  kt = BOLTZMANN * t;

  for (nion = 0; nion < nions; nion++)
    {

      if (ion[nion].nlevels > 0)
	//Calculate data on levels using a weighed BB assumption
	{
	  m = ion[nion].firstlevel;
	  m_ground = m;
	  //store ground state - in case energy neq 0(SS)
	  z = config[m].g;
	  //Makes explicit assumption that first level is ground

	  for (n = 1; n < ion[nion].nlevels; n++)
	    {
	      m++;
	      z +=
		weight * config[m].g *
		exp ((-config[m].ex + config[m_ground].ex) / kt);
	    }
	}
      else if (ion[nion].nlte > 0)
	//Calculate using "non-lte" levels
	{
	  m = ion[nion].first_nlte_level;
	  m_ground = m;
	  //store ground state - in case energy neq 0(SS)
	  z = config[m].g;
	  //This statement makes an explicit assumption that first level is ground

	  for (n = 1; n < ion[nion].nlte; n++)
	    {
	      m++;
	      z +=
		weight * config[m].g *
		exp ((-config[m].ex + config[m_ground].ex) / kt);
	    }
	}
      else
	{
	  z = ion[nion].g;
	}


      xplasma->partition[nion] = z;

    }

  levels (xplasma, mode);	//levels from reduced BB

  return (0);
}





/***********************************************************
                                       Southampton university

Synopsis:
	partition_functions_2 calculates the partition functions
	for a pair of states for a given temperature. This is needed
	to support the pairwise ioinzation state calculaions where 
	the saha equation is applied to a pair of states at a
	useful temperature, then corrected. Part of this is 
	a requirement to get the partition functions at that 
	temperature for just tow two states of interest. It is 
	wasteful to calculate all of the states at each temperature.

Arguments:

	xplasma - the cell we are interested in. Used to communicate back
		the new partition functions.
	xnion - the upper ion in the pair we are currently working on
	temp - the temperature we are using for the saha equation. In the
		original formulation of partition, this is got from xplasma,
		in this case we don't want to use the cell temperature.

Returns:

	changes the partition functions of just two ions we are currently working on
	
Description:



Notes:

	There is no need for the weight term in the calculation since this is only called when
		we are making the assumption that we are in LTE for the saha equation. 



History:
	2012Feb	nsh - began coding
	2012Sep	nsh - added weight as a calling value - this allows one to produce gs only with w=0, or LTE with w=1 or to produce a correction factor with W = the measured value	

**************************************************************/



int
partition_functions_2 (xplasma, xnion, temp, weight)
     PlasmaPtr xplasma;
     int xnion;
     double temp;
     double weight;
{
  int nion;
  double partition ();

  int n, m;
  int m_ground;
  double z, kt;


  /* Calculate the partition function for each ion in turn */
  kt = BOLTZMANN * temp;

  for (nion = xnion - 1; nion < xnion + 1; nion++)
    {

      if (ion[nion].nlevels > 0)
	//Calculate data on levels using a weighed BB assumption
	{
	  m = ion[nion].firstlevel;
	  m_ground = m;
	  //store ground state - in case energy neq 0(SS)
	  z = config[m].g;
	  //Makes explicit assumption that first level is ground

	  for (n = 1; n < ion[nion].nlevels; n++)
	    {
	      m++;
	      z +=
		weight * config[m].g *
		exp ((-config[m].ex + config[m_ground].ex) / kt);
	    }
	}
      else if (ion[nion].nlte > 0)
	//Calculate using "non-lte" levels
	{
	  m = ion[nion].first_nlte_level;
	  m_ground = m;
	  //store ground state - in case energy neq 0(SS)
	  z = config[m].g;
	  //This statement makes an explicit assumption that first level is ground

	  for (n = 1; n < ion[nion].nlte; n++)
	    {
	      m++;
	      z +=
		weight * config[m].g *
		exp ((-config[m].ex + config[m_ground].ex) / kt);
	    }
	}
      else
	{
	  z = ion[nion].g;
	}


      xplasma->partition[nion] = z;

    }






  return (0);
}
