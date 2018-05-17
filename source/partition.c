
/***********************************************************/
/** @file  partition.c
 * @author ksl,nsh
 * @date   April 2018
 *
 * @brief  This file includes calculations of  partition functions
 * for ions in plasma cells under a number of assumptions
 *
 * @bug It seems likely that this file and levels.c should be combined
 * into a single file.  There are very comments (by nsh) about eliminating
 * some of the modes that need investigation as well.  
 *
 ***********************************************************/



#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>


#include "atomic.h"
#include "python.h"


//OLD /***********************************************************
//OLD                                        Space Telescope Science Institute
//OLD 
//OLD Synopsis:
//OLD 	partition_functions calculates the partition functions
//OLD 	for a single cell in the grid
//OLD 
//OLD Arguments:
//OLD 
//OLD 
//OLD Returns:
//OLD 	
//OLD Description:
//OLD 
//OLD Mode here is identical to that used by nebular_concentrations, e.g
//OLD 
//OLD 0 = LTE using t_r
//OLD 1 = LTE using t_e
//OLD 2 = Lucy and Mazzali
//OLD 
//OLD 
//OLD 
//OLD Notes:
//OLD 
//OLD 	0800802 - This is a rewritten version of two routines, do_partitions
//OLD 	and partition.   It is still true that levels is almost identical
//OLD 	to some of the code here and therefore it is unclear why it is needed.
//OLD 	levels is always called from partition_functions
//OLD 
//OLD 
//OLD History:
//OLD 	080802	ksl	60b -- Brought levels routine into
//OLD 			do_partitions since they are always
//OLD 			called one after the other.  The routines
//OLD 			are almost identical and it remains
//OLD 			unclear why this was coded twice
//OLD 	080802	ksl	Combined two routines do_partitions
//OLD 			and partition into one in order to
//OLD 			make this routine very similar to
//OLD 			levels
//OLD 	080804	ksl	Removed hubeny calculation of g from
//OLD 			this code.  The advantage of the 
//OLD 			Hubeny calculation was that it had
//OLD 			a correction for density that we
//OLD 			no longer had, but it was not being
//OLD 			accessed by the code at all.
//OLD 
//OLD **************************************************************/




/**********************************************************/
/** 
 * @brief      calculates the partition functions
 * 	for a single cell in the grid (and then calls a routine that
 * 	calculates the level populations
 *
 * @param [in] PlasmaPtr  xplasma   A plasma cell in the wind
 * @param [in] int  mode   An integer which determines the way in which the partition funciton will be calculated
 * @return     Always returns 0
 *
 * @details
 *
 * The possibilites are:
 * * NEBULARMODE_TR 0        LTE using t_r
 * * NEBULARMODE_TE 1        LTE using t_e
 * * NEBULARMODE_ML93 2      ML93 using a nebular approximation correction to LTE
 * * NEBULARMODE_NLTE_SIM 3  // Non_LTE with SS modification (Probably could be removed)
 * * NEBULARMODE_LTE_GROUND 4        // A mode which forces all levels to the GS - this is used when we are modelling J_nu.
 *
 *
 * ### Notes ###
 * 0800802 - This is a rewritten version of two routines, do_partitions
 * 	and partition.   It is still true that levels is almost identical
 * 	to some of the code here and therefore it is unclear why it is needed.
 * 	levels is always called from partition_functions
 *
 **********************************************************/

int
partition_functions (xplasma, mode)
     PlasmaPtr xplasma;
     int mode;
{
  int nion;
  double partition ();
  double t, weight;

  int n, m;
  int m_ground;                 /* added by SS Jan 05 */
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
  else if (mode == NEBULARMODE_NLTE_SIM)        /*NSH 120912 This mode is more or less defunct. When the last vestigies of the mode 3 ionizetion scheme (the original sim PL correction) is removed, this can go too */
  {
    //Non LTE calculation with non BB radiation field. Use T_e to get partition functions, same as mode 1- 
    t = xplasma->t_e;
    weight = 1;
  }
  else if (mode == NEBULARMODE_LTE_GROUND)      /*This is used to set partition functions  to ground state only, used when
	  											we have a J_nu model rather than using dilute blackbodies.
	  											This is achieved by setting W to 0. At this point, 
                                                  the temperature is a moot point, so lest go with t_e, since this is only 
                                                  going to be called if we are doing a power law calculation */
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
        z += weight * config[m].g * exp ((-config[m].ex + config[m_ground].ex) / kt);
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
        z += weight * config[m].g * exp ((-config[m].ex + config[m_ground].ex) / kt);
      }
    }
    else
    {
      z = ion[nion].g;
    }


    xplasma->partition[nion] = z;

  }

  levels (xplasma, mode);       //levels from reduced BB

  return (0);
}





//OLD /***********************************************************
//OLD                                        Southampton university
//OLD 
//OLD Synopsis:
//OLD 	partition_functions_2 calculates the partition functions
//OLD 	for a pair of states for a given temperature. This is needed
//OLD 	to support the pairwise ioinization state calculations where 
//OLD 	the saha equation is applied to a pair of states at a
//OLD 	useful temperature, then corrected. Part of this is 
//OLD 	a requirement to get the partition functions at that 
//OLD 	temperature for just tow two states of interest. It is 
//OLD 	wasteful to calculate all of the states at each temperature.
//OLD 
//OLD Arguments:
//OLD 
//OLD 	xplasma - the cell we are interested in. Used to communicate back
//OLD 		the new partition functions.
//OLD 	xnion - the upper ion in the pair we are currently working on
//OLD 	temp - the temperature we are using for the saha equation. In the
//OLD 		original formulation of partition, this is got from xplasma,
//OLD 		in this case we don't want to use the cell temperature.
//OLD 
//OLD Returns:
//OLD 
//OLD 	changes the partition functions of just two ions we are currently working on
//OLD 	
//OLD Description:
//OLD 
//OLD 
//OLD 
//OLD Notes:
//OLD 
//OLD 	There is no need for the weight term in the calculation since this is only called when
//OLD 		we are making the assumption that we are in LTE for the saha equation. 
//OLD 
//OLD 
//OLD 
//OLD History:
//OLD 	2012Feb	nsh - began coding
//OLD 	2012Sep	nsh - added weight as a calling value - this allows one to produce gs 
//OLD                 only with w=0, or LTE with w=1 or to produce a correction factor with W = the measured value	
//OLD 
//OLD **************************************************************/




/**********************************************************/
/** 
 * @brief      calculates the partition functions
 * 	for a pair of ions at a given temperature. 
 *
 *
 * @param [in] PlasmaPtr  xplasma   - the cell we are interested in. Used to communicate back
 * @param [in] int  xnion   The ion number for the upper ion of a pair 
 * @param [in] double  temp  the temperature we are using for the saha equation. 
 * @param [in] double  weight The ratio of the intensity of the radiation field to that of a BB
 * @return     changes the partition functions of just two ions we are currently working on
 *
 * @details
 * 	This routine is used 
 * 	to support the pairwise ioinization state calculations where 
 * 	the saha equation is applied to a pair of states at a
 * 	temperature, then corrected. Part of this is 
 * 	a requirement to get the partition functions at that 
 * 	temperature for just two states of interest. It is 
 * 	wasteful to calculate all of the states at each temperature.
 *
 * 	The results are stored in xplasma->partition[nion]
 *
 * ### Notes ###
 * @bug According to the historical notes, there is no need for the weight term in the calculation 
 * since this is only called when
 * we are making the assumption that we are in LTE for the saha equation.
 *
 * One of the ways in which this routine differs from partition_functions is that
 * here the temperature and weight are passed to the routing directly instead of
 * being taken from the Plasma structure.
 *
 **********************************************************/

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
        z += weight * config[m].g * exp ((-config[m].ex + config[m_ground].ex) / kt);
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
        z += weight * config[m].g * exp ((-config[m].ex + config[m_ground].ex) / kt);
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
