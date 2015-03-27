/***********************************************************
                                       Southampton university

 Synopsis:
	test_dielectronic is a stand alone code to produce data relating to trying to understand dielectronic recombination in all its forms!
 
Arguments:		


Returns:
 
Description:	
	
Notes:

History:
	 
	0512 - nsh wrote it
 	


**************************************************************/



#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "atomic.h"
int
main (argc, argv)
     int argc;
     char *argv[];
{
  WindPtr w;
  PhotPtr p;
  int n, i, j, nion;
  int ntemp, nne, temp_start, tempdiv, ntmin;
  double ne[11], temp1, nh, xtemp, weight, xsaha, t_r, t_e, temp;
  double b, *partition, recomb_fudge, pi_fudge, gs_fudge, test, integral;
  double fthresh, fmax;
  double alpha_agn, distance, agn_ip, lum_agn, const_agn;
  int lstart, lmax;
  double temp_func ();
  FILE *fp_h;
  FILE *fp_he;
  FILE *fp_c;
  FILE *fp_n;
  FILE *fp_o;
  FILE *fp1_o;
  FILE *fp_fe;
  PlasmaPtr xplasma;
  partition = xplasma->partition;	/* Set the partition function array to that held for the cell */
  temp_start = 299;
  ntemp = 401;
  tempdiv = 100.0;
  nne = 1;
  nh = 1e10;
  weight = 0.01;
  xne = nh;

  strcpy (geo.atomic_filename, "atomic/standard70");
  printf ("atomic_filename=%s\n", geo.atomic_filename);
  get_atomic_data (geo.atomic_filename);

  printf (*, "We are alive!!\n");

  for (i = 300, i < 901, i++)
    {
      temp = pow (float (i) / 100);
      printf ("%f", temp);
      compute_dr_coeffs (temp);
      for (j = 2; j < 5; j++)
	{
	  printf (*, " %6.3e", dr_coeffs[n]);
	}
      printf (*, "\n");
    }





  return EXIT_SUCCESS;
}
