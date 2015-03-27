
/***********************************************************
                                       Space Telescope Science Institute

Synopsis:  This file contains routines that make use of atomic data to
	calculate basic atomic rates.  (THESE ROUTINES SHOULD NOT DEPEND
	ON A HAVING A WIND ARRAY!!!!
Arguments:


Returns:
 
Description:	 
	
        02jun   ksl     Created from existing routines in python.

**************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "atomic.h"



/***********************************************************
                                       Space Telescope Science Institute

 Synopsis: 
	double sigma_phot(x_ptr,freq)	calculates the photionization crossection due to the transition 
	associated with x_ptr at frequency freq
Arguments:

Returns:
	
Description:	 
	sigma_phot uses Verner et al.'s interpolation formulae for the photoionization crossection
	to calculate the bound free (or photoionization) optical depth.  The data must
	have been into the photoionization structures xphot with get_atomic_data and
	the densities of individual ions must have been calculated previously.

Notes:

History:
	98jul	ksl	Coded (actually moved from a subroutine called kappa_ds)

**************************************************************/

double
sigma_phot (x_ptr, freq)
     struct photoionization *x_ptr;
     double freq;
{
  double ft;
  double x, y;
  double f1, f2, f3;
  double xsection;

  ft = x_ptr->freq_t;		/* threshold frequency */

  if (ft < freq && freq < x_ptr->freq_max)
    {
      x = freq / x_ptr->freq0 - x_ptr->y0;
      y = sqrt (x * x + x_ptr->y1 * x_ptr->y1);

      /* This was line fixed by CK in 1998 Jul */
      f1 = (x - 1.0) * (x - 1.0) + x_ptr->yw * x_ptr->yw;

      f2 = pow (y, 0.5 * x_ptr->p - 5.5);
      f3 = pow (1.0 + sqrt (y / x_ptr->ya), -x_ptr->p);
      xsection = x_ptr->sigma * f1 * f2 * f3;	// the photoinization xsection

      return (xsection);
    }
  else
    return (0.0);
}



/* 
 *    a21 alculates and returns the Einstein A coefficient 
 *       History:
 *          98aug        ksl     Coded and debugged
 *             99jan        ksl Modified so would shortcircuit calculation if 
 *                called multiple times for same a
 *                 */
#define A21_CONSTANT 7.429297e-22	// 8 * PI * PI * E * E / (MELEC * C * C * C)

struct lines *a21_line_ptr;
double a21_a;

double
a21 (line_ptr)
     struct lines *line_ptr;
{
  double freq;

  if (a21_line_ptr != line_ptr)
    {
      freq = line_ptr->freq;
      a21_a =
	A21_CONSTANT * line_ptr->gl / line_ptr->gu * freq * freq *
	line_ptr->f;
      a21_line_ptr = line_ptr;
    }

  return (a21_a);
}
