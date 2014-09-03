/***********************************************************
                                       Space Telescope Science Institute

 Synopsis:
 
 	These are a series of routines which are designed to generate a random sampling
 	of a probability distribution using cumulative distribution functions.
 	
 	The main routines are
 	
		pdf_gen_from_func	(&pdf,&func,xmin,xmax,njumps,jump)		
			Generate a cumulative distrubution function (cdf) from a function
 	          
		pdf_gen_from_array(&pdf,x,y,n_xy,xmin,xmax,njumps,jump)		
			Generate a cumulative distribution function (cdf) from an
 			array of values

		pdf_get_rand(&pdf)				
			Generate a single sample from a cdf defined by either of the first
 			two routines 	
 	
 	Sometimes it is expensive computationally to recalculate the distribution function
 	every time you change the limits.   An example of this is the distribution function
 	for a bb where the distribution function is quite smooth (and can be represented in a 
 	dimensionless manner) but calculating the function and the distribution function
 	can be expensive.  The following two routines allow one to sample a part of a
 	pdf.  Note that pdf_gen-from_func or pdf_gen_from_array should have been
 	called previously!
 
		pdf_limit(&pdf,xmin,xmax)  sets limit1 and limit2 so that one can generate
 			distributions from a limited range within a previously created pdf.

		pdf_get_rand_limit(&pdf) will samples the pdf but limits the range of
			the return to lie between xmin and xmax from pdf_limit.

	
	There are some internal routines

		gen_array_from_func (func, xmin, xmax, pdfsteps) generates a cdf from a function
			on values linearly spaced between xmin and xmax.  pdfsteps gies the
			number of points.

 	There is also a simple diagnostic routine
		pdf_check(&pdf)

	which returns a nonzero number if certain checks of the pdf fail
	and a routine  

		pdf_to_file(&pdf,filename)				
 		  	to write a pdf structure to a file


								
Arguments for pdf_gen_from_func
	PdfPtr pdf;		A ptr to a pdf structure (see below)
	double (*func)()	The probablility density function to be integrated 
				to create the pdf
	double xmin,xmax;	The range over which the function is to be 
				integrated
	double jump[];		Specifix values of x, perhaps a Lyman edge, to be
				included in the pdf_structure.  These must be in
				ascending order but do not have to lie between xmin
				and xmax
	int njumps;		The size of jumps[]

Arguments for pdf_gen_from_array are same as above except
 	
	double x[],y[]		The one dimensional arrays containg the places x where the cumulative 
				distribution function density has been evaluated.  The array need not 
				be uniformly spaced. y is the value of the cdf at the point x 
	double d[]		The rate of change of the probability density at x, which is calculated
				from the CDF. -- Added 57i
	int n_xy		The size of x[] and y[]

Arguments for pdf_get_rand
	PdfPtr pdf		pdf (See below)	is a structure which contains the cumulative
				distribution function.

Arguments for pdf_limit
	double xmin,xmax	Here xmin and xmax will cause pdf[].limit1 and pdf[].limit2 to be set
				in such a way that when you call pdf_get_rand_limit the distribution
				will be sampled only between xmin and xmax.  Watch out though
				because if xmin and xmax lie outside of the original value of xmin
				and xmax then the distribution will not extend this far.							
Arguments for pdf_to_file
	PdfPtr pdf
	char filename[]		If you can't figure this out, you are in deep trouble!
				
Returns:
	All of the routines return 0 on successful completion, except pdf_get_rand which returns
	a random value between xmin and xmax.  Multiple calls to pdf_get_rand should regenerate
	the probability density function.
	
	If there is an error the routines usually send an error message to the screen and logfile
	and cause the program to stop!  Occasionally the program will send a nonzero return to
	indicate an error.
 
Description:

	The purpose of these routines is first to generate a cumulative distribution function
	and then allow one to sample that distribution function.  The procedure is either
	to call pdf_gen_from_func if you have a function whose probability distribution you
	want (as might be the case if you were generating bremsstrahlung photons and you knew
	the formula)  or pdf_gen_from_array (as might be the case if you were trying to
	generate Monte Carlo spectra from precalculated Kurucz models).
	
	The result of calling either of these routines will be to populate a structure of
	the form Pdf
		pdf->x[]  will contain NPDF values of x between xmin and xmax.  The values will be
			spaced so that it is almost but not quite true that if you generate a uniform
			number n between 0 and NPDF you could just read off the value pdf[n].x and
			that by doing this multiple times you could regenerate the distribution function
		pdf->y[]  will contain values between 0 and 1 and is the integral of the cumulative
			distribution function from xmin to pdf[].x
			
	Once a pdf has been generated, one samples the pdf by calls to pdf_get_rand which 
	creates a random number between 0 and 1, finds the elements in pdf which surround
	the random number and interpolates to return a value x. 
	
	It is possible to force specific values of x to appear in the pdf.  This is desirable
	if there are edges where the probability density changes rapidly, as for example
	at the Lyman edge in a stellar spectrum.  By using jumps properly you can assure
	that edges will be sharp in the regenerated distribution (since otherwse pdf_get_rand
	will linearly interpolate between two points on the opposite sides of the edge).

Notes:
	Must be compiled with log.c since it makes use of the "python" logging routines for
	errors.

	The fact that these routines start with pdf... is historical and represents a poor
	choice of nomenclature.  We reallly mean comulative distribution functions

	It would probably be desirable to modify this so that the size of the parallel 
	arrays being generated could be determined from inside the program.  In that case
	NPDF would have to become the maximum size.


	Error??: For Kurucz models, The probability distribution beyond the He edge  
	is often very peaked at the long wavelength end and rapidly changing.  In this 
	instance, some bizarre results can be obtained because NPDF is not enough to 
	properly populate the pdf.  One could increase NPDF, but this would slow 
	everything down.  This is an another argument that the whole routine pdf_gen_from_array
	neeeds to be rewritten to avoid using NPDF altogether.  One way to do this
	would be simply to use the values of x and y in creating an inital pdf.  There
	is however a question of resampling.
	

History:
	97aug	ksl 	Initial versions of these and similar routines coded as part of python project.
	98mar	ksl	Recoded all of the pdf routines so that jumps were incorporated and so
			that they all use the pdf structure. 
	98mar25	ksl	Fixed a problem with gen_from_func which occured when there were
			discontinues changes in the pdf which resulted in multiple values of
			the pdf being associated with the same probability.  Also, fixed I had
			not done proper accounting when there were specific values at which
			I want the pdf recorded. 
	98jul27	ksl	Fixed a problem with renormalization inf pdf_gen_from_func for the case when
			one wanted the pdf at specific values.  Discovered, and have fixed a problem
			in pdf_gen_from_array and pdf_gen_from_func associated with requireing the
			pdf at specific values.  Note: while I have tested this for pdf_gen_from func
			I have not really tested it for pdf_gen_from_array and this needs to be done.!!
			Added a routine pdf_check which looks for errors in a pdf_array.
	98jul27 ksl	Added check to assure that if one does put in points at specific values
                        that the probability density was not zero over this interval which would
                        produce two y values that are exactly the same.	
	00dec30	ksl	Fixed small error to change order of truth calculation in a while statement that
			was causing the program to fail on Proga's alpha
	06sep	ksl	57i -- Previously, one assumed the probability density was effectively
			uniform between the places where the cumulative distribution function
			was calculated.  Now we attempt to make a correction for this.
	06oct	ksl	58 -- After a lot of testing, I decided to baseline the version
			of pdf that makes a first order correction for non-uniformity of
			the probability distribution within intervals where the cdf is
			calculated.  The version here corresponds to 57ib.  I still don't understand
			some of the timing issues associated with this, as there are circumstances
			where this more excact calculation seems to take less time to exeucute
			in python than the old version.  The old version should still be in the
			directory (and is called pdf_uniform.c)
	10oct	ksl	Changed the way in wich pdf_gen_from function works so it will change the
			the number of points that are used in the array pdf_array in situations
			where the binning is too course.  In the process, eliminated pdf_init.  
			It was only called by one routine.  
 
**************************************************************/

#include <stdio.h>
#include <stdlib.h>

#include "atomic.h"
#include "python.h"

/*  The structure is defined in python.h.  Here for reference only */
//#define NPDF 200

//typedef struct Pdf {
//      double x[NPDF+1];
//      double y[NPDF+1];
//      double d(NPDF + 1];      
//      double limit1,limit2;
//      double norm;           //The scaling factor which would renormalize the pdf
//} *PdfPtr,pdf_dummy;


#define PDFSTEPS 10000		// This is the initial value of PDFSTEPS
int pdf_steps_current;		// This is the value of pdfsteps at this point in time
int init_pdf = 0;
double *pdf_array;

/* Generate a pdf structure from a function.  

Description:

This routine stores values of the function in  pdf_array

	02feb	ksl	Modified to handle a problem in which the
			pdf was not monotonic because of the 
			interaction between jumps and desired
			values of the pdf.  This is probably more
			like pdf_gen_from_array.  At this point
			it is not entirely clear to me that 
			pdf_gen_from function should not simple
			call pdf_gen_from array...but that is
			not currently the case 
	10oct	ksl	Modified (and reorganized) in an attempt to calculate a 
			finer grid of points in pdf_array, in situations
			where array was so course that we did not
			adequately sample the CDF.

*/
int
pdf_gen_from_func (pdf, func, xmin, xmax, njumps, jump)
     PdfPtr pdf;
     double (*func) (double);
     double xmin, xmax;
     double jump[];
     int njumps;
{
  double xstep;
  double y;
  int j, m, mm, n;
  int njump_min, njump_max;
  int icheck, pdfsteps;
  int pdf_check (), recalc_pdf_from_cdf ();
  double gen_array_from_func (), delta;

  njump_min = njump_max = 0;
  /* Check the input data before proceeding */
  if (xmax <= xmin)
    {
      Error ("pdf_gen_from_func: xmin %g <= xmax %g\n", xmin, xmax);
      exit (0);
    }

  if (njumps > 0)
    {
      for (j = 1; j < njumps; j++)
	{
	  if (jump[j] <= jump[j - 1])
	    {
	      Error
		("pdf_gen_from_func: jump[%d]=%g <=jump[%d]=%g out of order\n",
		 j, jump[j], j - 1, jump[j - 1]);
	      exit (0);
	    }
	}

      njump_min = 0;
      while (njump_min < njumps && jump[njump_min] <= xmin)
	njump_min++;
      njump_max = 0;
      while (njump_max < njumps && jump[njump_max] < xmax)
	njump_max++;
      /* So at this point njump_min will point to the first jump which is betewen xmin and
         xmax or it will equal to njumps in which case there were no jumps which were greater
         than xmin.

         Similarly njump_max will be the point to the first jump above xmax or if there are
         no jumps above xmax, then it will be njumps. */

      njumps = njump_max - njump_min;
    }

  /* OK, all the input data seems OK */

  /* Construct what is effectively is the definite integral from xmin to x. Note
     however that currently pdf_array[0] corresponds awkwardly to the integral
     from xmin to xstep.  

     There code increases the size of the array if the steps appear so large that they
     will affect the digitization.  

     Problems might occur if there are very large jumps in the function. 

     10oct - ksl - Incorprated into the code in attempting to create a power law.
   */

  delta = 1.;
  n = 0;
  pdfsteps = PDFSTEPS;

  while (n < 3)
    {
      delta = gen_array_from_func (func, xmin, xmax, pdfsteps);
      if (delta < 0.1 / NPDF)
	break;
      pdfsteps *= 10;
      n = n + 1;
    }

  xstep = (xmax - xmin) / pdfsteps;

/* So at this point pdf_array contains an unnormalized version of the CDF for the function */
  pdf->x[0] = xmin;
  pdf->y[0] = 0;

  n = 0;			//This is the position in pdf_array
  mm = 1;			//This is the index to a desired value of y, with no jumps
  j = njump_min;		//This refers to the jumps
  for (m = 1; m < NPDF; m++)
    {
      y = (float) mm / (NPDF - njumps);	// Desired value of y ignoring jumps

      while (pdf_array[n] < y && n < pdfsteps)	// Work one's way through pdf_array
	{
	  if (j < njump_max && jump[j] <= xmin + (n + 1) * xstep)
	    {
	      pdf->x[m] = xmin + (n + 1) * xstep;	//Not exactly jump but close
	      pdf->y[m] = pdf_array[n];
	      j++;		//increment the jump number
	      m++;		//increment the pdf structure number
	    }
	  n++;
	}

      /* So at this point pdf_array[n-1] < x and pdf_array[n]>x */
      pdf->x[m] = xmin + (n + 1) * xstep;
      pdf->y[m] = pdf_array[n];
      mm++;			// increment the number associated with the desired y ignoring jumps
      /* So pdf->y will contain numbers from 0 to 1 */
    }

  pdf->x[NPDF] = xmax;
  pdf->y[NPDF] = 1.0;
  pdf->norm = 1.;		/* pdf_gen_from array produces a properly nomalized cdf and so the
				           normalization is 1.  110629 ksl */

/* Calculate the gradients */
  recalc_pdf_from_cdf (pdf);	// 57ib 
  /* Check the pdf */
  if ((icheck = pdf_check (pdf)) != 0)
    {
      Error ("pdf_gen_from_function: error %d on pdf_check\n", icheck);
    }
  return (icheck);

}


/* 

gen_array_from_func


This is a routine which is called by pdf_gen_from_func which simply calculates the cumulative
distribution of the function in equally spaced steps between xmin and xmax.  The CDF is properly
normalized.


pdfsteps is the number of steps that are calculated.  The routine returns the largest change
in the cdf between any two points in the grid.

Notes:

10oct - ksl -This is rather brute force.  An alternative would have been to have increased the density of points
only in the region where it was needed.  This could be done by looking at the regions where the
sampling was poor in the initial calculation and addeding new points in in a new calculations.  But
this would require one to pass both the values of the CDF and the positions where the values were
taken.  This would be  not be hard but would require more extensive modifications to pdf_gen_from func than 
currently.  

11jun	ksl	The name of this function is rather misleading reflecting its history.  At one point
		it simple generated the value of the function at the various points in the grid.  That
		is no longer the case, and has led to confusion.  

History

10oct	ksl(python_69)	Coded to enable one to adaptively increase the density of points
11jun	ksl(69d)	This routine creates a properly normalized CDF.

*/

double
gen_array_from_func (func, xmin, xmax, pdfsteps)
     double (*func) (double);
     double xmin, xmax;
     int pdfsteps;
{

  double x, z, xstep;
  double sum;
  int m, n;
  int idiag;
  double delta;

  xstep = (xmax - xmin) / pdfsteps;
  idiag = 0;

  /* First, allocate an array for internal use, if that
     has not already been done */


  if (init_pdf == 0 || pdf_steps_current < pdfsteps)
    {

      if (pdf_array != NULL)
	free (pdf_array);

      if ((pdf_array = calloc (sizeof (x), pdfsteps)) == NULL)
	{
	  Error ("pdf: Could not allocate space for pdf_array\n");
	  exit (0);
	}
      init_pdf = 1;
      pdf_steps_current = pdfsteps;
    }


  for (n = 0; n < pdfsteps; n++)
    {
      x = xmin + (n + 0.5) * xstep;
      if ((z = (*func) (x)) < 0 || z > 1.e50 || sane_check (z))
	{
	  Error ("pdf_gen_from_func: probability density %g < 0 at %g\n", z,
		 x);
	}
      if (n == 0)
	pdf_array[0] = z;
      else
	pdf_array[n] = pdf_array[n - 1] + z;
      /* Check to see if the integral seems too large */
      if (pdf_array[n] > 1.0e100 && idiag == 0)
	{
	  idiag = 1;
	  Log ("pdf_gen_from_func: ZOWIE  n %d z %g pdf_array[n] %g x %g\n",
	       n, z, pdf_array[n], x);
	  for (m = 0; m < n; m += 100)
	    {
	      z = (*func) (x);
	      Log ("pdf_gen_from_func: zowie n %d x %g z %g pdf_array %g\n",
		   m, x = xmin + (0.5 * m) * xstep, z, pdf_array[m]);
	    }
	}
    }
  /* Thus, pdf_array is proportional to  the definite integral from xmin to x 
     (where x=xmin+(n+1)*xstep) */

  sum = pdf_array[pdfsteps - 1];

  if (sane_check (sum))
    {
      Error ("pdf_gen_from_func:sane_check Sum %f is NaN\n", sum);
    }

  /* Renormalize the array so that this really is a cumulative distribution function */
  delta = 0;
  for (n = 0; n < pdfsteps; n++)
    {
      pdf_array[n] /= sum;
      if (n > 0)
	{
	  x = pdf_array[n] - pdf_array[n - 1];
	  if (x > delta)
	    {
	      delta = x;
	    }
	}

    }


  return (delta);
}

// End of function


/*  
pdf_gen_from_array

The one dimensional arrays containg the places x where the probability
density y has been evaluated.  The array need not be uniformly spaced.  
The assumption made is that y is the probability at the point x and that 
one can linearly interpolate between points 

Notes:


	06sep -- There are differences between this and pdf_gen_from_func
	that look historical -- ksl


History:
	98jul	ck	Added check to see if pdf, really cdf, is all zeros
	98jul27 ksl	Changed the logic in the way the pdf array is created.  
			The old method was not correct.  
	02jul	ksl	Rewrote much of routine to try to handle situations in
			which the array implied a distribution which was very skewed 
			toward either end. Specifically, I eliminated the old approach
			that did a pseudo integration along the lines of the one
			in pdf_gen_from_func. The new approach should also speed the
			routine up significantly, and eliminate some of the other
			errors which have plagued earlier versions, e.g. non-monotonic
			pdfs and values just out of bounds.  (See the comments in
			python_43/pdf.c)
	04mar	ksl	A familiar problem of the pdf being out of order appeared 
			in pdf_gen_from_array when I increased the number of 
			elements in the pdf array.  I found an error in the way q 
			for jumps was calculated and added checks (and a shuffle)
			for situations where the cdf was out of order.
	090122	ksl	Increase size of PDF_ARRAY to reflect the fact that we are
			tending towards models with more points in them and updated 
			the error message associated with allowing too few points
			so that problems with this would be easier to update in
			future.
	1405	JM -- Increased PDF array for use with disk14 models
*/

#define PDF_ARRAY  28000

double pdf_x[PDF_ARRAY], pdf_y[PDF_ARRAY], pdf_z[PDF_ARRAY];
int pdf_n;

int
pdf_gen_from_array (pdf, x, y, n_xy, xmin, xmax, njumps, jump)
     PdfPtr pdf;
     double x[], y[];
     int n_xy;
     double xmin, xmax;
     int njumps;
     double jump[];
{
  int allzero;
  int m, n, nn, j;
  double sum, q, xx, yy;
  int njump_min, njump_max;
  double ysum;
  int echeck, pdf_check (), recalc_pdf_from_cdf ();


  /* Check the inputs */
  if (xmax < xmin)
    {
      Error ("pdf_gen_from_array: xmin %g <= xmax %g\n", xmin, xmax);
      return (-1);
    }

/* Determine which jumps are important */
  njump_min = njump_max = 0;
  if (njumps > 0)
    {
      for (j = 1; j < njumps; j++)
	{
	  if (jump[j] <= jump[j - 1])
	    {
	      Error
		("pdf_gen_from_array: jump[%d]=%g <=jump[%d]=%g out of order\n");
	      return (-1);
	    }
	}
      njump_min = 0;
      while (njump_min < njumps && jump[njump_min] <= xmin)
	njump_min++;
      njump_max = 0;
      while (njump_max < njumps && jump[njump_max] < xmax)
	njump_max++;
      /* So at this point njump_min will point to the first jump which is betewen xmin and
         xmax or it will equal to njumps in which case there were no jumps which were greater
         than xmin.

         Similarly njump_max will be the point to the first jump above xmax or if there are
         no jumps above xmax, then it will be njumps. */

      njumps = njump_max - njump_min;
    }
/* Finished inital processing of jumps */

  allzero = 0;
  for (n = 0; n < n_xy; n++)
    {
      if (y[n] < 0)
	{
	  Error
	    ("pdf_gen_from_array: probability density %g < 0 at element %d\n",
	     y[n], n);
	  return (-1);
	}
      if (y[n] > 0)
	{
	  allzero = 1;
	};
    }
  /* OK, all the input data seems OK */



/* Shuffle x and y into pdf_xx and pdf_yy allowing for xmin and xmax */

  if (xmax < x[0] || xmin > x[n_xy - 1] || allzero == 0)
    {				// These are special (probably nonsensical) cases
      pdf_x[0] = xmin;
      pdf_z[0] = 0.;
      pdf_x[1] = xmax;
      pdf_z[0] = 1.;
      sum = 1.0;
      pdf_n = 2;
      Error
	("pdf_gen_from_array: all y's were zero or xmin xmax out of range of array x-- returning uniform distribution %d\n",
	 allzero);

    }
  else
    {
      m = 0;
      while (x[m] < xmin)
	m++;			// Find the bottom boundary
      pdf_x[0] = xmin;
      if (m == 0)
	{
	  pdf_y[0] = y[0];	//Assume prob. density is constant outside array lims.
	}
      else
	{
	  q = (xmin - x[m - 1]) / (x[m] - x[m - 1]);
	  pdf_y[0] = q * y[m] + (1. - q) * y[m - 1];
	}
      pdf_n = 1;
      // Completed first element; now do those that are completely in the grid
      while (x[m] < xmax && m < n_xy)
	{
	  pdf_x[pdf_n] = x[m];
	  pdf_y[pdf_n] = y[m];
	  m++;
	  pdf_n++;
	  if (pdf_n > PDF_ARRAY)
	    {
	      Error
		("pdf_gen_from_array: pdf_n (%d) exceeded maximum array size PDF_ARRAY (%d) \n",
		 pdf_n, PDF_ARRAY);
	      Error
		("pdf_gen_from_array: n_xy %d xmin %f xmax %f njumps %d\n",
		 n_xy, xmin, xmax, njumps);
	      Error
		("pdf_gen_from_array: Consider increasing PDF_ARRAY to a value > n_xy + njumps, and recompiling\n");
	      exit (0);
	    }
	}
      // Now worry about the last element
      pdf_x[pdf_n] = xmax;
      if (m < n_xy - 1)
	{
	  q = (xmax - x[m]) / (x[m + 1] - x[m]);
	  pdf_y[pdf_n] = q * y[m + 1] + (1. - q) * y[m];
	}
      else
	{
	  pdf_y[pdf_n] = y[m - 1];	//Again assume constant prob. density outside lims
	}
      pdf_n++;

/* So at this point, have probability density in pdf_x, pdf_y for the points
 * specified by the input array but we want the cumulative distribution
 */

      pdf_z[0] = 0.;
      for (n = 1; n < pdf_n; n++)
	{
	  pdf_z[n] =
	    pdf_z[n - 1] + 0.5 * (pdf_y[n - 1] + pdf_y[n]) * (pdf_x[n] -
							      pdf_x[n - 1]);
	}
      sum = pdf_z[pdf_n - 1];

      for (n = 1; n < pdf_n; n++)
	pdf_z[n] /= sum;
/* So pdf_z contains a properly normalized cdf on the points specified by
   the input array, or more explicitly, at the points specied in the array
   pdf_x
*/

    }

  /* From this we construct the cumulative distribution function on our more uniform
     grid.  04March -- ksl -- There is a problem that is sometimes appearing that
     appears to be due to the jumps.  The basic problem arises because of jumps
     and the fact that q is calculated differently for jumps than for non-jump
     situations.  The effect is a non-monotonic pdf in some situations where one
     has to deal with jumps.  To avoid this problem, I have added a check and 
     shuffled the pdf when it would otherwise be monotonic.  It is not obvious
     to me that this is the most elegant way to deal with this problem.  
   */



  pdf->x[0] = xmin;
  pdf->y[0] = 0;

  j = njump_min;
  m = 0;
  nn = 1;
  for (n = 1; n < NPDF; n++)
    {
      ysum = ((double) nn) / (NPDF - njumps);	/* This is the target with no jumps */

      while (pdf_z[m] < ysum)
	{
	  while (j < njump_max && jump[j] <= pdf_x[m])
	    {
	      pdf->x[n] = jump[j];
	      q = (jump[j] - pdf_x[m - 1]) / (pdf_x[m] - pdf_x[m - 1]);
	      pdf->y[n] = q * pdf_z[m] + (1. - q) * pdf_z[m - 1];
/* Note that pdf_gen_from_array only produces a term close to the desired break */
	      j++;		// increment the jump number

	      if (pdf->x[n] < pdf->x[n - 1])
		{		// Then we need to shuffle the pdf
		  xx = pdf->x[n - 1];
		  yy = pdf->y[n - 1];
		  pdf->y[n - 1] = pdf->y[n];
		  pdf->x[n - 1] = pdf->x[n];
		  pdf->x[n] = xx;
		  pdf->y[n] = yy;
		}

	      n++;		// now increment n
	    }

	  m++;			//increment m if necessary

	}

      q = (ysum - pdf_z[m - 1]) / (pdf_z[m] - pdf_z[m - 1]);
      pdf->x[n] = q * pdf_x[m] + (1. - q) * pdf_x[m - 1];
      pdf->y[n] = ysum;


      nn++;
    }

  pdf->x[NPDF] = xmax;
  pdf->y[NPDF] = 1.0;
  pdf->norm = sum;		/* The normalizing factor that would convert the function we
				   have been given into a proper probability density function */
/* Calculate the gradients */
  recalc_pdf_from_cdf (pdf);	// 57ib 

  if ((echeck = pdf_check (pdf)) != 0)
    {
      Error ("pdf_gen_from_array: error %d on pdf_check\n", echeck);
    }
  return (echeck);
}



/* 
pdf_get_rand


Our searching routine assumes we can predict roughly where in the
structure the right values will be 

History:

	06sep	ksl	57i -- Modified to account for the fact
			that the probability density is not 
			uniform between intervals.
*/

double
pdf_get_rand (pdf)
     PdfPtr pdf;
{
  double x, r;
  int i, j;
  double q;
  double a, b, c, s[2];
  int xquadratic ();


/* Find the interval within which x lies */
  r = rand () / MAXRAND;	/* r must be slightly less than 1 */
  i = r * NPDF;			/* so i initially lies between 0 and NPDF-1 */

  while (pdf->y[i + 1] < r && i < NPDF - 1)
    i++;
  while (pdf->y[i] > r && i > 0)
    i--;

/* Now calculate a place within that interval */

  q = rand () / MAXRAND;

  a = 0.5 * (pdf->d[i + 1] - pdf->d[i]);
  b = pdf->d[i];
  c = (-0.5) * (pdf->d[i + 1] + pdf->d[i]) * q;

  if ((j = xquadratic (a, b, c, s)) < 0)
    {
      Error ("pdf_get_rand: %d\n", j);
    }
  else
    {
      q = s[j];
      if (q < 0 || q > 1)
	{
	  Error ("pdf_get_rand: q out of range %d  %f\n", j, q);
	}
    }

  x = pdf->x[i] * (1. - q) + pdf->x[i + 1] * q;


  if (!(pdf->x[0] <= x && x <= pdf->x[NPDF]))
    {
      Error ("pdf_get_rand: %g %d %g %g\n", r, i, q, x);
    }
  return (x);

}



/* 
We have defined the pdf in such a way that one can generate a random number 0 and
1, find the place in the y array that corresponds to this and map back onto the 
x values (which is what is returned. We now want to limit the range of the returned
x values, but we need to know in terms of the original array where thiese values
are and the generate a random number between some minimum (>0) and maximum (<1) whichi
will allow one to find the x that corresponds. The point of this routine is to
find that minimum and maximum value. Having done this one can call pdf_get_rand_limit


Error Some of this is not quite right.  It's intended to let us generate limits in
the CDF but we need x too

	06sep	ksl	57i -- Added the xlimits to the pd
*/

int
pdf_limit (pdf, xmin, xmax)
     PdfPtr pdf;
     double xmin, xmax;
{
  int i;
  double q;
  if (pdf->y[NPDF] != 1.0)
    {
      Error ("pdf_limit: pdf not defined!)");
      exit (0);
    }
  if (xmin >= pdf->x[NPDF])
    {
      Error ("pdf_limit: xmin %g > pdf->x[NPDF] %g\n", xmin, pdf->x[NPDF]);
//      exit (0);
    }
  if (xmax <= pdf->x[0])
    {
      Error ("pdf_limit: xmax %g < pdf->x[0] %g\n", xmax, pdf->x[0]);
      exit (0);
    }

/* Set the limits for the minimum */

  if (xmin <= pdf->x[0])
    {
      pdf->limit1 = 0;
      pdf->x1 = pdf->x[0];
    }
  else
    {
      pdf->x1 = xmin;
      i = 0;
      while (xmin > pdf->x[i])
	{
	  i++;
	}
      q = (xmin - pdf->x[i - 1]) / (pdf->x[i] - pdf->x[i - 1]);
      pdf->limit1 = pdf->y[i - 1] + q * (pdf->y[i] - pdf->y[i - 1]);
    }

/* Now set the limits for the maximum */

  if (xmax >= pdf->x[NPDF])
    {
      pdf->limit2 = 1.0;
      pdf->x2 = pdf->x[NPDF];
    }
  else
    {
      pdf->x2 = xmax;
      i = NPDF;
      while (xmax <= pdf->x[i])
	{
	  i--;
	}
      q = (xmax - pdf->x[i]) / (pdf->x[i + 1] - pdf->x[i]);
      pdf->limit2 = pdf->y[i] + q * (pdf->y[i + 1] - pdf->y[i]);
    }

  return (0);
}

/*
pdf_get_rand_limit (pdf)

History
	06sep	ksl	57h -- Modified to account for the fact
			that the probability density is not 
			uniform between intervals.

*/
double
pdf_get_rand_limit (pdf)
     PdfPtr pdf;
{
  double x, r;
  int i, j;
  double q;
  double a, b, c, s[2];
  int xquadratic ();

  r = rand () / MAXRAND;	/* r must be slightly less than 1 */
  r = r * pdf->limit2 + (1. - r) * pdf->limit1;

  i = r * NPDF;

  while (pdf->y[i + 1] < r && i < NPDF - 1)
    i++;
  while (pdf->y[i] > r && i > 0)
    i--;

  while (TRUE)
    {
      q = rand () / MAXRAND;

      a = 0.5 * (pdf->d[i + 1] - pdf->d[i]);
      b = pdf->d[i];
      c = (-0.5) * (pdf->d[i + 1] + pdf->d[i]) * q;

      if ((j = xquadratic (a, b, c, s)) < 0)
	{
	  Error ("pdf_get_rand: %d\n", j);
	}
      else
	{
	  q = s[j];
	  if (q < 0 || q > 1)
	    {
	      Error ("pdf_get_rand: q out of range %d  %f\n", j, q);
	    }
	}

      x = pdf->x[i] * (1. - q) + pdf->x[i + 1] * q;

      if (pdf->x1 < x && x < pdf->x2)
	break;
    }

  return (x);

}

/* 

Write the full structure of the cumulattive distribution function to
a file

110628	ksl	Added lines to write out all of the information contained
		in the pdf structure in an attempt to debug an error in
		anisowind
*/

int
pdf_to_file (pdf, filename)
     PdfPtr pdf;
     char filename[];
{
  FILE *fopen (), *fptr;
  int n;

  fptr = fopen (filename, "w");

  fprintf (fptr, "# limits (portion.to.sample)   %10.4g %10.4g\n",
	   pdf->limit1, pdf->limit2);
  fprintf (fptr, "# x1 x2  Range(to.be.returned) %10.4g %10.4g\n", pdf->x1,
	   pdf->x2);
  fprintf (fptr, "# norm   Scale.factor          %10.4g \n", pdf->norm);

  fprintf (fptr, "#x y  1-y d\n");
  for (n = 0; n <= NPDF; n++)
    fprintf (fptr, "%10.4g	%14.8g %14.8e  %10.4g\n", pdf->x[n],
	     pdf->y[n], 1. - pdf->y[n], pdf->d[n]);

  fclose (fptr);
  return (0);
}

/* This is a simple routine to check some basics of the
   cumulative distribution function.

   History:
   98oct        ksl     Error returns improved so they can be interpreted.
   Also the requirement for a monotonic distribution
   function has been relaxed, in the sense that two
   successive values in x and y can be the same.  This was
   allowed because there was a case in which I had
   required values of the pdf at specific points but
   in the actual event, the partial distribution function
   between these two points was everywhere zero and hence
   the cumulative distribution function had two y values
   the same.  This can cause an error in pdf_get_rand
   in principle.
 */

int
pdf_check (pdf)
     PdfPtr pdf;
{
  int n;
  double x, y;
  int hcheck, icheck, jcheck, kcheck, fcheck;

  hcheck = icheck = jcheck = kcheck = fcheck = 0;

  x = pdf->x[0];
  y = pdf->y[0];
  if (y != 0.0)
    {
      Error
	("pdf_check: cumulative distribution function should start at 0 not %e\n",
	 y);
      hcheck = 1;
    }
  if (pdf->y[NPDF] != 1.0)
    {
      Error
	("pdf_check: cumulative distribution function should end at 1 not %e\n",
	 pdf->y[NPDF - 1]);
      icheck = 1;
    }

  for (n = 1; n < NPDF + 1; n++)
    {
      // Note the equal sign here 
      if (x <= pdf->x[n])
	x = pdf->x[n];
      else
	{
	  jcheck = 1;
	  Error ("pdf_check: x problem n %d x %f pdf->x[n] %f\n", n, x,
		 pdf->x[n]);
	}
      // Note the equal sign here 
      if (y <= pdf->y[n])
	y = pdf->y[n];
      else
	{
	  kcheck = 1;
	  Error ("pdf_check: y problem n %d y %f pdf->y[n] %f\n", n, y,
		 pdf->y[n]);
	}
    }
  if (jcheck == 1)
    {
      Error ("pdf_check: x values in pdf should be monotonic\n");
    }
  if (kcheck == 1)
    {
      Error
	("pdf_check: cumulative distribution function should be monotonic\n");
    }

  if (hcheck != 0 || icheck != 0 || jcheck != 0 || kcheck != 0)
    {
      pdf_to_file (pdf, "pdf.diag");
      fcheck = hcheck * 1000 + icheck * 100 + jcheck * 10 + kcheck;
      Error ("pdf_check %d %d %d %d %d\n", fcheck, hcheck, icheck, jcheck,
	     kcheck);
    }

/* Next section is a dangerous attempt to repair the pdf  */

  if (hcheck != 0)
    pdf->y[0] = 0.0;
  if (icheck != 0)
    pdf->y[NPDF] = 1.0;
  if (jcheck != 0)
    {
      for (n = 0; n < NPDF; n++)
	{
	  if (pdf->x[n] >= pdf->x[n + 1])
	    pdf->x[n + 1] = pdf->x[n] + 1.e-20;
	}
    }
  if (kcheck != 0)
    {
      for (n = 0; n < NPDF; n++)
	{
	  if (pdf->y[n] >= pdf->y[n] + 1)
	    pdf->y[n + 1] = pdf->y[n] + 1.e-20;
	}
    }

  return (fcheck);
}


/***********************************************************
                Space Telescope Science Institute
                                                                                             
Synopsis:
	recalc_pdf_from_cdf(pdf)
                                                                                             
Arguments:
                                                                                             
Returns:
                                                                                             
Description:

	Calculate the pdf form a cdf 
	allow one to calculate a first order correction
	to a uniform distibution between the points
	where the gradient has been calculated.
                                                                                             
Notes:

	d is calculated by differencing the cdf.  Thus
	this routine should work both for cdf's calculated
	from functions and from arrays
                                                                                             
History:
	06sep	ksl	57i -- Initiated to try an improve
			the calculation of the distribuion
			within a pdf interval.
	06nov	ksl	58b: Fixed problem occuring when there
			were two points in cdf with same x
                                                                                             
**************************************************************/

int
recalc_pdf_from_cdf (pdf)
     PdfPtr pdf;
{
  int n;
  double dx1, dx2, dy1, dy2;

  for (n = 1; n < NPDF; n++)
    {
      dy1 = pdf->y[n] - pdf->y[n - 1];
      dx1 = pdf->x[n] - pdf->x[n - 1];

      dy2 = pdf->y[n + 1] - pdf->y[n];
      dx2 = pdf->x[n + 1] - pdf->x[n];

      if (dx1 != 0.0 && dx2 != 0.0)
	{
	  pdf->d[n] = 0.5 * (dy1 / dx1 + dy2 / dx2);
	}
      else if (dx1 != 0.0)
	{
	  pdf->d[n] = dy1 / dx1;
	}
      else if (dx2 != 0.0)
	{
	  pdf->d[n] = dy2 / dx2;
	}
      else
	{
	  pdf->d[n] = 0.0;
	  Error ("recalc_pdf_from_cdf: dx1 and dx2 both 0 at  %d %f\n", n,
		 pdf->x[n]);
	}

    }
  /* Fill in the ends */
  pdf->d[0] = pdf->d[1];
  pdf->d[NPDF] = pdf->d[NPDF - 1];

  return (0);
}
