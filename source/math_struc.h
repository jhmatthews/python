
#define NCDF 30000              //The default size for these arrays.  This needs to be greater than
                                //the size of any model that is read in, hence larger than NWAVE_EXTRACT in models.h
#define FUNC_CDF  2000          //The size for CDFs made from functional form CDFs
#define ARRAY_PDF 1000          //The size for PDFs to be turned into CDFs from arrays


/**
* This is the structure for storing cumulative distribution functions. The CDFs are
* generated from a function which is usually only proportional to the probability density
* function or from an array.  It is sometimes useful, e.g. in calculating the reweighting function to
* have access to the proper normalization.  
*
* x1 and x2 and limit1 and limit2 are set by cdf_limit, when one want to generate 
* random numbers with in a range of a pre-existing cdf.  x1 and x2 define the limits in
* x, and limit1 and limit2 define the limits in y.  Assuming these have been defined,
* then one uses cdf_gen_rand_limit to generate numbers between these two limits.
* limit1 and limit2 can be defined see cdf_limit and cdf_get_rand_limit so that one can
* use a cdf to generate numbers between these two limits.
*/
typedef struct Cdf
{
  double x[NCDF]; /**< Positions for which the CDF is calculated */
  double y[NCDF]; /**< The value of the CDF at x */
  double d[NCDF]; /**< The rate of change of the CDF at x */
  double limit1;  /**< LowerLimit (running from 0 to 1) that define a portion of the CDF to sample */
  double limit2;  /**< Upper Limit (running from 0 to 1) that define a portion of the CDF to sample */
  double x1;      /**< Lower limit if it exists on what is returned */
  double x2;      /**< Uppper limit if existst on what is returned */
  double norm;    /**< The scaling factor which would renormalize the CDF */
  int ncdf;       /**<  Size of this CDF */
}
 *CdfPtr, cdf_dummy;


 /**
   * A structure which defines a rotation matrix for
   * corrdinage sytstem transformations.
   */
struct basis
{
  double a[3][3];  /**< the rotation matrix */

};

extern struct Cdf cdf_vcos;
extern struct Cdf cdf_vdipole;



#define LINELENGTH     256
