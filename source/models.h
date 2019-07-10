

/* Note that there are some limits on how large NMODS*NWAVES can be.  When you exceed this limit, you will get a segmentation
 * fault on startup */

/* Note that there is some overlap between these definitions and those in python.h and so mixing these two
 * header files can be dangerous.  
 *
 * NDIM --> MOD_NDIM  -- 081026
 */

 /* 1405 JM -- Increased PDF array for use with disk14 models- also removed duplication of ncomps */
/* 1409 JM -- Increased LINELEN to 160 */

#define NWAVES  28000           // The maximum number of wavelength bins in models
#define NDIM	10              // The maximum number of free parameters
#define NCOMPS	10              //The maximum number of separate components
#define NPARS	10              //The maximum number of parameters in a component (not all variable)
#define NMODS   2000            //The total number of models read in in all components
#define LINELEN 160             //This should always be the same as LINELENGTH in python.h!


//#include      "pdf.h"

int ncomps;                     // The number of components that have been read.  
int nmods_tot;                  // The total number of models that have been read in 


/* This is the structure that describes an individual continuum model. 
 * mods is the set of all models that are read
 */
struct Model
{
  char name[LINELEN];
  double par[NPARS];
  double w[NWAVES];
  double f[NWAVES];
  int nwaves;
}
mods[NMODS];

/* There is one element of comp for each set of models of the same type, i.e. if
one reads in a list of WD atmosphers this will occupy one componenet here */

struct ModSum
{
  char name[LINELEN];
  int npars;                    // Number of parameters that describe the model
  int nmods, modstart, modstop; //number of models of this type and the index into the model struct
  double min[NPARS];            // The minimum and maximum for each "free" parameter of the model
  double max[NPARS];
  int nwaves;                   //All models in each comp should have same wavelengths;
  struct Model xmod;            //The current intepolated model of this type 
  struct Cdf xcdf;              //The current cumulative distribution function for this component
}
comp[NCOMPS];




/* In modsum, current[0] often refers to a normalization.  Therefore for parallism, a set of
models which really only have one parameter, e.g. T, will store that parmeter in .par[1] */
