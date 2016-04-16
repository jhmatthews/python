#ifdef MPI_ON
#include "mpi.h"
#endif 

int np_mpi_global;               /// Global variable which holds the number of MPI processes

int rank_global; 

// DEBUG is deprecated, see #111, #120
//#define DEBUG 				0	/* 0 means do not debug */
int verbosity;			/* verbosity level. 0 low, 10 is high */

/* the functions contained in log., rdpar.c and lineio.c are
   declare deparately from templates. This is because some functions
   only use log.h and don't use python.h due to repeated definitions */
#include "log.h"

/* In python_43 the assignment of the WindPtr size has been moved from a fixed
value determined by values in python.h to a values which are adjustable from
within python */
int ndim;			// Define the fundamental dimension of the grid
int mdim;
int NDIM, MDIM, NDIM2;
int NPLASMA;			//The number of cells with non-zero volume or the size of plasma structure

char basename[132];		// The root of the parameter file name being used by python

/* These are tunable parameters that control various aspects of python
 * and the assorted programs.  In some cases they affect the "care" with
 * which a calculation is made, and this can in principle affect the
 * speed of the program.  One day if may be desirable to change some
 * of these parameters with time in the program.  At present they are
 * simply collected here
 * 
 * */

double dfudge;			// This is the push-through distance
double DFUDGE;
#define VCHECK	1.e6		// The maximum allowable error in calculation of the velocity in calculate_ds




/* 57h -- Changed several defined variables to numbers to allow one to vary them 
in the process of running the code */
double SMAX_FRAC;		/* In translate_in_wind, a limit is placed on the maximum distance a
				   photon can travel in one step.  It is a fraction SMAX_FRAC of the
				   distance of the photon from the origin.  This had been hardwired to
				   0.1 for up to 57h.  Changing it to 0.5 speeds up the current version
				   of the code by as much as a factor of 2 for small sized grids.  This
				   had been introduced as part of the attempt to assure ourselves that
				   line shapes were calculated as accurately as possilble. 
				 */
double DENSITY_PHOT_MIN;	/* This constant is a minimum density for the purpose of calculating
				   photoionization heating and recombination cooling.  It is important that heating and cooling
				   be calculated self-consistently.  Program speed is somewhat sensitive 
				   to this parameter, at the 10% level if raised from 1e-3 to 1.  There is a 
				   trade-off since lower minima may give better results, especially for macro atoms. */

//#define SMAX_FRAC     0.1  
#define LDEN_MIN        1e-3	/* The minimum density required for a line to be conidered for scattering
				   or emission in calculate_ds and lum_lines */


/* End of "care factor" definition */


#define RMIN   				1.e9
#define RMAX   				3.e10
#define VWIND  				2.e8
#define MDOT  				1.e-9*MSOL/YR
#define TSTAR 				30000.	/*Sets a floor on the frequency range used to calculate
						   ionization balance in the wind.  Only used in python.c.  May 
						   be superfluous */
#define BETA  				1.0
#define KAPPA_CONT 			4.
#define EPSILON  			1.e-6	/* A general purpose fairly small number */
#define NSTAT 				10    // JM increased this to ten to allow for adiabatic
#define VMAX                		1.e9
#define TAU_MAX				20.	/* Sets an upper limit in extract on when
						   a photon can be assumed to be completely absorbed */

//#define SELECT_NEBULAR                1 /*non zero means to use the nebular approximation; 0 implies
//                                                              use LTE populations based on t_rad*/
#define DANG_LIVE_OR_DIE   2.0	/* If constructing photons from a live or die run of the code, the
				   angle over which photons will be accepted must be defined */


//#define NPHOT                                 10000
int NPHOT;			/* As of python_40, NPHOT must be defined in the main program using
				   python.h */

#define NWAVE  			       10000	//Increasing from 4000 to 10000 (SS June 04)
#define MAXSCAT 			50

/* Define the structures */


/* Geometry is an actual structure which exists at the beginning of th program
   It carries the variables which define the geometry.  Reasonable values of each of
   these should be defined before it is altered with inputs fromt eh terminal. 
   The geometry structure is used to transfer all of the information about a wind


 */

/* Definitions of spectral types, which are all negative because when
 * one reads a spectrum from a list of models these are numbered beginning
 * with zero, see the discussion in get_models.c   080518 - ksl - 60a
 */
#define SPECTYPE_BB      -1
#define SPECTYPE_UNIFORM -2
#define SPECTYPE_POW     -4
#define SPECTYPE_CL_TAB  -5
#define SPECTYPE_BREM    -6
#define SPECTYPE_NONE	 -3

/* Number of model_lists that one can have, should be the same as NCOMPS in models.h */
#define NCOMPS 	10
#define LINELENGTH 	160

struct geometry
{
/* 67 - ksl This section added to allow for restarting the program, and adds parameters used
 * in the calculation */

  int wcycle, pcycle;		/* The number of completed ionization and spectrum cycles */
  int wcycles, pcycles;		/* The number of ionization and spectrum cycles desired */

/* Begin description of the actual geometery */
  enum coord_type_enum
  	{	SPHERICAL 	=0,
		CYLIND 		=1,
		RTHETA		=2,
		CYLVAR      =3
  	} coord_type;

  int ndim, mdim;	/* The type of geometry and dimensionality of the wind array. 
				   0=1-d spherical, 1=cylindrical, 2 = spherical polar, 3=cylindrical
				   but the z coordinate changes with rho in an attempt to allow for
				   a vertically extended disk....
				   ndim is the dimensionality of the first dimension.  In the CV case
				   it is in the plane of the disk. Mdim is generally along the z axis
				 */
  int nplasma, nmacro;		/*The number of cells in the plasma and macro structures 08mar ksl */
  double rmax, rmax_sq;		/* The maximum distance to which a photon should be followed */
  double mstar, rstar, rstar_sq, tstar, gstar;	/* Basic parameters for the WD */
  double twind;			/* temperature of wind */
  double tmax;			/*NSH 120817 the maximim temperature of any element of the model - used to help estimate things for an exponential representation of the spectrum in a cell */
  enum system_type_enum 
  	{	SYSTEM_TYPE_STAR	  =0, 
  		SYSTEM_TYPE_BINARY	=1, 
  		SYSTEM_TYPE_AGN		  =2
  	} 	system_type;
  int disk_type;		/*0 --> no disk, 
				   1 --> a standard disk in xy plane, 
				   2 --> a vertically extended disk 
				   Note that this definition is new to Python52 */
  int disk_illum;		/*Treatment of effects of illumination on the disk. 
				   0--> Disk simply absorbs the radiation and it is lost
				   1--> Disk reradiates the radiation immediately via electron scattering
				   2--> Disk radiation is absorbed and changes the temperature of the disk for
				   future ionization cycles
				   3--> Disk illumination is treated in terms of an analytic approximation
				   04Aug ksl -- this parameter added for Python52
				 */
  int disk_tprofile;
  double disk_mdot;		/* mdot of  DISK */
  double diskrad, diskrad_sq;
  double disk_z0, disk_z1;	/* For vertically extended disk, z=disk_z0*(r/diskrad)**disk_z1 */
  int wind_type;		/*Basic prescription for wind(0=SV,1=speherical , 2 can imply old file */
  int log_linear;		/*0 -> the grid spacing will be logarithmic in x and z, 1-> linear */
  double xlog_scale, zlog_scale;	/* Scale factors for setting up a logarithmic grid, the [1,1] cell
					   will be located at xlog_scale,zlog_scale */
  int star_radiation, disk_radiation;	/* 1 means consider radiation from star, disk,  bl, and/or wind */
  int bl_radiation, wind_radiation, agn_radiation;
  int matom_radiation;		/* Added by SS Jun 2004: for use in macro atom computations of detailed spectra
				   - 1 means use emissivities for BOTH macro atom levels and kpkts. 0 means don't
				   (which is correct for the ionization cycles. */
  int ioniz_mode;		/* describes the type of ionization calculation which will
				   be carried out.  0=on the spot, 1=LTE, 2=fixed ionization
				   fractions,  3 means to recalculate the ionization structure 
				   based on the energy absorbed in the wind (mod_on_the_spot), 4
				   is a test.  It is currently set to do the same as 3, except
				   that ground state mulitpliciites are used instead of 
				   a partition function */
  int macro_ioniz_mode;		/* Added by SS Apr04 to control the use of macro atom populations and
				   ionization fractions. If it is set to 1 then macro atom populations
				   computed from estimators are used. If set to 0 then the macro atom
				   populations are computed as for minor ions. By default it is set to
				   0 initially in python.c and then set to 1 the first time that
				   Monte Carlo estimators are normalised. */
  int ioniz_or_extract;		/* Set to 1 (true) during ionization cycles, set to 0 (false) during calculation of
				   detailed spectrum.  Originally introduced by SS in July04 as he added
				   macro atoms.  Name changed by ksl (57h) since this variable can be used more
				   generally to speed up the extract portion of the calculation.
				 */
  int macro_simple;		/* Added by SS May04 for diagnostics. As default this is set to 0. A full
				   Macro Atom calculation is performed in that case. If it is set to 1 it means
				   that although Macro Atom data has been read in, all lines/continua are treated
				   using the simplified two-level approximation. Such a calculation should reproduce
				   the same results as pre-Macro Atom versions of the code. */
  int partition_mode;		/* Diagnostic to force the partition function to be calculated in
				   a specific way. */
  int line_mode;		/*0, the atomosphere is a completely absorbing and no photons
				   will be scattered.  In this mode, assuming the wind is a source
				   of emission, the emissivity will be the Einstein A coefficient
				   1, the atmosphere is completely scattering, there will be no
				   interchange of energy between the photons and the electrons
				   as a result of readiation transfer
				   2, then a simple single scattering approximation is applied in which
				   case the scattered flux is just  A21/(C21+A21). 
				   3, then radiation trapping is included as well.
				   6, If set to 6 initially, this switches on the macro atom stuff
				   and then behaves like 3. (SS)
				 */
  int scatter_mode;		/*The way in which scattering for resonance lines is treated 
				   0  isotropic
				   1  anisotropic
				   2  thermally broadened anisotropic
				 */
  int rt_mode;			/* radiative transfer mode (0=Sobolev,1=simple (used only by balance) */
  /* This IS now used by Python - set to 2 for Macro Atom method. Set to 1
     for non-Macro Atom methods (SS) */

  /* 71 - 111229  - ksl - These are the frequency bands used when calculating parameters like a power law slope
   * in limited regions.  Moved inside geo becuase we need to know this information in py_wind
   */
#define  NXBANDS 20		/* the maximum number of bands that can be defined */
  int nxfreq;			/* the number of bands actually used */
  double xfreq[NXBANDS + 1];	/* the band limits  */


  /* The spectral types are SPECTYPE_BB for bb, SPECTYPE_UNIFORM for a uniform spectral distribution, 
   * SPECTYPE_POW for a power law, 0 or more from a filelist.
   * A value of SPECTYPE_NONE indicates no emission is expected from this particular source */
  int star_ion_spectype, star_spectype;	/* The type of spectrum used to create the continuum
					   for the star in the ionization and final spectrum calculation */
  int disk_ion_spectype, disk_spectype;	/* The type of spectrum used to create the continuum
					   for the disk in the ionization and final spectrum calculation */
  int bl_ion_spectype, bl_spectype;	/* The type of spectrum used to create the continuum
					   for the bl in the ionization and final spectrum calculation */
  int agn_ion_spectype, agn_spectype;	/* The type of spectrum used to create the continuum
					   for the agn in the ionization and final spectrum calculation */
  char model_list[NCOMPS][LINELENGTH];	/* The file which contains the model names and the associated values for the model */

  /* Generic parameters for the wind */
  double wind_mdot, stellar_wind_mdot;	/* Mass loss rate in disk and stellar wind */
  double wind_rmin, wind_rmax;	/*Spherical extent of the wind */
  double wind_rho_min, wind_rho_max;	/*Min/Max rho for wind in disk plane */
  double wind_thetamin, wind_thetamax;	/*Angles defining inner and outer cones of wind, measured from disk plane */
  double mdot_norm;		/*A normalization factor used in SV wind, and Knigge wind */
  int adiabatic;		/*0-> Do not include adiabatic heating in calculating the cooling of the wind
				   1-> Use adiabatic heating in calculating the cooling of the wind
				 */
  int auger_ionization;		/*0 -> Do not include innershell photoionization /Auger effects; 1-> include them */
  /* Parameters defining Shlossman & Vitello Wind */
  double sv_lambda;		/* power law exponent describing from  what portion of disk wind is radiated */
  double sv_rmin, sv_rmax, sv_thetamin, sv_thetamax, sv_gamma;	/* parameters defining the goemetry of the wind */
  double sv_v_zero;		/* velocity at base of wind */
  double sv_r_scale, sv_alpha;	/* the scale length and power law exponent for the velocity law */
  double sv_v_infinity;		/* the factor by which the velocity at infinity exceeds the excape velocity */

  /* Paramater for the Elvis AGN wind - closely based on SV */
  double elvis_offset;		/*This is a vertical offset for a region where the
				   wind rises vertically from the disk */

  /* Parameters defining Knigge Wind */
  double kn_dratio;		/* parameter describing collimation of wind */
  double kn_lambda;		/* power law exponent describing from  what portion of disk wind is radiated */
//    double kn_rmin, kn_rmax, kn_thetamin, kn_thetamax, kn_gamma;      /* parameters defining the goemetry of the wind */
//    double kn_v_zero;         /* velocity at base of wind */
  double kn_r_scale, kn_alpha;	/* the scale length and power law exponent for the velocity law */
  double kn_v_infinity;		/* the factor by which the velocity at infinity exceeds the excape velocity */
  double kn_v_zero;		/* NSH 19/04/11 - Added in as the multiple of the sound speed to use as the initial velocity */

  /* Parameters describing Castor and Larmors spherical wind */
  double cl_v_zero, cl_v_infinity, cl_beta;	/* Power law exponent */
  double cl_rmin, cl_rmax;

  /* Parameters describing a spherical shell test wind */
  double shell_vmin, shell_vmax, shell_beta;
  double shell_rmin, shell_rmax;

  /*Parameters defining a corona in a ring above a disk */
  double corona_rmin, corona_rmax;	//the minimum and maximu radius of the corona

  double corona_base_density, corona_scale_height;	//the density at the base of the corona and the scale height

  double corona_vel_frac;	// the radial velocity of the corona in units of the keplerian velocity

/* The filling factior for the wind or corona */
  double fill;
/* Initial values for defining wind structure for a planar geometry.  These are currently only used by balance and this
   may not be the best approach generally and depending on where this ends up. Some consolidation is desirable */
  double pl_vol, pl_vmax;
  double pl_t_r, pl_t_e, pl_w;
  double pl_nh;

  double lum_tot, lum_star, lum_disk, lum_bl, lum_wind;	/* The total luminosities of the disk, star, bl, & wind 
							   are actually not used in a fundamental way in the program */
  double lum_agn;		/*The total luminosity of the AGN or point source at the center */

/* The next four variables added by nsh Apr 2012 to allow broken power law to match the cloudy table command */
  double agn_cltab_low;		//break at which the low frequency power law ends
  double agn_cltab_hi;		//break at which the high frequency power law cuts in
  double agn_cltab_low_alpha;	//photon index for the low frequency end
  double agn_cltab_hi_alpha;	//photon index for the high frequency end	


  double lum_ff, lum_fb, lum_lines;	/* The luminosity of the wind as a result of ff, fb, and line radiation */
  double lum_comp;		/*1108 NSH The luminosity of the wind as a result of compton cooling */
  double lum_di;		/* 1409 NSH The direct ionization luminosity */
  double lum_dr;		/*1109 NSH The luminosity of the wind due to dielectronic recombination */
  double lum_adiabatic;		/*1209 NSH The cooling of the wind due to adiabatic expansion */
  double heat_adiabatic;		/*1307 NSH The heating of the wind due to adiabatic heating - split out from lum_adiabatic to get an accurate idea of whether it is important */
  double f_tot, f_star, f_disk, f_bl, f_agn, f_wind;	/* The integrated specific L between a freq min and max which are
							   used to establish the fraction of photons of various types */

/* These variables are copies of the lum variables above, and are only calculated during ionization cycles
   This is a bugfix for JM130621, windsave bug */
  double lum_ff_ioniz, lum_fb_ioniz, lum_lines_ioniz;	
  double lum_comp_ioniz;
  double lum_di_ioniz;		/* 1409 NSH The direct ionization luminosity */		
  double lum_dr_ioniz;		
  double lum_adiabatic_ioniz;	
  double lum_wind_ioniz, lum_star_ioniz, lum_disk_ioniz, lum_bl_ioniz, lum_tot_ioniz;



  double f_matom, f_kpkt;	/*Added by SS Jun 2004 - to be used in computations of detailed spectra - the
				   energy emitted in the band via k-packets and macro atoms respectively. */

// The next set of parameters relate to the secondary
  double m_sec, q;		/* Mass of the secondary, mass ratio of system */
  double period;		/* Period of the systems in seconds */
  double a, l1, l2, phi;	/* Separation of primary and secondary, distance of l1 from primary,phi at l1 */
  double l1_from_m2, r2_far;	/* Distance to l1 from m2, distance to far side of secondary from primary */
  double r2_width;		/* Maximum width of Roche surface of secondary in the plane of orbit */

  double t_bl;			/*temperature of the boundary layer */
  double weight;		/*weight factor for photons/defined in define_phot */

// The next set of parameters relate to the central source of an AGN
  double brem_temp;       /*The temperature of a bremsstrahlung source */
  double brem_alpha;       /*The exponent of the nu term for a bremstrahlung source */

  double pl_low_cutoff;  /* accessible only in advanced mode- see #34. default to zero */

  double alpha_agn;		/*The power law index of a BH at the center of an AGN.  Note that the luminosity
				   of the agn is elsewhere in the structure
				 */
  double const_agn;		/*The constant for the Power law, there are lots of ways of defining the PL which is best? */
  double r_agn;			/* radius of the "photosphere" of the BH in the AGN.  */
  double d_agn;			/* the distance to the agn - only used in balance to calculate the ioinsation fraction */

// 70b - ksl - 110809  The next set of sources relate to a compton torus that is initially at least just related to AGN

  int compton_torus;		/* 0 if there is no Compton torus; 1 otherwise */
  double compton_torus_rmin;	/* The minimum radius of the torus */
  double compton_torus_rmax;	/* The maximum radius of the torus */
  double compton_torus_zheight;	/* The height of the torus. */
  double compton_torus_tau;	/* The optical depth through the torus at the height. */
  double compton_torus_te;	/* The initial temperature of the torus */

//70i - nsh 111007 - put lum_ioniz and n_ioniz into the geo structure. This will allow a simple estimate of ionisation parameter to be computed;

  double n_ioniz, lum_ioniz;

// The next set of parameters describe the input datafiles that are read
  char atomic_filename[132];	/* 54e -- The masterfile for the atomic data */
  char fixed_con_file[132];	/* 54e -- For fixed concentrations, the file specifying concentrations */

  //Added by SWM for reverberation mapping
  enum reverb_enum {
    REV_NONE=0, REV_PHOTON=1, REV_WIND=2
  } reverb;
  int reverb_path_bins, reverb_theta_bins; 
}
geo;


typedef struct plane /*SWM 10-10-14 - Switched to TypeDef */
{
  double x[3];			/* A position included in the plane (usally the "center" */
  double lmn[3];		/* A unit vector perpendicular to the plane (usually in the "positive" direction */
} plane_dummy, *PlanePtr;
plane_dummy plane_l1, plane_sec, plane_m2_far;	/* these all define planes which are perpendicular to the line of sight from the 
					   primary to the seconday */


/* Note that since we are interested in biconical flows, our definition of a cone is not exactly
 * what one might guess.  The cone is defined in the positive z direction but reflected through 
 * the xy plane.  
 * 56d -- Beginning with 56d, ksl has switched to a new definition of cones, that is intended to
 * make it possible to use ds_to_cone easier as part of different coordinate systems.  The new definition
 * is based on the intersection of the cone with the z axis rather than the intersection with
 * the disk plane.  At present both definitions are used in the program and therefore both shold
 * be defined.  Once the new definition is promulgated through the entire program, and verified
 * the old definitions can be elimiated.  05jul -- ksl
 */

typedef struct cone
{
  double z;			/* The place where the cone intersects the z axis (used after 56d) */
  double dzdr;			/* the slope (used after 56d) */
}
cone_dummy, *ConePtr;

ConePtr cones_rtheta;		/*A ptr to the cones that define the theta directions in rtheta coods */

struct cone windcone[2];	/* The cones that define the boundary of winds like SV or kwd */



#define NRINGS	301		/* The actual number of rings completely defined
				   is NRINGS-1 ... or from 0 to NRINGS-2.  This is
				   because you need an outer radius...but the rest
				   of this element is not filled in. */

struct xdisk
{
  double r[NRINGS];		/* The inner radius of an annulus */
  double t[NRINGS];		/* The temperature at the middle of the annulus */
  double g[NRINGS];		/* The gravity at the middle of the annulus */
  double v[NRINGS];		/* The velocity at the middle of the annulus */
  double heat[NRINGS];		/* The total energy flux of photons hitting each annulus */
  double ave_freq[NRINGS];	/* The flux weighted average of frequency of photons hiiting each annulus */
  double w[NRINGS];		/* The radiative weight of the photons that hit the disk */
  double t_hit[NRINGS];		/* The effective T of photons hitting the disk */
  int nphot[NRINGS];		/*The number of photons created in each annulus */
  int nhit[NRINGS];		/*The number of photons which hit each annulus */
}
disk, qdisk;			/* disk defines zones in the disk which in a specified frequency band emit equal amounts
				   of radiation. qdisk stores the amount of heating of the disk as a result of
				   illumination by the star or wind */

#define NBLMODEL 100

struct blmodel
{
  int n_blpts;
  float r[NBLMODEL];
  float t[NBLMODEL];
}
blmod;


/*
    SWN 6-2-15
    Wind paths is defined per cell and contains a binned array holding the spectrum of paths. Layers are
    For each frequency:
      For each path bin:
        What's the total fluxback of all these photons entering the cell?
*/
typedef struct wind_paths
{
  double* ad_freq_path_flux;  //Array[by frequency, then path] of total flux of photons with the given v&p
  int*    ai_freq_path_num;   //Array[by frequency, then path] of the number of photons in this bin
  double* ad_freq_flux;       //Array[by frequency] of total flux of photons with the given v
  int*    ai_freq_num;        //Array[by frequency] of the number of photons in this bin
  double  d_flux, d_path;     //Total flux, average path
  int     i_num;              //Number of photons hitting this cell
} wind_paths_dummy, *Wind_Paths_Ptr;

typedef struct path_data
{
  double* ad_path_bin;              //Array of bins for the path histograms
  double* ad_freq_bin;              //Array of bins for the frequency histograms
  int     i_path_bins, i_obs;       //Number of bins, number of observers
  int     i_theta_res;              //Number of angular bins when outputting observer paths
} path_data_dummy, *Path_Data_Ptr;
Path_Data_Ptr path_data;
Path_Data_Ptr g_path_data;

/* 	This structure defines the wind.  The structure w is allocated in the main
	routine.  The total size of the structure will be NDIM x MDIM, and the two
	dimenssions do not need to be the same.  The order of the
	cells in the structure is such that the as you increse the cell number by one
	z increases the fastest.

57+ -- The wind structure was reduced to contain just information about the geometry.  
Variables for wind cells that actually have volume in the wind are now in plasmamain, 
or macromain.  The wind structure still contains a volume element, which is the volume
of that cell in the wind, This is used in some cases to determine whether the cell
has any portion in the wind.  

Note that most of the macro atom information is in a separate structure.  This was
done to make it easier to control the size of the entire structure   06jul-ksl

 */
#define NIONIZ	5		/*The number of ions (normally H and He) for which one separately tracks ionization 
				   and recombinations */
#define LPDF   3		/*The number of bins into which the line luminosity is divided in the course pdf
				   created by total_line_emission 
				   Reduced from 10 to 3 by SS for testing data with few lines. */


/* 061104 -- 58b -- ksl -- Added definitions to characterize whether a cell is in the wind. */
/* 110810 -- ksl - these are assigned to w->inwind, and are used to help determine how photons 
that go through a cell are treated.  Most of the assignments are based on whether the
volume in the wind is zero, which is determined in cylind_volumes for the cylindrical wind
and in correpsonding reoutines elsewehere.  These routines are called from the routine define_wind.
W_IGNORE is currently set in define_wind itself.  The values of these variables are used
in translate_in_wind.

Note that where_in_wind has been modified to use some of the same returns.  In general, the idea
is that if a new component is to be added, it should be added with by with two varialles ALL in whatever
and PART in whatever, as n and n+1
*/

typedef struct wind
{
  int nwind;			/*A self-reference to this cell in the wind structure */
  int nplasma;			/*A cross refrence to the corresponding cell in the plasma structure */
  double x[3];			/*position of inner vertex of cell */
  double xcen[3];		/*position of the "center" of a cell (Added by ksl for 52a--04Aug) */
  double r, rcen;		/*radial location of cell (Used for spherical, spherical polar
				   coordinates. (Added by ksl for 52a --04Aug) */
  double theta, thetacen;	/*Angle of coordinate from z axis (Added by ksl for 52a -- 04Aug) */
  double dtheta,dr;    /* widths of bins, used in hydro import mode*/
  struct cone wcone;		/*56d -- cone structure that defines the bottom edge of the cell in 
				   CYLVAR coordinates */
  double v[3];			/*velocity at inner vertex of cell.  For 2d coordinate systems this
				   is defined in the xz plane */
  double v_grad[3][3];		/*velocity gradient tensor  at the inner vertex of the cell NEW */
  double div_v;			/*Divergence of v at center of cell */
  double dvds_ave;		/* Average value of dvds */
  double dvds_max, lmn[3];	/*The maximum value of dvds, and the direction in a cell in cylindrical coords */
  double vol;			/* valid volume of this cell (that is the volume of the cell that is considered
				   to be in the wind.  This differs from the volume in the Plasma structure
				   where the volume is the volume that is actually filled with material. */
  enum inwind_enum
  	{	W_IGNORE=-2, 	W_NOT_INWIND=-1, 
  		W_ALL_INWIND=0, W_PART_INWIND=1, 
  		W_ALL_INTORUS=2, W_PART_INTORUS=3
  	}	inwind;			
  Wind_Paths_Ptr paths;         // SWM 6-2-15 Path data struct for each cell
}
wind_dummy, *WindPtr;

WindPtr wmain;

/* 57+ - 06jun -- plasma is a new structure that contains information about the properties of the
plasma in regions of the geometry that are actually included n the wind 

	07jul	ksl	Added volume to the structure.  The value of this should be the same
			as the corresponding volume in the Wind structure, as we are still
			using the Wind volume for some tests of whether the a photon can
			interact in a cell.
*/

/* 70 - 1108 - Define wavelengths in which to record gross spectrum in a cell, see also xave_freq and xj in plasma structure */
/* ksl - It would probably make more sense to define these in the same ways that bands are done for the generation of photons, or to
 * do both the same way at least.  Choosing to do this in two different ways makes the program confusing. The structure that has
 * the photon generation is called xband */
/* 78 - 1407 - NSH - changed several elements (initially those of size nions) in the plasma array to by dynamically allocated.
They are now pointers in the array. */


typedef struct plasma
{
  int nwind;			/*A cross reference to the corresponding cell in the  wind structure */
  int nplasma;			/*A self reference to this  in the plasma structure */
  double ne;			/* electron density in the shell */
  double rho;			/*density at the center of the cell. For clumped models, this is rho of the clump */
  double vol;			/* volume of this cell (more specifically the volume  that is filled with material
				   which can differe from the valid volume of the cell due to clumping. */
  double *density;	/*The number density of a specific ion.  This needs to correspond
				   to the ion order obtained by get_atomic_data. 78 - changed to dynamic allocation*/
  double *partition;	/*The partition function for each  ion. 78 - changed to dynamic allocation */
  double levden[NLTE_LEVELS];	/*The number density (occupation number?) of a specific level */

  double *PWdenom;	/*The denominator in the pairwise ionization solver. Sicne this is computed at a temperature 
				   chosen from basic ioinzation proerties to be good for this ion, it should not change
				   very much from cycle to cycle - hence we shold be able to speed up the code by storing 
				   it and refering to it if the temperature has not changed much. 78 - changed to dynamic allocation */
  double *PWdtemp;	/*The temperature at which the pairwise denominator was calculated last. 78 - changed to dynamic allocation */
  double *PWnumer;	/* The numberator in the pairwise approach. When we are chasing the true density
				   by carying n_e - this value will not change, so we canspeed things up a lot
				   by not recomputing it!. 78 - changed to dynamic allocation */
  double *PWntemp;	/* The temperature at which the stored pairwise numerator was last computed at. This
				   is used in the BB version of the pairwise correction factor. 78 - changed to dynamic allocation */

  double kappa_ff_factor;	/* Multiplicative factor for calculating the FF heating for                                      a photon. */

  /* Two new objects in the structure to hol the number of line resonant scatters and the number of electron scatters */

  int nscat_es;
  int nscat_res;



  double recomb_simple[NTOP_PHOT];	/* "alpha_e - alpha" (in Leon's notation) for b-f processes in simple atoms. */

/* Begining of macro information */
  double kpkt_emiss;		/*This is the specific emissivity due to the conversion k-packet -> r-packet in the cell
				   in the frequency range that is required for the final spectral synthesis. (SS) */

  double kpkt_abs;		/* k-packet equivalent of matom_abs. (SS) */

  int kbf_use[NTOP_PHOT];	/* List of the indices of the photoionization processes to be used for kappa_bf. (SS) */
  int kbf_nuse;			/* Total number of photoionization processes to be used for kappa_bf. (SS) */

/* End of macro information */


  double t_r, t_r_old;		/*radiation temperature of cell */
  double t_e, t_e_old;		/*electron temperature of cell */
  double dt_e, dt_e_old;	/*How much t_e changed in the previous iteration */
  double heat_tot, heat_tot_old;	/* heating from all sources */
  double heat_lines, heat_ff;
  double heat_comp;		/* 1108 NSH The compton heating for the cell */
  double heat_ind_comp;		/* 1205 NSH The induced compton heatingfor the cell */
  double heat_lines_macro, heat_photo_macro;	/* bb and bf heating due to macro atoms. Subset of heat_lines 
						   and heat_photo. SS June 04. */
  double heat_photo, heat_z;	/*photoionization heating total and of metals */
  double heat_auger;       /* photoionization heating due to inner shell ionizations */  
  double w;			/*The dilution factor of the wind */
  int ntot;			/*Total number of photon passages */

  /* NSH 15/4/11 - added some counters to give a rough idea of where photons from various sources are ending up */
  /* NSH 111005  - changed counters to real variables, that allows us to take account of differening weights of photons */
  /* ksl - ???? The reason these are doubles if that what Nick did was to truly count the photons, but it is
   * not clera why that is a good idea.  I have converted them back to mean the number of packets in radiation.c.  It
   * would be simpler if this was an array rather than individual
   * variables
   */
  int ntot_star;
  int ntot_bl;
  int ntot_disk;		/* NSH 15/4/11 Added to count number of photons from the disk in the cell */
  int ntot_wind;
  int ntot_agn;			/* NSH 15/4/11 Added to count number of photons from the AGN in the cell */
  double mean_ds;		/* NSH 6/9/12 Added to allow a check that a thin shell is really optcially thin */
  int n_ds;			/* NSH 6/9/12 Added to allow the mean dsto be computed */
  int nrad;			/* Total number of photons radiated within the cell */
  int nioniz;			/* Total number of photons capable of ionizing H */
  double *ioniz, *recomb;	/* Number of ionizations and recombinations for each ion.
					   The sense is ionization from ion[n], and recombinations 
					   to each ion[n] . 78 - changed to dynamic allocation*/
  int *scatters;		/* 68b - The number of scatters in this cell for each ion. 78 - changed to dynamic allocation*/
  double *xscatters;	/* 68b - Diagnostic measure of energy scattered out of beam on extract. 78 - changed to dynamic allocation */
  double *heat_ion;	/* The amount of energy being transferred to the electron pool
				   sby this ion via photoionization. 78 - changed to dynamic allocation */
  double *lum_ion;	/* The amount of energy being released from the electron pool
				   by this ion via recombination. 78 - changed to dynamic allocation */
  double j, ave_freq, lum;	/*Respectively mean intensity, intensity_averaged frequency, 
				   luminosity and absorbed luminosity of shell */
  double xj[NXBANDS], xave_freq[NXBANDS];	/* 1108 NSH frequency limited versions of j and ave_freq */
  double fmin[NXBANDS];	/* the minimum freqneucy photon seen in a band - this is incremented during photon flight */
  double fmax[NXBANDS];	/* the maximum frequency photon seen in a band - this is incremented during photon flight */
  double fmin_mod[NXBANDS];	/* the minimum freqneucy that the model should be applied for */
  double fmax_mod[NXBANDS];	/* the maximum frequency that the model should be applied for */



  double j_direct,j_scatt; /* 1309 NSH mean intensity due to direct photons and scattered photons */
  double ip_direct,ip_scatt; /* 1309 NSH mean intensity due to direct photons and scattered photons */
  double xsd_freq[NXBANDS];	/*1208 NSH the standard deviation of the frequency in the band */
  int nxtot[NXBANDS];		/* 1108 NSH the total number of photon passages in frequency bands */
  double max_freq;		/*1208 NSH The maximum frequency photon seen in this cell */
  double lum_lines, lum_ff, lum_adiabatic;
  double lum_comp;		/* 1108 NSH The compton luminosity of the cell */
  double lum_di;		/* 1409 NSH The direct ionization luminosity */
  double lum_dr;		/* 1109 NSH The dielectronic recombination luminosity of the cell */
  double lum_fb, lum_z;		/*fb luminosity & fb of metals metals */
  double lum_rad, lum_rad_old;	/* The specfic radiative luminosity in frequencies defined by freqmin
				   and freqmax.  This will depend on the last call to total_emission */


  double lum_ioniz;
  double lum_lines_ioniz, lum_ff_ioniz, lum_adiabatic_ioniz;
  double lum_comp_ioniz;		/* 1108 NSH The compton luminosity of the cell */
  double lum_di_ioniz;		/* 1409 NSH The direct ionization luminosity */
  double lum_dr_ioniz;		/* 1109 NSH The dielectronic recombination luminosity of the cell */
  double lum_fb_ioniz, lum_z_ioniz;		/*fb luminosity & fb of metals metals */
  double lum_rad_ioniz;	/* The specfic radiative luminosity in frequencies defined by freqmin
				   and freqmax.  This will depend on the last call to total_emission */


  double dmo_dt[3];		/*Radiative force of wind */
  int npdf;			/* The number of points actually used in the luminosity pdf */
  int pdf_x[LPDF];		/* The line numbers of *line_ptr which form the boundaries the luminosity pdf */
  double pdf_y[LPDF];		/* Where the pdf is stored -- values between 0 and 1 */
  double gain;			/* The gain being used in interations of the structure */
  double converge_t_r, converge_t_e, converge_hc;	/* Three measures of whether the program believes the grid is converged.
							   The first wo  are the fraction changes in t_r, t_e between this and the last cycle. The third
							   number is the fraction between heating and cooling divided by the sum of the 2       */
  int trcheck, techeck, hccheck;	/* NSH the individual convergence checks used to calculate converge_whole.  Each of these values
					   is 0 if the fractional change or in the case of the last check error is less than a value, currently
					   set to 0.05.  ksl 111126   
NSH 130725 - this number is now also used to say if the cell is over temperature - it is set to 2 in this case   */
  int converge_whole, converging;	/* converge_whole is the sum of the indvidual convergence checks.  It is 0 if all of the
					   convergence checks indicated convergence.subroutine convergence feels point is converged, converging is an
					   indicator of whether the program thought the cell is on the way to convergence 0 implies converging */



  double gamma_inshl[NAUGER];	/*MC estimator that will record the inner shell ionization rate - very similar to macro atom-style estimators */
  /* 1108 Increase sim estimators to cover all of the bands */
  /* 1208 Add parameters for an exponential representation, and a switch to say which we prefer. */
  enum spec_mod_type_enum 
  	{  	
      SPEC_MOD_PL=1, 
  		SPEC_MOD_EXP=2,
      SPEC_MOD_FAIL=-1
	  } spec_mod_type[NXBANDS];	/* NSH 120817 A switch to say which type of representation we are using for this band in this cell. Negative means we have no useful representation, 0 means power law, 1 means exponential */

  double pl_alpha[NXBANDS];	/*Computed spectral index for a power law spectrum representing this cell NSH 120817 - changed name from sim_alpha to PL_alpha */
//  double pl_w[NXBANDS];		/*This is the computed weight of a PL spectrum in this cell - not the same as the dilution factor NSH 120817 - changed name from sim_w to pl_w */
  double pl_log_w[NXBANDS];    /* NSH 131106 - this is the log version of the power law weight. It is in an attempt to allow very large values of alpha to work with the PL spectral model to avoide NAN problems. The pl_w version can be deleted once testing is complete */


  double exp_temp[NXBANDS];	/*NSH 120817 - The effective temperature of an exponential representation of the radiation field in a cell */
  double exp_w[NXBANDS];	/*NSH 120817 - The prefector of an exponential representation of the radiation field in a cell */
  double sim_ip;		/*Ionisation parameter for the cell as defined in Sim etal 2010 */
  double ferland_ip;		/* IP calculaterd from equation 5.4 in hazy1 - assuming allphotons come from 0,0,0 and the wind is transparent */
  double ip;			/*NSH 111004 Ionization parameter calculated as number of photons over the lyman limit entering a cell, divided by the number density of hydrogen for the cell */
  double xi;			/*NSH 151109 Ionization parameter as defined by Taratr et al 1969 and described in Hazy. Its the ionizing flux over the number of hydrogen atoms */
  //int kpkt_rates_known;
  //COOLSTR kpkt_rates;
} plasma_dummy, *PlasmaPtr;

PlasmaPtr plasmamain;

/* A storage area for photons.  The idea is that it is sometimes time-consuming to create the
cumulative distribution function for a process, but trivial to create more than one photon 
of a particular type once one has the cdf,  This appears to be case for f fb photons.  But 
the same procedure could be used for line photons */

#define NSTORE 10
typedef struct photon_store
{
  int n;			/* This is the photon number that was last used */
  double t, f1, f2, freq[NSTORE];

} photon_store_dummy, *PhotStorePtr;

PhotStorePtr photstoremain;

/* 
060616 -- ksl -- 57g -- Modified the plamsa structure to add a new macro stucture that is only needed for macro
atoms. This structure contains most of the variables that were previously dominating the total size of the plasma
array. This effectively solves a problem with producing a huge wind_save file, as well as making the size of the
executable much smaller for the simple atom case.

0608 -- ksl -- 57h -- There is now and extended  discussion in gridwind.c about possible ways to restructure 
the macro structure in order to dynamically set the size of arrays like jbar and jbar_old.  For now, recompilation 
of the code is required.  

0803 -- ksl -- 60 -- The first index is the level, or config  number.  get_atomicdata assures that macro levels, 
if they exist have lower level numbers than other types of levels.  

0911 -- ksl - 68f -- The structure is allocated in a complicated fashion to minimize the total amount of space
taken up by the macro structure, particularly when it is written out to disk.  First the basic array structurre
is allocated (in calloc_macro) and then space for the various arrays contatined in the maccro pointer, like
jbar are allcoated in calloc_esimators.  
*/


typedef struct macro
{
  double *jbar;
  /* This will store the Sobolev mean intensity in transitions which is needed 
     for Macro Atom jumping probabilities. The indexing is by configuration (the 
     NLTE_LEVELS) and then by the upward bound-bound jumps from that level 
     (the NBBJUMPS) (SS) */

  double *jbar_old;

  double *gamma;
  /* This is similar to the jbar but for bound-free transitions. It records the 
     appropriate photoionisation rate co-efficient. (SS) */

  double *gamma_old;

  double *gamma_e;
  /* This is Leon's gamma_e: very similar to gamma but energy weighted. Needed
     for division of photoionisation energy into excitation and k-packets. (SS) */

  double *gamma_e_old;

  double *alpha_st;
  /* Same as gamma but for stimulated recombination rather than photoionisation. (SS)*/

  double *alpha_st_old;

  double *alpha_st_e;
  /* Same as gamma_e but for stimulated recombination rather than photoionisation. (SS) */

  double *alpha_st_e_old;

  double *recomb_sp;
  /* Spontaneous recombination. (SS) */

  double *recomb_sp_e;
  /* "e" version of the spontaneous recombination coefficient. (SS) */

  double *matom_emiss;
  /* This is the specific emissivity due to the de-activation of macro atoms in the cell
     in the frequency range that is required for the final spectral synthesis. (SS) */

  double *matom_abs;  
  /* This is the energy absorbed by the macro atom levels - recorded during the ionization 
     cycles and used to get matom_emiss (SS) */

  /* This portion of the macro structure  is not written out by windsave */
  int kpkt_rates_known;

  double *cooling_bf;
  double *cooling_bf_col;
  double *cooling_bb;

  /* set of cooling rate stores, which are calculated for each macro atom each cycle,
     and used to select destruction rates for kpkts */
  double cooling_normalisation;
  double cooling_bbtot, cooling_bftot, cooling_bf_coltot;
  double cooling_ff;
  double cooling_adiabatic;     // this is just lum_adiabatic / vol / ne

  // double lte_pops[NLEVELS_MACRO];  //a store of LTE level populations for the macro atom.  
  // double lte_pops_norm;

} macro_dummy, *MacroPtr;

MacroPtr macromain;

int xxxpdfwind;			// When 1, line luminosity calculates pdf

int size_Jbar_est, size_gamma_est, size_alpha_est;

#define TMAX_FACTOR			1.5	/*Factor by which t_e can exceed
						   t_r in order for absorbed to 
						   match emitted flux */
#define TMIN				2000.
/* ??? TMIN appears to be used both for the minimum temperature and for 
   calculating the fraction of recombinations that go to the ground state.  This
   looks like a problem ksl-98jul???? */

#define TMAX    5e8/*NSH 130725 - this is the maximum temperature permitted - this was introduced following problems with adaibatically heated cells increasing forever. The value was suggested by DP as a sensible compton teperature for the PK05/P05 Zeus models.*/


//These constants are used in the various routines which compute ionization state
#define SAHA 4.82907e15		/* 2* (2.*PI*MELEC*k)**1.5 / h**3  (Calculated in constants) */
#define MAXITERATIONS	200	//The number of loops to do to try to converge in ne
#define FRACTIONAL_ERROR 0.03	//The change in n_e which causes a break out of the loop for ne
#define THETAMAX	 1e4	//Used in initial calculation of n_e
#define MIN_TEMP	100.	//  ??? this is another minimum temperature - it is used as the minimum tempersture in (JM -- in what??)

// these definitions are for various ionization modes
#define IONMODE_ML93_FIXTE 0                  // Lucy Mazzali using existing t_e (no HC balance)
#define IONMODE_LTE 1                          // LTE using t_r
#define IONMODE_FIXED 2                       // Hardwired concentrations
#define IONMODE_ML93 3                        // Lucy Mazzali
#define IONMODE_PAIRWISE_ML93 6               // pairwise version of Lucy Mazzali 
#define IONMODE_PAIRWISE_SPECTRALMODEL 7      // pariwise modeled J_nu approach
#define IONMODE_MATRIX_BB 8	              // matrix solver BB model
#define IONMODE_MATRIX_SPECTRALMODEL 9        // matrix solver spectral model

// and the corresponding modes in nebular_concentrations
#define NEBULARMODE_TR 0                       // LTE using t_r
#define NEBULARMODE_TE 1                       // LTE using t_e
#define NEBULARMODE_ML93 2                     // ML93 using correction
#define NEBULARMODE_PAIRWISE_ML93 6            // pairwise ML93
#define NEBULARMODE_PAIRWISE_SPECTRALMODEL 7   // pairwise spectral models
#define NEBULARMODE_MATRIX_BB 8	               // matrix solver BB model
#define NEBULARMODE_MATRIX_SPECTRALMODEL 9     // matrix solver spectral model

#define NDIM_MAX 500                // maximum size of the grid in each dimension
double wind_x[NDIM_MAX], wind_z[NDIM_MAX];	/* These define the edges of the cells in the x and z directions */
double wind_midx[NDIM_MAX], wind_midz[NDIM_MAX];	/* These define the midpoints of the cells in the x and z directions */

/* Next two lines are for cyl_var coordinates.  They are used in locating the appropriate 
 * locating the appropriate cell, for example by cylvar_where_in_grid
 */

double wind_z_var[NDIM_MAX][NDIM_MAX];
double wind_midz_var[NDIM_MAX][NDIM_MAX];

typedef struct photon
{
  double x[3];			/* Vector containing position of packet */
  double lmn[3];		/*direction cosines of this packet */
  double freq, freq_orig;    /* current and original frequency of this packet */
  double w,w_orig;		       /* current and original weight of this packet */
  double tau;
  enum istat_enum 
  	{	
    P_INWIND           =0,	//in wind,
		P_SCAT             =1,	//in process of scattering,
		P_ESCAPE           =2,	//Escaped to reach the universe,
		P_HIT_STAR         =3,	//absorbed by photosphere of star,
		P_HIT_DISK         =7,	//Banged into disk
		P_ABSORB           =6,	//Photoabsorbed within wind
		P_TOO_MANY_SCATTERS=4,	//in wind after MAXSCAT scatters
		P_ERROR            =5,	//Too many calls to translate without something happening
		P_SEC              =8,	//Photon hit secondary
		P_ADIABATIC        =9 	//records that a photon created a kpkt which was destroyed by adiabatic cooling
  	} 	istat;					   /*status of photon.*/

  int nscat;			/*number of scatterings */
  int nres;			/*The line number in lin_ptr of last scatter or wind line creation */
  int nnscat;			/* Used for the thermal trapping model of
				   anisotropic scattering to carry the number of
				   scattering to "extract" when needed for wind
				   generated photons SS05. */
  int nrscat;			/* number of resonance scatterings */
  int grid;			/*grid position of the photon in the wind, if
				   the photon is in the wind.  If the photon is not
				   in the wind, then -1 implies inside the wind cone and  
				   -2 implies outside the wind */

  enum origin_enum
  	{	PTYPE_STAR=0, 
  		PTYPE_BL=1, 
  		PTYPE_DISK=2,
  		PTYPE_WIND=3,
  		PTYPE_AGN=4,
      PTYPE_STAR_MATOM=10,     
      PTYPE_BL_MATOM=11, 
      PTYPE_DISK_MATOM=12,
      PTYPE_WIND_MATOM=13,
      PTYPE_AGN_MATOM=14
  	} 	origin;				/* Where this photon originated.  If the photon has
		   					         scattered it's "origin" may be changed to "wind".*/
                      /* note that we add 10 to origin when processed by a macro-atom
                         which means we need these values in the enum list */
  int np;			/*NSH 13/4/11 - an internal pointer to the photon number so 
				   so we can write out details of where the photon goes */
  double path; /* SWM - Photon path length */

}
p_dummy, *PhotPtr;

/* minimum value for tau for p_escape_from_tau function- below this we 
   set to p_escape_ to 1 */
#define TAU_MIN 1e-6


/* 68b - ksl - This is a structure in which the history of a single photon bundle can be recorded
 * See phot_util   phot_hist().  It needs to be used carefully.  if phot_hist_on is true
 * then photon histories will be attempted.
 */

#define MAX_PHOT_HIST	1000
int n_phot_hist, phot_hist_on, phot_history_spectrum;
struct photon xphot_hist[MAX_PHOT_HIST];

struct basis
{
  double a[3][3];

}
basis_cartesian;


/* The next section defines the spectrum arrays.  The spectrum structure contains
   the selection criteria for each spectrum as well as the array in which to store the
   spectrum.  The first MSPEC spectra are reserved for the total generated spectrum,
   total emitted spectrum, the total spectrum of photons which are scattered and 
   the total spectrum of photons which are absorbed.  The remainder of the spectra pertain 
   to the spectrum as seen from various directions. Note that the space for the spectra 
   are allocated using a calloc statement in spectrum_init.  1409-ksl-A new spectrum
   was added.  This will be the first of the spectra.  It is simply the generated spectrum
   before passing through the wind.  It has the orginal weights as generated.  */



#define MSPEC                            7
int nspectra;			/* After create_spectrum, the number of elements allocated for s, or 
				   alternatively the number of spectra one has to work with.  Note that
				   general s[0],s[1] and s[2] are the escaping, scattered and absorbed photons,
				   while elements higher than this will contain spectra as seen by different observers */


int nscat[MAXSCAT + 1], nres[MAXSCAT + 1], nstat[NSTAT];

typedef struct spectrum
{
  char name[40];
  float freqmin, freqmax, dfreq;
  float lfreqmin, lfreqmax, ldfreq;	/* NSH 1302 - values for logarithmic spectra */
  double lmn[3];
  double mmax, mmin;		/* Used only in live or die situations, mmax=cos(angle-DANG_LIVE_OR_DIE)
				   and mmim=cos(angle+DANG_LIVE_OR_DIE).   In actually defining this
				   one has to worry about signs and exactly whether you are above
				   or below the plane */
  double renorm;		/* Factor used in Live or die option where it is 4*PI/domega;
				   1.0 otherwise */
  int nphot[NSTAT];		/* Used to help define what happened to photons when trying to construct this
				   particular spectrum */
  int nscat;			/* nscat>999 -> select everything
				   0 < nscat < MAXScat select only those photons which have
				   scattered nscat times, number of scattering photons in spectrum, if nscat is negative
				   then all photons with more than |nscat| are included.  */
  int top_bot;			/* 0 ->select photons regardless of location 
				   >0     -> select only photons whose "last" position is above the disk
				   <0    -> select only photons whose last position is below the disk */
  double x[3], r;		/* The position and radius of a special region from which to extract spectra  */
  double f[NWAVE];
  double lf[NWAVE];		/* a second array to hole the extracted spectrum in log units */
  double lfreq[NWAVE];		/* We need to hold what freqeuncy intervals our logarithmic spectrum has been taken over */
}
spectrum_dummy, *SpecPtr;

/*1409 - ksl - Replaced variable s with xxspec to avoid confusion in case some wanted to use a variable s */
SpecPtr xxspec;




/* Parameters used only by py_wind 
 * py_wind_projecti	0 -> simply print the various parameters without 
 * 			atempting toproject onto a yz plane
 * 			1 -> project onto a yz plane.
 */

int py_wind_min, py_wind_max, py_wind_delta, py_wind_project;
double *aaa;			// A pointer to an array used by py_wind

/* This is the structure needed for a cumulative distribution function. The CDFs are
generated from a function which is usually only proportional to the probability density
function.  It is sometimes useful, e.g. in calculating the reweighting function to
have access to the proper normalization.  Since the one needs the normalization to
properly create the CDF, this was added for python_43.2  */
#define NPDF 200

typedef struct Pdf
{
  double x[NPDF + 1];		/* Positions for which the probability density
				   is calculated */
  double y[NPDF + 1];		/* The value of the CDF at x */
  double d[NPDF + 1];		/* 57i -- the rate of change of the probability
				   density at x */
  double limit1, limit2;	/* Limits (running from 0 to 1) that define a portion
				   of the CDF to sample */
  double x1, x2;		/* limits if they exist on what is returned */
  double norm;			//The scaling factor which would renormalize the pdf
}
 *PdfPtr, pdf_dummy;


/* Variable used to allow something to be printed out the first few times
   an even occurs */
int itest, jtest;

char hubeny_list[132];		//Location of listing of files representing hubeny atmospheres




// Allow for a diagnostic file 

FILE *epltptr;			//TEST
// diag_on_off is deprecated see #111, #120
//int diag_on_off;		// on is non-zero  //TEST


/* These variables are stored or used by the routines for anisotropic scattering */
/* Allow for the transfer of tau info to scattering routine */

/* JM 1411 -- tau_x_dvds doesn't appear to be used anywhere, so I've 
   made it a local variable rather than global */  
//double tau_x_dvds;		//tau_x_dvds/dvds is the actual tau
//double tau_scatter_min;               //Set in subroutine scatter for use by extract

struct Pdf pdf_randwind_store[100];
PdfPtr pdf_randwind;
struct photon phot_randwind;

/* N.B. pdf_randwind and phot_randwind are used in the routine anisowind for 
as part of effort to incorporate anisotropic scattering in to python.  
Added for python_43.2 */


/* Provide generally for having arrays which descibe the 3 xyz axes. 
these are initialized in main  */

double x_axis[3];
double y_axis[3];
double z_axis[3];

/* These are structures associated with frequency limits/s used for photon
generation and for calculating heating and cooling */

#define NBANDS 20
struct xbands
{
  double f1[NBANDS], f2[NBANDS];
  double alpha[NBANDS];
  double pl_const[NBANDS];
  double min_fraction[NBANDS];
  double nat_fraction[NBANDS];	// The fraction of the accepted luminosity in this band
  double used_fraction[NBANDS];
  double flux[NBANDS];		//The "luminosity" within a band
  double weight[NBANDS];
  int nphot[NBANDS];
  int nbands;			// Actual number of bands in use
}
xband;


/* The next section contains the freebound structures that can be used for both the
 * specific emissivity of a free-bound transition, and for the recombination coefficient
 * assuming the array has been initialized, which can take a few minutes
*/

#define NTEMPS	60		// The number of temperatures which are stored in each fbstruct
				/* NSH this was increased from 30 to 60 to take account of 3 extra OOM 
				intemperature we wanted to have in fb */ 
#define NFB	10		// The maximum number of frequency intervals for which the fb emission is calculated

struct fbstruc
{
  double f1, f2;
  double emiss[NIONS][NTEMPS];
}
freebound[NFB];

double xnrecomb[NIONS][NTEMPS];	// There is only one set of recombination coefficients
double fb_t[NTEMPS];
int nfb;			// Actual number of freqency intervals calculated


//This is a new structure to contain the frequency range of the final spectrum
//During the ionization cycles, the emissivity due to k-packets and macro atom
//deactivations in this range will be computed and then used in the final spectral
//synthesis part of the code (SS June04).
struct emiss_range
{
  double fmin, fmax;		// min and max frequency required in the final spectrum
}
em_rnge;


#include "version.h"		/*54f -- Added so that version can be read directly */
#include "templates.h"
#include "recipes.h"

// 04apr ksl -- made kap_bf external so can be passed around variables
double kap_bf[NLEVELS];



// 12jun nsh - some commands to enable photon logging in given cells. There is also a pointer in the geo

FILE *pstatptr;			//NSH 120601 - pointer to a diagnostic file that will contain photon data for given cells
int cell_phot_stats;		//1=do  it, 0=dont do it
#define  NCSTAT 10		//The maximum number of cells we are going to log
int ncstat;			// the actual number we are going to log
int ncell_stats[NCSTAT];	//the numbers of the cells we are going to log


/* Added variables which count number of times two situations occur (See #91) */
int nerr_no_Jmodel;
int nerr_Jmodel_wrong_freq;



// advanced mode variables
struct advanced_modes
{
  /* these are all 0=off, 1=yes */
  int iadvanced;                // this is controlled by the -d flag, global mode control.
  int save_cell_stats;          // want to save photons statistics by cell
  int ispy;                     // want to use the ispy function
  int keep_ioncycle_windsaves;  // want to save wind file each ionization cycle
  int track_resonant_scatters;  // want to track resonant scatters
  int save_extract_photons;     // we want to save details on extracted photons
  int print_windrad_summary;    // we want to print the wind rad summary each cycle
  int adjust_grid;              // the user wants to adjust the grid scale
  int diag_on_off;              // extra diagnostics
  int use_debug;                // print out debug statements
  int print_dvds_info;          // print out information on the velocity gradients
  int keep_photoabs;            // keep photoabsorption in final spectrum
  int quit_after_inputs;        // quit after inputs read in, testing mode
  int fixed_temp;               // do not alter temperature from that set in the parameter file
  int zeus_connect;				// We are connecting to zeus, do not seek new temp and output a heating and cooling file
  int rand_seed_usetime;        // default random number seed is fixed, not based on time
}
modes;


FILE *optr;  //pointer to a diagnostic file that will contain dvds information



struct filenames
{
  char root[LINELENGTH];        // main rootname
  char windsave[LINELENGTH];    // wind save filename
  char specsave[LINELENGTH];    // spec save filename
  char diag[LINELENGTH];        // diag file
  char diagfolder[LINELENGTH];  // diag folder
  char old_windsave[LINELENGTH];// old windsave name
  char input[LINELENGTH];       // input name if creating new pf file
  char lspec[LINELENGTH];       // log_spec_tot fname
  char wspec[LINELENGTH];       // spectot fname
  char disk[LINELENGTH];        // disk diag fname
  char tprofile[LINELENGTH];    // non standard tprofile fname
  char phot[LINELENGTH];        // photfile e.g. python.phot
  char windrad[LINELENGTH];     // wind rad file
  char spec[LINELENGTH];        // .spec file
}
files;


        
#define NMAX_OPTIONS 20 

/* these two variables are used by xdefine_phot() in photon_gen.c 
   to set the mode for get_matom_f()in matom.c and tell it 
   whether it has already calculated the matom emissivities or not. */
#define CALCULATE_MATOM_EMISSIVITIES 0
#define USE_STORED_MATOM_EMISSIVITIES 1