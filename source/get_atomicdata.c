
/***********************************************************/
/** @file  get_atomicdata.c
 * @author ksl
 * @date   March, 2018
 *
 * @brief  Read in all of the atomic data for use with Python
 * and other similar programs
 *
 ***********************************************************/


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "atomic.h"
#include "log.h"
// If routines are added cproto > atomic_proto.h should be run
#include "atomic_proto.h"




#define LINELENGTH 400
#define MAXWORDS    20


/**********************************************************/
/**
 * @brief      generalized subroutine for reading atomic data
 *  	into a set of structures defined in "atomic.h"
 *
 *
 * @param [in] char  masterfile[]   The name of the "masterfile" which refers to other files which contain the data
 * @return     Always returns 0
 *
 * @details
 *
 * get_atomic_data reads in all the atomic data.  It also converts the data to cgs units unless
 * 	otherwise noted, e.g ionization potentials are converted to ergs.
 *
 *
 *
 * The masterfile is a list of other files, which contain atomic data for specific puruposes.
 *
 * All of the files are ascii.  The information is keyword based.
 *
 * The order of the data is important.  Elements should be defined before ions; ions
 * before levels, levels before  lines etc.  In most if not all cases, one can either define all
 * of the elements first, all of the ions, etc ... or the first element and its ions, the
 * second element and its ions etc. Commenting out an element in the datafile has the effect of
 * eliminating all of the data associated with that element, e.g ions, lines, etc.
 *
 * The program assumes that the ions for a given element are grouped together,
 * that the levels for a given element are grouped together, etc.
 *
 * If one wants to read both Topbase and VFKY photoionization x-sections, then the topbase
 * x-sections should be read first.  The routine will ignore data in the VFKY list if there
 * is a pre-existing Topbase data.  The routine will stop if you but a VFKY x-section first
 * on the assumption that Topbase data should trump VFKY data. (Making the program more flexible
 * would be difficult because you would have to drop various bits of Topbase data out of the
 * middle of the structure or mark it NG or something.)
 *
 * Similarly, as a general rule the macro information should be provided before the simple
 * ions are described, and the macro lines should be presented before simple lines.
 *
 * The only data which can be mixed up completely is the line array.
 *
 * Finally, get_atomic_data creates a pointer array to  lines, which has the lines in
 * frequency ascending order.   This is actually done by a small subroutine index_lines
 * which in turn calls a Numerical Recipes routine.
 *
 * ### Notes ###
 *
 * get_atomic data is intended to be stand-alone, that is one should be able to use it for routines
 * other than Python, e.g for another routine intended to calculate the ionization state of
 * a plasma in collisional equilibrium.
 *
 * To this end, the routines populate stuctures in atomic.h, which are not part of python.h, and 
 * one should avoid calling routines like Exit(0) that are very python centric.  It's important
 * that future modifications to get_atomic_data maintain this independence.
 *
 *
 *
 **********************************************************/

int
get_atomic_data (masterfile)
     char masterfile[];
{
  int match;
  FILE *fptr, *mptr;            //, fopen ();
  char aline[LINELENGTH];       //, *fgets ();
  char file[LINELENGTH], atomic_file[LINELENGTH];

  char word[LINELENGTH];
  int n, m, i, j;
  int n1, n2;                   //081115 nsh two new counters for DR - use new pointers to avoid any clashes!
  int nparam;                   //081115 nsh temperary holder for number of DR parameters
  double drp[MAX_DR_PARAMS];    //081115 nsh array to hold DR parameters prior to putting into structure
  double btrr[T_RR_PARAMS];     //0712 nsh array to hole badnell total RR params before putting into structure
  int ne, w;                    //081115 nsh new variables for DR variables
  int mstart, mstop;
  int nelem;
  double gl, gu;
  double el, eu;
  int qnum;
  double qqnum, ggg, gg;
  int istate, z, nion;
  int iistate, zz;
  int levl, levu;
  int in, il;                   //The levels used in inner shell data
  double q;
  double freq, f, exx, et, p;
  double the_ground_frac[20];
  char choice;
  int lineno;                   /* the line number in the file beginning with 1 */
  int simple_line_ignore[NIONS], cstren_no_line;
  int nwords;
  int nlte, nmax;
  int mflag;                    //flag to identify reading data for macro atoms
  int nconfigl, nconfigu;       //internal labels for configurations
  int islp, ilv, np;
  char configname[15];
  double e, rl;
  double xe[NCROSS], xx[NCROSS];
  // double a21 ();
  int nlines_simple;
  int nspline;
  double tmin;
  int nions_simple, nions_macro;
  int nlevels_simple;
  int ntop_phot_simple, ntop_phot_macro;
  int bb_max, bf_max;
  int lev_type;
  int nn;
//  double yield;
  double gstemp[BAD_GS_RR_PARAMS];      //Temporary storage for badnell resolved GS RR rates
  double temp[LINELENGTH];      //Temporary storage for data read in off a line this is enogh if every character on the
  char gsflag, drflag;          //Flags to say what part of data is being read in for DR and RR
  double gstmin, gstmax;        //The range of temperatures for which all ions have GS RR rates
  double gsqrdtemp, gfftemp, s1temp, s2temp, s3temp;    //Temporary storage for gaunt factors
  int n_elec_yield_tot;         //The number of inner shell cross sections with matching electron yield arrays
  int inner_no_e_yield;         //The number of inner shell cross sections with no yields
//  int n_fluor_yield_tot;        //The number of inner shell cross sections with matching fluorescent photon yield arrays
  double I, Ea;                 //The ionization energy and mean electron energy for electron yields
//  double energy;                //The energy of inner shell fluorescent photons
  int c_l, c_u;                 //Chianti level indicators
  double en, gf, hlt, sp;       //Parameters in collision strangth file
  int type;                     //used in collision strength

  /* define which files to read as data files */


/* Allocate structures for storage of data */

  if (ele != NULL)
  {
    free (ele);
  }
  ele = (ElemPtr) calloc (sizeof (ele_dummy), NELEMENTS);

  if (ele == NULL)
  {
    Error ("There is a problem in allocating memory for the element structure\n");
    exit (0);
  }
  else
  {
    Log_silent
      ("Allocated %10d bytes for each of %6d elements of   elements totaling %10.0f Mb \n",
       sizeof (ele_dummy), NELEMENTS, 1.e-6 * NELEMENTS * sizeof (ele_dummy));
  }


  if (ion != NULL)
  {
    free (ion);
  }
  ion = (IonPtr) calloc (sizeof (ion_dummy), NIONS);

  if (ion == NULL)
  {
    Error ("There is a problem in allocating memory for the ion structure\n");
    exit (0);
  }
  else
  {
    Log_silent
      ("Allocated %10d bytes for each of %6d elements of       ions totaling %10.1f Mb \n",
       sizeof (ion_dummy), NIONS, 1.e-6 * NIONS * sizeof (ion_dummy));
  }



  if (config != NULL)
  {
    free (config);
  }
  config = (ConfigPtr) calloc (sizeof (config_dummy), NLEVELS);

  if (config == NULL)
  {
    Error ("There is a problem in allocating memory for the config structure\n");
    exit (0);
  }
  else
  {
    Log_silent
      ("Allocated %10d bytes for each of %6d elements of     config totaling %10.1f Mb \n",
       sizeof (config_dummy), NLEVELS, 1.e-6 * NLEVELS * sizeof (config_dummy));
  }




  if (line != NULL)
  {
    free (line);
  }
  line = (LinePtr) calloc (sizeof (line_dummy), NLINES);

  if (line == NULL)
  {
    Error ("There is a problem in allocating memory for the line structure\n");
    exit (0);
  }
  else
  {
    Log_silent
      ("Allocated %10d bytes for each of %6d elements of       line totaling %10.1f Mb \n",
       sizeof (line_dummy), NLINES, 1.e-6 * NLINES * sizeof (line_dummy));
  }



  /* Initialize variables */



  phot_freq_min = VERY_BIG;
  inner_freq_min = VERY_BIG;

  for (n = 0; n < NELEMENTS; n++)
  {
    strcpy (ele[n].name, "none");
    ele[n].z = (-1);
    ele[n].abun = (-1);
    ele[n].firstion = (-1);
    ele[n].nions = (-1);
    ele[n].istate_max = (-1);
  }

  for (n = 0; n < NIONS; n++)
  {
    simple_line_ignore[n] = 0;  // diagnostic counter for how many lines ignored
    ion[n].z = (-1);
    ion[n].istate = (-1);
    ion[n].nelem = (-1);
    ion[n].ip = (-1);
    ion[n].g = (-1);
    ion[n].nmax = (-1);
    ion[n].firstlevel = (-1);
    ion[n].nlevels = (-1);
    ion[n].first_nlte_level = (-1);
    ion[n].first_levden = (-1);
    ion[n].nlte = (-1);
    ion[n].phot_info = (-1);
    ion[n].macro_info = (-1);   //Initialise - don't know if using Macro Atoms or not: set to -1 (SS)
    ion[n].ntop_first = 0;      // The fact that ntop_first and ntop  are initialized to 0 and not -1 is important
    ion[n].ntop_ground = 0;     //NSH 0312 initialize the new pointer for GS cross sections
    ion[n].ntop = 0;
    ion[n].nxphot = (-1);
    ion[n].lev_type = (-1);     // Initialise to indicate we don't know what types of configurations will be read
    ion[n].drflag = 0;          //Initialise to indicate as far as we know, there are no dielectronic recombination parameters associated with this ion.
    ion[n].total_rrflag = 0;    //Initialise to say this ion has no Badnell total recombination data
    ion[n].nxtotalrr = -1;      //Initialise the pointer into the bad_t_rr structure.
    ion[n].bad_gs_rr_t_flag = 0;        //Initialise to say this ion has no Badnell ground state recombination data
    ion[n].bad_gs_rr_r_flag = 0;        //Initialise to say this ion has no Badnell ground state recombination data
    ion[n].nxbadgsrr = -1;      //Initialise the pointer into the bad_gs_rr structure.
    ion[n].dere_di_flag = 0;    //Initialise to say this ion has no Dere DI rate data
    ion[n].nxderedi = -1;       //Initialise the pointer into the Dere DI rate structure
    ion[n].n_inner = 0;         //Initialise the pointer to say we have no inner shell ionization cross sections
    for (i = 0; i < N_INNER; i++)
      ion[n].nxinner[i] = -1;   //Inintialise the inner shell pointer array
  }

  nlevels = nxphot = nphot_total = ntop_phot = nauger = ndrecomb = n_inner_tot = 0;     //Added counter for DR//
  n_elec_yield_tot = 0;         //Counter for electron yield
  //  n_fluor_yield_tot = 0;     and fluorescent photon yields

  /*This initializes the top_phot array - it is used for all ionization processes so some elements
     are only used in some circumstances
   */

  for (n = 0; n < NLEVELS; n++)
  {
    phot_top[n].nlev = (-1);
    phot_top[n].uplev = (-1);
    phot_top[n].nion = (-1);    //the ion to which this cross section belongs
    phot_top[n].n_elec_yield = -1;      //pointer to the electron yield array (for inner shell)
//    phot_top[n].n_fluor_yield = -1;     //pointer to the fluorescent photon yield (for inner shell)
    phot_top[n].n = -1;         //pointer to shell (inner shell)
    phot_top[n].l = -1;         //pointer to l subshell (inner shell only)
    phot_top[n].z = (-1);       //atomic number
    phot_top[n].np = (-1);      //number of points in the fit
    phot_top[n].macro_info = (-1);      //Initialise - don't know if using Macro Atoms or not: set to -1 (SS)
    for (j = 0; j < NCROSS; j++)        //initialise the crooss sectiond
    {
      phot_top[n].freq[j] = (-1);
      phot_top[n].x[j] = (-1);
    }
    phot_top[n].f = (-1);       //last frequency
    phot_top[n].sigma = 0.0;    //last cross section
  }


  for (n = 0; n < NIONS * N_INNER; n++) //Initialise atomic arrasy with dimension NIONS*NINNER
  {
    inner_cross[n].nlev = (-1);
    inner_cross[n].uplev = (-1);
    inner_cross[n].nion = inner_elec_yield[n].nion = inner_fluor_yield[n].nion = (-1);
    inner_cross[n].n_elec_yield = -1;
//    inner_cross[n].n_fluor_yield = -1;
    inner_cross[n].n = inner_elec_yield[n].n = inner_fluor_yield[n].n = (-1);
    inner_cross[n].l = inner_elec_yield[n].l = inner_fluor_yield[n].l = (-1);
    inner_cross[n].z = inner_elec_yield[n].z = inner_fluor_yield[n].z = (-1);
    inner_elec_yield[n].I = inner_elec_yield[n].Ea = 0.0;
    inner_fluor_yield[n].freq = inner_fluor_yield[n].yield = 0.0;
    for (j = 0; j < 10; j++)
      inner_elec_yield[n].prob[j] = 0.0;
    inner_cross[n].np = (-1);
    inner_cross[n].macro_info = (-1);   //Initialise - don't know if using Macro Atoms or not: set to -1 (SS)
    for (j = 0; j < NCROSS; j++)
    {
      inner_cross[n].freq[j] = (-1);
      inner_cross[n].x[j] = (-1);
    }
    inner_cross[n].f = (-1);
    inner_cross[n].sigma = 0.0;
  }




  for (i = 0; i < NLEVELS; i++)
  {
    config[i].n_bbu_jump = 0;   // initialising the number of jumps from each level to 0. (SS)
    config[i].n_bbd_jump = 0;
    config[i].n_bfu_jump = 0;
    config[i].n_bfd_jump = 0;
  }

  for (n = 0; n < NLINES; n++)
  {
    line[n].freq = -1;
    line[n].f = 0;
    line[n].nion = -1;
    line[n].gl = line[n].gu = 0;
    line[n].el = line[n].eu = 0.0;
    line[n].macro_info = -1;
    line[n].coll_index = -999;
  }

/* The following lines initialise the dielectronic recombination structure */
  for (n = 0; n < NIONS; n++)
  {
    drecomb[n].nion = -1;
    drecomb[n].nparam = -1;     //the number of parameters - it varies from ion to ion
    for (n1 = 0; n1 < MAX_DR_PARAMS; n1++)
    {
      drecomb[n].c[n1] = 0.0;
      drecomb[n].e[n1] = 0.0;
    }
  }


/* The following lines initialise the badnell total radiative recombination rate structure */
  for (n = 0; n < NIONS; n++)
  {
    total_rr[n].nion = -1;
    for (n1 = 0; n1 < T_RR_PARAMS; n1++)
    {
      total_rr[n].params[n1] = 0.0;
    }

  }


/* The following lines initialise the badnell total radiative recombination rate structure */
  for (n = 0; n < NIONS; n++)
  {
    bad_gs_rr[n].nion = -1;
    for (n1 = 0; n1 < BAD_GS_RR_PARAMS; n1++)
    {
      bad_gs_rr[n].temps[n1] = 0.0;
      bad_gs_rr[n].rates[n1] = 0.0;
    }

  }
  gstmin = 0.0;
  gstmax = 1e99;


/* The following lines initialise the dere direct ionization rate struture */
  for (n = 0; n < NIONS; n++)
  {
    dere_di_rate[n].nion = -1;
    dere_di_rate[n].xi = 0.0;
    dere_di_rate[n].min_temp = 1e99;
    for (n1 = 0; n1 < DERE_DI_PARAMS; n1++)
    {
      dere_di_rate[n].temps[n1] = 0.0;
      dere_di_rate[n].rates[n1] = 0.0;
    }

  }




/* The following lines initialise the Sutherland gaunt factors */
  gaunt_n_gsqrd = 0;            //The number of sets of scaled temperatures we have data for
  for (n = 0; n < MAX_GAUNT_N_GSQRD; n++)
  {
    gaunt_total[n].log_gsqrd = 0.0;
    gaunt_total[n].gff = 0.0;
    gaunt_total[n].s1 = 0.0;
    gaunt_total[n].s2 = 0.0;
    gaunt_total[n].s3 = 0.0;
  }


/* The following lines initialise the collision strengths */
  n_coll_stren = 0;             //The number of data sets
  cstren_no_line = 0;           // counter to track how many times we don't find a matching line 
  for (n = 0; n < NLINES; n++)
  {
    coll_stren[n].n = -1;       //Internal index
    coll_stren[n].lower = -1;   //The lower energy level - this is in Chianti notation and is currently unused
    coll_stren[n].upper = -1;   //The upper energy level - this is in Chianti notation and is currently unused
    coll_stren[n].energy = 0.0; //The energy of the transition
    coll_stren[n].gf = 0.0;
    coll_stren[n].hi_t_lim = 0.0;       //The high temerature limit
    coll_stren[n].n_points = 0.0;       //The number of points in the splie fit
    coll_stren[n].type = -1;    //The type of fit, this defines how one computes the scaled temperature and scaled coll strength
    coll_stren[n].scaling_param = 0.0;  //The scaling parameter C used in the Burgess and Tully calculations
    for (n1 = 0; n1 < N_COLL_STREN_PTS; n1++)
    {
      coll_stren[n].sct[n1] = 0.0;      //The scaled temperature points in the fit
      coll_stren[n].scups[n1] = 0.0;
    }
  }



  choice = 'x';                 /* Start by assuming you cannot determine what kind of line it is */
  nelements = 0;
  nions = nions_simple = nions_macro = 0;
  nlevels = nlevels_simple = nlevels_macro = 0;
  ntop_phot_simple = ntop_phot_macro = 0;
  nlte_levels = 0;
  nlines = nlines_simple = nlines_macro = 0;
  lineno = 0;
  nxphot = 0;
/*mflag is set initially to 1, in order to establish that
macro-lines need to be read in before any "simple" lines.  This
is so that we can assure that the first lines in the line array
are macro-lines.   It is important to recognize that the lin_ptr
structure does not have this property! */
  mflag = 1;


/* Completed all initialization */

  /* OK now we can try to read in the data from the data files */

  if ((mptr = fopen (masterfile, "r")) == NULL)
  {
    Error ("Get_atomic_data: Could not find masterfile %s in current directory\n", masterfile);
    exit (1);
  }

  Log ("Get_atomic_data: Reading from masterfile %s\n", masterfile);

/* Open and read each line in the masterfile in turn */

  while (fgets (aline, LINELENGTH, mptr) != NULL)
  {
    if (sscanf (aline, "%s", file) == 1 && file[0] != '#')
    {

      /*
       * Open one of the files designated in the masterfile and begin to read it
       */

      if ((fptr = fopen (file, "r")) == NULL)
      {
        Error ("Get_atomic_data: Could not open %s in current directory\n", file);
        exit (1);
      }

      Log_silent ("Get_atomic_data: Reading data from %s\n", atomic_file);
      lineno = 1;

      /* Main loop for reading each data file line by line */

      while (fgets (aline, LINELENGTH, fptr) != NULL)
      {
        lineno++;

        strcpy (word, "");      /*For reasons which are not clear, word needs to be reinitialized every time to
                                   properly deal with blank lines */

        if (sscanf (aline, "%s", word) == 0 || strlen (word) == 0)
          choice = 'c';         /*It's a a blank line, treated like a comment */
        else if (strncmp (word, "!", 1) == 0)
          choice = 'c';
        else if (strncmp (word, "#", 1) == 0)
          choice = 'c';         /* It's a comment */
        else if (strncmp (word, "CSTREN", 6) == 0)      //collision strengths
          choice = 'C';
        else if (strncmp (word, "Element", 5) == 0)
          choice = 'e';
        else if (strncmp (word, "Ion", 3) == 0)
          choice = 'i';
        else if (strncmp (word, "LevTop", 6) == 0)
          choice = 'N';
        else if (strncmp (word, "LevMacro", 8) == 0)    // This indicated leves for a Macro Atom (SS)
          choice = 'N';
        else if (strncmp (word, "Level", 3) == 0)       // There are various records of this type
          choice = 'n';
        else if (strncmp (word, "Phot", 4) == 0)        // There are various records of this type
          choice = 'w';         // Macro Atom Phots are a subset of these (SS)
        else if (strncmp (word, "Line", 4) == 0)
          choice = 'r';
        else if (strncmp (word, "LinMacro", 8) == 0)    //This indicates lines for a Macro Atom (SS)
          choice = 'r';
        else if (strncmp (word, "Frac", 4) == 0)
          choice = 'f';         /*ground state fractions */
        else if (strncmp (word, "InnerVYS", 8) == 0)
          choice = 'I';         /*Its a set of inner shell photoionization cross sections */
        else if (strncmp (word, "DR_BADNL", 8) == 0)    /* It's a badnell type dielectronic recombination file */
          choice = 'D';
        else if (strncmp (word, "DR_SHULL", 8) == 0)    /*its a schull type dielectronic recombination */
          choice = 'S';
        else if (strncmp (word, "RR_BADNL", 8) == 0)    /*Its a badnell type line in the total RR file */
          choice = 'T';
        else if (strncmp (word, "DI_DERE", 7) == 0)     /*Its a data file giving direct ionization rates from Dere (2007) */
          choice = 'd';
        else if (strncmp (word, "RR_SHULL", 8) == 0)    /*Its a shull type line in the total RR file */
          choice = 's';
        else if (strncmp (word, "BAD_GS_RR", 9) == 0)   /*Its a badnell resolved ground state RR file */
          choice = 'G';
        else if (strncmp (word, "FF_GAUNT", 8) == 0)    /*Its a data file giving the temperature averaged gaunt factors from Sutherland (1998) */
          choice = 'g';
        else if (strncmp (word, "Kelecyield", 10) == 0) /*Electron yield from inner shell ionization fro Kaastra and Mewe */
          choice = 'K';
//        else if (strncmp (word, "Kphotyield", 10) == 0) /*Floruescent photon yield from IS ionization from Kaastra and Mewe */
//          choice = 'F';
        else if (strncmp (word, "*", 1) == 0);  /* It's a continuation so record type remains same */

        else
          choice = 'z';         /* Who knows what it is */


        switch (choice)
        {
/**
 * @section Elements
 *
 * A typical element has the following format
 *
 *  Element    6    C    8.56
 *
 * where 6 here refers to z of the elemnt, C is the name, and 8.56 is the abundance relative to H at 12
 *
 * */
        case 'e':
          if (sscanf (aline, "%*s %d %s %le", &ele[nelements].z, ele[nelements].name, &ele[nelements].abun) != 3)
          {
            Error ("Get_atomic_data: file %s line %d: Element line incorrectly formatted\n", file, lineno);
            Error ("Get_atomic_data: %s\n", aline);
            exit (0);
          }
          ele[nelements].abun = pow (10., ele[nelements].abun - 12.0);  /* Immediate replace by number density relative to H */
          nelements++;
          if (nelements > NELEMENTS)
          {
            Error ("getatomic_data: file %s line %d: More elements than allowed. Increase NELEMENTS in atomic.h\n", file, lineno);
            Error ("Get_atomic_data: %s\n", aline);
            exit (0);
          }
          break;


/**
 * @section Ions
 *
 * Ions are specified with the keyword Ion.  A typical entry for ions looks like
 *
 *  @verbatim
 *  IonV    C   6   4   2   64.49400 1000   10     1s^22s(2S_{1/2})
 *  @endverbatim
 *
 * where C is the elment name, 6 is the elemnt z, 4 is the ionization state (in conventional
 * astrophysics notation, 2 is the mulitplicity  of the ground state, and 1000 and 10 refer
 * to maximum number of allowed LTE and nLTE configurations.  The configuration is not actually
 * required and is for information only.
 *
 * @detail
 *
 * Note that keyword parsing is done only on the first 3 letters, so Ion and IonV are equivlalent
 *
 * The routine still supports an earlier format which consists of the first 6 "words" in which case
 * default values are assigned to for the maximum number of LTE/simple and nLTE/macro atom configurations.
 *
 * @bug Exactly what is meant by LTE here needs clarification.  Is there code here that should
 * be removed because we do not need the "earlier" format?
 *
 */

        case 'i':

          if ((nwords = sscanf (aline, "%*s %*s %d %d %le %le %d %d", &z, &istate, &gg, &p, &nmax, &nlte)) != 6)
          {
            Error ("get_atomic_data: file %s line %d: Ion istate line incorrectly formatted\n", file, lineno);
            Error ("Get_atomic_data: %s\n", aline);
            exit (0);
          }
// Now check that an element line for this ion has already been read
          n = 0;
          while (ele[n].z != z && n < nelements)
            n++;
          if (n == nelements)
          {

            Debug ("get_atomic_data: file %s line %d has ion for unknown element with z %d\n", file, lineno, z);
            break;
          }

// Now populate the ion structure

          if (nlte > 0)
          {                     // Then we want to consider some of these levels as non-lte
            ion[nions].first_levden = nlte_levels;      /* This is the index to into
                                                           the levden aray */
            ion[nions].n_lte_max = nlte;        //Reserve this many elements of levden
            nlte_levels += nlte;
            if (nlte_levels > NLTE_LEVELS)
            {
              Error ("get_atomic_data: nlte_levels (%d) > NLTE_LEVELS (%d)\n", nlte_levels, NLTE_LEVELS);
              exit (0);
            }

          }
          ion[nions].z = z;
          ion[nions].istate = istate;
          ion[nions].g = gg;
          ion[nions].ip = p * EV2ERGS;
          ion[nions].nmax = nmax;
/* Use the keyword IonM to classify the ion as a macro-ion (IonM) or not (simply Ion) */
          if (strncmp (word, "IonM", 4) == 0)
          {
            ion[nions].macro_info = 1;
            nions_macro++;
          }
          else
          {
            ion[nions].macro_info = 0;
            nions_simple++;
          }
          nions++;
          if (nions == NIONS)
          {
            Error
              ("getatomic_data: file %s line %d: %d ions is more than %d allowed. Increase NIONS in atomic.h\n",
               file, lineno, nions, NIONS);
            exit (0);
          }
          break;

/**
 * @section levels Levels or Configurations
 *
 * This version of get_atomicdata can read various types of configuration inputs, some of which should gradually be
 * exsized from the program.  Originally, the idea was that Levels would immediately follow an ion, and in that case
 * all that was need to describe a level, was a level number, a multiplicity, and an energy for the level.  An example
 * of this type of line follows:
 *
 *  @verbatim
 *  Ion	He	2	1	1	24.587  6 2
 *  Level 1 1 0.0
 *  Level 2 3 19.819
 *  @endverbatim
 *
 * It is a requirement that when an OLDSTYLE line is used, it must immediately follow the ion.
 *
 * The second type of LEVEL style was developed solely to allow the levels to be disjoint from
 * the ions.  Disjoint in this case implies that the level must follow the ion, but it does not
 * have to be immediately after the ion line.
 *
 * It was developed to allow one to use levels developed from the Kurucz
 * line list, and the new Verner level file.
 *
 * A typical KURUCZSTYLE record would look like:
 *
 *  @verbatim
 *  Comment-- Ion H  1
 *  \# There are 17 unique levels for 1 1
 *  Level   1   1   0    2   0.000000
 *  Level   1   1   1    2  10.200121
 *  Level   1   1   2    4  10.200166
 *  Level   1   1   3    2  12.089051
 *  @endverbatim
 * where the colums are z, istate (in conventianal notation), a unique level no,
 * the multiplicity of the level and the exitation energy of the level in eV
 *
 *
 *
 * The new Topbase style records (created by py_parse_levels) look like this:
 *
 *
 *  @verbatim
 *  ======================================================================================
 *  #      i NZ NE iSLP iLV iCONF                 E(RYD)      TE(RYD)   gi     EQN    RL(NS)
 *  # ======================================================================================
 *  # ======================================================================================
 *  LevTop  1  1  200  1 -13.605698   0.000000  2  1.0000 1.00e+21 () 1s
 *  LevTop  1  1  200  2  -3.401425  10.204273  2  2.0000 1.00e+21 () 2s
 *  LevTop  1  1  211  1  -3.401425  10.204273  6  2.0000 1.60e-09 () 2p
 *  LevTop  1  1  220  1  -1.511743  12.093955 10  3.0000 1.55e-08 () 3d
 *  LevTop  1  1  200  3  -1.511743  12.093955  2  3.0000 1.59e-07 () 3s
 *  LevTop  1  1  211  2  -1.511743  12.093955  6  3.0000 5.28e-09 () 3p
 *  LevTop  2  1  100  1 -24.310389   0.000000  1  0.7481 1.00e+21 1s2
 *  LevTop  2  1  300  1  -4.766334  19.544041  3  1.6895 1.00e+21 1s 2s
 *  LevTop  2  1  100  2  -3.941516  20.368818  1  1.8579 1.00e+21 1s 2s
 *  @endverbatim
 *
 * The details are:
 *
 * | Col | Name    | Description
 * | --- | ------- | -----------
 * | 2   |  NZ     | the atomic number of the level
 * | 3   |  NE     | the ionization state according to the usual astronomical convention.
 * |     |         | Note that Topbase uses the number of electrons but this is fixed in my reformatting program.
 * | 4   | iSLP    | A coded version of the angular momentum
 * | 5   | iLV     | The level number (starts with 1)
 * | 6   | E(Ryd)  | The energy in eV  relative to the continuum
 * | 7   | TE(Ryd) | The energy in ev  realative to the ground state
 * | 8   | gi      | The multiplicity
 * | 9   | eqn     | The equivalent quantum number  (Not necessarily an integer)
 * | 10  | rl	     | The radiative lifetime
 * | 11  | config  | A string showing the configuration
 *
 *
 * \subsection matoms Macro Atoms (SS)
 *
 * When the macro atom method is to be used the level data should look like:
 * e.g.
 *
 *  @verbatim
 *  #         z ion lvl ion_pot   ex_energy  g  rad_rate
 *  LevMacro  1  1  1 -13.605698   0.000000  2  1.00e+21 () n=1
 *  LevMacro  1  1  2  -3.401425  10.204273  8  1.60e-09 () n=2
 *  LevMacro  1  1  3  -1.511743  12.093955 18  1.00e-08 () n=3
 *  LevMacro  1  1  4  -0.850356  12.755342 32  1.00e-08 () n=4
 *  LevMacro  1  2  1   0.000000  13.605698  1  1.00e+21 () cnt
 *  @endverbatim
 *
 *
 * The details are:
 *
 * | Col | Name      | Description
 * |---- | --------- | ------------
 * | 2   | z         | Atomic number
 * | 3   | ion       | Ionisation state (usual astronomical convention: neutral=1)
 * | 4   | lvl       | An index for the level: this index must be assigned consistently in all the macro input files
 * |     |           | (i.e. also in photoionisation data). Each level index in a particular ion must be unique.
 * | 5   | ion_pot   | Ionisation potential (electron volts)
 * | 6   | ex_energy | Excitation energy (electron volts). A consistent definition of energy=0 must be used for all ions of the same element:
 * |     |           | ground state energy for the neutral species = 0.
 * | 7   | g         | Statistical weight
 * | 8   | rad_rate  | Radiative lifetime
 * | 9   |           | A string giving a name for the level (not used by code)
 *
 *
 * */


        case 'N':
/*
    It's a non-lte level, i.e. one for which we are going to calculate populations, at least for some number of these.
	For these, we have to set aside space in the levden array in the plasma structure.  This is used for topbase
	photoionization and macro atoms
*/

/* ?? ksl This mix and match situation may be too much.  We are storing both macro level densities and so-called
topbase level densities in some of the same arrays in python.  Leave for now, but it may be difficult to keep
the program working in both cases, and certainly mixed cases  04apr ksl  */

/* 080810 -- ksl -- 62 -- I have changed the way levels are created so that one can only read one type
 * of levels for each ion.  Note also that all of the confiruations for a single ion need to be read together.
 * It will be possible to read other types of records but one should not mix levels of different ions (This
 * last bit is not actually new.
 */

          if (strncmp (word, "LevTop", 6) == 0)
          {                     //Its a TOPBASESTYLE level
            sscanf (aline,
                    "%*s %d %d %d %d %le %le %le %le %le %15c \n", &zz, &iistate, &islp, &ilv, &e, &exx, &ggg, &qqnum, &rl, configname);
            istate = iistate;
            z = zz;
            gg = ggg;
            exx *= EV2ERGS;     // Convert energy above ground to ergs
            mflag = -1;         //record that this is a LevTop not LevMacro read
            lev_type = 2;       // It's a topbase record
          }

          else if (strncmp (word, "LevMacro", 8) == 0)
          {                     //It's a Macro Atom level (SS)
            sscanf (aline, "%*s %d %d %d %le %le %le %le %15c \n", &zz, &iistate, &ilv, &e, &exx, &ggg, &rl, configname);
            islp = -1;          //these indices are not going to be used so just leave
            qqnum = -1;         //them at -1
            mflag = 1;          //record Macro read
            lev_type = 1;       // It's a Macro record
            istate = iistate;
            z = zz;
            gg = ggg;
            exx *= EV2ERGS;     // Convert energy to ergs
          }
          else
          {
            Error ("get_atomic_data: file %s line %d: Level line incorrectly formatted\n", file, lineno);
            Error ("Get_atomic_data: %s\n", aline);
            exit (0);
            return (0);
          }
// Now check that the ion for this level is already known.  If not break out
          n = 0;
          while ((ion[n].z != z || ion[n].istate != istate) && n < nions)
            n++;
          if (n == nions)
          {

            Debug ("get_atomic_data: file %s line %d has level for unknown ion \n", file, lineno);
            break;
          }

          /* Check that for a non-lte (macro atom) that the configuration level is not greater than was allowed for in the
           * ion line
           */

          if (lev_type == 1 && ilv > ion[n].n_lte_max)
          {
            Error ("get_atomic_data: macro level %d ge %d for z %d  istate %d\n", ilv, ion[n].n_lte_max, ion[n].z, ion[n].istate);
            //exit(0);
            break;
          }

/*  So now we know that this level can be associated with an ion

Now either set the type of level that will be used for this ion or set it if
a level type has not been established
		   */

          if (ion[n].lev_type == (-1))
          {
            ion[n].lev_type = lev_type;
          }
          else if (ion[n].lev_type != lev_type)
          {
            break;
          }

/*
 Now check 1) if it was a LevMacro that there isn't already a LevTop (if there was then
 something has gone wrong in the input order).
 2) if it was a LevTop that there wasn't already a LevMacro.If there was then ignore the new LevTop data
 for this ion. (SS)
*/


// Next steps should never happen; we have added a more robust mechanism to prevent any kind of mix and match above
          if (ion[n].macro_info == 1 && mflag == -1)
          {                     //it is already flagged as macro atom - current read is for LevTop - don't use it (SS)
            Error ("Get_atomic_data: file %s  Ignoring LevTop data for ion %d - already using Macro Atom data\n", file, n);
            break;
          }
          if (ion[n].macro_info == 0 && mflag == 1)
          {                     //It is already flagged as simple atom and this is  before MacroAtom data - so ignore.  ksl
            Error ("Get_atomic_data: file %s  Trying to read MacroAtom data after LevTop data for ion %d. Not allowed\n", file, n);
            break;
          }


          // case where data will be used (SS)
          if (mflag == 1)
          {
            config[nlevels].macro_info = 1;

            /* Extra check added here to be sure that the level emissivities used in the
               detailed spectrum calculation won't get messed up. The next loop should
               never trigger and can probably be deleted but I just want to check it for now.
               SS June 04. */

            if (nlevels_macro != nlevels)
            {
              Error ("get_atomicdata: Simple level has appeared before macro level. Not allowed.\n");
              exit (0);
            }
            nlevels_macro++;

            if (nlevels_macro > NLEVELS_MACRO)
            {
              Error ("get_atomicdata: Too many macro atom levels. Increase NLEVELS_MACRO. Abort. \n");
              exit (0);
            }
          }
          else
          {
            config[nlevels].macro_info = 0;
            nlevels_simple++;
          }

          config[nlevels].z = z;
          config[nlevels].istate = istate;
          config[nlevels].isp = islp;
          config[nlevels].ilv = ilv;
          config[nlevels].nion = n;     //Internal index to ion structure
          config[nlevels].q_num = qqnum;
          config[nlevels].g = gg;
          config[nlevels].ex = exx;
          config[nlevels].rad_rate = rl;
          /* SS Aug 2005
             Previously, the line above set the rad_rate to 0 and is was never used.
             Now I'm setting it to the radiative lifetime of the level.
             It will now be used in the macro atom calculation - if the lifetime is
             set to be long (infinite) in the input data, the level is assumed to be
             collisional supported by the ground state (i.e. has the LTE excitation
             fraction relative to ground).
           */


          if (ion[n].n_lte_max > 0)
          {                     // Then this ion wants nlte levels
            if (ion[n].first_nlte_level < 0)
            {                   // Then this is the first one that has been found
              ion[n].first_nlte_level = nlevels;
              ion[n].nlte = 1;
              config[nlevels].nden = ion[n].first_levden;
            }
            else if (ion[n].n_lte_max > ion[n].nlte)
            {
              config[nlevels].nden = ion[n].first_levden + ion[n].nlte;
              ion[n].nlte++;
            }
            else
            {
              config[nlevels].nden = -1;
            }
          }
          else
          {
            config[nlevels].nden = -1;
          }


/* Now associate this config with the levden array where appropriate.  The -1 is because ion[].nlte
is already incremented

*/

          if (ion[n].firstlevel < 0)
          {
            ion[n].firstlevel = nlevels;
            ion[n].nlevels = 1;
          }
          else
            ion[n].nlevels++;



          nlevels++;

          if (nlevels > NLEVELS)
          {
            Error ("getatomic_data: file %s line %d: More energy levels than allowed. Increase NLEVELS in atomic.h\n", file, lineno);
            exit (0);
          }
          break;

        case 'n':              // Its an "LTE" level

          if (sscanf (aline, "%*s %d %d %d %le %le\n", &zz, &iistate, &qnum, &gg, &exx) == 5)   //IT's KURUCZSTYLE
          {
            istate = iistate;
            z = zz;
            exx *= EV2ERGS;
            qqnum = ilv = qnum;
            lev_type = 0;       // It's a Kurucz-style record

          }
          else                  // Read an OLDSTYLE level description
          if (sscanf (aline, "%*s  %d %le %le\n", &qnum, &gg, &exx) == 3)
          {
            exx *= EV2ERGS;
            qqnum = ilv = qnum;
            lev_type = -2;      // It's an old style record, one which is only here for backward compatibility
          }
          else
          {
            Error ("get_atomic_data: file %s line %d: Level line incorrectly formatted\n", file, lineno);
            Error ("Get_atomic_data: %s\n", aline);
            exit (0);
            return (0);
          }
/* Check whether the ion for this level is known.  If not, skip the level */

// Next section is identical already to case N
          n = 0;
          while ((ion[n].z != z || ion[n].istate != istate) && n < nions)
            n++;
          if (n == nions)
          {

            Debug ("get_atomic_data: file %s line %d has level for unknown ion \n", file, lineno);
            break;

          }

          /* Now either set the type of level that will be used for this ion or set it if
           * a level type has not been established
           */

          if (ion[n].lev_type == (-1))
          {
            ion[n].lev_type = lev_type;
          }
          else if (ion[n].lev_type != lev_type)
          {

            break;
          }
//  End section known to be idential to case N


/* Check whether this is a macro-ion.  If it is a macro-ion, but the level appears to be described as a
simple level (i.e without a keyword LeVMacro), then skip it, since a macro-ion has to have all the levels
described as macro-levels. */
          if (ion[n].macro_info == 1)
          {
            Error ("get_atomic_data: file %s line %d has simple level for ion[%d], which is a macro-ion\n", file, lineno, n);
            break;
          }
/* Check to prevent one from adding simple levels to an ionized that already has some nlte levels.  Note that
   an ion may have simple levels, i.e. levels with no entries in the plasma structure levden array, but this
   will only be the case if there are too many of this type of level.
*/
          if (ion[n].nlte > 0)
          {
            Error ("get_atomic_data:  file %s line %d has simple level for ion[%d], which has non_lte_levels\n", file, lineno, n);
            break;
          }

/*  Check whether we already have too many levels specified for this ion. If so, skip */
          if (ion[n].nmax == ion[n].nlevels)
          {

            Debug ("get_atomic_data: file %s line %d has level exceeding the number allowed for ion[%d]\n", file, lineno, n);

            break;
          }
//  So now we know that this level can be associated with an ion

          config[nlevels].z = z;
          config[nlevels].istate = istate;
          config[nlevels].isp = islp;
          config[nlevels].ilv = ilv;
          config[nlevels].nion = n;     //Internal index to ion structure
          config[nlevels].q_num = qqnum;
          config[nlevels].g = gg;
          config[nlevels].ex = exx;
          if (ion[n].firstlevel < 0)
          {
            ion[n].firstlevel = nlevels;
            ion[n].nlevels = 1;
          }
          else
            ion[n].nlevels++;


/* Now declare that this level has no corresponding element in the levden array which is part
   of the plasma stucture.  To do this set config[].ndent to -1
*/

          config[nlevels].nden = -1;

          config[nlevels].rad_rate = 0.0;       // ?? Set emission oscillator strength for the level to zero

          nlevels_simple++;
          nlevels++;
          if (nlevels > NLEVELS)
          {
            Error ("getatomic_data: file %s line %d: More energy levels than allowed. Increase NLEVELS in atomic.h\n", file, lineno);
            exit (0);
          }
          break;



/**
 * @section  Photoionization
 *
 * Until at least Oct 2001, Python used photoionization crossections from Verner, Ferland, Korista, and Yakolev (DFKY)
 * The routine sigma_phot(xptr, freq) calculates the crossection based on this.
 *
 *
 * The original topbase records look like this
 *
 * @verbatim
 *  ================================================
 *        I  NZ  NE  ISLP  ILV        E(RYD)      NP
 *  ================================================
 *        1   2   1   200    1  -4.00000E+00     101
 *   3.960000E+00 1.619E+00
 *   4.000000E+00 1.576E+00
 *   4.101260E+00 1.474E+00
 *   4.205084E+00 1.379E+00
 *   4.311536E+00 1.289E+00
 *   4.420684E+00 1.206E+00
 * @endverbatim
 *
 * They are converted to something that is more compatible with Topbase
 * by py_top_phot
 *
 * The new topbase style records look like this.
 *
 * @verbatim
 * PhotTopS  2  1 200    1    54.422791 101
 * PhotTop    54.422791 1.576e-18
 * PhotTop    66.508772 9.197e-19
 * PhotTop    78.594753 5.817e-19
 * PhotTop    90.680734 3.906e-19
 * PhotTop   102.766715 2.746e-19
 * PhotTop   114.852696 2.001e-19
 * PhotTop   126.938677 1.502e-19
 * PhotTop   139.024658 1.155e-19
 * PhotTop   151.110639 9.057e-20
 * @endverbatim
 *
 * The main changes from the original records are that energy levels
 * have been converted to eV, cross sections to cm**2, and the electron
 * number has been converted to conventional astronomical notiation
 * for the ionstate.
 *
 * ### Macro atoms (SS)
 *
 * For the Macro Atom method the input photoionisation data should look like
 * (one entry for every photoionisation process):
 *
 * @verbatim
 * z  ion ll ul   threshold  npts
 * PhotMacS  1  1  1  1    13.605698  50
 * PhotMac    13.605698 6.304e-18
 * PhotMac    16.627193 3.679e-18
 * PhotMac    19.648688 2.327e-18
 * PhotMac    22.670183 1.563e-18
 * PhotMac    25.691679 1.098e-18
 * PhotMac    28.713174 8.004e-19
 * PhotMac    31.734669 6.006e-19
 * PhotMac    34.756165 4.618e-19
 * @endverbatim
 *
 * Details for the first row of each entry:
 * | Col | Name  | Description
 * | --- | ----- | -----------
 * | 2   | z     | Atomic number
 * | 3   | ion   | Ionisation stage (astronomical convention): LOWER ION (i.e. the one that GETS IONISED)
 * |     |       | it is always assumed that the upper stage is this+1!
 * | 4   | ll    | index for the level in the LOWER ION (must match the index given to the level in input level data)
 * | 5   | ul    | index for the level in the UPPER ION (must match the index given to the level in the input level data)
 * | 6   | thresh| threshold energy of edge (in eV)
 * | 7   | npts  | number of data points in the following table
 *
 * * Then follows the table of energy and cross-section values (same as TopBase above)
 * * 04dec	ksl	Modified this section so that "unknown" level information is skipped so that
 *   		one need modify only the higher level elements_ions file
 */

        case 'w':
          if (strncmp (word, "PhotMacS", 8) == 0)
          {
            // It's a Macro atom entry - similar format to TOPBASE - see below (SS)
            sscanf (aline, "%*s %d %d %d %d %le %d\n", &z, &istate, &levl, &levu, &exx, &np);
            islp = -1;
            ilv = -1;

            for (n = 0; n < np; n++)
            {
              //Read the photo. records but do nothing with them until verifyina a valid level
              if (fgets (aline, LINELENGTH, fptr) == NULL)
              {
                Error ("Get_atomic_data: Problem reading topbase photoionization record\n");
                Error ("Get_atomic_data: %s\n", aline);
                exit (0);
              }
              sscanf (aline, "%*s %le %le", &xe[n], &xx[n]);
              lineno++;
            }

            // Locate upper state
            n = 0;
            while ((config[n].z != z || config[n].istate != (istate + 1)        //note that the upper config will (SS)
                    || config[n].ilv != levu) && n < nlevels)   //be the next ion up (istate +1) (SS)
              n++;
            if (n == nlevels)
            {
              Error_silent ("get_atomic_data: No configuration found to match upper state for phot. line %d\n", lineno);
              break;            //Need to match the configuration for macro atoms - break if not found.
            }


            // Locate lower state
            m = 0;
            while ((config[m].z != z || config[m].istate != istate      //Now searching for the lower
                    || config[m].ilv != levl) && m < nlevels)   //configuration (SS)
              m++;
            if (m == nlevels)
            {
              Error_silent ("get_atomic_data: No configuration found to match lower state for phot. line %d\n", lineno);
              break;            //Need to match the configuration for macro atoms - break if not found.
            }

            // Populate upper state info
            phot_top[ntop_phot].uplev = n;      //store the level in the upper ion (SS)
            config[n].bfd_jump[config[n].n_bfd_jump] = ntop_phot;       //record the line index as a downward bf Macro Atom jump (SS)
            phot_top[ntop_phot].down_index = config[n].n_bfd_jump;      //record jump index in the photoionization structure
            config[n].n_bfd_jump += 1;  //note that there is one more downwards bf jump available (SS)
            if (config[n].n_bfd_jump > NBFJUMPS)
            {
              Error ("get_atomic_data: Too many downward b-f jump for ion %d\n", config[n].istate);
              exit (0);
            }


            // Populate lower state info
            phot_top[ntop_phot].nlev = m;       //store lower configuration then find upper configuration(SS)
            config[m].bfu_jump[config[m].n_bfu_jump] = ntop_phot;       //record the line index as an upward bf Macro Atom jump (SS)
            phot_top[ntop_phot].up_index = config[m].n_bfu_jump;        //record the jump index in the photoionization structure
            config[m].n_bfu_jump += 1;  //note that there is one more upwards bf jump available (SS)
            if (config[m].n_bfu_jump > NBFJUMPS)
            {
              Error ("get_atomic_data: Too many upward b-f jump for ion %d\n", config[m].istate);
              exit (0);
            }


            phot_top[ntop_phot].nion = config[m].nion;
            phot_top[ntop_phot].z = z;
            phot_top[ntop_phot].istate = istate;
            phot_top[ntop_phot].np = np;
            phot_top[ntop_phot].nlast = -1;
            phot_top[ntop_phot].macro_info = 1;

            if (ion[config[m].nion].phot_info == -1)
            {
              ion[config[m].nion].phot_info = 1;        /* Mark this ion as using TOPBASE photo */
              ion[config[m].nion].ntop_first = ntop_phot;
            }

            /* next line sees if the topbase level just read in is the ground state -
               if it is, the ion structure element ntop_ground is set to that topbase level number
               note that m is the lower level here */
            if (m == config[ion[config[n].nion].first_nlte_level].ilv)
            {
              ion[config[n].nion].ntop_ground = ntop_phot;
            }

            ion[config[m].nion].ntop++;

            // Finish up this section by storing the photionization data properly

            for (n = 0; n < np; n++)
            {
              phot_top[ntop_phot].freq[n] = xe[n] * EV2ERGS / H;        // convert from eV to freqency
              phot_top[ntop_phot].x[n] = xx[n]; // leave cross sections in  CGS
            }
            if (phot_freq_min > phot_top[ntop_phot].freq[0])
              phot_freq_min = phot_top[ntop_phot].freq[0];


            ntop_phot_macro++;
            ntop_phot++;
            nphot_total++;

            if (nphot_total > NTOP_PHOT)
            {
              Error ("get_atomicdata: More macro photoionization cross sections that NTOP_PHOT (%d).  Increase in atomic.h\n", NTOP_PHOT);
              exit (0);
            }
            break;
          }



          else if (strncmp (word, "PhotTopS", 8) == 0)
          {
            // It's a TOPBASE style photoionization record, beginning with the summary record
            sscanf (aline, "%*s %d %d %d %d %le %d\n", &z, &istate, &islp, &ilv, &exx, &np);
            for (n = 0; n < np; n++)
            {                   //Read the topbase photoionization records
              if (fgets (aline, LINELENGTH, fptr) == NULL)
              {
                Error ("Get_atomic_data: Problem reading topbase photoionization record\n");
                Error ("Get_atomic_data: %s\n", aline);
                exit (0);
              }
              sscanf (aline, "%*s %le %le", &xe[n], &xx[n]);
              lineno++;

            }

            n = 0;


            /* additional check to assure that records were
             * only matched with levels whose density was being tracked in levden.  This
             * is now necesary since a change was made to use topbase levels for calculating
             * partition functions
             */

            while ((config[n].nden == -1
                    || config[n].z != z || config[n].istate != istate || config[n].isp != islp || config[n].ilv != ilv) && n < nlevels)
              n++;
            if (n == nlevels)
            {

              Debug ("No level found to match PhotTop data in file %s on line %d. Data ignored.\n", file, lineno);
              break;            // There was no pre-existing ion
            }
            if (ion[config[n].nion].macro_info == 0)    //this is not a macro atom level (SS)
            {
              phot_top[ntop_phot].nlev = n;     // level associated with this crossection.
              phot_top[ntop_phot].nion = config[n].nion;
              phot_top[ntop_phot].z = z;
              phot_top[ntop_phot].istate = istate;
              phot_top[ntop_phot].np = np;
              phot_top[ntop_phot].nlast = -1;
              phot_top[ntop_phot].macro_info = 0;

              /* next line sees if the topbase level just read in is the ground state -
                 if it is, the ion structure element ntop_ground is set to that topbase level number */
              if (islp == config[ion[config[n].nion].first_nlte_level].isp && ilv == config[ion[config[n].nion].first_nlte_level].ilv)
              {
                ion[config[n].nion].ntop_ground = ntop_phot;
              }


              if (ion[config[n].nion].phot_info == -1)
              {
                ion[config[n].nion].phot_info = 1;      /* Mark this ion as using TOPBASE photo */
                ion[config[n].nion].ntop_first = ntop_phot;

              }
              else if (ion[config[n].nion].phot_info == (0))
              {
                Error
                  ("Get_atomic_data: file %s VFKY and Topbase photoionization x-sections in wrong order for nion %d\n",
                   file, config[n].nion);
                Error ("             Read topbase x-sections before VFKY if using both types!!\n");
                exit (0);
              }
              ion[config[n].nion].ntop++;
              for (n = 0; n < np; n++)
              {
                phot_top[ntop_phot].freq[n] = xe[n] * EV2ERGS / H;      // convert from eV to freqency

                phot_top[ntop_phot].x[n] = xx[n];       // leave cross sections in  CGS
              }
              if (phot_freq_min > phot_top[ntop_phot].freq[0])
                phot_freq_min = phot_top[ntop_phot].freq[0];


              ntop_phot_simple++;
              ntop_phot++;
              nphot_total++;

              /* check to assure we did not exceed the allowed number of photoionization records */
              if (nphot_total > NTOP_PHOT)
              {
                Error
                  ("get_atomicdata: More TopBase photoionization cross sections that NTOP_PHOT (%d).  Increase in atomic.h\n", NTOP_PHOT);
                exit (0);
              }
            }
            else
            {
              Error
                ("Get_atomic_data: photoionisation data ignored since previously read Macro Atom input for the same ion. File: %s line: %d \n",
                 file, lineno);
            }
            break;
          }


          /* Check that there is an ion which has the same ionization state as this record
             otherwise it must be a VFKY style record and so read with that format */

          else if (strncmp (word, "PhotVfkyS", 8) == 0)
          {
            // It's a VFKY style photoionization record, beginning with the summary record
            sscanf (aline, "%*s %d %d %d %d %le %d\n", &z, &istate, &islp, &ilv, &exx, &np);
            for (n = 0; n < np; n++)
            {
              //Read the topbase photoionization records
              if (fgets (aline, LINELENGTH, fptr) == NULL)
              {
                Error ("Get_atomic_data: Problem reading Vfky photoionization record\n");
                Error ("Get_atomic_data: %s\n", aline);
                exit (0);
              }
              sscanf (aline, "%*s %le %le", &xe[n], &xx[n]);
              lineno++;

            }

            for (nion = 0; nion < nions; nion++)
            {
              if (ion[nion].z == z && ion[nion].istate == istate && ion[nion].macro_info != 1)
              {
                if (ion[nion].phot_info == -1)
                {
                  /* Then there is a match */
                  phot_top[nphot_total].nlev = ion[nion].firstlevel;    // ground state
                  phot_top[nphot_total].nion = nion;
                  phot_top[nphot_total].z = z;
                  phot_top[nphot_total].istate = istate;
                  phot_top[nphot_total].np = np;
                  phot_top[nphot_total].nlast = -1;
                  phot_top[nphot_total].macro_info = 0;

                  ion[nion].phot_info = 0;      /* Mark this ion as using VFKY photo */
                  ion[nion].nxphot = nphot_total;

                  for (n = 0; n < np; n++)
                  {
                    phot_top[nphot_total].freq[n] = xe[n] * EV2ERGS / H;        // convert from eV to freqency
                    phot_top[nphot_total].x[n] = xx[n]; // leave cross sections in  CGS
                  }
                  if (phot_freq_min > phot_top[ntop_phot].freq[0])
                    phot_freq_min = phot_top[ntop_phot].freq[0];
                  nxphot++;
                  nphot_total++;
                }

                else if (ion[nion].phot_info == 1 && ion[nion].macro_info != 1)
                  /* We already have a topbase cross section, but the VFKY
                     data is superior for the ground state, so we replace that data with the current data
                     JM 1508 -- don't do this with macro-atoms for the moment */
                {
                  phot_top[ion[nion].ntop_ground].nlev = ion[nion].firstlevel;  // ground state
                  phot_top[ion[nion].ntop_ground].nion = nion;
                  phot_top[ion[nion].ntop_ground].z = z;
                  phot_top[ion[nion].ntop_ground].istate = istate;
                  phot_top[ion[nion].ntop_ground].np = np;
                  phot_top[ion[nion].ntop_ground].nlast = -1;
                  phot_top[ion[nion].ntop_ground].macro_info = 0;
                  ion[nion].phot_info = 2;      //We mark this as having hybrid data - VFKY ground, TB excited, potentially VFKY innershell
                  for (n = 0; n < np; n++)
                  {
                    phot_top[ion[nion].ntop_ground].freq[n] = xe[n] * EV2ERGS / H;      // convert from eV to freqency
                    phot_top[ion[nion].ntop_ground].x[n] = xx[n];       // leave cross sections in  CGS
                  }
                  if (phot_freq_min > phot_top[ion[nion].ntop_ground].freq[0])
                    phot_freq_min = phot_top[ion[nion].ntop_ground].freq[0];
                  Debug
                    ("Get_atomic_data: file %s  Replacing ground state topbase photoionization for ion %d with VFKY photoionization\n",
                     file, nion);
                }
              }
            }

            if (nxphot > NIONS)
            {
              Error ("getatomic_data: file %s line %d: More photoionization edges than IONS.\n", file, lineno);
              exit (0);
            }
            if (nphot_total > NTOP_PHOT)
            {
              Error ("get_atomicdata: More photoionization cross sections that NTOP_PHOT (%d).  Increase in atomic.h\n", NTOP_PHOT);
              exit (0);
            }

            break;
          }
          else
          {
            Error ("get_atomic_data: file %s line %d: photoionization line incorrectly formatted\n", file, lineno);
            Log ("Make sure you are using the tabulated verner cross sections (photo_vfky_tabulated.data)\n");
            Error ("Get_atomic_data: %s\n", aline);
            exit (0);
          }

          /* Input inner shell cross section data */



        case 'I':
          if (sscanf (aline, "%*s %d %d %d %d %le %d\n", &z, &istate, &in, &il, &exx, &np) != 6)
          {
            Error ("Inner shell ionization data incorrectly formatted\n");
            Error ("Get_atomic_data: %s\n", aline);
            exit (0);
          }
          for (n = 0; n < np; n++)
          {
            //Read the topbase photoionization records
            if (fgets (aline, LINELENGTH, fptr) == NULL)
            {
              Error ("Get_atomic_data: Problem reading VY inner shell record\n");
              Error ("Get_atomic_data: %s\n", aline);
              exit (0);
            }
            sscanf (aline, "%*s %le %le", &xe[n], &xx[n]);
            lineno++;
          }
          for (nion = 0; nion < nions; nion++)
          {
            if (ion[nion].z == z && ion[nion].istate == istate && ion[nion].macro_info != 1)
            {
              /* Then there is a match */
              inner_cross[n_inner_tot].nlev = ion[nion].firstlevel;     //All these are for the ground state
              inner_cross[n_inner_tot].nion = nion;
              inner_cross[n_inner_tot].np = np;
              inner_cross[n_inner_tot].z = z;
              inner_cross[n_inner_tot].istate = istate;
              inner_cross[n_inner_tot].n = in;
              inner_cross[n_inner_tot].l = il;
              inner_cross[n_inner_tot].nlast = -1;
              ion[nion].n_inner++;      /*Increment the number of inner shells */
              ion[nion].nxinner[ion[nion].n_inner] = n_inner_tot;
              for (n = 0; n < np; n++)
              {
                inner_cross[n_inner_tot].freq[n] = xe[n] * EV2ERGS / H; // convert from eV to freqency
                inner_cross[n_inner_tot].x[n] = xx[n];  // leave cross sections in  CGS
              }
              if (inner_freq_min > inner_cross[n_inner_tot].freq[0])
                inner_freq_min = inner_cross[n_inner_tot].freq[0];
              n_inner_tot++;

            }
          }
          if (n_inner_tot > N_INNER * NIONS)
          {
            Error ("getatomic_data: file %s line %d: Inner edges than we have room for.\n", file, lineno);
            exit (0);
          }
          break;


          /*Input data for innershell ionization followed by
             Auger effect */
/**
 * @section Auger
 */
/*
		case 'A':
		  if (sscanf (aline,
			      "%*s %d %d %d %d %le %le %le %le %le %le %le",
			      &z, &istate, &nn, &nl, &yield, &arad, &etarad,
			      &adi, &t0di, &bdi, &t1di) != 11)
		    {
		      Error ("Auger input incorrectly formatted\n");
		      Error ("Get_atomic_data: %s\n", aline);
		      exit (0);
		    }
		  if (nauger < NAUGER)
		    {
		      if ((vptr =
			   fopen ("atomic/photo_verner.data", "r")) == NULL)
			{
			  Error
			    ("get_atomic data:  Could not open photo_verner.data\n");
			  exit (0);
			}
		      ion_index = -2;
		      target_index = -1;
		      for (n_verner = 0; n_verner < 1696; n_verner++)
			{
			  fscanf (vptr,
				  "%d %d %d %d %le %le %le %le %le %le\n",
				  &dumz, &dumistate, &dumnn, &dumnl, &dumE_th,
				  &dumE_0, &dumSigma, &dumya, &dumP, &dumyw);
			  if ((dumz == z)
			      && (dumistate == (dumz - istate + 1))
			      && (dumnn == nn) && (dumnl == nl))
			    {
			      // Now need to check that this ion is really in the data set and find which it is
			      ion_index = -1;
			      for (n = 0; n < nions; n++)
				{
				  if (ion[n].z == z
				      && ion[n].istate == istate)
				    {
				      ion_index = n;
				      augerion[nauger].nion = n;
				      augerion[nauger].yield = yield;
				      augerion[nauger].arad = arad;
				      augerion[nauger].etarad = etarad;
				      augerion[nauger].adi = adi;
				      augerion[nauger].t0di = t0di;
				      augerion[nauger].bdi = bdi;
				      augerion[nauger].t1di = t1di;

				      augerion[nauger].z = z;
				      augerion[nauger].istate = istate;
				      augerion[nauger].n = nn;
				      augerion[nauger].l = nl;
				      augerion[nauger].freq_t = dumE_th / HEV;
				      augerion[nauger].E_0 = dumE_0;
				      augerion[nauger].Sigma = dumSigma * 1.e-18;	//input in megabarns
				      augerion[nauger].ya = dumya;
				      augerion[nauger].yw = dumyw;
				      augerion[nauger].P = dumP;
				      nauger++;
				    }

				  //We also want the index for the
				     targe ion (i.e. the one that is
				     two ionization stages up from the
				     one found above
				  if (ion[n].z == z
				      && ion[n].istate == istate + 2)
				    {
				      target_index = n;
				    }

				}
			      if (ion_index == -1)
				{
				  Error
				    ("Get_atomic_data: Failed to find ion to match Auger input data. Ignoring.  %d %d %d %d\n",
				     dumz, dumistate, dumnn, dumnl);
				}
			      else
				{
				  augerion[nauger - 1].nion_target =
				    target_index;
				}
			    }
			}
		      fclose (vptr);
		      if (ion_index == -2)
			{
			  Error
			    ("Get_atomic_data: Failed to find source data to match Auger input data. Ignoring. %d %d %d %d\n",
			     z, istate, nn, nl);
			}
		      else
			{
			  Log
			    ("Matched Auger ionization input to Z %d istate %d going to Z %d istate %d.\n",
			     ion[augerion[nauger - 1].nion].z,
			     ion[augerion[nauger - 1].nion].istate,
			     ion[augerion[nauger - 1].nion_target].z,
			     ion[augerion[nauger - 1].nion_target].istate);
			}
		    }
		  else
		    {
		      Error
			("Get_atomic_data: NAUGER is filled up! Ignoring input data %d.\n",
			 nauger);
		    }
		  break;
*/


/**
 * @section Lines
 *
 * here are various types of line files which get_atomicdata must parse
 * This is the original format
 *
 * @verbatim
 * Line 1 1 1215.673584 0.139000  2  2
 * Line 1 1 1215.668213 0.277000  2  4
 * @endverbatim
 *
 * This is the first format developed for including excited state lines.  It adds the energy of the lower
 * and upper level
 *
 * @verbatim
 * Line  6  3 1174.933000    0.115845   3   5    6.496462   17.050273
 * Line  6  3 1175.263000    0.282488   1   3    6.493525   17.044369
 * Line  6  3 1175.590000    0.069804   3   3    6.496462   17.044369
 * @endverbatim
 *
 * Finally there is a new format developed to link lines to configurations.
 *
 * @verbatim
 * Line  1  1 1215.673584  0.139000   2   2     0.000000    10.200121    0    1
 * Line  1  1 1215.668213  0.277000   2   4     0.000000    10.200166    0    2
 * @endverbatim
 *
 *
 * ### Macro Atoms (SS)
 *
 * For the macro atom method the lines are input in the following format:
 *
 *
 * @verbatim
 * # z = element, ion= ionstage, f = osc. str., gl(gu) = stat. we. lower(upper) level
 * # el(eu) = energy lower(upper) level (eV), ll(lu) = lvl index lower(upper) level
 * #          z ion lambda      f         gl  gu    el          eu          ll   lu
 * LinMacro  1  1 1215.671021  0.416200   2   8     0.000000    10.200143    1    2
 * LinMacro  1  1 1025.722046  0.079100   2  18     0.000000    12.089062    1    3
 * LinMacro  1  1  972.540000  0.028990   2  32     0.000000    12.755342    1    4
 * LinMacro  1  1 6562.800000  0.640700   8  18    10.200143    12.089062    2    3
 * LinMacro  1  1 4861.320000  0.119300   8  32    10.200143    12.755342    2    4
 * LinMacro  1  1 18751.00000  0.842100  18  32    12.089062    12.755342    3    4
 * @endverbatim
 *
 * Details:
 *
 * | Col | Description
 * | --- | -----------
 * | 2   | Atomic Number
 * | 3   | Ionisation stage (astronomical convention)
 * | 4   | Wavelength of transition (AA)
 * | 5   | oscillator strength
 * | 6   | stat. weight for lower level of transition
 * | 7   | stat. weight for upper level of transition
 * | 8   | energy of lower level of transition (eV)
 * | 9   | energy of upper level of transition (eV)
 * | 10  | index for lower level of transition (MUST match the index assigned to this level in level data
 * | 11  | index for upper level of transition (MUST match the index assigned to this level in level data
 *
 *   04dec ksl -- I have modified the next section to try to conform to the "philosophy" of
 *   get_atomic_data which is that one does not need to modify the more detailed files, e.g
 *   the LinMacro file, if one has eliminated a particular level in the elements ion file.
 *   Basically this was accomplished by checking both the upper and lower level and breaking
 *   out if either was not accounted for.
*/
        case 'r':
          if (strncmp (word, "LinMacro", 8) == 0)
          {                     //It's a macro atoms line(SS)
            if (mflag != 1)
            {
              Error ("get_atomicdata: Can't read macro-line after some simple lines. Reorder the input files!\n");
              exit (0);
            }

            mflag = 1;          //flag to identify macro atom case (SS)
            nwords = sscanf (aline, "%*s %d %d %le %le %le %le %le %le %d %d", &z, &istate, &freq, &f, &gl, &gu, &el, &eu, &levl, &levu);
            if (nwords != 10)
            {
              Error ("get_atomic_data: file %s line %d: LinMacro line incorrectly formatted\n", file, lineno);
              Error ("Get_atomic_data: %s\n", aline);
              exit (0);
            }

            el = EV2ERGS * el;
            eu = EV2ERGS * eu;
            //need to identify the configurations associated with the upper and lower levels (SS)
            n = 0;
            while ((config[n].z != z || config[n].istate != istate || config[n].ilv != levl) && n < nlevels)
              n++;
            if (n == nlevels)
            {
              Error_silent ("Get_atomic_data: No configuration found to match lower level of line %d\n", lineno);
              break;
            }


            m = 0;
            while ((config[m].z != z || config[m].istate != istate || config[m].ilv != levu) && m < nlevels)
              m++;
            if (m == nlevels)
            {
              Error_silent ("Get_atomic_data: No configuration found to match upper level of line %d\n", lineno);
              break;
            }

            /* Now that we know this is a valid transition for the macro atom record the data */

            nconfigl = n;       //record lower configuration (SS)
            config[n].bbu_jump[config[n].n_bbu_jump] = nlines;  //record the line index as an upward bb Macro Atom jump(SS)
            line[nlines].down_index = config[n].n_bbu_jump;     //record the index for the jump in the line structure
            config[n].n_bbu_jump += 1;  //note that there is one more upwards jump available (SS)
            if (config[n].n_bbu_jump > NBBJUMPS)
            {
              Error ("get_atomic_data: Too many upward b-b jumps for ion %d\n", config[n].istate);
              exit (0);
            }

            nconfigu = m;       //record upper configuration (SS)
            config[m].bbd_jump[config[m].n_bbd_jump] = nlines;  //record the line index as a downward bb Macro Atom jump (SS)
            line[nlines].up_index = config[m].n_bbd_jump;       //record jump index in line structure
            config[m].n_bbd_jump += 1;  //note that there is one more downwards jump available (SS)
            if (config[m].n_bbd_jump > NBBJUMPS)
            {
              Error ("get_atomic_data: Too many downward b-b jumps for ion %d\n", config[m].istate);
              exit (0);
            }


          }
          else
          {                     //It's not a macro atom line (SS)
// It would have been better to define mflag = 0 since that is what we want to set
// macro_info to if it is an old-style line, but keep it this way for now.  ksl
            mflag = -1;         //a flag to mark this as not a macro atom case (SS)
            nconfigl = -1;
            nconfigu = -1;
            nwords = sscanf (aline, "%*s %d %2d %le %le %le %le %le %le %d %d", &z, &istate, &freq, &f, &gl, &gu, &el, &eu, &levl, &levu);
            if (nwords == 6)
            {
              el = 0.0;
              eu = H * C / (freq * 1e-8);       // Convert Angstroms to ergs
              levl = -1;
              levu = -1;

            }
            else if (nwords == 8)
            {                   // Then the file contains the energy levels of the transitions

              el = EV2ERGS * el;
              eu = EV2ERGS * eu;
              levl = -1;
              levu = -1;
            }
            else if (nwords == 10)
            {                   // Then the file contains energy levels and level numbers
              el = EV2ERGS * el;
              eu = EV2ERGS * eu;
            }

            else
            {
              Error ("get_atomic_data: file %s line %d: Resonance line incorrectly formatted\n", file, lineno);
              Error ("Get_atomic_data: %s\n", aline);
              exit (0);
            }
          }

          if (el > eu)
            Error ("get_atomic_data: file %s line %d : line has el (%f) > eu (%f)\n", file, lineno, el, eu);
          for (n = 0; n < nions; n++)
          {
            if (ion[n].z == z && ion[n].istate == istate)
            {                   /* Then there is a match */
              if (freq == 0 || f <= 0 || gl == 0 || gu == 0)
              {
                Error_silent ("getatomic_data: line input incomplete: %s\n", aline);
                break;
              }
              //
              //define macro atom case (SS)
/* XXXX  04 April ksl -- Right now have enforced a clean separation between macro-ions and simple-ions
but this is proably not what we want if we move all bf & fb transitions to macro-ion approach.  We
would like to have simple lines for macro-ions */
              if (ion[n].macro_info == 1 && mflag == -1)
              {
                /* count how many times this happens to report to user */
                simple_line_ignore[n] += 1;
                break;
              }

              if (ion[n].macro_info == -1 && mflag == 1)
              {
                Error ("Getatomic_data: Macro Atom line data supplied for ion %d\n but there is no suitable level data\n", n);
                exit (0);
              }
              line[nlines].nion = n;
              line[nlines].z = z;
              line[nlines].istate = istate;
              line[nlines].freq = C / (freq * 1e-8);    /* convert Angstroms to frequency */
              line[nlines].f = f;
              line[nlines].gl = gl;
              line[nlines].gu = gu;
              line[nlines].levl = levl;
              line[nlines].levu = levu;
              line[nlines].el = el;
              line[nlines].eu = eu;
              line[nlines].nconfigl = nconfigl;
              line[nlines].nconfigu = nconfigu;
              line[nlines].coll_index = -999;   //Tokick off with we assume there is no collisional strength data
              if (mflag == -1)
              {
                line[nlines].macro_info = 0;    // It's an old-style line`
                nlines_simple++;
              }
              else
              {
                line[nlines].macro_info = 1;    //It's a macro line
                nlines_macro++;
              }
              nlines++;
            }
          }
          if (nlines > NLINES)
          {
            Error ("getatomic_data: file %s line %d: More lines than allowed. Increase NLINES in atomic.h\n", file, lineno);
            exit (0);
          }
          break;

/** @section Ground state fractions
 */
        case 'f':
          if (sscanf
              (aline,
               "%*s %d %d %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le",
               &z, &istate, &the_ground_frac[0], &the_ground_frac[1],
               &the_ground_frac[2], &the_ground_frac[3],
               &the_ground_frac[4], &the_ground_frac[5],
               &the_ground_frac[6], &the_ground_frac[7],
               &the_ground_frac[8], &the_ground_frac[9],
               &the_ground_frac[10], &the_ground_frac[11],
               &the_ground_frac[12], &the_ground_frac[13],
               &the_ground_frac[14], &the_ground_frac[15],
               &the_ground_frac[16], &the_ground_frac[17], &the_ground_frac[18], &the_ground_frac[19]) != 22)
          {
            Error ("get_atomic_data: file %s line %d ground state fracs   frac table incorrectly formatted\n", file, lineno);
            Error ("Get_atomic_data: %s\n", aline);
            exit (0);
          }
          for (n = 0; n < nions; n++)
          {
            if (ion[n].z == z && ion[n].istate == istate)
            {                   /* Then there is a match */
              ground_frac[n].z = z;
              ground_frac[n].istate = istate;
              for (j = 0; j < 20; j++)
              {
                ground_frac[n].frac[j] = the_ground_frac[j];
              }
            }
          }
          break;



/** @section Dielectronic Recombination - type 1
 * This section reads in dielectronic recombination rates.
 * The data comes from Chianti, and there are two types of dielectronic data.
 * This is type 1 - which has fis decribed by Arnaud & Raymond (1992)
 * There are two fitting parameters, E and C which represent different
 * degrees in the fit - one sums over a set of parameters.
 *
 * A typical ion data would look like this
 *
 *
 @verbatim
 * DR_BADNL E 2 2 4.5560e+05 5.5520e+05 8.9820e+05 0.0000e+00 0.0000e+00 0.0000e+00 0.0000e+00 0.0000e+00
 * DR_BADNL C 2 2 5.9660e-04 1.6130e-04 -2.2230e-05 0.0000e+00 0.0000e+00 0.0000e+00 0.0000e+00 0.0000e+00
 @endverbatim
 *
 * where there are two lines for a given ion - this is for helium (z=2) 2 (istate=2)
 *
 */

        case 'D':              /* Dielectronic recombination data read in. */
          nparam = sscanf (aline, "%*s %s %d %d %le %le %le %le %le %le %le %le %le", &drflag, &z, &ne, &drp[0], &drp[1], &drp[2], &drp[3], &drp[4], &drp[5], &drp[6], &drp[7], &drp[8]);       //split and assign the line
          nparam -= 3;          //take 4 off the nparam to give the number of actual parameters
          if (nparam > 9 || nparam < 1) //     trap errors - not as robust as usual because there are a varaible number of parameters...
          {
            Error ("Something wrong with dielectronic recombination data\n", file, lineno);
            Error ("Get_atomic_data: %s\n", aline);
            exit (0);
          }

          istate = ne;          //         get the ionisation state we are recombining from

          for (n = 0; n < nions; n++)   //Loop over ions to find the correct place to put the data
          {
            if (ion[n].z == z && ion[n].istate == istate)       // this works out which ion we are dealing with
            {
              if (ion[n].drflag == 0)   //This is the first time we have dealt with this ion
              {
                drecomb[ndrecomb].nion = n;     //put the ion number into the DR structure
                drecomb[ndrecomb].nparam = nparam;      //Put the number of parameters we ware going to read in, into the DR structure so we know what to iterate over later
                ion[n].nxdrecomb = ndrecomb;    //put the number of the DR into the ion
                drecomb[ndrecomb].type = DRTYPE_BADNELL;        //define the type of data
                ndrecomb++;     //increment the counter of number of dielectronic recombination parameter sets
                ion[n].drflag++;        //increment the flag by 1. We will do this rather than simply setting it to 1 so we will get errors if we do this more than once....

              }
              if (drflag == 'E')        // this ion has no parameters, so it must be the first time through
              {

                n1 = ion[n].nxdrecomb;  //     Get the pointer to the correct bit of the recombination coefficient array. This should already be set from the first time through
                for (n2 = 0; n2 < nparam; n2++)
                {
                  drecomb[n1].e[n2] = drp[n2];  //we are getting e parameters
                }


              }
              else if (drflag == 'C')   //                  must be the second time though, so no need to read in all the other things
              {
                n1 = ion[n].nxdrecomb;  //     Get the pointer to the correct bit of the recombination coefficient array. This should already be set from the first time through
                for (n2 = 0; n2 < nparam; n2++)
                {
                  drecomb[n1].c[n2] = drp[n2];  //           we are getting e parameters
                }
              }

            }                   //close if statement that selects appropriate ion to add data to
          }                     //close loop over ions


          break;

/** @section Dielectronic Recombination - type 2
 * This section reads in type 2 dielectronic recombination rates.
 * The data comes from Chianti, and there are two types of dielectronic data.
 * This is type 2 - which has a fit decribed by Shull & van Steenberg (1982)
 * There are four fitting parameters.
 * A typical ion data would look like this
 *
 *
 @verbatim
 * DR_SHULL 17 2 2.6510e-03 2.9300e-02 2.4120e+05 4.9410e+05
 @endverbatim
 *
 * where there are two lines for a given ion - this is for z=17 (istate=2)
 *
 */



        case 'S':
          nparam = sscanf (aline, "%*s %d %d %le %le %le %le ", &z, &ne, &drp[0], &drp[1], &drp[2], &drp[3]);   //split and assign the line
          nparam -= 2;          //take 4 off the nparam to give the number of actual parameters
          if (nparam > 4 || nparam < 1) //     trap errors - not as robust as usual because there are a varaible number of parameters...
          {
            Error ("Something wrong with dielectronic recombination data\n", file, lineno);
            Error ("Get_atomic_data: %s\n", aline);
            exit (0);
          }

          istate = ne;          //         get the ionisation state we are recombining from

          for (n = 0; n < nions; n++)   //Loop over ions to find the correct place to put the data
          {
            if (ion[n].z == z && ion[n].istate == istate)       // this works out which ion we are dealing with
            {
              if (ion[n].drflag == 0)   //This is the first time we have dealt with this ion
              {
                drecomb[ndrecomb].nion = n;     //put the ion number into the DR structure
                drecomb[ndrecomb].nparam = nparam;      //Put the number of parameters we ware going to read in, into the DR structure so we know what to iterate over later
                ion[n].nxdrecomb = ndrecomb;    //put the number of the DR into the ion
                drecomb[ndrecomb].type = DRTYPE_SHULL;  //define the type of data
                ndrecomb++;     //increment the counter of number of dielectronic recombination parameter sets
                ion[n].drflag++;        //increment the flag by 1. We will do this rather than simply setting it to 1 so we will get errors if we do this more than once....

              }
              n1 = ion[n].nxdrecomb;    //     Get the pointer to the correct bit of the recombination coefficient array. This should already be set from the first time through
              for (n2 = 0; n2 < nparam; n2++)
              {
                drecomb[n1].shull[n2] = drp[n2];        //we are getting e parameters
              }
            }
          }
          break;

/**
 * @section total radiative Recombination rates from Chianti - type 1 and 2
 * !RR RATE COEFFICIENT FITS (C)20110412 N. R. BADNELL, DEPARTMENT OF PHYSICS, UNIVERSITY OF STRATHCLYDE, GLASGOW G4 0NG, UK.
 * !This is total radiative rate into all possible states or the recombined ion.
 * !This data was downloaded from http://amdpp.phys.strath.ac.uk/tamoc/DATA/RR/ on 13/July/2013
 * !It is described in http://adsabs.harvard.edu/abs/2006ApJS..167..334B
 * !Only modification to original file is to insert the BAD_T_RR label to the front of each line and comments
 * !There are metastable levels in here, M=1 means we are recombining from the ground state, which is what we tend to want.
 * Some of these have one four parameters, a fit introduced by Verner & Ferland (1996) and
 * some have six, which is from Gu (2003). A four parameter data type is shown below
 *
 * @verbatim
 * !  Z  N  M  W      A        B        T0         T1        C        T2
 * BAD_T_RR  1  0  1  1  8.318E-11  0.7472  2.965E+00  7.001E+05
 * BAD_T_RR  2  0  1  1  1.818E-10  0.7492  1.017E+01  2.786E+06
 * BAD_T_RR  3  0  1  1  2.867E-10  0.7493  2.108E+01  6.268E+06
 * BAD_T_RR  4  0  1  1  3.375E-10  0.7475  4.628E+01  1.121E+07
 * BAD_T_RR  5  0  1  1  4.647E-10  0.7484  6.142E+01  1.753E+07
 * @endverbatim
 * */

        case 'T':              /*Badnell type total raditive rate coefficients read in */

          nparam = sscanf (aline, "%*s %d %d %d %le %le %le %le %le %le", &z, &ne, &w, &btrr[0], &btrr[1], &btrr[2], &btrr[3], &btrr[4], &btrr[5]);     //split and assign the line
          nparam -= 3;          //take 4 off the nparam to give the number of actual parameters
          if (nparam > 6 || nparam < 1) //     trap errors - not as robust as usual because there are a varaible number of parameters...
          {
            Error ("Something wrong with badnell total RR data\n", file, lineno);
            Error ("Get_atomic_data: %s\n", aline);
            exit (0);
          }

          istate = ne;          //         get the traditional ionisation state
          for (n = 0; n < nions; n++)   //Loop over ions to find the correct place to put the data
          {
            if (ion[n].z == z && ion[n].istate == istate)       // this works out which ion we are dealing with
            {
              if (ion[n].total_rrflag == 0)     // this ion has no parameters, so it must be the first time through
              {
                total_rr[n_total_rr].nion = n;  //put the ion number into the bad_t_rr structure
                ion[n].nxtotalrr = n_total_rr;  /*put the number of the bad_t_rr into the ion
                                                   structure so we can go either way. */
                total_rr[n_total_rr].type = RRTYPE_BADNELL;
                for (n1 = 0; n1 < nparam; n1++)
                {
                  total_rr[n_total_rr].params[n1] = btrr[n1];   //we are getting  parameters
                }
                ion[n].total_rrflag++;  //increment the flag by 1. We will do this rather than simply setting it to 1 so we will get errors if we do this more than once....
                n_total_rr++;   //increment the counter of number of dielectronic recombination parameter sets
              }
              else if (ion[n].total_rrflag > 0) //       unexpected second line matching z and charge
              {
                Error ("More than one badnell total RR rate for ion %i\n", n);
                Error ("Get_atomic_data: %s\n", aline);
                exit (0);
              }
              else              //if flag is not a positive number, we have a problem
              {
                Error ("Total radiative recombination flag giving odd results\n");
                exit (0);
              }
            }                   //close if statement that selects appropriate ion to add data to
          }                     //close loop over ions



          break;

/**
 * @section Total radiative recombination - type 3
 * These are some rare ion recombination rates that have
 * a two parameter fit by Aldrovandi & Pequignot (1973)
 * They have data like
 * @verbatim
 * RR_SHULL 22 2 1.0965e-11 6.9901e-01
 * RR_SHULL 22 3 1.8600e-11 7.2800e-01
 * @endverbatim
 *
 * For erroneous historical reasons they are referred to as SHULL
 *
 */




        case 's':
          nparam = sscanf (aline, "%*s %d %d %le %le ", &z, &ne, &btrr[0], &btrr[1]);   //split and assign the line
          nparam -= 2;          //take 4 off the nparam to give the number of actual parameters
          if (nparam > 6 || nparam < 1) //     trap errors - not as robust as usual because there are a varaible number of parameters...
          {
            Error ("Something wrong with shull total RR data\n", file, lineno);
            Error ("Get_atomic_data: %s\n", aline);
            exit (0);
          }

          istate = ne;          //         get the traditional ionisation state
          for (n = 0; n < nions; n++)   //Loop over ions to find the correct place to put the data
          {
            if (ion[n].z == z && ion[n].istate == istate)       // this works out which ion we are dealing with
            {
              if (ion[n].total_rrflag == 0)     // this ion has no parameters, so it must be the first time through
              {
                total_rr[n_total_rr].nion = n;  //put the ion number into the bad_t_rr structure
                ion[n].nxtotalrr = n_total_rr;  /*put the number of the bad_t_rr into the ion
                                                   structure so we can go either way. */
                total_rr[n_total_rr].type = RRTYPE_SHULL;
                for (n1 = 0; n1 < nparam; n1++)
                {
                  total_rr[n_total_rr].params[n1] = btrr[n1];   //we are getting  parameters
                }
                ion[n].total_rrflag++;  //increment the flag by 1. We will do this rather than simply setting it to 1 so we will get errors if we do this more than once....
                n_total_rr++;   //increment the counter of number of dielectronic recombination parameter sets
              }
              else if (ion[n].total_rrflag > 0) //       unexpected second line matching z and charge
              {
                Error ("More than one total RR rate for ion %i\n", n);
                Error ("Get_atomic_data: %s\n", aline);
                exit (0);
              }
              else              //if flag is not a positive number, we have a problem
              {
                Error ("Total radiative recombination flag giving odd results\n");
                exit (0);
              }
            }                   //close if statement that selects appropriate ion to add data to
          }                     //close loop over ions
          break;


/**
 * @section Ground state recombination rate data
 * !RR RATE COEFFICIENT FITS (C)20110412 N. R. BADNELL, DEPARTMENT OF PHYSICS, UNIVERSITY OF STRATHCLYDE, GLASGOW G4 0NG, UK.
 * !This is resolved radiative rate from ground state to ground state.
 * !This data was downloaded from http://amdpp.phys.strath.ac.uk/tamoc/DATA/RR/ on 13/July/2013
 * !It is described in http://adsabs.harvard.edu/abs/2006ApJS..167..334B
 * !THe data comes as one file for each ion, with many rates. This file contains the first rate line in each file
 * ! which should relate to the ground state of the upper ion recombining into the ground state of the recombined ion.
 * !The data is tabulated in temperature, first line for each ion is temp, second is rate.
 * !As with other Badnell data, first number is z, second is number of remianing electrons.
 *
 * @verbatim
 * BAD_GS_RR T 1 0 1.00E+01 2.00E+01 5.00E+01 1.00E+02 2.00E+02 5.00E+02 1.00E+03 2.00E+03 5.00E+03 1.00E+04 2.00E+04 5.00E+04 1.00E+05 2.00E+05 5.00E+05 1.00E+06 2.00E+06 5.00E+06 1.00E+07
 * BAD_GS_RR R 1 0 5.21E-12 3.68E-12 2.33E-12 1.65E-12 1.16E-12 7.35E-13 5.18E-13 3.65E-13 2.28E-13 1.58E-13 1.08E-13 6.21E-14 3.88E-14 2.28E-14 1.01E-14 5.05E-15 2.36E-15 7.90E-16 3.27E-16
 * @endverbatim

 */

        case 'G':
          nparam = sscanf (aline, "%*s %s %d %d %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le", &gsflag, &z, &ne, &gstemp[0], &gstemp[1], &gstemp[2], &gstemp[3], &gstemp[4], &gstemp[5], &gstemp[6], &gstemp[7], &gstemp[8], &gstemp[9], &gstemp[10], &gstemp[11], &gstemp[12], &gstemp[13], &gstemp[14], &gstemp[15], &gstemp[16], &gstemp[17], &gstemp[18]);   //split and assign the line
          nparam -= 3;          //take 4 off the nparam to give the number of actual parameters
          if (nparam > 19 || nparam < 1)        //     trap errors - not as robust as usual because there are a varaible number of parameters...
          {
            Error ("Something wrong with badnell GS RR data\n");
            Error ("Get_atomic_data: %s\n", aline);
            exit (0);
          }
          istate = z - ne + 1;  //         get the traditional ionisation state
          for (n = 0; n < nions; n++)   //Loop over ions to find the correct place to put the data
          {
            if (ion[n].z == z && ion[n].istate == istate)       // this works out which ion we are dealing with
            {
              if (ion[n].bad_gs_rr_t_flag == 0 && ion[n].bad_gs_rr_r_flag == 0) //This is first set of this type of data for this ion
              {
                bad_gs_rr[n_bad_gs_rr].nion = n;        //put the ion number into the bad_t_rr structure
                ion[n].nxbadgsrr = n_bad_gs_rr; //put the number of the bad_t_rr into the ion structure so we can go either way.
                n_bad_gs_rr++;  //increment the counter of number of ground state RR
              }
              /*Now work out what type of line it is, and where it needs to go */
              if (gsflag == 'T')        //it is a temperature line
              {
                if (ion[n].bad_gs_rr_t_flag == 0)       //and we need a temp line for this ion
                {
                  if (gstemp[0] > gstmin)
                    gstmin = gstemp[0];
                  if (gstemp[18] < gstmax)
                    gstmax = gstemp[18];
                  ion[n].bad_gs_rr_t_flag = 1;  //set the flag
                  for (n1 = 0; n1 < nparam; n1++)
                  {
                    bad_gs_rr[ion[n].nxbadgsrr].temps[n1] = gstemp[n1];
                  }
                }
                else if (ion[n].bad_gs_rr_t_flag == 1)  //we already have a temp line for this ion
                {
                  Error ("More than one temp line for badnell GS RR rate for ion %i\n", n);
                  Error ("Get_atomic_data: %s\n", aline);
                  exit (0);
                }
                else            //some other odd thing had happened
                {
                  Error ("Get_atomic_data: %s\n", aline);
                  exit (0);
                }
              }
              else if (gsflag == 'R')   //it is a rate line
              {
                if (ion[n].bad_gs_rr_r_flag == 0)       //and we need a rate line for this ion
                {
                  ion[n].bad_gs_rr_r_flag = 1;  //set the flag
                  for (n1 = 0; n1 < nparam; n1++)
                  {
                    bad_gs_rr[ion[n].nxbadgsrr].rates[n1] = gstemp[n1];
                  }
                }
                else if (ion[n].bad_gs_rr_r_flag == 1)  //we already have a rate line for this ion
                {
                  Error ("More than one rate line for badnell GS RR rate for ion %i\n", n);
                  Error ("Get_atomic_data: %s\n", aline);
                  exit (0);
                }
                else            //some other odd thing had happened
                {
                  Error ("Get_atomic_data: %s\n", aline);
                  exit (0);
                }
              }
              else              //We have some problem with this line
              {
                Error ("Get_atomic_data: %s\n", aline);
                exit (0);
              }
            }                   //end of loop over dealing with data for a discovered ion
          }                     //end of loop over ions

          break;

/**
 * @section gaunt factor
 * The following are lines to read in temperature averaged free-free gaunt factors
 * from the data of Sutherland (1998). The atomic file is basically unchanged
 * from the data on the website, just with the top few lines commented out,
 * and a label prepended to each line
 * The data is a spline fit to the gaunt factor as a function of g - the reduced
 * temperature.  The file format is shown below; there are 45 lines
 * @verbatim
 * # log(g^2) <gff(g^2)>       s1           s2            s3
 * # ------------------------------------------------------------
 * FF_GAUNT -4.00 1.113883E+00  1.348000E-02  1.047100E-02 -1.854972E-03
 * FF_GAUNT -3.80 1.116983E+00  1.744580E-02  9.358014E-03  5.564917E-03
 * FF_GAUNT -3.60 1.120891E+00  2.185680E-02  1.269696E-02  4.970250E-03
 * FF_GAUNT -3.40 1.125810E+00  2.753201E-02  1.567911E-02  6.429140E-03
* @endverbatim

 */
        case 'g':
          nparam = sscanf (aline, "%*s %le %le %le %le %le", &gsqrdtemp, &gfftemp, &s1temp, &s2temp, &s3temp);  //split and assign the line
          if (nparam > 5 || nparam < 1) //     trap errors
          {
            Error ("Something wrong with sutherland gaunt data\n");
            Error ("Get_atomic_data: %s\n", aline);
            exit (0);
          }
          if (gaunt_n_gsqrd == 0 || gsqrdtemp > gaunt_total[gaunt_n_gsqrd - 1].log_gsqrd)       //We will use it if it's our first piece of data or is in order
          {
            gaunt_total[gaunt_n_gsqrd].log_gsqrd = gsqrdtemp;   //The scaled electron temperature squared for this array
            gaunt_total[gaunt_n_gsqrd].gff = gfftemp;
            gaunt_total[gaunt_n_gsqrd].s1 = s1temp;
            gaunt_total[gaunt_n_gsqrd].s2 = s2temp;
            gaunt_total[gaunt_n_gsqrd].s3 = s3temp;
            gaunt_n_gsqrd++;
          }
          else
          {
            Error ("Something wrong with gaunt data\n");
            Error ("Get_atomic_data %s\n", aline);
            exit (0);
          }



          break;
/**
 * @section direct (collisional) ionization data from Dere 07.
 * #Title: Ionization rate coefficients for elements H to Zn (Dere+, 2007)
 * #Table  J_A_A_466_771_table29:
 * #Title: Spline fits, multiplied by a factor of 10^6^, to the scaled ionization rate coefficients
 * #Column Z       (I2)    [1/30] Nuclear charge   [ucd=phys.atmol.element]
 * #Column Ion     (I2)    [1/30] Ion (spectroscopic notation)     [ucd=phys.atmol.ionStage]
 * #Column Nspl    (I2)    [15/20] Number of splines       [ucd=meta.number]
 * #Column I       (F9.3)  Ionization potential    [ucd=phys.energy;phys.atmol.ionization]
 * #Column Tmin    (E10.3) Minimum temperature     [ucd=phys.temperature]
 * #Column x1-X20      (F7.4)  Scaled temperature 1 (1)        [ucd=phys.temperature]
 * #Column rho1 -rho20   (F8.4)  ? Scaled rate coefficient 1 (2) [ucd=arith.rate;phys.atmol.collisional]
 */
        case 'd':
          nparam = sscanf (aline, "%*s %d %d %d %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le", &z, &istate, &nspline, &et, &tmin, &temp[0], &temp[1], &temp[2], &temp[3], &temp[4], &temp[5], &temp[6], &temp[7], &temp[8], &temp[9], &temp[10], &temp[11], &temp[12], &temp[13], &temp[14], &temp[15], &temp[16], &temp[17], &temp[18], &temp[19], &temp[20], &temp[21], &temp[22], &temp[23], &temp[24], &temp[25], &temp[26], &temp[27], &temp[28], &temp[29], &temp[30], &temp[31], &temp[32], &temp[33], &temp[34], &temp[35], &temp[36], &temp[37], &temp[38], &temp[39]);     //split and assign the line

          if (nparam != 5 + (nspline * 2))      //     trap errors
          {
            Error ("Something wrong with Dere DI data\n");
            Error ("Get_atomic_data: %s\n", aline);
            exit (0);
          }
          for (n = 0; n < nions; n++)   //Loop over ions to find the correct place to put the data
          {
            if (ion[n].z == z && ion[n].istate == istate)       // this works out which ion we are dealing with
            {
              if (ion[n].dere_di_flag == 0)     //This is first set of this type of data for this ion
              {
                ion[n].dere_di_flag = 1;
                dere_di_rate[n_dere_di_rate].nion = n;  //put the ion number into the dere_di_rate structure
                ion[n].nxderedi = n_dere_di_rate;       //put the number of the dere_di_rate into the ion structure so we can go either way.
                dere_di_rate[n_dere_di_rate].xi = et;
                dere_di_rate[n_dere_di_rate].min_temp = tmin;
                dere_di_rate[n_dere_di_rate].nspline = nspline;
                for (n1 = 0; n1 < nspline; n1++)
                {
                  dere_di_rate[n_dere_di_rate].temps[n1] = temp[n1];
                  dere_di_rate[n_dere_di_rate].rates[n1] = temp[n1 + nspline] * 1e-6;

                }
                n_dere_di_rate++;       //increment the counter of number of ground state RR
              }
              else
              {
                Error ("Get_atomic_data: More than one Dere DI rate for ion %i\n", n);
              }
            }
          }
          break;

/**
 * @section Electron yield - goes with auger ionization rates
 * #This is electron yield data from Kaastra & Mewe (1993) -1993A&AS...97..443K
 * #It is processed from data downloaded from http://vizier.cfa.harvard.edu/viz-bin/VizieR?-source=J/A+AS/97/443
 * #Via the python script kaastra_2_py.py. Example lines of file below
 * #
 * @verbatim
 * #Label z state n l IP mean_electron_energy Prob_of_1e Prob_of_2e Prob_of_3e ....
 * Kelecyield 4 1 1 0 1.1500e+02 9.280e+01 0.000e+00 1.000e+04 0.000e+00 0.000e+00 0.000e+00 0.000e+00 0.000e+00 0.000e+00 0.000e+00 0.000e+00
 * Kelecyield 5 1 1 0 1.9200e+02 1.639e+02 6.000e+00 9.994e+03 0.000e+00 0.000e+00 0.000e+00 0.000e+00 0.000e+00 0.000e+00 0.000e+00 0.000e+00
 * @endverbatim

*/

        case 'K':
          nparam =
            sscanf (aline,
                    "%*s %d %d %d %d %le %le %le %le %le %le %le %le %le %le %le %le",
                    &z, &istate, &in, &il, &I, &Ea, &temp[0],
                    &temp[1], &temp[2], &temp[3], &temp[4], &temp[5], &temp[6], &temp[7], &temp[8], &temp[9]);
          if (nparam != 16)
          {
            Error ("Something wrong with electron yield data\n");
            Error ("Get_atomic_data %s\n", aline);
            exit (0);
          }
          for (n = 0; n < n_inner_tot; n++)
          {
            if (inner_cross[n].z == z && inner_cross[n].istate == istate && inner_cross[n].n == in && inner_cross[n].l == il)
            {
              if (inner_cross[n].n_elec_yield == -1)    /*This is the first yield data for this vacancy */
              {
                inner_elec_yield[n_elec_yield_tot].nion = n;    /*This yield refers to this ion */
                inner_cross[n].n_elec_yield = n_elec_yield_tot;
                inner_elec_yield[n_elec_yield_tot].z = z;
                inner_elec_yield[n_elec_yield_tot].istate = istate;
                inner_elec_yield[n_elec_yield_tot].n = in;
                inner_elec_yield[n_elec_yield_tot].l = il;
                inner_elec_yield[n_elec_yield_tot].I = I * EV2ERGS;
                inner_elec_yield[n_elec_yield_tot].Ea = Ea * EV2ERGS;
                for (n1 = 0; n1 < 10; n1++)
                {
                  inner_elec_yield[n_elec_yield_tot].prob[n1] = temp[n1] / 10000.0;
                }
                n_elec_yield_tot++;
              }
              else
              {
                Error ("Get_atomic_data: more than one electron yield record for inner_cross %i z=%i istate=%i\n", n, z, istate);
              }
            }
          }
          break;
/**
 * @section Fluorescent photon yield from inner shell ionization - not currently used but read in.
		  now not read in - see #499
 * #This is fluorescent yield data from Kaastra & Mewe (1993) -1993A&AS...97..443K
 * #It is processed from data downloaded from http://vizier.cfa.harvard.edu/viz-bin/VizieR?-source=J/A+AS/97/443
 * #Via the python script kaastra_2_py.py.
 * Data format below
 * #
 * @verbatim
 * #Label z state n l photon_energy yield
 * Kphotyield 5 1 1 0 1.837e+02 6.000e-04
 * Kphotyield 5 1 1 0 1.690e+01 7.129e-01
 * @endverbatim
 
        case 'F':
          nparam = sscanf (aline, "%*s %d %d %d %d %le %le ", &z, &istate, &in, &il, &energy, &yield);
          if (nparam != 6)
          {
            Error ("Something wrong with fluorescent yield data\n");
            Error ("Get_atomic_data %s\n", aline);
            exit (0);
          }
          for (n = 0; n < n_inner_tot; n++)
          {
            if (inner_cross[n].z == z && inner_cross[n].istate == istate && inner_cross[n].n == in && inner_cross[n].l == il)
            {
              if (inner_cross[n].n_fluor_yield == -1)   //This is the first yield data for this vacancy 
              {
                inner_fluor_yield[n_fluor_yield_tot].nion = n;  //This yield refers to this ion 
                inner_cross[n].n_fluor_yield = n_fluor_yield_tot;
                inner_fluor_yield[n_fluor_yield_tot].z = z;
                inner_fluor_yield[n_fluor_yield_tot].istate = istate;
                inner_fluor_yield[n_fluor_yield_tot].n = in;
                inner_fluor_yield[n_fluor_yield_tot].l = il;
                inner_fluor_yield[n_fluor_yield_tot].freq = energy / HEV;
                inner_fluor_yield[n_fluor_yield_tot].yield = yield;
                n_fluor_yield_tot++;
              }
              else
              {
                Error ("Get_atomic_data: more than one fluorescent yield record for inner_cross %i\n", n);
              }
            }
          }
          break;
		  */

/**
 * @section Chianti Collision Strengths
 * The lines below read in collision strength data from Chianti (after Burgess and Tully). The original
 *		  is stored in .scups files in chianti. The python script searches for matches to the lines_linked_ver_2.py
 *		  data file, on the basis of energy and oscillator strength (and z and state). As a rather hamfisted
 *		  approach, all of the line data is stored on the same line as the collision strength data in
 *		  the file, so it is possible to make a very accurate match. It does mean quite a lot is read in
 *		  from a line, but most is thrown away. Currently the code below matches based on z, state, upper and
 *		  lower level numbers and oscillator strength
 * Each line is defned with a set of three lines
 * first line is CSTREN - this contains line data to try and match
 * Second two lines are collision strength data from Burgess and Tully 1992A&A...254..436B
 * These are generally 5 pojnt but up to 20 point spline fits of the T vs upsilon data. Typical lines are below
 * @verbatim
CSTREN Line  1  1 1215.673584  0.139000   2   2     0.000000    10.200121    0    1       1      3   7.500e-01   2.772e-01   1.478e+00    5    1   1.700e+00
SCT   0.000e+00   2.500e-01   5.000e-01   7.500e-01   1.000e+00
SCUPS    1.132e-01   2.708e-01   5.017e-01   8.519e-01   1.478e+00
 * @endverbatim



 *		  */

        case 'C':
          nparam =
            (sscanf
             (aline,
              "%*s %*s %d %2d %le %le %le %le %le %le %d %d %d %d %le %le %le %d %d %le",
              &z, &istate, &freq, &f, &gl, &gu, &el, &eu, &levl, &levu, &c_l, &c_u, &en, &gf, &hlt, &np, &type, &sp));
          if (nparam != 18)
          {
            Error ("Get_atomic_data: file %s line %d: Collision strength line incorrectly formatted\n", file, lineno);
            Error ("Get_atomic_data: %s\n", aline);
            exit (0);
          }
          match = 0;
          for (n = 0; n < nlines; n++)  //loop over all the lines we have read in - look for a match
          {
            if (line[n].z == z && line[n].istate == istate
                && line[n].levl == levl && line[n].levu == levu && line[n].gl == gl && line[n].gu == gu && line[n].f == f)
            {
              if (line[n].coll_index > -1)      //We already have a collision strength record from this line - throw an error and quit
              {
                Error ("Get_atomic_data More than one collision strength record for line %i\n", n);
                exit (0);
              }
              match = 1;
              coll_stren[n_coll_stren].n = n_coll_stren;
              coll_stren[n_coll_stren].lower = c_l;
              coll_stren[n_coll_stren].upper = c_u;
              coll_stren[n_coll_stren].energy = en;
              coll_stren[n_coll_stren].gf = gf;
              coll_stren[n_coll_stren].hi_t_lim = hlt;
              coll_stren[n_coll_stren].n_points = np;
              coll_stren[n_coll_stren].type = type;
              coll_stren[n_coll_stren].scaling_param = sp;

              line[n].coll_index = n_coll_stren;        //point the line to its matching collision strength

              //We now read in two lines of fitting data
              if (fgets (aline, LINELENGTH, fptr) == NULL)
              {
                Error ("Get_atomic_data: Problem reading collision strength record\n");
                Error ("Get_atomic_data: %s\n", aline);
                exit (0);
              }

              /* JM 1709 -- increased number of entries read up to max of 20 */
              nparam =
                sscanf (aline,
                        "%*s %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le",
                        &temp[0], &temp[1], &temp[2], &temp[3],
                        &temp[4], &temp[5], &temp[6], &temp[7],
                        &temp[8], &temp[9], &temp[10], &temp[11],
                        &temp[12], &temp[13], &temp[14], &temp[15], &temp[16], &temp[17], &temp[18], &temp[19]);

              for (nn = 0; nn < np; nn++)
              {
                coll_stren[n_coll_stren].sct[nn] = temp[nn];
              }
              if (fgets (aline, LINELENGTH, fptr) == NULL)
              {
                Error ("Get_atomic_data: Problem reading collision strength record\n");
                Error ("Get_atomic_data: %s\n", aline);
                exit (0);
              }

              nparam =
                sscanf (aline,
                        "%*s %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le",
                        &temp[0], &temp[1], &temp[2], &temp[3],
                        &temp[4], &temp[5], &temp[6], &temp[7],
                        &temp[8], &temp[9], &temp[10], &temp[11],
                        &temp[12], &temp[13], &temp[14], &temp[15], &temp[16], &temp[17], &temp[18], &temp[19]);

              for (nn = 0; nn < np; nn++)
              {
                coll_stren[n_coll_stren].scups[nn] = temp[nn];
              }
              n_coll_stren++;
            }
          }
          if (match == 0)       //Fix for an error where a line match isn't found - this then causes the next two lines to be skipped
          {
            fgets (aline, LINELENGTH, fptr);
            fgets (aline, LINELENGTH, fptr);
            cstren_no_line++;
          }
          break;

        case 'c':              /* It was a comment line so do nothing */
          break;
        case 'z':
        default:
          Error ("get_atomicdata: Could not interpret line %d in file %s: %s\n", lineno, file, aline);
          break;
        }

        strcpy (aline, "");
      }

      fclose (fptr);
    }
    /*End of do loop for reading a particular file of data */
  }

/* End of main do loop for reading all of the the data. At this point we can close
   the masterfile
 */

  fclose (mptr);
/* OK now summarize the data that has been read*/

  n_elec_yield_tot = 0;         //Reset this numnber, we are now going to use it to check we have yields for all inner shells
//  n_fluor_yield_tot = 0;        //Reset this numnber, we are now going to use it to check we have yields for all inner shells
  inner_no_e_yield = 0;

  for (n = 0; n < n_inner_tot; n++)
  {
    if (inner_cross[n].n_elec_yield != -1)
      n_elec_yield_tot++;
    else
      inner_no_e_yield++;
//      Error_silent ("get_atomicdata: No inner electron yield data for inner cross section %i\n", n);
//    if (inner_cross[n].n_fluor_yield != -1)
//    n_fluor_yield_tot++;

  }

  Log ("Data of %3d elements, %3d ions, %5d levels, %5d lines, and %5d topbase records\n", nelements, nions, nlevels, nlines, ntop_phot);
  Log
    ("Macro   %3d elements, %3d ions, %5d levels, %5d lines, and %5d topbase records\n",
     nelements, nions_macro, nlevels_macro, nlines_macro, ntop_phot_macro);
  Log
    ("Simple  %3d elements, %3d ions, %5d levels, %5d lines, and %5d topbase records\n",
     nelements, nions_simple, nlevels_simple, nlines_simple, ntop_phot_simple);
  Log ("We have read in %3d photoionization cross sections\n", nphot_total);
  Log ("                %3d are topbase \n", ntop_phot);
  Log ("                %3d are VFKY \n", nxphot);
  Log ("We have read in %5d   Chiantic collision strengths\n", n_coll_stren);   //1701 nsh collision strengths
  Log ("We have read in %3d Inner shell photoionization cross sections\n", n_inner_tot);        //110818 nsh added a reporting line about dielectronic recombination coefficients
  Log ("                %3d have matching electron yield data\n", n_elec_yield_tot);
//  Log ("                %3d have matching fluorescent yield data\n", n_fluor_yield_tot);

  Log ("We have read in %3d Dielectronic recombination coefficients\n", ndrecomb);      //110818 nsh added a reporting line about dielectronic recombination coefficients
  Log ("We have read in %3d Badnell totl Radiative rate coefficients\n", n_total_rr);
  Log ("We have read in %3d Badnell GS   Radiative rate coefficients over the temp range %e to %e\n", n_bad_gs_rr, gstmin, gstmax);
  Log ("We have read in %3d Scaled electron temperature frequency averaged gaunt factors\n", gaunt_n_gsqrd);
  Log ("The minimum frequency for photoionization is %8.2e\n", phot_freq_min);
  Log ("The minimum frequency for inner shell ionization is %8.2e\n", inner_freq_min);

  /* report ignored simple lines for macro-ions */
  for (n = 0; n < NIONS; n++)
  {
    if (simple_line_ignore[n] > 0)
      Error ("Ignored %d simple lines for macro-ion %d\n", simple_line_ignore[n], n);
  }
  /* report ignored collision strengths */
  if (cstren_no_line > 0)
    Error ("Ignored %d collision strengths with no matching line transition\n", cstren_no_line);
  if (inner_no_e_yield > 0)
    Error ("Ignoring %d inner shell cross sections because no matching yields\n", inner_no_e_yield);



/* Now begin a series of calculations with the data that has been read in in order
to prepare it for use by other programs*/

/* Calculate the conversion factor between rho and nh.  Note that our approximation here is
somewhat inexact  because it assumes all elements have equal neutrons and protons. */

  q = 1.;
  for (nelem = 1; nelem < nelements; nelem++)
  {
    q += ele[nelem].abun * 2. * ele[nelem].z;
  }

  q = MPROT * q;
  rho2nh = 1. / q;


/* Now find the first and last ion associated with each named element and
exit if there is an element with no ions */

  for (nelem = 0; nelem < nelements; nelem++)
  {
    n = 0;
    while (ion[n].z != ele[nelem].z && n < nions)
      n++;                      /* Find the first ion of that element for which there is data */

    ele[nelem].firstion = n;
    ion[n].nelem = nelem;

    /* find the highest ion stage and the number of ions */
    ele[nelem].istate_max = ion[n].istate;

    while (ion[n].z == ele[nelem].z && n < nions)
    {
      if (ele[nelem].istate_max < ion[n].istate)
        ele[nelem].istate_max = ion[n].istate;
      ion[n].nelem = nelem;
      n++;
    }
    ele[nelem].nions = n - ele[nelem].firstion;
    if (ele[nelem].firstion == nions)
    {
      ele[nelem].firstion = -1; /* There were no ions for this element */
// In principle, there might be a program which uses elements but not ions, but it seems unlikely,
// therefore stop if there is not at least one ion for each element
      Error ("Get_atomic_data: There were no ions for element %d %s\n", nelem, ele[nelem].name);
      exit (0);
    }
  }


/* Now attempt to associate lines with levels. If an association is found use, create
a total emission oscillator strength for the level....really ought to be radiative lifefime */


/* This next loop is connects levels to configurations for "simple atoms".
 For Macro Atoms files, the data files already contain this infomration and
 so the check is avoided by checking tte macro_info flag SS
*/



  for (n = 0; n < nlines; n++)
  {
    if (ion[line[n].nion].macro_info == 0)      // not a macro atom (SS)
    {
      mstart = ion[line[n].nion].firstlevel;
      mstop = mstart + ion[line[n].nion].nlevels;

      m = mstart;
      while (config[m].ilv != line[n].levl && m < mstop)
        m++;
      if (m < mstop)
        line[n].nconfigl = m;
      else
        line[n].nconfigl = -9999;

      m = mstart;
      while (config[m].ilv != line[n].levu && m < mstop)
        m++;
      if (m < mstop)
      {
        line[n].nconfigu = m;
        config[m].rad_rate += a21 (&line[n]);
      }
      else
        line[n].nconfigu = -9999;
    }
  }

/* Check that all of the macro_info variables are initialized to 1
or zero so that simple checks of true and false can be used for them */

  for (n = 0; n < nions; n++)
  {
    if (ion[n].macro_info == -1)
    {
      Error ("Ion %d for element %s and ion %d is of unknown type\n", n, ion[n].z, ion[n].istate);
      exit (0);
    }
  }

  for (n = 0; n < nlevels; n++)
  {
    if (config[n].macro_info == -1)
    {
      Error ("Level %d for element %s and ion %d is of unknown type\n", n, config[n].z, config[n].istate);
      exit (0);
    }
  }

  for (n = 0; n < nlines; n++)
  {
    if (line[n].macro_info == -1)
    {
      Error ("Level %d for element %s and ion %d is of unknown type\n", n, line[n].z, line[n].istate);
      exit (0);
    }
  }

  for (n = 0; n < nphot_total; n++)
  {
    if (phot_top[n].macro_info == -1)
    {
      Error ("Photoionization cross-section %d for element %s and ion %d is of unknown type\n", n, phot_top[n].z, phot_top[n].istate);
      exit (0);
    }
  }

/* Finally evaluate how close we are to limits set in the structures */

  Log ("get_atomic_data: Evaluation:  There are %6d elements     while %6d are currently allowed\n", nelements, NELEMENTS);
  Log ("get_atomic_data: Evaluation:  There are %6d ions         while %6d are currently allowed\n", nions, NIONS);
  Log ("get_atomic_data: Evaluation:  There are %6d levels       while %6d are currently allowed\n", nlevels, NLEVELS);
  Log ("get_atomic_data: Evaluation:  There are %6d lines        while %6d are currently allowed\n", nlines, NLINES);
  Log ("get_atomic_data: Evaluation:  There are %6d macro levels while %6d are currently allowed\n", nlevels_macro, NLEVELS_MACRO);

  bb_max = 0;
  bf_max = 0;

  for (i = 0; i < nlevels; i++)
  {
    if (bb_max < config[i].n_bbu_jump)
      bb_max = config[i].n_bbu_jump;
    if (bb_max < config[i].n_bbd_jump)
      bb_max = config[i].n_bbd_jump;
    if (bf_max < config[i].n_bfu_jump)
      bf_max = config[i].n_bfu_jump;
    if (bf_max < config[i].n_bfd_jump)
      bf_max = config[i].n_bfd_jump;
  }


  Log ("get_atomic_data: Evaluation:  The maximum value bb jumps is %d while %d are currently allowed\n", bb_max, NBBJUMPS);
  Log ("get_atomic_data: Evaluation:  The maximum value bf jumps is %d while %d are currently allowed\n", bb_max, NBBJUMPS);

/* Now, write the data to a file so you can check it later if you wish */
/* this is controlled by one of the -d flag modes, defined in atomic.h */
  if (write_atomicdata)
  {

    if ((fptr = fopen ("data.out", "w")) == NULL)
    {
      Error ("get_atomic data:  Could not open data.out\n");
      exit (0);
    }

    fprintf (fptr, "This file contains data which was read in by get_atomicdata\n");

    fprintf (fptr, "Data of %d elements, %d ions, and %d levels %d lines\n", nelements, nions, nlevels, nlines);


    /* Write the element array */
    fprintf (fptr, "Element Data:\n");
    for (nelem = 0; nelem < nelements; nelem++)
    {
      fprintf (fptr, "Element %2d %5s firstion %2d nions %2d\n", nelem, ele[nelem].name, ele[nelem].firstion, ele[nelem].nions);
    }

    /* Write the ion array */
    fprintf (fptr, "Ion data:\n");
    for (n = 0; n < nions; n++)
    {
      fprintf (fptr,
               "ion %3d z %3d istate %3d firstlevel %3d nlevels %3d potential %8.3g\n",
               n, ion[n].z, ion[n].istate, ion[n].firstlevel, ion[n].nlevels, ion[n].ip / EV2ERGS);
    }

    /* Write the excitation level data */
    fprintf (fptr, "Excitation levels: There are %d levels\n", nlevels);
    for (n = 0; n < nlevels; n++)
      fprintf (fptr, "n %3d q %.1f g %3.0f ex %8.3g\n", n, config[n].q_num, config[n].g, config[n].ex);

    /* Write the photoionization data  */
    fprintf (fptr, "Photoionization data: There are %d edges\n", ntop_phot + nxphot);
    for (n = 0; n < ntop_phot + nxphot; n++)
    {
      fprintf (fptr, "n %3d z %2d istate %3d sigma %8.2e freq[0] %8.2e\n",
               n, phot_top[n].z, phot_top[n].istate, phot_top[n].sigma, phot_top[n].freq[0]);
    }

    /* Write the resonance line data to the file */

    fprintf (fptr, "Line data: There are %d lines\n", nlines);

    for (n = 0; n < nlines; n++)
    {
      fprintf (fptr, "n %3d ion %3d freq %8.1e f %6.3f\n", n, line[n].nion, line[n].freq, line[n].f);
    }

    /* Write the ground fraction data to the file */
    fprintf (fptr, "Ground frac data (just first and last fracs here as a check):\n");

    for (n = 0; n < NIONS; n++)
    {
      fprintf (fptr, "%3d %3d %6.3f %6.3f\n", ground_frac[n].z, ground_frac[n].istate, ground_frac[n].frac[0], ground_frac[n].frac[19]);
    }

    fclose (fptr);
  }                             // end of if statement based on modes.write_atomicdata



  /* Finally create frequency ordered pointers to the various portions
   * of the atomic data
   */

  /* Index the lines */
  index_lines ();

/* Index the topbase photoionization structure by threshold freqeuncy */
  if (ntop_phot + nxphot > 0)
    index_phot_top ();
/* Index the topbase photoionization structure by threshold freqeuncy */
  if (n_inner_tot > 0)
    index_inner_cross ();


  check_xsections ();           // debug routine, only prints if verbosity > 4

  return (0);
}



/**********************************************************/
/**
 * @brief      sort the lines into frequency order
 *
 * @return     Always returns 0
 *
 * @details
 *
 * ### Notes ###
 * All use a Numerical recipes routine indexx for this
 *
 **********************************************************/

int
index_lines ()
{
  float *freqs, foo;
  int *index, ioo;
  int n;
  void indexx ();

  /* Allocate memory for some modestly large arrays */
  freqs = calloc (sizeof (foo), NLINES + 2);
  index = calloc (sizeof (ioo), NLINES + 2);

  freqs[0] = 0;
  for (n = 0; n < nlines; n++)
    freqs[n + 1] = line[n].freq;        /* So filled matrix
                                           elements run from 1 to nlines */

  indexx (nlines, freqs, index);        /* Note that this math recipes routine
                                           expects arrays to run from 1 to nlines inclusive */

  /* The for loop indices are complicated by the numerical recipes routine,
     which is a simple translation of a fortran routine.
     Specifically, index array elements 1 to nlines are now filled,
     and the numbers run from 1 to nlines, but the
     pointer array is only filled from elements 0 to nlines -1 */

  /* SS - adding quantity "where_in_list" to line structure so that it is easy to from emission
     in recombination line to correct place in line list. */


  for (n = 0; n < nlines; n++)
  {
    lin_ptr[n] = &line[index[n + 1] - 1];
    line[index[n + 1] - 1].where_in_list = n;
  }

  /* Free the memory for the arrays */
  free (freqs);
  free (index);

  return (0);
}


/**********************************************************/
/**
 * @brief      Index the topbase photoionzation crossections by frequency
 *
 *
 * @return     Always returns 0
 *
 * @details
 *
 * The results are stored in phot_top_ptr
 *
 * ### Notes ###
 * Adapted from index_lines as part to topbase
 *             addition to python
 *
 **********************************************************/

int
index_phot_top ()
{
  float *freqs, foo;
  int *index, ioo;
  int n;
  void indexx ();

  /* Allocate memory for some modestly large arrays */
  freqs = calloc (sizeof (foo), ntop_phot + nxphot + 2);
  index = calloc (sizeof (ioo), ntop_phot + nxphot + 2);

  freqs[0] = 0;
  for (n = 0; n < ntop_phot + nxphot; n++)
    freqs[n + 1] = phot_top[n].freq[0]; /* So filled matrix
                                           elements run from 1 to ntop_phot */

  indexx (ntop_phot + nxphot, freqs, index);    /* Note that this math recipes routine
                                                   expects arrays to run from 1 to nlines inclusive */

  /* The for loop indices are complicated by the numerical recipes routine,
     which is a simple translation of a fortran routine.
     Specifically, index array elements 1 to nlines are now filled,
     and the numbers run from 1 to nlines, but the
     pointer array is only filled from elements 0 to nlines -1 */

  for (n = 0; n < ntop_phot + nxphot; n++)
  {
    phot_top_ptr[n] = &phot_top[index[n + 1] - 1];
  }

  /* Free the memory for the arrays */
  free (freqs);
  free (index);

  return (0);

}


/**********************************************************/
/**
 * @brief      Index inner shell xsections in frequency order
 *
 * @return     Alwasy returns 0
 *
 * @details
 * The rusults are stored in inner_cross_ptr
 *
 * ### Notes ###
 * ??? NOTES ???
 *
 **********************************************************/

int
index_inner_cross ()
{
  float *freqs, foo;
  int *index, ioo;
  int n;
  void indexx ();

  /* Allocate memory for some modestly large arrays */
  freqs = calloc (sizeof (foo), n_inner_tot + 2);
  index = calloc (sizeof (ioo), n_inner_tot + 2);

  freqs[0] = 0;
  for (n = 0; n < n_inner_tot; n++)
    freqs[n + 1] = inner_cross[n].freq[0];      /* So filled matrix
                                                   elements run from 1 to ntop_phot */

  indexx (n_inner_tot, freqs, index);   /* Note that this math recipes routine
                                           expects arrays to run from 1 to nlines inclusive */

  /* The for loop indices are complicated by the numerical recipes routine,
     which is a simple translation of a fortran routine.
     Specifically, index array elements 1 to nlines are now filled,
     and the numbers run from 1 to nlines, but the
     pointer array is only filled from elements 0 to nlines -1 */

  for (n = 0; n < n_inner_tot; n++)
  {
    inner_cross_ptr[n] = &inner_cross[index[n + 1] - 1];
  }

  /* Free the memory for the arrays */
  free (freqs);
  free (index);

  return (0);

}


/* Numerical recipes routine used by index_lines which in turn is used by get_atomic_data */


/**********************************************************/
/**
 * @brief      Numerical recipes routine used by various routines to sort the atomic data
 *
 * @param [in] int  n   The dimension of the array that needs sorting
 * @param [in] float  arrin[]   The array to sort
 * @param [out] int  indx[]   The array that contains poionters to the sorted array
 * @return     N/A
 *
 * @details
 *
 * ### Notes ###
 * See nuerical recipes
 *
 * @bug Should probably replaces this with the equivalent gsl routine
 *
 **********************************************************/

void
indexx (n, arrin, indx)
     int n, indx[];
     float arrin[];
{
  int l, j, ir, indxt, i;
  float q;
/* NSH 1408 - This routine fails in the very odd circumstance that n=1 so we now do a little test here */
  if (n < 2)
  {
    Log_silent ("Nothing for indexx to do! Only one element\n");
    indx[0] = 0;
    indx[1] = 1;                /* NSH 1707 - We still need to populate the array */
    return;
  }

  for (j = 1; j <= n; j++)
    indx[j] = j;
  l = (n >> 1) + 1;
  ir = n;
  for (;;)
  {
    if (l > 1)
      q = arrin[(indxt = indx[--l])];
    else
    {
      q = arrin[(indxt = indx[ir])];
      indx[ir] = indx[1];
      if (--ir == 1)
      {
        indx[1] = indxt;
        return;
      }
    }
    i = l;
    j = l << 1;
    while (j <= ir)
    {
      if (j < ir && arrin[indx[j]] < arrin[indx[j + 1]])
        j++;
      if (q < arrin[indx[j]])
      {
        indx[i] = indx[j];
        j += (i = j);
      }
      else
        j = ir + 1;
    }
    indx[i] = indxt;
  }
}




/**********************************************************/
/**
 * @brief      uses freqmin and freqmax to set limits on which lines in the
 * frequency ordered list of lines need to be
 * considered when, for example,  calculating which lines have Sobolev surfaces
 * along a path. 
 *
 *
 * @param [in] double  freqmin   The minimum frequency we are interested in
 * @param [in] double  freqmax   The maximum frequency we are interested in
 * @return     the number of lines that are potentially in resonance.
 *
 * limit_lines sets the external variables nline_min and nline_max for use with
 * other routines
 *
 *
 * If limit_lines
 * 	returns 0 there are no lines of interest and one does not need to worry about any
 * 	resonances at this frequency.  If limit_lines returns a number greater than 0, then
 * 	the lines of interest are defined by nline_min and nline_max (inclusive) in atomic.h
 * 	nline_delt is also set which is the number of lines that are in the specified range
 *
 * @details
 * limit_lines  define the lines that are close to a given frequency.  The degree of closeness
 * 	is defined by v. The routine must be used in conjuction with get_atomic_data which
 * 	will have created an ordered list of the lines.
 *
 * ### Notes ###
 * Limit_lines needs to be used somewhat carefully.  Carefully means checking the
 * 	return value of limit_lines or equivalently nline_delt=0.  If nline_delt=0 then there
 * 	were no lines in the region of interest.  Assuming there were lines in thte range,
 * 	one must sum over lines from nline_min to nline_max inclusive.
 *
 * 	One might wonder why nline_max is not set to one larger than the last line which
 * 	is in range.  This is because depending on how the velocity is trending you may
 * 	want to sum from the highest frequency line to the lowest.
 *
 *
 **********************************************************/

int
limit_lines (freqmin, freqmax)
     double freqmin, freqmax;
{

  int nmin, nmax, n;
  double f;


  if (freqmin > lin_ptr[nlines - 1]->freq || freqmax < lin_ptr[0]->freq)
  {
    nline_min = 0;
    nline_max = 0;
    nline_delt = 0;
    return (0);
  }

  f = freqmin;

  nmin = 0;
  nmax = nlines - 1;
  n = (nmin + nmax) >> 1;       // Compute a midpoint >> is a bitwise right shift

  while (n != nmin)
  {
    if (lin_ptr[n]->freq < f)
      nmin = n;
    if (lin_ptr[n]->freq >= f)
      nmax = n;
    n = (nmin + nmax) >> 1;     // Compute a midpoint >> is a bitwise right shift
  }

  nline_min = nmin;

  f = freqmax;
  nmin = 0;
  nmax = nlines - 1;
  n = (nmin + nmax) >> 1;       // Compute a midpoint >> is a bitwise right shift

  while (n != nmin)
  {
    if (lin_ptr[n]->freq <= f)
      nmin = n;
    if (lin_ptr[n]->freq > f)
      nmax = n;
    n = (nmin + nmax) >> 1;     // Compute a midpoint >> is a bitwise right shift
  }

  nline_max = nmax;


  return (nline_delt = nline_max - nline_min + 1);
}




/**********************************************************/
/**
 * @brief      Perform sanity checks on xsection data
 *
 * @return     Alway returns 0; errors are written to the log file
 *
 * @details
 * Results are written to the log file
 *
 * ### Notes ###
 *
 * @bug  A careful look at this routine is warranted as get_atomic_data has changed over the years.
 *
 **********************************************************/

int
check_xsections ()
{
  int nion, n;

  for (n = 0; n < nphot_total; n++)
  {
    nion = phot_top[n].nion;
    if (ion[nion].phot_info == 1)
      Debug
        ("Topbase Ion %i Z %i istate %i nground %i ilv %i ntop %i f0 %8.4e IP %8.4e\n",
         nion, ion[nion].z, ion[nion].istate, ion[nion].ntop_ground, phot_top[n].nlev, ion[nion].ntop, phot_top[n].freq[0], ion[nion].ip);
    else if (ion[nion].phot_info == 0)
      Debug ("Vfky Ion %i Z %i istate %i nground %i f0 %8.4e IP %8.4e\n",
             nion, ion[nion].z, ion[nion].istate, ion[nion].nxphot, phot_top[n].freq[0], ion[nion].ip);

    /* some simple checks -- could be made more robust */
    if (ion[nion].n_lte_max == 0 && ion[nion].phot_info == 1)
    {
      Error
        ("get_atomicdata: not tracking levels for ion %i z %i istate %i, yet marked as topbase xsection!\n",
         nion, ion[nion].z, ion[nion].istate);
    }
    if (ion[nion].phot_info != 1 && ion[nion].macro_info)
    {
      Error
        ("get_atomicdata: macro atom but no topbase xsection! ion %i z %i istate %i, yet marked as topbase xsection!\n",
         nion, ion[nion].z, ion[nion].istate);
    }

    if (ion[nion].macro_info != phot_top[n].macro_info)
    {
      Error ("Macro info for ion doesn't match with xsection!! Could be serious.\n");
    }
  }

  return 0;
}




/*

   q21 calculates and returns the collisional de-excitation coefficient q21

   Notes:  This uses a simple approximation for the effective_collision_strength omega.
   ?? I am not sure where this came from at present and am indeed no sure that
   omega is omega(2->1) as it should be.  This should be carefully checked.??
   This has no been improved, and wherever possible we use tabulated collision
   stength data from Chianti (after Burgess and Tully 1992)

   c21=n_e * q21
   q12 = g_2/g_1 q21 exp(-h nu / kT )

   History:
   98aug        ksl     Recoded from other routines so that one always calculates q21 in
			the same way.
	01nov	ksl	Add tracking mechanism so if called with same conditions it
			returns without recalculation
	12oct	nsh	Added, then commented out approximate gaunt factor given in
			hazy 2.
	17jan	nsh Added code to use the actual collision strength date from chianti

 */


/// (8*PI)/(sqrt(3) *nu_1Rydberg
#define ECS_CONSTANT 4.773691e16

struct lines *q21_line_ptr;
double q21_a, q21_t_old;


/**********************************************************/
/**
 * @brief      calculates and returns the collisional de-excitation coefficient q21
 *
 * @param [in] struct lines *  line_ptr   A single line
 * @param [in] double  t   The temperature at wich to calculate the coefficient
 * @return    The collisional de-excitiation coefficient
 *
 * @details
 * This routine uses stored coefficients if we have them, or the Van Regemortor
 * approximation if we do not.
 *
 * ### Notes ###
 * the relevant paper to consult here is Van Regemorter 1962. We use an effective gaunt
 * factor to calculate collision strengths. There is one regime in which kt < hnu. For that
 * consult equation 4.20 and 4.21 of Hazy.
************************************************************/

double
q21 (line_ptr, t)
     struct lines *line_ptr;
     double t;
{
  double gaunt;
//OLD  double gbar;
  double omega;
  double u0;
  double upsilon ();


  if (q21_line_ptr != line_ptr || t != q21_t_old)
  {


    u0 = (BOLTZMANN * t) / (H * line_ptr->freq);

    if (line_ptr->istate == 1 && u0 < 2)        // neutrals at low energy. Used 2 to give continuous function.
      gaunt = u0 / 10.0;
    else                        // low energy electrons, positive ions
      gaunt = 0.2;



    if (line_ptr->coll_index < 0)       //if we do not have a collision strength for this line use the g-bar formulation
    {
      omega = ECS_CONSTANT * line_ptr->gl * gaunt * line_ptr->f / line_ptr->freq;
    }
    else                        //otherwise use the collision strength directly. NB what we call omega, most people including hazy call upsilon.
    {
      omega = upsilon (line_ptr->coll_index, u0);
//OLD      gbar = omega / ECS_CONSTANT / line_ptr->gl / line_ptr->f * line_ptr->freq;
    }


    q21_a = 8.629e-6 / (sqrt (t) * line_ptr->gu) * omega;
    q21_t_old = t;
  }

  return (q21_a);
}


/**********************************************************/
/**
 * @brief      Calculate the collisional excitation coefficient for a line
 *
 * @param [in] struct lines *  line_ptr   The line of interest
 * @param [in] double  t   The temperature of interest
 * @return     The collisional excitation coeffient
 *
 * @details
 * The routine calls q21 and uses detailed balance to derive q12
 *
 * ### Notes ###
 *
 **********************************************************/

double
q12 (line_ptr, t)
     struct lines *line_ptr;
     double t;
{
  double x;
  double q21 ();
  double exp ();

  x = line_ptr->gu / line_ptr->gl * q21 (line_ptr, t) * exp (-H_OVER_K * line_ptr->freq / t);

  return (x);
}


/*
   a21 alculates and returns the Einstein A coefficient
   History:
   98aug        ksl     Coded and debugged
   99jan        ksl Modified so would shortcircuit calculation if
   called multiple times for same a
 */
#define A21_CONSTANT 7.429297e-22       // 8 * PI * PI * E * E / (MELEC * C * C * C)

struct lines *a21_line_ptr;
double a21_a;


/**********************************************************/
/**
 * @brief      Calculate the Einstein A coefficient for aline
 *
 * @param [in out] struct lines *  line_ptr   The structure that describes a single line
 * @return     The Einstein A coefficient
 *
 * @details
 * Using the oscillator strecnth and the multiplicity (which are
 * read in) calculate A
 *
 * ### Notes ###
 *
 **********************************************************/

double
a21 (line_ptr)
     struct lines *line_ptr;
{
  double freq;

  if (a21_line_ptr != line_ptr)
  {
    freq = line_ptr->freq;
    a21_a = A21_CONSTANT * line_ptr->gl / line_ptr->gu * freq * freq * line_ptr->f;
    a21_line_ptr = line_ptr;
  }

  return (a21_a);
}


/**********************************************************/
/**
 * @brief      calculates the thermally averaged collision strength thermally excited line emission.
 *
 * @param [in out] int  n_coll   the index of the collision strength record we are working with
 * @param [in out] double  u0  - kT_e/hnu - where nu is the transition frequency for the line of interest and T_e is the electron temperature
 * @return     upsilon - the thermally averaged collision strength for a given line at a given temp
 *
 * @details
 *
 * ### Notes ###
 * It uses data extracted from Chianti stored in coll_stren.
 * The paper to consult is Burgess and Tully A&A 254,436 (1992).
 * u0 is the ratio of Boltzmans constant times the electron temperature in the cell
 * divided by Plancks constant times the frequency of the line under analysis.
 *
 *
 **********************************************************/

double
upsilon (n_coll, u0)
     int n_coll;
     double u0;
{
  double x;                     //The scaled temperature
  double y;                     //The scaled collision sterngth
  double upsilon;               //The actual collision strength
  int linterp ();

  /* first we compute x. This is the "reduced temperature" from
     Burgess & Tully 1992. */
  x = 0.0;
  if (coll_stren[n_coll].type == 1 || coll_stren[n_coll].type == 4)
  {
    x = 1. - (log (coll_stren[n_coll].scaling_param) / log (u0 + coll_stren[n_coll].scaling_param));
  }
  else if (coll_stren[n_coll].type == 2 || coll_stren[n_coll].type == 3)
  {
    x = u0 / (u0 + coll_stren[n_coll].scaling_param);
  }
  else
  {
    Error ("upsilon - coll_stren %i has no type %g\n", coll_stren[n_coll].type);
    exit (0);
  }


  /* we now compute y from the interpolation formulae
     y is the reduced upsilon from Burgess & Tully 1992. */
  linterp (x, coll_stren[n_coll].sct, coll_stren[n_coll].scups, coll_stren[n_coll].n_points, &y, 0);

  /*  now we extract upsilon from y  - there are four different parametrisations */

  upsilon = 0.0;
  if (coll_stren[n_coll].type == 1)
  {
    upsilon = y * (log (u0 + exp (1)));
  }
  else if (coll_stren[n_coll].type == 2)
  {
    upsilon = y;
  }
  else if (coll_stren[n_coll].type == 3)
  {
    upsilon = y / (u0 + 1);
  }
  else if (coll_stren[n_coll].type == 4)
  {
    upsilon = y * (log (u0 + coll_stren[n_coll].scaling_param));
  }
  else
  {
    Error ("upsilon - coll_stren %i has no type %g\n", coll_stren[n_coll].type);
    exit (0);
  }
  return (upsilon);
}






/**********************************************************/
/**
 * @brief Perform linear/logarithmic interpolation of an array
 *
 * @param [in] double  value   The value used for interpoaltion
 * @param [in] double  array[]   An array containing a set of asceding values
 * @param [in] int  npts   The size of the array
 * @param [out] int *  ival  The lower index to use in the interpolation
 * of the array which need to be used to interpolate on
 * @param [out] double *  f   The fraction of the upper point to use in the interpolation
 * @param [in] int  mode  A switch to choose linear(0)  or lograrithmic(1) interpolation
 * @return     Usually returns 0, but returns -1 if the input value is less
 * than the first elememnt in the array, and 1 if it is greater than the
 * last element int he array.  In either of these cases, the fractions are
 * set up only to access one array element
 *
 * @details
 *
 *
 * Typically one has two parallel arrays, one containing a set of values,
 * in the original case frequencies, and another containing some function
 * of those values.  This routine finds the fraction of the function at
 * two point in the data array to interpolate.
 *
 * The values for the input array can be interpolated linear or logarithmically.
 * that is the fraction that is returned are based on the values in the
 * array or the logarithm of them.
 *
 *
 * The routine uses bisection
 *
 *
 *
 * ### Notes ###
 *
 * fraction is a utility written to speed up the search for the bracketing
 * topbase photionization x-sections, but in principle it should work in
 * other situations within python (which at one time at least contributed
 * significantly to the runtime for Python).
 *
 * Today, fraction is called directly in the routines that are used to set up coordinate
 * systems, and indirectly through linterp (below) which should be inspected
 * to see how the logarithmic interpolation is supposed to work.
 *
 * The routine is similar to the numerical * recipes routine locate.
 * There may be a gsl routine as well.  The routine
 * should probably be replaced.
 *
 *
 **********************************************************/

int
fraction (value, array, npts, ival, f, mode)
     double array[];            // The array in we want to search
     int npts, *ival;           // ival is the lower point
     double value;              // The value we want to index
     double *f;                 // The fractional "distance" to the next point in the array
     int mode;                  // 0 = compute in linear space, 1=compute in log space
{
  int imin, imax, ihalf;

  if (value < array[0])
  {
    *ival = 0;
    *f = 0.0;
    return (-1);
  }

  imax = npts - 1;
  if (value > array[imax])
  {
    *ival = npts - 2;
    *f = 1.0;
    return (1);
  }



  imin = 0;

/* In what follows, there is a specific issue as to how to treat
the situation where the value is exactly on an array element. Right
now this is set to try to identify the array element below the one
on which the value sits, and to set the fraction to 1.  This was
to reflect the behavior of the search routine in where_in_grid. */

  while (imax - imin > 1)
  {
    ihalf = (imin + imax) >> 1; // Compute a midpoint >> is a bitwise right shift
    if (value > array[ihalf])
    {
      imin = ihalf;
    }
    else
      imax = ihalf;
  }

// So array[imin] just <= value

  if (mode == 0)
    *f = (value - array[imin]) / (array[imax] - array[imin]);   //linear interpolation
  else if (mode == 1)
    *f = (log (value) - log (array[imin])) / (log (array[imax]) - log (array[imin]));   //log interpolation
  else
  {
    Error ("Fraction - unknown mode %i\n", mode);
    exit (0);
    return (0);
  }

  *ival = imin;

  return (0);
}







/**********************************************************/
/**
 * @brief      Perform a linear interpolation on two parallel arrays, the first
 * of which contains a set of values to be interpolated and the second of
 * which has the function at those values
 *
 * @param [in] double  x   A value
 * @param [in] double  xarray[]   The array that is interpolated
 * @param [in] double  yarray[]   The array that contains a function of the values in xarray
 * @param [in] int  xdim   The length of the two arrays
 * @param [out] double *  y   The resulting interpolated value
 * @param [in] int  mode   A switch to choose linear(0) or "logarithmic" (1)
 * interpolation
 *
 * @return     The number of the array element that is used for the lower of the two
 * elements that are interpolated on.
 *
 * @details
 * Given a number x, and an array of x's in xarray, and functional
 * values y = f(x) in yarray and the dimension of xarray and yarray,
 * linterp calculates y, the linearly interpolated value of y(x). It
 * also returns the element in the xarray so that one can consider
 * not doing a search to find the right array element in certain
 * circumstances

 *
 * ### Notes ###
 *
 * For mode 0, the value that is returned is
 *
 * (1-f)*y[nelem]+f*[nelem+1)
 *
 * For mode 1, the value returned is
 * exp ((1. - f) * log (y[nelem]) + f * log (y[nelem + 1]))
 *
 *
 **********************************************************/

int
linterp (x, xarray, yarray, xdim, y, mode)
     double x;                  // The value that we wish to index i
     double xarray[], yarray[];
     int xdim;
     double *y;
     int mode;                  //0 = linear, 1 = log
{
  int nelem = 0;
  double frac;


  fraction (x, xarray, xdim, &nelem, &frac, mode);

  if (mode == 0)
    *y = (1. - frac) * yarray[nelem] + frac * yarray[nelem + 1];
  else if (mode == 1)
    *y = exp ((1. - frac) * log (yarray[nelem]) + frac * log (yarray[nelem + 1]));
  else
  {
    Error ("linterp - unknown mode %i\n", mode);
    exit (0);
  }

  return (nelem);

}
