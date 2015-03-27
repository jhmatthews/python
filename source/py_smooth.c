
/*                    Space Telescope Science Institute


   Synopsis:

   The program py_smooth  is designed to smooth a spectrum produced
   by python  

   It uses the rdpar file genspec.pf

   Description: 


   The program requires a series of input files to run whose names are input
   parameters.  All inputs come through a rdpar file

   lambda.min      Minimum wavelength
   lambda.max      Maximum wavelength
   resolution      Used to establish dlambda
   lambda.rescale  wavelength at which to rescale flux.  If lambda.rescale
   lies outside range defined by lambda min and max
   there will be no rescaling
   flux.rescale    Flux to which to rescale flux
   instrument.res  Optional convolution of the spectrum with a Gaussion
   to simulate.  If less than zero no convolution is
   performed; if > 0, then value is taken to be sigma
   of resolution
   out.file



   Notes:


   The way in which the gaussian smoothing is done is by brute
   force and could likely be improved significantly.  If NSPEC is increased
   then the read statement for the input spectra needs to be changed.


   History:
	97 May	ksl	Adapted from genspec
	97 Aug	ksl	Added routines to read and write files in "models" directory on the 
			MAC.  Modified the program so that it handles the possibility that
			there are fewer than NSPEC separate spectra.  
	97oct9	ksl	Implemented rescaling
	02feb	ksl	Eliminated the MAC bits
	04dec	ksl	Added a few lines to eliminate warnings when compiled
			with -O3
	090304	ksl	Added command line options to make it easier to do many spectra
			at once

   *********************************************************************** */


#include	<stdio.h>
#include	<math.h>
#include    <string.h>
#include    <stdlib.h>

#include "atomic.h"
#include "python.h"


#define NW		10000
#define NSPEC           20


typedef struct sss
{
  float freq, wav, x[NSPEC];
}
s_dummy, *SpPtr;


int
main (argc, argv)
     int argc;
     char *argv[];

{
  FILE *fp, *fout, *fopen ();
  char pffile[LINELENGTH];
  char root[LINELENGTH];

  float wavmin, wavmax, res, lambda, q;
  float wwavmin, wwavmax;
  float lrescale, frescale, xsigma;

  int nspec;			/* The actual number of spectra read in, not to be confused with NSPEC */
  int nwords;

  char cont_file[LINELENGTH], out_file[LINELENGTH];
  char line[LINELENGTH];

  int imax;			/* number of points in input spectrum */


  int npts, i, n, nn;
  float gauss ();
  float fwhm;

  SpPtr spec_in, spec_out;



  /* initialize a few variables to avoid warnings */
  fp = fout = NULL;
  nwords = 0;
  imax = 0;

  if ((spec_in = (SpPtr) calloc (sizeof (s_dummy), NW)) == NULL)
    {
      printf ("Could not allocate spec_in\n");
      exit (0);
    }
  if ((spec_out = (SpPtr) calloc (sizeof (s_dummy), NW)) == NULL)
    {
      printf ("Could not allocate spec_out\n");
      exit (0);
    }


/* Set initial values for everything */
  wavmin = 850;
  wavmax = 1850;
  res = 0.57;
  lrescale = 1350;
  frescale = 1e-11;
  fwhm = 3;
  strcpy (cont_file, "mod.spec");
  strcpy (out_file, "mod.spec_smo");


/* Now read the command line */

  if (argc == 1)
    strcpy (pffile, "py_smooth.pf");
  else
    {

      for (i = 1; i < argc; i++)
	{

	  if (strcmp (argv[i], "-h") == 0)
	    {
	      printf ("Usage: py_smooth [-s spec]  [file.pf]\n");
	      exit (0);
	    }
	  else if (strcmp (argv[i], "-s") == 0)
	    {
	      if (sscanf (argv[i + 1], "%s", cont_file) != 1)
		{
		  Error
		    ("py_smooth: Expected source spectrum  after -s switch\n");
		  exit (0);
		}
	      if (strcmp (out_file, "mod.spec_smo") == 0)
		{
		  strcpy (out_file, "construct");
		}
	      i++;

	    }
	  else if (strcmp (argv[i], "-o") == 0)
	    {
	      if (sscanf (argv[i + 1], "%s", out_file) != 1)
		{
		  Error
		    ("py_smooth: Expected output file  after -o switch\n");
		  exit (0);
		}
	      i++;

	    }
	  else if (strncmp (argv[i], "-", 1) == 0)
	    {
	      Error ("python: Unknown switch %s\n", argv[i]);
	      exit (0);
	    }



	}
      /* The last argument is always the pf file */
      strcpy (pffile, argv[argc - 1]);
      if (strstr (pffile, ".pf") == NULL)
	strcat (pffile, ".pf");
    }
  printf ("Opening parameter file %s\n", pffile);

  opar (pffile);


/* Now read the input data */


  rdflo ("lambda.min", &wavmin);
  rdflo ("lambda.max", &wavmax);
  rdflo ("delta.lambda", &res);

  rdflo ("lambda.rescale", &lrescale);
  rdflo ("flux.rescale", &frescale);
  rdflo ("instrument.fwhm", &fwhm);


  if (res > 1)
    printf ("Warning: Resolution is now bin width in Angstroms\n");

  xsigma = fwhm / (sqrt (8. * log (2.)));

  lambda = wavmin;
  for (npts = 0; npts < NW; npts++)
    {
      spec_out[npts].wav = lambda;
      spec_out[npts].freq = HC / (lambda * 1.e-8);
      if (lambda > wavmax)
	break;
      lambda += res;
    }

  printf ("Number of points in output spectrum %d\n", npts);


  if (strcmp (cont_file, "mod.spec") == 0)
    {
      rdstr ("cont.file", cont_file);
    }
  else
    {
      Log ("Using spectrum %s\n", cont_file);
    }


  if (strcmp (out_file, "construct") == 0)
    {
      strcpy (out_file, "");
      strcpy (out_file, cont_file);
      strcat (out_file, "_smo");
    }
  else if (strcmp (out_file, "mod.spec_smo") == 0)
    {
      strcpy (out_file, "");
      get_root (out_file, cont_file);
      strcat (out_file, ".spec_smo");
      rdstr ("out.file", out_file);
    }


  if (strncmp (root, "mod", 3) == 0)
    cpar ("mod_smooth.pf");
  else
    cpar ("py_smooth.pf");


/* Generate the continuum */

  if (strncmp (cont_file, "none", 4) != 0)
    {
      if ((fp = fopen (cont_file, "r")) == NULL)
	{
	  printf ("Failed to open %s\n", cont_file);
	  exit (0);
	}
      /* Now open the output file */
      fout = fopen (out_file, "w");


      imax = 0;
      while (fgets (line, LINELENGTH, fp) != NULL)
	{
	  if (line[0] != '#')
	    {
	      if ((nwords =
		   sscanf (line,
			   "%e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e",
			   &spec_in[imax].freq, &spec_in[imax].wav,
			   &spec_in[imax].x[0], &spec_in[imax].x[1],
			   &spec_in[imax].x[2], &spec_in[imax].x[3],
			   &spec_in[imax].x[4], &spec_in[imax].x[5],
			   &spec_in[imax].x[6], &spec_in[imax].x[7],
			   &spec_in[imax].x[8], &spec_in[imax].x[9],
			   &spec_in[imax].x[10], &spec_in[imax].x[11],
			   &spec_in[imax].x[12], &spec_in[imax].x[13],
			   &spec_in[imax].x[14], &spec_in[imax].x[15],
			   &spec_in[imax].x[16], &spec_in[imax].x[17],
			   &spec_in[imax].x[18], &spec_in[imax].x[19])) != 22)
		{
//                                              printf("Error: Incorrect number of spectra in input spectrum %d\n",nwords);
		  //                                              exit(0);
		}
	      imax++;
	      if (imax > NWAVE)
		{
		  printf ("Error: Reading input file imax > NWAVE\n");
		  exit (0);

		}
	    }
	  else
	    fprintf (fout, "%s", line);
	}
    }

/* Now convolve the spectrum with an instrumental response function */
  imax--;
  nspec = nwords - 2;
  fclose (fp);

  for (n = 0; n < npts; n++)
    {
      wwavmin = spec_out[n].wav;
      wwavmax = spec_out[n + 1].wav;

      for (i = 0; i < imax; i++)
	{
	  if ((q = gauss (wwavmin, wwavmax, spec_in[i].wav, xsigma)) > 0.0)
	    {
	      for (nn = 0; nn < nspec; nn++)
		{
		  spec_out[n].x[nn] +=
		    q * (spec_in[i].x[nn]) * (spec_in[i].wav -
					      spec_in[i + 1].wav) / (wwavmax -
								     wwavmin);
		}
	    }
	}

    }

/* Now rescale the spectrum if desired */
  printf ("OK checking %e %e %e\n", wavmin, lrescale, wavmax);
  if (wavmin < lrescale && lrescale < wavmax)
    {
      printf ("OK rescaling to frescale %e at lrescle %e\n", frescale,
	      lrescale);
      for (i = 0; i < nspec; i++)
	{
	  n = 0;
	  while (spec_out[n].wav < lrescale)
	    n++;
	  q = frescale / spec_out[n].x[i];
	  for (n = 0; n < npts; n++)
	    spec_out[n].x[i] *= q;
	}
    }
/* Now print out the spectrum */

  for (i = 0; i < npts; i++)
    {
      fprintf (fout, "%e %8.3f", spec_out[i].freq, spec_out[i].wav);
      for (nn = 0; nn < nspec; nn++)
	{
	  fprintf (fout, " %8.2e", spec_out[i].x[nn]);
	}
      fprintf (fout, "\n");

    }
  fclose (fout);

  return EXIT_SUCCESS;
}


#define SIGMA 4.

/* gauss returns the integral of a normalized gaussina in the 
   interval specified by xmin and xmax */

#define		ISTEPS	20
#define		LIMIT	4

float
gauss (xmin, xmax, x, sigma)
     float xmin, xmax, x, sigma;

{
  double exp ();
  double xx, delta;
  long i;
  float a;

/* Check to see whether integration is warranted */

  if (xmax < x - LIMIT * sigma)
    return (0.0);
  if (xmin > x + LIMIT * sigma)
    return (0.0);
  if (xmin < x - LIMIT * sigma && xmax > x + LIMIT * sigma)
    return (1.0);

/* Now do the integration */

  delta = (xmax - xmin) / (ISTEPS);

  a = 0.0;

  for (i = 0; i < ISTEPS; i++)
    {
      xx = (xmin + delta * (i + 0.5) - x) / sigma;
      a += exp ((-0.5) * xx * xx);
    }

  a /= 2.5066283 * sigma;
  a *= delta;

  return (a);

}
