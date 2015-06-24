/*
 * The subroutines in this file handle outputting spectrum array reverb data
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "atomic.h"
#include "python.h"

/***********************************************************
Synopsis:
	delay_dump(filename, p, f1, f2, nspec)
		Calculates the delay for the photons. Cloned from spectra.c

Arguments:
	filename		File to output to
	PhotPtr p		Array of photons (pass global p)
	int nspec		Observer spectrum to dump for

Returns:
	To file only

Description:
	Prep wipes clean the output file and writes the header, dump then
	builds up the output file in batches.

Notes:

History:
	14aug	-	Written by SWM
***********************************************************/
char	delay_dump_file[LINELENGTH];
int		delay_dump_bank_size = 65535, delay_dump_bank_curr = 0;
PhotPtr	delay_dump_bank;
int     *delay_dump_bank_ex;

double
delay_to_observer(PhotPtr pp)
{
	plane_dummy	observer;
	observer.x[0] = pp->lmn[0] * geo.rmax;
	observer.x[1] = pp->lmn[1] * geo.rmax;
	observer.x[2] = pp->lmn[2] * geo.rmax;
	observer.lmn[0] = pp->lmn[0];
	observer.lmn[1] = pp->lmn[1];
	observer.lmn[2] = pp->lmn[2];
	return (pp->path + ds_to_plane(&observer, pp));
}

int
delay_dump_prep(char filename[], int restart_stat, int i_rank)
{
	FILE *fopen(), *fptr;
	char string[LINELENGTH], c_file[LINELENGTH], c_rank[LINELENGTH];
	int	i;

	//Get output filename
		strcpy(c_file, filename);		//Copy filename to new string
		strcat(c_file, ".delay_dump");
	if (i_rank > 0) {
		sprintf(c_rank, "%i", i_rank);	//Write thread to string
		strcat(c_file, c_rank);			//Append thread to filename
	}
	strcpy(delay_dump_file, c_file);	//Store modified filename for later

	//Allocate and zero dump files and set extract status
	delay_dump_bank = (PhotPtr) calloc(sizeof(p_dummy), delay_dump_bank_size);
	delay_dump_bank_ex = (int *)calloc(sizeof(int), delay_dump_bank_size);
	for (i = 0; i < delay_dump_bank_size; i++)
		delay_dump_bank_ex[i] = 0;

	//Check whether we should continue
	if (restart_stat == 1) {
		Log("delay_dump_prep: Resume run, skipping writeout\n");
		return (0);
	}

	/* If this isn't a continue run, prep the output file */
	if ((fptr = fopen(delay_dump_file, "w")) == NULL) {
		Error("delay_dump_prep: Unable to open %s for writing\n",
		      delay_dump_file);
		exit(0);
	}

	if (i_rank > 0) 
	{
		fprintf(fptr, "# Delay dump file for slave process %d\n", i_rank);
	}
	else
	{
		// Construct and write a header string for the output file
		fprintf(fptr, "# Python Version %s\n", VERSION);
		get_time(string);
		fprintf(fptr, "# Date	%s\n#  \n", string);
		fprintf(fptr, "# \n# Freq      Wavelength  Weight   "
			" Last X     Last Y     Last Z    "
		    " Scatters   RScatter   Delay      Extracted  Spectrum   Origin   \n");
	}
	fclose(fptr);
	return (0);
}

/***********************************************************
Synopsis:
	delay_dump_finish()
		Outputs any leftover photons.

Arguments:
	None

Returns:
	To file only

Description:
	Dumps last remaining extracted photons to file and frees arrays.

Notes:

History:
	6/2/15	-	Written by SWM
***********************************************************/
int 
delay_dump_finish (void)
{
	if (delay_dump_bank_curr > 0) 
	{
		delay_dump(delay_dump_bank, delay_dump_bank_curr - 1, 1);
	}
	free(delay_dump_bank);
	free(delay_dump_bank_ex);
	return (0);
}

/***********************************************************
Synopsis:
	delay_dump_combine(iRanks)
		Finishes off delay dump for multithreaded runs

Arguments:
	int iRanks	Combines delay dump file from each thread

Returns:
	To file only

Description:
	Uses system commands to achieve this quickly

Notes:

History:
	6/2/15	-	Written by SWM
***********************************************************/
int
delay_dump_combine(int iRanks)
{
	FILE *fopen();
	char cCall     [LINELENGTH];

	//Yes this is done as a system call and won 't work on Windows machines. Lazy solution!
	sprintf(cCall, "cat %s[0-9]* >> %s", delay_dump_file, delay_dump_file);
	if (system(cCall) < 0)
		Error("delay_dump_combine: Error calling system command '%s'", cCall);
	return (0);
}

/***********************************************************
Synopsis:
	delay_dump(PhotPtr p, int np, int nspec, int iExtracted)
		Dumps array of photons provided to file

Arguments:
	PhotPtr p 		Array of photons to dump
	int np 			Length of array
	int iExtracted 	Whether or not this array is of extracted photons

Returns:
	To file only

Description:
	Dumps to this thread's delay_dump_file the passed photon
	array.

Notes:

History:
	6/2/15	-	Written by SWM
***********************************************************/
int
delay_dump(PhotPtr p, int np, int iExtracted)
{
	FILE *fopen(), *fptr;
	int	 nphot, mscat, mtopbot, i;
	double zangle, minpath;

	if(geo.system_type == SYSTEM_TYPE_STAR || geo.system_type == SYSTEM_TYPE_BINARY)
		minpath = geo.rmax - geo.rstar;	
	else if(geo.system_type == SYSTEM_TYPE_AGN)
		minpath = geo.r_agn;

	
	/*
	 * Open a file for writing the spectrum
	 */
	if ((fptr = fopen(delay_dump_file, "a")) == NULL) 
	{
		Error("delay_dump: Unable to reopen %s for writing\n",
		      delay_dump_file);
		exit(0);
	}
	for (nphot = 0; nphot < np; nphot++) {
		if (p[nphot].istat == P_ESCAPE && p[nphot].nscat > 0) {
			zangle = fabs(p[nphot].lmn[2]);
			/*
			 * Complicated if statement to allow one to choose
			 * whether to construct the spectrum from all photons
			 * or just from photons which have scattered a
			 * specific number of times.  01apr13--ksl-Modified
			 * if statement to change behavior on negative
			 * numbers to say that a negative number for mscat
			 * implies that you accept any photon with |mscat| or
			 * more scatters
			 */
			for (i = MSPEC; i < nspectra; i++) {
				if (((mscat = xxspec[i].nscat) > 999 || p[i].nscat == mscat
				|| (mscat < 0 && p[nphot].nscat >= (-mscat)))
				    && ((mtopbot = xxspec[i].top_bot) == 0
					|| (mtopbot * p[nphot].x[2]) > 0)) {
					if (xxspec[i].mmin < zangle && zangle < xxspec[i].mmax) {
						fprintf(fptr,
							"%10.5g %10.5g %10.5g %+10.5g %+10.5g %+10.5g %3d     %3d     %10.5g %5d %5d %5d\n",
							p[nphot].freq, C * 1e8 / p[nphot].freq, p[nphot].w, 
							p[nphot].x[0], p[nphot].x[1], p[nphot].x[2],
							p[nphot].nscat, p[nphot].nrscat,
							(delay_to_observer(&p[nphot]) - (geo.rmax-minpath)) / C,
							(iExtracted ? delay_dump_bank_ex[nphot] : 0),
							i - MSPEC, p[nphot].origin);
					}
				}
			}
		}
	}
	fclose(fptr);
	return (0);
}

/***********************************************************
Synopsis:
	delay_dump_single(PhotPtr pp, int extract_phot)
		Finishes off delay dump for multithreaded runs

Arguments:
	PhotPtr pp 			Photon to dump
	int extract_phot	Whether or not this was extracted

Returns:
	To internal arrays only

Description:
	Pushes the passed photon into the delay_dump_bank. If
	this fills it, dumps current bank to file and restarts it.

Notes:

History:
	6/2/15	-	Written by SWM
***********************************************************/
int
delay_dump_single(PhotPtr pp, int extract_phot)
{
	stuff_phot(pp, &delay_dump_bank[delay_dump_bank_curr]); 		//Bank single photon in temp array
		delay_dump_bank_ex[delay_dump_bank_curr] = extract_phot; 	//Record if it's extract photon
	if (delay_dump_bank_curr == delay_dump_bank_size - 1)			//If temp array is full
	{
		delay_dump(delay_dump_bank, delay_dump_bank_size, 1);
		delay_dump_bank_curr = 0;										//Dump to file, zero array position
		int		i;
		for (i = 0; i < delay_dump_bank_size; i++)
			delay_dump_bank_ex[i] = 0;									//Zero extract status of single array
	}
	else
	{
		delay_dump_bank_curr++;
	}
	return (0);
}

/***********************************************************
Synopsis:
	path_data_constructor(double r_rad_min, double r_rad_max,
						  int i_bins, int i_angles)
		Returns global data container for wind path info

Arguments:
	r_rad_min 		Minimum disk radius
	r_rad_max 		Maximum disk radius, used for making
					the path distance arrays
	i_bins 			Number of path distance bins
	i_angles 		Number of observer angles

Returns:
	Pointer to prepped global path data info.

Description:
	Sets up global data used by the individual cells for
	binning purposes.

Notes:

History:
	9/2/15	-	Written by SWM
***********************************************************/
Path_Data_Ptr
path_data_constructor(double r_rad_min, double r_rad_max, int i_path_bins,
		      int i_angles, double freqmin, double freqmax,
		      int i_theta_res)
{
	int		i;
	Path_Data_Ptr	data = (Path_Data_Ptr) calloc(sizeof(path_data_dummy), 1);

	if (data == NULL) 
	{
		Error("path_data_constructor: Could not allocate memory\n");
		exit(0);
	}
	data->i_theta_res = i_theta_res;
	data->i_obs = i_angles;
	data->i_path_bins = i_path_bins;
	data->ad_path_bin = (double *)calloc(sizeof(double), i_path_bins + 1);
	for (i = 0; i <= i_path_bins; i++) 
	{
		data->ad_path_bin[i] =
			r_rad_min + i * (r_rad_max * 5.0 - r_rad_min) / i_path_bins;
	}
	data->ad_freq_bin = (double *)calloc(sizeof(double), NWAVE);
	for (i = 0; i < NWAVE; i++) 
	{
		data->ad_freq_bin[i] = freqmin + i * (freqmax - freqmin) / (NWAVE - 1);
	}
	return (data);
}

int
path_data_init(double r_rad_min, double r_rad_max, int i_path_bins,
	       int i_angles, double r_freq_min, double r_freq_max,
	       int i_theta_res)
{
	g_path_data =
	(Path_Data_Ptr) path_data_constructor(r_rad_min, r_rad_max, i_path_bins,
					   i_angles, r_freq_min, r_freq_max,
					      i_theta_res);
	return (0);
}

/***********************************************************
Synopsis:
	wind_paths_constructor(WindPtr wind)
		Constructs a wind path record and returns it.

Arguments:
	Wind_Ptr wind	Wind cell the paths are being made for

Returns:
  	Pointer to finished wind paths

Description:
	Constructs a wind cell's path data and returns a pointer
	to it. Also calculates the distance to each extract
	viewpoint from the cell via a dummy photon.

Notes:

History:
	9/2/15	-	Written by SWM
***********************************************************/
Wind_Paths_Ptr
wind_paths_constructor(WindPtr wind)
{
	Wind_Paths_Ptr	paths =
	(Wind_Paths_Ptr) calloc(sizeof(wind_paths_dummy), 1);

	if (paths == NULL) 
	{
		Error
			("wind_paths_constructor: Could not allocate memory for cell %d\n",
			 wind->nwind);
		exit(0);
	}
	paths->ad_freq_flux = (double *)calloc(sizeof(double), NWAVE);
	paths->ai_freq_num = (int *)calloc(sizeof(int), NWAVE);
	paths->ad_freq_path_flux =
		(double *)calloc(sizeof(double), NWAVE * g_path_data->i_path_bins);
	paths->ai_freq_path_num =
		(int *)calloc(sizeof(int), NWAVE * g_path_data->i_path_bins);
	if (paths->ad_freq_path_flux == NULL || paths->ad_freq_flux == NULL
	 || paths->ai_freq_path_num == NULL || paths->ai_freq_num == NULL) 
	{
		Error
			("wind_paths_constructor: Could not allocate memory for cell %d bins\n",
			 wind->nwind);
		exit(0);
	}
	return (paths);
}

/***********************************************************
Synopsis:
	reverb_init(WindPtr wind, int nangles, path_bins, theta_bins,
				double freqmin, freqmax)
		Initialises module variables

Arguments:
	WindPtr wind 	Wind cell array, passed further down
	int nangles 	Number of angles used
	int path_bins 	Number of bins to sort photon paths by
	int theta_bins  For output only, number of angular bins
	double freqmin
	double freqmax 	Frequency range to record paths over

Returns:

Description:

Notes:

History:
	3/15	-	Written by SWM
***********************************************************/
int
reverb_init(WindPtr wind, int nangles, double freqmin, double freqmax)
{
	if (geo.reverb == REV_WIND) 
	{
		if (geo.wind_radiation == 0) 
		{
			Error
			("Wind radiation is off but wind-based path tracking is enabled!\n");
			exit(0);
		} 
		else 
		{
			Log ("Wind cell-based path tracking is enabled for frequency range %g-%g\n",
				 freqmin, freqmax);
			path_data_init(0.0, geo.rmax, geo.reverb_path_bins, nangles,
				   freqmin, freqmax, geo.reverb_theta_bins);
			wind_paths_init(wind);
		}
	}
	else if (geo.reverb == REV_PHOTON)
		Log("Photon-based path tracking is enabled.\n");


	return (0);
}

/***********************************************************
Synopsis:
	wind_paths_init(WindPtr wind)
		Initialises wind path structures

Arguments:
	WindPtr wind	Wind to initialise

Returns:

Description:
	Iterates over each wind cell, and sets up the wind path
	data therein. This should ideally be folded into the
	regular wind setup at some point.

Notes:

History:
	10/2/15	-	Written by SWM
***********************************************************/
int
wind_paths_init(WindPtr wind)
{
	int	i;

	for (i = 0; i < NDIM * MDIM; i++) 
	{
		wind[i].paths = (Wind_Paths_Ptr) wind_paths_constructor(&wind[i]);
	}
	return (0);
}

/***********************************************************
Synopsis:
	wind_paths_add_phot(Wind_Paths_Ptr paths, PhotPtr pp)
		Adds a photon's weight to a wind cell's delay array.

Arguments:
	PhotPtr pp 				Photon to dump
	Wind_Paths_Ptr paths	Paths (from wind) to add to

Returns:

Description:
	Adds photon to frequency and path bin as appropriate.
	As wind array is 2d, .5th dimension added via keeping
	two path sets for + and - x-axis positions). No full 3-d
	treatment required as Extract mode's viewpoint always
	lies along X-axis.

Notes:

History:
	9/2/15	-	Written by SWM
***********************************************************/
int
wind_paths_add_phot(WindPtr wind, PhotPtr pp)
{
	int		i         , j;

	for (i = 0; i < NWAVE - 1; i++) 
	{
		if (pp->freq >= g_path_data->ad_freq_bin[i] &&
		    pp->freq < g_path_data->ad_freq_bin[i + 1]) 
		{
			for (j = 0; j < g_path_data->i_path_bins; j++) 
			{
				if (pp->path >= g_path_data->ad_path_bin[j] &&
				    pp->path <= g_path_data->ad_path_bin[j + 1]) 
				{
					wind->paths->ad_freq_path_flux[i * g_path_data->i_path_bins + j] += pp->w;
					wind->paths->ai_freq_path_num[i * g_path_data->i_path_bins +  j]++;
				}
			}
		}
	}
	return (0);
}

/***********************************************************
Synopsis:
	wind_paths_gen_phot(Wind_Paths_Ptr paths, PhotPtr pp)
		Adds a path delay from the origin cell to photon

Arguments:
	PhotPtr pp 				Photon to set path of
	Wind_Paths_Ptr paths	Paths (from wind) to use

Returns:

Description:
	Picks a weighted random frequency bin, then weighted
	random path bin, then assigns the photon a path selected
	from within that bin (nonweighted uniform randomly).

Notes:
	May want to correlate photon energy with wind path in
	future!

History:
	26/3/15	-	Written by SWM
***********************************************************/
int
wind_paths_gen_phot(WindPtr wind, PhotPtr pp)
{
	double		r_rand , r_total;
	int		i_freq    , i_path;

	i_freq = -1;
	i_path = -1;

	r_total = 0.0;
	r_rand = wind->paths->d_flux * rand() / MAXRAND;
	while (r_rand < r_total) 
	{
		r_total += wind->paths->ad_freq_flux[++i_freq];
		if (i_freq >= NWAVE) 
		{
			Error
				("wind_paths_gen_phot: No path data in wind cell %d at %g %g %g\n",
			   wind->nwind, wind->x[0], wind->x[1], wind->x[2]);
			exit(0);
		}
	}

	r_total = 0;
	r_rand = wind->paths->d_flux * rand() / MAXRAND;
	while (r_rand < r_total) 
	{
		r_total += wind->paths->ad_freq_flux[i_freq * g_path_data->i_path_bins +
						  (++i_path)];
		if (i_path > g_path_data->i_path_bins)
		{
			Error("wind_paths_gen_phot: No path data in wind cell %d at %g %g %g\n",
			   wind->nwind, wind->x[0], wind->x[1], wind->x[2]);
			exit(0);
		}
	}

	pp->path = g_path_data->ad_freq_bin[i_path - 1] +
		(g_path_data->ad_freq_bin[i_path] -
		 g_path_data->ad_freq_bin[i_path - 1]) * (rand() / MAXRAND);
	return (0);
}


/***********************************************************
Synopsis:
	wind_paths_evaluate(WindPtr wind)
		Evaluates wind path details for a cycle

Arguments:
	WindPtr wind	Wind to evaluate

Returns:

Description:
	Iterates over each wind cell, recording the total flux
	in each bin in the cell's arrays for use later on, as
	well as making a simple 'average path' calculation.

Notes:

History:
	26/2/15	-	Written by SWM
***********************************************************/
int
wind_paths_single_evaluate(Wind_Paths_Ptr paths)
{
	int		i         , j;

	paths->d_flux = 0.0;
	paths->d_path = 0.0;
	paths->i_num = 0.0;
	for (i = 0; i < NWAVE - 1; i++) {
		paths->ad_freq_flux[i] = 0.0;
		paths->ai_freq_num[i] = 0;

		for (j = 0; j < g_path_data->i_path_bins; j++) {
			paths->ad_freq_flux[i] +=
				paths->ad_freq_path_flux[i * g_path_data->i_path_bins + j];
			paths->ai_freq_num[i] +=
				paths->ai_freq_path_num[i * g_path_data->i_path_bins + j];
			paths->d_path +=
				paths->ad_freq_path_flux[i * g_path_data->i_path_bins +
					 j] * (g_path_data->ad_path_bin[j] +
					       g_path_data->ad_path_bin[j +
								  1]) / 2.0;
		}

		paths->d_flux += paths->ad_freq_flux[i];
		paths->i_num += paths->ai_freq_num[i];
	}
	if (paths->d_flux > 0.0)
		paths->d_path /= paths->d_flux;
	return (0);
}

int
wind_paths_evaluate(WindPtr wind)
{
	int		i;
	for (i = 0; i < NDIM * MDIM; i++) {
		if (wind[i].inwind >= 0)
			wind_paths_single_evaluate(wind[i].paths);
	}
	return (0);
}


/***********************************************************
Synopsis:
	wind_paths_point_index(int i, j, k, i_top)
		Given trz index in wind, returns vtk file index

Arguments:
	int i, j, k 	Theta, radius and height of cell
	int i_top		Whether point is + or -ive height

Returns:
	int n 			Index for vtk poly

Description:

Notes:

History:
	26/2/15	-	Written by SWM
***********************************************************/
int
wind_paths_point_index(int i, int j, int k, int i_top)
{
	int		n;
	n = i * 2 * (g_path_data->i_theta_res + 1) * MDIM +
		j * 2 * (g_path_data->i_theta_res + 1) + k * 2 + i_top;
	return (n);

}

int
wind_paths_output(WindPtr wind, char c_file_in[])
{
	FILE           *fopen(), *fptr;
	char		c_file    [LINELENGTH];
	int		i         , j, k, n, i_obs, i_cells, i_points;
	double		r_theta, r_x, r_y, r_err;
	PhotPtr		p_test = calloc(sizeof(p_dummy), 1);

	//Get output filename
		strcpy(c_file, c_file_in);	//Copy filename to new string
		strcat(c_file, ".wind_paths.vtk");

		if ((fptr = fopen(c_file, "w")) == NULL) {
		Error("wind_paths_output: Unable to open %s for writing\n", c_file);
		exit(0);
	}
	Log("Outputting wind path information to file '%s'.\n", c_file);


	i_cells = 2 * (NDIM - 1) * (MDIM - 1) * g_path_data->i_theta_res;
	i_points = 2 * NDIM * MDIM * (g_path_data->i_theta_res + 1);

	fprintf(fptr, "# vtk DataFile Version 2.0\n");
	fprintf(fptr, "Wind file data\nASCII\n");
	fprintf(fptr, "DATASET UNSTRUCTURED_GRID\n");
	fprintf(fptr, "POINTS %d float\n", i_points);
	for (i = 0; i < NDIM; i++) {
		for (j = 0; j < MDIM; j++) {
			wind_ij_to_n(i, j, &n);

			for (k = 0; k <= g_path_data->i_theta_res; k++) {
				r_theta = k * (PI / (double)g_path_data->i_theta_res);
				r_x = wind[n].x[0] * cos(r_theta);
				r_y = wind[n].x[0] * sin(r_theta);
				fprintf(fptr, "%10.5g %10.5g %10.5g\n", r_x, r_y,
					wind[n].x[2]);
				fprintf(fptr, "%10.5g %10.5g %10.5g\n", r_x, r_y,
					-wind[n].x[2]);
			}
		}
	}
	fprintf(fptr, "\n");

	fprintf(fptr, "CELLS %d %d\n", i_cells, 9 * i_cells);
	for (i = 0; i < NDIM - 1; i++) {
		for (j = 0; j < MDIM - 1; j++) {
			for (k = 0; k < g_path_data->i_theta_res; k++) {
				fprintf(fptr, "8 %d %d %d %d %d %d %d %d\n",
					wind_paths_point_index(i, j, k, 1),
				     wind_paths_point_index(i, j, k + 1, 1),
				 wind_paths_point_index(i, j + 1, k + 1, 1),
				     wind_paths_point_index(i, j + 1, k, 1),
				     wind_paths_point_index(i + 1, j, k, 1),
				 wind_paths_point_index(i + 1, j, k + 1, 1),
					wind_paths_point_index(i + 1, j + 1, k + 1, 1),
				wind_paths_point_index(i + 1, j + 1, k, 1));
				fprintf(fptr, "8 %d %d %d %d %d %d %d %d\n",
					wind_paths_point_index(i, j, k, 0),
				     wind_paths_point_index(i, j, k + 1, 0),
				 wind_paths_point_index(i, j + 1, k + 1, 0),
				     wind_paths_point_index(i, j + 1, k, 0),
				     wind_paths_point_index(i + 1, j, k, 0),
				 wind_paths_point_index(i + 1, j, k + 1, 0),
					wind_paths_point_index(i + 1, j + 1, k + 1, 0),
				wind_paths_point_index(i + 1, j + 1, k, 0));
			}
		}
	}

	fprintf(fptr, "CELL_TYPES %d\n", i_cells);
	for (i = 0; i < i_cells; i++)
		fprintf(fptr, "12\n");
	fprintf(fptr, "\n");

	fprintf(fptr, "CELL_DATA %d\n", i_cells);

	fprintf(fptr, "SCALARS phot_count float 1\n");
	fprintf(fptr, "LOOKUP_TABLE default\n");
	for (i = 0; i < NDIM - 1; i++) {
		for (j = 0; j < MDIM - 1; j++) {
			wind_ij_to_n(i, j, &n);
			for (k = 0; k < g_path_data->i_theta_res; k++) {
				fprintf(fptr, "%d\n", wind[n].paths->i_num);
				fprintf(fptr, "%d\n", wind[n].paths->i_num);
			}
		}
	}

	fprintf(fptr, "SCALARS path_errors float 1\n");
	fprintf(fptr, "LOOKUP_TABLE default\n");
	for (i = 0; i < NDIM - 1; i++) {
		for (j = 0; j < MDIM - 1; j++) {
			wind_ij_to_n(i, j, &n);
			for (k = 0; k < g_path_data->i_theta_res; k++) {
				if (wind[n].paths->i_num > 0) {
					r_err = sqrt((double)wind[n].paths->i_num) /
						(double)wind[n].paths->i_num;
					fprintf(fptr, "%g\n", r_err);
					fprintf(fptr, "%g\n", r_err);
				} else {
					fprintf(fptr, "-1\n");
					fprintf(fptr, "-1\n");
				}

			}
		}
	}

	fprintf(fptr, "SCALARS path_rel_diff_from_direct float 1\n");
	fprintf(fptr, "LOOKUP_TABLE default\n");
	for (i = 0; i < NDIM - 1; i++) {
		for (j = 0; j < MDIM - 1; j++) {
			wind_ij_to_n(i, j, &n);
			for (k = 0; k < g_path_data->i_theta_res; k++) {
				if (wind[n].paths->i_num > 0) {
					double		f_diff;
					f_diff = wind[n].paths->d_path;
					f_diff -= (sqrt(wind[n].xcen[0] * wind[n].xcen[0] +
					 wind[n].xcen[1] * wind[n].xcen[1] +
					  wind[n].xcen[2] * wind[n].xcen[2])
						   - geo.rstar);
					f_diff = fabs(f_diff);
					f_diff /= wind[n].paths->d_path;

					fprintf(fptr, "%g\n", f_diff);
					fprintf(fptr, "%g\n", f_diff);
				} else {
					fprintf(fptr, "-1\n");
					fprintf(fptr, "-1\n");
				}
			}
		}
	}

	fprintf(fptr, "SCALARS path_average float 1\n");
	fprintf(fptr, "LOOKUP_TABLE default\n");
	for (i = 0; i < NDIM - 1; i++) {
		for (j = 0; j < MDIM - 1; j++) {
			wind_ij_to_n(i, j, &n);
			for (k = 0; k < g_path_data->i_theta_res; k++) {
				if (wind[n].paths->i_num > 0) {
					fprintf(fptr, "%g\n", wind[n].paths->d_path);
					fprintf(fptr, "%g\n", wind[n].paths->d_path);
				} else {
					fprintf(fptr, "-1\n");
					fprintf(fptr, "-1\n");
				}
			}
		}
	}

	for (i_obs = 0; i_obs < g_path_data->i_obs; i_obs++) {
		fprintf(fptr, "SCALARS path_%s float 1\n", xxspec[MSPEC + i_obs].name);
		fprintf(fptr, "LOOKUP_TABLE default\n");
		stuff_v(xxspec[MSPEC + i_obs].lmn, p_test->lmn);

		for (i = 0; i < NDIM - 1; i++) {
			for (j = 0; j < MDIM - 1; j++) {
				wind_ij_to_n(i, j, &n);
				for (k = 0; k < g_path_data->i_theta_res; k++) {
					if (wind[n].paths->i_num > 0) {
						r_theta = k * (PI / (double)g_path_data->i_theta_res);
						p_test->x[0] = wind[n].xcen[0] * cos(r_theta);
						p_test->x[1] = wind[n].xcen[0] * sin(r_theta);
						p_test->x[2] = wind[n].xcen[2];
						fprintf(fptr, "%g\n",
						     wind[n].paths->d_path +
						 delay_to_observer(p_test));
						p_test->x[2] = -wind[n].xcen[2];
						fprintf(fptr, "%g\n",
						     wind[n].paths->d_path +
						 delay_to_observer(p_test));
					} else {
						fprintf(fptr, "-1\n");
						fprintf(fptr, "-1\n");
					}
				}
			}
		}
	}

	free(p_test);
	return (0);
}
