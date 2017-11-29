/* bb.c */
double planck(double t, double freqmin, double freqmax);
double get_rand_pow(double x1, double x2, double alpha);
double get_rand_exp(double alpha_min, double alpha_max);
double integ_planck_d(double alphamin, double alphamax);
int init_integ_planck_d(void);
double planck_d(double alpha);
double emittance_bb(double freqmin, double freqmax, double t);
double check_fmax(double fmin, double fmax, double temp);
/* get_atomicdata.c */
int get_atomic_data(char masterfile[]);
int index_lines(void);
int index_phot_top(void);
int index_inner_cross(void);
int index_collisions(void);
void indexx(int n, float arrin[], int indx[]);
int limit_lines(double freqmin, double freqmax);
int check_xsections(void);
/* python.c */
int main(int argc, char *argv[]);
/* photon2d.c */
int translate(WindPtr w, PhotPtr pp, double tau_scat, double *tau, int *nres);
int translate_in_space(PhotPtr pp);
double ds_to_wind(PhotPtr pp);
int translate_in_wind(WindPtr w, PhotPtr p, double tau_scat, double *tau, int *nres);
int walls(PhotPtr p, PhotPtr pold, double *normal);
/* photon_gen.c */
int define_phot(PhotPtr p, double f1, double f2, long nphot_tot, int ioniz_or_final, int iwind, int freq_sampling);
double populate_bands(double f1, double f2, int ioniz_or_final, int iwind, struct xbands *band);
int xdefine_phot(double f1, double f2, int ioniz_or_final, int iwind);
int xmake_phot(PhotPtr p, double f1, double f2, int ioniz_or_final, int iwind, double weight, int iphot_start, int nphotons);
int star_init(double freqmin, double freqmax, int ioniz_or_final, double *f);
int photo_gen_star(PhotPtr p, double r, double t, double weight, double f1, double f2, int spectype, int istart, int nphot);
double disk_init(double rmin, double rmax, double m, double mdot, double freqmin, double freqmax, int ioniz_or_final, double *ftot);
int photo_gen_disk(PhotPtr p, double weight, double f1, double f2, int spectype, int istart, int nphot);
int phot_gen_sum(char filename[], char mode[]);
double bl_init(double lum_bl, double t_bl, double freqmin, double freqmax, int ioniz_or_final, double *f);
int photon_checks(PhotPtr p, double freqmin, double freqmax, char *comment);
/* parse.c */
int parse_command_line(int argc, char *argv[]);
int help(void);
/* saha.c */
int nebular_concentrations(PlasmaPtr xplasma, int mode);
int concentrations(PlasmaPtr xplasma, int mode);
int saha(PlasmaPtr xplasma, double ne, double t);
int lucy(PlasmaPtr xplasma);
int lucy_mazzali1(double nh, double t_r, double t_e, double www, int nelem, double ne, double density[], double xne, double newden[]);
int fix_concentrations(PlasmaPtr xplasma, int mode);
double get_ne(double density[]);
/* spectra.c */
int spectrum_init(double f1, double f2, int nangle, double angle[], double phase[], int scat_select[], int top_bot_select[], int select_extract, double rho_select[], double z_select[], double az_select[], double r_select[]);
int spectrum_create(PhotPtr p, double f1, double f2, int nangle, int select_extract);
int spectrum_summary(char filename[], char mode[], int nspecmin, int nspecmax, int select_spectype, double renorm, int loglin, int iwind);
int spectrum_restart_renormalise(int nangle);
/* wind2d.c */
int define_wind(void);
int where_in_grid(int ndom, double x[]);
int vwind_xyz(int ndom, PhotPtr p, double v[]);
int wind_div_v(WindPtr w);
double rho(WindPtr w, double x[]);
int mdot_wind(WindPtr w, double z, double rmax);
int get_random_location(int n, double x[]);
int zero_scatters(void);
int check_corners_inwind(int n);
int check_grid(void);
/* wind.c */
int where_in_wind(double x[], int *ndomain);
int wind_check(WindPtr www, int n);
double model_velocity(int ndom, double x[], double v[]);
int model_vgrad(int ndom, double x[], double v_grad[][3]);
double model_rho(int ndom, double x[]);
/* vvector.c */
double dot(double a[], double b[]);
double length(double a[]);
int renorm(double a[], double scalar);
int cross(double a[], double b[], double c[]);
int vmove(double u[], double lmn[], double s, double result[]);
int vsub(double u[], double v[], double result[]);
int vadd(double u[], double v[], double result[]);
int stuff_v(double vin[], double vout[]);
double dot_tensor_vec(double tensor[3][3], double vin[3], double vout[3]);
int project_from_xyz_cyl(double a[], double b[], double result[]);
int project_from_cyl_xyz(double a[], double b[], double result[]);
int create_basis(double u[], double v[], struct basis *basis_new);
int project_from(struct basis *basis_from, double v_in[], double v_out[]);
int project_to(struct basis *basis_from, double v_in[], double v_out[]);
int reorient(struct basis *basis_from, struct basis *basis_to, double v_from[], double v_to[]);
/* debug.c */
int DebugStr(char *string);
/* recipes.c */
double qromb(double (*func)(double), double a, double b, double eps);
double trapzd(double (*func)(double), double a, double b, int n);
void polint(double xa[], double ya[], int n, double x, double *y, double *dy);
double zbrent(double (*func)(double), double x1, double x2, double tol);
void spline(double x[], double y[], int n, double yp1, double ypn, double y2[]);
void splint(double xa[], double ya[], double y2a[], int n, double x, double *y);
double expint(int n, double x);
double *vector(int i, int j);
void free_vector(double *a, int i, int j);
double rtsafe(void (*funcd)(double, double *, double *), double x1, double x2, double xacc);
double golden(double ax, double bx, double cx, double (*f)(double), double tol, double *xmin);
/* trans_phot.c */
int trans_phot(WindPtr w, PhotPtr p, int iextract);
int trans_phot_single(WindPtr w, PhotPtr p, int iextract);
/* phot_util.c */
int stuff_phot(PhotPtr pin, PhotPtr pout);
int move_phot(PhotPtr pp, double ds);
int comp_phot(PhotPtr p1, PhotPtr p2);
int phot_hist(PhotPtr p, int iswitch);
int phot_history_summarize(void);
double ds_to_cone(ConePtr cc, struct photon *p);
double ds_to_sphere(double r, struct photon *p);
double ds_to_sphere2(double x[], double r, struct photon *p);
int quadratic(double a, double b, double c, double r[]);
double ds_to_plane(struct plane *pl, struct photon *p);
double ds_to_closest_approach(double x[], struct photon *p, double *impact_parameter);
/* resonate.c */
double calculate_ds(WindPtr w, PhotPtr p, double tau_scat, double *tau, int *nres, double smax, int *istat);
int select_continuum_scattering_process(double kap_cont, double kap_es, double kap_ff, PlasmaPtr xplasma);
double kappa_bf(PlasmaPtr xplasma, double freq, int macro_all);
int kbf_need(double fmin, double fmax);
double sobolev(WindPtr one, double x[], double den_ion, struct lines *lptr, double dvds);
int doppler(PhotPtr pin, PhotPtr pout, double v[], int nres);
int scatter(PhotPtr p, int *nres, int *nnscat);
/* radiation.c */
int radiation(PhotPtr p, double ds);
double kappa_ff(PlasmaPtr xplasma, double freq);
double sigma_phot(struct topbase_phot *x_ptr, double freq);
double sigma_phot_verner(struct innershell *x_ptr, double freq);
double den_config(PlasmaPtr xplasma, int nconf);
double pop_kappa_ff_array(void);
int update_banded_estimators(PlasmaPtr xplasma, PhotPtr p, double ds, double w_ave);
double mean_intensity(PlasmaPtr xplasma, double freq, int mode);
/* init.c */
int init_log_and_windsave(int restart_stat);
double setup_dfudge(void);
int setup_windcone(void);
int setup_created_files(void);
/* wind_updates2d.c */
int wind_update(WindPtr (w));
int wind_rad_init(void);
int wind_ip(void);
/* windsave.c */
int wind_save(char filename[]);
int wind_read(char filename[]);
int wind_complete(WindPtr w);
int spec_save(char filename[]);
int spec_read(char filename[]);
/* extract.c */
int extract(WindPtr w, PhotPtr p, int itype);
int extract_one(WindPtr w, PhotPtr pp, int itype, int nspec);
/* cdf.c */
int cdf_gen_from_func(CdfPtr cdf, double (*func)(double), double xmin, double xmax, int njumps, double jump[]);
double gen_array_from_func(double (*func)(double), double xmin, double xmax, int pdfsteps);
int cdf_gen_from_array(CdfPtr cdf, double x[], double y[], int n_xy, double xmin, double xmax);
double cdf_get_rand(CdfPtr cdf);
int cdf_limit(CdfPtr cdf, double xmin, double xmax);
double cdf_get_rand_limit(CdfPtr cdf);
int cdf_to_file(CdfPtr cdf, char filename[]);
int cdf_check(CdfPtr cdf);
int calc_cdf_gradient(CdfPtr cdf);
int cdf_array_fixup(double *x, double *y, int n_xy);
/* roche.c */
int binary_basics(void);
double ds_to_roche_2(PhotPtr p);
int hit_secondary(PhotPtr p);
double pillbox(PhotPtr p, double *smin, double *smax);
double phi(double s);
double dphi_ds(double s);
double d2phi_ds2(double s);
double roche_width(double x);
double roche2_width_max(void);
void roche(double s, double *value, double *derivative);
void roche_deriv(double s, double *value, double *derivative);
/* random.c */
int randvec(double a[], double r);
int randvcos(double lmn[], double north[]);
double vcos(double x);
/* stellar_wind.c */
int get_stellar_wind_params(int ndom);
double stellar_velocity(int ndom, double x[], double v[]);
double stellar_rho(int ndom, double x[]);
/* homologous.c */
int get_homologous_params(int ndom);
double homologous_velocity(int ndom, double x[], double v[]);
double homologous_rho(int ndom, double x[]);
/* hydro_import.c */
int get_hydro_wind_params(int ndom);
int get_hydro(int ndom);
double hydro_velocity(int ndom, double x[], double v[]);
double hydro_rho(double x[]);
double hydro_temp(double x[]);
int rtheta_make_hydro_grid(WindPtr w, int ndom);
int rtheta_hydro_volumes(int ndom, WindPtr w);
int hydro_frac(double coord, double coord_array[], int imax, int *cell1, int *cell2, double *frac);
double hydro_interp_value(double array[], int im, int ii, int jm, int jj, double f1, double f2);
int hydro_restart(int ndom);
/* corona.c */
int get_corona_params(int ndom);
double corona_velocity(int ndom, double x[], double v[]);
double corona_rho(int ndom, double x[]);
/* knigge.c */
int get_knigge_wind_params(int ndom);
double kn_velocity(int ndom, double x[], double v[]);
double kn_rho(int ndom, double x[]);
double kn_vzero(double r);
double kn_wind_mdot_integral(double r);
double kn_rho_zero(int ndom, double r);
/* disk.c */
double tdisk(double m, double mdot, double r);
double teff(double t, double x);
double gdisk(double mass, double mdot, double rmin);
double geff(double g0, double x);
double vdisk(double x[], double v[]);
double zdisk(double r);
double ds_to_disk(struct photon *p, int miss_return);
void disk_deriv(double s, double *value, double *derivative);
int qdisk_init(void);
int qdisk_save(char *diskfile, double ztot);
int read_non_standard_disk_profile(char *tprofile);
/* lines.c */
double total_line_emission(WindPtr one, double f1, double f2);
double lum_lines(WindPtr one, int nmin, int nmax);
int lum_pdf(PlasmaPtr xplasma, double lumlines);
double q21(struct lines *line_ptr, double t);
double q12(struct lines *line_ptr, double t);
double a21(struct lines *line_ptr);
double two_level_atom(struct lines *line_ptr, PlasmaPtr xplasma, double *d1, double *d2);
double line_nsigma(struct lines *line_ptr, PlasmaPtr xplasma);
double scattering_fraction(struct lines *line_ptr, PlasmaPtr xplasma);
double p_escape(struct lines *line_ptr, PlasmaPtr xplasma);
double p_escape_from_tau(double tau);
int line_heat(PlasmaPtr xplasma, PhotPtr pp, int nres);
double upsilon(int n_coll, double u0);
/* continuum.c */
double one_continuum(int spectype, double t, double g, double freqmin, double freqmax);
double emittance_continuum(int spectype, double freqmin, double freqmax, double t, double g);
/* emission.c */
double wind_luminosity(double f1, double f2);
double total_emission(WindPtr one, double f1, double f2);
int photo_gen_wind(PhotPtr p, double weight, double freqmin, double freqmax, int photstart, int nphot);
double one_line(WindPtr one, double freqmin, double freqmax, int *nres);
double total_free(WindPtr one, double t_e, double f1, double f2);
double ff(WindPtr one, double t_e, double freq);
double one_ff(WindPtr one, double f1, double f2);
double gaunt_ff(double gsquared);
/* cooling.c */
double cooling(PlasmaPtr xxxplasma, double t);
double xtotal_emission(WindPtr one, double f1, double f2);
double adiabatic_cooling(WindPtr one, double t);
double wind_cooling(double f1, double f2);
/* recomb.c */
double fb_topbase_partial(double freq);
double integ_fb(double t, double f1, double f2, int nion, int fb_choice, int mode);
double total_fb(WindPtr one, double t, double f1, double f2, int fb_choice, int mode);
double one_fb(WindPtr one, double f1, double f2);
int num_recomb(PlasmaPtr xplasma, double t_e, int mode);
double fb(PlasmaPtr xplasma, double t, double freq, int ion_choice, int fb_choice);
int init_freebound(double t1, double t2, double f1, double f2);
double get_nrecomb(double t, int nion, int mode);
double get_fb(double t, int nion, int narray, int fb_choice, int mode);
double xinteg_fb(double t, double f1, double f2, int nion, int fb_choice);
double xinteg_inner_fb(double t, double f1, double f2, int nion, int fb_choice);
double total_rrate(int nion, double T);
double gs_rrate(int nion, double T);
int sort_and_compress(double *array_in, double *array_out, int npts);
int compare_doubles(const void *a, const void *b);
/* diag.c */
int open_diagfile(void);
int get_extra_diagnostics(void);
int save_photon_stats(WindPtr one, PhotPtr p, double ds, double w_ave);
/* sv.c */
int get_sv_wind_params(int ndom);
double sv_velocity(double x[], double v[], int ndom);
double sv_rho(int ndom, double x[]);
double sv_find_wind_rzero(int ndom, double p[]);
int sv_zero_init(double p[]);
double sv_zero_r(double r);
double sv_theta_wind(int ndom, double r);
double sv_wind_mdot_integral(double r);
/* ionization.c */
int ion_abundances(PlasmaPtr xplasma, int mode);
int convergence(PlasmaPtr xplasma);
int check_convergence(void);
int one_shot(PlasmaPtr xplasma, int mode);
double calc_te(PlasmaPtr xplasma, double tmin, double tmax);
double zero_emit(double t);
/* ispy.c */
int ispy_init(char filename[], int icycle);
int ispy_close(void);
int ispy(PhotPtr p, int n);
/* levels.c */
int levels(PlasmaPtr xplasma, int mode);
/* gradv.c */
double dvwind_ds(PhotPtr p);
int dvds_ave(void);
/* reposition.c */
int reposition(PhotPtr p);
/* anisowind.c */
int randwind(PhotPtr p, double lmn[3], double north[3]);
double vrandwind(double x);
double reweightwind(PhotPtr p);
int make_cdf_randwind(double tau);
int randwind_thermal_trapping(PhotPtr p, int *nnscat);
/* util.c */
int fraction(double value, double array[], int npts, int *ival, double *f, int mode);
int linterp(double x, double xarray[], double yarray[], int xdim, double *y, int mode);
int coord_fraction(int ndom, int ichoice, double x[], int ii[], double frac[], int *nelem);
int where_in_2dcell(int ichoice, double x[], int n, double *fx, double *fz);
int wind_n_to_ij(int ndom, int n, int *i, int *j);
int wind_ij_to_n(int ndom, int i, int j, int *n);
int wind_x_to_n(double x[], int *n);
/* density.c */
double get_ion_density(int ndom, double x[], int nion);
/* detail.c */
int detailed_balance(PlasmaPtr xplasma, int nelem, double newden[]);
int rebalance(double rates_up[], double rates_down[], double fraction[], int ntot);
int wind_update_after_detailed_balance(PlasmaPtr xplasma, int nelem, double newden[]);
/* bands.c */
int bands_init(int imode, struct xbands *band);
int freqs_init(double freqmin, double freqmax);
/* time.c */
double timer(void);
int get_time(char curtime[]);
/* matom.c */
int matom(PhotPtr p, int *nres, int *escape);
double b12(struct lines *line_ptr);
double alpha_sp(struct topbase_phot *cont_ptr, PlasmaPtr xplasma, int ichoice);
double alpha_sp_integrand(double freq);
int kpkt(PhotPtr p, int *nres, int *escape);
int fake_matom_bb(PhotPtr p, int *nres, int *escape);
int fake_matom_bf(PhotPtr p, int *nres, int *escape);
int emit_matom(WindPtr w, PhotPtr p, int *nres, int upper);
double matom_emit_in_line_prob(WindPtr one, struct lines *line_ptr_emit);
/* estimators.c */
int bf_estimators_increment(WindPtr one, PhotPtr p, double ds);
int bb_estimators_increment(WindPtr one, PhotPtr p, double tau_sobolev, double dvds, int nn);
int mc_estimator_normalise(int n);
double total_fb_matoms(PlasmaPtr xplasma, double t_e, double f1, double f2);
double total_bb_cooling(PlasmaPtr xplasma, double t_e);
double macro_bb_heating(PlasmaPtr xplasma, double t_e);
double macro_bf_heating(PlasmaPtr xplasma, double t_e);
int bb_simple_heat(PlasmaPtr xplasma, PhotPtr p, double tau_sobolev, double dvds, int nn);
int check_stimulated_recomb(PlasmaPtr xplasma);
int get_dilute_estimators(PlasmaPtr xplasma);
double get_gamma(struct topbase_phot *cont_ptr, PlasmaPtr xplasma);
double gamma_integrand(double freq);
double get_gamma_e(struct topbase_phot *cont_ptr, PlasmaPtr xplasma);
double gamma_e_integrand(double freq);
double get_alpha_st(struct topbase_phot *cont_ptr, PlasmaPtr xplasma);
double alpha_st_integrand(double freq);
double get_alpha_st_e(struct topbase_phot *cont_ptr, PlasmaPtr xplasma);
double alpha_st_e_integrand(double freq);
/* wind_sum.c */
int xtemp_rad(WindPtr w);
/* yso.c */
int get_yso_wind_params(int ndom);
double yso_velocity(int ndom, double x[], double v[]);
double yso_rho(int ndom, double x[]);
/* cylindrical.c */
double cylind_ds_in_cell(PhotPtr p);
int cylind_make_grid(int ndom, WindPtr w);
int cylind_wind_complete(int ndom, WindPtr w);
int cylind_volumes(int ndom, WindPtr w);
int cylind_where_in_grid(int ndom, double x[]);
int cylind_get_random_location(int n, double x[]);
int cylind_extend_density(int ndom, WindPtr w);
int cylind_is_cell_in_wind(int n);
/* rtheta.c */
double rtheta_ds_in_cell(PhotPtr p);
int rtheta_make_grid(WindPtr w, int ndom);
int rtheta_make_cones(int ndom, WindPtr w);
int rtheta_wind_complete(int ndom, WindPtr w);
int rtheta_volumes(int ndom, WindPtr w);
int rtheta_where_in_grid(int ndom, double x[]);
int rtheta_get_random_location(int n, double x[]);
int rtheta_extend_density(int ndom, WindPtr w);
int rtheta_is_cell_in_wind(int n);
/* spherical.c */
double spherical_ds_in_cell(PhotPtr p);
int spherical_make_grid(WindPtr w, int ndom);
int spherical_wind_complete(int ndom, WindPtr w);
int spherical_volumes(int ndom, WindPtr w);
int spherical_where_in_grid(int ndom, double x[]);
int spherical_get_random_location(int n, double x[]);
int spherical_extend_density(int ndom, WindPtr w);
int shell_make_grid(WindPtr w, int ndom);
/* cylind_var.c */
double cylvar_ds_in_cell(PhotPtr p);
int cylvar_make_grid(WindPtr w, int ndom);
int cylvar_wind_complete(int ndom, WindPtr w);
int cylvar_volumes(int ndom, WindPtr w);
int cylvar_where_in_grid(int ndom, double x[], int ichoice, double *fx, double *fz);
int cylvar_get_random_location(int n, double x[]);
int cylvar_extend_density(int ndom, WindPtr w);
int cylvar_coord_fraction(int ndom, int ichoice, double x[], int ii[], double frac[], int *nelem);
/* bilinear.c */
int bilin(double x[], double x00[], double x01[], double x10[], double x11[], double *f, double *g);
int xquadratic(double a, double b, double c, double r[]);
/* gridwind.c */
int create_maps(int ichoice);
int calloc_wind(int nelem);
int calloc_plasma(int nelem);
int check_plasma(PlasmaPtr xplasma, char message[]);
int calloc_macro(int nelem);
int calloc_estimators(int nelem);
int calloc_dyn_plasma(int nelem);
/* partition.c */
int partition_functions(PlasmaPtr xplasma, int mode);
int partition_functions_2(PlasmaPtr xplasma, int xnion, double temp, double weight);
/* signal.c */
int xsignal(char *root, char *format, ...);
int xsignal_rm(char *root);
int set_max_time(char *root, double t);
int check_time(char *root);
/* auger_ionization.c */
int auger_ionization(PlasmaPtr xplasma);
/* agn.c */
double agn_init(double r, double lum, double alpha, double freqmin, double freqmax, int ioniz_or_final, double *f);
double emittance_pow(double freqmin, double freqmax, double lum, double alpha);
double emittance_bpow(double freqmin, double freqmax, double lum, double alpha);
int photo_gen_agn(PhotPtr p, double r, double alpha, double weight, double f1, double f2, int spectype, int istart, int nphot);
/* shell_wind.c */
int get_shell_wind_params(int ndom);
/* compton.c */
double kappa_comp(PlasmaPtr xplasma, double freq);
double kappa_ind_comp(PlasmaPtr xplasma, double freq);
double total_comp(WindPtr one, double t_e);
double klein_nishina(double nu);
int compton_dir(PhotPtr p, PlasmaPtr xplasma);
double compton_func(double f);
double sigma_compton_partial(double f, double x);
double alpha(double nu);
double beta(double nu);
double comp_cool_integrand(double nu);
/* torus.c */
double ds_to_cylinder(double rho, struct photon *p);
/* zeta.c */
double compute_zeta(double temp, int nion, int mode);
/* dielectronic.c */
int compute_dr_coeffs(double temp);
double total_dr(WindPtr one, double t_e);
/* spectral_estimators.c */
int spectral_estimators(PlasmaPtr xplasma);
double pl_alpha_func_log(double alpha);
double pl_logmean(double alpha, double lnumin, double lnumax);
double pl_log_w(double j, double alpha, double lnumin, double lnumax);
double pl_log_stddev(double alpha, double lnumin, double lnumax);
double exp_temp_func(double exp_temp);
double exp_mean(double exp_temp, double numin, double numax);
double exp_w(double j, double exp_temp, double numin, double numax);
double exp_stddev(double exp_temp, double numin, double numax);
/* variable_temperature.c */
int variable_temperature(PlasmaPtr xplasma, int mode);
double pi_correct(double xtemp, int nion, PlasmaPtr xplasma, int mode);
double temp_func(double solv_temp);
/* matom_diag.c */
int matom_emiss_report(void);
/* direct_ion.c */
int compute_di_coeffs(double T);
int compute_qrecomb_coeffs(double T);
double total_di(WindPtr one, double t_e);
double q_ioniz_dere(int nion, double t_e);
double q_ioniz(struct topbase_phot *cont_ptr, double electron_temperature);
double q_recomb_dere(struct topbase_phot *cont_ptr, double electron_temperature);
double q_recomb(struct topbase_phot *cont_ptr, double electron_temperature);
/* pi_rates.c */
double calc_pi_rate(int nion, PlasmaPtr xplasma, int mode, int type);
double tb_planck1(double freq);
double tb_logpow1(double freq);
double tb_exp1(double freq);
/* matrix_ion.c */
int matrix_ion_populations(PlasmaPtr xplasma, int mode);
int populate_ion_rate_matrix(PlasmaPtr xplasma, double rate_matrix[nions][nions], double pi_rates[nions], double inner_rates[n_inner_tot], double rr_rates[nions], double b_temp[nions], double xne, int xelem[nions]);
int solve_matrix(double *a_data, double *b_data, int nrows, double *x, int nplasma);
/* para_update.c */
int communicate_estimators_para(void);
int gather_spectra_para(int nspec_helper, int nspecs);
int communicate_matom_estimators_para(void);
/* setup.c */
double get_stellar_params(void);
double get_disk_params(void);
int get_bl_and_agn_params(double lstar);
int get_meta_params(void);
int get_standard_care_factors(void);
/* setup_domains.c */
int get_domain_params(int ndom);
int get_wind_params(int ndom);
int get_line_transfer_mode(void);
/* photo_gen_matom.c */
double get_kpkt_f(void);
double get_matom_f(int mode);
int photo_gen_kpkt(PhotPtr p, double weight, int photstart, int nphot);
int photo_gen_matom(PhotPtr p, double weight, int photstart, int nphot);
/* macro_gov.c */
int macro_gov(PhotPtr p, int *nres, int matom_or_kpkt, int *which_out);
int macro_pops(PlasmaPtr xplasma, double xne);
/* windsave2table_sub.c */
int do_windsave2table(char *root);
int create_master_table(int ndom, char rootname[]);
int create_heat_table(int ndom, char rootname[]);
int create_ion_table(int ndom, char rootname[], int iz);
double *get_ion(int ndom, int element, int istate, int iswitch);
double *get_one(int ndom, char variable_name[]);
/* import.c */
int import_wind(int ndom);
int import_make_grid(WindPtr w, int ndom);
double import_velocity(int ndom, double *x, double *v);
int get_import_wind_params(int ndom);
double import_rho(int ndom, double *x);
/* import_spherical.c */
int import_1d(int ndom, char *filename);
int spherical_make_grid_import(WindPtr w, int ndom);
double velocity_1d(int ndom, double *x, double *v);
double rho_1d(int ndom, double *x);
/* import_cylindrical.c */
int import_cylindrical(int ndom, char *filename);
int cylindrical_make_grid_import(WindPtr w, int ndom);
double velocity_cylindrical(int ndom, double *x, double *v);
double rho_cylindrical(int ndom, double *x);
/* import_rtheta.c */
int import_rtheta(int ndom, char *filename);
int rtheta_make_grid_import(WindPtr w, int ndom);
double velocity_rtheta(int ndom, double *x, double *v);
double rho_rtheta(int ndom, double *x);
/* reverb.c */
double delay_to_observer(PhotPtr pp);
int delay_dump_prep(int restart_stat);
int delay_dump_finish(void);
int delay_dump_combine(int i_ranks);
int delay_dump(PhotPtr p, int np);
int delay_dump_single(PhotPtr pp, int i_spec);
/* paths.c */
Wind_Paths_Ptr wind_paths_constructor(WindPtr wind);
int reverb_init(WindPtr wind);
int wind_paths_init(WindPtr wind);
int line_paths_add_phot(WindPtr wind, PhotPtr pp, int *nres);
int wind_paths_add_phot(WindPtr wind, PhotPtr pp);
int simple_paths_gen_phot(PhotPtr pp);
double r_draw_from_path_histogram(Wind_Paths_Ptr PathPtr);
int wind_paths_gen_phot(WindPtr wind, PhotPtr pp);
int line_paths_gen_phot(WindPtr wind, PhotPtr pp, int nres);
int wind_paths_evaluate_single(Wind_Paths_Ptr paths);
int wind_paths_evaluate(WindPtr wind, int i_rank);
int wind_paths_dump(WindPtr wind, int rank_global);
int wind_paths_output_dump(WindPtr wind, int i_rank);
int wind_paths_point_index(int i, int j, int k, int i_top, DomainPtr dom);
int wind_paths_sphere_point_index(int i, int j, int k);
int wind_paths_output_vtk(WindPtr wind, int ndom);
/* setup2.c */
int init_geo(void);
int get_spectype(int yesno, char *question, int *spectype);
int init_advanced_modes(void);
int init_observers(void);
PhotPtr init_photons(void);
int init_ionization(void);
/* run.c */
int calculate_ionization(int restart_stat);
int make_spectra(int restart_stat);
/* brem.c */
double emittance_brem(double freqmin, double freqmax, double lum, double t);
double integ_brem(double freq);
double brem_d(double alpha);
double get_rand_brem(double freqmin, double freqmax);
/* search_light.c */
int search_light_init(void);
int photo_gen_search_light(PhotPtr p, double r, double alpha, double weight, double f1, double f2, int spectype, int istart, int nphot);
/* synonyms.c */
int check_synonyms(char new_question[], char old_question[]);
/* py_wind_sub.c */
int zoom(int direction);
int overview(WindPtr w, char rootname[]);
int position_summary(WindPtr w);
int abs_summary(WindPtr w, char rootname[], int ochoice);
int adiabatic_cooling_summary(WindPtr w, char rootname[], int ochoice);
int lum_summary(WindPtr w, char rootname[], int ochoice);
int photo_summary(WindPtr w, char rootname[], int ochoice);
int recomb_summary(WindPtr w, char rootname[], int ochoice);
int electron_summary(WindPtr w, char rootname[], int ochoice);
int rho_summary(WindPtr w, char rootname[], int ochoice);
int plasma_cell(WindPtr w, char rootname[], int ochoice);
int freq_summary(WindPtr w, char rootname[], int ochoice);
int nphot_summary(WindPtr w, char rootname[], int ochoice);
int temp_summary(WindPtr w, char rootname[], int ochoice);
int temp_rad(WindPtr w, char rootname[], int ochoice);
int weight_summary(WindPtr w, char rootname[], int ochoice);
int velocity_summary(WindPtr w, char rootname[], int ochoice);
int mo_summary(WindPtr w, char rootname[], int ochoice);
int vol_summary(WindPtr w, char rootname[], int ochoice);
int wind_element(WindPtr w);
int tau_h_summary(WindPtr w, char rootname[], int ochoice);
int coolheat_summary(WindPtr w, char rootname[], int ochoice);
int complete_file_summary(WindPtr w, char root[], int ochoice);
int wind_reg_summary(WindPtr w, char rootname[], int ochoice);
int dvds_summary(WindPtr w, char rootname[], int ochoice);
int inner_shell_summary(WindPtr w, char rootname[], int ochoice);
int IP_summary(WindPtr w, char rootname[], int ochoice);
int alpha_summary(WindPtr w, char rootname[], int ochoice);
int J_summary(WindPtr w, char rootname[], int ochoice);
int J_scat_summary(WindPtr w, char rootname[], int ochoice);
int phot_split(WindPtr w, char rootname[], int ochoice);
int thompson(WindPtr w, char rootname[], int ochoice);
int nscat_split(WindPtr w, char rootname[], int ochoice);
int convergence_summary(WindPtr w, char rootname[], int ochoice);
int convergence_all(WindPtr w, char rootname[], int ochoice);
int model_bands(WindPtr w, char rootname[], int ochoice);
int heatcool_summary(WindPtr w, char rootname[], int ochoice);
int complete_physical_summary(WindPtr w, char rootname[], int ochoice);
int complete_ion_summary(WindPtr w, char rootname[], int ochoice);
double get_density_or_frac(PlasmaPtr xplasma, int element, int istate, int frac_choice);
int find_ion(int element, int istate);
int find_element(int element);
int get_los_dvds(WindPtr w, char rootname[], int ochoice);
int grid_summary(WindPtr w, char rootname[], int ochoice);
/* py_wind_ion.c */
int ion_summary(WindPtr w, int element, int istate, int iswitch, char rootname[], int ochoice);
int tau_ave_summary(WindPtr w, int element, int istate, double freq, char rootname[], int ochoice);
int line_summary(WindPtr w, char rootname[], int ochoice);
int total_emission_summary(WindPtr w, char rootname[], int ochoice);
int modify_te(WindPtr w, char rootname[], int ochoice);
int partial_measure_summary(WindPtr w, int element, int istate, char rootname[], int ochoice);
int collision_summary(WindPtr w, char rootname[], int ochoice);
/* py_wind_write.c */
int write_array(char filename[], int choice);
int display(char name[]);
/* py_wind_macro.c */
int xadiabatic_cooling_summary(WindPtr w, char rootname[], int ochoice);
int macro_summary(WindPtr w, char rootname[], int ochoice);
int ion_overview(int icell);
int config_overview(int n, int icell);
int depcoef_overview(int icell);
int copy_plasma(PlasmaPtr x1, PlasmaPtr x2);
int dealloc_copied_plasma(PlasmaPtr xcopy);
int depcoef_overview_specific(int version, int nconfig, WindPtr w, char rootname[], int ochoice);
int level_popsoverview(int nplasma, WindPtr w, char rootname[], int ochoice);
int level_emissoverview(int nlev, WindPtr w, char rootname[], int ochoice);
int level_escapeoverview(int nlev, WindPtr w, char rootname[], int ochoice);
int level_tauoverview(int nlev, WindPtr w, char rootname[], int ochoice);
/* py_wind.c */
int main(int argc, char *argv[]);
int one_choice(int choice, char *root, int ochoice);
int py_wind_help(void);
