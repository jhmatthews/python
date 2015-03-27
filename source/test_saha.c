/***********************************************************
                                       Space Telescope Science Institute

 Synopsis:
	Saha_inv is a quick program to generate lots of saha abundances for a range of temperatures and electron densities.
 
Arguments:		


Returns:
 
Description:	
	
Notes:

History:
	 
	0212 - nsh wrote it
 	
 	Look in Readme.c for more text concerning the early history of the program.

**************************************************************/



#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "atomic.h"
struct topbase_phot *xtop;
double g1, g2, xne, xip;
#define SAHA 4.82907e15
double xip, xxxne, qromb_temp;
#include "python.h"
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
  double alpha_agn, distance, agn_ip, lum_agn, const_agn, IPstart, IPstop, IP;
  int lstart, lmax;
  double xweight;
  double temp_func (), trr_rate;
  FILE *fp_h, *fp1_h;
  FILE *fp_he, *fp1_he;
  FILE *fp_c, *fp1_c;
  FILE *fp_n, *fp1_n;
  FILE *fp_o, *fp1_o;
  FILE *fp_fe;
  FILE *c_trr;
  FILE *c_prr;
  PlasmaPtr xplasma;
  partition = xplasma->partition;	/* Set the partition function array to that held for the cell */
  temp_start = 400;
  ntemp = 101;
  tempdiv = 100.0;
  nne = 1;
  nh = 1e10;
  weight = 0.0001;
  xne = nh;

  strcpy (geo.atomic_filename, "data/standard73");
  printf ("atomic_filename=%s\n", geo.atomic_filename);
  get_atomic_data (geo.atomic_filename);


/*xip=8.801388e-11 ;
xxxne=1e10;
/*for (i=3000;i<5001;i++)
	{
	temp=pow(10,(i/1000.0));
	test=temp_func(temp);
   	printf ("HHHHH i=%i,temp=%e,tempfunc=%e\n",i,temp,test);
	}*/
/*xip=5.139*EV2ERGS;
for (i=-6;i<21;i++)
	{
	xxxne=pow(10,i);
   	test=zbrent(temp_func,1e2,1e8,10);
	printf ("HHHHH xxxne=%e, temp=%e\n",xxxne,test);
	}
xip=13.6*EV2ERGS;
for (i=-6;i<21;i++)
	{
	xxxne=pow(10,i);
   	test=zbrent(temp_func,1e2,1e8,10);
	printf ("HHHHH xxxne=%e, temp=%e\n",xxxne,test);
	}
*/
  p = (PhotPtr) calloc (sizeof (p_dummy), NPHOT);
  NDIM = ndim = 1;
  MDIM = mdim = 1;
  NDIM2 = NDIM * MDIM;

  calloc_wind (NDIM2);
  w = wmain;

  NPLASMA = 1;
  calloc_plasma (NPLASMA);
  calloc_dyn_plasma (NPLASMA);

  plasmamain[0].rho = nh / rho2nh;
  plasmamain[0].ne = nh;

//  xplasma=plasmamain;
/*
fp_o=fopen("hydrogen_details.dat","w");
	xip=13.6*EV2ERGS;
	for (i=-6;i<21;i++)
		{
		xxxne=pow(10,i);
   		test=zbrent(temp_func,1e2,1e8,10);
		plasmamain[0].t_e=test;
		plasmamain[0].t_r=test;
		plasmamain[0].w=1;
		nion=1;
		partition_functions_2 (&plasmamain[0], nion, test);
		printf ("p1=%f p2=%f\n",plasmamain[0].partition[nion-1],plasmamain[0].partition[nion]);
		xsaha = SAHA * pow (test, 1.5);
	  	b = xsaha * plasmamain[0].partition[nion]
	    	* exp (-ion[nion - 1].ip / (BOLTZMANN * test)) / (xxxne *
							   plasmamain[0].partition[nion-1]);
      		ntmin = ion[nion-1].ntop_ground;
     		n = ntmin;
      		xtop = &phot_top[n];
      		fthresh = xtop->freq[0];
      		fmax = xtop->freq[xtop->np - 1];
		printf ("fthresh=%e, fmax=%e\n",fthresh,fmax);
		qromb_temp=test;  //The numerator is for the actual radiation temperature
		printf ("temp=%e\n",temp);
      		integral=qromb (tb_planck1, fthresh, fmax, 1.e-4);
		fprintf (fp_o,"xxxne %e temp %e saha %e integral %e\n",xxxne,test,b,integral);
		}
fclose(fp_o);

fp_o=fopen("Hydrogen_integral.dat","w");
	for (i=20;i<91;i++)
		{
		qromb_temp=pow(10,i/10.0);
      		integral=qromb (tb_planck1, fthresh, fmax, 1.e-4);
		fprintf (fp_o,"temp %e integral %e\n",temp,integral);
		}
fclose(fp_o);


  fp_o=fopen("atomic_details.dat","w");  
  for (j=0;j<nelements; j++)
	{
	fprintf(fp_o,"Element %i,%s,%i,%i,%i,%e\n",j,ele[j].name,ele[j].z,ele[j].firstion,ele[j].nions,ele[j].abun);
	} 
 for (j=0;j<nions; j++)
	{
	fprintf(fp_o,"Ion %i,%i,%i,%e,%f,%i,%i,%i,%i,%i,%i\n",j,ion[j].z,ion[j].istate,ion[j].ip,ion[j].g,ion[j].drflag,drecomb[ion[j].nxdrecomb].type,ion[j].bad_gs_rr_t_flag,ion[j].bad_gs_rr_r_flag,ion[j].total_rrflag,total_rr[ion[j].nxtotalrr].type);
	} 
  for (j=0;j<nlevels; j++)
	{
	fprintf(fp_o,"Level %i,%i,%i,%i,%i,%i,%i,%i,%i,%i,%i,%i,%f,%f,%e,%e,%i,%i,%i\n,",j,config[j].z,config[j].istate,config[j].nion,config[j].nden,config[j].isp,config[j].ilv,config[j].macro_info,config[j].n_bbu_jump,config[j].n_bbd_jump,config[j].n_bfu_jump,config[j].n_bfd_jump,config[j].g,config[j].q_num,config[j].ex,config[j].rad_rate,config[j].bbu_indx_first,config[j].bfu_indx_first,config[j].bfd_indx_first);

	}
 for (j=0;j<nxphot; j++)
	{
	fprintf(fp_o,"VFKY %i,%i,%i,%i,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e\n,",j,xphot[j].nion,xphot[j].z,xphot[j].istate,xphot[j].freq_t,xphot[j].freq_max,xphot[j].freq0,xphot[j].sigma,xphot[j].ya,xphot[j].p,xphot[j].yw,xphot[j].y0,xphot[j].y1,xphot[j].f_last,xphot[j].sigma_last);

	}
 for (j=0;j<ntop_phot; j++)
	{
	fprintf(fp_o,"topbase %i,%i,%i,%i,%i,%i,%i,%i,%i,%i,%i,%e,%e,%e,%e\n",j,phot_top[j].nlev,phot_top[j].uplev,phot_top[j].nion,phot_top[j].z,phot_top[j].istate,phot_top[j].np,phot_top[j].nlast,phot_top[j].macro_info,phot_top[j].down_index,phot_top[j].up_index,phot_top[j].freq[0],phot_top[j].x[0],phot_top[j].f,phot_top[j].sigma);

	}

  fclose(fp_o);
/*
  fp_o=fopen("O4details.dat","w");
  for (j=20;j<29;j++ ) 
	{		
	fprintf(fp_o," %6.3e",ion[j].ip);
	}
	fprintf(fp_o,"\n");
  for (j=20;j<29;j++ ) 
	{		
	fprintf(fp_o," %6.3f",ion[j].ip/EV2ERGS);
	}
	fprintf(fp_o,"\n");
	fprintf(fp_o,"Number of topbase cross sections");
  for (j=20;j<29;j++ ) 
	{		
	fprintf(fp_o," %i",ion[j].ntop);
	}
	fprintf(fp_o,"\n");
	fprintf(fp_o,"First topbase cross section");
  for (j=20;j<29;j++ ) 
	{		
	fprintf(fp_o," %i",ion[j].ntop_first);
	}
	fprintf(fp_o,"\n");
	fprintf(fp_o,"First nlte_level");
  for (j=20;j<29;j++ ) 
	{		
	fprintf(fp_o," %i",ion[j].first_nlte_level);
	}
	fprintf(fp_o,"\n");
	fprintf(fp_o,"Number of nlte levels");
  for (j=20;j<29;j++ ) 
	{		
	fprintf(fp_o," %i",ion[j].nlte);
	}
	fprintf(fp_o,"\n");
	fprintf(fp_o,"first nlte level threshold");
  for (j=20;j<29;j++ ) 
	{		
	fprintf(fp_o," %e",(config[ion[j].first_nlte_level].ex)/EV2ERGS);
	}
	fprintf(fp_o,"\n");
	fprintf(fp_o,"ground state topbase");
  for (j=20;j<29;j++ ) 
	{		
	fprintf(fp_o," %i",ion[j].ntop_ground);
	}
	fprintf(fp_o,"\n");
	fprintf(fp_o,"Topbase threshold");
  for (j=20;j<29;j++ ) 
	{		
	fprintf(fp_o," %6.3e",phot_top[ion[j].ntop_ground].freq[0]);
	}
	fprintf(fp_o,"\n");
	fprintf(fp_o,"Verner  threshold");
  for (j=20;j<29;j++ ) 
	{		
	fprintf(fp_o," %6.3e",xphot[ion[j].nxphot].freq_t);
	}
	fprintf(fp_o,"\n");
	fprintf(fp_o,"From xi threshold");
  for (j=20;j<29;j++ ) 
	{		
	fprintf(fp_o," %6.3e",(ion[j].ip/EV2ERGS)/HEV);
	}
	fprintf(fp_o,"\n");

   fprintf(fp_o,"%i\n",ntemp);

  plasmamain[0].t_e=t_e=t_r=pow(10,5);
  plasmamain[0].w=1.0;
  nion=23;



  temp1=temp_start;
	for (i=0;i<ntemp;i++)
		{
		temp1++;
		temp=pow(10,temp1/tempdiv);
  		partition_functions (&plasmamain[0],0);
		xsaha = SAHA * pow (temp, 1.5);
	  	b = xsaha * plasmamain[0].partition[nion]
	    	* exp (-ion[nion - 1].ip / (BOLTZMANN * temp)) / (xne *
			plasmamain[0].partition[nion-1]);
                recomb_fudge = sqrt (t_e / temp);
		pi_fudge = bb_correct_2 (temp, t_r, plasmamain[0].w, nion);
		gs_fudge = compute_zeta (t_e, nion-1, 1e14, 1e18, 1); 
		fprintf(fp_o,"%e %e %e %e %e \n",temp,b,recomb_fudge,pi_fudge,gs_fudge);
}

fclose(fp_o);
*/
/*
  fp_h=fopen("hydrogen_saha.out","w");
  fp1_h=fopen("hydrogen_partition.out","w");
  fp_he=fopen("helium_saha.out","w");
  fp1_he=fopen("helium_partition.out","w");
  fp_c=fopen("carbon_saha.out","w");
  fp1_c=fopen("carbon_partition.out","w");
  fp_n=fopen("nitrogen_saha.out","w");
  fp1_n=fopen("nitrogen_partition.out","w");
  fp_o=fopen("oxygen_saha.out","w");
  fp1_o=fopen("oxygen_partition.out","w");
	fprintf(fp_h,"%i %i\n",ntemp,nne);
	fprintf(fp_he,"%i %i\n",ntemp,nne);
	fprintf(fp_c,"%i %i\n",ntemp,nne);
	fprintf(fp_n,"%i %i\n",ntemp,nne);
	fprintf(fp_o,"%i %i\n",ntemp,nne);



	temp1=temp_start;
	for (i=0;i<ntemp;i++)
		{
		temp1++;
		temp=pow(10,temp1/tempdiv);
		temp=30000;
		plasmamain[0].t_e=temp;
  		plasmamain[0].ne=nh;
		plasmamain[0].t_r=temp;
		plasmamain[0].w=1;
		cardona_part_func (&plasmamain[0]);
  		saha (&plasmamain[0],  nh, temp);
		fprintf(fp_h,"%6.3e %6.3e %6.3e %6.3e\n",temp,plasmamain[0].ne,plasmamain[0].density[0],plasmamain[0].density[1]);
		fprintf(fp1_h,"%6.3e %6.3e %6.3e %6.3e\n",temp,plasmamain[0].ne,plasmamain[0].partition[0],plasmamain[0].partition[1]);
		fprintf(fp_he,"%e %e",temp,plasmamain[0].ne);		
		fprintf(fp1_he,"%e %e",temp,plasmamain[0].ne);		
		for (j=2;j<5;j++ ) 
			{		
			fprintf(fp_he," %6.3e",plasmamain[0].density[j]);
			fprintf(fp1_he," %6.3e",plasmamain[0].partition[j]);
			}
		fprintf(fp_he,"\n");
		fprintf(fp1_he,"\n");
		fprintf(fp_c,"%e %e",temp,plasmamain[0].ne);
		fprintf(fp1_c,"%e %e",temp,plasmamain[0].ne);				
		for (j=5;j<12;j++ ) 
			{		
			fprintf(fp_c," %6.3e",plasmamain[0].density[j]);
			fprintf(fp1_c," %6.3e",plasmamain[0].partition[j]);
			}
		fprintf(fp_c,"\n");
		fprintf(fp1_c,"\n");
		fprintf(fp_n,"%e %e",temp,plasmamain[0].ne);
		fprintf(fp1_n,"%e %e",temp,plasmamain[0].ne);				
		for (j=12;j<20;j++ ) 
			{		
			fprintf(fp_n," %6.3e",plasmamain[0].density[j]);
			fprintf(fp1_n," %6.3e",plasmamain[0].partition[j]);
			}
		fprintf(fp_n,"\n");
		fprintf(fp1_n,"\n");
		fprintf(fp_o,"%e %e",temp,plasmamain[0].ne);
		fprintf(fp1_o,"%e %e",temp,plasmamain[0].ne);		
		for (j=20;j<29;j++ ) 
			{		
			fprintf(fp_o," %6.3e",plasmamain[0].density[j]);
			fprintf(fp1_o," %6.3e",plasmamain[0].partition[j]);
			}
		fprintf(fp_o,"\n");
		fprintf(fp1_o,"\n");
		}

	
  fclose(fp_h);
  fclose(fp1_h);
  fclose(fp_he);
  fclose(fp1_he);
  fclose(fp_c);
  fclose(fp1_c);
  fclose(fp_n);
  fclose(fp1_n);
  fclose(fp_o);
  fclose(fp1_o);


/* fp_h=fopen("hydrogen_nc.out","w");
  fp_he=fopen("helium_nc.out","w");
  fp_c=fopen("carbon_nc.out","w");
  fp1_o=fopen("carbon_partition.out","w");
  fp_n=fopen("nitrogen_nc.out","w");
  fp_o=fopen("oxygen_nc.out","w");
	fprintf(fp_h,"%i %i\n",ntemp,nne);
	fprintf(fp_he,"%i %i\n",ntemp,nne);
	fprintf(fp_c,"%i %i\n",ntemp,nne);
	fprintf(fp_n,"%i %i\n",ntemp,nne);
	fprintf(fp_o,"%i %i\n",ntemp,nne);
 


	temp1=temp_start;
	for (i=0;i<ntemp;i++)
		{
		temp1++;
		temp=pow(10,temp1/tempdiv);
		temp=20000;
		plasmamain[0].t_e=temp;
		plasmamain[0].t_r=temp;
		plasmamain[0].w=weight;
  		nebular_concentrations (&plasmamain[0], 0);
		fprintf(fp_h,"%6.3e %6.3e %6.3e %6.3e\n",temp,plasmamain[0].ne,plasmamain[0].density[0],plasmamain[0].density[1]);
		fprintf(fp_he,"%e %e",temp,plasmamain[0].ne);		
		for (j=2;j<5;j++ ) 
			{		
			fprintf(fp_he," %6.3e",plasmamain[0].density[j]);
			}
		fprintf(fp_he,"\n");
		fprintf(fp_c,"%e %e",temp,plasmamain[0].ne);
		fprintf(fp1_o,"%e %e",temp,plasmamain[0].ne);		
		for (j=5;j<12;j++ ) 
			{		
			fprintf(fp_c," %6.3e",plasmamain[0].density[j]);
			fprintf(fp1_o," %6.3e",plasmamain[0].partition[j]);
			}
		fprintf(fp_c,"\n");
		fprintf(fp1_o,"\n");
		fprintf(fp_n,"%e %e",temp,plasmamain[0].ne);		
		for (j=12;j<20;j++ ) 
			{		
			fprintf(fp_n," %6.3e",plasmamain[0].density[j]);
			}
		fprintf(fp_n,"\n");
		fprintf(fp_o,"%e %e",temp,plasmamain[0].ne);		
		for (j=20;j<29;j++ ) 
			{		
			fprintf(fp_o," %6.3e",plasmamain[0].density[j]);
			}
		fprintf(fp_o,"\n");
		}

  fclose(fp_h);
  fclose(fp_he);
  fclose(fp_c);
  fclose(fp_n);
  fclose(fp_o);
*/

  fp_h=fopen("hydrogen_lm.out","w");
  fp_he=fopen("helium_lm.out","w");
  fp_c=fopen("carbon_lm.out","w");
  fp_n=fopen("nitrogen_lm.out","w");
  fp_o=fopen("oxygen_lm.out","w");
  fp_fe=fopen("iron_lm.out","w");
//  c_trr=fopen("carbon_trr.out","w");
	fprintf(fp_h,"%i %i\n",ntemp,nne);
	fprintf(fp_he,"%i %i\n",ntemp,nne);
	fprintf(fp_c,"%i %i\n",ntemp,nne);
	fprintf(fp_n,"%i %i\n",ntemp,nne);
	fprintf(fp_o,"%i %i\n",ntemp,nne);
	fprintf(fp_fe,"%i %i\n",ntemp,nne);
//	fprintf(c_trr,"%i %i\n",ntemp,nne);

	temp1=temp_start;
	for (i=0;i<ntemp;i++)
		{
		temp1=temp1+1;
		temp=pow(10,temp1/tempdiv);
		plasmamain[0].t_e=temp*0.9;
		plasmamain[0].t_r=temp;
		plasmamain[0].w=weight;
//		printf ("WEIGHT=%f\n",weight);
  		nebular_concentrations (&plasmamain[0], 2);
		fprintf(fp_h,"%6.3e %6.3e %6.3e %6.3e\n",temp,plasmamain[0].ne,plasmamain[0].density[0],plasmamain[0].density[1]);
		fprintf(fp_he,"%e %e",temp,plasmamain[0].ne);		
		for (j=2;j<5;j++ ) 
			{		
			fprintf(fp_he," %6.3e",plasmamain[0].density[j]);
			}
		fprintf(fp_he,"\n");
		fprintf(fp_c,"%e %e",temp,plasmamain[0].ne);
//		fprintf(c_trr,"%e %e",temp,plasmamain[0].ne);
		for (j=5;j<12;j++ ) 
			{
//			trr_rate=badnell_gs_rr (j,temp);
//			fprintf(c_trr," %6.3e",trr_rate);		
			fprintf(fp_c," %6.3e",plasmamain[0].density[j]);
			}
		fprintf(fp_c,"\n");
//		fprintf(c_trr," \n");		
		fprintf(fp_n,"%e %e",temp,plasmamain[0].ne);		
		for (j=12;j<20;j++ ) 
			{		
			fprintf(fp_n," %6.3e",plasmamain[0].density[j]);
			}
		fprintf(fp_n,"\n");
		fprintf(fp_o,"%e %e",temp,plasmamain[0].ne);
		for (j=20;j<29;j++ ) 
			{		
			fprintf(fp_o," %6.3e",plasmamain[0].density[j]);
			}
		fprintf(fp_o,"\n");
		fprintf(fp_fe,"%e %e",temp,plasmamain[0].ne);		
		for (j=151;j<178;j++ ) 
			{		
			fprintf(fp_fe," %6.3e",plasmamain[0].density[j]);
			}
		fprintf(fp_fe,"\n");
		}

  fclose(fp_h);
  fclose(fp_he);
  fclose(fp_c);
  fclose(fp_n);
  fclose(fp_o);
  fclose(fp_fe);
//  fclose(c_trr);

/*


  fp_h=fopen("hydrogen_vt.out","w");
  fp_he=fopen("helium_vt.out","w");
  fp_c=fopen("carbon_vt.out","w");
  fp_n=fopen("nitrogen_vt.out","w");
  fp_o=fopen("oxygen_vt.out","w");
  fp_fe=fopen("iron_vt.out","w");
	fprintf(fp_h,"%i %i\n",ntemp,nne);
	fprintf(fp_he,"%i %i\n",ntemp,nne);
	fprintf(fp_c,"%i %i\n",ntemp,nne);
	fprintf(fp_n,"%i %i\n",ntemp,nne);
	fprintf(fp_o,"%i %i\n",ntemp,nne);
	fprintf(fp_fe,"%i %i\n",ntemp,nne);

	temp1=temp_start;
	for (i=0;i<ntemp;i++)
		{
		temp1=temp1+1;
		temp=pow(10,temp1/tempdiv);
		printf ("Temperature = %e\n",temp);
		plasmamain[0].t_e=temp*0.9;
		plasmamain[0].t_r=temp;
		plasmamain[0].w=weight;
		printf("Going to the new code\n");
  		nebular_concentrations (&plasmamain[0], 6);
		fprintf(fp_h,"%6.3e %6.3e %6.3e %6.3e\n",temp,plasmamain[0].ne,plasmamain[0].density[0],plasmamain[0].density[1]);
		fprintf(fp_he,"%e %e",temp,plasmamain[0].ne);		
		for (j=2;j<5;j++ ) 
			{		
			fprintf(fp_he," %6.3e",plasmamain[0].density[j]);
			}
		fprintf(fp_he,"\n");
		fprintf(fp_c,"%e %e",temp,plasmamain[0].ne);		
		for (j=5;j<12;j++ ) 
			{		
			fprintf(fp_c," %6.3e",plasmamain[0].density[j]);
			}
		fprintf(fp_c,"\n");
		fprintf(fp_n,"%e %e",temp,plasmamain[0].ne);		
		for (j=12;j<20;j++ ) 
			{		
			fprintf(fp_n," %6.3e",plasmamain[0].density[j]);
			}
		fprintf(fp_n,"\n");
		fprintf(fp_o,"%e %e",temp,plasmamain[0].ne);
		for (j=20;j<29;j++ ) 
			{		
			fprintf(fp_o," %6.3e",plasmamain[0].density[j]);
			}
		fprintf(fp_o,"\n");
		fprintf(fp_fe,"%e %e",temp,plasmamain[0].ne);		
		for (j=151;j<178;j++ ) 
			{		
			fprintf(fp_fe," %6.3e",plasmamain[0].density[j]);
			}
		fprintf(fp_fe,"\n");
		}

  fclose(fp_h);
  fclose(fp_he);
  fclose(fp_c);
  fclose(fp_n);
  fclose(fp_o);
  fclose(fp_fe);


/*
alpha_agn=-1.5;
lstart=250;
lmax=450;
distance=1e15;
nh=1e10;

  fp_h=fopen("hydrogen_pl.out","w");
  fp_he=fopen("helium_pl.out","w");
  fp_c=fopen("carbon_pl.out","w");
  fp_n=fopen("nitrogen_pl.out","w");
  fp_o=fopen("oxygen_pl.out","w");
  fp_fe=fopen("iron_pl.out","w");
	fprintf(fp_h,"%i\n",lmax-lstart);
	fprintf(fp_he,"%i\n",lmax-lstart);
	fprintf(fp_c,"%i\n",lmax-lstart);
	fprintf(fp_n,"%i\n",lmax-lstart);
	fprintf(fp_o,"%i\n",lmax-lstart);
	fprintf(fp_fe,"%i\n",lmax-lstart);
 
		geo.nxfreq=1;
		geo.xfreq[0]=1e14;
		geo.xfreq[1]=1e18;
		plasmamain[0].sim_alpha[0]=alpha_agn;
		xband.nbands=1;
    		xband.f1[0]=1e14;
		xband.f2[0]=1e18;
		plasmamain[0].t_e=1e6;
		plasmamain[0].t_r=1e6;
  		plasmamain[0].rho=nh/rho2nh;
		plasmamain[0].ne=nh;

	for (i=lstart;i<lmax;i++)
		{
		lum_agn=pow(10,(i/10.0));
		const_agn = lum_agn / (((pow (10000/HEV, alpha_agn + 1.)) - pow (2000/HEV, alpha_agn + 1.0)) /
	   	(alpha_agn + 1.0));
		printf ("const_agn=%e\n",const_agn);
		agn_ip=const_agn*(((pow (50000/HEV, alpha_agn + 1.0)) - pow (100/HEV,alpha_agn + 1.0)) /  				(alpha_agn + 1.0));
		agn_ip /= (distance*distance);
		agn_ip /= nh;
		Log("i=%i,Ionisation Parameter=%e\n",i,(agn_ip));
		
	
		plasmamain[0].sim_w[0]=const_agn/((4.*PI)*(4.*PI*distance*distance));
		printf("weight=%e\n",plasmamain[0].sim_w[0]);
  		nebular_concentrations (&plasmamain[0], 7);
		fprintf(fp_h,"%6.3e %6.3e %6.3e %6.3e\n",agn_ip,plasmamain[0].ne,plasmamain[0].density[0],plasmamain[0].density[1]);
		fprintf(fp_he,"%e %e",agn_ip,plasmamain[0].ne);		
		for (j=2;j<5;j++ ) 
			{		
			fprintf(fp_he," %6.3e",plasmamain[0].density[j]);
			}
		fprintf(fp_he,"\n");
		fprintf(fp_c,"%e %e",agn_ip,plasmamain[0].ne);		
		for (j=5;j<12;j++ ) 
			{		
			fprintf(fp_c," %6.3e",plasmamain[0].density[j]);
			}
		fprintf(fp_c,"\n");
		fprintf(fp_n,"%e %e",agn_ip,plasmamain[0].ne);		
		for (j=12;j<20;j++ ) 
			{		
			fprintf(fp_n," %6.3e",plasmamain[0].density[j]);
			}
		fprintf(fp_n,"\n");
		fprintf(fp_o,"%e %e",agn_ip,plasmamain[0].ne);		
		for (j=20;j<29;j++ ) 
			{		
			fprintf(fp_o," %6.3e",plasmamain[0].density[j]);
			}
		fprintf(fp_o,"\n");
		fprintf(fp_fe,"%e %e",agn_ip,plasmamain[0].ne);		
		for (j=151;j<178;j++ ) 
			{		
			fprintf(fp_fe," %6.3e",plasmamain[0].density[j]);
			}
		fprintf(fp_fe,"\n");
		}  


  alpha_agn = -1.2;

  IPstart = 20;
  IPstop = 41;
  lstart = 250;
  lmax = 450;
  distance = 1e11;
  nh = 1e5;

  fp_h = fopen ("hydrogen_sim.out", "w");
  fp_he = fopen ("helium_sim.out", "w");
  fp_c = fopen ("carbon_sim.out", "w");
  fp_n = fopen ("nitrogen_sim.out", "w");
  fp_o = fopen ("oxygen_sim.out", "w");
  fp_fe = fopen ("iron_sim.out", "w");
  fprintf (fp_h, "%i\n", lmax - lstart);
  fprintf (fp_he, "%i\n", lmax - lstart);
  fprintf (fp_c, "%i\n", lmax - lstart);
  fprintf (fp_n, "%i\n", lmax - lstart);
  fprintf (fp_o, "%i\n", lmax - lstart);
  fprintf (fp_fe, "%i\n", lmax - lstart);

  geo.nxfreq = 1;
  geo.xfreq[0] = 1e14;
  geo.xfreq[1] = 1e20;
  plasmamain[0].pl_alpha[0] = alpha_agn;
  plasmamain[0].spec_mod_type[0] = SPEC_MOD_PL;
  xband.nbands = 1;
  xband.f1[0] = 1e14;
  xband.f2[0] = plasmamain[0].max_freq = 1e18;
  plasmamain[0].t_e = 1e6;
  plasmamain[0].t_r = 1e6;
  plasmamain[0].rho = nh / rho2nh;
  plasmamain[0].ne = nh;

  for (i = IPstart; i < IPstop; i++)
    {
      IP = pow (10.0, i / 10.0);
      lum_agn = (IP * distance * distance * nh);
      const_agn =
	lum_agn /
	(((pow (50000 / HEV, alpha_agn + 1.)) -
	  pow (100 / HEV, alpha_agn + 1.0)) / (alpha_agn + 1.0));
      printf ("IP=%f, lum_agn=%e, const_agn=%e\n", IP, lum_agn, const_agn);
      agn_ip =
	const_agn *
	(((pow (50000 / HEV, alpha_agn + 1.0)) -
	  pow (100 / HEV, alpha_agn + 1.0)) / (alpha_agn + 1.0));
      agn_ip /= (distance * distance);
      agn_ip /= nh;
      Log ("i=%i,Ionisation Parameter=%f\n", i, (agn_ip));


      xweight =
	const_agn / ((4. * PI) * (4. * PI * distance * distance));
      printf ("weight=%e\n", xweight);
      variable_temperature (&plasmamain[0], 7);
      fprintf (fp_h, "%6.3e %6.3e %6.3e %6.3e\n", agn_ip, plasmamain[0].ne,
	       plasmamain[0].density[0], plasmamain[0].density[1]);
      fprintf (fp_he, "%e %e", agn_ip, plasmamain[0].ne);
      for (j = 2; j < 5; j++)
	{
	  fprintf (fp_he, " %6.3e", plasmamain[0].density[j]);
	}
      fprintf (fp_he, "\n");
      fprintf (fp_c, "%e %e", agn_ip, plasmamain[0].ne);
      for (j = 5; j < 12; j++)
	{
	  fprintf (fp_c, " %6.3e", plasmamain[0].density[j]);
	}
      fprintf (fp_c, "\n");
      fprintf (fp_n, "%e %e", agn_ip, plasmamain[0].ne);
      for (j = 12; j < 20; j++)
	{
	  fprintf (fp_n, " %6.3e", plasmamain[0].density[j]);
	}
      fprintf (fp_n, "\n");
      fprintf (fp_o, "%e %e", agn_ip, plasmamain[0].ne);
      for (j = 20; j < 29; j++)
	{
	  fprintf (fp_o, " %6.3e", plasmamain[0].density[j]);
	}
      fprintf (fp_o, "\n");
      fprintf (fp_fe, "%e %e", agn_ip, plasmamain[0].ne);
      for (j = 151; j < 178; j++)
	{
	  fprintf (fp_fe, " %6.3e", plasmamain[0].density[j]);
	}
      fprintf (fp_fe, "\n");
    }


*/





  return EXIT_SUCCESS;
}
