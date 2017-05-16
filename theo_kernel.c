#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "cosmo_mad.h"

float main(int argc, char *argv[])
{
	int ii_;

	FILE *file_ndz_lens;
	FILE *file_ndz_sour;
	FILE *file_pofk;
	FILE *file_pars;

	file_ndz_lens = fopen(argv[1],"r");
	file_ndz_sour = fopen(argv[2],"r");
	file_pofk     = fopen(argv[3],"r");
	file_pars     = fopen(argv[8],"r");

	int nbins_lens = (int)atoi(argv[4]);
	int nbins_sour = (int)atoi(argv[5]);
	int nbins_pofk = (int)atoi(argv[6]);

	char z_c[200], phi_c[200];
	double *z_lens = (double *)malloc(nbins_lens * sizeof(double));
	double *z_sour = (double *)malloc(nbins_sour * sizeof(double));
	double *phi_lens = (double *)malloc(nbins_lens * sizeof(double));
	double *phi_sour = (double *)malloc(nbins_sour * sizeof(double));

	for(ii_=0; ii_<nbins_lens; ii_++)
	{
		fscanf(file_ndz_lens,"%s %s\n",z_c,phi_c);
		z_lens[ii_]   = (double)atof(z_c);
		phi_lens[ii_] = (double)atof(phi_c);
	}
	for(ii_=0; ii_<nbins_sour; ii_++)
	{

		fscanf(file_ndz_sour,"%s %s\n",z_c,phi_c);
		z_sour[ii_]   = (double)atof(z_c);
		phi_sour[ii_] = (double)atof(phi_c);
	}

	char k_c[200], Pk_c[200];
	double *k  = (double *)malloc(nbins_pofk * sizeof(double));
	double *Pk = (double *)malloc(nbins_pofk * sizeof(double));

	for(ii_=0; ii_<nbins_pofk; ii_++)
	{
		fscanf(file_pofk,"%s %s\n",k_c,Pk_c);
		k[ii_]  = (double)atof(k_c);
		Pk[ii_] = (double)atof(Pk_c);
	}

	double tht = (double)atof(argv[7]);
	double cth = cos(tht*M_PI/180.);
	
	double cosmopars[7];
	char pars_c[200],dummy_c[200];
	for(ii_=0; ii_<7; ii_++)
	{
		fscanf(file_pars,"%s %s\n",dummy_c,pars_c);
		cosmopars[ii_] = (double)atof(pars_c);
	}
	
	double OmM  = cosmopars[0];
	double OmL  = cosmopars[1];
	double Omb  = cosmopars[2];
	double w0   = cosmopars[3];
	double wa   = cosmopars[4];
	double h0   = cosmopars[5];
	double Tcmb = cosmopars[6];
	
	Csm_params *pars = csm_params_new();
	csm_unset_gsl_eh();
	csm_set_verbosity(0);
	csm_background_set(pars,OmM,OmL,Omb,w0,wa,h0,Tcmb);
	csm_set_linear_pk(pars,argv[3],-5,3,0.01,0.95,0.836);
	csm_set_nonlinear_pk(pars,argv[3]);

	int zl_, zs_;
	double ww=0.;
	for(zl_=0; zl_<nbins_lens-1; zl_++)
	{
		double rl = csm_radial_comoving_distance(pars,1./(1.+z_lens[zl_]));
		double gf = csm_growth_factor( pars,1./(1.+z_lens[zl_]) )/csm_growth_factor(pars,1);

		double w0 = 0.;
		for(ii_=0; ii_<nbins_pofk-1; ii_++)
		{
			w0 += (1./(2.*M_PI))*k[ii_]*(k[ii_+1]-k[ii_])*0.5*(Pk[ii_+1]+Pk[ii_])*
                                      sin(rl*k[ii_]*tht*M_PI/180.)/(rl*k[ii_]*tht*M_PI/180.);
		}
	
		double K = 0;
		for(zs_=zl_+1; zs_<nbins_sour-1; zs_++)
		{
			double rs = csm_radial_comoving_distance(pars,1./(1.+z_sour[zs_]));
			K += (z_sour[zs_+1]-z_sour[zs_])*phi_sour[zs_]*(rs-rl)/rs;
		}
		
		ww += (3./2.)*(OmM/(3000.*3000.))*h0*h0*rl*(1./(1+z_lens[zl_]))*(z_lens[zl_+1]-z_lens[zl_])*phi_lens[zl_]*K*w0*gf*gf;

	}
	printf("%f",ww);
	return ww;
}
