#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "cosmo_mad.h"

int main(int argc, char *argv[])
{
	FILE *f_out, *fndz1, *fndz2, *f_cov, *f_mic, *f_cl, *fpk;
	f_out = fopen(argv[3],"w");
	fndz1 = fopen(argv[1],"r");
	fndz2 = fopen(argv[2],"r");
	fpk   = fopen("Planck_2015_matterpower.dat","r");

	int th, ii, l, jj, zz1, zz2, zz3;
	double kPk[583], Pk[583];
	char kPk_c[200], Pk_c[583];

	for(ii=0; ii<583; ii++)
	{
		fscanf(fpk,"%s %s\n",kPk_c,Pk_c);
		kPk[ii] = (double)atof(kPk_c);
		Pk[ii]  = (double)atof(Pk_c);
	}

	Csm_params *pars = csm_params_new();
	csm_unset_gsl_eh();
	csm_set_verbosity(0);
	csm_background_set(pars,0.32,0.68,0.047,-1.,0.,0.68,2.725);
	csm_set_linear_pk(pars,"Planck_2015_matterpower.dat",-5,3,0.01,0.95,0.836);
	csm_set_nonlinear_pk(pars,"Planck_2015_matterpower.dat");
	//csm_set_xi_multipole_splines(pars);

	/*double theta[50];
	for( th=0; th<50; th++)
	{
		theta[th] = exp((th+0.5)/50*(log(1)-log(0.01))+log(0.01));
	}*/
	double theta[6] = {0.0147, 0.0316, 0.068, 0.146, 0.316, 0.681};
	//double theta[6] = {0.0077, 0.0188, 0.0454, 0.1099, 0.2659, 0.6430};

	char z_c[200], phi_c[200];
	double z1[200], z2[200], phi1[200], phi2[200];
	double integralz1 = 0., integralz2 = 0.;
	for(ii=0; ii<200; ii++)
	{
		fscanf(fndz1,"%s %s\n",z_c,phi_c);
		z1[ii]   = atof(z_c);
		phi1[ii] = atof(phi_c);

		fscanf(fndz2,"%s %s\n",z_c,phi_c);
		z2[ii]   = atof(z_c);
		phi2[ii] = atof(phi_c);
	}
	for(ii=0; ii<200; ii++)
	{
		printf("%d -- %f %f %f %f\n",ii,z1[ii],z2[ii],phi1[ii],phi2[ii]);
	}
	double w_cross[6],w_lens[6],w_back[6];
	for(th=0; th<6; th++)
	{
		double ww = 0.;

		double r0 = 3000.;
		double H0 = 68.00, Om0 = 0.31;
			
		double cth = cos(theta[th]*M_PI/180.);
		for(zz1=0; zz1<199; zz1++)
		{
			double r1 = csm_radial_comoving_distance(pars,1./(1.+z1[zz1]));
			double gf1 = csm_growth_factor( pars,1./(1.+z1[zz1]) )/csm_growth_factor(pars,1);

			double w0 = 0.;
			for(ii=0; ii<582; ii++)
			{
				w0 += (1./(2.*M_PI))*kPk[ii]*(kPk[ii+1]-kPk[ii])*0.5*(Pk[ii+1]+Pk[ii])*
				      sin(r1*kPk[ii]*theta[th]*M_PI/180.)/(r1*kPk[ii]*theta[th]*M_PI/180.);
			}

			double K = 0.;
			for(zz2=zz1+1; zz2<199; zz2++)
			{
				double r2 = csm_radial_comoving_distance(pars,1./(1.+z2[zz2]));
				K += (z2[zz2+1]-z2[zz2])*phi2[zz2]*(r2-r1)/r2;	
			}
			//printf("%f\t%f\t%f\n",r1,K,w0);
			ww +=  (1./(0.68*0.68))*(3.*Om0/(3000.*3000.))*(H0*H0/(100.*100.))*r1*(1./(1+z1[zz1]))*(z1[zz1+1]-z1[zz1])*phi1[zz1]*K*w0*gf1*gf1;
		}
		w_cross[th] = ww;
		fprintf(f_out,"%f\t%f\n",theta[th],ww);
		printf("%f\t%f\n",theta[th],ww);
	}
	return 0;
}
