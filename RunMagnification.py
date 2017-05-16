#!/usr/bin/python

import numpy
import subprocess

if __name__ == "__main__":
	theta_range = numpy.logspace(-2,1,20)
	cosmopars   = 'cosmo.pars'
	ndz_lens    = 'phi_lens_i.txt'
	ndz_sour    = 'phi_i2.txt'
	pdk_file    = 'Planck_2015_matterpower.dat'
	outfile     = 'test_out.ssv'
	

	nbins_lens  = str(len(numpy.loadtxt(ndz_lens)))
	nbins_sour  = str(len(numpy.loadtxt(ndz_sour)))
	nbins_pdk   = str(len(numpy.loadtxt(pdk_file)))

	ww = []
	for tht_ in theta_range:
		proc = subprocess.Popen(['./theo_kernel.exe',ndz_lens,ndz_sour,pdk_file,nbins_lens,nbins_sour,nbins_pdk,str(tht_),cosmopars],stdout=subprocess.PIPE)
		ww.append( float(proc.stdout.readline()) )

	
	theo = zip(theta_range,ww)
	numpy.savetxt(outfile,numpy.array(theo))
		

		
