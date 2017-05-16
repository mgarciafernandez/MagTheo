# MagTheo
Script to compute theoretical predictions for magnification

## Description of the files

* `cosmo.pars` file where to put the cosmological parameters.
* `./RunMagnification.py` script to call the C kernel that does the actual computation of the correlation function.
* `./theo_kernel.exe` binary that does the actual calculation of the correlation function
* `./theo_kernel.c` source code of the binary.
* `./Makefile_cosmomad` makefile to compile the C code (warning! put properly the PATH).
* `./cosmomad-0.10.tar.gz` cosmomad libraries (needed).

* `Planck_2015_matterpower.dat`, `phi_i2.txt` and `phi_lens.txt` are the input files that give as result `test_out.csv`. This is just for testing pruposes.

## Description of the pipeline.

The file `./theo_kernel.exe` does the actual calculation. It takes as input the n(z) filenames of the lens and source sample, the P(k) as given by CAMB, the number of bins of each file and the angle at which to compute the correlation function. It returns the value of the correlation function.

The file `./RunMagnification.py` simply call the previous binary. It can be customized at will. There are provided the names of the files.

## Format of the input files
The file `cosmo.pars` should not be changed (except the numbers). The ordering of the cosmological parameters must be what is on that file (and separated with a blanckspace form the value).

The n(z) files must be .ssv

The P(k) can be given directly by CAMB or .ssv.
