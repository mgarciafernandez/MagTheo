# MagTheo
Script to compute theoretical predictions for magnification

## Description of the files

* `cosmo.pars` file where to put the cosmological parameters.
* `./RunMagnification.py` script to call the C kernel that does the actual computation of the correlation function.
* `./theo_kernel.exe` binary that does the actual calculation of the correlation function
* `./theo_kernel.c` source code of the binary.
* `./Makefile_cosmomad` makefile to compile the C code (warning! put properly the PATH).
* `./cosmomad-0.10.tar.gz` cosmomad libraries (needed).
