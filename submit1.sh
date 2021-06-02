#!/bin/bash                                                                     
#SBATCH -N 1                                                                    
#SBATCH -n 1
#SBATCH -J nitrogen
#SBATCH --mem=10G

module load fftw

./azot ARGUMENT > stdout 2> stderr_ARGUMENT




