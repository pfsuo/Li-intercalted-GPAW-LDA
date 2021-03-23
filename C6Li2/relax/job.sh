#!/bin/bash
#SBATCH -J C6Li2
#SBATCH -n 24
#SBATCH -N 1
#SBATCH -o out
#SBATCH --time=2400:00:00
#SBATCH -x cu[09-18]

#export OMP_NUM_THREADS=16
gpaw -P 24 python relax.py
