#!/bin/bash
#SBATCH -J C6Li2_band
#SBATCH -N 1
#SBATCH -n 24
#SBATCH -o out
#SBATCH --time=72:00:00
#SBATCH -x cu[15-18]

gpaw -P 24 python band.py
gpaw python plot_band.py
