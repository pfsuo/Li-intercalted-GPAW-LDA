#!/bin/bash
#SBATCH -J C8Li_band
#SBATCH -N 1
#SBATCH -n 24
#SBATCH -o out
#SBATCH --time=4:00:00
#SBATCH -x cu[01,05]

#rm *.out *png
gpaw -P 24 python band.py
gpaw python plot_band.py
