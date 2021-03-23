#!/bin/bash
#SBATCH -J C8Li_eels
#SBATCH -N 1
#SBATCH -n 40
#SBATCH -o out
#SBATCH --time=2400:00:00
#SBATCH -w cu[16]

export OMP_NUM_THREADS=5
gpaw -P 8 python eels_calc.py
