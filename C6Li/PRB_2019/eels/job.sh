#!/bin/bash
#SBATCH -J PRB_2019
#SBATCH -N 1
#SBATCH -n 40
#SBATCH -o out
#SBATCH --time=2400:00:00
#SBATCH -w cu16

export OMP_NUM_THREADS=5
#rm -rf RPA ALDA
gpaw -P 8 python eels_calc.py
