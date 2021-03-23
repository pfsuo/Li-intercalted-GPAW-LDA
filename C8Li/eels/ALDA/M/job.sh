#!/bin/bash
#SBATCH -J C8LDA_M
#SBATCH -N 1
#SBATCH -n 52
#SBATCH -o out
#SBATCH --time=2400:00:00
#SBATCH -w cu[18]

ulimit -s unlimited
gpaw -P 6 python eels.py
