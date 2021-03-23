#!/bin/bash
#SBATCH -J C8LDA_A
#SBATCH -N 1
#SBATCH -n 40
#SBATCH -o out
#SBATCH --time=2400:00:00
#SBATCH -w cu[16]

ulimit -s unlimited
gpaw -P 6 python eels.py
