#!/bin/bash
#SBATCH -J C8RPA_K
#SBATCH -N 1
#SBATCH -n 40
#SBATCH -o out
#SBATCH --time=2400:00:00
#SBATCH -w cu16

ulimit -s unlimited
gpaw -P 4 python eels.py
