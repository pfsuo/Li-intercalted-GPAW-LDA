#!/bin/bash
#SBATCH -J C8RPA_M
#SBATCH -N 1
#SBATCH -n 52
#SBATCH -o out
#SBATCH --time=2400:00:00
#SBATCH -w cu18

gpaw -P 6 python eels.py
