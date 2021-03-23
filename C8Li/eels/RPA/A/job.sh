#!/bin/bash
#SBATCH -J C8RPA_A
#SBATCH -N 1
#SBATCH -n 40
#SBATCH -o out
#SBATCH --time=2400:00:00
#SBATCH -w cu16

gpaw -P 10 python eels.py
