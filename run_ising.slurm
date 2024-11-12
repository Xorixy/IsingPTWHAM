#!/bin/bash
#SBATCH --job-name=ising_pt_wham
#SBATCH -n 4

mpirun ./Ising

