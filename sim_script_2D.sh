#!/bin/bash
#SBATCH -p standard
#SBATCH -t 5-00:00:00
#SBATCH -c 1
#SBATCH -n 1
#SBATCH -a 1-100
#SBATCH --mem=6g
unset DISPLAY
module load matlab
matlab -r Ising_bluehive_sim_script_2D
