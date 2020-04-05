#!/bin/bash

#SBATCH --job-name=StochasticCRISPROutput_virus
#SBATCH --output=StochasticCRISPROutput_virus.out
#SBATCH --error=StochasticCRISPROutput_virus.err
#SBATCH --partition=broadwl
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1

python StochasticOutput_virus.py 10 15
