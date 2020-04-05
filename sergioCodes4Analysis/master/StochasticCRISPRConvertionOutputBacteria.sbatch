#!/bin/bash

#SBATCH --job-name=StochasticCRISPROutput_bact
#SBATCH --output=StochasticCRISPROutput_bact.out
#SBATCH --error=StochasticCRISPROutput_bact.err
#SBATCH --partition=broadwl
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1

python StochasticOutput_bact.py 10 15
