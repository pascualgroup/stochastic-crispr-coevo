#!/bin/bash

#SBATCH --job-name=StochasticCRISPRAbundanceOutput.py
#SBATCH --output=StochasticCRISPRAbundanceOutput.out
#SBATCH --error=StochasticCRISPRAbundanceOutput.err
#SBATCH --partition=broadwl
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1

python AbundanceStochasticOutput.py
