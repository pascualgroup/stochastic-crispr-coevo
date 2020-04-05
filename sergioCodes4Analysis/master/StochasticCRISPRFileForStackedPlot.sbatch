#!/bin/bash

#SBATCH --job-name=StochasticCRISPRAbundanceStackedPlot
#SBATCH --output=StochasticCRISPRAbundanceStackedPlot.out
#SBATCH --error=StochasticCRISPRAbundanceStackedPlot.err
#SBATCH --partition=broadwl
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1

python makeFileForAbundanceStackedPlot.py
