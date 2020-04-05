#!/bin/bash

#SBATCH --job-name=StochasticCRISPR
#SBATCH --output=StochasticCRISPR.out
#SBATCH --error=StochasticCRISPR.err
#SBATCH --partition=broadwl
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1

julia main.jl
