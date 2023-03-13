#!/bin/sh

#SBATCH --account=pi-pascualmm
#SBATCH --partition=broadwl
#SBATCH --job-name=crispr-test
#SBATCH --tasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2000m
#SBATCH --time=4:00:00
#SBATCH --chdir=/home/armun/crispr-sweep-6-9-2021/individual-test
#SBATCH --output=output.txt
#SBATCH --mail-user=armun@uchicago.edu
module purge
# Uncomment this to use the Midway-provided Julia:
module load julia
julia /home/armun/crispr-sweep-6-9-2021/main.jl parameters.json &> script_output.txt
