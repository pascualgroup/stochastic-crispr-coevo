#!/usr/bin/env python3

import os
import random
import json

SCRIPT_DIR = os.path.abspath(os.path.dirname(__file__))
JULIA_SCRIPT_PATH = os.path.abspath(
    os.path.join(SCRIPT_DIR, 'julia', 'main.jl')
)

RUNS_PATH = os.path.join(SCRIPT_DIR, 'runs')

N_REPLICATES = 30

PARAMETERS = {
    "t_final" : 6000.0,
    "t_output" : 1.0,
    
    "n_bstrains" : 1,
    "n_hosts_per_bstrain" : 100,
    "n_vstrains" : 1,
    "n_particles_per_vstrain" : 100,
    
    "n_protospacers" : 15,
    "u_n_spacers_max" : 10,
    "p_crispr_failure_prob" : 1e-5,
    "q_spacer_acquisition_prob" : 1e-5,
    "r_growth_rate" : 1,
    "K_carrying_capacity" : 10**(5.5),
    "beta_burst_size" : 50,
    "phi_adsorption_rate" : 1e-7,
    "m_viral_decay_rate" : 0.1,
    "mu_viral_mutation_rate" : 5e-7,
    "rho_c_density_cutoff" : 1,
    
    "d_death_rate" : 0,
    
    "g_immigration_rate" : 0
}


SBATCH_TEMPLATE = \
'''#!/bin/bash

#SBATCH --job-name={job_name}
#SBATCH --chdir={job_dir}
#SBATCH --account=pi-pascualmm
#SBATCH --output=simulation_stdout.txt
#SBATCH --error=simulation_stderr.txt
#SBATCH --partition=broadwl
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=20:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=armun@uchicago.edu

module purge
module load julia

julia {julia_script_path} parameters.json
'''

SUBMIT_HEADER = \
'''#!/bin/bash

module load julia

'''


STACKED_PLOTS = \
'''#!/usr/bin/env python3

import pandas as pd
import numpy as np
import scipy as sp
import math
import matplotlib.pyplot as plt
import sys
import os

import seaborn as sns
from scipy import stats

def StackedPlotDF(data,tp):
    
    stacked_plot = pd.DataFrame()
    ID = "bact"

    if (tp=="bact"):
        ID = "bstrain_id"
    elif (tp=="vir"):
        ID = "vstrain_id"
    
    
    strains = data[ID].unique()
    tl = len(data["t"].unique())
    
    stacked_plot["t"] = data["t"].unique()

    for s in strains:

        Abs = np.zeros(tl)
        
        tmp = data[data[ID]==s]
        
        for i in tmp["t"].values:
            abun = tmp[tmp["t"]==i]["abundance"].values
            Abs[int(i)] = abun
        
        stacked_plot[s] = Abs
        
    
    return stacked_plot




# ######################################
SCRIPT_PATH = os.path.abspath(os.path.dirname(__file__))

bact = pd.read_csv(os.path.join(SCRIPT_PATH, 'babundance.csv'), delimiter=',')
phage = pd.read_csv(os.path.join(SCRIPT_PATH, 'vabundance.csv'), delimiter=',')
#data = pd.read_csv(path+'time-series-data.txt', delimiter=' ')



stacked_plot = StackedPlotDF(bact,"bact")
Vstacked_plot = StackedPlotDF(phage,"vir")


#this is the relative path of a particular simulation
sim_dir = os.path.relpath(SCRIPT_PATH,   os.path.join(SCRIPT_PATH,'..','..'));

#plt.figure(figsize=(10, 10), dpi= 80)
stacked_plot.plot.area(x="t",stacked=True, legend=False, linewidth=0);
#plt.show()
plt.title('Microbial Abundances: ' + sim_dir )
plt.xlabel('Time t')
plt.ylabel('Abundances N_i')
plt.tight_layout()
plt.savefig(os.path.join(SCRIPT_PATH,'..','..','..','plots',sim_dir,'Microbe-Abundance_stacked_plot.png'),dpi=500)
plt.close()


#plt.figure(figsize=(10, 10), dpi= 80)
Vstacked_plot.plot.area(x="t",stacked=True, legend=False, linewidth=0);
#plt.show()
plt.title('Viral Abundances: ' + sim_dir )
plt.xlabel('Time t')
plt.ylabel('Abundances V_i')
plt.tight_layout()
plt.savefig(os.path.join(SCRIPT_PATH,'..','..','..','plots',sim_dir,'Virus-Abundance_stacked_plot.png'),dpi=500)
plt.close()
    
'''


SBATCH_TEMPLATE_PY = \
'''#!/bin/bash

#SBATCH --job-name={job_name}
#SBATCH --chdir={job_dir}
#SBATCH --account=pi-pascualmm
#SBATCH --output=plot_stdout.txt
#SBATCH --error=plot_stderr.txt
#SBATCH --partition=broadwl
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=20:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=armun@uchicago.edu

module purge
module load python

python {python_script_path}
'''

SUBMIT_HEADER_PY = \
'''#!/bin/bash

module load python

'''



def main():
    # Create root output directory
    os.makedirs(RUNS_PATH)
    
    # Seed RNG
    seed_rng = random.SystemRandom()
    
    # Set up all directories
    run_paths = []
    for n_protospacers in [15]:
        for u_n_spacers_max in [10]:
            for g_X in [0, 1, 10, 100]:
                for run_path in set_up_replicates(seed_rng, n_protospacers, u_n_spacers_max, g_X):
                    run_paths.append(run_path)
    
    # Create a shell script to submit everything
    with open(os.path.join(SCRIPT_DIR, 'submit.sh'), 'w') as f:
        f.write(SUBMIT_HEADER)
        
        for run_path in run_paths:
            f.write('cd {}\n'.format(run_path))
            f.write('sbatch run.sbatch\n')

              
    # Create a shell script to submit plot makers
    with open(os.path.join(SCRIPT_DIR, 'submitPlots.sh'), 'w') as f:
        f.write(SUBMIT_HEADER_PY)
             
        for run_path in run_paths:
            f.write('cd {}\n'.format(run_path))
            f.write('sbatch runPLOTS.sbatch\n')

def set_up_replicates(seed_rng, n_protospacers, u_n_spacers_max, g_immigration_rate_X):
    g_immigration_rate = g_immigration_rate_X * 1 # second factor holds the place of the order of magnitude
    
    # Create a parameters file for each replicate
    for i in range(N_REPLICATES):
        # Set up parameters for this replicate
        parameters = {
            **PARAMETERS,
            "n_protospacers" : n_protospacers,
            "u_n_spacers_max" : u_n_spacers_max,
            "g_immigration_rate" : g_immigration_rate,
            "rng_seed" : seed_rng.randrange(1, 2**31 - 1)
        }
        
        # Create subdirectory for replicate
        run_path = os.path.join(
            RUNS_PATH,
            'nps={0}-u={1}-g={2}'.format(
                n_protospacers,
                u_n_spacers_max,
                g_immigration_rate_X
            ),
            '{0}'.format(i)
        )
        
        os.makedirs(run_path)
        
        # Create plot directory for replicate. makeStackedPlot.py will save into this directory
        plots_path = os.path.join(
                 SCRIPT_DIR, 'plots',
                 'nps={0}-u={1}-g={2}'.format(
                     n_protospacers,
                     u_n_spacers_max,
                     g_immigration_rate_X
                 ),
                 '{0}'.format(i)
             )
        
        os.makedirs(plots_path)
        
        
        
        # Save parameters file for replicate
        params_path = os.path.join(run_path, 'parameters.json')
        with open(params_path, 'w') as f:
            json.dump(parameters, f, indent = 4)
        
        # Save sbatch file for replicate
        sbatch_path = os.path.join(run_path, 'run.sbatch')
        with open(sbatch_path, 'w') as f:
            f.write(SBATCH_TEMPLATE.format(
                job_name = 'SC-nps={0}-u={1}-{2}'.format(
                    n_protospacers, u_n_spacers_max, i
                ),
                job_dir = run_path,
                julia_script_path = JULIA_SCRIPT_PATH
            ))
    
        
        
        stackedplot_path = os.path.join(run_path,'makeStackedPlots.py')
        with open(stackedplot_path, 'w') as f:
            f.write(STACKED_PLOTS)
        
        sbatchPLOTS_path = os.path.join(run_path, 'runPLOTS.sbatch')
        with open(sbatchPLOTS_path, 'w') as f:
            f.write(SBATCH_TEMPLATE_PY.format(
                    job_name = 'SC-nps={0}-u={1}-{2}'.format(
                        n_protospacers, u_n_spacers_max, i
                    ) + '-PLOTS',
                    job_dir = run_path,
                    python_script_path = stackedplot_path
                ))
                
                
        
        yield run_path

if __name__ == '__main__':
    main()
