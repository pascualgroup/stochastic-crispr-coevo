#!/usr/bin/env python3

import os
import random
import json

SCRIPT_DIR = os.path.abspath(os.path.dirname(__file__))
JULIA_SCRIPT_PATH = os.path.abspath(
    os.path.join(SCRIPT_DIR, 'julia', 'main.jl')
)

RUNS_PATH = os.path.join(SCRIPT_DIR, 'runs')

N_REPLICATES = 1

PARAMETERS = {
    "t_final" : 5000.0,
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
    "mu_viral_mutation_rate" : 1-5e-7,
    "rho_c_density_cutoff" : 0.1,
    
    "d_death_rate" : 1e-18,
    
    "g_immigration_rate" : 1e-18
}


SBATCH_TEMPLATE = \
'''#!/bin/bash

#SBATCH --job-name={job_name}
#SBATCH --chdir={job_dir}
#SBATCH --account=pi-pascualmm
#SBATCH --output=stdout.txt
#SBATCH --error=stderr.txt
#SBATCH --partition=broadwl
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=02:00:00

module purge
module load julia

julia {julia_script_path} parameters.json
'''

SUBMIT_HEADER = \
'''#!/bin/bash

module load julia

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
            for g_X_1e18 in [0, 1, 2, 4]:
                for run_path in set_up_replicates(seed_rng, n_protospacers, u_n_spacers_max, g_X_1e18):
                    run_paths.append(run_path)
    
    # Create a shell script to submit everything
    with open(os.path.join(SCRIPT_DIR, 'submit.sh'), 'w') as f:
        f.write(SUBMIT_HEADER)
        
        for run_path in run_paths:
            f.write('cd {}\n'.format(run_path))
            f.write('sbatch run.sbatch\n')

def set_up_replicates(seed_rng, n_protospacers, u_n_spacers_max, g_immigration_rate_X_1e18):
    g_immigration_rate = g_immigration_rate_X_1e18 * 1e-18
    
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
                g_immigration_rate_X_1e18
            ),
            '{0}'.format(i)
        )
        os.makedirs(run_path)
        
        # Create plot directory for replicate
        plots_path = os.path.join(
                 SCRIPT_DIR, 'plots',
                 'nps={0}-u={1}-g={2}'.format(
                     n_protospacers,
                     u_n_spacers_max,
                     g_immigration_rate_X_1e18
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
        
        yield run_path

if __name__ == '__main__':
    main()
