#!/usr/bin/env python3

import os
import random
import json

SCRIPT_DIR = os.path.abspath(os.path.dirname(__file__))
JULIA_SCRIPT_PATH = os.path.abspath(
    os.path.join(SCRIPT_DIR, '..', 'julia', 'main.jl')
)

RUNS_PATH = os.path.join(SCRIPT_DIR, 'runs')

N_REPLICATES = 35

# This is the seed that generates seeds for each simulation. This is to compare simulations with and without immigration. This is not part of original code per Ed Baskerville
metaseed = 111
random.seed(metaseed)
seedList = [random.randrange(1, 2**63 - 1) for i in range(N_REPLICATES)]

RUN_PARAMETERS = {
    "t_final" : 5000.0,
    "t_output" : 1.0
}

INITIALIZATION_PARAMETERS = {
    "n_bstrains" : 1,
    "n_hosts_per_bstrain" : 100,
    "n_vstrains" : 1,
    "n_particles_per_vstrain" : 100
}

MODEL_PARAMETERS = {
    "n_protospacers" : 15,
    "u_n_spacers_max" : 10,
    "p_crispr_failure_prob" : 1e-5,
    "q_spacer_acquisition_prob" : 1e-5,
    "r_growth_rate" : 1,
    "K_carrying_capacity" : 10**(5.5),
    "beta_burst_size" : 50,
    "phi_adsorption_rate" : 1e-7,
    "m_viral_decay_rate" : 0.1,
    "mu1_viral_mutation_rate_coupled" : 1-5e-7,
    "mu2_viral_mutation_rate_decoupled" : 5e-7,
    "rho_c_density_cutoff" : 0.1,
    
    "d_death_rate" : 1e-18,
    
    #"decouple_viral_mutation" : false,
    
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
    
    # Set up all directories
    run_paths = []
    for n_protospacers in [15]:
        for u_n_spacers_max in [10]:
            for run_path in set_up_replicates(n_protospacers, u_n_spacers_max):
                run_paths.append(run_path)
    
    # Create a shell script to submit everything
    with open(os.path.join(SCRIPT_DIR, 'submit.sh'), 'w') as f:
        f.write(SUBMIT_HEADER)
        
        for run_path in run_paths:
            f.write('cd {}\n'.format(run_path))
            f.write('sbatch run.sbatch\n')

def set_up_replicates(n_protospacers, u_n_spacers_max):
    # Set up model parameters for these replicates
    model_parameters = {
        **MODEL_PARAMETERS,
        "n_protospacers" : n_protospacers,
        "u_n_spacers_max" : u_n_spacers_max
    }
    
    # Create a parameters file for each replicate
    for i in range(N_REPLICATES):
        # Set up run parameters (random seed) for each replicate
        run_parameters = {
            **RUN_PARAMETERS,
            "rng_seed" : seedList[i]
        }
        
        # Create subdirectory for replicate
        run_path = os.path.join(
            RUNS_PATH,
            'nps={0}-u={1}'.format(
                n_protospacers,
                u_n_spacers_max
            ),
            '{0}'.format(i)
        )
        os.makedirs(run_path)
        
        # Save parameters file for replicate
        params_path = os.path.join(run_path, 'parameters.json')
        with open(params_path, 'w') as f:
            json.dump({
                **run_parameters,
                **INITIALIZATION_PARAMETERS,
                **model_parameters
            }, f, indent = 4)
        
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
