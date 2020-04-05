#=
main:
- Julia version: 
- Author: ed
- Date: 2019-10-01
=#


using Logging
using Random
using Dates

# Uncomment this to see all debugging output:
# Logging.global_logger(
#     Logging.SimpleLogger(
#         stderr,
#         Logging.Debug
#     )
# )

push!(LOAD_PATH, ".")
using StochasticCrispr

function make_run_parameters() :: RunParameters
    p = RunParameters()
    
    p.t_final = 3000.0
    p.t_output = 1.0
    
#     p.rng_seed = 16332687649108451454   #set seed (it would be great taking it optionally as a command line argument)
    
    p
end

function make_initialization_parameters() :: InitializationParameters
    p = InitializationParameters()

    p.n_bstrains = 1
    p.n_hosts_per_bstrain = 100
    # p.n_spacers = 8
#     p.n_spacers = 10    
    
    p.n_vstrains = 1
    p.n_particles_per_vstrain = 100
    p.n_protospacers = 15

    p
end

function make_model_parameters() :: Parameters
    p = Parameters()

    #p.u_n_spacers_max = 8
    p.u_n_spacers_max = 10

#     p.v_n_protospacers_max = 10
    p.v_n_protospacers_max = 15

    p.p_crispr_failure_prob = 1e-5
    p.q_spacer_acquisition_prob = 1e-5
    p.r_growth_rate = 1
    p.K_carrying_capacity = 10^5.5
    p.beta_burst_size = 50
    p.phi_adsorption_rate = 1e-7
    p.m_viral_decay_rate = 0.1
#     p.mu_mutation_rate = 5e-7
    p.mu_mutation_rate = 1e-7
    p.rho_c_density_cutoff = 0.1

    p.d_death_rate = 0.05

    p
end

run_model(
    make_run_parameters(),
    make_initialization_parameters(),
    make_model_parameters()
)
