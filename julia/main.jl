#=
main:
- Julia version: 
- Author: ed
- Date: 2019-10-01
=#

push!(LOAD_PATH, ".")

using Logging
using StochasticCrispr
using Dates

function get_initialization_parameters() :: InitializationParameters
    p = InitializationParameters()

    p.n_bstrains = 1
    p.n_hosts_per_bstrain = 1
    p.n_spacers = 8


    modelparams = get_model_parameters()
    validate(modelparams)
    p.n_vstrains = 1
    p.n_particles_per_vstrain = 1
    p.n_protospacers = 10

    p
end

function get_model_parameters() :: Parameters
    p = Parameters()

    p.u_n_spacers_max = 8

    p.v_n_protospacers_max = 10

    p.p_crispr_failure_prob = 1e-5
    p.q_spacer_acquisition_prob = 1e-5
    p.r_growth_rate = 1
    p.K_carrying_capacity = 10^5.5
    p.beta_burst_size = 50
    p.phi_adsorption_rate = 1e-7
    p.m_viral_decay_rate = 0.1
    p.mu_mutation_rate = 5e-7
    p.rho_c_density_cutoff = 0.1

    p.d_death_rate = 0.05

    p
end


function main()
    logger = SimpleLogger(stderr, Logging.Info)
    global_logger(logger)

    wallclock_start = now()
    println("starting at:", wallclock_start)

    initparams = get_initialization_parameters()

    modelparams = get_model_parameters()
    validate(modelparams)

    state = State(initparams, modelparams)
    println("Initial spacers:\n", state.bstrains.spacers)
    println("Initial ppacers:\n", state.vstrains.pspacers)
    validate(state)

    sim = Simulator(modelparams, 0.0, state)
    @debug "top_level_event_rates" sim.top_level_event_rates

    println("sim.t: ", sim.t)
    t_final = 40.0
    n_events::UInt64 = 0
    while sim.t < t_final
        n_events += 1
        do_next_event!(sim, t_final)
        #println("sim.state.t: ", sim.t)
    end
    println("sim.t: ", sim.t)
    println(sim.state)
    println(sim.state.bstrains.total_abundance)
    println("total number of events: ", n_events)

    wallclock_end = now()
    println("ending at ", wallclock_end)
    elapsed_time = Dates.value(wallclock_end - wallclock_start) / 1000.0
    println("elapsed time: ", elapsed_time, " seconds")
 end

main()
