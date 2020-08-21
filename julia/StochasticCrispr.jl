module StochasticCrispr

# Temporary comment to test git

using Util
using Random
using Logging
using Distributions
using StatsBase
using DelimitedFiles
using Dates
using JSON2

export run_model
export RunParameters, InitializationParameters, Parameters
export State, Simulator
export validate, do_next_event!

mutable struct RunParameters
    "Simulation end time"
    t_final::Float64

    "Time between output events"
    t_output::Float64

    "Seed for random number generator"
    rng_seed::Union{UInt64, Nothing}

    "Enable output?"
    enable_output::Bool

    function RunParameters()
        p = new()
        p.rng_seed = nothing
        p.enable_output = true
        p
    end
end
JSON2.@format RunParameters noargs

function validate(p::RunParameters)
    @assert 0.0 < p.t_final <= 1000
    @assert 0.0 < p.t_output <= p.t_final
end

mutable struct InitializationParameters
    "Number of initial bacterial strains"
    n_bstrains::UInt64

    "Number of initial hosts per bacterial strain"
    n_hosts_per_bstrain::UInt64

    # "Number of initial spacers per bacterial strain"
    # n_spacers::UInt64

    "Number of initial virus strains"
    n_vstrains::UInt64

    "Number of initial particles per bacterial strain"
    n_particles_per_vstrain::UInt64

    "Number of initial protospacers per virus strain"
    n_protospacers::UInt64

    InitializationParameters() = new()
end
JSON2.@format InitializationParameters noargs

function validate(p::InitializationParameters)
    @assert p.n_bstrains >= 1
    @assert p.n_hosts_per_bstrain >= 1
    @assert p.n_vstrains >= 1
    @assert p.n_particles_per_vstrain >= 1
    @assert p.n_protospacers >= 1
end

mutable struct Parameters
    "Maximum number of spacers in a bacterial strain"
    u_n_spacers_max

    "Maximum number of protospacers in a virus strain"
    v_n_protospacers_max

    "CRIPSR failure probability"
    p_crispr_failure_prob::Float64

    "New spacer acquisition probability"
    q_spacer_acquisition_prob::Float64

    "Growth rate at 0 (1/h)"
    r_growth_rate::Float64

    "Carrying capacity (1/mL)"
    K_carrying_capacity::Float64

    "Burst size"
    beta_burst_size::UInt64

    "Adsorption rate"
    phi_adsorption_rate::Float64

    "Viral decay rate"
    m_viral_decay_rate::Float64

    "Mutation rate"
    mu_mutation_rate::Float64

    "Density cutoff: used to scale volume of system, and therefore discrete population sizes"
    rho_c_density_cutoff::Float64

    "Constant death rate (not in Childs model)"
    d_death_rate::Float64

    Parameters() = new()
end
JSON2.@format Parameters noargs

function validate(p::Parameters)
    @assert 1 <= p.u_n_spacers_max
    @assert 1 <= p.v_n_protospacers_max
    @assert 0.0 <= p.p_crispr_failure_prob <= 1.0
    @assert 0.0 <= p.q_spacer_acquisition_prob <= 1.0
    @assert 0.1 < p.r_growth_rate <= 10
    @assert 0.0 < p.K_carrying_capacity <= 1e7
    @assert 0 <= p.beta_burst_size <= 1000
    @assert 0.0 <= p.phi_adsorption_rate <= 1e-6
    @assert 0.0 <= p.m_viral_decay_rate <= 1.0
    @assert 0.0 <= p.mu_mutation_rate <= 1e-5
    @assert 0.05 <= p.rho_c_density_cutoff <= 1.0
    @assert 0.0 < p.d_death_rate
end

function open_csv(prefix, header...)
    filename = "$prefix.csv"
    if ispath(filename)
        error("$filename already exists. You should delete output, or run in a different directory.")
    else
        file = open(filename, "w")
        write_csv(file, header...)
        file
    end
end

function write_csv(file, row...)
    for i = 1:lastindex(row)
        print(file, row[i])
        if i < lastindex(row)
            print(file, ",")
        end
    end
    print(file, "\n")
end

function run_model(rp::RunParameters, ip::InitializationParameters, mp::Parameters)
    validate(rp)
    validate(ip)
    validate(mp)
    
    meta_file = open_csv("meta", "key", "value")

    # Use random seed if provided, or generate one
    rng_seed = rp.rng_seed === nothing ? rand(RandomDevice(), UInt64) : rp.rng_seed
    write_csv(meta_file, "rng_seed", rng_seed)

    # Record start time
    start_time = now()
    write_csv(meta_file, "start_time", start_time)

    # Initialize & validate model state
    state = State(ip, mp)
    validate(state)

    # Initialize simulator
    sim = Simulator(rp, mp, 0.0, state, MersenneTwister(rng_seed))
    
    # Initial output
    if rp.enable_output
        @info "initial output"
        write_periodic_output(sim)
        write_strains(
            state.bstrains.strain_file,
            state.bstrains.spacers_file,
            sim.t,
            state.bstrains.ids,
            state.bstrains.spacers
        )
        write_strains(
            state.vstrains.strain_file,
            state.vstrains.spacers_file,
            sim.t,
            state.vstrains.ids,
            state.vstrains.pspacers
        )
    end
    
    # Simulation loop
    t_next_output = 0.0
    
    while sim.t < rp.t_final
        # Simulate exactly until the next output time 
        t_next_output = min(rp.t_final, t_next_output + rp.t_output)
        
        @info "Beginning period: $(sim.t) to $(t_next_output)"
        
        while sim.t < t_next_output
            # Perform the next event.
            # If the next event time is computed to be
            # greater than t_next_output, no event will occur
            # and time will simply advance exactly to t_next_output.
            do_next_event!(sim, t_next_output)
        end

        @debug "event counts:" total=n_events, breakdown=sim.event_counts
        @debug "bstrains:" total_abund=state.bstrains.total_abundance abund=state.bstrains.abundance spacers=state.bstrains.spacers
        @debug "vstrains:" total_abund=state.vstrains.total_abundance abund=state.vstrains.abundance pspacers=state.vstrains.pspacers
        
        # Write periodic output
        @debug "rp.enable_output" rp.enable_output
        if rp.enable_output
            write_periodic_output(sim)
        end
        
        @assert sim.t == t_next_output
    end

    # Record end time and elapsed
    end_time = now()
    write_csv(meta_file, "end_time", end_time)
    elapsed_seconds = Dates.value(end_time - start_time) / 1000.0
    write_csv(meta_file, "elapsed_seconds", elapsed_seconds)

    close(meta_file)
end

mutable struct BStrains
    next_id::UInt64
    ids::Vector{UInt64}

    abundance::Vector{UInt64}
    total_abundance::UInt64

    spacers::Vector{Vector{UInt64}}
    
    strain_file::IOStream
    spacers_file::IOStream
    abundance_file::IOStream

    function BStrains(n_strains, n_hosts_per_strain)
        new(
            # next_id
            n_strains + 1,

            # ids
            1:n_strains,

            # abundance
            repeat([n_hosts_per_strain], n_strains),

            # total_abundance
            n_strains * n_hosts_per_strain,

            # spacers
            repeat([[]], n_strains),
            
            # strain_file
            open_csv("bstrains", "t_creation", "bstrain_id", "parent_bstrain_id", "infecting_vstrain_id"),
            
            # spacers_file
            open_csv("bspacers", "bstrain_id", "spacer_id"),
            
            # abundance_file
            open_csv("babundance", "t", "bstrain_id", "abundance")
        )
    end
end

function remove_bstrain!(b::BStrains, index)
    # This is only used when a bstrain has gone extinct
    @assert b.abundance[index] == 0

    @debug "Removing bstrain" id=b.ids[index] index=index

    swap_with_end_and_remove!(b.ids, index)
    swap_with_end_and_remove!(b.abundance, index)
    swap_with_end_and_remove!(b.spacers, index)
end

mutable struct VStrains
    next_id::UInt64
    ids::Vector{UInt64}

    abundance::Vector{UInt64}
    total_abundance::UInt64

    next_pspacer_id::UInt64
    pspacers::Vector{Vector{UInt64}}
    
    strain_file::IOStream
    spacers_file::IOStream
    abundance_file::IOStream

    function VStrains(n_strains, n_particles_per_strain, n_pspacers_init)
        next_id = n_strains + 1
        ids = Vector(1:n_strains)
        abundance = repeat([n_particles_per_strain], n_strains)
        total_abundance = n_strains * n_particles_per_strain
        next_pspacer_id = 1 + n_strains * n_pspacers_init
        pspacers = [
            Vector(1:n_pspacers_init) .+ repeat([n_pspacers_init * (i - 1)], n_pspacers_init)
            for i = 1:n_strains
        ]

        new(
            next_id,
            ids,
            abundance,
            total_abundance,
            next_pspacer_id,
            pspacers,
            open_csv("vstrains", "t_creation", "vstrain_id", "parent_vstrain_id", "infected_bstrain_id"),
            open_csv("vpspacers", "vstrain_id", "spacer_id"),
            open_csv("vabundance", "t", "vstrain_id", "abundance")
        )
    end
end

function remove_vstrain!(v::VStrains, index)
    # This is only used when a bstrain has gone extinct
    @assert v.abundance[index] == 0

    @debug "Removing vstrain" id=v.ids[index] index=index

    swap_with_end_and_remove!(v.ids, index)
    swap_with_end_and_remove!(v.abundance, index)
    swap_with_end_and_remove!(v.pspacers, index)
end

mutable struct State
    bstrains::BStrains
    vstrains::VStrains

    function State(
        n_bstrains, n_hosts_per_bstrain, # n_spacers_init,
        n_vstrains, n_particles_per_vstrain, n_pspacers_init
    )
        new(
            BStrains(n_bstrains, n_hosts_per_bstrain),
            VStrains(n_vstrains, n_particles_per_vstrain, n_pspacers_init)
        )
    end
end

function State(ip::InitializationParameters, mp::Parameters)
    State(
        ip.n_bstrains, ip.n_hosts_per_bstrain, # ip.n_spacers,
        ip.n_vstrains, ip.n_particles_per_vstrain, ip.n_protospacers
    )
end


const BACTERIAL_GROWTH = 1
const BACTERIAL_DEATH = 2
const VIRAL_DECAY = 3
const CONTACT = 4

const TOP_LEVEL_EVENTS = [
    BACTERIAL_GROWTH,
    BACTERIAL_DEATH,
    VIRAL_DECAY,
    CONTACT,
]

mutable struct Simulator
    runparams::RunParameters
    parameters::Parameters
    t::Float64
    state::State
    rng::MersenneTwister

    top_level_event_rates::Vector{Float64}
    event_counts::Vector{UInt64}

    summary_file::IOStream
end

function Simulator(runparams::RunParameters, params::Parameters, t_init::Float64, state_init::State, rng::MersenneTwister)
    sim = Simulator(
        runparams, params, t_init, state_init, rng,
        zeros(length(TOP_LEVEL_EVENTS)), zeros(length(TOP_LEVEL_EVENTS)),
        open_csv("summary", "t", "bacterial_abundance", "viral_abundance")
    )
    for event_id = TOP_LEVEL_EVENTS
        update_rate!(event_id, sim)
    end
    sim
end

function do_next_event!(sim::Simulator, t_max::Float64)
    p = sim.parameters
    s = sim.state

    @assert length(sim.top_level_event_rates) == length(TOP_LEVEL_EVENTS)
    R = sum(sim.top_level_event_rates)

    # Draw next event time using total rate
    t_next = sim.t + randexp(sim.rng) / R

    if t_next > t_max
        sim.t = t_max
    else
        @debug "event_rates:" sim.top_level_event_rates

        # Sample next top-level event proportional to event rate
        event_id = sample(sim.rng, TOP_LEVEL_EVENTS, Weights(sim.top_level_event_rates, R))
        sim.event_counts[event_id] += 1

        @debug "begin do_event()" event=event t=t_next
        do_event!(event_id, sim, t_next)
        update_all_rates!(sim)
        @debug "end do_event()"

        sim.t = t_next

        @debug "bstrains.total_abundance:" sim.state.bstrains.total_abundance
        @debug "VStrains.total_abundance:" sim.state.vstrains.total_abundance
    end
end

function write_periodic_output(sim)
    s = sim.state

    write_summary(sim.summary_file, sim.t, s)
    
    write_abundances(
        s.bstrains.abundance_file, sim.t, s.bstrains.ids, s.bstrains.abundance
    )
    write_abundances(
        s.vstrains.abundance_file, sim.t, s.vstrains.ids, s.vstrains.abundance
    )
    
    flush(s.bstrains.strain_file)
    flush(s.bstrains.spacers_file)
    flush(s.vstrains.strain_file)
    flush(s.vstrains.spacers_file)
end

function write_summary(file, t, state)
    write_csv(file, t, state.bstrains.total_abundance, state.vstrains.total_abundance)
    flush(file)
end

function write_strains(strain_file, spacers_file, t_creation, ids, spacers)
    @debug "write_strains" ids
    for i = 1:lastindex(ids)
        write_strain(strain_file, t_creation, ids[i], 0, 0)
        write_spacers(spacers_file, ids[i], spacers[i])
    end
    flush(strain_file)
    flush(spacers_file)
end

function write_strain(file, t_creation, id, parent_id, other_id)
    write_csv(file, t_creation, id, parent_id, other_id)
end

function write_spacers(file, id, spacers)
    @debug "write_spacers" id spacers
    for spacer_id = spacers
        write_csv(file, id, spacer_id)
    end
end

function write_abundances(file, t, ids, abundance)
    @debug "write_abundances" t ids abundance
    @debug "lastindex(ids)" lastindex(ids)
    for i = 1:lastindex(ids)
        @debug "lastindex(ids)" lastindex(ids)
        write_csv(file, t, ids[i], abundance[i])
    end
    flush(file)
end

function update_all_rates!(sim::Simulator)
    for i = 1:length(TOP_LEVEL_EVENTS)
        sim.top_level_event_rates[i] = get_rate(TOP_LEVEL_EVENTS[i], sim)
    end
end

function update_rate!(event_id, sim::Simulator)
    sim.top_level_event_rates[event_id] = get_rate(event_id, sim)
end


### EVENT DISPATCH

function get_rate(event_id, sim::Simulator)
  if event_id == BACTERIAL_GROWTH
    get_rate_bacterial_growth(sim)
  elseif event_id == BACTERIAL_DEATH
    get_rate_bacterial_death(sim)
  elseif event_id == VIRAL_DECAY
    get_rate_viral_decay(sim)
  elseif event_id == CONTACT
    get_rate_contact(sim)
  end
end

function do_event!(event_id, sim::Simulator, t::Float64)
  if event_id == BACTERIAL_GROWTH
    do_event_bacterial_growth!(sim, t)
  elseif event_id == BACTERIAL_DEATH
    do_event_bacterial_death!(sim, t)
  elseif event_id == VIRAL_DECAY
    do_event_viral_decay!(sim, t)
  elseif event_id == CONTACT
    do_event_contact!(sim, t)
  end
end

### BACTERIAL GROWTH EVENT

function get_rate_bacterial_growth(sim::Simulator)
    p = sim.parameters
    s = sim.state

    # Birth rate has a truncated logistic form:
    # B(N) = b0 * N * (1 - N / C) [N < C]
    # B(N) = 0 [N >= C]

    # To match the Childs model, we need the birth rate at N = 0 to
    # offset the death rate to yield a total growth rate of r:
    #
    # b0 = r + d
    r = p.r_growth_rate
    d = p.d_death_rate
    b0 = r + d

    # And we need the birth rate at N = K * V to similarly equal d:
    #
    # b0 * (1 - KV/C) = d
    # =>
    # C = KV / (1 - d / b0)

    KV = p.K_carrying_capacity / p.rho_c_density_cutoff
    C = KV / (1 - d / b0)

    # Assumes total_abundance is correct
    N = s.bstrains.total_abundance

    # Birth rate is truncated to be nonnegative:
    max(0, b0 * N * (1 - N / C))
end

function do_event_bacterial_growth!(sim::Simulator, t::Float64)
    p = sim.parameters
    s = sim.state
    rng = sim.rng

    N_vec = s.bstrains.abundance
    N = s.bstrains.total_abundance

    # Choose a strain proportional to abundance
    strain_index = sample_linear_integer_weights(rng, N_vec, N)

    # Update abundance and total abundance
    s.bstrains.abundance[strain_index] += 1
    s.bstrains.total_abundance += 1
end


### BACTERIAL DEATH EVENT ###

function get_rate_bacterial_death(sim::Simulator)
    p = sim.parameters
    s = sim.state

    N = s.bstrains.total_abundance
    d = p.d_death_rate

    d * N
end

function do_event_bacterial_death!(sim::Simulator, t::Float64)
    p = sim.parameters
    s = sim.state
    rng = sim.rng

    N_vec = s.bstrains.abundance
    N = s.bstrains.total_abundance

    # Choose a strain proportional to abundance.
    # This is OK since per-capita death rate is the same across all strains.
    strain_index = sample_linear_integer_weights(rng, N_vec, N)

    # Update abundance and total abundance
    s.bstrains.abundance[strain_index] -= 1
    s.bstrains.total_abundance -= 1

    # Remove extinct strain
    if s.bstrains.abundance[strain_index] == 0
        remove_bstrain!(s.bstrains, strain_index)
    end
end


### VIRAL DECAY EVENT ###

function get_rate_viral_decay(sim::Simulator)
    p = sim.parameters
    s = sim.state
    m = p.m_viral_decay_rate
    V = s.vstrains.total_abundance

    m * V
end

function do_event_viral_decay!(sim::Simulator, t::Float64)
    p = sim.parameters
    s = sim.state
    rng = sim.rng

    V_vec = s.vstrains.abundance
    V = s.vstrains.total_abundance

    # Choose a strain proportional to abundance.
    # This is OK since per-capita death rate is the same across all strains.
    strain_index = sample_linear_integer_weights(rng, V_vec, V)

    # Update abundance and total abundance
    s.vstrains.abundance[strain_index] -= 1
    s.vstrains.total_abundance -= 1

    # Remove extinct strain
    if s.vstrains.abundance[strain_index] == 0
        remove_vstrain!(s.vstrains, strain_index)
    end
end


### CONTACT EVENT ###

function get_rate_contact(sim::Simulator)
    p = sim.parameters
    s = sim.state

    phi = p.phi_adsorption_rate
    N = s.bstrains.total_abundance
    V = s.vstrains.total_abundance
    rho = p.rho_c_density_cutoff

    # phi * (N * rho) * (V * rho) * rho
    # = phi * (N / vol) * (V /vol) * vol
    # = contacts per unit time
    phi * N * V * rho
end

function do_event_contact!(sim::Simulator, t::Float64)
    rng = sim.rng
    params = sim.parameters
    s = sim.state

    N = s.bstrains.total_abundance
    V = s.vstrains.total_abundance

    N_vec = s.bstrains.abundance
    V_vec = s.vstrains.abundance

    # Choose bacterial strain and viral strain proportional to population size
    iB = sample_linear_integer_weights(rng, N_vec, N)
    jV = sample_linear_integer_weights(rng, V_vec, V)

    # Every contact results in the reduction of the viral population size by 1
    s.vstrains.abundance[jV] -= 1
    s.vstrains.total_abundance -= 1

    should_infect = false
    should_acquire_spacer = false
    if is_immune(s.bstrains.spacers[iB], s.vstrains.pspacers[jV])
        @debug "Immune" t
        # If immune, infect anyway with probability p
        if rand(rng) < params.p_crispr_failure_prob
            should_infect = true
        end
    else
        @debug "Not immune" t
        # If not immune, defend (and acquire spacer) with probability q
        if rand(rng) < params.q_spacer_acquisition_prob
            should_acquire_spacer = true
        else
            should_infect = true
        end
    end

    if should_infect
        infect!(sim, t, iB, jV)
    elseif should_acquire_spacer
        acquire_spacer!(sim, t, iB, jV)
    end
end

function infect!(sim::Simulator, t::Float64, iB, jV)
    rng = sim.rng
    params = sim.parameters
    s = sim.state

    mu = params.mu_mutation_rate
    beta = params.beta_burst_size
    
    @debug "Infecting!" t
    # Reduce bacterial population
    @assert s.bstrains.abundance[iB] > 0
    s.bstrains.abundance[iB] -= 1
    s.bstrains.total_abundance -= 1

    # Calculate number of mutations in each virus particle
    old_pspacers = s.vstrains.pspacers[jV]
    n_pspacers = length(old_pspacers)
    
    # The number of mutations for each new virus particle is binomially distributed.
    # n_mut is left with just the nonzero draws--that is, the number of
    # mutations for each *mutated* virus particle.
    n_mut = filter(x -> x > 0, rand(rng, Binomial(n_pspacers, mu), beta))
    n_with_mut = length(n_mut)

    @debug "mutations: " n_mut n_with_mut

    # Increment population: just the virus particles that don't have mutations
    s.vstrains.abundance[jV] += beta - n_with_mut
    s.vstrains.total_abundance += beta - n_with_mut

    # Perform mutations for mutated particles
    for i = 1:n_with_mut
        # Draw which loci are mutated using the previously drawn number of mutations,
        # and create new protospacers
        mut_loci = sample(rng, 1:n_pspacers, n_mut[i]; replace=false, ordered=false)
        new_pspacers = copy(old_pspacers)
        for locus = mut_loci
            new_pspacers[locus] = s.vstrains.next_pspacer_id
            s.vstrains.next_pspacer_id += 1
        end

        @debug "Mutating virus" t mut_loci old_pspacers new_pspacers

        @debug t "creating new viral strain"
        id = s.vstrains.next_id
        s.vstrains.next_id += 1
        push!(s.vstrains.ids, id)
        push!(s.vstrains.abundance, 1)
        s.vstrains.total_abundance += 1
        push!(s.vstrains.pspacers, new_pspacers)
    
        if sim.runparams.enable_output
            write_strain(s.vstrains.strain_file, t, id, s.vstrains.ids[jV], s.bstrains.ids[iB])
            write_spacers(s.vstrains.spacers_file, id, new_pspacers)
        end
    end
end

function acquire_spacer!(sim::Simulator, t::Float64, iB, jV)
    rng = sim.rng
    params = sim.parameters
    s = sim.state
    
    @debug "Acquiring spacer!" t

    # Choose among protospacers in the infecting strain not already acquired
    # (If all have already been acquired, don't do anything.)
    missing_spacers = setdiff(s.vstrains.pspacers[jV], s.bstrains.spacers[iB])
    if length(missing_spacers) > 0
        # Create new bacterial strain with modified spacers
        s.bstrains.abundance[iB] -= 1
    
        # Add spacer, dropping the oldest one if we're at capacity
        old_spacers = s.bstrains.spacers[iB]
    
        new_spacers = if length(old_spacers) == params.u_n_spacers_max
            old_spacers[2:length(old_spacers)]
        else
            copy(old_spacers)
        end
        new_spacers_from_old = copy(new_spacers)
    
        push!(new_spacers, rand(rng, missing_spacers))

        mutated_strain_index = findfirst(x -> x == new_spacers, s.bstrains.spacers)
        if mutated_strain_index === nothing
            @debug "Creating new bacterial strain"
        
            @debug "Old spacers:" old_spacers
            @debug "New spacers from old spacers:" new_spacers_from_old
            @debug "Missing spacers:" missing_spacers
            @debug "All new spacers:" new_spacers
        
            id = s.bstrains.next_id
            s.bstrains.next_id += 1
            push!(s.bstrains.ids, id)
            push!(s.bstrains.abundance, 1)
            push!(s.bstrains.spacers, new_spacers)
        
            if sim.runparams.enable_output
                write_strain(s.bstrains.strain_file, t, id, s.bstrains.ids[iB], s.vstrains.ids[jV])
                write_spacers(s.bstrains.spacers_file, id, new_spacers)
            end
        else
            @debug "using old bacterial strain" mutated_strain_index
            s.bstrains.abundance[mutated_strain_index] += 1
        end
    end
end

function is_immune(spacers::Vector{UInt64}, pspacers::Vector{UInt64})
    length(intersect(spacers, pspacers)) > 0
end


### VALIDATION TO ENSURE CODE CORRECTNESS IN THE FACE OF OPTIMIZATIONS ###

function validate(s::State)
    validate(s.bstrains)
    validate(s.vstrains)
end

function validate(b::BStrains)
    @assert b.total_abundance == sum(b.abundance)
    @assert b.next_id > maximum(b.ids)
    @assert length(b.abundance) == length(b.ids)
    @assert length(b.spacers) == length(b.ids)
end

function validate(v::VStrains)
    @assert v.total_abundance == sum(v.abundance)
    @assert v.next_id > maximum(v.ids)
    @assert length(v.abundance) == length(v.ids)
    @assert length(v.pspacers) == length(v.ids)
end

end # module StochasticCrispr


