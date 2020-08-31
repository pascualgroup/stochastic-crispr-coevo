using Random
using Logging
using Distributions
using StatsBase
using DelimitedFiles
using Dates
using JSON2


### PARAMETERS ###

mutable struct Parameters
    "Simulation end time"
    t_final::Union{Float64, Nothing}
    
    "Time between output events"
    t_output::Union{Float64, Nothing}
    
    "Seed for random number generator"
    rng_seed::Union{UInt64, Nothing}
    
    "Enable output?"
    enable_output::Union{Bool, Nothing}
    
    function RunParameters()
        p = new()
        p.rng_seed = nothing
        p.enable_output = true
        p
    end
    "Number of initial bacterial strains"
    n_bstrains::Union{UInt64, Nothing}
    
    "Number of initial hosts per bacterial strain"
    n_hosts_per_bstrain::Union{UInt64, Nothing}
    
    # "Number of initial spacers per bacterial strain"
    # n_spacers::UInt64
    
    "Number of initial virus strains"
    n_vstrains::Union{UInt64, Nothing}
    
    "Number of initial particles per bacterial strain"
    n_particles_per_vstrain::Union{UInt64, Nothing}
    
    "Number of initial protospacers per virus strain"
    n_protospacers::Union{UInt64, Nothing}
    
    InitializationParameters() = new()
    
    "Maximum number of spacers in a bacterial strain"
    u_n_spacers_max::Union{UInt64, Nothing}
    
    "CRIPSR failure probability"
    p_crispr_failure_prob::Union{Float64, Nothing}
    
    "New spacer acquisition probability"
    q_spacer_acquisition_prob::Union{Float64, Nothing}
    
    "Growth rate at 0 (1/h)"
    r_growth_rate::Union{Float64, Nothing}
    
    "Carrying capacity (1/mL)"
    K_carrying_capacity::Union{Float64, Nothing}
    
    "Burst size"
    beta_burst_size::Union{UInt64, Nothing}
    
    "Adsorption rate"
    phi_adsorption_rate::Union{Float64, Nothing}
    
    "Viral decay rate"
    m_viral_decay_rate::Union{Float64, Nothing}
    
    "Mutation rate for contact-coupled viral mutations"
    mu1_viral_mutation_rate_coupled::Union{Float64, Nothing}
    
    "Mutation rate for contact-decoupled viral mutations"
    mu2_viral_mutation_rate_decoupled::Union{Float64, Nothing}
    
    "Density cutoff: used to scale volume of system, and therefore discrete population sizes"
    rho_c_density_cutoff::Union{Float64, Nothing}
    
    "Constant death rate (not in Childs model)"
    d_death_rate::Union{Float64, Nothing}
    
    "If true, decouple viral mutation events from contact process as a separate event"
    decouple_viral_mutation::Union{Bool, Nothing}

    function Parameters()
        p = new()
        p.rng_seed = nothing
        p.enable_output = true
        p
    end
end
JSON2.@format Parameters noargs

function validate(p::Parameters)
    @assert p.t_final !== nothing
    @assert p.t_output !== nothing
    
    @assert p.n_bstrains !== nothing
    @assert p.n_hosts_per_bstrain !== nothing
    @assert p.n_vstrains !== nothing
    @assert p.n_particles_per_vstrain !== nothing
    @assert p.n_protospacers !== nothing
    
    @assert p.u_n_spacers_max !== nothing
    @assert p.p_crispr_failure_prob !== nothing
    @assert p.q_spacer_acquisition_prob !== nothing
    @assert p.r_growth_rate !== nothing
    @assert p.K_carrying_capacity !== nothing
    @assert p.beta_burst_size !== nothing
    @assert p.phi_adsorption_rate !== nothing
    @assert p.m_viral_decay_rate !== nothing
    @assert p.rho_c_density_cutoff !== nothing
    @assert p.d_death_rate !== nothing
    
    @assert p.decouple_viral_mutation !== nothing
    @assert (p.decouple_viral_mutation || p.mu1_viral_mutation_rate_coupled !== nothing)
    @assert (!p.decouple_viral_mutation || p.mu2_viral_mutation_rate_decoupled !== nothing)
end

function load_parameters_from_json(filename) :: Parameters
    str = open(f -> read(f, String), filename)
    params = JSON2.read(str, Parameters)
    validate(params)
    params
end


### EVENT CONSTANTS ###

const N_EVENTS = 5
const EVENTS = 1:N_EVENTS

const (
    BACTERIAL_GROWTH,
    BACTERIAL_DEATH,
    VIRAL_DECAY,
    CONTACT,
    VIRAL_MUTATION,
) = EVENTS


### SIMULATION STATE ###

mutable struct Strains
    next_id::UInt64
    ids::Vector{UInt64}

    abundance::Vector{UInt64}
    total_abundance::UInt64

    spacers::Vector{Vector{UInt64}}

    strain_file::IOStream
    spacers_file::IOStream
    abundance_file::IOStream
end

function make_bstrains(n_strains, n_hosts_per_strain)
    Strains(
        n_strains + 1,
        1:n_strains,
        repeat([n_hosts_per_strain], n_strains),
        n_strains * n_hosts_per_strain,
        repeat([[]], n_strains),
        open_csv("bstrains", "t_creation", "bstrain_id", "parent_bstrain_id", "infecting_vstrain_id"),
        open_csv("bspacers", "bstrain_id", "spacer_id"),
        open_csv("babundance", "t", "bstrain_id", "abundance")
    )
end

function make_vstrains(n_strains, n_particles_per_strain, n_pspacers_init)
    next_id = n_strains + 1
    ids = Vector(1:n_strains)
    abundance = repeat([n_particles_per_strain], n_strains)
    total_abundance = n_strains * n_particles_per_strain
    pspacers = [
        Vector(1:n_pspacers_init) .+ repeat([n_pspacers_init * (i - 1)], n_pspacers_init)
        for i = 1:n_strains
    ]

    Strains(
        next_id,
        ids,
        abundance,
        total_abundance,
        pspacers,
        open_csv("vstrains", "t_creation", "vstrain_id", "parent_vstrain_id", "infected_bstrain_id"),
        open_csv("vpspacers", "vstrain_id", "spacer_id"),
        open_csv("vabundance", "t", "vstrain_id", "abundance")
    )
end

function remove_strain!(strains, index)
    # This is only used when a strain has gone extinct
    @assert strains.abundance[index] == 0

    @debug "Removing strain" id=strains.ids[index] index=index

    swap_with_end_and_remove!(strains.ids, index)
    swap_with_end_and_remove!(strains.abundance, index)
    swap_with_end_and_remove!(strains.spacers, index)
end

mutable struct State
    bstrains::Strains
    vstrains::Strains
    next_pspacer_id::UInt64

    function State(
        n_bstrains, n_hosts_per_bstrain, # n_spacers_init,
        n_vstrains, n_particles_per_vstrain, n_pspacers_init
    )
        next_pspacer_id = 1 + n_vstrains * n_pspacers_init
        new(
            make_bstrains(n_bstrains, n_hosts_per_bstrain),
            make_vstrains(n_vstrains, n_particles_per_vstrain, n_pspacers_init),
            next_pspacer_id
        )
    end
end

function State(p::Parameters)
    State(
        p.n_bstrains, p.n_hosts_per_bstrain,
        p.n_vstrains, p.n_particles_per_vstrain, p.n_protospacers
    )
end

mutable struct Simulation
    params::Parameters
    t::Float64
    state::State
    rng::MersenneTwister

    event_rates::Vector{Float64}
    event_counts::Vector{UInt64}
    
    meta_file::IOStream
    summary_file::IOStream
    
    function Simulation(p::Parameters)
        meta_file = open_csv("meta", "key", "value")

        # Use random seed if provided, or generate one
        rng_seed = p.rng_seed === nothing ? UInt64(rand(RandomDevice(), UInt32)) : p.rng_seed
        p.rng_seed = rng_seed
        write_csv(meta_file, "rng_seed", rng_seed)
        
        # Save parameters as loaded
        write_json_to_file(p, "parameters_out.json")

        # Record start time
        start_time = now()
        write_csv(meta_file, "start_time", start_time)

        # Initialize & validate model state
        state = State(p)
        validate(state)

        sim = new(
            p, 0.0, state, MersenneTwister(rng_seed),
            zeros(length(EVENTS)), zeros(length(EVENTS)),
            meta_file,
            open_csv("summary", "t", "bacterial_abundance", "viral_abundance")
        )
        update_rates!(sim)
        sim
    end
end



### SIMULATION LOOP ###

function simulate(sim::Simulation)
    state = sim.state
    p = sim.params
    
    # Initial output
    if p.enable_output
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
            state.vstrains.spacers
        )
    end

    # Simulation loop
    t_next_output = 0.0

    while sim.t < p.t_final
        # Simulate exactly until the next output time
        t_next_output = min(p.t_final, t_next_output + p.t_output)

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
        @debug "vstrains:" total_abund=state.vstrains.total_abundance abund=state.vstrains.abundance pspacers=state.vstrains.spacers

        # Write periodic output
        @debug "p.enable_output" p.enable_output
        if p.enable_output
            write_periodic_output(sim)
        end

        @assert sim.t == t_next_output
    end

    # Record end time and elapsed
    end_time = now()
    write_csv(sim.meta_file, "end_time", end_time)
    elapsed_seconds = Dates.value(end_time - start_time) / 1000.0
    write_csv(sim.meta_file, "elapsed_seconds", elapsed_seconds)

    close(sim.meta_file)
end

function do_next_event!(sim::Simulation, t_max::Float64)
    p = sim.params
    s = sim.state

    @assert length(sim.event_rates) == length(EVENTS)
    R = sum(sim.event_rates)

    # Draw next event time using total rate
    t_next = sim.t + randexp(sim.rng) / R

    if t_next > t_max
        sim.t = t_max
    else
        @debug "event_rates:" sim.event_rates

        # Sample next top-level event proportional to event rate
        event_id = sample(sim.rng, EVENTS, Weights(sim.event_rates, R))
        sim.event_counts[event_id] += 1

        @debug "begin do_event()" event=event t=t_next
        do_event!(event_id, sim, t_next)
        update_rates!(sim)
        @debug "end do_event()"

        sim.t = t_next

        @debug "bstrains.total_abundance:" sim.state.bstrains.total_abundance
        @debug "VStrains.total_abundance:" sim.state.vstrains.total_abundance
    end
end

function update_rates!(sim::Simulation)
    for i = EVENTS
        sim.event_rates[i] = get_rate(i, sim)
    end
end


### EVENT DISPATCH ###

# This used to be done using Val/multiple dispatch,
# but this is easier to understand for Julia newbies.

function get_rate(event_id, sim::Simulation)
    if event_id == BACTERIAL_GROWTH
        get_rate_bacterial_growth(sim)
    elseif event_id == BACTERIAL_DEATH
        get_rate_bacterial_death(sim)
    elseif event_id == VIRAL_DECAY
        get_rate_viral_decay(sim)
    elseif event_id == CONTACT
        get_rate_contact(sim)
    elseif event_id == VIRAL_MUTATION
        get_rate_viral_mutation(sim)
    else
        error("unknown event")
    end
end

function do_event!(event_id, sim::Simulation, t::Float64)
    if event_id == BACTERIAL_GROWTH
        do_event_bacterial_growth!(sim, t)
    elseif event_id == BACTERIAL_DEATH
        do_event_bacterial_death!(sim, t)
    elseif event_id == VIRAL_DECAY
        do_event_viral_decay!(sim, t)
    elseif event_id == CONTACT
        do_event_contact!(sim, t)
    elseif event_id == VIRAL_MUTATION
        do_event_viral_mutation!(sim, t)
    else
        error("unknown event")
    end
end


### BACTERIAL GROWTH EVENT ###

function get_rate_bacterial_growth(sim::Simulation)
    p = sim.params
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

function do_event_bacterial_growth!(sim::Simulation, t::Float64)
    p = sim.params
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

function get_rate_bacterial_death(sim::Simulation)
    p = sim.params
    s = sim.state

    N = s.bstrains.total_abundance
    d = p.d_death_rate

    d * N
end

function do_event_bacterial_death!(sim::Simulation, t::Float64)
    p = sim.params
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
        remove_strain!(s.bstrains, strain_index)
    end
end


### VIRAL DECAY EVENT ###

function get_rate_viral_decay(sim::Simulation)
    p = sim.params
    s = sim.state
    m = p.m_viral_decay_rate
    V = s.vstrains.total_abundance

    m * V
end

function do_event_viral_decay!(sim::Simulation, t::Float64)
    p = sim.params
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
        remove_strain!(s.vstrains, strain_index)
    end
end


### CONTACT EVENT ###

function get_rate_contact(sim::Simulation)
    p = sim.params
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

function do_event_contact!(sim::Simulation, t::Float64)
    rng = sim.rng
    params = sim.params
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
    if is_immune(s.bstrains.spacers[iB], s.vstrains.spacers[jV])
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

function infect!(sim::Simulation, t::Float64, iB, jV)
    rng = sim.rng
    p = sim.params
    s = sim.state
    
    beta = p.beta_burst_size

    @debug "Infecting!" t
    # Reduce bacterial population
    @assert s.bstrains.abundance[iB] > 0
    s.bstrains.abundance[iB] -= 1
    s.bstrains.total_abundance -= 1

    # Calculate number of mutations in each virus particle
    old_pspacers = s.vstrains.spacers[jV]
    n_pspacers = length(old_pspacers)
    
    # Just adjust viral population upward with the burst size
    s.vstrains.abundance[jV] += beta
    s.vstrains.total_abundance += beta
    
    # Perform mutations if we're in coupled mutation mode
    if !p.decouple_viral_mutation
        mu1 = p.mu1_viral_mutation_rate_coupled
        
        # The number of mutations for each new virus particle is binomially distributed.
        # n_mut is left with just the nonzero draws--that is, the number of
        # mutations for each *mutated* virus particle.
        n_mut = filter(x -> x > 0, rand(rng, Binomial(n_pspacers, mu1), beta))
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
            mutate_virus!(sim, jV, mut_loci, iB)
        end
    end
end

function acquire_spacer!(sim::Simulation, t::Float64, iB, jV)
    rng = sim.rng
    p = sim.params
    s = sim.state

    @debug "Acquiring spacer!" t

    # Choose among protospacers in the infecting strain not already acquired
    # (If all have already been acquired, don't do anything.)
    missing_spacers = setdiff(s.vstrains.spacers[jV], s.bstrains.spacers[iB])
    if length(missing_spacers) > 0
        # Create new bacterial strain with modified spacers
        s.bstrains.abundance[iB] -= 1

        # Add spacer, dropping the oldest one if we're at capacity
        old_spacers = s.bstrains.spacers[iB]

        new_spacers = if length(old_spacers) == p.u_n_spacers_max
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

            if p.enable_output
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


### VIRAL MUTATION EVENT (used only if decouple_viral_mutation == true) ###

function get_rate_viral_mutation(sim::Simulation)
    p = sim.params
    if p.decouple_viral_mutation
        s = sim.state
        mu2 = p.mu2_viral_mutation_rate_decoupled
        V = s.vstrains.total_abundance
        
        # NB: mu2 rate is per-protospacer
        mu2 * p.n_protospacers * V
    else
        0.0
    end
end

function do_event_viral_mutation!(sim::Simulation, t::Float64)
    @assert sim.params.decouple_viral_mutation
    
    rng = sim.rng
    p = sim.params
    s = sim.state
    
    V = s.vstrains.total_abundance

    N_vec = s.bstrains.abundance
    V_vec = s.vstrains.abundance

    jV = sample_linear_integer_weights(rng, V_vec, V)

    old_pspacers = s.vstrains.spacers[jV]
    n_pspacers = length(old_pspacers)
    mut_loci = [rand(rng, 1:n_pspacers)]
    
    mutate_virus!(sim, jV, mut_loci, 0)
end


### VIRAL MUTATION FUNCTION (shared by contact-coupled and contact-decoupled mutation) ###

function mutate_virus!(sim, virus_id, mut_loci, contact_b_id)
    p = sim.params
    s = sim.state
    rng = sim.rng
    
    s.vstrains.abundance[virus_id] -= 1
    s.vstrains.total_abundance -= 1
    
    old_pspacers = s.vstrains.spacers[virus_id]
    new_pspacers = copy(old_pspacers)
    for locus = mut_loci
        new_pspacers[locus] = s.next_pspacer_id
        s.next_pspacer_id += 1
    end
    
    @debug "Mutating virus" t mut_loci old_pspacers new_pspacers
    
    id = s.vstrains.next_id
    s.vstrains.next_id += 1
    push!(s.vstrains.ids, id)
    push!(s.vstrains.abundance, 1)
    s.vstrains.total_abundance += 1
    push!(s.vstrains.spacers, new_pspacers)
    
    if p.enable_output
        write_strain(s.vstrains.strain_file, sim.t, id, s.vstrains.ids[virus_id], contact_b_id)
        write_spacers(s.vstrains.spacers_file, id, new_pspacers)
    end
end


### VALIDATION TO ENSURE CODE CORRECTNESS IN THE FACE OF OPTIMIZATIONS ###

function validate(s::State)
    validate(s.bstrains)
    validate(s.vstrains)
end

function validate(strains::Strains)
    @assert strains.total_abundance == sum(strains.abundance)
    @assert strains.next_id > maximum(strains.ids)
    @assert length(strains.abundance) == length(strains.ids)
    @assert length(strains.spacers) == length(strains.ids)
end


### RNG UTILITY FUNCTIONS ##

"""
Returns an index from 1:length(w) with probability proportional to w.

s must be precomputed to be sum(w).
"""
function sample_linear_integer_weights(rng::MersenneTwister, w::Vector{UInt64}, s::UInt64)
    i = rand(rng, 1:s)
    cs = 0
    for j = 1:(length(w) - 1)
        cs += w[j]
        if i <= cs
            return j
        end
    end
    length(w)
end

"""
Removes an item in the middle of an array that does not need to be kept ordered in constant time.

The item is replaced with the item at the end of the array, and then the item at the end of the
array is removed.
"""
function swap_with_end_and_remove!(a, index)
    if index != lastindex(a)
        setindex!(a, a[lastindex(a)], index)
    end
    pop!(a)
    nothing
end


### OUTPUT FUNCTIONS ###

function write_json_to_file(x, filename)
    if ispath(filename)
        error("$filename already exists. You should delete output, or run in a different directory.")
    else
        file = open(filename, "w")
        print(file, JSON2.write(x))
        file
    end
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

