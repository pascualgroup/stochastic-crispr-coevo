### EVENT CONSTANTS ###

const N_EVENTS = 5
const EVENTS = 1:N_EVENTS

const (
MICROBIAL_GROWTH,
MICROBIAL_DEATH,
VIRAL_DECAY,
CONTACT,
MICROBIAL_IMMIGRATION
) = EVENTS



############### SIMULATION LOOP ################

function simulate(sim::Simulation)
state = sim.state
p = sim.params
db = sim.db

execute(db, "BEGIN TRANSACTION")
# Initial output
if p.enable_output
    @info "initial output"
    write_periodic_output(sim)

    write_strains(sim, "bstrains", "bspacers", "bgrowthrates",
                    state.bstrains.ids, state.bstrains.spacers,
                    state.bstrains.growthalleles, state.bstrains.growthrates)
    write_strains(sim, "vstrains", "vpspacers", nothing,
                    state.vstrains.ids, state.vstrains.spacers,
                    nothing, nothing)
    # println("checkpoint1")
end

# Simulation loop
t_next_output = 0.0
while sim.t < p.t_final
    # Simulate exactly until the next output time
    t_next_output = min(p.t_final, t_next_output + p.t_output)

    beginTime = sim.t
    while sim.t < t_next_output # This continues until sum of event times meet sampling time
        do_next_event!(sim, t_next_output)
        # println("oh")
    end
    endTime = sim.t

    @info "Completed period: $(beginTime) to $(endTime)"

    @assert sim.t >= t_next_output
end
execute(db, "COMMIT")

# Record end time and elapsed
end_time = now()
elapsed_seconds = Dates.value(end_time - start_time) / 1000.0

execute(db, "BEGIN TRANSACTION")
execute(db,
"INSERT INTO meta VALUES (?,?)",
["end_time", Dates.format(end_time, "yyyy-mm-ddTHH:MM:SS")]
)

execute(db,
"INSERT INTO meta VALUES (?,?)",
["elapsed_seconds", elapsed_seconds]
)
execute(db, "COMMIT")

end




################## EVENT FUNCTIONS ##################

function do_next_event!(sim::Simulation, t_max::Float64)
    p = sim.params
    s = sim.state
    db = sim.db

    @assert length(sim.event_rates) == length(EVENTS)
    R = sum(sim.event_rates)

    # Draw next event time using total rate
    t_next = sim.t + randexp(sim.rng) / R

    if t_next > t_max
        if t_next > p.t_final

            sim.t = p.t_final

                # Write periodic output
                @debug "p.enable_output" p.enable_output
                if p.enable_output
                    write_periodic_output(sim)
                end

                @debug "bstrains:" total_abund=s.bstrains.total_abundance abund=s.bstrains.abundance spacers=s.bstrains.spacers
                @debug "vstrains:" total_abund=s.vstrains.total_abundance abund=s.vstrains.abundance pspacers=s.vstrains.spacers

        else
            sim.t = t_max

                # Write periodic output
                @debug "p.enable_output" p.enable_output
                if p.enable_output
                    write_periodic_output(sim)
                end

                @debug "bstrains:" total_abund=s.bstrains.total_abundance abund=s.bstrains.abundance spacers=s.bstrains.spacers
                @debug "vstrains:" total_abund=s.vstrains.total_abundance abund=s.vstrains.abundance pspacers=s.vstrains.spacers

            # Time advances to value greater than current output time

                # Sample next top-level event proportional to event rate
                total = sum(sim.event_counts)
                @debug "event counts:" total=total, breakdown=sim.event_counts
                @debug "event_rates:" rates=sim.event_rates

                event_id = sample(sim.rng, EVENTS, Weights(sim.event_rates, R))
                sim.event_counts[event_id] += 1

                @debug "begin do_event()" event_id=event_id t=t_next
                do_event!(event_id, sim, t_next)
                update_rates!(sim)
                @debug "end do_event()"

            sim.t = t_next

                @debug "bstrains:" total_abund=s.bstrains.total_abundance abund=s.bstrains.abundance spacers=s.bstrains.spacers
                @debug "vstrains:" total_abund=s.vstrains.total_abundance abund=s.vstrains.abundance pspacers=s.vstrains.spacers
        end
    else
                total = sum(sim.event_counts)
                @debug "event counts:" total=total, breakdown=sim.event_counts
                @debug "event_rates:" rates=sim.event_rates

                # Sample next top-level event proportional to event rate
                event_id = sample(sim.rng, EVENTS, Weights(sim.event_rates, R))
                sim.event_counts[event_id] += 1

                @debug "begin do_event()" event_id=event_id t=t_next
                do_event!(event_id, sim, t_next)
                update_rates!(sim)
                @debug "end do_event()"

            sim.t = t_next # Time advances

                @debug "bstrains:" total_abund=s.bstrains.total_abundance abund=s.bstrains.abundance spacers=s.bstrains.spacers
                @debug "vstrains:" total_abund=s.vstrains.total_abundance abund=s.vstrains.abundance pspacers=s.vstrains.spacers
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
    if event_id == MICROBIAL_GROWTH
        get_rate_microbial_growth(sim)
    elseif event_id == MICROBIAL_DEATH
        get_rate_microbial_death(sim)
    elseif event_id == VIRAL_DECAY
        get_rate_viral_decay(sim)
    elseif event_id == CONTACT
        get_rate_contact(sim)
    elseif event_id == MICROBIAL_IMMIGRATION
        get_rate_microbial_immigration(sim)
    else
        error("unknown event")
    end
end


function do_event!(event_id, sim::Simulation, t::Float64)
    if event_id == MICROBIAL_GROWTH
        do_event_microbial_growth!(sim, t)
    elseif event_id == MICROBIAL_DEATH
        do_event_microbial_death!(sim, t)
    elseif event_id == VIRAL_DECAY
        do_event_viral_decay!(sim, t)
    elseif event_id == CONTACT
        do_event_contact!(sim, t)
    elseif event_id == MICROBIAL_IMMIGRATION
        do_event_microbial_immigration!(sim, t)
    else
        error("unknown event")
    end
end


### MICROBIAL GROWTH EVENT ###

function get_rate_microbial_growth(sim::Simulation)
    p = sim.params
    s = sim.state
    C = p.microbe_carrying_capacity

    N_vec = s.bstrains.abundance
    R_vec = s.bstrains.growthrates
    # Assumes total_abundance is correct
    N = s.bstrains.total_abundance
    # Birth rate is truncated to be nonnegative:
    # if sim.t > 15
    #     println("R_vec: $(R_vec), N_vec: $(N_vec), spacers: $(s.bstrains.spacers)")
    # end
    max(0, sum((R_vec .* N_vec) * (1 - N / C)))
end


function do_event_microbial_growth!(sim::Simulation, t::Float64)
    p = sim.params
    s = sim.state
    rng = sim.rng

    N_vec = s.bstrains.abundance
    R_vec = s.bstrains.growthrates
    N = s.bstrains.total_abundance

    # Choose a strain proportional to abundance
    # strain_index = sample_linear_integer_weights(rng, N_vec, N)
    strain_index = sample(rng,Weights(N_vec.*R_vec))
    # Update abundance and total abundance
     s.bstrains.abundance[strain_index] += 1
    s.bstrains.total_abundance += 1
end


### MICROBIAL DEATH EVENT ###

function get_rate_microbial_death(sim::Simulation)
    p = sim.params
    s = sim.state

    N = s.bstrains.total_abundance
    d = p.microbe_death_rate

    d * N
end


function do_event_microbial_death!(sim::Simulation, t::Float64)
    p = sim.params
    s = sim.state
    rng = sim.rng

    N_vec = s.bstrains.abundance
    N = s.bstrains.total_abundance

    # Choose a strain proportional to abundance.
    # This is OK since per-capita death rate is the same across all strains.
    strain_index = sample_linear_integer_weights(rng, N_vec, N)

    # Update abundance and total abundance
    @assert s.bstrains.abundance[strain_index] > 0
    s.bstrains.abundance[strain_index] -= 1
    s.bstrains.total_abundance -= 1

    # Remove extinct strain but keep first index that is reserved
    ## for memoryless immigrants
    if s.bstrains.abundance[strain_index] == 0 && strain_index != 1
        remove_strain!(s.bstrains, strain_index)
    end

end

### MICROBIAL IMMIGRATION EVENT ###

function get_rate_microbial_immigration(sim::Simulation)
    p = sim.params
    s = sim.state
    g = p.microbe_immigration_rate
    C = p.microbe_carrying_capacity
    # Assumes total_abundance is correct
    N = s.bstrains.total_abundance

    # Immigration rate is truncated to be nonnegative:
    max(0, g * (1 - N / C))
end


function do_event_microbial_immigration!(sim::Simulation, t::Float64)
    s = sim.state
    s.bstrains.abundance[1] += 1
    s.bstrains.total_abundance += 1
end


### VIRAL DECAY EVENT ###

function get_rate_viral_decay(sim::Simulation)
    p = sim.params
    s = sim.state
    m = p.viral_decay_rate
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
    @assert s.vstrains.abundance[strain_index] > 0
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

    phi = p.adsorption_rate
    N = s.bstrains.total_abundance
    V = s.vstrains.total_abundance

    phi * N * V
end


function do_event_contact!(sim::Simulation, t::Float64)
    rng = sim.rng
    params = sim.params
    s = sim.state

    N = s.bstrains.total_abundance
    V = s.vstrains.total_abundance

    N_vec = s.bstrains.abundance
    V_vec = s.vstrains.abundance

    # Choose microbial strain and viral strain proportional to population size
    iB = sample_linear_integer_weights(rng, N_vec, N)
    jV = sample_linear_integer_weights(rng, V_vec, V)

    # Every contact results in the reduction of the viral population size by 1
    @assert s.vstrains.abundance[jV] > 0
    s.vstrains.abundance[jV] -= 1
    s.vstrains.total_abundance -= 1

    should_infect = false
    should_acquire_spacer = false
    if is_immune(s.bstrains.spacers[iB], s.vstrains.spacers[jV])
        @debug "Immune" t = sim.t
        # If immune, infect anyway with probability p
        if params.crispr_failure_prob > 0
            if rand(rng) < params.crispr_failure_prob
                should_infect = true
            end
        end
    else
        @debug "Not immune" t = sim.t
        # If not immune, defend (and acquire spacer) with probability q
        if params.spacer_acquisition_prob > 0
            if rand(rng) < params.spacer_acquisition_prob
                should_acquire_spacer = true
            else
                should_infect = true
            end
        else
            should_infect = true
        end
    end

    if should_infect
        infect!(sim, t, iB, jV)
    elseif should_acquire_spacer
        acquire_spacer!(sim, t, iB, jV)
    end

    # Remove extinct strain
    if s.vstrains.abundance[jV] == 0
        remove_strain!(s.vstrains, jV)
    end

    # Remove extinct strain but keep first index that is reserved
    ## for memoryless immigrants
    if s.bstrains.abundance[iB] == 0 && iB != 1
        remove_strain!(s.bstrains, iB)
    end

    @assert s.bstrains.abundance[1] >= 0
    ###############
end


function infect!(sim::Simulation, t::Float64, iB, jV)
    rng = sim.rng
    p = sim.params
    s = sim.state

    beta = p.viral_burst_size

    @debug "Infecting!" t = sim.t
    # Reduce microbial population
    @assert s.bstrains.abundance[iB] > 0
    s.bstrains.abundance[iB] -= 1
    s.bstrains.total_abundance -= 1

    # Calculate number of mutations in each virus particle
    old_pspacers = s.vstrains.spacers[jV]
    n_pspacers = length(old_pspacers)

    # Just adjust viral population upward with the burst size
    #s.vstrains.abundance[jV] += beta    This is redundant considering the code below
    #s.vstrains.total_abundance += beta

    # Perform mutations
    mu = p.viral_mutation_rate

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
        #loops through n_with_mut mutated viruses among the beta: makes and adds new identity and adds abundance
        # Draw which loci are mutated using the previously drawn number of mutations,
        # and create new protospacers
        mut_loci = sample(rng, 1:n_pspacers, n_mut[i]; replace=false, ordered=false)
        mutate_virus!(sim, jV, mut_loci, iB)
    end
end


function acquire_spacer!(sim::Simulation, t::Float64, iB, jV)
    rng = sim.rng
    p = sim.params
    s = sim.state

    @debug "Acquiring spacer!" t = sim.t
    # Choose among protospacers in the infecting strain not already acquired
    # (If all have already been acquired, don't do anything.)
    missing_spacers = setdiff(s.vstrains.spacers[jV], s.bstrains.spacers[iB])
    if length(missing_spacers) > 0
        # Create new microbial strain with modified spacers
        @assert s.bstrains.abundance[iB] > 0
        s.bstrains.abundance[iB] -= 1 #This does not need a change of total abundance alongside. We are just rearranging the identities
        # Add spacer, dropping the oldest one if we're at capacity
        old_spacers = s.bstrains.spacers[iB]
        new_spacers = if length(old_spacers) == p.n_spacers_max
            old_spacers[2:length(old_spacers)]        ### This removes thes first locus
        else                                          ### of the CRISPR cassette. Is this what original C code was?
            copy(old_spacers)
        end
        new_spacers_from_old = copy(new_spacers)
        push!(new_spacers, rand(rng, missing_spacers))
        mutatedSpacer_strain_index = findall(x -> x == new_spacers, s.bstrains.spacers)

        should_mutate = false
        if p.microbe_mutation_prob > 0
            if rand(rng) < p.microbe_mutation_prob
                should_mutate = true
            end
        end
        if should_mutate
            if p.directional
                if rand(rng) < 0.5
                    sign = -1
                else
                    sign = 1
                end
                newallele = s.bstrains.growthalleles[iB] + sign * p.allelic_change
                if newallele > p.max_allele || newallele < 2 * (p.center_allele) - p.max_allele
                    newallele = s.bstrains.growthalleles[iB] - sign * p.allelic_change
                end
            else
                newallele = rand(rng, Uniform(p.center_allele - p.max_allele, p.max_allele))
            end
        end
        if should_mutate
            allele = newallele
        else
            allele = s.bstrains.growthalleles[iB]
        end
        mutatedGrowth_strain_index = findall(x -> x == allele, s.bstrains.growthalleles)
        overlap = intersect(mutatedSpacer_strain_index, mutatedGrowth_strain_index)
        @assert length(overlap) < 2 "strain index: $(iB),
                                    strain id: $(s.bstrains.ids[iB]),
                                    spacers: $(old_spacers),
                                    directional: $(p.directional),
                                    allelic_change: $(p.allelic_change),
                                    allele_max: $(p.max_allele),
                                    newallele: $(allele),
                                    new_spacers: $(new_spacers),
                                    overlap: $(overlap),
                                    overlap strain ids: $(map(x->s.bstrains.ids[x],overlap)),
                                    spacers: $(map(x->s.bstrains.spacers[x],overlap)),
                                    growthalleles: $(map(x->s.bstrains.growthalleles[x],overlap))"

        mutated_strain_index = nothing
        if length(intersect(mutatedSpacer_strain_index, mutatedGrowth_strain_index)) == 1
            mutated_strain_index = intersect(mutatedSpacer_strain_index, mutatedGrowth_strain_index)
        end
        if mutated_strain_index == nothing
            newgrowthrate = fitness(allele, p.center_allele, p.max_allele, p.max_fitness, p.evofunction)
            @debug "Creating new microbial strain"
            @debug "Old spacers:" old_spacers
            @debug "New spacers from old spacers:" new_spacers_from_old
            @debug "Missing spacers:" missing_spacers
            @debug "All new spacers:" new_spacers
            @debug "Old growth rate:" s.bstrains.growthrates[strain_index]
            @debug "New growth rate:" newgrowthrate
            id = s.bstrains.next_id
            s.bstrains.next_id += 1
            push!(s.bstrains.ids, id) #Adds new strain identity
            push!(s.bstrains.abundance, 1) #Adds abundance of one to new strain identity
            push!(s.bstrains.spacers, new_spacers)
            push!(s.bstrains.growthalleles, allele)
            push!(s.bstrains.growthrates, newgrowthrate)
            # println(s.bstrains.growthrates)
            if p.enable_output
                write_strain(sim, "bstrains", id, s.bstrains.ids[iB], s.vstrains.ids[jV])
                write_spacers(sim, "bspacers", id, new_spacers)
                write_growth_rate(sim, "bgrowthrates", id, allele, newgrowthrate)
            end
        else
            @debug "using old microbial strain" mutated_strain_index
            s.bstrains.abundance[mutated_strain_index...] += 1
        end
    end
end


function is_immune(spacers::Vector{UInt64}, pspacers::Vector{UInt64})
    length(intersect(spacers, pspacers)) > 0
end


function mutate_virus!(sim, virus_id, mut_loci, contact_b_id)
    p = sim.params
    s = sim.state
    rng = sim.rng
    #s.vstrains.abundance[virus_id] -= 1 #This is redundant considering
    #s.vstrains.total_abundance -= 1    #how contact already removes one and only
                                        #adds beta - n_with_mut
    old_pspacers = s.vstrains.spacers[virus_id]
    new_pspacers = copy(old_pspacers)
    for locus = mut_loci
        new_pspacers[locus] = s.next_pspacer_id
        s.next_pspacer_id += 1
    end

    @debug "Mutating virus" t mut_loci old_pspacers new_pspacers

    id = s.vstrains.next_id
    s.vstrains.next_id += 1
    push!(s.vstrains.ids, id) #Adds new strain identity
    push!(s.vstrains.abundance, 1) #Adds abundance of one to new strain identity
    s.vstrains.total_abundance += 1
    push!(s.vstrains.spacers, new_pspacers)

    if p.enable_output
        write_strain(sim, "vstrains", id, s.vstrains.ids[virus_id], s.bstrains.ids[contact_b_id])
        write_spacers(sim, "vpspacers", id, new_pspacers)
    end
end
