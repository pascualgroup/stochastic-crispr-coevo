module StochasticCrispr

using Util
using Random
using Logging
using StatsBase
using DataStructures

export InitializationParameters, Parameters, validate, State, Simulator, do_next_event!

mutable struct InitializationParameters
    "Number of initial bacterial strains"
    n_bstrains::UInt64

    "Number of initial hosts per bacterial strain"
    n_hosts_per_bstrain::UInt64

    "Number of initial spacers per bacterial strain"
    n_spacers::UInt64

    "Number of initial virus strains"
    n_vstrains::UInt64

    "Number of initial particles per bacterial strain"
    n_particles_per_vstrain::UInt64

    "Number of initial protospacers per virus strain"
    n_protospacers::UInt64

    InitializationParameters() = new()
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

function validate(p::Parameters)
    @assert 0.0 <= p.p_crispr_failure_prob <= 1.0
    @assert 0.0 <= p.q_spacer_acquisition_prob <= 1.0
    @assert 0.1 < p.r_growth_rate <= 10
    @assert 0.0 < p.K_carrying_capacity <= 1e7
    @assert 0 <= p.beta_burst_size <= 1000
    @assert 0.0 <= p.phi_adsorption_rate <= 1e-6
    @assert 0.0 <= p.m_viral_decay_rate <= 1.0
    @assert 0.0 <= p.mu_mutation_rate <= 1e-5
    @assert 0.05 <= p.rho_c_density_cutoff <= 1.0
end

mutable struct BStrains
    next_id::UInt64
    ids::Vector{UInt64}

    abundance::Vector{UInt64}
    total_abundance::UInt64

    spacers::Vector{Vector{UInt64}}

    function BStrains(n_strains, n_hosts_per_strain, n_spacers_init, spacer_id_start)
        new(
            n_strains + 1, 1:n_strains,
            repeat([n_hosts_per_strain], n_strains), n_strains * n_hosts_per_strain,
            fill_id_buffers(n_strains, n_spacers_init, spacer_id_start)
        )
    end
end

mutable struct VStrains
    next_id::UInt64
    ids::Vector{UInt64}

    abundance::Vector{UInt64}
    total_abundance::UInt64

    next_pspacer_id::UInt64
    pspacers::Vector{Vector{UInt64}}

    function VStrains(n_strains, n_particles_per_strain, n_pspacers_init, pspacer_id_start)
        new(
            n_strains + 1, 1:n_strains,
            repeat([n_particles_per_strain], n_strains), n_strains * n_particles_per_strain,
            pspacer_id_start + n_strains * n_pspacers_init,
            fill_id_buffers(n_strains, n_pspacers_init, pspacer_id_start)
        )
    end
end

"""
Create n_buffers circular buffers with specified maximum capacity, filled with sequential ID

E.g.:

[1,2,3,4]
[5,6,7,8]
[9,10,11,12]
"""
function fill_id_buffers(n_buffers, n_init, id_start)
    @info "fill_id_buffers():" n_buffers n_init id_start
    buffers = Vector{Vector{UInt64}}()
    for i::UInt64 = 1:n_buffers
        push!(buffers, Vector{UInt64}())
        append!(
            buffers[i],
            Vector(1:n_init) .+ repeat([id_start - 1 + n_init * (i - 1)], n_init)
        )
    end
    buffers
end

mutable struct State
    bstrains::BStrains
    vstrains::VStrains

    function State(
        n_bstrains, n_hosts_per_bstrain, n_spacers_init,
        n_vstrains, n_particles_per_vstrain, n_pspacers_init
    )
        new(
            BStrains(
                n_bstrains, n_hosts_per_bstrain,
                n_spacers_init, 1
            ),
            VStrains(
                n_vstrains, n_particles_per_vstrain,
                n_pspacers_init, 1 + n_bstrains * n_spacers_init
            )
        )
    end
end

function State(ip::InitializationParameters, mp::Parameters)
    State(
        ip.n_bstrains, ip.n_hosts_per_bstrain, ip.n_spacers,
        ip.n_vstrains, ip.n_particles_per_vstrain, ip.n_protospacers
    )
end


"""
    List of top-level events. Each of these are instances of a different *type*, which allows
    dispatch to call a different function for each value.

    That is, Val(:BacterialGrowth) is an instance of type Val{:BacterialGrowth},
    and Val(:ViralExtinction) is an instance of type Val{:ViralExtinction}.

    Therefore, if the chosen event is Val(:BacterialGrowth), the do_event function that will be
    called is the one with signature

    function do_event(event::Val{:BacterialGrowth})

    See:

    https://docs.julialang.org/en/v1/manual/types/#Parametric-Types-1
    https://docs.julialang.org/en/v1/manual/types/#"Value-types"-1
"""
const TOP_LEVEL_EVENTS = [
    Val(:BacterialGrowth),
    Val(:BacterialDeath),
    Val(:ViralDecay),
    Val(:ViralMutation),
    Val(:Contact),
]

mutable struct Simulator
    parameters::Parameters
    t::Float64
    state::State
    rng::MersenneTwister

    top_level_event_rates::Vector{Float64}
    event_counts::Vector{UInt64}
end

function Simulator(params::Parameters, t_init::Float64, state_init::State, rng::MersenneTwister)
    sim = Simulator(
        params, t_init, state_init, rng,
        zeros(length(TOP_LEVEL_EVENTS)), zeros(length(TOP_LEVEL_EVENTS))
    )
    for e = TOP_LEVEL_EVENTS
        update_rate!(e, sim)
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
        event = sample(sim.rng, TOP_LEVEL_EVENTS, Weights(sim.top_level_event_rates, R))
        sim.event_counts[findfirst(x -> x == event, TOP_LEVEL_EVENTS)] += 1

        @debug "begin do_event()" event=event t=t_next
        do_event!(event, sim, t_next)
        update_all_rates!(sim)
        @debug "end do_event()"

        sim.t = t_next

        @debug "bstrains.total_abundance:" sim.state.bstrains.total_abundance
        @debug "VStrains.total_abundance:" sim.state.vstrains.total_abundance
    end
end

function update_all_rates!(sim::Simulator)
    for i = 1:length(TOP_LEVEL_EVENTS)
        sim.top_level_event_rates[i] = get_rate(TOP_LEVEL_EVENTS[i], sim)
    end
end

function update_rate!(e::Val, sim::Simulator)
    event_index = findfirst(x -> x == e, TOP_LEVEL_EVENTS)
    @debug "event_index:" event_index
    sim.top_level_event_rates[event_index] = get_rate(e, sim)
end


### BACTERIAL GROWTH EVENT

function get_rate(e::Val{:BacterialGrowth}, sim::Simulator)
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

function do_event!(e::Val{:BacterialGrowth}, sim::Simulator, t::Float64)
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

    # Update affected rates
#     update_rate!(Val(:BacterialGrowth), sim)
#     update_rate!(Val(:BacterialDeath), sim)
#     update_rate!(Val(:Contact), sim)
end


### BACTERIAL DEATH EVENT ###

function get_rate(e::Val{:BacterialDeath}, sim::Simulator)
    p = sim.parameters
    s = sim.state

    N = s.bstrains.total_abundance
    d = p.d_death_rate

    d * N
end

function do_event!(e::Val{:BacterialDeath}, sim::Simulator, t::Float64)
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

    # Update affected rates
#     update_rate!(Val(:BacterialGrowth), sim)
#     update_rate!(Val(:BacterialDeath), sim)
#     update_rate!(Val(:Contact), sim)
end


### VIRAL DECAY EVENT ###

function get_rate(e::Val{:ViralDecay}, sim::Simulator)
    p = sim.parameters
    s = sim.state
    m = p.m_viral_decay_rate
    V = s.vstrains.total_abundance

    m * V
end

function do_event!(e::Val{:ViralDecay}, sim::Simulator, t::Float64)
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
end


### VIRAL MUTATION EVENT ###

function get_rate(e::Val{:ViralMutation}, sim::Simulator)
    0.0
end

function do_event!(e::Val{:ViralMutation}, sim::Simulator, t::Float64)
    p = sim.parameters
    s = sim.state
end


### CONTACT EVENT ###

function get_rate(e::Val{:Contact}, sim::Simulator)
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

function do_event!(e::Val{:Contact}, sim::Simulator, t::Float64)
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
        @debug "Infecting!" t
        # Reduce bacterial population
        @assert s.bstrains.abundance[iB] > 0
        s.bstrains.abundance[iB] -= 1
        s.bstrains.total_abundance -= 1

        # Increase viral population
        s.vstrains.abundance[jV] += params.beta_burst_size
        s.vstrains.total_abundance += params.beta_burst_size
    elseif should_acquire_spacer
        @info "Acquiring spacer!" t
        @assert !should_infect

        # Create new bacterial strain with modified spacers
        s.bstrains.abundance[iB] -= 1

        new_spacers = s.bstrains.spacers[iB][2:length(s.bstrains.spacers[iB])]
        push!(new_spacers, rand(rng, s.vstrains.pspacers[jV]))
        @info "new_spacers" t new_spacers

        mutated_strain_index = findfirst(x -> x == new_spacers, s.bstrains.spacers)
        if mutated_strain_index === nothing
            @info "creating new strain"
            id = s.bstrains.next_id
            s.bstrains.next_id += 1
            push!(s.bstrains.ids, id)
            push!(s.bstrains.abundance, 1)
            push!(s.bstrains.spacers, new_spacers)
        else
            @info "using old strain" mutated_strain_index
            s.bstrains.abundance[mutated_strain_index] += 1
        end
    end

    # Update affected rates
#     update_rate!(Val(:BacterialGrowth), sim)
#     update_rate!(Val(:BacterialDeath), sim)
#     update_rate!(Val(:Contact), sim)
end

function is_immune(spacers::Vector{UInt64}, pspacers::Vector{UInt64})
    length(intersect(spacers, pspacers)) > 0
end


### VALIDATION TO ENSURE CODE CORRECTNESS IN THE FACE OF OPTIMIZATIONS ###

function validate(s::State)
    validate(s.bstrains)
    validate(s.vstrains)
#     if length(s.bstrains) > 0
#         @assert next_bstrain_id >= maximum(s.bstrains.id)
#     end

#     if length(s.vstrains) > 0
#         @assert next_vstrain_id >= maximum(s.vstrains.id)
#     end
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

# function validate(vs::VStrain)
# end

end # module StochasticCrispr


