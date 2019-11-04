module StochasticCrispr

using Util
using Random
using Logging
using StatsBase

export Parameters, validate, State, Simulator, do_next_event!

mutable struct Parameters
    "CRIPSR failure probability"
    p_crispr_failure_prob::Float64

    "New spacer acquisition probability"
    q_spacer_acquisition_prob::Float64

    "Growth rate at 0 (1/h)"
    r_growth_rate::Float64

    "Carrying capacity (1/mL)"
    K_carrying_capacity::Float64

    "Burst size"
    beta_burst_size::Float64

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
    @assert 0.0 <= p.beta_burst_size <= 1000
    @assert 0.0 <= p.phi_adsorption_rate <= 1e-6
    @assert 0.0 <= p.m_viral_decay_rate <= 1.0
    @assert 0.0 <= p.mu_mutation_rate <= 1e-5
    @assert 0.05 <= p.rho_c_density_cutoff <= 1.0
end

struct Spacer
end

struct PSpacer
end

mutable struct BStrains
    next_id::UInt64
    abundance::Vector{Int64}
    total_abundance::Int64

    function BStrains(n)
        new(n + 1, ones(n), n)
    end
end

# mutable struct BStrain
#     id::UInt64
#     abundance::Int64
#     spacers::Vector{Spacer}
#
#     BStrain(id, abundance) = new(id, abundance, [])
# end

# mutable struct VStrains
#     next_id::UInt64
#     bstrains::Vector{VStrain}
#
#     function BStrains()
#         new(1, [], 0)
#     end
# end

# mutable struct VStrain
#     id::UInt64
#     abundance::Int64
#     pspacers::Vector{PSpacer}
#
#     VStrain(id, abundance) = new(id, abundance, [])
# end

mutable struct State
    bstrains::BStrains
#     vstrains::VStrains

    State() = new(BStrains(2000))
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
    Val(:ViralDeath),
    Val(:ViralMutation),
    Val(:Infection),
]

mutable struct Simulator
    parameters::Parameters
    t::Float64
    state::State
    rng::MersenneTwister

    top_level_event_rates::Vector{Float64}
end

function Simulator(params::Parameters, t_init::Float64, state_init::State)
    sim = Simulator(
        params, t_init, state_init, MersenneTwister(), zeros(length(TOP_LEVEL_EVENTS))
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

        @debug "begin do_event()" event=event t=t_next
        do_event!(event, sim, t_next)
        @debug "end do_event()"

        sim.t = t_next

        @debug "bstrains.total_abundance:" sim.state.bstrains.total_abundance
#         @debug "VStrains.total_abundance:" sim.state.vstrains.total_abundance
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
    update_rate!(Val(:BacterialGrowth), sim)
    update_rate!(Val(:BacterialDeath), sim)
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
    update_rate!(Val(:BacterialGrowth), sim)
    update_rate!(Val(:BacterialDeath), sim)
end


### VIRAL DEATH EVENT ###

function get_rate(e::Val{:ViralDeath}, sim::Simulator)
    0.0
end

function do_event!(e::Val{:ViralDeath}, sim::Simulator, t::Float64)
    p = sim.parameters
    s = sim.state
end


### VIRAL MUTATION EVENT ###

function get_rate(e::Val{:ViralMutation}, sim::Simulator)
    0.0
end

function do_event!(e::Val{:ViralMutation}, sim::Simulator, t::Float64)
    p = sim.parameters
    s = sim.state
end


### INFECTION EVENT ###

function get_rate(e::Val{:Infection}, sim::Simulator)
    0.0
end

function do_event!(e::Val{:Infection}, sim::Simulator, t::Float64)
    p = sim.parameters
    s = sim.state
end

### VALIDATION TO ENSURE CODE CORRECTNESS IN THE FACE OF OPTIMIZATIONS ###

function validate(s::State)
    validate(s.bstrains)
#     validate(s.vstrains)
#     if length(s.bstrains) > 0
#         @assert next_bstrain_id >= maximum(s.bstrains.id)
#     end

#     if length(s.vstrains) > 0
#         @assert next_vstrain_id >= maximum(s.vstrains.id)
#     end
end

function validate(bss::BStrains)
#     if length(bss.bstrains) > 0
#         @assert bss.next_id > maximum(id.(bss.bstrains))
#     end
    @assert bss.total_abundance == sum(bss.abundance)
end

# function validate(bs::BStrain)
# end

# function validate(bs::VStrains)
#     if length(vs.bstrains) > 0
#         @assert vs.next_id > maximum(id.(vs.vstrains))
#     end
#
#     for bs in bss
#         validate(bs)
#     end
# end

# function validate(vs::VStrain)
# end

end # module StochasticCrispr


