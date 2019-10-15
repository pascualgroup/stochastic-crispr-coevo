module StochasticCrispr

using Random
using StatsBase

export Parameters, validate, State, Simulator, do_next_event!

mutable struct Parameters
    "Viral mutation rate per virus"
    beta::Float64

    Parameters() = new()
end

function validate(p::Parameters)
    @assert p.beta >= 0 && p.beta <= 100
end

struct Spacer
end

struct PSpacer
end

mutable struct BStrain
    id::UInt64
    abundance::Float64
    spacers::Vector{Spacer}
end

mutable struct VStrain
    id::UInt64
    abundance::Float64
    pspacers::Vector{PSpacer}
end

mutable struct State
    next_bstrain_id::UInt64
    bstrains::Vector{BStrain}

    next_vstrain_id::UInt64
    vstrains::Vector{VStrain}

    function State()
        new(1, [], 1, [])
    end
end

function id(x)
    x.id
end

function validate(s::State)
    if length(s.bstrains) > 0
        @assert next_bstrain_id >= maximum(id.(s.bstrains))
    end

    if length(s.vstrains) > 0
        @assert next_vstrain_id >= maximum(id.(s.vstrains))
    end
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
    Val(:ViralExtinction),
    Val(:BacterialExtinction),
    Val(:ViralMutation),
    Val(:Infection),
]

const TOP_LEVEL_EVENT_INDICES = map(i -> TOP_LEVEL_EVENTS[i], 1:length(TOP_LEVEL_EVENTS))

mutable struct Simulator
    parameters::Parameters
    t::Float64
    state::State
    rng::MersenneTwister

    top_level_event_rates::Vector{Float64}
end

function Simulator(params::Parameters, t_init::Float64, state_init::State)
    Simulator(params, t_init, state_init, MersenneTwister(), ones(length(TOP_LEVEL_EVENTS)))
end

function do_next_event!(sim::Simulator, t_max::Float64)
    p = sim.parameters
    s = sim.state

    @assert length(sim.top_level_event_rates) == length(TOP_LEVEL_EVENTS)
    R = sum(sim.top_level_event_rates)

    t_next = sim.t + randexp(sim.rng) / R

    if t_next > t_max
        sim.t = t_max
    else
        event = sample(sim.rng, TOP_LEVEL_EVENTS, Weights(sim.top_level_event_rates, R))
        #event = TOP_LEVEL_EVENTS[1]
        #println(t_next, event)
        do_event!(event, sim, t_next)
        sim.t = t_next
    end
end

function do_event!(e::Val{:BacterialGrowth}, sim::Simulator, t::Float64)
    p = sim.parameters
    s = sim.state
end

function do_event!(e::Val{:ViralExtinction}, sim::Simulator, t::Float64)
    p = sim.parameters
    s = sim.state
end

function do_event!(e::Val{:BacterialExtinction}, sim::Simulator, t::Float64)
    p = sim.parameters
    s = sim.state
end

function do_event!(e::Val{:ViralMutation}, sim::Simulator, t::Float64)
    p = sim.parameters
    s = sim.state
end

function do_event!(e::Val{:Infection}, sim::Simulator, t::Float64)
    p = sim.parameters
    s = sim.state
end

end
