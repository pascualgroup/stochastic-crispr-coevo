#!/usr/bin/env julia

println("(Julia compilation delay...)")

using Random
using Distributions
using SQLite
using DataFrames
using SQLite.DBInterface: execute
using StatsBase
using Parameters
using JSON

SAMPLE_PATH = ARGS[1]
sampleTime = parse(Float64, ARGS[2])
vstrainID = parse(Int64, ARGS[3])
vAbund0 = 1

# Ensure that UInt64 and UInt32 gets printed in decimal, not hex
Base.show(io::IO, x::UInt64) = Base.print(io, x)
Base.show(io::IO, y::UInt32) = Base.print(io, y)
##
rngDiv = MersenneTwister(9876543456)

## Define Paths ##
SCRIPT_PATH = abspath(dirname(PROGRAM_FILE))
include(joinpath(SCRIPT_PATH, "simulation", "src", "fitness.jl"))
dbOutputPath = joinpath(SCRIPT_PATH, "initial-conditions.sqlite") # cluster
##
if isfile(dbOutputPath)
    # error("initial-conditions.sqlite already exists; delete first")
    rm(dbOutputPath, force=true)
end
dbOutput = SQLite.DB(dbOutputPath)
dbSample = SQLite.DB(SAMPLE_PATH)

execute(dbOutput, "CREATE TABLE babundance (bstrain_id INTEGER, abundance INTEGER)")
execute(dbOutput, "CREATE TABLE bspacers (bstrain_id INTEGER, spacer_id INTEGER)")
execute(dbOutput, "CREATE TABLE vabundance (vstrain_id INTEGER, abundance INTEGER)")
execute(dbOutput, "CREATE TABLE vpspacers (vstrain_id INTEGER, spacer_id INTEGER)")
##

function make_sweep_params_script(SCRIPT_PATH)
    t_final = Float64(2000.0)
    t_output = Float64(1.0)
    rng_seed = UInt64(1) # this will be generated randomly upon calling generate-sweep.jl
    enable_output = true
    n_hosts_per_bstrain = UInt64(100) # pick any value
    n_particles_per_vstrain = UInt64(100) # pick any value
    n_protospacers = UInt64(15)
    n_spacers_max = UInt64(10)

    execute(dbOutput, "BEGIN TRANSACTION")
    execute(dbOutput, "ATTACH DATABASE '$(SAMPLE_PATH)' as dbSim")
    execute(dbOutput, "INSERT INTO babundance (bstrain_id, abundance) SELECT bstrain_id, abundance FROM dbSim.babundance WHERE t = $(sampleTime);")
    bstrainIDs = [bID for (bID,) in execute(dbSample, "SELECT bstrain_id FROM babundance WHERE t = $(sampleTime)")]
    n_bstrains = UInt64(length(bstrainIDs)) # pick any value
    n_vstrains = UInt64(1) # pick any value
    execute(
        dbOutput,
        "INSERT INTO bspacers (bstrain_id, spacer_id) SELECT bstrain_id, spacer_id 
          FROM dbSim.bspacers WHERE bstrain_id in ($(join(bstrainIDs,", ")));"
    )
    execute(dbOutput, "COMMIT")
    execute(dbOutput, "INSERT INTO vabundance VALUES ($(vstrainID),$(vAbund0))")
    pspacers = [sID for (sID,) in execute(dbSample, "SELECT spacer_id FROM vpspacers WHERE vstrain_id = $(vstrainID)")]
    for spacerID in pspacers
        execute(dbOutput, "INSERT INTO vpspacers VALUES (?,?)", (Int64(vstrainID), Int64(spacerID)))
    end
    ## Thease are parameters for actual dynamics
    crispr_failure_prob = Float64(0)
    spacer_acquisition_prob = Float64(1.9e-05)
    microbe_mutation_prob = Float64(0)
    microbe_carrying_capacity = Float64(sum([abund for (abund,) in execute(dbOutput,"SELECT abundance FROM babundance")]))
    viral_burst_size = UInt64(50)
    adsorption_rate = Float64(1e-07)
    viral_decay_rate = Float64(0.1)
    viral_mutation_rate = Float64(1.04e-06)
    microbe_growth_rate = Float64(1.0)
    microbe_death_rate = Float64(0.1)
    microbe_immigration_rate = Float64(0)
    ##
    numParams = []
    push!(numParams, length(n_bstrains))
    push!(numParams, length(n_vstrains))
    push!(numParams, length(microbe_immigration_rate))
    push!(numParams, length(microbe_growth_rate))
    push!(numParams, length(microbe_death_rate))
    push!(numParams, length(viral_mutation_rate))
    push!(numParams, length(viral_decay_rate))
    push!(numParams, length(adsorption_rate))
    push!(numParams, length(viral_burst_size))
    push!(numParams, length(microbe_carrying_capacity))
    push!(numParams, length(microbe_mutation_prob))
    push!(numParams, length(spacer_acquisition_prob))
    push!(numParams, length(crispr_failure_prob))
    push!(numParams, length(n_spacers_max))
    push!(numParams, length(n_protospacers))
    push!(numParams, length(n_particles_per_vstrain))
    push!(numParams, length(n_hosts_per_bstrain))
    push!(numParams, length(rng_seed))
    push!(numParams, length(t_output))
    push!(numParams, length(t_final))

    params = ParamSweep(;
        t_final=t_final,
        t_output=t_output,
        rng_seed=rng_seed,
        enable_output=enable_output,
        n_bstrains=n_bstrains,
        n_hosts_per_bstrain=n_hosts_per_bstrain,
        n_vstrains=n_vstrains,
        n_particles_per_vstrain=n_particles_per_vstrain,
        n_protospacers=n_protospacers,
        n_spacers_max=n_spacers_max,
        crispr_failure_prob=crispr_failure_prob,
        spacer_acquisition_prob=spacer_acquisition_prob,
        microbe_mutation_prob=microbe_mutation_prob,
        microbe_carrying_capacity=microbe_carrying_capacity,
        viral_burst_size=viral_burst_size,
        adsorption_rate=adsorption_rate,
        viral_decay_rate=viral_decay_rate,
        viral_mutation_rate=viral_mutation_rate,
        microbe_growth_rate=microbe_growth_rate,
        microbe_death_rate=microbe_death_rate,
        microbe_immigration_rate=microbe_immigration_rate
    )
    validate(params)
    params_json = pretty_json(params)
    open(joinpath(SCRIPT_PATH, "sweep-parameters.json"), "w") do f
        println(f, params_json)
    end
    @info "This sweep has $(prod(numParams)) parameter combinations.
    Multiply this value with the number of replicates intended to evaluate numbers of cores necessary."
end

function pretty_json(params)
    d = Dict(fn => getfield(params, fn) for fn in fieldnames(typeof(params)))
    io = IOBuffer()
    JSON.print(io, d, 2)
    String(take!(io))
end

@with_kw mutable struct ParamSweep
    "Simulation end time"
    t_final::Union{Float64,Array{Float64,1},Nothing}

    "Time between output events"
    t_output::Union{Float64,Array{Float64,1},Nothing}

    "Seed for random number generator"
    rng_seed::Union{UInt64,Array{UInt64,1},Nothing}

    "Enable output?"
    enable_output::Union{Bool,Nothing}

    "Number of initial bacterial strains"
    n_bstrains::Union{UInt64,Array{UInt64,1},Nothing}

    "Number of initial hosts per bacterial strain"
    n_hosts_per_bstrain::Union{UInt64,Array{UInt64,1},Nothing}

    "Number of initial virus strains"
    n_vstrains::Union{UInt64,Array{UInt64,1},Nothing}

    "Number of initial particles per bacterial strain"
    n_particles_per_vstrain::Union{UInt64,Array{UInt64,1},Nothing}

    "Number of initial protospacers per virus strain"
    n_protospacers::Union{UInt64,Array{UInt64,1},Nothing}

    "Maximum number of spacers in a bacterial strain"
    n_spacers_max::Union{UInt64,Array{UInt64,1},Nothing}

    "CRIPSR failure probability [p]"
    crispr_failure_prob::Union{Float64,Array{Float64,1},Nothing}

    "New spacer acquisition probability [q]"
    spacer_acquisition_prob::Union{Float64,Array{Float64,1},Nothing}

    "Microbial mutation probability [rho]"
    microbe_mutation_prob::Union{Float64,Array{Float64,1},Nothing}

    "Carrying capacity (1/mL) = [K]"
    microbe_carrying_capacity::Union{Float64,Array{Float64,1},Nothing}

    "Burst size [beta]"
    viral_burst_size::Union{UInt64,Array{UInt64,1},Nothing}

    "Adsorption rate [phi]"
    adsorption_rate::Union{Float64,Array{Float64,1},Nothing}

    "Viral decay rate [m]"
    viral_decay_rate::Union{Float64,Array{Float64,1},Nothing}

    "Mutation rate [mu]"
    viral_mutation_rate::Union{Float64,Array{Float64,1},Nothing}

    "Growth rate at 0 (1/h) = [r]"
    microbe_growth_rate::Union{Float64,Array{Float64,1},Nothing}

    "Constant death rate (not in Childs model) [d]"
    microbe_death_rate::Union{Float64,Array{Float64,1},Nothing}

    "Constant immigration rate (not in Childs model) [eta]"
    microbe_immigration_rate::Union{Float64,Array{Float64,1},Nothing}
end


function validate(p::ParamSweep)
    @assert p.t_final !== nothing
    @assert p.t_output !== nothing

    @assert p.n_bstrains !== nothing
    @assert p.n_hosts_per_bstrain !== nothing
    @assert p.n_vstrains !== nothing
    @assert p.n_particles_per_vstrain !== nothing
    @assert p.n_protospacers !== nothing
    @assert p.n_spacers_max !== nothing

    @assert p.crispr_failure_prob !== nothing
    @assert p.spacer_acquisition_prob !== nothing
    @assert p.microbe_mutation_prob !== nothing
    @assert p.microbe_carrying_capacity !== nothing
    @assert p.viral_burst_size !== nothing
    @assert p.adsorption_rate !== nothing
    @assert p.viral_decay_rate !== nothing
    @assert p.viral_mutation_rate !== nothing
    @assert p.microbe_growth_rate !== nothing
    @assert p.microbe_death_rate !== nothing
    @assert p.microbe_immigration_rate !== nothing
end


make_sweep_params_script(SCRIPT_PATH)
