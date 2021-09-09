println("(Annoying Julia compilation delay...)")

include("simulation/src/setup.jl")
include("simulation/src/structures.jl")
include("simulation/src/util.jl")

SCRIPT_PATH = abspath(dirname(PROGRAM_FILE))

function make_sweep_params_script()
    numParams = []
    t_final = Float64(2000.0)
    push!(numParams, length(t_final))
    t_output = Float64(1.0)
    push!(numParams, length(t_output))
    rng_seed = UInt64(1)
    push!(numParams, length(rng_seed))

    enable_output = true

    n_bstrains = UInt64(1)
    push!(numParams, length(n_bstrains))
    n_hosts_per_bstrain = UInt64(100)
    push!(numParams, length(n_hosts_per_bstrain))
    n_vstrains = UInt64(1)
    push!(numParams, length(n_vstrains))
    n_particles_per_vstrain = UInt64(100)
    push!(numParams, length(n_particles_per_vstrain))
    n_protospacers = [UInt64(5), UInt64(10), UInt64(15)]
    push!(numParams, length(n_protospacers))
    n_spacers_max = [UInt64(5), UInt64(10), UInt64(15)]
    push!(numParams, length(n_spacers_max))
    crispr_failure_prob = [Float64(5e-06), Float64(1e-05), Float64(1.5e-05), Float64(2e-05)]
    push!(numParams, length(crispr_failure_prob))
    spacer_acquisition_prob = [Float64(5e-06), Float64(1e-05), Float64(1.5e-05), Float64(2e-05)]
    push!(numParams, length(spacer_acquisition_prob))
    microbe_growth_rate = Float64(1.0)
    push!(numParams, length(microbe_growth_rate))
    microbe_carrying_capacity = Float64(316228.0) #~10^(5.5)
    push!(numParams, length(microbe_carrying_capacity))
    viral_burst_size = UInt64(50)
    push!(numParams, length(viral_burst_size))
    adsorption_rate = Float64(1e-07)
    push!(numParams, length(adsorption_rate))
    viral_decay_rate = Float64(0.1)
    push!(numParams, length(viral_decay_rate))
    viral_mutation_rate = [Float64(5e-07), Float64(7.5e-07), Float64(1e-06), Float64(1.75e-06), Float64(2e-06)]
    push!(numParams, length(viral_mutation_rate))
    microbe_death_rate = Float64(0)
    push!(numParams, length(microbe_death_rate))
    microbe_immigration_rate = Float64(0)
    push!(numParams, length(microbe_immigration_rate))


    params = ParamSweep(;
        t_final = t_final,
        t_output = t_output,
        rng_seed = rng_seed,
        enable_output = enable_output,
        n_bstrains = n_bstrains,
        n_hosts_per_bstrain = n_hosts_per_bstrain,
        n_vstrains = n_vstrains,
        n_particles_per_vstrain = n_particles_per_vstrain,
        n_protospacers = n_protospacers,
        n_spacers_max = n_spacers_max,
        crispr_failure_prob = crispr_failure_prob,
        spacer_acquisition_prob = spacer_acquisition_prob,
        microbe_growth_rate = microbe_growth_rate,
        microbe_carrying_capacity = microbe_carrying_capacity,
        viral_burst_size = viral_burst_size,
        adsorption_rate = adsorption_rate,
        viral_decay_rate = viral_decay_rate,
        viral_mutation_rate = viral_mutation_rate,
        microbe_death_rate = microbe_death_rate,
        microbe_immigration_rate = microbe_immigration_rate,
    )
    validate(params)
    params_json = pretty_json(params)
    open(joinpath(SCRIPT_PATH, "sweep-parameters.json"), "w") do f
        println(f, params_json)
    end
    println("This sweep has a maximum of $(prod(numParams)) parameter combinations.
    Multiply this value with the number of replicates intended to evaluate numbers of cores necessary.")
end

function pretty_json(params)
    d = Dict(fn => getfield(params, fn) for fn in fieldnames(typeof(params)))
    io = IOBuffer()
    JSON.print(io, d, 2)
    String(take!(io))
end

@with_kw mutable struct ParamSweep
    "Simulation end time"
    t_final::Union{Float64, Array{Float64,1}, Nothing}

    "Time between output events"
    t_output::Union{Float64, Array{Float64,1}, Nothing}

    "Seed for random number generator"
    rng_seed::Union{UInt64, Array{UInt64,1}, Nothing}

    "Enable output?"
    enable_output::Union{Bool, Nothing}

    "Number of initial bacterial strains"
    n_bstrains::Union{UInt64, Array{UInt64,1}, Nothing}

    "Number of initial hosts per bacterial strain"
    n_hosts_per_bstrain::Union{UInt64, Array{UInt64,1}, Nothing}

    # "Number of initial spacers per bacterial strain"
    # n_spacers::UInt64

    "Number of initial virus strains"
    n_vstrains::Union{UInt64, Array{UInt64,1}, Nothing}

    "Number of initial particles per bacterial strain"
    n_particles_per_vstrain::Union{UInt64, Array{UInt64,1}, Nothing}

    "Number of initial protospacers per virus strain"
    n_protospacers::Union{UInt64, Array{UInt64,1}, Nothing}

    #InitializationParameters() = new() # WHAT IS THIS????

    "Maximum number of spacers in a bacterial strain"
    n_spacers_max::Union{UInt64, Array{UInt64,1}, Nothing}

    "CRIPSR failure probability [p]"
    crispr_failure_prob::Union{Float64, Array{Float64,1}, Nothing}

    "New spacer acquisition probability [q]"
    spacer_acquisition_prob::Union{Float64, Array{Float64,1}, Nothing}

    "Growth rate at 0 (1/h) = [r]"
    microbe_growth_rate::Union{Float64, Array{Float64,1}, Nothing}

    "Carrying capacity (1/mL) = [K]"
    microbe_carrying_capacity::Union{Float64, Array{Float64,1}, Nothing}

    "Burst size [beta]"
    viral_burst_size::Union{UInt64, Array{UInt64,1}, Nothing}

    "Adsorption rate [phi]"
    adsorption_rate::Union{Float64, Array{Float64,1}, Nothing}

    "Viral decay rate [m]"
    viral_decay_rate::Union{Float64, Array{Float64,1}, Nothing}

    "Mutation rate [mu]"
    viral_mutation_rate::Union{Float64, Array{Float64,1}, Nothing}

    "Constant death rate (not in Childs model) [d]"
    microbe_death_rate::Union{Float64, Array{Float64,1}, Nothing}

    "Constant immigration rate (not in Childs model) [eta]"
    microbe_immigration_rate::Union{Float64, Array{Float64,1}, Nothing}
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
    @assert p.microbe_growth_rate !== nothing
    @assert p.microbe_carrying_capacity !== nothing
    @assert p.viral_burst_size !== nothing
    @assert p.adsorption_rate !== nothing
    @assert p.viral_decay_rate !== nothing
    @assert p.viral_mutation_rate !== nothing
    @assert p.microbe_death_rate !== nothing
    @assert p.microbe_immigration_rate !== nothing
end


make_sweep_params_script()
