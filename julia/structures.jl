using Random
using Distributions
using StatsBase
#using DelimitedFiles
using Dates
using Parameters
#using JSON2

using SQLite: DB, Stmt, bind!
using SQLite.DBInterface: execute

### PARAMETERS ###

@with_kw mutable struct Params
    "Simulation end time"
    t_final::Union{Float64, Nothing}

    "Time between output events"
    t_output::Union{Float64, Nothing}

    "Seed for random number generator"
    rng_seed::Union{UInt64, Nothing}

    "Enable output?"
    enable_output::Union{Bool, Nothing}

    # REMOVE COMMENTED BLOCK BELOW
    #function RunParameters()
        #p = new()
        #p.rng_seed = nothing
        #p.enable_output = true
        #p
    #end

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

    #InitializationParameters() = new() # WHAT IS THIS????

    "Maximum number of spacers in a bacterial strain"
    n_spacers_max::Union{UInt64, Nothing}

    "CRIPSR failure probability [p]"
    crispr_failure_prob::Union{Float64, Nothing}

    "New spacer acquisition probability [q]"
    spacer_acquisition_prob::Union{Float64, Nothing}

    "Growth rate at 0 (1/h) = [r]"
    microbe_growth_rate::Union{Float64, Nothing}

    "Carrying capacity (1/mL) = [K]"
    microbe_carrying_capacity::Union{Float64, Nothing}

    "Burst size [beta]"
    viral_burst_size::Union{UInt64, Nothing}

    "Adsorption rate [phi]"
    adsorption_rate::Union{Float64, Nothing}

    "Viral decay rate [m]"
    viral_decay_rate::Union{Float64, Nothing}

    "Mutation rate [mu]"
    viral_mutation_rate::Union{Float64, Nothing}

    "Constant death rate (not in Childs model) [d]"
    microbe_death_rate::Union{Float64, Nothing}

    "Constant immigration rate (not in Childs model) [eta]"
    microbe_immigration_rate::Union{Float64, Nothing}

    #function Params() #this function is important for proper function of JSON2
        #new()
        ##p = new()
        ##p.rng_seed = nothing # when feeding json script, with pre-defined entry,
        ###this function will not successfully change values
        ##p.enable_output = true # applies to this as well...
        ##p
    #end
end






### SIMULATION STATE ###





mutable struct Strains
    next_id::UInt64
    ids::Vector{UInt64}

    abundance::Vector{UInt64}
    total_abundance::UInt64

    spacers::Vector{Vector{UInt64}}

    #strain_file::IOStream
    #spacers_file::IOStream
    #abundance_file::IOStream
end




function make_bstrains(n_strains, n_hosts_per_strain)
    Strains(
        n_strains + 1,
        1:n_strains,
        repeat([n_hosts_per_strain], n_strains),
        n_strains * n_hosts_per_strain,
        repeat([[]], n_strains)#,
        #open_csv("bstrains", "t_creation", "bstrain_id", "parent_bstrain_id", "infecting_vstrain_id"),
        #open_csv("bspacers", "bstrain_id", "spacer_id"),
        #open_csv("babundance", "t", "bstrain_id", "abundance")
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
        pspacers#,
        #open_csv("vstrains", "t_creation", "vstrain_id", "parent_vstrain_id", "infected_bstrain_id"),
        #open_csv("vpspacers", "vstrain_id", "spacer_id"),
        #open_csv("vabundance", "t", "vstrain_id", "abundance")
    )
end







#########################################################
##################THIS SAVES NEW STRAINS#################

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





function State(p::Params)
    State(
        p.n_bstrains, p.n_hosts_per_bstrain,
        p.n_vstrains, p.n_particles_per_vstrain, p.n_protospacers
    )
end





#########################################################






mutable struct Simulation
    params::Params
    t::Float64
    state::State
    rng::MersenneTwister

    event_rates::Vector{Float64}
    event_counts::Vector{UInt64}

    db::DB

    #meta_file::IOStream
    #summary_file::IOStream

    function Simulation(p::Params)
        #meta_file = open_csv("meta", "key", "value")

        db = initialize_database()

        # Use random seed if provided, or generate one
        rng_seed = p.rng_seed === nothing ? UInt64(rand(RandomDevice(), UInt32)) : p.rng_seed
        p.rng_seed = rng_seed

        #write_csv(meta_file, "rng_seed", rng_seed)

        execute(db,
            "INSERT INTO meta VALUES (?,?)", ["rng_seed", Int64(rng_seed)]
        )

        # Save parameters as loaded
        #write_json_to_file(p, "parameters_out.json")

        #write_csv(meta_file, "start_time", start_time)

        execute(db,
            "INSERT INTO meta VALUES (?,?)",
            ["start_time", Dates.format(start_time, "yyyy-mm-ddTHH:MM:SS")]
        )

        # Initialize & validate model state
        state = State(p)
        validate(state)

        sim = new(
            p, 0.0, state, MersenneTwister(rng_seed),
            zeros(length(EVENTS)), zeros(length(EVENTS)),
            #meta_file,
            #open_csv("summary", "t", "microbial_abundance", "viral_abundance")
            db
        )
        update_rates!(sim)
        sim
    end
end
