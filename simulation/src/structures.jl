
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

    "Number of initial bacterial strains"
    n_bstrains::Union{UInt64, Nothing}

    "Number of initial hosts per bacterial strain"
    n_hosts_per_bstrain::Union{UInt64, Nothing}

    "Number of initial virus strains"
    n_vstrains::Union{UInt64, Nothing}

    "Number of initial particles per bacterial strain"
    n_particles_per_vstrain::Union{UInt64, Nothing}

    "Number of initial protospacers per virus strain"
    n_protospacers::Union{UInt64, Nothing}

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


function init_params(d_symb::Dict{Symbol,Any})
    Params(;
        t_final = d_symb[:t_final][1],

        t_output = d_symb[:t_output][1],

        rng_seed = d_symb[:rng_seed][1],

        enable_output = d_symb[:enable_output][1],

        n_bstrains = d_symb[:n_bstrains][1],

        n_hosts_per_bstrain = d_symb[:n_hosts_per_bstrain][1],

        n_vstrains = d_symb[:n_vstrains][1],

        n_particles_per_vstrain = d_symb[:n_particles_per_vstrain][1],

        n_protospacers = d_symb[:n_protospacers][1],

        n_spacers_max = d_symb[:n_spacers_max][1],

        crispr_failure_prob = d_symb[:crispr_failure_prob][1],

        spacer_acquisition_prob = d_symb[:spacer_acquisition_prob][1],

        microbe_growth_rate = d_symb[:spacer_acquisition_prob][1],

        microbe_carrying_capacity = d_symb[:microbe_carrying_capacity][1],

        viral_burst_size = d_symb[:viral_burst_size][1],

        adsorption_rate = d_symb[:adsorption_rate][1],

        viral_decay_rate = d_symb[:viral_decay_rate][1],

        viral_mutation_rate = d_symb[:viral_mutation_rate][1],

        microbe_death_rate = d_symb[:microbe_death_rate][1],

        microbe_immigration_rate = d_symb[:microbe_immigration_rate][1],
    )
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
        # println("$([p.n_bstrains, p.n_hosts_per_bstrain,
        # p.n_vstrains, p.n_particles_per_vstrain, p.n_protospacers])")
        # println("$(map(x->typeof(x),[p.n_bstrains, p.n_hosts_per_bstrain,
        # p.n_vstrains, p.n_particles_per_vstrain, p.n_protospacers]))")
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
