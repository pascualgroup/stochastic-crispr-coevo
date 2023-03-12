struct params
    adsorption_rate::Float64
    viral_mutation_rate::Float64
    n_protospacers::Float64
    spacer_acquisition_prob::Float64
    viral_burst_size::Float64
    viral_decay_rate::Float64
    function params(dbTempSim)
        (comboID,) = execute(dbTempSim, "SELECT combo_id FROM runs WHERE run_id = $(run_id)")
        comboID = comboID.combo_id
        (param,) = execute(dbTempSim, "SELECT adsorption_rate, viral_mutation_rate,
            n_protospacers, spacer_acquisition_prob, viral_burst_size, viral_decay_rate
            FROM param_combos WHERE combo_id = $(comboID)")
            new(
                param.adsorption_rate,
                param.viral_mutation_rate,
                param.n_protospacers,
                param.spacer_acquisition_prob,
                param.viral_burst_size,
                param.viral_decay_rate
                )
    end
end

mutable struct modHierarchy
    dbTri::DB
    dbMatch::DB
    dbSim::DB
    dbOutput::DB
    matchID::Int64
    matchtypes::Dict{Int64,Vector{Int64}}
    matches::Dict{Int64,Vector{Int64}}
    escapes::Dict{Int64,Vector{Int64}}
    lineages::Dict{Int64,Vector{Int64}}
    ispacers::Dict{Int64,Vector{Int64}}
    ibstrains::Dict{Int64,Vector{Int64}}
    function modHierarchy(dbTempTri::DB,dbTempMatch::DB,dbTempSim::DB,dbOutput::DB)
        new(
            dbTempTri,
            dbTempMatch,
            dbTempSim,
            dbOutput,
            Int64(0),
            Dict{Int64,Vector{Int64}}(),
            Dict{Int64,Vector{Int64}}(),
            Dict{Int64,Vector{Int64}}(),
            Dict{Int64,Vector{Int64}}(),
            Dict{Int64,Vector{Int64}}(),
            Dict{Int64,Vector{Int64}}()
            )
    end
end

mutable struct laststate
    q::Dict{Int64,Float64}
    spacers::Vector{Int64}
    function laststate()
        new(
            Dict{Int64,Float64}(),
            Vector{Int64}()
            )
    end
end

mutable struct intState
    dbTri::DB
    dbMatch::DB
    dbSim::DB
    dbOutput::DB
    time::Float64
    matches::Vector{Int64}
    escapes::Vector{Int64}
    sBiomass::Dict{Int64,Int64}
    iBiomass::Dict{Int64,Int64}
    function intState(matchStructure::modHierarchy,time::Float64)
        new(
            matchStructure.dbTri,
            matchStructure.dbMatch,
            matchStructure.dbSim,
            matchStructure.dbOutput,
            time,
            Vector{Int64}(),
            Vector{Int64}([0]),
            Dict{Int64,Int64}(),
            Dict{Int64,Int64}()
            )
    end
end

mutable struct intPcomponents
    dbSim::DB
    time::Float64
    sigdig::Int64
    birth::Dict{Int64,Float64}
    death::Dict{Int64,Float64}
    mut::Dict{Int64,Float64}
    lysis::Dict{Int64,Float64}
    prootmin::Dict{Int64,Float64}
    proots::Dict{Int64,Vector{Float64}}
    plambert::Dict{Int64,Float64}
    pactual::Dict{Int64,Float64}
    pintegrated::Dict{Int64,Float64}
    function intPcomponents(current::intState)
        new(
            current.dbSim,
            current.time,
            Int64(7),
            Dict{Int64,Float64}(),
            Dict{Int64,Float64}(),
            Dict{Int64,Float64}(),
            Dict{Int64,Float64}(),
            Dict{Int64,Float64}(),
            Dict{Int64,Vector{Float64}}(),
            Dict{Int64,Float64}(),
            Dict{Int64,Float64}(),
            Dict{Int64,Float64}()
            )
    end
end
