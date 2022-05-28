#!/usr/bin/env julia

println("(Julia compilation delay...)")

using SQLite
using DataFrames
using SQLite.DBInterface: execute
using SQLite: DB
using Combinatorics: combinations

run_id = ARGS[1]

## Define Paths ##
SCRIPT_PATH = abspath(dirname(PROGRAM_FILE))
#
dbSimPath = joinpath(SCRIPT_PATH,"..","..","..","simulation","sweep_db_gathered.sqlite") # cluster
dbOutputPath = joinpath("probability-of-emergence_output.sqlite") # cluster

# dbSimPath = joinpath("/Volumes/Yadgah","crispr-sweep-7-2-2022/isolates/runID3297-c66-r47/runID3297-c66-r47.sqlite") # local
# dbMatchPath = joinpath("/Volumes/Yadgah","crispr-sweep-7-2-2022/isolates/runID3297-c66-r47/matches_output.sqlite")
# # dbOutputPath = joinpath("/Volumes/Yadgah/crispr-sweep-7-2-2022/isolates/runID3297-c66-r47/probability-of-emergence_output.sqlite") # local

if isfile(dbOutputPath)
    error("probability-of-emergence_output.sqlite already exists; delete first")
end
##
dbSim = SQLite.DB(dbSimPath)
(cr,) = execute(dbSim,"SELECT combo_id,replicate FROM runs WHERE run_id = $(run_id)")
dbMatchPath = joinpath(SCRIPT_PATH,"..","..","isolates",
    "runID$(run_id)-c$(cr.combo_id)-r$(cr.replicate)","matches_output.sqlite")
if !isfile(dbMatchPath)
    error("matches_output.sqlite does not exist; compute first")
end
dbTempMatch = SQLite.DB(dbMatchPath)
dbOutput = SQLite.DB(dbOutputPath)

execute(dbOutput, "CREATE TABLE pExtinction (t REAL, match_id INTEGER, p REAL)")
execute(dbOutput, "CREATE TABLE pEmergenceProd (t REAL, p REAL)")
execute(dbOutput, "CREATE TABLE pEmergenceSum (t REAL, p REAL)")
execute(dbOutput, "CREATE TABLE pExtWeighted
    (t REAL, time_specific_match_id INTEGER, p REAL)")
execute(dbOutput, "CREATE TABLE Rweighted
    (t REAL, time_specific_match_id INTEGER, R REAL)")
execute(dbOutput, "CREATE TABLE Ravg (t REAL, R REAL)")
execute(dbOutput, "CREATE TABLE R0 (t REAL, R REAL)")
execute(dbOutput, "CREATE TABLE match_phenotypes
    (t REAL, time_specific_match_id INTEGER, phenotype INTEGER)")
execute(dbOutput, "CREATE TABLE match_phenotype_abundances
    (t REAL, time_specific_match_id INTEGER, abundance INTEGER, frequency INTEGER)")

dbTempSim = SQLite.DB()
execute(dbTempSim, "CREATE TABLE babundance (t REAL, bstrain_id, abundance INTEGER)")
execute(dbTempSim, "CREATE TABLE vabundance (t REAL, vstrain_id, abundance INTEGER)")
execute(dbTempSim, "CREATE TABLE bspacers (bstrain_id INTEGER, spacer_id INTEGER)")

execute(dbTempSim, "BEGIN TRANSACTION")
execute(dbTempSim,"ATTACH DATABASE '$(dbSimPath)' as dbSim")
execute(dbTempSim,"INSERT INTO babundance (t, bstrain_id, abundance)
    SELECT t, bstrain_id, abundance FROM dbSim.babundance WHERE run_id = $(run_id);")
execute(dbTempSim,"INSERT INTO vabundance (t, vstrain_id, abundance)
    SELECT t, vstrain_id, abundance FROM dbSim.vabundance WHERE run_id = $(run_id);")
execute(dbTempSim,"INSERT INTO bspacers (bstrain_id,spacer_id)
    SELECT bstrain_id,spacer_id FROM dbSim.bspacers WHERE run_id = $(run_id);")
execute(dbTempSim, "COMMIT")


execute(dbTempSim, "BEGIN TRANSACTION")
execute(dbTempSim, "CREATE INDEX babundance_index ON babundance (t,bstrain_id,abundance)")
execute(dbTempSim, "CREATE INDEX vabundance_index ON vabundance (t,vstrain_id,abundance)")
execute(dbTempSim, "CREATE INDEX bspacers_index ON bspacers (bstrain_id,spacer_id)")
execute(dbTempSim, "COMMIT")

struct params
    adsorption_rate::Float64
    viral_mutation_rate::Float64
    n_protospacers::Float64
    spacer_acquisition_prob::Float64
    viral_burst_size::Float64
    viral_decay_rate::Float64
    innovation_rate::Float64
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
                param.viral_decay_rate,
                1 - (1-param.viral_mutation_rate)^param.n_protospacers
                )
    end
end

mutable struct vPhenoInfo
    potential::Vector{Vector{Int64}}
    current::Vector{Vector{Int64}}
end
mutable struct bPhenoInfo
    potential::Dict{Vector{Int64},Vector{Int64}}
    current::Dict{Vector{Int64},Vector{Int64}}
end
mutable struct massInfo
    potential::Dict{Vector{Int64},Int64}
    current::Dict{Vector{Int64},Int64}
end

mutable struct matchHierarchy
    dbMatch::DB
    dbSim::DB
    dbOutput::DB
    time::Float64
    parameters::params
    matchtypes::Vector{Int64}
    phenotypes::vPhenoInfo
    strainclasses::Dict{Vector{Int64},Vector{Int64}}
    sBstrains::bPhenoInfo
    iBstrains::bPhenoInfo
    phenoBiomass::Dict{Vector{Int64},Int64}
    phenoFrequency::Dict{Vector{Int64},Float64}
    sBiomass::massInfo
    iBiomass::massInfo
    function matchHierarchy(dbTempMatch::DB,dbTempSim::DB,dbOutput::DB,t::Float64)
        new(
            dbTempMatch,
            dbTempSim,
            dbOutput,
            t,
            params(dbTempSim),
            Vector{Int64}(),
            vPhenoInfo(
                Vector{Vector{Int64}}(),Vector{Vector{Int64}}()
                ),
            Dict{Vector{Int64},Vector{Int64}}(),
            bPhenoInfo(
                Dict{Vector{Int64},Vector{Int64}}(),Dict{Vector{Int64},Vector{Int64}}()
                ),
            bPhenoInfo(
                Dict{Vector{Int64},Vector{Int64}}(),Dict{Vector{Int64},Vector{Int64}}()
                ),
                Dict{Vector{Int64},Int64}(),
                Dict{Vector{Int64},Float64}(),
            massInfo(
                Dict{Vector{Int64},Int64}(),Dict{Vector{Int64},Int64}()
                ),
            massInfo(
                Dict{Vector{Int64},Int64}(),Dict{Vector{Int64},Int64}()
                )
            )
    end
end

mutable struct extinction
    dbMatch::DB
    dbSim::DB
    dbOutput::DB
    R0::Float64
    birthRates::Dict{Vector{Int64},Float64}
    deathRates::Dict{Vector{Int64},Float64}
    birthMutRates::Dict{Vector{Int64},Float64}
    pointEscapeProb::Dict{Int64,Float64}
    probability::Dict{Vector{Int64},Float64}
    pWeight::Dict{Vector{Int64},Float64}
    pRoots::Dict{Vector{Int64},Vector{Float64}}
    Rweight::Dict{Vector{Int64},Float64}
    Ravg::Float64
    product::Float64
    sum::Float64

    function extinction(matchStructure::matchHierarchy)
        new(
            matchStructure.dbMatch,
            matchStructure.dbSim,
            matchStructure.dbOutput,
            0,
            Dict{Vector{Int64},Float64}(),
            Dict{Vector{Int64},Float64}(),
            Dict{Vector{Int64},Float64}(),
            Dict{Vector{Int64},Float64}(),
            Dict{Vector{Int64},Float64}(),
            Dict{Vector{Int64},Float64}(),
            Dict{Vector{Int64},Vector{Float64}}(),
            Dict{Vector{Int64},Float64}(),
            0,
            1,
            0
            )
    end
end

function structure(dbTempMatch,dbTempSim,dbOutput,time)
    matchStructure = matchHierarchy(dbTempMatch,dbTempSim,dbOutput,Float64(time))
    currentStructure!(matchStructure)
    potentialStructure!(matchStructure)
    phenoBiomass!(matchStructure)
    updateMatchesDB!(matchStructure)
    return matchStructure
end

function currentStructure!(matchStructure::matchHierarchy)
    time = matchStructure.time
    dbTempSim = matchStructure.dbSim
    dbTempMatch = matchStructure.dbMatch
    for (matchID,) in execute(dbTempTri,"SELECT DISTINCT match_id
            FROM vmatch_phenotypes WHERE t = $(time)
            ORDER BY vstrain_id")
        phenotype = [spacerID for (spacerID,) in
        execute(dbTempTri,"SELECT phenotype
        FROM vmatch_phenotypes
        WHERE t = $(time) AND match_id = $(matchID)
        ORDER BY phenotype")]
        # println("phenotype is $(phenotype)")
        identifyCurrentMatches!(matchStructure)
        push!(matchStructure.phenotypes.current,matchID)
        union!(matchStructure.matchtypes,length(phenotype))
    end
end


function identifyCurrentMatches!(
    matchStructure::matchHierarchy,phenotype::Vector{Int64},vstrain_id::Int64)
    time = matchStructure.time
    matchID = matchStructure.matchID
    dbTempSim = matchStructure.dbSim
    dbTempMatch = matchStructure.dbMatch

    matchStructure.sBstrains.current[matchID] = [Int64(strain)
        for (strain,) in execute(dbTempMatch, "SELECT bstrain_id
            FROM bstrain_to_vstrain_0matches
            WHERE t = $(time) AND vstrain_id = $(vstrain_id)")]
    matchStructure.sBstrains.potential[phenotype] =
        copy(matchStructure.sBstrains.current[phenotype])

    allBstrains = [Int64(strain)
        for (strain,) in execute(dbTempSim, "SELECT bstrain_id
            FROM babundance
            WHERE t = $(time)") if strain != 1]
    matchStructure.iBstrains.current[phenotype] =
        setdiff(allBstrains,matchStructure.sBstrains.current[phenotype])
    matchStructure.iBstrains.potential[phenotype] =
        copy(matchStructure.iBstrains.current[phenotype])

    if length(matchStructure.sBstrains.current[phenotype]) == 0
        matchStructure.sBiomass.current[phenotype] = 0
        matchStructure.sBiomass.potential[phenotype] =
            copy(matchStructure.sBiomass.current[phenotype])
    else
        matchStructure.sBiomass.current[phenotype] = sum([Int64(abund)
            for (abund,) in execute(dbTempSim, "SELECT abundance
                FROM babundance
                WHERE t = $(time) AND bstrain_id in
                ($(join(matchStructure.sBstrains.current[phenotype],", ")))")])
        matchStructure.sBiomass.potential[phenotype] =
            copy(matchStructure.sBiomass.current[phenotype])
    end
    if length(matchStructure.iBstrains.current[phenotype]) == 0
        matchStructure.iBiomass.current[phenotype] = 0
        matchStructure.iBiomass.potential[phenotype] =
            copy(matchStructure.sBiomass.current[phenotype])
    else
        matchStructure.iBiomass.current[phenotype] = sum([Int64(abund)
            for (abund,) in execute(dbTempSim, "SELECT abundance
                FROM babundance
                WHERE t = $(time) AND bstrain_id in
                ($(join(matchStructure.iBstrains.current[phenotype],", ")))")])
        matchStructure.iBiomass.potential[phenotype] =
            copy(matchStructure.iBiomass.current[phenotype])
    end
end

function potentialStructure!(matchStructure::matchHierarchy)
    time = matchStructure.time
    dbTempSim = matchStructure.dbSim
    dbTempMatch = matchStructure.dbMatch
    for phenotype in matchStructure.phenotypes.current
        # println("phenotype is $(phenotype)")
        if length(phenotype) == 1
            if !in(phenotype,matchStructure.phenotypes.potential)
                identifyPotentialMatches!(matchStructure,phenotype)
                push!(matchStructure.phenotypes.potential,phenotype)
            end
            continue
        end
        for k in collect(1:length(phenotype))
            # println("k is $(k)")
            for subPhenotype in collect(combinations(phenotype,k))
                subPheno = sort(subPhenotype)
                # println("subPhenotype is $(subPhenotype)")
                if !in(subPheno,matchStructure.phenotypes.potential)
                    identifyPotentialMatches!(matchStructure,subPheno)
                    push!(matchStructure.phenotypes.potential,subPheno)
                end
            end
        end
    end
end


# This just gets a (sub-)phenotype and checks for
# and logs susceptible and immune host strains
function identifyPotentialMatches!(
    matchStructure::matchHierarchy,subPhenotype::Vector{Int64})
    time = matchStructure.time
    dbTempSim = matchStructure.dbSim
    dbTempMatch = matchStructure.dbMatch

    matchStructure.iBstrains.potential[subPhenotype] = [Int64(strain)
        for (strain,) in execute(dbTempSim,
            "SELECT DISTINCT bstrain_id
            FROM bspacers
            WHERE spacer_id in ($(join(subPhenotype,", ")))")]
    allBstrains = [Int64(strain)
        for (strain,) in execute(dbTempSim, "SELECT bstrain_id
            FROM babundance
            WHERE t = $(time)")]
    matchStructure.sBstrains.potential[subPhenotype] =
        setdiff(allBstrains,matchStructure.iBstrains.potential[subPhenotype])
    if length(matchStructure.iBstrains.potential[subPhenotype]) == 0
        matchStructure.iBiomass.potential[subPhenotype] = 0
    else
        matchStructure.iBiomass.potential[subPhenotype] = sum([Int64(abund)
            for (abund,) in execute(dbTempSim, "SELECT abundance
                FROM babundance
                WHERE t = $(time) AND bstrain_id in
                ($(join(matchStructure.iBstrains.potential[subPhenotype],", ")))")])
    end
    if length(matchStructure.sBstrains.potential[subPhenotype]) == 0
        matchStructure.sBiomass.potential[subPhenotype] = 0
    else
        matchStructure.sBiomass.potential[subPhenotype] = sum([Int64(abund)
            for (abund,) in execute(dbTempSim, "SELECT abundance
                FROM babundance
                WHERE t = $(time) AND bstrain_id in
                ($(join(matchStructure.sBstrains.potential[subPhenotype],", ")))")])
    end
end


function phenoBiomass!(matchStructure::matchHierarchy)
    dbTempSim = matchStructure.dbSim
    time = matchStructure.time
    Vtotal = sum([abund for (abund,) in
                    execute(dbTempSim, "SELECT abundance
                        FROM vabundance WHERE t = $(time)")])
    for phenotype in matchStructure.phenotypes.current
        V = sum([Int64(abund) for (abund,) in execute(dbTempSim, "SELECT abundance
            FROM vabundance
            WHERE t = $(time) AND vstrain_id in
            ($(join(matchStructure.strainclasses[phenotype],", ")))")])
        matchStructure.phenoBiomass[phenotype] = V
        matchStructure.phenoFrequency[phenotype] = V/Vtotal
    end
    allVstrains = [strain for (strain,) in execute(dbTempSim,"SELECT vstrain_id
                    FROM vabundance WHERE t = $(matchStructure.time)")]
    match0 = setdiff(allVstrains,vcat(values(matchStructure.strainclasses)...))
    if length(match0) == 0
        matchStructure.phenoBiomass[[0]] = 0
        matchStructure.phenoFrequency[[0]] = 0
    else
        V = sum([Int64(abund) for (abund,) in execute(dbTempSim, "SELECT abundance
            FROM vabundance
            WHERE t = $(time) AND vstrain_id in
            ($(join(match0,", ")))")])
        matchStructure.phenoBiomass[[0]] = V
        matchStructure.phenoFrequency[[0]] = V/Vtotal
    end
end

function updateMatchesDB!(matchStructure::matchHierarchy)
    time = matchStructure.time
    dbTempSim = matchStructure.dbSim
    dbOutput = matchStructure.dbOutput
    match_id = 1
    execute(dbOutput, "INSERT INTO match_phenotype_abundances VALUES (?,?,?,?)",
        (time, 0,
        matchStructure.phenoBiomass[[0]],
        matchStructure.phenoFrequency[[0]]))
    for phenotype in matchStructure.phenotypes.current
        for spacer_id in phenotype
            execute(dbOutput, "INSERT INTO match_phenotypes VALUES (?,?,?)",
                (time, match_id,spacer_id))
        end
        execute(dbOutput, "INSERT INTO match_phenotype_abundances VALUES (?,?,?,?)",
            (time, match_id,
            matchStructure.phenoBiomass[phenotype],
            matchStructure.phenoFrequency[phenotype]))
        match_id += 1
    end
end

function computeR0!(matchStructure::matchHierarchy,Extinction::extinction)
    phi = matchStructure.parameters.adsorption_rate
    q = matchStructure.parameters.spacer_acquisition_prob
    beta = matchStructure.parameters.viral_burst_size
    d = matchStructure.parameters.viral_decay_rate
    N = sum([abund for (abund,) in
            execute(matchStructure.dbSim, "SELECT abundance
                FROM babundance WHERE t = $(matchStructure.time)")])
    # R0 = (phi*q*N+d)/(beta*phi*(1-q)*N)

    Extinction.birthRates[[0]] = beta*phi*(1-q)*N
    a = Extinction.birthRates[[0]]
    # println("a is $(a)")
    #
    Extinction.birthMutRates[[0]] = 0
    #
    Extinction.deathRates[[0]] = phi*q*N + d
    m = Extinction.deathRates[[0]]
    if a > 0
        Extinction.pRoots[[0]] =
            [1,m/a]
    else
        Extinction.pRoots[[0]] = [1]
    end
    println("pRoots of 0match is $(Extinction.pRoots[[0]])")
    Extinction.probability[[0]] =
        minimum(Extinction.pRoots[[0]])
    Extinction.R0 = 1/Extinction.probability[[0]]
    return Float64(Extinction.R0)
end


function computeExtinctionProbs(matchStructure::matchHierarchy)
    phi = matchStructure.parameters.adsorption_rate
    q = matchStructure.parameters.spacer_acquisition_prob
    beta = matchStructure.parameters.viral_burst_size
    d = matchStructure.parameters.viral_decay_rate
    mu = matchStructure.parameters.viral_mutation_rate
    N = sum([abund for (abund,) in
            execute(matchStructure.dbSim, "SELECT abundance
                FROM babundance WHERE t = $(matchStructure.time)")])
    # (V,) = execute(matchStructure.dbSim, "SELECT viral_abundance
    #         FROM summary WHERE t = $(matchStructure.time)")
    # V = V.viral_abundance
    Extinction = extinction(matchStructure)
    computeR0!(matchStructure,Extinction)
    R0 = Extinction.R0
    if length(matchStructure.matchtypes) == 0
        return Extinction
    end
    for matchtype in collect(1:maximum(matchStructure.matchtypes))
        # println("matchtype is $(matchtype)")
        Extinction.pointEscapeProb[matchtype] = mu #*(1-mu)^(matchtype+1)
        pE = Extinction.pointEscapeProb[matchtype]
        l = matchtype
        for phenotype in [pheno
                for pheno in matchStructure.phenotypes.potential
                    if length(pheno)==matchtype]
            # println("phenotype is $(phenotype)")
            if length(phenotype) == 1
                S = matchStructure.sBiomass.potential[phenotype]
                I = matchStructure.iBiomass.potential[phenotype]
                #
                Extinction.birthRates[phenotype] = beta*phi*(1-q)*(1-l*pE)*S
                a = round(Extinction.birthRates[phenotype],digits=6)
                #
                Extinction.birthMutRates[phenotype] = beta*phi*(1-q)*pE*S
                b = round(Extinction.birthMutRates[phenotype],digits=6)
                #
                Extinction.deathRates[phenotype] = phi*I + phi*q*S + d
                m = round(Extinction.deathRates[phenotype],digits=6)
                R = 1/R0
                R = round(R,digits=6)
                if a > 0
                    Extinction.pRoots[phenotype] =
                        [(a+b+m-b*R-sqrt(-4*a*m+(a+b+m-b*R)^2))/(2*a),
                            (a+b+m-b*R+sqrt(-4*a*m+(a+b+m-b*R)^2))/(2*a)]
                else
                    Extinction.pRoots[phenotype] =
                        [m/(b+m-b*R)]
                end
                # println("roots are $(Extinction.pRoots[phenotype])")
                # if isnan(minimum(Extinction.pRoots[phenotype]))
                #     Extinction.probability[phenotype] = 0
                # else
                    Extinction.probability[phenotype] =
                        minimum(Extinction.pRoots[phenotype])
                # end
            else
                S = matchStructure.sBiomass.potential[phenotype]
                I = matchStructure.iBiomass.potential[phenotype]
                #
                Extinction.birthRates[phenotype] = beta*phi*(1-q)*(1-l*pE)*S
                a = round(Extinction.birthRates[phenotype],digits=6)
                # println("a is $(a)")
                #
                Extinction.birthMutRates[phenotype] = beta*phi*(1-q)*pE*S
                b = round(Extinction.birthMutRates[phenotype],digits=6)
                #
                Extinction.deathRates[phenotype] = phi*I + phi*q*S + d
                m = round(Extinction.deathRates[phenotype],digits=6)
                R = 0
                for pheno in collect(combinations(phenotype,length(phenotype)-1))
                    # print("subPheno is $(pheno)")
                    R += Extinction.probability[sort(pheno)]
                end
                R = round(R,digits=6)
                # println("phenotype is $(phenotype)")
                # println("R is $(R)")
                # println("arg of sqrt is $(-4*a*m+(a+b+m-b*R)^2)")
                # println("a is $(a)")
                # println("b is $(b)")
                # println("m is $(m)")
                # println("S is $(S)")
                # println("I is $(I)")
                if a > 0
                    Extinction.pRoots[phenotype] =
                        [(a+b+m-b*R-sqrt(-4*a*m+(a+b+m-b*R)^2))/(2*a),
                            (a+b+m-b*R+sqrt(-4*a*m+(a+b+m-b*R)^2))/(2*a)]
                else
                    Extinction.pRoots[phenotype] =
                        [m/(b+m-b*R)]
                end

                # println("roots are $(Extinction.pRoots[phenotype])")
                # if isnan(minimum(Extinction.pRoots[phenotype]))
                #     Extinction.probability[phenotype] = 0
                # else
                    if minimum(Extinction.pRoots[phenotype]) > 1 &&
                            minimum(Extinction.pRoots[phenotype]) < 1.05
                        Extinction.probability[phenotype] = 1
                    else
                        Extinction.probability[phenotype] =
                            minimum(Extinction.pRoots[phenotype])
                    end
                # end
            end
        end
    end
    return Extinction
end

function updateProbsDB!(match::Bool,time::Float64,match_id::Int64,
        phenotype::Vector{Int64}, Extinction::extinction,T::Bool)
    dbOutput = Extinction.dbOutput
    if match
        execute(dbOutput, "INSERT INTO pExtinction VALUES (?,?,?)",
            (time, match_id, Extinction.probability[phenotype]))
        execute(dbOutput, "INSERT INTO Rweighted VALUES (?,?,?)",
            (time, match_id, Extinction.Rweight[phenotype]))
        execute(dbOutput, "INSERT INTO pExtWeighted VALUES (?,?,?)",
            (time, match_id,Extinction.pWeight[phenotype]))
        if T
            execute(dbOutput, "INSERT INTO Ravg VALUES (?,?)",
                (time, Extinction.Ravg))
            execute(dbOutput, "INSERT INTO R0 VALUES (?,?)",
                (time, Extinction.R0))
            execute(dbOutput, "INSERT INTO pExtWeighted VALUES (?,?,?)",
                (time, 0,Extinction.pWeight[[0]]))
            execute(dbOutput, "INSERT INTO Rweighted VALUES (?,?,?)",
                (time, 0, Extinction.Rweight[[0]]))
            execute(dbOutput, "INSERT INTO pEmergenceProd VALUES (?,?)",
                (time, 1-Extinction.product))
            execute(dbOutput, "INSERT INTO pEmergenceSum VALUES (?,?)",
                (time, 1-Extinction.sum))
        end
        return  Int64(1)
    else
        execute(dbOutput, "INSERT INTO pExtinction VALUES (?,?,?)",
            (time, 0, 1/Extinction.R0))
        execute(dbOutput, "INSERT INTO Ravg VALUES (?,?)",
            (time, Extinction.Ravg))
        execute(dbOutput, "INSERT INTO R0 VALUES (?,?)",
            (time, Extinction.R0))
        execute(dbOutput, "INSERT INTO pExtWeighted VALUES (?,?,?)",
            (time, 0,Extinction.pWeight[[0]]))
        execute(dbOutput, "INSERT INTO Rweighted VALUES (?,?,?)",
            (time, 0, Extinction.Rweight[[0]]))
        execute(dbOutput, "INSERT INTO pEmergenceProd VALUES (?,?)",
            (time, 1-Extinction.product))
        execute(dbOutput, "INSERT INTO pEmergenceSum VALUES (?,?)",
            (time, 1-Extinction.sum))
    end
end

function assembleProbs!(
        matchStructure::matchHierarchy,Extinction::extinction)
    dbTempSim = matchStructure.dbSim
    time = matchStructure.time
    match_id = 1
    last = false
    Extinction.Ravg = matchStructure.phenoFrequency[[0]]*Extinction.R0
    Extinction.Rweight[[0]] = matchStructure.phenoFrequency[[0]]*Extinction.R0
    Extinction.pWeight[[0]] = matchStructure.phenoFrequency[[0]]*(1/Extinction.R0)
    Extinction.product = (1/Extinction.R0)^(matchStructure.phenoBiomass[[0]])
    Extinction.sum = matchStructure.phenoFrequency[[0]]*(1/Extinction.R0)
    if length(matchStructure.phenotypes.current) == 0
        updateProbsDB!(false,time,match_id,Vector{Int64}(),Extinction,last)
        return
    end
    for phenotype in matchStructure.phenotypes.current
        if phenotype == matchStructure.phenotypes.current[end]
            last = true
        end
        V = matchStructure.phenoBiomass[phenotype]
        f = matchStructure.phenoFrequency[phenotype]
        if isnan(Extinction.probability[phenotype]) || Extinction.probability[phenotype] == 0
            println("Rweight of $(phenotype) is undefined; either Inf or 0*Inf!")
            Extinction.pWeight[phenotype] = 0
            Extinction.Rweight[phenotype] = NaN
            Extinction.Ravg = NaN
            return
        else
            Extinction.Rweight[phenotype] = f*1/(Extinction.probability[phenotype])
            Extinction.Ravg += Extinction.Rweight[phenotype]
        end
        Extinction.pWeight[phenotype] = f*(Extinction.probability[phenotype])
        Extinction.sum += Extinction.pWeight[phenotype]
        Extinction.product = (Extinction.product)*(Extinction.probability[phenotype])^V
        match_id += updateProbsDB!(true,time,match_id,phenotype,Extinction,last)
    end
end

function phenotypeStats!(matchStructure::matchHierarchy)
    Extinction = computeExtinctionProbs(matchStructure)
    assembleProbs!(matchStructure,Extinction)
    return Extinction
end

for (t,) in execute(dbTempSim,"SELECT DISTINCT t FROM vabundance")
    if t == 0
         continue
    end
    println("Computing structure and probabilities at time = $(t)")
    matchStructure = structure(dbTempMatch,dbTempSim,dbOutput,t)
    phenotypeStats!(matchStructure)
end


function createindices()
    println("(Creating run_id indices...)")
    db = SQLite.DB(dbOutputPath)
    execute(db, "BEGIN TRANSACTION")
    for (table_name,) in execute(
        db, "SELECT name FROM sqlite_schema
        WHERE type='table' ORDER BY name;")
        # cols = [info.name for info in execute(db,"PRAGMA table_info($(table_name))")]
        if in(table_name,["bmatch_phenotypes","vmatch_phenotypes"])
            execute(db, "CREATE INDEX $(table_name)_index ON $(table_name) (match_id)")
        end
        if in(table_name,["single_match_tripartite_networks"])
            execute(db, "CREATE INDEX $(table_name)_index ON $(table_name) (t)")
        end
        if in(table_name,["bmatch_phenotypes_singles","vmatch_phenotypes_singles",
                            "bmatches","vmatches","bmatches_abundance","vmatches_abundance"])
            execute(db, "CREATE INDEX $(table_name)_index ON $(table_name) (t, match_id)")
        end
        if in(table_name,["vmatches_susceptible_babundance","vmatches_susceptible_bstrains"])
            execute(db, "CREATE INDEX $(table_name)_index ON $(table_name) (t, vmatch_id)")
        end
    end
    execute(db, "COMMIT")
end

createindices()




println("Complete!")
