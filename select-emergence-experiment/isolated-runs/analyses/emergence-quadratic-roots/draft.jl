#!/usr/bin/env julia

println("(Julia compilation delay...)")

using SQLite
using DataFrames
using SQLite.DBInterface: execute
using SQLite: DB
using DataFrames
using Combinatorics: combinations

run_id = ARGS[1]

## Define Paths ##
SCRIPT_PATH = abspath(dirname(PROGRAM_FILE))
#
dbSimPath = joinpath(SCRIPT_PATH,"..","..","..","simulation","sweep_db_gathered.sqlite") # cluster
dbOutputPath = joinpath("probability-of-emergence_output.sqlite") # cluster

# dbSimPath = joinpath("/Volumes/Yadgah","crispr-sweep-7-2-2022/isolates/runID3297-c66-r47/runID3297-c66-r47.sqlite") # local
# dbMatchPath = joinpath("/Volumes/Yadgah","crispr-sweep-7-2-2022/isolates/runID3297-c66-r47/matches_output.sqlite") # local
# dbTriPath = joinpath("/Volumes/Yadgah","crispr-sweep-7-2-2022/isolates/runID3297-c66-r47/tripartite-networks_output.sqlite") # local
# dbOutputPath = joinpath("/Volumes/Yadgah/crispr-sweep-7-2-2022/isolates/runID3297-c66-r47/probability-of-emergence_output.sqlite") # local
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
dbTempTri = SQLite.DB(dbTriPath)
dbOutput = SQLite.DB(dbOutputPath)

execute(dbOutput, "CREATE TABLE pExtinction (t REAL, match_id INTEGER, p REAL)")
execute(dbOutput, "CREATE TABLE pEmergenceSum (t REAL, p REAL)")
execute(dbOutput, "CREATE TABLE pExtWeighted
    (t REAL, vmatch_id INTEGER, p REAL)")
execute(dbOutput, "CREATE TABLE R0 (t REAL, R REAL)")

dbTempSim = SQLite.DB()
# dbTempSim = SQLite.DB(dbSimPath)
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
execute(dbTempSim, "CREATE INDEX babundance_index ON babundance (t,bstrain_id)")
execute(dbTempSim, "CREATE INDEX vabundance_index ON vabundance (t,vstrain_id)")
execute(dbTempSim, "CREATE INDEX bspacers_index ON bspacers (bstrain_id,spacer_id)")
execute(dbTempSim, "COMMIT")

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

mutable struct hierarchy
    dbTri::DB
    dbMatch::DB
    dbSim::DB
    dbOutput::DB
    matchID::Int64
    matchtypes::Dict{Int64,Vector{Int64}}
    matches::Dict{Int64,Vector{Int64}}
    escapes::Dict{Int64,Vector{Int64}}
    lineages::Dict{Int64,Vector{Int64}}
    function hierarchy(dbTempTri::DB,dbTempMatch::DB,dbTempSim::DB,dbOutput::DB)
        new(
            dbTempTri,
            dbTempMatch,
            dbTempSim,
            dbOutput,
            Int64(0),
            Dict{Int64,Vector{Int64}}(),
            Dict{Int64,Vector{Int64}}(),
            Dict{Int64,Vector{Int64}}(),
            Dict{Int64,Vector{Int64}}()
            )
    end
end

mutable struct tstructure
    dbTri::DB
    dbMatch::DB
    dbSim::DB
    dbOutput::DB
    time::Float64
    matches::Vector{Int64}
    escapes::Vector{Int64}
    sBstrains::Dict{Int64,Vector{Int64}}
    iBstrains::Dict{Int64,Vector{Int64}}
    sBiomass::Dict{Int64,Int64}
    iBiomass::Dict{Int64,Int64}
    function tstructure(matchStructure::hierarchy,time::Float64)
        new(
            matchStructure.dbTri,
            matchStructure.dbMatch,
            matchStructure.dbSim,
            matchStructure.dbOutput,
            time,
            Vector{Int64}(),
            Vector{Int64}([0]),
            Dict{Int64,Vector{Int64}}(),
            Dict{Int64,Vector{Int64}}(),
            Dict{Int64,Int64}(),
            Dict{Int64,Int64}()
            )
    end
end

mutable struct invasion
    time::Float64
    birth::Dict{Int64,Float64}
    death::Dict{Int64,Float64}
    mut::Dict{Int64,Float64}
    probability::Dict{Int64,Float64}
    proots::Dict{Int64,Vector{Float64}}
    function invasion(current::tstructure)
        new(
            current.time,
            Int64(0),
            Dict{Int64,Float64}(),
            Dict{Int64,Float64}(),
            Dict{Int64,Float64}(),
            Dict{Int64,Float64}(),
            Dict{Int64,Vector{Float64}}()
            )
    end
end

function find0matches(dbMatchPath,dbTriPath)
    dbTempMatch = SQLite.DB(dbMatchPath)
    dbTempTri = SQLite.DB(dbTriPath)
    times = Vector{Float64}()
    allv0matches = Vector{Int64}()
    abundances = Vector{Int64}()
    for (t,) in execute(dbTempSim,"SELECT DISTINCT t FROM vabundance")
        v0matches = [vstrainID for (vstrainID,) in
            execute(dbTempMatch,"SELECT DISTINCT vstrain_id
                FROM bstrain_to_vstrain_0matches WHERE t = $(t)")]
        vmatches = [vstrainID for (vstrainID,) in
            execute(dbTempMatch,"SELECT DISTINCT vstrain_id
                FROM bstrain_to_vstrain_matches WHERE t = $(t)")]
        v0matches = setdiff(v0matches,vmatches)
        if length(v0matches) != 0
            push!(times,repeat([Float64(t)],length(v0matches))...)
            push!(allv0matches,v0matches...)
            V = sum([Int64(abund) for (abund,) in execute(dbTempSim, "SELECT abundance
                FROM vabundance
                WHERE t = $(t) AND vstrain_id in
                ($(join(v0matches,", ")))")])
            push!(abundances,V)
        end
    end
    v0 = DataFrame(t = times, vstrain_id = allv0matches)
    v0 |> SQLite.load!(dbTempTri,
                                "v0matches",ifnotexists=true)
    v0 = DataFrame(t = unique(times), vabundance = abundances)
    v0 |> SQLite.load!(dbTempTri,
                                "v0matches_abundances",ifnotexists=true)
    execute(dbTempTri, "CREATE INDEX v0matches_index ON v0matches (t)")
    execute(dbTempTri, "CREATE INDEX v0matches_abundances_index
                        ON v0matches_abundances (t)")
end

function matchLengths(dbTriPath)
    dbTempTri = SQLite.DB(dbTriPath)
    matchIDs = Vector{Int64}()
    matchLengths = Vector{Int64}()
    for (matchID,) in execute(dbTempTri,"SELECT DISTINCT match_id FROM vmatch_phenotypes")
        matchLength = length([spacerID for (spacerID,) in
            execute(dbTempTri,"SELECT phenotype
                FROM vmatch_phenotypes WHERE match_id = $(matchID)")])
        push!(matchIDs,matchID)
        push!(matchLengths,matchLength)
    end
    mLengths = DataFrame(match_id = matchIDs, match_length = matchLengths)
    mLengths |> SQLite.load!(dbTempTri,
                                "vmatch_lengths",ifnotexists=true)

    execute(dbTempTri, "CREATE INDEX vmatch_lengths_index
                        ON vmatch_lengths (match_id,match_length)")
end

function potentialStructure!(matchStructure::hierarchy)
    dbTempSim = matchStructure.dbSim
    dbTempMatch = matchStructure.dbMatch
    dbTempTri = matchStructure.dbTri
    matchIDs = Dict(Vector([Int64(spacerID) for (spacerID,) in
                    execute(dbTempTri,"SELECT phenotype
                    FROM vmatch_phenotypes
                    WHERE match_id = $(matchID)
                    ORDER BY phenotype")]) => Int64(matchID)
                        for matchID in
                            [matchID for (matchID,) in
                            execute(dbTempTri,"SELECT DISTINCT match_id
                            FROM vmatch_phenotypes ORDER BY match_id")]
                        )
    matchStructure.matchID = maximum(values(matchIDs)) + 1
    phenotypes = keys(matchIDs)
    lineages = Dict{Int64,Vector{Vector{Int64}}}()
    println("Generating all possible match phenotypes from...")
    for phenotype in phenotypes
        println("existing phenotype $(phenotype)")
        matchID = matchIDs[phenotype]
        mlength = length(phenotype)
        for k in collect(1:mlength-1)
            for subPhenotype in collect(combinations(phenotype,k))
                subPhenotype = sort(subPhenotype)
                if !in(subPhenotype,keys(matchIDs))
                    matchIDs[subPhenotype] = matchStructure.matchID
                    matchStructure.matchID += 1
                end
            end
        end
    end
    for phenotype in keys(matchIDs)
        println("Compiling single escapes for $(phenotype)...")
        if in(length(phenotype),keys(matchStructure.matchtypes))
            push!(matchStructure.matchtypes[length(phenotype)],matchIDs[phenotype])
        else
            matchStructure.matchtypes[length(phenotype)] = Vector{Int64}()
            push!(matchStructure.matchtypes[length(phenotype)],matchIDs[phenotype])
        end
        matchStructure.matches[matchIDs[phenotype]] = phenotype
        matchID = matchIDs[phenotype]

        if length(phenotype) == 1
            matchStructure.escapes[matchIDs[phenotype]] = [0]
        else
            matchStructure.escapes[matchID] = Vector{Int64}()
            for subPhenotype in collect(combinations(phenotype,length(phenotype)-1))
                push!(matchStructure.escapes[matchID],matchIDs[sort(subPhenotype)])
            end
        end
    end
    return matchIDs
end

function findEscapeLineages!(matchStructure::hierarchy,matchIDs::Dict{Vector{Int64},Int64})
    dbTempSim = matchStructure.dbSim
    dbTempMatch = matchStructure.dbMatch
    dbTempTri = matchStructure.dbTri
    numPhenos = length(keys(matchStructure.matches))
    for matchID in keys(matchStructure.matches)
        phenotype = matchStructure.matches[matchID]
        println("Compiling complete lineage for $(phenotype)...")
        matchStructure.lineages[matchID] = Vector{Int64}()
        mlength = length(phenotype)
        for k in collect(1:mlength-1)
            for subPhenotype in collect(combinations(phenotype,k))
                subPhenotype = sort(subPhenotype)
                push!(matchStructure.lineages[matchID],matchIDs[subPhenotype])
            end
        end
        numPhenos -= 1
        println("$(numPhenos) phenotypes left to search")
    end
end

function microbeBiomass!(matchStructure::hierarchy,current::tstructure)
    time = current.time
    dbTempSim = matchStructure.dbSim
    dbTempTri = matchStructure.dbTri
    current.matches = [Int64(matchID) for (matchID,) in
                                execute(dbTempTri,"SELECT DISTINCT match_id
                                FROM vmatches WHERE t = $(time)
                                ORDER BY match_id")]
    matchIDs = current.matches
    union!(current.escapes,matchIDs)
    for matchID in matchIDs
        union!(current.escapes,matchStructure.lineages[matchID])
    end
    tv0 = [t for (t,) in execute(dbTempTri,"SELECT t FROM v0matches_abundances")]
    if in(time,tv0)
        current.matches = vcat(Vector{Int64}([0]),current.matches)
    end
    bstrains = [Int64(strain)
                    for (strain,) in execute(dbTempSim, "SELECT bstrain_id
                        FROM babundance
                        WHERE t = $(time)")]
    spacerIDs = Vector{Int64}()
    for matchID in matchIDs
        union!(spacerIDs,matchStructure.matches[matchID])
    end
    for spacerID in spacerIDs
        current.iBstrains[spacerID] = [Int64(strain)
            for (strain,) in execute(dbTempSim,
                "SELECT DISTINCT bstrain_id
                FROM bspacers
                WHERE spacer_id = $(spacerID)")]
    end
    numEscapes = length(current.escapes)
    for matchID in current.escapes
        println("Computing abundances susceptible and immune to phenotype $(matchID)")
        if matchID == 0
            B = sum([Int64(abund) for (abund,)
                        in execute(dbTempSim, "SELECT abundance
                        FROM babundance WHERE t = $(time)")])
            current.sBiomass[0] = B
            current.iBiomass[0] = 0
            numEscapes -= 1
            println("$(numEscapes) phenotypes left")
            continue
        end
        immuneStrains = Vector{Int64}()
        for spacerID in matchStructure.matches[matchID]
            union!(immuneStrains,current.iBstrains[spacerID])
        end

        current.iBiomass[matchID] =
            sum([Int64(abund) for (abund,)
                in execute(dbTempSim, "SELECT abundance
                    FROM babundance
                    WHERE t = $(time) AND bstrain_id in
                    ($(join(immuneStrains,", ")))")])

        current.sBstrains[matchID] = setdiff(bstrains,immuneStrains)
        if length(current.sBstrains[matchID]) > 0
            current.sBiomass[matchID] =
                sum([Int64(abund) for (abund,)
                    in execute(dbTempSim, "SELECT abundance
                        FROM babundance
                        WHERE t = $(time) AND bstrain_id in
                        ($(join(current.sBstrains[matchID],", ")))")])
        else
            current.sBiomass[matchID] = 0
        end
        numEscapes -= 1
        println("$(numEscapes) phenotypes left")
    end
end


function computeR0!(extinction::invasion,parameters::params)
    time = extinction.time
    phi = parameters.adsorption_rate
    q = parameters.spacer_acquisition_prob
    beta = parameters.viral_burst_size
    d = parameters.viral_decay_rate
    N = sum([abund for (abund,) in
            execute(matchStructure.dbSim, "SELECT abundance
                FROM babundance WHERE t = $(time)")])
    extinction.birth[0] = beta*phi*(1-q)*N
    extinction.mut[0] = 0
    extinction.death[0] = phi*q*N + d
    a = extinction.birth[0]
    m = extinction.death[0]
    if a > 0
        extinction.proots[0] =
            Vector{Float64}([1,m/a])
    else
        extinction.proots[0] = Vector{Float64}([1])
    end
    println("proots of 0match is $(extinction.proots[0])")
    extinction.probability[0] =
        minimum(extinction.proots[0])
    # R0 = (phi*q*N+d)/(beta*phi*(1-q)*N)
    # extinction.R0 = 1/extinction.probability[0]
    # return Float64(extinction.R0)
end


function computeProbabilities(matchStructure::hierarchy,
                                current::tstructure,parameters::params)
    time = current.time
    phi = parameters.adsorption_rate
    q = parameters.spacer_acquisition_prob
    beta = parameters.viral_burst_size
    d = parameters.viral_decay_rate
    mu = parameters.viral_mutation_rate
    sigdig = 7
    N = sum([abund for (abund,) in
            execute(matchStructure.dbSim, "SELECT abundance
                FROM babundance WHERE t = $(time)")])

    extinction = invasion(current)
    computeR0!(extinction,parameters)

    if current.matches == [0]
        println("this works")
        maxLength = 0
    else
        maxLength = maximum([Int64(mL) for (mL,)
            in execute(dbTempTri, "SELECT match_length
                FROM vmatch_lengths WHERE match_id in
                ($(join(current.matches,", ")))")])
    end

    for matchtype in collect(1:maxLength)
        # println("matchtype is $(matchtype)")
        l = matchtype
        println("Computing extinction probability for match lengths of $(matchtype)...")
        matchIDs = intersect(matchStructure.matchtypes[matchtype],current.escapes)
        for matchID in matchIDs
            # println("phenotype is $(phenotype)")
            if length(matchStructure.matches[matchID]) == 1
                S = current.sBiomass[matchID]
                I = current.iBiomass[matchID]

                extinction.birth[matchID] = beta*phi*(1-q)*(1-l*mu)*S
                a = extinction.birth[matchID]
                # a = round(extinction.birth[matchID],digits=sigdig)

                extinction.mut[matchID] = beta*phi*(1-q)*l*mu*S
                b = extinction.mut[matchID]
                # b = round(extinction.mut[matchID],digits=sigdig)

                extinction.death[matchID] = phi*I + phi*q*S + d
                c = extinction.death[matchID]
                # c = round(extinction.death[matchID],digits=sigdig)

                Qsum = extinction.probability[0]
                # Qsum = round(Qsum,digits=6)

                if a !=  0.0
                    if -4*a*c+(a+b+c-b*Qsum)^2 < 0.0
                        extinction.proots[matchID] =
                        map(x->round(x,digits=sigdig), [(a+b+c-b*Qsum)/(2*a)])
                    println("proots of match $(matchID)
                        is $(extinction.proots[matchID])")
                    else
                        extinction.proots[matchID] =
                            map(x->round(x,digits=sigdig),
                                [(a+b+c-b*Qsum-sqrt(-4*a*c+(a+b+c-b*Qsum)^2))/(2*a),
                                (a+b+c-b*Qsum+sqrt(-4*a*c+(a+b+c-b*Qsum)^2))/(2*a)])
                        println("proots of match $(matchID)
                            is $(extinction.proots[matchID])")
                    end
                else
                    extinction.proots[matchID] =
                        map(x->round(x,digits=sigdig),[c/(b+c-b*Qsum)])
                    println("proots of match $(matchID)
                        is $(extinction.proots[matchID])")
                end

                extinction.probability[matchID] =
                    minimum(extinction.proots[matchID])
            else
                S = current.sBiomass[matchID]
                I = current.iBiomass[matchID]

                extinction.birth[matchID] = beta*phi*(1-q)*(1-l*mu)*S
                a = extinction.birth[matchID]
                # a = round(extinction.birth[matchID],digits=sigdig)

                extinction.mut[matchID] = beta*phi*(1-q)*l*mu*S
                b = extinction.mut[matchID]
                # b = round(extinction.mut[matchID],digits=sigdig)

                extinction.death[matchID] = phi*I + phi*q*S + d
                c = extinction.death[matchID]
                # c = round(extinction.death[matchID],digits=sigdig)

                Qsum = 0
                for escapeID in matchStructure.escapes[matchID]
                    Qsum += extinction.probability[escapeID]
                end
                # Qsum = round(Qsum,digits=sigdig)

                if a != 0.0
                    if -4*a*c+(a+b+c-b*Qsum)^2 < 0.0
                        extinction.proots[matchID] =
                        map(x->round(x,digits=sigdig), [(a+b+c-b*Qsum)/(2*a)])
                    println("proots of match $(matchID)
                        is $(extinction.proots[matchID])")
                    else
                        extinction.proots[matchID] =
                            map(x->round(x,digits=sigdig),
                                [(a+b+c-b*Qsum-sqrt(-4*a*c+(a+b+c-b*Qsum)^2))/(2*a),
                                (a+b+c-b*Qsum+sqrt(-4*a*c+(a+b+c-b*Qsum)^2))/(2*a)])
                        println("proots of match $(matchID)
                            is $(extinction.proots[matchID])")
                    end
                else
                    extinction.proots[matchID] =
                        map(x->round(x,digits=sigdig),[c/(b+c-b*Qsum)])
                    println("proots of match $(matchID)
                        is $(extinction.proots[matchID])")
                end

                if minimum(extinction.proots[matchID]) > 1
                    extinction.probability[matchID] = 1
                else
                    extinction.probability[matchID] =
                        minimum(extinction.proots[matchID])
                end
            end
        end
    end
    return extinction
end


function assemble(current::tstructure,extinction::invasion)
    time = current.time
    dbTempSim = current.dbSim
    dbTempTri = current.dbTri
    dbOutput = current.dbOutput

    pextinction = map(x->extinction.probability[x],current.matches)

    matchAbunds =  [Int64(abund) for (abund,)
        in execute(dbTempTri, "SELECT vabundance
            FROM vmatches_abundances
            WHERE t = $(time) AND match_id in
            ($(join(current.matches,", "))) ORDER BY match_id")]

    if in(0,current.matches)
        matchAbunds = vcat(sum([Int64(abund) for (abund,)
            in execute(dbTempTri, "SELECT vabundance
                FROM v0matches_abundances
                WHERE t = $(time)")]),matchAbunds)
    end

    V = sum([Int64(abund) for (abund,)
        in execute(dbTempSim, "SELECT abundance
            FROM vabundance
            WHERE t = $(time)")])

    DataFrame(  t = repeat([Float64(time)],length(current.matches)),
                vmatch_id = current.matches, p_extinction = pextinction,
                p_ext_weighted = matchAbunds.*pextinction*1/V,
                p_emerge_weighted = matchAbunds.*(ones(length(pextinction))-pextinction)*1/V,
                vfrequency = matchAbunds*1/V,
                vabundance = matchAbunds
                ) |> SQLite.load!(dbOutput, "vmatch_extinction",ifnotexists=true)

    return sum(matchAbunds.*pextinction*1/V)
end


function emergence(matchStructure::hierarchy,
                    current::tstructure,parameters::params)
    extinction = computeProbabilities(matchStructure,current,parameters)
    expectation = assemble(current,extinction)
    return expectation
end


function escape()
    matchStructure = hierarchy(dbTempTri,dbTempMatch,dbTempSim,dbOutput)
    parameters = params(matchStructure.dbSim)
    matchIDs = potentialStructure!(matchStructure)
    findEscapeLineages!(matchStructure,matchIDs)
    pExpectation = Vector{Float64}()
    times = Vector{Float64}()
    for (t,) in execute(dbTempSim,"SELECT DISTINCT t FROM vabundance")
        if t == 166
            return
        end
        push!(times,t)
        println("Computing structure and probabilities at time = $(t)")
        current = tstructure(matchStructure,Float64(t))
        microbeBiomass!(matchStructure,current)
        push!(pExpectation, emergence(matchStructure,current,parameters))
    end
    DataFrame(  t = times,
                p_emerge_expected = ones(length(times))-pExpectation,
                ) |> SQLite.load!(dbOutput, "vmatch_expectation",ifnotexists=true)
end

function createindices()
    println("(Creating run_id indices...)")
    db = SQLite.DB(dbOutputPath)
    execute(db, "BEGIN TRANSACTION")
    for (table_name,) in execute(
        db, "SELECT name FROM sqlite_schema
        WHERE type='table' ORDER BY name;")
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


find0matches(dbMatchPath,dbTriPath)
matchLengths(dbTriPath)

escape()
createindices()
println("Complete!")
