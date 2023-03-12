#!/usr/bin/env julia

println("(Julia compilation delay...)")

using SQLite
using DataFrames
using SQLite.DBInterface: execute
using SQLite: DB
using DataFrames
using Combinatorics: combinations
using LambertW
using SpecialFunctions

run_id = ARGS[1]

## Define Paths ##
SCRIPT_PATH = abspath(dirname(PROGRAM_FILE))
#
dbSimPath = joinpath(SCRIPT_PATH,"..","..","..","simulation","sweep_db_gathered.sqlite") # cluster
dbOutputPath = joinpath("emergence-lambert-root_output.sqlite") # cluster

# dbSimPath = joinpath("/Volumes/Yadgah","crispr-sweep-7-2-2022/isolates/runID3297-c66-r47/runID3297-c66-r47.sqlite") # local
# dbMatchPath = joinpath("/Volumes/Yadgah","crispr-sweep-7-2-2022/isolates/runID3297-c66-r47/matches_output.sqlite") # local
# dbTriPath = joinpath("/Volumes/Yadgah","crispr-sweep-7-2-2022/isolates/runID3297-c66-r47/tripartite-networks_output.sqlite") # local
# dbOutputPath = joinpath("/Volumes/Yadgah/crispr-sweep-7-2-2022/isolates/runID3297-c66-r47/emergence-lambert-root_output.sqlite") # local
if isfile(dbOutputPath)
    error("emergence-lambert-root.sqlite already exists; delete first")
end
##
dbTempSim = SQLite.DB(dbSimPath)
(cr,) = execute(dbTempSim,"SELECT combo_id,replicate FROM runs WHERE run_id = $(run_id)")
dbMatchPath = joinpath(SCRIPT_PATH,"..","..","isolates",
    "runID$(run_id)-c$(cr.combo_id)-r$(cr.replicate)","matches_output.sqlite")
if !isfile(dbMatchPath)
    error("matches_output.sqlite does not exist; compute first")
end
dbTriPath = joinpath(SCRIPT_PATH,"..","..","isolates",
    "runID$(run_id)-c$(cr.combo_id)-r$(cr.replicate)","tripartite-networks_output.sqlite")
if !isfile(dbMatchPath)
    error("tripartite-networks.sqlite does not exist; compute first")
end
dbTempMatch = SQLite.DB(dbMatchPath)
dbTempTri = SQLite.DB(dbTriPath)
dbOutput = SQLite.DB(dbOutputPath)

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
execute(dbTempSim, "CREATE INDEX babundance_index ON babundance (t,bstrain_id)")
execute(dbTempSim, "CREATE INDEX vabundance_index ON vabundance (t,vstrain_id)")
execute(dbTempSim, "CREATE INDEX bspacers_index ON bspacers (bstrain_id,spacer_id)")
execute(dbTempSim, "COMMIT")

include(joinpath(SCRIPT_PATH,"structures-integration.jl"))
include(joinpath(SCRIPT_PATH,"computeActualRoot.jl"))
include(joinpath(SCRIPT_PATH,"integrateExtinction.jl"))

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

function potentialStructure!(matchStructure::intHierarchy)
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
    # println("Generating all possible match phenotypes from...")
    # This loop gives all escape phenotypes (subphenotypes) of each existing
    # phenotype in the dynamics a match ID
    for phenotype in phenotypes
        # println("existing phenotype $(phenotype)")
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
    # This loop saves the lengths of all match types as a Dictionary (matchLength:[matchIDs]
    # also saves escape match IDs as a Dictionary (matchID0:[matchIDEscape]).
    # # Note that this does not save the actual phenotype
    for phenotype in keys(matchIDs)
        # println("Compiling single escapes for $(phenotype)...")
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
    findImmuneStrains!(matchStructure)
    # This saves all potential match phenotypes to tripartite database
    if !in("potential_vmatch_phenotypes",
                [table_name for (table_name,) in execute(dbTempTri,
                "SELECT name FROM sqlite_schema
                WHERE type='table' ORDER BY name;")])
        matchidlist = Vector{Int64}()
        phenolist = Vector{Int64}()
        for matchid in sort([keys(matchStructure.matches)...])
            append!(phenolist,matchStructure.matches[matchid])
            append!(matchidlist,
                fill(matchid,length(matchStructure.matches[matchid]))
                    )
        end
        DataFrame(match_id = matchidlist, phenotype = phenolist) |>
            SQLite.load!(dbTempTri,
                "potential_vmatch_phenotypes",ifnotexists=true)
        execute(dbTempTri, "CREATE INDEX potential_vmatch_phenotypes_index
                            ON potential_vmatch_phenotypes (match_id)")
    end
    if !in("potential_vmatch_single_escapes",
                [table_name for (table_name,) in execute(dbTempTri,
                "SELECT name FROM sqlite_schema
                WHERE type='table' ORDER BY name;")])
        matchidlist = Vector{Int64}()
        phenolist = Vector{Int64}()
        for matchid in sort([keys(matchStructure.escapes)...])
            append!(phenolist,sort([matchStructure.escapes[matchid]...]))
            append!(matchidlist,
                fill(matchid,length(matchStructure.escapes[matchid]))
                    )
        end
        DataFrame(match_id = matchidlist, escape_match_id = phenolist) |>
            SQLite.load!(dbTempTri,
                "potential_vmatch_single_escapes",ifnotexists=true)
        execute(dbTempTri, "CREATE INDEX potential_vmatch_single_escapes_index
                            ON potential_vmatch_single_escapes (match_id)")
    end
    return matchIDs
end
# This functipn saves all possible subphenotypes (not just escapes) as a Dictionary (matchID0:[matchIDsubphenotypes])
# this is what I call a "match lineage"
function findEscapeLineages!(matchStructure::intHierarchy,matchIDs::Dict{Vector{Int64},Int64})
    dbTempSim = matchStructure.dbSim
    dbTempMatch = matchStructure.dbMatch
    dbTempTri = matchStructure.dbTri
    numPhenos = length(keys(matchStructure.matches))
    for matchID in keys(matchStructure.matches)
        phenotype = matchStructure.matches[matchID]
        # println("Compiling complete lineage for $(phenotype)...")
        matchStructure.lineages[matchID] = Vector{Int64}()
        mlength = length(phenotype)
        for k in collect(1:mlength-1)
            for subPhenotype in collect(combinations(phenotype,k))
                subPhenotype = sort(subPhenotype)
                push!(matchStructure.lineages[matchID],matchIDs[subPhenotype])
            end
        end
        # numPhenos -= 1
        # println("$(numPhenos) phenotypes left to search")
    end
end

function findImmuneStrains!(matchStructure::intHierarchy)
    dbTempSim = matchStructure.dbSim
    dbTempTri = matchStructure.dbTri
    time = current.time

    bstrains = [Int64(strain)
                    for (strain,) in execute(dbTempSim, "SELECT bstrain_id
                        FROM babundance
                        WHERE t = $(time)")]
    spacerIDs = [Int64(spacerID) for (spacerID,) in execute(dbTempSim,
                "SELECT DISTINCT spacer_id
                FROM bspacers ORDER BY spacer_id")]
    for spacerID in spacerIDs
        matchStructure.ispacers[spacerID] =
                    [Int64(strainID) for (strainID,) in execute(dbTempSim,
                    "SELECT DISTINCT bstrain_id
                    FROM bspacers WHERE spacer_id = $(spacerID)
                    ORDER BY bstrain_id")]
    end
    for matchID in [keys(matchStructure.matches)...]
        for pheno in matchStructure.matches[matchID]
            union!(matchStructure.ibstrains[matchID],
                        matchStructure.ispacers[pheno])
        end
    end
end

function microbeBiomass!(matchStructure::intHierarchy,current::intState)
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
    spacerIDs = [Int64(spacerID) for (spacerID),) in execute(dbTempSim,
                "SELECT DISTINCT spacer_id
                FROM bspacers
                WHERE bstrain_id in
                ($(join(bstrains,", "))) ORDER BY spacer_id")]
    # println("current.iBstrains: $(current.iBstrains)")
    numEscapes = length(current.escapes)
    for matchID in [keys(matchStructure.matches)...]
        # println("Computing abundances susceptible and immune to phenotype $(matchID)")
        if matchID == 0
            B = sum([Int64(abund) for (abund,)
                        in execute(dbTempSim, "SELECT abundance
                        FROM babundance WHERE t = $(time)")])
            current.sBiomass[0] = B
            current.iBiomass[0] = 0
            # numEscapes -= 1
            # println("$(numEscapes) phenotypes left")
            continue
        end
        immuneStrains = intersect(matchStructure.ibstrains[matchID],bstrains)
        if length(immuneStrains) > 0
            current.iBiomass[matchID] =
                sum([Int64(abund) for (abund,)
                    in execute(dbTempSim, "SELECT abundance
                        FROM babundance
                        WHERE t = $(time) AND bstrain_id in
                        ($(join(immuneStrains,", ")))")])
        else
            current.iBiomass[matchID] = 0
        end

        sSstrains = setdiff(bstrains,immuneStrains)
        if length(sSstrains) > 0
            current.sBiomass[matchID] =
                sum([Int64(abund) for (abund,)
                    in execute(dbTempSim, "SELECT abundance
                        FROM babundance
                        WHERE t = $(time) AND bstrain_id in
                        ($(join(current.sBstrains[matchID],", ")))")])
        else
            current.sBiomass[matchID] = 0
        end
        # numEscapes -= 1
        # println("$(numEscapes) phenotypes left")
    end
    # sbiomasses = map(x->current.sBiomass[x],sort([current.escapes...]))
    # ibiomasses = map(x->current.iBiomass[x],sort([current.escapes...]))
    # DataFrame(t = fill(time,length(current.escapes)),
    #             match_id = sort([current.escapes...]),
    #             bsusceptible = sbiomasses,
    #             bimmune = ibiomasses
    #             ) |> SQLite.load!(dbTempTri,"potential_vmatches_babundances",
    #                     ifnotexists=true)
end

# This function assembles all pertinent data from objects and saves into database
function assemble(current::intState,extinction::intPcomponents)
    time = current.time
    dbTempSim = current.dbSim
    dbTempTri = current.dbTri
    dbOutput = current.dbOutput
    sigdig = extinction.sigdig

    pextinctionAct = map(x->extinction.pactual[x],current.matches)
    pextinctionInt = map(x->extinction.pintegrated[x],current.matches)
    lysis = map(x->extinction.lysis[x],current.matches)
    birth = map(x->extinction.birth[x],current.matches)
    death = map(x->extinction.death[x],current.matches)
    mut = map(x->extinction.mut[x],current.matches)

    # pextinctionEsc = map(x->extinction.prootmin[x],current.escapes)
    # birthEsc = map(x->extinction.birth[x],current.escapes)
    # deathEsc = map(x->extinction.death[x],current.escapes)
    # mutEsc = map(x->extinction.mut[x],current.escape)

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

    f = 1/V

    # DataFrame(  t = fill(Float64(time),length(current.matches)),
    #             vmatch_id = current.matches, p_extinction_root = pextinctionRoot,
    #             p_extinction_lambert = pextinctionLamb,
    #             lysis = lysis, birth = birth, death = death, mutation = mut,
    #             pactual_ext_weighted = matchAbunds.*pextinctionAct*f,
    #             pactual_emerge_weighted = matchAbunds.*(ones(length(pextinctionAct))-pextinctionAct)*f,
    #             pintegrated_ext_weighted = matchAbunds.*pextinctionInt*f,
    #             pintegrated_emerge_weighted = matchAbunds.*(ones(length(pextinctionInt))-pextinctionInt)*f,
    #             vfrequency = matchAbunds*f,
    #             vabundance = matchAbunds
    #             ) |> SQLite.load!(dbOutput, "existing_vmatch_extinction",ifnotexists=true)

    return sum(matchAbunds.*pextinctionAct*f),
            sum(matchAbunds.*pextinctionInt*f)
end

function emergence(matchStructure::intHierarchy,
                    current::state,pstate::probstate,parameters::params,dt::Float64)
    extinction = computeActualRoot!(matchStructure,current,parameters)
    integrateExtinction!(matchStructure,current,parameters,extinction,
                            pstate,dt)
    actualExpectation, intExpectation = assemble(current,extinction)
    return actualExpectation, intExpectation
end
# This function runs the whole analysis
function escape(pLambExpectation,pRootExpectation,pActualExpectation,pIntExpectation,matchStructure,matchIDs)
    # # This initializes the intHierarchy of matchtypes that can emerge from escape
    # matchStructure = intHierarchy(dbTempTri,dbTempMatch,dbTempSim,dbOutput)
    # # This calls and save the birth/death/mutation parameter rates of the models
    # parameters = params(matchStructure.dbSim)
    # matchIDs = potentialStructure!(matchStructure)
    # findEscapeLineages!(matchStructure,matchIDs)
    # pIntExpectation = Vector{Float64}()
    # pActualExpectation = Vector{Float64}()
    pstate = probstate()
    times = Vector{Float64}()
    t0 = 0
    for (t,) in execute(dbTempSim,"SELECT DISTINCT t FROM vabundance")
        if t == 150
            return
        end
        push!(times,t)
        println("Computing structure and probabilities at time = $(t)")
        current = intState(matchStructure,Float64(t))
        microbeBiomass!(matchStructure,current)
        # pinitial, intExpectation, lambertExpectation, rootExpectation =
        #         emergence(matchStructure,current,parameters)
        actualExpectation, integratedExpectation =
                emergence(matchStructure,current,pstate,parameters,t-t0)
        push!(pIntExpectation, integratedExpectation)
        push!(pActualExpectation, actualExpectation)
        t0 = t
    end
    # DataFrame(  t = times,
    #             p_emerge_actual_expected = ones(length(times))-pActualExpectation,
    #             p_emerge_integrated_expected = ones(length(times))-pIntExpectation
    #             ) |> SQLite.load!(dbOutput, "vmatch_expectation",ifnotexists=true)
end
# This checks for table names that indicate infection network exists in triparite database
# if not it modifies tripartite database to include infection network with associated abundances
# also modifies to include length of matches for each virus matchtype
function modifyTripartiteData()
    dbTempTri = SQLite.DB(dbTriPath)
    if !issubset(["v0matches","v0matches_abundances","vmatch_lengths"],[table_name for (table_name,) in execute(
        dbTempTri, "SELECT name FROM sqlite_schema
        WHERE type='table' ORDER BY name;")])

        find0matches(dbMatchPath,dbTriPath)
        matchLengths(dbTriPath)
    end
end

function createindices()
    println("(Creating run_id indices...)")
    db = SQLite.DB(dbOutputPath)
    execute(db, "BEGIN TRANSACTION")
    for (table_name,) in execute(
        db, "SELECT name FROM sqlite_schema
        WHERE type='table' ORDER BY name;")
        if in(table_name,["vmatch_expectation"])
            execute(db, "CREATE INDEX $(table_name)_index ON $(table_name) (t)")
        end
        if in(table_name,["existing_vmatch_extinction"])
            execute(db, "CREATE INDEX $(table_name)_index ON $(table_name) (t, vmatch_id)")
        end
    end
    execute(db, "COMMIT")
end


modifyTripartiteData()
escape()
createindices()
println("Complete!")
