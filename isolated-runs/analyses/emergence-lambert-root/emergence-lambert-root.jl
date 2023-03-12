#!/usr/bin/env julia

println("(Julia compilation delay...)")

using SQLite
using DataFrames
using SQLite.DBInterface: execute
using SQLite: DB
using DataFrames
using Combinatorics: combinations
using LambertW
using Roots

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
# dbOutputPath = joinpath("/Volumes/Yadgah/crispr-sweep-7-2-2022/isolates/emergence-lambert-root_output.sqlite") # local

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
# dbTempSim = SQLite.DB(dbSimPath) # local
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

include(joinpath(SCRIPT_PATH,"structures.jl"))
include(joinpath(SCRIPT_PATH,"computeExtinctionRoots.jl"))
include(joinpath(SCRIPT_PATH,"computeLambertRoot.jl"))

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
    dbOutput = matchStructure.dbOutput
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
        execute(dbOutput, "CREATE INDEX potential_vmatch_phenotypes_index
                            ON potential_vmatch_phenotypes (match_id)")
    end

    matchidlist = Vector{Int64}()
    escapelist = Vector{Int64}()
    for matchid in sort([keys(matchStructure.escapes)...])
        append!(escapelist,sort([matchStructure.escapes[matchid]...]))
        append!(matchidlist,
            fill(matchid,length(matchStructure.escapes[matchid]))
                )
    end
    DataFrame(match_id = matchidlist, escape_match_id = escapelist) |>
        SQLite.load!(dbOutput,
            "potential_vmatch_single_escapes",ifnotexists=true)
    execute(dbOutput, "CREATE INDEX potential_vmatch_single_escapes_index
                        ON potential_vmatch_single_escapes (match_id)")

    return matchIDs
end
# This functipn saves all possible subphenotypes (not just escapes) as a Dictionary (matchID0:[matchIDsubphenotypes])
# this is what I call a "match lineage"
function findEscapeLineages!(matchStructure::hierarchy,matchIDs::Dict{Vector{Int64},Int64})
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

function microbeBiomass!(matchStructure::hierarchy,current::state)
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

function computeR0!(extinction::pcomponents,parameters::params)
    time = extinction.time
    phi = parameters.adsorption_rate
    q = parameters.spacer_acquisition_prob
    beta = parameters.viral_burst_size
    d = parameters.viral_decay_rate
    sigdig = extinction.sigdig
    N = sum([abund for (abund,) in
            execute(extinction.dbSim, "SELECT abundance
                FROM babundance WHERE t = $(time)")])
    extinction.birth[0] = beta*phi*(1-q)*N
    extinction.mut[0] = 0
    extinction.death[0] = phi*q*N + d
    a = extinction.birth[0]
    m = extinction.death[0]
    if a > 0
        extinction.proots[0] =
            map(x->round(x,digits=sigdig),[1,m/a])
    else
        extinction.proots[0] = Vector{Float64}([1])
    end
    # println("proots of 0match is $(extinction.proots[0])")
    extinction.prootmin[0] =
        minimum(extinction.proots[0])
    # R0 = (phi*q*N+d)/(beta*phi*(1-q)*N)
    # extinction.R0 = 1/extinction.prootmin[0]
    # return Float64(extinction.R0)
end

# This function assembles all pertinent data from objects and saves into database
function assemble(current::state, extinction::pcomponents,escapes::Dict{Int64,Vector{Int64}})
    time = current.time
    dbTempSim = current.dbSim
    dbTempTri = current.dbTri
    dbOutput = current.dbOutput
    sigdig = extinction.sigdig

    pextinctionRoot = map(x -> extinction.prootmin[x], current.matches)
    pextinctionLamb = map(x -> extinction.plambert[x], current.matches)
    pextinctionAct = map(x -> extinction.pactual[x], current.matches)
    lysis = map(x -> extinction.lysis[x], current.matches)
    birth = map(x -> extinction.birth[x], current.matches)
    death = map(x -> extinction.death[x], current.matches)
    mut = map(x -> extinction.mut[x], current.matches)

    # matchAbunds = [Int64(abund) for (abund,)
    #                in
    #                execute(
    #     dbTempTri,
    #     "SELECT vabundance
    #     FROM vmatches_abundances
    #     WHERE t = $(time) AND match_id in
    #     ($(join(current.matches,", "))) ORDER BY match_id"
    #         )]

    #         if in(0, current.matches)
    #             matchAbunds = vcat(sum([Int64(abund) for (abund,)
    #                                     in
    #                                     execute(
    #                     dbTempTri,
    #                     "SELECT vabundance
    #     FROM v0matches_abundances
    #     WHERE t = $(time)"
    #                 )]), matchAbunds)
    #         end

    #         V = sum([Int64(abund) for (abund,)
    #                 in
    #                 execute(
    #             dbTempSim,
    #             "SELECT abundance
    #     FROM vabundance
    #     WHERE t = $(time)"
    #         )])

    # f = 1 / V

    DataFrame(t=fill(Float64(time), length(current.matches)),
        vmatch_id=current.matches,
        p_extinction_lambert=pextinctionLamb,
        p_extinction_actual=pextinctionAct,
        lysis=lysis, death=death, mutation=mut,
        p_extinction_root=pextinctionRoot, birth=birth#,
        # plambert_ext_weighted=matchAbunds .* pextinctionLamb * f,
        # plambert_emerge_weighted=matchAbunds .* (ones(length(pextinctionLamb)) - pextinctionLamb) * f,
        # pactual_ext_weighted=matchAbunds .* pextinctionAct * f,
        # pactual_emerge_weighted=matchAbunds .* (ones(length(pextinctionAct)) - pextinctionAct) * f,
        # vfrequency=matchAbunds * f,
        # vabundance=matchAbunds
    ) |> SQLite.load!(dbOutput, "existing_vmatch_extinction", ifnotexists=true)

    escapelist = Vector{Int64}()
    for matchid in current.matches
        if matchid == 0
            continue
        end
        append!(escapelist, sort([escapes[matchid]...]))
    end
    if length(escapelist) > 0
        unique!(escapelist)
        pextinctionRoot = map(x -> extinction.prootmin[x], escapelist)
        pextinctionLamb = map(x -> extinction.plambert[x], escapelist)
        pextinctionAct = map(x -> extinction.pactual[x], escapelist)
        lysis = map(x -> extinction.lysis[x], escapelist)
        birth = map(x -> extinction.birth[x], escapelist)
        death = map(x -> extinction.death[x], escapelist)
        mut = map(x -> extinction.mut[x], escapelist)

        DataFrame(t=fill(Float64(time), length(escapelist)),
            escape_vmatch_id=escapelist,
            p_extinction_lambert=pextinctionLamb,
            p_extinction_actual=pextinctionAct,
            lysis=lysis, death=death, mutation=mut,
            p_extinction_root=pextinctionRoot, birth=birth
        ) |> SQLite.load!(dbOutput, "single_escapes_vmatch_extinction", ifnotexists=true)
    end

    # DataFrame(  t = fill(Float64(time),length(current.escapes)),
    #             vmatch_id = current.escapes, p_extinction = pextinctionEsc,
    #             birth = birthEsc, death = deathEsc, mutation = mutEsc
    #             ) |> SQLite.load!(dbOutput, "potential_vmatch_extinction",ifnotexists=true)

    # return sum(matchAbunds .* pextinctionRoot * f),
    # sum(matchAbunds .* pextinctionLamb * f),
    # sum(matchAbunds .* pextinctionAct * f)
end

function emergence!(matchStructure::hierarchy,
                    current::state,parameters::params)
    extinction = computeExtinctionRoots(matchStructure,current,parameters)
    computeLambertRoot!(matchStructure,current,parameters,extinction)
    computeActualRoot!(matchStructure,current,parameters,extinction)
    # rootExpectation, lambertExpectation, actualExpectation = assemble(current,extinction)
    assemble(current,extinction,matchStructure.escapes)
    # return rootExpectation, lambertExpectation, actualExpectation
end
# This function runs the whole analysis
function escape()
    # This initializes the hierarchy of matchtypes that can emerge from escape
    matchStructure = hierarchy(dbTempTri,dbTempMatch,dbTempSim,dbOutput)
    # This calls and save the birth/death/mutation parameter rates of the models
    parameters = params(matchStructure.dbSim)
    matchIDs = potentialStructure!(matchStructure)
    findEscapeLineages!(matchStructure,matchIDs)
    # pActualExpectation = Vector{Float64}()
    # pLambExpectation = Vector{Float64}() 
    # pRootExpectation = Vector{Float64}()
    # times = Vector{Float64}()
    for (t,) in execute(dbTempSim,"SELECT DISTINCT t FROM vabundance")
        # if t == 150
        #     return
        # end
        # push!(times,t)
        println("Computing structure and probabilities at time = $(t)")
        current = state(matchStructure,Float64(t))
        microbeBiomass!(matchStructure,current)
        # pinitial, intExpectation, lambertExpectation, rootExpectation =
        #         emergence(matchStructure,current,parameters)
        # rootExpectation, lambertExpectation, actualExpectation =
        #         emergence(matchStructure,current,parameters)
        emergence!(matchStructure,current,parameters)
        # push!(pLambExpectation, lambertExpectation)
        # push!(pRootExpectation, rootExpectation)
        # push!(pActualExpectation, actualExpectation)
    end
    # DataFrame(  t = times,
    #             p_emerge_root_expected = ones(length(times))-pRootExpectation,
    #             p_emerge_lambert_expected = ones(length(times))-pLambExpectation,
    #             p_emerge_actual_expected = ones(length(times))-pActualExpectation
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
        # if in(table_name,["vmatch_expectation"])
        #     execute(db, "CREATE INDEX $(table_name)_index ON $(table_name) (t)")
        # end
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
