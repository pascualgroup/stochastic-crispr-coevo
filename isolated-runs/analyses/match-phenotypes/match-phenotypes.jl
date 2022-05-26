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
dbOutputPath = joinpath("match-phenotypes_output.sqlite") # cluster

# dbSimPath = joinpath("/Volumes/Yadgah","crispr-sweep-7-2-2022/isolates/runID3297-c66-r47/runID3297-c66-r47.sqlite") # local
# dbMatchPath = joinpath("/Volumes/Yadgah","crispr-sweep-7-2-2022/isolates/runID3297-c66-r47/matches_output.sqlite")
# # dbOutputPath = joinpath("/Volumes/Yadgah/crispr-sweep-7-2-2022/isolates/runID3297-c66-r47/match-phenotypes_output.sqlite") # local

if isfile(dbOutputPath)
    error("match-phenotypes_output.sqlite already exists; delete first")
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

execute(dbOutput, "CREATE TABLE bmatches
    (t REAL, match_id INTEGER, bstrain_id INTEGER)")
execute(dbOutput, "CREATE TABLE vmatches
    (t REAL, match_id INTEGER, vstrain_id INTEGER)")
execute(dbOutput, "CREATE TABLE bmatch_phenotypes
    (match_id INTEGER, phenotype INTEGER)")
execute(dbOutput, "CREATE TABLE vmatch_phenotypes
    (match_id INTEGER, phenotype INTEGER)")
execute(dbOutput, "CREATE TABLE bmatch_phenotypes_singles
    (t REAL, match_id INTEGER, spacer_id INTEGER)")
execute(dbOutput, "CREATE TABLE vmatch_phenotypes_singles
    (t REAL, match_id INTEGER, spacer_id INTEGER)")
execute(dbOutput, "CREATE TABLE tripartite_network
    (t REAL, vmatch_id INTEGER, bmatch_id INTEGER, spacer_id INTEGER)")
execute(dbOutput, "CREATE TABLE infection_network
    (t REAL, vmatch_id INTEGER, bstrain_id INTEGER)")
execute(dbOutput, "CREATE TABLE bmatches_abundances
    (t REAL, bmatch_id INTEGER,
    babundance INTEGER, bfrequency REAL,
    s_vabundance INTEGER, s_vfrequency REAL,
    i_vabundance INTEGER, i_vfrequency REAL)")
execute(dbOutput, "CREATE TABLE vmatches_abundances
    (t REAL, vmatch_id INTEGER,
    vabundance INTEGER, vfrequency REAL,
    s_babundance INTEGER, s_bfrequency REAL,
    i_babundance INTEGER, i_bfrequency REAL)")

dbTempSim = SQLite.DB()
execute(dbTempSim, "CREATE TABLE bspacers (bstrain_id INTEGER, spacer_id INTEGER)")

execute(dbTempSim, "BEGIN TRANSACTION")
execute(dbTempSim,"ATTACH DATABASE '$(dbSimPath)' as dbSim")
execute(dbTempSim,"INSERT INTO bspacers (bstrain_id,spacer_id)
    SELECT bstrain_id,spacer_id FROM dbSim.bspacers WHERE run_id = $(run_id);")
execute(dbTempSim, "COMMIT")


execute(dbTempSim, "BEGIN TRANSACTION")
execute(dbTempSim, "CREATE INDEX bspacers_index ON bspacers (bstrain_id,spacer_id)")
execute(dbTempSim, "COMMIT")



mutable struct tripartite
    dbMatch::DB
    dbSim::DB
    dbOutput::DB
    time::Float64
    vphenoID::Int64
    bphenoID::Int64

    vmatchtypes::Vector{Int64}
    vphenotypes::Vector{Vector{Int64}}
    vsinglematches::Dict{Vector{Int64},Vector{Int64}}
    vmatchIDs::Dict{Vector{Int64},Int64}
    vstrainclasses::Dict{Vector{Int64},Vector{Int64}}
    sBstrains::Dict{Vector{Int64},Vector{Int64}}
    iBstrains::Dict{Vector{Int64},Vector{Int64}}
    vphenoBiomass::Dict{Vector{Int64},Int64}
    vphenoFrequency::Dict{Vector{Int64},Float64}
    sBBiomass::Dict{Vector{Int64},Int64}
    iBBiomass::Dict{Vector{Int64},Int64}
    sBfrequency::Dict{Vector{Int64},Float64}()
    iBfrequency::Dict{Vector{Int64},Float64}()

    bmatchtypes::Vector{Int64}
    bphenotypes::Vector{Vector{Int64}}
    bsinglematches::Dict{Vector{Int64},Vector{Int64}}
    bmatchIDs::Dict{Vector{Int64},Int64}
    bstrainclasses::Dict{Vector{Int64},Vector{Int64}}
    sVstrains::Dict{Vector{Int64},Vector{Int64}}
    iVstrains::Dict{Vector{Int64},Vector{Int64}}
    bphenoBiomass::Dict{Vector{Int64},Int64}
    bphenoFrequency::Dict{Vector{Int64},Float64}
    sVBiomass::Dict{Vector{Int64},Int64}
    iVBiomass::Dict{Vector{Int64},Int64}
    sVfrequency::Dict{Vector{Int64},Float64}()
    iVfrequency::Dict{Vector{Int64},Float64}()
    function tripartite(dbTempMatch::DB,dbTempSim::DB,dbOutput::DB,t::Float64)
        new(
            dbTempMatch,
            dbTempSim,
            dbOutput,
            t,
            Int64(1),
            Int64(1),
            Vector{Int64}(),
            Vector{Vector{Int64}}(),
            Dict{Vector{Int64},Vector{Int64}}(),
            Dict{Vector{Int64},Int64}(),
            Dict{Vector{Int64},Vector{Int64}}(),
            Dict{Vector{Int64},Vector{Int64}}(),
            Dict{Vector{Int64},Vector{Int64}}(),
            Dict{Vector{Int64},Int64}(),
            Dict{Vector{Int64},Float64}(),
            Dict{Vector{Int64},Int64}(),
            Dict{Vector{Int64},Int64}(),
            Dict{Vector{Int64},Float64}(),
            Dict{Vector{Int64},Float64}(),
            Vector{Int64}(),
            Vector{Vector{Int64}}(),
            Dict{Vector{Int64},Vector{Int64}}(),
            Dict{Vector{Int64},Int64}(),
            Dict{Vector{Int64},Vector{Int64}}(),
            Dict{Vector{Int64},Vector{Int64}}(),
            Dict{Vector{Int64},Vector{Int64}}(),
            Dict{Vector{Int64},Int64}(),
            Dict{Vector{Int64},Float64}(),
            Dict{Vector{Int64},Int64}(),
            Dict{Vector{Int64},Int64}(),
            Dict{Vector{Int64},Float64}(),
            Dict{Vector{Int64},Float64}()
            )
    end
end


function structure!(matchStructure::tripartite)
    clearCurrentStructure!(matchStructure)
    currentStructure!(matchStructure)
    phenoBiomass!(matchStructure)
    return matchStructure
end

function currentStructure!(matchStructure::tripartite)
    time = matchStructure.time
    dbTempSim = matchStructure.dbSim
    dbTempMatch = matchStructure.dbMatch
    for (vstrain_id,) in execute(dbTempMatch,"SELECT DISTINCT vstrain_id
            FROM vstrain_matched_pspacers WHERE t = $(time)
            ORDER BY vstrain_id")
        if length([strain for (strain,) in
                execute(dbTempMatch,"SELECT bstrain_id
                    FROM bstrain_to_vstrain_0matches
                    WHERE t = $(time) AND vstrain_id = $(vstrain_id)")]) == 0
            continue
        end
        phenotype = sort([type for (type,) in
        execute(dbTempMatch,"SELECT matched_pspacer_id
        FROM vstrain_matched_pspacers
        WHERE t = $(time) AND vstrain_id = $(vstrain_id)
        ORDER BY matched_pspacer_id")])
        # println("phenotype is $(phenotype)")
        if in(phenotype,matchStructure.vphenotypes)
            push!(matchStructure.vstrainclasses[phenotype],Int64(vstrain_id))
        else
            push!(matchStructure.vphenotypes,phenotype)
            matchStructure.vstrainclasses[phenotype] = [Int64(vstrain_id)]
            identifyCurrentMatches!(matchStructure,phenotype,vstrain_id,true)
            matchStructure.vphenotypes[phenotype] = matchStructure.vphenoID
            union!(matchStructure.vmatchtypes,length(phenotype))
            matchStructure.vphenoID += Int64(1)
        end

    end

    for (bstrain_id,) in execute(dbTempMatch,"SELECT DISTINCT bstrain_id
            FROM bstrain_matched_pspacers WHERE t = $(time)
            ORDER BY bstrain_id")
        if length([strain for (strain,) in
                execute(dbTempMatch,"SELECT vstrain_id
                    FROM bstrain_to_vstrain_matches
                    WHERE t = $(time) AND match_length = 1
                    AND bstrain_id = $(bstrain_id)")]) == 0
            continue
        end
        phenotype = sort([type for (type,) in
        execute(dbTempMatch,"SELECT matched_pspacer_id
        FROM bstrain_matched_spacers
        WHERE t = $(time) AND bstrain_id = $(bstrain_id)
        ORDER BY matched_pspacer_id")])
        # println("phenotype is $(phenotype)")
        if in(phenotype,matchStructure.bphenotypes)
            push!(matchStructure.bstrainclasses[phenotype],Int64(bstrain_id))
        else
            push!(matchStructure.bphenotypes,phenotype)
            matchStructure.bstrainclasses[phenotype] = [Int64(bstrain_id)]
            identifyCurrentMatches!(matchStructure,phenotype,bstrain_id,false)
            matchStructure.bphenotypes[phenotype] = matchStructure.bphenoID
            union!(matchStructure.bmatchtypes,length(phenotype))
            matchStructure.bphenoID += Int64(1)
        end
    end
end

function clearCurrentStructure!(matchStructure::tripartite)
    matchStructure.vstrainclasses = Dict{Vector{Int64},Vector{Int64}}()
    matchStructure.vsinglematches = Dict{Vector{Int64},Vector{Int64}}()
    matchStructure.sBstrains =  Dict{Vector{Int64},Vector{Int64}}()
    matchStructure.iBstrains = Dict{Vector{Int64},Vector{Int64}}()
    matchStructure.vphenoBiomass = Dict{Vector{Int64},Int64}()
    matchStructure.vphenoFrequency = Dict{Vector{Int64},Float64}()
    matchStructure.sBBiomass = Dict{Vector{Int64},Int64}()
    matchStructure.iBBiomass = Dict{Vector{Int64},Int64}()
    matchStructure.sBfrequency = Dict{Vector{Int64},Float64}()
    matchStructure.iBfrequency = Dict{Vector{Int64},Float64}()

    matchStructure.bstrainclasses = Dict{Vector{Int64},Vector{Int64}}()
    matchStructure.bsinglematches = Dict{Vector{Int64},Vector{Int64}}()
    matchStructure.sVstrains = Dict{Vector{Int64},Vector{Int64}}()
    matchStructure.iVstrains = Dict{Vector{Int64},Vector{Int64}}()
    matchStructure.bphenoBiomass = Dict{Vector{Int64},Int64}()
    matchStructure.bphenoFrequency = Dict{Vector{Int64},Float64}()
    matchStructure.sVBiomass = Dict{Vector{Int64},Int64}()
    matchStructure.iVBiomass = Dict{Vector{Int64},Int64}()
    matchStructure.sVfrequency = Dict{Vector{Int64},Float64}()
    matchStructure.iVfrequency = Dict{Vector{Int64},Float64}()
end

function identifyCurrentMatches!(
    matchStructure::tripartite,phenotype::Vector{Int64},strain_id::Int64,vb::Bool)
    time = matchStructure.time
    dbTempSim = matchStructure.dbSim
    dbTempMatch = matchStructure.dbMatch
    if vb
        Btotal = sum([abund for (abund,) in
                        execute(dbTempSim, "SELECT abundance
                            FROM babundance WHERE t = $(time)")])
        matchStructure.sBstrains[phenotype] = [Int64(strain)
            for (strain,) in execute(dbTempMatch, "SELECT bstrain_id
                FROM bstrain_to_vstrain_0matches
                WHERE t = $(time) AND vstrain_id = $(strain_id)")]
        matchStructure.iBstrains[phenotype] = [Int64(strain)
            for (strain,) in execute(dbTempMatch, "SELECT bstrain_id
                FROM bstrain_to_vstrain_matches
                WHERE t = $(time)
                AND match_length = 1
                AND vstrain_id = $(strain_id)")]
        timeIDs = [matchID
            for (matchID,) in execute(dbTempMatch, "SELECT time_specific_match_id
                FROM bstrain_to_vstrain_matches
                WHERE t = $(time)
                AND match_length = 1
                AND vstrain_id = $(strain_id)")]
        matchStructure.vsinglematches[phenotype] = [spacerID
            for (spacerID,) in execute(dbTempMatch, "SELECT spacer_id
                FROM matches_spacers
                WHERE t = $(time)
                AND time_specific_match_id
                in ($(join(timeIDs,", ")))
                ORDER BY spacer_id")]

        # allBstrains = [Int64(strain)
        #     for (strain,) in execute(dbTempSim, "SELECT bstrain_id
        #         FROM babundance
        #         WHERE t = $(time)") if strain != 1]
        # matchStructure.iBstrains[phenotype] =
        #     setdiff(allBstrains,matchStructure.sBstrains[phenotype])


        matchStructure.sBBiomass[phenotype] = sum([Int64(abund)
            for (abund,) in execute(dbTempSim, "SELECT abundance
                FROM babundance
                WHERE t = $(time) AND bstrain_id in
                ($(join(matchStructure.sBstrains[phenotype],", ")))")])
        matchStructure.sBfrequency[phenotype] =
            matchStructure.sBBiomass[phenotype]/Btotal

        if length(matchStructure.iBstrains[phenotype]) == 0
            matchStructure.iBBiomass[phenotype] = 0
            matchStructure.iBfrequency[phenotype] = 0
        else
            matchStructure.iBBiomass[phenotype] = sum([Int64(abund)
                for (abund,) in execute(dbTempSim, "SELECT abundance
                    FROM babundance
                    WHERE t = $(time) AND bstrain_id in
                    ($(join(matchStructure.iBstrains[phenotype],", ")))")])
            matchStructure.iBfrequency[phenotype] =
                matchStructure.iBBiomass[phenotype]/Btotal
        end
    else
        Vtotal = sum([abund for (abund,) in
                        execute(dbTempSim, "SELECT abundance
                            FROM vabundance WHERE t = $(time)")])
        matchStructure.sVstrains[phenotype] = [Int64(strain)
            for (strain,) in execute(dbTempMatch, "SELECT vstrain_id
                FROM bstrain_to_vstrain_0matches
                WHERE t = $(time) AND bstrain_id = $(strain_id)")]
        matchStructure.iVstrains[phenotype] = [Int64(strain)
            for (strain,) in execute(dbTempMatch, "SELECT vstrain_id
                FROM bstrain_to_vstrain_matches
                WHERE t = $(time)
                AND match_length = 1
                AND bstrain_id = $(strain_id)")]
        timeIDs = [matchID
            for (matchID,) in execute(dbTempMatch, "SELECT time_specific_match_id
                FROM bstrain_to_vstrain_matches
                WHERE t = $(time)
                AND match_length = 1
                AND bstrain_id = $(strain_id)")]
        matchStructure.bsinglematches[phenotype] = [spacerID
            for (spacerID,) in execute(dbTempMatch, "SELECT spacer_id
                FROM matches_spacers
                WHERE t = $(time)
                AND time_specific_match_id
                in ($(join(timeIDs,", ")))
                ORDER BY spacer_id")]

        # allVstrains = [Int64(strain)
        #     for (strain,) in execute(dbTempSim, "SELECT vstrain_id
        #         FROM vabundance
        #         WHERE t = $(time)")]
        # matchStructure.iVstrains[phenotype] =
        #     setdiff(allVstrains,matchStructure.sVstrains[phenotype])

        if length(matchStructure.sVstrains[phenotype]) == 0
            matchStructure.sVBiomass[phenotype] = 0
            matchStructure.sVfrequency[phenotype] = 0
        else
            matchStructure.sVBiomass[phenotype] = sum([Int64(abund)
                for (abund,) in execute(dbTempSim, "SELECT abundance
                    FROM vabundance
                    WHERE t = $(time) AND vstrain_id in
                    ($(join(matchStructure.sVstrains[phenotype],", ")))")])
            matchStructure.sVfrequency[phenotype] =
                matchStructure.sVBiomass[phenotype]/Vtotal
        end

        if length(matchStructure.iVstrains[phenotype]) == 0
            matchStructure.iVBiomass[phenotype] = 0
            matchStructure.iVfrequency[phenotype] = 0
        else
            matchStructure.iVBiomass[phenotype] = sum([Int64(abund)
                for (abund,) in execute(dbTempSim, "SELECT abundance
                    FROM vabundance
                    WHERE t = $(time) AND vstrain_id in
                    ($(join(matchStructure.iVstrains[phenotype],", ")))")])
            matchStructure.iVfrequency[phenotype] =
                matchStructure.iVBiomass[phenotype]/Vtotal
        end
    end
end


function phenoBiomass!(matchStructure::tripartite)
    dbTempSim = matchStructure.dbSim
    time = matchStructure.time
    Vtotal = sum([abund for (abund,) in
                    execute(dbTempSim, "SELECT abundance
                        FROM vabundance WHERE t = $(time)")])
    for phenotype in keys(matchStructure.vstrainclasses)
        V = sum([Int64(abund) for (abund,) in execute(dbTempSim, "SELECT abundance
            FROM vabundance
            WHERE t = $(time) AND vstrain_id in
            ($(join(matchStructure.vstrainclasses[phenotype],", ")))")])
        matchStructure.vphenoBiomass[phenotype] = V
        matchStructure.vphenoFrequency[phenotype] = V/Vtotal
    end
    allVstrains = [strain for (strain,) in execute(dbTempSim,"SELECT vstrain_id
                    FROM vabundance WHERE t = $(matchStructure.time)")]
    match0 = setdiff(allVstrains,vcat(values(matchStructure.vstrainclasses)...))
    if length(match0) == 0
        matchStructure.vphenoBiomass[[0]] = 0
        matchStructure.vphenoFrequency[[0]] = 0
    else
        V = sum([Int64(abund) for (abund,) in execute(dbTempSim, "SELECT abundance
            FROM vabundance
            WHERE t = $(time) AND vstrain_id in
            ($(join(match0,", ")))")])
        matchStructure.vphenoBiomass[[0]] = V
        matchStructure.vphenoFrequency[[0]] = V/Vtotal
    end


    Btotal = sum([abund for (abund,) in
                    execute(dbTempSim, "SELECT abundance
                        FROM babundance WHERE t = $(time)")])
    for phenotype in keys(matchStructure.bstrainclasses)
        B = sum([Int64(abund) for (abund,) in execute(dbTempSim, "SELECT abundance
            FROM babundance
            WHERE t = $(time) AND bstrain_id in
            ($(join(matchStructure.bstrainclasses[phenotype],", ")))")])
        matchStructure.bphenoBiomass[phenotype] = B
        matchStructure.bphenoFrequency[phenotype] = B/Btotal
    end
    allBstrains = [strain for (strain,) in execute(dbTempSim,"SELECT bstrain_id
                    FROM babundance WHERE t = $(matchStructure.time)")]
    match0 = setdiff(allBstrains,vcat(values(matchStructure.bstrainclasses)...))
    if length(match0) == 0
        matchStructure.bphenoBiomass[[0]] = 0
        matchStructure.bphenoFrequency[[0]] = 0
    else
        B = sum([Int64(abund) for (abund,) in execute(dbTempSim, "SELECT abundance
            FROM babundance
            WHERE t = $(time) AND bstrain_id in
            ($(join(match0,", ")))")])
        matchStructure.bphenoBiomass[[0]] = B
        matchStructure.bphenoFrequency[[0]] = B/Btotal
    end
end

function updateMatchesDB(matchStructure::tripartite)
    time = matchStructure.time
    dbOutput = matchStructure.dbOutput

    for phenotype in keys(matchStructure.vstrainclasses)
        matchID = matchStructure.vmatchIDs[phenotype]
        for spacerID in phenotype
            execute(dbOutput, "INSERT INTO vmatch_phenotypes VALUES (?,?)",
                (matchID,spacerID))
        end
        for vstrainID in matchStructure.vstrainclasses[phenotype]
            execute(dbOutput, "INSERT INTO vmatches VALUES (?,?)",
                (time, matchID,vstrainID))
        end
        execute(dbOutput, "INSERT INTO vmatches_abundances VALUES
            (?,?,?,?,?,?)",
            (time, matchID,
            matchStructure.vphenoBiomass[phenotype],
            matchStructure.vphenoFrequency[phenotype],
            matchStructure.sBBiomass[phenotype],
            matchStructure.sBfrequency[phenotype],
            matchStructure.iBBiomass[phenotype],
            matchStructure.iBfrequency[phenotype]))
    end

    for phenotype in keys(matchStructure.bstrainclasses)
        matchID = matchStructure.bmatchIDs[phenotype]
        for spacerID in phenotype
            execute(dbOutput, "INSERT INTO bmatch_phenotypes VALUES (?,?)",
                (matchID,spacerID))
        end
        for bstrainID in matchStructure.bstrainclasses[phenotype]
            execute(dbOutput, "INSERT INTO bmatches VALUES (?,?)",
                (time, matchID,bstrainID))
        end
        execute(dbOutput, "INSERT INTO bmatches_abundances VALUES
            (?,?,?,?,?,?)",
            (time, matchID,
            matchStructure.bphenoBiomass[phenotype],
            matchStructure.bphenoFrequency[phenotype],
            matchStructure.sVBiomass[phenotype],
            matchStructure.sVfrequency[phenotype],
            matchStructure.iVBiomass[phenotype],
            matchStructure.iVfrequency[phenotype]))
    end



end

function logPhenotypesDB(matchStructure::tripartite)
    dbOutput = matchStructure.dbOutput

    for phenotype in matchStructure.vphenotypes
        matchID = matchStructure.vmatchIDs[phenotype]
        for spacerID in phenotype
            execute(dbOutput, "INSERT INTO vmatch_phenotypes VALUES (?,?)",
                (Int64(matchID),Int64(spacerID)))
        end
    end
    for phenotype in matchStructure.bphenotypes
        matchID = matchStructure.bmatchIDs[phenotype]
        for spacerID in phenotype
            execute(dbOutput, "INSERT INTO bmatch_phenotypes VALUES (?,?)",
                (Int64(matchID),Int64(spacerID)))
        end
    end
end

function logTripartiteDB(matchStructure::tripartite)
    time = matchStructure.time
    dbOutput = matchStructure.dbOutput

    for phenotype in keys(matchStructure.vstrainclasses)
        if length(matchStructure.vsinglematches[phenotype]) == 0
            continue
        end
        vmatchID = matchStructure.vmatchIDs[phenotype]
        for spacerID in matchStructure.vsinglematches[phenotype]
            execute(dbOutput, "INSERT INTO vmatch_phenotypes_singles VALUES (?,?,?)",
                (time,Int64(vmatchID),spacerID))
        end
    end

    for phenotype in keys(matchStructure.bstrainclasses)
        if length(matchStructure.bsinglematches[phenotype]) == 0
            continue
        end
        bmatchID = matchStructure.bmatchIDs[phenotype]
        for spacerID in matchStructure.bsinglematches[phenotype]
            execute(dbOutput, "INSERT INTO bmatch_phenotypes_singles VALUES (?,?,?)",
                (time,Int64(bmatchID),spacerID))
        end
    end

    vspacers = [spacerID for (spacerID,) in
                    execute(dbOutput, "SELECT DISTINCT spacer_id
                    FROM vmatch_phenotypes_singles
                    WHERE t = $(time)
                    ORDER BY spacer_id")]
    bspacers = [spacerID for (spacerID,) in
                    execute(dbOutput, "SELECT DISTINCT spacer_id
                    FROM bmatch_phenotypes_singles
                    WHERE t = $(time)
                    ORDER BY spacer_id")]

    @assert vspacers == bspacers

    for spacerID in execute(dbOutput, "SELECT DISTINCT spacer_id
                        FROM vmatch_phenotypes_singles
                        ORDER BY spacer_id")
        for (vmatchID,) in execute(dbOutput, "SELECT DISTINCT vmatch_id
                        FROM vmatch_phenotypes_singles
                        WHERE t = $(time)
                        AND spacer_id = $(spacerID)
                        ORDER BY vmatch_id")
            for (bmatchID,) in execute(dbOutput, "SELECT DISTINCT bmatch_id
                            FROM bmatch_phenotypes_singles
                            WHERE t = $(time)
                            AND spacer_id = $(spacerID)
                            ORDER BY bmatch_id")
                execute(dbOutput, "INSERT INTO tripartite_networks VALUES (?,?,?,?)",
                    (time,Int64(vmatchID),Int64(bmatchID),spacerID))
            end
        end
    end
end



matchStructure = tripartite(dbTempMatch,dbTempSim,dbOutput,Float64(0))
for (t,) in execute(dbTempSim,"SELECT DISTINCT t FROM vabundance")
    if t == 0
         continue
    end
    println("Computing tripartite structure at time = $(t)")
    matchStructure.time = t
    structure!(matchStructure)
    # updateMatchesDB!(matchStructure)
    logTripartiteDB(matchStructure)
end
logPhenotypesDB(matchStructure)



function createindices()
    println("(Creating run_id indices...)")
    db = SQLite.DB(dbOutputPath)
    execute(db, "BEGIN TRANSACTION")
    for (table_name,) in execute(
        db, "SELECT name FROM sqlite_schema
        WHERE type='table' ORDER BY name;")
        execute(db, "CREATE INDEX $(table_name)_index ON $(table_name) (t)")
    end
    execute(db, "COMMIT")
end
createindices()




println("Complete!")
