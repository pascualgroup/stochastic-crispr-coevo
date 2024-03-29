#!/usr/bin/env julia

### The full match phenotypes for the viruses needed to be classified so that
### we can use the data on probability of emergence
# Full match phenotypes of micobes are also classified, however with the
# condition that they have at least one single match with a viral strain.
##
println("(Julia compilation delay...)")

using SQLite
using DataFrames
using DataFramesMeta
using SQLite.DBInterface: execute
using SQLite: DB

run_id = ARGS[1]

## Define Paths ##
SCRIPT_PATH = abspath(dirname(PROGRAM_FILE))
#
dbSimPath = joinpath(SCRIPT_PATH, "..", "..", "..", "simulation", "sweep_db_gathered.sqlite") # cluster
dbOutputPath = joinpath("tripartite-networks_output.sqlite") # cluster

# dbSimPath = joinpath("/Volumes/Yadgah","crispr-sweep-7-2-2022/isolates/runID3297-c66-r47/runID3297-c66-r47.sqlite") # local
# dbMatchPath = joinpath("/Volumes/Yadgah","crispr-sweep-7-2-2022/isolates/runID3297-c66-r47/matches_output.sqlite") # local
# dbOutputPath = joinpath("/Volumes/Yadgah/crispr-sweep-7-2-2022/isolates/runID3297-c66-r47/tripartite-networks_output.sqlite") # local
# dbSimPath = joinpath("/Users/armun/Dropbox/Current/Projects/microbe-virus-crispr/stochastic-crispr/repository-2/isolated-runs/isolates/runID3297-c66-r47/runID3297-c66-r47.sqlite") # local
# dbMatchPath = joinpath("/Users/armun/Dropbox/Current/Projects/microbe-virus-crispr/stochastic-crispr/repository-2/isolated-runs/isolates/runID3297-c66-r47/matches_output.sqlite") # local
# dbOutputPath = joinpath("/Users/armun/Dropbox/Current/Projects/microbe-virus-crispr/stochastic-crispr/repository-2/tripartite-networks_output.sqlite") # local


if isfile(dbOutputPath)
    error("tripartite-networks_output.sqlite already exists; delete first")
end
##
dbTempSim = SQLite.DB(dbSimPath)
(cr,) = execute(dbTempSim, "SELECT combo_id,replicate FROM runs WHERE run_id = $(run_id)")
if length(ARGS) > 1 && ARGS[end] == "combo"
    if isfile(joinpath(SCRIPT_PATH, "..", "..", "..", "data-analysis",
        "gathered-analyses", "matches", "matches.sqlite"))
        dbMatchPath = joinpath(SCRIPT_PATH, "..", "..", "..", "data-analysis",
            "gathered-analyses", "matches", "matches.sqlite")
    elseif isfile(joinpath(SCRIPT_PATH, "..", "..", "isolates",
        "comboID$(cr[1])", "matchesC$(cr[1]).sqlite"))
        dbMatchPath = joinpath(SCRIPT_PATH, "..", "..", "isolates",
            "comboID$(cr[1])", "matchesC$(cr[1]).sqlite")
    else
        error("neither /data-analysis/gathered-analyses/matches/matches.sqlite exists \n
            nor does /isolated-runs/isolates/comboID$(cr[1])/matchesC$(cr[1]).sqlite; \n
            compute one of them first")
    end
end
if length(ARGS) == 1
    dbMatchPath = joinpath(SCRIPT_PATH, "..", "..", "isolates",
        "runID$(run_id)-c$(cr.combo_id)-r$(cr.replicate)", "matches_output.sqlite")
    if !isfile(dbMatchPath)
        error("matches_output.sqlite does not exist; compute first")
    end

end
dbTempMatch = SQLite.DB(dbMatchPath)
dbOutput = SQLite.DB(dbOutputPath)
dbTempSim = SQLite.DB()
# dbTempSim = SQLite.DB(dbSimPath)
execute(dbTempSim, "CREATE TABLE babundance (t REAL, bstrain_id, abundance INTEGER)")
execute(dbTempSim, "CREATE TABLE vabundance (t REAL, vstrain_id, abundance INTEGER)")
execute(dbTempSim, "BEGIN TRANSACTION")
execute(dbTempSim, "ATTACH DATABASE '$(dbSimPath)' as dbSim")
execute(dbTempSim,"INSERT INTO babundance (t, bstrain_id, abundance) 
                    SELECT t, bstrain_id, abundance 
                    FROM dbSim.babundance WHERE run_id = $(run_id);")
execute(dbTempSim,"INSERT INTO vabundance (t, vstrain_id, abundance) 
                    SELECT t, vstrain_id, abundance 
                    FROM dbSim.vabundance WHERE run_id = $(run_id);")
execute(dbTempSim, "COMMIT")
execute(dbTempSim, "BEGIN TRANSACTION")
execute(dbTempSim, "CREATE INDEX babundance_index ON babundance (t,bstrain_id,abundance)")
execute(dbTempSim, "CREATE INDEX vabundance_index ON vabundance (t,vstrain_id,abundance)")
execute(dbTempSim, "COMMIT")

mutable struct tripartite
    dbMatch::DB
    dbSim::DB
    dbOutput::DB
    time::Float64
    strainID::Int64
    vphenoID::Int64
    bphenoID::Int64
    phenotype::Vector{Int64}
    vb::Bool

    vmatchtypes::Vector{Int64}
    vphenotypes::Vector{Vector{Int64}}
    vsinglematches::Dict{Vector{Int64},Vector{Int64}}
    vmatchIDs::Dict{Vector{Int64},Int64}
    vstrainclasses::Dict{Vector{Int64},Vector{Int64}}
    vphenoBiomass::Dict{Vector{Int64},Int64}
    sbstrains::Dict{Vector{Int64},Vector{Int64}}
    sBiomass::Dict{Vector{Int64},Int64}
    ibstrains::Dict{Vector{Int64},Vector{Int64}}
    iBiomass::Dict{Vector{Int64},Int64}

    bmatchtypes::Vector{Int64}
    bphenotypes::Vector{Vector{Int64}}
    bsinglematches::Dict{Vector{Int64},Vector{Int64}}
    bmatchIDs::Dict{Vector{Int64},Int64}
    bstrainclasses::Dict{Vector{Int64},Vector{Int64}}
    bphenoBiomass::Dict{Vector{Int64},Int64}

    function tripartite(dbTempMatch::DB, dbTempSim::DB, dbOutput::DB, t::Float64)
        new(
            dbTempMatch,
            dbTempSim,
            dbOutput,
            t,
            Int64(0),
            Int64(1),
            Int64(1),
            Vector{Int64}(),
            true,
            Vector{Int64}(),
            Vector{Vector{Int64}}(),
            Dict{Vector{Int64},Vector{Int64}}(),
            Dict{Vector{Int64},Int64}(),
            Dict{Vector{Int64},Vector{Int64}}(),
            Dict{Vector{Int64},Int64}(),
            Dict{Vector{Int64},Vector{Int64}}(),
            Dict{Vector{Int64},Int64}(),
            Dict{Vector{Int64},Vector{Int64}}(),
            Dict{Vector{Int64},Int64}(),
            Vector{Int64}(),
            Vector{Vector{Int64}}(),
            Dict{Vector{Int64},Vector{Int64}}(),
            Dict{Vector{Int64},Int64}(),
            Dict{Vector{Int64},Vector{Int64}}(),
            Dict{Vector{Int64},Int64}()
        )
    end
end

if length(ARGS) > 1 && ARGS[end] == "combo"
    include(joinpath(SCRIPT_PATH, "combo-functions.jl"))
else
    include(joinpath(SCRIPT_PATH, "functions.jl"))
end

function structure!(matchStructure::tripartite)
    clearCurrentStructure!(matchStructure)
    currentStructure!(matchStructure)
    phenoBiomass!(matchStructure)
    return matchStructure
end

function newPhenotype!(matchStructure::tripartite)
    vb = matchStructure.vb
    phenotype = matchStructure.phenotype
    if vb
        push!(matchStructure.vphenotypes, phenotype)
        matchStructure.vstrainclasses[phenotype] = [matchStructure.strainID]
        matchStructure.vmatchIDs[phenotype] = matchStructure.vphenoID
        union!(matchStructure.vmatchtypes, length(phenotype))
        matchStructure.vphenoID += Int64(1)
    else
        push!(matchStructure.bphenotypes, phenotype)
        matchStructure.bstrainclasses[phenotype] = [matchStructure.strainID]
        matchStructure.bmatchIDs[phenotype] = matchStructure.bphenoID
        union!(matchStructure.bmatchtypes, length(phenotype))
        matchStructure.bphenoID += Int64(1)
    end
    identifyCurrentMatches!(matchStructure)
end

function newCurrent!(matchStructure::tripartite)
    vb = matchStructure.vb
    phenotype = matchStructure.phenotype
    if vb
        matchStructure.vstrainclasses[phenotype] = [matchStructure.strainID]
    else
        matchStructure.bstrainclasses[phenotype] = [matchStructure.strainID]
    end
    identifyCurrentMatches!(matchStructure)
end

function clearCurrentStructure!(matchStructure::tripartite)
    matchStructure.vb = true
    matchStructure.vstrainclasses = Dict{Vector{Int64},Vector{Int64}}()
    matchStructure.vsinglematches = Dict{Vector{Int64},Vector{Int64}}()
    matchStructure.vphenoBiomass = Dict{Vector{Int64},Int64}()
    matchStructure.sbstrains = Dict{Vector{Int64},Vector{Int64}}()
    matchStructure.sBiomass = Dict{Vector{Int64},Int64}()
    matchStructure.ibstrains = Dict{Vector{Int64},Vector{Int64}}()
    matchStructure.iBiomass = Dict{Vector{Int64},Int64}()
    matchStructure.bstrainclasses = Dict{Vector{Int64},Vector{Int64}}()
    matchStructure.bsinglematches = Dict{Vector{Int64},Vector{Int64}}()
    matchStructure.bphenoBiomass = Dict{Vector{Int64},Int64}()
end

function phenoBiomass!(matchStructure::tripartite)
    dbTempSim = matchStructure.dbSim
    time = matchStructure.time
    for phenotype in keys(matchStructure.vstrainclasses)
        V = sum([Int64(abund) for (abund,) in execute(
            dbTempSim,
            "SELECT abundance
FROM vabundance
WHERE t = $(time) AND vstrain_id in
($(join(matchStructure.vstrainclasses[phenotype],", ")))"
        )])
        matchStructure.vphenoBiomass[phenotype] = V
    end
    for phenotype in keys(matchStructure.bstrainclasses)
        B = sum([Int64(abund) for (abund,) in execute(
            dbTempSim,
            "SELECT abundance
FROM babundance WHERE t = $(time) AND bstrain_id in 
($(join(matchStructure.bstrainclasses[phenotype],", ")))"
        )])
        matchStructure.bphenoBiomass[phenotype] = B
    end
end

function updateMatchesDB(matchStructure::tripartite)
    time = matchStructure.time
    dbOutput = matchStructure.dbOutput
    matchIDs = Vector{Int64}()
    strainIDs = Vector{Int64}()
    biomass = Vector{Int64}()
    smatchIDsB = Vector{Int64}()
    sbstrainIDs = Vector{Int64}()
    imatchIDsB = Vector{Int64}()
    ibstrainIDs = Vector{Int64}()
    sBiomass = Vector{Int64}()
    iBiomass = Vector{Int64}()
    for phenotype in keys(matchStructure.vstrainclasses)
        matchID = matchStructure.vmatchIDs[phenotype]
        strains = matchStructure.vstrainclasses[phenotype]
        matchIDs = vcat(matchIDs,
            repeat([matchID], length(strains)))
        strainIDs = vcat(strainIDs, strains)

        biomass = vcat(biomass, matchStructure.vphenoBiomass[phenotype])

        sbstrains = matchStructure.sbstrains[phenotype]
        smatchIDsB = vcat(smatchIDsB,
            repeat([matchID], length(sbstrains)))
        sbstrainIDs = vcat(sbstrainIDs, sbstrains)
        sBiomass = vcat(sBiomass, matchStructure.sBiomass[phenotype])

        ibstrains = matchStructure.ibstrains[phenotype]
        imatchIDsB = vcat(imatchIDsB,
            repeat([matchID], length(ibstrains)))
        ibstrainIDs = vcat(ibstrainIDs, ibstrains)
        iBiomass = vcat(iBiomass, matchStructure.iBiomass[phenotype])
    end

    matches = DataFrame(t=repeat([Float64(time)], length(matchIDs)),
        match_id=matchIDs, vstrain_id=strainIDs)
    matches |> SQLite.load!(dbOutput,
        "vmatches", ifnotexists=true)
    unique!(matchIDs)
    matches = DataFrame(t=repeat([Float64(time)], length(matchIDs)),
        match_id=matchIDs, vabundance=biomass,
        bsusceptible=sBiomass, bimmune=iBiomass)
    matches |> SQLite.load!(dbOutput,
        "vmatches_abundances", ifnotexists=true)
    matches = DataFrame(t=repeat([Float64(time)], length(smatchIDsB)),
        vmatch_id=smatchIDsB, bstrain_id=sbstrainIDs)
    matches |> SQLite.load!(dbOutput,
        "vmatches_susceptible_bstrains", ifnotexists=true)
    matches = DataFrame(t=repeat([Float64(time)], length(imatchIDsB)),
        vmatch_id=imatchIDsB, bstrain_id=ibstrainIDs)
    matches |> SQLite.load!(dbOutput,
        "vmatches_immune_bstrains", ifnotexists=true)


    matchIDs = Vector{Int64}()
    strainIDs = Vector{Int64}()
    biomass = Vector{Int64}()
    for phenotype in keys(matchStructure.bstrainclasses)
        matchID = matchStructure.bmatchIDs[phenotype]
        strains = matchStructure.bstrainclasses[phenotype]
        matchIDs = vcat(matchIDs,
            repeat([matchID], length(strains)))
        strainIDs = vcat(strainIDs, strains)
        biomass = vcat(biomass, matchStructure.bphenoBiomass[phenotype])
    end

    matches = DataFrame(t=repeat([Float64(time)], length(matchIDs)),
        match_id=matchIDs, bstrain_id=strainIDs)
    matches |> SQLite.load!(dbOutput,
        "bmatches", ifnotexists=true)
    unique!(matchIDs)
    matches = DataFrame(t=repeat([Float64(time)], length(matchIDs)),
        match_id=matchIDs, babundance=biomass)
    matches |> SQLite.load!(dbOutput,
        "bmatches_abundances", ifnotexists=true)
end

function logPhenotypesDB(matchStructure::tripartite)
    dbOutput = matchStructure.dbOutput
    matchIDs = Vector{Int64}()
    spacerIDs = Vector{Int64}()
    for phenotype in matchStructure.vphenotypes
        matchID = matchStructure.vmatchIDs[phenotype]
        matchIDs = vcat(matchIDs,
            repeat([matchID], length(phenotype)))
        spacerIDs = vcat(spacerIDs, phenotype)
    end
    phenotypes = DataFrame(match_id=matchIDs, phenotype=spacerIDs)
    phenotypes |> SQLite.load!(dbOutput,
        "vmatch_phenotypes", ifnotexists=true)
    matchIDs = Vector{Int64}()
    spacerIDs = Vector{Int64}()
    for phenotype in matchStructure.bphenotypes
        matchID = matchStructure.bmatchIDs[phenotype]
        matchIDs = vcat(matchIDs,
            repeat([matchID], length(phenotype)))
        spacerIDs = vcat(spacerIDs, phenotype)
    end
    phenotypes = DataFrame(match_id=matchIDs, phenotype=spacerIDs)
    phenotypes |> SQLite.load!(dbOutput,
        "bmatch_phenotypes", ifnotexists=true)
end

function tripartiteNetwork!(matchStructure::tripartite)
    for (t,) in execute(dbTempSim, "SELECT DISTINCT t FROM vabundance")
        if t == 0
            continue
        end
        println("Computing tripartite structure at time = $(t)")
        matchStructure.time = t
        structure!(matchStructure)
        updateMatchesDB(matchStructure)
        logTripartiteDB(matchStructure)
    end
    logPhenotypesDB(matchStructure)
end

function createindices()
    println("(Creating run_id indices...)")
    db = SQLite.DB(dbOutputPath)
    execute(db, "BEGIN TRANSACTION")
    for (table_name,) in execute(
        db, "SELECT name FROM sqlite_schema
        WHERE type='table' ORDER BY name;")
        # cols = [info.name for info in execute(db,"PRAGMA table_info($(table_name))")]
        if in(table_name, ["bmatch_phenotypes", "vmatch_phenotypes"])
            execute(db, "CREATE INDEX $(table_name)_index ON $(table_name) (match_id)")
        end
        if in(table_name, ["single_spacer_matches"])
            execute(db, "CREATE INDEX $(table_name)_index ON $(table_name) (t, spacer_id)")
        end
        if in(table_name, ["single_spacer_match_diversity"])
            execute(db, "CREATE INDEX $(table_name)_index ON $(table_name) (t)")
        end
        if in(table_name, ["single_match_tripartite_networks"])
            execute(db, "CREATE INDEX $(table_name)_index ON $(table_name) (t,vstrain_id)")
        end
        if in(table_name, ["vmatches_susceptible_bstrains", "vmatches_immune_bstrains"])
            execute(db, "CREATE INDEX $(table_name)_index ON $(table_name) (t, vmatch_id)")
        end
    end
    execute(db, "COMMIT")
end

matchStructure = tripartite(dbTempMatch, dbTempSim, dbOutput, Float64(0))
tripartiteNetwork!(matchStructure)
createindices()
println("Complete!")





# function logTripartiteDB(matchStructure::tripartite)
#     time = matchStructure.time
#     dbOutput = matchStructure.dbOutput
#     for phenotype in keys(matchStructure.vsinglematches) # this is not a bug
#         matches = matchStructure.vsinglematches[phenotype]
#         vmatchID = matchStructure.vmatchIDs[phenotype]
#         matchIDs = vcat(matchIDs,
#             repeat([vmatchID], length(matches)))
#         spacerIDs = vcat(spacerIDs, matches)
#     end
#     vsingles = DataFrame(t=repeat([Float64(time)], length(matchIDs)),
#         match_id=matchIDs, spacer_id=spacerIDs)
#     vsingles |> SQLite.load!(dbOutput,
#         "vmatch_phenotypes_singles", ifnotexists=true)
#     vsingles = vsingles[:, Not(:t)]
#     rename!(vsingles, :match_id => :vmatch_id)

#     matchIDs = Vector{Int64}()
#     spacerIDs = Vector{Int64}()
#     for phenotype in keys(matchStructure.bstrainclasses) # this is not a bug
#         matches = matchStructure.bsinglematches[phenotype]
#         bmatchID = matchStructure.bmatchIDs[phenotype]
#         matchIDs = vcat(matchIDs,
#             repeat([bmatchID], length(matches)))
#         spacerIDs = vcat(spacerIDs, matches)
#     end
#     bsingles = DataFrame(t=repeat([Float64(time)], length(matchIDs)),
#         match_id=matchIDs, spacer_id=spacerIDs)
#     bsingles |> SQLite.load!(dbOutput,
#         "bmatch_phenotypes_singles", ifnotexists=true)
#     bsingles = bsingles[:, Not(:t)]
#     rename!(bsingles, :match_id => :bmatch_id)

#     tripartiteDF = innerjoin(vsingles, bsingles, on=:spacer_id)
#     pureSingles = nonunique(DataFrame(vmatch_id=tripartiteDF.vmatch_id,
#         bmatch_id=tripartiteDF.bmatch_id))
#     tripartiteDF = tripartiteDF[pureSingles.==false, :]
#     insertcols!(tripartiteDF, 1, :t => repeat([Float64(time)], size(tripartiteDF)[1]))
#     tripartiteDF |> SQLite.load!(dbOutput,
#         "single_match_tripartite_networks", ifnotexists=true)
# end
