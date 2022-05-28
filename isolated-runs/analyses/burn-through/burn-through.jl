#!/usr/bin/env julia

println("(Julia compilation delay...)")

using SQLite
using DataFrames
using SQLite.DBInterface: execute

run_id = ARGS[1]

## Define Paths ##
SCRIPT_PATH = abspath(dirname(PROGRAM_FILE))
#
dbSimPath = joinpath(SCRIPT_PATH,"..","..","..","simulation","sweep_db_gathered.sqlite") # cluster
dbOutputPath = joinpath("burn-through_output.sqlite") # cluster

# dbSimPath = joinpath("/Volumes/Yadgah","crispr-sweep-7-2-2022/isolates/runID3297-c66-r47/runID3297-c66-r47.sqlite") # local
# dbMatchPath = joinpath("/Volumes/Yadgah","crispr-sweep-7-2-2022/isolates/runID3297-c66-r47/matches_output.sqlite") # local
# dbMatchDivPath = joinpath("/Volumes/Yadgah","crispr-sweep-7-2-2022/isolates/runID3297-c66-r47/match-diversity_output.sqlite") # local
# dbProbPath = joinpath("/Volumes/Yadgah/crispr-sweep-7-2-2022/isolates/runID3297-c66-r47/probability-of-emergence_output.sqlite") # local
# dbOutputPath = joinpath("/Volumes/Yadgah/crispr-sweep-7-2-2022/isolates/runID3297-c66-r47/burn-through_output.sqlite") #

if isfile(dbOutputPath)
    error("burn-through_output.sqlite already exists; delete first")
end
##
dbTempSim = SQLite.DB(dbSimPath)
(cr,) = execute(dbTempSim,"SELECT combo_id,replicate FROM runs WHERE run_id = $(run_id)")
dbMatchPath = joinpath(SCRIPT_PATH,"..","..","isolates",
    "runID$(run_id)-c$(cr.combo_id)-r$(cr.replicate)","matches_output.sqlite")
if !isfile(dbMatchPath)
    error("matches_output.sqlite does not exist; compute first")
end
dbMatchDivPath = joinpath(SCRIPT_PATH,"..","..","isolates",
    "runID$(run_id)-c$(cr.combo_id)-r$(cr.replicate)","match-diversity_output.sqlite")
if !isfile(dbMatchDivPath)
    error("match-diversity_output.sqlite does not exist; compute first")
end
dbProbPath = joinpath(SCRIPT_PATH,"..","..","isolates",
    "runID$(run_id)-c$(cr.combo_id)-r$(cr.replicate)","probability-of-emergence_output.sqlite")
if !isfile(dbProbPath)
    error("probability-of-emergence.sqlite does not exist; compute first")
end
dbTriPath = joinpath(SCRIPT_PATH,"..","..","isolates",
    "runID$(run_id)-c$(cr.combo_id)-r$(cr.replicate)","tripartite-networks_output.sqlite")
if !isfile(dbProbPath)
    error("tripartite-networks_output.sqlite does not exist; compute first")
end
dbTempMatch = SQLite.DB(dbMatchPath)
dbTempMatchDiv = SQLite.DB(dbMatchDivPath)
dbTempProb = SQLite.DB(dbProbPath)
dbOutput = SQLite.DB(dbOutputPath)

execute(dbOutput, "CREATE TABLE burn_through_probability (t REAL, p REAL)")

dbTempSim = SQLite.DB()
execute(dbTempSim, "CREATE TABLE babundance (t REAL, bstrain_id, abundance INTEGER)")
execute(dbTempSim, "CREATE TABLE vabundance (t REAL, vstrain_id, abundance INTEGER)")

execute(dbTempSim, "BEGIN TRANSACTION")
execute(dbTempSim,"ATTACH DATABASE '$(dbSimPath)' as dbSim")
execute(dbTempSim,"INSERT INTO babundance (t, bstrain_id, abundance)
    SELECT t, bstrain_id, abundance FROM dbSim.babundance WHERE run_id = $(run_id);")
execute(dbTempSim,"INSERT INTO vabundance (t, vstrain_id, abundance)
    SELECT t, vstrain_id, abundance FROM dbSim.vabundance WHERE run_id = $(run_id);")
execute(dbTempSim, "COMMIT")


execute(dbTempSim, "BEGIN TRANSACTION")
execute(dbTempSim, "CREATE INDEX babundance_index ON babundance (t,bstrain_id,abundance)")
execute(dbTempSim, "CREATE INDEX vabundance_index ON vabundance (t,vstrain_id,abundance)")
execute(dbTempSim, "COMMIT")


function pEmergeSpacers()
    for (t,) in execute(dbTempProb,"SELECT DISTINCT t FROM match_phenotypes")
        println("Computing burn-through probability at t = $(t)")
        for (matchID,) in execute(dbTempMatch, "SELECT DISTINCT vmatch_id
            FROM singe_match_tripartite_networks
            WHERE t = $(t) ORDER BY vmatch_id")
            spacers = [spacerID for (spacerID,) in execute(dbTempMatch,
                "SELECT DISTINCT spacer_id
                FROM singe_match_tripartite_networks
                WHERE t = $(t) AND vmatch_id = $(matchID)")]

            (p,) = execute(dbTempProb,"SELECT p
            FROM pExtinction
            WHERE t = $(t) AND match_id = $(matchClassID)")
        end
        bstrains = [bstrain_id for (bstrain_id,) in
            execute(dbTempMatch, "SELECT bstrain_id FROM bstrain_to_vstrain_matches
            WHERE t = $(t) AND match_length = 1 AND vstrain_id in ($(join(vstrains0,", ")))
            ORDER BY time_specific_match_id")]
        expFreq = 0

        for bstrainID in bstrains
            match_ids = [match_id for (match_id,) in
                execute(dbTempMatch, "SELECT time_specific_match_id FROM bstrain_to_vstrain_matches
                WHERE t = $(t) AND match_length = 1 AND bstrain_id = $(bstrainID)
                ORDER BY time_specific_match_id")]
            spacer_ids = [spacer_ids for (spacer_ids,) in
                execute(dbTempMatch, "SELECT DISTINCT spacer_id FROM matches_spacers
                WHERE t = $(t) AND time_specific_match_id in ($(join(match_ids,", ")))
                ORDER BY time_specific_match_id")]
            matchClassIDs = [matchID for (matchID,) in execute(dbTempProb,"SELECT DISTINCT
                    time_specific_match_id FROM match_phenotypes
                    WHERE t = $(t) AND phenotype in ($(join(spacer_ids,", ")))
                    ORDER BY time_specific_match_id")]
            pBstrain = 0
            for matchClassID in matchClassIDs # this is missing those who are not matched? kind of..
                pheno = [matchID for (matchID,) in execute(dbTempProb,"SELECT
                        phenotype FROM match_phenotypes
                        WHERE t = $(t) AND time_specific_match_id = $(matchClassID)
                        ORDER BY phenotype")]
                (f,) = execute(dbTempProb,"SELECT frequency
                FROM match_phenotype_abundances
                WHERE t = $(t) AND time_specific_match_id = $(matchClassID)")
                (p,) = execute(dbTempProb,"SELECT p
                FROM pExtinction
                WHERE t = $(t) AND match_id = $(matchClassID)")
                pBstrain += f.frequency*(1-p.p)*1/length(pheno) # viral match frequency needs to be normalized!!!!!
            end
            bAbund = sum([abund for (abund,) in
                    execute(dbTempSim,"SELECT abundance
                    FROM babundance
                    WHERE t = $(t) AND bstrain_id = $(bstrainID)")])
            # this is missing those who are already susceptible, but is included in p_ext...
            bTotal = sum([abund for (abund,) in
                    execute(dbTempSim,"SELECT abundance
                    FROM babundance WHERE t = $(t)")])

            expFreq += bAbund/bTotal*pBstrain
        end
        execute(dbOutput, "INSERT INTO burn_through_probability
            VALUES (?,?)",(t,expFreq))
    end
end


pEmergeSpacers()

println("Complete!")

# pSpacer = 0
# for (match_id,) in execute(dbTempProb,"SELECT DISTINCT
#         time_specific_match_id FROM match_phenotypes
#         WHERE phenotype = $(spacerID)
#         ORDER BY time_specific_match_id")
#     phenotype = [spacer for (spacer,) in
#     execute(dbTempProb,"SELECT phenotype FROM match_phenotypes
#     WHERE t = $(t) AND time_specific_match_id = $(match_id)")]
#     if !issubset(phenotype,spacers)
#         continue
#     end
#     (f,) = execute(dbTempProb,"SELECT frequency
#     FROM match_phenotype_abundances
#     WHERE t = $(t) AND time_specific_match_id = $(match_id)")
#     (p,) = execute(dbTempProb,"SELECT p
#     FROM pExtinction
#     WHERE t = $(t) AND time_specific_match_id = $(match_id)")
#     pSpacer += f.frequency*(1-p.p)
# end
# bStrainTotal = sum([abund for (abund,) in
#             execute(dbTempSim,"SELECT abundance
#             FROM babundance WHERE t = $(t)
#             AND bstrain_id in ($(join(bstrains,", ")))")])
# bTotal = sum([abund for (abund,) in
#             execute(dbTempSim,"SELECT abundance
#             FROM babundance WHERE t = $(t)")])
# execute(dbTempProb, "INSERT INTO pEmergeSpacers VALUES (?,?,?,?)",
#     (t,spacerID,pSpacer,bStrainTotal/bTotal))
