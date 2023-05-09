#!/usr/bin/env julia

println("(Julia compilation delay...)")

using SQLite
using DataFrames
using SQLite.DBInterface: execute
using SQLite: DB
using DataFrames

run_id = ARGS[1]
## Define Paths ##
SCRIPT_PATH = abspath(dirname(PROGRAM_FILE))
#
dbSimPath = joinpath(SCRIPT_PATH, "..", "..", "..", "simulation", "sweep_db_gathered.sqlite") # cluster
dbOutputPath = joinpath("experiment_output.sqlite") # cluster
# dbSimPath = joinpath("//Volumes/Yadgah/crispr-sweep-7-2-2022/emergence-experiment/simulation/sweep_db_gathered.sqlite") # local
# dbOutputPath = joinpath("/Volumes/Yadgah/experiment_output.sqlite") # local

if isfile(dbOutputPath)
    error("experiment.sqlite already exists; delete first")
end
##
dbTempSim = SQLite.DB(dbSimPath)
dbOutput = SQLite.DB(dbOutputPath)
vstrainID = [vID for (vID,) in execute(dbTempSim, "SELECT DISTINCT vstrain_id 
                FROM vabundance WHERE t = 0")]
function identifyInitialMatches(vstrainID)
    matches = Int64[]
    match_length = Int64[]
    vspacers = [sID for (sID,) in
                execute(
        dbTempSim,
        "SELECT DISTINCT spacer_id FROM vpspacers 
        WHERE vstrain_id = $(vstrainID)")]
    for (bstrainID,) in execute(
        dbTempSim,
        "SELECT DISTINCT bstrain_id FROM babundance
        WHERE t = 0 AND abundance != 0 ORDER BY bstrain_id")
        println("bstrain: $(bstrainID)")
        bspacers = [sID for (sID,) in
                    execute(dbTempSim,"SELECT DISTINCT spacer_id FROM bspacers 
                    WHERE bstrain_id = $(bstrainID)")]
        matchLength = length(intersect(bspacers, vspacers))
        append!(matches, Int64(bstrainID))
        append!(match_length, Int64(matchLength))
    end
    df = DataFrame(vstrain_id=vstrainID * ones(Int64, length(matches)), bstrain_id=matches, match_length=match_length)
    df |> SQLite.load!(dbOutput,
        "initial_matches", ifnotexists=true)
    return df
end
execute(dbOutput, "CREATE TABLE emergence (emergence INTEGER, run_id INTEGER)")
initialstrains = identifyInitialMatches(vstrainID)
# initialstrains = matches[matches.match_length.>0, :][:, :bstrain_id]
initialstrains = matches[matches.match_length.>=0, :][:, :bstrain_id]
for (run_id,) in execute(dbTempSim,"SELECT run_id FROM runs")
    println("runID: $(run_id)")
    maxT = maximum([t for (t,) in execute(
        dbTempSim, "SELECT DISTINCT t FROM babundance 
        WHERE run_id = $(run_id) AND abundance != 0")])
    finalstrains = [bID for (bID,) in execute(
        dbTempSim, "SELECT bstrain_id FROM babundance 
        WHERE run_id = $(run_id) AND abundance != 0 AND t = $(maxT)")]
    remaining = intersect(initialstrains, finalstrains)
    if length(remaining) < length(initialstrains)
        execute(dbOutput,"INSERT INTO emergence VALUES (?,?)",(1,run_id))
    else
        execute(dbOutput,"INSERT INTO emergence VALUES (?,?)",(0,run_id))
    end
end

escape()
println("Complete!")

zeros = matches[matches.match_length.==0, :][:, :bstrain_id]