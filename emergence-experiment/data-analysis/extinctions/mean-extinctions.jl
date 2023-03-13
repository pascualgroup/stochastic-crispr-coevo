#!/usr/bin/env julia

println("(Julia compilation delay...)")

using SQLite
import SQLite.DBInterface.execute
using DataFrames
using StatsBase
using Statistics

combo_id = ARGS[1]
analysisType = "extinctions"

## Define Paths ##
SCRIPT_PATH = abspath(dirname(PROGRAM_FILE))
dbAnalysisPath = joinpath(SCRIPT_PATH,"..","gathered-analyses",analysisType,"$(analysisType).sqlite") # cluster
dbSimInfoPath = joinpath(SCRIPT_PATH,"..","..","simulation","sweep_db.sqlite") # cluster
dbOutputPath = joinpath("mean-$(analysisType)_output.sqlite") # cluster

#dbAnalysisPath = joinpath("/Volumes/Yadgah/$(analysisType).sqlite") # local
#dbSimInfoPath = joinpath("/Volumes/Yadgah/sweep_db.sqlite") # local
#dbOutputPath = joinpath("/Volumes/Yadgah/mean-$(analysisType)_output.sqlite") # local

if isfile(dbOutputPath)
    error("mean-$(analysisType)_output.sqlite already exists; delete first")
end # cluster
##

dbAnalysis = SQLite.DB(dbAnalysisPath)
dbSimInfo = SQLite.DB(dbSimInfoPath)
dbOutput = SQLite.DB(dbOutputPath)

dbAnalysisTemp = SQLite.DB()
execute(dbAnalysisTemp, "CREATE TABLE extinction_occurrence(microbes INTEGER, viruses INTEGER)")
execute(dbAnalysisTemp, "CREATE TABLE simulation_end_time(microbe_end_time REAL, virus_end_time REAL)")

run_ids = ["$(run_id)" for (run_id,) in execute(dbAnalysis, "SELECT run_id FROM runs WHERE combo_id = $(combo_id) ORDER BY run_id")]
run_idStmt = join(run_ids,", ")
const numReplicates = length(run_ids)

execute(dbAnalysisTemp, "BEGIN TRANSACTION")
execute(dbAnalysisTemp,"ATTACH DATABASE '$(dbAnalysisPath)' as dbAnalysis")
execute(dbAnalysisTemp,
"INSERT INTO extinction_occurrence(microbes, viruses)
SELECT microbes, viruses FROM dbAnalysis.extinction_occurrence WHERE run_id in ($(run_idStmt)) ORDER BY run_id"
)
execute(dbAnalysisTemp, "CREATE INDEX occurrences ON extinction_occurrence (microbes,viruses)")
execute(dbAnalysisTemp,
"INSERT INTO simulation_end_time(microbe_end_time, virus_end_time)
SELECT microbe_end_time, virus_end_time FROM dbAnalysis.simulation_end_time WHERE run_id in ($(run_idStmt)) ORDER BY run_id"
)
execute(dbAnalysisTemp, "CREATE INDEX end_times ON simulation_end_time (microbe_end_time, virus_end_time)")
execute(dbAnalysisTemp, "COMMIT")

function meanExtinctions()
    execute(dbOutput, "CREATE TABLE mean_extinction_occurrences(microbes REAL, viruses REAL, num_microbial_extinction_replicates INTEGER,
    num_viral_extinction_replicates INTEGER, total_replicates_of_combo INTEGER)")

    execute(dbOutput, "CREATE TABLE mean_simulation_end_times(microbe_end_time REAL, virus_end_time REAL)")

    #for (table_name,) in execute(
        #dbAnalysisTemp, "SELECT name FROM sqlite_schema WHERE type='table' ORDER BY name;"
        #)
        #tableCols = ["$(table_info.name)" for table_info in execute(dbAnalysisTemp,"PRAGMA table_info($(table_name))")]
        #tableColsType = ["$(table_info.name) $(table_info.type)" for table_info in execute(dbAnalysisTemp,"PRAGMA table_info($(table_name))")]
        #numCols = length(tableCols)
        #colStmt = join(tableCols,", ")
        #colTypeStmt = join(tableColsType,", ")

        analyses = DataFrame(execute(dbAnalysisTemp, "SELECT microbes,viruses FROM extinction_occurrence"))
        analyses2 = DataFrame(execute(dbAnalysisTemp, "SELECT microbe_end_time, virus_end_time  FROM simulation_end_time"))

        #meanVals = []
        #for i in 1:numCols
            #analysesVec = [analyses[j][i] for j in 1:numReplicates]
            #meanVal = sum(analysesVec)/numReplicates
            #push!(meanVals,"$(meanVal)")
        #end
        #valSpaces = ["?" for i in 1:numCols]
        #valSpaces = join(valSpaces,", ")
        #meanVals = join(meanVals,", ")
        execute(dbOutput, "BEGIN TRANSACTION")
        execute(dbOutput,"INSERT INTO mean_extinction_occurrences VALUES (?,?,?,?,?)",(mean(analyses[:,:microbes]),mean(analyses[:,:viruses]),
        sum(analyses[:,:microbes]),sum(analyses[:,:viruses]),numReplicates))
        execute(dbOutput,"INSERT INTO mean_simulation_end_times VALUES ($(mean(analyses2[:,:microbe_end_time])),$(mean(analyses2[:,:virus_end_time])))")
        execute(dbOutput, "COMMIT")
    #end

end

meanExtinctions()

println("Complete!")
