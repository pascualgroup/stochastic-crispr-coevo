#!/usr/bin/env julia

println("(Annoying Julia compilation delay...)")

using SQLite
import SQLite.DBInterface.execute
using DataFrames
using StatsBase
using Statistics

combo_id = ARGS[1]
analysisType = "walls-shannon"
divType = "shannon"

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
execute(dbAnalysisTemp, "CREATE TABLE microbial_peakwall_count(num_walls INTEGER, run_id INTEGER)")
execute(dbAnalysisTemp, "CREATE TABLE microbial_peakwall_durations(wall_number INTEGER, duration REAL, run_id INTEGER)")
execute(dbAnalysisTemp, "CREATE TABLE threshold_values (upper_threshold INTEGER, lower_threshold INTEGER, $(divType)_threshold REAL)")

run_ids = ["$(run_id)" for (run_id,) in execute(dbAnalysis, "SELECT run_id FROM runs WHERE combo_id = $(combo_id)")]
run_idStmt = join(run_ids,", ")
const numReplicates = length(run_ids)

execute(dbAnalysisTemp, "BEGIN TRANSACTION")
execute(dbAnalysisTemp,"ATTACH DATABASE '$(dbAnalysisPath)' as dbAnalysis")
execute(dbAnalysisTemp,
"INSERT INTO microbial_peakwall_count(num_walls, run_id)
SELECT num_walls, run_id FROM dbAnalysis.microbial_peakwall_count WHERE run_id in ($(run_idStmt))"
)
execute(dbAnalysisTemp,
"INSERT INTO microbial_peakwall_durations(wall_number, duration, run_id)
SELECT wall_number, duration, run_id FROM dbAnalysis.microbial_peakwall_durations WHERE run_id in ($(run_idStmt))"
)
execute(dbAnalysisTemp,
"INSERT INTO threshold_values(upper_threshold, lower_threshold, $(divType)_threshold)
SELECT upper_threshold, lower_threshold, $(divType)_threshold FROM dbAnalysis.threshold_values WHERE run_id in ($(run_ids[1]))"
)
execute(dbAnalysisTemp, "CREATE INDEX wall_counts ON microbial_peakwall_count (num_walls,run_id)")
execute(dbAnalysisTemp, "CREATE INDEX durations ON microbial_peakwall_durations (wall_number,duration,run_id)")
execute(dbAnalysisTemp, "COMMIT")

function wallstats(dbAnalysisTemp)
    execute(dbOutput, "
    CREATE TABLE microbial_mean_walldurations
    (num_walls INTEGER, mean_wallduration REAL, mean_longest_wallduration REAL, num_replicates INTEGER,
    total_replicates_of_combo INTEGER, proportion_of_total_replicates REAL)
    ")
    execute(dbOutput, "
    CREATE TABLE microbial_most_frequent_num_walls
    (most_frequent_num_walls INTEGER, expected_wallduration_of_most_frequents REAL,
    num_replicates_with_most_frequent INTEGER,total_replicates_of_combo INTEGER,
    proportion_of_total_replicates REAL)
    ")
    execute(dbOutput, "
    CREATE TABLE microbial_wall_statistics
    (wall_occurrence REAL, num_replicates_with_walls INTEGER, total_replicates_of_combo INTEGER,
    expected_num_walls REAL, std_num_walls REAL, expected_num_walls_5percentThreshold REAL, std_num_walls_5percentThreshold REAL,
    expected_num_walls_per_replicate_with_walls REAL, expected_num_walls_per_replicate_with_walls_5percentThreshold REAL,
    avg_of_most_frequent_num_walls REAL, mean_wallduration_of_avg_most_frequents REAL)
    ")
    execute(dbOutput, "
    CREATE TABLE threshold_values
    (upper_threshold INTEGER, lower_threshold INTEGER, $(divType)_threshold REAL)
    ")

    execute(dbOutput, "BEGIN TRANSACTION")
    for (numWalls,) in execute(dbAnalysisTemp, "SELECT DISTINCT num_walls FROM microbial_peakwall_count ORDER BY num_walls")
        if numWalls > 0
            reps = ["$(run_id)" for (run_id,) in execute(dbAnalysisTemp, "SELECT run_id FROM microbial_peakwall_count WHERE num_walls = $(numWalls)")]
            numReps = length(reps)
            repStmt = join(reps,", ")
            wallDuration = DataFrame(execute(dbAnalysisTemp,
            "SELECT wall_number, duration, run_id FROM microbial_peakwall_durations WHERE run_id in ($(repStmt))"))
            durations = wallDuration[wallDuration.wall_number.>0,:][:,:duration]
            longests = [maximum(wallDuration[(wallDuration.wall_number.>0) .& (wallDuration.run_id .== run_id),:][:,:duration])
                for run_id in unique(wallDuration[wallDuration.wall_number.>0,:][:,:run_id])]
            execute(dbOutput,"INSERT INTO microbial_mean_walldurations VALUES (?,?,?,?,?,?)",
            (numWalls,mean(durations),mean(longests),numReps,numReplicates,numReps/numReplicates)
            )
        else
            reps = ["$(run_id)" for (run_id,) in execute(dbAnalysisTemp, "SELECT run_id FROM microbial_peakwall_count WHERE num_walls = $(numWalls)")]
            numReps = length(reps)
            repStmt = join(reps,", ")
            execute(dbOutput,"INSERT INTO microbial_mean_walldurations VALUES (?,?,?,?,?,?)",
            (numWalls,0,0,numReps,numReplicates,numReps/numReplicates)
            )
        end
    end
    execute(dbOutput, "COMMIT")

    numWalls = DataFrame(execute(dbAnalysisTemp, "SELECT num_walls,run_id FROM microbial_peakwall_count"))
    mostFrequents = modes(numWalls.num_walls)
    countWalls = countmap(numWalls.num_walls)

    execute(dbOutput, "BEGIN TRANSACTION")
    for i in 1:lastindex(mostFrequents)
        (meanDuration,) = execute(dbOutput,"SELECT mean_wallduration FROM microbial_mean_walldurations WHERE num_walls = $(mostFrequents[i])")
        execute(dbOutput, "INSERT INTO microbial_most_frequent_num_walls
        VALUES (?,?,?,?,?)",
        (mostFrequents[i],meanDuration.mean_wallduration,
        countWalls[mostFrequents[i]],numReplicates,countWalls[mostFrequents[i]]/numReplicates)
        )
    end
    execute(dbOutput, "COMMIT")

    numReps = length(numWalls[(numWalls.num_walls .> 0),:][:,:num_walls])
    if numReps > 0
        wall_occurrence = numReps/numReplicates
        avgNumWallsPerRepWithWall = sum(numWalls[(numWalls.num_walls .> 0),:][:,:num_walls])/numReps
        numWallsThresh = filter(x->countWalls[x]/numReplicates>0.05,numWalls.num_walls)
        durationMostFrequents = [duration for (duration,) in execute(dbOutput,
        "SELECT expected_wallduration_of_most_frequents FROM microbial_most_frequent_num_walls")
        ]
        execute(dbOutput, "BEGIN TRANSACTION")
        execute(dbOutput, "
        INSERT INTO microbial_wall_statistics
        VALUES (?,?,?,?,?,?,?,?,?,?,?)",
        (wall_occurrence, numReps, numReplicates,
        mean(numWalls.num_walls), std(numWalls.num_walls), mean(numWallsThresh), std(numWallsThresh),
        avgNumWallsPerRepWithWall,mean(numWallsThresh[numWallsThresh .> 0]),mean(mostFrequents),mean(durationMostFrequents))
        )
        execute(dbOutput, "COMMIT")
    else
        execute(dbOutput, "BEGIN TRANSACTION")
        execute(dbOutput, "
        INSERT INTO microbial_wall_statistics
        VALUES (?,?,?,?,?,?,?,?,?,?,?)",
        (0, numReps, numReplicates, 0, 0, 0, 0, 0, 0, 0, 0))
        execute(dbOutput, "COMMIT")
    end


    (threshold,) = execute(dbAnalysisTemp,"SELECT * FROM threshold_values")

    execute(dbOutput, "BEGIN TRANSACTION")
    execute(dbOutput, "
    INSERT INTO threshold_values
    VALUES (?,?,?)",
    (threshold[1],threshold[2],threshold[3])
    )
    execute(dbOutput, "COMMIT")

end

wallstats(dbAnalysisTemp)

println("Complete!")

# for (table_name,) in execute(
#     dbAnalysisTemp, "SELECT name FROM sqlite_schema WHERE type='table' ORDER BY name;"
#     )
#     if table_name == "microbial_peakwall_count"
#         tableCols = ["$(table_info.name)" for table_info in execute(dbAnalysisTemp,"PRAGMA table_info($(table_name))")]
#         tableColsType = ["$(table_info.name) $(table_info.type)" for table_info in execute(dbAnalysisTemp,"PRAGMA table_info($(table_name))")]
#         numCols = length(tableCols)
#         colStmt = join(tableCols,", ")
#         colTypeStmt = join(tableColsType,", ")
#         execute(dbOutput, "CREATE TABLE $(table_name) ($(colTypeStmt...))")
#
#         analyses = execute(dbAnalysisTemp, "SELECT $(colStmt...) FROM $(table)")
#         meanVals = []
#         for i in 1:numCols
#             analysesVec = [analyses[j][i] for j in 1:numReplicates]
#             meanVal = sum(analysesVec)/numReplicates
#             push!(meanVals,meanVal)
#         end
#         valSpaces = ["?" for 1:numCols]
#         valSpaces = join(valSpaces,", ")
#         execute(dbOutput, "BEGIN TRANSACTION")
#         execute(dbOutput,"INSERT INTO $(table_name) VALUES ($(valSpaces))",(meanVals...))
#         execute(dbOutput, "COMMIT")
# end
