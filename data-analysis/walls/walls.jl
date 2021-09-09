#!/usr/bin/env julia

println("(Annoying Julia compilation delay...)")

using SQLite
using DataFrames
using SQLite.DBInterface: execute

run_id = ARGS[1] # cluster & local
#run_id = 1 # local

upperThreshold = ARGS[2]
lowerThreshold = ARGS[3]

thresholdVals =  DataFrame(upper_threshold = upperThreshold,lower_threshold = lowerThreshold)

#upperThreshold = 3e5; # for testing
#lowerThreshold = 1.2e5; # for testing

SCRIPT_PATH = abspath(dirname(PROGRAM_FILE))

#dbSimPath = joinpath(SCRIPT_PATH,"..","..","simulation","sweep_db_gathered.sqlite") # cluster
dbSimPath = joinpath("/Volumes/Yadgah/sweep_db_gathered.sqlite") # local machine

dbSim = SQLite.DB(dbSimPath)

series = DataFrame(execute(dbSim, "SELECT t, microbial_abundance FROM summary WHERE run_id = $(run_id)"))

#dbOutputPath = joinpath("walls_output.sqlite") # cluster
#run(`cd`) # local machine
dbOutputPath = joinpath("/Volumes/Yadgah/walls_output.sqlite") # local machine
if isfile(dbOutputPath)
    error("$(dbOutputPath) already exists")
end

dbOutput = SQLite.DB(dbOutputPath)

#dbRichPath = joinpath(SCRIPT_PATH,"..","gathered-analyses","richness","richness.sqlite") # cluster
#if !isfile(joinpath(dbRichPath) # cluster
    #error("richness.sqlite is missing; compute richness first") # cluster
#end # cluster
dbRichPath = joinpath("/Volumes/Yadgah/richness.sqlite") # local machine

##The following commands ar enot necessary as sqlite.load! automatically creates the table
##execute(dbOutput, "CREATE TABLE microbial_wall_series (t REAL, wall_presence INTEGER)")
##execute(dbOutput, "CREATE TABLE microbial_wall_durations (t REAL, wall_presence INTEGER)")
##execute(dbOutput, "CREATE TABLE microbial_num_walls (t REAL, wall_presence INTEGER)")
##

# Create temporary database that is a copy of the main database at the run_id value of the script's argument
dbTempRich = SQLite.DB()
execute(dbTempRich, "CREATE TABLE richness (t REAL, vrichness INTEGER, brichness INTEGER)")

execute(dbTempRich, "BEGIN TRANSACTION")
execute(dbTempRich,"ATTACH DATABASE '$(dbRichPath)' as dbRich")
execute(dbTempRich,"INSERT INTO richness(t, vrichness, brichness) SELECT t, vrichness, brichness FROM dbRich.richness WHERE run_id = $(run_id);")
execute(dbTempRich, "COMMIT")

execute(dbTempRich, "BEGIN TRANSACTION")
execute(dbTempRich, "CREATE INDEX richness_index ON richness (t, vrichness, brichness)")
execute(dbTempRich, "COMMIT")


richness = DataFrame(execute(dbTempRich, "SELECT t, vrichness, brichness FROM richness"))


# Identify crossings of thresholds
coarseSeries = series[(series.microbial_abundance .>= upperThreshold) .| (series.microbial_abundance .<= lowerThreshold),:]
upper = 2*(series.microbial_abundance .>= upperThreshold); # 2 signifies values above upper threshold
lower = series.microbial_abundance .<= lowerThreshold; # 1 signifies values below lower threshold, 0 for values between lower and upper threshold,
coarseSeries = upper + lower

# Initialize data frames
wallSeries = DataFrame(t=series.t, wall_presence=Array{Int64,1}(zeros(length(series.t))))
wallDurations =  DataFrame(begin_row = Int64[],end_row = Int64[], begin_t = Float64[], end_t = Float64[],
wall_number=Int64[], duration=Float64[])
numWalls =  DataFrame(num_walls = Int64[])

# This function counts number of walls, logs time series of walls their respective durations
function wallCount(coarseSeries,wallSeries)
    #let
        walls = 0
        check = 0
        j1,j2,j3,j4 = 0,0,0,0
        for i in 1:length(coarseSeries)
            println("Scanning through time series point $(i)")
            if i == length(coarseSeries)
                return walls
            end
            if coarseSeries[i] == 1 && coarseSeries[i+1] == 0
                check = 1
                j1 = i
            end
            if coarseSeries[i] == 0 && coarseSeries[i+1] == 2
                if check == 1
                    check = 2
                    j2 = i
                end
            end
            if coarseSeries[i] == 2 && coarseSeries[i+1] == 0
                if check == 2
                    check = 3
                    j3 = i
                end
            end
            if coarseSeries[i] == 0 && coarseSeries[i+1] == 1
                if check == 3
                    j4 = i
                    walls = walls + 1

                    wallSeries.wall_presence[j1:j4] = Array{Int64,1}(ones(length(j1:j4)))

                    wall_duration = wallSeries.t[j4] - wallSeries.t[j1]
                    push!(wallDurations,[j1 j4 wallSeries.t[j1] wallSeries.t[j4] walls wall_duration])

                    check = 0
                end
            end
        end
    return walls
    #end
end

walls = wallCount(coarseSeries,wallSeries);

push!(numWalls,[Int64(walls)])

thresholdVals |> SQLite.load!(dbOutput,"threshold_values",ifnotexists=true)
wallSeries |> SQLite.load!(dbOutput,"microbial_wall_series",ifnotexists=true)
wallDurations |> SQLite.load!(dbOutput,"microbial_wall_durations",ifnotexists=true)
numWalls |> SQLite.load!(dbOutput,"microbial_num_walls",ifnotexists=true)
