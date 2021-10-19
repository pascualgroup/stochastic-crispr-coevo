#!/usr/bin/env julia

println("(Annoying Julia compilation delay...)")

using SQLite
using DataFrames
using SQLite.DBInterface: execute

run_id = ARGS[1]

upperThreshold = Int64(parse(Float64,ARGS[2]))
lowerThreshold = Int64(parse(Float64,ARGS[3]))
hill2Threshold = parse(Float64,ARGS[4])

SCRIPT_PATH = abspath(dirname(PROGRAM_FILE))

dbOutputPath = joinpath("walls_output.sqlite") # cluster
#dbOutputPath = joinpath("/Volumes/Yadgah/walls_output.sqlite") # local
if isfile(dbOutputPath)
    error("$(dbOutputPath) already exists")
end
# Create database for output data
dbOutput = SQLite.DB(dbOutputPath)

thresholdVals =  DataFrame(upper_threshold = upperThreshold,lower_threshold = lowerThreshold, hill2_threshold = hill2Threshold)
thresholdVals |> SQLite.load!(dbOutput,"threshold_values",ifnotexists=true)

# This function counts number of peaks, logs time series of peaks their respective durations
function peakwallCount(upperThreshold,lowerThreshold,hill2Threshold)

    ## Define Paths ##
    SCRIPT_PATH = abspath(dirname(PROGRAM_FILE))

    dbSimPath = joinpath(SCRIPT_PATH,"..","..","simulation","sweep_db_gathered.sqlite") # cluster
    #run(`cd`) # local
    #dbSimPath = joinpath("/Volumes/Yadgah/sweep_db_gathered.sqlite") # local

    #dbShanPath = joinpath("/Volumes/Yadgah/shannon.sqlite") # local
    dbShanPath = joinpath(SCRIPT_PATH,"..","gathered-analyses","shannon","shannon.sqlite") # cluster
    ##

    # Connect to time series data of abundances
    dbSim = SQLite.DB(dbSimPath)

    # Save as data frame
    series = DataFrame(execute(dbSim, "SELECT t, microbial_abundance FROM summary WHERE run_id = $(run_id)"))

    # Identify crossings of thresholds
    coarseSeries = series[(series.microbial_abundance .>= upperThreshold) .| (series.microbial_abundance .<= lowerThreshold),:]
    upper = 2*(series.microbial_abundance .>= upperThreshold); # 2 signifies values above upper threshold
    lower = series.microbial_abundance .<= lowerThreshold; # 1 signifies values below lower threshold, 0 for values between lower and upper threshold,
    coarseSeries = upper + lower

    # Initialize each table of output database as data frame
    peakSeries = DataFrame(t=series.t, peak_presence=Array{Int64,1}(zeros(length(series.t))))
    peakDurations =  DataFrame(begin_row = Int64[],end_row = Int64[], begin_t = Float64[], end_t = Float64[],
    peak_number=Int64[], wall_number=Int64[], duration=Float64[])
    wallSeries = DataFrame(t=series.t, wall_presence=Array{Int64,1}(zeros(length(series.t))), bhill2 = Array{Float64,1}(zeros(length(series.t))), vhill2 = Array{Float64,1}(zeros(length(series.t))))

    # Create temporary database that is a copy of the main database at the run_id value of the script's argument
    #dbTemp = SQLite.DB("/Volumes/Yadgah/timeSeries$(run_id).sqlite") # local
    dbTemp = SQLite.DB()
    execute(dbTemp, "CREATE TABLE hill_no2 (t REAL, vhill2 REAL, bhill2 REAL)")

    execute(dbTemp, "BEGIN TRANSACTION")
    execute(dbTemp,"ATTACH DATABASE '$(dbShanPath)' as dbShan")
    execute(dbTemp,"INSERT INTO hill_no2(t, vhill2, bhill2) SELECT t, vhill2, bhill2 FROM dbShan.hill_no2 WHERE run_id = $(run_id);")
    execute(dbTemp, "COMMIT")

    execute(dbTemp, "BEGIN TRANSACTION")
    execute(dbTemp, "CREATE INDEX hill2_index ON hill_no2 (t, vhill2, bhill2)")
    execute(dbTemp, "COMMIT")

    peaks = 0
    walls = 0
    check = 0
    j1,j2,j3,j4 = 0,0,0,0

    for i in 1:length(coarseSeries)
        println("Scanning through time series point $(i)")
        if i == length(coarseSeries)
            return peaks, walls, peakSeries, peakDurations, wallSeries
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
                peaks = peaks + 1

                peakSeries.peak_presence[j1:j4] = Array{Int64,1}(ones(length(j1:j4)))
                peak_duration = peakSeries.t[j4] - peakSeries.t[j1]

                timeStmt = map(x->string(" OR t = $(x)"), peakSeries.t[(j1+1):j4])
                timeStmt = append!([string("t = $(peakSeries.t[j1])")],timeStmt)
                bhill2 = [hill.bhill2 for hill in execute(dbTemp, "SELECT bhill2 FROM hill_no2 WHERE $(timeStmt...) ORDER BY t")]
                vhill2 = [hill.vhill2 for hill in execute(dbTemp, "SELECT vhill2 FROM hill_no2 WHERE $(timeStmt...) ORDER BY t")]

                if maximum(bhill2) >= hill2Threshold
                    walls = walls + 1
                    wallSeries.wall_presence[j1:j4] = Array{Int64,1}(ones(length(j1:j4)))
                    wallSeries.bhill2[j1:j4] = bhill2[:,:]
                    wallSeries.vhill2[j1:j4] = vhill2[:,:]
                    push!(peakDurations,[j1 j4 peakSeries.t[j1] peakSeries.t[j4] peaks walls peak_duration])
                else
                    push!(peakDurations,[j1 j4 peakSeries.t[j1] peakSeries.t[j4] peaks 0 peak_duration])
                end

                check = 0
            end
        end
    end
    return peaks, walls, peakSeries, peakDurations, wallSeries
#end
end

peaks, walls, peakSeries, peakDurations, wallSeries = peakwallCount(upperThreshold,lowerThreshold,hill2Threshold);

count =  DataFrame(num_peaks = peaks, num_walls = walls)

peakSeries |> SQLite.load!(dbOutput,"microbial_peak_series",ifnotexists=true)
peakDurations |> SQLite.load!(dbOutput,"microbial_peakwall_durations",ifnotexists=true)
count |> SQLite.load!(dbOutput,"microbial_peakwall_count",ifnotexists=true)
wallSeries |> SQLite.load!(dbOutput,"microbial_wall_series",ifnotexists=true)

dbSimInfoPath = joinpath(SCRIPT_PATH,"..","..","simulation","sweep_db.sqlite") # cluster
#run(`cd`) # local
#dbSimInfoPath = joinpath("/Volumes/Yadgah/sweep_db.sqlite") # local
dbSimInfo = SQLite.DB(dbSimInfoPath)
tableNamesTypes = ["$(table_info.name) $(table_info.type)" for table_info in execute(dbSimInfo,"PRAGMA table_info(param_combos)")]
tableNamesTypes = join(tableNamesTypes,", ")

execute(dbOutput, "CREATE TABLE runs (run_id INTEGER, combo_id INTEGER, replicate INTEGER)")
execute(dbOutput, "CREATE TABLE param_combos ($(tableNamesTypes...))")

tableNames = ["$(table_info.name)" for table_info in execute(dbSimInfo,"PRAGMA table_info(param_combos)")]
tableNames = join(tableNames,", ")
execute(dbOutput, "BEGIN TRANSACTION")
execute(dbOutput,"ATTACH DATABASE '$(dbSimInfoPath)' as dbSimInfo")
execute(dbOutput,"INSERT INTO param_combos($(tableNames)) SELECT * FROM dbSimInfo.param_combos")
execute(dbOutput,"INSERT INTO runs (run_id, combo_id, replicate) SELECT run_id, combo_id, replicate FROM dbSimInfo.runs")
execute(dbOutput, "COMMIT")

println("Complete!")
