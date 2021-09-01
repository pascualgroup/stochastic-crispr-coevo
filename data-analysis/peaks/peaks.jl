#!/usr/bin/env julia

println("(Annoying Julia compilation delay...)")

using SQLite
using DataFrames
using SQLite.DBInterface: execute

run_id = ARGS[1] # cluster

upperThreshold = ARGS[2]
lowerThreshold = ARGS[3]

#upperThreshold = 3e5; # for testing
#lowerThreshold = 1.2e5; # for testing

SCRIPT_PATH = abspath(dirname(PROGRAM_FILE))

#dbSimPath = joinpath(SCRIPT_PATH,"..","..","simulation","sweep_db_gathered.sqlite") # cluster
dbSimPath = joinpath("/Volumes/Yadgah/sweep_db_gathered.sqlite") # local machine

dbSim = SQLite.DB(dbSimPath)

series = DataFrame(execute(dbSim, "SELECT t, microbial_abundance FROM summary WHERE run_id = $(run_id)"))

#dbOutputPath = joinpath("peaks_output.sqlite") # cluster
#run(`cd`) # local machine
dbOutputPath = joinpath("/Volumes/Yadgah/peaks_output.sqlite") # local machine
if isfile(dbOutputPath)
    error(dbOutputPath)
end

dbOutput = SQLite.DB(dbOutputPath)

#dbRichPath = joinpath(SCRIPT_PATH,"..","gathered-analyses","richness","richness.sqlite") # cluster
#if !isfile(joinpath(dbRichPath) # cluster
    #error("richness.sqlite is missing; compute richness first") # cluster
#end # cluster
dbRichPath = joinpath("/Volumes/Yadgah/richness.sqlite") # local machine

##The following commands ar enot necessary as sqlite.load! automatically creates the table
##execute(dbOutput, "CREATE TABLE microbial_peak_series (t REAL, peak_presence INTEGER)")
##execute(dbOutput, "CREATE TABLE microbial_peak_durations (t REAL, peak_presence INTEGER)")
##execute(dbOutput, "CREATE TABLE microbial_num_peaks (t REAL, peak_presence INTEGER)")
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
peakSeries = DataFrame(t=series.t, peak_presence=Array{Int64,1}(zeros(length(series.t))))
peakDurations =  DataFrame(begin_row = Int64[],end_row = Int64[], begin_t = Float64[], end_t = Float64[],
peak_number=Int64[], duration=Float64[])
numpeaks =  DataFrame(num_peaks = Int64[])

# This function counts number of peaks, logs time series of peaks their respective durations
function peakCount(coarseSeries,peakSeries)
    #let
        peaks = 0
        check = 0
        j1,j2,j3,j4 = 0,0,0,0
        for i in 1:length(coarseSeries)
            println("Scanning through time series point $(i)")
            if i == length(coarseSeries)
                return peaks
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
                    push!(peakDurations,[j1 j4 peakSeries.t[j1] peakSeries.t[j4] peaks peak_duration])

                    check = 0
                end
            end
        end
    return peaks
    #end
end

peaks = peakCount(coarseSeries,peakSeries);

push!(numpeaks,[Int64(peaks)])

peakSeries |> SQLite.load!(dbOutput,"microbial_peak_series",ifnotexists=true)
peakDurations |> SQLite.load!(dbOutput,"microbial_peak_durations",ifnotexists=true)
numpeaks |> SQLite.load!(dbOutput,"microbial_num_peaks",ifnotexists=true)