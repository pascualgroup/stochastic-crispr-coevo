#!/usr/bin/env julia

println("(Annoying Julia compilation delay...)")

using SQLite
using DataFrames
using SQLite.DBInterface: execute

run_id = ARGS[1] # cluster & local

upperThreshold = Int64(parse(Float64,ARGS[2]))
lowerThreshold = Int64(parse(Float64,ARGS[3]))

#dbOutputPath = joinpath("peaks_output.sqlite") # cluster
dbOutputPath = joinpath("/Volumes/Yadgah/peaks_output.sqlite") # local
if isfile(dbOutputPath)
    error("$(dbOutputPath) already exists")
end
dbOutput = SQLite.DB(dbOutputPath)

thresholdVals =  DataFrame(upper_threshold = upperThreshold,lower_threshold = lowerThreshold)
thresholdVals |> SQLite.load!(dbOutput,"threshold_values",ifnotexists=true)

#upperThreshold = 3e5; # for testing
#lowerThreshold = 1.2e5; # for testing


# This function counts number of peaks, logs time series of peaks their respective durations
function peakCount(upperThreshold,lowerThreshold)
    ## Define Paths ##
    SCRIPT_PATH = abspath(dirname(PROGRAM_FILE))

    #dbSimPath = joinpath(SCRIPT_PATH,"..","..","simulation","sweep_db_gathered.sqlite") # cluster
    dbSimPath = joinpath("/Volumes/Yadgah/sweep_db_gathered.sqlite") # local
    ##

    dbSim = SQLite.DB(dbSimPath)

    series = DataFrame(execute(dbSim, "SELECT t, microbial_abundance FROM summary WHERE run_id = $(run_id)"))

    # Identify crossings of thresholds
    coarseSeries = series[(series.microbial_abundance .>= upperThreshold) .| (series.microbial_abundance .<= lowerThreshold),:]
    upper = 2*(series.microbial_abundance .>= upperThreshold); # 2 signifies values above upper threshold
    lower = series.microbial_abundance .<= lowerThreshold; # 1 signifies values below lower threshold, 0 for values between lower and upper threshold,
    coarseSeries = upper + lower

    # Initialize data frames
    peakSeries = DataFrame(t=series.t, peak_presence=Array{Int64,1}(zeros(length(series.t))))
    peakDurations =  DataFrame(begin_row = Int64[],end_row = Int64[], begin_t = Float64[], end_t = Float64[],
    peak_number=Int64[], duration=Float64[])

    peaks = 0
    check = 0
    j1,j2,j3,j4 = 0,0,0,0
    for i in 1:length(coarseSeries)
        println("Scanning through time series point $(i)")
        if i == length(coarseSeries)
            return peaks, peakSeries, peakDurations
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
    return peaks, peakSeries, peakDurations
end

peaks, peakSeries, peakDurations = peakCount(upperThreshold,lowerThreshold);

numpeaks =  DataFrame(num_peaks = peaks)

peakSeries |> SQLite.load!(dbOutput,"microbial_peak_series",ifnotexists=true)
peakDurations |> SQLite.load!(dbOutput,"microbial_peak_durations",ifnotexists=true)
numpeaks |> SQLite.load!(dbOutput,"microbial_num_peaks",ifnotexists=true)
