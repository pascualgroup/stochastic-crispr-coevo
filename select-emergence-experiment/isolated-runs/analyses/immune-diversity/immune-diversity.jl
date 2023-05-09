#!/usr/bin/env julia

println("(Julia compilation delay...)")

using SQLite
using DataFrames
using SQLite.DBInterface: execute

run_id = ARGS[1]

uPercent = parse(Float64, ARGS[2]) / 100
lPercent = parse(Float64, ARGS[3]) / 100
shannonThreshold = parse(Float64, ARGS[4])

##
SCRIPT_PATH = abspath(dirname(PROGRAM_FILE))
dbSimPath = joinpath(SCRIPT_PATH, "..", "..", "..", "simulation", "sweep_db_gathered.sqlite") # cluster
dbOutputPath = joinpath("walls-extinctions_output.sqlite") # cluster
# dbSimPath = joinpath("/Volumes/Yadgah","crispr-sweep-7-2-2022/isolates/runID3297-c66-r47/runID3297-c66-r47.sqlite") # local
# dbOutputPath = joinpath("/Volumes/Yadgah/crispr-sweep-7-2-2022/isolates/runID3297-c66-r47/tripartite-networks_output.sqlite") # local

if isfile(dbOutputPath)
    error("walls-extinctions_output.sqlite already exists; delete first")
end
##
dbTempSim = SQLite.DB(dbSimPath)
(cr,) = execute(dbTempSim, "SELECT combo_id,replicate FROM runs WHERE run_id = $(run_id)")
if length(ARGS) > 1 && ARGS[end] == "combo"
    if isfile(joinpath(SCRIPT_PATH, "..", "..", "..", "data-analysis",
        "gathered-analyses", "shannon", "shannon.sqlite"))
        dbShanPath = joinpath(SCRIPT_PATH, "..", "..", "..", "data-analysis",
            "gathered-analyses", "shannon", "shannon.sqlite")
    elseif isfile(joinpath(SCRIPT_PATH, "..", "..", "isolates",
        "comboID$(cr[1])", "shannonC$(cr[1]).sqlite"))
        dbShanPath = joinpath(SCRIPT_PATH, "..", "..", "isolates",
            "comboID$(cr[1])", "shannonC$(cr[1]).sqlite")
    else
        error("neither /data-analysis/gathered-analyses/shannon/shannon.sqlite exists \n
            nor does /isolated-runs/isolates/comboID$(cr[1])/shannonC$(cr[1]).sqlite; \n
            compute one of them first")
    end
end
if length(ARGS) == 1
    dbShanPath = joinpath(SCRIPT_PATH, "..", "..", "isolates",
        "runID$(run_id)-c$(cr.combo_id)-r$(cr.replicate)", "shannon_output.sqlite")
    if !isfile(dbShanPath)
        error("shannon_output.sqlite does not exist; compute first")
    end

end
dbTempShan = SQLite.DB(dbShanPath)
dbTempTri = SQLite.DB(dbTriPath)
dbOutput = SQLite.DB(dbOutputPath)
# dbTempSim = SQLite.DB(dbSimPath) # local
dbTempSim = SQLite.DB()
execute(dbTempSim, "CREATE TABLE babundance (t REAL, bstrain_id, abundance INTEGER)")
execute(dbTempSim, "CREATE TABLE vabundance (t REAL, vstrain_id, abundance INTEGER)")
execute(dbTempSim, "BEGIN TRANSACTION")
execute(dbTempSim, "ATTACH DATABASE '$(dbSimPath)' as dbSim")
execute(
    dbTempSim,
    "INSERT INTO babundance (t, bstrain_id, abundance)
SELECT t, bstrain_id, abundance FROM dbSim.babundance WHERE run_id = $(run_id);"
)
execute(
    dbTempSim,
    "INSERT INTO vabundance (t, vstrain_id, abundance)
SELECT t, vstrain_id, abundance FROM dbSim.vabundance WHERE run_id = $(run_id);"
)
execute(dbTempSim, "COMMIT")
execute(dbTempSim, "BEGIN TRANSACTION")
execute(dbTempSim, "CREATE INDEX babundance_index ON babundance (t,bstrain_id,abundance)")
execute(dbTempSim, "CREATE INDEX vabundance_index ON vabundance (t,vstrain_id,abundance)")
execute(dbTempSim, "COMMIT")




# This function counts number of peaks, logs time series of peaks their respective durations
function peakwallCount(uPercent, lPercent, shannonThreshold, dbOutput, dbSim)
    (comboID,) = execute(dbSim, "SELECT combo_id FROM runs WHERE run_id = $(run_id)")
    comboID = comboID.combo_id
    (cCapacity,) = execute(dbSim, "SELECT microbe_carrying_capacity FROM param_combos WHERE combo_id = $(comboID)")
    cCapacity = cCapacity.microbe_carrying_capacity

    # Save as data frame
    series = DataFrame(execute(dbSim, "SELECT t, microbial_abundance FROM summary WHERE run_id = $(run_id)"))

    upperThreshold = floor(cCapacity * uPercent)
    lowerThreshold = floor(cCapacity * lPercent)

    thresholdVals = DataFrame(carrying_capacity=cCapacity,
        upper_percent=uPercent, lower_percent=lPercent,
        upper_threshold=upperThreshold, lower_threshold=lowerThreshold,
        shannon_threshold=shannonThreshold)
    thresholdVals |> SQLite.load!(dbOutput, "threshold_values", ifnotexists=true)

    # Identify crossings of thresholds
    coarseSeries = series[(series.microbial_abundance.>=upperThreshold).|(series.microbial_abundance.<=lowerThreshold), :]
    upper = 2 * (series.microbial_abundance .>= upperThreshold) # 2 signifies values above upper threshold
    lower = series.microbial_abundance .<= lowerThreshold # 1 signifies values below lower threshold, 0 for values between lower and upper threshold,
    coarseSeries = upper + lower

    # Initialize each table of output database as data frame
    peakSeries = DataFrame(t=series.t, peak_presence=Array{Int64,1}(zeros(length(series.t))))
    peakDurations = DataFrame(begin_row=Int64[], end_row=Int64[], begin_t=Float64[], end_t=Float64[],
        peak_number=Int64[], wall_number=Int64[], duration=Float64[])
    wallSeries = DataFrame(t=series.t, wall_presence=Array{Int64,1}(zeros(length(series.t))), bshannon=Array{Float64,1}(zeros(length(series.t))), vshannon=Array{Float64,1}(zeros(length(series.t))))

    # Create temporary database that is a copy of the main database at the run_id value of the script's argument
    #dbTemp = SQLite.DB("/Volumes/Yadgah/timeSeries$(run_id).sqlite") # local
    dbTemp = SQLite.DB()
    execute(dbTemp, "CREATE TABLE shannon_diversity (t REAL, vshannon REAL, bshannon REAL)")

    execute(dbTemp, "BEGIN TRANSACTION")
    execute(dbTemp, "ATTACH DATABASE '$(dbShanPath)' as dbShan")
    execute(dbTemp, "INSERT INTO shannon_diversity(t, vshannon, bshannon) SELECT t, vshannon, bshannon FROM dbShan.shannon_diversity WHERE run_id = $(run_id);")
    # execute(dbTemp,"INSERT INTO shannon_diversity(t, vshannon, bshannon) SELECT t, vshannon, bshannon FROM dbShan.shannon_diversity;")
    execute(dbTemp, "COMMIT")

    execute(dbTemp, "BEGIN TRANSACTION")
    execute(dbTemp, "CREATE INDEX shannon_index ON shannon_diversity (t, vshannon, bshannon)")
    execute(dbTemp, "COMMIT")

    peaks = 0
    walls = 0
    check = 0
    j1, j2, j3, j4 = 0, 0, 0, 0

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

                timeStmt = map(x -> string(" OR t = $(x)"), peakSeries.t[(j1+1):j4])
                timeStmt = append!([string("t = $(peakSeries.t[j1])")], timeStmt)
                bshannon = [shannon.bshannon for shannon in execute(dbTemp, "SELECT bshannon FROM shannon_diversity WHERE $(timeStmt...) ORDER BY t")]
                vshannon = [shannon.vshannon for shannon in execute(dbTemp, "SELECT vshannon FROM shannon_diversity WHERE $(timeStmt...) ORDER BY t")]

                if maximum(bshannon) >= shannonThreshold
                    walls = walls + 1
                    wallSeries.wall_presence[j1:j4] = Array{Int64,1}(ones(length(j1:j4)))
                    wallSeries.bshannon[j1:j4] = bshannon[:, :]
                    wallSeries.vshannon[j1:j4] = vshannon[:, :]
                    push!(peakDurations, [j1 j4 peakSeries.t[j1] peakSeries.t[j4] peaks walls peak_duration])
                else
                    push!(peakDurations, [j1 j4 peakSeries.t[j1] peakSeries.t[j4] peaks 0 peak_duration])
                end

                check = 0
            end
        end
    end
    return peaks, walls, peakSeries, peakDurations, wallSeries
    #end
end
function extinction()
    microbesDF = DataFrame(execute(dbTemp, "SELECT t,microbial_abundance FROM summary ORDER BY t"))
    virusesDF = DataFrame(execute(dbTemp, "SELECT t,viral_abundance FROM summary ORDER by t"))
    microbesDF = microbesDF[(microbesDF.t.!=0), :]
    virusesDF = virusesDF[(virusesDF.t.!=0), :]
    microbeExt = 0
    virusExt = 0
    simEndTime = 0

    if issubset(0, microbesDF[:, :microbial_abundance])
        microbeExt = 1
        mSimEndTime = maximum(microbesDF[(microbesDF.microbial_abundance.!=0), :].t) + 1
        println(microbeExt)
    else
        microbeExt = 0
        mSimEndTime = maximum(microbesDF.t)
        println(microbeExt)
    end

    if issubset(0, virusesDF[:, :viral_abundance])
        if issubset(0, microbesDF[:, :microbial_abundance])
            virusExt = 0
            vSimEndTime = maximum(virusesDF.t)
        else
            virusExt = 1
            vSimEndTime = maximum(virusesDF[(virusesDF.viral_abundance.!=0), :].t) + 1
        end
        println(virusExt)
    else
        virusExt = 0
        vSimEndTime = maximum(virusesDF.t)
        println(virusExt)
    end

    execute(dbOutput, "BEGIN TRANSACTION")
    execute(dbOutput, "INSERT INTO extinction_occurrence VALUES (?,?)", (microbeExt, virusExt))
    execute(dbOutput, "INSERT INTO simulation_end_time VALUES ($(mSimEndTime),$(vSimEndTime))")
    execute(dbOutput, "COMMIT")
end

peaks, walls, peakSeries, peakDurations, wallSeries = peakwallCount(uPercent, lPercent, shannonThreshold, dbOutput, dbSim);
count = DataFrame(num_peaks=peaks, num_walls=walls)

peakSeries |> SQLite.load!(dbOutput, "microbial_peak_series", ifnotexists=true)
peakDurations |> SQLite.load!(dbOutput, "microbial_peakwall_durations", ifnotexists=true)
count |> SQLite.load!(dbOutput, "microbial_peakwall_count", ifnotexists=true)
wallSeries |> SQLite.load!(dbOutput, "microbial_wall_series", ifnotexists=true)

println("Processing extinction occurrences of run $(run_id)")
extinction()
println("Complete!")
