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
dbOutputPath = joinpath("walls-diversity_output.sqlite") # cluster
# dbSimPath = joinpath("/Volumes/Yadgah","crispr-sweep-7-2-2022/isolates/runID3297-c66-r47/runID3297-c66-r47.sqlite") # local
# dbTriPath = joinpath("/Volumes/Yadgah/comboID66/tripartite-networksC66.sqlite") # local
# dbShanPath = joinpath("/Volumes/Yadgah/comboID66/shannonC66.sqlite") # local
# dbMatchPath = joinpath("/Volumes/Yadgah/comboID66/matchesC66.sqlite") # local
# dbOutputPath = joinpath("/Volumes/Yadgah/comboID66/walls-diversity.sqlite") # local

if isfile(dbOutputPath)
    error("walls-diversity_output.sqlite already exists; delete first")
end
##
dbTempSim = SQLite.DB(dbSimPath)
(comboID,) = execute(dbTempSim, "SELECT combo_id FROM runs WHERE run_id = $(run_id)")
comboID = comboID.combo_id
(cCapacity,) = execute(dbTempSim, "SELECT microbe_carrying_capacity FROM param_combos WHERE combo_id = $(comboID)")
cCapacity = cCapacity.microbe_carrying_capacity
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
    if isfile(joinpath(SCRIPT_PATH, "..", "..", "..", "data-analysis",
        "gathered-analyses", "tripartite-networks", "tripartite-networks.sqlite"))
        dbMatchPath = joinpath(SCRIPT_PATH, "..", "..", "..", "data-analysis",
            "gathered-analyses", "matches", "matches.sqlite")
        dbTriPath = joinpath(SCRIPT_PATH, "..", "..", "..", "data-analysis",
            "gathered-analyses", "tripartite-networks", "tripartite-networks.sqlite")
    elseif isfile(joinpath(SCRIPT_PATH, "..", "..", "isolates",
        "comboID$(cr[1])", "tripartite-networksC$(cr[1]).sqlite"))
        dbMatchPath = joinpath(SCRIPT_PATH, "..", "..", "isolates",
            "comboID$(cr[1])", "matchesC$(cr[1]).sqlite")
        dbTriPath = joinpath(SCRIPT_PATH, "..", "..", "isolates",
            "comboID$(cr[1])", "tripartite-networksC$(cr[1]).sqlite")
    else
        error("Neither /data-analysis/gathered-analyses/tripartite-networks/tripartite-networks.sqlite exists \n
            nor does /isolated-runs/isolates/comboID$(cr[1])/tripartite-networksC$(cr[1]).sqlite; \n
            compute one of them first. \n
            Note that `matches` analysis needs to be run before `tripartite-networks` analysis.")
    end
end
if length(ARGS) == 1
    dbShanPath = joinpath(SCRIPT_PATH, "..", "..", "isolates",
        "runID$(run_id)-c$(cr.combo_id)-r$(cr.replicate)", "shannon_output.sqlite")
    if !isfile(dbShanPath)
        error("shannon_output.sqlite does not exist; compute first")
    end
    dbMatchPath = joinpath(SCRIPT_PATH, "..", "..", "isolates",
        "runID$(run_id)-c$(cr.combo_id)-r$(cr.replicate)", "matches_output.sqlite")
    dbTriPath = joinpath(SCRIPT_PATH, "..", "..", "isolates",
        "runID$(run_id)-c$(cr.combo_id)-r$(cr.replicate)", "tripartite-networks_output.sqlite")
    if !isfile(dbMatchPath)
        error("tripartite-networks_output.sqlite does not exist; compute first. \n 
        Note that `matches` analysis needs to be run before `tripartite-networks` analysis.")
    end
end
dbTempMatch = SQLite.DB(dbMatchPath)
dbTempShan = SQLite.DB(dbShanPath)
dbTempTri = SQLite.DB(dbTriPath)
dbOutput = SQLite.DB(dbOutputPath)
execute(dbOutput, "CREATE TABLE extinction_occurrence (microbes INTEGER, viruses INTEGER)")
execute(dbOutput, "CREATE TABLE simulation_end_time (microbe_end_time REAL, virus_end_time REAL)")
# dbTempSim = SQLite.DB(dbSimPath) # local
dbTempSim = SQLite.DB()
execute(dbTempSim, "CREATE TABLE summary (t REAL, microbial_abundance INTEGER, viral_abundance INTEGER)")
execute(dbTempSim, "CREATE TABLE babundance (t REAL, bstrain_id INTEGER, abundance INTEGER)")
execute(dbTempSim, "BEGIN TRANSACTION")
execute(dbTempSim, "ATTACH DATABASE '$(dbSimPath)' as dbSim")
execute(dbTempSim,"INSERT INTO summary (t, microbial_abundance, viral_abundance) 
                    SELECT t, microbial_abundance, viral_abundance 
                    FROM dbSim.summary WHERE run_id = $(run_id);")
execute(dbTempSim,"INSERT INTO babundance (t, bstrain_id, abundance) 
                    SELECT t, bstrain_id, abundance 
                    FROM dbSim.babundance WHERE run_id = $(run_id);")
execute(dbTempSim, "COMMIT")
execute(dbTempSim, "BEGIN TRANSACTION")
execute(dbTempSim, "CREATE INDEX summary_index ON summary (t)")
execute(dbTempSim, "COMMIT")

if length(ARGS) > 1 && ARGS[end] == "combo"
    include(joinpath(SCRIPT_PATH, "combo-functions.jl"))
else
    include(joinpath(SCRIPT_PATH, "functions.jl"))
end

function extinction()
    microbesDF = DataFrame(execute(dbTempSim, "SELECT t,microbial_abundance FROM summary ORDER BY t"))
    virusesDF = DataFrame(execute(dbTempSim, "SELECT t,viral_abundance FROM summary ORDER by t"))
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

peaks, walls, peakSeries, peakDurations, wallSeries = peakwallCount(uPercent, lPercent, shannonThreshold, dbOutput, dbTempSim);
count = DataFrame(num_peaks=peaks, num_walls=walls)

peakSeries |> SQLite.load!(dbOutput, "microbial_peak_series", ifnotexists=true)
peakDurations |> SQLite.load!(dbOutput, "microbial_peakwall_durations", ifnotexists=true)
count |> SQLite.load!(dbOutput, "microbial_peakwall_count", ifnotexists=true)
wallSeries |> SQLite.load!(dbOutput, "microbial_wall_series", ifnotexists=true)

println("Processing extinction occurrences of run $(run_id)")
extinction()
diversityToExtinction()
println("Complete!")
