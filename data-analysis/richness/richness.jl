#!/usr/bin/env julia

println("(Annoying Julia compilation delay...)")

using SQLite
using DataFrames
using SQLite.DBInterface: execute


run_id = ARGS[1] # cluster & local
#run_id = 1 # local

SCRIPT_PATH = abspath(dirname(PROGRAM_FILE)) # cluster

dbSimPath = joinpath(SCRIPT_PATH,"..","..","simulation","sweep_db_gathered.sqlite") # cluster
#dbSimPath = joinpath("/Volumes/Yadgah/sweep_db_gathered.sqlite") # local machine

dbSim = SQLite.DB(dbSimPath)
(combo_id,) = execute(dbSim,"SELECT combo_id FROM runs WHERE run_id = ?",(run_id,))
combo_id = combo_id.combo_id
(CARRYING_CAP,) = execute(dbSim,"SELECT microbe_carrying_capacity FROM param_combos WHERE combo_id = ?",(combo_id,))
CARRYING_CAP = CARRYING_CAP.microbe_carrying_capacity

if isfile("richness_output.sqlite") # cluster
    error("richness_output.sqlite already exists; delete first") # cluster
end # cluster
dbOutput = SQLite.DB("richness_output.sqlite") # cluster
#dbOutput = SQLite.DB("/Volumes/Yadgah/richness_output.sqlite") # local machine
execute(dbOutput, "CREATE TABLE richness (t REAL, vrichness INTEGER, brichness INTEGER)")

# Create temporary database that is a copy of the main database at the run_id value of the script's argument
#dbTemp = SQLite.DB("/Volumes/Yadgah/timeSeries$(run_id).sqlite") # local
dbTemp = SQLite.DB() # cluster
execute(dbTemp, "CREATE TABLE summary (t REAL, microbial_abundance INTEGER)")
execute(dbTemp, "CREATE TABLE babundance (t REAL, bstrain_id INTEGER, abundance INTEGER)")
execute(dbTemp, "CREATE TABLE vabundance (t REAL, vstrain_id INTEGER, abundance INTEGER)")

execute(dbTemp, "BEGIN TRANSACTION")
execute(dbTemp,"ATTACH DATABASE '$(dbSimPath)' as dbSim")
execute(dbTemp,"INSERT INTO summary(t, microbial_abundance) SELECT t, microbial_abundance FROM dbSim.summary WHERE run_id = $(run_id);")
execute(dbTemp,"INSERT INTO babundance (t, bstrain_id, abundance) SELECT t, bstrain_id, abundance FROM dbSim.babundance WHERE run_id = $(run_id);")
execute(dbTemp,"INSERT INTO vabundance (t, vstrain_id, abundance) SELECT t, vstrain_id, abundance FROM dbSim.vabundance WHERE run_id = $(run_id);")
execute(dbTemp, "COMMIT")

execute(dbTemp, "BEGIN TRANSACTION")
execute(dbTemp, "CREATE INDEX summary_index ON summary (t,microbial_abundance)")
execute(dbTemp, "CREATE INDEX bstrain_index ON babundance (t,bstrain_id)")
execute(dbTemp, "CREATE INDEX vstrain_index ON vabundance (t,vstrain_id)")
execute(dbTemp, "COMMIT")

println("Processing richness of run $(run_id)")
function richness(dbTemp,dbOutput)
    for (time,) in execute(dbTemp, "SELECT DISTINCT t FROM summary")
        #time = time + 400 # for testing
        println("Computing richness at time $(time)")
        bstrains = [strain.bstrain_id for strain in execute(dbTemp, "SELECT DISTINCT bstrain_id FROM babundance WHERE t = ?", (time,))]
        brichness = length(bstrains)

        #bstrains = DataFrame(execute(dbTemp, "SELECT DISTINCT bstrain_id FROM babundance WHERE t = ?", (time,)))
        #brichness = length(bstrains.bstrain_id)

        (babund,) = execute(dbTemp, "SELECT microbial_abundance FROM summary WHERE t = ?", (time,))
        #THE EXPRESSION ABOVE DOESN'T SAVE AS A NUMBER,
        #IT ONLY SAVES AS A NUMBER WHEN A LOOP ARGUMENT
        babund = babund.microbial_abundance
        #println(brichness)

        vstrains = [strain.vstrain_id for strain in execute(dbTemp, "SELECT DISTINCT vstrain_id FROM vabundance WHERE t = ?", (time,))]
        vrichness = length(vstrains)

        #vstrains = DataFrame(execute(dbTemp, "SELECT DISTINCT vstrain_id FROM vabundance WHERE t = ?", (time,)))
        #vrichness = length(vstrains.vstrain_id)
        #println([brichness,vrichness])

        execute(dbOutput, "BEGIN TRANSACTION")
        execute(dbOutput, "INSERT INTO richness VALUES (?,?,?)", (time,vrichness,brichness))
        execute(dbOutput, "COMMIT")

        if vrichness == 0 && babund >= CARRYING_CAP
            println("Virus has gone extinct! Also, microbes are at carrying capacity. Truncating data and projecting rest of computation...")

            timeData = [summary.t for summary in execute(dbTemp, "SELECT DISTINCT t FROM summary")]
            timeDataTrunc = timeData[timeData.>time]
            timesLeft = length(timeDataTrunc)
            tempDF = DataFrame(t = Array{Float64,1}(timeDataTrunc),vrichness = Array{Int64,1}(zeros(timesLeft)),brichness = Array{Int64,1}(brichness*ones(timesLeft)))

            #timeData = DataFrame(execute(dbTemp, "SELECT DISTINCT t FROM summary"))
            #timeDataTrunc = timeData[check.t.>time,:]
            #timesLeft = nrow(timeDataTrunc)
            #tempDF = DataFrame(t = timeDataTrunc[:,"t"],vrichness = Array{Int64,1}(zeros(timesLeft)),brichness = Array{Int64,1}(babund*ones(timesLeft)))

            tempDF |> SQLite.load!(dbOutput,"richness",ifnotexists=true)

            return
        end

    end
end

richness(dbTemp,dbOutput)

println("Complete!")
