#!/usr/bin/env julia

println("(Annoying Julia compilation delay...)")

using SQLite
using DataFrames
using SQLite.DBInterface: execute


run_id = ARGS[1]

## Define Paths ##
SCRIPT_PATH = abspath(dirname(PROGRAM_FILE))

dbSimPath = joinpath(SCRIPT_PATH,"..","..","simulation","sweep_db_gathered.sqlite") # cluster
#dbSimPath = joinpath("/Volumes/Yadgah/sweep_db_gathered.sqlite") # local

dbSimInfoPath = joinpath(SCRIPT_PATH,"..","..","simulation","sweep_db.sqlite") # cluster
#dbSimInfoPath = joinpath("/Volumes/Yadgah/sweep_db.sqlite") # local

if isfile("shannon_output.sqlite")
    error("shannon_output.sqlite already exists; delete first")
end # cluster
dbOutput = SQLite.DB("shannon_output.sqlite") # cluster
#dbOutput = SQLite.DB("/Volumes/Yadgah/shannon_output.sqlite") # local
##

dbSimInfo = SQLite.DB(dbSimInfoPath)
(combo_id,) = execute(dbSimInfo,"SELECT combo_id FROM runs WHERE run_id = ?",(run_id,))
combo_id = combo_id.combo_id
(CARRYING_CAP,) = execute(dbSimInfo,"SELECT microbe_carrying_capacity FROM param_combos WHERE combo_id = ?",(combo_id,))
CARRYING_CAP = CARRYING_CAP.microbe_carrying_capacity

execute(dbOutput, "CREATE TABLE hill_no2 (t REAL, vhill2 REAL, bhill2 REAL)")

# Create temporary database that is a copy of the main database at the run_id value of the script's argument
#dbTemp = SQLite.DB("/Volumes/Yadgah/timeSeries$(run_id).sqlite") # local
dbTemp = SQLite.DB()
execute(dbTemp, "CREATE TABLE summary (t REAL, microbial_abundance INTEGER, viral_abundance INTEGER)")
execute(dbTemp, "CREATE TABLE babundance (t REAL, bstrain_id INTEGER, abundance INTEGER)")
execute(dbTemp, "CREATE TABLE vabundance (t REAL, vstrain_id INTEGER, abundance INTEGER)")

execute(dbTemp, "BEGIN TRANSACTION")
execute(dbTemp,"ATTACH DATABASE '$(dbSimPath)' as dbSim")
execute(dbTemp,"INSERT INTO summary(t, microbial_abundance,viral_abundance) SELECT t, microbial_abundance,viral_abundance FROM dbSim.summary WHERE run_id = $(run_id);")
execute(dbTemp,"INSERT INTO babundance (t, bstrain_id, abundance) SELECT t, bstrain_id, abundance FROM dbSim.babundance WHERE run_id = $(run_id);")
execute(dbTemp,"INSERT INTO vabundance (t, vstrain_id, abundance) SELECT t, vstrain_id, abundance FROM dbSim.vabundance WHERE run_id = $(run_id);")
execute(dbTemp, "COMMIT")

execute(dbTemp, "BEGIN TRANSACTION")
execute(dbTemp, "CREATE INDEX summary_index ON summary (t,microbial_abundance,viral_abundance)")
execute(dbTemp, "CREATE INDEX bstrain_index ON babundance (t,bstrain_id)")
execute(dbTemp, "CREATE INDEX vstrain_index ON vabundance (t,vstrain_id)")
execute(dbTemp, "COMMIT")

println("Processing shannon entropies of run $(run_id)")
function shannon(dbTemp,dbOutput)
    for (time,) in execute(dbTemp, "SELECT DISTINCT t FROM summary")
        #time = time + 400 # for testing
        println("Computing shannon entropy at time $(time)")

        (totAbunds,) = execute(dbTemp, "SELECT microbial_abundance,viral_abundance FROM summary WHERE t = ?", (time,))
        #THE EXPRESSION ABOVE DOESN'T SAVE AS A NUMBER,
        #IT ONLY SAVES AS A NUMBER WHEN A LOOP ARGUMENT
        btotal = totAbunds.microbial_abundance
        vtotal = totAbunds.viral_abundance

        bsubAbunds = [bsubAbunds.abundance for bsubAbunds in execute(dbTemp, "SELECT abundance FROM babundance WHERE t = ?", (time,))]
        brelAbunds = bsubAbunds./btotal
        brelAbunds = brelAbunds[brelAbunds.>=1e-200]

        # This is assuming diversity will not reach +Inf
        #if true in isnan.(brelAbunds.*(log.(brelAbunds)))
            #bhill2 = 0
        #else
            bhill2 = exp(-1*sum(brelAbunds.*(log.(brelAbunds))))
        #end

        vsubAbunds = [vsubAbunds.abundance for vsubAbunds in execute(dbTemp, "SELECT abundance FROM vabundance WHERE t = ?", (time,))]
        vrelAbunds = vsubAbunds./vtotal
        vrelAbunds = vrelAbunds[vrelAbunds.>=1e-200]

        # This is assuming diversity will not reach +Inf
        #if true in isnan.(vrelAbunds.*(log.(vrelAbunds)))
            #vhill2 = 0
        #else
            vhill2 = exp(-1*sum(vrelAbunds.*(log.(vrelAbunds))))
        #end

        execute(dbOutput, "BEGIN TRANSACTION")
        execute(dbOutput, "INSERT INTO hill_no2 VALUES (?,?,?)", (time,vhill2,bhill2))
        execute(dbOutput, "COMMIT")

        if vtotal == 0 && btotal >= CARRYING_CAP
            println("Virus has gone extinct! Also, microbes are at carrying capacity. Truncating data and projecting rest of computation...")

            timeData = [summary.t for summary in execute(dbTemp, "SELECT DISTINCT t FROM summary")]
            timeDataTrunc = timeData[timeData.>time]
            timesLeft = length(timeDataTrunc)

            bsubAbunds = [bsubAbunds.abundance for bsubAbunds in execute(dbTemp, "SELECT abundance FROM babundance WHERE t = ?", (time,))]
            brelAbunds = bsubAbunds./btotal
            brelAbunds = brelAbunds[brelAbunds.>=1e-200]

            # This is assuming diversity will not reach +Inf
            #if true in isnan.(brelAbunds.*(log.(brelAbunds)))
                #bhill2 = 0
            #else
                bhill2 = exp(-1*sum(brelAbunds.*(log.(brelAbunds))))
            #end


            tempDF = DataFrame(t = Array{Float64,1}(timeDataTrunc),vhill2 = Array{Float64,1}(zeros(timesLeft)),bhill2 = Array{Float64,1}(bhill2*ones(timesLeft)))
            tempDF |> SQLite.load!(dbOutput,"hill_no2",ifnotexists=true)

            return
        end
    end
end

shannon(dbTemp,dbOutput)

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
