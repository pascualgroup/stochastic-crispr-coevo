#!/usr/bin/env julia

println("(Annoying Julia compilation delay...)")

using SQLite
using DataFrames
using SQLite.DBInterface: execute


## Define Paths ##
SCRIPT_PATH = abspath(dirname(PROGRAM_FILE)) # cluster

#DBOutputPath = joinpath(SCRIPT_PATH,"..","gathered-analyses","phases","evophases.sqlite") # cluster
dbOutputPath = joinpath("/Volumes/Yadgah/evophases.sqlite") # local

if isfile(dbOutputPath)
    error("evophases.sqlite already exists; please move or delete")
end

#mkpath(joinpath(SCRIPT_PATH,"..","gathered-analyses","phases")) # cluster

# I don't think I need this; maybe delete
#dbSimPath = joinpath(SCRIPT_PATH,"..","..","simulation","sweep_db_gathered.sqlite") # cluster
#dbSimPath = joinpath("/Volumes/Yadgah/sweep_db_gathered.sqlite") # local

#dbSimInfoPath = joinpath(SCRIPT_PATH,"..","..","simulation","sweep_db.sqlite") # cluster
dbSimInfoPath = joinpath("/Volumes/Yadgah/sweep_db.sqlite") # local

#dbWallPath = joinpath(SCRIPT_PATH,"..","gathered-analyses","walls","walls.sqlite") # cluster
dbWallPath = joinpath("/Volumes/Yadgah/walls.sqlite") # local

#dbExtinctPath = joinpath(SCRIPT_PATH,"..","gathered-analyses","extinctions","extinctions.sqlite") # cluster # Activate once extinction.jl script is made
#dbExtinctPath = joinpath("/Volumes/Yadgah/extinctions.sqlite") # local # Activate once extinction.jl script is made
##

#dbExtinct = SQLite.DB(dbExtinctPath) # Activate once extinction.jl script is made
dbOutput = SQLite.DB(dbOutputPath)
dbSimInfo = SQLite.DB(dbSimInfoPath)
dbWall = SQLite.DB(dbWallPath)

timeStmt = map(x->string(" OR t = $(x)"), peakSeries.t[(j1+1):j4])
timeStmt = append!([string("t = $(peakSeries.t[j1])")],timeStmt)
bhill2 = [hill.bhill2 for hill in execute(dbTemp, "SELECT bhill2 FROM hill_no2 WHERE $(timeStmt...) ORDER BY t")]
vhill2 = [hill.vhill2 for hill in execute(dbTemp, "SELECT vhill2 FROM hill_no2 WHERE $(timeStmt...) ORDER BY t")]



execute(dbOutput, "CREATE TABLE average_wall_presence (combo_id INTEGER, wall_presence REAL)")
execute(dbOutput, "CREATE TABLE average_num_walls (combo_id INTEGER, num_walls REAL)")
execute(dbOutput, "CREATE TABLE average_wall_duration (combo_id INTEGER, duration REAL)")
#execute(dbOutput, "CREATE TABLE average_extinction_presence (combo_id INTEGER, viral_extinction_presence REAL, microbial_extinction_presence REAL)")


# YOU ARE HERE!!
for (combo_id,) in execute(dbSimInfo,"SELECT combo_id FROM param_combos")
    (runrep,) = execute(dbSimInfo, "SELECT DISTINCT run_id,replicate FROM runs WHERE combo_id = ? ORDER BY replicate", (combo_id,))

    execute(dbWall, "SELECT FROM summary WHERE t = ?", (time,))

end


(CARRYING_CAP,) = execute(dbSimInfo,"SELECT microbe_carrying_capacity FROM param_combos WHERE combo_id = ?",(combo_id,))
CARRYING_CAP = CARRYING_CAP.microbe_carrying_capacity

execute(dbOutput, "CREATE TABLE richness (t REAL, vrichness INTEGER, brichness INTEGER)")

# Create temporary database that is a copy of the main database at the run_id value of the script's argument
#dbTemp = SQLite.DB("/Volumes/Yadgah/timeSeries$(run_id).sqlite") # local

dbTemp = SQLite.DB()
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
