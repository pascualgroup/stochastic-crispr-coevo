#!/usr/bin/env julia

println("(Julia compilation delay...)")

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
dbOutputPath = joinpath("extinctions_output.sqlite") # cluster
#dbOutputPath = joinpath("/Volumes/Yadgah/extinctions_output.sqlite") # local

if isfile(dbOutputPath)
    error("extinctions_output.sqlite already exists; delete first")
end # cluster
dbOutput = SQLite.DB(dbOutputPath)
##

execute(dbOutput, "CREATE TABLE extinction_occurrence (microbes INTEGER, viruses INTEGER)")
execute(dbOutput, "CREATE TABLE simulation_end_time (microbe_end_time REAL, virus_end_time REAL)")

# Create temporary database that is a copy of the main database at the run_id value of the script's argument
#dbTemp = SQLite.DB("/Volumes/Yadgah/timeSeries$(run_id).sqlite") # local
dbTemp = SQLite.DB()
execute(dbTemp, "CREATE TABLE summary (t REAL, microbial_abundance INTEGER, viral_abundance INTEGER)")

execute(dbTemp, "BEGIN TRANSACTION")
execute(dbTemp,"ATTACH DATABASE '$(dbSimPath)' as dbSim")
execute(dbTemp,"INSERT INTO summary(t, microbial_abundance,viral_abundance) SELECT t, microbial_abundance,viral_abundance FROM dbSim.summary WHERE run_id = $(run_id);")
execute(dbTemp, "COMMIT")

execute(dbTemp, "BEGIN TRANSACTION")
execute(dbTemp, "CREATE INDEX summary_index ON summary (t,microbial_abundance,viral_abundance)")
execute(dbTemp, "COMMIT")

println("Processing extinction occurrences of run $(run_id)")
function extinction()
    microbesDF = DataFrame(execute(dbTemp, "SELECT t,microbial_abundance FROM summary ORDER BY t"))
    virusesDF = DataFrame(execute(dbTemp, "SELECT t,viral_abundance FROM summary ORDER by t"))
    microbesDF = microbesDF[(microbesDF.t .!=0),:]
    virusesDF = virusesDF[(virusesDF.t .!=0),:]
    microbeExt = 0
    virusExt = 0
    simEndTime = 0

    if issubset(0,microbesDF[:,:microbial_abundance])
        microbeExt = 1
        mSimEndTime = maximum(microbesDF[(microbesDF.microbial_abundance .!=0),:].t) + 1
        println(microbeExt)
    else
        microbeExt = 0
        mSimEndTime = maximum(microbesDF.t)
        println(microbeExt)
    end

    if issubset(0,virusesDF[:,:viral_abundance])
        if issubset(0,microbesDF[:,:microbial_abundance])
            virusExt = 0
            vSimEndTime = maximum(virusesDF.t)
        else
            virusExt = 1
            vSimEndTime = maximum(virusesDF[(virusesDF.viral_abundance .!=0),:].t) + 1
        end
        println(virusExt)
    else
        virusExt = 0
        vSimEndTime = maximum(virusesDF.t)
        println(virusExt)
    end

    execute(dbOutput, "BEGIN TRANSACTION")
    execute(dbOutput, "INSERT INTO extinction_occurrence VALUES (?,?)", (microbeExt,virusExt))
    execute(dbOutput, "INSERT INTO simulation_end_time VALUES ($(mSimEndTime),$(vSimEndTime))")
    execute(dbOutput, "COMMIT")
    
    dbSim = SQLite.DB(dbSimPath)
    tableNames = [table for (table,) in 
                execute(dbSim,"SELECT name FROM sqlite_master WHERE type='table' 
                AND name in ('bextinctions', 'vextinctions');")]
    if length(tableNames) > 0
        execute(dbOutput, "ATTACH DATABASE '$(dbSimPath)' as dbSim")
        for table in tableNames
            if table == "vextinctions"
                strainID = "vstrain_id"
                execute(dbOutput, "CREATE TABLE $(table) ($(strainID) INTEGER, t_extinction REAL)")
                execute(dbOutput, "BEGIN TRANSACTION")
                execute(dbOutput, "INSERT INTO $(table)($(strainID), t_extinction) SELECT $(strainID), t_extinction FROM dbSim.$(table) WHERE run_id = $(run_id);")
                execute(dbOutput, "COMMIT")
            end
            if table == "bextinctions"
                strainID = "bstrain_id"
                execute(dbOutput, "CREATE TABLE $(table) ($(strainID) INTEGER, t_extinction REAL)")
                execute(dbOutput, "BEGIN TRANSACTION")
                execute(dbOutput, "INSERT INTO $(table)($(strainID), t_extinction) SELECT $(strainID), t_extinction FROM dbSim.$(table) WHERE run_id = $(run_id);")
                execute(dbOutput, "COMMIT")
            end
        end
    end
end

extinction()

println("Complete!")
