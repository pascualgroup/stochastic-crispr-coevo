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

if isfile("simpson_output.sqlite")
    error("simpson_output.sqlite already exists; delete first")
end # cluster
dbOutput = SQLite.DB("simpson_output.sqlite") # cluster
#dbOutput = SQLite.DB("/Volumes/Yadgah/simpson_output.sqlite") # local
##

dbSimInfo = SQLite.DB(dbSimInfoPath)
(combo_id,) = execute(dbSimInfo,"SELECT combo_id FROM runs WHERE run_id = ?",(run_id,))
combo_id = combo_id.combo_id
(CARRYING_CAP,) = execute(dbSimInfo,"SELECT microbe_carrying_capacity FROM param_combos WHERE combo_id = ?",(combo_id,))
CARRYING_CAP = CARRYING_CAP.microbe_carrying_capacity

execute(dbOutput, "CREATE TABLE simpson_diversity (t REAL, vsimpson REAL, bsimpson REAL)")

# Create temporary database that is a copy of the main database at the run_id value of the script's argument
#dbTemp = SQLite.DB("/Volumes/Yadgah/timeSeries$(run_id).sqlite") # local
dbTemp = SQLite.DB() # cluster
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

println("Processing inverse simpson index of run $(run_id)")
function simpson(dbTemp,dbOutput)
    for (time,) in execute(dbTemp, "SELECT DISTINCT t FROM summary")
        #time = time + 400 # for testing
        println("Computing inverse simpson index at time $(time)")

        (totAbunds,) = execute(dbTemp, "SELECT microbial_abundance,viral_abundance FROM summary WHERE t = ?", (time,))

        btotal = totAbunds.microbial_abundance
        vtotal = totAbunds.viral_abundance
        #println(bsimpson)

        bsubAbunds = [bsubAbunds.abundance for bsubAbunds in execute(dbTemp, "SELECT abundance FROM babundance WHERE t = ?", (time,))]
        brelAbunds = bsubAbunds./btotal
        brelAbunds = brelAbunds[brelAbunds.>=1e-200]

        if btotal == 0
            bsimpson = 0
        else
            bsimpson = 1/sum(brelAbunds.*brelAbunds)
        end

        vsubAbunds = [vsubAbunds.abundance for vsubAbunds in execute(dbTemp, "SELECT abundance FROM vabundance WHERE t = ?", (time,))]
        vrelAbunds = vsubAbunds./vtotal
        vrelAbunds = vrelAbunds[vrelAbunds.>=1e-200]

        if vtotal == 0
            vsimpson = 0
        else
            vsimpson = 1/sum(vrelAbunds.*vrelAbunds)
        end

        execute(dbOutput, "BEGIN TRANSACTION")
        execute(dbOutput, "INSERT INTO simpson_diversity VALUES (?,?,?)", (time,vsimpson,bsimpson))
        execute(dbOutput, "COMMIT")
    end
end

simpson(dbTemp,dbOutput)

println("Complete!")
