#!/usr/bin/env julia

## This script extracts a very particular part of the full tripartite structure.
# Full match phenotypes for viruses are classified, with the condition that
# they have a susceptible group of hosts
### The full match phenotypes for the viruses needed to be classified so that
### we can use the data on probability of emergence
# Full match phenotypes of micobes are also classified, however with the
# condition that they have at least one single match with a viral strain.
##
println("(Julia compilation delay...)")

using SQLite
using DataFrames
using SQLite.DBInterface: execute
using SQLite: DB

run_id = ARGS[1]

## Define Paths ##
SCRIPT_PATH = abspath(dirname(PROGRAM_FILE))
#
dbSimPath = joinpath(SCRIPT_PATH, "..", "..", "..", "simulation", "sweep_db_gathered.sqlite") # cluster
dbOutputPath = joinpath("shannon_output.sqlite") # cluster

# dbSimPath = joinpath("/Volumes/Yadgah","crispr-sweep-7-2-2022/isolates/runID3297-c66-r47/runID3297-c66-r47.sqlite") # local
# dbMatchPath = joinpath("/Volumes/Yadgah","crispr-sweep-7-2-2022/isolates/runID3297-c66-r47/matches_output.sqlite") # local
# dbOutputPath = joinpath("/Volumes/Yadgah/crispr-sweep-7-2-2022/isolates/runID3297-c66-r47/tripartite-networks_output.sqlite") # local

if isfile(dbOutputPath)
    error("shannon_output.sqlite already exists; delete first")
end
##
dbOutput = SQLite.DB(dbOutputPath)
dbTempSim = SQLite.DB()
# dbTempSim = SQLite.DB(dbSimPath)
execute(dbOutput, "CREATE TABLE shannon_diversity (t REAL, vshannon REAL, bshannon REAL)")
dbTempSim = SQLite.DB()
execute(dbTempSim, "CREATE TABLE summary (t REAL, microbial_abundance INTEGER, viral_abundance INTEGER)")
execute(dbTempSim, "CREATE TABLE babundance (t REAL, bstrain_id INTEGER, abundance INTEGER)")
execute(dbTempSim, "CREATE TABLE vabundance (t REAL, vstrain_id INTEGER, abundance INTEGER)")

execute(dbTempSim, "BEGIN TRANSACTION")
execute(dbTempSim, "ATTACH DATABASE '$(dbSimPath)' as dbSim")
execute(dbTempSim, "INSERT INTO summary(t, microbial_abundance,viral_abundance) SELECT t, microbial_abundance,viral_abundance FROM dbSim.summary WHERE run_id = $(run_id);")
execute(dbTempSim, "INSERT INTO babundance (t, bstrain_id, abundance) SELECT t, bstrain_id, abundance FROM dbSim.babundance WHERE run_id = $(run_id);")
execute(dbTempSim, "INSERT INTO vabundance (t, vstrain_id, abundance) SELECT t, vstrain_id, abundance FROM dbSim.vabundance WHERE run_id = $(run_id);")
execute(dbTempSim, "COMMIT")

execute(dbTempSim, "BEGIN TRANSACTION")
execute(dbTempSim, "CREATE INDEX summary_index ON summary (t,microbial_abundance,viral_abundance)")
execute(dbTempSim, "CREATE INDEX bstrain_index ON babundance (t,bstrain_id)")
execute(dbTempSim, "CREATE INDEX vstrain_index ON vabundance (t,vstrain_id)")
execute(dbTempSim, "COMMIT")

println("Processing shannon entropies of run $(run_id)")
function shannon(dbTempSim, dbOutput)
    for (time,) in execute(dbTempSim, "SELECT DISTINCT t FROM summary")
        #time = time + 400 # for testing
        println("Computing shannon entropy at time $(time)")

        (totAbunds,) = execute(dbTempSim, "SELECT microbial_abundance,viral_abundance FROM summary WHERE t = ?", (time,))
        btotal = totAbunds.microbial_abundance
        vtotal = totAbunds.viral_abundance

        bsubAbunds = [bsubAbunds.abundance for bsubAbunds in execute(dbTempSim, "SELECT abundance FROM babundance WHERE t = ?", (time,))]
        brelAbunds = bsubAbunds ./ btotal
        brelAbunds = brelAbunds[brelAbunds.>=1e-200]

        if btotal == 0
            bshannon = 0
        else
            bshannon = exp(-1 * sum(brelAbunds .* (log.(brelAbunds))))
        end

        vsubAbunds = [vsubAbunds.abundance for vsubAbunds in execute(dbTempSim, "SELECT abundance FROM vabundance WHERE t = ?", (time,))]
        vrelAbunds = vsubAbunds ./ vtotal
        vrelAbunds = vrelAbunds[vrelAbunds.>=1e-200]

        if vtotal == 0
            vshannon = 0
        else
            vshannon = exp(-1 * sum(vrelAbunds .* (log.(vrelAbunds))))
        end

        execute(dbOutput, "BEGIN TRANSACTION")
        execute(dbOutput, "INSERT INTO shannon_diversity VALUES (?,?,?)", (time, vshannon, bshannon))
        execute(dbOutput, "COMMIT")
    end
end

shannon(dbTempSim, dbOutput)

println("Complete!")

