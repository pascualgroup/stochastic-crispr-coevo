#!/usr/bin/env julia

println("(Julia compilation delay...)")

using SQLite
using DataFrames
using SQLite.DBInterface: execute

run_id = ARGS[1]

## Define Paths ##
SCRIPT_PATH = abspath(dirname(PROGRAM_FILE))

dbSimPath = joinpath(SCRIPT_PATH,"..","..","simulation","sweep_db_gathered.sqlite") # cluster
dbMatchPath = #######
dbOutputPath = joinpath("match-diversity_output.sqlite") # cluster

# # dbSimPath = joinpath("/Volumes/Yadgah/sweep_db_gathered.sqlite") # local
# dbSimPath = joinpath("/Volumes/Yadgah","run_id1455_combo73_replicate15.sqlite") # local
# dbOutputPath = joinpath("/Volumes/Yadgah/matches_output.sqlite") # local
# dbSimPath = joinpath("/Volumes/Yadgah/crispr-sweep-7-2-2022/isolated-runs/
# dbMatchPath = #######
# run_id3297_combo66_replicate47/run_id3297_combo66_replicate47.sqlite")

if isfile(dbOutputPath)
    error("match-diversity_output.sqlite already exists; delete first")
end # cluster
##
