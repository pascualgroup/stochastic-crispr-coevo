#!/usr/bin/env julia

using SQLite
using DataFrames
using SQLite.DBInterface: execute

dbAnalysis = SQLite.DB(joinpath("..","data_analysis.sqlite"))
dbSim = SQLite.DB(joinpath("..","..","simulation","sweep_db_gathered.sqlite"))

run_id = ARGS[1]

println("Processing: Richness of Run $(run_id)")
for (time,) in execute(dbSim, "SELECT DISTINCT t FROM summary WHERE run_id =?",(run_id,))

    #bstrains = [strain.bstrain_id for strain in execute(db, "SELECT DISTINCT bstrain_id FROM babundance WHERE t = ?", (time,))]
    #brichness = sum(bstrains)
    bstrains = DataFrame(execute(dbSim, "SELECT DISTINCT bstrain_id FROM babundance WHERE t = ?", (time,)))
    brichness = sum(bstrains.bstrain_id)
    #println(brichness)

    #vstrains = [strain.vstrain_id for strain in execute(db, "SELECT DISTINCT vstrain_id FROM vabundance WHERE t = ?", (time,))]
    #vrichness = sum(vstrains)
    vstrains = DataFrame(execute(dbSim, "SELECT DISTINCT vstrain_id FROM vabundance WHERE t = ?", (time,)))
    vrichness = sum(vstrains.vstrain_id)
    #println([brichness,vrichness])

    execute(dbAnalysis, "BEGIN TRANSACTION")
    execute(dbAnalysis, "INSERT INTO richness VALUES (?,?,?,?)", (run_id,time,vrichness,brichness))
    execute(dbAnalysis, "COMMIT")
    #return
end
println("Complete!")
