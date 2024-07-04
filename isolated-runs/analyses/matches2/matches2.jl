#!/usr/bin/env julia

println("(Julia compilation delay...)")

using SQLite
using DataFrames
using SQLite.DBInterface: execute
# using BenchmarkTools


run_id = ARGS[1]
tInc = parse(Float64,ARGS[2])

## Define Paths ##
SCRIPT_PATH = abspath(dirname(PROGRAM_FILE))

dbSimPath = joinpath(SCRIPT_PATH,"..","..","..","simulation","sweep_db_gathered.sqlite") # cluster
dbOutputPath = joinpath("matches_output.sqlite") # cluster

# dbSimPath = joinpath("/Volumes/Yadgah/sweep-6-5-2024/sweep_db_gathered.sqlite") # local
# dbOutputPath = joinpath("/Volumes/Yadgah/matches_output.sqlite") # local

if isfile(dbOutputPath)
    rm(dbOutputPath, force=true)
    # error("matches_output.sqlite already exists; delete first")
end # cluster
##

dbSim = SQLite.DB(dbSimPath)
dbOutput = SQLite.DB(dbOutputPath)

execute(dbOutput, "CREATE TABLE bstrain_to_vstrain_matches (t REAL, bstrain_id INTEGER, vstrain_id INTEGER,
bstrain_specific_match_id INTEGER, match_length INTEGER)")
execute(dbOutput, "CREATE TABLE matches_spacers (bstrain_id INTEGER, bstrain_specific_match_id INTEGER, match_type INTEGER, spacer_id INTEGER)")
execute(dbOutput, "CREATE TABLE bstrain_to_vstrain_0matches (t REAL, bstrain_id INTEGER, vstrain_id INTEGER,
bstrain_specific_0match_id INTEGER)")

# Create temporary database that is a copy of the main database at the run_id value of the script's argument
# tmpPath = joinpath(ENV["SLURM_TMPDIR"], "matches-$(run_id).sqlite")
# tmpPath = joinpath("/Volumes/Yadgah/test.sqlite") # local
println("this is the temp path: $(tmpPath)")
# rm(tmpPath, force=true)
# dbTemp = SQLite.DB(tmpPath)
dbTemp = SQLite.DB()
execute(dbTemp, "CREATE TABLE babundance (t REAL, bstrain_id, abundance INTEGER)")
execute(dbTemp, "CREATE TABLE vabundance (t REAL, vstrain_id INTEGER)")
execute(dbTemp, "CREATE TABLE bspacers (bstrain_id INTEGER, spacer_id INTEGER)")
execute(dbTemp, "CREATE TABLE vpspacers (vstrain_id INTEGER, spacer_id INTEGER)")


execute(dbTemp, "BEGIN TRANSACTION")
execute(dbTemp,"ATTACH DATABASE '$(dbSimPath)' as dbSim")
execute(dbTemp,"INSERT INTO babundance (t, bstrain_id, abundance) SELECT t, bstrain_id, abundance FROM dbSim.babundance WHERE run_id = $(run_id);")
execute(dbTemp,"INSERT INTO bspacers (bstrain_id,spacer_id) SELECT bstrain_id,spacer_id FROM dbSim.bspacers WHERE run_id = $(run_id);")
execute(dbTemp, "COMMIT")


execute(dbTemp, "BEGIN TRANSACTION")
execute(dbTemp, "CREATE INDEX babundance_index ON babundance (t,bstrain_id,abundance)")
execute(dbTemp, "CREATE INDEX bspacers_index ON bspacers (bstrain_id,spacer_id)")
execute(dbTemp, "COMMIT")


println("Identifying matches of run $(run_id)")
function identifyMatches()
    times = [t for (t,) in execute(dbTemp, "SELECT DISTINCT t FROM babundance")]
    tMax = maximum(times)
    times = collect(0:tInc:tMax)
    vstrainsDF = DataFrame(execute(dbSim,"SELECT DISTINCT t, vstrain_id FROM vabundance 
                                    WHERE run_id = $(run_id) 
                                    AND t in ($(join(times,", ")))"))
    # vpspacers = Dict((vstrain_id, 
    #                         [spacer_id for (spacer_id,) in
    #                             execute(dbTemp, "SELECT spacer_id FROM 
    #                                             vpspacers 
    #                                             WHERE vstrain_id = $(vstrain_id)")]
    #                         ) 
    #                         for vstrain_id in vstrainsDF[:,:vstrain_id])
    vpspacersDF = DataFrame(execute(dbSim,"SELECT DISTINCT vstrain_id, spacer_id FROM vpspacers 
                                WHERE run_id = $(run_id)
                                AND vstrain_id in ($(join(unique(vstrainsDF.vstrain_id),", ")))"))
    println("last bstrainID is $(maximum([bstrain_id for (bstrain_id,) in execute(dbTemp, "SELECT DISTINCT bstrain_id FROM babundance
                            WHERE t in ($(join(times,", "))) ORDER BY bstrain_id")]))")
    for (bstrain_id,) in execute(dbTemp, "SELECT DISTINCT bstrain_id FROM babundance
                            WHERE t in ($(join(times,", "))) ORDER BY bstrain_id")
        println("bstrainID: $(bstrain_id)")
        strainTimes = [t for (t,) in execute(dbTemp, "SELECT DISTINCT t FROM babundance
                            WHERE t in ($(join(times,", "))) 
                            AND bstrain_id = $(bstrain_id) 
                            AND abundance > 0 
                            ORDER BY t")]
        strainTimes = strainTimes[strainTimes .<= maximum(vstrainsDF.t)]
        if length(strainTimes) == 0
            continue
        end
        vstrains = copy(vstrainsDF[map(x->in(x,strainTimes),vstrainsDF.t),:])
        vstrains = innerjoin(vstrains,vpspacersDF,on=[:vstrain_id])
        spacers = [spacer_id for (spacer_id,) in
            execute(dbTemp, "SELECT spacer_id FROM bspacers WHERE bstrain_id = $(bstrain_id)")]
        vmatches = vstrains[map(x -> in(x, spacers), vstrains.spacer_id), :]
        vmatches = unique(vmatches[:, :vstrain_id])
        if length(vmatches) == 0
            vstrains = unique(select(vstrains,Not(:spacer_id)))
            vstrains[!,:bstrain_id] = repeat([bstrain_id],length(vstrains[:,:vstrain_id]))
            vstrains[!,:bstrain_specific_0match_id] = copy(vstrains[:,:vstrain_id])
            ordered = sort(unique(vstrains[:,:vstrain_id]))
            matchIDs = Dict((ordered[i],i) for i in  collect(1:length(ordered)))
            replace!(vstrains.bstrain_specific_0match_id,matchIDs...)
            vstrains[!,[:t,:bstrain_id,:vstrain_id,:bstrain_specific_0match_id]] |> SQLite.load!(dbOutput, "bstrain_to_vstrain_0matches", ifnotexists=true)
            continue
        end
        #
        vstrains = vstrains[map(x -> in(x, spacers), vstrains.spacer_id), :]
        network = DataFrame(vstrain_id=vmatches, bstrain_id=repeat([bstrain_id],length(vmatches)))
        network[!, :bstrain_specific_match_id] = copy(network[:,:vstrain_id])
        ordered = sort(unique(network[:,:vstrain_id]))
        matchIDs = Dict((ordered[i],i) for i in  collect(1:length(ordered)))
        replace!(network.bstrain_specific_match_id,matchIDs...)
        network[!, :match_length] = copy(network[:, :vstrain_id])
        match_lengths = Dict((vstrain_id,
                length(intersect(spacers,unique(vstrains[vstrains.vstrain_id.==vstrain_id,:spacer_id]))))
                for vstrain_id in vmatches)     
        replace!(network.match_length,match_lengths...)
        network = innerjoin(network,unique(select(vstrains,Not(:spacer_id))),on=[:vstrain_id])
        network[!,[:t,:bstrain_id,:vstrain_id,:bstrain_specific_match_id,:match_length]] |> 
            SQLite.load!(dbOutput, "bstrain_to_vstrain_matches", ifnotexists=true)
        network = innerjoin(unique(network[:, [:bstrain_id, :vstrain_id, :match_length, :bstrain_specific_match_id]]),
                                unique(vstrains[:,[:vstrain_id,:spacer_id]]),on=[:vstrain_id])
        rename!(network, [:match_length => :match_type])
        network[:, [:bstrain_id, :bstrain_specific_match_id, :match_type, :spacer_id]] |>
            SQLite.load!(dbOutput, "matches_spacers", ifnotexists=true)
        #
        vstrains = copy(vstrainsDF[map(x -> in(x, strainTimes), vstrainsDF.t), :])
        vstrains = vstrains[map(x -> !in(x, vmatches), vstrains.vstrain_id), :]
        if length(vstrains[:,:vstrain_id]) > 0 
            vstrains[!, :bstrain_id] = repeat([bstrain_id], length(vstrains[:, :vstrain_id]))
            vstrains[!,:bstrain_specific_0match_id] = copy(vstrains[:,:vstrain_id])
            ordered = sort(unique(vstrains[:,:vstrain_id]))
            matchIDs = Dict((ordered[i],i) for i in  collect(1:length(ordered)))
            replace!(vstrains.bstrain_specific_0match_id,matchIDs...)
            vstrains[!, [:t, :bstrain_id, :vstrain_id, :bstrain_specific_0match_id]] |> SQLite.load!(dbOutput, "bstrain_to_vstrain_0matches", ifnotexists=true)
            continue
        end
    end
end


identifyMatches()

function createindices()
    println("(Creating run_id indices...)")
    db = SQLite.DB(dbOutputPath)
    execute(db, "BEGIN TRANSACTION")
    for (table_name,) in execute(
        db, "SELECT name FROM sqlite_schema
        WHERE type='table' ORDER BY name;")
        # cols = [info.name for info in execute(db,"PRAGMA table_info($(table_name))")]
        if in(table_name,["bstrain_to_vstrain_matches","bstrain_to_vstrain_0matches"])
            execute(db, "CREATE INDEX $(table_name)_index ON $(table_name) (t, bstrain_specific_match_id, bstrain_id)")
        end
        if in(table_name,["matches_spacers"])
            execute(db, "CREATE INDEX $(table_name)_index ON $(table_name) (t, bstrain_specific_match_id, bstrain_id, match_type)")
        end
    end
    execute(db, "COMMIT")
end
createindices()

println("Complete!")
