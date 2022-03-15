#!/usr/bin/env julia

println("(Annoying Julia compilation delay...)")

using SQLite
using DataFrames
using SQLite.DBInterface: execute


run_id = ARGS[1]

## Define Paths ##
SCRIPT_PATH = abspath(dirname(PROGRAM_FILE))

dbSimPath = joinpath(SCRIPT_PATH,"..","..","..","simulation","sweep_db_gathered.sqlite") # cluster
dbOutputPath = joinpath("clade-abundances_output.sqlite") # cluster

# dbSimPath = joinpath("/Volumes/Yadgah/sweep_db_gathered.sqlite") # local
# dbSimPath = joinpath("/Volumes/Yadgah/crispr-sweep-19-1-2022/simulation/sweep_db_gathered.sqlite") # local
# dbSimPath = joinpath("/Volumes/Yadgah","run_id1455_combo73_replicate15.sqlite") # local
# dbSimPath = joinpath("/Volumes/Yadgah","run_id1343_combo68_replicate3.sqlite") # local
# dbOutputPath = joinpath("/Volumes/Yadgah/clade-abundances_output.sqlite") # local

##

dbSim = SQLite.DB(dbSimPath)
dbOutput = SQLite.DB(dbOutputPath)

execute(dbOutput, "CREATE TABLE clade_babundances (t REAL, clade_id INTEGER, abundance INTEGER)")
execute(dbOutput, "CREATE TABLE babundances (t REAL, clade_id INTEGER, bstrain_id INTEGER, abundance INTEGER)")
execute(dbOutput, "CREATE TABLE clade_vabundances (t REAL, clade_id INTEGER, abundance INTEGER)")
execute(dbOutput, "CREATE TABLE vabundances (t REAL, clade_id INTEGER, vstrain_id INTEGER, abundance INTEGER)")

# Create temporary database that is a copy of the main database at the run_id value of the script's argument
dbTemp = SQLite.DB()
#dbTemp = SQLite.DB("/Volumes/Yadgah/test.sqlite") # local
execute(dbTemp, "CREATE TABLE babundance (t REAL, bstrain_id INTEGER, abundance INTEGER)")
execute(dbTemp, "CREATE TABLE vabundance (t REAL, vstrain_id INTEGER, abundance INTEGER)")
execute(dbTemp, "CREATE TABLE bstrains (bstrain_id INTEGER, parent_bstrain_id INTEGER, infecting_vstrain_id INTEGER)")
execute(dbTemp, "CREATE TABLE vstrains (vstrain_id INTEGER, parent_vstrain_id INTEGER, infected_bstrain_id INTEGER)")


execute(dbTemp, "BEGIN TRANSACTION")
execute(dbTemp,"ATTACH DATABASE '$(dbSimPath)' as dbSim")
execute(dbTemp,"INSERT INTO babundance (t, bstrain_id, abundance) SELECT t, bstrain_id, abundance FROM dbSim.babundance WHERE run_id = $(run_id);")
execute(dbTemp,"INSERT INTO vabundance (t, vstrain_id, abundance) SELECT t, vstrain_id, abundance FROM dbSim.vabundance WHERE run_id = $(run_id);")
execute(dbTemp,"INSERT INTO bstrains (bstrain_id,parent_bstrain_id,infecting_vstrain_id)
SELECT bstrain_id,parent_bstrain_id,infecting_vstrain_id
FROM dbSim.bstrains WHERE run_id = $(run_id);")
execute(dbTemp,"INSERT INTO vstrains (vstrain_id,parent_vstrain_id,infected_bstrain_id)
SELECT vstrain_id,parent_vstrain_id,infected_bstrain_id
FROM dbSim.vstrains WHERE run_id = $(run_id);")
execute(dbTemp, "COMMIT")


execute(dbTemp, "BEGIN TRANSACTION")
execute(dbTemp, "CREATE INDEX babundance_index ON babundance (t,bstrain_id,abundance)")
execute(dbTemp, "CREATE INDEX vabundance_index ON vabundance (t,vstrain_id)")
execute(dbTemp, "CREATE INDEX bstrains_index ON bstrains (bstrain_id,parent_bstrain_id,infecting_vstrain_id)")
execute(dbTemp, "CREATE INDEX vstrains_index ON vstrains (vstrain_id,parent_vstrain_id,infected_bstrain_id)")
execute(dbTemp, "COMMIT")




function microbeCladeAbundances()
    lineages = Dict{Int64,Array{Int64,1}}()
    for (MRCA_id,) in execute(dbTemp,"SELECT bstrain_id FROM bstrains WHERE parent_bstrain_id = 1 ORDER BY bstrain_id")
        lineages[MRCA_id] = searchLineage(MRCA_id,"microbe")
    end
        lineages[1] = [1]
    for (time,) in execute(dbTemp, "SELECT DISTINCT t FROM babundance ORDER BY t")
        println("Tracking microbial clades at time $(time)")
        # if time > 200.0
        #     return
        # end
        strains = [bstrain_id for (bstrain_id,) in
        execute(dbTemp, "SELECT bstrain_id FROM babundance WHERE t = $(time)")]
        if length(strains) == 1 && strains... == 1
            (abundance,) = execute(dbTemp, "SELECT abundance FROM babundance WHERE t = $(time)
            AND bstrain_id = 1")
            execute(dbOutput, "INSERT INTO clade_babundances VALUES (?,?,?)",(time,1,abundance.abundance))
            execute(dbOutput, "INSERT INTO babundances VALUES (?,?,?,?)",(time,1,1,abundance.abundance))
            continue
        end
        clades = [MRCA_id for (MRCA_id,) in execute(dbTemp,"SELECT bstrain_id FROM bstrains
            WHERE parent_bstrain_id = 1 ORDER BY bstrain_id")]
        for MRCA_id in [1;clades]
            descendents = intersect(lineages[MRCA_id],strains)
            if length(descendents) == 0
                continue
            end
            abundances = [abundance for (abundance,) in
            execute(dbTemp, "SELECT abundance FROM babundance WHERE t = $(time)
            AND bstrain_id in ($(join(descendents,", ")))")]
            execute(dbOutput, "INSERT INTO clade_babundances VALUES (?,?,?)",(time,MRCA_id,sum(abundances)))
            for desc_id in descendents
                (abundance,) = execute(dbTemp, "SELECT abundance FROM babundance WHERE t = $(time)
                AND bstrain_id = $(desc_id)")
                execute(dbOutput, "INSERT INTO babundances VALUES (?,?,?,?)",(time,MRCA_id,desc_id,abundance.abundance))
            end
        end
    end
end

function virusCladeAbundances()
    lineages = Dict{Int64,Array{Int64,1}}()
    for (MRCA_id,) in execute(dbTemp,"SELECT vstrain_id FROM vstrains WHERE parent_vstrain_id = 1 ORDER BY vstrain_id")
        lineages[MRCA_id] = searchLineage(MRCA_id,"virus")
    end
        lineages[1] = [1]
    for (time,) in execute(dbTemp, "SELECT DISTINCT t FROM vabundance ORDER BY t")
        println("Tracking viral clades at time $(time)")
        # if time > 200.0
        #     return
        # end
        strains = [vstrain_id for (vstrain_id,) in
        execute(dbTemp, "SELECT vstrain_id FROM vabundance WHERE t = $(time)")]
        if length(strains) == 1 && strains... == 1
            (abundance,) = execute(dbTemp, "SELECT abundance FROM vabundance WHERE t = $(time)
            AND vstrain_id = 1")
            execute(dbOutput, "INSERT INTO clade_vabundances VALUES (?,?,?)",(time,1,abundance.abundance))
            execute(dbOutput, "INSERT INTO vabundances VALUES (?,?,?,?)",(time,1,1,abundance.abundance))
            continue
        end
        clades = [MRCA_id for (MRCA_id,) in execute(dbTemp,"SELECT vstrain_id
        FROM vstrains WHERE parent_vstrain_id = 1 ORDER BY vstrain_id")]
        for MRCA_id in [1;clades]
            descendents = intersect(lineages[MRCA_id],strains)
            if length(descendents) == 0
                continue
            end
            abundances = [abundance for (abundance,) in
            execute(dbTemp, "SELECT abundance FROM vabundance WHERE t = $(time)
            AND vstrain_id in ($(join(descendents,", ")))")]
            execute(dbOutput, "INSERT INTO clade_vabundances VALUES (?,?,?)",(time,MRCA_id,sum(abundances)))
            for desc_id in descendents
                (abundance,) = execute(dbTemp, "SELECT abundance FROM vabundance WHERE t = $(time)
                AND vstrain_id = $(desc_id)")
                execute(dbOutput, "INSERT INTO vabundances VALUES (?,?,?,?)",(time,MRCA_id,desc_id,abundance.abundance))
            end
        end
    end
end


function searchLineage(MRCA_id::Int64,strainType::String)
    if strainType == "virus"
        type = "v"
    end
    if strainType == "microbe"
        type = "b"
    end
    numDescendents = 1
    lineage = Int64[MRCA_id]
    ancestors = Int64[MRCA_id]
    descendents = Int64[]
    while numDescendents > 0
        for id in ancestors
            desc_ids = [id for (id,) in
            execute(dbTemp,"SELECT $(type)strain_id FROM $(type)strains
            WHERE parent_$(type)strain_id in ($(id)) ORDER BY $(type)strain_id")]
            append!(descendents,desc_ids)
            lineage = append!(lineage,desc_ids)
            unique!(lineage)
        end
        numDescendents = length(descendents)
        ancestors = descendents[:]
        descendents = Int64[]
    end
    return sort(lineage)
end

microbeCladeAbundances()
virusCladeAbundances()

println("Complete!")
