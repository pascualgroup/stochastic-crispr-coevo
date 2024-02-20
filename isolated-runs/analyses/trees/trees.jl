#!/usr/bin/env julia

println("(Julia compilation delay...)")

using SQLite
using DataFrames
using SQLite.DBInterface: execute


run_id = ARGS[1]

## Define Paths ##
SCRIPT_PATH = abspath(dirname(PROGRAM_FILE))

dbSimPath = joinpath(SCRIPT_PATH,"..","..","..","simulation","sweep_db_gathered.sqlite") # cluster
dbOutputPath = joinpath("trees_output.sqlite") # cluster

# dbSimPath = joinpath("/Volumes/Yadgah/sweep_db_gathered.sqlite") # local
# dbSimPath = joinpath("/Volumes/Yadgah","run_id1455_combo73_replicate15.sqlite") # local
# dbSimPath = joinpath("/Volumes/Yadgah", "sylvain-martin-collab/12_MOI3/isolates/runID4209-c15-r9/runID4209-c15-r9.sqlite") # local 
# dbOutputPath = joinpath("/Volumes/Yadgah/trees_output.sqlite") # local

##

dbOutput = SQLite.DB(dbOutputPath)

execute(dbOutput, "CREATE TABLE tree_babundance (t REAL, tree_bstrain_id INTEGER, abundance INTEGER)")
execute(dbOutput, "CREATE TABLE tree_vabundance (t REAL, tree_vstrain_id INTEGER, abundance INTEGER)")
execute(dbOutput, "CREATE TABLE bstrain_creation_extinction (bstrain_id INTEGER, t_creation REAL, t_extinction REAL)")
execute(dbOutput, "CREATE TABLE tree_bstrain_creation_extinction
(tree_bstrain_id INTEGER, t_creation REAL, t_extinction REAL,
tree_parent_bstrain_id INTEGER, tree_infecting_vstrain_id)")
execute(dbOutput, "CREATE TABLE vstrain_creation_extinction (vstrain_id INTEGER, t_creation REAL, t_extinction REAL)")
execute(dbOutput, "CREATE TABLE tree_vstrain_creation_extinction
(tree_vstrain_id INTEGER, t_creation REAL, t_extinction REAL,
tree_parent_vstrain_id INTEGER, tree_infected_bstrain_id INTEGER)")

# Create temporary database that is a copy of the main database at the run_id value of the script's argument
dbTemp = SQLite.DB()
#dbTemp = SQLite.DB("/Volumes/Yadgah/test.sqlite") # local
# execute(dbTemp, "CREATE TABLE babundance (t REAL, bstrain_id INTEGER, abundance INTEGER)")
# execute(dbTemp, "CREATE TABLE vabundance (t REAL, vstrain_id INTEGER, abundance INTEGER)")
execute(dbTemp, "CREATE TABLE bstrains (t_creation REAL, bstrain_id INTEGER, parent_bstrain_id INTEGER, infecting_vstrain_id INTEGER)")
execute(dbTemp, "CREATE TABLE vstrains (t_creation REAL, vstrain_id INTEGER, parent_vstrain_id INTEGER, infected_bstrain_id INTEGER)")
execute(dbTemp, "CREATE TABLE babundance (t REAL, bstrain_id INTEGER, abundance INTEGER)")
execute(dbTemp, "CREATE TABLE vabundance (t REAL, vstrain_id INTEGER, abundance INTEGER)")


execute(dbTemp, "BEGIN TRANSACTION")
execute(dbTemp,"ATTACH DATABASE '$(dbSimPath)' as dbSim")
# execute(dbTemp,"INSERT INTO babundance (t, bstrain_id, abundance) SELECT t, bstrain_id, abundance FROM dbSim.babundance WHERE run_id = $(run_id);")
# execute(dbTemp,"INSERT INTO vabundance (t, vstrain_id, abundance) SELECT t, vstrain_id, abundance FROM dbSim.vabundance WHERE run_id = $(run_id);")
execute(dbTemp,"INSERT INTO bstrains (t_creation, bstrain_id,parent_bstrain_id,infecting_vstrain_id)
SELECT t_creation, bstrain_id,parent_bstrain_id,infecting_vstrain_id
FROM dbSim.bstrains WHERE run_id = $(run_id);")
execute(dbTemp,"INSERT INTO vstrains (t_creation, vstrain_id,parent_vstrain_id,infected_bstrain_id)
SELECT t_creation, vstrain_id,parent_vstrain_id,infected_bstrain_id
FROM dbSim.vstrains WHERE run_id = $(run_id);")
execute(dbTemp,"INSERT INTO babundance (t, bstrain_id, abundance)
SELECT t, bstrain_id, abundance FROM dbSim.babundance WHERE run_id = $(run_id);")
execute(dbTemp,"INSERT INTO vabundance (t, vstrain_id, abundance)
SELECT t, vstrain_id, abundance FROM dbSim.vabundance WHERE run_id = $(run_id);")
execute(dbTemp, "COMMIT")


execute(dbTemp, "BEGIN TRANSACTION")
execute(dbTemp, "CREATE INDEX babundance_index ON babundance (t,bstrain_id,abundance)")
execute(dbTemp, "CREATE INDEX vabundance_index ON vabundance (t,vstrain_id, abundance)")
execute(dbTemp, "CREATE INDEX bstrains_index ON bstrains (t_creation, bstrain_id,parent_bstrain_id,infecting_vstrain_id)")
execute(dbTemp, "CREATE INDEX vstrains_index ON vstrains (t_creation, vstrain_id,parent_vstrain_id,infected_bstrain_id)")
execute(dbTemp, "COMMIT")


function orderTreeStrains(strainType::String)
    if strainType == "virus"
        s = "v"
        MRCA_ids = [MRCA_id for (MRCA_id,) in execute(dbTemp,"SELECT $(s)strain_id
        FROM $(s)strains WHERE parent_$(s)strain_id = 1 ORDER BY $(s)strain_id")]
        strainTreeOrder = DataFrame(vstrain_id = Int64[0,1], tree_vstrain_id = Int64[0,1])

        lineages = Dict{Int64,Array{Int64,1}}()
        for (strain_id,) in execute(dbTemp,"SELECT DISTINCT $(s)strain_id FROM $(s)strains
            ORDER BY t_creation")
            if strain_id == 1
                continue
            end
            lineages[strain_id] = searchLineage(strain_id,strainType)
        end
        while length(MRCA_ids) > 0
            println("Finding positions for MRCA: $(MRCA_ids)")
            nextGeneration!(strainTreeOrder,MRCA_ids,lineages,strainType)
            MRCA_ids = [MRCA_id for (MRCA_id,) in execute(dbTemp,"SELECT $(s)strain_id
            FROM $(s)strains WHERE parent_$(s)strain_id in ($(join(MRCA_ids,", ")))
            ORDER BY $(s)strain_id")]
        end
        strainTreeOrder |> SQLite.load!(dbOutput,"tree_vstrain_order",ifnotexists=true)
    end
    if strainType == "microbe"
        s = "b"
        MRCA_ids = [MRCA_id for (MRCA_id,) in execute(dbTemp,"SELECT $(s)strain_id
        FROM $(s)strains WHERE parent_$(s)strain_id = 1 ORDER BY $(s)strain_id")]
        strainTreeOrder = DataFrame(bstrain_id = Int64[0,1], tree_bstrain_id = Int64[0,1])

        lineages = Dict{Int64,Array{Int64,1}}()
        for (strain_id,) in execute(dbTemp,"SELECT DISTINCT $(s)strain_id FROM $(s)strains
            ORDER BY t_creation")
            if strain_id == 1
                continue
            end
            lineages[strain_id] = searchLineage(strain_id,strainType)
        end
        while length(MRCA_ids) > 0
            println("Finding positions for MRCA: $(MRCA_ids)")
            nextGeneration!(strainTreeOrder,MRCA_ids,lineages,strainType)
            MRCA_ids = [MRCA_id for (MRCA_id,) in execute(dbTemp,"SELECT $(s)strain_id
            FROM $(s)strains WHERE parent_$(s)strain_id in ($(join(MRCA_ids,", ")))
            ORDER BY $(s)strain_id")]
        end
        strainTreeOrder |> SQLite.load!(dbOutput,"tree_bstrain_order",ifnotexists=true)
    end
end


function nextGeneration!(strainTreeOrder::DataFrame,
    MRCA_ids::Array{Int64,1},lineages::Dict,strainType::String)
    # MRCA_ids = sort([Int64(key) for key in keys(lineages)])
    if strainType == "virus"
        s = "v"
        parentStrainIDs = Int64[]
        for MRCA_id in MRCA_ids
            (strain,) = execute(dbTemp, "SELECT parent_$(s)strain_id FROM $(s)strains
            WHERE $(s)strain_id = $(MRCA_id) ORDER BY parent_$(s)strain_id")
            union!(parentStrainIDs,strain.parent_vstrain_id)
        end
        sort!(parentStrainIDs)
        println("Parents of MRCA are: $(parentStrainIDs)")
        treePositions!(strainTreeOrder,parentStrainIDs,lineages,strainType)
    end
    if strainType == "microbe"
        s = "b"
        parentStrainIDs = Int64[]
        for MRCA_id in MRCA_ids
            (strain,) = execute(dbTemp, "SELECT parent_$(s)strain_id FROM $(s)strains
            WHERE $(s)strain_id = $(MRCA_id) ORDER BY parent_$(s)strain_id")
            union!(parentStrainIDs,strain.parent_bstrain_id)
        end
        sort!(parentStrainIDs)
        println("Parents of MRCA are: $(parentStrainIDs)")
        treePositions!(strainTreeOrder,parentStrainIDs,lineages,strainType)
    end
end

function treePositions!(strainTreeOrder::DataFrame,
    parentStrainIDs::Array{Int64,1},lineages::Dict,strainType::String)
    if strainType == "virus"
        s = "v"
        for parentID in parentStrainIDs
            println("Locating descendants of strain $(parentID)")
            subMRCA_ids = [strain for (strain,) in execute(dbTemp, "SELECT $(s)strain_id FROM $(s)strains
            WHERE parent_$(s)strain_id = $(parentID) ORDER BY $(s)strain_id")]
            firstTreeStrainID = strainTreeOrder[strainTreeOrder.vstrain_id .== parentID,:].tree_vstrain_id[1]
            for subMRCA in subMRCA_ids
                if subMRCA == last(subMRCA_ids)
                    push!(strainTreeOrder,[subMRCA, firstTreeStrainID + 1])
                    continue
                end
                treeStrainID = 1 + firstTreeStrainID + sum([length(lineages[nestMRCA_id])
                for nestMRCA_id in subMRCA_ids if nestMRCA_id > subMRCA])
                push!(strainTreeOrder,[subMRCA, treeStrainID])
            end
        end
    end
    if strainType == "microbe"
        s = "b"
        for parentID in parentStrainIDs
            println("Locating descendants of strain $(parentID)")
            subMRCA_ids = [strain for (strain,) in execute(dbTemp, "SELECT $(s)strain_id FROM $(s)strains
            WHERE parent_$(s)strain_id = $(parentID) ORDER BY $(s)strain_id")]
            firstTreeStrainID = strainTreeOrder[strainTreeOrder.bstrain_id .== parentID,:].tree_bstrain_id[1]
            for subMRCA in subMRCA_ids
                if subMRCA == last(subMRCA_ids)
                    push!(strainTreeOrder,[subMRCA, firstTreeStrainID + 1])
                    continue
                end
                treeStrainID = 1 + firstTreeStrainID + sum([length(lineages[nestMRCA_id])
                for nestMRCA_id in subMRCA_ids if nestMRCA_id > subMRCA])
                push!(strainTreeOrder,[subMRCA, treeStrainID])
            end
        end
    end
end

function newAbundanceTimeSeries(strainType::String)
    if strainType == "virus"
        s = "v"
        for (time,) in execute(dbTemp,"SELECT DISTINCT t
            FROM $(s)abundance ORDER BY t")
            println("Changing strain IDs at time = $(time)")
            for (strain,) in execute(dbTemp,"SELECT $(s)strain_id FROM $(s)abundance
                WHERE t = $(time) ORDER BY $(s)strain_id")
                (new,) = execute(dbOutput, "SELECT tree_$(s)strain_id
                FROM tree_$(s)strain_order WHERE $(s)strain_id = $(strain)")
                (abundance,) = execute(dbTemp, "SELECT abundance
                FROM $(s)abundance WHERE $(s)strain_id = $(strain) AND t = $(time)")
                execute(dbOutput, "INSERT INTO tree_$(s)abundance VALUES (?,?,?)",
                (time,new.tree_vstrain_id,abundance.abundance))
            end
        end
    end
    if strainType == "microbe"
        s = "b"
        for (time,) in execute(dbTemp,"SELECT DISTINCT t
            FROM $(s)abundance ORDER BY t")
            println("Changing strain IDs at time = $(time)")
            for (strain,) in execute(dbTemp,"SELECT $(s)strain_id FROM $(s)abundance
                WHERE t = $(time) ORDER BY $(s)strain_id")
                (new,) = execute(dbOutput, "SELECT tree_$(s)strain_id
                FROM tree_$(s)strain_order WHERE $(s)strain_id = $(strain)")
                (abundance,) = execute(dbTemp, "SELECT abundance
                FROM $(s)abundance WHERE $(s)strain_id = $(strain) AND t = $(time)")
                execute(dbOutput, "INSERT INTO tree_$(s)abundance VALUES (?,?,?)",
                (time,new.tree_bstrain_id,abundance.abundance))
            end
        end
    end
end

function findExtinctionTimes(strainType::String)
    if strainType == "virus"
        s = "v"
        strainsBefore = Int64[1]
        endTime = maximum([time for (time,) in execute(dbTemp,"SELECT DISTINCT t
            FROM $(s)abundance ORDER BY t")])
        for (time,) in execute(dbTemp,"SELECT DISTINCT t
            FROM $(s)abundance ORDER BY t")
            println("Searching for strain extinctions at time = $(time)")
            strains = [strain for (strain,) in execute(dbTemp,"SELECT $(s)strain_id FROM $(s)abundance
                WHERE t = $(time)")]
            extinct = setdiff(strainsBefore,strains)
            if length(extinct) > 0
                for strain in extinct
                    (t,) = execute(dbTemp, "SELECT t_creation FROM $(s)strains
                    WHERE $(s)strain_id = $(strain)")
                    execute(dbOutput,"INSERT INTO $(s)strain_creation_extinction
                    VALUES (?,?,?)",(strain,t.t_creation,time-1))

                    (pstrain,) = execute(dbTemp, "SELECT parent_$(s)strain_id FROM $(s)strains
                    WHERE $(s)strain_id = $(strain)")
                    (ibstrain,) = execute(dbTemp, "SELECT infected_bstrain_id FROM $(s)strains
                    WHERE $(s)strain_id = $(strain)")

                    (pstrain,) = execute(dbOutput, "SELECT tree_$(s)strain_id
                    FROM tree_$(s)strain_order WHERE $(s)strain_id = $(pstrain.parent_vstrain_id)")
                    (ibstrain,) = execute(dbOutput, "SELECT tree_bstrain_id
                    FROM tree_bstrain_order WHERE bstrain_id = $(ibstrain.infected_bstrain_id)")

                    (tree,) = execute(dbOutput, "SELECT tree_$(s)strain_id
                    FROM tree_$(s)strain_order WHERE $(s)strain_id = $(strain)")
                    execute(dbOutput,"INSERT INTO tree_$(s)strain_creation_extinction
                    VALUES (?,?,?,?,?)",(tree.tree_vstrain_id,t.t_creation,time-1,
                    pstrain.tree_vstrain_id,ibstrain.tree_bstrain_id))
                end
            end
            if time == endTime
                for strain in strains
                    (t,) = execute(dbTemp, "SELECT t_creation FROM $(s)strains
                    WHERE $(s)strain_id = $(strain)")
                    execute(dbOutput,"INSERT INTO $(s)strain_creation_extinction
                    VALUES (?,?,?)",(strain,t.t_creation,time))

                    (pstrain,) = execute(dbTemp, "SELECT parent_$(s)strain_id FROM $(s)strains
                    WHERE $(s)strain_id = $(strain)")
                    (ibstrain,) = execute(dbTemp, "SELECT infected_bstrain_id FROM $(s)strains
                    WHERE $(s)strain_id = $(strain)")

                    (pstrain,) = execute(dbOutput, "SELECT tree_$(s)strain_id
                    FROM tree_$(s)strain_order WHERE $(s)strain_id = $(pstrain.parent_vstrain_id)")
                    (ibstrain,) = execute(dbOutput, "SELECT tree_bstrain_id
                    FROM tree_bstrain_order WHERE bstrain_id = $(ibstrain.infected_bstrain_id)")

                    (tree,) = execute(dbOutput, "SELECT tree_$(s)strain_id
                    FROM tree_$(s)strain_order WHERE $(s)strain_id = $(strain)")
                    execute(dbOutput,"INSERT INTO tree_$(s)strain_creation_extinction
                    VALUES (?,?,?,?,?)",(tree.tree_vstrain_id,t.t_creation,time,
                    pstrain.tree_vstrain_id,ibstrain.tree_bstrain_id))
                end
            end
            strainsBefore = strains[:]
        end
    end
    if strainType == "microbe"
        s = "b"
        strainsBefore = Int64[1]
        endTime = maximum([time for (time,) in execute(dbTemp,"SELECT DISTINCT t
            FROM $(s)abundance ORDER BY t")])
        naiveAlive = true # CHANGE THIS
        # WHEN YOU INCLUDE A BETTER VERSION OF IMMIGRATION
        for (time,) in execute(dbTemp,"SELECT DISTINCT t
            FROM $(s)abundance ORDER BY t")
            println("Searching for strain extinctions at time = $(time)")
            strains = [strain for (strain,) in execute(dbTemp,"SELECT $(s)strain_id FROM $(s)abundance
                WHERE t = $(time)")]
            extinct = setdiff(strainsBefore,strains)
            (naive,) = execute(dbTemp,"SELECT abundance FROM $(s)abundance
            WHERE $(s)strain_id = 1 AND t = $(time)") # CHANGE THIS
            # WHEN YOU INCLUDE A BETTER VERSION OF IMMIGRATION
            if naive.abundance == 0 && naiveAlive # CHANGE THIS
            # WHEN YOU INCLUDE A BETTER VERSION OF IMMIGRATION
                (t,) = execute(dbTemp, "SELECT t_creation FROM $(s)strains
                WHERE $(s)strain_id = 1")
                execute(dbOutput,"INSERT INTO $(s)strain_creation_extinction
                VALUES (?,?,?)",(1,t.t_creation,time-1))

                (pstrain,) = execute(dbTemp, "SELECT parent_$(s)strain_id FROM $(s)strains
                WHERE $(s)strain_id = 1")
                (ivstrain,) = execute(dbTemp, "SELECT infecting_vstrain_id FROM $(s)strains
                WHERE $(s)strain_id = 1")

                (pstrain,) = execute(dbOutput, "SELECT tree_$(s)strain_id
                FROM tree_$(s)strain_order WHERE $(s)strain_id = $(pstrain.parent_bstrain_id)")
                (ivstrain,) = execute(dbOutput, "SELECT tree_vstrain_id
                FROM tree_vstrain_order WHERE vstrain_id = $(ivstrain.infecting_vstrain_id)")

                execute(dbOutput,"INSERT INTO tree_$(s)strain_creation_extinction
                VALUES (?,?,?,?,?)",(1,t.t_creation,time-1,
                pstrain.tree_bstrain_id,ivstrain.tree_vstrain_id))
                naiveAlive = false
            end
            if length(extinct) > 0
                for strain in extinct
                    (t,) = execute(dbTemp, "SELECT t_creation FROM $(s)strains
                    WHERE $(s)strain_id = $(strain)")
                    execute(dbOutput,"INSERT INTO $(s)strain_creation_extinction
                    VALUES (?,?,?)",(strain,t.t_creation,time-1))

                    (pstrain,) = execute(dbTemp, "SELECT parent_$(s)strain_id FROM $(s)strains
                    WHERE $(s)strain_id = $(strain)")
                    (ivstrain,) = execute(dbTemp, "SELECT infecting_vstrain_id FROM $(s)strains
                    WHERE $(s)strain_id = $(strain)")

                    (pstrain,) = execute(dbOutput, "SELECT tree_$(s)strain_id
                    FROM tree_$(s)strain_order WHERE $(s)strain_id = $(pstrain.parent_bstrain_id)")
                    (ivstrain,) = execute(dbOutput, "SELECT tree_vstrain_id
                    FROM tree_vstrain_order WHERE vstrain_id = $(ivstrain.infecting_vstrain_id)")

                    (tree,) = execute(dbOutput, "SELECT tree_$(s)strain_id
                    FROM tree_$(s)strain_order WHERE $(s)strain_id = $(strain)")
                    execute(dbOutput,"INSERT INTO tree_$(s)strain_creation_extinction
                    VALUES (?,?,?,?,?)",(tree.tree_bstrain_id,t.t_creation,time-1,
                    pstrain.tree_bstrain_id,ivstrain.tree_vstrain_id))
                end
            end
            if time == endTime
                strains = [strain for (strain,abundance)
                in execute(dbTemp,"SELECT $(s)strain_id, abundance FROM $(s)abundance
                    WHERE t = $(time)") if abundance != 0] # CHANGE THIS
                    # WHEN YOU INCLUDE A BETTER VERSION OF IMMIGRATION
                for strain in strains
                    (t,) = execute(dbTemp, "SELECT t_creation FROM $(s)strains
                    WHERE $(s)strain_id = $(strain)")
                    execute(dbOutput,"INSERT INTO $(s)strain_creation_extinction
                    VALUES (?,?,?)",(strain,t.t_creation,time))

                    (pstrain,) = execute(dbTemp, "SELECT parent_$(s)strain_id FROM $(s)strains
                    WHERE $(s)strain_id = $(strain)")
                    (ivstrain,) = execute(dbTemp, "SELECT infecting_vstrain_id FROM $(s)strains
                    WHERE $(s)strain_id = $(strain)")

                    (pstrain,) = execute(dbOutput, "SELECT tree_$(s)strain_id
                    FROM tree_$(s)strain_order WHERE $(s)strain_id = $(pstrain.parent_bstrain_id)")
                    (ivstrain,) = execute(dbOutput, "SELECT tree_vstrain_id
                    FROM tree_vstrain_order WHERE vstrain_id = $(ivstrain.infecting_vstrain_id)")

                    (tree,) = execute(dbOutput, "SELECT tree_$(s)strain_id
                    FROM tree_$(s)strain_order WHERE $(s)strain_id = $(strain)")
                    execute(dbOutput,"INSERT INTO tree_$(s)strain_creation_extinction
                    VALUES (?,?,?,?,?)",(tree.tree_bstrain_id,t.t_creation,time,
                    pstrain.tree_bstrain_id,ivstrain.tree_vstrain_id))
                end
            end
            strainsBefore = strains[:]
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

orderTreeStrains("microbe")
newAbundanceTimeSeries("microbe")
orderTreeStrains("virus")
newAbundanceTimeSeries("virus")
findExtinctionTimes("microbe")
findExtinctionTimes("virus")


function createindices()
    println("(Creating run_id indices...)")
    db = SQLite.DB(dbOutputPath)
    execute(db, "BEGIN TRANSACTION")
    for (table_name,) in execute(
        db, "SELECT name FROM sqlite_schema
        WHERE type='table' ORDER BY name;")
        # cols = [info.name for info in execute(db,"PRAGMA table_info($(table_name))")]
        if in(table_name,["tree_bstrain_order"])
            execute(db, "CREATE INDEX $(table_name)_index ON $(table_name) (bstrain_id)")
        end
        if in(table_name,["tree_vstrain_order"])
            execute(db, "CREATE INDEX $(table_name)_index ON $(table_name) (vstrain_id)")
        end
        if in(table_name,["tree_babundance","tree_vabundance"])
            execute(db, "CREATE INDEX $(table_name)_index ON $(table_name) (t)")
        end
    end
    execute(db, "COMMIT")
end
createindices()


println("Complete!")
