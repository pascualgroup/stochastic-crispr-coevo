#!/usr/bin/env julia

println("(Annoying Julia compilation delay...)")

using SQLite
using DataFrames
using SQLite.DBInterface: execute


run_id = ARGS[1]

## Define Paths ##
SCRIPT_PATH = abspath(dirname(PROGRAM_FILE))

#dbSimPath = joinpath(SCRIPT_PATH,"..","..","simulation","sweep_db_gathered.sqlite") # cluster
#dbOutputPath = joinpath("memory_output.sqlite") # cluster

dbSimPath = joinpath("/Volumes/Yadgah/sweep_db_gathered.sqlite") # local
dbSimPath = joinpath("/Volumes/Yadgah","run_id1455_combo73_replicate15.sqlite") # local
dbOutputPath = joinpath("/Volumes/Yadgah/memory_output.sqlite") # local

if isfile(dbOutputPath)
    error("memory_output.sqlite already exists; delete first")
end # cluster
##

dbSim = SQLite.DB(dbSimPath)
dbOutput = SQLite.DB(dbOutputPath)

execute(dbOutput, "CREATE TABLE bstrain_memory_size (t REAL, bstrain_id INTEGER,
num_spacers INTEGER)")
execute(dbOutput, "CREATE TABLE bstrain_memory_abundances (t REAL,
memory_size INTEGER, abundance INTEGER)")

execute(dbOutput, "CREATE TABLE bstrain_overlap (t REAL, bstrain_id1 INTEGER, bstrain_id2 INTEGER, overlap_id)")

execute(dbOutput, "CREATE TABLE overlap_ids (t REAL, bstrain_id1 INTEGER, bstrain_id2 INTEGER, overlap_id)")

# Create temporary database that is a copy of the main database at the run_id value of the script's argument
dbTemp = SQLite.DB()
#dbTemp = SQLite.DB("/Volumes/Yadgah/test.sqlite") # local
execute(dbTemp, "CREATE TABLE summary (t REAL, microbial_abundance INTEGER, viral_abundance INTEGER)")
execute(dbTemp, "CREATE TABLE babundance (t REAL, bstrain_id INTEGER, abundance INTEGER)")
execute(dbTemp, "CREATE TABLE vabundance (t REAL, vstrain_id INTEGER, abundance INTEGER)")
execute(dbTemp, "CREATE TABLE bspacers (bstrain_id INTEGER, spacer_id INTEGER)")
execute(dbTemp, "CREATE TABLE vpspacers (vstrain_id INTEGER, spacer_id INTEGER)")


execute(dbTemp, "BEGIN TRANSACTION")
execute(dbTemp,"ATTACH DATABASE '$(dbSimPath)' as dbSim")
execute(dbTemp,"INSERT INTO summary(t, microbial_abundance,viral_abundance) SELECT t, microbial_abundance,viral_abundance FROM dbSim.summary WHERE run_id = $(run_id);")
execute(dbTemp,"INSERT INTO babundance (t, bstrain_id, abundance) SELECT t, bstrain_id, abundance FROM dbSim.babundance WHERE run_id = $(run_id);")
execute(dbTemp,"INSERT INTO vabundance (t, vstrain_id, abundance) SELECT t, vstrain_id, abundance FROM dbSim.vabundance WHERE run_id = $(run_id);")
execute(dbTemp,"INSERT INTO bspacers (bstrain_id,spacer_id) SELECT bstrain_id,spacer_id FROM dbSim.bspacers WHERE run_id = $(run_id);")
execute(dbTemp,"INSERT INTO vpspacers (vstrain_id,spacer_id) SELECT vstrain_id, spacer_id FROM dbSim.vpspacers WHERE run_id = $(run_id);")
execute(dbTemp, "COMMIT")


execute(dbTemp, "BEGIN TRANSACTION")
execute(dbTemp, "CREATE INDEX summary_index ON summary (t,microbial_abundance,viral_abundance)")
execute(dbTemp, "CREATE INDEX babundance_index ON babundance (t,bstrain_id,abundance)")
execute(dbTemp, "CREATE INDEX vabundance_index ON vabundance (t,vstrain_id)")
execute(dbTemp, "CREATE INDEX bspacers_index ON bspacers (bstrain_id,spacer_id)")
execute(dbTemp, "CREATE INDEX vspacers_index ON vpspacers (vstrain_id,spacer_id)")
execute(dbTemp, "COMMIT")


println("Processing match diversity of run $(run_id)")
function identifyMatches()
    for (time,) in execute(dbTemp, "SELECT t FROM summary ORDER BY t")
        println("Computing matches at time $(time)")
        if time > 30.0
            return
        end
        match_id = 1
        (btotal,) = execute(dbTemp, "SELECT microbial_abundance FROM summary WHERE t = $(time)")
        btotal = btotal.microbial_abundance
        (vtotal,) = execute(dbTemp, "SELECT viral_abundance FROM summary WHERE t = $(time)")
        vtotal = vtotal.viral_abundance

        for (bstrain_id,) in execute(dbTemp, "SELECT bstrain_id FROM babundance
            WHERE t = $(time) ORDER BY bstrain_id")
            matchedSpacers = Int64[]
            spacerMatches = [spacer_id for (spacer_id,) in
            execute(dbTemp, "SELECT spacer_id FROM bspacers WHERE bstrain_id = $(bstrain_id)")]

            for (vstrain_id,) in execute(dbTemp, "SELECT vstrain_id FROM vabundance
                WHERE t = $(time) ORDER BY vstrain_id")
                pspacers = [pspacer_id for (pspacer_id,) in
                execute(dbTemp, "SELECT spacer_id FROM vpspacers WHERE vstrain_id = $(vstrain_id)")]

                if bstrain_id == 1 # naive hosts
                    match_type = 0
                else
                    intersect!(spacerMatches,pspacers)
                    append!(matchedSpacers,spacerMatches)
                    unique!(matchedSpacers)
                    match_type = length(spacerMatches) # number of matches
                end
                if match_type == 0
                    continue
                end
                execute(dbOutput, "BEGIN TRANSACTION")
                for spacer_id in spacerMatches
                    execute(dbOutput, "INSERT INTO match_ids VALUES (?,?,?)",(time, match_id, spacer_id))
                end
                execute(dbOutput, "INSERT INTO bstrain_to_vstrain_matches VALUES (?,?,?,?,?)",
                (time, bstrain_id, vstrain_id, match_id, match_type))
                execute(dbOutput, "COMMIT")

                match_id += 1
            end
        end
    end
end

function bstrainTypeAbundances()
    for (time,) in execute(dbTemp, "SELECT t FROM summary ORDER BY t")
        println("Computing microbial match abundances at time $(time)")

        for (bstrain_id,) in execute(dbOutput, "SELECT DISTINCT bstrain_id FROM bstrain_to_vstrain_matches WHERE t = $(time)")
            println("bstrain_id: $(bstrain_id)")
            match_ids = [match_id for (match_id,) in
            execute(dbOutput, "SELECT match_id FROM bstrain_to_vstrain_matches WHERE bstrain_id = $(bstrain_id) AND t= $(time)")]
            spacers = unique([spacer_id for (spacer_id,) in
            execute(dbOutput, "SELECT spacer_id FROM matches
            WHERE match_id in ($(join(match_ids,", ")))")])

            for spacer_id in spacers
                execute(dbOutput, "INSERT INTO bstrain_matched_memories VALUES (?,?,?)",
                (time,bstrain_id,spacer_id))
            end

            (abundance,) = execute(dbTemp, "SELECT abundance FROM babundance WHERE bstrain_id = $(bstrain_id) AND t= $(time)")
            execute(dbOutput, "INSERT INTO bstrain_matches VALUES (?,?,?,?)",
            (time,bstrain_id,length(spacers),abundance.abundance))
        end
        for (match_type,) in execute(dbOutput,"SELECT num_matches FROM bstrain_matches ORDER BY num_matches")
            abundances = [abundance for (abundance,) in execute(dbOutput, "SELECT abundance
            FROM bstrain_matches WHERE num_matches = $(match_type) WHERE t = $(time)")]
            execute(dbOutput, "INSERT INTO bstrain_type_abundances VALUES (?,?,?,?)",
            (time,match_type,length(abundances),sum(abundances)))
        end
    end
end

function vstrainTypeAbundances()
    println("Computing viral match abundances at time $(time)")
    for (vstrain_id,) in execute(dbOutput, "SELECT DISTINCT vstrain_id FROM bstrain_to_vstrain_matches")
        match_ids = [match_id for (match_id,) in
        execute(dbOutput, "SELECT match_id FROM bstrain_to_vstrain_matches WHERE vstrain_id = $(vstrain_id)")]
        spacers = unique([spacer_id for (spacer_id,) in
        execute(dbOutput, "SELECT spacer_id FROM bstrain_to_vstrain_matches
        WHERE match_id in ($(join(match_ids,", ")))")])

        for spacer_id in spacers
            execute(dbOutput, "INSERT INTO vstrain_matched_repertoires VALUES (?,?,?)",
            (time,vstrain_id,spacer_id))
        end

        (abundance,) = execute(dbTemp, "SELECT abundance FROM vabundance WHERE vstrain_id = $(vstrain_id)")
        execute(dbOutput, "INSERT INTO vstrain_matches VALUES (?,?,?,?)",
        (time,vstrain_id,length(spacers),abundance.abundance))
    end
    for (match_type,) in execute(dbOutput,"SELECT num_matches FROM vstrain_matches ORDER BY num_matches")
        abundances = [abundance for (abundance,) in execute(dbOutput, "SELECT abundance
        FROM vstrain_matches WHERE num_matches = $(match_type)")]
        execute(dbOutput, "INSERT INTO vstrain_type_abundances VALUES (?,?,?,?)",
        (time,match_type,length(abundances),sum(abundances)))
    end
end

identifyMatches()
bstrainTypeAbundances()
vstrainTypeAbundances()


# function matchTypeAbundances(time,btotal,vtotal)
#     totalMatches = length([match for (match,) in execute(dbOutput, "SELECT DISTINCT match_id FROM matches WHERE t = $(time)")])
#     matchDivTotal = 0
#     for (match_type,) in execute(dbOutput,"SELECT DISTINCT match_type FROM bstrain_to_vstrain_matches WHERE t = $(time) ORDER BY match_type")
#         bstrains = [bstrain for (bstrain,) in execute(dbOutput, "SELECT DISTINCT bstrain_id FROM bstrain_to_vstrain_matches WHERE t = $(time)
#         AND match_type = $(match_type)")]
#         bstrainAbundances = [abundance for (abundance,) in execute(dbTemp, "SELECT abundance FROM babundance WHERE t = $(time)
#         AND bstrain_id in ($(join(bstrains,", ")))")]
#         bstrainTotal = sum(bstrainAbundances)
#         vstrains = [vstrain for (vstrain,) in execute(dbOutput, "SELECT DISTINCT vstrain_id FROM bstrain_to_vstrain_matches WHERE t = $(time)
#         AND match_type = $(match_type)")]
#         vstrainAbundances = [abundance for (abundance,) in execute(dbTemp, "SELECT abundance FROM vabundance WHERE t = $(time)
#         AND vstrain_id in ($(join(vstrains,", ")))")]
#         vstrainTotal = sum(vstrainAbundances)
#
#
#         (first,) = execute(dbOutput, "SELECT DISTINCT match_id FROM matches WHERE t = $(time)
#             AND match_type = $(match_type) ORDER BY match_id")
#         match_ref = [spacer for (spacer,) in
#             execute(dbOutput, "SELECT spacer_id FROM matches WHERE t = $(time) AND match_id = $(first.match_id)
#             ORDER BY spacer_id")]
#             matchDiv = 1
#             matchDivTotal += 1
#             #if time == 19.0
#                 #println("reference matches: $(match_ref)")
#             #end
#         for (match_id,) in execute(dbOutput, "SELECT DISTINCT match_id FROM matches WHERE t = $(time)
#             AND match_type = $(match_type) ORDER BY match_id")
#             matches = [spacer for (spacer,) in
#                 execute(dbOutput, "SELECT spacer_id FROM matches WHERE t = $(time) AND match_id = $(match_id)
#                 ORDER BY spacer_id")]
#                 #if time == 19.0
#                     #println("matches: $(matches)")
#                     #println("intersection: $(intersect(matches,match_ref))")
#                     #println("length: $(length(intersect(matches,match_ref)))")
#                     #println("match_type: $(match_type)")
#                     #println("truth: $(length(intersect(matches,match_ref)) != match_type)")
#                 #end
#
#             if length(intersect(matches,match_ref)) !== match_type
#                 matchDiv += 1
#                 matchDivTotal += 1
#                 match_ref = [spacer for (spacer,) in
#                     execute(dbOutput, "SELECT spacer_id FROM matches WHERE t = $(time) AND match_id = $(match_id)
#                     ORDER BY spacer_id")]
#             end
#         end
#
#         allSpacers = [spacer for spacer in execute(dbOutput, "SELECT spacer_id FROM matches
#         WHERE t = $(time) AND match_type = $(match_type)")]
#
#         match_ids = [match_id for match_id in execute(dbOutput, "SELECT DISTINCT match_id FROM matches WHERE t = $(time)
#             AND match_type = $(match_type)")]
#         execute(dbOutput, "INSERT INTO matchtype_abundances VALUES (?,?,?,?,?,?,?,?,?,?,?,?)",
#         (time, match_type, length(match_ids), length(match_ids)/totalMatches, matchDiv,
#         length(unique(allSpacers)), bstrainTotal/btotal, vstrainTotal/vtotal,
#         bstrainTotal, vstrainTotal, length(bstrains), length(vstrains)))
#     end
#
#     execute(dbOutput, "INSERT INTO total_abundances VALUES (?,?,?,?,?)",(time, btotal, vtotal, totalMatches, matchDivTotal))
# end




println("Complete!")
