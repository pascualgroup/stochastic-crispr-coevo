#!/usr/bin/env julia

println("(Annoying Julia compilation delay...)")

using SQLite
using DataFrames
using SQLite.DBInterface: execute
# using BenchmarkTools


run_id = ARGS[1]

## Define Paths ##
SCRIPT_PATH = abspath(dirname(PROGRAM_FILE))

dbSimPath = joinpath(SCRIPT_PATH,"..","..","..","simulation","sweep_db_gathered.sqlite") # cluster
dbOutputPath = joinpath("matches_output.sqlite") # cluster

# # dbSimPath = joinpath("/Volumes/Yadgah/sweep_db_gathered.sqlite") # local
# dbSimPath = joinpath("/Volumes/Yadgah","run_id1455_combo73_replicate15.sqlite") # local
# dbOutputPath = joinpath("/Volumes/Yadgah/matches_output.sqlite") # local

if isfile(dbOutputPath)
    error("matches_output.sqlite already exists; delete first")
end # cluster
##

dbSim = SQLite.DB(dbSimPath)
# dbTemp = SQLite.DB(dbSimPath)
dbOutput = SQLite.DB(dbOutputPath)

execute(dbOutput, "CREATE TABLE bstrain_to_vstrain_matches (t REAL, bstrain_id INTEGER, vstrain_id INTEGER,
time_specific_match_id INTEGER, match_length INTEGER)")

execute(dbOutput, "CREATE TABLE matches_spacers (t REAL, time_specific_match_id INTEGER, match_type INTEGER, spacer_id INTEGER)")

execute(dbOutput, "CREATE TABLE bstrain_matched_spacers
(t REAL, bstrain_id INTEGER, matched_spacer_id INTEGER)")

execute(dbOutput, "CREATE TABLE vstrain_matched_pspacers
(t REAL, vstrain_id INTEGER, matched_pspacer_id INTEGER)")

execute(dbOutput, "CREATE TABLE bstrain_num_matched_spacers (t REAL, bstrain_id INTEGER,
num_matched_spacers INTEGER)")

execute(dbOutput, "CREATE TABLE vstrain_num_matched_pspacers (t REAL, vstrain_id INTEGER,
num_matched_pspacers INTEGER)")

execute(dbOutput, "CREATE TABLE bstrain_to_vstrain_0matches (t REAL, bstrain_id INTEGER, vstrain_id INTEGER,
time_specific_0match_id INTEGER)")

execute(dbOutput, "CREATE TABLE bstrain_num_total_spacers (bstrain_id INTEGER, num_total_spacers INTEGER)")

# Create temporary database that is a copy of the main database at the run_id value of the script's argument
# @btime begin
# @time begin
dbTemp = SQLite.DB()
# dbTemp = SQLite.DB("/Volumes/Yadgah/test.sqlite") # local
execute(dbTemp, "CREATE TABLE babundance (t REAL, bstrain_id INTEGER)")
execute(dbTemp, "CREATE TABLE bstrain1Abundance (t REAL, abundance INTEGER)")
execute(dbTemp, "CREATE TABLE vabundance (t REAL, vstrain_id INTEGER)")
execute(dbTemp, "CREATE TABLE bspacers (bstrain_id INTEGER, spacer_id INTEGER)")
execute(dbTemp, "CREATE TABLE vpspacers (vstrain_id INTEGER, spacer_id INTEGER)")


execute(dbTemp, "BEGIN TRANSACTION")
execute(dbTemp,"ATTACH DATABASE '$(dbSimPath)' as dbSim")
execute(dbTemp,"INSERT INTO babundance (t, bstrain_id) SELECT t, bstrain_id FROM dbSim.babundance WHERE run_id = $(run_id);")
execute(dbTemp,"INSERT INTO bstrain1Abundance (t, abundance) SELECT t, abundance FROM dbSim.babundance
WHERE run_id = $(run_id) AND bstrain_id = 1;")
execute(dbTemp,"INSERT INTO vabundance (t, vstrain_id) SELECT t, vstrain_id FROM dbSim.vabundance WHERE run_id = $(run_id);")
execute(dbTemp,"INSERT INTO bspacers (bstrain_id,spacer_id) SELECT bstrain_id,spacer_id FROM dbSim.bspacers WHERE run_id = $(run_id);")
execute(dbTemp,"INSERT INTO vpspacers (vstrain_id,spacer_id) SELECT vstrain_id, spacer_id FROM dbSim.vpspacers WHERE run_id = $(run_id);")
execute(dbTemp, "COMMIT")


execute(dbTemp, "BEGIN TRANSACTION")
execute(dbTemp, "CREATE INDEX babundance_index ON babundance (t,bstrain_id)")
execute(dbTemp, "CREATE INDEX bstrain1_index ON bstrain1Abundance (t)")
execute(dbTemp, "CREATE INDEX vabundance_index ON vabundance (t,vstrain_id)")
execute(dbTemp, "CREATE INDEX bspacers_index ON bspacers (bstrain_id,spacer_id)")
execute(dbTemp, "CREATE INDEX vspacers_index ON vpspacers (vstrain_id,spacer_id)")
execute(dbTemp, "COMMIT")


println("Identifying matches of run $(run_id)")
function identifyMatches()
    for (time,) in execute(dbTemp, "SELECT DISTINCT t FROM babundance ORDER BY t")
        println("Searching (non-)matches at time $(time)")
        match_id = 1
        nomatch_id = 1
        # if time > 60.0
        #     return
        # end
        vstrainMatches = Dict((vstrain_id, Int64[]) for (vstrain_id,) in
        execute(dbTemp, "SELECT vstrain_id FROM vabundance
            WHERE t = $(time) ORDER BY vstrain_id"))
        for (bstrain_id,) in execute(dbTemp, "SELECT bstrain_id FROM babundance
            WHERE t = $(time) ORDER BY bstrain_id")
            if bstrain_id == 1
                (abund1,) = execute(dbTemp, "SELECT abundance FROM bstrain1Abundance
                    WHERE t = $(time)")
                if abund1.abundance > 0
                    execute(dbOutput, "BEGIN TRANSACTION")
                    execute(dbOutput, "INSERT INTO bstrain_num_matched_spacers VALUES (?,?,?)",
                    (time, bstrain_id, 0))
                    for (vstrain_id,) in execute(dbTemp, "SELECT vstrain_id FROM vabundance
                        WHERE t = $(time) ORDER BY vstrain_id")
                        execute(dbOutput, "INSERT INTO bstrain_to_vstrain_0matches VALUES (?,?,?,?)",
                        (time, bstrain_id, vstrain_id, nomatch_id))
                        nomatch_id += 1
                    end
                    execute(dbOutput, "COMMIT")
                end
                continue
            end
            spacers = [spacer_id for (spacer_id,) in
            execute(dbTemp, "SELECT spacer_id FROM bspacers WHERE bstrain_id = $(bstrain_id)")]
            matchedSpacers, match_id, vstrainMatches, nomatch_id =
            bstrainTOvstrain(time,bstrain_id,spacers,match_id,vstrainMatches,nomatch_id)
            numMatched = length(matchedSpacers)
            execute(dbOutput, "BEGIN TRANSACTION")
            execute(dbOutput, "INSERT INTO bstrain_num_matched_spacers
            VALUES (?,?,?)",(time, bstrain_id, numMatched))
            if numMatched > 0
                for spacer_id in matchedSpacers
                    execute(dbOutput, "INSERT INTO bstrain_matched_spacers
                    VALUES (?,?,?)", (time, bstrain_id, spacer_id))
                end
            end
            execute(dbOutput, "COMMIT")
        end
        identify_vstrainMatches(time,vstrainMatches)
    end
end

function bstrainTOvstrain(time,bstrain_id,spacers,match_id,vstrainMatches::Dict,nomatch_id)
    matchedSpacers = Int64[]
    for (vstrain_id,) in execute(dbTemp, "SELECT vstrain_id FROM vabundance
        WHERE t = $(time) ORDER BY vstrain_id")

        pspacers = [pspacer_id for (pspacer_id,) in
        execute(dbTemp, "SELECT spacer_id FROM vpspacers WHERE vstrain_id = $(vstrain_id)")]
        matches = intersect(spacers,pspacers)
        union!(matchedSpacers,matches)
        matchLength = length(matches) # number of matches

        if matchLength == 0
            execute(dbOutput, "INSERT INTO bstrain_to_vstrain_0matches VALUES (?,?,?,?)",
            (time, bstrain_id, vstrain_id, nomatch_id))
            nomatch_id += 1
            continue
        end
        union!(vstrainMatches[vstrain_id],matches)

        execute(dbOutput, "BEGIN TRANSACTION")
        for spacer_id in matches
            execute(dbOutput, "INSERT INTO matches_spacers VALUES (?,?,?,?)",(time, match_id, length(matches), spacer_id))
        end
        execute(dbOutput, "INSERT INTO bstrain_to_vstrain_matches VALUES (?,?,?,?,?)",
        (time, bstrain_id, vstrain_id, match_id, matchLength))
        execute(dbOutput, "COMMIT")

        match_id += 1
    end
    return matchedSpacers, match_id, vstrainMatches, nomatch_id
end

function identify_vstrainMatches(time,vstrainMatches::Dict)
    execute(dbOutput, "BEGIN TRANSACTION")
    for (vstrain_id,) in execute(dbTemp, "SELECT vstrain_id FROM vabundance
        WHERE t = $(time) ORDER BY vstrain_id")
        if !issubset(vstrain_id,[keys(vstrainMatches)...])
            execute(dbOutput, "INSERT INTO vstrain_num_matched_pspacers VALUES (?,?,?)",
            (time, vstrain_id, 0))
            continue
        else
            execute(dbOutput, "INSERT INTO vstrain_num_matched_pspacers VALUES (?,?,?)",
            (time, vstrain_id, length(vstrainMatches[vstrain_id])))
        end
        for pspacer_id in vstrainMatches[vstrain_id]
            execute(dbOutput, "INSERT INTO vstrain_matched_pspacers
            VALUES (?,?,?)", (time, vstrain_id, pspacer_id))
        end
    end
    execute(dbOutput, "COMMIT")
end

function numTotalSpacers()
    println("Counting total number of spacers of each bstrain")
    for (bstrain_id,) in execute(dbTemp, "SELECT DISTINCT bstrain_id FROM babundance
        ORDER BY bstrain_id")
        spacers = [spacer_id for (spacer_id,) in
        execute(dbTemp, "SELECT spacer_id FROM bspacers WHERE bstrain_id = $(bstrain_id)")]
        execute(dbOutput, "INSERT INTO bstrain_num_total_spacers VALUES (?,?)",(bstrain_id, length(spacers)))
    end
end

identifyMatches()
numTotalSpacers()

# end # @btime
# end # @time

println("Complete!")
