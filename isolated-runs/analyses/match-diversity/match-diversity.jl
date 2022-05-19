#!/usr/bin/env julia

println("(Julia compilation delay...)")

using SQLite
using DataFrames
using SQLite.DBInterface: execute

run_id = ARGS[1]

## Define Paths ##
SCRIPT_PATH = abspath(dirname(PROGRAM_FILE))

dbSimPath = joinpath(SCRIPT_PATH,"..","..","..","simulation","sweep_db_gathered.sqlite") # cluster
dbMatchPath = #######
dbOutputPath = joinpath("match-diversity_output.sqlite") # cluster

# # dbSimPath = joinpath("/Volumes/Yadgah/sweep_db_gathered.sqlite") # local
# dbSimPath = joinpath("/Volumes/Yadgah","run_id1455_combo73_replicate15.sqlite") # local
# dbSimPath = joinpath("/Volumes/Yadgah","crispr-sweep-7-2-2022/isolates/runID3297-c66-r47/runID3297-c66-r47.sqlite") # local
# dbMatchPath = joinpath("/Volumes/Yadgah","crispr-sweep-7-2-2022/isolates/runID3297-c66-r47/matches_output.sqlite")
# dbOutputPath = joinpath("/Volumes/Yadgah/crispr-sweep-7-2-2022/isolates/runID3297-c66-r47/match-diversity_output.sqlite") # local

if isfile(dbOutputPath)
    error("match-diversity_output.sqlite already exists; delete first")
end # cluster
##

dbSim = SQLite.DB(dbSimPath)
# dbTempSim = SQLite.DB(dbSimPath)
# dbTempMatch = SQLite.DB(dbMatchPath)
dbOutput = SQLite.DB(dbOutputPath)

execute(dbOutput, "CREATE TABLE single_locus_escape_match_diversity (t REAL,
virus_shannon_diversity REAL, microbe_shannon_diversity REAL)")
# execute(dbOutput, "CREATE TABLE single_locus_escape_matches (t REAL, vstrain_id INTEGER, bstrain_id INTEGER, spacer_id INTEGER)")

execute(dbOutput, "CREATE TABLE single_locus_escape_match_abundances (t REAL,
spacer_id INTEGER, vabundance INTEGER, babundance INTEGER)")

execute(dbOutput, "CREATE TABLE vMatchAbundances (t REAL, match0 INTEGER,
match1 INTEGER, match2 INTEGER, match3 INTEGER)")
execute(dbOutput, "CREATE TABLE bMatchAbundances
(t REAL, match0 INTEGER, match1 INTEGER,
match2 INTEGER, match3 INTEGER)")


dbTempSim = SQLite.DB()
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

## FINISH THIS PART
dbTempMatch = SQLite.DB()
##


function singleLocusMatches()
    for (time,) in execute(dbTempSim, "SELECT DISTINCT t FROM babundance ORDER BY t")
        println("Identifying single locus escape matches at time $(time)")
        escapeMatches = let
            vstrains0 = [vstrain_id for (vstrain_id,) in
                execute(dbTempMatch, "SELECT DISTINCT vstrain_id FROM bstrain_to_vstrain_0matches
                WHERE t = $(time) ORDER BY vstrain_id")]
                if length(vstrains0) == 0
                    0
                else
                    vstrains1 = [vstrain_id for (vstrain_id,) in
                        execute(dbTempMatch, "SELECT vstrain_id FROM bstrain_to_vstrain_matches
                        WHERE t = $(time) AND match_length = 1 AND vstrain_id in ($(join(vstrains0,", ")))
                        ORDER BY time_specific_match_id")]
                        if length(vstrains1) == 0
                            0
                        else
                            bstrains1 = [bstrain_id for (bstrain_id,) in
                                execute(dbTempMatch, "SELECT bstrain_id FROM bstrain_to_vstrain_matches
                                WHERE t = $(time) AND match_length = 1 AND vstrain_id in ($(join(vstrains0,", ")))
                                ORDER BY time_specific_match_id")]
                            match_ids = [match_id for (match_id,) in
                                execute(dbTempMatch, "SELECT time_specific_match_id FROM bstrain_to_vstrain_matches
                                WHERE t = $(time) AND match_length = 1 AND vstrain_id in ($(join(vstrains0,", ")))
                                ORDER BY time_specific_match_id")]
                            spacer_ids = [spacer_ids for (spacer_ids,) in
                                execute(dbTempMatch, "SELECT spacer_id FROM matches_spacers
                                WHERE t = $(time) AND time_specific_match_id in ($(join(match_ids,", ")))
                                ORDER BY time_specific_match_id")]
                            DataFrame(t=time*ones(length(vstrains1)),
                            vstrain_id = vstrains1, bstrain_id = bstrains1, spacer_id = spacer_ids)
                        end
                end
            end
        if typeof(escapeMatches) == DataFrame
            escapeMatches |> SQLite.load!(dbOutput,"single_locus_escape_matches",ifnotexists=true)
            vAbunds = []
            bAbunds = []
            for spacer in unique(escapeMatches.spacer_id)
                bstrains1 = unique(escapeMatches[escapeMatches.spacer_id .== spacer,:].bstrain_id)
                vstrains1 = unique(escapeMatches[escapeMatches.spacer_id .== spacer,:].vstrain_id)
                vAbund = sum([abund for (abund,) in
                    execute(dbTempSim, "SELECT abundance FROM vabundance
                    WHERE t = $(time) AND vstrain_id in
                    ($(join(vstrains1,", ")))")])
                append!(vAbunds,vAbund)
                bAbund = sum([abund for (abund,) in
                    execute(dbTempSim, "SELECT abundance FROM babundance
                    WHERE t = $(time) AND bstrain_id in
                    ($(join(bstrains1,", ")))")])
                append!(bAbunds,bAbund)
                execute(dbOutput, "INSERT INTO single_locus_escape_match_abundances VALUES (?,?,?,?)",
                (time, Int(spacer), vAbund, bAbund))
            end
            vRelAbunds = vAbunds/sum(vAbunds)
            bRelAbunds = bAbunds/sum(bAbunds)
            vLociShannon = exp(-1*sum(vRelAbunds.*(log.(vRelAbunds))))
            bLociShannon = exp(-1*sum(bRelAbunds.*(log.(bRelAbunds))))
            execute(dbOutput, "INSERT INTO single_locus_escape_match_diversity VALUES (?,?,?)",
            (time, vLociShannon, bLociShannon))
        end


    end
end

function matchAbundances()
    for (time,) in execute(dbTempSim, "SELECT DISTINCT t FROM babundance ORDER BY t")
        println("Computing match abundances at time $(time)")
        v0Abund  = let
            vstrains0 = [vstrain_id for (vstrain_id,) in
                execute(dbTempMatch, "SELECT DISTINCT vstrain_id FROM bstrain_to_vstrain_0matches
                WHERE t = $(time) ORDER BY vstrain_id")]
            if length(vstrains0) > 0
                sum([abund for (abund,) in
                    execute(dbTempSim, "SELECT abundance FROM vabundance
                    WHERE t = $(time) AND vstrain_id in ($(join(vstrains0,", ")))")])
            else
                0
            end
        end
        v1Abund  = let
            vstrains1 = [vstrain_id for (vstrain_id,) in
                execute(dbTempMatch, "SELECT DISTINCT vstrain_id FROM bstrain_to_vstrain_matches
                WHERE t = $(time) AND match_length = 1
                ORDER BY time_specific_match_id")]
            if length(vstrains1) > 0
                sum([abund for (abund,) in
                    execute(dbTempSim, "SELECT abundance FROM vabundance
                    WHERE t = $(time) AND vstrain_id in ($(join(vstrains1,", ")))")])
            else
                0
            end
        end
        v2Abund  = let
            vstrains2 = [vstrain_id for (vstrain_id,) in
                execute(dbTempMatch, "SELECT DISTINCT vstrain_id FROM bstrain_to_vstrain_matches
                WHERE t = $(time) AND match_length = 2
                ORDER BY time_specific_match_id")]
            if length(vstrains2) > 0
                sum([abund for (abund,) in
                    execute(dbTempSim, "SELECT abundance FROM vabundance
                    WHERE t = $(time) AND vstrain_id in ($(join(vstrains2,", ")))")])
            else
                0
            end
        end
        v3Abund  = let
            vstrains3 = [vstrain_id for (vstrain_id,) in
                execute(dbTempMatch, "SELECT DISTINCT vstrain_id FROM bstrain_to_vstrain_matches
                WHERE t = $(time) AND match_length = 3
                ORDER BY time_specific_match_id")]
            if length(vstrains3) > 0
                sum([abund for (abund,) in
                    execute(dbTempSim, "SELECT abundance FROM vabundance
                    WHERE t = $(time) AND vstrain_id in ($(join(vstrains3,", ")))")])
            else
                0
            end
        end
        execute(dbOutput, "INSERT INTO vMatchAbundances VALUES (?,?,?,?,?)",
        (time, v0Abund, v1Abund, v2Abund, v3Abund))


        b0Abund  = let
            bstrains0 = [bstrain_id for (bstrain_id,) in
                execute(dbTempMatch, "SELECT DISTINCT bstrain_id FROM bstrain_to_vstrain_0matches
                WHERE t = $(time) ORDER BY bstrain_id")]
            if length(bstrains0) > 0
                sum([abund for (abund,) in
                    execute(dbTempSim, "SELECT abundance FROM babundance
                    WHERE t = $(time) AND bstrain_id in ($(join(bstrains0,", ")))")])
            else
                0
            end
        end
        b1Abund  = let
            bstrains1 = [bstrain_id for (bstrain_id,) in
                execute(dbTempMatch, "SELECT DISTINCT bstrain_id FROM bstrain_to_vstrain_matches
                WHERE t = $(time) AND match_length = 1
                ORDER BY time_specific_match_id")]
            if length(bstrains1) > 0
                sum([abund for (abund,) in
                    execute(dbTempSim, "SELECT abundance FROM babundance
                    WHERE t = $(time) AND bstrain_id in ($(join(bstrains1,", ")))")])
            else
                0
            end
        end
        b2Abund  = let
            bstrains2 = [bstrain_id for (bstrain_id,) in
                execute(dbTempMatch, "SELECT DISTINCT bstrain_id FROM bstrain_to_vstrain_matches
                WHERE t = $(time) AND match_length = 2
                ORDER BY time_specific_match_id")]
            if length(bstrains2) > 0
                sum([abund for (abund,) in
                    execute(dbTempSim, "SELECT abundance FROM babundance
                    WHERE t = $(time) AND bstrain_id in ($(join(bstrains2,", ")))")])
            else
                0
            end
        end
        b3Abund  = let
            bstrains3 = [bstrain_id for (bstrain_id,) in
                execute(dbTempMatch, "SELECT DISTINCT bstrain_id FROM bstrain_to_vstrain_matches
                WHERE t = $(time) AND match_length = 3
                ORDER BY time_specific_match_id")]
            if length(bstrains3) > 0
                sum([abund for (abund,) in
                    execute(dbTempSim, "SELECT abundance FROM babundance
                    WHERE t = $(time) AND bstrain_id in ($(join(bstrains3,", ")))")])
            else
                0
            end
        end
        execute(dbOutput, "INSERT INTO bMatchAbundances VALUES (?,?,?,?,?)",
        (time, b0Abund, b1Abund, b2Abund, b3Abund))
    end
end


matchAbundances()
singleLocusMatches()


println("Complete!")
