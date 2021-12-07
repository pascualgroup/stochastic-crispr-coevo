#!/usr/bin/env julia

println("(Annoying Julia compilation delay...)")

using SQLite
using DataFrames
using SQLite.DBInterface: execute


run_id = ARGS[1]

## Define Paths ##
SCRIPT_PATH = abspath(dirname(PROGRAM_FILE))

dbSimPath = joinpath(SCRIPT_PATH,"..","..","simulation","sweep_db_gathered.sqlite") # cluster
dbMatchPath =  joinpath(SCRIPT_PATH,"..","gathered-analyses","matches","matches.sqlite") # local
dbOutputPath = joinpath("match-abundances_output.sqlite") # cluster

# # dbSimPath = joinpath("/Volumes/Yadgah/sweep_db_gathered.sqlite") # local
# dbSimPath = joinpath("/Volumes/Yadgah","run_id1455_combo73_replicate15.sqlite") # local
# dbMatchPath =  joinpath("/Volumes/Yadgah/matches.sqlite") # local
# dbOutputPath = joinpath("/Volumes/Yadgah/match-abundances_output.sqlite") # local

if isfile(dbOutputPath)
    error("match-abundances_output.sqlite already exists; delete first")
end # cluster
##

dbSim = SQLite.DB(dbSimPath)
dbOutput = SQLite.DB(dbOutputPath)

execute(dbOutput, "CREATE TABLE total_matched_bstrain_abundances (t REAL, bstrain_match_length INTEGER,
num_bstrains INTEGER, bfrequency REAL, babundance INTEGER,
num_viral0matches INTEGER, vfrequency_0matches REAL, vabundance_0matches INTEGER,
num_vstrains_matched INTEGER, vfrequency_matched REAL, vabundance_matched INTEGER)")

# execute(dbOutput, "CREATE TABLE match_type_abundances (t REAL, match_type INTEGER,
# num_match_ids INTEGER, match_type_richness INTEGER,
# matched_bfrequency REAL, matched_vfrequency REAL,
# matched_babundance INTEGER, matched_vabundance INTEGER,
# num_bstrains_matched INTEGER, num_vstrains_matched INTEGER)")

# execute(dbOutput, "CREATE TABLE matched_bstrain_abundances (t REAL, bstrain_match_length INTEGER,
# bstrain_id INTEGER, bfrequency REAL, babundance INTEGER,
# num_viral0matches INTEGER, vfrequency_0matches REAL, vabundance_0matches INTEGER,
# num_vstrains_matched INTEGER, vfrequency_matched REAL, vabundance_matched INTEGER)")

# execute(dbOutput, "CREATE TABLE matched_vstrain_abundances (t REAL, vstrain_match_length INTEGER,
# vstrain_id INTEGER, vfrequency REAL, vabundance INTEGER,
# num_microbial0matches INTEGER, bfrequency_0matches REAL, babundance_0matches INTEGER,
# num_bstrains_matched INTEGER, bfrequency_matched REAL, babundance_matched INTEGER)")

# execute(dbOutput, "CREATE TABLE total_matched_vstrain_abundances (t REAL, vstrain_match_length INTEGER,
# num_vstrains INTEGER, vfrequency REAL, vabundance INTEGER,
# num_microbial0matches INTEGER, bfrequency_0matches REAL, babundance_0matches INTEGER,
# num_bstrains_matched INTEGER, bfrequency_matched REAL, babundance_matched INTEGER)")

# Create temporary database that is a copy of the main database at the run_id value of the script's argument
dbTempSim = SQLite.DB()
#dbTempSim = SQLite.DB("/Volumes/Yadgah/testSim.sqlite") # local
execute(dbTempSim, "CREATE TABLE summary (t REAL, microbial_abundance INTEGER, viral_abundance INTEGER)")
execute(dbTempSim, "CREATE TABLE babundance (t REAL, bstrain_id INTEGER, abundance INTEGER)")
execute(dbTempSim, "CREATE TABLE vabundance (t REAL, vstrain_id INTEGER, abundance INTEGER)")

execute(dbTempSim, "BEGIN TRANSACTION")
execute(dbTempSim,"ATTACH DATABASE '$(dbSimPath)' as dbSim")
execute(dbTempSim,"INSERT INTO summary(t, microbial_abundance,viral_abundance) SELECT t, microbial_abundance,viral_abundance FROM dbSim.summary WHERE run_id = $(run_id);")
execute(dbTempSim,"INSERT INTO babundance (t, bstrain_id, abundance) SELECT t, bstrain_id, abundance FROM dbSim.babundance WHERE run_id = $(run_id);")
execute(dbTempSim,"INSERT INTO vabundance (t, vstrain_id, abundance) SELECT t, vstrain_id, abundance FROM dbSim.vabundance WHERE run_id = $(run_id);")
execute(dbTempSim, "COMMIT")


execute(dbTempSim, "BEGIN TRANSACTION")
execute(dbTempSim, "CREATE INDEX summary_index ON summary (t,microbial_abundance,viral_abundance)")
execute(dbTempSim, "CREATE INDEX babundance_index ON babundance (t,bstrain_id,abundance)")
execute(dbTempSim, "CREATE INDEX vabundance_index ON vabundance (t,vstrain_id)")
execute(dbTempSim, "COMMIT")

dbTempMatches = SQLite.DB()
#dbTempMatches = SQLite.DB("/Volumes/Yadgah/testMatches.sqlite") # local
execute(dbTempMatches, "CREATE TABLE bstrain_num_matched_spacers (t REAL, bstrain_id INTEGER, num_matched_spacers INTEGER)")
execute(dbTempMatches, "CREATE TABLE bstrain_to_vstrain_matches (t REAL, bstrain_id INTEGER, vstrain_id INTEGER, match_length INTEGER)")

execute(dbTempMatches, "BEGIN TRANSACTION")
execute(dbTempMatches,"ATTACH DATABASE '$(dbMatchPath)' as dbMatch")
execute(dbTempMatches,"INSERT INTO bstrain_num_matched_spacers (t, bstrain_id, num_matched_spacers)
SELECT t, bstrain_id, num_matched_spacers FROM dbMatch.bstrain_num_matched_spacers WHERE run_id = $(run_id);")
execute(dbTempMatches,"INSERT INTO bstrain_to_vstrain_matches (t, bstrain_id, vstrain_id, match_length)
SELECT t, bstrain_id, vstrain_id, match_length FROM dbMatch.bstrain_to_vstrain_matches WHERE run_id = $(run_id);")
execute(dbTempMatches, "COMMIT")


execute(dbTempMatches, "BEGIN TRANSACTION")
execute(dbTempMatches, "CREATE INDEX strain_to_strain_index ON bstrain_to_vstrain_matches (t, bstrain_id, vstrain_id, match_length)")
execute(dbTempMatches, "CREATE INDEX vabundance_index ON bstrain_num_matched_spacers (t, bstrain_id, num_matched_spacers)")
execute(dbTempMatches, "COMMIT")

println("Processing bstrain match abundances of run $(run_id)")

# function matchedBstrains()
#     for (time,) in execute(dbTempSim, "SELECT t FROM summary ORDER BY t")
#         println("Searching matched bstrains at time $(time)")
#         (totalBabundance,) = execute(dbTempSim, "SELECT microbial_abundance
#         FROM summary WHERE t = $(time)")
#         (totalVabundance,) = execute(dbTempSim, "SELECT viral_abundance
#         FROM summary WHERE t = $(time)")
#         totalBabundance = totalBabundance.microbial_abundance
#         totalVabundance = totalVabundance.viral_abundance
#         for (numMatched,) in execute(dbTempMatches,"SELECT DISTINCT num_matched_spacers
#             FROM bstrain_num_matched_spacers WHERE t = $(time)")
#             strains = [strain for (strain,) in execute(dbTempMatches,"SELECT bstrain_id
#             FROM bstrain_num_matched_spacers WHERE t = $(time)
#             AND num_matched_spacers = $(numMatched)")]
#             subAbundance = sum([abundance for (abundance,) in execute(dbTempMatches,
#             "SELECT abundance FROM babundance WHERE bstrain_id in ($(join(strains,", ")))
#             AND t = $(time)")])
#
#             execute(dbOutput,"INSERT INTO matched_bstrain_abundances
#             VALUES (?,?,?,?,?,?,?)",(time, numMatched,subAbundance/totalBabundance,
#             subAbundance, num_vstrains_matched, vfrequency_matched,
#             vabundance_matched))
#
#         end
#
#     end
# end
#
# function matchedVstrains()
#     for (time,) in execute(dbTempSim, "SELECT t FROM summary ORDER BY t")
#         println("Searching matched vstrains at time $(time)")
#         (totalAbundance,) = execute(dbTempSim, "SELECT microbial_abundance
#         FROM summary WHERE t = $(time)")
#         (totalVabundance,) = execute(dbTempSim, "SELECT viral_abundance
#         FROM summary WHERE t = $(time)")
#         totalBabundance = totalBabundance.microbial_abundance
#         totalVabundance = totalVabundance.viral_abundance
#         for (numMatched,) in execute(dbMatches,"SELECT DISTINCT num_matched_spacers
#             FROM bstrain_num_matched_spacers WHERE t = $(time)")
#             strains = [strain for (strain,) in execute(dbTempMatches,"SELECT bstrain_id
#             FROM bstrain_num_matched_spacers WHERE t = $(time)
#             AND num_matched_spacers = $(numMatched)")]
#             subAbundance = sum([abundance for (abundance,) in execute(dbTempMatches,
#             "SELECT abundance FROM babundance WHERE bstrain_id in ($(join(strains,", ")))
#             AND t = $(time)")])
#             execute(dbOutput,"INSERT INTO matched_bstrain_abundances
#             VALUES (?,?,?,?,?,?,?,?,?,?,?)",(time, numMatched,subAbundance/totalAbundance,
#             subAbundance, num_vstrains_matched, vfrequency_matched,
#             vabundance_matched))
#
#         end
#
#     end
# end

function total_matchedBstrains()
    for (time,) in execute(dbTempSim, "SELECT t FROM summary ORDER BY t")
        # if time > 30.0
        #     return
        # end
        println("Searching total matched bstrains at time $(time)")
        (totalBabundance,) = execute(dbTempSim, "SELECT microbial_abundance
        FROM summary WHERE t = $(time)")
        (totalVabundance,) = execute(dbTempSim, "SELECT viral_abundance
        FROM summary WHERE t = $(time)")
        totalBabundance = totalBabundance.microbial_abundance
        totalVabundance = totalVabundance.viral_abundance
        for (numMatched,) in execute(dbTempMatches,"SELECT DISTINCT num_matched_spacers
            FROM bstrain_num_matched_spacers WHERE t = $(time)")
            strainsB = [strain for (strain,) in execute(dbTempMatches,"SELECT bstrain_id
            FROM bstrain_num_matched_spacers WHERE t = $(time)
            AND num_matched_spacers = $(numMatched)")]
            if length(strainsB) == 0
                numVstrains = length([strain for (strain,) in execute(dbTempSim,
                "SELECT DISTINCT vstrain_id FROM vabaundance WHERE t = $(time)")])
                execute(dbOutput,"INSERT INTO total_matched_bstrain_abundances
                VALUES (?,?,?,?,?,?,?,?,?,?,?)",(time, numMatched,
                0, 0, 0,
                0, 0, 0,
                0, 0, 0))
                continue
            end

            subAbundanceB = sum([abundance for (abundance,) in execute(dbTempSim,
            "SELECT abundance FROM babundance WHERE bstrain_id in ($(join(strainsB,", ")))
            AND t = $(time)")])

            strainsVAll = [strain for (strain,) in execute(dbTempSim,
            "SELECT DISTINCT vstrain_id FROM vabundance WHERE t = $(time)")]

            if numMatched == 0 && subAbundanceB == 0
                execute(dbOutput,"INSERT INTO total_matched_bstrain_abundances
                VALUES (?,?,?,?,?,?,?,?,?,?,?)",(time, numMatched,
                0, 0, 0,
                0, 0, 0,
                0, 0, 0))
                continue
            end

            if numMatched == 0
                execute(dbOutput,"INSERT INTO total_matched_bstrain_abundances
                VALUES (?,?,?,?,?,?,?,?,?,?,?)",(time, numMatched,
                length(strainsB), subAbundanceB/totalBabundance, subAbundanceB,
                length(strainsVAll), 1, totalVabundance,
                0, 0, 0))
                continue
            end

            strainsVMatched = [strain for (strain,) in execute(dbTempMatches,
            "SELECT DISTINCT vstrain_id FROM bstrain_to_vstrain_matches
            WHERE bstrain_id in ($(join(strainsB,", "))) AND t = $(time)")]

            strainsV0 = setdiff(strainsVAll,strainsVMatched)

            subAbundanceVMatched = sum([abundance for (abundance,) in execute(dbTempSim,
            "SELECT abundance FROM vabundance WHERE vstrain_id in ($(join(strainsVMatched,", ")))
            AND t = $(time)")])

            if length(strainsV0) == 0
                execute(dbOutput,"INSERT INTO total_matched_bstrain_abundances
                VALUES (?,?,?,?,?,?,?,?,?,?,?)",(time, numMatched,
                length(strainsB),subAbundanceB/totalBabundance, subAbundanceB,
                0, 0, 0,
                length(strainsVMatched), subAbundanceVMatched/totalVabundance, subAbundanceVMatched))
                continue
            end

            subAbundanceV0 = sum([abundance for (abundance,) in execute(dbTempSim,
            "SELECT abundance FROM vabundance WHERE vstrain_id in ($(join(strainsV0,", ")))
            AND t = $(time)")])

            execute(dbOutput,"INSERT INTO total_matched_bstrain_abundances
            VALUES (?,?,?,?,?,?,?,?,?,?,?)",(time, numMatched,
            length(strainsB),subAbundanceB/totalBabundance, subAbundanceB,
            length(strainsV0), subAbundanceV0/totalVabundance, subAbundanceV0,
            length(strainsVMatched), subAbundanceVMatched/totalVabundance, subAbundanceVMatched))

        end

    end
end

# function total_matchedVstrains()
#     for (time,) in execute(dbTempSim, "SELECT t FROM summary ORDER BY t")
#         println("Searching total matched vstrains at time $(time)")
#         (totalBabundance,) = execute(dbTempSim, "SELECT microbial_abundance
#         FROM summary WHERE t = $(time)")
#         (totalVabundance,) = execute(dbTempSim, "SELECT viral_abundance
#         FROM summary WHERE t = $(time)")
#         totalBabundance = totalBabundance.microbial_abundance
#         totalVabundance = totalVabundance.viral_abundance
#         for (numMatched,) in execute(dbMatches,"SELECT DISTINCT num_matched_spacers
#             FROM bstrain_num_matched_spacers WHERE t = $(time)")
#             strains = [strain for (strain,) in execute(dbTempMatches,"SELECT bstrain_id
#             FROM bstrain_num_matched_spacers WHERE t = $(time)
#             AND num_matched_spacers = $(numMatched)")]
#             if length(strains) == 0
#                 execute(dbOutput,"INSERT INTO matched_bstrain_abundances
#                 VALUES (?,?,?,?,?,?,?)",(time, numMatched, 0, 0, 0, 0, 0))
#                 continue
#             end
#
#             subAbundance = sum([abundance for (abundance,) in execute(dbTempMatches,
#             "SELECT abundance FROM babundance WHERE bstrain_id in ($(join(strains,", ")))
#             AND t = $(time)")])
#
#             strainsV = [strain for (strain,) in execute(dbTempMatches,
#             "SELECT vstrain_id FROM bstrain_to_vstrain_matches
#             WHERE bstrain_id in ($(join(strains,", "))) AND t = $(time)")]
#
#             subAbundanceV = sum([abundance for (abundance,) in execute(dbTempMatches,
#             "SELECT abundance FROM vabundance WHERE vstrain_id in ($(join(strains2,", ")))
#             AND t = $(time)")])
#
#             execute(dbOutput,"INSERT INTO matched_bstrain_abundances
#             VALUES (?,?,?,?,?,?,?)",(time, numMatched,subAbundance/totalAbundance,
#             subAbundance, length(strainsV), subAbundanceV/totalVabundance,
#             subAbundanceV))
#
#         end
#
#     end
# end

# function matchAbundances()
#     for (time,) in execute(dbTempSim, "SELECT DISTINCT t FROM summary ORDER BY t")
#         # println("time is $(time)")
#         if time > 30.0
#             return
#         end
#         (btotal,) = execute(dbTempSim, "SELECT microbial_abundance FROM summary WHERE t = $(time)")
#         btotal = btotal.microbial_abundance
#         (vtotal,) = execute(dbTempSim, "SELECT viral_abundance FROM summary WHERE t = $(time)")
#         vtotal = vtotal.viral_abundance
#         totalMatches = length([match for (match,) in execute(dbOutput, "SELECT DISTINCT match_id FROM match_types WHERE t = $(time)")])
#         matchDivTotal = 0
#         for (match_type,) in execute(dbOutput,"SELECT DISTINCT match_type FROM strain_matches WHERE t = $(time) ORDER BY match_type")
#             bstrains = [bstrain for (bstrain,) in execute(dbOutput, "SELECT DISTINCT bstrain_id FROM strain_matches WHERE t = $(time)
#             AND match_type = $(match_type)")]
#             bstrainAbundances = [abundance for (abundance,) in execute(dbTempSim, "SELECT abundance FROM babundance WHERE t = $(time)
#             AND bstrain_id in ($(join(bstrains,", ")))")]
#             bstrainTotal = sum(bstrainAbundances)
#             vstrains = [vstrain for (vstrain,) in execute(dbOutput, "SELECT DISTINCT vstrain_id FROM strain_matches WHERE t = $(time)
#             AND match_type = $(match_type)")]
#             vstrainAbundances = [abundance for (abundance,) in execute(dbTempSim, "SELECT abundance FROM vabundance WHERE t = $(time)
#             AND vstrain_id in ($(join(vstrains,", ")))")]
#             vstrainTotal = sum(vstrainAbundances)
#
#
#             (first,) = execute(dbOutput, "SELECT DISTINCT match_id FROM match_types WHERE t = $(time)
#                 AND match_type = $(match_type) ORDER BY match_id")
#             match_ref = [spacer for (spacer,) in
#                 execute(dbOutput, "SELECT spacer_id FROM match_types WHERE t = $(time) AND match_id = $(first.match_id)
#                 ORDER BY spacer_id")]
#                 matchDiv = 1
#                 matchDivTotal += 1
#                 #if time == 19.0
#                     #println("reference matches: $(match_ref)")
#                 #end
#             for (match_id,) in execute(dbOutput, "SELECT DISTINCT match_id FROM match_types WHERE t = $(time)
#                 AND match_type = $(match_type) ORDER BY match_id")
#                 matches = [spacer for (spacer,) in
#                     execute(dbOutput, "SELECT spacer_id FROM match_types WHERE t = $(time) AND match_id = $(match_id)
#                     ORDER BY spacer_id")]
#                     #if time == 19.0
#                         #println("matches: $(matches)")
#                         #println("intersection: $(intersect(matches,match_ref))")
#                         #println("length: $(length(intersect(matches,match_ref)))")
#                         #println("match_type: $(match_type)")
#                         #println("truth: $(length(intersect(matches,match_ref)) != match_type)")
#                     #end
#
#                 if length(intersect(matches,match_ref)) !== match_type
#                     matchDiv += 1
#                     matchDivTotal += 1
#                     match_ref = [spacer for (spacer,) in
#                         execute(dbOutput, "SELECT spacer_id FROM match_types WHERE t = $(time) AND match_id = $(match_id)
#                         ORDER BY spacer_id")]
#                 end
#             end
#
#             allSpacers = [spacer for spacer in execute(dbOutput, "SELECT spacer_id FROM match_types
#             WHERE t = $(time) AND match_type = $(match_type)")]
#
#             match_ids = [match_id for match_id in execute(dbOutput, "SELECT DISTINCT match_id FROM match_types WHERE t = $(time)
#                 AND match_type = $(match_type)")]
#             execute(dbOutput, "INSERT INTO match_abundances VALUES (?,?,?,?,?,?,?,?,?,?,?,?)",
#             (time, match_type, length(match_ids), length(match_ids)/totalMatches, matchDiv,
#             length(unique(allSpacers)), bstrainTotal/btotal, vstrainTotal/vtotal,
#             bstrainTotal, vstrainTotal, length(bstrains), length(vstrains)))
#         end
#
#         execute(dbOutput, "INSERT INTO total_abundances VALUES (?,?,?,?,?)",(time, btotal, vtotal, totalMatches, matchDivTotal))
#     end
# end


total_matchedBstrains()

# matchedBstrains()
#
# matchedVstrains()
# total_matchedVstrains()
# # matchAbundances()

println("Complete!")
