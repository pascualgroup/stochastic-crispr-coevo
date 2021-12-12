#!/usr/bin/env julia

println("(Annoying Julia compilation delay...)")

using SQLite
using DataFrames
using SQLite.DBInterface: execute


run_id = ARGS[1]

## Define Paths ##
SCRIPT_PATH = abspath(dirname(PROGRAM_FILE))

dbSimPath = joinpath(SCRIPT_PATH,"..","..","simulation","sweep_db_gathered.sqlite") # cluster
dbMatchPath =  joinpath(SCRIPT_PATH,"..","gathered-analyses","matches","matches.sqlite") # cluster
dbCladePath =  joinpath(SCRIPT_PATH,"..","gathered-analyses","clade-abundances","clade-abundances.sqlite") # cluster
dbOutputPath = joinpath("match-abundances_output.sqlite") # cluster

# # dbSimPath = joinpath("/Volumes/Yadgah/sweep_db_gathered.sqlite") # local
# dbSimPath = joinpath("/Volumes/Yadgah","run_id1455_combo73_replicate15.sqlite") # local
# dbMatchPath =  joinpath("/Volumes/Yadgah/matches.sqlite") # local
# dbCladePath = joinpath("/Volumes/Yadgah/clade-abundances.sqlite") # local
# dbOutputPath = joinpath("/Volumes/Yadgah/match-abundances_output.sqlite") # local

if isfile(dbOutputPath)
    error("match-abundances_output.sqlite already exists; delete first")
end # cluster
##

# dbSim = SQLite.DB(dbSimPath)
dbOutput = SQLite.DB(dbOutputPath)

execute(dbOutput, "CREATE TABLE total_matched_bstrain_abundances (t REAL, bstrain_match_length INTEGER,
num_bstrains INTEGER, bfrequency REAL, babundance INTEGER,
num_vstrains_matched INTEGER, vfrequency_matched REAL, vabundance_matched INTEGER)")

execute(dbOutput, "CREATE TABLE bclades_match_type_abundances (t REAL, bclade_id INTEGER,
match_type INTEGER, matched_bfrequency REAL, matched_vfrequency REAL,
matched_babundance INTEGER, matched_vabundance INTEGER,
num_bstrains_matched INTEGER, num_vstrains_matched INTEGER)")

execute(dbOutput, "CREATE TABLE bmatch_type_abundances (t REAL, bmatch_type INTEGER,
matched_bfrequency REAL, matched_to_vfrequency REAL,
matched_babundance INTEGER, matched_to_vabundance INTEGER,
num_bstrains_matched INTEGER, num_vstrains_matched_to INTEGER)")

execute(dbOutput, "CREATE TABLE vmatch_type_abundances (t REAL, vmatch_type INTEGER,
matched_vfrequency REAL, matched_to_bfrequency REAL,
matched_vabundance INTEGER, matched_to_babundance INTEGER,
num_vstrains_matched INTEGER, num_bstrains_matched_to INTEGER)")

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

dbTempClades = SQLite.DB()
#dbTempClades = SQLite.DB("/Volumes/Yadgah/testClades.sqlite") # local
execute(dbTempClades, "CREATE TABLE babundances (t REAL, clade_id INTEGER, bstrain_id INTEGER, abundance INTEGER)")
execute(dbTempClades, "BEGIN TRANSACTION")
execute(dbTempClades,"ATTACH DATABASE '$(dbCladePath)' as dbClade")
execute(dbTempClades,"INSERT INTO babundances (t, clade_id, bstrain_id, abundance)
SELECT t, clade_id, bstrain_id, abundance FROM dbClade.babundances WHERE run_id = $(run_id);")
execute(dbTempClades, "COMMIT")
execute(dbTempClades, "BEGIN TRANSACTION")
execute(dbTempClades, "CREATE INDEX bclade_index ON babundances (t, clade_id, bstrain_id, abundance)")
execute(dbTempClades, "COMMIT")


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
                VALUES (?,?,?,?,?,?,?,?)",(time, numMatched,
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
                VALUES (?,?,?,?,?,?,?,?)",(time, numMatched,
                0, 0, 0,
                0, 0, 0))
                continue
            end

            if numMatched == 0
                execute(dbOutput,"INSERT INTO total_matched_bstrain_abundances
                VALUES (?,?,?,?,?,?,?,?)",(time, numMatched,
                length(strainsB), subAbundanceB/totalBabundance, subAbundanceB,
                0, 0, 0))
                continue
            end

            strainsVMatched = [strain for (strain,) in execute(dbTempMatches,
            "SELECT DISTINCT vstrain_id FROM bstrain_to_vstrain_matches
            WHERE bstrain_id in ($(join(strainsB,", "))) AND t = $(time)")]

            subAbundanceVMatched = sum([abundance for (abundance,) in execute(dbTempSim,
            "SELECT abundance FROM vabundance WHERE vstrain_id in ($(join(strainsVMatched,", ")))
            AND t = $(time)")])

            execute(dbOutput,"INSERT INTO total_matched_bstrain_abundances
            VALUES (?,?,?,?,?,?,?,?)",(time, numMatched,
            length(strainsB),subAbundanceB/totalBabundance, subAbundanceB,
            length(strainsVMatched), subAbundanceVMatched/totalVabundance, subAbundanceVMatched))

        end

    end
end

function matchAbundances()
    for (time,) in execute(dbTempSim, "SELECT DISTINCT t FROM summary ORDER BY t")
        println("Searching match types at time $(time)")
        # if time > 30.0
        #     return
        # end
        (btotal,) = execute(dbTempSim, "SELECT microbial_abundance FROM summary WHERE t = $(time)")
        btotal = btotal.microbial_abundance
        (vtotal,) = execute(dbTempSim, "SELECT viral_abundance FROM summary WHERE t = $(time)")
        vtotal = vtotal.viral_abundance
        vstrains = [vstrain for (vstrain,) in execute(dbTempSim,
        "SELECT DISTINCT vstrain_id FROM vabundance WHERE t = $(time)")]
        vstrainsMatched = [vstrain for (vstrain,) in execute(dbTempMatches,
        "SELECT DISTINCT vstrain_id FROM bstrain_to_vstrain_matches
        WHERE t = $(time)")]
        setdiff!(vstrains,vstrainsMatched)
        println(vstrains)
        zeroMatches(time,btotal,vtotal,vstrains)
        matchTypes(time,btotal,vtotal)
    end
end

function zeroMatches(time,btotal,vtotal,vstrains)
    totalNumVstrainsMatched = length(vstrains)
    vstrainAbundances = [abundance for (abundance,) in execute(dbTempSim,
    "SELECT abundance FROM vabundance WHERE t = $(time)
    AND vstrain_id in ($(join(vstrains,", ")))")]
    if totalNumVstrainsMatched == 0
        totalvstrainAbundances = 0
    else
        totalvstrainAbundances = sum(vstrainAbundances)
    end

    totalbstrainAbundances = 0
    totalNumBstrainsMatched = 0
    for (clade_id,) in execute(dbTempClades,"SELECT DISTINCT clade_id
        FROM babundances WHERE t = $(time) ORDER BY clade_id")
        bstrains = [bstrain for (bstrain,) in execute(dbTempClades,
        "SELECT DISTINCT bstrain_id FROM babundances WHERE t = $(time)
        AND clade_id = $(clade_id)")]
        bstrainsMatched = [bstrain for (bstrain,) in execute(dbTempMatches,
        "SELECT DISTINCT bstrain_id FROM bstrain_to_vstrain_matches
        WHERE t = $(time)")]
        setdiff!(bstrains,bstrainsMatched)
        numBstrains = length(bstrains)
        if numBstrains == 0
            execute(dbOutput, "INSERT INTO bclades_match_type_abundances
            VALUES (?,?,?,?,?,?,?,?,?)",
            (time, clade_id, 0,
            0, totalvstrainAbundances/vtotal,
            0, totalvstrainAbundances,
            numBstrains, totalNumVstrainsMatched))
            continue
        end
        totalNumBstrainsMatched += numBstrains
        bstrainAbundances = [abundance for (abundance,) in execute(dbTempSim,
        "SELECT abundance FROM babundance WHERE t = $(time)
        AND bstrain_id in ($(join(bstrains,", ")))")]
        bstrainTotal = sum(bstrainAbundances)
        totalbstrainAbundances += bstrainTotal

        execute(dbOutput, "INSERT INTO bclades_match_type_abundances
        VALUES (?,?,?,?,?,?,?,?,?)",
        (time, clade_id, 0,
        bstrainTotal/btotal, totalvstrainAbundances/vtotal,
        bstrainTotal, totalvstrainAbundances,
        numBstrains, totalNumVstrainsMatched))
    end
    execute(dbOutput, "INSERT INTO match_type_abundances
    VALUES (?,?,?,?,?,?,?,?)",
    (time, 0,
    totalbstrainAbundances/btotal, totalvstrainAbundances/vtotal,
    totalbstrainAbundances, totalvstrainAbundances,
    totalNumBstrainsMatched, totalNumVstrainsMatched))
end

function matchTypes(time,btotal,vtotal)
    for (match_type,) in execute(dbTempMatches,"SELECT DISTINCT match_length
        FROM bstrain_to_vstrain_matches WHERE t = $(time) ORDER BY match_length")
        totalbstrainAbundances = 0
        totalNumBstrainsMatched = 0
        allVstrains = []
        allBstrains = []
        for (clade_id,) in execute(dbTempClades,"SELECT DISTINCT clade_id
            FROM babundances WHERE t = $(time) ORDER BY clade_id")
            bstrains = [bstrain for (bstrain,) in execute(dbTempClades,
            "SELECT DISTINCT bstrain_id FROM babundances WHERE t = $(time)
            AND clade_id = $(clade_id)")]
            bstrainsMatchType = [bstrain for (bstrain,) in execute(dbTempMatches,
            "SELECT DISTINCT bstrain_id FROM bstrain_to_vstrain_matches
            WHERE t = $(time)
            AND match_length = $(match_type)")]
            intersect!(bstrains,bstrainsMatchType)
            if length(bstrains) == 0
                continue
            end
            numBstrains = length(bstrains)
            totalNumBstrainsMatched += numBstrains
            bstrainAbundances = [abundance for (abundance,) in
            execute(dbTempSim, "SELECT abundance FROM babundance
            WHERE t = $(time)
            AND bstrain_id in ($(join(bstrains,", ")))")]
            bstrainTotal = sum(bstrainAbundances)
            totalbstrainAbundances += bstrainTotal

            vstrains = [vstrain for (vstrain,) in execute(dbTempMatches,
            "SELECT DISTINCT vstrain_id FROM bstrain_to_vstrain_matches
            WHERE t = $(time)
            AND match_length = $(match_type)
            AND bstrain_id in ($(join(bstrains,", ")))")]
            union!(allVstrains,vstrains)
            numVstrains = length(vstrains)
            totalNumVstrainsMatched += numVstrains
            vstrainAbundances = [abundance for (abundance,) in
            execute(dbTempSim, "SELECT abundance FROM vabundance
            WHERE t = $(time)
            AND vstrain_id in ($(join(vstrains,", ")))")]
            vstrainTotal = sum(vstrainAbundances)

            execute(dbOutput, "INSERT INTO bclades_match_type_abundances
            VALUES (?,?,?,?,?,?,?,?,?)",
            (time, clade_id, match_type,
            bstrainTotal/btotal, vstrainTotal/vtotal,
            bstrainTotal, vstrainTotal,
            numBstrains, numVstrains))
        end

        totalvstrainAbundances = [abundances for (abundances,)
        in execute(dbTempSim, "SELECT abundance FROM vabundance WHERE t = $(time)
        AND vstrain_id in ($(join(allVstrains,", ")))")]

        execute(dbOutput, "INSERT INTO bmatch_type_abundances
        VALUES (?,?,?,?,?,?,?,?)",
        (time, match_type,
        totalbstrainAbundances/btotal, totalvstrainAbundances/vtotal,
        totalbstrainAbundances, totalvstrainAbundances,
        totalNumBstrainsMatched, length(allVstrains)))

        vstrains = [vstrain for (vstrain,) in execute(dbTempSim,
        "SELECT DISTINCT vstrain_id FROM vabundances WHERE t = $(time)")]
        vstrainsMatchType = [vstrain for (vstrain,) in execute(dbTempMatches,
        "SELECT DISTINCT vstrain_id FROM bstrain_to_vstrain_matches
        WHERE t = $(time)
        AND match_length = $(match_type)")]
        intersect!(vstrains,vstrainsMatchType)

        if length(vstrains) == 0
            ?
        end
        execute(dbOutput, "INSERT INTO match_type_abundances
        VALUES (?,?,?,?,?,?,?,?)",
        (time, match_type,
        totalbstrainAbundances/btotal, totalvstrainAbundances/vtotal,
        totalbstrainAbundances, totalvstrainAbundances,
        totalNumBstrainsMatched, totalNumVstrainsMatched))
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

# total_matchedBstrains()
matchAbundances()

# matchedBstrains()
#
# matchedVstrains()
# total_matchedVstrains()


println("Complete!")
