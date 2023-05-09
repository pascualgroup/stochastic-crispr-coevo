#!/usr/bin/env julia

println("(Annoying Julia compilation delay...)")

using SQLite
using DataFrames
using SQLite.DBInterface: execute


run_id = ARGS[1]

## Define Paths ##
SCRIPT_PATH = abspath(dirname(PROGRAM_FILE))

#dbSimPath = joinpath(SCRIPT_PATH,"..","..","simulation","sweep_db_gathered.sqlite") # cluster
#dbSimInfoPath = joinpath(SCRIPT_PATH,"..","..","simulation","sweep_db.sqlite") # cluster
#dbOutputPath = joinpath("match-diversity_output.sqlite") # cluster

dbSimPath = joinpath("/Volumes/Yadgah/sweep_db_gathered.sqlite") # local
dbSimInfoPath = joinpath("/Volumes/Yadgah/sweep_db.sqlite") # local
dbOutputPath = joinpath("/Volumes/Yadgah/match-diversity_output.sqlite") # local

if isfile(dbOutputPath)
    error("match-diversity_output.sqlite already exists; delete first")
end # cluster
##

dbSimInfo = SQLite.DB(dbSimInfoPath)
dbSim = SQLite.DB(dbSimPath)
dbOutput = SQLite.DB(dbOutputPath)

execute(dbOutput, "CREATE TABLE total_abundances (t REAL, babundance INTEGER, vabundance INTEGER,
total_num_matches_w_redundancy INTEGER, total_match_richness INTEGER)")
execute(dbOutput, "CREATE TABLE strain_matches (t REAL, bstrain_id INTEGER, vstrain_id INTEGER,
match_id INTEGER, match_type INTEGER,
bstrain_frequency REAL, vstrain_frequency REAL,
bstrain_abundance INTEGER, vstrain_abundance INTEGER)")
execute(dbOutput, "CREATE TABLE match_types (t REAL, match_id INTEGER, spacer_id INTEGER, match_type INTEGER)")
execute(dbOutput, "CREATE TABLE IDI (t REAL, IDI REAL)")
execute(dbOutput, "CREATE TABLE PDI (t REAL, PDI REAL, maxPDI REAL, DRI REAL, PDI_degree1 REAL,
PDI_degree2 REAL, PDI_degree3 REAL, PDI_degree4 REAL, PDI_degree5 REAL)")
execute(dbOutput, "CREATE TABLE HVI (t REAL, HVI REAL)")


execute(dbOutput, "CREATE TABLE match_abundances (t REAL, match_type INTEGER,
num_match_ids INTEGER, proportion_of_total_matches_w_redundancy REAL, match_type_richness INTEGER, num_unique_spacers INTEGER,
matched_bfrequency REAL, matched_vfrequency REAL,
matched_babundance INTEGER, matched_vabundance INTEGER,
num_bstrains_matched INTEGER, num_vstrains_matched INTEGER)")

execute(dbOutput, "CREATE TABLE spacer_match_frequencies (t REAL, spacer_id INTEGER,
bfrequency_with_spacer_id REAL, vfrequency_with_spacer_id REAL,
babundance_with_spacer_id INTEGER, vabundance_with_spacer_id INTEGER)")


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
function analyzeMatches()
    for (time,) in execute(dbTemp, "SELECT DISTINCT t FROM summary ORDER BY t")
        # println("time is $(time)")
        if time > 30.0
            return
        end
        # println("Identifying matches and computing immunity indices at time $(time)")
        IDI = 0
        PDI = 0
        HVI = 0
        DRI = 0
        PDI_degree1 = 0 # should the order number be automated???
        PDI_degree2 = 0
        PDI_degree3 = 0
        PDI_degree4 = 0
        PDI_degree5 = 0
        match_id = 1
        (btotal,) = execute(dbTemp, "SELECT microbial_abundance FROM summary WHERE t = $(time)")
        btotal = btotal.microbial_abundance
        (vtotal,) = execute(dbTemp, "SELECT viral_abundance FROM summary WHERE t = $(time)")
        vtotal = vtotal.viral_abundance

        spacerFrequencies = DataFrame(t = Float64[], spacer_id = Int64[],
        bfrequency_with_spacer_id = Float64[], vfrequency_with_spacer_id = Float64[],
        babundance_with_spacer_id = Int64[], vabundance_with_spacer_id = Int64[])

        maxbAbundance = maximum([abundance for (abundance,) in
        execute(dbTemp, "SELECT abundance FROM babundance WHERE t = $(time) ORDER BY bstrain_id")])

        for (bstrain_id,) in execute(dbTemp, "SELECT bstrain_id FROM babundance WHERE t = $(time) ORDER BY bstrain_id")
            (babundance,) = execute(dbTemp, "SELECT abundance FROM babundance WHERE t = $(time)
            AND bstrain_id = $(bstrain_id)")
            babundance = babundance.abundance

            #bipartite
            for (vstrain_id,) in execute(dbTemp, "SELECT vstrain_id FROM vabundance WHERE t = $(time) ORDER BY vstrain_id")
                (vabundance,) = execute(dbTemp, "SELECT abundance FROM vabundance WHERE t = $(time)
                AND vstrain_id = $(vstrain_id)")
                vabundance = vabundance.abundance
                spacerMatches = [spacer_id for (spacer_id,) in
                execute(dbTemp, "SELECT spacer_id FROM bspacers WHERE bstrain_id = $(bstrain_id)")]
                pspacers = [pspacer_id for (pspacer_id,) in
                execute(dbTemp, "SELECT spacer_id FROM vpspacers WHERE vstrain_id = $(vstrain_id)")]

                if bstrain_id == 1 # naive hosts
                    match_type = 0
                else
                    intersect!(spacerMatches,pspacers)
                    match_type = length(spacerMatches) # number of matches
                end
                if match_type == 0
                    # println("Skipping iteration vstrain_id = $(vstrain_id)")
                    HVI += (babundance/btotal)*(vabundance/vtotal)
                    continue
                end
                # println("Skip check: iteration vstrain_id = $(vstrain_id)")

                execute(dbOutput, "BEGIN TRANSACTION")
                for spacer_id in spacerMatches
                    execute(dbOutput, "INSERT INTO match_types VALUES (?,?,?,?)",(time, match_id, spacer_id, match_type))
                end
                execute(dbOutput, "INSERT INTO strain_matches VALUES (?,?,?,?,?,?,?,?,?)",
                (time, bstrain_id, vstrain_id, match_id, match_type,
                babundance/btotal,vabundance/vtotal,babundance,vabundance))
                execute(dbOutput, "COMMIT")

                match_id += 1
                IDI += (babundance/btotal)*(vabundance/vtotal)*match_type

                #tripartite
                dPDI, dDRI, dPDI_degree1, dPDI_degree2, dPDI_degree3, dPDI_degree4, dPDI_degree5 =
                compute_deltaPDI(time,bstrain_id,vstrain_id,babundance,vabundance,btotal,vtotal,spacerMatches,pspacers,maxbAbundance)
                PDI += dPDI
                DRI += dDRI
                PDI_degree1 += dPDI_degree1
                PDI_degree2 += dPDI_degree2
                PDI_degree3 += dPDI_degree3
                PDI_degree4 += dPDI_degree4
                PDI_degree5 += dPDI_degree5
            end
        end


        indSpacers = []
        for (match_id,) in execute(dbOutput, "SELECT match_id FROM strain_matches WHERE t = $(time)")
            matches = [spacer for (spacer,) in
                execute(dbOutput, "SELECT spacer_id FROM match_types WHERE t = $(time) AND match_id = $(match_id)")]
            push!(indSpacers,matches...)
            unique!(indSpacers)
        end

        execute(dbOutput, "BEGIN TRANSACTION")
        for spacer_id in indSpacers
            babundance = 0
            vabundance = 0
            # I CAN MAKE THIS MORE EFFICIENT
            for (bstrain_id,) in execute(dbOutput, "SELECT DISTINCT bstrain_id FROM strain_matches WHERE t = $(time)")
                spacers = [spacer for (spacer,) in
                execute(dbTemp, "SELECT spacer_id FROM bspacers WHERE bstrain_id = $(bstrain_id)")]
                if issubset(spacer_id,spacers)
                    (bStrainAbundance,) = execute(dbTemp, "SELECT abundance FROM babundance WHERE t = $(time)
                    AND bstrain_id = $(bstrain_id)")
                    babundance += bStrainAbundance.abundance
                end
            end
            for (vstrain_id,) in execute(dbOutput, "SELECT DISTINCT vstrain_id FROM strain_matches WHERE t = $(time)")
                pspacers = [pspacer for (pspacer,) in
                execute(dbTemp, "SELECT spacer_id FROM vpspacers WHERE vstrain_id = $(vstrain_id)")]
                if issubset(spacer_id,pspacers)
                    (vStrainAbundance,) = execute(dbTemp, "SELECT abundance FROM vabundance WHERE t = $(time)
                    AND vstrain_id = $(vstrain_id)")
                    vabundance += vStrainAbundance.abundance
                end
            end
            execute(dbOutput, "INSERT INTO spacer_match_frequencies VALUES (?,?,?,?,?,?)",
            (time, spacer_id, babundance/btotal, vabundance/vtotal, babundance,vabundance))
        end
        execute(dbOutput, "COMMIT")

        #spacerFrequencies |> SQLite.load!(dbOutput,"spacer_match_frequencies",ifnotexists=true)

        #updateSingleSpacerFrequencies!(time,spacerFrequencies,spacer_id,babundance,vabundance,btotal,vtotal)

        PDI = PDI/2
        DRI = DRI/2
        PDI_degree1 = PDI_degree1/2
        PDI_degree2 = PDI_degree2/2
        PDI_degree3 = PDI_degree3/2
        PDI_degree4 = PDI_degree4/2
        PDI_degree5 = PDI_degree5/2

        num_bstrains = length([bstrain_id for (bstrain_id,)
        in execute(dbTemp, "SELECT bstrain_id FROM babundance WHERE t = $(time)")])

        execute(dbOutput, "BEGIN TRANSACTION")
        execute(dbOutput, "INSERT INTO IDI VALUES (?,?)",(time, IDI))
        execute(dbOutput, "INSERT INTO PDI VALUES (?,?,?,?,?,?,?,?,?)",(time, PDI, 1-1/num_bstrains,
        DRI, PDI_degree1, PDI_degree2, PDI_degree3, PDI_degree4, PDI_degree5))
        execute(dbOutput, "INSERT INTO HVI VALUES (?,?)",(time, HVI))
        execute(dbOutput, "COMMIT")

        # println("Computing match diversities and abundances at time $(time)")
        matchAbundances(time,btotal,vtotal)
    end

end


function matchAbundances(time,btotal,vtotal)
    totalMatches = length([match for (match,) in execute(dbOutput, "SELECT DISTINCT match_id FROM match_types WHERE t = $(time)")])
    matchDivTotal = 0
    for (match_type,) in execute(dbOutput,"SELECT DISTINCT match_type FROM strain_matches WHERE t = $(time) ORDER BY match_type")
        bstrains = [bstrain for (bstrain,) in execute(dbOutput, "SELECT DISTINCT bstrain_id FROM strain_matches WHERE t = $(time)
        AND match_type = $(match_type)")]
        bstrainAbundances = [abundance for (abundance,) in execute(dbTemp, "SELECT abundance FROM babundance WHERE t = $(time)
        AND bstrain_id in ($(join(bstrains,", ")))")]
        bstrainTotal = sum(bstrainAbundances)
        vstrains = [vstrain for (vstrain,) in execute(dbOutput, "SELECT DISTINCT vstrain_id FROM strain_matches WHERE t = $(time)
        AND match_type = $(match_type)")]
        vstrainAbundances = [abundance for (abundance,) in execute(dbTemp, "SELECT abundance FROM vabundance WHERE t = $(time)
        AND vstrain_id in ($(join(vstrains,", ")))")]
        vstrainTotal = sum(vstrainAbundances)


        (first,) = execute(dbOutput, "SELECT DISTINCT match_id FROM match_types WHERE t = $(time)
            AND match_type = $(match_type) ORDER BY match_id")
        match_ref = [spacer for (spacer,) in
            execute(dbOutput, "SELECT spacer_id FROM match_types WHERE t = $(time) AND match_id = $(first.match_id)
            ORDER BY spacer_id")]
            matchDiv = 1
            matchDivTotal += 1
            #if time == 19.0
                #println("reference matches: $(match_ref)")
            #end
        for (match_id,) in execute(dbOutput, "SELECT DISTINCT match_id FROM match_types WHERE t = $(time)
            AND match_type = $(match_type) ORDER BY match_id")
            matches = [spacer for (spacer,) in
                execute(dbOutput, "SELECT spacer_id FROM match_types WHERE t = $(time) AND match_id = $(match_id)
                ORDER BY spacer_id")]
                #if time == 19.0
                    #println("matches: $(matches)")
                    #println("intersection: $(intersect(matches,match_ref))")
                    #println("length: $(length(intersect(matches,match_ref)))")
                    #println("match_type: $(match_type)")
                    #println("truth: $(length(intersect(matches,match_ref)) != match_type)")
                #end

            if length(intersect(matches,match_ref)) !== match_type
                matchDiv += 1
                matchDivTotal += 1
                match_ref = [spacer for (spacer,) in
                    execute(dbOutput, "SELECT spacer_id FROM match_types WHERE t = $(time) AND match_id = $(match_id)
                    ORDER BY spacer_id")]
            end
        end

        allSpacers = [spacer for spacer in execute(dbOutput, "SELECT spacer_id FROM match_types
        WHERE t = $(time) AND match_type = $(match_type)")]

        match_ids = [match_id for match_id in execute(dbOutput, "SELECT DISTINCT match_id FROM match_types WHERE t = $(time)
            AND match_type = $(match_type)")]
        execute(dbOutput, "INSERT INTO match_abundances VALUES (?,?,?,?,?,?,?,?,?,?,?,?)",
        (time, match_type, length(match_ids), length(match_ids)/totalMatches, matchDiv,
        length(unique(allSpacers)), bstrainTotal/btotal, vstrainTotal/vtotal,
        bstrainTotal, vstrainTotal, length(bstrains), length(vstrains)))
    end

    execute(dbOutput, "INSERT INTO total_abundances VALUES (?,?,?,?,?)",(time, btotal, vtotal, totalMatches, matchDivTotal))
end


function updateSingleSpacerFrequencies!(time,spacerFrequencies,spacer_id,babundance,vabundance,btotal,vtotal)
    if issubset(spacer_id,spacerFrequencies.spacer_id)
        @. spacerFrequencies.babundance_with_spacer_id =
        ifelse((spacerFrequencies.spacer_id == spacer_id) & (spacerFrequencies.t == time),
        spacerFrequencies.babundance_with_spacer_id + babundance,
        spacerFrequencies.babundance_with_spacer_id)

        @. spacerFrequencies.bfrequency_with_spacer_id =
        ifelse((spacerFrequencies.spacer_id == spacer_id) & (spacerFrequencies.t == time),
        (spacerFrequencies.babundance_with_spacer_id)/btotal,
        spacerFrequencies.bfrequency_with_spacer_id)

        @. spacerFrequencies.vabundance_with_spacer_id =
        ifelse((spacerFrequencies.spacer_id == spacer_id) & (spacerFrequencies.t == time),
        spacerFrequencies.vabundance_with_spacer_id + vabundance,
        spacerFrequencies.vabundance_with_spacer_id)

        @. spacerFrequencies.vfrequency_with_spacer_id =
        ifelse((spacerFrequencies.spacer_id == spacer_id) & (spacerFrequencies.t == time),
        (spacerFrequencies.vabundance_with_spacer_id)/vtotal,
        spacerFrequencies.vfrequency_with_spacer_id)
    else
        push!(spacerFrequencies,[time,spacer_id, babundance/btotal, vabundance/vtotal, babundance, vabundance])
    end
end


function compute_deltaPDI(time,bstrain_id1,vstrain_id,babundance1,vabundance,btotal,vtotal,spacerMatches1,pspacers,maxbAbundance)
    dDRI, dPDI, dPDI_degree1, dPDI_degree2, dPDI_degree3, dPDI_degree4, dPDI_degree5 = 0,0,0,0,0,0,0
    for (bstrain_id2,) in execute(dbTemp, "SELECT bstrain_id FROM babundance WHERE t = $(time)")
        ## YOU ARE HERE!
        if bstrain_id2 == bstrain_id1
            # println("Skipping iteration bstrain2 = $(bstrain_id2), because bstrain1 = $(bstrain_id1)")
            continue
        end
        # println("Skip check: iteration bstrain2 = $(bstrain_id2)")

        (babundance2,) = execute(dbTemp, "SELECT abundance FROM babundance WHERE t = $(time)
        AND bstrain_id = $(bstrain_id2)")
        babundance2 = babundance2.abundance
        spacers2 = [spacer_id for (spacer_id,) in
        execute(dbTemp, "SELECT spacer_id FROM bspacers WHERE bstrain_id = $(bstrain_id2)")]
        spacerMatches2 = intersect(spacers2,pspacers)
        match_type2 = length(spacerMatches2)

        if match_type2 == 0
            # println("Skipping iteration bstrain2 = $(bstrain_id2), because bstrain2 no match w vstrain = $(vstrain_id)")
            continue
        end

        # println("Skip check: iteration for bstrain2 = $(bstrain_id2) and vstrain = $(vstrain_id) match")

        in1not2 = length(setdiff(spacerMatches1,spacerMatches2))

        in2not1 = length(setdiff(spacerMatches2,spacerMatches1))

        disjointness = in1not2 + in2not1

        if in1not2 == 0 && in2not1 == 0
            dDRI += (1-abs(babundance1-babundance2)/maxbAbundance)*(babundance1/btotal)*(babundance2/btotal)*(vabundance/vtotal)
        end

        if disjointness > 0
            dPDI += (1-abs(babundance1-babundance2)/maxbAbundance)*(babundance1/btotal)*(babundance2/btotal)*(vabundance/vtotal)
        end

        if disjointness == 1
            dPDI_degree1 += (1-abs(babundance1-babundance2)/maxbAbundance)*(babundance1/btotal)*(babundance2/btotal)*(vabundance/vtotal)
        end

        if disjointness == 2
            dPDI_degree2 += (1-abs(babundance1-babundance2)/maxbAbundance)*(babundance1/btotal)*(babundance2/btotal)*(vabundance/vtotal)
        end

        if disjointness == 3
            dPDI_degree3 += (1-abs(babundance1-babundance2)/maxbAbundance)*(babundance1/btotal)*(babundance2/btotal)*(vabundance/vtotal)
        end

        if disjointness == 4
            dPDI_degree4 += (1-abs(babundance1-babundance2)/maxbAbundance)*(babundance1/btotal)*(babundance2/btotal)*(vabundance/vtotal)
        end

        if disjointness == 5
            dPDI_degree5 += (1-abs(babundance1-babundance2)/maxbAbundance)*(babundance1/btotal)*(babundance2/btotal)*(vabundance/vtotal)
        end

        # if in1not2 == 0 && in2not1 == 0
        #     dDRI += (1-abs(babundance1-babundance2)/maximum([babundance1,babundance2]))
        #     *(babundance1/btotal)*(babundance2/btotal)*(vabundance/vtotal)
        # end
        #
        # if disjointness > 0
        #     dPDI += (1-abs(babundance1-babundance2)/maximum([babundance1,babundance2]))
        #     *(babundance1/btotal)*(babundance2/btotal)*(vabundance/vtotal)
        # end
        #
        # if disjointness == 1
        #     dPDI_degree1 += (1-abs(babundance1-babundance2)/maximum([babundance1,babundance2]))
        #     *(babundance1/btotal)*(babundance2/btotal)*(vabundance/vtotal)
        # end
        #
        # if disjointness == 2
        #     dPDI_degree2 += (1-abs(babundance1-babundance2)/maximum([babundance1,babundance2]))
        #     *(babundance1/btotal)*(babundance2/btotal)*(vabundance/vtotal)
        # end
        #
        # if disjointness == 3
        #     dPDI_degree3 += (1-abs(babundance1-babundance2)/maximum([babundance1,babundance2]))
        #     *(babundance1/btotal)*(babundance2/btotal)*(vabundance/vtotal)
        # end
        #
        # if disjointness == 4
        #     dPDI_degree4 += (1-abs(babundance1-babundance2)/maximum([babundance1,babundance2]))
        #     *(babundance1/btotal)*(babundance2/btotal)*(vabundance/vtotal)
        # end
        #
        # if disjointness == 5
        #     dPDI_degree5 += (1-abs(babundance1-babundance2)/maximum([babundance1,babundance2]))
        #     *(babundance1/btotal)*(babundance2/btotal)*(vabundance/vtotal)
        # end
    end

    return dPDI, dDRI, dPDI_degree1, dPDI_degree2, dPDI_degree3, dPDI_degree4, dPDI_degree5
end

analyzeMatches()

println("Complete!")
