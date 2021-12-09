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
#dbOutputPath = joinpath("DI_output.sqlite") # cluster

dbSimPath = joinpath("/Volumes/Yadgah/sweep_db_gathered.sqlite") # local
dbSimInfoPath = joinpath("/Volumes/Yadgah/sweep_db.sqlite") # local
dbOutputPath = joinpath("/Volumes/Yadgah/DI_output.sqlite") # local

if isfile(dbOutputPath)
    error("DI_output.sqlite already exists; delete first")
end # cluster
##

dbSimInfo = SQLite.DB(dbSimInfoPath)
dbSim = SQLite.DB(dbSimPath)
dbOutput = SQLite.DB(dbOutputPath)


execute(dbOutput, "CREATE TABLE IDI (t REAL, IDI REAL)")
execute(dbOutput, "CREATE TABLE PDI (t REAL, PDI REAL, maxPDI REAL, DRI REAL, PDI_degree1 REAL,
PDI_degree2 REAL, PDI_degree3 REAL, PDI_degree4 REAL, PDI_degree5 REAL)")
execute(dbOutput, "CREATE TABLE HVI (t REAL, HVI REAL)")


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
function computeDI()
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
                    HVI += (babundance/btotal)*(vabundance/vtotal)
                    continue
                end

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


        # DO THIS IN A MORE EFFICIENT WAY
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
    end

    return dPDI, dDRI, dPDI_degree1, dPDI_degree2, dPDI_degree3, dPDI_degree4, dPDI_degree5
end

computeDI()

println("Complete!")
