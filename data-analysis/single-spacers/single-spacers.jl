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
#dbMatchesPath = joinpath(SCRIPT_PATH,"..","gathered-analyses","match-diversity","match-diversity.sqlite") # local
#dbOutputPath = joinpath("match-diversity_output.sqlite") # cluster

dbSimPath = joinpath("/Volumes/Yadgah/sweep_db_gathered.sqlite") # local
dbSimInfoPath = joinpath("/Volumes/Yadgah/sweep_db.sqlite") # local
dbMatchesPath = joinpath("/Volumes/Yadgah/match-diversity.sqlite") # local
dbOutputPath = joinpath("/Volumes/Yadgah/single-spacers_output.sqlite") # local

if isfile(dbOutputPath)
    error("single_spacers_output.sqlite already exists; delete first")
end # cluster
##

dbSimInfo = SQLite.DB(dbSimInfoPath)
dbSim = SQLite.DB(dbSimPath)
dbMatches = SQLite.DB(dbMatchesPath)
dbOutput = SQLite.DB(dbOutputPath)


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
function singleSpacers()
    for (time,) in execute(dbTemp, "SELECT DISTINCT t FROM summary ORDER BY t")
        # println("time is $(time)")
        if time > 30.0
            return
        end
        (btotal,) = execute(dbTemp, "SELECT microbial_abundance FROM summary WHERE t = $(time)")
        btotal = btotal.microbial_abundance
        (vtotal,) = execute(dbTemp, "SELECT viral_abundance FROM summary WHERE t = $(time)")
        vtotal = vtotal.viral_abundance

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
            for (bstrain_id,) in execute(dbMatches, "SELECT DISTINCT bstrain_id FROM strain_matches WHERE t = $(time)")
                spacers = [spacer for (spacer,) in
                execute(dbTemp, "SELECT spacer_id FROM bspacers WHERE bstrain_id = $(bstrain_id)")]
                if issubset(spacer_id,spacers)
                    (bStrainAbundance,) = execute(dbTemp, "SELECT abundance FROM babundance WHERE t = $(time)
                    AND bstrain_id = $(bstrain_id)")
                    babundance += bStrainAbundance.abundance
                end
            end
            for (vstrain_id,) in execute(dbMatches, "SELECT DISTINCT vstrain_id FROM strain_matches WHERE t = $(time)")
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

    end

end

singleSpacers()

println("Complete!")
