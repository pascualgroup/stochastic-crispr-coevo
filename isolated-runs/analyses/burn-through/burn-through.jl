#!/usr/bin/env julia

println("(Julia compilation delay...)")

using SQLite
using DataFrames
using SQLite.DBInterface: execute
using Statistics

run_id = ARGS[1]

## Define Paths ##
SCRIPT_PATH = abspath(dirname(PROGRAM_FILE))
#
dbSimPath = joinpath(SCRIPT_PATH,"..","..","..","simulation","sweep_db_gathered.sqlite") # cluster
dbOutputPath = joinpath("burn-through_output.sqlite") # cluster

# dbSimPath = joinpath("/Volumes/Yadgah","crispr-sweep-7-2-2022/isolates/runID3297-c66-r47/runID3297-c66-r47.sqlite") # local
# dbMatchDivPath = joinpath("/Volumes/Yadgah","crispr-sweep-7-2-2022/isolates/runID3297-c66-r47/match-diversity_output.sqlite") # local
# dbMatchPath = joinpath("/Volumes/Yadgah","crispr-sweep-7-2-2022/isolates/runID3297-c66-r47/matches_output.sqlite") # local
# dbProbPath = joinpath("/Volumes/Yadgah/crispr-sweep-7-2-2022/isolates/runID3297-c66-r47/probability-of-emergence_output.sqlite") # local
# dbTriPath = joinpath("/Volumes/Yadgah/crispr-sweep-7-2-2022/isolates/runID3297-c66-r47/tripartite-networks_output.sqlite") # local
# dbOutputPath = joinpath("/Volumes/Yadgah/crispr-sweep-7-2-2022/isolates/runID3297-c66-r47/burn-through_output.sqlite") # local

if isfile(dbOutputPath)
    error("burn-through_output.sqlite already exists; delete first")
end
##
dbTempSim = SQLite.DB(dbSimPath)
(cr,) = execute(dbTempSim,"SELECT combo_id,replicate FROM runs WHERE run_id = $(run_id)")
dbMatchPath = joinpath(SCRIPT_PATH,"..","..","isolates",
    "runID$(run_id)-c$(cr.combo_id)-r$(cr.replicate)","matches_output.sqlite")
if !isfile(dbMatchPath)
    error("matches_output.sqlite does not exist; compute first")
end
dbMatchDivPath = joinpath(SCRIPT_PATH,"..","..","isolates",
    "runID$(run_id)-c$(cr.combo_id)-r$(cr.replicate)","match-diversity_output.sqlite")
if !isfile(dbMatchDivPath)
    error("match-diversity_output.sqlite does not exist; compute first")
end
dbProbPath = joinpath(SCRIPT_PATH,"..","..","isolates",
    "runID$(run_id)-c$(cr.combo_id)-r$(cr.replicate)","probability-of-emergence_output.sqlite")
if !isfile(dbProbPath)
    error("probability-of-emergence.sqlite does not exist; compute first")
end
dbTriPath = joinpath(SCRIPT_PATH,"..","..","isolates",
    "runID$(run_id)-c$(cr.combo_id)-r$(cr.replicate)","tripartite-networks_output.sqlite")
if !isfile(dbProbPath)
    error("tripartite-networks_output.sqlite does not exist; compute first")
end
dbTempMatch = SQLite.DB(dbMatchPath)
dbTempMatchDiv = SQLite.DB(dbMatchDivPath)
dbTempProb = SQLite.DB(dbProbPath)
dbTempTri = SQLite.DB(dbTriPath)
dbOutput = SQLite.DB(dbOutputPath)

execute(dbOutput, "CREATE TABLE burn_through_probability (t REAL, p REAL)")

dbTempSim = SQLite.DB()
execute(dbTempSim, "CREATE TABLE summary (t REAL, viral_abundance INTEGER, microbial_abundance INTEGER)")

execute(dbTempSim, "BEGIN TRANSACTION")
execute(dbTempSim,"ATTACH DATABASE '$(dbSimPath)' as dbSim")
execute(dbTempSim,"INSERT INTO summary (t, bstrain_id, abundance)
    SELECT t, bstrain_id, abundance FROM dbSim.babundance WHERE run_id = $(run_id);")
execute(dbTempSim,"INSERT INTO vabundance (t, vstrain_id, abundance)
    SELECT t, vstrain_id, abundance FROM dbSim.vabundance WHERE run_id = $(run_id);")
execute(dbTempSim, "COMMIT")


execute(dbTempSim, "BEGIN TRANSACTION")
execute(dbTempSim, "CREATE INDEX babundance_index ON babundance (t,bstrain_id,abundance)")
execute(dbTempSim, "CREATE INDEX vabundance_index ON vabundance (t,vstrain_id,abundance)")
execute(dbTempSim, "COMMIT")


# function pEmergeSpacers()
    vtriAbundance = Vector{Int64}()
    btriAbundance = Vector{Int64}()
    sBtriFrequency = Vector{Float64}()
    sbAbundance = Vector{Int64}()
    pburn = Vector{Float64}()
    times = Vector{Float64}()
    for (t,) in execute(dbTempTri,"SELECT DISTINCT t FROM single_match_tripartite_networks")
        println("Computing burn-through probability at t = $(t)")
        push!(times,t)
        (V,) = execute(dbTempSim,"SELECT viral_abundance
                        FROM summary WHERE t = $(t)")
        V = V.viral_abundance
        (B,) = execute(dbTempSim,"SELECT microbial_abundance
                        FROM summary WHERE t = $(t)")
        B = B.microbial_abundance
        vmatchIDs = [vmatchID for (vmatchID,) in
            execute(dbTempTri, "SELECT DISTINCT vmatch_id
            FROM single_match_tripartite_networks
            WHERE t = $(t) ORDER BY vmatch_id")]
        bmatchIDs = [match_id for (match_id,) in
            execute(dbTempTri, "SELECT DISTINCT bmatch_id
            FROM single_match_tripartite_networks
            WHERE t = $(t)")]

        Vtri = sum([abund for (abund,) in
            execute(dbTempTri, "SELECT vabundance
                FROM vmatches_abundances
                WHERE t = $(t)
                AND match_id in ($(join(vmatchIDs,", ")))")])
        push!(vtriAbundance,Vtri)
        sBtri = mean(1/Float64(B)*[Float64(abund) for (abund,) in
            execute(dbTempTri, "SELECT bsusceptible
                FROM vmatches_abundances
                WHERE t = $(t)
                AND match_id in ($(join(vmatchIDs,", ")))")])
        push!(sBtriFrequency,sBtri)
        Btri = sum([abund for (abund,) in
            execute(dbTempTri, "SELECT babundance
                FROM bmatches_abundances
                WHERE t = $(t)
                AND match_id in ($(join(bmatchIDs,", ")))")])
        push!(btriAbundance,Btri)


        v0strains = [strainID for (strainID,)
                        in execute(dbTempMatch,"SELECT DISTINCT vstrain_id
                        FROM bstrain_to_vstrain_matches
                        WHERE t = $(t) AND match_length = 1")]
        b0strains = [strainID for (strainID,) in
                        execute(dbTempMatch, "SELECT DISTINCT bstrain_id
                            FROM bstrain_to_vstrain_0matches
                            WHERE t = $(t)
                            AND vstrain_id in ($(join(v0strains,", ")))")]
        sbAbund = sum([strainID for (strainID,) in
                        execute(dbTempSim, "SELECT abundance
                            FROM babundance
                            WHERE t = $(t)
                            AND bstrain_id in ($(join(b0strains,", ")))")])
        push!(sbAbundance,sbAbund)

        pvburn = 0
        for vmatchID in vmatchIDs
            numSpacers = length([spacer_id for (spacer_id,) in
                execute(dbTempTri, "SELECT DISTINCT spacer_id
                FROM single_match_tripartite_networks
                WHERE t = $(t) AND vmatch_id = $(vmatchID)
                ORDER BY spacer_id")])

            bmatchIDs = [match_id for (match_id,) in
                execute(dbTempTri, "SELECT DISTINCT bmatch_id
                FROM single_match_tripartite_networks
                WHERE t = $(t) AND vmatch_id = $(vmatchID)")]

            (subV,) = execute(dbTempTri,"SELECT vabundance
                        FROM vmatches_abundances
                        WHERE t = $(t) AND match_id = $(vmatchID)")
            subV = subV.vabundance

            subB = sum([abund for (abund,) in
                execute(dbTempTri, "SELECT babundance
                    FROM bmatches_abundances
                    WHERE t = $(t)
                    AND match_id in ($(join(bmatchIDs,", ")))")])

            (pext,) = execute(dbTempProb,"SELECT p_extinction
                                    FROM vmatch_extinction
                                    WHERE t = $(t) AND vmatch_id = $(vmatchID)")
            pext = pext.p_extinction
            pemerge = 1 - pext
            pvburn += subV/V*subB/B*1/numSpacers*pemerge
        end
        push!(pburn,pvburn)

        # DataFrame(  t = times,
        #             vabundance = vtriAbundance,
        #             babundance = btriAbundance,
        #             ) |> SQLite.load!(dbOutput, "tripartite_abundances",ifnotexists=true)
        # DataFrame(  t = times,
        #             p_burn = pburn,
        #             ) |> SQLite.load!(dbOutput, "burn_through",ifnotexists=true)
    end
# end


pEmergeSpacers()

println("Complete!")

pburn

p = [p for (p,) in execute(dbTempProb,"SELECT p_emerge_expected FROM vmatch_expectation")]
times = [t for (t,) in execute(dbTempProb,"SELECT t FROM vmatch_expectation")]

plot(times,sbAbundance)


plot(times,vtriAbundance)




savefig("/Volumes/Yadgah/sbabund.pdf")
