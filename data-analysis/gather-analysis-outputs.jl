#!/usr/bin/env julia

"""
This script gathers the output generated by the sweep in `generate-sweep.jl`
into a single SQLite database.
This is probably more than you want, but it shows generally how to consolidate
things.
"""

println("(Annoying Julia compilation delay...)")

using SQLite
import SQLite.DBInterface.execute

analysisType = ARGS[1]
analysisDir = "$(analysisType)"

# Get relevant paths and cd to the script path.
SCRIPT_PATH = abspath(dirname(PROGRAM_FILE))
#cd(SCRIPT_PATH)

function main()
    # Source database containing experiment metadata
    #if !ispath("sweep_db.sqlite")
        #error("`sweep_db.sqlite` does not exist; please run the sweep first")
    #end

    # New database
    #if ispath("sweep_db_gathered.sqlite")
        #error("`sweep_db_gathered.sqlite` exists; please move or delete")
    #end

    # Copy source database as starting point
    #cp("sweep_db.sqlite", "sweep_db_gathered.sqlite")

    an_dir = joinpath("gathered-analyses")
    if !ispath(an_dir)
        mkpath(an_dir)
    end

    # Connect or create data analysis databsase
    dbAnalysis = SQLite.DB(joinpath(an_dir,"$(analysisType).sqlite"))

    if !isfile("$(analysisType)jobs.sqlite") # cluster
        error("$(analysisType)jobs.sqlite is missing; please run analysis first") # cluster
    end # cluster
    run_pairs = let
        db = SQLite.DB(joinpath(analysisDir,"$(analysisType)jobs.sqlite"))
        [(run_id, run_dir) for (run_id, run_dir) in execute(db, "SELECT run_id, run_dir FROM job_runs")]
    end

    is_first = true
    table_names = []
    for (run_id, run_dir) in run_pairs
        rundb_path = joinpath(analysisDir,run_dir, "$(analysisType)_output.sqlite")

        if !ispath(rundb_path)
            println("Missing: $(rundb_path)")
            continue
        else
            println("Processing: $(rundb_path)")
        end

        # Start up new connection; weird things seem to happen with attached databases otherwise
        db = SQLite.DB(joinpath(an_dir,"$(analysisType).sqlite"))
        execute(db, "BEGIN TRANSACTION")

        # Connect to output database as sub-db inside connection
        execute(db, "ATTACH DATABASE '$(joinpath(analysisDir,run_dir, "$(analysisType)_output.sqlite"))' as rundb")
        # If this is the first run, initialize tables
        if is_first
            for (table_name, sql) in execute(
                db, "SELECT tbl_name, sql FROM rundb.sqlite_master WHERE type = 'table' AND name NOT LIKE 'sqlite_%';"
            )
                #if table_name == "meta"
                    #continue
                #end

                # Use original SQL to create table with the same columns
                execute(db, sql)

                # Add a run_id column
                execute(db, "ALTER TABLE $(table_name) ADD COLUMN run_id INTEGER")

                push!(table_names, table_name)
            end

            is_first = false
        end

        # Copy all data from single-run DB into consolidated DB
        #execute(db, "INSERT INTO run_meta SELECT ?, * FROM rundb.meta", (run_id,))
        for table_name in table_names
            execute(db, "INSERT INTO $(table_name) SELECT *, ? FROM rundb.$(table_name)", (run_id,))
        end

        execute(db, "COMMIT")
    end
end

main()
