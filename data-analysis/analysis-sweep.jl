#!/usr/bin/env julia

"""
The purpose of this file is to illustrate how to do a parameter sweep in Julia.
This script loops through parameter combinations, and replicates with different
random seeds, and generates files necessary to perform runs on a local machine
or on a SLURM cluster.
To use this script for an experiment, you should copy this directory to a new
location, modify the parameter sweeps, and modify the relative paths to
`preamble.jl` and `ROOT_PATH` below to be correct, and the run
`./generate-runs.jl`
When the experiment is complete, you can collect all output files into a single
SQLite database via
`./gather-output.jl`
For each run, it creates a directory, `runs/c<combo_id>/r<replicate>`, and adds
entries to a SQLite database of run information, to make it easy to identify
runs and collate output.
It also divides runs into jobs suitable for execution on a single cluster node
or local machine. The runs are specified as lines in the job's `runs.txt`
file, and the job is specified in a `job.sbatch` file, which can be run directly
as a shell script or submitted to a SLURM cluster.
Each job uses the script `varmodel3/runmany.jl` to run a single-node, multi-core
queue of runs, with one run running on each core at any time.
This script also generates a script `submit_jobs.sh`, which submits every job to
SLURM at once.
Runs are divided into at most `N_JOBS_MAX` jobs that make use of at most
`N_CORES_PER_JOB_MAX` for the cluster node's local queue.
This allows you to work within limits set by your cluster administrator.
If you have no limits, you should set `N_JOBS_MAX` to a very large number,
and set `N_CORES_PER_JOB_MAX = 1`, so that the cluster can dynamically
balance runs across cluster nodes as the experiment runs.
To modify configuration settings for SLURM jobs, edit the template string in
the `generate_jobs()` function.
"""

# for other data change "richness" folder (lines ...), "make-richness-plots" (line 48)!!!!
# make sure data type has its own directory, analogous to "richness" in the
# data analysis folder

println("(Annoying Julia compilation delay...)")

include(joinpath("..","simulation","src","setup.jl"))

using SQLite

analysisType = ARGS[1]
analysisDir = "$(analysisType)"

SCRIPT_PATH = abspath(dirname(PROGRAM_FILE))
ROOT_RUN_SCRIPT = joinpath(SCRIPT_PATH,analysisDir,"$(analysisType).py")
ROOT_RUNMANY_SCRIPT = joinpath(SCRIPT_PATH,"src", "runmany.jl")
cd(SCRIPT_PATH)



# Number of replicates for each parameter combination
const N_REPLICATES = 30

# Number of SLURM jobs to generate
const N_JOBS_MAX = 100
const N_CORES_PER_JOB_MAX = 14 # Half a node (14) is easier to get scheduled than a whole one

function main()
    # Root run directory
    if !ispath(joinpath("..","simulation","runs"))
        error("`runs` does not exist; please simulate first.")
    end

    # Root job directory
    if !ispath(joinpath("..","simulation","jobs"))
        error("`jobs` does not exist; please simulate first.")
    end

    if isfile(joinpath("data_analysis.sqlite"))
        @info "\n\n\n\nNotice: `data_analysis.sqlite` already exists.\n
        Make sure other data is updated and is in accordance.\n\n\n\n"
    end
    # Connect or create data analysis databsase
    dbAnalysis = SQLite.DB(joinpath("data_analysis.sqlite"))

    # Create table for type of analysis. This will write over existing data
    execute(dbAnalysis, "DROP TABLE IF EXISTS $(analysisType)")
    execute(dbAnalysis, "CREATE TABLE $(analysisType) (run_id INTEGER, t REAL, microbe REAL, virus REAL)")

    # Connect to simulation data
    dbSim = SQLite.DB(joinpath("..","simulation","sweep_db.sqlite"))

    # Create little database that corresponds analysis runs to jobIDs for troubleshooting
    dbTempJobs = SQLite.DB(joinpath(analysisDir,"$(analysisType)jobs.sqlite"))
    execute(dbTempJobs, "DROP TABLE IF EXISTS jobs")
    execute(dbTempJobs, "CREATE TABLE jobs (job_id INTEGER, job_dir TEXT)")
    execute(dbTempJobs, "DROP TABLE IF EXISTS job_runs")
    execute(dbTempJobs, "CREATE TABLE job_runs (job_id INTEGER, run_id INTEGER, run_dir TEXT)")

    generate_analysis_runs(dbSim)
    generate_analysis_jobs(dbSim,dbTempJobs)

end

function generate_analysis_runs(dbSim) # This function generates the directories
    # for the individual parameter sets and corresponding replicates. It also
    # generates shell scripts for each run and corresponding parameter file.

    # Loop through parameter combinations and replicates, generating a run directory
    # `runs/c<combo_id>/r<replicate>` for each one.
    combo_id = 1
    run_id = 1
    for viral_mutation_rate in (1.0e-06, 2.0e-06)
        for spacer_acquisition_prob in (1e-06, 1e-05)
            for crispr_failure_prob in (1e-06, 1e-05)
                println("Processing analysis scripts for c$(combo_id): viral_mutation_rate = $(viral_mutation_rate),
                    spacer_acquisition_prob = $(spacer_acquisition_prob),
                    crispr_failure_prob = $(crispr_failure_prob)"
                )

                for replicate in 1:N_REPLICATES

                    run_dir = joinpath(analysisDir,"runs", "c$(combo_id)", "r$(replicate)")
                    if ispath(run_dir)
                        error("Please delete `$(analysisType)/runs`.")
                    end
                    mkpath(run_dir)


                    # Generate shell script to perform a single run
                    run_script = joinpath(run_dir, "run.sh")
                    open(run_script, "w") do f
                        print(f, """
                        #!/bin/sh
                        cd `dirname \$0`
                        julia $(ROOT_RUN_SCRIPT) $(run_id) &> output.txt
                        """)
                    end
                    run(`chmod +x $(run_script)`) # Make run script executable

                    run_id += 1
                end
                combo_id += 1
            end
        end
    end
end

function generate_analysis_jobs(dbSim,dbTempJobs)
    println("Assigning analysis runs to jobs...")

    # Assign runs to jobs (round-robin)
    job_id = 1
    for (run_id, run_dir) in execute(dbSim, "SELECT run_id, run_dir FROM runs ORDER BY replicate, combo_id")
        execute(dbTempJobs, "BEGIN TRANSACTION")
        execute(dbTempJobs, "INSERT INTO job_runs VALUES (?,?,?)", (job_id, run_id, run_dir))
        execute(dbTempJobs, "COMMIT")
        # Mod-increment job ID
        job_id = (job_id % N_JOBS_MAX) + 1
    end

    # Create job directories containing job scripts and script to submit all jobs
    submit_file = open("submit_analysis_jobs.sh", "w")
    println(submit_file, """
    #!/bin/sh
    cd `dirname \$0`
    """)

    for (job_id,) in execute(dbTempJobs, "SELECT DISTINCT job_id FROM job_runs ORDER BY job_id")

        job_dir = joinpath(analysisDir,"jobs", "$(job_id)")
        if !ispath(job_dir)
            mkpath(job_dir)
        else
            error("Please delete `$(analysisType)/jobs` and `$(analysisType)/runs`.")
        end

        # Get all run directories for this job
        run_dirs = [run_dir for (run_dir,) in execute(dbTempJobs,
            """
            SELECT run_dir FROM job_runs
            WHERE job_id = ?
            """,
            (job_id,)
        )]

        n_cores = min(length(run_dirs), N_CORES_PER_JOB_MAX)

        # Write out list of runs
        open(joinpath(job_dir, "runs.txt"), "w") do f
            for run_dir in run_dirs
                run_script = joinpath(SCRIPT_PATH, analysisDir, run_dir, "run.sh")
                println(f, run_script)
            end
        end

        # Create job sbatch file
        job_sbatch = joinpath(job_dir, "job.sbatch")
        open(job_sbatch, "w") do f
            print(f, """
            #!/bin/sh
            #SBATCH --account=pi-pascualmm
            #SBATCH --partition=broadwl
            #SBATCH --job-name=crispr-$(analysisType)-$(job_id)
            #SBATCH --tasks=1
            #SBATCH --cpus-per-task=$(n_cores)
            #SBATCH --mem-per-cpu=2000m
            #SBATCH --time=4:00:00
            #SBATCH --chdir=$(joinpath(SCRIPT_PATH, job_dir))
            #SBATCH --output=output.txt
            #SBATCH --mail-user=armun@uchicago.edu
            module purge
            # Uncomment this to use the Midway-provided Julia:
            module load julia
            julia $(ROOT_RUNMANY_SCRIPT) $(n_cores) runs.txt
            """) # runs.txt is for parallel processing
        end
        run(`chmod +x $(job_sbatch)`) # Make run script executable (for local testing)

        execute(dbTempJobs, "BEGIN TRANSACTION")
        execute(dbTempJobs, "INSERT INTO jobs VALUES (?,?)", (job_id, job_dir))
        execute(dbTempJobs, "COMMIT")

        println(submit_file, "sbatch $(job_sbatch)")
    end
    close(submit_file)
    run(`chmod +x submit_analysis_jobs.sh`) # Make submit script executable
end


main()
