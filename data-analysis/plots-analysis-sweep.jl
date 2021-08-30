#!/usr/bin/env julia

println("(Annoying Julia compilation delay...)")

include(joinpath("..","simulation","src","setup.jl"))

using SQLite

analysisType = ARGS[1]
analysisDir = "$(analysisType)"

SCRIPT_PATH = abspath(dirname(PROGRAM_FILE))
ROOT_RUN_SCRIPT = joinpath(SCRIPT_PATH,analysisDir,"make-$(analysisType)-plots.py")
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
        error("`../simulation/runs` does not exist; please simulate time series first.")
    end

    # Root job directory
    if !ispath(joinpath("..","simulation","jobs"))
        error("`../simulation/jobs` does not exist; please simulate time series first.")
    end

    # Connect to simulation data
    dbSim = SQLite.DB(joinpath("..","simulation","sweep_db.sqlite"))

    # Use the little database that corresponds analysis runs to jobIDs for troubleshooting
    if !isfile(joinpath(analysisDir,"$(analysisType)jobs.sqlite"))
        error("`$(analysisType)jobs.sqlite` in /$(analysisType) is missing. Please analyze simulation data first.")
    end
    dbTempJobs = SQLite.DB(joinpath(analysisDir,"$(analysisType)jobs.sqlite"))
    execute(dbTempJobs, "DROP TABLE IF EXISTS plot_jobs")
    execute(dbTempJobs, "CREATE TABLE plot_jobs (job_id INTEGER, job_dir TEXT)")
    execute(dbTempJobs, "DROP TABLE IF EXISTS plot_job_runs")
    execute(dbTempJobs, "CREATE TABLE plot_job_runs (job_id INTEGER, run_id INTEGER, run_dir TEXT)")

    generate_plot_runs(dbSim)
    generate_plot_jobs(dbSim,dbTempJobs)

end

function generate_plot_runs(dbSim) # This function generates the directories
    # for the individual parameter sets and corresponding replicates. It also
    # generates shell scripts for each run and corresponding parameter file.

    # Loop through parameter combinations and replicates, generating a run directory
    # `runs/c<combo_id>/r<replicate>` for each one.
    combo_id = 1
    run_id = 1
    for viral_mutation_rate in (1.0e-06, 2.0e-06)
        for spacer_acquisition_prob in (1e-06, 1e-05)
            for crispr_failure_prob in (1e-06, 1e-05)
                println("Processing plot scripts for c$(combo_id): viral_mutation_rate = $(viral_mutation_rate),
                    spacer_acquisition_prob = $(spacer_acquisition_prob),
                    crispr_failure_prob = $(crispr_failure_prob)"
                )

                for replicate in 1:N_REPLICATES

                    run_dir = joinpath(analysisDir,"runs", "c$(combo_id)", "r$(replicate)")
                    if !ispath(run_dir)
                        error("Please analyze data first; `$(analysisType)/runs` is missing.")
                    end

                    plot_dir = joinpath("gathered-analyses",analysisDir,"plots", "c$(combo_id)", "r$(replicate)")
                    if ispath(plot_dir)
                        error("Please delete `gathered-analyses/$(analysisType)/plots`.")
                    end
                    mkpath(plot_dir)

                    # Generate shell script to perform a single run
                    run_script = joinpath(run_dir, "runplotmaker.sh")
                    open(run_script, "w") do f
                        print(f, """
                        #!/bin/sh
                        cd `dirname \$0`
                        python $(ROOT_RUN_SCRIPT) $(run_id) &> plot_output.txt
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

function generate_plot_jobs(dbSim,dbTempJobs)
    println("Assigning plot runs to jobs...")

    # Assign runs to jobs (round-robin)
    job_id = 1
    for (run_id, run_dir) in execute(dbSim, "SELECT run_id, run_dir FROM runs ORDER BY replicate, combo_id")
        execute(dbTempJobs, "BEGIN TRANSACTION")
        execute(dbTempJobs, "INSERT INTO plot_job_runs VALUES (?,?,?)", (job_id, run_id, run_dir))
        execute(dbTempJobs, "COMMIT")
        # Mod-increment job ID
        job_id = (job_id % N_JOBS_MAX) + 1
    end

    # Create job directories containing job scripts and script to submit all jobs
    submit_file = open("submit_plot_jobs.sh", "w")
    println(submit_file, """
    #!/bin/sh
    cd `dirname \$0`
    """)

    for (job_id,) in execute(dbTempJobs, "SELECT DISTINCT job_id FROM plot_job_runs ORDER BY job_id")

        job_dir = joinpath(analysisDir,"plotjobs", "$(job_id)")
        if !ispath(job_dir)
            mkpath(job_dir)
        else
            error("Please delete `$(analysisType)/plotjobs` and `gathered-analyses/$(analysisType)/plots`.")
        end

        # Get all run directories for this job
        run_dirs = [run_dir for (run_dir,) in execute(dbTempJobs,
            """
            SELECT run_dir FROM plot_job_runs
            WHERE job_id = ?
            """,
            (job_id,)
        )]


        n_cores = min(length(run_dirs), N_CORES_PER_JOB_MAX) #!


        # Write out list of runs
        open(joinpath(job_dir, "plot_runs.txt"), "w") do f
            for run_dir in run_dirs
                run_script = joinpath(SCRIPT_PATH, analysisDir, run_dir, "runplotmaker.sh")
                println(f, run_script)
            end
        end

        # Create job sbatch file
        job_sbatch = joinpath(job_dir, "plotjob.sbatch")
        open(job_sbatch, "w") do f
            print(f, """
            #!/bin/sh
            #SBATCH --account=pi-pascualmm
            #SBATCH --partition=broadwl
            #SBATCH --job-name=crispr-$(analysisType)-plots-$(job_id)
            #SBATCH --tasks=1
            #SBATCH --cpus-per-task=$(n_cores)
            #SBATCH --mem-per-cpu=6500m
            #SBATCH --time=4:00:00
            #SBATCH --chdir=$(joinpath(SCRIPT_PATH, job_dir))
            #SBATCH --output=plot_output.txt
            #SBATCH --mail-user=armun@uchicago.edu
            module purge
            # Uncomment this to use the Midway-provided Julia:
            module load julia
            julia $(ROOT_RUNMANY_SCRIPT) $(n_cores) plot_runs.txt
            """) # runs.txt is for parallel processing
        end
        run(`chmod +x $(job_sbatch)`) # Make run script executable (for local testing)

        execute(dbTempJobs, "BEGIN TRANSACTION")
        execute(dbTempJobs, "INSERT INTO plot_jobs VALUES (?,?)", (job_id, job_dir))
        execute(dbTempJobs, "COMMIT")

        println(submit_file, "sbatch $(job_sbatch)")
    end
    close(submit_file)
    run(`chmod +x submit_plot_jobs.sh`) # Make submit script executable
end


main()
