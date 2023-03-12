#!/usr/bin/env julia

println("(Annoying Julia compilation delay...)")

include(joinpath("..","simulation","src","setup.jl"))

using SQLite

analysisType = ARGS[1]

SCRIPT_PATH = abspath(dirname(PROGRAM_FILE))
ROOT_RUN_SCRIPT = joinpath(SCRIPT_PATH,analysisType,"mean-$(analysisType).jl")
ROOT_RUNMANY_SCRIPT = joinpath(SCRIPT_PATH,"src", "runmany.jl")
cd(SCRIPT_PATH)

# Number of SLURM jobs to generate
const N_JOBS_MAX = 100
const N_CORES_PER_JOB_MAX = 28 # Half a node (14) is easier to get scheduled than a whole one
const mem_per_cpu = 2000 # in MB 100MB = 1 GB


function main()
    # Root run directory
    if ispath(joinpath(SCRIPT_PATH,analysisType,"mean","runs"))
        error("Please move or delete `/$(analysisType)/mean/runs`.")
    end

    if ispath(joinpath(SCRIPT_PATH,analysisType,"mean","jobs"))
        error("Please move or delete `/$(analysisType)/mean/jobs`.")
    end

    if !isfile(joinpath(SCRIPT_PATH,"gathered-analyses",analysisType,"$(analysisType).sqlite"))
        error("`/gathered-analyses/$(analysisType).sqlite` is missing; please analyze data for individual runs first.")
    end

    # Connect to simulation data
    dbSim = SQLite.DB(joinpath(SCRIPT_PATH,"..","simulation","sweep_db.sqlite"))

    # Create little database that corresponds analysis runs to jobIDs for troubleshooting
    dbTempJobs = SQLite.DB(joinpath(analysisType,"mean$(analysisType)jobs.sqlite"))
    execute(dbTempJobs, "DROP TABLE IF EXISTS jobs")
    execute(dbTempJobs, "CREATE TABLE jobs (job_id INTEGER, job_dir TEXT)")
    execute(dbTempJobs, "DROP TABLE IF EXISTS job_runs")
    execute(dbTempJobs, "CREATE TABLE job_runs (job_id INTEGER, combo_id INTEGER, run_dir TEXT)")

    numSubmits = generate_mean_runs(dbSim)
    generate_mean_jobs(dbSim,dbTempJobs,numSubmits)
end

function generate_mean_runs(dbSim::DB) # This function generates the directories
    # for the individual parameter sets and corresponding replicates. It also
    # generates shell scripts for each run and corresponding parameter file.

    # Loop through parameter combinations and replicates, generating a run directory
    # `runs/c<combo_id>/r<replicate>` for each one.
    run_count = 0
    println("Processing analysis script for each run")
    for (combo_id,) in execute(dbSim, "SELECT combo_id FROM param_combos ORDER BY combo_id")
        #println("Processing analysis script for combination $(combo_id)/replicate $(replicate)"
        #) # local
        run_dir = joinpath(analysisType,"mean","runs", "c$(combo_id)")
        @assert !ispath(run_dir)
        mkpath(run_dir)

        # Generate shell script to perform a single run
        run_script = joinpath(run_dir, "run.sh")
        open(run_script, "w") do f
            print(f, """
            #!/bin/sh
            cd `dirname \$0`
            julia $(ROOT_RUN_SCRIPT) $(combo_id) &> output.txt
            """)
        end
        run(`chmod +x $(run_script)`) # Make run script executable
        run_count += 1
    end
    return numSubmits = Int64(ceil(run_count/(N_JOBS_MAX*N_CORES_PER_JOB_MAX)))
end

function generate_mean_jobs(dbSim::DB,dbTempJobs::DB,numSubmits::Int64)
    println("Assigning analysis runs to jobs...")

    # Assign runs to jobs (round-robin)
    job_id = 1
    n_cores_count = 0
    execute(dbTempJobs, "BEGIN TRANSACTION")

    for (combo_id,) in execute(dbSim, "SELECT combo_id FROM param_combos ORDER BY combo_id")
        run_dir = joinpath("mean","runs", "c$(combo_id)")
        execute(dbTempJobs, "INSERT INTO job_runs VALUES (?,?,?)", (job_id, combo_id, run_dir))
        # Mod-increment job ID
        job_id = mod(combo_id,N_JOBS_MAX*numSubmits) + 1
    end

    submitScripts = IOStream[]
    for script in 1:numSubmits
        push!(submitScripts,open("$(script)_mean-analysis-submit-jobs.sh", "w"))
        println(submitScripts[script], """
        #!/bin/sh
        cd `dirname \$0`
        """)
    end

    for (job_id,) in execute(dbTempJobs, "SELECT DISTINCT job_id FROM job_runs ORDER BY job_id")

        job_dir = joinpath(SCRIPT_PATH,analysisType,"mean","jobs","$(job_id)")
        @assert !ispath(job_dir)
        mkpath(job_dir)

        # Get all run directories for this job
        run_dirs = [run_dir for (run_dir,) in execute(dbTempJobs,
            """
            SELECT run_dir FROM job_runs
            WHERE job_id = ?
            """,
            (job_id,)
        )]

        n_cores = min(length(run_dirs), N_CORES_PER_JOB_MAX)

        if n_cores > n_cores_count
            n_cores_count = n_cores
        end

        # Write out list of runs
        open(joinpath(job_dir, "runs.txt"), "w") do f
            for run_dir in run_dirs
                run_script = joinpath(SCRIPT_PATH, analysisType, run_dir, "run.sh")
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
            #SBATCH --job-name=crispr-mean-$(analysisType)-$(job_id)
            #SBATCH --tasks=1
            #SBATCH --cpus-per-task=$(n_cores)
            #SBATCH --mem-per-cpu=$(mem_per_cpu)m
            #SBATCH --time=1-12:00:00
            #SBATCH --chdir=$(joinpath(SCRIPT_PATH, job_dir))
            #SBATCH --output=output.txt
            #SBATCH --mail-user=armun@uchicago.edu
            # Uncomment this to use the Midway-provided Julia:
            module load julia
            julia $(ROOT_RUNMANY_SCRIPT) $(n_cores) runs.txt
            """) # runs.txt is for parallel processing
        end
        run(`chmod +x $(job_sbatch)`) # Make run script executable (for local testing)

        execute(dbTempJobs, "INSERT INTO jobs VALUES (?,?)", (job_id, job_dir))

        submitScript = mod(job_id,numSubmits) + 1
        println(submitScripts[submitScript], "sbatch $(job_sbatch)")
    end
    execute(dbTempJobs, "COMMIT")
    map(close,submitScripts)

    #run(`chmod +x submit_analysis_jobs.sh`) # Make submit script executable
    @info "
    Sweep will be submitted via $(numSubmits) `mean-analysis-submit-jobs.sh` script(s).
    Each `mean-analysis_submit_jobs.sh` script submits $(N_JOBS_MAX) jobs.
    Each job will use $(n_cores_count) cpus (cores) at most, where each cpu will use $(mem_per_cpu/1000)GB.
    Each job therefore will use at most $(n_cores_count*mem_per_cpu/1000)GB of memory in total.
    "
end

main()
