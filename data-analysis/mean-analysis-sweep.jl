#!/usr/bin/env julia

println("(Annoying Julia compilation delay...)")

include(joinpath("..","simulation","src","setup.jl"))

using SQLite

analysisType = ARGS[1]

SCRIPT_PATH = abspath(dirname(PROGRAM_FILE))

ROOT_COMBO_SCRIPT = joinpath(SCRIPT_PATH,"$(analysisType)","mean-$(analysisType).jl")
ROOT_RUNMANY_SCRIPT = joinpath(SCRIPT_PATH,"src", "runmany.jl")
cd(SCRIPT_PATH)

# Number of SLURM jobs to generate
const N_JOBS_MAX = 100
const N_CORES_PER_JOB_MAX = 20 # Half a node (14) is easier to get scheduled than a whole one
const mem_per_cpu = 2000 # in MB 100MB = 1 GB

function main()
    # Root run directory

    if !ispath(joinpath(SCRIPT_PATH,"gathered-analyses","$(analysisType)","$(analysisType).sqlite"))
        error("`$(analysisType).sqlite` is missing; please simulate and gather $(analysisType) analysis first.")
    end

    if ispath(joinpath(SCRIPT_PATH,analysisType,"mean-runs"))
        error("Please move or delete `/$(analysisType)/mean-runs`.")
    end

    if ispath(joinpath(SCRIPT_PATH,analysisType,"mean-jobs"))
        error("Please move or delete `/$(analysisType)/mean-jobs`.")
    end

    # Connect to simulation data
    dbSim = SQLite.DB(joinpath(SCRIPT_PATH,"..","simulation","sweep_db.sqlite"))

    # Create little database that corresponds analysis runs to jobIDs for troubleshooting
    dbTempJobs = SQLite.DB(joinpath(analysisType,"mean$(analysisType)jobs.sqlite"))
    execute(dbTempJobs, "DROP TABLE IF EXISTS jobs")
    execute(dbTempJobs, "CREATE TABLE jobs (job_id INTEGER, job_dir TEXT)")
    execute(dbTempJobs, "DROP TABLE IF EXISTS job_combos")
    execute(dbTempJobs, "CREATE TABLE job_combos (job_id INTEGER, combo_id INTEGER, combo_dir TEXT)")

    numSubmits = generate_analysis_runs(dbSim)
    generate_analysis_jobs(dbSim,dbTempJobs,numSubmits)
end

function generate_analysis_runs(dbSim::DB) # This function generates the directories
    # for the individual parameter sets and corresponding replicates. It also
    # generates shell scripts for each run and corresponding parameter file.

    # Loop through parameter combinations and replicates, generating a run directory
    # `runs/c<combo_id>/r<replicate>` for each one.
    run_count = 0
    println("Processing analysis script for each combination")
    for (combo_id,) in execute(dbSim, "SELECT combo_id FROM param_combos")
        #println("Processing analysis script for combination $(combo_id)/replicate $(replicate)"
        #) # local
        combo_dir = joinpath(analysisType,"mean-runs", "c$(combo_id)")
        @assert !ispath(combo_dir)
        mkpath(combo_dir)

        # Generate shell script to perform a single run
        combo_script = joinpath(combo_dir, "run.sh")
        open(combo_script, "w") do f
            print(f, """
            #!/bin/sh
            cd `dirname \$0`
            julia $(ROOT_COMBO_SCRIPT) $(combo_id) &> output.txt
            """)
        end
        run(`chmod +x $(combo_script)`) # Make run script executable
        run_count += 1
    end
    return numSubmits = Int64(ceil(run_count/(N_JOBS_MAX*N_CORES_PER_JOB_MAX)))
end

function generate_analysis_jobs(dbSim::DB,dbTempJobs::DB,numSubmits::Int64)
    println("Assigning analysis runs to jobs...")

    # Assign runs to jobs (round-robin)
    job_id = 1
    n_cores_count = 0

    execute(dbTempJobs, "BEGIN TRANSACTION")

    for (combo_id,) in execute(dbSim, "SELECT combo_id FROM param_combos")
        combo_dir = joinpath("mean-runs", "c$(combo_id)")
        execute(dbTempJobs, "INSERT INTO job_combos VALUES (?,?,?)", (job_id, combo_id, combo_dir))
        # Mod-increment job ID
        job_id = mod(job_id,N_JOBS_MAX*numSubmits) + 1
    end

    submitScripts = IOStream[]
    for script in 1:numSubmits
        push!(submitScripts,open("$(script)_mean-analysis-submit-jobs.sh", "w"))
        println(submitScripts[script], """
        #!/bin/sh
        cd `dirname \$0`
        """)
    end

    for (job_id,) in execute(dbTempJobs, "SELECT DISTINCT job_id FROM job_combos ORDER BY job_id")

        job_dir = joinpath(SCRIPT_PATH,analysisType,"mean-jobs", "$(job_id)")
        @assert !ispath(job_dir)
        mkpath(job_dir)

        # Get all run directories for this job
        combo_dirs = [combo_dir for (combo_dir,) in execute(dbTempJobs,
            """
            SELECT combo_dir FROM job_combos
            WHERE job_id = ?
            """,
            (job_id,)
        )]

        n_cores = min(length(combo_dirs), N_CORES_PER_JOB_MAX)

        if n_cores > n_cores_count
            n_cores_count = n_cores
        end

        # Write out list of runs
        open(joinpath(job_dir, "runs.txt"), "w") do f
            for combo_dir in combo_dirs
                combo_script = joinpath(SCRIPT_PATH, analysisType, combo_dir, "run.sh")
                println(f, combo_script)
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
            #SBATCH --time=4:00:00
            #SBATCH --chdir=$(joinpath(SCRIPT_PATH, job_dir))
            #SBATCH --output=output.txt
            #SBATCH --mail-user=armun@uchicago.edu
            module purge
            # Uncomment this to use the Midway-provided Julia:
            module load julia
            julia $(ROOT_RUNMANY_SCRIPT) $(n_cores) runs.txt
            """)
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
    Sweep will be submitted via $(numSubmits) `analysis_submit_jobs.sh` script(s).
    Each `analysis_submit_jobs.sh` script submits $(N_JOBS_MAX) jobs.
    Each job will use $(n_cores_count) cpus (cores) at most, where each cpu will use $(mem_per_cpu/1000)GB.
    Each job therefore will use at most $(n_cores_count*mem_per_cpu/1000)GB of memory in total.
    "
end


main()