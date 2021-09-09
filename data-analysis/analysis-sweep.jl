#!/usr/bin/env julia

println("(Annoying Julia compilation delay...)")

include(joinpath("..","simulation","src","setup.jl"))

using SQLite

analysisType = ARGS[1]
analysisDir = "$(analysisType)"

if analysisType == "peaks" && length(ARGS) < 2
    error("`peaks` analysis needs two arguments: upperThreshold lowerThreshold")
elseif analysisType == "walls" && length(ARGS) < 2
    error("`walls` analysis needs two arguments: upperThreshold lowerThreshold")
end

if analysisType == "walls" && !isfile(joinpath(SCRIPT_PATH,"richness","richness.sqlite"))
    error("`/richness/richness.sqlite` is missing; please analyze richness first.")
end

SCRIPT_PATH = abspath(dirname(PROGRAM_FILE))
ROOT_RUN_SCRIPT = joinpath(SCRIPT_PATH,analysisDir,"$(analysisType).jl")
ROOT_RUNMANY_SCRIPT = joinpath(SCRIPT_PATH,"src", "runmany.jl")
cd(SCRIPT_PATH)

# Number of SLURM jobs to generate
const N_JOBS_MAX = 100
const N_CORES_PER_JOB_MAX = 20 # Half a node (14) is easier to get scheduled than a whole one
const mem_per_cpu = 2000 # in MB 100MB = 1 GB

function main()
    # Root run directory
    if ispath(joinpath(SCRIPT_PATH,analysisDir,"runs"))
        error("Please move or delete `/$(analysisType)/runs`.")
    end

    if ispath(joinpath(SCRIPT_PATH,analysisDir,"jobs"))
        error("Please move or delete `/$(analysisType)/jobs`.")
    end

    if !ispath(joinpath(SCRIPT_PATH,"..","simulation","runs"))
        error("`/../simulation/runs` is missing; please simulate time series first.")
    end

    # Root job directory
    if !ispath(joinpath(SCRIPT_PATH,"..","simulation","jobs"))
        @info "Note that `/../simulation/jobs` is missing."
    end

    # Connect to simulation data
    dbSim = SQLite.DB(joinpath(SCRIPT_PATH,"..","simulation","sweep_db.sqlite"))

    # Create little database that corresponds analysis runs to jobIDs for troubleshooting
    dbTempJobs = SQLite.DB(joinpath(analysisDir,"$(analysisType)jobs.sqlite"))
    execute(dbTempJobs, "DROP TABLE IF EXISTS jobs")
    execute(dbTempJobs, "CREATE TABLE jobs (job_id INTEGER, job_dir TEXT)")
    execute(dbTempJobs, "DROP TABLE IF EXISTS job_runs")
    execute(dbTempJobs, "CREATE TABLE job_runs (job_id INTEGER, run_id INTEGER, run_dir TEXT)")

    numSubmits = generate_analysis_runs(dbSim)
    generate_analysis_jobs(dbSim,dbTempJobs,numSubmits)
end

function generate_analysis_runs(dbSim::DB) # This function generates the directories
    # for the individual parameter sets and corresponding replicates. It also
    # generates shell scripts for each run and corresponding parameter file.

    # Loop through parameter combinations and replicates, generating a run directory
    # `runs/c<combo_id>/r<replicate>` for each one.
    run_count = 0
    println("Processing analysis script for each run")
    for (run_id, combo_id, replicate) in execute(dbSim, "SELECT run_id,combo_id,replicate FROM runs")
        #println("Processing analysis script for combination $(combo_id)/replicate $(replicate)"
        #) # local
        run_dir = joinpath(analysisDir,"runs", "c$(combo_id)", "r$(replicate)")
        @assert !ispath(run_dir)
        mkpath(run_dir)

        argString = map(x->string("$(x) "), ARGS) # the space after $(x) is important
        popfirst!(argString)

        # Generate shell script to perform a single run
        run_script = joinpath(run_dir, "run.sh")
        open(run_script, "w") do f
            print(f, """
            #!/bin/sh
            cd `dirname \$0`
            julia $(ROOT_RUN_SCRIPT) $(run_id) $(argString...) &> output.txt
            """)
        end
        run(`chmod +x $(run_script)`) # Make run script executable
        run_count += 1
    end
    return numSubmits = Int64(ceil(run_count/(N_JOBS_MAX*N_CORES_PER_JOB_MAX)))
end

function generate_analysis_jobs(dbSim::DB,dbTempJobs::DB,numSubmits::Int64)
    println("Assigning analysis runs to jobs...")

    # Assign runs to jobs (round-robin)
    job_id = 1
    job_count = 0
    n_cores_count = 0

    execute(dbTempJobs, "BEGIN TRANSACTION")

    for (run_id, run_dir) in execute(dbSim, "SELECT run_id, run_dir FROM runs ORDER BY replicate, combo_id")
        execute(dbTempJobs, "INSERT INTO job_runs VALUES (?,?,?)", (job_id, run_id, run_dir))
        # Mod-increment job ID
        job_id = (job_id % N_JOBS_MAX*numSubmits) + 1
        if job_id > job_count
            job_count = job_id
        end
    end

    submitScripts = IOStream[]
    for script in 1:numSubmits
        push!(submitScripts,open("$(script)-analysis_submit_jobs.sh", "w"))
        println(submitScripts[script], """
        #!/bin/sh
        cd `dirname \$0`
        """)
    end

    for (job_id,) in execute(dbTempJobs, "SELECT DISTINCT job_id FROM job_runs ORDER BY job_id")

        job_dir = joinpath(SCRIPT_PATH,analysisDir,"jobs", "$(job_id)")
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
            #SBATCH --mem-per-cpu=$(mem_per_cpu)m
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

        execute(dbTempJobs, "INSERT INTO jobs VALUES (?,?)", (job_id, job_dir))

        submitScript = (job_id % numSubmits) + 1
        println(submitScripts[submitScript], "sbatch $(job_sbatch)")
    end
    execute(dbTempJobs, "COMMIT")
    map(close,submitScripts)

    #run(`chmod +x submit_analysis_jobs.sh`) # Make submit script executable
    @info "
    Sweep will be submitted via $(numSubmits) `analysis_submit_jobs.sh` script(s).
    Each `analysis_submit_jobs.sh` script submits $(job_count) jobs.
    Each job will use $(n_cores_count) cpus (cores) at most, where each cpu will use $(mem_per_cpu/1000)GB.
    Each job therefore will use at most $(n_cores_count*mem_per_cpu/1000)GB of memory in total.
    "
end


main()
