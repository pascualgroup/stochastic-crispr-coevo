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

# Number of SLURM jobs to generate
const N_JOBS_MAX = 100
const N_CORES_PER_JOB_MAX = 20 # Half a node (14) is easier to get scheduled than a whole one
const mem_per_cpu = 8500 # in MB 100MB = 1 GB

function main()
    # Root run directory
    if ispath(joinpath(SCRIPT_PATH,"gathered-analyses",analysisDir,"plots"))
        error("Please move or delete `/gathered-analyses/$(analysisType)/plots`.")
    end

    if ispath(joinpath(SCRIPT_PATH,analysisDir,"plotjobs"))
        error("Please move or delete `/$(analysisType)/plotjobs`.")
    end

    if !ispath(joinpath(SCRIPT_PATH,"..","simulation","runs"))
        @info "Note that `/../simulation/runs` is missing."
    end

    # Root job directory
    if !ispath(joinpath(SCRIPT_PATH,"..","simulation","jobs"))
        @info "Note that `/../simulation/jobs` is missing."
    end

    if !ispath(joinpath(SCRIPT_PATH,analysisDir,"runs"))
        error("`/$(analysisType)/runs` is missing; please analyze data first")
    end

    if !ispath(joinpath(SCRIPT_PATH,analysisDir,"jobs"))
        @info "Note that `/$(analysisType)/jobs` is missing."
    end

    # Use the little database that corresponds analysis runs to jobIDs for troubleshooting
    if !isfile(joinpath(SCRIPT_PATH,analysisDir,"$(analysisType)jobs.sqlite"))
        error("`/$(analysisType)/$(analysisType)jobs.sqlite` is missing; please analyze simulation data first.")
    end

    # Connect to simulation data
    dbSim = SQLite.DB(joinpath(SCRIPT_PATH,"..","simulation","sweep_db.sqlite"))

    dbTempJobs = SQLite.DB(joinpath(SCRIPT_PATH,analysisDir,"$(analysisType)jobs.sqlite"))
    execute(dbTempJobs, "DROP TABLE IF EXISTS plot_jobs")
    execute(dbTempJobs, "CREATE TABLE plot_jobs (job_id INTEGER, job_dir TEXT)")
    execute(dbTempJobs, "DROP TABLE IF EXISTS plot_job_runs")
    execute(dbTempJobs, "CREATE TABLE plot_job_runs (job_id INTEGER, run_id INTEGER, run_dir TEXT)")

    numSubmits = generate_plot_runs(dbSim::DB)
    generate_plot_jobs(dbSim,dbTempJobs,numSubmits)
end

function generate_plot_runs(dbSim::DB) # This function generates the directories
    # for the individual parameter sets and corresponding replicates. It also
    # generates shell scripts for each run and corresponding parameter file.

    # Loop through parameter combinations and replicates, generating a run directory
    # `runs/c<combo_id>/r<replicate>` for each one.
    run_count = 0
    println("Processing analysis plot script for each run")
    for (run_id, combo_id, replicate) in execute(dbSim, "SELECT run_id,combo_id,replicate FROM runs")
        #println("Processing analysis plot script for combination $(combo_id)/replicate $(replicate)"
        #) # local
        run_dir = joinpath(SCRIPT_PATH,analysisDir,"runs", "c$(combo_id)", "r$(replicate)")
        @assert ispath(run_dir)

        plot_dir = joinpath(SCRIPT_PATH,"gathered-analyses",analysisDir,"plots", "c$(combo_id)", "r$(replicate)")
        @assert !ispath(plot_dir)
        mkpath(plot_dir)

        # Generate shell script to perform a single run
        run_script = joinpath(run_dir, "runplotmaker.sh")
        open(run_script, "w") do f
            print(f, """
            #!/bin/sh
            cd `dirname \$0`
            module load python/anaconda-2021.05
            python $(ROOT_RUN_SCRIPT) $(run_id) &> plot_output.txt
            """)
        end
        run(`chmod +x $(run_script)`) # Make run script executable
        run_count += 1
    end
    return numSubmits = Int64(ceil(run_count/(N_JOBS_MAX*N_CORES_PER_JOB_MAX)))
end

function generate_plot_jobs(dbSim::DB,dbTempJobs::DB,numSubmits::Int64)
    println("Assigning analysis plot runs to jobs...")

    # Assign runs to jobs (round-robin)
    job_id = 1
    job_count = 0
    n_cores_count = 0

    execute(dbTempJobs, "BEGIN TRANSACTION")

    for (run_id, run_dir) in execute(dbSim, "SELECT run_id, run_dir FROM runs ORDER BY replicate, combo_id")
        execute(dbTempJobs, "INSERT INTO plot_job_runs VALUES (?,?,?)", (job_id, run_id, run_dir))
        # Mod-increment job ID
        job_id = (job_id % N_JOBS_MAX*numSubmits) + 1
        if job_id > job_count
            job_count = job_id
        end
    end

    submitScripts = IOStream[]
    for script in 1:numSubmits
        push!(submitScripts,open("$(script)-plot-analysis_submit_jobs.sh", "w"))
        println(submitScripts[script], """
        #!/bin/sh
        cd `dirname \$0`
        """)
    end

    for (job_id,) in execute(dbTempJobs, "SELECT DISTINCT job_id FROM plot_job_runs ORDER BY job_id")

        job_dir = joinpath(SCRIPT_PATH,analysisDir,"plotjobs", "$(job_id)")
        @assert !ispath(job_dir)
        mkpath(job_dir)

        # Get all run directories for this job
        run_dirs = [run_dir for (run_dir,) in execute(dbTempJobs,
            """
            SELECT run_dir FROM plot_job_runs
            WHERE job_id = ?
            """,
            (job_id,)
        )]

        n_cores = min(length(run_dirs), N_CORES_PER_JOB_MAX)

        if n_cores > n_cores_count
            n_cores_count = n_cores
        end

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
            #SBATCH --job-name=crispr-plots-$(analysisType)-$(job_id)
            #SBATCH --tasks=1
            #SBATCH --cpus-per-task=$(n_cores)
            #SBATCH --mem-per-cpu=$(mem_per_cpu)m
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

        execute(dbTempJobs, "INSERT INTO plot_jobs VALUES (?,?)", (job_id, job_dir))

        submitScript = (job_id % numSubmits) + 1
        println(submitScripts[submitScript], "sbatch $(job_sbatch)")
    end
    execute(dbTempJobs, "COMMIT")
    map(close,submitScripts)

    #run(`chmod +x submit_plot_jobs.sh`) # Make submit script executable
    @info "
    Sweep will be submitted via $(numSubmits) `plot-analysis_submit_jobs.sh` script(s).
    Each `plot-analysis_submit_jobs.sh` script submits $(job_count) jobs.
    Each job will use $(n_cores_count) cpus (cores) at most, where each cpu will use $(mem_per_cpu/1000)GB.
    Each job therefore will use at most $(n_cores_count*mem_per_cpu/1000)GB of memory in total.
    "
end


main()
