#!/usr/bin/env julia

println("(Julia compilation delay...)")

include(joinpath("..","simulation","src","setup.jl"))

using SQLite

analysisType = ARGS[1]
analysisDir = "$(analysisType)"

SCRIPT_PATH = abspath(dirname(PROGRAM_FILE))
ROOT_PATH = joinpath(SCRIPT_PATH)
# ROOT_PATH = joinpath("/scratch/cgsb/pascual/armun/crispr/vary-fitness/15_MOI3/data-analysis")
dbSimPath = joinpath(SCRIPT_PATH, "..", "simulation", "sweep_db.sqlite")
# dbSimPath = joinpath("/Volumes/Yadgah/15_MOI3/sweep_db.sqlite")

# if analysisType == "peaks" && length(ARGS) < 2
#     error("`peaks` analysis needs two arguments: upper threshold, lower threshold")
# elseif analysisType == "walls-shannon" && length(ARGS) < 3
#     error("`walls` analysis needs three arguments: upper percent threshold, lower percent threshold, diversity threshold")
# end

# if analysisType == "walls-shannon" && !isfile(joinpath(SCRIPT_PATH,"shannon","shannon.sqlite"))
#     error("`/shannon/shannon.sqlite` is missing; please analyze shannon first.")
# end

# if analysisType == "walls-shannon" && parse(Float64,ARGS[2]) <= parse(Float64,ARGS[3])
#     error("upper percent threshold must be greater than lower percent threshold")
# end

# if analysisType == "walls-shannon" && parse(Float64,ARGS[2]) > 100
#     error("upper percent cannot be greater than 100%")
# end

# if analysisType == "match-diversity" && !isfile(joinpath(SCRIPT_PATH,"matches","matches.sqlite"))
#     error("`/matches/matches.sqlite` is missing; please analyze matches first.")
# end

ROOT_RUN_SCRIPT = joinpath(ROOT_PATH,analysisDir,"$(analysisType).jl")
ROOT_RUNMANY_SCRIPT = joinpath(ROOT_PATH,"src", "runmany.jl")
cd(SCRIPT_PATH)

# Number of SLURM jobs to generate
const N_JOBS_MAX = 1000
const N_PROCESSES = 15 # Half a node (14) is easier to get scheduled than a whole one
const N_CORES_PER_JOB_MAX = 15
const mem_per_cpu = 6000 # in MB 100MB = 1 GB

function main()
    # Root run directory
    # if ispath(joinpath(SCRIPT_PATH,analysisDir,"runs"))
    #     error("Please move or delete `/$(analysisType)/runs`.")
    # end

    # if ispath(joinpath(SCRIPT_PATH,analysisDir,"jobs"))
    #     error("Please move or delete `/$(analysisType)/jobs`.")
    # end

    # if !ispath(joinpath(SCRIPT_PATH,"..","simulation","runs"))
    #     error("`/../simulation/runs` is missing; please simulate time series first.")
    # end

    # Root job directory
    if !ispath(joinpath(SCRIPT_PATH,"..","simulation","jobs"))
        @info "Note that `/../simulation/jobs` is missing."
    end

    # Connect to simulation data
    dbSim = SQLite.DB(dbSimPath)

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
        #if !issubset(run_id,[1343, 1344, 1345])
            #continue
        #end
        # if run_id % 2 != 0
        #     continue
        # end
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
            /share/apps/julia/1.6.1/bin/julia $(ROOT_RUN_SCRIPT) $(run_id) $(argString...) &> output.txt
            """)
        end
        run(`chmod +x $(run_script)`) # Make run script executable
        run_count += 1
    end
    return numSubmits = Int64(ceil(run_count/(N_JOBS_MAX*N_PROCESSES)))
end

function generate_analysis_jobs(dbSim::DB,dbTempJobs::DB,numSubmits::Int64)
    println("Assigning analysis runs to jobs...")

    # Assign runs to jobs (round-robin)
    job_id = 1
    # n_cores_count = 0

    execute(dbTempJobs, "BEGIN TRANSACTION")

    for (run_id, run_dir) in execute(dbSim, "SELECT run_id, run_dir FROM runs ORDER BY replicate, combo_id")
        #if !issubset(run_id,[1343, 1344, 1345, 1455, 1550, 1160, 1657, 2950, 2954])
            #continue
        #end
        # if run_id % 2 != 0
        #     continue
        # end
        execute(dbTempJobs, "INSERT INTO job_runs VALUES (?,?,?)", (job_id, run_id, run_dir))
        # Mod-increment job ID
        job_id = mod(job_id,N_JOBS_MAX*numSubmits) + 1
    end

    submitScripts = IOStream[]
    for script in 1:numSubmits
        push!(submitScripts,open("$(script)_analysis-submit-jobs.sh", "w"))
        println(submitScripts[script], """
        #!/bin/sh
        cd `dirname \$0`
        """)
    end

    for (job_id,) in execute(dbTempJobs, "SELECT DISTINCT job_id FROM job_runs ORDER BY job_id")

        job_dir = joinpath(analysisDir,"jobs", "$(job_id)")
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

        if N_CORES_PER_JOB_MAX === nothing
            n_cores = min(length(run_dirs), N_PROCESSES)
        else
            n_cores = min(length(run_dirs), N_CORES_PER_JOB_MAX)
        end

        # if n_cores > n_cores_count
        #     n_cores_count = n_cores
        # end

        # Write out list of runs
        open(joinpath(job_dir, "runs.txt"), "w") do f
            for run_dir in run_dirs
                run_script = joinpath(ROOT_PATH, analysisDir, run_dir, "run.sh")
                println(f, run_script)
            end
        end

        # Create job sbatch file
        job_sbatch = joinpath(job_dir, "job.sbatch")
        open(job_sbatch, "w") do f
            print(f, """
            #!/bin/sh
            #SBATCH --nodes=1
            #SBATCH --ntasks-per-node=$(n_cores)
            #SBATCH --cpus-per-task=1
            #SBATCH --job-name=$(job_id)$(analysisType)
            #SBATCH --mem-per-cpu=$(mem_per_cpu)m
            #SBATCH --time=1-12:00:00
            #SBATCH --chdir=$(joinpath(ROOT_PATH, job_dir))
            #SBATCH --output=output.txt
            #SBATCH --mail-type=END,FAIL,REQUEUE
            #SBATCH --mail-user=al8784@nyu.edu
            module purge
            module load julia/1.6.1
            ~/julia/my-julia $(ROOT_RUNMANY_SCRIPT) $(n_cores) runs.txt
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
    Sweep will be submitted via $(numSubmits) `analysis-submit-jobs.sh` script(s).
    Each `analysis_submit_jobs.sh` script submits $(N_JOBS_MAX) jobs.
    Each job will use $(n_cores) cpus (cores) at most, where each cpu will use $(mem_per_cpu/1000)GB.
    Each job therefore will use at most $(n_cores*mem_per_cpu/1000)GB of memory in total.
    "
end


main()
