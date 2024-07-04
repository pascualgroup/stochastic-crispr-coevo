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

println("(Julia compilation delay...)")

include("src/setup.jl")
include("src/structures.jl")
include("src/util.jl")

using SQLite
using DataFrames

SCRIPT_PATH = abspath(dirname(PROGRAM_FILE))
ROOT_PATH = joinpath("/scratch/cgsb/pascual/armun/crispr/vary-q-mu/sweep-6-5-2024_50reps/simulation")
ROOT_RUN_SCRIPT = joinpath(ROOT_PATH, "src", "main-sweep.jl")
ROOT_RUNMANY_SCRIPT = joinpath(ROOT_PATH, "src", "runmany.jl")
cd(SCRIPT_PATH)

# Number of replicates for each parameter combination
const N_REPLICATES = 50

# Number of SLURM jobs to generate
const N_JOBS_MAX = 500
const N_PROCESSES = 16 # Half a node (14) is easier to get scheduled than a whole one
const N_CORES_PER_JOB_MAX = 16
const mem_per_cpu = 8000 # in MB 100MB = 1 GB

db = SQLite.DB(joinpath(SCRIPT_PATH,"sweep_db.sqlite")) # the function of this database
# is to log run and job ids of individual simulation directory names
json_str = read(joinpath(SCRIPT_PATH,"..","sweep-parameters.json"), String)
paramStringTrunc = JSON.parse(json_str)
# Deleting any key will designate the parameter as a base (unchanging parameter)
# with the exception of the seed, which is generated for every individual replicate
delete!(paramStringTrunc, "rng_seed")
delete!(paramStringTrunc, "enable_output")

function main()
    # Root run directory
    if ispath(joinpath(SCRIPT_PATH,"runs"))
        error("`/simulation/runs` already exists; please move or delete.")
    end
    mkdir("runs")

    # Root job directory
    if ispath(joinpath(SCRIPT_PATH,"jobs"))
        error("`/simulation/jobs` already exists; please move or delete.")
    end
    mkdir("jobs")

    # Database of experiment information
    if isfile(joinpath(SCRIPT_PATH,"sweep_db.sqlite)"))
        error("`/simulation/sweep_db.sqlite` already exists; please move or delete")
    end

    paramKeys = ["$(k) REAL" for (k,v) in paramStringTrunc]
    paramColumns = join(paramKeys,", ")

    execute(db, "CREATE TABLE meta (key, value)")
    execute(db, "CREATE TABLE param_combos (combo_id INTEGER, $(paramColumns))")
    execute(db, "CREATE TABLE runs (run_id INTEGER, combo_id INTEGER, replicate INTEGER, rng_seed INTEGER, run_dir TEXT, params TEXT)")
    execute(db, "CREATE TABLE jobs (job_id INTEGER, job_dir TEXT)")
    execute(db, "CREATE TABLE job_runs (job_id INTEGER, run_id INTEGER)")

    numCombos = generate_runs(db)
    generate_jobs(db,numCombos)
end

function generate_runs(db::DB) # This function generates the directories
    # for the individual parameter sets and corresponding replicates. It also
    # generates shell scripts for each run and corresponding parameter file.

    # System random device used to generate seeds
    seed_rng = RandomDevice()

    json_str = read(joinpath(SCRIPT_PATH,"..","sweep-parameters.json"), String)
    paramString = JSON.parse(json_str)
    paramSymb = Dict((Symbol(k), v) for (k, v) in paramString)
    base_params = init_params(paramSymb)
    validate(base_params)

    # execute(db, "BEGIN TRANSACTION")
    execute(db, "INSERT INTO meta VALUES (?, ?)", ("base_params", pretty_json(base_params)))
    # execute(db, "COMMIT")

    paramKeys = [k for (k,) in paramStringTrunc]
    paramVals = [v for (nothing, v) in paramStringTrunc]
    combos = collect(Base.product(paramVals...))

    combosMod = let
        paramValsMod = [v for (nothing, v) in paramStringTrunc] 
        collect(Base.product(paramValsMod...))
    end
    
    # execute(db, "BEGIN TRANSACTION")
    run_id = 1
    execute(db, "BEGIN TRANSACTION")
    for combo_id in collect(1:length(combos))
        println("Processing parameter combination $(combo_id)")
        paramSpace = ["?" for i in 1:(length(paramVals)+1)]
        paramSpaces = join(paramSpace,", ")
        execute(db, "INSERT INTO param_combos VALUES ($(paramSpaces))",
        (combo_id,combosMod[combo_id]...)
        )

        paramSymbTrunc = Dict(map(Symbol,paramKeys).=>combos[combo_id])

        for replicate in 1:N_REPLICATES
            # t0 = now()
            rng_seed = rand(seed_rng, 1:typemax(Int64))
            params = Params(base_params;rng_seed = rng_seed,
                        paramSymbTrunc...
                        )

            run_dir = joinpath("runs", "c$(combo_id)", "r$(replicate)")
            @assert !ispath(run_dir)
            mkpath(run_dir)
            # tp = now()
            # Generate parameters file
            params_json = pretty_json(params)
            open(joinpath(run_dir, "parameters.json"), "w") do f
                println(f, params_json)
            end
            # println("$(replicate):json takes $(Dates.value(tp - t0) / 1000.0)")

            # tscript = now()
            # Generate shell script to perform a single run
            run_script = joinpath(run_dir, "run.sh")
            open(run_script, "w") do f
                print(f, """
                #!/bin/sh
                cd `dirname \$0`
                /share/apps/julia/1.6.1/bin/julia $(ROOT_RUN_SCRIPT) parameters.json &> output.txt
                """)
            end
            # println("$(replicate):script takes $(Dates.value(tscript - t0) / 1000.0)")
            run(`chmod +x $(run_script)`) # Make run script executable

            # Save all run info (including redundant stuff for reference) into DB
            execute(db, "INSERT INTO runs VALUES (?, ?, ?, ?, ?, ?)", (run_id, combo_id, replicate, rng_seed, run_dir, params_json))
            run_id += 1
        end
    end
    execute(db, "COMMIT")
    # execute(db, "COMMIT")
    return Int64(length(combos))
end

function generate_jobs(db::DB,numCombos::Int64)
    println("Assigning runs to jobs...")

    numSubmits = Int64(ceil(numCombos*N_REPLICATES/(N_JOBS_MAX*N_PROCESSES)))

    # Assign runs to jobs (round-robin)
    job_id = 1
    n_cores = 0

    # execute(db, "BEGIN TRANSACTION")
    execute(db, "BEGIN TRANSACTION")
    for (run_id,) in execute(db, "SELECT run_id FROM runs ORDER BY combo_id, replicate")
        execute(db, "INSERT INTO job_runs VALUES (?,?)", (job_id, run_id))
        # Mod-increment job ID
        job_id = mod(job_id,N_JOBS_MAX*numSubmits) + 1
    end

    submitScripts = IOStream[]
    for script in 1:numSubmits
        push!(submitScripts,open("$(script)_submit-jobs.sh", "w"))
        println(submitScripts[script], """
        #!/bin/sh
        cd `dirname \$0`
        """)
    end

    # Create job directories containing job scripts and script to submit all jobs
    for (job_id,) in execute(db, "SELECT DISTINCT job_id FROM job_runs ORDER BY job_id")
        job_dir = joinpath("jobs", "$(job_id)")
        mkpath(job_dir) # this is the directory of a "job". Here a inside a job id folder there is
        # a corresponding runs.txt file that is opened with runmany.jl in order to parallelize individual runs

        # Get all run directories for this job
        run_dirs = [run_dir for (run_dir,) in execute(db,
            """
            SELECT run_dir FROM job_runs, runs
            WHERE job_runs.job_id = ?
            AND runs.run_id = job_runs.run_id
            """,
            (job_id,)
        )]

        if N_CORES_PER_JOB_MAX === nothing
            n_cores = min(length(run_dirs), N_PROCESSES)
        else
            n_cores = min(length(run_dirs), N_CORES_PER_JOB_MAX)
        end
        # Write out list of runs
        open(joinpath(job_dir, "runs.txt"), "w") do f
            for run_dir in run_dirs
                run_script = joinpath(ROOT_PATH, run_dir, "run.sh")
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
            #SBATCH --job-name=$(job_id)crispr
            #SBATCH --mem-per-cpu=$(mem_per_cpu)m
            #SBATCH --time=1-12:00:00
            #SBATCH --chdir=$(joinpath(ROOT_PATH, job_dir))
            #SBATCH --output=output.txt
            #SBATCH --mail-type=FAIL
            #SBATCH --mail-user=al8784@nyu.edu
            module purge
            module load julia/1.6.1
            ~/julia/my-julia $(ROOT_RUNMANY_SCRIPT) $(n_cores) runs.txt
            """) # runs.txt is for parallel processing
        end
        run(`chmod +x $(job_sbatch)`) # Make run script executable (for local testing)

        execute(db, "INSERT INTO jobs VALUES (?,?)", (job_id, job_dir,))

        submitScript = mod(job_id,numSubmits) + 1
        println(submitScripts[submitScript], "sbatch $(job_sbatch)")
    end
    execute(db, "COMMIT")
    map(close,submitScripts)

    #run(`chmod +x submit_jobs.sh`) # Make submit script executable
    @info "
    Sweep will be submitted via $(numSubmits) `submit_jobs.sh` script(s).
    Each `submit_jobs.sh` script submits $(N_JOBS_MAX) jobs.
    Each job will use $(n_cores) cpus at most, where each cpu will use $(mem_per_cpu/1000)GB.
    Each job therefore will use at most $(n_cores*mem_per_cpu/1000)GB of memory in total.
    "
end

function pretty_json(params)
    d = Dict(fn => getfield(params, fn) for fn in fieldnames(typeof(params)))
    io = IOBuffer()
    JSON.print(io, d, 2)
    String(take!(io))
end

main()
