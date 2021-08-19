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

println("(Annoying Julia compilation delay...)")

using Random
using SQLite
import SQLite.DBInterface.execute
using DelimitedFiles

# Get relevant paths and cd to the script path.
# NB: use actual relative locations of varmodel3 root relative to your script.
include("../../preamble.jl")
SCRIPT_PATH = abspath(dirname(PROGRAM_FILE))
ROOT_PATH = abspath(joinpath(SCRIPT_PATH, "..", "..")) #CHANGE THESE DIRECTORIES ACCORDINGLY
ROOT_RUN_SCRIPT = joinpath(ROOT_PATH, "run.jl")
ROOT_RUNMANY_SCRIPT = joinpath(ROOT_PATH, "runmany.jl")
cd(SCRIPT_PATH)

# Number of replicates for each parameter combination
const N_REPLICATES = 2

# Number of SLURM jobs to generate
const N_JOBS_MAX = 100
const N_CORES_PER_JOB_MAX = 14 # Half a node, easier to get scheduled than a whole one

function main()
    # Root run directory
    if ispath("runs")
        error("`runs` already exists; please move or delete.")
    end
    mkdir("runs")

    # Root job directory
    if ispath("jobs")
        error("`jobs` already exists; please move or delete.")
    end
    mkdir("jobs")

    # Database of experiment information
    if ispath("sweep_db.sqlite")
        error("`sweep_db.sqlite` already exists; please move or delete")
    end
    db = SQLite.DB(joinpath("sweep_db.sqlite")) # the function of this database
    # is to log run and job ids of individual simulation directory names
    execute(db, "CREATE TABLE meta (key, value)")
    execute(db, "CREATE TABLE param_combos (combo_id INTEGER, mutation_rate REAL, transmissibility REAL)")
    execute(db, "CREATE TABLE runs (run_id INTEGER, combo_id INTEGER, replicate INTEGER, rng_seed INTEGER, run_dir TEXT, params TEXT)")
    execute(db, "CREATE TABLE jobs (job_id INTEGER, job_dir TEXT)")
    execute(db, "CREATE TABLE job_runs (job_id INTEGER, run_id INTEGER)")

    generate_runs(db)
    generate_jobs(db)
end

function generate_runs(db) # This function generates the directories
    # for the individual parameter sets and corresponding replicates. It also
    # generates shell scripts for each run and corresponding parameter file.

    # System random device used to generate seeds
    seed_rng = RandomDevice()

    # Base parameter set, copied/modified for each combination/replicate
    base_params = init_base_params()
    validate(base_params)
    execute(db, "INSERT INTO meta VALUES (?, ?)", ("base_params", pretty_json(base_params)))

    # Loop through parameter combinations and replicates, generating a run directory
    # `runs/c<combo_id>/r<replicate>` for each one.
    combo_id = 1
    run_id = 1
    for mutation_rate in (0.5e-8, 1.0e-8)
        for transmissibility in (0.25, 0.5)
            println("Processing c$(combo_id): mutation_rate = $(mutation_rate), transmissibility = $(transmissibility)")

            execute(db, "INSERT INTO param_combos VALUES (?, ?, ?)", (combo_id, mutation_rate, transmissibility))

            for replicate in 1:N_REPLICATES
                rng_seed = rand(seed_rng, 1:typemax(Int64))
                params = Params(
                    base_params;
                    rng_seed = rng_seed,
                    mutation_rate = mutation_rate,
                    transmissibility = transmissibility
                ) # <-- via JSON"1". this params function either needs to be adapted for CRISPR code,
                # or I need to adapt CRISPR code to adapt to functuon

                run_dir = joinpath("runs", "c$(combo_id)", "r$(replicate)")
                @assert !ispath(run_dir)
                mkpath(run_dir)

                # Generate parameters file
                params_json = pretty_json(params)
                open(joinpath(run_dir, "parameters.json"), "w") do f
                    println(f, params_json)
                end


                # I think this needs a "joinpath" code block for model.jl,
                # structures.jl, output.jl, util.jl, model.jl??????? ##########


                # Generate shell script to perform a single run
                run_script = joinpath(run_dir, "run.sh")
                open(run_script, "w") do f
                    print(f, """
                    #!/bin/sh
                    cd `dirname \$0`
                    julia $(ROOT_RUN_SCRIPT) parameters.json &> output.txt
                    """)# what does \$0 mean???? # ROOT_RUN_SCRIPT = run.jl which is analogous to main.jl
                    # Double check what runs.jl looks like. Don't forget to put "include()" scripts into a "preamble.jl"
                    # what is --check-bounds=no -03???
                    # this creates a shell script for each simulation... this is not an sbatch
                    #EACH OF THESE IN AN INDIVIDUAL DIRECTORY IS SAVED LINE-BY-LINE in RUNS.TXTTTT.
                    # what does &> mean???
                    #check julia --check-bounds=no -O3 $(ROOT_RUN_SCRIPT) parameters.json &> output.txt on desktop
                end
                run(`chmod +x $(run_script)`) # Make run script executable
                # WHY NOT 777 but rather +x? TEST THIS ON MIDWAY

                # Save all run info (including redundant stuff for reference) into DB
                execute(db, "INSERT INTO runs VALUES (?, ?, ?, ?, ?, ?)", (run_id, combo_id, replicate, rng_seed, run_dir, params_json))

                run_id += 1
            end
            combo_id += 1
        end
    end
end

function generate_jobs(db)
    println("Assigning runs to jobs...")

    # Assign runs to jobs (round-robin)
    job_id = 1
    for (run_id, run_dir) in execute(db, "SELECT run_id, run_dir FROM runs ORDER BY replicate, combo_id")
        execute(db, "INSERT INTO job_runs VALUES (?,?)", (job_id, run_id))

        # Mod-increment job ID
        job_id = (job_id % N_JOBS_MAX) + 1
    end

    # Create job directories containing job scripts and script to submit all jobs
    submit_file = open("submit_jobs.sh", "w")
    println(submit_file, """
    #!/bin/sh
    cd `dirname \$0`
    """)
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
        n_cores = min(length(run_dirs), N_CORES_PER_JOB_MAX)

        # Write out list of runs
        open(joinpath(job_dir, "runs.txt"), "w") do f
            for run_dir in run_dirs
                run_script = joinpath(SCRIPT_PATH, run_dir, "run.sh")
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
            #SBATCH --job-name=var-$(job_id)
            #SBATCH --tasks=1
            #SBATCH --cpus-per-task=$(n_cores)
            #SBATCH --mem-per-cpu=2000m
            #SBATCH --time=4:00:00
            #SBATCH --chdir=$(joinpath(SCRIPT_PATH, job_dir))
            #SBATCH --output=output.txt
            module purge
            # Uncomment this to use the Midway-provided Julia:
            module load julia
            julia $(ROOT_RUNMANY_SCRIPT) $(n_cores) runs.txt
            """) # runs.txt is for parallel processing
            # WHAT DOES OUTPUT.TXT LOOK LIKE??
        end
        run(`chmod +x $(job_sbatch)`) # Make run script executable (for local testing)

        execute(db, "INSERT INTO jobs VALUES (?,?)", (job_id, job_dir,))
        println(submit_file, "sbatch $(job_sbatch)")
    end
    close(submit_file)
    run(`chmod +x submit_jobs.sh`) # Make submit script executable
end

function pretty_json(params) #LEARN WHAT THIS FUNCTION IS DOING...
    d = Dict(fn => getfield(params, fn) for fn in fieldnames(typeof(params)))
    io = IOBuffer()
    JSON.print(io, d, 2)
    String(take!(io))
end

function init_base_params()
    Params(
        t_final = 2000.0,

        t_output = 1.0,

        rng_seed = nothing,

        enable_output = true,

        n_bstrains = 1,

        n_hosts_per_bstrain = 100,

        n_vstrains = 1,

        n_particles_per_vstrain = 100,

        n_protospacers = 15,

        u_n_spacers_max = 10,

        p_crispr_failure_prob = 1e-05,

        q_spacer_acquisition_prob = 1e-05,

        r_growth_rate = 1,

        K_carrying_capacity = 316227.7660168379,

        beta_burst_size = 50,

        phi_adsorption_rate = 1e-01,

        m_viral_decay_rate = 0.1,

        mu_viral_mutation_rate = 1e-06,

        d_death_rate = 0,

        g_immigration_rate = 1,
    )
end

main()
