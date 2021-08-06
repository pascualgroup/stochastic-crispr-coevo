using Dates
using JSON2

using SQLite: DB, Stmt, bind!
using SQLite.DBInterface: execute

### OUTPUT FUNCTIONS ###

function initialize_database()
    if isfile("output.sqlite")
        error("output.sqlite already exists; delete first")
    end
    db = DB("output.sqlite")

    execute(db, "CREATE TABLE meta (key, value);")
    execute(db, """
        CREATE TABLE summary (
            t REAL,
            bacterial_abundance INTEGER,
            viral_abundance INTEGER
        );
    """)

    execute(db, """
        CREATE TABLE bstrains (
            t_creation REAL,
            bstrain_id INTEGER,
            parent_bstrain_id INTEGER,
            infecting_vstrain_id INTEGER
        )
    """)
    execute(db, """
        CREATE TABLE bspacers (
            bstrain_id INTEGER,
            spacer_id INTEGER
        )
    """)
    execute(db, """
        CREATE TABLE babundance (
            t REAL,
            bstrain_id INTEGER,
            abundance INTEGER
        )
    """)

    execute(db, """
        CREATE TABLE vstrains (
            t_creation REAL,
            vstrain_id INTEGER,
            parent_vstrain_id INTEGER,
            infected_bstrain_id INTEGER
        )
    """)
    execute(db, """
        CREATE TABLE vpspacers (
            vstrain_id INTEGER,
            spacer_id INTEGER
        )
    """)
    execute(db, """
        CREATE TABLE vabundance (
            t REAL,
            vstrain_id INTEGER,
            abundance INTEGER
        )
    """)

    db
end

function write_json_to_file(x, filename)
    if ispath(filename)
        error("$filename already exists. You should delete output, or run in a different directory.")
    else
        file = open(filename, "w")
        print(file, JSON2.write(x))
        file
    end
end

#function open_csv(prefix, header...)
    #filename = "$prefix.csv"
    #if ispath(filename)
        #error("$filename already exists. You should delete output, or run in a different directory.")
    #else
        #file = open(filename, "w")
        #write_csv(file, header...)
        #file
    #end
#end

#function write_csv(file, row...)
    #for i = 1:lastindex(row)
        #print(file, row[i])
        #if i < lastindex(row)
            #print(file, ",")
        #end
    #end
    #print(file, "\n")
#end

function write_periodic_output(sim) # this only writes the summary
    #and the abundances of the different strains
    #s = sim.state

    #write_summary(sim.summary_file, sim.t, s)

    #write_abundances(
        #s.bstrains.abundance_file, sim.t, s.bstrains.ids, s.bstrains.abundance
    #)
    #write_abundances(
        #s.vstrains.abundance_file, sim.t, s.vstrains.ids, s.vstrains.abundance
    #)

    #flush(s.bstrains.strain_file)
    #flush(s.bstrains.spacers_file)
    #flush(s.vstrains.strain_file)
    #flush(s.vstrains.spacers_file)


    write_summary(sim)
    write_abundances(sim, "babundance", s.bstrains.ids, s.bstrains.abundance)
    write_abundances(sim, "vabundance", s.vstrains.ids, s.vstrains.abundance)
end

function write_summary(sim)

    db = sim.db
    t = sim.t
    state = sim.state
    babundance = state.bstrains.total_abundance
    vabundance = state.vstrains.total_abundance

    execute(db,
        "INSERT INTO summary VALUES (?,?,?)",
        [t, babundance, vabundance]
    )
end

function write_strains(sim, strains_table_name, spacers_table_name, ids, spacers)
    #function write_strains(strain_file, spacers_file, t_creation, ids, spacers)
    ## where does t_creation go??
    @debug "write_strains" ids
    for i = 1:lastindex(ids)
        write_strain(sim, strains_table_name, ids[i], 0, 0)
        write_spacers(sim, spacers_table_name, ids[i], spacers[i])
    end
end

function write_strain(sim, strains_table_name, id, parent_id, other_id) #where does t_creation go?
    #function write_strain(file, t_creation, id, parent_id, other_id) #only initialization of simulation
    # and events such as "acquire spacer" and "mutate virus" use this write new strains.
    #write_periodic_output is used to update abundances thereafter.
    # THIS FUNCTION ONLY INSERTS a new strain into the file

    execute(sim.db,
        "INSERT INTO $strains_table_name VALUES (?,?,?,?)",
        [sim.t, id, parent_id, other_id]
    )
end

function write_spacers(sim, spacers_table_name, id, spacers)
    @debug "write_spacers" id spacers
    stmt = Stmt(sim.db, "INSERT INTO $spacers_table_name VALUES (?,?)")
    for spacer_id = spacers
        execute(stmt, [Int64(id), Int64(spacer_id)])
    end
end

function write_abundances(sim, table_name, ids, abundance)
    @debug "write_abundances" t ids abundance
    @debug "lastindex(ids)" lastindex(ids)
    stmt = Stmt(sim.db, "INSERT INTO $table_name VALUES (?,?,?)")
    t = sim.t
    for i = 1:lastindex(ids)
        @debug "lastindex(ids)" lastindex(ids)
        execute(stmt, [t, ids[i], abundance[i]])
    end
    #flush(file)
end






## Output Functions from var model ##








"""
This file defines types and functions for database output.
Currently, it is shared across model variants, but in the future may need to be
broken up to be variant-specific.
The `write_output!()` function defined here calls three functions defined by the
code for individual model variants:
* `write_summary()`: writes out summary data
* `write_host_samples()`: writes out sampled hosts and infections
* `write_gene_strain_counts()`: writes out gene and strain counts
There is a sophisticated Julia system called Tables.jl that would allow the code
here to be less repetitive. With the goal of not introducing too many layers
of abstraction, this code uses SQLite more directly; a future Julia-oriented
maintainer may want to convert the code.
"""

using SQLite: DB, Stmt
import SQLite.DBInterface.execute

"""
Database type encapsulating SQLite.DB and prepared insert statements for tables.
Prepared insert statements are used for the sake of performance; they allow
SQLite to avoid parsing and compiling the statements with every inserted row.
"""
struct VarModelDB
    db::DB
    meta::Stmt
    summary::Stmt
    gene_strain_counts::Stmt
    sampled_hosts::Stmt
    sampled_infections::Stmt
    sampled_infection_genes::Stmt
end

struct CRISPRModelDB
    db::DB
    meta::Stmt
    summary::Stmt
    gene_strain_counts::Stmt
    total_microbe_abundance::Stmt
    total_virus_abundance::Stmt
    sampled_infection_genes::Stmt
end

"""
Type encapsulating various summary statistics gathered between summary periods.
"""
@with_kw mutable struct SummaryStats
    start_datetime::DateTime
    n_events::Int = 0
    n_bites::Int = 0
    n_infected_bites::Int = 0
    n_infected_bites_with_space::Int = 0
    n_transmitting_bites::Int = 0
    n_transmissions::Int = 0
end

"""
SummaryStats constructor.
"""
function SummaryStats()
    SummaryStats(start_datetime = now())
end

"""
Reset counts and time in a SummaryStats object.
Called after summary statistics are written out in model variant-specific
`write_summary()` function.
"""
function reset!(stats::SummaryStats, start_datetime)
    stats.start_datetime = start_datetime
    stats.n_bites = 0
    stats.n_infected_bites = 0
    stats.n_infected_bites_with_space = 0
    stats.n_transmitting_bites = 0
    stats.n_transmissions = 0
end

"""
Pass commands issued to a `VarModelDB` on to the underlying `SQLite.DB`.
"""
function execute(db::VarModelDB, cmd)
    execute(db.db, cmd)
end

"""
Initialize database.
Note: if the number of columns in a table is modified, the corresponding
`make_insert_statement()` call must be updated to match the number of columns.
(This could be automated if it seems worth it.)
"""
function initialize_database()
    if isfile(P.output_db_filename)
        error("$(P.output_db_filename) already exists; delete first")
    end
    db = DB(P.output_db_filename)

    # Note: all of this could be done with the Julia Tables package, which
    # might make it less cumbersome, but also introduces conceptual overhead.

    execute(db, "CREATE TABLE meta (key, value);")

    execute(db, """
        CREATE TABLE summary (
            time INTEGER,
            n_infections_liver INTEGER,
            n_infections_active INTEGER,
            n_infections INTEGER,
            n_infected_liver INTEGER,
            n_infected_active INTEGER,
            n_infected INTEGER,
            n_bites INTEGER,
            n_infected_bites INTEGER,
            n_infected_bites_with_space INTEGER,
            n_transmitting_bites INTEGER,
            n_transmissions INTEGER,
            exec_time INTEGER
        );
    """)

    execute(db, """
        CREATE TABLE gene_strain_counts (
            time INTEGER,
            n_circulating_genes INTEGER,
            n_circulating_strains INTEGER
        );
    """)

    execute(db, """
        CREATE TABLE sampled_hosts (
            time INTEGER,
            id INTEGER,
            birth_time REAL,
            death_time REAL,
            n_infections_liver INTEGER,
            n_infections_active INTEGER
        )
    """)

    execute(db, """
        CREATE TABLE sampled_infections (
            time INTEGER,
            host_id INTEGER,
            infection_id INTEGER,
            infection_time REAL,
            strain_id INTEGER,
            expression_index INTEGER
        );
    """)

    allele_columns = join(["allele_id_$(i) INTEGER" for i in 1:P.n_loci], ", ")
    execute(db, """
        CREATE TABLE sampled_infection_genes(
            infection_id INTEGER,
            expression_index INTEGER,
            $(allele_columns),
            UNIQUE(infection_id, expression_index) ON CONFLICT IGNORE
        );
    """)

    VarModelDB(
        db,
        make_insert_statement(db, "meta", 2),
        make_insert_statement(db, "summary", 13),
        make_insert_statement(db, "gene_strain_counts", 3),
        make_insert_statement(db, "sampled_hosts", 6),
        make_insert_statement(db, "sampled_infections", 6),
        make_insert_statement(db, "sampled_infection_genes", 2 + P.n_loci)
    )
end

"""
Construct an prepared insert statement for a particular table.
The statement covers all columns in the table, and `n_columns` must match the
actual number of columns in the table.
"""
function make_insert_statement(db, table_name, n_columns)
    qmarks = join(repeat(["?"], n_columns), ",")
    Stmt(
        db,
        "INSERT INTO $(table_name) VALUES ($(qmarks))"
    )
end


### OUTPUT FUNCTIONS ###

function write_output!(db, t, s, stats)
    if P.t_burnin !== missing && t < P.t_burnin
        return
    end

    if t % minimum(
        (P.summary_period, P.host_sampling_period, P.gene_strain_count_period,)
    ) == 0
        println("t = $(t)")

#         println("write_output!($(t))")

        execute(db, "BEGIN TRANSACTION")

        if t % P.summary_period == 0
            write_summary(db, t, s, stats)
        end

        if t % P.host_sampling_period == 0
            write_host_samples(db, t, s)
        end

        if t % P.gene_strain_count_period == 0
            write_gene_strain_counts(db, t, s)
        end

        execute(db, "COMMIT")
    end
end

"""
Write output to `summary` table
"""
function write_summary(db, t, s, stats)
#     println("write_summary($(t))")

    # Compute number of infections (liver, active, and both) with a simple tally/sum.
    n_infections_liver = sum(length(host.liver_infections) for host in s.hosts)
    n_infections_active = sum(length(host.active_infections) for host in s.hosts)
    n_infections = n_infections_liver + n_infections_active

    # Compute number of individuals with infections (liver, active, or either).
    n_infected_liver = sum(length(host.liver_infections) > 0 for host in s.hosts)
    n_infected_active = sum(length(host.active_infections) > 0 for host in s.hosts)
    n_infected = sum(
        length(host.liver_infections) > 0 || length(host.active_infections) > 0
        for host in s.hosts
    )

    # Compute elapsed time in seconds.
    next_datetime = now()
    exec_time = Dates.value(next_datetime - stats.start_datetime) / 1000.0

    # Write to summary table.
    execute(db.summary, (
        t,
        n_infections_liver,
        n_infections_active,
        n_infections,
        n_infected_liver,
        n_infected_active,
        n_infected,
        stats.n_bites,
        stats.n_infected_bites,
        stats.n_infected_bites_with_space,
        stats.n_transmitting_bites,
        stats.n_transmissions,
        exec_time
    ))

    # Reset counters and elapsed time.
    reset!(stats, next_datetime)
end

"""
Write output for periodically sampled hosts.
"""
function write_host_samples(db, t, s)
#     println("write_host_samples($(t))")

    # Sample `host_sample_size` hosts randomly (without replacement).
    hosts = sample(s.hosts, P.host_sample_size, replace = false)

    # For each host, write out birth/death time and each infection.
    for host in hosts
        execute(
            db.sampled_hosts,
            (
                t, Int64(host.id), host.t_birth, host.t_death,
                length(host.liver_infections), length(host.active_infections)
            )
        )

        for infection in host.liver_infections
            write_infection(db, t, host, infection)
        end

        for infection in host.active_infections
            write_infection(db, t, host, infection)
        end
    end
end

"""
Write output for a single infection from a sampled host.
Infections are written to the `sampled_infections` table, and the genes for the
infections are written to `sampled_infection_genes`.
`expression_index` indicates which gene is currently expressed. Liver-stage
infections are indicated using a SQLite `NULL` (Julia `missing`) for
`expression_index`.
For each infection, the `sampled_infection_genes` table contains one row for
each expression index, with the final columns containing the allele IDs for each
locus.
"""
function write_infection(db, t, host, infection)
#     println("write_infection($(t), $(host.id), $(infection.id))")
    execute(
        db.sampled_infections,
        (
            t, Int64(host.id), Int64(infection.id), infection.t_infection, Int64(infection.strain_id),
            infection.expression_index == 0 ? missing : Int64(infection.expression_index)
        )
    )
    for i in 1:P.n_genes_per_strain
        execute(db.sampled_infection_genes, vcat([infection.id, i], infection.genes[:,i]))
    end
end

"""
Write gene and strain counts to the `gene_strain_counts` table.
Counts are not maintained dynamically during the simulation; this function
simply scans all host infections and assembles sets of all genes and all
strains.
"""
function write_gene_strain_counts(db, t, s)
#     println("write_gene_strain_counts($(t))")
    genes::Set{Gene} = Set()
    strains::BitSet = BitSet()

    for host in s.hosts
        count_genes_and_strains!(genes, strains, host.liver_infections)
        count_genes_and_strains!(genes, strains, host.active_infections)
    end

    execute(db.gene_strain_counts, (t, length(genes), length(strains)))
end

"""
Count genes and strains for a particular list of infections in a particular host.
This function simply adds genes and strains from each infection to corresponding
sets.
"""
function count_genes_and_strains!(genes, strains, infections)
#     println("count_circulating_genes_and_strains()")
    for infection in infections
        for i in 1:P.n_genes_per_strain
            push!(genes, infection.genes[:,i])
        end
        push!(strains, infection.strain_id)
    end
end
