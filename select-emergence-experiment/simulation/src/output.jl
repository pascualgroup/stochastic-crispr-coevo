### OUTPUT FUNCTIONS ###

function initialize_database()
    print(cd(pwd))
    if isfile("output.sqlite")
        error("output.sqlite already exists; delete first")
    end
    db = DB("output.sqlite")

    execute(db, "CREATE TABLE meta (key, value);")
    execute(db, """
        CREATE TABLE summary (
            t REAL,
            microbial_abundance INTEGER,
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
        CREATE TABLE bextinctions (
            bstrain_id INTEGER,
            t_extinction REAL
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
        CREATE TABLE vextinctions (
            vstrain_id INTEGER,
            t_extinction REAL
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

function write_periodic_output(sim) # this only writes the summary
    s = sim.state
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

function write_strains(sim, strains_table_name, spacers_table_name,
                        other_table_name, ids, spacers, parameter1, parameter2)
    @debug "write_strains" ids
    for i = 1:lastindex(ids)
        write_strain(sim, strains_table_name, ids[i], 0, 0)
        write_spacers(sim, spacers_table_name, ids[i], spacers[i])
    end
    if other_table_name == "bgrowthrates"
        for i = 1:lastindex(ids)
            write_growth_rate(sim, other_table_name, ids[i], parameter1[i], parameter2[i])
        end
    end
end

function write_strain(sim, strains_table_name, id, parent_id, other_id)
    execute(sim.db,
        "INSERT INTO $strains_table_name VALUES (?,?,?,?)",
        [sim.t, id, parent_id, other_id]
    )
end

function write_extinction(sim, strains_table_name, id)
    execute(sim.db,
        "INSERT INTO $strains_table_name VALUES (?,?)",
        [id, sim.t]
    )
end

function write_spacers(sim, spacers_table_name, id, spacers)
    @debug "write_spacers" id spacers
    stmt = Stmt(sim.db, "INSERT INTO $spacers_table_name VALUES (?,?)")
    for spacer_id = spacers
        execute(stmt, [Int64(id), Int64(spacer_id)])
    end
end

function write_growth_rate(sim, table_name, id, lAllele, gRate)
    @debug "write_growth_rates" id gRate
    stmt = Stmt(sim.db, "INSERT INTO $table_name VALUES (?,?,?)")
    execute(stmt, [Int64(id), Float64(lAllele), Float64(gRate)])
end

function write_abundances(sim, table_name, ids, abundance)
    t = sim.t
    @debug "write_abundances" t ids abundance
    @debug "lastindex(ids)" lastindex(ids)
    stmt = Stmt(sim.db, "INSERT INTO $table_name VALUES (?,?,?)")
    for i = 1:lastindex(ids)
        @debug "lastindex(ids)" lastindex(ids)
        execute(stmt, [t, ids[i], abundance[i]])
    end
end
