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

#function write_json_to_file(x, filename)
    #if ispath(filename)
        #error("$filename already exists. You should delete output, or run in a different directory.")
    #else
        #file = open(filename, "w")
        #print(file, JSON2.write(x))
        #file
    #end
#end

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
    s = sim.state

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

function write_extinction(sim, strains_table_name, id)
    execute(sim.db,
        "INSERT INTO $strains_table_name VALUES (?,?)",
        [id, sim.t]
    )
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
end
