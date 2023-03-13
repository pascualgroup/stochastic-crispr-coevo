function potentialStructure!(matchStructure::hierarchy)
    dbTempSim = matchStructure.dbSim
    dbTempMatch = matchStructure.dbMatch
    dbTempTri = matchStructure.dbTri
    matchIDs = Dict(Vector([Int64(spacerID) for (spacerID,) in
                    execute(dbTempTri,"SELECT phenotype
                    FROM vmatch_phenotypes
                    WHERE match_id = $(matchID)
                    ORDER BY phenotype")]) => Int64(matchID)
                        for matchID in
                            [matchID for (matchID,) in
                            execute(dbTempTri,"SELECT DISTINCT match_id
                            FROM vmatch_phenotypes ORDER BY match_id")]
                        )
    matchStructure.matchID = maximum(values(matchIDs)) + 1
    phenotypes = keys(matchIDs)
    println("Generating all possible match phenotypes from...")
    # This loop gives all escape phenotypes (subphenotypes) of each existing
    # phenotype in the dynamics a match ID
    for phenotype in phenotypes
        println("existing phenotype $(phenotype)")
        matchID = matchIDs[phenotype]
        mlength = length(phenotype)
        for k in collect(1:mlength-1)
            for subPhenotype in collect(combinations(phenotype,k))
                subPhenotype = sort(subPhenotype)
                if !in(subPhenotype,keys(matchIDs))
                    matchIDs[subPhenotype] = matchStructure.matchID
                    matchStructure.matchID += 1
                end
            end
        end
    end
    # This loop saves the lengths of all match types as a Dictionary (matchLength:[matchIDs]
    # also saves escape match IDs as a Dictionary (matchID0:[matchIDEscape]).
    # # Note that this does not save the actual phenotype
    for phenotype in keys(matchIDs)
        println("Compiling single escapes for $(phenotype)...")
        if in(length(phenotype),keys(matchStructure.matchtypes))
            push!(matchStructure.matchtypes[length(phenotype)],matchIDs[phenotype])
        else
            matchStructure.matchtypes[length(phenotype)] = Vector{Int64}()
            push!(matchStructure.matchtypes[length(phenotype)],matchIDs[phenotype])
        end
        matchStructure.matches[matchIDs[phenotype]] = phenotype
        matchID = matchIDs[phenotype]

        if length(phenotype) == 1
            matchStructure.escapes[matchIDs[phenotype]] = [0]
        else
            matchStructure.escapes[matchID] = Vector{Int64}()
            for subPhenotype in collect(combinations(phenotype,length(phenotype)-1))
                push!(matchStructure.escapes[matchID],matchIDs[sort(subPhenotype)])
            end
        end
    end
    # This saves all potential match phenotypes to tripartite database
    if !in("potential_vmatch_phenotypes",
                [table_name for (table_name,) in execute(dbTempTri,
                "SELECT name FROM sqlite_schema
                WHERE type='table' ORDER BY name;")])
        matchidlist = Vector{Int64}()
        phenolist = Vector{Int64}()
        for matchid in sort([keys(matchStructure.matches)...])
            append!(phenolist,matchStructure.matches[matchid])
            append!(matchidlist,
                fill(matchid,length(matchStructure.matches[matchid]))
                    )
        end
        DataFrame(match_id = matchidlist, phenotype = phenolist) |>
            SQLite.load!(dbTempTri,
                "potential_vmatch_phenotypes",ifnotexists=true)
        execute(dbTempTri, "CREATE INDEX potential_vmatch_phenotypes_index
                            ON potential_vmatch_phenotypes (match_id)")
    end
    if !in("potential_vmatch_single_escapes",
                [table_name for (table_name,) in execute(dbTempTri,
                "SELECT name FROM sqlite_schema
                WHERE type='table' ORDER BY name;")])
        matchidlist = Vector{Int64}()
        phenolist = Vector{Int64}()
        for matchid in sort([keys(matchStructure.escapes)...])
            append!(phenolist,sort([matchStructure.escapes[matchid]...]))
            append!(matchidlist,
                fill(matchid,length(matchStructure.escapes[matchid]))
                    )
        end
        DataFrame(match_id = matchidlist, escape_match_id = phenolist) |>
            SQLite.load!(dbTempTri,
                "potential_vmatch_single_escapes",ifnotexists=true)
        execute(dbTempTri, "CREATE INDEX potential_vmatch_single_escapes_index
                            ON potential_vmatch_single_escapes (match_id)")
    end
    return matchIDs
end
# This functipn saves all possible subphenotypes (not just escapes) as a Dictionary (matchID0:[matchIDsubphenotypes])
# this is what I call a "match lineage"
function findEscapeLineages!(matchStructure::hierarchy,matchIDs::Dict{Vector{Int64},Int64})
    dbTempSim = matchStructure.dbSim
    dbTempMatch = matchStructure.dbMatch
    dbTempTri = matchStructure.dbTri
    numPhenos = length(keys(matchStructure.matches))
    for matchID in keys(matchStructure.matches)
        phenotype = matchStructure.matches[matchID]
        println("Compiling complete lineage for $(phenotype)...")
        matchStructure.lineages[matchID] = Vector{Int64}()
        mlength = length(phenotype)
        for k in collect(1:mlength-1)
            for subPhenotype in collect(combinations(phenotype,k))
                subPhenotype = sort(subPhenotype)
                push!(matchStructure.lineages[matchID],matchIDs[subPhenotype])
            end
        end
        # numPhenos -= 1
        # println("$(numPhenos) phenotypes left to search")
    end
end

function microbeBiomass!(matchStructure::hierarchy,current::state)
    time = current.time
    dbTempSim = matchStructure.dbSim
    dbTempTri = matchStructure.dbTri
    current.matches = [Int64(matchID) for (matchID,) in
                                execute(dbTempTri,"SELECT DISTINCT match_id
                                FROM vmatches WHERE t = $(time)
                                ORDER BY match_id")]
    matchIDs = current.matches
    union!(current.escapes,matchIDs)
    for matchID in matchIDs
        union!(current.escapes,matchStructure.lineages[matchID])
    end
    tv0 = [t for (t,) in execute(dbTempTri,"SELECT t FROM v0matches_abundances")]
    if in(time,tv0)
        current.matches = vcat(Vector{Int64}([0]),current.matches)
    end
    bstrains = [Int64(strain)
                    for (strain,) in execute(dbTempSim, "SELECT bstrain_id
                        FROM babundance
                        WHERE t = $(time)")]
    spacerIDs = Vector{Int64}()
    for matchID in matchIDs
        union!(spacerIDs,matchStructure.matches[matchID])
    end
    for spacerID in spacerIDs
        current.iBstrains[spacerID] = [Int64(strain)
            for (strain,) in execute(dbTempSim,
                "SELECT DISTINCT bstrain_id
                FROM bspacers
                WHERE spacer_id = $(spacerID)")]
    end
    numEscapes = length(current.escapes)
    for matchID in current.escapes
        # println("Computing abundances susceptible and immune to phenotype $(matchID)")
        if matchID == 0
            B = sum([Int64(abund) for (abund,)
                        in execute(dbTempSim, "SELECT abundance
                        FROM babundance WHERE t = $(time)")])
            current.sBiomass[0] = B
            current.iBiomass[0] = 0
            # numEscapes -= 1
            # println("$(numEscapes) phenotypes left")
            continue
        end
        immuneStrains = Vector{Int64}()
        for spacerID in matchStructure.matches[matchID]
            union!(immuneStrains,current.iBstrains[spacerID])
        end

        current.iBiomass[matchID] =
            sum([Int64(abund) for (abund,)
                in execute(dbTempSim, "SELECT abundance
                    FROM babundance
                    WHERE t = $(time) AND bstrain_id in
                    ($(join(immuneStrains,", ")))")])

        current.sBstrains[matchID] = setdiff(bstrains,immuneStrains)
        if length(current.sBstrains[matchID]) > 0
            current.sBiomass[matchID] =
                sum([Int64(abund) for (abund,)
                    in execute(dbTempSim, "SELECT abundance
                        FROM babundance
                        WHERE t = $(time) AND bstrain_id in
                        ($(join(current.sBstrains[matchID],", ")))")])
        else
            current.sBiomass[matchID] = 0
        end
        # numEscapes -= 1
        # println("$(numEscapes) phenotypes left")
    end
    # sbiomasses = map(x->current.sBiomass[x],sort([current.escapes...]))
    # ibiomasses = map(x->current.iBiomass[x],sort([current.escapes...]))
    # DataFrame(t = fill(time,length(current.escapes)),
    #             match_id = sort([current.escapes...]),
    #             bsusceptible = sbiomasses,
    #             bimmune = ibiomasses
    #             ) |> SQLite.load!(dbTempTri,"potential_vmatches_babundances",
    #                     ifnotexists=true)
end

function computeR0!(extinction::pcomponents,parameters::params)
    time = extinction.time
    phi = parameters.adsorption_rate
    q = parameters.spacer_acquisition_prob
    beta = parameters.viral_burst_size
    d = parameters.viral_decay_rate
    sigdig = extinction.sigdig
    N = sum([abund for (abund,) in
            execute(extinction.dbSim, "SELECT abundance
                FROM babundance WHERE t = $(time)")])
    extinction.birth[0] = beta*phi*(1-q)*N
    extinction.mut[0] = 0
    extinction.death[0] = phi*q*N + d
    a = extinction.birth[0]
    m = extinction.death[0]
    if a > 0
        extinction.proots[0] =
            map(x->round(x,digits=sigdig),[1,m/a])
    else
        extinction.proots[0] = Vector{Float64}([1])
    end
    # println("proots of 0match is $(extinction.proots[0])")
    extinction.prootmin[0] =
        minimum(extinction.proots[0])
    # R0 = (phi*q*N+d)/(beta*phi*(1-q)*N)
    # extinction.R0 = 1/extinction.prootmin[0]
    # return Float64(extinction.R0)
end

# This function assembles all pertinent data from objects and saves into database
function assemble(current::state,extinction::pcomponents)
    time = current.time
    dbTempSim = current.dbSim
    dbTempTri = current.dbTri
    dbOutput = current.dbOutput
    sigdig = extinction.sigdig

    pextinctionRoot = map(x->extinction.prootmin[x],current.matches)
    pextinctionLamb = map(x->extinction.plambert[x],current.matches)
    pextinctionInt = map(x->extinction.pintegrated[x],current.matches)
    birth = map(x->extinction.birth[x],current.matches)
    death = map(x->extinction.death[x],current.matches)
    mut = map(x->extinction.mut[x],current.matches)

    # pextinctionEsc = map(x->extinction.prootmin[x],current.escapes)
    # birthEsc = map(x->extinction.birth[x],current.escapes)
    # deathEsc = map(x->extinction.death[x],current.escapes)
    # mutEsc = map(x->extinction.mut[x],current.escape)

    matchAbunds =  [Int64(abund) for (abund,)
        in execute(dbTempTri, "SELECT vabundance
            FROM vmatches_abundances
            WHERE t = $(time) AND match_id in
            ($(join(current.matches,", "))) ORDER BY match_id")]

    if in(0,current.matches)
        matchAbunds = vcat(sum([Int64(abund) for (abund,)
            in execute(dbTempTri, "SELECT vabundance
                FROM v0matches_abundances
                WHERE t = $(time)")]),matchAbunds)
    end

    V = sum([Int64(abund) for (abund,)
        in execute(dbTempSim, "SELECT abundance
            FROM vabundance
            WHERE t = $(time)")])

    f = 1/V
    return sum(matchAbunds.*pextinctionRoot*f), sum(matchAbunds.*pextinctionLamb*f), sum(matchAbunds.*pextinctionInt*f)
end

function emergence(matchStructure::hierarchy,
                    current::state,parameters::params,
                    init::initial,dt::Float64)
    extinction = computeExtinctionRoots(matchStructure,current,parameters)
    computeLambertRoot!(matchStructure,current,parameters,extinction)
    integrateExtinction!(matchStructure,current,parameters,extinction,
                            init,dt)
    # intExpectation, lambertExpectation, rootExpectation = assemble(current,extinction)
    rootExpectation, lambertExpectation, intExpectation = assemble(current,extinction)
    # return extinction.pintegrated, intExpectation, lambertExpectation, rootExpectation
    return rootExpectation, lambertExpectation, intExpectation
end
# This function runs the whole analysis
function escape(matchStructure, parameters,matchIDs)
    # This initializes the hierarchy of matchtypes that can emerge from escape
    # matchStructure = hierarchy(dbTempTri,dbTempMatch,dbTempSim,dbOutput)
    # # This calls and save the birth/death/mutation parameter rates of the models
    # parameters = params(matchStructure.dbSim)
    # matchIDs = potentialStructure!(matchStructure)
    # findEscapeLineages!(matchStructure,matchIDs)
    pIntExpectation = Vector{Float64}()
    pLambExpectation = Vector{Float64}()
    pRootExpectation = Vector{Float64}()
    times = Vector{Float64}()
    init = initial()
    t0 = 0
    for (t,) in execute(dbTempSim,"SELECT DISTINCT t FROM vabundance")
        if t == 200
            return pRootExpectation, pLambExpectation, pIntExpectation
        end
        push!(times,t)
        println("Computing structure and probabilities at time = $(t)")
        current = state(matchStructure,Float64(t))
        microbeBiomass!(matchStructure,current)
        # pinitial, intExpectation, lambertExpectation, rootExpectation =
        #         emergence(matchStructure,current,parameters)
        rootExpectation, lambertExpectation, intExpectation =
                emergence(matchStructure,current,parameters,
                            init,Float64(t-t0))
        push!(pIntExpectation, intExpectation)
        push!(pLambExpectation, lambertExpectation)
        push!(pRootExpectation, rootExpectation)
        t0 = t
    end
    return pRootExpectation, pLambExpectation, pIntExpectation
end

plot(ones(length(pRootExpectation))-pRootExpectation)

plot([2*ones(length(c2))-c2,ones(length(c3))-c3])

plot([(ones(length(c1))-c1)/(1-minimum(c1)),
(ones(length(c2))-c2)/(1-minimum(c2))])

escape()
