function extinctionEqn(A,M,beta,Qprod,Q)
    -A*Qprod*Q^beta - M + Q*(A+M)
end


function computeActualRoot!(matchStructure::modHierarchy,
                                current::intState,pstate::laststate,
                                parameters::params)
    dbTempTri = matchStructure.dbTri
    time = current.time
    phi = parameters.adsorption_rate
    q = parameters.spacer_acquisition_prob
    beta = parameters.viral_burst_size
    d = parameters.viral_decay_rate
    mu = parameters.viral_mutation_rate
    N = sum([abund for (abund,) in
            execute(matchStructure.dbSim, "SELECT abundance
                FROM babundance WHERE t = $(time)")])
    extinction = intPcomponents(current)
    sigdig = extinction.sigdig
    extinction.lysis[0] = phi*(1-q)*N
    extinction.death[0] = phi*q*N + d
    a = extinction.lysis[0]
    m = extinction.death[0]
    root = fzero(x->extinctionEqn(a,m,beta,1,x), 10^(-16))
    extinction.pactual[0] = root
    if time == 0.0
        pstate.q[0] = root
        extinction.pintegrated[0] = root
    end

    maxLength = maximum([Int64(mL) for (mL,)
        in execute(dbTempTri, "SELECT DISTINCT match_length
            FROM vmatch_lengths")])

    for matchtype in collect(1:maxLength)
        # println("matchtype is $(matchtype)")
        l = matchtype
        matchIDs = intersect(matchStructure.matchtypes[matchtype],current.escapes)
        for matchID in matchIDs
            # println("phenotype is $(phenotype)")
            S = current.sBiomass[matchID]
            I = current.iBiomass[matchID]
            # println("S = $(S), I = $(I)")
            extinction.lysis[matchID] = phi*(1-q)*S
            a = extinction.lysis[matchID]
            # a = round(extinction.lysis[matchID],digits=sigdig)

            extinction.death[matchID] = phi*I + phi*q*S + d
            m = extinction.death[matchID]
            # c = round(extinction.death[matchID],digits=sigdig)
            # println("a = $(a), m = $(m)")
            if matchtype == 1
                Qprod = extinction.pactual[0]
            else
                Qprod = 1
                for escapeID in matchStructure.escapes[matchID]
                    if !in(escapeID,[keys(pstate.q)...])
                        pstate.q[escapeID] = extinction.plambert[escapeID]
                    end
                    x = extinction.pactual[escapeID]
                    Qprod = Qprod*x^(beta*mu)
                end
            end
            # println("Qprod = $(Qprod)")
            # Qprod = round(Qprod,digits=6)

            root = fzero(x->extinctionEqn(a,m,beta*(1-mu*l),Qprod,x), 10^(-16))
            extinction.pactual[matchID] = root
            if time == 0.0
                pstate.q[matchID] = root
                extinction.pintegrated[matchID] = root
            end
            # println("plambert of match $(matchID)
            #     is $(extinction.plambert[matchID])")
        end
    end
    return extinction
end














function computeActualRoot!(matchStructure::modHierarchy,
                                current::intState,pstate::laststate,
                                parameters::params)
    dbTempTri = matchStructure.dbTri
    time = current.time
    phi = parameters.adsorption_rate
    q = parameters.spacer_acquisition_prob
    beta = parameters.viral_burst_size
    d = parameters.viral_decay_rate
    mu = parameters.viral_mutation_rate
    N = sum([abund for (abund,) in
            execute(matchStructure.dbSim, "SELECT abundance
                FROM babundance WHERE t = $(time)")])
    extinction = intPcomponents(current)
    sigdig = extinction.sigdig
    extinction.lysis[0] = phi*(1-q)*N
    extinction.death[0] = phi*q*N + d
    a = extinction.lysis[0]
    m = extinction.death[0]
    root = fzero(x->extinctionEqn(a,m,beta,1,x), 10^(-16))
    extinction.pactual[0] = root
    if time == 0.0
        pstate.q[0] = root
        if 0 in current.matches
            extinction.pintegrated[0] = root
            # println("this works: $(extinction.pintegrated[0])")
        end
    end


    maxLength = maximum([Int64(mL) for (mL,)
        in execute(dbTempTri, "SELECT DISTINCT match_length
            FROM vmatch_lengths")])

    for matchtype in collect(1:maxLength)
        # println("matchtype is $(matchtype)")
        l = matchtype
        for matchID in matchStructure.matchtypes[matchtype]
            # println("phenotype is $(phenotype)")
            S = current.sBiomass[matchID]
            I = current.iBiomass[matchID]
            # println("S = $(S), I = $(I)")
            extinction.lysis[matchID] = phi*(1-q)*S
            a = extinction.lysis[matchID]
            # a = round(extinction.lysis[matchID],digits=sigdig)

            extinction.death[matchID] = phi*I + phi*q*S + d
            m = extinction.death[matchID]
            # c = round(extinction.death[matchID],digits=sigdig)
            # println("a = $(a), m = $(m)")
            if matchtype == 1
                Qprod = extinction.pactual[0]
            else
                Qprod = 1
                for escapeID in matchStructure.escapes[matchID]
                    x = extinction.pactual[escapeID]
                    Qprod = Qprod*x^(beta*mu)
                end
            end
            # println("Qprod = $(Qprod)")
            # Qprod = round(Qprod,digits=6)

            root = fzero(x->extinctionEqn(a,m,beta*(1-mu*l),Qprod,x), 10^(-16))
            extinction.pactual[matchID] = root
            if time == 0.0
                pstate.q[matchID] = root
                extinction.pintegrated[matchID] = root
            end
            # println("plambert of match $(matchID)
            #     is $(extinction.plambert[matchID])")
        end
    end
    return extinction
end
