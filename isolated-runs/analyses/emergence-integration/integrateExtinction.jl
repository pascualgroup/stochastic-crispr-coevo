function dQ(A,M,beta,Qprod,Q0)
    return -A*Qprod*Q0^beta - M + Q0*(A+M)
end


function integrateExtinction!(matchStructure::modHierarchy,
                                current::intState,parameters::params,
                                extinction::intPcomponents,
                                pstate::laststate,dt::Float64)
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
    # println("pstate is $(pstate.q)")
    if time == 0.0
        return
    end
    a = extinction.lysis[0]
    m = extinction.death[0]
    q0 = pstate.q[0] + dQ(a,m,beta,1,pstate.q[0])*dt
    qInt = pstate.q[0] +
        dt/2*(dQ(a,m,beta,1,pstate.q[0]) + dQ(a,m,beta,1,q0))
    pstate.q[0] = qInt

    if 0 in current.matches
        extinction.pintegrated[0] = qInt
    end

    if current.matches == [0]
        maxLength = 0
    else
        maxLength = maximum([Int64(mL) for (mL,)
            in execute(dbTempTri, "SELECT match_length
                FROM vmatch_lengths WHERE match_id in
                ($(join(current.matches,", ")))")])
    end

    for matchtype in collect(1:maxLength)
        # println("matchtype is $(matchtype)")
        l = matchtype
        # println("Computing Integrated probability for match lengths of $(matchtype)...")
        matchIDs = intersect(matchStructure.matchtypes[matchtype],current.escapes)
        for matchID in matchIDs

            S = current.sBiomass[matchID]
            I = current.iBiomass[matchID]

            extinction.lysis[matchID] = phi*(1-q)*S
            a = extinction.lysis[matchID]
            # a = round(extinction.lysis[matchID],digits=sigdig)

            extinction.death[matchID] = phi*I + phi*q*S + d
            m = extinction.death[matchID]
            # c = round(extinction.death[matchID],digits=sigdig)

            if l == 1
                Qprod = pstate.q[0]^(beta*mu)
            else
                Qprod = 1
                for escapeID in matchStructure.escapes[matchID]
                    if !in(matchStructure.escapes[escapeID],[keys(pstate.q)...])

                    Qprod = Qprod*pstate.q[escapeID]^(beta*mu)
                end
            end

            q0 = pstate.q[matchID] + dQ(a,m,beta*(1-mu*l),Qprod,pstate.q[matchID])*dt
            qInt = pstate.q[matchID] +
                dt/2*(dQ(a,m,beta*(1-mu*l),Qprod,pstate.q[matchID]) + dQ(a,m,beta*(1-mu*l),Qprod,q0))

            pstate.q[matchID] = qInt
            extinction.pintegrated[matchID] = qInt
        end
    end
    return
end





function integrateExtinction!(matchStructure::modHierarchy,
                                current::intState,parameters::params,
                                extinction::intPcomponents,
                                pstate::probstate,dt::Float64)
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
    # println("pstate is $(pstate.q)")
    if time == 0.0
        return
    end
    a = extinction.lysis[0]
    m = extinction.death[0]
    q0 = pstate.q[0] + dQ(a,m,beta,1,pstate.q[0])*dt
    qInt = pstate.q[0] +
        dt/2*(dQ(a,m,beta,1,pstate.q[0]) + dQ(a,m,beta,1,q0))
    pstate.q[0] = qInt

    if 0 in current.matches
        extinction.pintegrated[0] = qInt
    end

    if current.matches == [0]
        maxLength = 0
    else
        maxLength = maximum([Int64(mL) for (mL,)
            in execute(dbTempTri, "SELECT match_length
                FROM vmatch_lengths WHERE match_id in
                ($(join(current.matches,", ")))")])
    end

    for matchtype in collect(1:maxLength)
        # println("matchtype is $(matchtype)")
        l = matchtype
        # println("Computing Integrated probability for match lengths of $(matchtype)...")
        for matchID in matchStructure.matchtypes[matchtype]

            S = current.sBiomass[matchID]
            I = current.iBiomass[matchID]

            extinction.lysis[matchID] = phi*(1-q)*S
            a = extinction.lysis[matchID]
            # a = round(extinction.lysis[matchID],digits=sigdig)

            extinction.death[matchID] = phi*I + phi*q*S + d
            m = extinction.death[matchID]
            # c = round(extinction.death[matchID],digits=sigdig)

            if l == 1
                Qprod = pstate.q[0]^(beta*mu)
            else
                Qprod = 1
                for escapeID in matchStructure.escapes[matchID]
                    Qprod = Qprod*pstate.q[escapeID]^(beta*mu)
                end
            end

            q0 = pstate.q[matchID] + dQ(a,m,beta*(1-mu*l),Qprod,pstate.q[matchID])*dt
            qInt = pstate.q[matchID] +
                dt/2*(dQ(a,m,beta*(1-mu*l),Qprod,pstate.q[matchID]) + dQ(a,m,beta*(1-mu*l),Qprod,q0))

            pstate.q[matchID] = qInt
            extinction.pintegrated[matchID] = qInt
        end
    end
    return
end
