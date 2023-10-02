function computeLambertRoot!(matchStructure::hierarchy,
                                current::state,parameters::params,
                                extinction::pcomponents)
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
    sigdig = extinction.sigdig
    extinction.lysis[0] = phi*(1-q)*N
    a = extinction.birth[0]
    m = extinction.death[0]
    lambArg = exp((q-1)*N*beta*phi/(d+N*phi))*(q-1)*N*beta*phi/(d+N*phi)
    extinction.plambert[0] =
        (d+q*N*phi)/(d+N*phi) - lambertw(lambArg)/beta

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
        # println("Computing Lambert probability for match lengths of $(matchtype)...")
        matchIDs = intersect(matchStructure.matchtypes[matchtype],current.escapes)
        for matchID in matchIDs
            # println("phenotype is $(phenotype)")
            S = current.sBiomass[matchID]
            I = current.iBiomass[matchID]

            extinction.lysis[matchID] = phi*(1-q)*S
            a = extinction.lysis[matchID]
            # a = round(extinction.lysis[matchID],digits=sigdig)

            extinction.death[matchID] = phi*I + phi*q*S + d
            c = extinction.death[matchID]
            # c = round(extinction.death[matchID],digits=sigdig)

            Qprod = 1
            for escapeID in matchStructure.escapes[matchID]
                x = extinction.plambert[escapeID]
                Qprod = Qprod*exp(-1*beta*mu*(1-x))
            end
            # Qprod = round(Qprod,digits=6)

            lambArg = exp((q-1)*S*beta*phi*(mu*l-1)/(d+N*phi))*
                        Qprod*S*beta*phi*(q-1)*(1-mu*l)*phi/(d+N*phi)
            extinction.plambert[matchID] =
                (d+(q*S+I)*phi)/(d+N*phi) + lambertw(lambArg)/(beta*(mu*l-1))
            # println("plambert of match $(matchID)
            #     is $(extinction.plambert[matchID])")
        end
    end
    return extinction
end

function extinctionEqn(A,M,beta,Qprod,Q)
    -A*Qprod*Q^beta - M + Q*(A+M)
end


function computeActualRoot!(matchStructure::hierarchy,
                                current::state,parameters::params,
                                extinction::pcomponents)
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
    sigdig = extinction.sigdig
    a = extinction.lysis[0]
    m = extinction.death[0]
    extinction.pactual[0] = fzero(x->extinctionEqn(a,m,beta,1,x), 10^(-16))

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
        # println("Computing Lambert probability for match lengths of $(matchtype)...")
        matchIDs = intersect(matchStructure.matchtypes[matchtype],current.escapes)
        for matchID in matchIDs
            # println("phenotype is $(phenotype)")
            S = current.sBiomass[matchID]
            I = current.iBiomass[matchID]

            extinction.lysis[matchID] = phi*(1-q)*S
            a = extinction.lysis[matchID]
            # a = round(extinction.lysis[matchID],digits=sigdig)

            extinction.death[matchID] = phi*I + phi*q*S + d
            m = extinction.death[matchID]
            # c = round(extinction.death[matchID],digits=sigdig)

            if l == 1
                Qprod = extinction.pactual[0]
            else
                Qprod = 1
                for escapeID in matchStructure.escapes[matchID]
                    x = extinction.pactual[escapeID]
                    Qprod = Qprod*x^(beta*mu)
                end
            end
            # Qprod = round(Qprod,digits=6)

            extinction.pactual[matchID] =
                fzero(x->extinctionEqn(a,m,beta*(1-mu*l),Qprod,x), 10^(-16))
            # println("plambert of match $(matchID)
            #     is $(extinction.plambert[matchID])")
        end
    end
    return extinction
end
