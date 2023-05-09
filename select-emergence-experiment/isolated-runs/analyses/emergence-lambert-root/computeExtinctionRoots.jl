function computeExtinctionRoots(matchStructure::hierarchy,
                                current::state,parameters::params)
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

    extinction = pcomponents(current)
    sigdig = extinction.sigdig
    computeR0!(extinction,parameters)

    if current.matches == [0]
        # println("this works")
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
        # println("Computing extinction probability for match lengths of $(matchtype)...")
        matchIDs = intersect(matchStructure.matchtypes[matchtype],current.escapes)
        for matchID in matchIDs
            # println("phenotype is $(phenotype)")
            if length(matchStructure.matches[matchID]) == 1
                S = current.sBiomass[matchID]
                I = current.iBiomass[matchID]

                extinction.birth[matchID] = beta*phi*(1-q)*(1-l*mu)*S
                a = extinction.birth[matchID]
                # a = round(extinction.birth[matchID],digits=sigdig)

                extinction.mut[matchID] = beta*phi*(1-q)*mu*S
                b = extinction.mut[matchID]
                # extinction.mut[matchID] = beta*phi*(1-q)*l*mu*S
                # b = extinction.mut[matchID]
                # b = round(extinction.mut[matchID],digits=sigdig)

                extinction.death[matchID] = phi*I + phi*q*S + d
                c = extinction.death[matchID]
                # c = round(extinction.death[matchID],digits=sigdig)

                Qsum = extinction.prootmin[0]
                # Qsum = round(Qsum,digits=6)

                if a !=  0.0
                    if -4*a*c+(a+b+c-b*Qsum)^2 < 0.0
                        extinction.proots[matchID] =
                        map(x->round(x,digits=sigdig), [(a+b+c-b*Qsum)/(2*a)])
                    # println("proots of match $(matchID)
                    #     is $(extinction.proots[matchID])")
                        @assert imag(extinction.proots[matchID][1]) == 0 "root is imaginary!"
                    else
                        extinction.proots[matchID] =
                            map(x->round(x,digits=sigdig),
                                [(a+b+c-b*Qsum-sqrt(-4*a*c+(a+b+c-b*Qsum)^2))/(2*a),
                                (a+b+c-b*Qsum+sqrt(-4*a*c+(a+b+c-b*Qsum)^2))/(2*a)])
                        # println("proots of match $(matchID)
                        #     is $(extinction.proots[matchID])")
                        @assert imag(extinction.proots[matchID][1]) == 0 "root is imaginary!"
                        @assert imag(extinction.proots[matchID][2]) == 0 "root is imaginary!"
                    end
                else
                    extinction.proots[matchID] =
                        map(x->round(x,digits=sigdig),[c/(b+c-b*Qsum)])
                    # println("proots of match $(matchID)
                    #     is $(extinction.proots[matchID])")
                    @assert imag(extinction.proots[matchID][1]) == 0 "root is imaginary!"
                end

                extinction.prootmin[matchID] =
                    minimum(extinction.proots[matchID])
            else
                S = current.sBiomass[matchID]
                I = current.iBiomass[matchID]

                extinction.birth[matchID] = beta*phi*(1-q)*(1-l*mu)*S
                a = extinction.birth[matchID]

                extinction.mut[matchID] = beta*phi*(1-q)*mu*S
                b = extinction.mut[matchID]

                extinction.death[matchID] = phi*I + phi*q*S + d
                c = extinction.death[matchID]

                Qsum = 0
                for escapeID in matchStructure.escapes[matchID]
                    Qsum += extinction.prootmin[escapeID]
                end

                if a != 0.0
                    if -4*a*c+(a+b+c-b*Qsum)^2 < 0.0
                        extinction.proots[matchID] =
                        map(x->round(x,digits=sigdig), [(a+b+c-b*Qsum)/(2*a)])
                    # println("proots of match $(matchID)
                    #     is $(extinction.proots[matchID])")
                        @assert imag(extinction.proots[matchID][1]) == 0 "root is imaginary!"
                    else
                        extinction.proots[matchID] =
                            map(x->round(x,digits=sigdig),
                                [(a+b+c-b*Qsum-sqrt(-4*a*c+(a+b+c-b*Qsum)^2))/(2*a),
                                (a+b+c-b*Qsum+sqrt(-4*a*c+(a+b+c-b*Qsum)^2))/(2*a)])
                        # println("proots of match $(matchID)
                        #     is $(extinction.proots[matchID])")
                        @assert imag(extinction.proots[matchID][1]) == 0 "root is imaginary!"
                        @assert imag(extinction.proots[matchID][2]) == 0 "root is imaginary!"
                    end
                else
                    extinction.proots[matchID] =
                        map(x->round(x,digits=sigdig),[c/(b+c-b*Qsum)])
                    # println("proots of match $(matchID)
                    #     is $(extinction.proots[matchID])")
                    @assert imag(extinction.proots[matchID][1]) == 0 "root is imaginary!"
                end



                # THIS IS BAD PRACTICE but I know what the solutions look like analytically
                # I.E. minimum roots that are greater than one, are so by a very small amount
                # due to floating point errors
                if minimum(extinction.proots[matchID]) > 1
                    extinction.prootmin[matchID] = 1
                else
                    extinction.prootmin[matchID] =
                        minimum(extinction.proots[matchID])
                end
            end
        end
    end
    return extinction
end
