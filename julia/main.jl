#=
main:
- Julia version: 
- Author: ed
- Date: 2019-10-01
=#

push!(LOAD_PATH, ".")

using StochasticCrispr

using Dates


function main()
    wallclock_start = now()
    print("starting at ", wallclock_start, "\n")

    params = Parameters()
    params.beta = 1.0
    validate(params)

    state = State()
    validate(state)

    sim = Simulator(params, 0.0, state)

    println("sim.t: ", sim.t)
    t_final = 10000000.0
    while sim.t < t_final
        do_next_event!(sim, t_final)
        #println("sim.state.t: ", sim.t)
    end
    println("sim.t: ", sim.t)
    println(sim.state)

    wallclock_end = now()
    println("ending at ", wallclock_end)
    elapsed_time = Dates.value(wallclock_end - wallclock_start) / 1000.0
    println("elapsed time: ", elapsed_time, " seconds")
 end

main()
