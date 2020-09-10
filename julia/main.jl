#!/usr/bin/env julia

import Pkg

Pkg.add("JSON2")
Pkg.add("StatsBase")
Pkg.add("JSON2")
Pkg.add("Distributions")


using Logging
using Random
using Dates
using JSON2

# Uncomment this to see all debugging output:
# Logging.global_logger(
#     Logging.SimpleLogger(
#         stderr,
#         Logging.Debug
#     )
# )

# Ensure that UInt64 gets printed in decimal, not hex
Base.show(io::IO, x::UInt64) = Base.print(io, x)

include("model.jl")

function main()
    params = load_parameters_from_json(ARGS[1])
    sim = Simulation(params)

    simulate(sim)
end

main()
