#!/usr/bin/env julia

using Logging

# Uncomment this to see all debugging output:
# Logging.global_logger(
#     Logging.SimpleLogger(
#         stderr,
#         Logging.Debug
#     )
# )

# Ensure that UInt64 and UInt32 gets printed in decimal, not hex
Base.show(io::IO, x::UInt64) = Base.print(io, x)
Base.show(io::IO, y::UInt32) = Base.print(io, y)

include("model.jl")

function main()
    params = load_parameters_from_json(ARGS[1])
    sim = Simulation(params)

    simulate(sim)
end

main()
