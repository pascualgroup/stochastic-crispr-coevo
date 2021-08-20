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

include("setup.jl")
include("structures.jl")
include("output.jl")
include("util.jl")
include("model.jl")






const P = let
    params_filename = if length(ARGS) == 0
        "parameters.json"
    elseif length(ARGS) == 1
        ARGS[1]
    else
        error("Usage: <path-to>/run.jl [parameters.json]")
    end

    json_str = read(params_filename, String)

    d_str = JSON.parse(json_str)
    d_symb = Dict((Symbol(k), v) for (k, v) in d_str)
    Params(; d_symb...)
end


# Run simulation
function main(P::Params)
    #P = load_parameters_from_json(ARGS[1]) # parameters
    sim = Simulation(P)
    simulate(sim)
end

# Record start time
start_time = now()

main(P)
