#!/usr/bin/env julia

#=
main:
- Julia version: 
- Author: ed
- Date: 2019-10-01
=#


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


# Make sure that we can find StochasticCrispr
push!(LOAD_PATH, abspath(dirname(PROGRAM_FILE)))

using StochasticCrispr

mutable struct AllParameters
    run_parameters::RunParameters
    initialization_parameters::InitializationParameters
    model_parameters::Parameters

    function AllParameters()
        new(RunParameters(), InitializationParameters(), Parameters())
    end
end
JSON2.@format AllParameters noargs

function load_parameters_from_json(filename) :: AllParameters
    str = open(f -> read(f, String), filename)
    JSON2.read(str, AllParameters)
end

function main()
    all_params = load_parameters_from_json(ARGS[1])
    validate(all_params.run_parameters)
    validate(all_params.initialization_parameters)
    validate(all_params.model_parameters)

    run_model(
        all_params.run_parameters,
        all_params.initialization_parameters,
        all_params.model_parameters
    )
end

main()