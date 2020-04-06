# Stochastic CRISPR model implementation

## Running a single simulation

There are three categories of parameters, organized into separate structs in the code:

1. Run parameters (`RunParameters` struct): start and stop time of the simulation, and random seed (`rng_seed`)
2. Initialization parameters (`InitializationParameters` struct): parameters governing the initial population
3. Model parameters (`Parameters` struct): parameters governing the dynamical process

In the simplest case, a single simulation can be run by constructing a new script to run the model with specific parameters.
An example of this is in `example-simple-run.jl`.

More likely, you'll want to do many runs with different parameters.
Rather than creating a new Julia source file for each run, you can construct a parameters file in JSON format for each run.

You can run a single simulation using `run.jl`, with a single argument giving the location of the parameters file, like this:

```sh
julia /path/to/stochastic-crispr/julia/run.jl parameters.json
```

An example parameters file is in `example-parameters.json`.

## Running a parameter sweep

To sweep across many parameter values, you'll want to write a script to generate a working directory for each simulation, containing a parameters file for that simulation and perhaps also a script that submits the simulation to the cluster.

A full example of this is provided in the `example-sweep` directory; see `README.md` in that directory for documentation.
