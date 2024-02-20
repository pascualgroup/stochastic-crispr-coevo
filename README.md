# Stochastic CRISPR model implementation

## Running a single simulation

You can run a single simulation using `main-sweep.jl`, with a single argument giving the location of the parameters file in JSON format, like this:

```sh
julia /simulation/src/main-sweep.jl parameters.json
```
An example parameters file is in `simulation/src/parameters.json`.

## Running a parameter sweep

To sweep across many parameter values, you'll want to write a script to generate a working directory for each simulation, containing a parameters file for that simulation and perhaps also a script that submits the simulation to the cluster.

An example of this script is `/simulation/generate-sweep.jl`, which calls a JSON file designating the parameter intervals (`sweep-parameters.json`). This script is written for a cluster operating with SLURM. Number of processes, cores and jobs can be customized. Run
```sh
julia /simulation/generate-sweep.jl
```


Working directories for each simulation and sbatch script will be generated, along with a shell script that submits jobs.

## Model structure and simulation algorithm

The model is a stochastic, continuous-time, discrete-event simulation (Gillespie-style) that steps from event to event, with a random amount of time, real-valued, between events.

## Model state

Model state is represented by the `State` struct, which contains two fields, `bstrains`, of type `BStrains`, and `vstrains`, of type `VStrains`.

The `BStrains` and `VStrains` structs are very similar.
Conceptually, they represent a collection of bacterial/viral strains, and although it would make conceptual sense to just have `bstrains::Vector{BStrain}`, for efficiency of memory layout they are implemented using parallel arrays.

## Event loop

Each step follows this algorithm:

1. The time to the next event is computed by drawing a waiting time from a Poisson process whose rate is the total rate of all possible events in the system.
2. The type of event is chosen proportional to the total rate of all events of that type (e.g., all possible bacterial growth events).
3. For each event type, the specific event is chosen by sampling proportional to the rates of the individual events. For example, for bacterial growth events, since bacterial cells grow at the same rate, the growth rate of particular strain is proportional to the abundance of that strain, and therefore a strain is chosen proportional to its abundance.
4. The event is executed.

There are four events:

1. Microbial growth: `const MICROBIAL_GROWTH = 1`
2. Microbial death: `const MICROBIAL_DEATH = 2`
3. Viral decay: `const VIRAL_DECAY = 3`
4. Contact: `const CONTACT = 4`
4. Bacterial immigration: `const MICROBIAL_IMMIGRATION = 5`

### Contact Events

1. Randomly select a bacterial strain and a viral strain proportional to population size
2. If immune, infect with probability p_crispr_failure_prob
3. If not immune, defend successfully & acquire spacer with prob q_spacer_acquisition_prob

Infect:

1. Draw the number of mutations for each of `beta` newly created virus particles; the number is `Binomial(n_pspacers, mu)`.
2. Unmutated virus particles just increase the abundance of the existing viral strain.
3. Mutated virus particles acquire newly created protospacers at randomly selected loci.

## Generating Figures

Please download analysis data (DOI 10.6084/m9.figshare.25254565), and put contents of the folder into a directory named `generated-data`. This directory should be in the same containing directory as `isolated-runs`, `simulation`, etc.

### Figure 1
Run the following
```sh
python /isolated-runs/tree-figure.py
```

### Figure 2
Run the following
```sh
python /isolated-runs/analyses-figures.py host-strain-collapses
```

### Protospacer Shannon Diversity in Figure 3
Run the following
```sh
python /isolated-runs/analyses-figures.py protospacer-div
```

### Tripartite graphs in Figure 3
Run the following
```sh
python /isolated-runs/triartitegraph.py
```

### Tripartite graphs in Figure 3
Run the following
```sh
python /isolated-runs/triartitegraph.py
```

### Figures 4 and S5
Run the following
```sh
python /isolated-runs/analyses-figures.py R
```

### Figures 3, S7 and S8
Run the following
```sh
python /isolated-runs/make-phase-diagrams.py
```

### Figure S1
Run the following
```sh
python /isolated-runs/analyses-figures.py diversity
```

### Figure S2
Run the following
```sh
python /isolated-runs/analyses-figures.py viralheatmap
```

### Figure S3
Run the following
```sh
python /isolated-runs/analyses-figures.py escapeheatmap
```

### Figure S4
Run the following
```sh
python /isolated-runs/analyses-figures.py escapeheatmap
```
