# Example parameter sweep

The script `example-sweep.py` shows how to generate working directories for many runs by sweeping across multiple sets of parameter values.

The script sweeps over two parameters, `n_protospacers` and `u_n_protospacers_max`, and generates a directory hierarchy in `runs` that looks like this:

```
nps=10-u=8/
    0/
        parameters.json
        run.sbatch
    1/
    ...
    19/
```

where the directory `nps=10-u=8/0` corresponds to replicate `0`, with a unique random seed, for parameter combination `n_protospacers = 10` and `u_n_protospacers_max = 8`.
Each run directory contains a file `parameters.json`, which is used to run the model, and a SLURM job script `run.sbatch`, which looks like this:

```sh
#!/bin/bash

#SBATCH --job-name=SC-nps=10-u=8-0
#SBATCH --output=StochasticCRISPR.out
#SBATCH --error=StochasticCRISPR.err
#SBATCH --partition=broadwl
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1

julia /Users/ed/uchicago/crispr/crispr-model/julia/run.jl parameters.json
```

These job scripts can be submitted individually, of course, or run directly using `bash`, e.g., on your local machine.

But you'll also notice a file `submit.sh`, which simply consists of a whole bunch of SLURM submission commands:

```sh
#!/bin/bash

module load julia

cd /Users/ed/uchicago/crispr/crispr-model/example-sweep/runs/nps=10-u=8/0
sbatch run.sbatch
cd /Users/ed/uchicago/crispr/crispr-model/example-sweep/runs/nps=10-u=8/1
sbatch run.sbatch
cd /Users/ed/uchicago/crispr/crispr-model/example-sweep/runs/nps=10-u=8/2
sbatch run.sbatch
...
```
