
#!/bin/sh

module load julia
 
max=$1  # number of simulations
 
for i in $(seq 1 $max)
do
        cp -r master sim_$i
        cd sim_$i
        sbatch StochasticCRISPR.sbatch
        cd ..
done

echo "OutputConvertions Submited!!"

