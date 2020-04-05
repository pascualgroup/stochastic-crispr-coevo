
#!/bin/sh


module load python

#max=$1
max=$(ls -d sim_* | wc -l)

for i in $(seq 1 $max)
do
        cd sim_$i
        sbatch StochasticCRISPRFileForStackedPlot.sbatch
        sbatch StochasticCRISPRConvertionAbundances.sbatch
        sbatch StochasticCRISPRConvertionOutputBacteria.sbatch
        sbatch StochasticCRISPRConvertionOutputVirus.sbatch        
        cd ..
done

echo "OutputConvertions Submited!!"
 


