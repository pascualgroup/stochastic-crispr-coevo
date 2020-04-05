# coding: utf-8


import pandas as pd
import numpy as np
import scipy as sp
import math
import matplotlib.pyplot as plt
import sys
import os


####################################################################


bact = pd.read_csv('babundance.csv', delimiter=',')
phage = pd.read_csv('vabundance.csv', delimiter=',')


bact['timesOfRecord'] = bact['t']
bact = bact.rename(columns={"t": "time", "bstrain_id": "label", "abundance": "Bdensity"})
tmpdf = bact['timesOfRecord']
bact.drop(labels=['timesOfRecord'], axis=1, inplace = True)
bact.insert(0, 'timesOfRecord', tmpdf)


phage['timesOfRecord'] = phage['t']
phage = phage.rename(columns={"t": "time", "vstrain_id": "label", "abundance": "Pdensity"})
tmpdf = phage['timesOfRecord']
phage.drop(labels=['timesOfRecord'], axis=1, inplace = True)
phage.insert(0, 'timesOfRecord', tmpdf)



bact.to_csv(r'mu1e-7_initialDiffDp1_S10P15_R-'+'seed'+'_Bacteria-abundance.txt',sep=' ', index = False)
phage.to_csv(r'mu1e-7_initialDiffDp1_S10P15_R-'+'seed'+'_Phage-abundance.txt',sep=' ', index = False)







