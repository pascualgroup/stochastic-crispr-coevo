
# coding: utf-8

import pandas as pd
import numpy as np
import scipy as sp
import math
import matplotlib.pyplot as plt
import sys
import os
# import mpld3


import seaborn as sns
from scipy import stats



def StackedPlotDF(data,tp):
    
    stacked_plot = pd.DataFrame()
    ID = "bact"

    if (tp=="bact"):
        ID = "bstrain_id"
    elif (tp=="vir"):
        ID = "vstrain_id"
    
    
    #strains = data['bstrain_id'].unique()
    strains = data[ID].unique()
    tl = len(data["t"].unique())
    
    stacked_plot["t"] = data["t"].unique()

    for s in strains:

        Abs = np.zeros(tl)
        
        #tmp = data[data["bstrain_id"]==s]
        tmp = data[data[ID]==s]
        
        for i in tmp["t"].values:
            abun = tmp[tmp["t"]==i]["abundance"].values
            Abs[int(i)] = abun
        
        stacked_plot[s] = Abs 
        
    
    return stacked_plot




# ######################################

SCRIPT_DIR = os.path.abspath(os.path.dirname(__file__))

bact = pd.read_csv(os.path.join(SCRIPT_DIR, 'babundance.csv'), delimiter=',')
phage = pd.read_csv(os.path.join(SCRIPT_DIR, 'vabundance.csv'), delimiter=',')
#data = pd.read_csv(path+'time-series-data.txt', delimiter=' ')



stacked_plot = StackedPlotDF(bact,"bact")
Vstacked_plot = StackedPlotDF(phage,"vir")

stacked_plot.to_csv(r'Bacteria_Abundance-StackedPlot-DataFrame.csv')
Vstacked_plot.to_csv(r'Virus_Abundance-StackedPlot-DataFrame.csv')

