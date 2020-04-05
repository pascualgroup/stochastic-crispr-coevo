
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


stacked_plot = pd.read_csv('Bacteria_Abundance-StackedPlot-DataFrame.csv', delimiter=',')
Vstacked_plot = pd.read_csv('Virus_Abundance-StackedPlot-DataFrame.csv', delimiter=',')


#plt.figure(figsize=(10, 10), dpi= 80)
stacked_plot.plot.area(stacked=True, legend=False, linewidth=0);
#plt.show()
plt.tight_layout()  
plt.savefig("Bacteria-Abundance_stacked_plot.png",dpi=500)
plt.close()   


#plt.figure(figsize=(10, 10), dpi= 80)
Vstacked_plot.plot.area(stacked=True, legend=False, linewidth=0);
#plt.show()
plt.tight_layout()  
plt.savefig("Virus-Abundance_stacked_plot.png",dpi=500)
plt.close()   
