# coding: utf-8


import pandas as pd
import numpy as np
import scipy as sp
import math
import matplotlib.pyplot as plt
import sys
import os


##############################################


def StochasticOutputBact(bactDF,spacersDF,l):

    bactTMP = bactDF.copy()
    #SP = np.zeros(l)
    
    for i in range(l):
        bactTMP["Spacer"+str(i+1)] = ""
        
    for i in bactTMP["bstrain_id"].unique():
        
        SP = np.zeros(l)
        sp = spacersDF[spacersDF['bstrain_id']==i]['spacer_id'].values
        for j in range(len(sp)):
            SP[j] = sp[j]

        tmp = bactTMP[bactTMP["bstrain_id"]==i]
        tmp.loc[:,'Spacer1':'Spacer'+str(l)] = SP

        bactTMP[bactTMP["bstrain_id"]==i] = tmp[tmp["bstrain_id"]==i]
        
    
    return bactTMP


##############################################


bact = pd.read_csv('babundance.csv', delimiter=',')
spacers = pd.read_csv('bspacers.csv', delimiter=',')


#Lsp = 10
#Lpt = 15

Lsp = int(sys.argv[1])
Lpt = int(sys.argv[2])



BactOutputDF = StochasticOutputBact(bact,spacers,Lsp)

BactOutputDF['timesOfRecord'] = BactOutputDF['t']
BactOutputDF = BactOutputDF.rename(columns={"t": "time", "bstrain_id": "label", "abundance": "density"})
tmpdf = BactOutputDF['timesOfRecord']
BactOutputDF.drop(labels=['timesOfRecord'], axis=1, inplace = True)
BactOutputDF.insert(0, 'timesOfRecord', tmpdf)


BactOutputDF.to_csv(r'mu1e-7_initialDiffDp1_S'+str(Lsp)+'P'+str(Lpt)+'_R-'+'seed'+'_data-bact.csv',sep=' ', index = False)

