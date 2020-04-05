# coding: utf-8


import pandas as pd
import numpy as np
import scipy as sp
import math
import matplotlib.pyplot as plt
import sys
import os


##############################################

def StochasticOutput(phageDF,protospacersDF):

    phageTMP = phageDF.copy()

    l = protospacersDF[protospacersDF['vstrain_id']==1]['spacer_id'].values
    for j in range(len(l)):
        phageTMP["Protospacer"+str(j+1)] = ""

    for i in phageTMP["vstrain_id"].unique():

        Pt = protospacersDF[protospacersDF['vstrain_id']==i]['spacer_id'].values

        tmp = phageTMP[phageTMP["vstrain_id"]==i]
        tmp.loc[:,'Protospacer1':'Protospacer'+str(len(l))] = Pt

        phageTMP[phageTMP["vstrain_id"]==i] = tmp[tmp["vstrain_id"]==i]
    
    return phageTMP


##############################################


phage = pd.read_csv('vabundance.csv', delimiter=',')
protospacers = pd.read_csv('vpspacers.csv', delimiter=',')


#Lsp = 10
Lsp = int(sys.argv[1])

Lpt = len(protospacers[protospacers['vstrain_id']==1]['spacer_id'].values)


VirOutputDF = StochasticOutput(phage,protospacers)

VirOutputDF['timesOfRecord'] = VirOutputDF['t']
VirOutputDF = VirOutputDF.rename(columns={"t": "time", "vstrain_id": "label", "abundance": "density"})
tmpdf = VirOutputDF['timesOfRecord']
VirOutputDF.drop(labels=['timesOfRecord'], axis=1, inplace = True)
VirOutputDF.insert(0, 'timesOfRecord', tmpdf)


VirOutputDF.to_csv(r'mu1e-7_initialDiffDp1_S'+str(Lsp)+'P'+str(Lpt)+'_R-'+'seed'+'_data-phage.csv', sep=' ', index = False)

