#!/usr/bin/env python3

import pandas as pd
import numpy as np
import scipy as sp
import math
import matplotlib.pyplot as plt
import sys
import os
import seaborn as sns
from scipy import stats
import sqlite3
import matplotlib.ticker as ticker
import matplotlib.cm as cm
from matplotlib.colors import Normalize
import matplotlib.colors as mc
import colorsys
from matplotlib.collections import LineCollection

run_id = sys.argv[1]

resolve = 500
imgType = "pdf" #or "png"

SCRIPT_PATH = os.path.abspath(os.path.dirname(__file__)) # cluster

DBSIM_PATH = os.path.join(SCRIPT_PATH,'..','..','..','simulation','sweep_db_gathered.sqlite') # cluster
# DBSIM_PATH = os.path.join('/Volumes','Yadgah','run_id1455_combo73_replicate15.sqlite')
# DBSIM_PATH = os.path.join('/Volumes','Yadgah','crispr-sweep-7-2-2022/isolates/runID3297-c66-r47/runID3297-c66-r47.sqlite')
DBPROB_PATH = os.path.join('/Volumes','Yadgah','probability-of-emergence_output.sqlite')

conSim = sqlite3.connect(DBSIM_PATH)
curSim = conSim.cursor()
# Designating plot path from simulation data
ID = curSim.execute('SELECT combo_id,replicate FROM runs WHERE run_id = {}'.format(run_id)).fetchall()
combo_id = ID[0][0]
replicate = ID[0][1]
RUN_DIR = os.path.join('runID{0}-c{1}-r{2}'.format(run_id,combo_id,replicate))

conProb = sqlite3.connect(DBPROB_PATH)
curProb = conProb.cursor()

print('SQLite Query: microbial abundance time series data')
microbe_stacked = pd.read_sql_query("SELECT t,bstrain_id,abundance FROM babundance WHERE run_id = {}".format(run_id), conSim)
microbe_stacked = microbe_stacked.pivot(index='t',columns='bstrain_id',values='abundance')
print('SQLite Query: viral abundance time series data')
virus_stacked = pd.read_sql_query("SELECT t,vstrain_id,abundance FROM vabundance WHERE run_id = {}".format(run_id), conSim)
virus_stacked = virus_stacked.pivot(index='t',columns='vstrain_id',values='abundance')

virusTotal = pd.read_sql_query("SELECT t,viral_abundance FROM summary WHERE run_id = {}".format(run_id), conSim)

if len(microbe_stacked.index) > len(virus_stacked.index):
    microbe_stacked.drop(index=list(set(microbe_stacked.index)-set(virus_stacked.index)),inplace=True)

pal = sns.color_palette("tab20b")

Ravg = pd.read_sql_query("SELECT t,R \
    FROM Ravg", conProb)

fig, ax = plt.subplots(3,sharex=True)
axes = [ax[0], ax[1], ax[2]]
microbe_stacked.plot.area(ax = axes[0],stacked=True,legend=False, linewidth=0,color=pal)
axes[1].plot(Ravg.t,Ravg.R)
axes[1].axhline(y=1, color='r', linestyle='-')
axes[2].plot(virusTotal.t,virusTotal.viral_abundance,color='tab:orange')
fig.tight_layout()
fig.savefig(os.path.join(SCRIPT_PATH,'..','..','isolates',RUN_DIR,'microbe-virus-stacked-abundances.{}'.format(imgType)),dpi=resolve)
fig.close()


pEmerge = pd.read_sql_query("SELECT t,p \
    FROM pEmergenceSum", conProb)
pal = sns.color_palette("tab20b")
fig, ax = plt.subplots(2,sharex=True)
axes = [ax[0], ax[1]]
microbe_stacked.plot.area(ax = axes[0],stacked=True,legend=False, linewidth=0,color=pal)
axes[1].plot(pEmerge.t,pEmerge.p)
fig.tight_layout()
fig.savefig(os.path.join(SCRIPT_PATH,'..','..','isolates',RUN_DIR,'microbe-virus-stacked-abundances.{}'.format(imgType)),dpi=resolve)
fig.close()


Rweighted = pd.read_sql_query("SELECT t,time_specific_match_id,R \
    FROM Rweighted ORDER BY t,time_specific_match_id", conProb)
pExt = pd.read_sql_query("SELECT t,match_id,p \
    FROM pExtinction", conProb)
R0 = pd.read_sql_query("SELECT t,R \
    FROM R0", conProb)
matchIDs = pd.read_sql_query("SELECT t,time_specific_match_id,phenotype \
    FROM match_phenotypes", conProb)

matchPhenos = []
phenoIDs = []
match_id = 1
newIDs = []
for t in np.unique(Rweighted.t):
    ids = np.unique(Rweighted[Rweighted.t==t].time_specific_match_id)
    ids.sort()
    for id in ids:
        if id == 0:
            newIDs.append(0)
            continue

        matchPheno = matchIDs[(matchIDs.t==t) & (matchIDs.time_specific_match_id == id)].phenotype.values
        matchPheno.sort()
        matchPheno = list(matchPheno)
        if matchPheno in matchPhenos:
            phenoID = matchPhenos.index(matchPheno) + 1
            newIDs.append(phenoID)
            continue
        else:
            matchPhenos.append(matchPheno)
            phenoIDs.append(match_id)
            newIDs.append(match_id)
            match_id += 1

Rweighted['pheno_id'] = newIDs



pExt['pheno_id'] = newIDs





print('Complete!')
