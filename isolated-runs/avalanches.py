import pandas as pd
import numpy as np
import scipy as sp
import math
import matplotlib.pyplot as plt
from matplotlib import ticker
import sys
import os
import seaborn as sns
from scipy import stats
import sqlite3
from random import sample
import matplotlib.cm as cm
from matplotlib.colors import Normalize


DBAVALANCHE_PATH = os.path.join('/Volumes/Yadgah/avalanches/avalanches.sqlite')
DBSIM_PATH = os.path.join('/Volumes/Yadgah/avalanches/sweep_db.sqlite')



conSim = sqlite3.connect(DBSIM_PATH)
curSim = conSim.cursor()
conAval = sqlite3.connect(DBAVALANCHE_PATH)
curAval = conAval.cursor()



epiTimeColumn = 'normalized_distance_to_wall_end_max'
viralMutProb = 1.68e-6
# spacerAcqProb = 2.5*10**(-5)
# comboSpace = pd.read_sql_query(
#     "SELECT combo_id \
# FROM param_combos WHERE viral_mutation_rate = {0} AND spacer_acquisition_prob = {1} \
#         ORDER BY combo_id"
#     .format(viralMutProb, spacerAcqProb),conSim)
comboSpace = pd.read_sql_query(
    "SELECT combo_id, viral_mutation_rate, spacer_acquisition_prob \
        FROM param_combos WHERE viral_mutation_rate = {0} \
        ORDER BY spacer_acquisition_prob"
    .format(viralMutProb),conSim)
curCombos = []
conCombos = []
for cID in comboSpace['combo_id']:
    DBCOMBO_PATH = os.path.join('/Volumes/Yadgah/avalanches/simulations/comboID-{0}.sqlite'.format(cID))
    con = sqlite3.connect(DBCOMBO_PATH)
    conCombos.append(con)
    curCombos.append(con.cursor())

avalancheSizeAvg = pd.DataFrame({'index' : [],'size':[], 'combo_id':[]})
for cIdx in range(0,len(comboSpace['combo_id'])):
    df = pd.read_sql_query(
        "SELECT * FROM avalanche_averages",conCombos[cIdx])
    df = pd.concat([df, pd.DataFrame({'combo_id':len(df)*list([comboSpace['combo_id'][cIdx]])})],\
                    axis=1)
    avalancheSizeAvg = pd.concat([avalancheSizeAvg,df], ignore_index=True)
    
runs = pd.read_sql_query(
    "SELECT combo_id, run_id \
FROM runs WHERE combo_id in ({0}) \
        ORDER BY combo_id"
    .format(', '.join(map(str, comboSpace['combo_id']))),conSim)
avalanches = pd.read_sql_query(
    "SELECT * \
FROM avalanche_index_distances WHERE run_id in ({0}) \
        ORDER BY run_id"
    .format(', '.join(map(str, runs['run_id']))),conAval)\
    .merge(runs,on=['run_id'])
avalanchesAll = pd.read_sql_query(
    "SELECT * \
FROM avalanche_all_index_distances WHERE run_id in ({0}) \
        ORDER BY run_id"
    .format(', '.join(map(str, runs['run_id']))),conAval)\
    .merge(runs,on=['run_id'])


# avalanchesTrunc = avalanches[avalanches.wall_number == wallNum]
avalanchesMeanAll = avalanchesAll.groupby(['combo_id','wall_number', 'avalanche_index'])\
    .agg(distance=(epiTimeColumn, 'mean'),
            std=(epiTimeColumn, 'std'),
            n=(epiTimeColumn, 'size'))\
            .reset_index()
avalanchesMean = avalanches.groupby(['combo_id','wall_number', 'avalanche_index'])\
    .agg(distance=(epiTimeColumn, 'mean'),
            std=(epiTimeColumn, 'std'),
            n=(epiTimeColumn, 'size'))\
            .reset_index()
wallNum = 1
fig, ax = plt.subplots(1, figsize=(5, 5))
ax = avalanchesMeanAll[avalanchesMeanAll.wall_number == wallNum].plot.scatter(x='avalanche_index',
                      y='distance', c='combo_id',colormap='tab10')
# ax.set_xscale('log', base=10)
# ax.set_yscale('log', base=10)
# ax.set_ylim(-0.5, 1.5)
plt.show()

###

avalanchesMeanAll = avalanchesAll.groupby(['combo_id','avalanche_index'])\
    .agg(distance=(epiTimeColumn, 'mean'),
            std=(epiTimeColumn, 'std'),
            n=(epiTimeColumn, 'size'))\
            .reset_index()
avalanchesMean = avalanches.groupby(['combo_id','avalanche_index'])\
    .agg(distance=(epiTimeColumn, 'mean'),
            std=(epiTimeColumn, 'std'),
            n=(epiTimeColumn, 'size'))\
            .reset_index()

comboMap = cm.get_cmap('jet').reversed()
norm = Normalize(
    vmin=float(0), vmax=float(len(comboSpace['combo_id'])))

fig, ax = plt.subplots(1, figsize=(5, 5))
for i in range(0,len(comboSpace))
    ax.plot(avalanchesMeanAll['avalanche_index'],
            avalanchesMeanAll['distance'],label='spacer acq. prob. {}'.format(),
            color=comboMap(norm(float(i))), alpha=0.3)
# ax.set_xscale('log', base=10)
# ax.set_yscale('log', base=10)
# ax.set_ylim(-0.5, 1.5)
plt.show()




fig, ax = plt.subplots(1, figsize=(5, 5))
ax = avalanchesMeanAll[(avalanchesMeanAll.wall_number == wallNum)&(avalanchesMeanAll.combo_id == comboSpace['combo_id'][3])].plot.scatter(x='avalanche_index',
                      y='distance', c='combo_id',colormap='tab10')
# ax.set_xscale('log', base=10)
# ax.set_yscale('log', base=10)
# ax.set_ylim(-0.5, 1.5)
plt.show()


fig, ax = plt.subplots(1, figsize=(5, 5))
ax = avalancheSizeAvg[avalancheSizeAvg.combo_id==comboSpace['combo_id'][3]].plot.scatter(x='index',
                      y='size', c='combo_id',colormap='tab10')
ax.set_xscale('log', base=10)
ax.set_yscale('log', base=10)
# ax.set_xlim(1, 100)
plt.show()


avalanchesMean = avalanches[avalanches.wall_number == 2].groupby(['combo_id','wall_number', 'avalanche_index'])\
    .agg(distance=(epiTimeColumn, 'mean'),
            std=(epiTimeColumn, 'std'),
            n=(epiTimeColumn, 'size'))\
            .reset_index()
fig, ax = plt.subplots(1, figsize=(5, 5))
ax = avalanchesMean.plot.scatter(x='avalanche_index',
                      y='distance', c='combo_id',colormap='tab10')
# ax.set_xscale('log', base=10)
# ax.set_yscale('log', base=10)
# ax.set_xlim(1, 100)
plt.show()
