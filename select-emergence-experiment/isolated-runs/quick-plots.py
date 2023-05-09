#!/usr/bin/env python3

import pandas as pd
import numpy as np
import scipy as sp
import math
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import ticker
import sys
import os
import seaborn as sns
from scipy import stats
import sqlite3
from scipy.interpolate import griddata
from scipy.interpolate import interp1d
from mpl_toolkits.axes_grid1 import make_axes_locatable

dir = 'sylvain-martin-collab/trial-3/gathered-analyses'
DBWALLS_PATH = os.path.join('/Volumes/Yadgah/{0}/walls-shannon/mean-walls-shannon.sqlite'.format(dir))
DBEXT_PATH = os.path.join('/Volumes/Yadgah/{0}/extinctions/mean-extinctions.sqlite'.format(dir))
conWalls = sqlite3.connect(DBWALLS_PATH)
curWalls = conWalls.cursor()
comboSpace = pd.read_sql_query(
"SELECT combo_id, max_allele \
FROM param_combos ORDER BY combo_id", conWalls)
print('SQLite Query: mean microbial wall occurences per parameter combination')
wallOccurrences = pd.read_sql_query("SELECT wall_occurrence, combo_id FROM microbial_wall_statistics ORDER BY combo_id", conWalls)
wallOccurrences = wallOccurrences.merge(comboSpace,on=['combo_id']).drop(columns=['combo_id'])\
                    .sort_values(by=['max_allele'])
numWalls = pd.read_sql_query("SELECT expected_num_walls, combo_id FROM microbial_wall_statistics ORDER BY combo_id", conWalls)
numWalls = numWalls.merge(comboSpace,on=['combo_id']).drop(columns=['combo_id'])\
                    .sort_values(by=['max_allele'])
conExt = sqlite3.connect(DBEXT_PATH)
curExt = conExt.cursor()
extOccurrences = pd.read_sql_query("SELECT microbes,viruses,combo_id FROM mean_extinction_occurrences ORDER BY combo_id", conExt)
extOccurrences = extOccurrences.merge(comboSpace,on=['combo_id']).drop(columns=['combo_id'])\
                    .sort_values(by=['max_allele'])
simEndTimes = pd.read_sql_query("SELECT microbe_end_time,virus_end_time,combo_id FROM mean_simulation_end_times ORDER BY combo_id", conExt)
simEndTimes  = simEndTimes .merge(comboSpace,on=['combo_id']).drop(columns=['combo_id'])\
                    .sort_values(by=['max_allele'])

fig, ax = plt.subplots(4,sharex=True)
ax[0].bar(wallOccurrences['max_allele'], wallOccurrences['wall_occurrence'], color ='darkred',
        width = 2)
ax[0].set_xticks([30,35,45,60,80,105,135])
ax[0].set_xticklabels([30,35,45,60,80,105,135])
ax[1].bar(numWalls['max_allele'], numWalls['expected_num_walls'], color ='darkorange',
        width = 2)
ax[1].set_xticks([30,35,45,60,80,105,135])
ax[1].set_xticklabels([30,35,45,60,80,105,135])
ax[2].bar(extOccurrences['max_allele'], extOccurrences['viruses'], color ='darkgreen',
        width = 2)
ax[2].set_xticks([30,35,45,60,80,105,135])
ax[2].set_xticklabels([30,35,45,60,80,105,135])
ax[3].bar(simEndTimes['max_allele'], simEndTimes['virus_end_time'], color ='darkblue',
        width = 2)
ax[3].set_xticks([30,35,45,60,80,105,135])
ax[3].set_xticklabels([30,35,45,60,80,105,135])
ax[0].set_ylabel(ylabel ='Expected Wall Occurrence',labelpad=15,fontsize=7)
ax[1].set_ylabel(ylabel ='Expected Number of Walls',labelpad=15,fontsize=7)
ax[2].set_ylabel(ylabel ='Expected Extinction Occurrence',labelpad=15,fontsize=7)
ax[3].set_ylabel(ylabel ='Expected Time to Viral Extinction',labelpad=15,fontsize=7)
ax[3].set_xlabel(xlabel = 'Max Allele',fontsize=7)


def fgm(a,centerA,maxA,maxFit):
    minA = 2*centerA-maxA
    scale = -1*maxFit/((centerA-minA)*(centerA - maxA))
    return -1*(a-minA)*(a-maxA)*scale



a = list(np.linspace(-30, 30, 100))
centerA = 0
maxFit = 1
for maxA in [30,35,45,60,80,105,135]:
    fig, ax = plt.subplots(1,sharex=True)
    ax.plot(a,list(map(lambda x: fgm(x,centerA,maxA,maxFit),a)),label='max allele: {0}'.format(maxA),\
            color='darkblue',linewidth=1.5)
    ax.set_ylim(0,1)
    ax.legend()
    fig.savefig(os.path.join('/Users/armun/Desktop/fgm{0}.png'.format(maxA)),dpi=500)

# fig, ax = plt.subplots(1,sharex=True)
# ax.plot([0,1e-5,1e-3,1e-2],[7/30,7/30,8/30,11/30],color='darkred',linewidth=2)
# ax.set_xticks([0,1e-5,1e-3,1e-2])
# ax.set_xticklabels([0,1e-5,1e-3,1e-2])
# ax.set_ylabel(ylabel ='Proportion of Realizations',labelpad=15,fontsize=7)
#
# fig, ax = plt.subplots(1,sharex=True)
# ax.plot([0,1e-5,1e-3,1e-2],[5/30,6/30,4/30,1/30],color='darkblue',linewidth=2)
# ax.set_xticks([0,1e-5,1e-3,1e-2])
# ax.set_xticklabels([0,1e-5,1e-3,1e-2])
# ax.set_ylabel(ylabel ='Proportion of Realizations',labelpad=15,fontsize=7)
