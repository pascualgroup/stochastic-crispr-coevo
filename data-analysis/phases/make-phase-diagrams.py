#!/usr/bin/env python3

import pandas as pd
import numpy as np
import scipy as sp
import math
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')
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

analysisType = "walls-shannon"

SCRIPT_PATH = os.path.abspath(os.path.dirname(__file__)) # cluster

DBSIM_PATH = os.path.join(SCRIPT_PATH,'..','..','simulation','sweep_db_gathered.sqlite') # cluster
#home_dir = os.system("cd ~") # local
#DBSIM_PATH = os.path.join('/Volumes','Yadgah','sweep_db_gathered.sqlite') # local

# DBRICH_PATH = os.path.join(SCRIPT_PATH,'..','gathered-analyses','richness','richness.sqlite') # cluster
# DBRICH_PATH = os.path.join('/Volumes','Yadgah','richness.sqlite') # local

DBSHAN_PATH = os.path.join(SCRIPT_PATH,'..','gathered-analyses','shannon','shannon.sqlite') # cluster
#DBSHAN_PATH = os.path.join('/Volumes','Yadgah','shannon.sqlite') # local

DBWALLS_PATH = os.path.join(SCRIPT_PATH,'..','gathered-analyses','{}','mean-{}.sqlite'.format(analysisType)) # cluster
# DBWALLS_PATH = os.path.join('/Volumes','Yadgah','mean-{}.sqlite'.format(analysisType)) # local
# DBWALLS_PATH = os.path.join('/Volumes','Yadgah','mean-{}_output.sqlite'.format(analysisType)) # local. combo_id fixed; for testing
# DBWALLS_PATH = os.path.join('mean-{}.sqlite'.format(analysisType)) # local

DBEXT_PATH = os.path.join(SCRIPT_PATH,'..','gathered-analyses','{}','mean-{}.sqlite'.format('extinctions')) # cluster
#DBEXT_PATH = os.path.join('/Volumes','Yadgah','mean-{}.sqlite'.format('extinctions')) # local
#DBEXT_PATH = os.path.join('/Volumes','Yadgah','mean-{}_output.sqlite'.format('extinctions')) # local. run_id fixed; for testing
# DBEXT_PATH = os.path.join('mean-{}.sqlite'.format('extinctions')) # local

conWalls = sqlite3.connect(DBWALLS_PATH)
curWalls = conWalls.cursor()
print('SQLite Query: mean microbial wall occurences per parameter combination')
extOccurrences = pd.read_sql_query("SELECT wall_occurrence, combo_id FROM microbial_wall_statistics ORDER BY combo_id", conWalls)
numWalls = pd.read_sql_query("SELECT expected_num_walls, combo_id FROM microbial_wall_statistics ORDER BY combo_id", conWalls)
meanWallDurations = pd.read_sql_query("SELECT mean_wallduration, num_replicates, total_replicates_of_combo, combo_id \
FROM microbial_mean_walldurations ORDER BY combo_id", conWalls)


conExt = sqlite3.connect(DBEXT_PATH)
curExt = conExt.cursor()
extOccurrences = pd.read_sql_query("SELECT microbes,viruses,combo_id FROM mean_extinction_occurrences ORDER BY combo_id", conExt)
simEndTimes = pd.read_sql_query("SELECT time,combo_id FROM mean_simulation_end_times ORDER BY combo_id", conExt)
# simEndTimes = pd.read_sql_query("SELECT microbe_end_time,virus_end_time,combo_id FROM mean_simulation_end_times ORDER BY combo_id", conExt)


expectedWallDurations = pd.DataFrame(columns = ['expected_wallduration', 'combo_id'])
for i in np.unique(meanWallDurations.combo_id):
    mwd = (meanWallDurations[meanWallDurations['combo_id']==i].num_replicates/ \
    meanWallDurations[meanWallDurations['combo_id']==i].total_replicates_of_combo)
    mwd = mwd*meanWallDurations[meanWallDurations['combo_id']==i].mean_wallduration
    expectedWallDurations = expectedWallDurations.append({'expected_wallduration' : sum(mwd), 'combo_id' : i},
                    ignore_index = True)



# make this more general for other parameters!!!!
print('SQLite Query: parameter space spanning combinations')
comboSpace = pd.read_sql_query(
"SELECT combo_id, crispr_failure_prob, spacer_acquisition_prob, viral_mutation_rate \
FROM param_combos ORDER BY combo_id", conWalls)

g = 15 # viral repertoire length
comboSpace['virion_mutation_prob'] = 1- (1- comboSpace.viral_mutation_rate)**g

minX = 1e-06
maxX = 2.8e-05
minY = 5.9e-07
maxY = 2.76e-5
grid_x, grid_y = np.mgrid[minX:maxX:20j, minY:maxY:20j]

grid_wallOccurrences = griddata((comboSpace.spacer_acquisition_prob.values,comboSpace.virion_mutation_prob.values),
wallOccurrences.wall_occurrence.values, (grid_x, grid_y), method='nearest')

grid_numWalls = griddata((comboSpace.spacer_acquisition_prob.values,comboSpace.virion_mutation_prob.values),
numWalls.expected_num_walls.values, (grid_x, grid_y), method='nearest')

grid_vExt = griddata((comboSpace.spacer_acquisition_prob.values,comboSpace.virion_mutation_prob.values),
extOccurrences.viruses.values, (grid_x, grid_y), method='nearest')

grid_mExt = griddata((comboSpace.spacer_acquisition_prob.values,comboSpace.virion_mutation_prob.values),
extOccurrences.microbes.values, (grid_x, grid_y), method='nearest')

grid_expectedDuration = griddata((comboSpace.spacer_acquisition_prob.values,comboSpace.virion_mutation_prob.values),
expectedWallDurations.expected_wallduration.values, (grid_x, grid_y), method='nearest')

grid_vEndTimes = griddata((comboSpace.spacer_acquisition_prob.values,comboSpace.virion_mutation_prob.values),
simEndTimes.time.values, (grid_x, grid_y), method='nearest')

# grid_vEndTimes = griddata((comboSpace.spacer_acquisition_prob.values,comboSpace.virion_mutation_prob.values),
# simEndTimes.virus_end_time.values, (grid_x, grid_y), method='nearest')
#
# grid_mEndTimes = griddata((comboSpace.spacer_acquisition_prob.values,comboSpace.virion_mutation_prob.values),
# simEndTimes.microbe_end_time.values, (grid_x, grid_y), method='nearest')

orig_map=plt.cm.get_cmap('viridis')
# reversing the original colormap using reversed() function
reversed_map = orig_map.reversed()
lowerY = 2.0e-06
upperY = 2.8e-05


fig = plt.figure(figsize=(10,10))
ax = fig.add_subplot(111)
im = ax.imshow(grid_wallOccurrences.T, extent=(minX,maxX,minY,maxY), aspect=1, origin='lower',cmap = reversed_map)
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="2%", pad=0.2)
fig.colorbar(im, cax)
ax.set_xlabel('spacer acquisition probabiliity: q', labelpad = 20)
ax.set_ylabel(r'mutant virion probability:   $1-(1-\mu)^g$', labelpad = 20)
ax.set_ylim(lowerY,uppperY)
ax.set_title('Probability of Wall Occurence')
# fig.savefig(os.path.join('/Volumes/Yadgah','test.png'),dpi=500)

fig.savefig(os.path.join('q-mu-phase-diagram_wall-occurrences.png'),dpi=500)


fig = plt.figure(figsize=(10,10))
ax = fig.add_subplot(111)
im = ax.imshow(grid_numWalls.T, extent=(minX,maxX,minY,maxY), aspect=1, origin='lower',cmap = reversed_map)
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="2%", pad=0.2)
fig.colorbar(im, cax)
ax.set_xlabel('spacer acquisition probabiliity: q', labelpad = 20)
ax.set_ylabel(r'mutant virion probability:   $1-(1-\mu)^g$', labelpad = 20)
ax.set_ylim(lowerY,uppperY)
ax.set_title('Expected Number of Walls')
# fig.savefig(os.path.join('/Volumes/Yadgah','test.png'),dpi=500)

fig.savefig(os.path.join('q-mu-phase-diagram_expected-num-walls.png'),dpi=500)


fig = plt.figure(figsize=(10,10))
ax = fig.add_subplot(111)
im = ax.imshow(grid_expectedDuration.T, extent=(minX,maxX,minY,maxY), aspect=1, origin='lower',cmap = reversed_map)
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="2%", pad=0.2)
fig.colorbar(im, cax)
ax.set_xlabel('spacer acquisition probabiliity: q', labelpad = 20)
ax.set_ylabel(r'mutant virion probability:   $1-(1-\mu)^g$', labelpad = 20)
ax.set_ylim(lowerY,uppperY)
ax.set_title('Expected Duration of Walls')
# fig.savefig(os.path.join('/Volumes/Yadgah','test.png'),dpi=500)

fig.savefig(os.path.join('q-mu-phase-diagram_expected-wall-durations.png'),dpi=500)


fig = plt.figure(figsize=(10,10))
ax = fig.add_subplot(111)
im = ax.imshow(grid_vExt.T, extent=(minX,maxX,minY,maxY), aspect=1, origin='lower',cmap = reversed_map)
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="2%", pad=0.2)
fig.colorbar(im, cax)
ax.set_xlabel('spacer acquisition probabiliity: q', labelpad = 20)
ax.set_ylabel(r'mutant virion probability:   $1-(1-\mu)^g$', labelpad = 20)
ax.set_ylim(lowerY,uppperY)
ax.set_title('Probability of Viral Extinction')
# fig.savefig(os.path.join('/Volumes/Yadgah','test.png'),dpi=500)

fig.savefig(os.path.join('q-mu-phase-diagram_virus-ext-occurrences.png'),dpi=500)

fig = plt.figure(figsize=(10,10))
ax = fig.add_subplot(111)
im = ax.imshow(grid_vEndTimes.T, extent=(minX,maxX,minY,maxY), aspect=1, origin='lower',cmap = reversed_map)
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="2%", pad=0.2)
fig.colorbar(im, cax)
ax.set_xlabel('spacer acquisition probabiliity: q', labelpad = 20)
ax.set_ylabel(r'mutant virion probability:   $1-(1-\mu)^g$', labelpad = 20)
ax.set_ylim(lowerY,uppperY)
ax.set_title('Expected Time Until Viral Extinction')
# fig.savefig(os.path.join('/Volumes/Yadgah','test.png'),dpi=500)

fig.savefig(os.path.join('q-mu-phase-diagram_virus-sim-end-times.png'),dpi=500)









# wallOccurrences = wallOccurrences[(wallOccurrences.combo_id >32)]
# comboSpace = comboSpace[(comboSpace.combo_id >32)]
# comboSpace = comboSpace[(comboSpace.crispr_failure_prob == 1e-05)]
# wallOccurrences = wallOccurrences[wallOccurrences['combo_id'].isin(comboSpace.combo_id.values)]
