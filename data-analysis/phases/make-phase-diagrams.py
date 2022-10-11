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

if len(sys.argv[1]) != 0:
    sweepDate = sys.argv[1]
sweepDate = '7-2-2022'
analysisType = "walls-shannon"
xaxis = 'microbe_carrying_capacity'
xvar = 'K'
xaxlabel = 'Microbial Carrying Capacity: K'
yaxis = 'viral_burst_size'
yvar = 'beta'
yaxlabel = r'Viral Burst Size:   $\beta$'
resolve = 500
interpMethod ='cubic'


xaxis = 'spacer_acquisition_prob'
xvar = 'q'
xaxlabel = 'Spacer Acquisition Probability: q'
yaxis = 'viral_mutation_rate'
yvar = 'mu'
yaxlabel = r'Mutant Virion Probability:   $1 - (1- \mu)^g$'
resolve = 500
interpMethod ='cubic'


SCRIPT_PATH = os.path.abspath(os.path.dirname(__file__)) # cluster

DBSIM_PATH = os.path.join(SCRIPT_PATH,'..','..','simulation','sweep_db_gathered.sqlite') # cluster
#DBSIM_PATH = os.path.join('/Volumes','Yadgah','sweep_db_gathered.sqlite') # local

# # DBRICH_PATH = os.path.join(SCRIPT_PATH,'..','gathered-analyses','richness','richness.sqlite') # cluster
# DBSHAN_PATH = os.path.join(SCRIPT_PATH,'..','gathered-analyses','shannon','shannon.sqlite') # cluster
# DBWALLS_PATH = os.path.join(SCRIPT_PATH,'..','gathered-analyses','{}','mean-{}.sqlite'.format(analysisType)) # cluster
# DBEXT_PATH = os.path.join(SCRIPT_PATH,'..','gathered-analyses','{}','mean-{}.sqlite'.format('extinctions')) # cluster

# DBRICH_PATH = os.path.join('/Volumes','Yadgah','richness.sqlite') # local
# DBSHAN_PATH = os.path.join('/Volumes','Yadgah','shannon.sqlite') # local
# DBWALLS_PATH = os.path.join('/Volumes','Yadgah','mean-{}.sqlite'.format(analysisType)) # local
# DBWALLS_PATH = os.path.join('/Volumes','Yadgah','mean-{}_output.sqlite'.format(analysisType)) # local. combo_id fixed; for testing
# DBWALLS_PATH = os.path.join('mean-{}.sqlite'.format(analysisType)) # local
DBWALLS_PATH = os.path.join('/Volumes/Yadgah/crispr-sweep-{0}/mean-walls-shannon.sqlite'.format(sweepDate))
#DBEXT_PATH = os.path.join('/Volumes','Yadgah','mean-{}.sqlite'.format('extinctions')) # local
#DBEXT_PATH = os.path.join('/Volumes','Yadgah','mean-{}_output.sqlite'.format('extinctions')) # local. run_id fixed; for testing
# DBEXT_PATH = os.path.join('mean-{}.sqlite'.format('extinctions')) # local
DBEXT_PATH = os.path.join('/Volumes/Yadgah/crispr-sweep-{0}/mean-extinctions.sqlite'.format(sweepDate))


# DBPROB_PATH = os.path.join('/Volumes','Yadgah',dir,'probability-of-emergence_output.sqlite') # local
# DBMATCH_PATH = os.path.join('/Volumes','Yadgah',dir,'matches_output.sqlite') # local
# DBTREE_PATH = os.path.join('/Volumes','Yadgah',dir,'trees_output.sqlite') # local
# DBTRI_PATH = os.path.join('/Volumes','Yadgah',dir,'tripartite-networks_output.sqlite') # local
# PLOT_PATH = os.path.join('/Volumes','Yadgah') # local

conWalls = sqlite3.connect(DBWALLS_PATH)
curWalls = conWalls.cursor()
print('SQLite Query: mean microbial wall occurences per parameter combination')
wallOccurrences = pd.read_sql_query("SELECT wall_occurrence, combo_id FROM microbial_wall_statistics ORDER BY combo_id", conWalls)
numWalls = pd.read_sql_query("SELECT expected_num_walls, combo_id FROM microbial_wall_statistics ORDER BY combo_id", conWalls)
meanWallDurations = pd.read_sql_query("SELECT mean_wallduration, num_replicates, total_replicates_of_combo, combo_id \
FROM microbial_mean_walldurations ORDER BY combo_id", conWalls)


conExt = sqlite3.connect(DBEXT_PATH)
curExt = conExt.cursor()
extOccurrences = pd.read_sql_query("SELECT microbes,viruses,combo_id FROM mean_extinction_occurrences ORDER BY combo_id", conExt)
simEndTimes = pd.read_sql_query("SELECT microbe_end_time,virus_end_time,combo_id FROM mean_simulation_end_times ORDER BY combo_id", conExt)
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
"SELECT combo_id, {0}, {1} \
FROM param_combos ORDER BY combo_id".format(xaxis,yaxis), conWalls)
yaxis = 'virion_mutation_prob'
g = 15 # viral repertoire length
comboSpace['virion_mutation_prob'] = 1- (1- comboSpace.viral_mutation_rate)**g
minX = min(comboSpace[xaxis])
maxX = max(comboSpace[xaxis])
minY = min(comboSpace[yaxis])
maxY = max(comboSpace[yaxis])
grid_x, grid_y = np.mgrid[minX:maxX:500j, minY:maxY:500j]

grid_wallOccurrences = griddata((comboSpace[xaxis].values,comboSpace[yaxis].values),
wallOccurrences.wall_occurrence.values, (grid_x, grid_y), method=interpMethod)
grid_vExt = griddata((comboSpace[xaxis].values,comboSpace[yaxis].values),
extOccurrences.viruses.values, (grid_x, grid_y), method=interpMethod)
grid_bExt = griddata((comboSpace[xaxis].values,comboSpace[yaxis].values),
extOccurrences.microbes.values, (grid_x, grid_y), method=interpMethod)

grid_numWalls = griddata((comboSpace[xaxis].values,comboSpace[yaxis].values),
numWalls.expected_num_walls.values, (grid_x, grid_y), method=interpMethod)
grid_expectedDuration = griddata((comboSpace[xaxis].values,comboSpace[yaxis].values),
expectedWallDurations.expected_wallduration.values, (grid_x, grid_y), method=interpMethod)
grid_vEndTimes = griddata((comboSpace[xaxis].values,comboSpace[yaxis].values),
simEndTimes.virus_end_time.values, (grid_x, grid_y), method=interpMethod)
grid_bEndTimes = griddata((comboSpace[xaxis].values,comboSpace[yaxis].values),
simEndTimes.microbe_end_time.values, (grid_x, grid_y), method=interpMethod)


orig_map=plt.cm.get_cmap('viridis')
# reversing the original colormap using reversed() function
reversed_map = orig_map.reversed()
reversed_map = 'jet'


lowerY = minY
upperY = maxY
fig = plt.figure(figsize=(10,10))
ax = fig.add_subplot(111)
im = ax.imshow(grid_wallOccurrences.T, extent=[minX,maxX,minY,maxY], aspect = 'auto', origin='lower',cmap = reversed_map)
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="2%", pad=0.2)
cbar = fig.colorbar(im, cax)
cbar.ax.tick_params(labelsize=13)
ax.set_xlabel(xaxlabel, labelpad = 20, fontsize = 20)
ax.set_ylabel(yaxlabel, labelpad = 20, fontsize = 20)
ax.tick_params(axis='x', labelsize= 13)
ax.tick_params(axis='y', labelsize= 13)
text = ax.xaxis.get_offset_text() # Get the text object
text.set_size(13) # # Set the size
text = ax.yaxis.get_offset_text() # Get the text object
text.set_size(13) # # Set the size
ax.set_ylim(lowerY,upperY)
ax.set_title('Probability of Sustained Host Control', pad = 30, fontsize = 30)
fig.savefig(os.path.join('/Volumes/Yadgah/{0}-{1}-phase-diagram_wall-occurrences.png'.format(xvar,yvar)),dpi=resolve)
fig.savefig(os.path.join('{0}-{1}-phase-diagram_wall-occurrences.png'.format(xvar,yvar)),dpi=resolve)


fig = plt.figure(figsize=(10,10))
ax = fig.add_subplot(111)
im = ax.imshow(grid_vExt.T, extent=(minX,maxX,minY,maxY), aspect='auto', origin='lower',cmap = reversed_map)
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="2%", pad=0.2)
cbar = fig.colorbar(im, cax)
cbar.ax.tick_params(labelsize=13)
ax.set_xlabel(xaxlabel, labelpad = 20, fontsize = 20)
ax.set_ylabel(yaxlabel, labelpad = 20, fontsize = 20)
ax.tick_params(axis='x', labelsize= 13)
ax.tick_params(axis='y', labelsize= 13)
text = ax.xaxis.get_offset_text() # Get the text object
text.set_size(13) # # Set the size
text = ax.yaxis.get_offset_text() # Get the text object
text.set_size(13) # # Set the size
ax.set_ylim(lowerY,upperY)
ax.set_title('Probability of Viral Extinction', pad=30, fontsize = 30)
fig.savefig(os.path.join('/Volumes/Yadgah/{0}-{1}-phase-diagram_virus-ext-occurrences.png'.format(xvar,yvar)),dpi=resolve)
fig.savefig(os.path.join('{0}-{1}-phase-diagram_virus-ext-occurrences.png'.format(xvar,yvar)),dpi=resolve)


fig = plt.figure(figsize=(10,10))
ax = fig.add_subplot(111)
im = ax.imshow(grid_bExt.T, extent=(minX,maxX,minY,maxY), aspect='auto', origin='lower',cmap = reversed_map)
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="2%", pad=0.2)
cbar = fig.colorbar(im, cax)
cbar.ax.tick_params(labelsize=13)
ax.set_xlabel(xaxlabel, labelpad = 20, fontsize = 20)
ax.set_ylabel(yaxlabel, labelpad = 20, fontsize = 20)
ax.tick_params(axis='x', labelsize= 13)
ax.tick_params(axis='y', labelsize= 13)
text = ax.xaxis.get_offset_text() # Get the text object
text.set_size(13) # # Set the size
text = ax.yaxis.get_offset_text() # Get the text object
text.set_size(13) # # Set the size
ax.set_ylim(lowerY,upperY)
ax.set_title('Probability of Microbial Extinction', pad=30, fontsize = 30)
fig.savefig(os.path.join('/Volumes/Yadgah/{0}-{1}-phase-diagram_microbe-ext-occurrences.png'.format(xvar,yvar)),dpi=resolve)
fig.savefig(os.path.join('{0}-{1}-phase-diagram_microbe-ext-occurrences.png'.format(xvar,yvar)),dpi=resolve)




fig = plt.figure(figsize=(10,10))
ax = fig.add_subplot(111)
im = ax.imshow(grid_numWalls.T, extent=[minX,maxX,minY,maxY], aspect='auto', origin='lower',cmap = reversed_map)
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="2%", pad=0.2)
fig.colorbar(im, cax)
ax.set_xlabel(xaxlabel, labelpad = 20)
ax.set_ylabel(yaxlabel, labelpad = 20)
ax.set_ylim(lowerY,upperY)
ax.set_title('Expected Number of Walls', pad = 20)
# fig.savefig(os.path.join('/Volumes/Yadgah','test.png'),dpi=resolve)
fig.savefig(os.path.join('{0}-{1}-phase-diagram_expected-num-walls.png'.format(xvar,yvar)),dpi=resolve)


fig = plt.figure(figsize=(10,10))
ax = fig.add_subplot(111)
im = ax.imshow(grid_expectedDuration.T, extent=(minX,maxX,minY,maxY), aspect='auto', origin='lower',cmap = reversed_map)
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="2%", pad=0.2)
fig.colorbar(im, cax)
ax.set_xlabel(xaxlabel, labelpad = 20)
ax.set_ylabel(yaxlabel, labelpad = 20)
ax.set_ylim(lowerY,upperY)
ax.set_title('Expected Duration of Walls', pad = 20)
# fig.savefig(os.path.join('/Volumes/Yadgah','test.png'),dpi=resolve)
fig.savefig(os.path.join('{0}-{1}-phase-diagram_expected-wall-durations.png'.format(xvar,yvar)),dpi=resolve)


fig = plt.figure(figsize=(10,10))
ax = fig.add_subplot(111)
im = ax.imshow(grid_vEndTimes.T, extent=(minX,maxX,minY,maxY), aspect='auto', origin='lower',cmap = reversed_map)
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="2%", pad=0.2)
fig.colorbar(im, cax)
ax.set_xlabel(xaxlabel, labelpad = 20)
ax.set_ylabel(yaxlabel, labelpad = 20)
ax.set_ylim(lowerY,upperY)
ax.set_title('Mean Viral Simulation End Time', pad=20)
# fig.savefig(os.path.join('/Volumes/Yadgah','test.png'),dpi=resolve)
fig.savefig(os.path.join('{0}-{1}-phase-diagram_virus-sim-end-times.png'.format(xvar,yvar)),dpi=resolve)

fig = plt.figure(figsize=(10,10))
ax = fig.add_subplot(111)
im = ax.imshow(grid_bEndTimes.T, extent=(minX,maxX,minY,maxY), aspect='auto', origin='lower',cmap = reversed_map)
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="2%", pad=0.2)
fig.colorbar(im, cax)
ax.set_xlabel(xaxlabel, labelpad = 20)
ax.set_ylabel(yaxlabel, labelpad = 20)
ax.set_ylim(lowerY,upperY)
ax.set_title('Mean Microbial Simulation End Time', pad=20)
# fig.savefig(os.path.join('/Volumes/Yadgah','test.png'),dpi=resolve)
fig.savefig(os.path.join('{0}-{1}-phase-diagram_microbe-sim-end-times.png'.format(xvar,yvar)),dpi=resolve)













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
grid_expectedDuration = griddata((comboSpace.spacer_acquisition_prob.values,comboSpace.virion_mutation_prob.values),
expectedWallDurations.expected_wallduration.values, (grid_x, grid_y), method='nearest')

grid_vExt = griddata((comboSpace.spacer_acquisition_prob.values,comboSpace.virion_mutation_prob.values),
extOccurrences.viruses.values, (grid_x, grid_y), method='nearest')
grid_mExt = griddata((comboSpace.spacer_acquisition_prob.values,comboSpace.virion_mutation_prob.values),
extOccurrences.microbes.values, (grid_x, grid_y), method='nearest')
grid_vEndTimes = griddata((comboSpace.spacer_acquisition_prob.values,comboSpace.virion_mutation_prob.values),
simEndTimes.time.values, (grid_x, grid_y), method='nearest')

# grid_vEndTimes = griddata((comboSpace.spacer_acquisition_prob.values,comboSpace.virion_mutation_prob.values),
# simEndTimes.virus_end_time.values, (grid_x, grid_y), method='nearest')
#
# grid_mEndTimes = griddata((comboSpace.spacer_acquisition_prob.values,comboSpace.virion_mutation_prob.values),
# simEndTimes.microbe_end_time.values, (grid_x, grid_y), method='nearest')

orig_map=plt.cm.get_cmap('viridis')
# reversing the original colormap using reversed() function
reversed_map = orig_map
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
