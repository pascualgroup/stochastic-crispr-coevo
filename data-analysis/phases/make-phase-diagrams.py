#!/usr/bin/env python3

import pandas as pd
import numpy as np
import scipy as sp
import math
import matplotlib.pyplot as plt
# import matplotlib
# matplotlib.use('Agg')
# from mpl_toolkits.mplot3d import Axes3D
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
xaxis = 'spacer_acquisition_prob'
xvar = 'q'
xaxlabel = 'q'
yaxis = 'viral_mutation_rate'
yvar = 'mu'
yaxlabel = r'$\mu$'
resolve = 500

#####
sweepDate = '24-2-2022'
analysisType = "walls-shannon"
xaxis = 'microbe_carrying_capacity'
xvar = 'K'
xaxlabel = r'$\phi$$K/d$'
yaxis = 'viral_burst_size'
yvar = 'beta'
yaxlabel = r'$\beta$'
resolve = 500





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
DBEXTMEAN_PATH = os.path.join('/Volumes/Yadgah/crispr-sweep-{0}/mean-extinctions.sqlite'.format(sweepDate))
DBEXT_PATH = os.path.join('/Volumes/Yadgah/crispr-sweep-{0}/extinctions.sqlite'.format(sweepDate))


conWalls = sqlite3.connect(DBWALLS_PATH)
curWalls = conWalls.cursor()
print('SQLite Query: mean microbial wall occurences per parameter combination')
wallOccurrences = pd.read_sql_query("SELECT wall_occurrence, combo_id FROM microbial_wall_statistics ORDER BY combo_id", conWalls)
numWalls = pd.read_sql_query("SELECT expected_num_walls, combo_id FROM microbial_wall_statistics ORDER BY combo_id", conWalls)
meanWallDurations = pd.read_sql_query("SELECT num_walls, mean_wallduration, num_replicates, total_replicates_of_combo, combo_id \
FROM microbial_mean_walldurations ORDER BY combo_id", conWalls)
meanWallDurations['mean_wallduration'] = meanWallDurations['mean_wallduration'] * (meanWallDurations['num_replicates']/meanWallDurations['total_replicates_of_combo'])
meanWallDurations = meanWallDurations.groupby(['combo_id']).agg(mean_wallduration=('mean_wallduration', 'sum')).reset_index()


conMExt = sqlite3.connect(DBEXTMEAN_PATH)
curMExt = conExt.cursor()
conExt = sqlite3.connect(DBEXT_PATH)
curExt = conExt.cursor()
extOccurrences = pd.read_sql_query("SELECT microbes,viruses,run_id FROM extinction_occurrence ORDER BY run_id", conExt)\
                .merge(\
                    pd.read_sql_query(
                        "SELECT run_id,combo_id FROM runs", conExt), on = ['run_id']\
                )
extOccurrences = extOccurrences.groupby(['combo_id']).agg(
    num_replicates=('run_id', 'size')).reset_index()\
    .merge(extOccurrences,on=['combo_id'])
extOccurrences = extOccurrences[extOccurrences.microbes==0]
extOccurrences = extOccurrences.groupby(['combo_id','num_replicates']).agg(
    viruses=('viruses', 'sum')).reset_index()
extOccurrences['viruses'] = extOccurrences['viruses']/extOccurrences['num_replicates']


extOccurrencesM = pd.read_sql_query("SELECT microbes,run_id FROM extinction_occurrence ORDER BY run_id", conExt)\
                .merge(\
                    pd.read_sql_query(
                        "SELECT run_id,combo_id FROM runs", conExt), on = ['run_id']\
                )
extOccurrencesM = extOccurrencesM.groupby(['combo_id']).agg(
    num_replicates=('run_id', 'size')).reset_index()\
    .merge(extOccurrencesM,on=['combo_id'])
extOccurrencesM = extOccurrencesM.groupby(['combo_id','num_replicates']).agg(
    microbes=('microbes', 'sum')).reset_index()
extOccurrencesM['microbes'] = extOccurrencesM['microbes']/extOccurrencesM['num_replicates']
extOccurrences = extOccurrences.drop(columns=['num_replicates']).merge(extOccurrencesM,on=['combo_id'])

simEndTimes = pd.read_sql_query("SELECT microbe_end_time,virus_end_time,combo_id FROM mean_simulation_end_times ORDER BY combo_id", conMExt)
# simEndTimes = pd.read_sql_query("SELECT time,combo_id FROM mean_simulation_end_times ORDER BY combo_id", conMExt)\
#                 .rename(columns={'time':'virus_end_time'}) ## for sweepDate = '7-2-2022'
# simEndTimes = pd.read_sql_query("SELECT microbe_end_time,virus_end_time,combo_id FROM mean_simulation_end_times ORDER BY combo_id", conExt)


# make this more general for other parameters!!!!
print('SQLite Query: parameter space spanning combinations')
####
####
comboSpace = pd.read_sql_query(
"SELECT combo_id, {0}, {1} \
FROM param_combos ORDER BY combo_id".format(xaxis,yaxis), conWalls)
yaxis = 'virion_mutation_prob'
g = 15  # viral repertoire length
comboSpace[yaxis] = comboSpace.viral_mutation_rate
###
###
comboSpace = pd.read_sql_query(
"SELECT combo_id, {0}, {1}, {2}, {3} \
FROM param_combos ORDER BY combo_id".format(xaxis, yaxis,'adsorption_rate','viral_decay_rate'), conWalls)
xaxis = 'maximal_adsoprtion_ratio'
comboSpace[xaxis] = comboSpace.adsorption_rate*\
        comboSpace.microbe_carrying_capacity/comboSpace.viral_decay_rate
###
###
minX = min(comboSpace[xaxis])
maxX = max(comboSpace[xaxis])
minY = min(comboSpace[yaxis])
maxY = max(comboSpace[yaxis])
grid_x, grid_y = np.mgrid[minX:maxX:10j, minY:maxY:10j]
interpMethod = 'linear'

# interpMethod ='nearest'
grid_wallOccurrences = griddata((comboSpace[xaxis].values,comboSpace[yaxis].values),
wallOccurrences.wall_occurrence.values, (grid_x, grid_y), method=interpMethod)
grid_vExt = griddata((comboSpace[xaxis].values,comboSpace[yaxis].values),
extOccurrences.viruses.values, (grid_x, grid_y), method=interpMethod)
grid_bExt = griddata((comboSpace[xaxis].values,comboSpace[yaxis].values),
extOccurrences.microbes.values, (grid_x, grid_y), method=interpMethod)

grid_numWalls = griddata((comboSpace[xaxis].values,comboSpace[yaxis].values),
numWalls.expected_num_walls.values, (grid_x, grid_y), method=interpMethod)
grid_expectedDuration = griddata((comboSpace[xaxis].values,comboSpace[yaxis].values),
meanWallDurations.mean_wallduration.values, (grid_x, grid_y), method=interpMethod)
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
cbar.ax.tick_params(labelsize=20)
ax.set_xlabel(xaxlabel, labelpad = 20, fontsize = 25)
ax.set_ylabel(yaxlabel, labelpad = 20, fontsize = 25)
ax.tick_params(axis='x', labelsize= 20)
ax.tick_params(axis='y', labelsize= 20)
text = ax.xaxis.get_offset_text() # Get the text object
text.set_size(20) # # Set the size
text = ax.yaxis.get_offset_text() # Get the text object
text.set_size(20) # # Set the size
ax.set_ylim(lowerY,upperY)
ax.set_title('Probability of Alternating Dynamics', pad = 30, fontsize = 20)
# fig.savefig(os.path.join('/Volumes/Yadgah/{0}-{1}-phase-diagram_wall-occurrences.png'.format(xvar,yvar)),dpi=resolve)
fig.savefig(os.path.join('/Users/armun/Desktop/{0}-{1}-phase-diagram_wall-occurrences.pdf'.format(xvar,yvar)),dpi=resolve)

# fig.savefig(os.path.join('/Users/armun/Desktop/check.pdf'.format(xvar,yvar)),dpi=resolve)


fig = plt.figure(figsize=(10,10))
ax = fig.add_subplot(111)
im = ax.imshow(grid_vExt.T, extent=(minX,maxX,minY,maxY), aspect='auto', origin='lower',cmap = reversed_map)
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="2%", pad=0.2)
cbar = fig.colorbar(im, cax)
cbar.ax.tick_params(labelsize=20)
ax.set_xlabel(xaxlabel, labelpad = 20, fontsize=25)
ax.set_ylabel(yaxlabel, labelpad = 20, fontsize=25)
ax.tick_params(axis='x', labelsize= 20)
ax.tick_params(axis='y', labelsize= 20)
text = ax.xaxis.get_offset_text() # Get the text object
text.set_size(20) # # Set the size
text = ax.yaxis.get_offset_text() # Get the text object
text.set_size(20) # # Set the size
ax.set_ylim(lowerY,upperY)
ax.set_title(r'Probability of Total Viral Extinction Before $t = 2000$', pad=30, fontsize = 20)
# fig.savefig(os.path.join('/Volumes/Yadgah/{0}-{1}-phase-diagram_virus-ext-occurrences.png'.format(xvar,yvar)),dpi=resolve)
fig.savefig(os.path.join('/Users/armun/Desktop/{0}-{1}-phase-diagram_virus-ext-occurrences.pdf'.format(xvar,yvar)),dpi=resolve)


fig = plt.figure(figsize=(10,10))
ax = fig.add_subplot(111)
im = ax.imshow(grid_bExt.T, extent=(minX,maxX,minY,maxY), aspect='auto', origin='lower',cmap = reversed_map)
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="2%", pad=0.2)
if sweepDate == '7-2-2022':
    cbar = fig.colorbar(im, cax, ticks=[0])
else:
    cbar = fig.colorbar(im, cax)
cbar.ax.tick_params(labelsize=20)
ax.set_xlabel(xaxlabel, labelpad = 20, fontsize=25)
ax.set_ylabel(yaxlabel, labelpad = 20, fontsize=25)
ax.tick_params(axis='x', labelsize= 20)
ax.tick_params(axis='y', labelsize= 20)
text = ax.xaxis.get_offset_text() # Get the text object
text.set_size(20) # # Set the size
text = ax.yaxis.get_offset_text() # Get the text object
text.set_size(20) # # Set the size
ax.set_ylim(lowerY,upperY)
ax.set_title(r'Probability of Total Microbial Extinction Before $t = 2000$', pad=30, fontsize = 20)
# fig.savefig(os.path.join('/Volumes/Yadgah/{0}-{1}-phase-diagram_microbe-ext-occurrences.pdf'.format(xvar,yvar)),dpi=resolve)
fig.savefig(os.path.join('/Users/armun/Desktop/{0}-{1}-phase-diagram_microbe-ext-occurrences.pdf'.format(xvar,yvar)),dpi=resolve)

fig = plt.figure(figsize=(10,10))
ax = fig.add_subplot(111)
im = ax.imshow(grid_numWalls.T, extent=[minX,maxX,minY,maxY], aspect='auto', origin='lower',cmap = reversed_map)
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="2%", pad=0.2)
cbar = fig.colorbar(im, cax)
cbar.ax.tick_params(labelsize=20)
ax.set_xlabel(xaxlabel, labelpad=20, fontsize=25)
ax.set_ylabel(yaxlabel, labelpad=20, fontsize=25)
ax.tick_params(axis='x', labelsize=20)
ax.tick_params(axis='y', labelsize=20)
text = ax.xaxis.get_offset_text()  # Get the text object
text.set_size(20)  # Set the size
text = ax.yaxis.get_offset_text()  # Get the text object
text.set_size(20)  # Set the size
ax.set_ylim(lowerY, upperY)
ax.set_title('Expected Number of SHC Periods', pad=30, fontsize=20)
# fig.savefig(os.path.join('/Volumes/Yadgah','test.png'),dpi=resolve)
fig.savefig(os.path.join('/Users/armun/Desktop/{0}-{1}-phase-diagram_expected-num-walls.pdf'.format(xvar,yvar)),dpi=resolve)


fig = plt.figure(figsize=(10,10))
ax = fig.add_subplot(111)
im = ax.imshow(grid_expectedDuration.T, extent=(minX,maxX,minY,maxY), aspect='auto', origin='lower',cmap = reversed_map)
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="2%", pad=0.2)
cbar = fig.colorbar(im, cax)
cbar.ax.tick_params(labelsize=20)
ax.set_xlabel(xaxlabel, labelpad=20, fontsize=25)
ax.set_ylabel(yaxlabel, labelpad=20, fontsize=25)
ax.tick_params(axis='x', labelsize=20)
ax.tick_params(axis='y', labelsize=20)
text = ax.xaxis.get_offset_text()  # Get the text object
text.set_size(20)  # Set the size
text = ax.yaxis.get_offset_text()  # Get the text object
text.set_size(20)  # Set the size
ax.set_ylim(lowerY, upperY)
ax.set_title('Expected Duration of SHC Periods', pad=30, fontsize=20)
# fig.savefig(os.path.join('/Volumes/Yadgah','test.png'),dpi=resolve)
fig.savefig(os.path.join('/Users/armun/Desktop/{0}-{1}-phase-diagram_expected-wall-durations.pdf'.format(xvar,yvar)),dpi=resolve)


fig = plt.figure(figsize=(10,10))
ax = fig.add_subplot(111)
im = ax.imshow(grid_vEndTimes.T, extent=(minX,maxX,minY,maxY), aspect='auto', origin='lower',cmap = reversed_map)
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="2%", pad=0.2)
cbar = fig.colorbar(im, cax)
cbar.ax.tick_params(labelsize=20)
ax.set_xlabel(xaxlabel, labelpad=20, fontsize=25)
ax.set_ylabel(yaxlabel, labelpad=20, fontsize=25)
ax.tick_params(axis='x', labelsize=20)
ax.tick_params(axis='y', labelsize=20)
text = ax.xaxis.get_offset_text()  # Get the text object
text.set_size(20)  # Set the size
text = ax.yaxis.get_offset_text()  # Get the text object
text.set_size(20)  # Set the size
ax.set_ylim(lowerY, upperY)
ax.set_title('Mean Time to Total Viral Extinction', pad=30, fontsize=20)
# fig.savefig(os.path.join('/Volumes/Yadgah','test.png'),dpi=resolve)
fig.savefig(os.path.join('/Users/armun/Desktop/{0}-{1}-phase-diagram_virus-sim-end-times.pdf'.format(xvar,yvar)),dpi=resolve)

fig = plt.figure(figsize=(10,10))
ax = fig.add_subplot(111)
im = ax.imshow(grid_bEndTimes.T, extent=(minX,maxX,minY,maxY), aspect='auto', origin='lower',cmap = reversed_map)
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="2%", pad=0.2)
cbar = fig.colorbar(im, cax)
if sweepDate == '24-2-2022':
    cbarticklabels = cbar.ax.get_yticks()
    cbarticklabels = ['{}'.format(int(i)) for i in cbarticklabels]
    cbarticklabels[-1] = r'$(2000,\infty)$'
    cbar.ax.set_yticklabels(cbarticklabels)
cbar.ax.tick_params(labelsize=20)
cbarticks = cbar.ax.get_yticklabels()
cbarticks[-1].set_fontsize(13)
ax.set_xlabel(xaxlabel, labelpad=20, fontsize=25)
ax.set_ylabel(yaxlabel, labelpad=18, fontsize=25)
ax.tick_params(axis='x', labelsize=20)
ax.tick_params(axis='y', labelsize=20)
text = ax.xaxis.get_offset_text()  # Get the text object
text.set_size(20)  # Set the size
text = ax.yaxis.get_offset_text()  # Get the text object
text.set_size(20)  # Set the size
ax.set_ylim(lowerY, upperY)
ax.set_title('Mean Time to Total Microbial Extinction', pad=30, fontsize=20)
# fig.savefig(os.path.join('/Volumes/Yadgah','test.png'),dpi=resolve)
fig.savefig(os.path.join('/Users/armun/Desktop/{0}-{1}-phase-diagram_microbe-sim-end-times.pdf'.format(xvar,yvar)),dpi=resolve)

################
################
################
################
################

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

################
################
################
################
################

DBWALLS_PATH = os.path.join('/Volumes/Yadgah/vary-q-K/mean-walls-shannon.sqlite')
DBEXTMEAN_PATH = os.path.join('/Volumes/Yadgah/vary-q-K/mean-extinctions.sqlite')
DBEXT_PATH = os.path.join('/Volumes/Yadgah/vary-q-K/extinctions.sqlite')
xaxis = 'spacer_acquisition_prob'
xvar = 'q'
xaxlabel = 'q'
yaxis = 'microbe_carrying_capacity'
yvar = 'K'
yaxlabel = r'$\phi$$K/d$'
resolve = 500

conWalls = sqlite3.connect(DBWALLS_PATH)
curWalls = conWalls.cursor()
print('SQLite Query: mean microbial wall occurences per parameter combination')
wallOccurrences = pd.read_sql_query("SELECT wall_occurrence, combo_id FROM microbial_wall_statistics ORDER BY combo_id", conWalls)
numWalls = pd.read_sql_query("SELECT expected_num_walls, combo_id FROM microbial_wall_statistics ORDER BY combo_id", conWalls)
meanWallDurations = pd.read_sql_query("SELECT num_walls, mean_wallduration, num_replicates, total_replicates_of_combo, combo_id \
FROM microbial_mean_walldurations ORDER BY combo_id", conWalls)
meanWallDurations['mean_wallduration'] = meanWallDurations['mean_wallduration'] * (meanWallDurations['num_replicates']/meanWallDurations['total_replicates_of_combo'])
meanWallDurations = meanWallDurations.groupby(['combo_id']).agg(mean_wallduration=('mean_wallduration', 'sum')).reset_index()


conMExt = sqlite3.connect(DBEXTMEAN_PATH)
curMExt = conExt.cursor()
conExt = sqlite3.connect(DBEXT_PATH)
curExt = conExt.cursor()
extOccurrences = pd.read_sql_query("SELECT microbes,viruses,run_id FROM extinction_occurrence ORDER BY run_id", conExt)\
                .merge(\
                    pd.read_sql_query(
                        "SELECT run_id,combo_id FROM runs", conExt), on = ['run_id']\
                )
extOccurrences = extOccurrences.groupby(['combo_id']).agg(
    num_replicates=('run_id', 'size')).reset_index()\
    .merge(extOccurrences,on=['combo_id'])
extOccurrences = extOccurrences[extOccurrences.microbes==0]
extOccurrences = extOccurrences.groupby(['combo_id','num_replicates']).agg(
    viruses=('viruses', 'sum')).reset_index()
extOccurrences['viruses'] = extOccurrences['viruses']/extOccurrences['num_replicates']


extOccurrencesM = pd.read_sql_query("SELECT microbes,run_id FROM extinction_occurrence ORDER BY run_id", conExt)\
                .merge(\
                    pd.read_sql_query(
                        "SELECT run_id,combo_id FROM runs", conExt), on = ['run_id']\
                )
extOccurrencesM = extOccurrencesM.groupby(['combo_id']).agg(
    num_replicates=('run_id', 'size')).reset_index()\
    .merge(extOccurrencesM,on=['combo_id'])
extOccurrencesM = extOccurrencesM.groupby(['combo_id','num_replicates']).agg(
    microbes=('microbes', 'sum')).reset_index()
extOccurrencesM['microbes'] = extOccurrencesM['microbes']/extOccurrencesM['num_replicates']
extOccurrences = extOccurrences.drop(columns=['num_replicates']).merge(extOccurrencesM,on=['combo_id'])

simEndTimes = pd.read_sql_query("SELECT microbe_end_time,virus_end_time,combo_id FROM mean_simulation_end_times ORDER BY combo_id", conMExt)
# simEndTimes = pd.read_sql_query("SELECT time,combo_id FROM mean_simulation_end_times ORDER BY combo_id", conMExt)\
#                 .rename(columns={'time':'virus_end_time'})
# simEndTimes = pd.read_sql_query("SELECT microbe_end_time,virus_end_time,combo_id FROM mean_simulation_end_times ORDER BY combo_id", conExt)


# make this more general for other parameters!!!!
print('SQLite Query: parameter space spanning combinations')
###
###
comboSpace = pd.read_sql_query(
"SELECT combo_id, {0}, {1}, {2}, {3} \
FROM param_combos ORDER BY combo_id".format(xaxis, yaxis,'adsorption_rate','viral_decay_rate'), conWalls)
yaxis = 'maximal_adsoprtion_ratio'
comboSpace[yaxis] = comboSpace.adsorption_rate*\
        comboSpace.microbe_carrying_capacity/comboSpace.viral_decay_rate
###
###
minX = min(comboSpace[xaxis])
maxX = max(comboSpace[xaxis])
minY = min(comboSpace[yaxis])
maxY = max(comboSpace[yaxis])
grid_x, grid_y = np.mgrid[minX:maxX:11j, minY:maxY:11j]
interpMethod = 'linear'

# interpMethod ='nearest'
grid_wallOccurrences = griddata((comboSpace[xaxis].values,comboSpace[yaxis].values),
wallOccurrences.wall_occurrence.values, (grid_x, grid_y), method=interpMethod)
grid_vExt = griddata((comboSpace[xaxis].values,comboSpace[yaxis].values),
extOccurrences.viruses.values, (grid_x, grid_y), method=interpMethod)
grid_bExt = griddata((comboSpace[xaxis].values,comboSpace[yaxis].values),
extOccurrences.microbes.values, (grid_x, grid_y), method=interpMethod)

grid_numWalls = griddata((comboSpace[xaxis].values,comboSpace[yaxis].values),
numWalls.expected_num_walls.values, (grid_x, grid_y), method=interpMethod)
grid_expectedDuration = griddata((comboSpace[xaxis].values,comboSpace[yaxis].values),
meanWallDurations.mean_wallduration.values, (grid_x, grid_y), method=interpMethod)
grid_vEndTimes = griddata((comboSpace[xaxis].values,comboSpace[yaxis].values),
simEndTimes.virus_end_time.values, (grid_x, grid_y), method=interpMethod)
grid_bEndTimes = griddata((comboSpace[xaxis].values,comboSpace[yaxis].values),
simEndTimes.microbe_end_time.values, (grid_x, grid_y), method=interpMethod)


reversed_map = 'jet'
lowerY = minY
upperY = maxY

fig = plt.figure(figsize=(10,10))
ax = fig.add_subplot(111)
im = ax.imshow(grid_wallOccurrences.T, extent=[minX,maxX,minY,maxY], aspect = 'auto', origin='lower',cmap = reversed_map)
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="2%", pad=0.2)
cbar = fig.colorbar(im, cax)
cbar.ax.tick_params(labelsize=20)
ax.set_xlabel(xaxlabel, labelpad = 20, fontsize = 25)
ax.set_ylabel(yaxlabel, labelpad = 20, fontsize = 25)
ax.tick_params(axis='x', labelsize= 20)
ax.tick_params(axis='y', labelsize= 20)
text = ax.xaxis.get_offset_text() # Get the text object
text.set_size(20) # # Set the size
text = ax.yaxis.get_offset_text() # Get the text object
text.set_size(20) # # Set the size
ax.set_ylim(lowerY,upperY)
ax.set_title('Probability of Alternating Dynamics', pad = 30, fontsize = 20)
# fig.savefig(os.path.join('/Volumes/Yadgah/{0}-{1}-phase-diagram_wall-occurrences.png'.format(xvar,yvar)),dpi=resolve)
fig.savefig(os.path.join('/Users/armun/Desktop/{0}-{1}-phase-diagram_wall-occurrences.pdf'.format(xvar,yvar)),dpi=resolve)

# fig.savefig(os.path.join('/Users/armun/Desktop/check.pdf'.format(xvar,yvar)),dpi=resolve)


fig = plt.figure(figsize=(10,10))
ax = fig.add_subplot(111)
im = ax.imshow(grid_vExt.T, extent=(minX,maxX,minY,maxY), aspect='auto', origin='lower',cmap = reversed_map)
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="2%", pad=0.2)
cbar = fig.colorbar(im, cax)
cbar.ax.tick_params(labelsize=20)
ax.set_xlabel(xaxlabel, labelpad = 20, fontsize=25)
ax.set_ylabel(yaxlabel, labelpad = 20, fontsize=25)
ax.tick_params(axis='x', labelsize= 20)
ax.tick_params(axis='y', labelsize= 20)
text = ax.xaxis.get_offset_text() # Get the text object
text.set_size(20) # # Set the size
text = ax.yaxis.get_offset_text() # Get the text object
text.set_size(20) # # Set the size
ax.set_ylim(lowerY,upperY)
ax.set_title(r'Probability of Total Viral Extinction Before $t = 2000$', pad=30, fontsize = 20)
# fig.savefig(os.path.join('/Volumes/Yadgah/{0}-{1}-phase-diagram_virus-ext-occurrences.png'.format(xvar,yvar)),dpi=resolve)
fig.savefig(os.path.join('/Users/armun/Desktop/{0}-{1}-phase-diagram_virus-ext-occurrences.pdf'.format(xvar,yvar)),dpi=resolve)


fig = plt.figure(figsize=(10,10))
ax = fig.add_subplot(111)
im = ax.imshow(grid_bExt.T, extent=(minX,maxX,minY,maxY), aspect='auto', origin='lower',cmap = reversed_map)
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="2%", pad=0.2)
cbar = fig.colorbar(im, cax)
cbar.ax.tick_params(labelsize=20)
cbarticks = cbar.ax.get_yticklabels()
cbarticks[-1].set_fontsize(13)
cbar.ax.tick_params(labelsize=20)
ax.set_xlabel(xaxlabel, labelpad = 20, fontsize=25)
ax.set_ylabel(yaxlabel, labelpad = 20, fontsize=25)
ax.tick_params(axis='x', labelsize= 20)
ax.tick_params(axis='y', labelsize= 20)
text = ax.xaxis.get_offset_text() # Get the text object
text.set_size(20) # # Set the size
text = ax.yaxis.get_offset_text() # Get the text object
text.set_size(20) # # Set the size
ax.set_ylim(lowerY,upperY)
ax.set_title(r'Probability of Total Microbial Extinction Before $t = 2000$', pad=30, fontsize = 20)
# fig.savefig(os.path.join('/Volumes/Yadgah/{0}-{1}-phase-diagram_microbe-ext-occurrences.pdf'.format(xvar,yvar)),dpi=resolve)
fig.savefig(os.path.join('/Users/armun/Desktop/{0}-{1}-phase-diagram_microbe-ext-occurrences.pdf'.format(xvar,yvar)),dpi=resolve)

fig = plt.figure(figsize=(10,10))
ax = fig.add_subplot(111)
im = ax.imshow(grid_numWalls.T, extent=[minX,maxX,minY,maxY], aspect='auto', origin='lower',cmap = reversed_map)
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="2%", pad=0.2)
cbar = fig.colorbar(im, cax)
cbar.ax.tick_params(labelsize=20)
ax.set_xlabel(xaxlabel, labelpad=20, fontsize=25)
ax.set_ylabel(yaxlabel, labelpad=20, fontsize=25)
ax.tick_params(axis='x', labelsize=20)
ax.tick_params(axis='y', labelsize=20)
text = ax.xaxis.get_offset_text()  # Get the text object
text.set_size(20)  # Set the size
text = ax.yaxis.get_offset_text()  # Get the text object
text.set_size(20)  # Set the size
ax.set_ylim(lowerY, upperY)
ax.set_title('Expected Number of SHC Periods', pad=30, fontsize=20)
# fig.savefig(os.path.join('/Volumes/Yadgah','test.png'),dpi=resolve)
fig.savefig(os.path.join('/Users/armun/Desktop/{0}-{1}-phase-diagram_expected-num-walls.pdf'.format(xvar,yvar)),dpi=resolve)


fig = plt.figure(figsize=(10,10))
ax = fig.add_subplot(111)
im = ax.imshow(grid_expectedDuration.T, extent=(minX,maxX,minY,maxY), aspect='auto', origin='lower',cmap = reversed_map)
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="2%", pad=0.2)
cbar = fig.colorbar(im, cax)
cbar.ax.tick_params(labelsize=20)
ax.set_xlabel(xaxlabel, labelpad=20, fontsize=25)
ax.set_ylabel(yaxlabel, labelpad=20, fontsize=25)
ax.tick_params(axis='x', labelsize=20)
ax.tick_params(axis='y', labelsize=20)
text = ax.xaxis.get_offset_text()  # Get the text object
text.set_size(20)  # Set the size
text = ax.yaxis.get_offset_text()  # Get the text object
text.set_size(20)  # Set the size
ax.set_ylim(lowerY, upperY)
ax.set_title('Expected Duration of SHC Periods', pad=30, fontsize=20)
# fig.savefig(os.path.join('/Volumes/Yadgah','test.png'),dpi=resolve)
fig.savefig(os.path.join('/Users/armun/Desktop/{0}-{1}-phase-diagram_expected-wall-durations.pdf'.format(xvar,yvar)),dpi=resolve)


fig = plt.figure(figsize=(10,10))
ax = fig.add_subplot(111)
im = ax.imshow(grid_vEndTimes.T, extent=(minX,maxX,minY,maxY), aspect='auto', origin='lower',cmap = reversed_map)
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="2%", pad=0.2)
cbar = fig.colorbar(im, cax)
cbar.ax.tick_params(labelsize=20)
ax.set_xlabel(xaxlabel, labelpad=20, fontsize=25)
ax.set_ylabel(yaxlabel, labelpad=20, fontsize=25)
ax.tick_params(axis='x', labelsize=20)
ax.tick_params(axis='y', labelsize=20)
text = ax.xaxis.get_offset_text()  # Get the text object
text.set_size(20)  # Set the size
text = ax.yaxis.get_offset_text()  # Get the text object
text.set_size(20)  # Set the size
ax.set_ylim(lowerY, upperY)
ax.set_title('Mean Time to Total Viral Extinction', pad=30, fontsize=20)
# fig.savefig(os.path.join('/Volumes/Yadgah','test.png'),dpi=resolve)
fig.savefig(os.path.join('/Users/armun/Desktop/{0}-{1}-phase-diagram_virus-sim-end-times.pdf'.format(xvar,yvar)),dpi=resolve)

fig = plt.figure(figsize=(10,10))
ax = fig.add_subplot(111)
im = ax.imshow(grid_bEndTimes.T, extent=(minX,maxX,minY,maxY), aspect='auto', origin='lower',cmap = reversed_map)
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="2%", pad=0.2)
cbar = fig.colorbar(im, cax)
# if sweepDate == '24-2-2022':
cbarticklabels = cbar.ax.get_yticks()
cbarticklabels = ['{}'.format(int(i)) for i in cbarticklabels]
cbarticklabels[-1] = r'$(2000,\infty)$'
cbar.ax.set_yticklabels(cbarticklabels)
cbar.ax.tick_params(labelsize=20)
cbarticks = cbar.ax.get_yticklabels()
cbarticks[-1].set_fontsize(13)
ax.set_xlabel(xaxlabel, labelpad=20, fontsize=25)
ax.set_ylabel(yaxlabel, labelpad=18, fontsize=25)
ax.tick_params(axis='x', labelsize=20)
ax.tick_params(axis='y', labelsize=20)
text = ax.xaxis.get_offset_text()  # Get the text object
text.set_size(20)  # Set the size
text = ax.yaxis.get_offset_text()  # Get the text object
text.set_size(20)  # Set the size
ax.set_ylim(lowerY, upperY)
ax.set_title('Mean Time to Total Microbial Extinction', pad=30, fontsize=20)
# fig.savefig(os.path.join('/Volumes/Yadgah','test.png'),dpi=resolve)
fig.savefig(os.path.join('/Users/armun/Desktop/{0}-{1}-phase-diagram_microbe-sim-end-times.pdf'.format(xvar,yvar)),dpi=resolve)