#!/usr/bin/env python3

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
from scipy.interpolate import griddata
from scipy.interpolate import interp1d
from mpl_toolkits.axes_grid1 import make_axes_locatable



PLOT_PATH = os.path.join('/Volumes/Yadgah')
resolve = 500

for sweepDate in ['7-2-2022','24-2-2022']:

    if sweepDate == '7-2-2022':
        xaxis = 'spacer_acquisition_prob'
        xvar = 'q'
        xaxlabel = 'q'
        yaxis = 'viral_mutation_rate'
        yvar = 'mu'
        yaxlabel = r'$\mu$'
        
    if sweepDate == '24-2-2022':
        analysisType = "walls-shannon"
        xaxis = 'microbe_carrying_capacity'
        xvar = 'K'
        xaxlabel = r'$\phi$$K/d$'
        yaxis = 'viral_burst_size'
        yvar = 'beta'
        yaxlabel = r'$\beta$'
    DBWALLS_PATH = os.path.join('/Volumes/Yadgah/crispr-sweep-{0}/mean-walls-shannon.sqlite'.format(sweepDate))
    DBEXTMEAN_PATH = os.path.join('/Volumes/Yadgah/crispr-sweep-{0}/mean-extinctions.sqlite'.format(sweepDate))
    DBEXT_PATH = os.path.join('/Volumes/Yadgah/crispr-sweep-{0}/extinctions.sqlite'.format(sweepDate))

    conWalls = sqlite3.connect(DBWALLS_PATH)
    curWalls = conWalls.cursor()
    print('SQLite Query: mean microbial wall occurences per parameter combination')
    wallOccurrences = pd.read_sql_query("SELECT wall_occurrence, combo_id FROM microbial_wall_statistics ORDER BY combo_id", conWalls)
    numWalls = pd.read_sql_query("SELECT expected_num_walls, combo_id FROM microbial_wall_statistics ORDER BY combo_id", conWalls)
    meanWallDurations = pd.read_sql_query("SELECT num_walls, mean_wallduration, num_replicates, total_replicates, combo_id \
    FROM microbial_mean_walldurations ORDER BY combo_id", conWalls)
    meanWallDurations['mean_wallduration'] = meanWallDurations['mean_wallduration'] * (meanWallDurations['num_replicates']/meanWallDurations['total_replicates'])
    meanWallDurations = meanWallDurations.groupby(['combo_id']).agg(mean_wallduration=('mean_wallduration', 'sum')).reset_index()


    conMExt = sqlite3.connect(DBEXTMEAN_PATH)
    curMExt = conMExt.cursor()
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

    if sweepDate == '24-2-2022':
        simEndTimes = pd.read_sql_query("SELECT microbe_end_time,virus_end_time,combo_id FROM mean_simulation_end_times ORDER BY combo_id", conMExt)
    if sweepDate == '7-2-2022':
        simEndTimes = pd.read_sql_query("SELECT time,combo_id FROM mean_simulation_end_times ORDER BY combo_id", conMExt)\
                        .rename(columns={'time':'virus_end_time'}) ## for sweepDate = '7-2-2022'
    # simEndTimes = pd.read_sql_query("SELECT microbe_end_time,virus_end_time,combo_id FROM mean_simulation_end_times ORDER BY combo_id", conExt)

    if sweepDate == '7-2-2022':
        comboSpace = pd.read_sql_query(
        "SELECT combo_id, {0}, {1} \
        FROM param_combos ORDER BY combo_id".format(xaxis,yaxis), conWalls)
        yaxis = 'virion_mutation_rate'
        comboSpace[yaxis] = comboSpace.viral_mutation_rate
    ###
    if sweepDate == '24-2-2022':
        comboSpace = pd.read_sql_query(
        "SELECT combo_id, {0}, {1}, {2}, {3} \
        FROM param_combos ORDER BY combo_id".format(xaxis, yaxis,'adsorption_rate','viral_decay_rate'), conWalls)
        xaxis = 'maximal_adsorption_ratio'
        comboSpace[xaxis] = comboSpace.adsorption_rate*\
                comboSpace.microbe_carrying_capacity/comboSpace.viral_decay_rate
    ###
    minX = min(comboSpace[xaxis])
    maxX = max(comboSpace[xaxis])
    minY = min(comboSpace[yaxis])
    maxY = max(comboSpace[yaxis])
    grid_x, grid_y = np.mgrid[minX:maxX:10j, minY:maxY:10j]
    interpMethod = 'linear'
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
    if sweepDate == '24-2-2022':
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
    ax.set_xlabel('{}'.format(xaxlabel), labelpad = 20, fontsize = 25)
    ax.set_ylabel('{}'.format(yaxlabel), labelpad = 20, fontsize = 25)
    ax.tick_params(axis='x', labelsize= 20)
    ax.tick_params(axis='y', labelsize= 20)
    text = ax.xaxis.get_offset_text() # Get the text object
    text.set_size(20) # # Set the size
    text = ax.yaxis.get_offset_text() # Get the text object
    text.set_size(20) # # Set the size
    ax.set_ylim(lowerY,upperY)
    ax.set_title('Probability of Alternating Dynamics', pad = 30, fontsize = 20)
    fig.savefig(os.path.join(PLOT_PATH,'{0}-{1}-phase-diagram_wall-occurrences.pdf'.format(xvar,yvar)),dpi=resolve)
    plt.close('all')


    fig = plt.figure(figsize=(10,10))
    ax = fig.add_subplot(111)
    im = ax.imshow(grid_vExt.T, extent=(minX,maxX,minY,maxY), aspect='auto', origin='lower',cmap = reversed_map)
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="2%", pad=0.2)
    cbar = fig.colorbar(im, cax)
    cbar.ax.tick_params(labelsize=20)
    ax.set_xlabel(xaxlabel, labelpad = 20, fontsize=25)
    ax.set_ylabel(yaxlabel, labelpad=20, fontsize=25)
    ax.tick_params(axis='x', labelsize= 20)
    ax.tick_params(axis='y', labelsize= 20)
    text = ax.xaxis.get_offset_text() # Get the text object
    text.set_size(20) # # Set the size
    text = ax.yaxis.get_offset_text() # Get the text object
    text.set_size(20) # # Set the size
    ax.set_ylim(lowerY,upperY)
    ax.set_title(r'Probability of Total Viral Extinction Before $t = 2000$', pad=30, fontsize = 20)
    fig.savefig(os.path.join(PLOT_PATH,'{0}-{1}-phase-diagram_virus-ext-occurrences.pdf'.format(xvar,yvar)),dpi=resolve)
    plt.close('all')


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
    fig.savefig(os.path.join(PLOT_PATH, '{0}-{1}-phase-diagram_microbe-ext-occurrences.pdf'.format(xvar,yvar)),dpi=resolve)
    plt.close('all')

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
    fig.savefig(os.path.join(PLOT_PATH, '{0}-{1}-phase-diagram_expected-num-walls.pdf'.format(xvar,yvar)),dpi=resolve)
    plt.close('all')


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
    fig.savefig(os.path.join(PLOT_PATH,'{0}-{1}-phase-diagram_expected-wall-durations.pdf'.format(xvar,yvar)),dpi=resolve)
    plt.close('all')


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
    fig.savefig(os.path.join(PLOT_PATH,'{0}-{1}-phase-diagram_virus-sim-end-times.pdf'.format(xvar,yvar)),dpi=resolve)
    plt.close('all')

    if sweepDate == '24-2-2022':
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
        fig.savefig(os.path.join(PLOT_PATH,'{0}-{1}-phase-diagram_microbe-sim-end-times.pdf'.format(xvar,yvar)),dpi=resolve)
        plt.close('all')

########################################################################
# spacer acqusition probability q vs. maximal adsorption ratio phi*K/d #
########################################################################
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
meanWallDurations = pd.read_sql_query("SELECT num_walls, mean_wallduration, num_replicates, total_replicates, combo_id \
FROM microbial_mean_walldurations ORDER BY combo_id", conWalls)
meanWallDurations['mean_wallduration'] = meanWallDurations['mean_wallduration'] * (meanWallDurations['num_replicates']/meanWallDurations['total_replicates'])
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
fig.savefig(os.path.join(PLOT_PATH,'{0}-{1}-phase-diagram_wall-occurrences.pdf'.format(xvar,yvar)),dpi=resolve)
plt.close('all')


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
fig.savefig(os.path.join(PLOT_PATH,'{0}-{1}-phase-diagram_virus-ext-occurrences.pdf'.format(xvar,yvar)),dpi=resolve)
plt.close('all')


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
fig.savefig(os.path.join(PLOT_PATH,'{0}-{1}-phase-diagram_microbe-ext-occurrences.pdf'.format(xvar,yvar)),dpi=resolve)
plt.close('all')

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
fig.savefig(os.path.join(PLOT_PATH,'{0}-{1}-phase-diagram_expected-num-walls.pdf'.format(xvar,yvar)),dpi=resolve)
plt.close('all')

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
fig.savefig(os.path.join(PLOT_PATH,'{0}-{1}-phase-diagram_expected-wall-durations.pdf'.format(xvar,yvar)),dpi=resolve)
plt.close('all')

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
fig.savefig(os.path.join(PLOT_PATH,'{0}-{1}-phase-diagram_virus-sim-end-times.pdf'.format(xvar,yvar)),dpi=resolve)
plt.close('all')

fig = plt.figure(figsize=(10,10))
ax = fig.add_subplot(111)
im = ax.imshow(grid_bEndTimes.T, extent=(minX,maxX,minY,maxY), aspect='auto', origin='lower',cmap = reversed_map)
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="2%", pad=0.2)
cbar = fig.colorbar(im, cax)
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
fig.savefig(os.path.join(PLOT_PATH,'{0}-{1}-phase-diagram_microbe-sim-end-times.pdf'.format(xvar,yvar)),dpi=resolve)
plt.close('all')


########################################################################
#### spacer acqusition probability q vs. viral repertoire size g #######
########################################################################
DBWALLS_PATH = os.path.join('/Volumes/Yadgah/g-sweep/17-5-2023/mean-walls-shannon.sqlite')
DBEXTMEAN_PATH = os.path.join('/Volumes/Yadgah/g-sweep/17-5-2023/mean-extinctions.sqlite')
DBEXT_PATH = os.path.join('/Volumes/Yadgah/g-sweep/17-5-2023/extinctions.sqlite')
conWalls = sqlite3.connect(DBWALLS_PATH)
curWalls = conWalls.cursor()
conMExt = sqlite3.connect(DBEXTMEAN_PATH)
curMExt = conExt.cursor()
conExt = sqlite3.connect(DBEXT_PATH)
curExt = conExt.cursor()


# make this more general for other parameters!!!!
print('SQLite Query: parameter space spanning combinations')
fixed = "viral_mutation_rate"
fixedValue = 1.04e-06
# fixedValue = 7.8e-07
# fixedValue = 4.0e-08
#
# xaxis = 'viral_mutation_rate'
# xvar = 'mu'
# xaxlabel = r'$\mu$'
#
xaxis = 'spacer_acquisition_prob'
xvar = 'q'
xaxlabel = 'q'
#
yaxis = 'n_protospacers'
yvar = 'g'
yaxlabel = 'g'
resolve = 500
####
####
comboSpace = pd.read_sql_query(
"SELECT combo_id, {0}, {1} \
FROM param_combos WHERE {2} = {3} ORDER BY combo_id".format(xaxis,yaxis,fixed,fixedValue), conWalls)
comboSpace[yaxis] = comboSpace.n_protospacers
# comboSpace[yaxis] = 1 - (1 - comboSpace.viral_mutation_rate)**fixedValue
###
wallOccurrences = pd.read_sql_query("SELECT wall_occurrence, combo_id FROM microbial_wall_statistics \
                                    WHERE combo_id in ({0}) ORDER BY combo_id".format(', '.join(map(str, comboSpace['combo_id']))), conWalls)
numWalls = pd.read_sql_query("SELECT expected_num_walls, combo_id FROM microbial_wall_statistics \
                             WHERE combo_id in ({0}) ORDER BY combo_id".format(', '.join(map(str, comboSpace['combo_id']))), conWalls)
meanWallDurations = pd.read_sql_query("SELECT num_walls, mean_wallduration, num_replicates, total_replicates, combo_id \
                                        FROM microbial_mean_walldurations \
                                        WHERE combo_id in ({0}) \
                                        ORDER BY combo_id".format(', '.join(map(str, comboSpace['combo_id']))), conWalls)
meanWallDurations['mean_wallduration'] = meanWallDurations['mean_wallduration'] * (meanWallDurations['num_replicates']/meanWallDurations['total_replicates'])
meanWallDurations = meanWallDurations.groupby(['combo_id']).agg(mean_wallduration=('mean_wallduration', 'sum')).reset_index()
#######
#######
extOccurrences = pd.read_sql_query("SELECT microbes,viruses,run_id FROM extinction_occurrence \
                                   ORDER BY run_id", conExt)\
                   .merge(
    pd.read_sql_query(
        "SELECT run_id,combo_id FROM runs \
            WHERE combo_id in ({0})".format(', '.join(map(str, comboSpace['combo_id'])))\
            , conExt), on=['run_id']
                        )
extOccurrences = extOccurrences.groupby(['combo_id']).agg(
    num_replicates=('run_id', 'size')).reset_index()\
    .merge(extOccurrences,on=['combo_id'])
extOccurrences = extOccurrences[extOccurrences.microbes==0]
extOccurrences = extOccurrences.groupby(['combo_id','num_replicates']).agg(
    viruses=('viruses', 'sum')).reset_index()
extOccurrences['viruses'] = extOccurrences['viruses']/extOccurrences['num_replicates']


extOccurrencesM = pd.read_sql_query("SELECT microbes,run_id FROM extinction_occurrence \
                                    ORDER BY run_id", conExt)\
                    .merge(
    pd.read_sql_query(
        "SELECT run_id,combo_id FROM runs \
            WHERE combo_id in ({0})".format(', '.join(map(str, comboSpace['combo_id'])))\
            , conExt), on=['run_id']
                        )
extOccurrencesM = extOccurrencesM.groupby(['combo_id']).agg(
    num_replicates=('run_id', 'size')).reset_index()\
    .merge(extOccurrencesM,on=['combo_id'])
extOccurrencesM = extOccurrencesM.groupby(['combo_id','num_replicates']).agg(
    microbes=('microbes', 'sum')).reset_index()
extOccurrencesM['microbes'] = extOccurrencesM['microbes']/extOccurrencesM['num_replicates']
extOccurrences = extOccurrences.drop(columns=['num_replicates']).merge(extOccurrencesM,on=['combo_id'])
##
simEndTimes = pd.read_sql_query("SELECT microbe_end_time,virus_end_time,run_id \
                                FROM simulation_end_time \
                                ORDER BY run_id", conExt)\
    .merge(
    pd.read_sql_query(
        "SELECT run_id,combo_id FROM runs \
            WHERE combo_id in ({0})".format(', '.join(map(str, comboSpace['combo_id']))), conExt), on=['run_id']
)
simEndTimes = simEndTimes.groupby(['combo_id']).agg(
    microbe_end_time=('microbe_end_time', 'mean'), virus_end_time=('virus_end_time', 'mean')).reset_index()



# externa hard drive
DBWALLS_PATH = os.path.join('/Volumes/Yadgah/g-sweep/25-5-2023/mean-walls-shannon.sqlite')
DBEXTMEAN_PATH = os.path.join('/Volumes/Yadgah/g-sweep/25-5-2023/mean-extinctions.sqlite')
DBEXT_PATH = os.path.join('/Volumes/Yadgah/g-sweep/25-5-2023/extinctions.sqlite')
conWalls = sqlite3.connect(DBWALLS_PATH)
curWalls = conWalls.cursor()
conMExt = sqlite3.connect(DBEXTMEAN_PATH)
curMExt = conExt.cursor()
conExt = sqlite3.connect(DBEXT_PATH)
curExt = conExt.cursor()


comboSpace2 = pd.read_sql_query(
"SELECT combo_id, {0}, {1} \
FROM param_combos WHERE {2} = {3} ORDER BY combo_id".format(xaxis,yaxis,fixed,fixedValue), conExt)
comboSpace2[yaxis] = comboSpace2.n_protospacers


wallOccurrences2 = pd.read_sql_query("SELECT wall_occurrence, combo_id FROM microbial_wall_statistics \
                                    WHERE combo_id in ({0}) ORDER BY combo_id".format(', '.join(map(str, comboSpace2['combo_id']))), conWalls)
numWalls2 = pd.read_sql_query("SELECT expected_num_walls, combo_id FROM microbial_wall_statistics \
                             WHERE combo_id in ({0}) ORDER BY combo_id".format(', '.join(map(str, comboSpace2['combo_id']))), conWalls)
meanWallDurations2 = pd.read_sql_query("SELECT num_walls, mean_wallduration, num_replicates, total_replicates, combo_id \
                                        FROM microbial_mean_walldurations \
                                        WHERE combo_id in ({0}) \
                                        ORDER BY combo_id".format(', '.join(map(str, comboSpace2['combo_id']))), conWalls)
meanWallDurations2['mean_wallduration'] = meanWallDurations2['mean_wallduration'] * \
    (meanWallDurations2['num_replicates'] /
     meanWallDurations2['total_replicates'])
meanWallDurations2 = meanWallDurations2.groupby(['combo_id']).agg(
    mean_wallduration=('mean_wallduration', 'sum')).reset_index()
#######
#######
extOccurrences2 = pd.read_sql_query("SELECT microbes,viruses,run_id FROM extinction_occurrence \
                                   ORDER BY run_id".format(', '.join(map(str, comboSpace2['combo_id']))), conExt)\
    .merge(
    pd.read_sql_query(
        "SELECT run_id,combo_id FROM runs \
            WHERE combo_id in ({0})".format(', '.join(map(str, comboSpace2['combo_id'])))\
            , conExt), on=['run_id']
)
extOccurrences2 = extOccurrences2.groupby(['combo_id']).agg(
    num_replicates=('run_id', 'size')).reset_index()\
    .merge(extOccurrences2, on=['combo_id'])
extOccurrences2 = extOccurrences2[extOccurrences2.microbes == 0]
extOccurrences2 = extOccurrences2.groupby(['combo_id', 'num_replicates']).agg(
    viruses=('viruses', 'sum')).reset_index()
extOccurrences2['viruses'] = extOccurrences2['viruses'] / \
    extOccurrences2['num_replicates']


extOccurrencesM2 = pd.read_sql_query("SELECT microbes,run_id FROM extinction_occurrence \
                                    ORDER BY run_id", conExt)\
    .merge(
    pd.read_sql_query(
        "SELECT run_id,combo_id FROM runs \
            WHERE combo_id in ({0})".format(', '.join(map(str, comboSpace2['combo_id'])))\
            , conExt), on=['run_id']
)
extOccurrencesM2 = extOccurrencesM2.groupby(['combo_id']).agg(
    num_replicates=('run_id', 'size')).reset_index()\
    .merge(extOccurrencesM2, on=['combo_id'])
extOccurrencesM2 = extOccurrencesM2.groupby(['combo_id', 'num_replicates']).agg(
    microbes=('microbes', 'sum')).reset_index()
extOccurrencesM2['microbes'] = extOccurrencesM2['microbes'] / \
    extOccurrencesM2['num_replicates']
extOccurrences2 = extOccurrences2.drop(
    columns=['num_replicates']).merge(extOccurrencesM2, on=['combo_id'])
##
simEndTimes2 = pd.read_sql_query("SELECT microbe_end_time,virus_end_time,run_id \
                                FROM simulation_end_time \
                                ORDER BY run_id", conExt)\
    .merge(
    pd.read_sql_query(
        "SELECT run_id,combo_id FROM runs \
            WHERE combo_id in ({0})".format(', '.join(map(str, comboSpace2['combo_id']))), conExt), on=['run_id']
)
simEndTimes2 = simEndTimes2.groupby(['combo_id']).agg(
    microbe_end_time=('microbe_end_time', 'mean'), virus_end_time=('virus_end_time', 'mean')).reset_index()



comboSpace['combo_id'] = [*range(1, len(comboSpace['combo_id'])+1, 1)]
wallOccurrences['combo_id'] = [*range(1, len(comboSpace['combo_id'])+1, 1)]
numWalls['combo_id'] = [*range(1, len(comboSpace['combo_id'])+1, 1)]
meanWallDurations['combo_id'] = [*range(1, len(comboSpace['combo_id'])+1, 1)]
extOccurrences['combo_id'] = [*range(1, len(comboSpace['combo_id'])+1, 1)]
simEndTimes['combo_id'] = [*range(1, len(comboSpace['combo_id'])+1, 1)]
#
comboSpace2['combo_id'] = [
    *range(comboSpace['combo_id'].values[-1]+1, \
           comboSpace['combo_id'].values[-1] + 1 + len(comboSpace2['combo_id']), 1)]
wallOccurrences2['combo_id'] = comboSpace2['combo_id'].values
numWalls2['combo_id'] = comboSpace2['combo_id'].values
meanWallDurations2['combo_id'] = comboSpace2['combo_id'].values
extOccurrences2['combo_id'] = comboSpace2['combo_id'].values
simEndTimes2['combo_id'] = comboSpace2['combo_id'].values
#
comboSpace = pd.concat([comboSpace,comboSpace2]).reset_index(drop=True)
wallOccurrences  = pd.concat(
    [wallOccurrences, wallOccurrences2]).reset_index(drop=True)
numWalls = pd.concat(
    [numWalls, numWalls2]).reset_index(drop=True)
meanWallDurations= pd.concat(
    [meanWallDurations, meanWallDurations2]).reset_index(drop=True)
extOccurrences = pd.concat(
    [extOccurrences, extOccurrences2]).reset_index(drop=True)
simEndTimes = pd.concat(
    [simEndTimes, simEndTimes2]).reset_index(drop=True)


minX = min(comboSpace[xaxis])
maxX = max(comboSpace[xaxis])
minY = min(comboSpace[yaxis])
maxY = max(comboSpace[yaxis])
y = np.array(sorted(np.unique(comboSpace['n_protospacers'].values)))
s = len(y)
grid_x, grid_y = np.mgrid[minX:maxX:11j, minY:maxY:150j]
# grid_y = y.copy()
# for i in range(0,len(grid_x)-1):
#     grid_y = np.vstack([grid_y,y])
interpMethod = 'linear'
# interpMethod ='nearest'
grid_wallOccurrences = griddata((comboSpace[xaxis].values,comboSpace[yaxis].values),
wallOccurrences.wall_occurrence.values, (grid_x, grid_y), method=interpMethod, rescale=True)
grid_numWalls = griddata((comboSpace[xaxis].values,comboSpace[yaxis].values),
                         numWalls.expected_num_walls.values, (grid_x, grid_y), method=interpMethod, rescale=True)
grid_expectedDuration = griddata((comboSpace[xaxis].values,comboSpace[yaxis].values),
meanWallDurations.mean_wallduration.values, (grid_x, grid_y), method=interpMethod, rescale=True)

grid_vExt = griddata((comboSpace[xaxis].values,comboSpace[yaxis].values),
extOccurrences.viruses.values, (grid_x, grid_y), method=interpMethod, rescale=True)
grid_bExt = griddata((comboSpace[xaxis].values,comboSpace[yaxis].values),
extOccurrences.microbes.values, (grid_x, grid_y), method=interpMethod, rescale=True)
grid_vEndTimes = griddata((comboSpace[xaxis].values,comboSpace[yaxis].values),
simEndTimes.virus_end_time.values, (grid_x, grid_y), method=interpMethod, rescale=True)
grid_bEndTimes = griddata((comboSpace[xaxis].values,comboSpace[yaxis].values),
simEndTimes.microbe_end_time.values, (grid_x, grid_y), method=interpMethod, rescale=True)


reversed_map = 'jet'
lowerY = minY
upperY = maxY


fig = plt.figure(figsize=(10, 10))
ax = fig.add_subplot(111)
im = ax.imshow(grid_wallOccurrences.T, extent=[
               minX, maxX, minY, maxY], aspect='auto', origin='lower', cmap=reversed_map)
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="2%", pad=0.2)
cbar = fig.colorbar(im, cax)
cbar.ax.tick_params(labelsize=20)
ax.set_xlabel(xaxlabel, labelpad=20, fontsize=25)
ax.set_ylabel(yaxlabel, labelpad=20, fontsize=25)
yticks = sorted(np.unique(comboSpace['n_protospacers'].values))
ax.set_yticks(yticks)
ax.set_yticklabels(['{}'.format(int(i)) for i in yticks])
ax.tick_params(axis='x', labelsize=20)
ax.tick_params(axis='y', labelsize=20)
yticks = ax.get_yticklabels()
yticks[0].set_fontsize(7)
yticks[1].set_fontsize(7)
yticks[2].set_fontsize(15)
text = ax.xaxis.get_offset_text()  # Get the text object
text.set_size(20)  # Set the size
ax.set_title('Probability of Alternating Dynamics', pad=30, fontsize=20)
fig.savefig(os.path.join(PLOT_PATH,'{0}-{1}-phase-diagram_wall-occurrences_{2}_{3}.pdf'.format(
    xvar, yvar, fixed, fixedValue)), dpi=resolve)
plt.close('all')


fig = plt.figure(figsize=(10, 10))
ax = fig.add_subplot(111)
im = ax.imshow(grid_vExt.T, extent=(minX, maxX, minY, maxY),
               aspect='auto', origin='lower', cmap=reversed_map)
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="2%", pad=0.2)
cbar = fig.colorbar(im, cax)
cbar.ax.tick_params(labelsize=20)
ax.set_xlabel(xaxlabel, labelpad=20, fontsize=25)
ax.set_ylabel(yaxlabel, labelpad=20, fontsize=25)
yticks = sorted(np.unique(comboSpace['n_protospacers'].values))
ax.set_yticks(yticks)
ax.set_yticklabels(['{}'.format(int(i)) for i in yticks])
ax.tick_params(axis='x', labelsize=20)
ax.tick_params(axis='y', labelsize=20)
yticks = ax.get_yticklabels()
yticks[0].set_fontsize(7)
yticks[1].set_fontsize(7)
yticks[2].set_fontsize(15)
text = ax.xaxis.get_offset_text()  # Get the text object
text.set_size(20)  # Set the size
text = ax.yaxis.get_offset_text()  # Get the text object
text.set_size(20)  # Set the size
ax.set_title(
    r'Probability of Total Viral Extinction Before $t = 2000$', pad=30, fontsize=20)
fig.savefig(os.path.join(PLOT_PATH,'{0}-{1}-phase-diagram_virus-ext-occurrences_{2}_{3}.pdf'.format(
    xvar, yvar, fixed, fixedValue)), dpi=resolve)
plt.close('all')


fig = plt.figure(figsize=(10, 10))
ax = fig.add_subplot(111)
im = ax.imshow(grid_bExt.T, extent=(minX, maxX, minY, maxY),
               aspect='auto', origin='lower', cmap=reversed_map)
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="2%", pad=0.2)
cbar = fig.colorbar(im, cax, ticks=[0])
cbar.ax.tick_params(labelsize=20)
ax.set_xlabel(xaxlabel, labelpad=20, fontsize=25)
ax.set_ylabel(yaxlabel, labelpad=20, fontsize=25)
yticks = sorted(np.unique(comboSpace['n_protospacers'].values))
ax.set_yticks(yticks)
ax.set_yticklabels(['{}'.format(int(i)) for i in yticks])
ax.tick_params(axis='x', labelsize=20)
ax.tick_params(axis='y', labelsize=20)
yticks = ax.get_yticklabels()
yticks[0].set_fontsize(7)
yticks[1].set_fontsize(7)
yticks[2].set_fontsize(15)
text = ax.xaxis.get_offset_text()  # Get the text object
text.set_size(20)  # Set the size
text = ax.yaxis.get_offset_text()  # Get the text object
text.set_size(20)  # Set the size
ax.set_title(
    r'Probability of Total Microbial Extinction Before $t = 2000$', pad=30, fontsize=20)
fig.savefig(os.path.join(PLOT_PATH,'{0}-{1}-phase-diagram_microbe-ext-occurrences_{2}_{3}.pdf'.format(
    xvar, yvar, fixed, fixedValue)), dpi=resolve)
plt.close('all')

fig = plt.figure(figsize=(10, 10))
ax = fig.add_subplot(111)
im = ax.imshow(grid_numWalls.T, extent=[
               minX, maxX, minY, maxY], aspect='auto', origin='lower', cmap=reversed_map)
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="2%", pad=0.2)
cbar = fig.colorbar(im, cax)
cbar.ax.tick_params(labelsize=20)
ax.set_xlabel(xaxlabel, labelpad=20, fontsize=25)
ax.set_ylabel(yaxlabel, labelpad=20, fontsize=25)
yticks = sorted(np.unique(comboSpace['n_protospacers'].values))
ax.set_yticks(yticks)
ax.set_yticklabels(['{}'.format(int(i)) for i in yticks])
ax.tick_params(axis='x', labelsize=20)
ax.tick_params(axis='y', labelsize=20)
yticks = ax.get_yticklabels()
yticks[0].set_fontsize(7)
yticks[1].set_fontsize(7)
yticks[2].set_fontsize(15)
text = ax.xaxis.get_offset_text()  # Get the text object
text.set_size(20)  # Set the size
text = ax.yaxis.get_offset_text()  # Get the text object
text.set_size(20)  # Set the size
ax.set_title('Expected Number of SHC Periods', pad=30, fontsize=20)
fig.savefig(os.path.join(PLOT_PATH,'{0}-{1}-phase-diagram_expected-num-walls_{2}_{3}.pdf'.format(
    xvar, yvar, fixed, fixedValue)), dpi=resolve)
plt.close('all')


fig = plt.figure(figsize=(10, 10))
ax = fig.add_subplot(111)
im = ax.imshow(grid_expectedDuration.T, extent=(
    minX, maxX, minY, maxY), aspect='auto', origin='lower', cmap=reversed_map)
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="2%", pad=0.2)
cbar = fig.colorbar(im, cax)
cbar.ax.tick_params(labelsize=20)
ax.set_xlabel(xaxlabel, labelpad=20, fontsize=25)
ax.set_ylabel(yaxlabel, labelpad=20, fontsize=25)
yticks = sorted(np.unique(comboSpace['n_protospacers'].values))
ax.set_yticks(yticks)
ax.set_yticklabels(['{}'.format(int(i)) for i in yticks])
ax.tick_params(axis='x', labelsize=20)
ax.tick_params(axis='y', labelsize=20)
yticks = ax.get_yticklabels()
yticks[0].set_fontsize(7)
yticks[1].set_fontsize(7)
yticks[2].set_fontsize(15)
text = ax.xaxis.get_offset_text()  # Get the text object
text.set_size(20)  # Set the size
text = ax.yaxis.get_offset_text()  # Get the text object
text.set_size(20)  # Set the size
# ax.set_ylim(lowerY, upperY)
ax.set_title('Expected Duration of SHC Periods', pad=30, fontsize=20)
fig.savefig(os.path.join(PLOT_PATH,'{0}-{1}-phase-diagram_expected-wall-durations_{2}_{3}.pdf'.format(
    xvar, yvar, fixed, fixedValue)), dpi=resolve)
plt.close('all')



fig = plt.figure(figsize=(10, 10))
ax = fig.add_subplot(111)
im = ax.imshow(grid_vEndTimes.T, extent=(minX, maxX, minY, maxY),
               aspect='auto', origin='lower', cmap=reversed_map)
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="2%", pad=0.2)
cbar = fig.colorbar(im, cax)
cbarticklabels = cbar.ax.get_yticks()
cbarticklabels = ['{}'.format(int(i)) for i in cbarticklabels]
cbarticklabels[-1] = r'$(2000,\infty)$'
cbar.ax.set_yticklabels(cbarticklabels)
cbar.ax.tick_params(labelsize=20)
cbarticks = cbar.ax.get_yticklabels()
cbarticks[-1].set_fontsize(13)
ax.set_xlabel(xaxlabel, labelpad=20, fontsize=25)
ax.set_ylabel(yaxlabel, labelpad=20, fontsize=25)
yticks = sorted(np.unique(comboSpace['n_protospacers'].values))
ax.set_yticks(yticks)
ax.set_yticklabels(['{}'.format(int(i)) for i in yticks])
ax.tick_params(axis='x', labelsize=20)
ax.tick_params(axis='y', labelsize=20)
yticks = ax.get_yticklabels()
yticks[0].set_fontsize(7)
yticks[1].set_fontsize(7)
yticks[2].set_fontsize(15)
text = ax.xaxis.get_offset_text()  # Get the text object
text.set_size(20)  # Set the size
text = ax.yaxis.get_offset_text()  # Get the text object
text.set_size(20)  # Set the size
# ax.set_ylim(lowerY, upperY)
ax.set_title('Mean Time to Total Viral Extinction', pad=30, fontsize=20)
fig.savefig(os.path.join(PLOT_PATH,'{0}-{1}-phase-diagram_virus-sim-end-times._{2}_{3}.pdf'.format(
    xvar, yvar, fixed, fixedValue)), dpi=resolve)
plt.close('all')

fig = plt.figure(figsize=(10, 10))
ax = fig.add_subplot(111)
im = ax.imshow(grid_bEndTimes.T, extent=(minX, maxX, minY, maxY),
               aspect='auto', origin='lower', cmap=reversed_map)
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
yticks = sorted(np.unique(comboSpace['n_protospacers'].values))
ax.set_yticks(yticks)
ax.set_yticklabels(['{}'.format(int(i)) for i in yticks])
ax.tick_params(axis='x', labelsize=20)
ax.tick_params(axis='y', labelsize=20)
yticks = ax.get_yticklabels()
yticks[0].set_fontsize(7)
yticks[1].set_fontsize(7)
yticks[2].set_fontsize(15)
text = ax.xaxis.get_offset_text()  # Get the text object
text.set_size(20)  # Set the size
text = ax.yaxis.get_offset_text()  # Get the text object
text.set_size(20)  # Set the size
# ax.set_ylim(lowerY, upperY)
ax.set_title('Mean Time to Total Microbial Extinction', pad=30, fontsize=20)
fig.savefig(os.path.join(PLOT_PATH,'{0}-{1}-phase-diagram_microbe-sim-end-times_{2}_{3}.pdf'.format(
    xvar, yvar, fixed, fixedValue)), dpi=resolve)
plt.close('all')
