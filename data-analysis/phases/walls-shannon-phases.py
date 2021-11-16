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

analysisType = "walls-shannon"

SCRIPT_PATH = os.path.abspath(os.path.dirname(__file__)) # cluster

DBSIM_PATH = os.path.join(SCRIPT_PATH,'..','..','simulation','sweep_db_gathered.sqlite') # cluster
#home_dir = os.system("cd ~") # local
#DBSIM_PATH = os.path.join('/Volumes','Yadgah','sweep_db_gathered.sqlite') # local

DBRICH_PATH = os.path.join(SCRIPT_PATH,'..','gathered-analyses','richness','richness.sqlite') # cluster
#DBRICH_PATH = os.path.join('/Volumes','Yadgah','richness.sqlite') # local

DBSHAN_PATH = os.path.join(SCRIPT_PATH,'..','gathered-analyses','shannon','shannon.sqlite') # cluster
#DBSHAN_PATH = os.path.join('/Volumes','Yadgah','shannon.sqlite') # local

DBWALLS_PATH = os.path.join(SCRIPT_PATH,'..','gathered-analyses','{}','mean-{}.sqlite'.format(analysisType)) # cluster
#DBWALLS_PATH = os.path.join('/Volumes','Yadgah','mean-{}.sqlite'.format(analysisType)) # local
#DBWALLS_PATH = os.path.join('/Volumes','Yadgah','mean-{}_output.sqlite'.format(analysisType)) # local. run_id fixed; for testing

DBEXT_PATH = os.path.join(SCRIPT_PATH,'..','gathered-analyses','{}','mean-{}.sqlite'.format('extinctions')) # cluster
#DBEXT_PATH = os.path.join('/Volumes','Yadgah','mean-{}.sqlite'.format('extinctions')) # local
#DBEXT_PATH = os.path.join('/Volumes','Yadgah','mean-{}_output.sqlite'.format('extinctions')) # local. run_id fixed; for testing

conWalls = sqlite3.connect(DBWALLS_PATH)
curWalls = conWalls.cursor()
print('SQLite Query: mean microbial wall occurences per parameter combination')
wallOccurrences = pd.read_sql_query("SELECT wall_occurrence FROM microbial_wall_statistics ORDER BY combo_id", conWalls)
# make this more general for other parameters!!!!
print('SQLite Query: parameter space spanning combinations')
comboSpace = pd.read_sql_query(
"SELECT combo_id, crispr_failure_prob, spacer_acquisition_prob, viral_mutation_rate \
FROM param_combos ORDER BY combo_id", conWalls)
comboSpaceTrunc = comboSpace[(comboSpace.combo_id < 2)|(comboSpace.combo_id >30)]

orig_map=plt.cm.get_cmap('viridis')
# reversing the original colormap using reversed() function
reversed_map = orig_map.reversed()
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
phases = ax.scatter(comboSpaceTrunc['crispr_failure_prob'], comboSpaceTrunc['spacer_acquisition_prob'],
comboSpaceTrunc['viral_mutation_rate'], c=wallOccurrences['wall_occurrence'],
lw=0, s=20, cmap = reversed_map)
#ax.ticklabel_format(style='sci',scilimits=(0,0))
#ax.xaxis.set_major_formatter(ticker.ScalarFormatter(useMathText=True))
fig.canvas.draw()
offsetx = ax.xaxis.get_offset_text().get_text()
offsety = ax.yaxis.get_offset_text().get_text()
offsetz = ax.zaxis.get_offset_text().get_text()
#ax.xaxis.set_major_formatter()
#ax.yaxis.set_major_formatter()
#ax.zaxis.set_major_formatter(ticker.ScalarFormatter(useOffset=False))
ax.set_xlabel(r'CRISPR Failure Probability $/$ {}'.format(offsetx),labelpad=2,fontsize=7)
ax.set_ylabel(r'Spacer Acquisition Probability $/$ {}'.format(offsety),labelpad=6,fontsize=7)
ax.set_zlabel(r'Viral Mutation Rate $/$ {}'.format(offsetz),labelpad=6,fontsize=7)
ax.tick_params(labelsize=5,pad=0)
ax.minorticks_off()
ax.xaxis.offsetText.set_visible(False)
ax.yaxis.offsetText.set_visible(False)
ax.zaxis.offsetText.set_visible(False)
ax.set_title("Wall Occurrences")
fig.colorbar(phases, ax=ax, shrink = 0.3, location ='bottom',pad=0.1)
elev = 15
azim = 170
ax.view_init(elev, azim)
fig.savefig(os.path.join('/Volumes/Yadgah','phasediagram_wall-occurrences-2.png'),dpi=500)
elev = 0
azim = 170
ax.view_init(elev, azim)
ax.xaxis.set_ticks([])
fig.savefig(os.path.join('/Volumes/Yadgah','phasediagram_wall-occurrences-1.png'),dpi=500)
#plt.show()


expectedNumWalls = pd.read_sql_query("SELECT expected_num_walls FROM microbial_wall_statistics ORDER BY combo_id", conWalls)
orig_map=plt.cm.get_cmap('viridis')
# reversing the original colormap using reversed() function
reversed_map = orig_map.reversed()
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
phases = ax.scatter(comboSpaceTrunc['crispr_failure_prob'], comboSpaceTrunc['spacer_acquisition_prob'],
comboSpaceTrunc['viral_mutation_rate'], c=expectedNumWalls['expected_num_walls'],
lw=0, s=20, cmap = reversed_map)
#ax.ticklabel_format(style='sci',scilimits=(0,0))
#ax.xaxis.set_major_formatter(ticker.ScalarFormatter(useMathText=True))
fig.canvas.draw()
offsetx = ax.xaxis.get_offset_text().get_text()
offsety = ax.yaxis.get_offset_text().get_text()
offsetz = ax.zaxis.get_offset_text().get_text()
#ax.xaxis.set_major_formatter()
#ax.yaxis.set_major_formatter()
#ax.zaxis.set_major_formatter(ticker.ScalarFormatter(useOffset=False))
ax.set_xlabel(r'CRISPR Failure Probability $/$ {}'.format(offsetx),labelpad=2,fontsize=7)
ax.set_ylabel(r'Spacer Acquisition Probability $/$ {}'.format(offsety),labelpad=6,fontsize=7)
ax.set_zlabel(r'Viral Mutation Rate $/$ {}'.format(offsetz),labelpad=6,fontsize=7)
ax.tick_params(labelsize=5,pad=0)
ax.minorticks_off()
ax.xaxis.offsetText.set_visible(False)
ax.yaxis.offsetText.set_visible(False)
ax.zaxis.offsetText.set_visible(False)
ax.set_title("Expected Number of Walls")
fig.colorbar(phases, ax=ax, shrink = 0.3, location ='bottom',pad=0.1)
elev = 15
azim = 170
ax.view_init(elev, azim)
fig.savefig(os.path.join('/Volumes/Yadgah','phasediagram_expected-num-walls-2.png'),dpi=500)
elev = 0
azim = 170
ax.view_init(elev, azim)
ax.xaxis.set_ticks([])
fig.savefig(os.path.join('/Volumes/Yadgah','phasediagram_expected-num-walls-1.png'),dpi=500)



expectedNumWalls = pd.read_sql_query("SELECT expected_num_walls_5percentThreshold FROM microbial_wall_statistics ORDER BY combo_id", conWalls)
orig_map=plt.cm.get_cmap('viridis')
# reversing the original colormap using reversed() function
reversed_map = orig_map.reversed()
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
phases = ax.scatter(comboSpaceTrunc['crispr_failure_prob'], comboSpaceTrunc['spacer_acquisition_prob'],
comboSpaceTrunc['viral_mutation_rate'], c=expectedNumWalls['expected_num_walls_5percentThreshold'],
lw=0, s=20, cmap = reversed_map)
#ax.ticklabel_format(style='sci',scilimits=(0,0))
#ax.xaxis.set_major_formatter(ticker.ScalarFormatter(useMathText=True))
fig.canvas.draw()
offsetx = ax.xaxis.get_offset_text().get_text()
offsety = ax.yaxis.get_offset_text().get_text()
offsetz = ax.zaxis.get_offset_text().get_text()
#ax.xaxis.set_major_formatter()
#ax.yaxis.set_major_formatter()
#ax.zaxis.set_major_formatter(ticker.ScalarFormatter(useOffset=False))
ax.set_xlabel(r'CRISPR Failure Probability $/$ {}'.format(offsetx),labelpad=2,fontsize=7)
ax.set_ylabel(r'Spacer Acquisition Probability $/$ {}'.format(offsety),labelpad=6,fontsize=7)
ax.set_zlabel(r'Viral Mutation Rate $/$ {}'.format(offsetz),labelpad=6,fontsize=7)
ax.tick_params(labelsize=5,pad=0)
ax.minorticks_off()
ax.xaxis.offsetText.set_visible(False)
ax.yaxis.offsetText.set_visible(False)
ax.zaxis.offsetText.set_visible(False)
ax.set_title(r'Expected Number of Walls (5$\%$ threshold)')
fig.colorbar(phases, ax=ax, shrink = 0.3, location ='bottom',pad=0.1)
elev = 15
azim = 170
ax.view_init(elev, azim)
fig.savefig(os.path.join('/Volumes/Yadgah','phasediagram_expected-num-walls-5percent-2.png'),dpi=500)
elev = 0
azim = 170
ax.view_init(elev, azim)
ax.xaxis.set_ticks([])
fig.savefig(os.path.join('/Volumes/Yadgah','phasediagram_expected-num-walls-5percent-1.png'),dpi=500)

wallDurations = pd.read_sql_query("SELECT mean_wallduration_of_avg_most_frequents FROM microbial_wall_statistics ORDER BY combo_id", conWalls)
orig_map=plt.cm.get_cmap('viridis')
# reversing the original colormap using reversed() function
reversed_map = orig_map.reversed()
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
phases = ax.scatter(comboSpaceTrunc['crispr_failure_prob'], comboSpaceTrunc['spacer_acquisition_prob'],
comboSpaceTrunc['viral_mutation_rate'], c=wallDurations['mean_wallduration_of_avg_most_frequents'],
lw=0, s=20, cmap = reversed_map)
#ax.ticklabel_format(style='sci',scilimits=(0,0))
#ax.xaxis.set_major_formatter(ticker.ScalarFormatter(useMathText=True))
fig.canvas.draw()
offsetx = ax.xaxis.get_offset_text().get_text()
offsety = ax.yaxis.get_offset_text().get_text()
offsetz = ax.zaxis.get_offset_text().get_text()
#ax.xaxis.set_major_formatter()
#ax.yaxis.set_major_formatter()
#ax.zaxis.set_major_formatter(ticker.ScalarFormatter(useOffset=False))
ax.set_xlabel(r'CRISPR Failure Probability $/$ {}'.format(offsetx),labelpad=2,fontsize=7)
ax.set_ylabel(r'Spacer Acquisition Probability $/$ {}'.format(offsety),labelpad=6,fontsize=7)
ax.set_zlabel(r'Viral Mutation Rate $/$ {}'.format(offsetz),labelpad=6,fontsize=7)
ax.tick_params(labelsize=5,pad=0)
ax.minorticks_off()
ax.xaxis.offsetText.set_visible(False)
ax.yaxis.offsetText.set_visible(False)
ax.zaxis.offsetText.set_visible(False)
ax.set_title('Mean Duration of Longest Wall (of Most Frequent No. of Walls)')
fig.colorbar(phases, ax=ax, shrink = 0.3, location ='bottom',pad=0.1)
elev = 15
azim = 170
ax.view_init(elev, azim)
fig.savefig(os.path.join('/Volumes/Yadgah','phasediagram_mean-wallduration-most-frequent-2.png'),dpi=500)
elev = 0
azim = 170
ax.view_init(elev, azim)
ax.xaxis.set_ticks([])
fig.savefig(os.path.join('/Volumes/Yadgah','phasediagram_mean-wallduration-most-frequent-1.png'),dpi=500)

mostFrequent = pd.read_sql_query("SELECT avg_of_most_frequent_num_walls FROM microbial_wall_statistics ORDER BY combo_id", conWalls)
orig_map=plt.cm.get_cmap('viridis')
# reversing the original colormap using reversed() function
reversed_map = orig_map.reversed()
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
phases = ax.scatter(comboSpaceTrunc['crispr_failure_prob'], comboSpaceTrunc['spacer_acquisition_prob'],
comboSpaceTrunc['viral_mutation_rate'], c=mostFrequent['avg_of_most_frequent_num_walls'],
lw=0, s=20, cmap = reversed_map)
#ax.ticklabel_format(style='sci',scilimits=(0,0))
#ax.xaxis.set_major_formatter(ticker.ScalarFormatter(useMathText=True))
fig.canvas.draw()
offsetx = ax.xaxis.get_offset_text().get_text()
offsety = ax.yaxis.get_offset_text().get_text()
offsetz = ax.zaxis.get_offset_text().get_text()
#ax.xaxis.set_major_formatter()
#ax.yaxis.set_major_formatter()
#ax.zaxis.set_major_formatter(ticker.ScalarFormatter(useOffset=False))
ax.set_xlabel(r'CRISPR Failure Probability $/$ {}'.format(offsetx),labelpad=2,fontsize=7)
ax.set_ylabel(r'Spacer Acquisition Probability $/$ {}'.format(offsety),labelpad=6,fontsize=7)
ax.set_zlabel(r'Viral Mutation Rate $/$ {}'.format(offsetz),labelpad=6,fontsize=7)
ax.tick_params(labelsize=5,pad=0)
ax.minorticks_off()
ax.xaxis.offsetText.set_visible(False)
ax.yaxis.offsetText.set_visible(False)
ax.zaxis.offsetText.set_visible(False)
ax.set_title('Most Frequent Number of Walls')
fig.colorbar(phases, ax=ax, shrink = 0.3, location ='bottom',pad=0.1)
elev = 15
azim = 170
ax.view_init(elev, azim)
fig.savefig(os.path.join('/Volumes/Yadgah','phasediagram_most-frequent-num-walls-2.png'),dpi=500)
elev = 0
azim = 170
ax.view_init(elev, azim)
ax.xaxis.set_ticks([])
fig.savefig(os.path.join('/Volumes/Yadgah','phasediagram_most-frequent-num-walls-1.png'),dpi=500)
