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

run_id = sys.argv[1]

upperThreshold = sys.argv[2]
lowerThreshold = sys.argv[3]

#upperThreshold = 3e5; # for testing
#lowerThreshold = 1.2e5; # for testing

SCRIPT_PATH = os.path.abspath(os.path.dirname(__file__)) # cluster

DBSIM_PATH = os.path.join(SCRIPT_PATH,'..','..','simulation','sweep_db_gathered.sqlite') # cluster
#home_dir = os.system("cd ~") # local
#DBSIM_PATH = os.path.join('/Volumes','Yadgah','sweep_db_gathered.sqlite') # local

conSim = sqlite3.connect(DBSIM_PATH)
curSim = conSim.cursor()
print('SQLite Query: microbial abundance time series data')
microbeSim = pd.read_sql_query("SELECT t,bstrain_id,abundance FROM babundance WHERE run_id = {}".format(run_id), conSim)
microbe_stacked = microbeSim.pivot(index='t',columns='bstrain_id',values='abundance')
print('SQLite Query: viral abundance time series data')
virusSim = pd.read_sql_query("SELECT t,vstrain_id,abundance FROM vabundance WHERE run_id = {}".format(run_id), conSim)
virus_stacked = virusSim.pivot(index='t',columns='vstrain_id',values='abundance')

# Designating plot path from simulation data
ID = curSim.execute('SELECT combo_id,replicate FROM runs WHERE run_id = {}'.format(run_id)).fetchall()
combo_id = ID[0][0]
replicate = ID[0][1]
PLOT_PATH = os.path.join(SCRIPT_PATH,'..','gathered-analyses','walls','plots','c{}'.format(combo_id),'r{}'.format(replicate)) # cluster
#PLOT_PATH = os.path.abspath(os.path.dirname(__file__)) # local
#PLOT_PATH = os.getcwd() # local. run_id fixed; for testing


DBRICH_PATH = os.path.join(SCRIPT_PATH,'..','gathered-analyses','richness','richness.sqlite') # cluster
#DBRICH_PATH = os.path.join('/Volumes','Yadgah','richness.sqlite') # local

conRich = sqlite3.connect(DBRICH_PATH)
curRich = conRich.cursor()
print('SQLite Query: virus strain richness data')
virusRichness = pd.read_sql_query("SELECT t, vrichness FROM richness WHERE run_id = {}".format(run_id), conRich)
print('SQLite Query: microbe immune richness data')
microbeRichness = pd.read_sql_query("SELECT t,brichness FROM richness WHERE run_id = {}".format(run_id), conRich)


DBWALLS_PATH = os.path.join(SCRIPT_PATH,'..','gathered-analyses','walls','walls.sqlite') # cluster
#DBWALLS_PATH = os.path.join('/Volumes','Yadgah','walls.sqlite') # local
#DBWALLS_PATH = os.path.join('/Volumes','Yadgah','walls_output.sqlite') # local. run_id fixed; for testing
conWalls = sqlite3.connect(DBWALLS_PATH)
curWalls = conWalls.cursor()

#print('SQLite Query: microbial walls data')
#microbeWalls = pd.read_sql_query("SELECT t,wall_presence FROM microbial_wall_series", conWalls) # run_id fixed; for testing
#numWalls = pd.read_sql_query("SELECT num_walls FROM microbial_num_walls", conWalls) # run_id fixed; for testing
#wallDurations = pd.read_sql_query("SELECT begin_row,end_row,begin_t,end_t,wall_number,duration FROM microbial_wall_durations", conWalls) # run_id fixed; for testing

print('SQLite Query: microbial walls data')
microbeWalls = pd.read_sql_query("SELECT t,wall_presence FROM microbial_wall_series WHERE run_id = {}".format(run_id), conWalls)
numWalls = pd.read_sql_query("SELECT num_walls FROM microbial_num_walls WHERE run_id = {}".format(run_id), conWalls)
wallDurations = pd.read_sql_query("SELECT begin_row,end_row,begin_t,end_t,wall_number,duration FROM microbial_wall_durations WHERE run_id = {}".format(run_id), conWalls)


print('Compiling microbial wall plot')
pal = sns.color_palette("tab20b")
fig, ax = plt.subplots()
#fig = plt.figure()
fig.suptitle('Microbial Immune Abundances (run{0}-c{1}-r{2})'.format(run_id,combo_id,replicate))

axes = [ax, ax.twinx()]
#axes[1].set_frame_on(False)
#axes[1].patch.set_visible(False)
#axes[1].set_xticks([])

microbe_stacked.plot.area(ax = axes[0],stacked=True,legend=False, linewidth=0,color=pal)
axes[0].set_ylabel(ylabel ='N_i',labelpad=15,fontsize=7)
axes[0].set_xlabel(xlabel = 'Time t',fontsize=7)
axes[0].ticklabel_format(style='sci',scilimits=(0,0))
axes[0].axhline(upperThreshold, color='k', lw=.5, linestyle='dotted')
axes[0].axhline(lowerThreshold, color='k', lw=.5, linestyle='dotted')
axes[1].set_yticks([])
axes[1].fill_between(microbeWalls["t"],0,1, where=microbeWalls["wall_presence"]>0,color="grey", alpha=0.5, linewidth=0,transform=axes[0].get_xaxis_transform())
fig.tight_layout()
fig.savefig(os.path.join(PLOT_PATH,'microbe-walls.png'),dpi=500)

print('Compiling microbial wall and richness plot')
fig, ax = plt.subplots(2,sharex=True)
#fig = plt.figure()
fig.suptitle('(run{0}-c{1}-r{2})'.format(run_id,combo_id,replicate))

axes = [ax[0], ax[0].twinx(), ax[1], ax[1].twinx()]
#axes[1].set_frame_on(False)
#axes[1].patch.set_visible(False)
#axes[1].set_xticks([])

microbe_stacked.plot.area(ax = axes[0],stacked=True,legend=False, linewidth=0,color=pal)
axes[0].set_ylabel(ylabel ='Microbial Immune Abundances N_i',labelpad=15,fontsize=7)
axes[0].set_xlabel(xlabel = 'Time t',fontsize=7)
axes[0].ticklabel_format(style='sci',scilimits=(0,0))
axes[0].axhline(upperThreshold, color='k', lw=.5, linestyle='dotted')
axes[0].axhline(lowerThreshold, color='k', lw=.5, linestyle='dotted')
axes[1].set_yticks([])
axes[1].fill_between(microbeWalls["t"],0,1, where=microbeWalls["wall_presence"]>0,interpolate=True, color="grey", alpha=0.5, linewidth=0,transform=axes[0].get_xaxis_transform())

microbeRichness.plot(x='t',xlabel = 'Time t',ax = axes[2],legend=False,color=pal[0],linewidth=0.75)
axes[2].set_ylabel(ylabel ='Microbial Immune Richness Sn',labelpad=15,fontsize=7)
axes[2].set_xlabel(xlabel = 'Time t',fontsize=7)
axes[2].ticklabel_format(style='sci',scilimits=(0,0))
axes[3].set_yticks([])
axes[3].fill_between(microbeWalls["t"],0,1, where=microbeWalls["wall_presence"]>0,interpolate=True, color="grey", alpha=0.5, linewidth=0,transform=axes[2].get_xaxis_transform())
fig.tight_layout()
fig.savefig(os.path.join(PLOT_PATH,'microbe-walls-richness.png'),dpi=500)

print('Compiling microbial wall plot projection on viral abundance and richness')
fig, ax = plt.subplots(2,sharex=True)
#fig = plt.figure()
fig.suptitle('(run{0}-c{1}-r{2})'.format(run_id,combo_id,replicate))
axes = [ax[0], ax[0].twinx(), ax[1], ax[1].twinx()]

virus_stacked.plot.area(ax = axes[0],stacked=True,legend=False, linewidth=0,color=pal)
axes[0].set_ylabel(ylabel ='Viral Strain Abundances V_i',labelpad=15,fontsize=7)
axes[0].set_xlabel(xlabel = 'Time t',fontsize=7)
axes[0].ticklabel_format(style='sci',scilimits=(0,0))
axes[1].set_yticks([])
axes[1].fill_between(microbeWalls["t"],0,1, where=microbeWalls["wall_presence"]>0,interpolate=True, color="grey", alpha=0.5, linewidth=0,transform=axes[0].get_xaxis_transform())

virusRichness.plot(x='t',xlabel = 'Time t',ax = axes[2],legend=False,color=pal[0],linewidth=0.75)
axes[2].set_ylabel(ylabel ='Viral Strain Richness Sv',labelpad=15,fontsize=7)
axes[2].set_xlabel(xlabel = 'Time t',fontsize=7)
axes[2].ticklabel_format(style='sci',scilimits=(0,0))
axes[3].set_yticks([])
axes[3].fill_between(microbeWalls["t"],0,1, where=microbeWalls["wall_presence"]>0,interpolate=True, color="grey", alpha=0.5, linewidth=0,transform=axes[2].get_xaxis_transform())
fig.tight_layout()
fig.savefig(os.path.join(PLOT_PATH,'microbe-walls-on-virus-richness.png'),dpi=500)




print('Complete!')
