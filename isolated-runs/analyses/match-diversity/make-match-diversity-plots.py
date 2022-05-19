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

run_id = sys.argv[1]

resolve = 500
imgType = "pdf" #or "png"

SCRIPT_PATH = os.path.abspath(os.path.dirname(__file__)) # cluster

DBSIM_PATH = os.path.join(SCRIPT_PATH,'..','..','..','simulation','sweep_db_gathered.sqlite') # cluster
# DBSIM_PATH = os.path.join('/Volumes','Yadgah','crispr-sweep-7-2-2022/isolates/runID3297-c66-r47/runID3297-c66-r47.sqlite')

conSim = sqlite3.connect(DBSIM_PATH)
curSim = conSim.cursor()
# Designating plot path from simulation data
ID = curSim.execute('SELECT combo_id,replicate FROM runs WHERE run_id = {}'.format(run_id)).fetchall()
combo_id = ID[0][0]
replicate = ID[0][1]
RUN_DIR = os.path.join('runID{0}-c{1}-r{2}'.format(run_id,combo_id,replicate))

DBMDIV_PATH = os.path.join(SCRIPT_PATH,'..','..','isolates',RUN_DIR,'match-diversity_output.sqlite') # cluster
# DBMDIV_PATH = os.path.join('/Volumes','Yadgah','crispr-sweep-7-2-2022/isolates/runID3297-c66-r47','match-diversity_output.sqlite') # local. run_id fixed; for testing

conMDiv = sqlite3.connect(DBMDIV_PATH)
curMDiv = conMDiv.cursor()


print('SQLite Query: microbial abundance time series data')
microbe_stacked = pd.read_sql_query("SELECT t,bstrain_id,abundance FROM babundance WHERE run_id = {}".format(run_id), conSim)
microbe_stacked = microbe_stacked.pivot(index='t',columns='bstrain_id',values='abundance')
print('SQLite Query: viral abundance time series data')
virus_stacked = pd.read_sql_query("SELECT t,vstrain_id,abundance FROM vabundance WHERE run_id = {}".format(run_id), conSim)
virus_stacked = virus_stacked.pivot(index='t',columns='vstrain_id',values='abundance')

if len(microbe_stacked.index) > len(virus_stacked.index):
    microbe_stacked.drop(index=list(set(microbe_stacked.index)-set(virus_stacked.index)),inplace=True)


print('SQLite Query: shannon match diversity')
vMatchDiv = pd.read_sql_query("SELECT t,virus_shannon_diversity \
    FROM single_locus_escape_match_diversity", conMDiv)
bMatchDiv = pd.read_sql_query("SELECT t,microbe_shannon_diversity \
    FROM single_locus_escape_match_diversity", conMDiv)
print('SQLite Query: single locus abundances')
vLocusAbunds = pd.read_sql_query("SELECT t,spacer_id,vabundance \
FROM single_locus_escape_match_abundances", conMDiv)
bLocusAbunds = pd.read_sql_query("SELECT t,spacer_id,babundance \
FROM single_locus_escape_match_abundances", conMDiv)
vLocusAbunds = vLocusAbunds.pivot(index='t',columns='spacer_id',values='vabundance')
bLocusAbunds = bLocusAbunds.pivot(index='t',columns='spacer_id',values='babundance')
print('SQLite Query: match abundances')
vMatchAbunds = pd.read_sql_query("SELECT t,match0,match1,match2,match3 \
FROM vMatchAbundances", conMDiv)
bMatchAbunds = pd.read_sql_query("SELECT t,match0,match1,match2,match3 \
FROM bMatchAbundances", conMDiv)

if len(vLocusAbunds.index) > len(virus_stacked.index):
    vLocusAbunds.drop(index=list(set(vLocusAbunds.index)-set(virus_stacked.index)),inplace=True)

if len(bLocusAbunds.index) > len(virus_stacked.index):
    bLocusAbunds.drop(index=list(set(bLocusAbunds.index)-set(virus_stacked.index)),inplace=True)

if len(vMatchAbunds.t) > len(virus_stacked.index):
    vMatchAbunds = vMatchAbunds.loc[vMatchAbunds.t.isin(virus_stacked.index.values)]

if len(bMatchAbunds.t) > len(virus_stacked.index):
    bMatchAbunds = bMatchAbunds.loc[bMatchAbunds.t.isin(virus_stacked.index.values)]


pal = sns.color_palette("tab20b")
# fig, ax = plt.subplots(1)
#fig = plt.figure()
fig, ax = plt.subplots(2,sharex=True)
#fig = plt.figure()
fig.suptitle('Strain Abundances (run{0}-c{1}-r{2})'.format(run_id,combo_id,replicate))
axes = [ax[0], ax[1]]
microbe_stacked.plot.area(ax = axes[0],stacked=True,legend=False, linewidth=0,color=pal)
axes[1].plot(bMatchDiv.t,bMatchDiv.microbe_shannon_diversity)
fig, ax = plt.subplots(1)
ax.plot(vMatchDiv.t,vMatchDiv.virus_shannon_diversity)

fig, ax = plt.subplots(1)
vLocusAbunds.plot.area(ax = ax,stacked=True,legend=True, linewidth=0,color=sns.color_palette("tab20c"))
plt.show()



fig, ax = plt.subplots(2,sharex=True)
axes = [ax[0], ax[1]]
microbe_stacked.plot.area(ax = axes[0],stacked=True,legend=False, linewidth=0,color=pal)
axes[1].plot(bMatchAbunds.t,bMatchAbunds.match0)
fig, ax = plt.subplots(2,sharex=True)
axes = [ax[0], ax[1]]
microbe_stacked.plot.area(ax = axes[0],stacked=True,legend=False, linewidth=0,color=pal)
axes[1].plot(vMatchAbunds.t,vMatchAbunds.match0, linestyle='solid', label = '0-matches')
axes[1].plot(vMatchAbunds.t,vMatchAbunds.match1, linestyle=(0,(5,1)), label = '1-matches')
axes[1].plot(vMatchAbunds.t,vMatchAbunds.match2, linestyle=(0,(5,5)), label = '2-matches')
axes[1].plot(vMatchAbunds.t,vMatchAbunds.match3, linestyle=(0,(5,10)), label = '3-matches')
axes[1].legend(loc = 'upper left')

fig, ax = plt.subplots(2,sharex=True)
axes = [ax[0], ax[1]]
microbe_stacked.plot.area(ax = axes[0],stacked=True,legend=False, linewidth=0,color=pal)
axes[1].plot(bMatchAbunds.t,bMatchAbunds.match0, linestyle='solid', label = '0-matches')
axes[1].plot(bMatchAbunds.t,bMatchAbunds.match1, linestyle='solid', label = '1-matches')
axes[1].plot(bMatchAbunds.t,bMatchAbunds.match2, linestyle='solid', label = '2-matches')
axes[1].plot(bMatchAbunds.t,bMatchAbunds.match3, linestyle='solid', label = '3-matches')
axes[1].legend(loc = 'upper left')


pal = sns.color_palette("tab20b")
print('Compiling microbial-virus stacked time series plots')
fig, ax = plt.subplots(2,sharex=True)
#fig = plt.figure()
fig.suptitle('Strain Abundances (run{0}-c{1}-r{2})'.format(run_id,combo_id,replicate))
axes = [ax[0], ax[1]]
microbe_stacked.plot.area(ax = axes[0],stacked=True,legend=False, linewidth=0,color=pal)
axes[0].set_ylabel(ylabel ='Microbial Immune Abundances N_i',labelpad=15,fontsize=7)
axes[0].set_xlabel(xlabel = 'Time t',fontsize=7)
axes[0].ticklabel_format(axis = 'y',style='sci',scilimits=(0,0))
axes[0].xaxis.set_minor_locator(ticker.MultipleLocator(25))
lim = axes[0].get_ylim()
axes[0].set_ylim(0,lim[1])
virus_stacked.plot.area(ax = axes[1],stacked=True,legend=False, linewidth=0,color=pal)
axes[1].set_ylabel(ylabel ='Viral Strain Abundances V_i',labelpad=15,fontsize=7)
axes[1].set_xlabel(xlabel = 'Time t',fontsize=7)
axes[1].ticklabel_format(axis = 'y',style='sci',scilimits=(0,0))
axes[1].xaxis.set_minor_locator(ticker.MultipleLocator(25))
lim = axes[1].get_ylim()
axes[1].set_ylim(0,lim[1])

fig.tight_layout()
fig.savefig(os.path.join(SCRIPT_PATH,'..','..','isolates',RUN_DIR,'microbe-virus-stacked-abundances.{}'.format(imgType)),dpi=resolve)
# fig.savefig(os.path.join('/Volumes/Yadgah','microbe-virus-stacked-abundances.{}'.format(imgType)),dpi=resolve)
plt.close(fig)
