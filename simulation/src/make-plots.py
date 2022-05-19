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

run_id = sys.argv[1]

if len(sys.argv) == 2:
    threshold = 0
else:
    threshold = float(sys.argv[2])/100

SCRIPT_PATH = os.path.abspath(os.path.dirname(__file__))
DBSIM_PATH = os.path.join(SCRIPT_PATH,'..','sweep_db_gathered.sqlite')
# DBSIM_PATH = os.path.join('/Volumes','Yadgah','run_id1455_combo73_replicate15.sqlite')
# DBSIM_PATH = os.path.join('/Volumes','Yadgah','crispr-sweep-7-2-2022/isolates/runID1723-c35-r23/runID1723-c35-r23.sqlite')
conSim = sqlite3.connect(DBSIM_PATH)
curSim = conSim.cursor()

ID = curSim.execute('SELECT combo_id,replicate FROM runs WHERE run_id = {}'.format(run_id)).fetchall()
combo_id = ID[0][0]
replicate = ID[0][1]

PLOT_PATH = os.path.join(SCRIPT_PATH,'..', 'plots','c{}'.format(combo_id),'r{}'.format(replicate))

print('SQLite Query: microbe abundance data')
microbe_stacked = pd.read_sql_query("SELECT t,bstrain_id,abundance \
    FROM babundance WHERE run_id = {}".format(run_id), conSim)
microbeThreshAbund = pd.read_sql_query("SELECT t, microbial_abundance \
    FROM summary WHERE run_id = {}".format(run_id), conSim)
microbeThreshAbund['microbial_abundance'] = \
    threshold*microbeThreshAbund['microbial_abundance']
microbeThreshAbund = \
    [microbeThreshAbund[microbeThreshAbund.t == t].microbial_abundance.values[0]
    for t in microbe_stacked.t]
microbe_stacked = microbe_stacked[microbe_stacked.abundance >= microbeThreshAbund]
microbe_stacked = microbe_stacked.pivot(index='t',columns='bstrain_id',values='abundance')

print('SQLite Query: virus abundance data')
virus_stacked = pd.read_sql_query("SELECT t,vstrain_id,abundance \
    FROM vabundance WHERE run_id = {}".format(run_id), conSim)
virusThreshAbund = pd.read_sql_query("SELECT t, viral_abundance \
    FROM summary WHERE run_id = {}".format(run_id), conSim)
virusThreshAbund['viral_abundance'] = \
    threshold*virusThreshAbund['viral_abundance']
virusThreshAbund = \
    [virusThreshAbund[virusThreshAbund.t == t].viral_abundance.values[0]
    for t in virus_stacked.t]
virus_stacked = virus_stacked[virus_stacked.abundance >= virusThreshAbund]
virus_stacked = virus_stacked.pivot(index='t',columns='vstrain_id',values='abundance')

if len(microbe_stacked.index) > len(virus_stacked.index):
    microbe_stacked.drop(index=list(set(microbe_stacked.index)-set(virus_stacked.index)),inplace=True)

pal = sns.color_palette("tab20b")

print('Compiling microbial strain time series plot')
fig, ax = plt.subplots(1)
microbe_stacked.plot.area(ax=ax, stacked=True, legend=False, linewidth=0,color=pal)
ax.set_ylabel(ylabel ='Microbial Strain Abundances',labelpad=15,fontsize=7)
ax.set_xlabel(xlabel = 'Time t',fontsize=7)
ax.ticklabel_format(axis = 'y',style='sci',scilimits=(0,0))
ax.xaxis.set_minor_locator(ticker.MultipleLocator(25))
lim = ax.get_ylim()
ax.set_ylim(0,lim[1])
plt.tight_layout()
plt.savefig(os.path.join(PLOT_PATH,'microbe-strain-abundances.png'),dpi=500)
plt.close(fig)

print('Compiling viral strain time series plot')
fig, ax = plt.subplots(1)
virus_stacked.plot.area(ax=ax, stacked=True, legend=False, linewidth=0,color=pal)
ax.set_ylabel(ylabel ='Viral Strain Abundances',labelpad=15,fontsize=7)
ax.set_xlabel(xlabel = 'Time t',fontsize=7)
ax.ticklabel_format(axis = 'y',style='sci',scilimits=(0,0))
ax.xaxis.set_minor_locator(ticker.MultipleLocator(25))
lim = ax.get_ylim()
ax.set_ylim(0,lim[1])
plt.tight_layout()
plt.savefig(os.path.join(PLOT_PATH,'virus-strain-abundances.png'),dpi=500)
plt.close(fig)

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
fig.savefig(os.path.join(PLOT_PATH,'microbe-virus-stacked-abundances.png'),dpi=500)
plt.close(fig)

print('Complete!')
