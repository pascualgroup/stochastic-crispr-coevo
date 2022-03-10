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

threshold = sys.argv[2]/100

#con = sqlite3.connect('/Volumes/Yadgah/sweep_db_gathered.sqlite')
SCRIPT_PATH = os.path.abspath(os.path.dirname(__file__))
DB_PATH = os.path.join(SCRIPT_PATH,'..','sweep_db_gathered.sqlite')
# DB_PATH = os.path.join('/Volumes','Yadgah','run_id1455_combo73_replicate15.sqlite')
con = sqlite3.connect(DB_PATH)
cur = con.cursor()

ID = cur.execute('SELECT combo_id,replicate FROM runs WHERE run_id = {}'.format(run_id)).fetchall()
#ID = pd.read_sql_query("SELECT combo_id,replicate FROM runs WHERE run_id = {}".format(run_id), con)

combo_id = ID[0][0]
replicate = ID[0][1]

PLOT_PATH = os.path.join(SCRIPT_PATH,'..', 'plots','c{}'.format(combo_id),'r{}'.format(replicate))

#PLOT_PATH = os.path.abspath(os.path.dirname(__file__))

print('SQLite Query: microbe abundance data')
microbe_stacked = pd.read_sql_query("SELECT t,bstrain_id,abundance FROM babundance WHERE run_id = {}".format(run_id), con)
microbe_stacked = microbe_stacked.pivot(index='t',columns='bstrain_id',values='abundance')
microbe_stacked = microbe_stacked.fillna(0)
thresholdAbundances = threshold*microbe_stacked.sum(axis=1).values
dfSums = pd.DataFrame(list(zip(*[thresholdAbundances for i in range(len(microbe_stacked.columns))])),
    index=microbe_stacked.index, columns=microbe_stacked.columns)
microbe_stacked = microbe_stacked.where(microbe_stacked >= dfSums, None)

print('SQLite Query: virus abundance data')
virus_stacked = pd.read_sql_query("SELECT t,vstrain_id,abundance FROM vabundance WHERE run_id = {}".format(run_id), con)
virus_stacked = virus_stacked.pivot(index='t',columns='vstrain_id',values='abundance')
virus_stacked = virus_stacked.fillna(0)
thresholdAbundances = threshold*virus_stacked.sum(axis=1).values
dfSums = pd.DataFrame(list(zip(*[thresholdAbundances for i in range(len(virus_stacked.columns))])),
    index=virus_stacked.index, columns=virus_stacked.columns)
virus_stacked = virus_stacked.where(virus_stacked >= dfSums, None)


if microbe['t'][microbe['t'].size-1] not in virus_stacked.index:
    microbe_stacked.drop(microbe['t'][microbe['t'].size-1],inplace=True)

pal = sns.color_palette("tab20b")

print('Compiling microbial strain time series plot')
fig, ax = plt.subplots(1)
microbe_stacked.plot.area(ax=ax, stacked=True, legend=False, linewidth=0,color=pal)
ax.set_ylabel(ylabel ='Microbial Strain Abundances',labelpad=15,fontsize=7)
ax.set_xlabel(xlabel = 'Time t',fontsize=7)
ax.ticklabel_format(axis = 'y',style='sci',scilimits=(0,0))
ax.xaxis.set_minor_locator(ticker.MultipleLocator(50))
lim = ax.get_ylim()
ax.set_ylim(0,lim[1])
plt.tight_layout()
plt.savefig(os.path.join(PLOT_PATH,'microbe-strain-abundances.png'),dpi=500)

print('Compiling viral strain time series plot')
fig, ax = plt.subplots(1)
virus_stacked.plot.area(ax=ax, stacked=True, legend=False, linewidth=0,color=pal)
ax.set_ylabel(ylabel ='Viral Strain Abundances',labelpad=15,fontsize=7)
ax.set_xlabel(xlabel = 'Time t',fontsize=7)
ax.ticklabel_format(axis = 'y',style='sci',scilimits=(0,0))
ax.xaxis.set_minor_locator(ticker.MultipleLocator(50))
lim = ax.get_ylim()
ax.set_ylim(0,lim[1])
plt.tight_layout()
plt.savefig(os.path.join(PLOT_PATH,'virus-strain-abundances.png'),dpi=500)

print('Compiling microbial-virus stacked time series plots')
fig, ax = plt.subplots(2,sharex=True)
#fig = plt.figure()
fig.suptitle('Strain Abundances (run{0}-c{1}-r{2})'.format(run_id,combo_id,replicate))
axes = [ax[0], ax[1]]
microbe_stacked.plot.area(ax = axes[0],stacked=True,legend=False, linewidth=0,color=pal)
axes[0].set_ylabel(ylabel ='Microbial Immune Abundances N_i',labelpad=15,fontsize=7)
axes[0].set_xlabel(xlabel = 'Time t',fontsize=7)
axes[0].ticklabel_format(axis = 'y',style='sci',scilimits=(0,0))
axes[0].xaxis.set_minor_locator(ticker.MultipleLocator(50))
lim = axes[0].get_ylim()
axes[0].set_ylim(0,lim[1])
virus_stacked.plot.area(ax = axes[1],stacked=True,legend=False, linewidth=0,color=pal)
axes[1].set_ylabel(ylabel ='Viral Strain Abundances V_i',labelpad=15,fontsize=7)
axes[1].set_xlabel(xlabel = 'Time t',fontsize=7)
axes[1].ticklabel_format(axis = 'y',style='sci',scilimits=(0,0))
axes[1].xaxis.set_minor_locator(ticker.MultipleLocator(50))
lim = axes[1].get_ylim()
axes[1].set_ylim(0,lim[1])

fig.tight_layout()
fig.savefig(os.path.join(PLOT_PATH,'microbe-virus-stacked-abundances.png'),dpi=500)

print('Complete!')
