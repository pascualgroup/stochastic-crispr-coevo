#!/usr/bin/env python3

import pandas as pd
import numpy as np
import math
import matplotlib.pyplot as plt
import sys
import os
import seaborn as sns
import sqlite3
import matplotlib.ticker as ticker
import matplotlib.gridspec as grid_spec
# import plotly.graph_objects as go
# import plotly.express as px
# from mpl_toolkits.mplot3d import Axes3D
import itertools



run_id = sys.argv[1]

resolve = 500
imgType = "png" #or "png"
vabundthreshold = .05

SCRIPT_PATH = os.path.abspath(os.path.dirname(__file__)) # cluster
dir = 'crispr-sweep-7-2-2022/isolates/runID3297-c66-r47'
run = 'runID3297-c66-r47'

DBSIM_PATH = os.path.join(SCRIPT_PATH,'..','..','..','simulation','sweep_db_gathered.sqlite') # cluster
# DBSIM_PATH = os.path.join('/Volumes','Yadgah',dir,'{}.sqlite'.format(run))

conSim = sqlite3.connect(DBSIM_PATH)
curSim = conSim.cursor()
# Designating plot path from simulation data
ID = curSim.execute('SELECT combo_id,replicate FROM runs WHERE run_id = {}'.format(run_id)).fetchall()
combo_id = ID[0][0]
replicate = ID[0][1]
CC = curSim.execute('SELECT microbe_carrying_capacity FROM param_combos WHERE combo_id = {}'.format(combo_id)).fetchall()
CC = CC[0][0]
RUN_DIR = os.path.join('runID{0}-c{1}-r{2}'.format(run_id,combo_id,replicate))

DBPROB_PATH = os.path.join(SCRIPT_PATH,'..','..','isolates',RUN_DIR,'probability-of-emergence_output.sqlite') # cluster
DBMATCH_PATH = os.path.join(SCRIPT_PATH,'..','..','isolates',RUN_DIR,'matches_output.sqlite') # cluster
DBTREE_PATH = os.path.join(SCRIPT_PATH,'..','..','isolates',RUN_DIR,'trees_output.sqlite') # cluster
DBTRI_PATH = os.path.join(SCRIPT_PATH,'..','..','isolates',RUN_DIR,'tripartite-networks_output.sqlite') # cluster
PLOT_PATH = os.path.join(SCRIPT_PATH,'..','..','isolates',RUN_DIR)
# DBPROB_PATH = os.path.join('/Volumes','Yadgah',dir,'probability-of-emergence_output.sqlite') # local
# DBPROB_PATH = os.path.join('/Volumes','Yadgah',dir,'emergence-lambert-root_output.sqlite')
# DBMATCH_PATH = os.path.join('/Volumes','Yadgah',dir,'matches_output.sqlite') # local
# DBTREE_PATH = os.path.join('/Volumes','Yadgah',dir,'trees_output.sqlite') # local
# DBTRI_PATH = os.path.join('/Volumes','Yadgah',dir,'tripartite-networks_output.sqlite') # local
# PLOT_PATH = os.path.join('/Volumes','Yadgah') # local

conProb = sqlite3.connect(DBPROB_PATH)
curProb = conProb.cursor()
conMatch = sqlite3.connect(DBMATCH_PATH)
curMatch = conMatch.cursor()

if os.path.exists(DBTRI_PATH) & os.path.exists(DBTREE_PATH):
    conTri = sqlite3.connect(DBTRI_PATH)
    curTri = conTri.cursor()
    conTree = sqlite3.connect(DBTREE_PATH)
    curTree = conTree.cursor()

print('SQLite Query: viral abundance time series data')
virus_stacked = pd.read_sql_query("SELECT t,vstrain_id,abundance FROM vabundance WHERE run_id = {}".format(run_id), conSim)
virus_total = pd.read_sql_query("SELECT t,viral_abundance FROM summary WHERE run_id = {}".format(run_id), conSim)\
.rename(columns={"viral_abundance": "vtotal"})
virus_total = virus_total[virus_total.t <= max(virus_stacked.t)]
# THIS IS FOR THE 3D DISTRIBUTIONS....
vMaxIDs = virus_stacked.set_index('t').groupby(['vstrain_id']).agg(t = ('abundance','idxmax'),\
                            vmaxAbund = ('abundance','max')).reset_index()
vMaxIDs = virus_total.merge(vMaxIDs,on=['t'])
vMaxIDs['vtotal'] = vabundthreshold*np.array(vMaxIDs['vtotal'])
keepStrains = list(vMaxIDs[vMaxIDs['vmaxAbund']>vMaxIDs['vtotal']]['vstrain_id'].values)
virus_stacked = virus_stacked[[(i in keepStrains) for i in virus_stacked.vstrain_id]].reset_index(drop=True)

print('SQLite Query: microbial abundance time series data')
microbe_stacked = pd.read_sql_query("SELECT t,bstrain_id,abundance FROM babundance WHERE run_id = {}".format(run_id), conSim)
microbe_total = pd.read_sql_query("SELECT t,microbial_abundance FROM summary WHERE run_id = {}".format(run_id), conSim)\
.rename(columns={"microbial_abundance": "btotal"})
microbe_stacked = microbe_stacked[microbe_stacked.t <= max(virus_stacked.t)]
microbe_total = microbe_total[microbe_total.t <= max(virus_stacked.t)]
microbe_stacked = microbe_stacked.pivot(index='t',columns='bstrain_id',values='abundance')
virus_stacked = virus_stacked.pivot(index='t',columns='vstrain_id',values='abundance')



pal = sns.color_palette("tab20c")
fig, ax = plt.subplots(2,sharex=True)
microbe_stacked.plot.area(ax = ax[0],stacked=True,legend=False, linewidth=0,color=pal,sort_columns=True)
ax[0].set_ylabel(ylabel ='Microbial Strain Abundances',labelpad=15,fontsize=15)
ax[0].set_xlabel(xlabel = 'Time t',fontsize=15)
ax[0].ticklabel_format(axis = 'y',style='sci',scilimits=(0,0))
ax[0].xaxis.set_minor_locator(ticker.MultipleLocator(25))
ax[0].tick_params(axis='x', labelsize= 10)
ax[0].tick_params(axis='y', labelsize= 10)
lim = ax[0].get_ylim()
ax[0].set_ylim(0,lim[1])
virus_stacked.plot.area(ax = ax[1],stacked=True,legend=False, linewidth=0,color=pal,sort_columns=True)
ax[1].set_ylabel(ylabel ='Viral Strain Abundances',labelpad=15,fontsize=15)
ax[1].set_xlabel(xlabel = 'Time t',fontsize=15)
ax[1].ticklabel_format(axis = 'y',style='sci',scilimits=(0,0))
ax[1].xaxis.set_minor_locator(ticker.MultipleLocator(25))
ax[1].tick_params(axis='x', labelsize= 10)
ax[1].tick_params(axis='y', labelsize= 10)
lim = ax[1].get_ylim()
ax[1].set_ylim(0,lim[1])


pEmergeExpected = pd.read_sql_query("SELECT t, p_emerge_weighted \
    FROM vmatch_extinction ORDER BY t", conProb).groupby(['t'])\
    .agg(p_exp=('p_emerge_weighted','sum')).reset_index()
pal = sns.color_palette("tab20c")
fig, ax = plt.subplots(2,sharex=True)
axes = [ax[0], ax[0].twinx(), ax[1], ax[1].twinx()]
microbe_stacked.plot.area(ax = axes[0],stacked=True,legend=False, linewidth=0,color=pal,sort_columns=True)
axes[0].set_ylabel(ylabel ='Microbial Strain Abundances',labelpad=15,fontsize=15)
axes[0].set_xlabel(xlabel = 'Time t',fontsize=15)
axes[0].ticklabel_format(axis = 'y',style='sci',scilimits=(0,0))
axes[0].xaxis.set_minor_locator(ticker.MultipleLocator(25))
axes[0].tick_params(axis='x', labelsize= 10)
axes[0].tick_params(axis='y', labelsize= 10)
lim = axes[0].get_ylim()
axes[0].set_ylim(0,lim[1])
axes[1].plot(virus_total['t'],virus_total['vtotal'],linewidth=0,color='grey')
axes[1].fill_between(virus_total['t'],virus_total['vtotal'], color='grey',alpha=0.6)
lim = axes[1].get_ylim()
axes[1].set_ylim(0,lim[1])
axes[2].plot(pEmergeExpected['t'],pEmergeExpected['p_exp'],linewidth=1,color='darkblue')
axes[2].set_ylabel(ylabel ='Probability of Viral Strain Emergence',labelpad=15,fontsize=15)
axes[2].set_xlabel(xlabel = 'Time t',fontsize=15)
axes[2].tick_params(axis='x', labelsize= 10)
axes[2].tick_params(axis='y', labelsize= 10)
axes[3].plot(virus_total['t'],virus_total['vtotal'],linewidth=0,color='grey')
axes[3].fill_between(virus_total['t'],virus_total['vtotal'], color='grey',alpha=0.3)
fig.tight_layout()


pEmergeExpected = pd.read_sql_query("SELECT t, pactual_emerge_weighted \
    FROM existing_vmatch_extinction ORDER BY t", conProb).groupby(['t'])\
    .agg(p_exp=('pactual_emerge_weighted','sum')).reset_index()
pal = sns.color_palette("tab20c")
fig, ax = plt.subplots(2,sharex=True)
axes = [ax[0], ax[0].twinx(), ax[1], ax[1].twinx()]
microbe_stacked.plot.area(ax = axes[0],stacked=True,legend=False, linewidth=0,color=pal,sort_columns=True)
axes[0].set_ylabel(ylabel ='Microbial Strain Abundances',labelpad=15,fontsize=15)
axes[0].set_xlabel(xlabel = 'Time t',fontsize=15)
axes[0].ticklabel_format(axis = 'y',style='sci',scilimits=(0,0))
axes[0].xaxis.set_minor_locator(ticker.MultipleLocator(25))
axes[0].tick_params(axis='x', labelsize= 10)
axes[0].tick_params(axis='y', labelsize= 10)
lim = axes[0].get_ylim()
axes[0].set_ylim(0,lim[1])
axes[1].plot(virus_total['t'],virus_total['vtotal'],linewidth=0,color='grey')
axes[1].fill_between(virus_total['t'],virus_total['vtotal'], color='grey',alpha=0.6)
lim = axes[1].get_ylim()
axes[1].set_ylim(0,lim[1])
axes[2].plot(pEmergeExpected['t'],pEmergeExpected['p_exp'],linewidth=1,color='darkblue')
axes[2].set_ylabel(ylabel ='Probability of Viral Strain Emergence',labelpad=15,fontsize=15)
axes[2].set_xlabel(xlabel = 'Time t',fontsize=15)
axes[2].tick_params(axis='x', labelsize= 10)
axes[2].tick_params(axis='y', labelsize= 10)
axes[3].plot(virus_total['t'],virus_total['vtotal'],linewidth=0,color='grey')
axes[3].fill_between(virus_total['t'],virus_total['vtotal'], color='grey',alpha=0.3)
fig.tight_layout()


pEmergeExpected = pd.read_sql_query("SELECT t, plambert_emerge_weighted \
    FROM existing_vmatch_extinction ORDER BY t", conProb).groupby(['t'])\
    .agg(p_exp=('plambert_emerge_weighted','sum')).reset_index()
pal = sns.color_palette("tab20c")
fig, ax = plt.subplots(2,sharex=True)
axes = [ax[0], ax[0].twinx(), ax[1], ax[1].twinx()]
microbe_stacked.plot.area(ax = axes[0],stacked=True,legend=False, linewidth=0,color=pal,sort_columns=True)
axes[0].set_ylabel(ylabel ='Microbial Strain Abundances',labelpad=15,fontsize=15)
axes[0].set_xlabel(xlabel = 'Time t',fontsize=15)
axes[0].ticklabel_format(axis = 'y',style='sci',scilimits=(0,0))
axes[0].xaxis.set_minor_locator(ticker.MultipleLocator(25))
axes[0].tick_params(axis='x', labelsize= 10)
axes[0].tick_params(axis='y', labelsize= 10)
lim = axes[0].get_ylim()
axes[0].set_ylim(0,lim[1])
axes[1].plot(virus_total['t'],virus_total['vtotal'],linewidth=0,color='grey')
axes[1].fill_between(virus_total['t'],virus_total['vtotal'], color='grey',alpha=0.6)
lim = axes[1].get_ylim()
axes[1].set_ylim(0,lim[1])
axes[2].plot(pEmergeExpected['t'],pEmergeExpected['p_exp'],linewidth=1,color='darkblue')
axes[2].set_ylabel(ylabel ='Probability of Viral Strain Emergence',labelpad=15,fontsize=15)
axes[2].set_xlabel(xlabel = 'Time t',fontsize=15)
axes[2].tick_params(axis='x', labelsize= 10)
axes[2].tick_params(axis='y', labelsize= 10)
axes[3].plot(virus_total['t'],virus_total['vtotal'],linewidth=0,color='grey')
axes[3].fill_between(virus_total['t'],virus_total['vtotal'], color='grey',alpha=0.3)
fig.tight_layout()



pEmergeExpected = pd.read_sql_query("SELECT t, vmatch_id, p_extinction_root, vfrequency \
    FROM existing_vmatch_extinction ORDER BY t", conProb)\
    .merge(frequencies,on=['t','vmatch_id']).rename(columns={"vfrequency":"vfreq"})
pEmergeExpected['p_extinction_root'] = 1- np.array(pEmergeExpected['p_extinction_root'])
pEmergeExpected=pEmergeExpected.rename(columns={"p_extinction_root":"proot_emerge_weighted"})
pEmergeExpected['proot_emerge_weighted'] = \
    pEmergeExpected['proot_emerge_weighted']*pEmergeExpected['vfreq']
pEmergeExpected = pEmergeExpected.groupby(['t'])\
    .agg(p_exp=('proot_emerge_weighted','sum')).reset_index()
pal = sns.color_palette("tab20c")
fig, ax = plt.subplots(2,sharex=True)
axes = [ax[0], ax[0].twinx(), ax[1], ax[1].twinx()]
microbe_stacked.plot.area(ax = axes[0],stacked=True,legend=False, linewidth=0,color=pal,sort_columns=True)
axes[0].set_ylabel(ylabel ='Microbial Strain Abundances',labelpad=15,fontsize=15)
axes[0].set_xlabel(xlabel = 'Time t',fontsize=15)
axes[0].ticklabel_format(axis = 'y',style='sci',scilimits=(0,0))
axes[0].xaxis.set_minor_locator(ticker.MultipleLocator(25))
axes[0].tick_params(axis='x', labelsize= 10)
axes[0].tick_params(axis='y', labelsize= 10)
lim = axes[0].get_ylim()
axes[0].set_ylim(0,lim[1])
axes[1].plot(virus_total['t'],virus_total['vtotal'],linewidth=0,color='grey')
axes[1].fill_between(virus_total['t'],virus_total['vtotal'], color='grey',alpha=0.6)
lim = axes[1].get_ylim()
axes[1].set_ylim(0,lim[1])
axes[2].plot(pEmergeExpected['t'],pEmergeExpected['p_exp'],linewidth=1,color='darkblue')
axes[2].set_ylabel(ylabel ='Probability of Viral Strain Emergence',labelpad=15,fontsize=15)
axes[2].set_xlabel(xlabel = 'Time t',fontsize=15)
axes[2].tick_params(axis='x', labelsize= 10)
axes[2].tick_params(axis='y', labelsize= 10)
axes[3].plot(virus_total['t'],virus_total['vtotal'],linewidth=0,color='grey')
axes[3].fill_between(virus_total['t'],virus_total['vtotal'], color='grey',alpha=0.3)
fig.tight_layout()




# THIS IS MULTIPLYING WITH MICROBE FREQUENCY BUT NORMALIZING WITH VIRUS & MICROBE
frequencies = pd.read_sql_query("SELECT t,match_id,bsusceptible \
        FROM vmatches_abundances", conTri)\
        .rename(columns={"match_id":"vmatch_id"})\
        .merge(microbe_total,on = ['t'])\
        .rename(columns={"bsusceptible":"bfreq"})
frequencies['bfreq'] = frequencies['bfreq']/frequencies['btotal']
norm = frequencies.groupby(['t']).agg(norm=('bfreq','sum')).reset_index()
frequencies = frequencies.merge(norm,on=['t'])
frequencies['bfreq'] = frequencies['bfreq']/frequencies['norm']
frequencies = frequencies.drop(columns=['btotal','norm'])
pEmergeExpected = pd.read_sql_query("SELECT t, vmatch_id, p_extinction_lambert, vfrequency \
    FROM existing_vmatch_extinction ORDER BY t", conProb)\
    .merge(frequencies,on=['t','vmatch_id']).rename(columns={"vfrequency":"vfreq"})
pEmergeExpected['vfreq'] = pEmergeExpected['vfreq']/pEmergeExpected['bfreq']
# n = pEmergeExpected.groupby(['t']).agg(n=('vfreq','size')).reset_index()
# pEmergeExpected = pEmergeExpected.merge(n,on=['t'])
# pEmergeExpected['vfreq'] = pEmergeExpected['vfreq']/pEmergeExpected['n']
# pEmergeExpected['vfreq'] = pEmergeExpected['bfreq']
# norm = pEmergeExpected.groupby(['t']).agg(norm=('vfreq','sum')).reset_index()
# pEmergeExpected = pEmergeExpected.merge(norm,on=['t'])
# pEmergeExpected['vfreq'] = pEmergeExpected['vfreq']/pEmergeExpected['norm']
pEmergeExpected['p_extinction_lambert'] = 1- np.array(pEmergeExpected['p_extinction_lambert'])
# pEmergeExpected['p_extinction_lambert'] = np.ones(len(pEmergeExpected['p_extinction_lambert']))
pEmergeExpected=pEmergeExpected.rename(columns={"p_extinction_lambert":"plambert_emerge_weighted"})
pEmergeExpected['plambert_emerge_weighted'] = \
    pEmergeExpected['plambert_emerge_weighted']*pEmergeExpected['vfreq']
# pEmergeExpected['plambert_emerge_weighted'] = 1/np.array(pEmergeExpected['plambert_emerge_weighted'])
pEmergeExpected = pEmergeExpected.groupby(['t'])\
    .agg(p_exp=('plambert_emerge_weighted','sum')).reset_index()
# pEmergeExpected['p_exp'] = 1/np.array(pEmergeExpected['p_exp'])
pal = sns.color_palette("tab20c")
fig, ax = plt.subplots(2,sharex=True)
axes = [ax[0], ax[0].twinx(), ax[1], ax[1].twinx()]
microbe_stacked.plot.area(ax = axes[0],stacked=True,legend=False, linewidth=0,color=pal,sort_columns=True)
axes[0].set_ylabel(ylabel ='Microbial Strain Abundances',labelpad=15,fontsize=15)
axes[0].set_xlabel(xlabel = 'Time t',fontsize=15)
axes[0].ticklabel_format(axis = 'y',style='sci',scilimits=(0,0))
axes[0].xaxis.set_minor_locator(ticker.MultipleLocator(25))
axes[0].tick_params(axis='x', labelsize= 10)
axes[0].tick_params(axis='y', labelsize= 10)
lim = axes[0].get_ylim()
axes[0].set_ylim(0,lim[1])
axes[1].plot(virus_total['t'],virus_total['vtotal'],linewidth=0,color='grey')
axes[1].fill_between(virus_total['t'],virus_total['vtotal'], color='grey',alpha=0.6)
lim = axes[1].get_ylim()
axes[1].set_ylim(0,lim[1])
axes[2].plot(pEmergeExpected['t'],pEmergeExpected['p_exp'],linewidth=1,color='darkblue')
axes[2].set_ylabel(ylabel ='Probability of Viral Strain Emergence',labelpad=15,fontsize=15)
axes[2].set_xlabel(xlabel = 'Time t',fontsize=15)
axes[2].tick_params(axis='x', labelsize= 10)
axes[2].tick_params(axis='y', labelsize= 10)
axes[3].plot(virus_total['t'],virus_total['vtotal'],linewidth=0,color='grey')
axes[3].fill_between(virus_total['t'],virus_total['vtotal'], color='grey',alpha=0.3)
fig.tight_layout()

# THIS IS MULTIPLYING WITH MICROBE FREQUENCY BUT NORMALIZING WITH VIRUS & MICROBE & p-ext
frequencies = pd.read_sql_query("SELECT t,match_id,bsusceptible \
        FROM vmatches_abundances", conTri)\
        .rename(columns={"match_id":"vmatch_id"})\
        .merge(microbe_total,on = ['t'])\
        .rename(columns={"bsusceptible":"bfreq"})
frequencies['bfreq'] = frequencies['bfreq']/frequencies['btotal']
norm = frequencies.groupby(['t']).agg(norm=('bfreq','sum')).reset_index()
frequencies = frequencies.merge(norm,on=['t'])
frequencies['bfreq'] = frequencies['bfreq']/frequencies['norm']
frequencies = frequencies.drop(columns=['btotal','norm'])
pEmergeExpected = pd.read_sql_query("SELECT t, vmatch_id, p_extinction_lambert, vfrequency \
    FROM existing_vmatch_extinction ORDER BY t", conProb)\
    .merge(frequencies,on=['t','vmatch_id']).rename(columns={"vfrequency":"vfreq"})
pEmergeExpected['vfreq'] = pEmergeExpected['vfreq']*pEmergeExpected['bfreq']
norm = pEmergeExpected.groupby(['t']).agg(norm=('vfreq','sum')).reset_index()
pEmergeExpected = pEmergeExpected.merge(norm,on=['t'])
# pEmergeExpected['vfreq'] = pEmergeExpected['vfreq']/pEmergeExpected['norm']
pEmergeExpected = pEmergeExpected.drop(columns=['norm'])
pEmergeExpected['p_extinction_lambert'] = 1- np.array(pEmergeExpected['p_extinction_lambert'])
norm = pEmergeExpected.groupby(['t']).agg(norm=('p_extinction_lambert','sum')).reset_index()
pEmergeExpected = pEmergeExpected.merge(norm,on=['t'])
pEmergeExpected['p_extinction_lambert'] = \
    pEmergeExpected['p_extinction_lambert']/pEmergeExpected['norm']
pEmergeExpected = pEmergeExpected.drop(columns=['norm'])
pEmergeExpected=pEmergeExpected.rename(columns={"p_extinction_lambert":"plambert_emerge_weighted"})
pEmergeExpected['plambert_emerge_weighted'] = \
    pEmergeExpected['plambert_emerge_weighted']*pEmergeExpected['vfreq']
pEmergeExpected = pEmergeExpected.groupby(['t'])\
    .agg(p_exp=('plambert_emerge_weighted','sum')).reset_index()
pal = sns.color_palette("tab20c")
fig, ax = plt.subplots(2,sharex=True)
axes = [ax[0], ax[0].twinx(), ax[1], ax[1].twinx()]
microbe_stacked.plot.area(ax = axes[0],stacked=True,legend=False, linewidth=0,color=pal,sort_columns=True)
axes[0].set_ylabel(ylabel ='Microbial Strain Abundances',labelpad=15,fontsize=15)
axes[0].set_xlabel(xlabel = 'Time t',fontsize=15)
axes[0].ticklabel_format(axis = 'y',style='sci',scilimits=(0,0))
axes[0].xaxis.set_minor_locator(ticker.MultipleLocator(25))
axes[0].tick_params(axis='x', labelsize= 10)
axes[0].tick_params(axis='y', labelsize= 10)
lim = axes[0].get_ylim()
axes[0].set_ylim(0,lim[1])
axes[1].plot(virus_total['t'],virus_total['vtotal'],linewidth=0,color='grey')
axes[1].fill_between(virus_total['t'],virus_total['vtotal'], color='grey',alpha=0.6)
lim = axes[1].get_ylim()
axes[1].set_ylim(0,lim[1])
axes[2].plot(pEmergeExpected['t'],pEmergeExpected['p_exp'],linewidth=1,color='darkblue')
axes[2].set_ylabel(ylabel ='Probability of Viral Strain Emergence',labelpad=15,fontsize=15)
axes[2].set_xlabel(xlabel = 'Time t',fontsize=15)
axes[2].tick_params(axis='x', labelsize= 10)
axes[2].tick_params(axis='y', labelsize= 10)
axes[3].plot(virus_total['t'],virus_total['vtotal'],linewidth=0,color='grey')
axes[3].fill_between(virus_total['t'],virus_total['vtotal'], color='grey',alpha=0.3)
fig.tight_layout()



# THIS IS WITH DIVIDING THE TRIPARTITE MICROBE FREQUENCY!
frequencies = pd.read_sql_query("SELECT t,match_id,bsusceptible \
        FROM vmatches_abundances", conTri)\
        .rename(columns={"match_id":"vmatch_id"})\
        .merge(microbe_total,on = ['t'])\
        .rename(columns={"bsusceptible":"bfreq"})
frequencies['bfreq'] = frequencies['bfreq']/frequencies['btotal']
norm = frequencies.groupby(['t']).agg(norm=('bfreq','sum')).reset_index()
frequencies = frequencies.merge(norm,on=['t'])
frequencies['bfreq'] = frequencies['bfreq']/frequencies['norm']
frequencies = frequencies.drop(columns=['btotal','norm'])
pEmergeExpected = pd.read_sql_query("SELECT t, vmatch_id, plambert_emerge_weighted \
    FROM existing_vmatch_extinction ORDER BY t", conProb)\
    .merge(frequencies,on=['t','vmatch_id'])
pEmergeExpected['plambert_emerge_weighted'] = pEmergeExpected['plambert_emerge_weighted']/pEmergeExpected['bfreq']
pEmergeExpected = pEmergeExpected.groupby(['t'])\
    .agg(p_exp=('plambert_emerge_weighted','sum')).reset_index()
pal = sns.color_palette("tab20c")
fig, ax = plt.subplots(2,sharex=True)
axes = [ax[0], ax[0].twinx(), ax[1], ax[1].twinx()]
microbe_stacked.plot.area(ax = axes[0],stacked=True,legend=False, linewidth=0,color=pal,sort_columns=True)
axes[0].set_ylabel(ylabel ='Microbial Strain Abundances',labelpad=15,fontsize=15)
axes[0].set_xlabel(xlabel = 'Time t',fontsize=15)
axes[0].ticklabel_format(axis = 'y',style='sci',scilimits=(0,0))
axes[0].xaxis.set_minor_locator(ticker.MultipleLocator(25))
axes[0].tick_params(axis='x', labelsize= 10)
axes[0].tick_params(axis='y', labelsize= 10)
lim = axes[0].get_ylim()
axes[0].set_ylim(0,lim[1])
axes[1].plot(virus_total['t'],virus_total['vtotal'],linewidth=0,color='grey')
axes[1].fill_between(virus_total['t'],virus_total['vtotal'], color='grey',alpha=0.6)
lim = axes[1].get_ylim()
axes[1].set_ylim(0,lim[1])
axes[2].plot(pEmergeExpected['t'],pEmergeExpected['p_exp'],linewidth=1,color='darkblue')
axes[2].set_ylabel(ylabel ='Probability of Viral Strain Emergence',labelpad=15,fontsize=15)
axes[2].set_xlabel(xlabel = 'Time t',fontsize=15)
axes[2].tick_params(axis='x', labelsize= 10)
axes[2].tick_params(axis='y', labelsize= 10)
axes[3].plot(virus_total['t'],virus_total['vtotal'],linewidth=0,color='grey')
axes[3].fill_between(virus_total['t'],virus_total['vtotal'], color='grey',alpha=0.3)
fig.tight_layout()


# THIS IS WITH DIVIDING THE TRIPARTITE MICROBE FREQUENCY!
frequencies = pd.read_sql_query("SELECT t,match_id,bsusceptible \
        FROM vmatches_abundances", conTri)\
        .rename(columns={"match_id":"vmatch_id"})\
        .merge(microbe_total,on = ['t'])\
        .rename(columns={"bsusceptible":"bfreq"})
frequencies['bfreq'] = frequencies['bfreq']/frequencies['btotal']
norm = frequencies.groupby(['t']).agg(norm=('bfreq','sum')).reset_index()
frequencies = frequencies.merge(norm,on=['t'])
frequencies['bfreq'] = frequencies['bfreq']/frequencies['norm']
frequencies = frequencies.drop(columns=['btotal','norm'])
pEmergeExpected = pd.read_sql_query("SELECT t, vmatch_id, pactual_emerge_weighted \
    FROM existing_vmatch_extinction ORDER BY t", conProb)\
    .merge(frequencies,on=['t','vmatch_id'])
pEmergeExpected['pactual_emerge_weighted'] = pEmergeExpected['pactual_emerge_weighted']/pEmergeExpected['bfreq']
pEmergeExpected = pEmergeExpected.groupby(['t'])\
    .agg(p_exp=('pactual_emerge_weighted','sum')).reset_index()
pal = sns.color_palette("tab20c")
fig, ax = plt.subplots(2,sharex=True)
axes = [ax[0], ax[0].twinx(), ax[1], ax[1].twinx()]
microbe_stacked.plot.area(ax = axes[0],stacked=True,legend=False, linewidth=0,color=pal,sort_columns=True)
axes[0].set_ylabel(ylabel ='Microbial Strain Abundances',labelpad=15,fontsize=15)
axes[0].set_xlabel(xlabel = 'Time t',fontsize=15)
axes[0].ticklabel_format(axis = 'y',style='sci',scilimits=(0,0))
axes[0].xaxis.set_minor_locator(ticker.MultipleLocator(25))
axes[0].tick_params(axis='x', labelsize= 10)
axes[0].tick_params(axis='y', labelsize= 10)
lim = axes[0].get_ylim()
axes[0].set_ylim(0,lim[1])
axes[1].plot(virus_total['t'],virus_total['vtotal'],linewidth=0,color='grey')
axes[1].fill_between(virus_total['t'],virus_total['vtotal'], color='grey',alpha=0.6)
lim = axes[1].get_ylim()
axes[1].set_ylim(0,lim[1])
axes[2].plot(pEmergeExpected['t'],pEmergeExpected['p_exp'],linewidth=1,color='darkblue')
axes[2].set_ylabel(ylabel ='Probability of Viral Strain Emergence',labelpad=15,fontsize=15)
axes[2].set_xlabel(xlabel = 'Time t',fontsize=15)
axes[2].tick_params(axis='x', labelsize= 10)
axes[2].tick_params(axis='y', labelsize= 10)
axes[3].plot(virus_total['t'],virus_total['vtotal'],linewidth=0,color='grey')
axes[3].fill_between(virus_total['t'],virus_total['vtotal'], color='grey',alpha=0.3)
fig.tight_layout()


frequencies = pd.read_sql_query("SELECT t,match_id,bsusceptible \
        FROM vmatches_abundances", conTri)\
        .rename(columns={"match_id":"vmatch_id"})\
        .merge(microbe_total,on = ['t'])\
        .rename(columns={"bsusceptible":"bfreq"})
frequencies['bfreq'] = frequencies['bfreq']/frequencies['btotal']
norm = frequencies.groupby(['t']).agg(norm=('bfreq','sum')).reset_index()
frequencies = frequencies.merge(norm,on=['t'])
frequencies['bfreq'] = frequencies['bfreq']/frequencies['norm']
frequencies = frequencies.drop(columns=['btotal','norm'])
pEmergeExpected = pd.read_sql_query("SELECT t, vmatch_id, p_extinction_lambert, vfrequency \
    FROM existing_vmatch_extinction ORDER BY t", conProb)\
    .merge(frequencies,on=['t','vmatch_id']).rename(columns={"vfrequency":"vfreq"})
pEmergeExpected['vfreq'] = pEmergeExpected['vfreq']/pEmergeExpected['bfreq']
norm = pEmergeExpected.groupby(['t']).agg(norm=('vfreq','sum')).reset_index()
pEmergeExpected = pEmergeExpected.merge(norm,on=['t'])
pEmergeExpected['vfreq'] = pEmergeExpected['vfreq']/pEmergeExpected['norm']
pEmergeExpected = pEmergeExpected.drop(columns=['norm'])
pEmergeExpected['p_extinction_lambert'] = 1- np.array(pEmergeExpected['p_extinction_lambert'])
pEmergeExpected=pEmergeExpected.rename(columns={"p_extinction_lambert":"plambert_emerge_weighted"})
pEmergeExpected['plambert_emerge_weighted'] = \
    pEmergeExpected['plambert_emerge_weighted']*pEmergeExpected['vfreq']
pEmergeExpected = pEmergeExpected.groupby(['t'])\
    .agg(p_exp=('plambert_emerge_weighted','sum')).reset_index()
# pEmergeExpected = pEmergeExpected.merge(microbe_total,on=['t'])
# pEmergeExpected['p_exp'] = pEmergeExpected['p_exp']/pEmergeExpected['btotal']
pal = sns.color_palette("tab20c")
fig, ax = plt.subplots(2,sharex=True)
axes = [ax[0], ax[0].twinx(), ax[1], ax[1].twinx()]
microbe_stacked.plot.area(ax = axes[0],stacked=True,legend=False, linewidth=0,color=pal,sort_columns=True)
axes[0].set_ylabel(ylabel ='Microbial Strain Abundances',labelpad=15,fontsize=15)
axes[0].set_xlabel(xlabel = 'Time t',fontsize=15)
axes[0].ticklabel_format(axis = 'y',style='sci',scilimits=(0,0))
axes[0].xaxis.set_minor_locator(ticker.MultipleLocator(25))
axes[0].tick_params(axis='x', labelsize= 10)
axes[0].tick_params(axis='y', labelsize= 10)
lim = axes[0].get_ylim()
axes[0].set_ylim(0,lim[1])
axes[1].plot(virus_total['t'],virus_total['vtotal'],linewidth=0,color='grey')
axes[1].fill_between(virus_total['t'],virus_total['vtotal'], color='grey',alpha=0.6)
lim = axes[1].get_ylim()
axes[1].set_ylim(0,lim[1])
axes[2].plot(pEmergeExpected['t'],pEmergeExpected['p_exp'],linewidth=1,color='darkblue')
axes[2].set_ylabel(ylabel ='Probability of Viral Strain Emergence',labelpad=15,fontsize=15)
axes[2].set_xlabel(xlabel = 'Time t',fontsize=15)
axes[2].tick_params(axis='x', labelsize= 10)
axes[2].tick_params(axis='y', labelsize= 10)
axes[3].plot(virus_total['t'],virus_total['vtotal'],linewidth=0,color='grey')
axes[3].fill_between(virus_total['t'],virus_total['vtotal'], color='grey',alpha=0.3)
fig.tight_layout()






print('Ordering spacers by time acquired...')
bstrainTimes = pd.read_sql_query("SELECT t_creation,bstrain_id \
    FROM bstrains WHERE run_id = {}".format(run_id), conSim)
bstrainExactTimes = pd.read_sql_query("SELECT t_creation,bstrain_id \
    FROM bstrains WHERE run_id = {}".format(run_id), conSim)
bstrainTimes.t_creation = list(map(math.ceil,bstrainTimes.t_creation))
bstrainTimes = bstrainTimes.rename(columns={"t_creation": "t"})
bmatches = pd.read_sql_query("SELECT t, match_id,bstrain_id \
    FROM bmatches", conTri)
bmatchTimes = pd.merge(bstrainTimes, bmatches, on=["t","bstrain_id"])
bmatchTimes = bmatchTimes.drop(columns=['t'])
bmatchTimesExact = pd.merge(bstrainExactTimes, bmatchTimes, on=["bstrain_id"])
bmatchPhenos = pd.read_sql_query("SELECT match_id, phenotype \
    FROM bmatch_phenotypes WHERE match_id in ({0})"
    .format(', '.join(map(str,np.unique(bmatchTimesExact.match_id)))), conTri)
spacersOrdered = []
for matchID in bmatchTimesExact.match_id:
    phenotype = bmatchPhenos[bmatchPhenos.match_id == matchID].phenotype.values
    newSpacer = list(set(phenotype) - set(spacersOrdered))
    if len(newSpacer) > 1:
        # print('oops')
        break
    spacersOrdered.extend(newSpacer)
    # print(phenotype)


newSpacerIDs = dict(zip(spacersOrdered,range(1,len(spacersOrdered)+1)))


print('SQLite Query: tripartite network data')
tripartiteNets = pd.read_sql_query("SELECT t, vmatch_id, spacer_id, bmatch_id \
    FROM single_match_tripartite_networks", conTri)
# vmatchspacers = pd.read_sql_query("SELECT t, vmatch_id, spacer_id \
#     FROM single_match_tripartite_networks", conTri)
# numSpacers = vmatchspacers.groupby(['t','vmatch_id']).agg(
#      numSpacers = ('spacer_id','size')
#      ).reset_index()
microbe_stacked = microbe_stacked.pivot(index='t',columns='bstrain_id',values='abundance')
vmatchstrain = pd.read_sql_query("SELECT t, match_id, vstrain_id \
    FROM vmatches", conTri)\
    .rename(columns={"match_id":"vmatch_id"})
pEmergence = pd.read_sql_query("SELECT t, vmatch_id, p_extinction \
    FROM vmatch_extinction ORDER BY t, vmatch_id", conProb)
pEmergence['p_emergence'] = 1 - np.array(pEmergence['p_extinction'])
pEmergence = pEmergence.drop(columns=['p_extinction'])
pEmergence = pEmergence.merge(vmatchstrain,on=['t','vmatch_id'])
bmatchstrain = pd.read_sql_query("SELECT t, match_id, bstrain_id \
    FROM bmatches", conTri)\
    .rename(columns={"match_id":"bmatch_id"})

v0bmatches = pd.read_sql_query("SELECT t, bstrain_id, vstrain_id \
                FROM bstrain_to_vstrain_0matches", conMatch)


###############################################################################
###############################################################################
###############################################################################
# V INFECTIOUS TO SPACER IDs WITH P_EMERGENCE AND TRIMICROBES AND NO VIRAL ABUNDANCE
spacerHeatPB_noV = tripartiteNets.drop(columns=['vmatch_id']).drop_duplicates()\
    .merge(bmatchstrain,on=['t','bmatch_id']).drop(columns=['bmatch_id'])\
    .merge(v0bmatches,on=['t','bstrain_id']).drop(columns=['bstrain_id'])\
    .drop_duplicates().merge(pEmergence,on=['t','vstrain_id'])
norm = spacerHeatPB_noV.groupby(['t']).agg(norm=('p_emergence','sum')).reset_index()
spacerHeatPB_noV = spacerHeatPB_noV.merge(norm,on=['t'])
spacerHeatPB_noV['p_emergence'] =  \
    [spacerHeatPB_noV['p_emergence'][i]/spacerHeatPB_noV['norm'][i] if spacerHeatPB_noV['norm'][i] != 0 else 0 \
    for i in range(0,len(spacerHeatPB_noV['p_emergence']))]
spacerHeatPB_noV = spacerHeatPB_noV.groupby(['t','spacer_id'])\
            .agg(p_spacer=('p_emergence','sum')).reset_index()
bfreq = tripartiteNets.drop(columns=['vmatch_id']).drop_duplicates()\
.merge(\
pd.read_sql_query("SELECT t,match_id,babundance \
    FROM bmatches_abundances",conTri).rename(columns={"match_id":"bmatch_id"}),\
on=['t','bmatch_id']).groupby(['t','spacer_id']).agg(bfreq=('babundance','sum'))\
.reset_index()
bfreq = bfreq.groupby(['t']).agg(btotal=('bfreq','sum')).reset_index()\
        .merge(bfreq,on=['t'])
bfreq['bfreq'] = bfreq['bfreq']/bfreq['btotal']
bfreq = bfreq.drop(columns=['btotal'])
spacerHeatPB_noV = spacerHeatPB_noV.merge(bfreq,on=['t','spacer_id'])
spacerHeatPB_noV['exp_freq'] = spacerHeatPB_noV['p_spacer']*spacerHeatPB_noV['bfreq']
expBFreq = spacerHeatPB_noV.drop(columns=['p_spacer','bfreq']).groupby(['t'])\
            .agg(exp_freq=('exp_freq','sum')).reset_index()
expectationB = expBFreq.groupby(['t']).agg(exp_bfreq=('exp_freq','sum')).reset_index()
shannon = []
for t in sorted(spacerHeatPB_noV['t'].unique()):
    div = list(spacerHeatPB_noV[(spacerHeatPB_noV['t']==t) & (spacerHeatPB_noV['exp_freq']!=0)]['exp_freq'])
    shannon.append(np.exp(-1*sum(np.array(div)*np.array(list(map(np.log,div))))))


expBFreq['shanDivBurn'] = shannon[:]
shannon = []
for t in sorted(spacerHeatPB_noV['t'].unique()):
    div = list(bfreq[(bfreq['t']==t)]['bfreq'])
    shannon.append(np.exp(-1*sum(np.array(div)*np.array(list(map(np.log,div))))))


expBFreq['shanDiv'] = shannon[:]

expBFreq['shanDivRatio'] = np.array(expBFreq['shanDivBurn'])/np.array(expBFreq['shanDiv'])


# spacerHeatPB_noV = spacerHeatPB_noV[spacerHeatPB_noV.exp_freq !=0 ]
# spacerHeatPB_noV['exp_freq'] = list(map(np.log,spacerHeatPB_noV.exp_freq))
spacerHeatPB_noV = spacerHeatPB_noV.replace({"spacer_id": newSpacerIDs})\
            .pivot_table(index='spacer_id', columns='t', values='exp_freq').fillna(0)

# SPACER DIVERSITY PLOT
pal = sns.color_palette("tab20c")
fig, ax = plt.subplots(2,sharex=True)
axes = [ax[0], ax[0].twinx(), ax[1], ax[1].twinx()]
microbe_stacked.plot.area(ax = axes[0],stacked=True,legend=False, linewidth=0,color=pal,sort_columns=True)
axes[0].set_ylabel(ylabel ='Microbial Strain Abundances',labelpad=15,fontsize=15)
axes[0].set_xlabel(xlabel = 'Time t',fontsize=15)
axes[0].ticklabel_format(axis = 'y',style='sci',scilimits=(0,0))
axes[0].xaxis.set_minor_locator(ticker.MultipleLocator(25))
axes[0].tick_params(axis='x', labelsize= 10)
axes[0].tick_params(axis='y', labelsize= 10)
lim = axes[0].get_ylim()
axes[0].set_ylim(0,lim[1])
axes[1].plot(virus_total['t'],virus_total['vtotal'],linewidth=0,color='grey')
axes[1].fill_between(virus_total['t'],virus_total['vtotal'], color='grey',alpha=0.6)
lim = axes[1].get_ylim()
axes[1].set_ylim(0,lim[1])
axes[2].plot(expBFreq['t'],expBFreq['shanDiv'],linewidth=1,color='darkblue')
axes[2].set_ylabel(ylabel ='Spacer Shannon Diversity',labelpad=15,fontsize=15)
axes[2].set_xlabel(xlabel = 'Time t',fontsize=15)
axes[2].tick_params(axis='x', labelsize= 10)
axes[2].tick_params(axis='y', labelsize= 10)
axes[3].plot(virus_total['t'],virus_total['vtotal'],linewidth=0,color='grey')
axes[3].fill_between(virus_total['t'],virus_total['vtotal'], color='grey',alpha=0.3)
fig.tight_layout()


# TARGETED SPACER DIVERSITY PLOT
pal = sns.color_palette("tab20c")
fig, ax = plt.subplots(2,sharex=True)
axes = [ax[0], ax[0].twinx(), ax[1], ax[1].twinx()]
microbe_stacked.plot.area(ax = axes[0],stacked=True,legend=False, linewidth=0,color=pal,sort_columns=True)
axes[0].set_ylabel(ylabel ='Microbial Strain Abundances',labelpad=15,fontsize=15)
axes[0].set_xlabel(xlabel = 'Time t',fontsize=15)
axes[0].ticklabel_format(axis = 'y',style='sci',scilimits=(0,0))
axes[0].xaxis.set_minor_locator(ticker.MultipleLocator(25))
axes[0].tick_params(axis='x', labelsize= 10)
axes[0].tick_params(axis='y', labelsize= 10)
lim = axes[0].get_ylim()
axes[0].set_ylim(0,lim[1])
axes[1].plot(virus_total['t'],virus_total['vtotal'],linewidth=0,color='grey')
axes[1].fill_between(virus_total['t'],virus_total['vtotal'], color='grey',alpha=0.6)
lim = axes[1].get_ylim()
axes[1].set_ylim(0,lim[1])
axes[2].plot(expBFreq['t'],expBFreq['shanDivRatio'],linewidth=1,color='darkblue')
axes[2].set_ylabel(ylabel ='Proportion of Spacers Targeted',labelpad=15,fontsize=15)
axes[2].set_xlabel(xlabel = 'Time t',fontsize=15)
axes[2].tick_params(axis='x', labelsize= 10)
axes[2].tick_params(axis='y', labelsize= 10)
axes[3].plot(virus_total['t'],virus_total['vtotal'],linewidth=0,color='grey')
axes[3].fill_between(virus_total['t'],virus_total['vtotal'], color='grey',alpha=0.3)
fig.tight_layout()

# EXPECTATION OF B_FREQUENCY PLOT
pal = sns.color_palette("tab20c")
fig, ax = plt.subplots(2,sharex=True)
axes = [ax[0], ax[0].twinx(), ax[1], ax[1].twinx()]
microbe_stacked.plot.area(ax = axes[0],stacked=True,legend=False, linewidth=0,color=pal,sort_columns=True)
axes[0].set_ylabel(ylabel ='Microbial Strain Abundances',labelpad=15,fontsize=15)
axes[0].set_xlabel(xlabel = 'Time t',fontsize=15)
axes[0].ticklabel_format(axis = 'y',style='sci',scilimits=(0,0))
axes[0].xaxis.set_minor_locator(ticker.MultipleLocator(25))
axes[0].tick_params(axis='x', labelsize= 10)
axes[0].tick_params(axis='y', labelsize= 10)
lim = axes[0].get_ylim()
axes[0].set_ylim(0,lim[1])
axes[1].plot(virus_total['t'],virus_total['vtotal'],linewidth=0,color='grey')
axes[1].fill_between(virus_total['t'],virus_total['vtotal'], color='grey',alpha=0.6)
lim = axes[1].get_ylim()
axes[1].set_ylim(0,lim[1])
axes[2].plot(expectationB['t'],expectationB['exp_bfreq'],linewidth=1,color='darkblue')
axes[2].set_ylabel(ylabel ='Expected Microbial Relative Abundance',labelpad=15,fontsize=15)
axes[2].set_xlabel(xlabel = 'Time t',fontsize=15)
axes[2].tick_params(axis='x', labelsize= 10)
axes[2].tick_params(axis='y', labelsize= 10)
axes[3].plot(virus_total['t'],virus_total['vtotal'],linewidth=0,color='grey')
axes[3].fill_between(virus_total['t'],virus_total['vtotal'], color='grey',alpha=0.3)
fig.tight_layout()



# EVERYTHING COMPARED
pal = sns.color_palette("tab20c")
fig, ax = plt.subplots(5,sharex=True)
axes = [ax[0], ax[0].twinx(), ax[1], ax[2], ax[3], ax[4], ax[4].twinx(),ax[1].twinx(),ax[2].twinx(),ax[3].twinx()]
microbe_stacked.plot.area(ax = axes[0],stacked=True,legend=False, linewidth=0,color=pal,sort_columns=True)
axes[0].set_ylabel(ylabel ='Microbial Strain Abundances',labelpad=15,fontsize=7)
axes[0].ticklabel_format(axis = 'y',style='sci',scilimits=(0,0))
axes[0].xaxis.set_minor_locator(ticker.MultipleLocator(25))
lim = axes[0].get_ylim()
axes[0].set_ylim(0,lim[1])
axes[1].plot(virus_total[virus_total.t<=max(microbe_stacked.index)]['t'],virus_total[virus_total.t<=max(microbe_stacked.index)]['vtotal'],linewidth=0,color='grey')
axes[1].fill_between(virus_total[virus_total.t<=max(microbe_stacked.index)]['t'],virus_total[virus_total.t<=max(microbe_stacked.index)]['vtotal'], color='grey',alpha=0.6)
lim = axes[1].get_ylim()
axes[1].set_ylim(0,lim[1])
axes[2].plot(expBFreq['t'],expBFreq['exp_freq'], color = 'darkorange')
axes[2].set_ylabel(ylabel ='Expected Burn Frequency',labelpad=15,fontsize=7)
# axes[3].plot(expBFreq['t'],expBFreq['displaced'])
axes[3].plot(expBFreq['t'],expBFreq['shanDivBurn'], color = 'darkgreen')
axes[3].set_ylabel(ylabel ='Targeted Spacer Diversity',labelpad=15,fontsize=7)
axes[4].plot(expBFreq['t'],expBFreq['shanDiv'], color='darkred')
axes[4].set_ylabel(ylabel ='Spacer Shannon Diversity',labelpad=15,fontsize=7)
axes[5].plot(expBFreq['t'],expBFreq['shanDivRatio'],color='darkblue')
axes[5].set_ylabel(ylabel ='Proportion of Spacers Targeted',labelpad=15,fontsize=7)
axes[6].plot(virus_total[virus_total.t<=max(microbe_stacked.index)]['t'],virus_total[virus_total.t<=max(microbe_stacked.index)]['vtotal'],linewidth=0,color='grey')
axes[6].fill_between(virus_total[virus_total.t<=max(microbe_stacked.index)]['t'],virus_total[virus_total.t<=max(microbe_stacked.index)]['vtotal'], color='grey',alpha=0.5)
lim = axes[1].get_ylim()
axes[6].set_ylim(0,lim[1])
axes[5].set_xlabel(xlabel = 'Time t',fontsize=7)
axes[7].plot(virus_total[virus_total.t<=max(microbe_stacked.index)]['t'],virus_total[virus_total.t<=max(microbe_stacked.index)]['vtotal'],linewidth=0,color='grey')
axes[7].fill_between(virus_total[virus_total.t<=max(microbe_stacked.index)]['t'],virus_total[virus_total.t<=max(microbe_stacked.index)]['vtotal'], color='grey',alpha=0.5)
axes[8].plot(virus_total[virus_total.t<=max(microbe_stacked.index)]['t'],virus_total[virus_total.t<=max(microbe_stacked.index)]['vtotal'],linewidth=0,color='grey')
axes[8].fill_between(virus_total[virus_total.t<=max(microbe_stacked.index)]['t'],virus_total[virus_total.t<=max(microbe_stacked.index)]['vtotal'], color='grey',alpha=0.5)
axes[9].plot(virus_total[virus_total.t<=max(microbe_stacked.index)]['t'],virus_total[virus_total.t<=max(microbe_stacked.index)]['vtotal'],linewidth=0,color='grey')
axes[9].fill_between(virus_total[virus_total.t<=max(microbe_stacked.index)]['t'],virus_total[virus_total.t<=max(microbe_stacked.index)]['vtotal'], color='grey',alpha=0.5)
fig.tight_layout()




# SPACER HEAT MAP WITH D_TARGET
pal = sns.color_palette("tab20c")
fig, ax = plt.subplots(3,sharex=True, gridspec_kw={'height_ratios': [1,4,1]})
axes = [ax[0], ax[0].twinx(), ax[1], ax[2], ax[2].twinx()]#, ax[2], ax[3]]
microbe_stacked.plot.area(ax = axes[0],stacked=True,legend=False, linewidth=0,color=pal,sort_columns=True)
axes[0].set_ylabel(ylabel ='Microbial Strain Abundances',labelpad=15,fontsize=7)
axes[0].set_xlabel(xlabel = 'Time t',fontsize=7)
axes[0].ticklabel_format(axis = 'y',style='sci',scilimits=(0,0))
axes[0].xaxis.set_minor_locator(ticker.MultipleLocator(25))
lim = axes[0].get_ylim()
axes[0].set_ylim(0,lim[1])
axes[1].plot(virus_total[virus_total.t<=max(microbe_stacked.index)]['t'],virus_total[virus_total.t<=max(microbe_stacked.index)]['vtotal'],linewidth=0,color='grey')
axes[1].fill_between(virus_total[virus_total.t<=max(microbe_stacked.index)]['t'],virus_total[virus_total.t<=max(microbe_stacked.index)]['vtotal'], color='grey',alpha=0.6)
lim = axes[1].get_ylim()
axes[1].set_ylim(0,lim[1])
axes[2].imshow(spacerHeatPB_noV,cmap='jet', aspect='auto',extent=[0, max(microbe_stacked.index), max(spacerHeatPB_noV.index), min(spacerHeatPB_noV.index)-1])
axes[2].set_yticks(np.arange(0.5,max(spacerHeatPB_noV.index)+0.5,1))
axes[2].set_yticklabels(spacersOrdered,fontsize=4)
axes[2].set_ylabel(ylabel ='Spacer IDs',labelpad=15,fontsize=10)
fig.tight_layout()

# axes[1].grid(which='major', color='w', linestyle='-', linewidth=.03)
# ax.set_yticks(np.arange(2, 10, 1), minor=True)
# axes[1].yaxis.grid(True)
# # axes[2].plot(expBFreq['t'],expBFreq['exp_freq'])
axes[3].plot(expBFreq['t'],expBFreq['shanDivRatio'],color='darkblue')
axes[3].set_ylabel(ylabel ='Proportion of Spacers Targeted',labelpad=15,fontsize=7)
axes[4].plot(virus_total[virus_total.t<=max(microbe_stacked.index)]['t'],virus_total[virus_total.t<=max(microbe_stacked.index)]['vtotal'],linewidth=0,color='grey')
axes[4].fill_between(virus_total[virus_total.t<=max(microbe_stacked.index)]['t'],virus_total[virus_total.t<=max(microbe_stacked.index)]['vtotal'], color='grey',alpha=0.5)



# SPACER HEAT MAP WITH EXPECTATION OF B_FREQUENCY
pal = sns.color_palette("tab20c")
fig, ax = plt.subplots(3,sharex=True, gridspec_kw={'height_ratios': [1,4,1]})
axes = [ax[0], ax[0].twinx(), ax[1], ax[2], ax[2].twinx()]#, ax[2], ax[3]]
microbe_stacked.plot.area(ax = axes[0],stacked=True,legend=False, linewidth=0,color=pal,sort_columns=True)
axes[0].set_ylabel(ylabel ='Microbial Strain Abundances',labelpad=15,fontsize=7)
axes[0].set_xlabel(xlabel = 'Time t',fontsize=7)
axes[0].ticklabel_format(axis = 'y',style='sci',scilimits=(0,0))
axes[0].xaxis.set_minor_locator(ticker.MultipleLocator(25))
lim = axes[0].get_ylim()
axes[0].set_ylim(0,lim[1])
axes[1].plot(virus_total[virus_total.t<=max(microbe_stacked.index)]['t'],virus_total[virus_total.t<=max(microbe_stacked.index)]['vtotal'],linewidth=0,color='grey')
axes[1].fill_between(virus_total[virus_total.t<=max(microbe_stacked.index)]['t'],virus_total[virus_total.t<=max(microbe_stacked.index)]['vtotal'], color='grey',alpha=0.6)
lim = axes[1].get_ylim()
axes[1].set_ylim(0,lim[1])
axes[2].imshow(spacerHeatPB_noV,cmap='jet', aspect='auto',extent=[0, max(microbe_stacked.index), max(spacerHeatPB_noV.index), min(spacerHeatPB_noV.index)-1])
axes[2].set_yticks(np.arange(0.5,max(spacerHeatPB_noV.index)+0.5,1))
axes[2].set_yticklabels(spacersOrdered,fontsize=4)
axes[2].set_ylabel(ylabel ='Spacer IDs',labelpad=15,fontsize=10)
fig.tight_layout()

# axes[1].grid(which='major', color='w', linestyle='-', linewidth=.03)
# ax.set_yticks(np.arange(2, 10, 1), minor=True)
# axes[1].yaxis.grid(True)
# # axes[2].plot(expBFreq['t'],expBFreq['exp_freq'])
axes[3].plot(expectationB['t'],expectationB['exp_bfreq'],color='darkblue')
axes[3].set_ylabel(ylabel ='Expected Microbial Relative Abundance,labelpad=15,fontsize=7)
axes[4].plot(virus_total[virus_total.t<=max(microbe_stacked.index)]['t'],virus_total[virus_total.t<=max(microbe_stacked.index)]['vtotal'],linewidth=0,color='grey')
axes[4].fill_between(virus_total[virus_total.t<=max(microbe_stacked.index)]['t'],virus_total[virus_total.t<=max(microbe_stacked.index)]['vtotal'], color='grey',alpha=0.5)









###############################################################################
###############################################################################
###############################################################################
# V INFECTIOUS TO SPACER IDs WITH P_EMERGENCE AND TRIMICROBES WITHHHHH VIRAL ABUNDANCE
spacerHeatPB_noV = tripartiteNets.drop(columns=['vmatch_id']).drop_duplicates()\
    .merge(bmatchstrain,on=['t','bmatch_id']).drop(columns=['bmatch_id'])\
    .merge(v0bmatches,on=['t','bstrain_id']).drop(columns=['bstrain_id'])\
    .drop_duplicates().merge(pEmergence,on=['t','vstrain_id'])
norm = spacerHeatPB_noV.groupby(['t']).agg(norm=('p_emergence','sum')).reset_index()
spacerHeatPB_noV = spacerHeatPB_noV.merge(norm,on=['t'])
spacerHeatPB_noV['p_emergence'] =  \
    [spacerHeatPB_noV['p_emergence'][i]/spacerHeatPB_noV['norm'][i] if spacerHeatPB_noV['norm'][i] != 0 else 0 \
    for i in range(0,len(spacerHeatPB_noV['p_emergence']))]
spacerHeatPB_noV = spacerHeatPB_noV.groupby(['t','spacer_id'])\
            .agg(p_spacer=('p_emergence','sum')).reset_index()
bfreq = tripartiteNets.drop(columns=['vmatch_id']).drop_duplicates()\
.merge(\
pd.read_sql_query("SELECT t,match_id,babundance \
    FROM bmatches_abundances",conTri).rename(columns={"match_id":"bmatch_id"}),\
on=['t','bmatch_id']).groupby(['t','spacer_id']).agg(bfreq=('babundance','sum'))\
.reset_index()
bfreq = bfreq.groupby(['t']).agg(btotal=('bfreq','sum')).reset_index()\
        .merge(bfreq,on=['t'])
bfreq['bfreq'] = bfreq['bfreq']/bfreq['btotal']
bfreq = bfreq.drop(columns=['btotal'])
spacerHeatPB_noV = spacerHeatPB_noV.merge(bfreq,on=['t','spacer_id'])
spacerHeatPB_noV['exp_freq'] = spacerHeatPB_noV['p_spacer']*spacerHeatPB_noV['bfreq']
expBFreq = spacerHeatPB_noV.drop(columns=['p_spacer','bfreq']).groupby(['t'])\
            .agg(exp_freq=('exp_freq','sum')).reset_index()
shannon = []
for t in sorted(spacerHeatPB_noV['t'].unique()):
    div = list(spacerHeatPB_noV[(spacerHeatPB_noV['t']==t) & (spacerHeatPB_noV['exp_freq']!=0)]['exp_freq'])
    shannon.append(np.exp(-1*sum(np.array(div)*np.array(list(map(np.log,div))))))


expBFreq['shanDivBurn'] = shannon[:]
shannon = []
for t in sorted(spacerHeatPB_noV['t'].unique()):
    div = list(bfreq[(bfreq['t']==t)]['bfreq'])
    shannon.append(np.exp(-1*sum(np.array(div)*np.array(list(map(np.log,div))))))


expBFreq['shanDiv'] = shannon[:]

expBFreq['shanDivRatio'] = np.array(expBFreq['shanDivBurn'])/np.array(expBFreq['shanDiv'])


# spacerHeatPB_noV = spacerHeatPB_noV[spacerHeatPB_noV.exp_freq !=0 ]
# spacerHeatPB_noV['exp_freq'] = list(map(np.log,spacerHeatPB_noV.exp_freq))
spacerHeatPB_noV = spacerHeatPB_noV.replace({"spacer_id": newSpacerIDs})\
            .pivot_table(index='spacer_id', columns='t', values='exp_freq').fillna(0)

# SPACER DIVERSITY PLOT
pal = sns.color_palette("tab20c")
fig, ax = plt.subplots(2,sharex=True)
axes = [ax[0], ax[0].twinx(), ax[1], ax[1].twinx()]
microbe_stacked.plot.area(ax = axes[0],stacked=True,legend=False, linewidth=0,color=pal,sort_columns=True)
axes[0].set_ylabel(ylabel ='Microbial Strain Abundances',labelpad=15,fontsize=15)
axes[0].set_xlabel(xlabel = 'Time t',fontsize=15)
axes[0].ticklabel_format(axis = 'y',style='sci',scilimits=(0,0))
axes[0].xaxis.set_minor_locator(ticker.MultipleLocator(25))
axes[0].tick_params(axis='x', labelsize= 10)
axes[0].tick_params(axis='y', labelsize= 10)
lim = axes[0].get_ylim()
axes[0].set_ylim(0,lim[1])
axes[1].plot(virus_total['t'],virus_total['vtotal'],linewidth=0,color='grey')
axes[1].fill_between(virus_total['t'],virus_total['vtotal'], color='grey',alpha=0.6)
lim = axes[1].get_ylim()
axes[1].set_ylim(0,lim[1])
axes[2].plot(expBFreq['t'],expBFreq['shanDiv'],linewidth=1,color='darkblue')
axes[2].set_ylabel(ylabel ='Spacer Shannon Diversity',labelpad=15,fontsize=15)
axes[2].set_xlabel(xlabel = 'Time t',fontsize=15)
axes[2].tick_params(axis='x', labelsize= 10)
axes[2].tick_params(axis='y', labelsize= 10)
axes[3].plot(virus_total['t'],virus_total['vtotal'],linewidth=0,color='grey')
axes[3].fill_between(virus_total['t'],virus_total['vtotal'], color='grey',alpha=0.3)
fig.tight_layout()


# TARGETED SPACER DIVERSITY PLOT
pal = sns.color_palette("tab20c")
fig, ax = plt.subplots(2,sharex=True)
axes = [ax[0], ax[0].twinx(), ax[1], ax[1].twinx()]
microbe_stacked.plot.area(ax = axes[0],stacked=True,legend=False, linewidth=0,color=pal,sort_columns=True)
axes[0].set_ylabel(ylabel ='Microbial Strain Abundances',labelpad=15,fontsize=15)
axes[0].set_xlabel(xlabel = 'Time t',fontsize=15)
axes[0].ticklabel_format(axis = 'y',style='sci',scilimits=(0,0))
axes[0].xaxis.set_minor_locator(ticker.MultipleLocator(25))
axes[0].tick_params(axis='x', labelsize= 10)
axes[0].tick_params(axis='y', labelsize= 10)
lim = axes[0].get_ylim()
axes[0].set_ylim(0,lim[1])
axes[1].plot(virus_total['t'],virus_total['vtotal'],linewidth=0,color='grey')
axes[1].fill_between(virus_total['t'],virus_total['vtotal'], color='grey',alpha=0.6)
lim = axes[1].get_ylim()
axes[1].set_ylim(0,lim[1])
axes[2].plot(expBFreq['t'],expBFreq['shanDivRatio'],linewidth=1,color='darkblue')
axes[2].set_ylabel(ylabel ='Proportion of Spacers Targeted',labelpad=15,fontsize=15)
axes[2].set_xlabel(xlabel = 'Time t',fontsize=15)
axes[2].tick_params(axis='x', labelsize= 10)
axes[2].tick_params(axis='y', labelsize= 10)
axes[3].plot(virus_total['t'],virus_total['vtotal'],linewidth=0,color='grey')
axes[3].fill_between(virus_total['t'],virus_total['vtotal'], color='grey',alpha=0.3)
fig.tight_layout()



# EVERYTHING COMPARED
pal = sns.color_palette("tab20c")
fig, ax = plt.subplots(5,sharex=True)
axes = [ax[0], ax[0].twinx(), ax[1], ax[2], ax[3], ax[4], ax[4].twinx(),ax[1].twinx(),ax[2].twinx(),ax[3].twinx()]
microbe_stacked.plot.area(ax = axes[0],stacked=True,legend=False, linewidth=0,color=pal,sort_columns=True)
axes[0].set_ylabel(ylabel ='Microbial Strain Abundances',labelpad=15,fontsize=7)
axes[0].ticklabel_format(axis = 'y',style='sci',scilimits=(0,0))
axes[0].xaxis.set_minor_locator(ticker.MultipleLocator(25))
lim = axes[0].get_ylim()
axes[0].set_ylim(0,lim[1])
axes[1].plot(virus_total[virus_total.t<=max(microbe_stacked.index)]['t'],virus_total[virus_total.t<=max(microbe_stacked.index)]['vtotal'],linewidth=0,color='grey')
axes[1].fill_between(virus_total[virus_total.t<=max(microbe_stacked.index)]['t'],virus_total[virus_total.t<=max(microbe_stacked.index)]['vtotal'], color='grey',alpha=0.6)
lim = axes[1].get_ylim()
axes[1].set_ylim(0,lim[1])
axes[2].plot(expBFreq['t'],expBFreq['exp_freq'], color = 'darkorange')
axes[2].set_ylabel(ylabel ='Expected Burn Frequency',labelpad=15,fontsize=7)
# axes[3].plot(expBFreq['t'],expBFreq['displaced'])
axes[3].plot(expBFreq['t'],expBFreq['shanDivBurn'], color = 'darkgreen')
axes[3].set_ylabel(ylabel ='Targeted Spacer Diversity',labelpad=15,fontsize=7)
axes[4].plot(expBFreq['t'],expBFreq['shanDiv'], color='darkred')
axes[4].set_ylabel(ylabel ='Spacer Shannon Diversity',labelpad=15,fontsize=7)
axes[5].plot(expBFreq['t'],expBFreq['shanDivRatio'],color='darkblue')
axes[5].set_ylabel(ylabel ='Proportion of Spacers Targeted',labelpad=15,fontsize=7)
axes[6].plot(virus_total[virus_total.t<=max(microbe_stacked.index)]['t'],virus_total[virus_total.t<=max(microbe_stacked.index)]['vtotal'],linewidth=0,color='grey')
axes[6].fill_between(virus_total[virus_total.t<=max(microbe_stacked.index)]['t'],virus_total[virus_total.t<=max(microbe_stacked.index)]['vtotal'], color='grey',alpha=0.5)
lim = axes[1].get_ylim()
axes[6].set_ylim(0,lim[1])
axes[5].set_xlabel(xlabel = 'Time t',fontsize=7)
axes[7].plot(virus_total[virus_total.t<=max(microbe_stacked.index)]['t'],virus_total[virus_total.t<=max(microbe_stacked.index)]['vtotal'],linewidth=0,color='grey')
axes[7].fill_between(virus_total[virus_total.t<=max(microbe_stacked.index)]['t'],virus_total[virus_total.t<=max(microbe_stacked.index)]['vtotal'], color='grey',alpha=0.5)
axes[8].plot(virus_total[virus_total.t<=max(microbe_stacked.index)]['t'],virus_total[virus_total.t<=max(microbe_stacked.index)]['vtotal'],linewidth=0,color='grey')
axes[8].fill_between(virus_total[virus_total.t<=max(microbe_stacked.index)]['t'],virus_total[virus_total.t<=max(microbe_stacked.index)]['vtotal'], color='grey',alpha=0.5)
axes[9].plot(virus_total[virus_total.t<=max(microbe_stacked.index)]['t'],virus_total[virus_total.t<=max(microbe_stacked.index)]['vtotal'],linewidth=0,color='grey')
axes[9].fill_between(virus_total[virus_total.t<=max(microbe_stacked.index)]['t'],virus_total[virus_total.t<=max(microbe_stacked.index)]['vtotal'], color='grey',alpha=0.5)
fig.tight_layout()




# SPACER HEAT MAP
pal = sns.color_palette("tab20c")
fig, ax = plt.subplots(3,sharex=True, gridspec_kw={'height_ratios': [1,4,1]})
axes = [ax[0], ax[0].twinx(), ax[1], ax[2], ax[2].twinx()]#, ax[2], ax[3]]
microbe_stacked.plot.area(ax = axes[0],stacked=True,legend=False, linewidth=0,color=pal,sort_columns=True)
axes[0].set_ylabel(ylabel ='Microbial Strain Abundances',labelpad=15,fontsize=7)
axes[0].set_xlabel(xlabel = 'Time t',fontsize=7)
axes[0].ticklabel_format(axis = 'y',style='sci',scilimits=(0,0))
axes[0].xaxis.set_minor_locator(ticker.MultipleLocator(25))
lim = axes[0].get_ylim()
axes[0].set_ylim(0,lim[1])
axes[1].plot(virus_total[virus_total.t<=max(microbe_stacked.index)]['t'],virus_total[virus_total.t<=max(microbe_stacked.index)]['vtotal'],linewidth=0,color='grey')
axes[1].fill_between(virus_total[virus_total.t<=max(microbe_stacked.index)]['t'],virus_total[virus_total.t<=max(microbe_stacked.index)]['vtotal'], color='grey',alpha=0.6)
lim = axes[1].get_ylim()
axes[1].set_ylim(0,lim[1])
axes[2].imshow(spacerHeatPB_noV,cmap='jet', aspect='auto',extent=[0, max(microbe_stacked.index), max(spacerHeatPB_noV.index), min(spacerHeatPB_noV.index)-1])
axes[2].set_yticks(np.arange(0.5,max(spacerHeatPB_noV.index)+0.5,1))
axes[2].set_yticklabels(spacersOrdered,fontsize=4)
axes[2].set_ylabel(ylabel ='Spacer IDs',labelpad=15,fontsize=10)
fig.tight_layout()

# axes[1].grid(which='major', color='w', linestyle='-', linewidth=.03)
# ax.set_yticks(np.arange(2, 10, 1), minor=True)
# axes[1].yaxis.grid(True)
# # axes[2].plot(expBFreq['t'],expBFreq['exp_freq'])
axes[3].plot(expBFreq['t'],expBFreq['shanDivRatio'],color='darkblue')
axes[3].set_ylabel(ylabel ='Proportion of Spacers Targeted',labelpad=15,fontsize=7)
axes[4].plot(virus_total[virus_total.t<=max(microbe_stacked.index)]['t'],virus_total[virus_total.t<=max(microbe_stacked.index)]['vtotal'],linewidth=0,color='grey')
axes[4].fill_between(virus_total[virus_total.t<=max(microbe_stacked.index)]['t'],virus_total[virus_total.t<=max(microbe_stacked.index)]['vtotal'], color='grey',alpha=0.5)


###############################################################################
###############################################################################
###############################################################################















# THIS IS FOR BURN WITH RESPECT TO STRAIN ID, NOT SPACER ID: LOG P_EMERGENCE AND TRIMICROBES
PB_noV = tripartiteNets.drop(columns=['vmatch_id','spacer_id']).drop_duplicates()\
    .merge(bmatchstrain,on=['t','bmatch_id']).drop(columns=['bmatch_id'])\
    .merge(v0bmatches,on=['t','bstrain_id']) \
    .drop_duplicates().merge(pEmergence,on=['t','vstrain_id'])\
    .drop(columns=['vmatch_id'])
norm = PB_noV.groupby(['t']).agg(norm=('p_emergence','sum')).reset_index()
PB_noV = PB_noV.merge(norm,on=['t'])
PB_noV['p_emergence'] =  \
    [PB_noV['p_emergence'][i]/PB_noV['norm'][i] if PB_noV['norm'][i] != 0 else 0 \
    for i in range(0,len(PB_noV['p_emergence']))]
PB_noV = PB_noV.groupby(['t','bstrain_id'])\
            .agg(p_bstrainBurn=('p_emergence','sum')).reset_index()

bfreq = tripartiteNets.drop(columns=['vmatch_id']).drop_duplicates()\
.merge(\
pd.read_sql_query("SELECT t,match_id,babundance \
    FROM bmatches_abundances",conTri).rename(columns={"match_id":"bmatch_id"}),\
on=['t','bmatch_id'])\
.merge(bmatchstrain,on=['t','bmatch_id']).drop(columns=['bmatch_id','spacer_id'])\
.drop_duplicates().rename(columns={"babundance":"bfreq"})

bfreq = bfreq.groupby(['t']).agg(btotal=('bfreq','sum')).reset_index()\
        .merge(bfreq,on=['t'])

bfreq['bfreq'] = bfreq['bfreq']/bfreq['btotal']
bfreq['bburn'] = bfreq['bfreq']/CC
bfreq = bfreq.drop(columns=['btotal'])

PB_noV = PB_noV.merge(bfreq,on=['t','bstrain_id'])
PB_noV['exp_freq'] = PB_noV['p_bstrainBurn']*PB_noV['bfreq']

expBFreq = PB_noV.groupby(['t']).agg(exp_freq=('exp_freq','sum')).reset_index()
shannon = []
for t in sorted(PB_noV['t'].unique()):
    div = list(PB_noV[(PB_noV['t']==t) & (PB_noV['exp_freq']!=0)]['exp_freq'])
    shannon.append(np.exp(-1*sum(np.array(div)*np.array(list(map(np.log,div))))))


expBFreq['shanDivBurn'] = shannon[:]

shannon = []
for t in sorted(PB_noV['t'].unique()):
    div = list(bfreq[(bfreq['t']==t)]['bfreq'])
    shannon.append(np.exp(-1*sum(np.array(div)*np.array(list(map(np.log,div))))))

expBFreq['shanDiv'] = shannon[:]

expBFreq['shanDivRatio'] = np.array(expBFreq['shanDivBurn'])/np.array(expBFreq['shanDiv'])

# spacerHeatPB_noV = spacerHeatPB_noV[spacerHeatPB_noV.exp_freq !=0 ]
# spacerHeatPB_noV['exp_freq'] = list(map(np.log,spacerHeatPB_noV.exp_freq))
PB_noV = PB_noV.pivot_table(index='bstrain_id', columns='t', values='exp_freq').fillna(0)

pal = sns.color_palette("tab20b")
fig, ax = plt.subplots(5,sharex=True)
axes = [ax[0], ax[0].twinx(), ax[1], ax[2], ax[3], ax[4], ax[4].twinx()]
microbe_stacked.plot.area(ax = axes[0],stacked=True,legend=False, linewidth=0,color=pal,sort_columns=True)
axes[0].set_ylabel(ylabel ='Microbial Strain Abundances',labelpad=15,fontsize=7)
axes[0].ticklabel_format(axis = 'y',style='sci',scilimits=(0,0))
axes[0].xaxis.set_minor_locator(ticker.MultipleLocator(25))
lim = axes[0].get_ylim()
axes[0].set_ylim(0,lim[1])
axes[1].plot(virus_total[virus_total.t<=max(microbe_stacked.index)]['t'],virus_total[virus_total.t<=max(microbe_stacked.index)]['vtotal'],linewidth=0,color='grey')
axes[1].fill_between(virus_total[virus_total.t<=max(microbe_stacked.index)]['t'],virus_total[virus_total.t<=max(microbe_stacked.index)]['vtotal'], color='grey',alpha=0.6)
lim = axes[1].get_ylim()
axes[1].set_ylim(0,lim[1])
axes[2].plot(expBFreq['t'],expBFreq['exp_freq'])
axes[2].set_ylabel(ylabel ='Expected Burn Frequency',labelpad=15,fontsize=7)
# axes[3].plot(expBFreq['t'],expBFreq['displaced'])
axes[3].plot(expBFreq['t'],expBFreq['shanDivBurn'])
axes[3].set_ylabel(ylabel ='Expected Spacer Burn Diversity',labelpad=15,fontsize=7)
axes[4].plot(expBFreq['t'],expBFreq['shanDiv'])
axes[4].set_ylabel(ylabel ='Spacer Shannon Diversity',labelpad=15,fontsize=7)
axes[5].plot(expBFreq['t'],expBFreq['shanDivRatio'])
axes[5].set_ylabel(ylabel ='R_burn',labelpad=15,fontsize=7)
axes[6].plot(virus_total[virus_total.t<=max(microbe_stacked.index)]['t'],virus_total[virus_total.t<=max(microbe_stacked.index)]['vtotal'],linewidth=0,color='grey')
axes[6].fill_between(virus_total[virus_total.t<=max(microbe_stacked.index)]['t'],virus_total[virus_total.t<=max(microbe_stacked.index)]['vtotal'], color='grey',alpha=0.6)
lim = axes[1].get_ylim()
axes[6].set_ylim(0,lim[1])
axes[5].set_xlabel(xlabel = 'Time t',fontsize=7)
fig.tight_layout()


























# V INFECTIOUS TO SPACER IDs
spacerHeat = tripartiteNets.drop(columns=['vmatch_id']).drop_duplicates()\
    .merge(bmatchstrain,on=['t','bmatch_id']).drop(columns=['bmatch_id'])\
    .merge(v0bmatches,on=['t','bstrain_id']).drop(columns=['bstrain_id'])\
    .drop_duplicates().merge(\
    pd.read_sql_query("SELECT t,vstrain_id,abundance \
        FROM vabundance WHERE run_id = {}".format(run_id), conSim)\
    ,on=['t','vstrain_id']).groupby(['t','spacer_id'])\
    .agg(vtotal=('abundance','sum')).reset_index()\
    .replace({"spacer_id": newSpacerIDs})\
    .pivot_table(index='spacer_id', columns='t', values='vtotal').fillna(0)
pal = sns.color_palette("tab20b")
fig, ax = plt.subplots(2,sharex=True)
axes = [ax[0], ax[1]]
microbe_stacked.plot.area(ax = axes[0],stacked=True,legend=False, linewidth=0,color=pal,sort_columns=True)
axes[0].set_ylabel(ylabel ='Microbial Strain Abundances',labelpad=15,fontsize=7)
axes[0].set_xlabel(xlabel = 'Time t',fontsize=7)
axes[0].ticklabel_format(axis = 'y',style='sci',scilimits=(0,0))
axes[0].xaxis.set_minor_locator(ticker.MultipleLocator(25))
lim = axes[0].get_ylim()
axes[0].set_ylim(0,lim[1])
axes[1].imshow(spacerHeat,cmap='viridis', aspect='auto')




# V INFECTIOUS TO SPACER IDs WITH P_EMERGENCE
spacerHeatP = tripartiteNets.drop(columns=['vmatch_id']).drop_duplicates()\
    .merge(bmatchstrain,on=['t','bmatch_id']).drop(columns=['bmatch_id'])\
    .merge(v0bmatches,on=['t','bstrain_id']).drop(columns=['bstrain_id'])\
    .drop_duplicates().merge(\
    pd.read_sql_query("SELECT t,vstrain_id,abundance \
        FROM vabundance WHERE run_id = {}".format(run_id), conSim)\
    ,on=['t','vstrain_id']).merge(pEmergence,on=['t','vstrain_id'])
spacerHeatP['p_emerge_abund'] = spacerHeatP['p_emergence']*spacerHeatP['abundance']
spacerHeatP = spacerHeatP.merge(virus_total,on=['t'])
spacerHeatP['p_emerge_abund'] = spacerHeatP['p_emerge_abund']/spacerHeatP['vtotal']
norm = spacerHeatP.groupby(['t']).agg(norm=('p_emerge_abund','sum')).reset_index()
spacerHeatP = spacerHeatP.merge(norm,on=['t'])
spacerHeatP['p_emerge_abund'] =  \
    [spacerHeatP['p_emerge_abund'][i]/spacerHeatP['norm'][i] if spacerHeatP['norm'][i] != 0 else 0 \
    for i in range(0,len(spacerHeatP['p_emerge_abund']))]
spacerHeatP = spacerHeatP.groupby(['t','spacer_id'])\
            .agg(p_spacer=('p_emerge_abund','sum')).reset_index()\
            .replace({"spacer_id": newSpacerIDs})\
            .pivot_table(index='spacer_id', columns='t', values='p_spacer').fillna(0)
pal = sns.color_palette("tab20b")
fig, ax = plt.subplots(2,sharex=True)
axes = [ax[0], ax[1]]
microbe_stacked.plot.area(ax = axes[0],stacked=True,legend=False, linewidth=0,color=pal,sort_columns=True)
axes[0].set_ylabel(ylabel ='Microbial Strain Abundances',labelpad=15,fontsize=7)
axes[0].set_xlabel(xlabel = 'Time t',fontsize=7)
axes[0].ticklabel_format(axis = 'y',style='sci',scilimits=(0,0))
axes[0].xaxis.set_minor_locator(ticker.MultipleLocator(25))
lim = axes[0].get_ylim()
axes[0].set_ylim(0,lim[1])
axes[1].imshow(spacerHeatP,cmap='viridis', aspect='auto')





# V INFECTIOUS TO SPACER IDs WITH P_EMERGENCE AND TRIMICROBES
spacerHeatPB = tripartiteNets.drop(columns=['vmatch_id']).drop_duplicates()\
    .merge(bmatchstrain,on=['t','bmatch_id']).drop(columns=['bmatch_id'])\
    .merge(v0bmatches,on=['t','bstrain_id']).drop(columns=['bstrain_id'])\
    .drop_duplicates().merge(\
    pd.read_sql_query("SELECT t,vstrain_id,abundance \
        FROM vabundance WHERE run_id = {}".format(run_id), conSim)\
    ,on=['t','vstrain_id']).merge(pEmergence,on=['t','vstrain_id'])
spacerHeatPB['p_emerge_abund'] = spacerHeatPB['p_emergence']*spacerHeatPB['abundance']
spacerHeatPB = spacerHeatPB.merge(virus_total,on=['t'])
spacerHeatPB['p_emerge_abund'] = spacerHeatPB['p_emerge_abund']/spacerHeatPB['vtotal']
norm = spacerHeatPB.groupby(['t']).agg(norm=('p_emerge_abund','sum')).reset_index()
spacerHeatPB = spacerHeatPB.merge(norm,on=['t'])
spacerHeatPB['p_emerge_abund'] =  \
    [spacerHeatPB['p_emerge_abund'][i]/spacerHeatPB['norm'][i] if spacerHeatPB['norm'][i] != 0 else 0 \
    for i in range(0,len(spacerHeatPB['p_emerge_abund']))]
spacerHeatPB = spacerHeatPB.groupby(['t','spacer_id'])\
            .agg(p_spacer=('p_emerge_abund','sum')).reset_index()
bfreq = tripartiteNets.drop(columns=['vmatch_id']).drop_duplicates()\
.merge(\
pd.read_sql_query("SELECT t,match_id,babundance \
    FROM bmatches_abundances",conTri).rename(columns={"match_id":"bmatch_id"}),\
on=['t','bmatch_id']).groupby(['t','spacer_id']).agg(bfreq=('babundance','sum'))\
.reset_index()
bfreq = bfreq.groupby(['t']).agg(btotal=('bfreq','sum')).reset_index()\
        .merge(bfreq,on=['t'])
bfreq['bfreq'] = bfreq['bfreq']/bfreq['btotal']
bfreq = bfreq.drop(columns=['btotal'])
spacerHeatPB = spacerHeatPB.merge(bfreq,on=['t','spacer_id'])
spacerHeatPB['exp_freq'] = spacerHeatPB['p_spacer']*spacerHeatPB['bfreq']
expBFreq = spacerHeatPB.drop(columns=['p_spacer','bfreq']).groupby(['t'])\
            .agg(exp_freq=('exp_freq','sum')).reset_index()
shannon = []
for t in sorted(spacerHeatPB['t'].unique()):
    div = list(spacerHeatPB[(spacerHeatPB['t']==t) & (spacerHeatPB['exp_freq']!=0)]['exp_freq'])
    shannon.append(np.exp(-1*sum(np.array(div)*np.array(list(map(np.log,div))))))

expBFreq['shanDiv'] = shannon[:]

# spacerHeatPB = spacerHeatPB[spacerHeatPB.exp_freq !=0 ]
# spacerHeatPB['exp_freq'] = list(map(np.log,spacerHeatPB.exp_freq))
spacerHeatPB = spacerHeatPB.replace({"spacer_id": newSpacerIDs})\
            .pivot_table(index='spacer_id', columns='t', values='exp_freq').fillna(0)

fig, ax = plt.subplots(2,sharex=True)
axes = [ax[0], ax[0].twinx(), ax[1]]#, ax[2], ax[3]]
microbe_stacked.plot.area(ax = axes[0],stacked=True,legend=False, linewidth=0,color=pal,sort_columns=True)
axes[0].set_ylabel(ylabel ='Microbial Strain Abundances',labelpad=15,fontsize=7)
axes[0].set_xlabel(xlabel = 'Time t',fontsize=7)
axes[0].ticklabel_format(axis = 'y',style='sci',scilimits=(0,0))
axes[0].xaxis.set_minor_locator(ticker.MultipleLocator(25))
lim = axes[0].get_ylim()
axes[0].set_ylim(0,lim[1])
axes[1].plot(virus_total[virus_total.t<=max(microbe_stacked.index)]['t'],virus_total[virus_total.t<=max(microbe_stacked.index)]['vtotal'],linewidth=0,color='grey')
axes[1].fill_between(virus_total['t'],virus_total['vtotal'], color='grey',alpha=0.6)
lim = axes[1].get_ylim()
axes[1].set_ylim(0,lim[1])
axes[2].imshow(spacerHeatPB,cmap='jet', aspect='auto',extent=[0, max(microbe_stacked.index), max(spacerHeatPB.index), min(spacerHeatPB.index)-1])
axes[2].set_yticks(np.arange(0.5,max(spacerHeatPB.index)+0.5,1))
axes[2].set_yticklabels(spacersOrdered,fontsize=4)
fig.tight_layout()
# axes[2].plot(expBFreq['t'],expBFreq['exp_freq'])
# axes[3].plot(expBFreq['t'],expBFreq['shanDiv'])


# JUST TRIMICROBES
bfreq = tripartiteNets.drop(columns=['vmatch_id']).drop_duplicates()\
.merge(\
pd.read_sql_query("SELECT t,match_id,babundance \
    FROM bmatches_abundances",conTri).rename(columns={"match_id":"bmatch_id"}),\
on=['t','bmatch_id']).groupby(['t','spacer_id']).agg(bfreq=('babundance','sum'))\
.reset_index()
bfreq = bfreq.groupby(['t']).agg(btotal=('bfreq','sum')).reset_index()\
        .merge(bfreq,on=['t'])
bfreq['bfreq'] = bfreq['bfreq']/bfreq['btotal']
shannon = []
times = []
for t in sorted(bfreq['t'].unique()):
    times.append(t)
    div = list(bfreq[(bfreq['t']==t) & (bfreq['bfreq']!=0)]['bfreq'])
    shannon.append(np.exp(-1*sum(np.array(div)*np.array(list(map(np.log,div))))))



bfreq = bfreq.drop(columns=['btotal'])\
    .replace({"spacer_id": newSpacerIDs})\
    .pivot_table(index='spacer_id', columns='t', values='bfreq').fillna(0)

pal = sns.color_palette("tab20b")
fig, ax = plt.subplots(3,sharex=True)
axes = [ax[0], ax[0].twinx(), ax[1], ax[2]]#, ax[2], ax[3]]
microbe_stacked.plot.area(ax = axes[0],stacked=True,legend=False, linewidth=0,color=pal,sort_columns=True)
axes[0].set_ylabel(ylabel ='Microbial Strain Abundances',labelpad=15,fontsize=7)
axes[0].set_xlabel(xlabel = 'Time t',fontsize=7)
axes[0].ticklabel_format(axis = 'y',style='sci',scilimits=(0,0))
axes[0].xaxis.set_minor_locator(ticker.MultipleLocator(25))
lim = axes[0].get_ylim()
axes[0].set_ylim(0,lim[1])
axes[1].plot(virus_total[virus_total.t<=max(microbe_stacked.index)]['t'],virus_total[virus_total.t<=max(microbe_stacked.index)]['vtotal'],linewidth=0,color='grey')
axes[1].fill_between(virus_total[virus_total.t<=max(microbe_stacked.index)]['t'],virus_total[virus_total.t<=max(microbe_stacked.index)]['vtotal'], color='grey',alpha=0.6)
lim = axes[1].get_ylim()
axes[1].set_ylim(0,lim[1])
axes[2].imshow(bfreq,cmap='jet', aspect='auto',extent=[0, max(bfreq.columns), max(bfreq.index), min(bfreq.index)-1])
axes[2].set_yticks(np.arange(0.5,max(bfreq.index)+0.5,1))
axes[2].set_yticklabels(spacersOrdered,fontsize=4)
axes[3].plot(times,shannon)
fig.tight_layout()









# pal = sns.color_palette("tab20b")
# vmatchstrain = pd.read_sql_query("SELECT t, match_id, vstrain_id \
#     FROM vmatches WHERE vstrain_id in ({0})"\
#     .format(', '.join(map(str,keepStrains))), conTri)\
#     .rename(columns={"match_id":"vmatch_id"})
# pEmergence = pd.read_sql_query("SELECT t, vmatch_id, p_extinction \
#     FROM vmatch_extinction ORDER BY t, vmatch_id", conProb)
# pEmergence['p_emergence'] = 1 - np.array(pEmergence['p_extinction'])
# pEmergence = pEmergence.drop(columns=['p_extinction'])
# pEmergence = pEmergence.merge(vmatchstrain,on=['t','vmatch_id'])
#
# # for t in pEmergence.t.unique()
# # for list(itertools.product([True,False], repeat=3)):
# #     pEmergence[(pEmergence.t == t) & ]
#
#
# # THIS PLOTS THE PREDICTION OF EXPECTED BURN FREQUENCY
# pEmergeVstrains = pEmergence.drop(columns=['vmatch_id'])
# pEmergeVstrains = pEmergeVstrains.merge(\
#                 pd.read_sql_query("SELECT t, bstrain_id, vstrain_id \
#                 FROM bstrain_to_vstrain_0matches", conMatch),on=['t','vstrain_id'])\
#                 .merge(microbe_stacked,on=['t','bstrain_id'])\
#                 .groupby(['t','p_emergence','vstrain_id']).agg(babundance=('abundance','sum'))\
#                 .reset_index()\
#                 .merge(microbe_total,on=['t'])
# pEmergeVstrains['burnNet'] = np.array(pEmergeVstrains['babundance'])/np.array(pEmergeVstrains['btotal'])
# pEmergeVstrains['burn'] = (CC - np.array(pEmergeVstrains['btotal'])\
#                             + np.array(pEmergeVstrains['babundance']))/CC
# pNorm = pEmergeVstrains.groupby(['t']).agg(pNorm=('p_emergence','sum')).reset_index()\
#         .merge(pEmergeVstrains[['t']],on=['t'])
# pNorm = pNorm['pNorm'].values
# pEmergeVstrains['p_emergence_normalized'] =  \
#     [pEmergeVstrains['p_emergence'][i]/pNorm[i] if pNorm[i] != 0 else 0 \
#     for i in range(0,len(pNorm))]
# pEmergeVstrains['burnProb'] = pEmergeVstrains['burn']* pEmergeVstrains['p_emergence_normalized']
# pEmergeVstrains['burnProbNet'] = pEmergeVstrains['burnNet']* pEmergeVstrains['p_emergence_normalized']
# pBurn = pEmergeVstrains.groupby(['t']).agg(expBurn=('burnProb','sum'),expBurnNet=('burnProbNet','sum')).reset_index()
#
#
# fig, ax = plt.subplots(2,sharex=True)
# axes = [ax[0], ax[1]]
# microbe_stacked = microbe_stacked.pivot(index='t',columns='bstrain_id',values='abundance')
# microbe_stacked.plot.area(ax = axes[0],stacked=True,legend=False, linewidth=0,color=pal,sort_columns=True)
# axes[0].set_ylabel(ylabel ='Microbial Strain Abundances',labelpad=15,fontsize=7)
# axes[0].set_xlabel(xlabel = 'Time t',fontsize=7)
# axes[0].ticklabel_format(axis = 'y',style='sci',scilimits=(0,0))
# axes[0].xaxis.set_minor_locator(ticker.MultipleLocator(25))
# lim = axes[0].get_ylim()
# axes[0].set_ylim(0,lim[1])
# # axes[0].set_xlim(0,350.0)
# # spacerHeat = pextWeighted.replace({"spacer_id": newSpacerIDs})
# # # spacerHeat['p_spacer'] = list(map(np.log,spacerHeat.p_spacer))
# # spacerHeat = spacerHeat.pivot_table(index='spacer_id', columns='t', values='p_spacer')
# # axes[1].imshow(spacerHeat,cmap='viridis', aspect='auto')
# # axes[1].set_xlim(0,350.0)
# microbe_total['freq_displaced'] = (CC-np.array(microbe_total['btotal']))/CC
# axes[1].plot(microbe_total['t'],microbe_total['freq_displaced'],linewidth=0.75)
# axes[1].plot(pBurn['t'],pBurn['expBurn'],linewidth=0.75)
# axes[1].plot(pBurn['t'],pBurn['expBurnNet'],linewidth=0.75)






vmatchstrain = pd.read_sql_query("SELECT t, match_id, vstrain_id \
    FROM vmatches WHERE vstrain_id in ({0})"\
    .format(', '.join(map(str,keepStrains))), conTri)\
    .rename(columns={"match_id":"vmatch_id"})
pEmergence = pd.read_sql_query("SELECT t, vmatch_id, p_extinction \
    FROM vmatch_extinction ORDER BY t, vmatch_id", conProb)
pEmergence['p_emergence'] = 1 - np.array(pEmergence['p_extinction'])
pEmergence = pEmergence.drop(columns=['p_extinction'])
pEmergence = pEmergence.merge(vmatchstrain,on=['t','vmatch_id'])

vtreeOrder = pd.read_sql_query("SELECT vstrain_id, tree_vstrain_id\
    FROM tree_vstrain_order WHERE vstrain_id in ({0})"\
    .format(', '.join(map(str,keepStrains))), conTree)\
    .rename(columns={"match_id":"vmatch_id"})

pEmergence = pEmergence.merge(vtreeOrder,on=['vstrain_id'])\
            .sort_values(by=['t','tree_vstrain_id']).reset_index(drop=True)
pEmergence['order'] = pEmergence['tree_vstrain_id'][:]
strainDict = {}
newID = 1
for id in sorted(vtreeOrder['tree_vstrain_id'].unique())[::-1]:
    strainDict[id] = newID
    newID += 10


pEmergence = pEmergence.replace({"order": strainDict})


fig = plt.figure()
ax1 = fig.add_subplot(111, projection='3d')
pEmergeTrunc = pEmergence[(pEmergence.t > 650) & (pEmergence.t < 800)]
pEmergeTrunc =  pEmergeTrunc[(np.array(pEmergeTrunc.t)%4) == 0]
pEmergeTrunc['order'] = pEmergeTrunc['tree_vstrain_id'][:]
strainDict = {}
# newID = 0
# for id in sorted(pEmergeTrunc['tree_vstrain_id'].unique())[::-1]: # THIS IS FOR OLDERST TO THE RIGHT
#     strainDict[id] = newID
#     newID += 1
# pEmergeTrunc = pEmergeTrunc.replace({"order": strainDict}).sort_values(by=['order'])
newID = 0
for id in sorted(pEmergeTrunc['tree_vstrain_id'].unique())[::-1]: # THIS IS FOR OLDEST TO THE LEFT
    strainDict[id] = newID
    newID += 10


pEmergeTrunc = pEmergeTrunc.replace({"order": strainDict}).sort_values(by=['t','order'])
colors = plt.cm.jet(pEmergeTrunc['p_emergence'])
dx = np.ones(len(pEmergeTrunc['t']))*3
dy = np.ones(len(pEmergeTrunc['order']))*3
ax1.bar3d(pEmergeTrunc['t'], pEmergeTrunc['order'], np.zeros(pEmergeTrunc['t'].shape),\
dx, dy, pEmergeTrunc['p_emergence'], color=colors)
ax1.set_yticklabels(pEmergeTrunc['vstrain_id'].unique())
plt.show()



fig = plt.figure(figsize=(10,10))
ax_objs = []
pEmergeTrunc = pEmergence[(pEmergence.t < 600)]
pEmergeTrunc =  pEmergeTrunc[(np.array(pEmergeTrunc.t)%4) == 0]
pEmergeTrunc['order'] = pEmergeTrunc['tree_vstrain_id'][:]
strainDict = {}
# newID = 0
# for id in sorted(pEmergeTrunc['tree_vstrain_id'].unique())[::-1]: # THIS IS FOR OLDERST TO THE RIGHT
#     strainDict[id] = newID
#     newID += 1
# pEmergeTrunc = pEmergeTrunc.replace({"order": strainDict}).sort_values(by=['order'])
newID = 0
for id in sorted(pEmergeTrunc['tree_vstrain_id'].unique()): # THIS IS FOR OLDEST TO THE LEFT
    strainDict[id] = newID
    newID += 1


pEmergeTrunc = pEmergeTrunc.replace({"order": strainDict}).sort_values(by=['t','order'])


times = list(pEmergeTrunc['t'].unique())
gs = (grid_spec.GridSpec(len(times),1))
j = 0
width = 1
textMin = -width/2 -1
textH = 0.05
for t in times:
    # creating new axes object and appending to ax_objs
    if j == 0:
        ax_objs.append(fig.add_subplot(gs[j:j+1, 0:]))
    else:
        ax_objs.append(fig.add_subplot(gs[j:j+1, 0:],sharex=ax_objs[0]))
    # plotting the distribution
    ax_objs[-1].bar(pEmergeTrunc[pEmergeTrunc.t == t]['order'].values,\
    pEmergeTrunc[pEmergeTrunc.t == t]['p_emergence'].values,\
    linewidth=1,edgecolor='black',align='center',width=width,\
    color= plt.cm.jet(pEmergeTrunc[pEmergeTrunc.t == t]['p_emergence']),alpha=0.75)
    # setting uniform x and y lims
    # ax_objs[-1].set_xlim(textMin, max(pEmergeTrunc['order'])+ width/2)
    ax_objs[-1].set_ylim(0,1.2)
    # make background transparent
    rect = ax_objs[-1].patch
    rect.set_alpha(0)
    # remove borders, axis ticks, and labels
    ax_objs[-1].set_yticklabels([])
    ax_objs[-1].set_ylabel('')
    if t == times[-1]:
        ax_objs[-1].set_xlabel("vstrainID", labelpad=10,fontsize=10,fontweight="bold")
    else:
        ax_objs[-1].set_xticklabels([])
        ax_objs[-1].xaxis.set_ticks_position('none')
    ax_objs[-1].yaxis.set_ticks_position('none')
    # spines = ["top","right","left","bottom"]
    spines = ["top","right","left"]
    for s in spines:
        ax_objs[-1].spines[s].set_visible(False)
    ax_objs[-1].text(textMin,textH,'t = {}'.format(t),fontsize=9,ha="center")
    j += 1

gs.update(hspace= -0.7)

# plt.tight_layout()
plt.show()






# 3D BAR CHART SUPPORTED ON MATCH IDs
# vmatchTimes = pd.read_sql_query("SELECT DISTINCT t, match_id FROM vmatches", conTri)
# vmatchTimes = vmatchTimes.rename(columns={"match_id": "vmatch_id"})
#
# matchIDsOrdered = {}
# orderID = 1
# for t in vmatchTimes['t']:
#     for matchID in vmatchTimes[vmatchTimes.t == t]['vmatch_id']:
#         if len(matchIDsOrdered.keys()) == 0:
#             matchIDsOrdered[int(matchID)] = int(orderID)
#             orderID += 1
#         else:
#             if matchID not in matchIDsOrdered.keys():
#                 matchIDsOrdered[int(matchID)] = int(orderID)
#                 orderID += 1
#
# vmatchTimes['order'] = vmatchTimes['vmatch_id'][:]
# vmatchTimes.replace({"order": matchIDsOrdered})
# vmatchTimes = vmatchTimes.sort_values(by=['t', 'order'])
# pEmergence['order'] = pEmergence['vmatch_id'][:]
# pEmergence.replace({"order": matchIDsOrdered})
# pEmergence = pEmergence.sort_values(by=['t', 'order'])
# fig = plt.figure()
# ax1 = fig.add_subplot(111, projection='3d')
# pEmergeTrunc = pEmergence[(pEmergence.t > 0) & (pEmergence.vmatch_id == 647)]
# colors = plt.cm.jet(pEmergeTrunc['p_emergence'])
# dx = np.ones(len(pEmergeTrunc['t']))
# dy = np.ones(len(pEmergeTrunc['vmatch_id']))*10
# ax1.bar3d(pEmergeTrunc['t'], pEmergeTrunc['vmatch_id'], np.zeros(pEmergeTrunc['t'].shape),\
# dx, dy, pEmergeTrunc['p_emergence'], color=colors)
# plt.show()




print('SQLite Query: viral abundance time series data')
virus_stacked = pd.read_sql_query("SELECT t,vstrain_id,abundance FROM vabundance WHERE run_id = {}".format(run_id), conSim)
virus_total = pd.read_sql_query("SELECT t,viral_abundance FROM summary WHERE run_id = {}".format(run_id), conSim)\
.rename(columns={"viral_abundance": "vtotal"})
virus_total = virus_total[virus_total.t <= max(virus_stacked.t)]
# THIS IS FOR THE 3D DISTRIBUTIONS....
vMaxIDs = virus_stacked.set_index('t').groupby(['vstrain_id']).agg(t = ('abundance','idxmax'),\
                            vmaxAbund = ('abundance','max')).reset_index()
vMaxIDs = virus_total.merge(vMaxIDs,on=['t'])
vMaxIDs['vtotal'] = vabundthreshold*np.array(vMaxIDs['vtotal'])
keepStrains = list(vMaxIDs[vMaxIDs['vmaxAbund']>vMaxIDs['vtotal']]['vstrain_id'].values)
virus_stacked = virus_stacked[[(i in keepStrains) for i in virus_stacked.vstrain_id]].reset_index(drop=True)

print('SQLite Query: microbial abundance time series data')
microbe_stacked = pd.read_sql_query("SELECT t,bstrain_id,abundance FROM babundance WHERE run_id = {}".format(run_id), conSim)
microbe_total = pd.read_sql_query("SELECT t,microbial_abundance FROM summary WHERE run_id = {}".format(run_id), conSim)\
.rename(columns={"microbial_abundance": "btotal"})
microbe_stacked = microbe_stacked[microbe_stacked.t <= max(virus_stacked.t)]
microbe_total = microbe_total[microbe_total.t <= max(virus_stacked.t)]

shannon = []
for t in sorted(microbe_total['t'].unique()):
    div = list(microbe_stacked[(microbe_stacked['t']==t) & (microbe_stacked['abundance']!=0)]['abundance'])
    div = np.array(div)/microbe_total[microbe_total['t']==t]['btotal'].values[0]
    shannon.append(np.exp(-1*sum(np.array(div)*np.array(list(map(np.log,div))))))

richness = []
for t in sorted(microbe_total['t'].unique()):
    div = list(microbe_stacked[(microbe_stacked['t']==t) & (microbe_stacked['abundance']!=0)]['abundance'])
    richness.append(len(div))

microbe_stacked = microbe_stacked.pivot(index='t',columns='bstrain_id',values='abundance')

pal = sns.color_palette("tab20c")
fig, ax = plt.subplots(2,sharex=True)
axes = [ax[0], ax[0].twinx(), ax[1], ax[1].twinx()]
microbe_stacked.plot.area(ax = axes[0],stacked=True,legend=False, linewidth=0,color=pal,sort_columns=True)
axes[0].set_ylabel(ylabel ='Microbial Strain Abundances',labelpad=15,fontsize=15)
axes[0].set_xlabel(xlabel = 'Time t',fontsize=15)
axes[0].ticklabel_format(axis = 'y',style='sci',scilimits=(0,0))
axes[0].xaxis.set_minor_locator(ticker.MultipleLocator(25))
axes[0].tick_params(axis='x', labelsize= 10)
axes[0].tick_params(axis='y', labelsize= 10)
lim = axes[0].get_ylim()
axes[0].set_ylim(0,lim[1])
axes[1].plot(virus_total['t'],virus_total['vtotal'],linewidth=0,color='grey')
axes[1].fill_between(virus_total['t'],virus_total['vtotal'], color='grey',alpha=0.6)
lim = axes[1].get_ylim()
axes[1].set_ylim(0,lim[1])
axes[2].plot(microbe_total['t'],shannon,linewidth=1,color='darkblue')
axes[2].set_ylabel(ylabel ='Strain Shannon Diversity',labelpad=15,fontsize=15)
axes[2].set_xlabel(xlabel = 'Time t',fontsize=15)
axes[2].tick_params(axis='x', labelsize= 10)
axes[2].tick_params(axis='y', labelsize= 10)
axes[3].plot(virus_total['t'],virus_total['vtotal'],linewidth=0,color='grey')
axes[3].fill_between(virus_total['t'],virus_total['vtotal'], color='grey',alpha=0.3)
fig.tight_layout()


pal = sns.color_palette("tab20c")
fig, ax = plt.subplots(2,sharex=True)
axes = [ax[0], ax[0].twinx(), ax[1], ax[1].twinx()]
microbe_stacked.plot.area(ax = axes[0],stacked=True,legend=False, linewidth=0,color=pal,sort_columns=True)
axes[0].set_ylabel(ylabel ='Microbial Strain Abundances',labelpad=15,fontsize=15)
axes[0].set_xlabel(xlabel = 'Time t',fontsize=15)
axes[0].ticklabel_format(axis = 'y',style='sci',scilimits=(0,0))
axes[0].xaxis.set_minor_locator(ticker.MultipleLocator(25))
axes[0].tick_params(axis='x', labelsize= 10)
axes[0].tick_params(axis='y', labelsize= 10)
lim = axes[0].get_ylim()
axes[0].set_ylim(0,lim[1])
axes[1].plot(virus_total['t'],virus_total['vtotal'],linewidth=0,color='grey')
axes[1].fill_between(virus_total['t'],virus_total['vtotal'], color='grey',alpha=0.6)
lim = axes[1].get_ylim()
axes[1].set_ylim(0,lim[1])
axes[2].plot(microbe_total['t'],np.array(richness) - np.array(shannon),linewidth=1,color='darkblue')
axes[2].set_ylabel(ylabel ='Strain Diversity Differential',labelpad=15,fontsize=15)
axes[2].set_xlabel(xlabel = 'Time t',fontsize=15)
axes[2].tick_params(axis='x', labelsize= 10)
axes[2].tick_params(axis='y', labelsize= 10)
axes[3].plot(virus_total['t'],virus_total['vtotal'],linewidth=0,color='grey')
axes[3].fill_between(virus_total['t'],virus_total['vtotal'], color='grey',alpha=0.3)
fig.tight_layout()


pal = sns.color_palette("tab20c")
fig, ax = plt.subplots(3,sharex=True)
axes = [ax[0], ax[0].twinx(), ax[1],ax[1].twinx(), ax[2], ax[2].twinx()]
microbe_stacked.plot.area(ax = axes[0],stacked=True,legend=False, linewidth=0,color=pal,sort_columns=True)
axes[0].set_ylabel(ylabel ='Microbial Strain Abundances',labelpad=15,fontsize=15)
axes[0].set_xlabel(xlabel = 'Time t',fontsize=15)
axes[0].ticklabel_format(axis = 'y',style='sci',scilimits=(0,0))
axes[0].xaxis.set_minor_locator(ticker.MultipleLocator(25))
axes[0].tick_params(axis='x', labelsize= 10)
axes[0].tick_params(axis='y', labelsize= 10)
lim = axes[0].get_ylim()
axes[0].set_ylim(0,lim[1])
axes[1].plot(virus_total['t'],virus_total['vtotal'],linewidth=0,color='grey')
axes[1].fill_between(virus_total['t'],virus_total['vtotal'], color='grey',alpha=0.6)
lim = axes[1].get_ylim()
axes[1].set_ylim(0,lim[1])

axes[2].plot(microbe_total['t'],richness,linewidth=1,color='darkblue')
axes[2].set_ylabel(ylabel ='Strain Richness',labelpad=15,fontsize=15)
axes[2].set_xlabel(xlabel = 'Time t',fontsize=15)
axes[2].tick_params(axis='x', labelsize= 10)
axes[2].tick_params(axis='y', labelsize= 10)
axes[3].plot(virus_total['t'],virus_total['vtotal'],linewidth=0,color='grey')
axes[3].fill_between(virus_total['t'],virus_total['vtotal'], color='grey',alpha=0.3)

axes[4].plot(microbe_total['t'],shannon,linewidth=1,color='darkblue')
axes[4].set_ylabel(ylabel ='Strain Shannon Diversity',labelpad=15,fontsize=15)
axes[4].set_xlabel(xlabel = 'Time t',fontsize=15)
axes[4].tick_params(axis='x', labelsize= 10)
axes[4].tick_params(axis='y', labelsize= 10)
axes[5].plot(virus_total['t'],virus_total['vtotal'],linewidth=0,color='grey')
axes[5].fill_between(virus_total['t'],virus_total['vtotal'], color='grey',alpha=0.3)
fig.tight_layout()








print('Complete!')
