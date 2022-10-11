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
from scipy.signal import find_peaks, peak_widths
from scipy.interpolate import interp1d, CubicSpline
import scipy.integrate as integrate
import sqlite3
import matplotlib.ticker as ticker
import colorsys
from mpl_toolkits.axes_grid1 import make_axes_locatable
import math

run_id = sys.argv[1]
resolve = 500
imgTypes = ["pdf"]

SCRIPT_PATH = os.path.abspath(os.path.dirname(__file__)) # cluster
dir = 'crispr-sweep-7-2-2022/isolates/runID3297-c66-r47' # local
run = 'runID3297-c66-r47' # local

DBSIM_PATH = os.path.join(SCRIPT_PATH,'..','..','..','sweep_db_gathered.sqlite') # cluster
# DBSIM_PATH = os.path.join('/Volumes','Yadgah',dir,'{}.sqlite'.format(run)) # local

conSim = sqlite3.connect(DBSIM_PATH)
curSim = conSim.cursor()
ID = curSim.execute('SELECT combo_id,replicate FROM runs WHERE run_id = {}'.format(run_id)).fetchall()
combo_id = ID[0][0]
replicate = ID[0][1]
CC = curSim.execute('SELECT microbe_carrying_capacity FROM param_combos WHERE combo_id = {}'.format(combo_id)).fetchall()
CC = CC[0][0]
RUN_DIR = os.path.join('runID{0}-c{1}-r{2}'.format(run_id,combo_id,replicate))

DBMATCH_PATH = os.path.join(SCRIPT_PATH,'..','..','isolates',,RUN_DIR,'matches_output.sqlite') # cluster
DBTRI_PATH = os.path.join(SCRIPT_PATH,'..','..','isolates',RUN_DIR,'tripartite-networks_output.sqlite') # cluster
DBPROB_PATH = os.path.join(SCRIPT_PATH,'..','..','isolates',RUN_DIR,'probability-of-emergence_output.sqlite') # cluster
# DBMATCH_PATH = os.path.join('/Volumes','Yadgah',dir,'matches_output.sqlite') # local
# DBTRI_PATH = os.path.join('/Volumes','Yadgah',dir,'tripartite-networks_output.sqlite') # local
# DBPROB_PATH = os.path.join('/Volumes','Yadgah',dir,'probability-of-emergence_output.sqlite') # local

conMatch = sqlite3.connect(DBMATCH_PATH)
curMatch = conMatch.cursor()
conTri = sqlite3.connect(DBTRI_PATH)
curTri = conTri.cursor()
conProb = sqlite3.connect(DBPROB_PATH)
curProb = conProb.cursor()

print('SQLite Query: microbial abundance time series data')
microbe_stacked = pd.read_sql_query("SELECT t,bstrain_id,abundance FROM babundance WHERE run_id = {}".format(run_id), conSim)
microbe_stacked = microbe_stacked.pivot(index='t',columns='bstrain_id',values='abundance')
print('SQLite Query: viral abundance time series data')
virus_stacked = pd.read_sql_query("SELECT t,vstrain_id,abundance FROM vabundance WHERE run_id = {}".format(run_id), conSim)
virus_stacked = virus_stacked.pivot(index='t',columns='vstrain_id',values='abundance')
if len(microbe_stacked.index) > len(virus_stacked.index):
    microbe_stacked.drop(index=list(set(microbe_stacked.index)-set(virus_stacked.index)),inplace=True)

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
        print('oops')
        break
    spacersOrdered.extend(newSpacer)
    print(phenotype)

newSpacerIDs = dict(zip(spacersOrdered,range(1,len(spacersOrdered)+1)))

print('SQLite Query: tripartite network data')
tripartiteNets = pd.read_sql_query("SELECT t, vmatch_id, spacer_id, bmatch_id \
    FROM single_match_tripartite_networks", conTri)
vmatchspacers = pd.read_sql_query("SELECT t, vmatch_id, spacer_id \
    FROM single_match_tripartite_networks", conTri)
numSpacers = vmatchspacers.groupby(['t','vmatch_id']).agg(
     numSpacers = ('spacer_id','size')
     ).reset_index()

# gather p_extinction
pextinction = pd.DataFrame()
for t in np.unique(tripartiteNets.t):
    vmatchIDs = np.unique(tripartiteNets[tripartiteNets.t == t].vmatch_id)
    pext = pd.read_sql_query("SELECT t, vmatch_id, p_extinction \
        FROM vmatch_extinction WHERE t = {0} \
        AND vmatch_id in ({1}) ORDER BY t, vmatch_id"
        .format(t,', '.join(map(str,vmatchIDs))), conProb)
    pextinction = pd.concat([pextinction, pext], ignore_index=True, sort=False)

# ## p_emergence from p_ext & normalize
# pextWeighted = pextinction[:]
# pextWeighted['p_extinction_perSpacer'] = 1 - pextWeighted['p_extinction']
# pextWeighted = pextWeighted.drop(columns=['p_extinction'])
# pnormal = pextWeighted.groupby(['t']) \
#                 .agg(p_normal=('p_extinction_perSpacer','sum')).reset_index()
# pextWeighted['p_extinction_perSpacer'] = pextWeighted['p_extinction_perSpacer']/pextWeighted.merge(pnormal, on =['t'])['p_normal']
# pextWeighted = pextWeighted.fillna(0)

# p_emergence from p_ext, weight by num of spacers, & normalize
vMatchProbs = pd.DataFrame()
vMatchProbs = pd.concat([vMatchProbs, pextinction])
vMatchProbs['p_emergence'] = list(1 - np.array(vMatchProbs['p_extinction']))
vMatchProbs = pd.merge(vMatchProbs, numSpacers, on=["t","vmatch_id"])
vMatchProbs['p_emergence_spacer'] = \
    np.array(vMatchProbs['p_emergence'])/np.array(vMatchProbs['numSpacers'])
vMatchProbs = vMatchProbs.drop(columns = ['numSpacers'])
pnormal = vMatchProbs.groupby(['t']) \
                .agg(p_normal=('p_emergence_spacer','sum')).reset_index()
vMatchProbs['p_emergence_spacer_normalized'] = vMatchProbs['p_emergence_spacer']/vMatchProbs.merge(pnormal, on =['t'])['p_normal']
vMatchProbs = vMatchProbs.fillna(0)


bmatchAbundances = pd.read_sql_query("SELECT t, match_id, babundance \
    FROM bmatches_abundances WHERE match_id in ({0})"
    .format(', '.join(map(str,np.unique(tripartiteNets.bmatch_id)))), conTri)\
    .rename(columns={"match_id": "bmatch_id"})\
    .merge(tripartiteNets,on=["t","bmatch_id"]).drop(columns=['spacer_id'])\
    .drop_duplicates()

# bmatchAbundances = bmatchAbundances.drop(columns=['vmatch_id']).drop_duplicates()\
#                     .groupby(['t']).agg(bmatchTotal=('babundance','sum')).reset_index()\
#                     .merge(bmatchAbundances,on=['t'])
#
# bmatchAbundances['burn_freq'] = bmatchAbundances['bmatchTotal'] - bmatchAbundances['babundance']
# bmatchAbundances['burn_freq'] = (CC-np.array(bmatchAbundances['burn_freq']))/CC

bmatchAbundances =  pd.read_sql_query("SELECT t,microbial_abundance FROM summary \
WHERE run_id = {0}".format(run_id),conSim).rename(columns={"microbial_abundance": "burn_freq"})\
.merge(bmatchAbundances, on=['t'])

# bmatchAbundances['burn_freq'] = (CC-np.array(bmatchAbundances['babundance']))/CC

bmatchAbundances['burn_freq'] = (CC-np.array(bmatchAbundances['burn_freq'])\
                                + np.array(bmatchAbundances['babundance']))/CC


burnFreq = vMatchProbs.drop(columns=['p_extinction','p_emergence','p_emergence_spacer'])\
            .merge(bmatchAbundances,on=['t','vmatch_id']).drop(columns=['babundance','vmatch_id'])\
            .groupby(['t','bmatch_id','burn_freq']).agg(p_emerge=('p_emergence_spacer_normalized','sum')).reset_index()

pburnNormal = burnFreq.groupby(['t']).agg(pburnNormal=('p_emerge','sum')).reset_index()\
                .merge(burnFreq,on=['t']).drop(columns=['bmatch_id','burn_freq','p_emerge'])
burnFreq['p_emerge'] = [i/j if j>0 else 0 for (i,j) in \
zip(np.array(burnFreq['p_emerge']),np.array(pburnNormal['pburnNormal']))]
burnFreq['burn_freq'] = np.array(burnFreq['burn_freq'])*np.array(burnFreq['p_emerge'])
burnFreq = burnFreq.drop(columns=['p_emerge','bmatch_id']).groupby(['t'])\
            .agg(exp_burn=('burn_freq','sum')).reset_index()

freqDisplaced =  pd.read_sql_query("SELECT t,microbial_abundance FROM summary \
WHERE run_id = {0}".format(run_id),conSim).rename(columns={"microbial_abundance": "bTotal"})\

freqDisplaced['freq_displaced'] = (CC-np.array(freqDisplaced['bTotal']))/CC

freqDisplaced = freqDisplaced[freqDisplaced['t'] <= max(virus_stacked.index)]

# grid_kws = {"height_ratios": (1, 1, 1), "hspace": (.1, .1)}
pal = sns.color_palette("tab20b")
# fig, ax = plt.subplots(2,sharex=True)
# axes = [ax[0], ax[1]]
fig, ax = plt.subplots(2,sharex=True)
axes = [ax[0], ax[1]]

microbe_stacked.plot.area(ax = axes[0],stacked=True,legend=False, linewidth=0,color=pal,sort_columns=True)
axes[0].set_ylabel(ylabel ='Microbial Strain Abundances',labelpad=15,fontsize=7)
axes[0].set_xlabel(xlabel = 'Time t',fontsize=7)
axes[0].ticklabel_format(axis = 'y',style='sci',scilimits=(0,0))
axes[0].xaxis.set_minor_locator(ticker.MultipleLocator(25))
lim = axes[0].get_ylim()
axes[0].set_ylim(0,lim[1])
# axes[0].set_xlim(0,350.0)
# spacerHeat = pextWeighted.replace({"spacer_id": newSpacerIDs})
# # spacerHeat['p_spacer'] = list(map(np.log,spacerHeat.p_spacer))
# spacerHeat = spacerHeat.pivot_table(index='spacer_id', columns='t', values='p_spacer')
# axes[1].imshow(spacerHeat,cmap='viridis', aspect='auto')
# axes[1].set_xlim(0,350.0)
burnFreq.plot(x='t',ax = axes[1],legend=False,linewidth=0.75)
axes[1].plot(freqDisplaced['t'],freqDisplaced['freq_displaced'],linewidth=0.75)













bmatchAbundances = pd.read_sql_query("SELECT t, match_id, babundance \
    FROM bmatches_abundances WHERE match_id in ({0})"
    .format(', '.join(map(str,np.unique(tripartiteNets.bmatch_id)))), conTri)
bmatchAbundances = bmatchAbundances.rename(columns={"match_id": "bmatch_id"})
bmatchAbundances = pd.merge(tripartiteNets, bmatchAbundances, on=["t","bmatch_id"])
bmatchAbundances = bmatchAbundances.drop(['vmatch_id'],1)
bmatchAbundances = bmatchAbundances.drop_duplicates()

# totalB = bmatchAbundances.drop(['spacer_id'],1)
# totalB = totalB.drop_duplicates()
totalB = bmatchAbundances.groupby(['t']).agg(
     bTotal = ('babundance','sum')
     ).reset_index()
bmatchAbundances = pd.merge(bmatchAbundances, totalB, on=["t"])
bmatchAbundances['bfrequency'] = np.array(bmatchAbundances['babundance'])*1./np.array(bmatchAbundances['bTotal'])






pBurnSpacers = bmatchAbundances.groupby(['t','spacer_id']).agg(
     p_burn = ('bfrequency','sum')
     ).reset_index()

pBurnSpacers = pd.merge(pBurnSpacers, pExtinctionSpacers, on=["t","spacer_id"])
pBurnSpacers['p_burn'] = np.array(pBurnSpacers.p_burn)*np.array(pBurnSpacers.p_spacer)
pBurnSpacers = pBurnSpacers.drop(['p_spacer'],1)



# this is p*B/Btotal
pExtinctionSpacers = pd.merge(pextinction, tripartiteNets, on=["t","vmatch_id"])
pExtinctionSpacers = pExtinctionSpacers.drop(['bmatch_id'],1)
pExtinctionSpacers = pExtinctionSpacers.drop_duplicates()
pExtinctionSpacers = pExtinctionSpacers.groupby(['t','spacer_id']).agg(
     p_spacer = ('p_extinction_perSpacer','sum')
     ).reset_index()
pBurnSpacers = bmatchAbundances.groupby(['t','spacer_id']).agg(
     p_burn = ('bfrequency','sum')
     ).reset_index()
pBurnSpacers['p_burn'] = np.array(pBurnSpacers.p_burn)*np.array(pExtinctionSpacers.p_spacer)
# grid_kws = {"height_ratios": (1, 1, 1), "hspace": (.1, .1)}
pal = sns.color_palette("tab20b")
fig, ax = plt.subplots(3,sharex=True)
axes = [ax[0], ax[1], ax[2]]
microbe_stacked.plot.area(ax = axes[0],stacked=True,legend=False, linewidth=0,color=pal,sort_columns=True)
axes[0].set_ylabel(ylabel ='Microbial Strain Abundances',labelpad=15,fontsize=7)
axes[0].set_xlabel(xlabel = 'Time t',fontsize=7)
axes[0].ticklabel_format(axis = 'y',style='sci',scilimits=(0,0))
axes[0].xaxis.set_minor_locator(ticker.MultipleLocator(25))
lim = axes[0].get_ylim()
axes[0].set_ylim(0,lim[1])
# axes[0].set_xlim(0,350.0)
spacerHeat = pBurnSpacers.replace({"spacer_id": newSpacerIDs})
# spacerHeat['p_burn'] = list(map(np.log,spacerHeat.p_burn))
spacerHeat = spacerHeat.pivot_table(index='spacer_id', columns='t', values='p_burn')
axes[1].imshow(spacerHeat,cmap='RdPu', aspect='auto')
# axes[1].set_xlim(0,350.0)
pBurnAvg = pBurnSpacers.groupby(['t']).agg(
     p_burn_exp = ('p_burn','sum')
     ).reset_index()
pBurnAvg.plot(x='t',ax = axes[2],legend=False,color=pal[0],linewidth=0.75)



vmatchAbundances = pd.read_sql_query("SELECT t, match_id, vabundance \
    FROM vmatches_abundances WHERE match_id in ({0})"
    .format(', '.join(map(str,np.unique(tripartiteNets.vmatch_id)))), conTri)
vmatchAbundances = vmatchAbundances.rename(columns={"match_id": "vmatch_id"})
vmatchAbundances= pd.merge(vmatchAbundances, tripartiteNets, on=["t","vmatch_id"])
vmatchAbundances = vmatchAbundances.drop(['bmatch_id'],1)
totalV = vmatchAbundances.groupby(['t']).agg(
     vTotal = ('vabundance','sum')
     ).reset_index()
pBurnSpacers2 = vmatchAbundances.drop_duplicates().groupby(['t','spacer_id']).agg(
     p_burn = ('vabundance','sum')
     ).reset_index()
pBurnSpacers2 = pd.merge(pBurnSpacers2, totalV, on=["t"])
pBurnSpacers2['p_burn'] = np.array(pBurnSpacers2['p_burn'])/np.array(pBurnSpacers2['vTotal'])
pBurnSpacers2 = pBurnSpacers2.drop(['vTotal'],1)
pBurnSpacers2['p_burn'] = np.array(pBurnSpacers2['p_burn'])*np.array(pBurnSpacers['p_burn'])
grid_kws = {"height_ratios": (1, 1), "hspace": (.1)}
pal = sns.color_palette("tab20b")
fig, ax = plt.subplots(2,sharex=True, gridspec_kw=grid_kws)
axes = [ax[0], ax[1]]
microbe_stacked.plot.area(ax = axes[0],stacked=True,legend=False, linewidth=0,color=pal,sort_columns=True)
axes[0].set_ylabel(ylabel ='Microbial Strain Abundances',labelpad=15,fontsize=7)
axes[0].set_xlabel(xlabel = 'Time t',fontsize=7)
axes[0].ticklabel_format(axis = 'y',style='sci',scilimits=(0,0))
axes[0].xaxis.set_minor_locator(ticker.MultipleLocator(25))
lim = axes[0].get_ylim()
axes[0].set_ylim(0,lim[1])
# axes[0].set_xlim(0,350.0)
spacerHeat = pBurnSpacers2.replace({"spacer_id": newSpacerIDs})
# spacerHeat['p_burn'] = list(map(np.log,spacerHeat.p_burn))
spacerHeat = spacerHeat.pivot_table(index='spacer_id', columns='t', values='p_burn')
axes[1].imshow(spacerHeat,cmap='viridis', aspect='auto')






grid_kws = {"height_ratios": (1, 1), "hspace": (.1)}
pal = sns.color_palette("tab20b")
fig, ax = plt.subplots(2,sharex=True, gridspec_kw=grid_kws)
axes = [ax[0], ax[1]]
microbe_stacked.plot.area(ax = axes[0],stacked=True,legend=False, linewidth=0,color=pal,sort_columns=True)
axes[0].set_ylabel(ylabel ='Microbial Strain Abundances',labelpad=15,fontsize=7)
axes[0].set_xlabel(xlabel = 'Time t',fontsize=7)
axes[0].ticklabel_format(axis = 'y',style='sci',scilimits=(0,0))
axes[0].xaxis.set_minor_locator(ticker.MultipleLocator(25))
lim = axes[0].get_ylim()
axes[0].set_ylim(0,lim[1])
# axes[0].set_xlim(0,350.0)
spacerHeat = pExtinctionSpacers.replace({"spacer_id": newSpacerIDs})
spacerHeat = spacerHeat.pivot_table(index='spacer_id', columns='t', values='p_spacer')
axes[1].imshow(spacerHeat,cmap='viridis', aspect='auto')




print('Complete!')
