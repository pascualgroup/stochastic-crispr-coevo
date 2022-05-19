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
from matplotlib.collections import LineCollection

run_id = sys.argv[1]

if len(sys.argv[2]) == 0:
    lifeTimeThreshold = 0
else:
    lifeTimeThreshold = float(sys.argv[2])

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

DBCLADE_PATH = os.path.join(SCRIPT_PATH,'..','..','isolates',RUN_DIR,'clade-abundances_output.sqlite') # cluster
# DBCLADE_PATH = os.path.join('/Volumes','Yadgah','clade-abundances_output.sqlite') # local. run_id fixed; for testing
# DBCLADE_PATH = os.path.join('/Volumes','Yadgah','crispr-sweep-7-2-2022/isolates/runID3297-c66-r47','clade-abundances_output.sqlite') # local. run_id fixed; for testing
DBTREE_PATH = os.path.join(SCRIPT_PATH,'..','..','isolates',RUN_DIR,'trees_output.sqlite') # cluster
# DBTREE_PATH = os.path.join('/Volumes','Yadgah','trees_output.sqlite') # local. run_id fixed; for testing
# DBTREE_PATH = os.path.join('/Volumes','Yadgah','crispr-sweep-7-2-2022/isolates/runID3297-c66-r47','trees_output.sqlite') # local. run_id fixed; for testing

print('SQLite Query: microbial abundance time series data')
microbeSim = pd.read_sql_query("SELECT t,bstrain_id,abundance FROM babundance WHERE run_id = {}".format(run_id), conSim)
microbe_stacked = microbeSim.pivot(index='t',columns='bstrain_id',values='abundance')
print('SQLite Query: viral abundance time series data')
virusSim = pd.read_sql_query("SELECT t,vstrain_id,abundance FROM vabundance WHERE run_id = {}".format(run_id), conSim)
virus_stacked = virusSim.pivot(index='t',columns='vstrain_id',values='abundance')

if (microbeSim['t'][microbeSim['t'].size-1] not in virus_stacked.index):
    microbe_stacked.drop(microbeSim['t'][microbeSim['t'].size-1],inplace=True)

if (virusSim['t'][virusSim['t'].size-1] not in microbe_stacked.index):
    virus_stacked.drop(virusSim['t'][virusSim['t'].size-1],inplace=True)


conClade = sqlite3.connect(DBCLADE_PATH)
curClade = conClade.cursor()
conTree = sqlite3.connect(DBTREE_PATH)
curTree = conTree.cursor()
print('SQLite Query: clade data')
microbeClades = pd.read_sql_query("SELECT DISTINCT clade_id, bstrain_id \
FROM babundances", conClade)
microbeCladeIDs = pd.read_sql_query("SELECT DISTINCT clade_id \
FROM babundances", conClade)
microbeCladeAbundances = pd.read_sql_query("SELECT t, clade_id, abundance \
FROM clade_babundances", conClade)
virusClades = pd.read_sql_query("SELECT DISTINCT clade_id, vstrain_id \
FROM vabundances", conClade)
virusCladeIDs = pd.read_sql_query("SELECT DISTINCT clade_id \
FROM vabundances", conClade)
virusCladeAbundances = pd.read_sql_query("SELECT t, clade_id, abundance \
FROM clade_vabundances", conClade)

print('SQLite Query: tree data')
bStrainTimes = pd.read_sql_query(
"SELECT tree_bstrain_id, t_creation, t_extinction, tree_parent_bstrain_id \
FROM tree_bstrain_creation_extinction", conTree)
bTreeAbundances = pd.read_sql_query(
"SELECT t, tree_bstrain_id, abundance \
FROM tree_babundance", conTree)
vStrainTimes = pd.read_sql_query(
"SELECT tree_vstrain_id, t_creation, t_extinction, tree_parent_vstrain_id \
FROM tree_vstrain_creation_extinction", conTree)
vTreeAbundances = pd.read_sql_query(
"SELECT t, tree_vstrain_id, abundance \
FROM tree_vabundance", conTree)

if bTreeAbundances['t'][bTreeAbundances['t'].size-1] not in vTreeAbundances.t.values:
    bTreeAbundances = bTreeAbundances[bTreeAbundances.t != bTreeAbundances['t'][bTreeAbundances['t'].size-1]]

if vTreeAbundances['t'][vTreeAbundances['t'].size-1] not in bTreeAbundances.t.values:
    vTreeAbundances = vTreeAbundances[vTreeAbundances.t != vTreeAbundances['t'][vTreeAbundances['t'].size-1]]


print('Compiling microbial tree plots')
# Get a color map
Mcmap = cm.get_cmap('turbo')
# Get normalize function (takes data in range [vmin, vmax] -> [0, 1])
Mnorm = Normalize(vmin=1, vmax=len(microbeCladeIDs))
microbeColorDict = {}
cladeColorDict = {}
cladeDict = {}
colorInd = 0
for clade_id in microbeCladeIDs.values:
    cladeColorDict[clade_id[0]] = colorInd
    cladeDict[clade_id[0]] = []
    colorInd += 1
for strain in microbeClades.bstrain_id.values:
    clade = microbeClades[microbeClades['bstrain_id']==strain]['clade_id'].values[0]
    cladeDict[clade] = np.append(cladeDict[clade],strain)
columnOrder = []
for clade_id in microbeCladeIDs.values[::-1]:
    columnOrder = np.append(columnOrder,cladeDict[clade_id[0]])
columnOrder = columnOrder.astype(int)
microbe_stacked = microbe_stacked[columnOrder]
for strain in microbe_stacked.columns.values:
    clade = microbeClades[microbeClades['bstrain_id']==strain]['clade_id'].values[0]
    microbeColorDict[strain] = Mcmap(Mnorm(np.arange(1, len(microbeCladeIDs)+1, 1)))[cladeColorDict[clade]]

fig1, ax1 = plt.subplots(1,sharex=True,figsize=(20,20))
ax1.set_xlim(0,np.min([np.max(vTreeAbundances.t.values),np.max(bTreeAbundances.t.values)]))
ax1.set_ylim(0,np.max(bStrainTimes.tree_bstrain_id)+1)
ax1.set_yticks([])
hlinecMicrobe = []
vlinecMicrobe = []
hcolorsMicrobe = []
vcolorsMicrobe = []
maxTickSizeMicrobe = 500
maxAbundanceMicrobe = np.max(bTreeAbundances.abundance.values)
markerIncMicrobe = maxTickSizeMicrobe/maxAbundanceMicrobe
markerColorsMicrobe = []
for clade_id in np.unique(microbeClades.clade_id):
    strainsOfClade = curClade.execute("SELECT DISTINCT bstrain_id FROM babundances \
    WHERE clade_id = {}".format(clade_id))
    strainsOfClade = [i[0] for i in strainsOfClade.fetchall()]
    treeStrainsOfClade = curTree.execute("SELECT DISTINCT tree_bstrain_id FROM tree_bstrain_order \
    WHERE bstrain_id in ({}) ORDER BY tree_bstrain_id".format(', '.join(map(str,strainsOfClade))))
    treeStrainsOfClade = [i[0] for i in treeStrainsOfClade.fetchall()]
    ax1.scatter(bTreeAbundances[bTreeAbundances['tree_bstrain_id'].isin(treeStrainsOfClade)]['t'],
    bTreeAbundances[bTreeAbundances['tree_bstrain_id'].isin(treeStrainsOfClade)]['tree_bstrain_id'],
    lw=.5,
    s=bTreeAbundances[bTreeAbundances['tree_bstrain_id'].isin(treeStrainsOfClade)]['abundance'].values*markerIncMicrobe,
    c = np.array([Mcmap(Mnorm(np.arange(1, len(microbeCladeIDs)+1, 1)))[cladeColorDict[clade_id]]]),marker='|')
    for strain in treeStrainsOfClade:
        # print(strain)
        tCreate = bStrainTimes[bStrainTimes.tree_bstrain_id == strain].t_creation.values[0]
        tExtinct = bStrainTimes[bStrainTimes.tree_bstrain_id == strain].t_extinction.values[0]
        parent = bStrainTimes[bStrainTimes.tree_bstrain_id == strain].tree_parent_bstrain_id.values[0]
        hlinecMicrobe.append([[tCreate, strain],[tExtinct, strain]])
        vlinecMicrobe.append([[tCreate, parent],[tCreate, strain]])
        hcolorsMicrobe.append(Mcmap(Mnorm(np.arange(1, len(microbeCladeIDs)+1, 1)))[cladeColorDict[clade_id]])
        vcolorsMicrobe.append(Mcmap(Mnorm(np.arange(1, len(microbeCladeIDs)+1, 1)))[cladeColorDict[clade_id]])
        markerColorsMicrobe.append(Mcmap(Mnorm(np.arange(1, len(microbeCladeIDs)+1, 1)))[cladeColorDict[clade_id]])




strainLineages = LineCollection(hlinecMicrobe, linestyles='solid', colors=hcolorsMicrobe,linewidths=(0.35))
creationLines = LineCollection(vlinecMicrobe, linestyles='solid', colors=vcolorsMicrobe, linewidths=(0.35))
ax1.add_collection(strainLineages)
ax1.add_collection(creationLines)
fig1.savefig(os.path.join(SCRIPT_PATH,'..','..','isolates',RUN_DIR,'microbe-tree.{}'.format(imgType)),dpi=resolve)
plt.close(fig1)
# fig1.savefig(os.path.join('/Volumes/Yadgah','microbe-tree.{}'.format(imgType)),dpi=resolve)
# plt.show()



print('Compiling stacked microbe clade abundances and tree plots')
fig, ax = plt.subplots(2,sharex=True, figsize=(20,20), gridspec_kw={'height_ratios': [1, 3]})
fig.suptitle('Strain Abundances (run{0}-c{1}-r{2})'.format(run_id,combo_id,replicate))
axes = [ax[0], ax[1]]
microbe_stacked.plot.area(ax = axes[0],stacked=True,legend=False, linewidth=0,color=microbeColorDict,sort_columns=True)
microbe_stacked.plot(stacked=True, ax=axes[0], legend=False, color='black',sort_columns=True,linewidth=.1)
axes[0].set_ylabel(ylabel ='Microbial Strain Abundances',labelpad=15,fontsize=7)
axes[0].set_xlabel(xlabel = 'Time t',fontsize=7)
axes[0].ticklabel_format(axis = 'y',style='sci',scilimits=(0,0))
axes[0].xaxis.set_minor_locator(ticker.MultipleLocator(25))
lim = axes[0].get_ylim()
axes[0].set_ylim(0,lim[1])
axes[0].set_xlim(0,np.min([np.max(vTreeAbundances.t.values),np.max(bTreeAbundances.t.values)]))
axes[1].set_xlabel(xlabel = 'Time t',fontsize=7)
axes[1].set_ylim(0,np.max(bStrainTimes.tree_bstrain_id)+1)
axes[1].set_xlim(0,np.min([np.max(vTreeAbundances.t.values),np.max(bTreeAbundances.t.values)]))
axes[1].xaxis.set_minor_locator(ticker.MultipleLocator(25))
axes[1].set_yticks([])
hlinecMicrobe = []
vlinecMicrobe = []
hcolorsMicrobe = []
vcolorsMicrobe = []
maxTickSizeMicrobe = 500
maxAbundanceMicrobe = np.max(bTreeAbundances.abundance.values)
markerIncMicrobe = maxTickSizeMicrobe/maxAbundanceMicrobe
markerColorsMicrobe = []
for clade_id in np.unique(microbeClades.clade_id):
    strainsOfClade = curClade.execute("SELECT DISTINCT bstrain_id FROM babundances \
    WHERE clade_id = {}".format(clade_id))
    strainsOfClade = [i[0] for i in strainsOfClade.fetchall()]
    treeStrainsOfClade = curTree.execute("SELECT DISTINCT tree_bstrain_id FROM tree_bstrain_order \
    WHERE bstrain_id in ({}) ORDER BY tree_bstrain_id".format(', '.join(map(str,strainsOfClade))))
    treeStrainsOfClade = [i[0] for i in treeStrainsOfClade.fetchall()]
    axes[1].scatter(bTreeAbundances[bTreeAbundances['tree_bstrain_id'].isin(treeStrainsOfClade)]['t'],
    bTreeAbundances[bTreeAbundances['tree_bstrain_id'].isin(treeStrainsOfClade)]['tree_bstrain_id'],
    lw=.5,
    s=bTreeAbundances[bTreeAbundances['tree_bstrain_id'].isin(treeStrainsOfClade)]['abundance'].values*markerIncMicrobe,
    c = np.array([Mcmap(Mnorm(np.arange(1, len(microbeCladeIDs)+1, 1)))[cladeColorDict[clade_id]]]),marker='|')
    for strain in treeStrainsOfClade:
        tCreate = bStrainTimes[bStrainTimes.tree_bstrain_id == strain].t_creation.values[0]
        tExtinct = bStrainTimes[bStrainTimes.tree_bstrain_id == strain].t_extinction.values[0]
        parent = bStrainTimes[bStrainTimes.tree_bstrain_id == strain].tree_parent_bstrain_id.values[0]
        hlinecMicrobe.append([[tCreate, strain],[tExtinct, strain]])
        vlinecMicrobe.append([[tCreate, parent],[tCreate, strain]])
        hcolorsMicrobe.append(Mcmap(Mnorm(np.arange(1, len(microbeCladeIDs)+1, 1)))[cladeColorDict[clade_id]])
        vcolorsMicrobe.append(Mcmap(Mnorm(np.arange(1, len(microbeCladeIDs)+1, 1)))[cladeColorDict[clade_id]])
        markerColorsMicrobe.append(Mcmap(Mnorm(np.arange(1, len(microbeCladeIDs)+1, 1)))[cladeColorDict[clade_id]])

strainLineages = LineCollection(hlinecMicrobe, linestyles='solid', colors=hcolorsMicrobe,linewidths=(0.35))
creationLines = LineCollection(vlinecMicrobe, linestyles='solid', colors=vcolorsMicrobe, linewidths=(0.35))
axes[1].add_collection(strainLineages)
axes[1].add_collection(creationLines)

fig.tight_layout()
fig.savefig(os.path.join(SCRIPT_PATH,'..','..','isolates',RUN_DIR,'microbe-clades-abundances-and-trees.{}'.format(imgType)),dpi=resolve)
# fig.savefig(os.path.join('/Volumes/Yadgah','microbe-clades-abundances-and-trees.{}'.format(imgType)),dpi=resolve)
# plt.show()
plt.close(fig)


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
fig.savefig(os.path.join(SCRIPT_PATH,'..','..','isolates',RUN_DIR,'microbe-strain-abundances.{}'.format(imgType)),dpi=resolve)
plt.close(fig)



Vcmap = cm.get_cmap('turbo')
# Get normalize function (takes data in range [vmin, vmax] -> [0, 1])
Vnorm = Normalize(vmin=1, vmax=len(virusCladeIDs))
virusColorDict = {}
cladeColorDict = {}
cladeDict = {}
colorInd = 0
for clade_id in virusCladeIDs.values:
    # print(id[0])
    cladeColorDict[clade_id[0]] = colorInd
    cladeDict[clade_id[0]] = []
    colorInd += 1

for strain in virusClades.vstrain_id.values:
    clade = virusClades[virusClades['vstrain_id']==strain]['clade_id'].values[0]
    cladeDict[clade] = np.append(cladeDict[clade],strain)

columnOrder = []
for id in virusCladeIDs.values[::-1]:
    columnOrder = np.append(columnOrder,cladeDict[id[0]])

columnOrder = columnOrder.astype(int)
virus_stacked = virus_stacked[columnOrder]

for strain in virus_stacked.columns.values:
    clade = virusClades[virusClades['vstrain_id']==strain]['clade_id'].values[0]
    # print('strain: {}, clade: {}'.format(strain,clade))
    virusColorDict[strain] = Vcmap(Vnorm(np.arange(1, len(virusCladeIDs)+1, 1)))[cladeColorDict[clade]]

fig2, ax2 = plt.subplots(1,sharex=True,figsize=(20,20))
ax2.set_xlim(0,np.min([np.max(vTreeAbundances.t.values),np.max(bTreeAbundances.t.values)]))
ax2.set_ylim(0,np.max(vStrainTimes.tree_vstrain_id)+1)
hlinecVirus = []
vlinecVirus = []
hcolorsVirus = []
vcolorsVirus = []
maxTickSizeVirus = 500
maxAbundanceVirus = np.max(vTreeAbundances.abundance.values)
markerIncVirus = maxTickSizeVirus/maxAbundanceVirus
markerColorsVirus = []
for clade_id in np.unique(virusClades.clade_id):
    strainsOfClade = curClade.execute("SELECT DISTINCT vstrain_id FROM vabundances \
    WHERE clade_id = {}".format(clade_id))
    strainsOfClade = [i[0] for i in strainsOfClade.fetchall()]
    treeStrainsOfCladeAll = curTree.execute("SELECT DISTINCT tree_vstrain_id FROM tree_vstrain_order \
    WHERE vstrain_id in ({}) ORDER BY tree_vstrain_id".format(', '.join(map(str,strainsOfClade))))
    treeStrainsOfClade = []
    for (i,) in treeStrainsOfCladeAll:
        tCreate = vStrainTimes[vStrainTimes.tree_vstrain_id == i].t_creation.values[0]
        tExtinct = vStrainTimes[vStrainTimes.tree_vstrain_id == i].t_extinction.values[0]
        if tExtinct - tCreate > lifeTimeThreshold:
            treeStrainsOfClade.append(i)
    ax2.scatter(vTreeAbundances[vTreeAbundances['tree_vstrain_id'].isin(treeStrainsOfClade)]['t'],
    vTreeAbundances[vTreeAbundances['tree_vstrain_id'].isin(treeStrainsOfClade)]['tree_vstrain_id'],
    lw=.5,
    s=vTreeAbundances[vTreeAbundances['tree_vstrain_id'].isin(treeStrainsOfClade)]['abundance'].values*markerIncVirus,
    c = np.array([Vcmap(Vnorm(np.arange(1, len(virusCladeIDs)+1, 1)))[cladeColorDict[clade_id]]]),marker='|')
    for strain in treeStrainsOfClade:
        tCreate = vStrainTimes[vStrainTimes.tree_vstrain_id == strain].t_creation.values[0]
        tExtinct = vStrainTimes[vStrainTimes.tree_vstrain_id == strain].t_extinction.values[0]
        parent = vStrainTimes[vStrainTimes.tree_vstrain_id == strain].tree_parent_vstrain_id.values[0]
        hlinecVirus.append([[tCreate, strain],[tExtinct, strain]])
        vlinecVirus.append([[tCreate, parent],[tCreate, strain]])
        hcolorsVirus.append(Vcmap(Vnorm(np.arange(1, len(virusCladeIDs)+1, 1)))[cladeColorDict[clade_id]])
        vcolorsVirus.append(Vcmap(Vnorm(np.arange(1, len(virusCladeIDs)+1, 1)))[cladeColorDict[clade_id]])
        markerColorsVirus.append(Vcmap(Vnorm(np.arange(1, len(virusCladeIDs)+1, 1)))[cladeColorDict[clade_id]])


strainLineages = LineCollection(hlinecVirus, linestyles='solid', colors=hcolorsVirus,linewidths=(0.35))
creationLines = LineCollection(vlinecVirus, linestyles='solid', colors=vcolorsVirus, linewidths=(0.35))
ax2.add_collection(strainLineages)
ax2.add_collection(creationLines)
fig2.savefig(os.path.join(SCRIPT_PATH,'..','..','isolates',RUN_DIR,'virus-tree.{}'.format(imgType)),dpi=resolve)
# fig2.savefig(os.path.join('/Volumes/Yadgah','virus-tree.{}'.format(imgType)),dpi=resolve)
# plt.show()
plt.close(fig2)


print('Compiling stacked virus clade abundances and tree plots')
fig, ax = plt.subplots(2,sharex=True, figsize=(20,20), gridspec_kw={'height_ratios': [1, 3]})
fig.suptitle('Strain Abundances (run{0}-c{1}-r{2})'.format(run_id,combo_id,replicate))
axes = [ax[0], ax[1]]
virus_stacked.plot.area(ax = axes[0],stacked=True,legend=False, linewidth=0,color=virusColorDict,sort_columns=True)
virus_stacked.plot(stacked=True, ax=axes[0], legend=False, color='black',sort_columns=True,linewidth=.1)
axes[0].set_ylabel(ylabel ='Viral Strain Abundances',labelpad=15,fontsize=7)
axes[0].set_xlabel(xlabel = 'Time t',fontsize=7)
axes[0].ticklabel_format(axis = 'y',style='sci',scilimits=(0,0))
axes[0].xaxis.set_minor_locator(ticker.MultipleLocator(25))
lim = axes[0].get_ylim()
axes[0].set_ylim(0,lim[1])
axes[0].set_xlim(0,np.min([np.max(vTreeAbundances.t.values),np.max(bTreeAbundances.t.values)]))
axes[1].set_xlabel(xlabel = 'Time t',fontsize=7)
axes[1].set_ylim(0,np.max(vStrainTimes.tree_vstrain_id)+1)
axes[1].set_xlim(0,np.min([np.max(vTreeAbundances.t.values),np.max(bTreeAbundances.t.values)]))
axes[1].xaxis.set_minor_locator(ticker.MultipleLocator(25))
axes[1].set_yticks([])
hlinecVirus = []
vlinecVirus = []
hcolorsVirus = []
vcolorsVirus = []
maxTickSizeVirus = 500
maxAbundanceVirus = np.max(vTreeAbundances.abundance.values)
markerIncVirus = maxTickSizeVirus/maxAbundanceVirus
markerColorsVirus = []
for clade_id in np.unique(virusClades.clade_id):
    strainsOfClade = curClade.execute("SELECT DISTINCT vstrain_id FROM vabundances \
    WHERE clade_id = {}".format(clade_id))
    strainsOfClade = [i[0] for i in strainsOfClade.fetchall()]
    treeStrainsOfCladeAll = curTree.execute("SELECT DISTINCT tree_vstrain_id FROM tree_vstrain_order \
    WHERE vstrain_id in ({}) ORDER BY tree_vstrain_id".format(', '.join(map(str,strainsOfClade))))
    treeStrainsOfClade = []
    for (i,) in treeStrainsOfCladeAll:
        tCreate = vStrainTimes[vStrainTimes.tree_vstrain_id == i].t_creation.values[0]
        tExtinct = vStrainTimes[vStrainTimes.tree_vstrain_id == i].t_extinction.values[0]
        if tExtinct - tCreate > lifeTimeThreshold:
            treeStrainsOfClade.append(i)
    axes[1].scatter(vTreeAbundances[vTreeAbundances['tree_vstrain_id'].isin(treeStrainsOfClade)]['t'],
    vTreeAbundances[vTreeAbundances['tree_vstrain_id'].isin(treeStrainsOfClade)]['tree_vstrain_id'],
    lw=.5,
    s=vTreeAbundances[vTreeAbundances['tree_vstrain_id'].isin(treeStrainsOfClade)]['abundance'].values*markerIncVirus,
    c = np.array([Vcmap(Vnorm(np.arange(1, len(virusCladeIDs)+1, 1)))[cladeColorDict[clade_id]]]),marker='|')
    for strain in treeStrainsOfClade:
        tCreate = vStrainTimes[vStrainTimes.tree_vstrain_id == strain].t_creation.values[0]
        tExtinct = vStrainTimes[vStrainTimes.tree_vstrain_id == strain].t_extinction.values[0]
        parent = vStrainTimes[vStrainTimes.tree_vstrain_id == strain].tree_parent_vstrain_id.values[0]
        hlinecVirus.append([[tCreate, strain],[tExtinct, strain]])
        vlinecVirus.append([[tCreate, parent],[tCreate, strain]])
        hcolorsVirus.append(Vcmap(Vnorm(np.arange(1, len(virusCladeIDs)+1, 1)))[cladeColorDict[clade_id]])
        vcolorsVirus.append(Vcmap(Vnorm(np.arange(1, len(virusCladeIDs)+1, 1)))[cladeColorDict[clade_id]])
        markerColorsVirus.append(Vcmap(Vnorm(np.arange(1, len(virusCladeIDs)+1, 1)))[cladeColorDict[clade_id]])


strainLineages = LineCollection(hlinecVirus, linestyles='solid', colors=hcolorsVirus,linewidths=(0.35))
creationLines = LineCollection(vlinecVirus, linestyles='solid', colors=vcolorsVirus, linewidths=(0.35))
axes[1].add_collection(strainLineages)
axes[1].add_collection(creationLines)

fig.tight_layout()
fig.savefig(os.path.join(SCRIPT_PATH,'..','..','isolates',RUN_DIR,'virus-clades-abundances-and-trees.{}'.format(imgType)),dpi=resolve)
# fig.savefig(os.path.join('/Volumes/Yadgah','virus-clades-abundances-and-trees.{}'.format(imgType)),dpi=resolve)
# plt.show()
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
fig.savefig(os.path.join(SCRIPT_PATH,'..','..','isolates',RUN_DIR,'virus-strain-abundances.{}'.format(imgType)),dpi=resolve)
# fig.savefig(os.path.join('/Volumes/Yadgah','virus-strain-abundances.{}'.format(imgType)),dpi=resolve)
plt.close(fig)

print('Compiling stacked microbial and viral clade abundance plots')
fig, ax = plt.subplots(2,sharex=True)
fig.suptitle('Strain Abundances (run{0}-c{1}-r{2})'.format(run_id,combo_id,replicate))
axes = [ax[0], ax[1]]
microbe_stacked.plot.area(ax = axes[0],stacked=True,legend=False, linewidth=0,color=microbeColorDict,sort_columns=True)
microbe_stacked.plot(stacked=True, ax=axes[0], legend=False, color='black',sort_columns=True,linewidth=.1)
axes[0].set_ylabel(ylabel ='Microbial Strain Abundances',labelpad=15,fontsize=7)
axes[0].set_xlabel(xlabel = 'Time t',fontsize=7)
axes[0].ticklabel_format(axis = 'y',style='sci',scilimits=(0,0))
axes[0].xaxis.set_minor_locator(ticker.MultipleLocator(25))
lim = axes[0].get_ylim()
axes[0].set_ylim(0,lim[1])
virus_stacked.plot.area(ax = axes[1],stacked=True,legend=False, linewidth=0,color=virusColorDict,sort_columns=True)
virus_stacked.plot(stacked=True, ax=axes[1], legend=False, color='black',sort_columns=True,linewidth=.1)
axes[1].set_ylabel(ylabel ='Viral Strain Abundances',labelpad=15,fontsize=7)
axes[1].set_xlabel(xlabel = 'Time t',fontsize=7)
axes[1].ticklabel_format(axis = 'y',style='sci',scilimits=(0,0))
axes[1].xaxis.set_minor_locator(ticker.MultipleLocator(25))
lim = axes[1].get_ylim()
axes[1].set_ylim(0,lim[1])
fig.tight_layout()
fig.savefig(os.path.join(SCRIPT_PATH,'..','..','isolates',RUN_DIR,'microbe-virus-clades-stacked-abundances.{}'.format(imgType)),dpi=resolve)
# fig.savefig(os.path.join('/Volumes/Yadgah','microbe-virus-clades-stacked-abundances.{}'.format(imgType)),dpi=resolve)
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
fig.savefig(os.path.join(SCRIPT_PATH,'..','..','isolates',RUN_DIR,'microbe-virus-stacked-abundances.{}'.format(imgType)),dpi=resolve)
# fig.savefig(os.path.join('/Volumes/Yadgah','microbe-virus-stacked-abundances.{}'.format(imgType)),dpi=resolve)
plt.close(fig)

print('Complete!')
