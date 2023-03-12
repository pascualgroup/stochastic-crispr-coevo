#!/usr/bin/env python3

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import sqlite3
import matplotlib.cm as cm
from matplotlib.colors import Normalize
import matplotlib.colors as mc
from matplotlib.collections import LineCollection
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.ticker as ticker

resolve = 500
imgTypes = ["png"]
treeImgTypes = ["png"]
graphImgTypes = ["png"]
cladetreefigxy = (20,20)
cladetreehratio = [1,3]
maxAbundTickSize = 500
colorpalClades = 'turbo'
colorpalSpacers = 'tab20b'
lifeTimeThreshold = 75
vabundthreshold = 0.04
hyperAnalyze = True
html = False
overlay = True
sSpacing = 2
vSpacing = 2
bSpacing = 2

dir = 'crispr-sweep-7-2-2022/isolates/runID3297-c66-r47'
run = 'runID3297-c66-r47'

DBSIM_PATH = os.path.join('/Volumes','Yadgah',dir,'{}.sqlite'.format(run))

conSim = sqlite3.connect(DBSIM_PATH)
curSim = conSim.cursor()
ID = curSim.execute('SELECT combo_id,replicate FROM runs WHERE run_id = {}'.format(run_id)).fetchall()
combo_id = ID[0][0]
replicate = ID[0][1]
RUN_DIR = os.path.join('runID{0}-c{1}-r{2}'.format(run_id,combo_id,replicate))

DBMATCH_PATH = os.path.join('/Volumes','Yadgah',dir,'matches_output.sqlite') # local
DBCLADE_PATH = os.path.join('/Volumes','Yadgah',dir,'clade-abundances_output.sqlite') # local
DBTREE_PATH = os.path.join('/Volumes','Yadgah',dir,'trees_output.sqlite') # local
DBTRI_PATH = os.path.join('/Volumes','Yadgah',dir,'tripartite-networks_output.sqlite') # local
PLOT_PATH = os.path.join('/Volumes','Yadgah')

def virusCladeTreePlot(run_id,DBSIM_PATH,DBCLADE_PATH,DBTREE_PATH,colorpalClades,maxTickSizeVirus,figxy,hratio,lifeTimeThreshold,vabundthreshold):
    conSim = sqlite3.connect(DBSIM_PATH)
    curSim = conSim.cursor()
    ID = curSim.execute('SELECT combo_id,replicate FROM runs WHERE run_id = {}'.format(run_id)).fetchall()
    combo_id = ID[0][0]
    replicate = ID[0][1]
    RUN_DIR = os.path.join('runID{0}-c{1}-r{2}'.format(run_id,combo_id,replicate))
    virus_stacked = pd.read_sql_query("SELECT t,vstrain_id,abundance FROM vabundance WHERE run_id = {}".format(run_id), conSim)
    virus_total = pd.read_sql_query("SELECT t,viral_abundance FROM summary WHERE run_id = {}".format(run_id), conSim)\
    .rename(columns={"viral_abundance": "vtotal"})
    vMaxIDs = virus_stacked.set_index('t').groupby(['vstrain_id']).agg(t = ('abundance','idxmax'),\
                                vmaxAbund = ('abundance','max')).reset_index()
    vMaxIDs = virus_total.merge(vMaxIDs,on=['t'])
    vMaxIDs['vtotal'] = vabundthreshold*np.array(vMaxIDs['vtotal'])
    keepStrains = list(vMaxIDs[vMaxIDs['vmaxAbund']>vMaxIDs['vtotal']]['vstrain_id'].values)
    virus_stacked = virus_stacked[[(i in keepStrains) for i in virus_stacked.vstrain_id]].reset_index(drop=True)
    virus_stacked = virus_stacked.pivot(index='t',columns='vstrain_id',values='abundance')
    conClade = sqlite3.connect(DBCLADE_PATH)
    curClade = conClade.cursor()
    conTree = sqlite3.connect(DBTREE_PATH)
    curTree = conTree.cursor()
    print('SQLite Query: clade data')
    virusClades = pd.read_sql_query("SELECT DISTINCT clade_id, vstrain_id \
    FROM vabundances WHERE vstrain in ({})".format(', '.join(map(str,keepStrains))), conClade)
    virusCladeIDs = pd.read_sql_query("SELECT DISTINCT clade_id \
    FROM vabundances WHERE vstrain in ({})".format(', '.join(map(str,keepStrains))), conClade)
    print('SQLite Query: tree data')
    ######
    keepTreeStrains = pd.read_sql_query("SELECT tree_vstrain_id\
        FROM tree_vstrain_order WHERE vstrain_id in ({0})"\
        .format(', '.join(map(str,keepStrains))), conTree)['tree_vstrain_id']
        .values
    ######
    strainDict={}
    newID = 1
    for id in sorted(keepTreeStrains):
        strainDict[id] = newID
        newID += 1
    vStrainTimes = pd.read_sql_query(
    "SELECT tree_vstrain_id, t_creation, t_extinction, tree_parent_vstrain_id \
    FROM tree_vstrain_creation_extinction \
    WHERE tree_vstrain_id = {0}
    AND tree_parent_vstrain_id = {0}".format(', '.join(map(str,keepTreeStrains))), conTree)
    vTreeAbundances = pd.read_sql_query(
    "SELECT t, tree_vstrain_id, abundance \
    WHERE tree_vstrain_id = {0} \
    FROM tree_vabundance".format(', '.join(map(str,keepTreeStrains))), conTree)
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
        if(set(cladeDict[id[0]]).issubset(virus_stacked.columns.values)):
            columnOrder = np.append(columnOrder,cladeDict[id[0]])
        else:
            strainsTrunc = list(set(cladeDict[id[0]])\
                            .intersection(virus_stacked.columns.values))
            columnOrder = np.append(columnOrder,strainsTrunc)
    columnOrder = columnOrder.astype(int)
    virus_stacked = virus_stacked[columnOrder]
    for strain in virus_stacked.columns.values:
        clade = virusClades[virusClades['vstrain_id']==strain]['clade_id'].values[0]
        # print('strain: {}, clade: {}'.format(strain,clade))
        virusColorDict[strain] = Vcmap(Vnorm(np.arange(1, len(virusCladeIDs)+1, 1)))[cladeColorDict[clade]]
    print('Compiling stacked virus clade abundances and tree plots')
    fig, ax = plt.subplots(2,sharex=True, figsize=figxy, gridspec_kw={'height_ratios': hratio})
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
    axes[0].set_xlim(0,np.max(vTreeAbundances.t.values))
    axes[1].set_xlabel(xlabel = 'Time t',fontsize=7)
    axes[1].set_ylim(0,np.max(vStrainTimes.tree_vstrain_id)+1)
    axes[1].set_xlim(0,np.max(vTreeAbundances.t.values))
    axes[1].xaxis.set_minor_locator(ticker.MultipleLocator(25))
    axes[1].set_yticks([])
    hlinecVirus = []
    vlinecVirus = []
    hcolorsVirus = []
    vcolorsVirus = []
    maxAbundanceVirus = np.max(vTreeAbundances.abundance.values)
    markerIncVirus = maxTickSizeVirus/maxAbundanceVirus
    markerColorsVirus = []
    treeStrainsAll = []
    for clade_id in np.unique(virusClades.clade_id):
        strainsOfClade = curClade.execute("SELECT DISTINCT vstrain_id FROM vabundances \
        WHERE clade_id = {0} AND vstrain_id in ({1})".format(clade_id,', '.join(map(str,keepStrains))))
        strainsOfClade = [i[0] for i in strainsOfClade.fetchall()]
        treeStrainsOfCladeAll = curTree.execute("SELECT DISTINCT tree_vstrain_id FROM tree_vstrain_order \
        WHERE vstrain_id in ({}) ORDER BY tree_vstrain_id".format(', '.join(map(str,strainsOfClade))))
        treeStrainsOfClade = []
        for (i,) in treeStrainsOfCladeAll:
            tCreate = vStrainTimes[vStrainTimes.tree_vstrain_id == i].t_creation.values[0]
            tExtinct = vStrainTimes[vStrainTimes.tree_vstrain_id == i].t_extinction.values[0]
            # if tExtinct - tCreate > lifeTimeThreshold:
            treeStrainsOfClade.append(i)
            treeStrainsAll.append(i)
        axes[1].scatter(vTreeAbundances[vTreeAbundances['tree_vstrain_id'].isin(treeStrainsOfClade)]['t'],
        vTreeAbundances[vTreeAbundances['tree_vstrain_id'].isin(treeStrainsOfClade)]['tree_vstrain_id'],
        lw=.5,
        s=vTreeAbundances[vTreeAbundances['tree_vstrain_id'].isin(treeStrainsOfClade)]['abundance'].values*markerIncVirus,
        c = np.array([Vcmap(Vnorm(np.arange(1, len(virusCladeIDs)+1, 1)))[cladeColorDict[clade_id]]]),marker='|')
        for strain in treeStrainsOfClade:
            tCreate = vStrainTimes[vStrainTimes.tree_vstrain_id == strain].t_creation.values[0]
            tExtinct = vStrainTimes[vStrainTimes.tree_vstrain_id == strain].t_extinction.values[0]
            parent = vStrainTimes[vStrainTimes.tree_vstrain_id == strain].tree_parent_vstrain_id.values[0]
            strain = strainDict[strain]
            parent = strainDict[parent]
            hlinecVirus.append([[tCreate, strain],[tExtinct, strain]])
            vlinecVirus.append([[tCreate, parent],[tCreate, strain]])
            hcolorsVirus.append(Vcmap(Vnorm(np.arange(1, len(virusCladeIDs)+1, 1)))[cladeColorDict[clade_id]])
            vcolorsVirus.append(Vcmap(Vnorm(np.arange(1, len(virusCladeIDs)+1, 1)))[cladeColorDict[clade_id]])
            markerColorsVirus.append(Vcmap(Vnorm(np.arange(1, len(virusCladeIDs)+1, 1)))[cladeColorDict[clade_id]])
    strainLineages = LineCollection(hlinecVirus, linestyles='solid', colors=hcolorsVirus,linewidths=(0.35))
    creationLines = LineCollection(vlinecVirus, linestyles='solid', colors=vcolorsVirus, linewidths=(0.35))
    axes[1].add_collection(strainLineages)
    axes[1].add_collection(creationLines)
    return virusColorDict, cladeColorDict, cladeDict, treeStrainsAll, fig, axes

def microbeCladeTreePlot(run_id,DBSIM_PATH,DBCLADE_PATH,DBTREE_PATH,colorpalClades,maxTickSizeMicrobe,figxy,hratio):
    conSim = sqlite3.connect(DBSIM_PATH)
    curSim = conSim.cursor()
    ID = curSim.execute('SELECT combo_id,replicate FROM runs WHERE run_id = {}'.format(run_id)).fetchall()
    combo_id = ID[0][0]
    replicate = ID[0][1]
    RUN_DIR = os.path.join('runID{0}-c{1}-r{2}'.format(run_id,combo_id,replicate))
    microbe_stacked = pd.read_sql_query("SELECT t,bstrain_id,abundance FROM babundance WHERE run_id = {}".format(run_id), conSim)
    print('SQLite Query: viral abundance time series data')
    virus_stacked = pd.read_sql_query("SELECT t FROM vabundance WHERE run_id = {}".format(run_id), conSim)
    microbe_stacked = microbe_stacked[microbe_stacked.t <= max(virus_stacked.t)]
    microbe_stacked = microbe_stacked.pivot(index='t',columns='bstrain_id',values='abundance')
    conClade = sqlite3.connect(DBCLADE_PATH)
    curClade = conClade.cursor()
    conTree = sqlite3.connect(DBTREE_PATH)
    curTree = conTree.cursor()
    print('SQLite Query: clade data')
    microbeClades = pd.read_sql_query("SELECT DISTINCT clade_id, bstrain_id \
    FROM babundances", conClade)
    microbeCladeIDs = pd.read_sql_query("SELECT DISTINCT clade_id \
    FROM babundances", conClade)
    print('SQLite Query: tree data')
    bStrainTimes = pd.read_sql_query(
    "SELECT tree_bstrain_id, t_creation, t_extinction, tree_parent_bstrain_id \
    FROM tree_bstrain_creation_extinction", conTree)
    bTreeAbundances = pd.read_sql_query(
    "SELECT t, tree_bstrain_id, abundance \
    FROM tree_babundance", conTree)
    bTreeAbundances = bTreeAbundances[bTreeAbundances.t <= max(virus_stacked.t)]
    Mcmap = cm.get_cmap('{}'.format(colorpalClades)) #usually cmap = 'turbo'
    # Get normalize function (takes data in range [vmin, vmax] -> [0, 1])
    Mnorm = Normalize(vmin=1, vmax=len(microbeCladeIDs))
    microbeColorDict = {}
    cladeColorDict = {}
    cladeDict = {}
    colorInd = 0
    for clade_id in microbeCladeIDs.values:
        # print(id[0])
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
    print('Compiling stacked microbe clade abundances and tree plots')
    fig, ax = plt.subplots(2,sharex=True, figsize=figxy, gridspec_kw={'height_ratios': hratio})
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
    axes[0].set_xlim(0,np.max(bTreeAbundances.t.values))
    axes[1].set_xlabel(xlabel = 'Time t',fontsize=7)
    axes[1].set_ylim(0,np.max(bStrainTimes.tree_bstrain_id)+1)
    axes[0].set_xlim(0,np.max(bTreeAbundances.t.values))
    axes[1].xaxis.set_minor_locator(ticker.MultipleLocator(25))
    axes[1].set_yticks([])
    hlinecMicrobe = []
    vlinecMicrobe = []
    hcolorsMicrobe = []
    vcolorsMicrobe = []
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
        lw=.8,
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
    return microbeColorDict, cladeColorDict, cladeDict, fig, axes

def susceptibleCladeTreePlot(run_id,DBSIM_PATH,DBCLADE_PATH,DBTREE_PATH,DBMATCH_PATH,colorpalClades,maxTickSizeMicrobe,figxy,hratio):
    conSim = sqlite3.connect(DBSIM_PATH)
    curSim = conSim.cursor()
    ID = curSim.execute('SELECT combo_id,replicate FROM runs WHERE run_id = {}'.format(run_id)).fetchall()
    combo_id = ID[0][0]
    replicate = ID[0][1]
    RUN_DIR = os.path.join('runID{0}-c{1}-r{2}'.format(run_id,combo_id,replicate))
    microbe_stacked = pd.read_sql_query("SELECT t,bstrain_id,abundance FROM babundance WHERE run_id = {}".format(run_id), conSim)
    print('SQLite Query: viral abundance time series data')
    virus_stacked = pd.read_sql_query("SELECT t FROM vabundance WHERE run_id = {}".format(run_id), conSim)
    microbe_stacked = microbe_stacked[microbe_stacked.t <= max(virus_stacked.t)]
    microbe_stacked = microbe_stacked.pivot(index='t',columns='bstrain_id',values='abundance')
    conClade = sqlite3.connect(DBCLADE_PATH)
    curClade = conClade.cursor()
    conTree = sqlite3.connect(DBTREE_PATH)
    curTree = conTree.cursor()
    conMatch = sqlite3.connect(DBMATCH_PATH)
    curMatch = conMatch.cursor()
    print('SQLite Query: clade data')
    microbeClades = pd.read_sql_query("SELECT DISTINCT clade_id, bstrain_id \
    FROM babundances", conClade)
    microbeCladeIDs = pd.read_sql_query("SELECT DISTINCT clade_id \
    FROM babundances", conClade)
    print('SQLite Query: tree data')
    bStrainTimes = pd.read_sql_query(
    "SELECT tree_bstrain_id, t_creation, t_extinction, tree_parent_bstrain_id \
    FROM tree_bstrain_creation_extinction", conTree)
    bstrain0vstrain = pd.read_sql_query(
    "SELECT t, bstrain_id, vstrain_id \
    FROM bstrain_to_vstrain_0matches", conMatch)
    vAbunds = pd.read_sql_query("SELECT t,vstrain_id,abundance FROM vabundance WHERE run_id = {}".format(run_id), conSim)
    bstrain0vstrain = bstrain0vstrain.merge(vAbunds,on=['t','vstrain_id'])\
    .groupby(['t','bstrain_id']).agg(abundance=('abundance','sum')).reset_index()
    bTreeAbundances = pd.read_sql_query(
    "SELECT t, tree_bstrain_id, abundance \
    FROM tree_babundance", conTree)
    bTreeAbundances = bTreeAbundances[bTreeAbundances.t <= max(virus_stacked.t)]
    Mcmap = cm.get_cmap('{}'.format(colorpalClades)) #usually cmap = 'turbo'
    # Get normalize function (takes data in range [vmin, vmax] -> [0, 1])
    Mnorm = Normalize(vmin=1, vmax=len(microbeCladeIDs))
    microbeColorDict = {}
    cladeColorDict = {}
    cladeDict = {}
    colorInd = 0
    for clade_id in microbeCladeIDs.values:
        # print(id[0])
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
    print('Compiling stacked microbe clade abundances and tree plots')
    fig, ax = plt.subplots(2,sharex=True, figsize=figxy, gridspec_kw={'height_ratios': hratio})
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
    axes[0].set_xlim(0,np.max(bTreeAbundances.t.values))
    axes[1].set_xlabel(xlabel = 'Time t',fontsize=7)
    axes[1].set_ylim(0,np.max(bStrainTimes.tree_bstrain_id)+1)
    axes[0].set_xlim(0,np.max(bTreeAbundances.t.values))
    axes[1].xaxis.set_minor_locator(ticker.MultipleLocator(25))
    axes[1].set_yticks([])
    hlinecMicrobe = []
    vlinecMicrobe = []
    hcolorsMicrobe = []
    vcolorsMicrobe = []
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
        lw=.8,
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
    return fig, axes
