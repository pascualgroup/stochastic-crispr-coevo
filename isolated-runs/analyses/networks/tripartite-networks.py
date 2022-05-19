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
import matplotlib.ticker as ticker
import matplotlib.cm as cm
from matplotlib.colors import Normalize
import matplotlib.colors as mc
from matplotlib.patches import Rectangle
import colorsys
from matplotlib.collections import LineCollection
from mpl_toolkits.axes_grid1 import make_axes_locatable

run_id = sys.argv[1]
print('Compiling plots and networks for run ID: {0}'.format(sys.argv[1]))
print('is this  what you are looking for: {}'.format(len(sys.argv)))
times = [float(sys.argv[i]) for i in list(range(2,len(sys.argv)-1))]
print('...for times: {0}...'.format(', '.join(map(str,times))))
threshold = float(sys.argv[-1])/100
print('...and with with {0}% abundances removed...'.format(float(sys.argv[-1])))
resolve = 500
imgTypes = ["pdf","png"]
cladetreefigxy = (20,20)
cladetreehratio = [1,3]
maxAbundTickSize = 500
colorpalClades = 'turbo'
colorpalSpacers = 'tab20b'

SCRIPT_PATH = os.path.abspath(os.path.dirname(__file__)) # cluster

dir = 'crispr-sweep-7-2-2022/isolates/runID3297-c66-r47'
run = 'runID3297-c66-r47'

# DBSIM_PATH = os.path.join(SCRIPT_PATH,'..','..','..','sweep_db_gathered.sqlite')
# DBSIM_PATH = os.path.join('/Volumes','Yadgah','crispr-sweep-7-2-2022/isolates/runID1723-c35-r23/runID1723-c35-r23.sqlite')
DBSIM_PATH = os.path.join('/Volumes','Yadgah',dir,'{}.sqlite'.format(run))

conSim = sqlite3.connect(DBSIM_PATH)
curSim = conSim.cursor()
ID = curSim.execute('SELECT combo_id,replicate FROM runs WHERE run_id = {}'.format(run_id)).fetchall()
combo_id = ID[0][0]
replicate = ID[0][1]
RUN_DIR = os.path.join('runID{0}-c{1}-r{2}'.format(run_id,combo_id,replicate))

# DBMATCH_PATH = os.path.join(SCRIPT_PATH,'..','..','isolates',,RUN_DIR,'matches_output.sqlite') # cluster
DBMATCH_PATH = os.path.join('/Volumes','Yadgah',dir,'matches_output.sqlite') # local
# DBCLADE_PATH = os.path.join(SCRIPT_PATH,'..','..','isolates',RUN_DIR,'clade-abundances_output.sqlite') # cluster
DBCLADE_PATH = os.path.join('/Volumes','Yadgah',dir,'clade-abundances_output.sqlite') # local. run_id fixed; for testing
# DBTREE_PATH = os.path.join(SCRIPT_PATH,'..','..','isolates',RUN_DIR,'trees_output.sqlite') # cluster
DBTREE_PATH = os.path.join('/Volumes','Yadgah',dir,'trees_output.sqlite') # local. run_id fixed; for testing


conMatch = sqlite3.connect(DBMATCH_PATH)
curMatch = conMatch.cursor()
conClade = sqlite3.connect(DBCLADE_PATH)
curClade = conClade.cursor()
conTree = sqlite3.connect(DBTREE_PATH)
curTree = conTree.cursor()


print('SQLite Query: microbial abundance time series data')
microbe_stacked = pd.read_sql_query("SELECT t,bstrain_id,abundance FROM babundance WHERE run_id = {}".format(run_id), conSim)
microbe_stacked = microbe_stacked.pivot(index='t',columns='bstrain_id',values='abundance')
print('SQLite Query: viral abundance time series data')
virus_stacked = pd.read_sql_query("SELECT t,vstrain_id,abundance FROM vabundance WHERE run_id = {}".format(run_id), conSim)
virus_stacked = virus_stacked.pivot(index='t',columns='vstrain_id',values='abundance')
if len(microbe_stacked.index) > len(virus_stacked.index):
    microbe_stacked.drop(index=list(set(microbe_stacked.index)-set(virus_stacked.index)),inplace=True)


print('SQLite Query: clade data')
microbeClades = pd.read_sql_query("SELECT DISTINCT clade_id, bstrain_id \
FROM babundances", conClade)
microbeCladeIDs = pd.read_sql_query("SELECT DISTINCT clade_id \
FROM babundances", conClade)
virusClades = pd.read_sql_query("SELECT DISTINCT clade_id, vstrain_id \
FROM vabundances", conClade)
virusCladeIDs = pd.read_sql_query("SELECT DISTINCT clade_id \
FROM vabundances", conClade)

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

def microbeCladeTreePlot(microbe_stacked,colorpalClades,maxTickSizeMicrobe,figxy,hratio):
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
microbeColorDict, cladeColorDict, cladeDict, fig, axes = microbeCladeTreePlot(microbe_stacked,colorpalClades,maxAbundTickSize,cladetreefigxy,cladetreehratio)

for time in times:
    print('Compiling tripartite (virus-spacer-microbe) networks at time {0}'.format(time))
    vThresholdAbundance = threshold*sum([abund[0] for abund in
        curSim.execute("SELECT abundance FROM vabundance \
        WHERE t = {0}".format(time)).fetchall()])
    bThresholdAbundance = threshold*sum([abund[0] for abund
        in curSim.execute("SELECT abundance FROM babundance \
        WHERE t = {0}".format(time)).fetchall()])
    vstrainThreshold = [vstrain_id[0] for vstrain_id in
        curSim.execute("SELECT vstrain_id, abundance FROM vabundance \
        WHERE t = {0}".format(time)) if vstrain_id[1] >= vThresholdAbundance]
    bstrainThreshold = [bstrain_id[0] for bstrain_id in
        curSim.execute("SELECT bstrain_id, abundance FROM babundance \
        WHERE t = {0}".format(time)) if bstrain_id[1] >= bThresholdAbundance]
    vstrains0=[vstrain_id[0] for vstrain_id in
        curMatch.execute("SELECT DISTINCT vstrain_id FROM bstrain_to_vstrain_0matches \
        WHERE t = {0} AND vstrain_id in ({1})".format(time,', '.join(map(str,vstrainThreshold))))]
    # phylogenetically order virus strains
    vstrains0 = [vstrain_id[0] for vstrain_id in
        curTree.execute("SELECT vstrain_id FROM tree_vstrain_order \
        WHERE vstrain_id in ({}) \
        ORDER BY tree_vstrain_id".format(', '.join(map(str,vstrains0))))]
    nonMatchDF = pd.read_sql_query(
        "SELECT vstrain_id, bstrain_id \
        FROM bstrain_to_vstrain_0matches WHERE t = {0} \
        AND vstrain_id in ({1}) AND bstrain_id in ({2})".format(time,
        ', '.join(map(str,vstrainThreshold)),', '.join(map(str,bstrainThreshold))), conMatch)
    nonMatchDF['spacer_id'] = len(nonMatchDF) * [0]
    singleMatchDF = pd.read_sql_query(
        "SELECT vstrain_id, bstrain_id \
        FROM bstrain_to_vstrain_matches WHERE t = {0} AND match_length = 1 \
        AND vstrain_id in ({1}) AND bstrain_id in ({2}) \
        ORDER BY time_specific_match_id".format(time,
        ', '.join(map(str,vstrains0)),', '.join(map(str,bstrainThreshold))), conMatch)

    matchIDs = [match_id[0] for match_id in
        curMatch.execute("SELECT time_specific_match_id FROM bstrain_to_vstrain_matches \
        WHERE t = {0} AND match_length = 1 \
        AND vstrain_id in ({1}) AND bstrain_id in ({2}) \
        ORDER BY time_specific_match_id".format(time,
        ', '.join(map(str,vstrains0)),', '.join(map(str,bstrainThreshold))))]
    spacerIDs = [spacer_id[0] for spacer_id in
        curMatch.execute("SELECT spacer_id FROM matches_spacers \
        WHERE t = {0} \
        AND time_specific_match_id in ({1}) \
        ORDER BY time_specific_match_id".format(time,
        ', '.join(map(str,matchIDs))))]
    spacerDict = {i:j for (i,j) in zip(np.unique(spacerIDs),list(range(1,len(np.unique(spacerIDs))+1)))}
    orderedSpacerIDs = [spacerDict[i] for i in spacerIDs]
    singleMatchDF['spacer_id'] = orderedSpacerIDs
    spacerIDs = list(np.unique(spacerIDs))
    spacerIDs.insert(0,0)
    Mcmap = cm.get_cmap('{}'.format(colorpalSpacers))
    spacerCMAP = mc.LinearSegmentedColormap.from_list(
        'spacers', [Mcmap(i) for i in range(len(spacerIDs))], len(spacerIDs))
    # doubleMatchDF = pd.read_sql_query(
    #     "SELECT vstrain_id, bstrain_id, match_length \
    #     FROM bstrain_to_vstrain_matches WHERE t = {0} \
    #     AND vstrain_id in ({1}) AND bstrain_id in ({2}) \
    #     AND match_length = 2\
    #     ORDER BY match_length".format(time,
    #     ', '.join(map(str,vstrains0)),', '.join(map(str,bstrainThreshold))), conMatch)

    matchNetwork = pd.concat([nonMatchDF,singleMatchDF])
    matchNetwork = matchNetwork.pivot(index='vstrain_id',columns='bstrain_id',values='spacer_id')
    # singleMatchNetwork = singleMatchDF.pivot(index='vstrain_id',columns='bstrain_id',values='spacer_id')

    bstrains=[bstrain_id[0] for bstrain_id in
        curMatch.execute("SELECT DISTINCT bstrain_id FROM bstrain_to_vstrain_0matches \
        WHERE t = {0} \
        AND bstrain_id in ({1}) AND vstrain_id in ({2})"
        .format(time,', '.join(map(str,bstrainThreshold)), ', '.join(map(str,vstrains0))))]
    bstrains = list(set().union(bstrains,
        [bstrain_id[0] for bstrain_id in
            curMatch.execute("SELECT DISTINCT bstrain_id FROM bstrain_to_vstrain_matches \
            WHERE t = {0} \
            AND bstrain_id in ({1}) AND vstrain_id in ({2})\
            AND match_length = 1"
            .format(time,', '.join(map(str,bstrainThreshold)),', '.join(map(str,vstrains0))))]
        ))
    # print(matchNetwork)
    # print(sort(bstrains))
    matchNetwork = matchNetwork.reindex(index = sorted(vstrains0))
    matchNetwork = matchNetwork[sorted(bstrains)]

    fig2, ax2 = plt.subplots(1,sharex=True)
    # cbar_kws=dict(pad=0.01,shrink=0.86, drawedges=True, ticks=range(int(maxfold)+1))
    ax2 = sns.heatmap(matchNetwork,cmap=spacerCMAP,mask=matchNetwork.isnull(),
    xticklabels=True, yticklabels=True)
    colorbar = ax2.collections[0].colorbar
    colorbar.set_ticks(list(
    np.linspace(.5*(len(spacerIDs)-1)/len(spacerIDs),(len(spacerIDs)-0.5)
    *(len(spacerIDs)-1)/len(spacerIDs),len(spacerIDs)
    )
    ))
    colorbar.set_ticklabels(spacerIDs)
    colorbar.set_label("Spacer ID")
    ax2.set_xlabel(xlabel ='Microbe Strain ID',labelpad=10)
    ax2.set_ylabel(ylabel = 'Virus Strain ID',labelpad=10)
    ax2.set_title('runID{0}-c{1}-r{2}\nTripartite Match Network at t = {3} (Time-ordered)'.format(run_id,combo_id,replicate,time),pad=20,fontsize=10)
    # ax2 = sns.heatmap(matchNetwork,cmap=spacerCMAP, xticklabels=True, yticklabels=True)
    ax2.set_yticklabels(ax2.get_ymajorticklabels(), fontsize = 2)
    ax2.set_xticklabels(ax2.get_xmajorticklabels(), fontsize = 5)
    ax2.yaxis.set_tick_params(width=.25)
    ax2.yaxis.set_tick_params(width=.25)
    # fig2.savefig(os.path.join(SCRIPT_PATH,'..','..','isolates',RUN_DIR,'tripartite-match-network-time-ordered-time{0}.{1}'.format(int(np.floor(time)),imgType)),dpi=resolve)
    for imgType in imgTypes:
        fig2.savefig(os.path.join('/Volumes/Yadgah','tripartite-match-network-time-ordered-time{0}.{1}'.format(int(np.floor(time)),imgType)),dpi=resolve)
    plt.close(fig2)

    # phylogenetically order microbe strains
    bstrains = [bstrain_id[0] for bstrain_id in
        curTree.execute("SELECT bstrain_id FROM tree_bstrain_order \
        WHERE bstrain_id in ({}) \
        ORDER BY tree_bstrain_id".format(', '.join(map(str,bstrains))))]
    treebstrains = {tree:actual for (tree,actual) in
    zip([bstrain_id[0]
            for bstrain_id in
                curTree.execute("SELECT tree_bstrain_id FROM tree_bstrain_order \
                WHERE bstrain_id in ({}) \
                ORDER BY tree_bstrain_id".format(', '.join(map(str,bstrains))))],bstrains
                )
        }
    axes[0].axvline(x=time, color='black', ls=':', lw=.7)
    axes[1].axvline(x=time, color='black', ls=':', lw=.7)
    for treestrain in treebstrains.keys():
        axes[1].text(time, treestrain, '{}'.format(treebstrains[treestrain]), fontsize=10, fontweight='bold')

    matchNetwork = matchNetwork.reindex(index = vstrains0[::-1])
    matchNetwork = matchNetwork[bstrains]

    fig2, ax2 = plt.subplots(1,sharex=True)
    # cbar_kws=dict(pad=0.01,shrink=0.86, drawedges=True, ticks=range(int(maxfold)+1))
    ax2 = sns.heatmap(matchNetwork,cmap=spacerCMAP,mask=matchNetwork.isnull(),
    xticklabels=True, yticklabels=True)
    colorbar = ax2.collections[0].colorbar
    colorbar.set_ticks(list(
    np.linspace(.5*(len(spacerIDs)-1)/len(spacerIDs),(len(spacerIDs)-0.5)
    *(len(spacerIDs)-1)/len(spacerIDs),len(spacerIDs)
    )
    ))
    colorbar.set_ticklabels(spacerIDs)
    colorbar.set_label("Spacer ID")
    ax2.set_xlabel(xlabel ='Microbe Strain ID',labelpad=10)
    ax2.set_ylabel(ylabel = 'Virus Strain ID',labelpad=10)
    ax2.set_title('runID{0}-c{1}-r{2}\nTripartite Match Network at t = {3} (Phylogenetically Ordered)'.format(run_id,combo_id,replicate,time),pad=20,fontsize=10)
    # ax2 = sns.heatmap(matchNetwork,cmap=spacerCMAP, xticklabels=True, yticklabels=True)
    ax2.set_yticklabels(ax2.get_ymajorticklabels(), fontsize = 2)
    ax2.set_xticklabels(ax2.get_xmajorticklabels(), fontsize = 5)
    ax2.yaxis.set_tick_params(width=.25)
    ax2.yaxis.set_tick_params(width=.25)
    McmapClades = cm.get_cmap('{}'.format(colorpalClades))
    Mnorm = Normalize(vmin=1, vmax=len(microbeCladeIDs))
    totalnumstrains = 0
    for clade_id in np.unique(microbeClades.clade_id)[::-1]:
        substrains = set.intersection(set(cladeDict[clade_id]),bstrains)
        if len(substrains) == 0:
            continue
        if totalnumstrains == 0:
            ax2.add_patch(Rectangle((0, 0), len(substrains), len(vstrains0), fill=False,
            edgecolor=McmapClades(Mnorm(np.arange(1, len(microbeCladeIDs)+1, 1)))[cladeColorDict[clade_id]],
            lw=3))
        else:
            ax2.add_patch(Rectangle((totalnumstrains+.1, 0), len(substrains)-.1, len(vstrains0), fill=False,
            edgecolor=McmapClades(Mnorm(np.arange(1, len(microbeCladeIDs)+1, 1)))[cladeColorDict[clade_id]],
            lw=3))
        totalnumstrains  += len(substrains)

    # fig2.savefig(os.path.join(SCRIPT_PATH,'..','..','isolates',RUN_DIR,'tripartite-match-network-phylo-ordered-time{0}.{1}'.format(int(np.floor(time)),imgType)),dpi=resolve)
    for imgType in imgTypes:
        fig2.savefig(os.path.join('/Volumes/Yadgah','tripartite-match-network-phylo-ordered-time{0}.{1}'.format(int(np.floor(time)),imgType)),dpi=resolve)
    plt.close(fig2)

    vstrainsAbundOrder=[vstrain_id[0] for vstrain_id in
        curSim.execute("SELECT vstrain_id FROM vabundance \
        WHERE t = {0} AND vstrain_id in ({1}) \
        ORDER BY abundance".format(time,', '.join(map(str,vstrains0))))]
    bstrainsAbundOrder=[bstrain_id[0] for bstrain_id in
        curSim.execute("SELECT bstrain_id FROM babundance \
        WHERE t = {0} AND bstrain_id in ({1}) \
        ORDER BY abundance".format(time,', '.join(map(str,bstrains))))]
    matchNetwork = matchNetwork.reindex(index = vstrainsAbundOrder[::-1])
    matchNetwork = matchNetwork[bstrainsAbundOrder]
    fig2, ax2 = plt.subplots(1,sharex=True)
    # cbar_kws=dict(pad=0.01,shrink=0.86, drawedges=True, ticks=range(int(maxfold)+1))
    ax2 = sns.heatmap(matchNetwork,cmap=spacerCMAP,mask=matchNetwork.isnull(),
    xticklabels=True, yticklabels=True)
    colorbar = ax2.collections[0].colorbar
    colorbar.set_ticks(list(
    np.linspace(.5*(len(spacerIDs)-1)/len(spacerIDs),(len(spacerIDs)-0.5)
    *(len(spacerIDs)-1)/len(spacerIDs),len(spacerIDs)
    )
    ))
    colorbar.set_ticklabels(spacerIDs)
    colorbar.set_label("Spacer ID")
    ax2.set_xlabel(xlabel ='Microbe Strain ID',labelpad=10)
    ax2.set_ylabel(ylabel = 'Virus Strain ID',labelpad=10)
    ax2.set_title('runID{0}-c{1}-r{2}\nTripartite Match Network at t = {3} (Abundance-ordered)'.format(run_id,combo_id,replicate,time),pad=20,fontsize=10)
    # ax2 = sns.heatmap(matchNetwork,cmap=spacerCMAP, xticklabels=True, yticklabels=True)
    ax2.set_yticklabels(ax2.get_ymajorticklabels(), fontsize = 2)
    ax2.set_xticklabels(ax2.get_xmajorticklabels(), fontsize = 5)
    ax2.yaxis.set_tick_params(width=.25)
    ax2.yaxis.set_tick_params(width=.25)
    # fig2.savefig(os.path.join(SCRIPT_PATH,'..','..','isolates',RUN_DIR,'tripartite-match-network-abundance-ordered-time{0}.{1}'.format(int(np.floor(time)),imgType)),dpi=resolve)
    for imgType in imgTypes:
        fig2.savefig(os.path.join('/Volumes/Yadgah','tripartite-match-network-abundance-ordered-time{0}.{1}'.format(int(np.floor(time)),imgType)),dpi=resolve)
    plt.close(fig2)

    pairwisebAbundDict = {}
    for bstrain in matchNetwork.columns.values:
        vstrains = matchNetwork[bstrain][matchNetwork[bstrain] == 0].index.values
        if len(vstrains) == 0:
            bAbund = curSim.execute("SELECT abundance FROM babundance \
                WHERE t = {0} AND bstrain_id in ({1})".format(time,bstrain)).fetchall()
            bAbund = bAbund[0][0]
            pairwisebAbundDict[bstrain] = [0, bAbund]
            continue
        vAbunds=[abund[0] for abund in
            curSim.execute("SELECT abundance FROM vabundance \
            WHERE t = {0} AND vstrain_id in ({1})".format(time,', '.join(map(str,vstrains))))]
        bAbund = curSim.execute("SELECT abundance FROM babundance \
            WHERE t = {0} AND bstrain_id in ({1})".format(time,bstrain)).fetchall()
        bAbund = bAbund[0][0]
        pairwisebAbundDict[bstrain] = [bAbund*sum(vAbunds), bAbund]

    pairwiseAbund = pd.DataFrame.from_dict(pairwisebAbundDict,orient='index',columns=['pairwiseAbundance', 'bAbundance'])
    bstrainPairwiseOrder = pairwiseAbund.sort_values(
        by=['pairwiseAbundance', 'bAbundance']).index.values

    pairwisevAbundDict = {}
    for vstrain in matchNetwork.index.values:
        bstrains = matchNetwork.loc[:,(matchNetwork.loc[[vstrain]] == 0).any()].loc[[vstrain]].columns.values
        if len(bstrains) == 0:
            vAbund = curSim.execute("SELECT abundance FROM vabundance \
                WHERE t = {0} AND vstrain_id in ({1})".format(time,vstrain)).fetchall()
            vAbund = vAbund[0][0]
            pairwisevAbundDict[vstrain] = [0, vAbund]
            continue
        bAbunds=[abund[0] for abund in
            curSim.execute("SELECT abundance FROM babundance \
            WHERE t = {0} AND bstrain_id in ({1})".format(time,', '.join(map(str,bstrains))))]
        vAbund = curSim.execute("SELECT abundance FROM vabundance \
            WHERE t = {0} AND vstrain_id in ({1})".format(time,vstrain)).fetchall()
        vAbund = vAbund[0][0]
        pairwisevAbundDict[vstrain] = [vAbund*sum(bAbunds), vAbund]

    pairwiseAbund = pd.DataFrame.from_dict(pairwisevAbundDict,orient='index',columns=['pairwiseAbundance', 'vAbundance'])
    vstrainPairwiseOrder = pairwiseAbund.sort_values(
        by=['pairwiseAbundance', 'vAbundance']).index.values
    matchNetwork = matchNetwork.reindex(index = vstrainPairwiseOrder[::-1])
    matchNetwork = matchNetwork[bstrainPairwiseOrder]

    fig2, ax2 = plt.subplots(1,sharex=True)
    # cbar_kws=dict(pad=0.01,shrink=0.86, drawedges=True, ticks=range(int(maxfold)+1))
    ax2 = sns.heatmap(matchNetwork,cmap=spacerCMAP,mask=matchNetwork.isnull(),
    xticklabels=True, yticklabels=True)
    colorbar = ax2.collections[0].colorbar
    colorbar.set_ticks(list(
    np.linspace(.5*(len(spacerIDs)-1)/len(spacerIDs),(len(spacerIDs)-0.5)
    *(len(spacerIDs)-1)/len(spacerIDs),len(spacerIDs)
    )
    ))
    colorbar.set_ticklabels(spacerIDs)
    colorbar.set_label("Spacer ID")
    ax2.set_xlabel(xlabel ='Microbe Strain ID',labelpad=10)
    ax2.set_ylabel(ylabel = 'Virus Strain ID',labelpad=10)
    ax2.set_title('runID{0}-c{1}-r{2}\nTripartite Match Network at t = {3} (0-match Abundance-ordered)'.format(run_id,combo_id,replicate,time),pad=20,fontsize=10)
    # ax2 = sns.heatmap(matchNetwork,cmap=spacerCMAP, xticklabels=True, yticklabels=True)
    ax2.set_yticklabels(ax2.get_ymajorticklabels(), fontsize = 2)
    ax2.set_xticklabels(ax2.get_xmajorticklabels(), fontsize = 5)
    ax2.yaxis.set_tick_params(width=.25)
    ax2.yaxis.set_tick_params(width=.25)
    # fig2.savefig(os.path.join(SCRIPT_PATH,'..','..','isolates',RUN_DIR,'tripartite-match-network-abundance-ordered-time{0}.{1}'.format(int(np.floor(time)),imgType)),dpi=resolve)
    for imgType in imgTypes:
        fig2.savefig(os.path.join('/Volumes/Yadgah','tripartite-match-network-0match-abundance-ordered-time{0}.{1}'.format(int(np.floor(time)),imgType)),dpi=resolve)
    plt.close(fig2)

    pairwisebAbundDict = {}
    for bstrain in matchNetwork.columns.values:
        vstrains = matchNetwork[bstrain].index.values
        if len(vstrains) == 0:
            bAbund = curSim.execute("SELECT abundance FROM babundance \
                WHERE t = {0} AND bstrain_id in ({1})".format(time,bstrain)).fetchall()
            bAbund = bAbund[0][0]
            pairwisebAbundDict[bstrain] = [0, bAbund]
            continue
        vAbunds=[abund[0] for abund in
            curSim.execute("SELECT abundance FROM vabundance \
            WHERE t = {0} AND vstrain_id in ({1})".format(time,', '.join(map(str,vstrains))))]
        bAbund = curSim.execute("SELECT abundance FROM babundance \
            WHERE t = {0} AND bstrain_id in ({1})".format(time,bstrain)).fetchall()
        bAbund = bAbund[0][0]
        pairwisebAbundDict[bstrain] = [bAbund*sum(vAbunds), bAbund]

    pairwiseAbund = pd.DataFrame.from_dict(pairwisebAbundDict,orient='index',columns=['pairwiseAbundance', 'bAbundance'])
    bstrainPairwiseOrder = pairwiseAbund.sort_values(
        by=['pairwiseAbundance', 'bAbundance']).index.values

    pairwisevAbundDict = {}
    for vstrain in matchNetwork.index.values:
        bstrains = matchNetwork.loc[[vstrain]].columns.values
        if len(bstrains) == 0:
            vAbund = curSim.execute("SELECT abundance FROM vabundance \
                WHERE t = {0} AND vstrain_id in ({1})".format(time,vstrain)).fetchall()
            vAbund = vAbund[0][0]
            pairwisevAbundDict[vstrain] = [0, vAbund]
            continue
        bAbunds=[abund[0] for abund in
            curSim.execute("SELECT abundance FROM babundance \
            WHERE t = {0} AND bstrain_id in ({1})".format(time,', '.join(map(str,bstrains))))]
        vAbund = curSim.execute("SELECT abundance FROM vabundance \
            WHERE t = {0} AND vstrain_id in ({1})".format(time,vstrain)).fetchall()
        vAbund = vAbund[0][0]
        pairwisevAbundDict[vstrain] = [vAbund*sum(bAbunds), vAbund]

    pairwiseAbund = pd.DataFrame.from_dict(pairwisevAbundDict,orient='index',columns=['pairwiseAbundance', 'vAbundance'])
    vstrainPairwiseOrder = pairwiseAbund.sort_values(
        by=['pairwiseAbundance', 'vAbundance']).index.values
    matchNetwork = matchNetwork.reindex(index = vstrainPairwiseOrder[::-1])
    matchNetwork = matchNetwork[bstrainPairwiseOrder]

    fig2, ax2 = plt.subplots(1,sharex=True)
    # cbar_kws=dict(pad=0.01,shrink=0.86, drawedges=True, ticks=range(int(maxfold)+1))
    ax2 = sns.heatmap(matchNetwork,cmap=spacerCMAP,mask=matchNetwork.isnull(),
    xticklabels=True, yticklabels=True)
    colorbar = ax2.collections[0].colorbar
    colorbar.set_ticks(list(
    np.linspace(.5*(len(spacerIDs)-1)/len(spacerIDs),(len(spacerIDs)-0.5)
    *(len(spacerIDs)-1)/len(spacerIDs),len(spacerIDs)
    )
    ))
    colorbar.set_ticklabels(spacerIDs)
    colorbar.set_label("Spacer ID")
    ax2.set_xlabel(xlabel ='Microbe Strain ID',labelpad=10)
    ax2.set_ylabel(ylabel = 'Virus Strain ID',labelpad=10)
    ax2.set_title('runID{0}-c{1}-r{2}\nTripartite Match Network at t = {3} (Pairwise Abundance-ordered)'.format(run_id,combo_id,replicate,time),pad=20,fontsize=10)
    # ax2 = sns.heatmap(matchNetwork,cmap=spacerCMAP, xticklabels=True, yticklabels=True)
    ax2.set_yticklabels(ax2.get_ymajorticklabels(), fontsize = 2)
    ax2.set_xticklabels(ax2.get_xmajorticklabels(), fontsize = 5)
    ax2.yaxis.set_tick_params(width=.25)
    ax2.yaxis.set_tick_params(width=.25)
    # fig2.savefig(os.path.join(SCRIPT_PATH,'..','..','isolates',RUN_DIR,'tripartite-match-network-abundance-ordered-time{0}.{1}'.format(int(np.floor(time)),imgType)),dpi=resolve)
    for imgType in imgTypes:
        fig2.savefig(os.path.join('/Volumes/Yadgah','tripartite-match-network-pairwise-abundance-ordered-time{0}.{1}'.format(int(np.floor(time)),imgType)),dpi=resolve)
    plt.close(fig2)



fig.tight_layout()
for imgType in imgTypes:
    fig.savefig(os.path.join('/Volumes/Yadgah','microbe-clades-abundances-and-trees.{0}'.format(imgType)),dpi=resolve)
# fig.savefig(os.path.join(SCRIPT_PATH,'..','..','isolates',RUN_DIR,'microbe-clades-abundances-and-trees.{0}'.format(imgType)),dpi=resolve)













# Microbe-to-microbe network
# vstrains0=[vstrain_id[0] for vstrain_id in
#     curMatch.execute("SELECT DISTINCT vstrain_id FROM bstrain_to_vstrain_0matches \
#     WHERE t = {0}".format(time))]
# nonMatchDF = pd.read_sql_query(
#     "SELECT vstrain_id, bstrain_id \
#     FROM bstrain_to_vstrain_0matches WHERE t = {0}".format(time), conMatch)
# singleMatchDF = pd.read_sql_query(
#     "SELECT vstrain_id, bstrain_id \
#     FROM bstrain_to_vstrain_matches WHERE t = {0} AND match_length = 1 \
#     AND vstrain_id in ({1}) \
#     ORDER BY time_specific_match_id".format(time,', '.join(map(str,vstrains0))),
#     conMatch)
#
# tripartiteMVM = pd.DataFrame(columns = ['bstrainID_0', 'bstrainID_1', 'vstrainID'])
# vfreqDict = {}
# mmDict = {'bstrainID_0':[],'bstrainID_1':[]}
# vtotal = sum([abund[0] for abund in
#     curSim.execute("SELECT abundance FROM vabundance \
#     WHERE t = {0}".format(time)).fetchall()])
# btotal = sum([abund[0] for abund
#     in curSim.execute("SELECT abundance FROM babundance \
#     WHERE t = {0}".format(time)).fetchall()])
# for bstrain1 in np.unique(singleMatchDF.bstrain_id.values):
#     for vstrain1 in singleMatchDF[singleMatchDF.bstrain_id == bstrain1].vstrain_id.values:
#         for bstrain0 in nonMatchDF[nonMatchDF.vstrain_id == vstrain1].bstrain_id.values:
#             if (bstrain0 not in mmDict['bstrainID_0']):
#                 mmDict['bstrainID_0'] = np.append(mmDict['bstrainID_0'],int(bstrain0))
#                 mmDict['bstrainID_1'] = np.append(mmDict['bstrainID_1'],int(bstrain1))
#             if (vstrain1 not in vfreqDict.keys()):
#                 vrel = curSim.execute("SELECT abundance FROM vabundance \
#                     WHERE t = {0} AND vstrain_id = {1}".format(time,vstrain1)).fetchall()
#                 vfreqDict[vstrain1] = vrel[0][0]/vtotal
#             # mvmAbundDict[bstrain1] = np.append(mmDict[bstrain1],tuple((int(bstrain0),vrel[0][0]/vtotal)))
#             tripartiteMVM.loc[tripartiteMVM.shape[0]] = [int(bstrain0),int(bstrain1),int(vstrain1)]
# edge = []
# for bstrain1 in np.unique(tripartiteMVM.bstrainID_1.values):
#     for bstrain0 in np.unique(tripartiteMVM[tripartiteMVM['bstrainID_1'] == bstrain1]
#     .bstrainID_0.values):
#         edge = np.append(
#             edge,sum([vfreqDict[vstrainID] for vstrainID in
#             tripartiteMVM[(tripartiteMVM['bstrainID_1'] == bstrain1)
#             & (tripartiteMVM['bstrainID_0'] == bstrain0)]['vstrainID'].values])
#             )
#
# tripartiteMVM['edge'] = edge


print('Complete!')
