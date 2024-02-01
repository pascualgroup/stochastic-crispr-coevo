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
from mpl_toolkits.mplot3d.art3d import Line3DCollection
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.ticker as ticker
import treefunctions as tree
import generateNetworks as networks
import timestamps as stamps
from scipy.stats import gennorm
from scipy.stats import powerlaw

run = sys.argv[1]

# if  sys.argv[2] ==  'all':
#     times = [float(sys.argv[i]) for i in list(range(5,len(sys.argv)))]
#     print('Compiling abundance plots, trees and networks for {0}'.format(sys.argv[1]))
#     print('...for times: {0}...'.format(', '.join(map(str,times))))
if  sys.argv[2] ==  'networks':
    times = [float(sys.argv[i]) for i in list(range(3,len(sys.argv)))]
    print('Compiling networks for {0}'.format(sys.argv[1]))
    print('...for times: {0}...'.format(', '.join(map(str,times))))
if  sys.argv[2] ==  'trees':
    print('Compiling trees for {0}'.format(sys.argv[1]))
    if len(sys.argv) > 5:
        times = [float(sys.argv[i]) for i in list(range(5,len(sys.argv)))]
        print('...for times: {0}...'.format(', '.join(map(str,times))))
if  sys.argv[2] ==  'abundances':
    print('Compiling abundance plots for {0}'.format(sys.argv[1]))

resolve = 500
imgTypes = ["pdf"]
treeImgTypes = ["pdf"]
# graphImgTypes = ["pdf"]
figxy = (15,10) # setting for tree abundance figure
hratio = [1,3] # setting for tree abundance figure
maxticksize = 100 # setting for abundances on individual branches of tree
treepalette = 'turbo' # Color palette for tree: use contiguous color palette
colorpalSpacers = 'tab20b' # Color palette for spacer in networks: use discrete color palette
babundthreshold  = float(sys.argv[3])/100 # this is %/100 of total population size
vabundthreshold  = float(sys.argv[4])/100 # this is %/100 of total population size
stacked = True
# that a strain has to surpass at its maximum value
# to be included in visualizations
# abundthreshold = 0.05
hyperAnalyze = False
overlay = True
sSpacing = 8
vSpacing = 8
bSpacing = 6

SCRIPT_PATH = os.path.abspath(os.path.dirname(__file__)) # cluster

# SQLITE path to simulation data
DBSIM_PATH = os.path.join(SCRIPT_PATH,'isolates',run,'{}.sqlite'.format(run))

# Call SQLITE simulation database and some related info
conSim = sqlite3.connect(DBSIM_PATH)
curSim = conSim.cursor()
run_id = curSim.execute('SELECT DISTINCT run_id FROM summary').fetchall()
run_id = run_id[0][0]
ID = curSim.execute('SELECT combo_id,replicate FROM runs WHERE run_id = {}'.format(run_id)).fetchall()
combo_id = ID[0][0]
replicate = ID[0][1]

# SQLITE paths, beware of manipulating these. Directories are structured accordingly
DBPROB_PATH = os.path.join(SCRIPT_PATH,'isolates',run,'emergence-lambert-root_output.sqlite') # cluster
DBMATCH_PATH = os.path.join(SCRIPT_PATH,'isolates',run,'matches_output.sqlite') # cluster
DBCLADE_PATH = os.path.join(SCRIPT_PATH,'isolates',run,'clade-abundances_output.sqlite') # cluster
DBTREE_PATH = os.path.join(SCRIPT_PATH,'isolates',run,'trees_output.sqlite') # cluster
DBTRI_PATH = os.path.join(SCRIPT_PATH,'isolates',run,'tripartite-networks_output.sqlite') # cluster
PLOT_PATH = os.path.join(SCRIPT_PATH,'isolates',run,'outputs')
dir = 'crispr-sweep-7-2-2022/isolates/runID3297-c66-r47'
run = 'runID3297-c66-r47'
# DBSIM_PATH = os.path.join('isolates', run, '{}.sqlite'.format(run))
# DBPROB_PATH = os.path.join('isolates',run,'emergence-lambert-root_output.sqlite') # cluster
# DBMATCH_PATH = os.path.join('isolates',run,'matches_output.sqlite') # cluster
# DBCLADE_PATH = os.path.join('isolates',run,'clade-abundances_output.sqlite') # cluster
# DBTREE_PATH = os.path.join('isolates',run,'trees_output.sqlite') # cluster
# DBTRI_PATH = os.path.join('isolates',run,'tripartite-networks_output.sqlite') # cluster
# PLOT_PATH = os.path.join('isolates',run,'outputs')
DBPROB_PATH = os.path.join('/Volumes','Yadgah',dir,'emergence-lambert-root_output.sqlite') # local
DBMATCH_PATH = os.path.join('/Volumes','Yadgah',dir,'matches_output.sqlite') # local
DBMATCHc_PATH = os.path.join('/Volumes','Yadgah','crispr-sweep-7-2-2022/isolated-runs/isolates/comboID66','matchesC66.sqlite') # local
DBTREE_PATH = os.path.join('/Volumes','Yadgah',dir,'trees_output.sqlite') # local
DBTRI_PATH = os.path.join('/Volumes','Yadgah',dir,'tripartite-networks_output.sqlite') # local
PLOT_PATH = os.path.join('/Volumes','Yadgah') # local
DBSIM_PATH = os.path.join('/Volumes','Yadgah',dir,'{}.sqlite'.format(run))
conSim = sqlite3.connect(DBSIM_PATH)
curSim = conSim.cursor()
run_id = curSim.execute('SELECT DISTINCT run_id FROM summary').fetchall()
run_id = run_id[0][0]
ID = curSim.execute('SELECT combo_id,replicate FROM runs WHERE run_id = {}'.format(run_id)).fetchall()
combo_id = ID[0][0]
replicate = ID[0][1]

if not os.path.isdir(PLOT_PATH):
    os.makedirs (PLOT_PATH)

# Call SQLITE analysis databases
conProb = sqlite3.connect(DBPROB_PATH)
curProb = conProb.cursor()
conMatch = sqlite3.connect(DBMATCH_PATH)
curMatch = conMatch.cursor()
conTree = sqlite3.connect(DBTREE_PATH)
curTree = conTree.cursor()
conTri = sqlite3.connect(DBTRI_PATH)
curTri = conTri.cursor()

# Call dataframes of total abundances
microbe_total = pd.read_sql_query("SELECT t,microbial_abundance FROM summary WHERE run_id = {}".format(run_id), conSim)\
.rename(columns={"microbial_abundance": "Microbial Abundance"})
virus_total = pd.read_sql_query("SELECT t,viral_abundance FROM summary WHERE run_id = {}".format(run_id), conSim)\
.rename(columns={"viral_abundance": "Viral Abundance"})
virus_total=virus_total[virus_total["Viral Abundance"]>0]
microbe_total = microbe_total[microbe_total.t <= max(virus_total.t)]


if sys.argv[2] == 'networks':
    vstrainmatchIDs, bstrainmatchIDs = networks.generate(run_id,combo_id,replicate,\
    conSim,conMatch,conTri,times,microbe_total,virus_total,\
    colorpalSpacers,imgTypes,PLOT_PATH,resolve)

if sys.argv[2] == 'trees':
    # Generate tree figures and treeID and colorID info
    vAbunds,vKeepTreeStrainsDF, vSpeciesColorDict, _, _, figV, axesV = \
    tree.speciesTreePlot(run_id,'virus', DBSIM_PATH,DBTREE_PATH,\
    treepalette,maxticksize,figxy,hratio,vabundthreshold,stacked,overlay)
    bAbunds,bKeepTreeStrainsDF, bSpeciesColorDict, _, _, figB, axesB = \
    tree.speciesTreePlot(run_id,'microbe', DBSIM_PATH,DBTREE_PATH,\
    treepalette,maxticksize,figxy,hratio,babundthreshold,stacked,overlay)
    if len(sys.argv) > 5:
        stamps.treetimes(run_id,curSim,conTri,axesV,axesB,times,False,\
        vAbunds,bAbunds,vKeepTreeStrainsDF,bKeepTreeStrainsDF)
    for imgType in treeImgTypes:
        if overlay:
            figB.tight_layout()
            figB.savefig(os.path.join(PLOT_PATH,'microbe-abundances-tree-time-slices_virus-overlay.{0}'\
                .format(imgType)),dpi=resolve)
            figV.tight_layout()
            figV.savefig(os.path.join(PLOT_PATH,'virus-abundances-tree-time-slices_microbe-overlay.{0}'\
                .format(imgType)),dpi=resolve)
        else:
            figB.tight_layout()
            figB.savefig(os.path.join(PLOT_PATH,'microbe-abundances-tree-time-slices.{0}'\
                .format(imgType)),dpi=resolve)
            figV.tight_layout()
            figV.savefig(os.path.join(PLOT_PATH,'virus-abundances-tree-time-slices.{0}'\
                .format(imgType)),dpi=resolve)

    bAbunds, bKeepTreeStrainsDF, bSpeciesColorDict, _, _, figB, axesB, figB2, axesB2 = \
    speciesTreeDiv2Plot(run_id,'microbe', DBSIM_PATH,DBTREE_PATH,\
    treepalette,100,figxy,[1,3],babundthreshold,350,750)

if sys.argv[2] == 'tree-diversity':
    # Generate tree figures and treeID and colorID info
    vAbunds,vKeepTreeStrainsDF, vSpeciesColorDict, _, _, figV, axesV = \
    tree.speciesTreeDivPlot(run_id,'virus', DBSIM_PATH,DBTREE_PATH,\
    treepalette,100,figxy,[1,3,1],vabundthreshold,stacked,overlay)
    bAbunds,bKeepTreeStrainsDF, bSpeciesColorDict, _, _, figB, axesB = \
    tree.speciesTreeDivPlot(run_id,'microbe', DBSIM_PATH,DBTREE_PATH,\
    treepalette,100,figxy,[1,3,1],babundthreshold,stacked,overlay)
    bAbunds, bKeepTreeStrainsDF, bSpeciesColorDict, _, _, figB, axesB, figB2, axesB2 = \
    tree.speciesTreeDiv2Plot(run_id,'microbe', DBSIM_PATH,DBTREE_PATH,\
    treepalette,100,figxy,[1,3],babundthreshold,350,750)
    for imgType in treeImgTypes:
        if overlay:
            figB.tight_layout()
            figB.savefig(os.path.join(PLOT_PATH,'microbe-tree-diversity_virus-overlay.{0}'\
                .format(imgType)),dpi=resolve)
            figV.tight_layout()
            figV.savefig(os.path.join(PLOT_PATH,'virus-tree-diversity_microbe-overlay.{0}'\
                .format(imgType)),dpi=resolve)
        else:
            figB.tight_layout()
            figB.savefig(os.path.join(PLOT_PATH,'microbe-tree-diversity.{0}'\
                .format(imgType)),dpi=resolve)
            figV.tight_layout()
            figV.savefig(os.path.join(PLOT_PATH,'virus-tree-diversity.{0}'\
                .format(imgType)),dpi=resolve)

if sys.argv[2] == 'abundances':
    bAbunds, _, bSpeciesColorDict, _, _, _, _  = \
    tree.speciesTreePlot(run_id,'microbe', DBSIM_PATH,DBTREE_PATH,\
    treepalette,maxticksize,figxy,hratio,babundthreshold,stacked,overlay)
    vAbunds, _, vSpeciesColorDict, _, _, _, _ = \
    tree.speciesTreePlot(run_id,'virus', DBSIM_PATH,DBTREE_PATH,\
    treepalette,maxticksize,figxy,hratio,vabundthreshold,stacked,overlay)
    bAbunds = bAbunds[bAbunds.t <= max(vAbunds.t)]
    bAbunds = bAbunds[bAbunds['abundance']>0]
    # This is to set abundances to log scale
    # bAbunds['abundance'] = np.log(np.array(bAbunds['abundance']))
    # vAbunds['abundance'] = np.log(np.array(vAbunds['abundance']))
    microbe_stacked = bAbunds.pivot(index='t',columns='tree_bstrain_id',values='abundance')
    virus_stacked = vAbunds.pivot(index='t',columns='tree_vstrain_id',values='abundance')
    fig, ax = plt.subplots(1)
    microbe_stacked.plot.area(ax = ax,stacked=True,legend=False, linewidth=0,color=bSpeciesColorDict,sort_columns=True)
    microbe_stacked.plot(stacked=True, ax=ax, legend=False, color='black',sort_columns=True,linewidth=.1)
    ax.set_ylabel(ylabel = 'Host Abundance',labelpad=15,fontsize=15)
    ax.set_xlabel(xlabel = 'Time t',fontsize=15)
    ax.ticklabel_format(axis = 'y',style='sci',scilimits=(0,0))
    ax.xaxis.set_minor_locator(ticker.MultipleLocator(25))
    lim = ax.get_ylim()
    ax.set_ylim(0,lim[1])
    fig.savefig(os.path.join(PLOT_PATH,'microbe-abundances.png'),dpi=resolve)
    plt.close(fig)
    fig, ax = plt.subplots(1)
    virus_stacked.plot.area(ax = ax,stacked=True,legend=False, linewidth=0,color=vSpeciesColorDict,sort_columns=True)
    virus_stacked.plot(stacked=True, ax=ax, legend=False, color='black',sort_columns=True,linewidth=.1)
    ax.set_ylabel(ylabel = 'Viral Abundance',labelpad=15,fontsize=15)
    ax.set_xlabel(xlabel = 'Time t',fontsize=15)
    ax.ticklabel_format(axis = 'y',style='sci',scilimits=(0,0))
    ax.xaxis.set_minor_locator(ticker.MultipleLocator(25))
    lim = ax.get_ylim()
    ax.set_ylim(0,lim[1])
    fig.savefig(os.path.join(PLOT_PATH,'virus-abundances.png'),dpi=resolve)
    plt.close(fig)
    fig, ax = plt.subplots(2,sharex=True,figsize=(20,8))
    axes = [ax[0], ax[1]]
    axes[0].set_ylabel(ylabel ='Host\nAbundance',labelpad=20,fontsize=20)
    axes[0].tick_params(axis='y', labelsize=20)
    axes[0].ticklabel_format(axis = 'y',style='sci',scilimits=(0,0))
    axes[0].yaxis.get_offset_text().set_fontsize(20)
    axes[0].xaxis.set_minor_locator(ticker.MultipleLocator(25))
    axes[1].set_ylabel(ylabel ='Viral\nAbundance',labelpad=20,fontsize=20)
    axes[1].tick_params(axis='y', labelsize=20)
    axes[1].ticklabel_format(axis = 'y',style='sci',scilimits=(0,0))
    axes[0].yaxis.get_offset_text().set_fontsize(20)
    axes[1].tick_params(axis='x', labelsize=20)
    axes[1].xaxis.set_minor_locator(ticker.MultipleLocator(25))
    microbe_stacked.plot.area(ax = axes[0],stacked=True,legend=False, linewidth=0,color=bSpeciesColorDict,sort_columns=True)
    # microbe_stacked.plot(stacked=True, ax=axes[0], legend=False, color='white',sort_columns=True,linewidth=.01)
    virus_stacked.plot.area(ax = axes[1],stacked=True,legend=False, linewidth=0,color=vSpeciesColorDict,sort_columns=True)
    # virus_stacked.plot(stacked=True, ax=axes[1], legend=False, color='white',sort_columns=True,linewidth=.01)
    axes[1].set_xlabel(xlabel="Time t", labelpad=20, fontsize=20)
    axes[0].margins(x=0)
    axes[1].margins(x=0)
    fig.tight_layout()
    fig.savefig(os.path.join(PLOT_PATH,'microbe-virus-stacked-abundances.png'),dpi=resolve)


if sys.argv[2] == 'viralheatmap':
    bAbunds, bKeepTreeStrainsDF, bSpeciesColorDict, _, _, _, _ = \
        tree.speciesTreePlot(run_id, 'microbe',   DBSIM_PATH,  DBTREE_PATH,
                            treepalette, maxticksize, figxy, hratio, babundthreshold, stacked, overlay)
    vAbunds, vKeepTreeStrainsDF, vSpeciesColorDict, hlinecVirus, vlinecVirus, _, _ = \
        tree.speciesTreePlot(run_id, 'virus',   DBSIM_PATH,  DBTREE_PATH,
                            treepalette, maxticksize, figxy, hratio, vabundthreshold, stacked, overlay)
    bAbunds = bAbunds[bAbunds.t <= max(vAbunds.t)]
    bAbunds = bAbunds[bAbunds['abundance'] > 0]
    bstrain0vstrain = pd.read_sql_query(
        "SELECT t, bstrain_id, vstrain_id \
    FROM bstrain_to_vstrain_0matches", conMatch)
    bstrain0vstrain = bstrain0vstrain.merge(vKeepTreeStrainsDF,
                                            on=['vstrain_id']).drop(columns=['vstrain_id', 'tree_vstrain_id'])
    microbe_stacked = pd.read_sql_query("SELECT t,bstrain_id,abundance FROM babundance \
                        WHERE run_id = {}".format(run_id), conSim)
    bstrain0vstrain = bstrain0vstrain.merge(
        microbe_stacked, on=['t', 'bstrain_id'])
    bstrain0vstrain = bstrain0vstrain.groupby(['t', 'new_tree_vstrain_id'])\
        .agg(susceptible=('abundance', 'sum')).reset_index()
    extraTimes = list({*np.arange(min(bstrain0vstrain['t']), max(bstrain0vstrain['t'])+1, 1)}
                    - {*np.array(bstrain0vstrain['t'])})
    for t in extraTimes:
        lastTime = max(bstrain0vstrain[bstrain0vstrain.t < t]['t'])
        new = bstrain0vstrain[bstrain0vstrain.t == lastTime].copy()
        new['t'] = np.ones(len(new['t']))*t
        bstrain0vstrain = pd.concat([bstrain0vstrain, new])
    bstrain0vstrain = bstrain0vstrain.sort_values(by=['t', 'new_tree_vstrain_id'])
    bAbunds = bAbunds.merge(bstrain0vstrain[['t']].drop_duplicates(), on=['t'])
    vAbunds = vAbunds.merge(bstrain0vstrain[['t']].drop_duplicates(), on=['t'])
    virus_total = virus_total.merge(
        bstrain0vstrain[['t']].drop_duplicates(), on=['t'])
    # bstrain0vstrain['susceptible'] = np.log(np.array(bstrain0vstrain['susceptible']))
    treeHeat = bstrain0vstrain.pivot_table(
        index='new_tree_vstrain_id', columns='t', values='susceptible').fillna(0)
    # treeHeat = treeHeat.reindex(index=treeHeat.index[::-1])
    microbe_stacked = bAbunds.pivot(
        index='t', columns='tree_bstrain_id', values='abundance')
    # virus_stacked = vAbunds.pivot(index='t',columns='tree_vstrain_id',values='abundance')
    fig, ax = plt.subplots(2, sharex=True, gridspec_kw={'height_ratios': [1, 6]})
    axes = [ax[0], ax[0].twinx(), ax[1]]  # , ax[2], ax[3]]
    microbe_stacked.plot.area(ax=axes[0], stacked=True, legend=False,
                            linewidth=1.5, color=bSpeciesColorDict, sort_columns=True)
    axes[0].set_ylabel(ylabel='Host\nAbundances',
                    labelpad=15, fontsize=7)
    axes[0].set_xlabel(xlabel='Time t', fontsize=7)
    axes[0].ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
    axes[0].xaxis.set_minor_locator(ticker.MultipleLocator(25))
    lim = axes[0].get_ylim()
    axes[0].set_ylim(0, lim[1])
    axes[1].plot(virus_total['t'], virus_total['Viral Abundance'],
                linewidth=0, color='grey')
    axes[1].fill_between(
        virus_total['t'], virus_total['Viral Abundance'], color='grey', alpha=0.6)
    lim = axes[1].get_ylim()
    axes[1].set_ylim(0, lim[1])
    strainLineages = LineCollection(
        np.array(hlinecVirus) - np.array([[0, 0.5], [0, 0.5]]), linestyles='solid', color='white', linewidths=(.3))
    creationLines = LineCollection(
        np.array(vlinecVirus) - np.array([[0, 0.5], [0, 0.5]]), linestyles='solid', colors='white', linewidths=(.3))
    axes[2].add_collection(strainLineages)
    axes[2].add_collection(creationLines)
    vKeepTreeStrainsDF = vKeepTreeStrainsDF.sort_values(by=['new_tree_vstrain_id'])
    strainsOrdered = vKeepTreeStrainsDF['vstrain_id']
    strainsOrdered = list(strainsOrdered[1::])
    axes[2].pcolor(treeHeat, cmap='turbo')
    divider = make_axes_locatable(axes[2])
    maxAbundanceSpecies = np.max(vAbunds.abundance.values)
    markerIncSpecies = 100/maxAbundanceSpecies
    axes[2].scatter(vAbunds['t'], np.array(vAbunds['tree_vstrain_id'])-0.5, lw=1.5,
                    s=vAbunds['abundance'].values*markerIncSpecies,
                    color='white', marker='|', alpha=0.25)
    axes[0].set_xlim(0, max(treeHeat.columns))
    axes[1].set_xlim(0, max(treeHeat.columns))
    axes[2].set_xlim(0, max(treeHeat.columns))
    axes[2].set_yticks(
        np.arange(.5, len(treeHeat.index), 1))
    axes[2].set_yticklabels(strainsOrdered, fontsize=6)
    # axes[2].set_yticks([])
    # axes[2].set_yticklabels([], fontsize=6)
    axes[2].set_ylabel(ylabel='Viral Strains', labelpad=15, fontsize=7)
    axes[2].margins(x=0)
    axes[2].margins(y=0)
    fig.tight_layout()
    fig.savefig(os.path.join(PLOT_PATH,'viral-heat-map.png'),dpi=resolve)
    ########

if sys.argv[2] == 'protospacer-div':
    bAbunds, bKeepTreeStrainsDF, bSpeciesColorDict, _, _, _, _ = \
        tree.speciesTreePlot(run_id, 'microbe',   DBSIM_PATH,  DBTREE_PATH,
                                treepalette, maxticksize, figxy, hratio, babundthreshold,stacked,overlay)
    vAbunds, vKeepTreeStrainsDF, vSpeciesColorDict, hlinecVirus, vlinecVirus, _, _ = \
        tree.speciesTreePlot(run_id, 'virus',   DBSIM_PATH,  DBTREE_PATH,
                                treepalette, maxticksize, figxy, hratio, vabundthreshold, stacked, overlay)
    bAbunds = bAbunds[bAbunds.t <= max(vAbunds.t)]
    bAbunds = bAbunds[bAbunds['abundance'] > 0]
    tripartiteNetwork = pd.read_sql_query("SELECT t, bstrain_id, vstrain_id, time_specific_match_id \
                            FROM bstrain_to_vstrain_matches WHERE match_length = 1", conMatch)
    spacerMatches = pd.read_sql_query("SELECT t, time_specific_match_id, spacer_id \
                            FROM matches_spacers", conMatch)
    tripartiteNetwork = tripartiteNetwork.merge(spacerMatches, on=['t', 'time_specific_match_id'])\
        .drop(columns=['time_specific_match_id'])
    microbe_stacked = pd.read_sql_query("SELECT t,bstrain_id,abundance FROM babundance \
                        WHERE run_id = {}".format(run_id), conSim)
    tripartiteNetwork = tripartiteNetwork.merge(microbe_stacked, on=['t', 'bstrain_id'])\
        .drop_duplicates()
    spacers = tripartiteNetwork[['t', 'spacer_id', 'bstrain_id', 'abundance']]\
        .drop_duplicates().groupby(['t', 'spacer_id'])\
        .agg(abundance=('abundance', 'sum')).reset_index()
    spacers = spacers.groupby(['t'])\
        .agg(total=('abundance', 'sum')).reset_index()\
        .merge(spacers, on=['t'])
    spacers['freq'] = spacers['abundance']/spacers['total']
    spacers['div'] = np.exp(-1*np.array(spacers['freq'])
                            * np.log(np.array(spacers['freq'])))
    spacers = spacers.groupby(['t'])\
        .agg(div=('div', 'prod')).reset_index()
    microbe_stacked = bAbunds.pivot(
        index='t', columns='tree_bstrain_id', values='abundance')
    fig, ax = plt.subplots(1,sharex=True, figsize=(20,5))
    axes = [ax, ax.twinx(), ax.twinx()]
    microbe_stacked.plot.area(ax=axes[0], stacked=True, legend=False,
                        linewidth=0, cmap='Purples', sort_columns=True, alpha = 0.25)
    axes[1].plot(virus_total['t'], virus_total['Viral Abundance'],
                linewidth=0, color='grey')
    axes[1].fill_between(
        virus_total['t'], virus_total['Viral Abundance'], color='grey', alpha=0.6)
    axes[0].set_yticklabels([])
    axes[0].set_yticks([])
    axes[0].margins(x=0)
    axes[0].xaxis.set_minor_locator(ticker.MultipleLocator(25))
    lim = axes[0].get_ylim()
    axes[0].set_ylim(0, lim[1])
    axes[1].set_yticklabels([])
    axes[1].set_yticks([])
    axes[1].margins(x=0)
    axes[1].xaxis.set_minor_locator(ticker.MultipleLocator(25))
    lim = axes[1].get_ylim()
    axes[1].set_ylim(0, lim[1])
    ###
    axes[2].plot(spacers['t'].values,spacers['div'].values,color='darkblue',linewidth=1.5,label = r'$D_p$')
    axes[2].yaxis.tick_left()
    axes[2].yaxis.set_label_position("left")
    axes[2].set_ylabel(ylabel=''.join(['Single Spacer Match\n',
                        r"Shannon Diversity"]), labelpad=10, fontsize=20)
    axes[2].legend(loc='upper right', fontsize=20)
    axes[2].xaxis.set_minor_locator(ticker.MultipleLocator(25))
    lim = axes[2].get_ylim()
    axes[2].set_ylim(0, lim[1])
    axes[2].margins(x=0)
    axes[2].margins(y=0)
    axes[0].tick_params(axis='x', labelsize=20)
    axes[2].tick_params(axis='y', labelsize=20)
    axes[0].set_xlabel(xlabel='Time t', fontsize=20, labelpad=10)
    axes[1].set_xlabel(xlabel='Time t', fontsize=20, labelpad=10)
    axes[2].set_xlabel(xlabel='Time t', fontsize=20, labelpad=10)
    fig.tight_layout()
    fig.savefig(os.path.join(PLOT_PATH, 'protospacer-div.pdf'), dpi=resolve)










if sys.argv[2] == 'escapeheatmap':
    bAbunds, bKeepTreeStrainsDF, bSpeciesColorDict, _, _, _, _ = \
        tree.speciesTreePlot(run_id, 'microbe',   DBSIM_PATH,  DBTREE_PATH,
                                treepalette, maxticksize, figxy, hratio, babundthreshold,stacked,overlay)
    vAbunds, vKeepTreeStrainsDF, vSpeciesColorDict, hlinecVirus, vlinecVirus, _, _ = \
        tree.mcspeciesTreePlot(run_id, 'virus',   DBSIM_PATH,  DBTREE_PATH,
                                treepalette, maxticksize, figxy, hratio, vabundthreshold, stacked, overlay)
    bAbunds = bAbunds[bAbunds.t <= max(vAbunds.t)]
    bAbunds = bAbunds[bAbunds['abundance'] > 0]
    tripartiteNetwork = pd.read_sql_query("SELECT t, bstrain_id, vstrain_id, time_specific_match_id \
                            FROM bstrain_to_vstrain_matches WHERE match_length = 1", conMatch)
    spacerMatches = pd.read_sql_query("SELECT t, time_specific_match_id, spacer_id \
                            FROM matches_spacers", conMatch)
    tripartiteNetwork = tripartiteNetwork.merge(spacerMatches, on=['t', 'time_specific_match_id'])\
        .drop(columns=['time_specific_match_id'])
    microbe_stacked = pd.read_sql_query("SELECT t,bstrain_id,abundance FROM babundance \
                        WHERE run_id = {}".format(run_id), conSim)
    tripartiteNetwork = tripartiteNetwork.merge(vKeepTreeStrainsDF, on=['vstrain_id'])\
        .drop(columns=['vstrain_id', 'tree_vstrain_id'])\
        .merge(microbe_stacked, on=['t', 'bstrain_id'])\
        .drop_duplicates()
    spacers = tripartiteNetwork[['t', 'spacer_id', 'bstrain_id', 'abundance']]\
        .drop_duplicates().groupby(['t', 'spacer_id'])\
        .agg(abundance=('abundance', 'sum')).reset_index()
    spacers = spacers.groupby(['t'])\
        .agg(total=('abundance', 'sum')).reset_index()\
        .merge(spacers, on=['t'])
    spacers['freq'] = spacers['abundance']/spacers['total']
    spacers['div'] = np.exp(-1*np.array(spacers['freq'])
                            * np.log(np.array(spacers['freq'])))
    tripartiteNetwork = tripartiteNetwork.drop(columns=['bstrain_id', 'abundance'])\
        .drop_duplicates()\
        .merge(spacers[['t', 'spacer_id', 'div']],
                on=['t', 'spacer_id'])
    tripartiteNetwork = tripartiteNetwork.groupby(['t', 'new_tree_vstrain_id'])\
        .agg(matches=('div', 'prod')).reset_index()
    spacers = spacers.groupby(['t'])\
        .agg(div=('div', 'prod')).reset_index()
    extraTimes = list(
        np.arange(min(vAbunds['t']), min(tripartiteNetwork['t']), 1))
    extraStrains = list({*np.array(vAbunds['tree_vstrain_id'])}
                        - {*np.array(tripartiteNetwork['new_tree_vstrain_id'])})
    if len(extraStrains) > 0:
        t0 = len(extraStrains)*[extraTimes[0]]
        strainsbegin = pd.DataFrame(
            {'t': t0,
                'new_tree_vstrain_id': extraStrains,
                'matches': len(extraStrains)*[0]})
    else:
        strainsbegin = pd.DataFrame(
            {'t': extraTimes,
                'new_tree_vstrain_id': len(extraTimes)*[tripartiteNetwork['new_tree_vstrain_id'][0]],
                'matches': len(extraTimes)*[0]})
    tripartiteNetwork = pd.concat([strainsbegin, tripartiteNetwork])
    extraTimes = list({*np.arange(min(vAbunds['t']), max(vAbunds['t'])+1, 1)}
                        - {*np.array(tripartiteNetwork['t'])})
    for t in sorted(extraTimes):
        lastTime = max(tripartiteNetwork[tripartiteNetwork.t < t]['t'])
        new = tripartiteNetwork[tripartiteNetwork.t == lastTime].copy()
        new['t'] = np.ones(len(new['t']))*t
        tripartiteNetwork = pd.concat([tripartiteNetwork, new])
    tripartiteNetwork = tripartiteNetwork.sort_values(
        by=['t', 'new_tree_vstrain_id'])
    treeHeat = tripartiteNetwork.pivot_table(
        index='new_tree_vstrain_id', columns='t', values='matches').fillna(0)
    # treeHeat = treeHeat.reindex(index=treeHeat.index[::-1])
    microbe_stacked = bAbunds.pivot(
        index='t', columns='tree_bstrain_id', values='abundance')
    # virus_stacked = vAbunds.pivot(index='t',columns='tree_vstrain_id',values='abundance')
    fig, ax = plt.subplots(2, sharex=True, figsize=(10,10),gridspec_kw={'height_ratios': [1, 3]})
    axes = [ax[0], ax[0].twinx(), ax[0].twinx(), ax[1]]
    microbe_stacked.plot.area(ax=axes[0], stacked=True, legend=False,
                        linewidth=0, cmap='Purples', sort_columns=True, alpha = 0.25)
    axes[1].plot(virus_total['t'], virus_total['Viral Abundance'],
                linewidth=0, color='grey')
    axes[1].fill_between(
        virus_total['t'], virus_total['Viral Abundance'], color='grey', alpha=0.6)
    axes[0].set_yticklabels([])
    axes[0].set_yticks([])
    axes[0].margins(x=0)
    axes[0].xaxis.set_minor_locator(ticker.MultipleLocator(25))
    lim = axes[0].get_ylim()
    axes[0].set_ylim(0, lim[1])
    axes[1].set_yticklabels([])
    axes[1].set_yticks([])
    axes[1].margins(x=0)
    axes[1].xaxis.set_minor_locator(ticker.MultipleLocator(25))
    lim = axes[1].get_ylim()
    axes[1].set_ylim(0, lim[1])
    ###
    axes[2].plot(spacers['t'],spacers['div'],color='darkblue',linewidth=1.5)
    axes[2].yaxis.tick_left()
    axes[2].yaxis.set_label_position("left")
    axes[2].set_ylabel(ylabel=''.join(['Singly-matched protospacer\n',
                        r"Shannon Diversity ($D_p$)"]), labelpad=15, fontsize=10)
    axes[2].xaxis.set_minor_locator(ticker.MultipleLocator(25))
    lim = axes[2].get_ylim()
    axes[2].set_ylim(0, lim[1])
    ###
    strainLineages = LineCollection(
        np.array(hlinecVirus) - np.array([[0, 0.5], [0, 0.5]]), linestyles='solid', color='white', linewidths=(.3))
    creationLines = LineCollection(
        np.array(vlinecVirus) - np.array([[0, 0.5], [0, 0.5]]), linestyles='solid', colors='white', linewidths=(.3))
    axes[3].add_collection(strainLineages)
    axes[3].add_collection(creationLines)
    vKeepTreeStrainsDF = vKeepTreeStrainsDF.sort_values(by=['new_tree_vstrain_id'])
    strainsOrdered = vKeepTreeStrainsDF['vstrain_id']
    strainsOrdered = list(strainsOrdered[1::])
    # axes[3].pcolor(np.round(treeHeat), cmap='turbo')
    axes[3].pcolor(treeHeat, cmap='turbo', snap=True, linewidth=0,rasterized=True)

    divider = make_axes_locatable(axes[3])
    maxAbundanceSpecies = np.max(vAbunds.abundance.values)
    markerIncSpecies = 100/maxAbundanceSpecies
    axes[3].scatter(vAbunds['t'], np.array(vAbunds['tree_vstrain_id'])-0.5, lw=1.5,
                    s=vAbunds['abundance'].values*markerIncSpecies,
                    color='white', marker='|', alpha=0.6)
    axes[0].set_xlim(0, max(treeHeat.columns))
    axes[1].set_xlim(0, max(treeHeat.columns))
    axes[2].set_xlim(0, max(treeHeat.columns))
    axes[3].set_xlim(0, max(treeHeat.columns))
    axes[3].set_yticks(
        np.arange(.5, len(treeHeat.index), 1))
    axes[3].set_yticklabels(strainsOrdered, fontsize=6)
    # axes[2].set_yticks([])
    # axes[2].set_yticklabels([], fontsize=6)
    axes[3].set_ylabel(ylabel='Viral Strain ID', labelpad=15, fontsize=10)
    axes[3].margins(x=0)
    axes[3].margins(y=0)
    axes[3].set_xlabel(xlabel='Time t', fontsize=10)
    fig.tight_layout()
    fig.savefig(os.path.join(PLOT_PATH, 'viral-escape-heat-map.pdf'), dpi=resolve)
    ####
    left = 200
    right = 600
    top = 15
    triNet = tripartiteNetwork[(tripartiteNetwork.matches>0)\
            &(tripartiteNetwork.t>=left)&(tripartiteNetwork.t<=right)\
                &(tripartiteNetwork.new_tree_vstrain_id <= top)]
    top = max(triNet.new_tree_vstrain_id)
    bottom = min(triNet.new_tree_vstrain_id)
    # vline3d, hline3d = tree.truncateTree(left,right,top,bottom,vlinecVirus,hlinecVirus)
    vline3d, hline3d = tree.truncateTree(left,right,top,bottom,vlinecVirus,hlinecVirus)
    zeros = len(np.array(vline3d))*2
    vline3d = np.concatenate((vline3d, np.reshape(
        np.array([np.zeros(zeros)]), (int(zeros/2), 2, 1))), axis=2)
    zeros = len(np.array(hline3d))*2
    # hline3d = np.concatenate((hline3d, np.reshape(
    #     np.array([np.zeros(zeros)]), (int(zeros/2), 2, 1))), axis=2)
    dx = len(triNet.t)*[1]  # Width of each bar
    dy = len(triNet.t)*[1]  # Depth of each bar
    dz = list(triNet.matches)
    xpos = list(np.array(triNet.t)-np.array(dx)/2)
    ypos = list(np.array(triNet.new_tree_vstrain_id)-np.array(dy)/2)  # y coordinates of each bar
    zpos = len(triNet.t)*[0]  # z coordinates of each bar      # Height of each bar
    strainsOrdered3d = vKeepTreeStrainsDF[\
                        vKeepTreeStrainsDF.new_tree_vstrain_id\
                        .isin(np.arange(bottom, top+1, 1))]\
                        .sort_values(by=['new_tree_vstrain_id'])
    strainsOrdered3d = list(strainsOrdered3d['vstrain_id'])
    ###
    cmap = cm.get_cmap('turbo')
    norm = Normalize(vmin=float(0), vmax=float(max(tripartiteNetwork.matches)))
    colors = cmap(norm(dz))
    fig = plt.figure(figsize=(20, 10))
    ax = fig.add_subplot(111, projection='3d')
    ax.bar3d(xpos, ypos, zpos, dx, dy, dz, color=colors,shade=True,alpha=1)
    creationLines = Line3DCollection(
        vline3d + np.array([[0, 0.25, 0], [0, 0, 0]]),
            linestyles='solid', colors='white', linewidths=(1.25))
    # strainLineages = Line3DCollection(
    #     hline3d + np.array([[0, 0, 0], [0, 0, 0]]),
    #     linestyles='solid', colors='white', linewidths=(1.25))
    # ax.add_collection3d(strainLineages,zs=0,zdir='z')
    ax.add_collection3d(creationLines, zs=0, zdir='z')
    # ax.set_yticks(np.arange(.5, len(treeHeat.index), 1))
    ax.set_yticks(np.arange(bottom, top+1, 1))
    ax.set_yticklabels(strainsOrdered3d, fontsize=7)
    ax.w_zaxis.set_pane_color(cmap(0))
    ax.xaxis._axinfo["grid"].update({"linewidth": 0.05})
    ax.yaxis._axinfo["grid"].update({"linewidth": 0.05})
    ax.xaxis.set_minor_locator(ticker.MultipleLocator(25))
    ax.axes.set_xlim3d(left=left, right=right)
    ax.axes.set_ylim3d(bottom=bottom, top=top)
    lim = ax.get_zlim()
    ax.set_zlim3d(bottom=0, top=lim[1])
    ax.set_xlabel(xlabel="Time t", labelpad=15, fontsize=10)
    ax.set_ylabel(ylabel = "Viral Strain ID", labelpad=15, fontsize=10)
    ax.zaxis.set_rotate_label(False)
    ax.set_zlabel(zlabel=''.join(['Single-match\nShannon Diversity\n',
                                    r"for Viral Strain $i$ ($D_i$)"]), labelpad=30, fontsize=10, rotation=0)
    ax.set_box_aspect((1/2+np.sqrt(5)/2, 1, 1))
    # ax.view_init(elev=25, azim=135)
    ax.view_init(elev=20, azim=118)
    ###
    left = 16
    right = 200
    bottom = 28
    top = 39
    triNet = tripartiteNetwork[(tripartiteNetwork.matches>0)\
            &(tripartiteNetwork.t>=left)&(tripartiteNetwork.t<=right)\
                &(tripartiteNetwork.new_tree_vstrain_id <= top)\
        & (tripartiteNetwork.new_tree_vstrain_id >= bottom)]
    top = max(triNet.new_tree_vstrain_id)
    bottom = min(triNet.new_tree_vstrain_id)
    # vline3d, hline3d = tree.truncateTree(left,right,top,bottom,vlinecVirus,hlinecVirus)
    vline3d, hline3d = truncateTree(left,right,top,bottom,vlinecVirus,hlinecVirus)
    zeros = len(np.array(vline3d))*2
    vline3d = np.concatenate((vline3d, np.reshape(
        np.array([np.zeros(zeros)]), (int(zeros/2), 2, 1))), axis=2)
    zeros = len(np.array(hline3d))*2
    hline3d = np.concatenate((hline3d, np.reshape(
        np.array([np.zeros(zeros)]), (int(zeros/2), 2, 1))), axis=2)
    dx = len(triNet.t)*[1]  # Width of each bar
    dy = len(triNet.t)*[1]  # Depth of each bar
    dz = list(triNet.matches)
    xpos = list(np.array(triNet.t)-np.array(dx)/2)
    ypos = list(np.array(triNet.new_tree_vstrain_id)-np.array(dy)/2)  # y coordinates of each bar
    zpos = len(triNet.t)*[0]  # z coordinates of each bar      # Height of each bar
    strainsOrdered3d = vKeepTreeStrainsDF[\
                        vKeepTreeStrainsDF.new_tree_vstrain_id\
                        .isin(np.arange(bottom, top+1, 1))]\
                        .sort_values(by=['new_tree_vstrain_id'])
    strainsOrdered3d = list(strainsOrdered3d['vstrain_id'])
    colors = cmap(norm(dz))
    ###
    fig = plt.figure(figsize=(20, 10))
    ax = fig.add_subplot(111, projection='3d')
    ax.bar3d(xpos, ypos, zpos, dx, dy, dz, color=colors,shade=True,alpha=1)
    creationLines = Line3DCollection(
        vline3d + np.array([[0, 0.25, 0], [0, 0, 0]]),
            linestyles='solid', colors='white', linewidths=(1.25))  
    ax.add_collection3d(creationLines, zs=0, zdir='z')
    # ax.set_yticks(np.arange(.5, len(treeHeat.index), 1))
    ax.set_yticks(np.arange(bottom, top+1, 1))
    ax.set_yticklabels(strainsOrdered3d, fontsize=7)
    ax.w_zaxis.set_pane_color(cmap(0))
    ax.xaxis._axinfo["grid"].update({"linewidth": 0.05})
    ax.yaxis._axinfo["grid"].update({"linewidth": 0.05})
    ax.xaxis.set_minor_locator(ticker.MultipleLocator(25))
    ax.axes.set_xlim3d(left=left, right=right)
    ax.axes.set_ylim3d(bottom=bottom, top=top)
    lim = ax.get_zlim()
    ax.set_zlim3d(bottom=0, top=lim[1])
    ax.set_xlabel(xlabel="Time t", labelpad=15, fontsize=10)
    ax.set_ylabel(ylabel = "Viral Strain ID", labelpad=15, fontsize=10)
    ax.zaxis.set_rotate_label(False)
    ax.set_zlabel(zlabel= "Shannon Diversity of\nSingly-matched\nProtospacers", labelpad=30, fontsize=10,rotation=0)
    ax.set_box_aspect((1/2+np.sqrt(5)/2, 1, 1))
    ax.view_init(elev=20, azim=118)
    plt.close('all')
    sc = cm.ScalarMappable(cmap=cmap, norm=norm)
    sc.set_array([])
    # fig.subplots_adjust(right=0.8)
    # put colorbar at desire position
    fig = plt.figure(figsize=(10,10))
    cbar_ax = fig.add_axes([0.85, 0.1, 0.02, 0.8])
    fig.colorbar(sc,cax=cbar_ax,ax=ax)


if sys.argv[2] == 'diversity':
    conMatchC = sqlite3.connect(DBMATCHc_PATH)
    curMatchC = conMatchC.cursor()
    bAbunds, bKeepTreeStrainsDF, bSpeciesColorDict, _, _, _, _ = \
    tree.speciesTreePlot(run_id, 'microbe', DBSIM_PATH, DBTREE_PATH,
    treepalette, maxticksize, figxy, hratio,babundthreshold,stacked,overlay)
    vAbunds, vKeepTreeStrainsDF, vSpeciesColorDict, hlinecVirus, vlinecVirus, _, _ = \
    tree.speciesTreePlot(run_id,'virus', DBSIM_PATH,DBTREE_PATH,\
    treepalette,maxticksize,figxy,hratio,vabundthreshold,stacked,overlay)
    microbe_stacked = pd.read_sql_query("SELECT t,bstrain_id,abundance FROM babundance \
                        WHERE run_id = {}".format(run_id), conSim)
    bmatchspacers = pd.read_sql_query("SELECT t, bstrain_id, time_specific_match_id \
                                FROM bstrain_to_vstrain_matches WHERE match_length = 1", conMatch)
    spacerMatches = pd.read_sql_query("SELECT t, time_specific_match_id, spacer_id\
                                FROM matches_spacers", conMatch)
    bmatchspacers = bmatchspacers.merge(spacerMatches, on=['t','time_specific_match_id'])\
        .drop(columns=['time_specific_match_id']).drop_duplicates()
    microbe_stacked = microbe_stacked[microbe_stacked.t <= max(vAbunds.t)]
    bAbunds = bAbunds[bAbunds.t <= max(vAbunds.t)]
    bfreq = microbe_stacked.groupby(['t','bstrain_id']).agg(bfreq=('abundance','sum'))\
                .reset_index()
    bfreq = bfreq.groupby(['t']).agg(btotal=('bfreq','sum')).reset_index()\
            .merge(bfreq,on=['t'])
    bfreq['bfreq'] = bfreq['bfreq']/bfreq['btotal']
    bfreq = bfreq.drop(columns=['btotal'])
    bfreq= bfreq[bfreq.bfreq!=0]
    straintime = sorted(bfreq['t'].unique())
    spacerArrays = pd.read_sql_query(
    "SELECT bstrain_id, num_total_spacers FROM bstrain_num_total_spacers WHERE run_id = {0}".format(run_id), conMatchC)
    spacerArrays = spacerArrays.merge(bfreq, on=['bstrain_id'])
    spacerArrays['num_total_spacers'] = spacerArrays['num_total_spacers'] * spacerArrays['bfreq']
    spacerArrays = spacerArrays.groupby(['t']).agg(
        length=('num_total_spacers', 'sum')).reset_index()
    strainCreation = pd.read_sql_query(
    "SELECT t_creation, bstrain_id FROM bstrains WHERE run_id = {0}".format(run_id), conSim)
    strainCreation = strainCreation.merge(bfreq, on=['bstrain_id'])
    strainCreation['t_creation'] = strainCreation['t_creation'] * \
        strainCreation['bfreq']
    strainCreation = strainCreation.groupby(['t']).agg(
        t_creation=('t_creation', 'sum')).reset_index()
    shannonStrain = []
    for t in sorted(bfreq['t'].unique()):
        div = list(bfreq[(bfreq['t'] == t)]['bfreq'])
        shannonStrain.append(np.exp(-1*sum(np.array(div)*np.array(list(map(np.log,div))))))
    richnessStrain = []
    for t in sorted(bfreq['t'].unique()):
        div = list(bfreq[(bfreq['t'] == t) & (
            bfreq['bfreq'] != 0)]['bfreq'])
        richnessStrain.append(len(div))
    bfreq = bmatchspacers\
    .merge(microbe_stacked,on=['t','bstrain_id'])\
        .groupby(['t','spacer_id']).agg(bfreq=('abundance','sum'))\
    .reset_index()
    bfreq = bfreq.groupby(['t']).agg(btotal=('bfreq','sum')).reset_index()\
            .merge(bfreq,on=['t'])
    bfreq['bfreq'] = bfreq['bfreq']/bfreq['btotal']
    bfreq = bfreq.drop(columns=['btotal'])
    bfreq = bfreq[bfreq.bfreq != 0]
    spacertime = sorted(bfreq['t'].unique())
    shannonSpacer = []
    for t in sorted(bfreq['t'].unique()):
        div = list(bfreq[(bfreq['t']==t)]['bfreq'])
        shannonSpacer.append(np.exp(-1*sum(np.array(div)*np.array(list(map(np.log,div))))))
    richnessSpacer = []
    for t in sorted(bfreq['t'].unique()):
        div = list(bfreq[(bfreq['t']==t)]['bfreq'])
        richnessSpacer.append(len(div))
    spacerCreation = pd.read_sql_query(
    "SELECT t_creation, bstrain_id FROM bstrains WHERE run_id = {0}".format(run_id), conSim)
    spacerCreation = spacerCreation.merge(pd.read_sql_query(
        "SELECT bstrain_id, spacer_id FROM bspacers WHERE run_id = {0}".format(run_id), conSim))
    spacerCreation = spacerCreation.groupby(['spacer_id']).agg(
        t_creation=('t_creation', 'min')).reset_index()\
        .merge(bfreq,on=['spacer_id']).drop_duplicates()
    spacerCreation['t_creation'] = spacerCreation['t_creation'] * spacerCreation['bfreq']
    spacerCreation = spacerCreation.groupby(['t']).agg(
        t_creation=('t_creation', 'sum')).reset_index()
    fig, ax = plt.subplots(7,sharex=True)
    axes = [ax[0], ax[0].twinx(), ax[0].twinx(), ax[4], ax[4].twinx(), 
            ax[4].twinx(), ax[3], ax[3].twinx(), ax[1], ax[1].twinx(), 
            ax[2], ax[2].twinx(),
            ax[5], ax[5].twinx(),
            ax[6], ax[6].twinx()]
    microbe_stacked = bAbunds.pivot(index='t',columns='tree_bstrain_id',values='abundance')
    microbe_stacked.plot.area(ax = axes[0],stacked=True,legend=False, linewidth=0,color=bSpeciesColorDict,sort_columns=True,alpha=0.25)
    microbe_stacked.plot(stacked=True, ax=axes[0], legend=False, color='white',sort_columns=True,linewidth=.1)
    axes[0].set_yticklabels([])
    axes[0].set_yticks([])
    axes[0].margins(x=0)
    microbe_stacked.plot.area(ax = axes[3],stacked=True,legend=False, linewidth=0,color=bSpeciesColorDict,sort_columns=True,alpha=0.25)
    microbe_stacked.plot(stacked=True, ax=axes[3], legend=False, color='white',sort_columns=True,linewidth=.1)
    axes[3].set_yticklabels([])
    axes[3].set_yticks([])
    axes[3].margins(x=0)
    axes[0].set_xlabel(xlabel = 'Time t',fontsize=7)
    axes[0].xaxis.set_minor_locator(ticker.MultipleLocator(25))
    axes[1].plot(straintime,richnessStrain,color='darkred',linewidth=1)
    axes[1].yaxis.tick_right()
    axes[1].set_ylabel(ylabel ='Richness',labelpad=15,fontsize=7,rotation=270)
    axes[2].plot(straintime,shannonStrain,color='darkblue',linewidth=1)
    axes[2].yaxis.tick_left()
    axes[2].yaxis.set_label_position("left")
    axes[2].set_ylabel(ylabel ='Host Strain\nShannon Diversity',labelpad=15,fontsize=7)
    axes[4].plot(spacertime,richnessSpacer,color='darkred',linewidth=1)
    axes[4].yaxis.tick_right()
    axes[4].set_ylabel(ylabel ='Richness',labelpad=15,fontsize=7,rotation=270)
    axes[5].plot(spacertime,shannonSpacer,color='darkblue',linewidth=1)
    axes[5].yaxis.tick_left()
    axes[5].yaxis.set_label_position("left")
    axes[5].set_ylabel(ylabel ='Single Spacer\nShannon Diversity',labelpad=15,fontsize=7)
    microbe_stacked.plot.area(ax = axes[6],stacked=True,legend=False, linewidth=0,color=bSpeciesColorDict,sort_columns=True,alpha=0.25)
    microbe_stacked.plot(stacked=True, ax=axes[6], legend=False, color='white',sort_columns=True,linewidth=.1)
    axes[6].set_yticklabels([])
    axes[6].set_yticks([])
    axes[6].margins(x=0)
    axes[6].set_xlabel(xlabel = 'Time t',fontsize=7)
    axes[7].plot(spacerArrays['t'],spacerArrays['length'],color='darkgreen',linewidth=1)
    axes[7].yaxis.tick_left()
    axes[7].yaxis.set_label_position("left")
    axes[7].set_ylabel(ylabel ='Expected Host\nSpacer Array Length',labelpad=15,fontsize=7)
    microbe_stacked.plot.area(ax = axes[8],stacked=True,legend=False, linewidth=0,color=bSpeciesColorDict,sort_columns=True,alpha=0.25)
    microbe_stacked.plot(stacked=True, ax=axes[8], legend=False, color='white',sort_columns=True,linewidth=.1)
    axes[8].set_yticklabels([])
    axes[8].set_yticks([])
    axes[8].margins(x=0)
    axes[8].set_xlabel(xlabel = 'Time t',fontsize=7)
    axes[8].xaxis.set_minor_locator(ticker.MultipleLocator(25))
    axes[9].plot(straintime,np.array(shannonStrain)/np.array(richnessStrain),color='purple',linewidth=1)
    axes[9].yaxis.tick_left()
    axes[9].yaxis.set_label_position("left")
    axes[9].set_ylabel(ylabel ='Host Strain\nShannon/Richness',labelpad=15,fontsize=7)
    microbe_stacked.plot.area(ax = axes[10],stacked=True,legend=False, linewidth=0,color=bSpeciesColorDict,sort_columns=True,alpha=0.25)
    microbe_stacked.plot(stacked=True, ax=axes[10], legend=False, color='white',sort_columns=True,linewidth=.1)
    axes[10].set_yticklabels([])
    axes[10].set_yticks([])
    axes[10].margins(x=0)
    axes[10].set_xlabel(xlabel = 'Time t',fontsize=7)
    axes[10].xaxis.set_minor_locator(ticker.MultipleLocator(25))
    axes[11].plot(strainCreation['t'],strainCreation['t_creation'],color='darkorange',linewidth=1)
    axes[11].yaxis.tick_left()
    axes[11].yaxis.set_label_position("left")
    axes[11].set_ylabel(ylabel ='Expected Host Strain\nCreation Time',labelpad=15,fontsize=7)
    microbe_stacked.plot.area(ax = axes[12],stacked=True,legend=False, linewidth=0,color=bSpeciesColorDict,sort_columns=True,alpha=0.25)
    microbe_stacked.plot(stacked=True, ax=axes[12], legend=False, color='white',sort_columns=True,linewidth=.1)
    axes[12].set_yticklabels([])
    axes[12].set_yticks([])
    axes[12].margins(x=0)
    axes[12].set_xlabel(xlabel = 'Time t',fontsize=7)
    axes[12].xaxis.set_minor_locator(ticker.MultipleLocator(25))
    axes[13].plot(spacertime, np.array(shannonSpacer)/np.array(richnessSpacer), color='purple', linewidth=1)
    axes[13].yaxis.tick_left()
    axes[13].yaxis.set_label_position("left")
    axes[13].set_ylabel(ylabel='Single Spacer\nShannon/Richness', labelpad=15, fontsize=7)
    microbe_stacked.plot.area(ax = axes[14],stacked=True,legend=False, linewidth=0,color=bSpeciesColorDict,sort_columns=True,alpha=0.25)
    microbe_stacked.plot(stacked=True, ax=axes[14], legend=False, color='white',sort_columns=True,linewidth=.1)
    axes[14].set_yticklabels([])
    axes[14].set_yticks([])
    axes[14].margins(x=0)
    axes[14].set_xlabel(xlabel = 'Time t',fontsize=7)
    axes[14].xaxis.set_minor_locator(ticker.MultipleLocator(25))
    axes[15].plot(spacerCreation['t'],spacerCreation['t_creation'],color='darkorange',linewidth=1)
    axes[15].yaxis.tick_left()
    axes[15].yaxis.set_label_position("left")
    axes[15].set_ylabel(ylabel='Expected Spacer\nAcquisition Time', labelpad=15, fontsize=7)
    fig.tight_layout()
    fig.savefig(os.path.join(PLOT_PATH,'diversity.pdf'),dpi=resolve)

if sys.argv[2] == 'creationtime':
    bAbunds, bKeepTreeStrainsDF, bSpeciesColorDict, _, _, _, _ = \
    tree.speciesTreePlot(run_id, 'microbe', DBSIM_PATH, DBTREE_PATH,
    treepalette, maxticksize, figxy, hratio,babundthreshold,stacked,overlay)
    vAbunds, vKeepTreeStrainsDF, vSpeciesColorDict, hlinecVirus, vlinecVirus, _, _ = \
    tree.speciesTreePlot(run_id,'virus', DBSIM_PATH,DBTREE_PATH,\
    treepalette,maxticksize,figxy,hratio,vabundthreshold,stacked,overlay)
    microbe_stacked = pd.read_sql_query("SELECT t,bstrain_id,abundance FROM babundance \
                        WHERE run_id = {}".format(run_id), conSim)
    microbe_stacked = microbe_stacked[microbe_stacked.t <= max(vAbunds.t)]
    bAbunds = bAbunds[bAbunds.t <= max(vAbunds.t)]
    bfreq = microbe_stacked.groupby(['t','bstrain_id']).agg(bfreq=('abundance','sum'))\
                .reset_index()
    bfreq = bfreq.groupby(['t']).agg(btotal=('bfreq','sum')).reset_index()\
            .merge(bfreq,on=['t'])
    bfreq['bfreq'] = bfreq['bfreq']/bfreq['btotal']
    bfreq = bfreq.drop(columns=['btotal'])
    bfreq= bfreq[bfreq.bfreq!=0]
    creationTimes = pd.read_sql_query("SELECT t_creation,bstrain_id FROM bstrains \
                        WHERE run_id = {}".format(run_id), conSim)\
                            .merge(bfreq,on=['bstrain_id'])
    creationTimes['exp'] = creationTimes['bfreq']*creationTimes['t_creation']
    creationTimes = creationTimes.groupby(['t']).agg(exp=('exp', 'sum'))\
                .reset_index()
    fig, ax = plt.subplots(1, sharex=True)
    ax.plot(creationTimes['t'],creationTimes['exp'])
   

    fig, ax = plt.subplots(2,sharex=True)
    axes = [ax[0], ax[0].twinx(), ax[0].twinx()]
    microbe_stacked = bAbunds.pivot(index='t',columns='tree_bstrain_id',values='abundance')
    microbe_stacked.plot.area(ax = axes[0],stacked=True,legend=False, linewidth=0,color=bSpeciesColorDict,sort_columns=True,alpha=0.25)
    microbe_stacked.plot(stacked=True, ax=axes[0], legend=False, color='white',sort_columns=True,linewidth=.1)
    axes[0].set_yticklabels([])
    axes[0].set_yticks([])
    axes[0].margins(x=0)
    microbe_stacked.plot.area(ax = axes[3],stacked=True,legend=False, linewidth=0,color=bSpeciesColorDict,sort_columns=True,alpha=0.25)
    microbe_stacked.plot(stacked=True, ax=axes[3], legend=False, color='white',sort_columns=True,linewidth=.1)
    axes[0].set_xlabel(xlabel = 'Time t',fontsize=7)
    axes[0].xaxis.set_minor_locator(ticker.MultipleLocator(25))
    axes[1].plot(straintime,richnessStrain,color='darkred',linewidth=1)
    axes[1].yaxis.tick_right()
    axes[1].set_ylabel(ylabel ='Richness',labelpad=15,fontsize=7)
    axes[2].plot(straintime,shannonStrain,color='darkblue',linewidth=1)
    axes[2].yaxis.tick_left()
    axes[2].yaxis.set_label_position("left")
    axes[2].set_ylabel(ylabel ='Host Immune Strain\nShannon Diversity',labelpad=15,fontsize=7)
    fig.savefig(os.path.join(PLOT_PATH,'host-strain-creation-time.png'),dpi=resolve)
    plt.close('all')
    fig, ax = plt.subplots(2,sharex=True)
    axes = [ax[0], ax[0].twinx(), ax[0].twinx()]
    microbe_stacked.plot.area(ax = axes[0],stacked=True,legend=False, linewidth=0,color=bSpeciesColorDict,sort_columns=True,alpha=0.25)
    microbe_stacked.plot(stacked=True, ax=axes[0], legend=False, color='white',sort_columns=True,linewidth=.1)
    axes[0].set_yticklabels([])
    axes[0].set_yticks([])
    axes[0].margins(x=0)
    microbe_stacked.plot.area(ax = axes[3],stacked=True,legend=False, linewidth=0,color=bSpeciesColorDict,sort_columns=True,alpha=0.25)
    microbe_stacked.plot(stacked=True, ax=axes[3], legend=False, color='white',sort_columns=True,linewidth=.1)
    axes[0].set_xlabel(xlabel = 'Time t',fontsize=7)
    axes[0].xaxis.set_minor_locator(ticker.MultipleLocator(25))
    axes[1].plot(spacertime,richnessSpacer,color='darkred',linewidth=1)
    axes[1].yaxis.tick_right()
    axes[1].set_ylabel(ylabel ='Richness',labelpad=15,fontsize=7)
    axes[2].plot(spacertime,shannonSpacer,color='darkblue',linewidth=1)
    axes[2].yaxis.tick_left()
    axes[2].yaxis.set_label_position("left")
    axes[2].set_ylabel(ylabel ='Single-match Spacer\nShannon Diversity',labelpad=15,fontsize=7)
    fig.tight_layout()
    fig.savefig(os.path.join(PLOT_PATH,'spacer-creation-time.png'),dpi=resolve)




if sys.argv[2] == 'emergence':
    bAbunds, bKeepTreeStrainsDF, bSpeciesColorDict, _, _, _, _  = \
    tree.speciesTreePlot(run_id,'microbe',DBSIM_PATH,DBTREE_PATH,\
    treepalette,maxticksize,figxy,hratio,babundthreshold,stacked,overlay)
    vAbunds, vKeepTreeStrainsDF, vSpeciesColorDict, _, _, _, _ = \
    tree.speciesTreePlot(run_id,'virus', DBSIM_PATH,DBTREE_PATH,\
    treepalette,maxticksize,figxy,hratio,vabundthreshold,stacked,overlay)
    p = curSim.execute('SELECT viral_mutation_rate,viral_burst_size,microbe_carrying_capacity, \
                        spacer_acquisition_prob, adsorption_rate, viral_decay_rate \
                        FROM param_combos WHERE combo_id = {}'.format(combo_id)).fetchall()
    mu = p[0][0]
    beta = p[0][1]
    K = p[0][2]
    q = p[0][3]
    phi = p[0][4]
    d = p[0][5]
    microbe_total = pd.read_sql_query("SELECT t,microbial_abundance FROM summary WHERE run_id = {}".format(run_id), conSim)\
    .rename(columns={"microbial_abundance": "btotal"})
    virus_total = pd.read_sql_query("SELECT t,viral_abundance FROM summary WHERE run_id = {}".format(run_id), conSim)\
    .rename(columns={"viral_abundance": "vtotal"})
    virus_total=virus_total[virus_total["vtotal"]>0]
    microbe_stacked = pd.read_sql_query("SELECT t,bstrain_id,abundance FROM babundance \
                        WHERE run_id = {}".format(run_id), conSim)
    virus_stacked = pd.read_sql_query("SELECT t,vstrain_id,abundance FROM vabundance \
                        WHERE run_id = {}".format(run_id), conSim)
    bstrain0vstrains = pd.read_sql_query("SELECT t, bstrain_id, vstrain_id\
        FROM bstrain_to_vstrain_0matches", conMatch)
    vmatches = pd.read_sql_query("SELECT t, match_id, vstrain_id\
        FROM vmatches", conTri).rename(columns={'match_id':'vmatch_id'})
    v0matches = pd.read_sql_query("SELECT t, vstrain_id \
        FROM v0matches", conTri)
    v0matches['vmatch_id'] = len(v0matches['vstrain_id'])*[0]
    vmatches = pd.concat([v0matches,vmatches])\
                .sort_values(by=['t','vmatch_id','vstrain_id'])
    tripartiteNets = pd.read_sql_query("SELECT t, bstrain_id, vstrain_id, time_specific_match_id \
                FROM bstrain_to_vstrain_matches WHERE match_length = 1", conMatch)#\
        # .merge(bstrain0vstrains[['t', 'vstrain_id']], on=['t', 'vstrain_id'])
    vAbunds = vAbunds.rename(columns={'tree_vstrain_id':'new_tree_vstrain_id'})\
        .merge(vKeepTreeStrainsDF,on=['new_tree_vstrain_id'])[['t','abundance','vstrain_id']]
    bstrain0vstrains = bstrain0vstrains.merge(vmatches,on=['t','vstrain_id'])\
                        .merge(microbe_stacked,on=['t','bstrain_id'])\
                        .groupby(['t','vmatch_id'])\
                        .agg(babundance=('abundance','sum')).reset_index()\
                        .merge(microbe_total,on=['t'])
    #####
    ####
    bstrain0vstrains['bfreq'] = \
        np.array(bstrain0vstrains['babundance'])\
        /np.array(bstrain0vstrains['btotal'])
    bstrain0vstrains = bstrain0vstrains.drop(columns=['babundance','btotal'])
    vmatches = vmatches.merge(virus_stacked,on=['t','vstrain_id'])\
            .groupby(['t','vmatch_id'])\
            .agg(vabundance=('abundance','sum')).reset_index()\
            .merge(virus_total,on=['t'])
    vmatches['vfreq'] = np.array(vmatches['vabundance'])/np.array(vmatches['vtotal'])
    bstrain0vstrains = bstrain0vstrains\
                        .merge(vmatches[['t','vmatch_id','vfreq','vabundance']],\
                            on=['t','vmatch_id'])
    ######
    bstrain0vstrains['weights'] = \
        np.array(bstrain0vstrains['vfreq'])*np.array(bstrain0vstrains['bfreq'])
    # bstrain0vstrains = bstrain0vstrains[bstrain0vstrains.weights>0]
    # bstrain0vstrains['weights'] = \
    #     1/np.array(bstrain0vstrains['vfreq'])
    # bstrain0vstrains['weights'] = \
    #     np.array(bstrain0vstrains['vfreq'])
    #####
    # bstrain0vstrains = bstrain0vstrains[bstrain0vstrains.vabundance <= 100]
    #
    norm = bstrain0vstrains.groupby(['t']).agg(norm=('weights','sum')).reset_index()
    bstrain0vstrains = bstrain0vstrains.merge(norm,on=['t'])
    bstrain0vstrains['weights'] = np.array(bstrain0vstrains['weights'])\
                                            /np.array(bstrain0vstrains['norm'])
    bstrain0vstrains = pd.read_sql_query("SELECT t, vmatch_id, p_extinction_lambert \
        FROM existing_vmatch_extinction ORDER BY t", conProb)\
        .merge(bstrain0vstrains,on=['t','vmatch_id'])
    ##
    fig, ax = plt.subplots(1)
    vmatches = pd.read_sql_query("SELECT t, match_id, vstrain_id\
        FROM vmatches", conTri).rename(columns={'match_id': 'vmatch_id'})
    bstrain0escapes = bstrain0vstrains[['t', 'vmatch_id', 'p_extinction_lambert']]\
        .merge(vmatches,on=['t','vmatch_id'])\
            .drop(columns=['vmatch_id']).merge(vAbunds[['t','vstrain_id']],on=['t','vstrain_id'])
    bstrain0escapes['pemerge'] = 1 - bstrain0escapes['p_extinction_lambert']
    bstrain0escapes.pivot(
        index='t', columns='vstrain_id', values='pemerge').plot(ax=ax, stacked=False, legend=False,
                        linewidth=1.5, cmap='turbo', sort_columns=True)
    ##
    bstrain0vstrains['pExp'] = \
        (1-(np.array(bstrain0vstrains['p_extinction_lambert']))**1)\
        * np.array(bstrain0vstrains['weights'])
    pEmergeExpected = bstrain0vstrains.groupby(['t'])\
        .agg(p_exp=('pExp', 'sum')).reset_index()
    #
    escapeProbs = pd.read_sql_query("SELECT t, escape_vmatch_id, p_extinction_lambert \
        FROM single_escapes_vmatch_extinction ORDER BY t", conProb)
    spacerMatches = pd.read_sql_query("SELECT t, time_specific_match_id, spacer_id \
            FROM matches_spacers", conMatch)
    vmatches = pd.read_sql_query("SELECT t, match_id, vstrain_id\
            FROM vmatches", conTri).rename(columns={'match_id': 'vmatch_id'})
    tripartiteNets = tripartiteNets.merge(spacerMatches, on=['t', 'time_specific_match_id'])\
            .drop(columns=['time_specific_match_id','bstrain_id','spacer_id'])\
                .drop_duplicates().merge(vmatches,on=['t','vstrain_id'])\
                    .drop(columns=['vstrain_id']).drop_duplicates()
    escapeProbs = pd.read_sql_query("SELECT match_id, escape_match_id \
        FROM potential_vmatch_single_escapes", conProb)\
            .rename(columns={'match_id':'vmatch_id','escape_match_id':'escape_vmatch_id'})\
                .merge(escapeProbs,on=['escape_vmatch_id'])\
                .merge(tripartiteNets,on=['t','vmatch_id'])
    escapeProbs = escapeProbs.groupby(['t', 'vmatch_id'])\
        .agg(length=('p_extinction_lambert', 'size')).reset_index()\
            .merge(escapeProbs,on=['t','vmatch_id'])
    tt = sorted(np.unique(escapeProbs.t))
    plt.close('all')
    for t in tt:
        escapes = set(escapeProbs[escapeProbs.t==t]['escape_vmatch_id'])
        matchIDs = set(escapeProbs[escapeProbs.t == t]['vmatch_id'])
        overlap = list(matchIDs.intersection(escapes))
        if len(overlap) > 0:
            # print('t:{0},overlap:{1}'.format(t,overlap))
            for o in overlap:
                escapeProbs = escapeProbs[~(np.array(escapeProbs.t==t)*np.array(escapeProbs.escape_vmatch_id==o))]
    ###
    min1 = escapeProbs.groupby(['t', 'vmatch_id'])\
        .agg(min=('p_extinction_lambert', 'min')).reset_index()\
        .merge(escapeProbs,on=['t', 'vmatch_id'])
    min1rows = list(np.array(min1['p_extinction_lambert']) ==
                    np.array(min1['min']))  # keep 1st rank
    min1 = min1[min1rows].drop_duplicates()  # keep 1st rank
    #
    min2 = escapeProbs[['t', 'vmatch_id','p_extinction_lambert']]\
        .merge(min1[['t', 'vmatch_id', 'min']].drop_duplicates(),
            on=['t', 'vmatch_id'])
    min1rows = list(np.array(min2['p_extinction_lambert']) !=
                    np.array(min2['min']))  # remove 1st rank
    min2 = min2[min1rows].drop(columns=['min'])  # remove 1st rank
    min2 = min2.groupby(['t', 'vmatch_id'])\
        .agg(min=('p_extinction_lambert', 'min')).reset_index()\
        .merge(min2, on=['t', 'vmatch_id'])
    min2rows = list(np.array(min2['p_extinction_lambert']) ==
                    np.array(min2['min']))  # keep 2nd rank
    min2 = min2[min2rows].drop_duplicates()  # keep 2nd rank
    #
    min3 = escapeProbs[['t', 'vmatch_id', 'p_extinction_lambert']]\
        .merge(min1[['t', 'vmatch_id', 'min']],
            on=['t', 'vmatch_id'])
    min1rows = list(np.array(min3['p_extinction_lambert']) !=
                    np.array(min3['min']))  # remove 1st rank
    min3 = min3[min1rows].drop(columns=['min'])  # remove 1st rank
    min3 = min3.merge(min2[['t', 'vmatch_id', 'min']], on=['t', 'vmatch_id'])
    min2rows = list(np.array(min3['p_extinction_lambert']) !=
                    np.array(min3['min']))  # remove 2nd rank
    min3 = min3[min2rows].drop(columns=['min'])  # remove 2nd rank
    min3 = min3.groupby(['t', 'vmatch_id'])\
        .agg(min=('p_extinction_lambert', 'min')).reset_index()\
        .merge(min3, on=['t', 'vmatch_id'])
    min3rows = list(np.array(min3['p_extinction_lambert']) ==
                    np.array(min3['min']))  # keep 3rd rank
    min3 = min3[min3rows].drop_duplicates()  # keep 3rd rank  
    # at least...
    minAtLeast2 = pd.concat([min1, min2]).sort_values(by=['t', 'vmatch_id'])
    minAtLeast3 = pd.concat([minAtLeast2, min3]).sort_values(
        by=['t', 'vmatch_id'])
    min1 = bstrain0vstrains[['t', 'vmatch_id', 'weights']]\
        .drop_duplicates().merge(min1, on=['t', 'vmatch_id'])
    min2 = bstrain0vstrains[['t', 'vmatch_id', 'weights']]\
        .drop_duplicates().merge(min2, on=['t', 'vmatch_id'])
    min3 = bstrain0vstrains[['t', 'vmatch_id', 'weights']]\
        .drop_duplicates().merge(min3, on=['t', 'vmatch_id'])
    minAtLeast2 = bstrain0vstrains[['t', 'vmatch_id', 'weights']]\
        .drop_duplicates().merge(minAtLeast2, on=['t', 'vmatch_id'])
    minAtLeast3 = bstrain0vstrains[['t', 'vmatch_id', 'weights']]\
        .drop_duplicates().merge(minAtLeast3, on=['t', 'vmatch_id'])
    ###
    norm = min1.groupby(['t']).agg(
        norm=('weights', 'sum')).reset_index()
    min1 = min1.merge(norm, on=['t'])
    min1['weights'] = ((np.array(min1['weights'])
                                    / np.array(min1['norm'])))
    min1['expmin'] = np.array(min1['min']) *\
        np.array(min1['weights'])
    norm = min2.groupby(['t']).agg(
        norm=('weights', 'sum')).reset_index()
    min2 = min2.merge(norm, on=['t'])
    min2['weights'] = ((np.array(min2['weights'])
                                    / np.array(min2['norm'])))
    min2['expmin'] = np.array(min2['min']) *\
        np.array(min2['weights'])
    norm = min3.groupby(['t']).agg(
        norm=('weights', 'sum')).reset_index()
    min3 = min3.merge(norm, on=['t'])
    min3['weights'] = ((np.array(min3['weights'])
                                    / np.array(min3['norm'])))
    min3['expmin'] = np.array(min3['min']) *\
        np.array(min3['weights'])
    norm = minAtLeast2.groupby(['t']).agg(
        norm=('weights', 'sum')).reset_index()
    minAtLeast2 = minAtLeast2.merge(norm, on=['t'])
    minAtLeast2['weights'] = \
        ((np.array(minAtLeast2['weights'])
        / np.array(minAtLeast2['norm'])))
    minAtLeast2['expmin'] = np.array(minAtLeast2['min']) *\
        np.array(minAtLeast2['weights'])
    norm = minAtLeast3.groupby(['t']).agg(
        norm=('weights', 'sum')).reset_index()
    minAtLeast3 = minAtLeast3.merge(norm, on=['t'])
    minAtLeast3['weights'] = \
        ((np.array(minAtLeast3['weights'])
        / np.array(minAtLeast3['norm'])))
    minAtLeast3['expmin'] = np.array(minAtLeast3['min']) *\
        np.array(minAtLeast3['weights'])
    #
    expmin1 = min1[['t', 'vmatch_id', 'expmin']]\
        .groupby(['t']).agg(exp=('expmin', 'sum')).reset_index()
    expmin2 = min2[['t', 'vmatch_id', 'expmin']]\
        .groupby(['t']).agg(exp=('expmin', 'sum')).reset_index()
    expmin3 = min3[['t', 'vmatch_id', 'expmin']]\
        .groupby(['t']).agg(exp=('expmin', 'sum')).reset_index()
    expminAtLeast2 = minAtLeast2[['t', 'vmatch_id', 'expmin']]\
        .groupby(['t']).agg(exp=('expmin', 'sum')).reset_index()
    expminAtLeast3 = minAtLeast3[['t', 'vmatch_id', 'expmin']]\
        .groupby(['t']).agg(exp=('expmin', 'sum')).reset_index()
    ##
    ##
    fig, ax = plt.subplots(2,sharex=True,figsize=(10,5))
    axes = [ax[0], ax[0].twinx(), ax[0].twinx(), ax[1]]
    microbe_stacked = bAbunds.pivot(index='t',columns='tree_bstrain_id',values='abundance')
    microbe_stacked.plot.area(ax = axes[0],stacked=True,legend=False, linewidth=0,color=bSpeciesColorDict,sort_columns=True,alpha=0.25)
    microbe_stacked.plot(stacked=True, ax=axes[0], legend=False, color='white',sort_columns=True,linewidth=.1)
    axes[1].plot(virus_total['t'], virus_total['vtotal'],
                        linewidth=0, color='grey')
    axes[1].fill_between(virus_total['t'], virus_total['vtotal'], color='grey', alpha=0.75)
    axes[2].plot(sorted(pEmergeExpected['t'].unique()),np.array(pEmergeExpected['p_exp']),\
                    color='darkblue',linewidth=1.5)
    pemergetrunc = pEmergeExpected.merge(expmin1, on=['t'])
    axes[2].plot(sorted(pemergetrunc['t'].unique()),
            1-np.array(pemergetrunc['exp']),
            color='darkred', linewidth=1.5, label='Rank 1')
    lim = axes[0].get_ylim()
    axes[0].set_ylim(2,lim[1])
    lim = axes[1].get_ylim()
    axes[1].set_ylim(0,lim[1])
    lim = axes[2].get_ylim()
    axes[2].set_ylim(0, lim[1])
    axes[0].set_yticklabels([])
    axes[0].set_yticks([])
    axes[1].set_yticklabels([])
    axes[1].set_yticks([])
    axes[2].yaxis.tick_left()
    axes[2].yaxis.set_label_position("left")
    axes[2].set_ylabel(ylabel ='Probability of\nEmergence',labelpad=15,fontsize=10)
    axes[3].plot(sorted(pEmergeExpected['t'].unique()), np.array(pEmergeExpected['p_exp']),
                color='darkblue', linewidth=1.5)
    # axes[3].plot(sorted(expmin1['t'].unique()), 1-np.array(expmin1['exp']),
    #              color='darkred', linewidth=1.5, label='Rank 1 Escape')
    # pemergetrunc = pEmergeExpected.merge(expmin1, on=['t'])
    # axes[3].plot(sorted(pemergetrunc['t'].unique()),
    #             1-np.array(pemergetrunc['exp']) - (pemergetrunc['p_exp']),
    #             color='darkred', linewidth=1.5, label='Rank 1')
    # pemergetrunc = pEmergeExpected.merge(expmin2, on=['t'])
    # axes[3].plot(sorted(pemergetrunc['t'].unique()),
    #             1-np.array(pemergetrunc['exp']) - (pemergetrunc['p_exp']),
    #             color='darkorange', linewidth=1.5, label='Rank 2')
    # pemergetrunc = pEmergeExpected.merge(expmin3, on=['t'])
    # axes[3].plot(sorted(pemergetrunc['t'].unique()),
    #             1-np.array(pemergetrunc['exp']) - (pemergetrunc['p_exp']),
    #             color='darkgreen', linewidth=1.5, label='Rank 3')
    ######
    # pemergetrunc = pEmergeExpected.merge(expmin1, on=['t'])
    # axes[3].plot(sorted(pemergetrunc['t'].unique()),
    #             1-np.array(pemergetrunc['exp']),
    #             color='darkred', linewidth=1.5, label='Rank 1')
    # pemergetrunc = pEmergeExpected.merge(expmin2, on=['t'])
    # axes[3].plot(sorted(pemergetrunc['t'].unique()),
    #             1-np.array(pemergetrunc['exp']),
    #             color='darkorange', linewidth=1.5, label='Rank 2')
    # pemergetrunc = pEmergeExpected.merge(expmin3, on=['t'])
    # axes[3].plot(sorted(pemergetrunc['t'].unique()),
    #             1-np.array(pemergetrunc['exp']),
    #             color='darkgreen', linewidth=1.5, label='Rank 3')
    ######
    avg = bstrain0vstrains[['t', 'vmatch_id', 'weights']]\
        .drop_duplicates().merge(escapeProbs, on=['t', 'vmatch_id'])
    norm = avg.groupby(['t']).agg(
        norm=('weights', 'sum')).reset_index()
    avg = avg.merge(norm, on=['t'])
    avg['weights'] = ((np.array(avg['weights'])
                        / np.array(avg['norm'])))
    avg['exp'] = np.array(avg['p_extinction_lambert'])*\
        np.array(avg['weights'])
    avg = avg[['t', 'vmatch_id', 'exp']]\
        .groupby(['t']).agg(exp=('exp', 'sum')).reset_index()
    axes[3].plot(sorted(avg['t'].unique()),
                    1-np.array(avg['exp']),
                    color='darkred', linewidth=1.5, label='Rank 1')
    ##
    axes[3].legend()
    lim = axes[3].get_ylim()
    axes[3].set_ylim(0, lim[1])
    axes[3].yaxis.tick_left()
    axes[3].yaxis.set_label_position("left")
    axes[3].set_ylabel(ylabel ='Emergence Differential\nupon Escape',labelpad=15,fontsize=10)
    axes[3].set_xlabel(xlabel = 'Time t',fontsize=10)
    axes[0].xaxis.set_minor_locator(ticker.MultipleLocator(25))
    axes[0].margins(x=0)
    axes[1].margins(x=0)
    axes[2].margins(x=0)
    axes[3].margins(x=0)
    fig.tight_layout()
    ##
    fig, ax = plt.subplots(2,sharex=True,figsize=(10,5))
    axes = [ax[0], ax[0].twinx(), ax[0].twinx(), ax[1]]
    microbe_stacked = bAbunds.pivot(index='t',columns='tree_bstrain_id',values='abundance')
    microbe_stacked.plot.area(ax = axes[0],stacked=True,legend=False, linewidth=0,color=bSpeciesColorDict,sort_columns=True,alpha=0.25)
    microbe_stacked.plot(stacked=True, ax=axes[0], legend=False, color='white',sort_columns=True,linewidth=.1)
    axes[1].plot(virus_total['t'], virus_total['vtotal'],
                        linewidth=0, color='grey')
    axes[1].fill_between(virus_total['t'], virus_total['vtotal'], color='grey', alpha=0.75)
    axes[2].plot(sorted(pEmergeExpected['t'].unique()),np.array(pEmergeExpected['p_exp']),\
                    color='darkblue',linewidth=1.5)
    lim = axes[0].get_ylim()
    axes[0].set_ylim(2,lim[1])
    lim = axes[1].get_ylim()
    axes[1].set_ylim(0,lim[1])
    lim = axes[2].get_ylim()
    axes[2].set_ylim(0, lim[1])
    axes[0].set_yticklabels([])
    axes[0].set_yticks([])
    axes[1].set_yticklabels([])
    axes[1].set_yticks([])
    axes[2].yaxis.tick_left()
    axes[2].yaxis.set_label_position("left")
    axes[2].set_ylabel(ylabel ='Probability of\nEmergence',labelpad=15,fontsize=10)
    axes[3].plot(sorted(pEmergeExpected['t'].unique()), np.array(pEmergeExpected['p_exp']),
                color='darkblue', linewidth=1.5)
    # axes[3].plot(sorted(expmin1['t'].unique()), 1-np.array(expmin1['exp']),
    #              color='darkred', linewidth=1.5, label='Rank 1 Escape')
    # pemergetrunc = pEmergeExpected.merge(expmin1, on=['t'])
    # axes[3].plot(sorted(pemergetrunc['t'].unique()),
    #             1-np.array(pemergetrunc['exp']) - (pemergetrunc['p_exp']),
    #             color='darkred', linewidth=1.5, label='Rank 1')
    axes[3].plot(sorted(pemergetrunc['t'].unique()),
                    1-np.array(pemergetrunc['exp']),
                    color='darkred', linewidth=1.5, label='Rank 1')
    axes[3].legend()
    lim = axes[3].get_ylim()
    axes[3].set_ylim(0, lim[1])
    axes[3].yaxis.tick_left()
    axes[3].yaxis.set_label_position("left")
    axes[3].set_ylabel(ylabel ='Emergence Differential\nupon Escape',labelpad=15,fontsize=10)
    axes[3].set_xlabel(xlabel = 'Time t',fontsize=10)
    axes[0].xaxis.set_minor_locator(ticker.MultipleLocator(25))
    axes[0].margins(x=0)
    axes[1].margins(x=0)
    axes[2].margins(x=0)
    axes[3].margins(x=0)
    fig.tight_layout()
    ###
    fig, ax = plt.subplots(1,figsize=(10,5))
    axes = [ax, ax.twinx(), ax.twinx()]
    microbe_stacked = bAbunds.pivot(index='t',columns='tree_bstrain_id',values='abundance')
    microbe_stacked.plot.area(ax = axes[0],stacked=True,legend=False, linewidth=0,color=bSpeciesColorDict,sort_columns=True,alpha=0.25)
    microbe_stacked.plot(stacked=True, ax=axes[0], legend=False, color='white',sort_columns=True,linewidth=.1)
    axes[1].plot(virus_total['t'], virus_total['vtotal'],
                        linewidth=0, color='grey')
    axes[1].fill_between(virus_total['t'], virus_total['vtotal'], color='grey', alpha=0.75)
    axes[2].plot(sorted(pEmergeExpected['t'].unique()),np.array(pEmergeExpected['p_exp']),\
                    color='darkblue',linewidth=1.5)
    lim = axes[0].get_ylim()
    axes[0].set_ylim(2,lim[1])
    lim = axes[1].get_ylim()
    axes[1].set_ylim(0,lim[1])
    lim = axes[2].get_ylim()
    axes[2].set_ylim(0, lim[1])
    axes[0].set_yticklabels([])
    axes[0].set_yticks([])
    axes[1].set_yticklabels([])
    axes[1].set_yticks([])
    axes[2].yaxis.tick_left()
    axes[2].yaxis.set_label_position("left")
    axes[2].set_ylabel(ylabel ='Probability of\nEmergenc (100 cutoff)',labelpad=15,fontsize=10)
    axes[0].xaxis.set_minor_locator(ticker.MultipleLocator(25))
    axes[0].margins(x=0)
    axes[1].margins(x=0)
    axes[2].margins(x=0)
    fig.tight_layout()


if sys.argv[2] == 'emergenceheatmap':
    bAbunds, bKeepTreeStrainsDF, bSpeciesColorDict, _, _, _, _ = \
        tree.speciesTreePlot(run_id, 'microbe',   DBSIM_PATH,  DBTREE_PATH,
                                treepalette, maxticksize, figxy, hratio, babundthreshold, stacked, overlay)
    vAbunds, vKeepTreeStrainsDF, vSpeciesColorDict, hlinecVirus, vlinecVirus, _, _ = \
        tree.speciesTreePlot(run_id, 'virus',   DBSIM_PATH,  DBTREE_PATH,
                                treepalette, maxticksize, figxy, hratio, vabundthreshold, stacked, overlay)
    bAbunds = bAbunds[bAbunds.t <= max(vAbunds.t)]
    bAbunds = bAbunds[bAbunds['abundance'] > 0]
    p = curSim.execute('SELECT viral_mutation_rate,viral_burst_size,microbe_carrying_capacity, \
                    spacer_acquisition_prob, adsorption_rate, viral_decay_rate \
                    FROM param_combos WHERE combo_id = {}'.format(combo_id)).fetchall()
    mu = p[0][0]
    beta = p[0][1]
    K = p[0][2]
    q = p[0][3]
    phi = p[0][4]
    d = p[0][5]
    microbe_total = pd.read_sql_query("SELECT t,microbial_abundance FROM summary WHERE run_id = {}".format(run_id), conSim)\
    .rename(columns={"microbial_abundance": "btotal"})
    virus_total = pd.read_sql_query("SELECT t,viral_abundance FROM summary WHERE run_id = {}".format(run_id), conSim)\
    .rename(columns={"viral_abundance": "vtotal"})
    virus_total=virus_total[virus_total["vtotal"]>0]
    microbe_stacked = pd.read_sql_query("SELECT t,bstrain_id,abundance FROM babundance \
                        WHERE run_id = {}".format(run_id), conSim)
    virus_stacked = pd.read_sql_query("SELECT t,vstrain_id,abundance FROM vabundance \
                        WHERE run_id = {}".format(run_id), conSim)
    vmatches = pd.read_sql_query("SELECT t, match_id, vstrain_id\
        FROM vmatches", conTri).rename(columns={'match_id':'vmatch_id'})
    v0matches = pd.read_sql_query("SELECT t, vstrain_id \
        FROM v0matches", conTri)
    v0matches['vmatch_id'] = len(v0matches['vstrain_id'])*[0]
    vmatches = pd.concat([v0matches,vmatches])\
                .sort_values(by=['t','vmatch_id','vstrain_id'])\
                .merge(vKeepTreeStrainsDF,on=['vstrain_id'])
    vAbunds = vAbunds.rename(columns={'tree_vstrain_id':'new_tree_vstrain_id'})
    pExt = pd.read_sql_query("SELECT t, vmatch_id, p_extinction_lambert \
        FROM existing_vmatch_extinction ORDER BY t", conProb)
    pExt['p_emerge'] = 1- pExt['p_extinction_lambert']
    ######
    tripartiteNetwork = pd.read_sql_query("SELECT t, bstrain_id, vstrain_id, time_specific_match_id \
                            FROM bstrain_to_vstrain_matches WHERE match_length = 1", conMatch)
    spacerMatches = pd.read_sql_query("SELECT t, time_specific_match_id, spacer_id \
                            FROM matches_spacers", conMatch)
    tripartiteNetwork = tripartiteNetwork.merge(spacerMatches, on=['t', 'time_specific_match_id'])\
        .drop(columns=['time_specific_match_id'])
    microbe_stacked = pd.read_sql_query("SELECT t,bstrain_id,abundance FROM babundance \
                        WHERE run_id = {}".format(run_id), conSim)
    tripartiteNetwork = tripartiteNetwork.merge(vKeepTreeStrainsDF, on=['vstrain_id'])\
        .drop(columns=['vstrain_id', 'tree_vstrain_id'])\
        .merge(microbe_stacked, on=['t', 'bstrain_id'])\
        .drop_duplicates()
    spacers = tripartiteNetwork[['t', 'spacer_id', 'bstrain_id', 'abundance']]\
        .drop_duplicates().groupby(['t', 'spacer_id'])\
        .agg(abundance=('abundance', 'sum')).reset_index()
    spacers = spacers.groupby(['t'])\
        .agg(total=('abundance', 'sum')).reset_index()\
        .merge(spacers, on=['t'])
    spacers['freq'] = spacers['abundance']/spacers['total']
    spacers['div'] = np.exp(-1*np.array(spacers['freq'])
                            * np.log(np.array(spacers['freq'])))
    tripartiteNetwork = tripartiteNetwork.drop(columns=['bstrain_id', 'abundance'])\
        .drop_duplicates()\
        .merge(spacers[['t', 'spacer_id', 'div']],
                on=['t', 'spacer_id'])
    tripartiteNetwork = tripartiteNetwork.groupby(['t', 'new_tree_vstrain_id'])\
        .agg(matches=('div', 'prod')).reset_index()
    spacers = spacers.groupby(['t'])\
        .agg(div=('div', 'prod')).reset_index()
    ######
    tripartiteNetwork = vAbunds.merge(vmatches,on=['t','new_tree_vstrain_id'])\
                    .merge(pExt,on=['t','vmatch_id']).drop(columns=['vmatch_id'])
    tripartiteNetwork = tripartiteNetwork[['t','new_tree_vstrain_id','p_emerge']]
    ####
    extraTimes = list({*np.arange(min(vAbunds['t']), max(vAbunds['t'])+1, 1)}
                        - {*np.array(tripartiteNetwork['t'])})
    for t in sorted(extraTimes):
        lastTime = max(tripartiteNetwork[tripartiteNetwork.t < t]['t'])
        new = tripartiteNetwork[tripartiteNetwork.t == lastTime].copy()
        new['t'] = np.ones(len(new['t']))*t
        tripartiteNetwork = pd.concat([tripartiteNetwork, new])

    tripartiteNetwork = tripartiteNetwork.sort_values(
        by=['t', 'new_tree_vstrain_id'])
    treeHeat = tripartiteNetwork.pivot_table(
        index='new_tree_vstrain_id', columns='t', values='p_emerge').fillna(0)
    # treeHeat = treeHeat.reindex(index=treeHeat.index[::-1])
    microbe_stacked = bAbunds.pivot(
        index='t', columns='tree_bstrain_id', values='abundance')
    # virus_stacked = vAbunds.pivot(index='t',columns='tree_vstrain_id',values='abundance')
    fig, ax = plt.subplots(2, sharex=True, figsize=(
        10, 10), gridspec_kw={'height_ratios': [1, 3]})
    axes = [ax[0], ax[0].twinx(), ax[0].twinx(), ax[1]]
    microbe_stacked.plot.area(ax=axes[0], stacked=True, legend=False,
                                linewidth=0, cmap='Purples', sort_columns=True, alpha=0.25)
    axes[1].plot(virus_total['t'], virus_total['vtotal'],
                    linewidth=0, color='grey')
    axes[1].fill_between(
        virus_total['t'], virus_total['vtotal'], color='grey', alpha=0.6)
    axes[0].set_yticklabels([])
    axes[0].set_yticks([])
    axes[0].margins(x=0)
    axes[0].xaxis.set_minor_locator(ticker.MultipleLocator(25))
    lim = axes[0].get_ylim()
    axes[0].set_ylim(0, lim[1])
    axes[1].set_yticklabels([])
    axes[1].set_yticks([])
    axes[1].margins(x=0)
    axes[1].xaxis.set_minor_locator(ticker.MultipleLocator(25))
    lim = axes[1].get_ylim()
    axes[1].set_ylim(0, lim[1])
    ###
    axes[2].plot(spacers['t'], spacers['div'], color='darkblue', linewidth=1.5)
    axes[2].yaxis.tick_left()
    axes[2].yaxis.set_label_position("left")
    axes[2].set_ylabel(
        ylabel='Shannon Diversity of\nSingly-matched Protospacers', labelpad=15, fontsize=10)
    axes[2].xaxis.set_minor_locator(ticker.MultipleLocator(25))
    lim = axes[2].get_ylim()
    axes[2].set_ylim(0, lim[1])
    ###
    strainLineages = LineCollection(
        np.array(hlinecVirus) - np.array([[0, 0.5], [0, 0.5]]), linestyles='solid', color='white', linewidths=(.3))
    creationLines = LineCollection(
        np.array(vlinecVirus) - np.array([[0, 0.5], [0, 0.5]]), linestyles='solid', colors='white', linewidths=(.3))
    axes[3].add_collection(strainLineages)
    axes[3].add_collection(creationLines)
    vKeepTreeStrainsDF = vKeepTreeStrainsDF.sort_values(
        by=['new_tree_vstrain_id'])
    strainsOrdered = vKeepTreeStrainsDF['vstrain_id']
    strainsOrdered = list(strainsOrdered[1::])
    # axes[3].pcolor(np.round(treeHeat), cmap='turbo')
    axes[3].pcolor(treeHeat, cmap='turbo', snap=True,
                    linewidth=0, rasterized=True)

    divider = make_axes_locatable(axes[3])
    maxAbundanceSpecies = np.max(vAbunds.abundance.values)
    markerIncSpecies = 100/maxAbundanceSpecies
    axes[3].scatter(vAbunds['t'], np.array(vAbunds['new_tree_vstrain_id'])-0.5, lw=1.5,
                    s=vAbunds['abundance'].values*markerIncSpecies,
                    color='white', marker='|', alpha=0.6)
    axes[0].set_xlim(0, max(treeHeat.columns))
    axes[1].set_xlim(0, max(treeHeat.columns))
    axes[2].set_xlim(0, max(treeHeat.columns))
    axes[3].set_xlim(0, max(treeHeat.columns))
    axes[3].set_xlabel(xlabel='Time t', fontsize=10)
    axes[3].set_yticks(
        np.arange(.5, len(treeHeat.index), 1))
    axes[3].set_yticklabels(strainsOrdered, fontsize=6)
    # axes[2].set_yticks([])
    # axes[2].set_yticklabels([], fontsize=6)
    axes[3].set_ylabel(ylabel='Viral Strain ID', labelpad=15, fontsize=10)
    axes[3].margins(x=0)
    axes[3].margins(y=0)
    fig.tight_layout()
    fig.savefig(os.path.join(
        PLOT_PATH, 'viral-emergence-heat-map.pdf'), dpi=resolve)
    ####
    left = 200
    right = 600
    top = 15
    triNet = tripartiteNetwork[(tripartiteNetwork.matches > 0)
                               & (tripartiteNetwork.t >= left) & (tripartiteNetwork.t <= right)
                               & (tripartiteNetwork.new_tree_vstrain_id <= top)]
    top = max(triNet.new_tree_vstrain_id)
    bottom = min(triNet.new_tree_vstrain_id)
    # vline3d, hline3d = tree.truncateTree(left,right,top,bottom,vlinecVirus,hlinecVirus)
    vline3d, hline3d = truncateTree(
        left, right, top, bottom, vlinecVirus, hlinecVirus)
    zeros = len(np.array(vline3d))*2
    vline3d = np.concatenate((vline3d, np.reshape(
        np.array([np.zeros(zeros)]), (int(zeros/2), 2, 1))), axis=2)
    zeros = len(np.array(hline3d))*2
    # hline3d = np.concatenate((hline3d, np.reshape(
    #     np.array([np.zeros(zeros)]), (int(zeros/2), 2, 1))), axis=2)
    dx = len(triNet.t)*[1]  # Width of each bar
    dy = len(triNet.t)*[1]  # Depth of each bar
    dz = list(triNet.matches)
    xpos = list(np.array(triNet.t)-np.array(dx)/2)
    ypos = list(np.array(triNet.new_tree_vstrain_id) -
                np.array(dy)/2)  # y coordinates of each bar
    # z coordinates of each bar      # Height of each bar
    zpos = len(triNet.t)*[0]
    strainsOrdered3d = vKeepTreeStrainsDF[
        vKeepTreeStrainsDF.new_tree_vstrain_id
        .isin(np.arange(bottom, top+1, 1))]\
        .sort_values(by=['new_tree_vstrain_id'])
    strainsOrdered3d = list(strainsOrdered3d['vstrain_id'])
    ###
    cmap = cm.get_cmap('turbo')
    norm = Normalize(vmin=float(0), vmax=float(max(tripartiteNetwork.matches)))
    colors = cmap(norm(dz))
    fig = plt.figure(figsize=(20, 10))
    ax = fig.add_subplot(111, projection='3d')
    ax.bar3d(xpos, ypos, zpos, dx, dy, dz, color=colors, shade=True, alpha=1)
    creationLines = Line3DCollection(
        vline3d + np.array([[0, 0.25, 0], [0, 0, 0]]),
        linestyles='solid', colors='white', linewidths=(1.25))
    # strainLineages = Line3DCollection(
    #     hline3d + np.array([[0, 0, 0], [0, 0, 0]]),
    #     linestyles='solid', colors='white', linewidths=(1.25))
    # ax.add_collection3d(strainLineages,zs=0,zdir='z')
    ax.add_collection3d(creationLines, zs=0, zdir='z')
    # ax.set_yticks(np.arange(.5, len(treeHeat.index), 1))
    ax.set_yticks(np.arange(bottom, top+1, 1))
    ax.set_yticklabels(strainsOrdered3d, fontsize=7)
    ax.w_zaxis.set_pane_color(cmap(0))
    ax.xaxis._axinfo["grid"].update({"linewidth": 0.05})
    ax.yaxis._axinfo["grid"].update({"linewidth": 0.05})
    ax.xaxis.set_minor_locator(ticker.MultipleLocator(25))
    ax.axes.set_xlim3d(left=left, right=right)
    ax.axes.set_ylim3d(bottom=bottom, top=top)
    lim = ax.get_zlim()
    ax.set_zlim3d(bottom=0, top=lim[1])
    ax.set_xlabel(xlabel="Time t", labelpad=15, fontsize=10)
    ax.set_ylabel(ylabel="Viral Strain ID", labelpad=15, fontsize=10)
    ax.zaxis.set_rotate_label(False)
    ax.set_zlabel(zlabel="Shannon Diversity of\nSingly-matched\nProtospacers",
                  labelpad=30, fontsize=10, rotation=0)
    ax.set_box_aspect((1/2+np.sqrt(5)/2, 1, 1))
    # ax.view_init(elev=25, azim=135)
    ax.view_init(elev=0, azim=118)
    ###
    left = 16
    right = 200
    bottom = 28
    top = 39
    triNet = tripartiteNetwork[(tripartiteNetwork.matches > 0)
                               & (tripartiteNetwork.t >= left) & (tripartiteNetwork.t <= right)
                               & (tripartiteNetwork.new_tree_vstrain_id <= top)
                               & (tripartiteNetwork.new_tree_vstrain_id >= bottom)]
    top = max(triNet.new_tree_vstrain_id)
    bottom = min(triNet.new_tree_vstrain_id)
    # vline3d, hline3d = tree.truncateTree(left,right,top,bottom,vlinecVirus,hlinecVirus)
    vline3d, hline3d = tree.truncateTree(
        left, right, top, bottom, vlinecVirus, hlinecVirus)
    zeros = len(np.array(vline3d))*2
    vline3d = np.concatenate((vline3d, np.reshape(
        np.array([np.zeros(zeros)]), (int(zeros/2), 2, 1))), axis=2)
    zeros = len(np.array(hline3d))*2
    hline3d = np.concatenate((hline3d, np.reshape(
        np.array([np.zeros(zeros)]), (int(zeros/2), 2, 1))), axis=2)
    dx = len(triNet.t)*[1]  # Width of each bar
    dy = len(triNet.t)*[1]  # Depth of each bar
    dz = list(triNet.matches)
    xpos = list(np.array(triNet.t)-np.array(dx)/2)
    ypos = list(np.array(triNet.new_tree_vstrain_id) -
                np.array(dy)/2)  # y coordinates of each bar
    # z coordinates of each bar      # Height of each bar
    zpos = len(triNet.t)*[0]
    strainsOrdered3d = vKeepTreeStrainsDF[
        vKeepTreeStrainsDF.new_tree_vstrain_id
        .isin(np.arange(bottom, top+1, 1))]\
        .sort_values(by=['new_tree_vstrain_id'])
    strainsOrdered3d = list(strainsOrdered3d['vstrain_id'])
    colors = cmap(norm(dz))
    ###
    fig = plt.figure(figsize=(20, 10))
    ax = fig.add_subplot(111, projection='3d')
    ax.bar3d(xpos, ypos, zpos, dx, dy, dz, color=colors, shade=True, alpha=1)
    creationLines = Line3DCollection(
        vline3d + np.array([[0, 0.25, 0], [0, 0, 0]]),
        linestyles='solid', colors='white', linewidths=(1.25))
    ax.add_collection3d(creationLines, zs=0, zdir='z')
    # ax.set_yticks(np.arange(.5, len(treeHeat.index), 1))
    ax.set_yticks(np.arange(bottom, top+1, 1))
    ax.set_yticklabels(strainsOrdered3d, fontsize=7)
    ax.w_zaxis.set_pane_color(cmap(0))
    ax.xaxis._axinfo["grid"].update({"linewidth": 0.05})
    ax.yaxis._axinfo["grid"].update({"linewidth": 0.05})
    ax.xaxis.set_minor_locator(ticker.MultipleLocator(25))
    ax.axes.set_xlim3d(left=left, right=right)
    ax.axes.set_ylim3d(bottom=bottom, top=top)
    lim = ax.get_zlim()
    ax.set_zlim3d(bottom=0, top=lim[1])
    ax.set_xlabel(xlabel="Time t", labelpad=15, fontsize=10)
    ax.set_ylabel(ylabel="Viral Strain ID", labelpad=15, fontsize=10)
    ax.zaxis.set_rotate_label(False)
    ax.set_zlabel(zlabel="Shannon Diversity of\nSingly-matched\nProtospacers",
                  labelpad=30, fontsize=10, rotation=0)
    ax.set_box_aspect((1/2+np.sqrt(5)/2, 1, 1))
    ax.view_init(elev=20, azim=118)
    plt.close('all')
    sc = cm.ScalarMappable(cmap=cmap, norm=norm)
    sc.set_array([])
    # fig.subplots_adjust(right=0.8)
    # put colorbar at desire position
    fig = plt.figure(figsize=(10, 10))
    cbar_ax = fig.add_axes([0.85, 0.1, 0.02, 0.8])
    fig.colorbar(sc, cax=cbar_ax, ax=ax)


if sys.argv[2] == 'R':
    bAbunds, bKeepTreeStrainsDF, bSpeciesColorDict, _, _, _, _  = \
    tree.speciesTreePlot(run_id,'microbe',DBSIM_PATH,DBTREE_PATH,\
    treepalette,maxticksize,figxy,hratio,babundthreshold,stacked,overlay)
    p = curSim.execute('SELECT viral_mutation_rate,viral_burst_size,microbe_carrying_capacity, \
                        spacer_acquisition_prob, adsorption_rate, viral_decay_rate \
                        FROM param_combos WHERE combo_id = {}'.format(combo_id)).fetchall()
    mu = p[0][0]
    beta = p[0][1]
    K = p[0][2]
    q = p[0][3]
    phi = p[0][4]
    d = p[0][5]
    vmatches = pd.read_sql_query("SELECT t, match_id, vstrain_id\
        FROM vmatches", conTri).rename(columns={'match_id':'vmatch_id'})
    v0matches = pd.read_sql_query("SELECT t, vstrain_id \
        FROM v0matches", conTri)
    v0matches['vmatch_id'] = len(v0matches['vstrain_id'])*[0]
    vmatches = pd.concat([v0matches,vmatches])\
                .sort_values(by=['t','vmatch_id','vstrain_id'])
    virus_stacked = pd.read_sql_query("SELECT t,vstrain_id,abundance FROM vabundance \
                    WHERE run_id = {}".format(run_id), conSim)
    microbe_total = pd.read_sql_query("SELECT t,microbial_abundance FROM summary \
                WHERE run_id = {}".format(run_id), conSim)\
                    .rename(columns={"microbial_abundance": "btotal"})
    virus_total = pd.read_sql_query("SELECT t,viral_abundance FROM summary \
                WHERE run_id = {}".format(run_id), conSim)\
                .rename(columns={"viral_abundance": "vtotal"})
    virus_total = virus_total[virus_total.vtotal>0]
    vfreq = vmatches.merge(virus_stacked,on=['t','vstrain_id'])\
                    .merge(virus_total,on=['t'])
    vfreq['vfreq'] = np.array(vfreq['abundance'])/np.array(vfreq['vtotal'])
    vfreq = vfreq.groupby(['t','vmatch_id','vtotal']).agg(
        vfreq=('vfreq', 'sum'),vabundance=('abundance','sum')).reset_index()
    ##
    vmatchExt = pd.read_sql_query("SELECT t, vmatch_id, p_extinction_lambert, lysis, death \
    FROM existing_vmatch_extinction", conProb).merge(vfreq[['t','vmatch_id','vfreq']],on=['t','vmatch_id'])
    vmatchExt['R0'] = np.array(vmatchExt['lysis'])/np.array(vmatchExt['death'])
    vmatchExt['death'] = vmatchExt['death'] - d
    vmatchExt = vmatchExt.merge(microbe_total,on=['t'])
    # vmatchExt['weights'] = np.array(vmatchExt['vfreq'])*\
    #     (np.array(vmatchExt['lysis'])+np.array(vmatchExt['death']))/np.array(vmatchExt['btotal'])
    vmatchExt['weights'] = np.array(vmatchExt['vfreq'])*(np.array(vmatchExt['lysis']))/np.array(vmatchExt['btotal'])
            #+np.array(vmatchExt['death']))\
    norm = vmatchExt.groupby(['t']).agg(
        norm=('weights', 'sum')).reset_index()
    vmatchExt = vmatchExt.merge(norm, on=['t'])
    vmatchExt['weights'] = np.array(
        vmatchExt['weights'])/np.array(vmatchExt['norm'])
    expR0 = vmatchExt[['t','vmatch_id','weights','R0']].copy()
    expR0['exp'] =np.array(expR0['weights'])*np.array(expR0['R0'])
    expR0 = expR0.groupby(['t']).agg(exp=('exp', 'sum')).reset_index()
    ##
    vescapeExt = pd.read_sql_query("SELECT t, escape_vmatch_id, lysis, death \
    FROM single_escapes_vmatch_extinction", conProb)
    vescapeExt['Resc'] = np.array(vescapeExt['lysis'])/np.array(vescapeExt['death'])
    escapes = pd.read_sql_query("SELECT match_id, escape_match_id \
    FROM potential_vmatch_single_escapes", conProb)\
        .rename(columns={'match_id':'vmatch_id','escape_match_id':'escape_vmatch_id'})
    vescapeExt = vmatchExt[['t','vmatch_id','p_extinction_lambert','R0','weights']]\
                    .merge(escapes, on=['vmatch_id'])\
                    .merge(vescapeExt[['t','escape_vmatch_id','Resc']],\
                            on=['t','escape_vmatch_id'])\
                    .sort_values(by=['t','vmatch_id'])
    norm = vescapeExt[['t','vmatch_id','weights']]\
            .drop_duplicates().groupby(['t']).agg(
            norm=('weights', 'sum')).reset_index()
    vescapeExt = vescapeExt.merge(norm, on=['t'])
    vescapeExt['weights'] = np.array(
        vescapeExt['weights'])/np.array(vescapeExt['norm'])
    vescapeExt = vescapeExt.drop(columns=['norm']).sort_values(by=['t','vmatch_id'])
    #############
    rank1 = vescapeExt.groupby(['t', 'vmatch_id'])\
        .agg(rank1=('Resc', 'max')).reset_index()\
            .merge(vescapeExt,on=['t','vmatch_id'])
    rows = list(np.array(rank1['Resc']) ==
                    np.array(rank1['rank1']))  # keep 1st rank
    rank1 = rank1[rows][['t','vmatch_id','weights','rank1']]\
                .drop_duplicates()  # keep 1st rank
    #
    rank2 = vescapeExt[['t', 'vmatch_id', 'weights','Resc']]\
        .merge(rank1[['t', 'vmatch_id', 'rank1']],
            on=['t', 'vmatch_id'])
    rows = list(np.array(rank2['Resc']) !=
                    np.array(rank2['rank1']))  # remove 1st rank
    rank2 = rank2[rows][['t','vmatch_id','weights','Resc']]
    rank2 = rank2.groupby(['t', 'vmatch_id'])\
        .agg(rank2=('Resc', 'max')).reset_index()\
        .merge(rank2, on=['t', 'vmatch_id'])
    rows = list(np.array(rank2['Resc']) ==
                    np.array(rank2['rank2']))  # keep 2nd rank
    rank2 = rank2[rows][['t','vmatch_id','weights','rank2']]\
            .drop_duplicates()  # keep 2nd rank
    #
    rank3 = vescapeExt[['t', 'vmatch_id', 'weights', 'Resc']]\
        .merge(rank1[['t', 'vmatch_id', 'rank1']],
            on=['t', 'vmatch_id'])
    rows = list(np.array(rank3['Resc']) !=
                    np.array(rank3['rank1']))  # remove 1st rank
    rank3 = rank3[rows][['t','vmatch_id','weights','Resc']]  # remove 1st rank
    rank3 = rank3.merge(rank2[['t', 'vmatch_id', 'rank2']], on=['t', 'vmatch_id'])
    rows = list(np.array(rank3['Resc']) !=
                    np.array(rank3['rank2']))  # remove 2nd rank
    rank3 = rank3[rows][['t','vmatch_id','weights','Resc']]  # remove 2nd rank
    rank3 = rank3.groupby(['t', 'vmatch_id'])\
        .agg(rank3=('Resc', 'max')).reset_index()\
        .merge(rank3, on=['t', 'vmatch_id'])
    rows = list(np.array(rank3['Resc']) ==
                    np.array(rank3['rank3']))  # keep 3rd rank
    rank3 = rank3[rows][['t','vmatch_id','weights','rank3']]\
            .drop_duplicates()  # keep 3rd rank
    ###
    norm = rank1.groupby(['t']).agg(
        norm=('weights', 'sum')).reset_index()
    rank1 = rank1.merge(norm, on=['t'])
    rank1['weights'] = np.array(
        rank1['weights'])/np.array(rank1['norm'])
    rank1 = rank1.drop(columns=['norm'])
    expRank1 = rank1.copy()
    expRank1['exp'] = np.array(expRank1['weights'])*np.array(expRank1['rank1'])
    expRank1 = expRank1.groupby(['t']).agg(
        exp=('exp', 'sum')).reset_index()
    #
    norm = rank2.groupby(['t']).agg(
        norm=('weights', 'sum')).reset_index()
    rank2 = rank2.merge(norm, on=['t'])
    rank2['weights'] = np.array(
        rank2['weights'])/np.array(rank2['norm'])
    rank2 = rank2.drop(columns=['norm'])
    expRank2 = rank2.copy()
    expRank2['exp'] = np.array(expRank2['weights'])*np.array(expRank2['rank2'])
    expRank2 = expRank2.groupby(['t']).agg(
        exp=('exp', 'sum')).reset_index()
    #
    norm = rank3.groupby(['t']).agg(
        norm=('weights', 'sum')).reset_index()
    rank3 = rank3.merge(norm, on=['t'])
    rank3['weights'] = np.array(
        rank3['weights'])/np.array(rank3['norm'])
    rank3 = rank3.drop(columns=['norm'])
    expRank3 = rank3.copy()
    expRank3['exp'] = np.array(expRank3['weights'])*np.array(expRank3['rank3'])
    expRank3 = expRank3.groupby(['t']).agg(
        exp=('exp', 'sum')).reset_index()
    ############
    ############
    microbe_stacked = pd.read_sql_query("SELECT t,bstrain_id,abundance FROM babundance \
                        WHERE run_id = {}".format(run_id), conSim)
    bstrain0vstrains = pd.read_sql_query("SELECT t, bstrain_id, vstrain_id\
        FROM bstrain_to_vstrain_0matches", conMatch)
    vmatches = pd.read_sql_query("SELECT t, match_id, vstrain_id\
        FROM vmatches", conTri).rename(columns={'match_id':'vmatch_id'})
    v0matches = pd.read_sql_query("SELECT t, vstrain_id \
        FROM v0matches", conTri)
    v0matches['vmatch_id'] = len(v0matches['vstrain_id'])*[0]
    vmatches = pd.concat([v0matches,vmatches])\
                .sort_values(by=['t','vmatch_id','vstrain_id'])
    bstrain0vstrains = bstrain0vstrains.merge(vmatches,on=['t','vstrain_id'])\
                        .merge(microbe_stacked,on=['t','bstrain_id'])\
                        .groupby(['t','vmatch_id'])\
                        .agg(babundance=('abundance','sum')).reset_index()\
                        .merge(microbe_total,on=['t'])
    bstrain0vstrains['bfreq'] = \
        np.array(bstrain0vstrains['babundance'])\
        /np.array(bstrain0vstrains['btotal'])
    bstrain0vstrains = bstrain0vstrains.drop(columns=['babundance','btotal'])
    vmatches = vmatches.merge(virus_stacked,on=['t','vstrain_id'])\
            .groupby(['t','vmatch_id'])\
            .agg(vabundance=('abundance','sum')).reset_index()\
            .merge(virus_total,on=['t'])
    vmatches['vfreq'] = np.array(vmatches['vabundance'])/np.array(vmatches['vtotal'])
    bstrain0vstrains = bstrain0vstrains\
                        .merge(vmatches[['t','vmatch_id','vfreq','vabundance']],\
                            on=['t','vmatch_id'])
    bstrain0vstrains['weights'] = \
        np.array(bstrain0vstrains['vfreq'])*np.array(bstrain0vstrains['bfreq'])
    norm = bstrain0vstrains.groupby(['t']).agg(norm=('weights','sum')).reset_index()
    bstrain0vstrains = bstrain0vstrains.merge(norm,on=['t'])
    bstrain0vstrains['weights'] = np.array(bstrain0vstrains['weights'])\
                                            /np.array(bstrain0vstrains['norm'])
    bstrain0vstrains = pd.read_sql_query("SELECT t, vmatch_id, p_extinction_lambert \
        FROM existing_vmatch_extinction ORDER BY t", conProb)\
        .merge(bstrain0vstrains,on=['t','vmatch_id'])
    bstrain0vstrains['exp'] = \
        (1-(np.array(bstrain0vstrains['p_extinction_lambert']))**1)\
        * np.array(bstrain0vstrains['weights'])
    pemerge = bstrain0vstrains.groupby(['t'])\
        .agg(exp=('exp', 'sum')).reset_index()
    #######
    ######
    fig, ax = plt.subplots(2,sharex=True,figsize=(20,8))
    axes = [ax[0], ax[0].twinx(), ax[0].twinx(), ax[1], ax[1].twinx()]
    microbe_stacked = bAbunds[bAbunds.t<=max(virus_total['t'])].pivot(index='t',columns='tree_bstrain_id',values='abundance')
    microbe_stacked.plot.area(ax = axes[0],stacked=True,legend=False, linewidth=0,cmap='Purples',alpha=0.25,sort_columns=True)
    # microbe_stacked.plot(stacked=True, ax=axes[0], legend=False, color='white',sort_columns=True,linewidth=)
    axes[1].plot(virus_total['t'],virus_total['vtotal'],linewidth=0,color='grey')
    axes[1].fill_between(virus_total['t'],virus_total['vtotal'], color='grey',alpha=0.6)
    axes[2].plot(pemerge['t'], np.array(pemerge['exp']), color='darkblue',
                 linewidth=1.5, label=r'$\langle P^*\rangle$')
    axes[3].plot(virus_total['t'],virus_total['vtotal'],linewidth=0,color='grey')
    axes[3].fill_between(virus_total['t'],virus_total['vtotal'], color='grey',alpha=0.6)
    axes[4].plot(expR0['t'],np.ones(len(expR0['t'])), color='grey',linewidth=1, linestyle='dashed')
    axes[4].plot(expR0['t'],beta*np.array(expR0['exp']), color='darkred',linewidth=1.5,label=r'$\langle R_0 \rangle$')
    minT = min(expRank1['t'])
    axes[4].plot(expRank1['t'], beta*np.array(expRank1['exp']) - beta*np.array(expR0[expR0.t>=minT]['exp']), 
                 color='darkgreen', linewidth=1.5, label=r'$\langle R_1 - R_0 \rangle$')
    axes[2].legend(loc='upper right', fontsize=20)
    axes[4].legend(loc='upper right',fontsize=20)
    axes[0].set_yticklabels([])
    axes[0].set_yticks([])
    axes[1].set_yticklabels([])
    axes[1].set_yticks([])
    axes[3].set_yticklabels([])
    axes[3].set_yticks([])
    axes[2].xaxis.set_minor_locator(ticker.MultipleLocator(25))
    axes[4].xaxis.set_minor_locator(ticker.MultipleLocator(25))
    axes[2].yaxis.tick_left()
    axes[2].yaxis.set_label_position("left")
    axes[2].set_ylabel(ylabel=''.join(['Expected Probability of\n',r"Emergence"]),\
                       fontsize=20, labelpad=20)
    axes[4].yaxis.tick_left()
    lim = axes[0].get_ylim()
    axes[0].set_ylim(0,lim[1])
    lim = axes[1].get_ylim()
    axes[1].set_ylim(0,lim[1])
    lim = axes[2].get_ylim()
    axes[2].set_ylim(0, lim[1])
    lim = axes[3].get_ylim()
    axes[3].set_ylim(0, lim[1])
    lim = axes[4].get_ylim()
    axes[4].set_ylim(0, lim[1])
    axes[0].margins(x=0)
    axes[1].margins(x=0)
    axes[2].margins(x=0)
    axes[3].margins(x=0)
    axes[4].margins(x=0)
    axes[4].set_xlabel(xlabel = 'Time t',fontsize=20, labelpad=20)
    axes[3].set_xlabel(xlabel = 'Time t',fontsize=20, labelpad=20)
    axes[0].tick_params(axis='x', labelsize=20)
    axes[3].tick_params(axis='x', labelsize=20)
    axes[2].tick_params(axis='x', labelsize=20)
    axes[4].tick_params(axis='x', labelsize=20)
    axes[2].tick_params(axis='x', labelsize=20)
    axes[4].tick_params(axis='y', labelsize=20)
    axes[2].tick_params(axis='y', labelsize=20)
    fig.tight_layout()
    #####
    #####
    fig, ax = plt.subplots(1,sharex=True,figsize=figxy)
    axes = [ax, ax.twinx(), ax.twinx()]
    microbe_stacked = bAbunds[bAbunds.t<=max(virus_total['t'])].pivot(index='t',columns='tree_bstrain_id',values='abundance')
    microbe_stacked.plot.area(ax = axes[0],stacked=True,legend=False, linewidth=0,cmap='Purples',alpha=0.25,sort_columns=True)
    # microbe_stacked.plot(stacked=True, ax=axes[0], legend=False, color='white',sort_columns=True,linewidth=)
    axes[1].plot(virus_total['t'],virus_total['vtotal'],linewidth=0,color='grey')
    axes[1].fill_between(virus_total['t'],virus_total['vtotal'], color='grey',alpha=0.6)
    axes[2].plot(pemerge['t'],np.array(pemerge['exp']), color='darkblue',linewidth=1.5)
    axes[0].set_yticklabels([])
    axes[0].set_yticks([])
    axes[1].set_yticklabels([])
    axes[1].set_yticks([])
    axes[2].yaxis.tick_left()
    axes[2].yaxis.set_label_position("left")
    axes[2].xaxis.set_minor_locator(ticker.MultipleLocator(25))
    axes[2].set_ylabel(ylabel=''.join(['Expected Probability of\n',r"Emergence $\langle P^*(t)\rangle$"]), fontsize=10)
    lim = axes[0].get_ylim()
    axes[0].set_ylim(0, lim[1])
    lim = axes[1].get_ylim()
    axes[1].set_ylim(0, lim[1])
    lim = axes[2].get_ylim()
    axes[2].set_ylim(0, lim[1])
    axes[0].margins(x=0)
    axes[1].margins(x=0)
    axes[2].margins(x=0)
    axes[0].set_xlabel(xlabel = 'Time t',fontsize=10)
    fig.tight_layout()
    #
    fig, ax = plt.subplots(1, sharex=True, figsize=figxy)
    axes = [ax, ax.twinx()]
    axes[0].plot(virus_total[(virus_total.t<=775)&(virus_total.t>=250)]['t'],\
        virus_total[(virus_total.t<=775)&(virus_total.t>=250)]['vtotal'],linewidth=0,color='grey')
    axes[0].fill_between(virus_total[(virus_total.t<=775)&(virus_total.t>=250)]['t'],\
        virus_total[(virus_total.t<=775)&(virus_total.t>=250)]['vtotal'], color='grey',alpha=0.6)
    axes[1].plot(expR0[(expR0.t<=775)&(expR0.t>=250)]['t'],\
        np.ones(len(expR0[(expR0.t<=775)&(expR0.t>=250)]['t'])), \
            color='grey',linewidth=1, linestyle='dashed')
    axes[1].plot(expR0[(expR0.t<=775)&(expR0.t>=250)]['t'],\
        beta*np.array(expR0[(expR0.t<=775)&(expR0.t>=250)]['exp']), \
            color='black',linewidth=1.5,label=r'$R_0$')
    axes[1].plot(expRank1[(expRank1.t<=775)&(expRank1.t>=250)]['t'], \
                    beta*np.array(expRank1[(expRank1.t <= 775) & (expRank1.t >= 250)]['exp'])\
                        - beta*np.array(expR0[(expR0.t<=775)&(expR0.t>=250)]['exp']),\
                    color='red', linewidth=1.5, label=r'$\rho_1$')
    axes[1].plot(expRank2[(expRank2.t<=775)&(expRank2.t>=250)]['t'],\
        beta*np.array(expRank2[(expRank2.t<=775)&(expRank2.t>=250)]['exp'])\
            - beta*np.array(expR0[(expR0.t<=775)&(expR0.t>=250)]['exp']),color='green',linewidth=1.5, label=r'$\rho_2$')
    axes[1].plot(expRank3[(expRank3.t<=775)&(expRank3.t>=250)]['t'],\
    beta*np.array(expRank3[(expRank3.t<=775)&(expRank3.t>=250)]['exp'])\
        - beta*np.array(expR0[(expR0.t<=775)&(expR0.t>=250)]['exp']), 
            color='blue',linewidth=1.5, label=r'$\rho_3$')
    axes[1].legend(loc='upper left')
    axes[0].set_yticklabels([])
    axes[0].set_yticks([])
    axes[1].yaxis.tick_left()
    axes[1].yaxis.set_label_position("left")
    axes[1].xaxis.set_minor_locator(ticker.MultipleLocator(25))
    axes[1].yaxis.tick_left()
    axes[1].set_xlabel(xlabel = 'Time t',fontsize=10)
    lim = axes[0].get_ylim()
    axes[0].set_ylim(0, lim[1])
    lim = axes[1].get_ylim()
    axes[1].set_ylim(0, lim[1])
    axes[0].margins(x=0)
    axes[1].margins(x=0)
    axes[0].set_xlabel(xlabel = 'Time t',fontsize=10)
    fig.tight_layout()
    fig.savefig(os.path.join(PLOT_PATH,'ranks.pdf'),dpi=resolve)



print('Complete!')
