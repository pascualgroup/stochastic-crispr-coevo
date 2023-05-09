#!/usr/bin/env python3

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import sqlite3
import seaborn as sns
import matplotlib.cm as cm
from matplotlib.colors import Normalize
import matplotlib.colors as mc
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.ticker as ticker
import treefunctions as tree
import generateNetworks as networks
import timestamps as stamps

run = sys.argv[1]

if  sys.argv[2] ==  'all':
    times = [float(sys.argv[i]) for i in list(range(5,len(sys.argv)))]
    print('Compiling abundance plots, trees and networks for {0}'.format(sys.argv[1]))
    print('...for times: {0}...'.format(', '.join(map(str,times))))
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
imgTypes = ["png"]
treeImgTypes = ["png"]
# graphImgTypes = ["pdf"]
figxy = (10,20) # setting for tree abundance figure
hratio = [1,3] # setting for tree abundance figure
maxticksize = 500 # setting for abundances on individual branches of tree
treepalette = 'Spectral' # Color palette for tree: use contiguous color palette
colorjumpV = 0 # this parameter determines whether there should be a jump in the gradient between each clade
# this should equal 0 or greater, no negative values.
colorjumpB = 0 # same as above but for microbes
colorpalSpacers = 'tab20b' # Color palette for spacer in networks: use discrete color palette
babundthreshold  = float(sys.argv[3])/100 # this is %/100 of total population size
vabundthreshold  = float(sys.argv[4])/100 # this is %/100 of total population size
# that a strain has to surpass at its maximum value
# to be included in visualizations
# abundthreshold = 0.05
# graph = float(sys.argv[3])
# justTree = float(sys.argv[4])
hyperAnalyze = True
overlay = True
sSpacing = 2
vSpacing = 2
bSpacing = 2

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
DBMATCH_PATH = os.path.join(SCRIPT_PATH,'isolates',run,'matches_output.sqlite') # cluster
# DBCLADE_PATH = os.path.join(SCRIPT_PATH,'isolates',run,'clade-abundances_output.sqlite') # cluster
# DBTREE_PATH = os.path.join(SCRIPT_PATH,'isolates',run,'trees_output.sqlite') # cluster
DBTRI_PATH = os.path.join(SCRIPT_PATH,'isolates',run,'tripartite-networks_output.sqlite') # cluster
PLOT_PATH = os.path.join(SCRIPT_PATH,'isolates',run,'outputs')
# dir = 'crispr-sweep-7-2-2022/isolates/runID3297-c66-r47'
# run = 'runID3297-c66-r47'
# dir = 'sylvain-martin-collab/runID146-c5-r6'
# run = 'runID146-c5-r6'
# # DBCLADE_PATH = os.path.join('/Volumes','Yadgah',dir,'clade-abundances_output.sqlite')
# DBMATCH_PATH = os.path.join('/Volumes','Yadgah',dir,'matches_output.sqlite') # local
# # DBTREE_PATH = os.path.join('/Volumes','Yadgah',dir,'trees_output.sqlite') # local
# DBTRI_PATH = os.path.join('/Volumes','Yadgah',dir,'tripartite-networks_output.sqlite') # local
# PLOT_PATH = os.path.join('/Volumes','Yadgah') # local
# DBSIM_PATH = os.path.join('/Volumes','Yadgah',dir,'{}.sqlite'.format(run))
# conSim = sqlite3.connect(DBSIM_PATH)
# curSim = conSim.cursor()
# run_id = curSim.execute('SELECT DISTINCT run_id FROM summary').fetchall()
# run_id = run_id[0][0]
# ID = curSim.execute('SELECT combo_id,replicate FROM runs WHERE run_id = {}'.format(run_id)).fetchall()
# combo_id = ID[0][0]
# replicate = ID[0][1]

if not os.path.isdir(PLOT_PATH):
    os.makedirs (PLOT_PATH)

# Call SQLITE analysis databases
conMatch = sqlite3.connect(DBMATCH_PATH)
curMatch = conMatch.cursor()
# conTree = sqlite3.connect(DBTREE_PATH)
# curTree = conTree.cursor()
conTri = sqlite3.connect(DBTRI_PATH)
curTri = conTri.cursor()

# Call dataframes of total abundances
microbe_total = pd.read_sql_query("SELECT t,microbial_abundance FROM summary WHERE run_id = {}".format(run_id), conSim)\
.rename(columns={"microbial_abundance": "Microbial Abundance"})
virus_total = pd.read_sql_query("SELECT t,viral_abundance FROM summary WHERE run_id = {}".format(run_id), conSim)\
.rename(columns={"viral_abundance": "Viral Abundance"})
virus_total=virus_total[virus_total["Viral Abundance"]>0]
microbe_total = microbe_total[microbe_total.t <= max(virus_total.t)]

if sys.argv[2] == 'all':
    vstrainmatchIDs, bstrainmatchIDs = networks.generate(run_id,combo_id,replicate,\
    conSim,conMatch,conTri,times,microbe_total,virus_total,\
    colorpalSpacers,imgTypes,PLOT_PATH,resolve)
    # Generate tree figures and treeID and colorID info
    vAbunds,vKeepTreeStrainsDF, vSpeciesColorDict, vCladeColorDict, _, _, figV, axesV = \
    tree.speciesCladeTreePlot(run_id,'virus',colorjumpV,DBSIM_PATH,DBCLADE_PATH,DBTREE_PATH,\
    treepalette,maxticksize,figxy,hratio,vabundthreshold)
    bAbunds,bKeepTreeStrainsDF, bSpeciesColorDict, bCladeColorDict, _, _, figB, axesB = \
    tree.speciesCladeTreePlot(run_id,'microbe',colorjumpB,DBSIM_PATH,DBCLADE_PATH,DBTREE_PATH,\
    treepalette,maxticksize,figxy,hratio,babundthreshold)
    stamps.treetimes(run_id, curSim,conTri,axesV,axesB,times,hyperAnalyze,\
    vAbunds,bAbunds,vKeepTreeStrainsDF,bKeepTreeStrainsDF)
    for imgType in treeImgTypes:
        if overlay:
            axesB.append(axesB[0].twinx())
            axesB[2].plot(virus_total['t'],virus_total['Viral Abundance'],linewidth=0,color='grey')
            axesB[2].fill_between(virus_total['t'],virus_total['Viral Abundance'], color='grey',alpha=0.6)
            lim = axesB[2].get_ylim()
            axesB[2].set_ylim(0,lim[1])
            figB.tight_layout()
            figB.savefig(os.path.join(PLOT_PATH,'microbe-clades-abundances-tree-time-slices_virus-overlay.{0}'\
                .format(imgType)),dpi=resolve)
            axesV.append(axesV[0].twinx())
            axesV[2].plot(microbe_total['t'],microbe_total['Microbial Abundance'],linewidth=0,color='grey')
            axesV[2].fill_between(microbe_total['t'],microbe_total['Microbial Abundance'], color='grey',alpha=0.4)
            lim = axesV[2].get_ylim()
            axesV[2].set_ylim(0,lim[1])
            figV.tight_layout()
            figV.savefig(os.path.join(PLOT_PATH,'virus-clades-abundances-tree-time-slices_microbe-overlay.{0}'\
                .format(imgType)),dpi=resolve)
        else:
            figB.tight_layout()
            figB.savefig(os.path.join(PLOT_PATH,'microbe-clades-abundances-tree-time-slices.{0}'\
                .format(imgType)),dpi=resolve)
            figV.tight_layout()
            figV.savefig(os.path.join(PLOT_PATH,'virus-clades-abundances-tree-time-slices.{0}'\
                .format(imgType)),dpi=resolve)



if sys.argv[2] == 'networks':
    vstrainmatchIDs, bstrainmatchIDs = networks.generate(run_id,combo_id,replicate,\
    conSim,conMatch,conTri,times,microbe_total,virus_total,\
    colorpalSpacers,imgTypes,PLOT_PATH,resolve)


if sys.argv[2] == 'trees':
    # Generate tree figures and treeID and colorID info
    vAbunds,vKeepTreeStrainsDF, vSpeciesColorDict, vCladeColorDict, _, _, figV, axesV = \
    tree.speciesCladeTreePlot(run_id,'virus',colorjumpV,DBSIM_PATH,DBCLADE_PATH,DBTREE_PATH,\
    treepalette,maxticksize,figxy,hratio,vabundthreshold)
    bAbunds,bKeepTreeStrainsDF, bSpeciesColorDict, bCladeColorDict, _, _, figB, axesB = \
    tree.speciesCladeTreePlot(run_id,'microbe',colorjumpB,DBSIM_PATH,DBCLADE_PATH,DBTREE_PATH,\
    treepalette,maxticksize,figxy,hratio,babundthreshold)
    if len(sys.argv) > 5:
        stamps.treetimes(run_id,curSim,conTri,axesV,axesB,times,False,\
        vAbunds,bAbunds,vKeepTreeStrainsDF,bKeepTreeStrainsDF)
    for imgType in treeImgTypes:
        if overlay:
            axesB.append(axesB[0].twinx())
            axesB[2].plot(virus_total['t'],virus_total['Viral Abundance'],linewidth=0,color='grey')
            axesB[2].fill_between(virus_total['t'],virus_total['Viral Abundance'], color='grey',alpha=0.6)
            lim = axesB[2].get_ylim()
            axesB[2].set_ylim(0,lim[1])
            figB.tight_layout()
            figB.savefig(os.path.join(PLOT_PATH,'microbe-clades-abundances-tree-time-slices_virus-overlay.{0}'\
                .format(imgType)),dpi=resolve)
            axesV.append(axesV[0].twinx())
            axesV[2].plot(microbe_total['t'],microbe_total['Microbial Abundance'],linewidth=0,color='grey')
            axesV[2].fill_between(microbe_total['t'],microbe_total['Microbial Abundance'], color='grey',alpha=0.4)
            lim = axesV[2].get_ylim()
            axesV[2].set_ylim(0,lim[1])
            figV.tight_layout()
            figV.savefig(os.path.join(PLOT_PATH,'virus-clades-abundances-tree-time-slices_microbe-overlay.{0}'\
                .format(imgType)),dpi=resolve)
        else:
            figB.tight_layout()
            figB.savefig(os.path.join(PLOT_PATH,'microbe-clades-abundances-tree-time-slices.{0}'\
                .format(imgType)),dpi=resolve)
            figV.tight_layout()
            figV.savefig(os.path.join(PLOT_PATH,'virus-clades-abundances-tree-time-slices.{0}'\
                .format(imgType)),dpi=resolve)


if sys.argv[2] == 'abundances':
    bAbunds, _, bSpeciesColorDict, _, _, _, _, _  = \
    tree.speciesCladeTreePlot(run_id,'microbe',colorjumpB,DBSIM_PATH,DBCLADE_PATH,DBTREE_PATH,\
    treepalette,maxticksize,figxy,hratio,babundthreshold)
    vAbunds, _, vSpeciesColorDict, _, _, _, _, _ = \
    tree.speciesCladeTreePlot(run_id,'virus',colorjumpV,DBSIM_PATH,DBCLADE_PATH,DBTREE_PATH,\
    treepalette,maxticksize,figxy,hratio,vabundthreshold)
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
    ax.set_ylabel(ylabel = 'Microbial Abundance',labelpad=15,fontsize=15)
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
    fig, ax = plt.subplots(2,sharex=True)
    axes = [ax[0], ax[1]]
    axes[0].set_ylabel(ylabel ='Microbial Abundance',labelpad=15,fontsize=7)
    axes[0].set_xlabel(xlabel = 'Time t',fontsize=7)
    axes[0].ticklabel_format(axis = 'y',style='sci',scilimits=(0,0))
    axes[0].xaxis.set_minor_locator(ticker.MultipleLocator(25))
    axes[1].set_ylabel(ylabel ='Viral Abundance',labelpad=15,fontsize=7)
    axes[1].set_xlabel(xlabel = 'Time t',fontsize=7)
    axes[1].ticklabel_format(axis = 'y',style='sci',scilimits=(0,0))
    axes[1].xaxis.set_minor_locator(ticker.MultipleLocator(25))
    microbe_stacked.plot.area(ax = axes[0],stacked=True,legend=False, linewidth=0,color=bSpeciesColorDict,sort_columns=True)
    microbe_stacked.plot(stacked=True, ax=axes[0], legend=False, color='black',sort_columns=True,linewidth=.1)
    virus_stacked.plot.area(ax = axes[1],stacked=True,legend=False, linewidth=0,color=vSpeciesColorDict,sort_columns=True)
    virus_stacked.plot(stacked=True, ax=axes[1], legend=False, color='black',sort_columns=True,linewidth=.1)
    fig.tight_layout()
    fig.savefig(os.path.join(PLOT_PATH,'microbe-virus-stacked-abundances.png'),dpi=resolve)


if sys.argv[2] == 'viralheatmap':
    bAbunds, bKeepTreeStrainsDF, bSpeciesColorDict, _, _, _, _, _  = \
    tree.speciesCladeTreePlot(run_id,'microbe',colorjumpB,DBSIM_PATH,DBCLADE_PATH,DBTREE_PATH,\
    treepalette,maxticksize,figxy,hratio,babundthreshold)
    vAbunds, vKeepTreeStrainsDF, vSpeciesColorDict, _, hlinecVirus, vlinecVirus, _, _ = \
    tree.speciesCladeTreePlot(run_id,'virus',colorjumpV,DBSIM_PATH,DBCLADE_PATH,DBTREE_PATH,\
    treepalette,maxticksize,figxy,hratio,vabundthreshold)
    # bAbunds, bKeepTreeStrainsDF, bSpeciesColorDict, _, _, _, _, _  = \
    # speciesCladeTreePlot(run_id,'microbe',colorjumpB,DBSIM_PATH,DBCLADE_PATH,DBTREE_PATH,\
    # treepalette,maxticksize,figxy,hratio,babundthreshold)
    # vAbunds, vKeepTreeStrainsDF, vSpeciesColorDict, _, hlinecVirus, vlinecVirus, _, _ = \
    # speciesCladeTreePlot(run_id,'virus',colorjumpV,DBSIM_PATH,DBCLADE_PATH,DBTREE_PATH,\
    # treepalette,maxticksize,figxy,hratio,vabundthreshold)
    bAbunds = bAbunds[bAbunds.t <= max(vAbunds.t)]
    bAbunds = bAbunds[bAbunds['abundance']>0]
    bstrain0vstrain = pd.read_sql_query(
    "SELECT t, bstrain_id, vstrain_id \
    FROM bstrain_to_vstrain_0matches", conMatch)
    bstrain0vstrain = bstrain0vstrain.merge(vKeepTreeStrainsDF,\
            on=['vstrain_id']).drop(columns=['vstrain_id','tree_vstrain_id'])
    microbe_stacked = pd.read_sql_query("SELECT t,bstrain_id,abundance FROM babundance \
                        WHERE run_id = {}".format(run_id), conSim)
    bstrain0vstrain = bstrain0vstrain.merge(microbeStacked,on=['t','bstrain_id'])
    bstrain0vstrain = bstrain0vstrain.groupby(['t','new_tree_vstrain_id'])\
                        .agg(susceptible=('abundance','sum')).reset_index()

    treeHeat = bstrain0vstrain.pivot_table(index='new_tree_vstrain_id', columns='t', values='susceptible').fillna(0)
    # treeHeat = treeHeat.reindex(index=treeHeat.index[::-1])
    microbe_stacked = bAbunds.pivot(index='t',columns='tree_bstrain_id',values='abundance')
    # virus_stacked = vAbunds.pivot(index='t',columns='tree_vstrain_id',values='abundance')
    fig, ax = plt.subplots(2,sharex=True, gridspec_kw={'height_ratios': [1,6]})
    axes = [ax[0], ax[0].twinx(), ax[1]]#, ax[2], ax[3]]
    microbe_stacked.plot.area(ax = axes[0],stacked=True,legend=False, linewidth=0,color=bSpeciesColorDict,sort_columns=True)
    axes[0].set_ylabel(ylabel ='Microbial Strain Abundances',labelpad=15,fontsize=7)
    axes[0].set_xlabel(xlabel = 'Time t',fontsize=7)
    axes[0].ticklabel_format(axis = 'y',style='sci',scilimits=(0,0))
    axes[0].xaxis.set_minor_locator(ticker.MultipleLocator(25))
    lim = axes[0].get_ylim()
    axes[0].set_ylim(0,lim[1])
    axes[1].plot(virus_total['t'],virus_total['Viral Abundance'],linewidth=0,color='grey')
    axes[1].fill_between(virus_total['t'],virus_total['Viral Abundance'], color='grey',alpha=0.6)
    lim = axes[1].get_ylim()
    axes[1].set_ylim(0,lim[1])
    axes[2].imshow(treeHeat,cmap='turbo', aspect='auto',origin='lower',\
    extent=[-0.5, max(microbe_stacked.index)+0.5, \
        min(treeHeat.index)-0.5, max(treeHeat.index)+0.5],\
    interpolation=None)
    strainLineages = LineCollection(hlinecVirus, linestyles='solid', color='white',linewidths=(.3))
    creationLines = LineCollection(vlinecVirus, linestyles='solid', colors='white', linewidths=(.3))
    axes[2].add_collection(strainLineages)
    axes[2].add_collection(creationLines)
    vKeepTreeStrainsDF = vKeepTreeStrainsDF.sort_values(by=['new_tree_vstrain_id'])
    strainsOrdered = vKeepTreeStrainsDF['vstrain_id']
    strainsOrdered = list(strainsOrdered[1::])
    # axes[2].set_xticks(np.arange(0,max(microbe_stacked.index)+0.5,1))
    axes[2].set_yticks(np.arange(1,max(treeHeat.index)+0.5,1))
    axes[2].set_yticklabels(strainsOrdered,fontsize=4)
    axes[2].set_ylabel(ylabel ='Viral Strains',labelpad=15,fontsize=10)
    fig.tight_layout()
    fig.savefig(os.path.join(PLOT_PATH,'viral-heat-map.png'),dpi=resolve)

if sys.argv[2] == 'diversity':
    # bAbunds, bKeepTreeStrainsDF, bSpeciesColorDict, _, _, _, _, _  = \
    # tree.speciesCladeTreePlot(run_id,'microbe',colorjumpB,DBSIM_PATH,DBCLADE_PATH,DBTREE_PATH,\
    # treepalette,maxticksize,figxy,hratio,babundthreshold)
    # vAbunds, vKeepTreeStrainsDF, vSpeciesColorDict, _, hlinecVirus, vlinecVirus, _, _ = \
    # tree.speciesCladeTreePlot(run_id,'virus',colorjumpV,DBSIM_PATH,DBCLADE_PATH,DBTREE_PATH,\
    # treepalette,maxticksize,figxy,hratio,vabundthreshold)
    #
    # bAbunds, bKeepTreeStrainsDF, bSpeciesColorDict, _, _, _, _, _  = \
    # speciesCladeTreePlot(run_id,'microbe',colorjumpB,DBSIM_PATH,DBCLADE_PATH,DBTREE_PATH,\
    # treepalette,maxticksize,figxy,hratio,babundthreshold)
    # vAbunds, vKeepTreeStrainsDF, vSpeciesColorDict, _, hlinecVirus, vlinecVirus, _, _ = \
    # speciesCladeTreePlot(run_id,'virus',colorjumpV,DBSIM_PATH,DBCLADE_PATH,DBTREE_PATH,\
    # treepalette,maxticksize,figxy,hratio,vabundthreshold)
    vtotal = pd.read_sql_query("SELECT t, viral_abundance \
        FROM summary WHERE run_id = {}".format(run_id), conSim)\
        .rename(columns={'viral_abundance': 'vtotal'})
    vtotal = vtotal[vtotal['vtotal'] > 0]
    t = max(vtotal['t'].values)
    print('SQLite Query: microbe abundance data')
    bAbunds = pd.read_sql_query("SELECT t,bstrain_id,abundance \
        FROM babundance WHERE run_id = {}".format(run_id), conSim)
    bAbunds = bAbunds[bAbunds.t <= t]
    btotal = pd.read_sql_query("SELECT t, microbial_abundance \
        FROM summary WHERE run_id = {}".format(run_id), conSim)\
        .rename(columns={'microbial_abundance': 'btotal'})
    btotal = btotal[btotal['t'] <= t]
    maxIDs =bAbunds.set_index('t').groupby(['bstrain_id']).agg(t = ('abundance','idxmax'),\
                                maxAbund = ('abundance','max')).reset_index()
    maxIDs = btotal.merge(maxIDs,on=['t'])
    maxIDs['btotal'] = babundthreshold*np.array(maxIDs['btotal'])
    keepStrains = list(maxIDs[maxIDs['maxAbund']>maxIDs['btotal']]['bstrain_id'].values)
    bAbunds = bAbunds[[(i in keepStrains) for i in bAbunds['bstrain_id']]].reset_index(drop=True)
    vAbunds = pd.read_sql_query("SELECT t,vstrain_id,abundance \
        FROM vabundance WHERE run_id = {}".format(run_id), conSim)
    maxIDs = vAbunds.set_index('t').groupby(['vstrain_id']).agg(t = ('abundance','idxmax'),\
                                maxAbund = ('abundance','max')).reset_index()
    maxIDs = vtotal.merge(maxIDs,on=['t'])
    maxIDs['vtotal'] = vabundthreshold*np.array(maxIDs['vtotal'])
    keepStrains = list(maxIDs[maxIDs['maxAbund']>maxIDs['vtotal']]['vstrain_id'].values)
    vAbunds = vAbunds[[(i in keepStrains) for i in vAbunds['vstrain_id']]].reset_index(drop=True)
    bAbunds = bAbunds.rename(columns={'bstrain_id':'tree_bstrain_id'})
    vAbunds = vAbunds.rename(columns={'vstrain_id':'tree_vstrain_id'})
    bSpeciesColorDict = sns.color_palette("tab20b")
    vSpeciesColorDict = sns.color_palette("tab20b")
    #####
    ##### This is richness and shannon diversity of hosts
    ####
    microbe_stacked = pd.read_sql_query("SELECT t,bstrain_id,abundance FROM babundance \
                        WHERE run_id = {}".format(run_id), conSim)
    microbe_stacked = microbe_stacked[microbe_stacked.t<=max(virus_total.t)]
    shannonStrain = []
    for t in virus_total['t']:
        div = list(microbe_stacked[(microbe_stacked['t']==t) & (microbe_stacked['abundance']!=0)]['abundance'])
        div = np.array(div)/microbe_total[microbe_total['t']==t]['Microbial Abundance'].values[0]
        shannonStrain.append(np.exp(-1*sum(np.array(div)*np.array(list(map(np.log,div))))))
    richnessStrain = []
    for t in virus_total['t']:
        div = list(microbe_stacked[(microbe_stacked['t']==t) & (microbe_stacked['abundance']!=0)]['abundance'])
        richnessStrain.append(len(div))
    fig, ax = plt.subplots(2,sharex=True)
    axes = [ax[0], ax[0].twinx(), ax[1], ax[1].twinx(), ax[1].twinx()]
    microbe_stacked = bAbunds[bAbunds.t<=max(virus_total['t'])].pivot(index='t',columns='tree_bstrain_id',values='abundance')
    microbe_stacked.plot.area(ax = axes[0],stacked=True,legend=False, linewidth=0,color=bSpeciesColorDict,sort_columns=True)
    microbe_stacked.plot(stacked=True, ax=axes[0], legend=False, color='black',sort_columns=True,linewidth=.1)
    axes[0].set_ylabel(ylabel ='Host Abundance',labelpad=15,fontsize=7)
    axes[0].set_xlabel(xlabel = 'Time t',fontsize=7)
    axes[0].ticklabel_format(axis = 'y',style='sci',scilimits=(0,0))
    axes[0].xaxis.set_minor_locator(ticker.MultipleLocator(25))
    axes[1].plot(virus_total['t'],virus_total['Viral Abundance'],linewidth=0,color='grey')
    axes[1].fill_between(virus_total['t'],virus_total['Viral Abundance'], color='grey',alpha=0.6)
    axes[0].margins(x=0)
    axes[1].margins(x=0)
    lim = axes[1].get_ylim()
    axes[1].set_ylim(0,lim[1])
    axes[2].fill_between(virus_total['t'],virus_total['Viral Abundance'], color='grey',alpha=0.6)
    lim = axes[2].get_ylim()
    axes[2].set_ylim(0,lim[1])
    axes[2].set_yticklabels([])
    axes[2].set_yticks([])
    axes[2].margins(x=0)
    axes[3].margins(x=0)
    axes[4].margins(x=0)
    axes[3].plot(virus_total['t'],richnessStrain,color='darkred',linewidth=1)
    axes[3].yaxis.tick_right()
    axes[3].set_ylabel(ylabel ='Richness',labelpad=15,fontsize=7)
    axes[4].plot(virus_total['t'],shannonStrain,color='darkblue',linewidth=1)
    axes[4].yaxis.tick_left()
    axes[4].yaxis.set_label_position("left")
    axes[4].set_ylabel(ylabel ='Host Immune Strain\nShannon Diversity',labelpad=15,fontsize=7)
    fig.tight_layout()
    fig.savefig(os.path.join(PLOT_PATH,'microbe-strain-richness-shannon-stacked.png'),dpi=resolve)
    plt.close(fig)
    #####
    ##### This is richness and shannon diversity of viruses
    #####
    virus_stacked = pd.read_sql_query("SELECT t,vstrain_id,abundance FROM vabundance \
                        WHERE run_id = {}".format(run_id), conSim)
    shannonStrainV = []
    for t in virus_total['t']:
        div = list(virus_stacked[(virus_stacked['t']==t) & (virus_stacked['abundance']!=0)]['abundance'])
        div = np.array(div)/virus_total[virus_total['t']==t]['Viral Abundance'].values[0]
        shannonStrainV.append(np.exp(-1*sum(np.array(div)*np.array(list(map(np.log,div))))))
    richnessStrainV = []
    for t in virus_total['t']:
        div = list(virus_stacked[(virus_stacked['t']==t) & (virus_stacked['abundance']!=0)]['abundance'])
        richnessStrainV.append(len(div))
    fig, ax = plt.subplots(2,sharex=True)
    axes = [ax[0], ax[0].twinx(), ax[1], ax[1].twinx(), ax[1].twinx()]
    axes[0].plot(microbe_total['t'],microbe_total['Microbial Abundance'],linewidth=0,color='grey')
    axes[0].fill_between(microbe_total['t'],microbe_total['Microbial Abundance'], color='grey',alpha=0.6)
    axes[0].set_yticklabels([])
    axes[0].set_yticks([])
    virus_stacked = vAbunds[vAbunds.t<=max(virus_total['t'])].pivot(index='t',columns='tree_vstrain_id',values='abundance')
    virus_stacked.plot.area(ax = axes[1],stacked=True,legend=False, linewidth=0,color=vSpeciesColorDict,sort_columns=True)
    virus_stacked.plot(stacked=True, ax=axes[1], legend=False, color='black',sort_columns=True,linewidth=.1)
    axes[1].yaxis.tick_left()
    axes[1].yaxis.set_label_position("left")
    axes[1].set_ylabel(ylabel ='Viral Abundance',labelpad=15,fontsize=7)
    axes[1].set_xlabel(xlabel = 'Time t',fontsize=7)
    axes[1].ticklabel_format(axis = 'y',style='sci',scilimits=(0,0))
    axes[1].xaxis.set_minor_locator(ticker.MultipleLocator(25))
    axes[0].margins(x=0)
    axes[1].margins(x=0)
    lim = axes[0].get_ylim()
    axes[0].set_ylim(0,lim[1])
    lim = axes[1].get_ylim()
    axes[1].set_ylim(0,lim[1])
    axes[2].fill_between(microbe_total['t'],microbe_total['Microbial Abundance'], color='grey',alpha=0.6)
    lim = axes[2].get_ylim()
    axes[2].set_ylim(0,lim[1])
    axes[2].set_yticklabels([])
    axes[2].set_yticks([])
    axes[2].margins(x=0)
    axes[3].margins(x=0)
    axes[4].margins(x=0)
    axes[3].plot(virus_total['t'],richnessStrainV,color='darkred',linewidth=1)
    axes[3].yaxis.tick_right()
    axes[3].set_ylabel(ylabel ='Richness',labelpad=15,fontsize=7)
    axes[4].plot(virus_total['t'],shannonStrainV,color='darkblue',linewidth=1)
    axes[4].yaxis.tick_left()
    axes[4].yaxis.set_label_position("left")
    axes[4].set_ylabel(ylabel ='Viral Strain\nShannon Diversity',labelpad=15,fontsize=7)
    fig.savefig(os.path.join(PLOT_PATH,'virus-strain-richness-shannon-stacked.png'),dpi=resolve)
    plt.close(fig)
    #####
    ##### This is proportion of diversity in single-match tripartite network
    ####
    microbe_stacked = pd.read_sql_query("SELECT t,bstrain_id,abundance FROM babundance \
                        WHERE run_id = {}".format(run_id), conSim)
    microbe_stacked = microbe_stacked[microbe_stacked.t<=max(virus_total.t)]
    tripartiteNets = pd.read_sql_query("SELECT t, vmatch_id, spacer_id, bmatch_id \
        FROM single_match_tripartite_networks", conTri)\
        .drop(columns=['vmatch_id','spacer_id']).drop_duplicates()
    bmatchstrain = pd.read_sql_query("SELECT t, match_id, bstrain_id \
        FROM bmatches", conTri)\
        .rename(columns={"match_id":"bmatch_id"})
    tripartiteNets = tripartiteNets.merge(bmatchstrain,on=['t','bmatch_id'])\
                        .drop(columns=['bmatch_id']).drop_duplicates()
    bfreq = tripartiteNets\
    .merge(microbe_stacked,on=['t','bstrain_id'])\
    .rename(columns={'abundance':'bfreq'})\
    .merge(microbe_total,on=['t']).reset_index(drop=True)
    bfreq['bfreq'] = bfreq['bfreq']/bfreq['Microbial Abundance']
    bfreq = bfreq.drop(columns=['Microbial Abundance'])
    tripartiteNets = tripartiteNets.groupby(['t']).agg(numStrains=('bstrain_id','size'))\
                        .reset_index()
    shannonStrain = []
    for t in sorted(tripartiteNets['t'].unique()):
        div = list(microbe_stacked[(microbe_stacked['t']==t) & (microbe_stacked['abundance']!=0)]['abundance'])
        div = np.array(div)/microbe_total[microbe_total['t']==t]['Microbial Abundance'].values[0]
        shannonStrain.append(np.exp(-1*sum(np.array(div)*np.array(list(map(np.log,div))))))
    shannonStrain2 = []
    for t in sorted(tripartiteNets['t'].unique()):
        div = bfreq[bfreq['t']==t]['bfreq']
        shannonStrain2.append(np.exp(-1*sum(np.array(div)*np.array(list(map(np.log,div))))))
    tripartiteNets['pShannon'] = np.array(shannonStrain2)/np.array(shannonStrain)
    richnessStrain = []
    for t in sorted(tripartiteNets['t'].unique()):
        div = list(microbe_stacked[(microbe_stacked['t']==t) & (microbe_stacked['abundance']!=0)]['abundance'])
        richnessStrain.append(len(div))
    tripartiteNets['pRichness'] = np.array(tripartiteNets['numStrains'])/\
                                    np.array(richnessStrain)
    fig, ax = plt.subplots(2,sharex=True)
    axes = [ax[0], ax[0].twinx(), ax[1], ax[1].twinx(), ax[1].twinx()]
    microbe_stacked = bAbunds[bAbunds.t<=max(virus_total['t'])].pivot(index='t',columns='tree_bstrain_id',values='abundance')
    microbe_stacked.plot.area(ax = axes[0],stacked=True,legend=False, linewidth=0,color=bSpeciesColorDict,sort_columns=True)
    microbe_stacked.plot(stacked=True, ax=axes[0], legend=False, color='black',sort_columns=True,linewidth=.1)
    axes[0].set_ylabel(ylabel ='Host Abundance',labelpad=15,fontsize=7)
    axes[0].set_xlabel(xlabel = 'Time t',fontsize=7)
    axes[0].ticklabel_format(axis = 'y',style='sci',scilimits=(0,0))
    axes[0].xaxis.set_minor_locator(ticker.MultipleLocator(25))
    axes[1].plot(virus_total['t'],virus_total['Viral Abundance'],linewidth=0,color='grey')
    axes[1].fill_between(virus_total['t'],virus_total['Viral Abundance'], color='grey',alpha=0.6)
    axes[0].margins(x=0)
    axes[1].margins(x=0)
    lim = axes[1].get_ylim()
    axes[1].set_ylim(0,lim[1])
    axes[2].fill_between(virus_total['t'],virus_total['Viral Abundance'], color='grey',alpha=0.6)
    lim = axes[2].get_ylim()
    axes[2].set_ylim(0,lim[1])
    axes[2].set_yticklabels([])
    axes[2].set_yticks([])
    axes[2].margins(x=0)
    axes[3].margins(x=0)
    axes[4].margins(x=0)
    axes[3].plot(sorted(tripartiteNets['t'].unique()),tripartiteNets['pRichness'],color='darkred',linewidth=1)
    axes[3].yaxis.tick_right()
    axes[3].set_ylabel(ylabel ='Proportion of Richness',labelpad=15,fontsize=7)
    axes[4].plot(sorted(tripartiteNets['t'].unique()),tripartiteNets['pShannon'],color='darkblue',linewidth=1)
    axes[4].yaxis.tick_left()
    axes[4].yaxis.set_label_position("left")
    axes[4].set_ylabel(ylabel ='Proportion of\nHost Strain Shannon Diversity\nBelonging to Single-Match Tripartite Network',\
                        labelpad=15,fontsize=7)
    fig.tight_layout()
    fig.savefig(os.path.join(PLOT_PATH,'microbe-single-match-tripartite-proportions.png'),dpi=resolve)
    plt.close(fig)
    ##### this is the tripartite host strain and spacer diversity
    microbe_stacked = pd.read_sql_query("SELECT t,bstrain_id,abundance FROM babundance \
                        WHERE run_id = {}".format(run_id), conSim)
    tripartiteNets = pd.read_sql_query("SELECT t \
        FROM single_match_tripartite_networks", conTri)
    microbe_stacked = microbe_stacked[microbe_stacked.t <= max(tripartiteNets.t)]
    bAbunds = bAbunds[bAbunds.t <= max(tripartiteNets.t)]
    bmatchspacers = pd.read_sql_query("SELECT t, bmatch_id \
        FROM single_match_tripartite_networks", conTri).drop_duplicates()
    bmatchspacers = bmatchspacers.merge(\
    pd.read_sql_query("SELECT match_id, bstrain_id \
        FROM bmatches", conTri).rename(columns={"match_id":"bmatch_id"})\
        ).drop(columns=['bmatch_id']).drop_duplicates()
    bfreq = bmatchspacers\
    .merge(microbe_stacked,on=['t','bstrain_id']).groupby(['t','bstrain_id']).agg(bfreq=('abundance','sum'))\
    .reset_index()
    bfreq = bfreq.groupby(['t']).agg(btotal=('bfreq','sum')).reset_index()\
            .merge(bfreq,on=['t'])
    bfreq['bfreq'] = bfreq['bfreq']/bfreq['btotal']
    bfreq = bfreq.drop(columns=['btotal'])
    bmatchspacers = bmatchspacers.merge(bfreq,on=['t','bstrain_id'])

    shannonStrain = []
    for t in sorted(tripartiteNets['t'].unique()):
        div = list(bmatchspacers[(bmatchspacers['t']==t)]['bfreq'])
        shannonStrain.append(np.exp(-1*sum(np.array(div)*np.array(list(map(np.log,div))))))
    richnessStrain = []
    for t in sorted(tripartiteNets['t'].unique()):
        div = list(microbe_stacked[(microbe_stacked['t']==t) & (microbe_stacked['abundance']!=0)]['abundance'])
        richnessStrain.append(len(div))
    bmatchspacers = pd.read_sql_query("SELECT t, spacer_id, bmatch_id \
    FROM single_match_tripartite_networks", conTri).drop_duplicates()
    bfreq = bmatchspacers\
    .merge(\
    pd.read_sql_query("SELECT t,match_id,babundance \
    FROM bmatches_abundances",conTri).rename(columns={"match_id":"bmatch_id"}),\
    on=['t','bmatch_id']).groupby(['t','spacer_id']).agg(bfreq=('babundance','sum'))\
    .reset_index()
    bfreq = bfreq.groupby(['t']).agg(btotal=('bfreq','sum')).reset_index()\
    .merge(bfreq,on=['t'])
    bfreq['bfreq'] = bfreq['bfreq']/bfreq['btotal']
    bfreq = bfreq.drop(columns=['btotal'])
    bmatchspacers = bmatchspacers.merge(bfreq,on=['t','spacer_id'])\
                .drop(columns=['bmatch_id']).drop_duplicates()
    shannonSpacer = []
    for t in sorted(tripartiteNets['t'].unique()):
        div = list(bmatchspacers[(bmatchspacers['t']==t)]['bfreq'])
        shannonSpacer.append(np.exp(-1*sum(np.array(div)*np.array(list(map(np.log,div))))))
    richnessSpacer = []
    for t in sorted(tripartiteNets['t'].unique()):
        div = list(bmatchspacers[(bmatchspacers['t']==t)]['bfreq'])
        richnessSpacer.append(len(div))
    fig, ax = plt.subplots(2,sharex=True)
    axes = [ax[0], ax[0].twinx(), ax[0].twinx(), ax[1], ax[1].twinx(), ax[1].twinx()]
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
    axes[1].plot(sorted(tripartiteNets['t'].unique()),richnessStrain,color='darkred',linewidth=1)
    axes[1].yaxis.tick_right()
    axes[1].set_ylabel(ylabel ='Richness',labelpad=15,fontsize=7)
    axes[1].set_xlabel(xlabel = 'Time t',fontsize=7)
    axes[2].plot(sorted(tripartiteNets['t'].unique()),shannonStrain,color='darkblue',linewidth=1)
    axes[2].yaxis.tick_left()
    axes[2].set_xlabel(xlabel = 'Time t',fontsize=7)
    axes[2].yaxis.set_label_position("left")
    axes[2].set_ylabel(ylabel ='Host Immune Strain\nShannon Diversity',labelpad=15,fontsize=7)
    axes[4].plot(sorted(tripartiteNets['t'].unique()),richnessSpacer,color='darkred',linewidth=1)
    axes[4].yaxis.tick_right()
    axes[4].set_ylabel(ylabel ='Richness',labelpad=15,fontsize=7)
    axes[4].set_xlabel(xlabel = 'Time t',fontsize=7)
    axes[5].plot(sorted(tripartiteNets['t'].unique()),shannonSpacer,color='darkblue',linewidth=1)
    axes[5].yaxis.tick_left()
    axes[5].yaxis.set_label_position("left")
    axes[5].set_ylabel(ylabel ='Single-match Spacer\nShannon Diversity',labelpad=15,fontsize=7)
    axes[5].set_xlabel(xlabel = 'Time t',fontsize=7)
    fig.tight_layout()
    fig.savefig(os.path.join(PLOT_PATH,'tripartite-host-strain-spacer-richness-shannon.png'),dpi=resolve)


if sys.argv[2] == 'escape':
    





if sys.argv[2] == 'meanfitness':
    vtotal = pd.read_sql_query("SELECT t, viral_abundance \
        FROM summary WHERE run_id = {}".format(run_id), conSim)\
        .rename(columns={'viral_abundance': 'vtotal'})
    vtotal = vtotal[vtotal['vtotal'] > 0]
    t = max(vtotal['t'].values)
    print('SQLite Query: microbe abundance data')
    bAbunds = pd.read_sql_query("SELECT t,bstrain_id,abundance \
        FROM babundance WHERE run_id = {}".format(run_id), conSim)
    bAbunds = bAbunds[bAbunds.t <= t]
    btotal = pd.read_sql_query("SELECT t, microbial_abundance \
        FROM summary WHERE run_id = {}".format(run_id), conSim)\
        .rename(columns={'microbial_abundance': 'btotal'})
    btotal = btotal[btotal['t'] <= t]
    maxIDs =bAbunds.set_index('t').groupby(['bstrain_id']).agg(t = ('abundance','idxmax'),\
                                maxAbund = ('abundance','max')).reset_index()
    maxIDs = btotal.merge(maxIDs,on=['t'])
    maxIDs['btotal'] = babundthreshold*np.array(maxIDs['btotal'])
    keepStrains = list(maxIDs[maxIDs['maxAbund']>maxIDs['btotal']]['bstrain_id'].values)
    bAbunds = bAbunds[[(i in keepStrains) for i in bAbunds['bstrain_id']]].reset_index(drop=True)
    vAbunds = pd.read_sql_query("SELECT t,vstrain_id,abundance \
        FROM vabundance WHERE run_id = {}".format(run_id), conSim)
    maxIDs = vAbunds.set_index('t').groupby(['vstrain_id']).agg(t = ('abundance','idxmax'),\
                                maxAbund = ('abundance','max')).reset_index()
    maxIDs = vtotal.merge(maxIDs,on=['t'])
    maxIDs['vtotal'] = vabundthreshold*np.array(maxIDs['vtotal'])
    keepStrains = list(maxIDs[maxIDs['maxAbund']>maxIDs['vtotal']]['vstrain_id'].values)
    vAbunds = vAbunds[[(i in keepStrains) for i in vAbunds['vstrain_id']]].reset_index(drop=True)
    bAbunds = bAbunds.rename(columns={'bstrain_id':'tree_bstrain_id'})
    vAbunds = bAbunds.rename(columns={'vstrain_id':'tree_vstrain_id'})
    bSpeciesColorDict = sns.color_palette("tab20b")
    vSpeciesColorDict = sns.color_palette("tab20b")
    microbe_stacked = pd.read_sql_query("SELECT t,bstrain_id,abundance FROM babundance \
                        WHERE run_id = {}".format(run_id), conSim)
    microbe_stacked = microbe_stacked[microbe_stacked.t<=max(virus_total.t)]
    virus_stacked = pd.read_sql_query("SELECT t,vstrain_id,abundance FROM vabundance \
                        WHERE run_id = {}".format(run_id), conSim)
    bfreq = microbe_stacked.merge(microbe_total, on=['t'])
    vfreq = virus_stacked.merge(virus_total, on=['t'])
    bfreq['bfreq'] = bfreq['abundance']/bfreq['Microbial Abundance']
    vfreq['vfreq'] = vfreq['abundance']/vfreq['Viral Abundance']
    bfreq = bfreq.drop(columns=['abundance','Microbial Abundance'])
    vfreq = vfreq.drop(columns=['abundance','Viral Abundance'])
    zeromatch = pd.read_sql_query("SELECT t, vstrain_id, bstrain_id\
        FROM bstrain_to_vstrain_0matches", conMatch).drop_duplicates()
    maxzeromatch = pd.read_sql_query("SELECT t, vstrain_id, bstrain_id\
        FROM bstrain_to_vstrain_0matches", conMatch).drop_duplicates()
    onematch = pd.read_sql_query("SELECT t, vstrain_id, bstrain_id\
        FROM bstrain_to_vstrain_matches WHERE match_length in (1)", conMatch)\
        .drop_duplicates()
    doublematch = pd.read_sql_query("SELECT t, vstrain_id, bstrain_id\
        FROM bstrain_to_vstrain_matches WHERE match_length in (2)", conMatch)\
        .drop_duplicates()
    triplematch = pd.read_sql_query("SELECT t, vstrain_id, bstrain_id\
        FROM bstrain_to_vstrain_matches WHERE match_length in (3)", conMatch)\
        .drop_duplicates()
    maxzeromatch = maxzeromatch.merge(vfreq,on=['t','vstrain_id'])\
                    .merge(bfreq,on=['t','bstrain_id'])
    maxzeromatch = maxzeromatch.groupby(['t','vstrain_id','vfreq'])\
                        .agg(bfreq=('bfreq','sum')).reset_index()
    maxzeromatch['meanfitness'] = maxzeromatch['vfreq']*maxzeromatch['bfreq']
    maxzeromatch = maxzeromatch.groupby(['t'])\
                        .agg(meanfitness=('meanfitness','sum')).reset_index()
    maxonematch = pd.concat([zeromatch, onematch], sort=False)
    maxonematch.sort_values(by=['t', 'vstrain_id'])
    maxonematch = maxonematch.merge(vfreq,on=['t','vstrain_id'])\
                    .merge(bfreq,on=['t','bstrain_id'])
    maxonematch = maxonematch.groupby(['t','vstrain_id','vfreq'])\
                        .agg(bfreq=('bfreq','sum')).reset_index()
    maxonematch['meanfitness'] = maxonematch['vfreq']*maxonematch['bfreq']
    maxonematch = maxonematch.groupby(['t'])\
                        .agg(meanfitness=('meanfitness','sum')).reset_index()
    maxtwomatch = pd.concat([zeromatch, onematch, doublematch], sort=False)
    maxtwomatch.sort_values(by=['t', 'vstrain_id'])
    maxtwomatch = maxtwomatch.merge(vfreq,on=['t','vstrain_id'])\
                    .merge(bfreq,on=['t','bstrain_id'])
    maxtwomatch = maxtwomatch.groupby(['t','vstrain_id','vfreq'])\
                        .agg(bfreq=('bfreq','sum')).reset_index()
    maxtwomatch['meanfitness'] = maxtwomatch['vfreq']*maxtwomatch['bfreq']
    maxtwomatch = maxtwomatch.groupby(['t'])\
                        .agg(meanfitness=('meanfitness','sum')).reset_index()
    maxthreematch = pd.concat([zeromatch, onematch, doublematch, triplematch], sort=False)
    maxthreematch.sort_values(by=['t', 'vstrain_id'])
    maxthreematch = maxthreematch.merge(vfreq,on=['t','vstrain_id'])\
                    .merge(bfreq,on=['t','bstrain_id'])
    maxthreematch = maxthreematch.groupby(['t','vstrain_id','vfreq'])\
                        .agg(bfreq=('bfreq','sum')).reset_index()
    maxthreematch['meanfitness'] = maxthreematch['vfreq']*maxthreematch['bfreq']
    maxthreematch = maxthreematch.groupby(['t'])\
                        .agg(meanfitness=('meanfitness','sum')).reset_index()
    fig, ax = plt.subplots(2,sharex=True)
    axes = [ax[0], ax[0].twinx(), ax[1], ax[1].twinx()]
    microbe_stacked = bAbunds[bAbunds.t<=max(virus_total['t'])].pivot(index='t',columns='tree_bstrain_id',values='abundance')
    microbe_stacked.plot.area(ax = axes[0],stacked=True,legend=False, linewidth=0,color=bSpeciesColorDict,sort_columns=True)
    microbe_stacked.plot(stacked=True, ax=axes[0], legend=False, color='black',sort_columns=True,linewidth=.1)
    axes[0].set_ylabel(ylabel ='Host Abundance',labelpad=15,fontsize=7)
    axes[0].set_xlabel(xlabel = 'Time t',fontsize=7)
    axes[0].ticklabel_format(axis = 'y',style='sci',scilimits=(0,0))
    axes[0].xaxis.set_minor_locator(ticker.MultipleLocator(25))
    axes[1].plot(virus_total['t'],virus_total['Viral Abundance'],linewidth=0,color='grey')
    axes[1].fill_between(virus_total['t'],virus_total['Viral Abundance'], color='grey',alpha=0.6)
    axes[0].margins(x=0)
    axes[1].margins(x=0)
    lim = axes[1].get_ylim()
    axes[1].set_ylim(0,lim[1])
    axes[2].fill_between(virus_total['t'],virus_total['Viral Abundance'], color='grey',alpha=0.6)
    lim = axes[2].get_ylim()
    axes[2].set_ylim(0,lim[1])
    axes[2].set_yticklabels([])
    axes[2].set_yticks([])
    axes[2].margins(x=0)
    axes[3].margins(x=0)
    lim = axes[3].get_ylim()
    axes[3].set_ylim(0,lim[1])
    axes[3].plot(maxzeromatch['t'],maxzeromatch['meanfitness'],label='0 match max',linewidth=1)
    axes[3].plot(maxonematch['t'],maxonematch['meanfitness'],label='1 match max',linewidth=1)
    axes[3].plot(maxtwomatch['t'],maxtwomatch['meanfitness'],label='2 match max',linewidth=1)
    axes[3].plot(maxthreematch['t'],maxthreematch['meanfitness'],label='3 match max',linewidth=1)
    axes[3].yaxis.tick_left()
    axes[3].yaxis.set_label_position("left")
    axes[3].set_ylabel(ylabel ='Mean Fitness',labelpad=15,fontsize=7)
    axes[3].legend()
    fig.tight_layout()
    fig.savefig(os.path.join(PLOT_PATH,'mean-fitness.png'),dpi=resolve)

vtotal = pd.read_sql_query("SELECT t, viral_abundance \
    FROM summary WHERE run_id = {}".format(run_id), conSim)\
    .rename(columns={'viral_abundance': 'vtotal'})
vtotal = vtotal[vtotal['vtotal'] > 0]
t = max(vtotal['t'].values)
print('SQLite Query: microbe abundance data')
bAbunds = pd.read_sql_query("SELECT t,bstrain_id,abundance \
    FROM babundance WHERE run_id = {}".format(run_id), conSim)
bAbunds = bAbunds[bAbunds.t <= t]
btotal = pd.read_sql_query("SELECT t, microbial_abundance \
    FROM summary WHERE run_id = {}".format(run_id), conSim)\
    .rename(columns={'microbial_abundance': 'btotal'})
btotal = btotal[btotal['t'] <= t]
maxIDs =bAbunds.set_index('t').groupby(['bstrain_id']).agg(t = ('abundance','idxmax'),\
                            maxAbund = ('abundance','max')).reset_index()
maxIDs = btotal.merge(maxIDs,on=['t'])
maxIDs['btotal'] = babundthreshold*np.array(maxIDs['btotal'])
keepStrains = list(maxIDs[maxIDs['maxAbund']>maxIDs['btotal']]['bstrain_id'].values)
bAbunds = bAbunds[[(i in keepStrains) for i in bAbunds['bstrain_id']]].reset_index(drop=True)
vAbunds = pd.read_sql_query("SELECT t,vstrain_id,abundance \
    FROM vabundance WHERE run_id = {}".format(run_id), conSim)
maxIDs = vAbunds.set_index('t').groupby(['vstrain_id']).agg(t = ('abundance','idxmax'),\
                            maxAbund = ('abundance','max')).reset_index()
maxIDs = vtotal.merge(maxIDs,on=['t'])
maxIDs['vtotal'] = vabundthreshold*np.array(maxIDs['vtotal'])
keepStrains = list(maxIDs[maxIDs['maxAbund']>maxIDs['vtotal']]['vstrain_id'].values)
vAbunds = vAbunds[[(i in keepStrains) for i in vAbunds['vstrain_id']]].reset_index(drop=True)
bAbunds = bAbunds.rename(columns={'bstrain_id':'tree_bstrain_id'})
vAbunds = bAbunds.rename(columns={'vstrain_id':'tree_vstrain_id'})
bSpeciesColorDict = sns.color_palette("tab20b")
vSpeciesColorDict = sns.color_palette("tab20b")
microbe_stacked = pd.read_sql_query("SELECT t,bstrain_id,abundance FROM babundance \
                    WHERE run_id = {}".format(run_id), conSim)
microbe_stacked = microbe_stacked[microbe_stacked.t<=max(virus_total.t)]
virus_stacked = pd.read_sql_query("SELECT t,vstrain_id,abundance FROM vabundance \
                    WHERE run_id = {}".format(run_id), conSim)
bfreq = microbe_stacked.merge(microbe_total, on=['t'])
bfreq['bfreq'] = bfreq['abundance']/bfreq['Microbial Abundance']
bfreq = bfreq.drop(columns=['abundance','Microbial Abundance'])
bgrowthrates = pd.read_sql_query("SELECT bstrain_id,growth_rate FROM bgrowthrates \
                    WHERE run_id = {}".format(run_id), conSim)
bgrowthrates = bgrowthrates.merge(bfreq,on=['bstrain_id'])
bgrowthrates['mean'] = bgrowthrates['bfreq']*bgrowthrates['growth_rate']
bgrowthrates =bgrowthrates.groupby(['t']).agg(meangrowthrate = ('mean','sum'))\
                .reset_index()
fig, ax = plt.subplots(2,sharex=True)
axes = [ax[0], ax[0].twinx(), ax[1], ax[1].twinx()]
microbe_stacked = bAbunds[bAbunds.t<=max(virus_total['t'])].pivot(index='t',columns='tree_bstrain_id',values='abundance')
microbe_stacked.plot.area(ax = axes[0],stacked=True,legend=False, linewidth=0,color=bSpeciesColorDict,sort_columns=True)
microbe_stacked.plot(stacked=True, ax=axes[0], legend=False, color='black',sort_columns=True,linewidth=.1)
axes[0].set_ylabel(ylabel ='Host Abundance',labelpad=15,fontsize=7)
axes[0].set_xlabel(xlabel = 'Time t',fontsize=7)
axes[0].ticklabel_format(axis = 'y',style='sci',scilimits=(0,0))
axes[0].xaxis.set_minor_locator(ticker.MultipleLocator(25))
axes[1].plot(virus_total['t'],virus_total['Viral Abundance'],linewidth=0,color='grey')
axes[1].fill_between(virus_total['t'],virus_total['Viral Abundance'], color='grey',alpha=0.6)
axes[0].margins(x=0)
axes[1].margins(x=0)
lim = axes[1].get_ylim()
axes[1].set_ylim(0,lim[1])
axes[2].fill_between(virus_total['t'],virus_total['Viral Abundance'], color='grey',alpha=0.6)
lim = axes[2].get_ylim()
axes[2].set_ylim(0,lim[1])
axes[2].set_yticklabels([])
axes[2].set_yticks([])
axes[2].margins(x=0)
axes[3].margins(x=0)
lim = axes[3].get_ylim()
axes[3].set_ylim(0,lim[1])
axes[3].plot(bgrowthrates['t'],bgrowthrates['meangrowthrate'],linewidth=1)
axes[3].yaxis.tick_left()
axes[3].yaxis.set_label_position("left")
axes[3].set_ylabel(ylabel ='Mean Intrinsic Fitness',labelpad=15,fontsize=7)
fig.tight_layout()
fig.savefig(os.path.join(PLOT_PATH,'mean-intrinsic-fitness.png'),dpi=resolve)









print('Complete!')


            # ##### This is proportion of diversity in higher-order tripartite network
            # microbe_stacked = pd.read_sql_query("SELECT t,bstrain_id,abundance FROM babundance \
            #                     WHERE run_id = {}".format(run_id), conSim)
            # microbe_stacked = microbe_stacked[microbe_stacked.t<=max(virus_total.t)]
            # susceptible = pd.read_sql_query("SELECT t, bstrain_id \
            #     FROM bstrain_to_vstrain_0matches", conMatch).drop_duplicates()
            # tripartiteNets = pd.read_sql_query("SELECT t, bstrain_id, match_length \
            #     FROM bstrain_to_vstrain_matches WHERE match_length = 1", conMatch)\
            #     .drop(columns=['match_length']).drop_duplicates()
            # tripartiteNets2 = pd.read_sql_query("SELECT t, bstrain_id, match_length \
            #     FROM bstrain_to_vstrain_matches WHERE match_length = 2", conMatch)\
            #     .drop(columns=['match_length']).drop_duplicates()
            # tripartiteNets3 = pd.read_sql_query("SELECT t, bstrain_id, match_length \
            #     FROM bstrain_to_vstrain_matches WHERE match_length = 2", conMatch)\
            #     .drop(columns=['match_length']).drop_duplicates()
            # # zero match
            # bfreq = susceptible\
            # .merge(microbe_stacked,on=['t','bstrain_id'])\
            # .rename(columns={'abundance':'bfreq'})\
            # .merge(microbe_total,on=['t']).reset_index(drop=True)
            # bfreq['bfreq'] = bfreq['bfreq']/bfreq['Microbial Abundance']
            # bfreq = bfreq.drop(columns=['Microbial Abundance'])
            # susceptible = susceptible.groupby(['t']).agg(numStrains=('bstrain_id','size'))\
            #                     .reset_index()
            # shannonStrain = []
            # for t in sorted(susceptible['t'].unique()):
            #     div = list(microbe_stacked[(microbe_stacked['t']==t) & (microbe_stacked['abundance']!=0)]['abundance'])
            #     div = np.array(div)/microbe_total[microbe_total['t']==t]['Microbial Abundance'].values[0]
            #     shannonStrain.append(np.exp(-1*sum(np.array(div)*np.array(list(map(np.log,div))))))
            #
            #
            # shannonStrain2 = []
            # for t in sorted(susceptible['t'].unique()):
            #     div = bfreq[bfreq['t']==t]['bfreq']
            #     shannonStrain2.append(np.exp(-1*sum(np.array(div)*np.array(list(map(np.log,div))))))
            #
            #
            # susceptible['pShannon'] = np.array(shannonStrain2)/np.array(shannonStrain)
            # richnessStrain = []
            # for t in sorted(susceptible['t'].unique()):
            #     div = list(microbe_stacked[(microbe_stacked['t']==t) & (microbe_stacked['abundance']!=0)]['abundance'])
            #     richnessStrain.append(len(div))
            #
            #
            # susceptible['pRichness'] = np.array(susceptible['numStrains'])/\
            #                                 np.array(richnessStrain)
            # # double match
            # bfreq = tripartiteNets2\
            # .merge(microbe_stacked,on=['t','bstrain_id'])\
            # .rename(columns={'abundance':'bfreq'})\
            # .merge(microbe_total,on=['t']).reset_index(drop=True)
            # bfreq['bfreq'] = bfreq['bfreq']/bfreq['Microbial Abundance']
            # bfreq = bfreq.drop(columns=['Microbial Abundance'])
            # tripartiteNets2 = tripartiteNets2.groupby(['t']).agg(numStrains=('bstrain_id','size'))\
            #                     .reset_index()
            # timeAll.append(list(sorted(tripartiteNets2['t'].unique())))
            # shannonStrain = []
            # for t in sorted(tripartiteNets2['t'].unique()):
            #     div = list(microbe_stacked[(microbe_stacked['t']==t) & (microbe_stacked['abundance']!=0)]['abundance'])
            #     div = np.array(div)/microbe_total[microbe_total['t']==t]['Microbial Abundance'].values[0]
            #     shannonStrain.append(np.exp(-1*sum(np.array(div)*np.array(list(map(np.log,div))))))
            #
            #
            # shannonStrain2 = []
            # for t in sorted(tripartiteNets2['t'].unique()):
            #     div = bfreq[bfreq['t']==t]['bfreq']
            #     shannonStrain2.append(np.exp(-1*sum(np.array(div)*np.array(list(map(np.log,div))))))
            #
            #
            # tripartiteNets2['pShannon'] = np.array(shannonStrain2)/np.array(shannonStrain)
            # richnessStrain = []
            # for t in sorted(tripartiteNets2['t'].unique()):
            #     div = list(microbe_stacked[(microbe_stacked['t']==t) & (microbe_stacked['abundance']!=0)]['abundance'])
            #     richnessStrain.append(len(div))
            #
            #
            # tripartiteNets2['pRichness'] = np.array(tripartiteNets2['numStrains'])/\
            #                                 np.array(richnessStrain)
            # # triple match
            # bfreq = tripartiteNets3\
            # .merge(microbe_stacked,on=['t','bstrain_id'])\
            # .rename(columns={'abundance':'bfreq'})\
            # .merge(microbe_total,on=['t']).reset_index(drop=True)
            # bfreq['bfreq'] = bfreq['bfreq']/bfreq['Microbial Abundance']
            # bfreq = bfreq.drop(columns=['Microbial Abundance'])
            # tripartiteNets3 = tripartiteNets3.groupby(['t']).agg(numStrains=('bstrain_id','size'))\
            #                     .reset_index()
            # timeAll.append(list(sorted(tripartiteNets3['t'].unique())))
            # shannonStrain = []
            # for t in sorted(tripartiteNets3['t'].unique()):
            #     div = list(microbe_stacked[(microbe_stacked['t']==t) & (microbe_stacked['abundance']!=0)]['abundance'])
            #     div = np.array(div)/microbe_total[microbe_total['t']==t]['Microbial Abundance'].values[0]
            #     shannonStrain.append(np.exp(-1*sum(np.array(div)*np.array(list(map(np.log,div))))))
            #
            #
            # shannonStrain3 = []
            # for t in sorted(tripartiteNets3['t'].unique()):
            #     div = bfreq[bfreq['t']==t]['bfreq']
            #     shannonStrain3.append(np.exp(-1*sum(np.array(div)*np.array(list(map(np.log,div))))))
            #
            #
            # tripartiteNets3['pShannon'] = np.array(shannonStrain2)/np.array(shannonStrain)
            # richnessStrain = []
            # for t in sorted(tripartiteNets3['t'].unique()):
            #     div = list(microbe_stacked[(microbe_stacked['t']==t) & (microbe_stacked['abundance']!=0)]['abundance'])
            #     richnessStrain.append(len(div))
            #
            #
            # tripartiteNets3['pRichness'] = np.array(tripartiteNets3['numStrains'])/\
            #                                 np.array(richnessStrain)
            # fig1, ax1 = plt.subplots(2,sharex=True)
            # fig2, ax2 = plt.subplots(2,sharex=True)
            # fig2, ax3 = plt.subplots(2,sharex=True)
            # microbe_stacked = bAbunds[bAbunds.t<=max(virus_total['t'])].pivot(index='t',columns='tree_bstrain_id',values='abundance')
            # microbe_stacked.plot.area(ax = axes[0],stacked=True,legend=False, linewidth=0,color=bSpeciesColorDict,sort_columns=True)
            # microbe_stacked.plot(stacked=True, ax=axes[0], legend=False, color='black',sort_columns=True,linewidth=.1)
            # axes[0].set_ylabel(ylabel ='Host Abundance',labelpad=15,fontsize=7)
            # axes[0].set_xlabel(xlabel = 'Time t',fontsize=7)
            # axes[0].ticklabel_format(axis = 'y',style='sci',scilimits=(0,0))
            # axes[0].xaxis.set_minor_locator(ticker.MultipleLocator(25))
            # axes[1].plot(virus_total['t'],virus_total['Viral Abundance'],linewidth=0,color='grey')
            # axes[1].fill_between(virus_total['t'],virus_total['Viral Abundance'], color='grey',alpha=0.6)
            # axes[0].margins(x=0)
            # axes[1].margins(x=0)
            # lim = axes[1].get_ylim()
            # axes[1].set_ylim(0,lim[1])
            # axes[2].fill_between(virus_total['t'],virus_total['Viral Abundance'], color='grey',alpha=0.6)
            # lim = axes[2].get_ylim()
            # axes[2].set_ylim(0,lim[1])
            # axes[2].set_yticklabels([])
            # axes[2].set_yticks([])
            # axes[2].margins(x=0)
            # axes[3].margins(x=0)
            # axes[4].margins(x=0)
            # axes[3].plot(sorted(tripartiteNets['t'].unique()),tripartiteNets['pRichness'],color='darkred',linewidth=1)
            # axes[3].yaxis.tick_right()
            # axes[3].set_ylabel(ylabel ='Proportion of Richness',labelpad=15,fontsize=7)
            # axes[4].plot(sorted(tripartiteNets['t'].unique()),tripartiteNets['pShannon'],color='darkblue',linewidth=1)
            # axes[4].yaxis.tick_left()
            # axes[4].yaxis.set_label_position("left")
            # axes[4].set_ylabel(ylabel ='Proportion of\nHost Strain Shannon Diversity\nBelonging to Single-Match Tripartite Network',\
            #                     labelpad=15,fontsize=7)
            # fig.savefig(os.path.join(PLOT_PATH,'microbe-single-match-tripartite-proportions.png'),dpi=resolve)
            # plt.close(fig)

        ##### this is the tripartite strain and spacer diversity
