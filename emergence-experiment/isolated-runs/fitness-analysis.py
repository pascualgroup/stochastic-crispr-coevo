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

run = sys.argv[1]

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
figxy = (15,15) # setting for tree abundance figure
hratio = [1,3] # setting for tree abundance figure
maxticksize = 400 # setting for abundances on individual branches of tree
treepalette = 'turbo' # Color palette for tree: use contiguous color palette
babundthreshold  = float(sys.argv[3])/100 # this is %/100 of total population size
vabundthreshold  = float(sys.argv[4])/100 # this is %/100 of total population size
overlay = True

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
DBTREE_PATH = os.path.join(SCRIPT_PATH,'isolates',run,'trees_output.sqlite') # cluster
DBTRI_PATH = os.path.join(SCRIPT_PATH,'isolates',run,'tripartite-networks_output.sqlite') # cluster
PLOT_PATH = os.path.join(SCRIPT_PATH,'isolates',run,'outputs')
# dir = 'crispr-sweep-7-2-2022/isolates/runID3297-c66-r47'
# run = 'runID3297-c66-r47'
# DBMATCH_PATH = os.path.join('/Volumes','Yadgah',dir,'matches_output.sqlite') # local
# DBTREE_PATH = os.path.join('/Volumes','Yadgah',dir,'trees_output.sqlite') # local
# DBTRI_PATH = os.path.join('/Volumes','Yadgah',dir,'tripartite-networks_output.sqlite') # local
# PLOT_PATH = os.path.join('/Volumes','Yadgah') # local
# DBSIM_PATH = os.path.join('/Volumes','Yadgah',dir,'{}.sqlite'.format(run))
# DBTREE_PATH = os.path.join('/Users/armun/Dropbox/Current/Projects/microbe-virus-crispr/stochastic-crispr/trees_output.sqlite')
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

if sys.argv[2] == 'trees':
    # Generate tree figures and treeID and colorID info
    vAbunds,vKeepTreeStrainsDF, vSpeciesColorDict, _, _, figV, axesV = \
    tree.speciesTreePlot(run_id,'virus',DBSIM_PATH,DBTREE_PATH,\
    treepalette,maxticksize,figxy,hratio,vabundthreshold)
    bAbunds,bKeepTreeStrainsDF, bSpeciesColorDict, _, _, figB, axesB = \
    tree.speciesTreePlot(run_id,'microbe',DBSIM_PATH,DBTREE_PATH,\
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
    bAbunds, _, bSpeciesColorDict, _, _, _, _  = \
    tree.speciesTreePlot(run_id,'microbe',DBSIM_PATH,DBTREE_PATH,\
    treepalette,maxticksize,figxy,hratio,babundthreshold)
    vAbunds, _, vSpeciesColorDict, _, _, _, _ = \
    tree.speciesTreePlot(run_id,'virus',DBSIM_PATH,DBTREE_PATH,\
    treepalette,maxticksize,figxy,hratio,vabundthreshold)
    bAbunds = bAbunds[bAbunds.t <= max(vAbunds.t)]
    bAbunds = bAbunds[bAbunds['abundance']>0]
    # This is to set abundances to log scale
    bAbunds['abundance'] = np.log(np.array(bAbunds['abundance']))
    vAbunds['abundance'] = np.log(np.array(vAbunds['abundance']))
    microbe_stacked = bAbunds.pivot(index='t',columns='tree_bstrain_id',values='abundance')
    virus_stacked = vAbunds.pivot(index='t',columns='tree_vstrain_id',values='abundance')
    fig, ax = plt.subplots(1)
    microbe_stacked.plot.area(ax = ax,stacked=True,legend=False, linewidth=0,color=bSpeciesColorDict,sort_columns=True)
    microbe_stacked.plot(stacked=True, ax=ax, legend=False, color='white',sort_columns=True,linewidth=.1)
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
    virus_stacked.plot(stacked=True, ax=ax, legend=False, color='white',sort_columns=True,linewidth=.1)
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
    microbe_stacked.plot(stacked=True, ax=axes[0], legend=False, color='white',sort_columns=True,linewidth=.1)
    virus_stacked.plot.area(ax = axes[1],stacked=True,legend=False, linewidth=0,color=vSpeciesColorDict,sort_columns=True)
    virus_stacked.plot(stacked=True, ax=axes[1], legend=False, color='white',sort_columns=True,linewidth=.1)
    fig.tight_layout()
    fig.savefig(os.path.join(PLOT_PATH,'microbe-virus-stacked-abundances.png'),dpi=resolve)

if sys.argv[2] == 'diversity':
    bAbunds, _, bSpeciesColorDict, _, _, _, _  = \
    tree.speciesTreePlot(run_id,'microbe',DBSIM_PATH,DBTREE_PATH,\
    treepalette,maxticksize,figxy,hratio,babundthreshold)
    vAbunds, _, vSpeciesColorDict, _, _, _, _ = \
    tree.speciesTreePlot(run_id,'virus',DBSIM_PATH,DBTREE_PATH,\
    treepalette,maxticksize,figxy,hratio,vabundthreshold)
    bAbunds = bAbunds[bAbunds.t <= max(vAbunds.t)]
    bAbunds = bAbunds[bAbunds['abundance']>0]
    # This is to set abundances to log scale
    # bAbunds['abundance'] = np.log(np.array(bAbunds['abundance']))
    # vAbunds['abundance'] = np.log(np.array(vAbunds['abundance']))
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
    microbe_stacked = bAbunds.pivot(index='t',columns='tree_bstrain_id',values='abundance')
    microbe_stacked.plot.area(ax = axes[0],stacked=True,legend=False, linewidth=0,color=bSpeciesColorDict,sort_columns=True)
    microbe_stacked.plot(stacked=True, ax=axes[0], legend=False, color='white',sort_columns=True,linewidth=.1)
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
    virus_stacked.plot(stacked=True, ax=axes[1], legend=False, color='white',sort_columns=True,linewidth=.1)
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
    microbe_stacked.plot(stacked=True, ax=axes[0], legend=False, color='white',sort_columns=True,linewidth=.1)
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
    ###
    ## This is clade diversity
    ###
    microbe_stacked = pd.read_sql_query("SELECT t,bstrain_id,abundance FROM babundance \
                        WHERE run_id = {}".format(run_id), conSim)
    microbe_stacked = microbe_stacked[microbe_stacked.abundance > 0]
    microbe_stacked = microbe_stacked[microbe_stacked.t<=max(virus_total.t)]
    bstrainclades = pd.read_sql_query("SELECT t,bstrain_id,lineage FROM MRCA_bstrains", conTree)
    bstrainclades = bstrainclades.merge(microbe_stacked,on=['t','bstrain_id'])
    bstrainclades = bstrainclades.groupby(['t','lineage'])\
                        .agg(lineageAbundance=('abundance','sum')).reset_index()
    btotal = bstrainclades.groupby(['t'])\
                        .agg(btotal=('lineageAbundance','sum')).reset_index()
    bstrainclades = bstrainclades.merge(btotal,on=['t'])
    bstrainclades['lineageFreq'] = np.array(bstrainclades['lineageAbundance'])\
                                    /np.array(bstrainclades['btotal'])
    shannonLineage = []
    for t in sorted(bstrainclades['t'].unique()):
        div = list(bstrainclades[(bstrainclades['t']==t)]['lineageFreq'])
        shannonLineage.append(np.exp(-1*sum(np.array(div)*np.array(list(map(np.log,div))))))
    richnessLineage = []
    for t in virus_total['t']:
        div = list(bstrainclades[(bstrainclades['t']==t)]['lineageFreq'])
        richnessLineage.append(len(div))
    fig, ax = plt.subplots(2,sharex=True)
    axes = [ax[0], ax[0].twinx(), ax[1], ax[1].twinx(), ax[1].twinx()]
    microbe_stacked = bAbunds.pivot(index='t',columns='tree_bstrain_id',values='abundance')
    microbe_stacked.plot.area(ax = axes[0],stacked=True,legend=False, linewidth=0,color=bSpeciesColorDict,sort_columns=True)
    microbe_stacked.plot(stacked=True, ax=axes[0], legend=False, color='white',sort_columns=True,linewidth=.1)
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
    axes[3].plot(virus_total['t'],richnessLineage,color='darkred',linewidth=1)
    axes[3].yaxis.tick_right()
    axes[3].set_ylabel(ylabel ='Host Clade Richness',labelpad=15,fontsize=7)
    axes[4].plot(virus_total['t'],shannonLineage,color='darkblue',linewidth=1)
    axes[4].yaxis.tick_left()
    axes[4].yaxis.set_label_position("left")
    axes[4].set_ylabel(ylabel ='Host Clade \nShannon Diversity',labelpad=15,fontsize=7)
    fig.tight_layout()
    fig.savefig(os.path.join(PLOT_PATH,'microbe-clade-richness-shannon-stacked.png'),dpi=resolve)
    plt.close(fig)




if sys.argv[2] == 'meanviralfitness':
    bAbunds, _, bSpeciesColorDict, _, _, _, _  = \
    tree.speciesTreePlot(run_id,'microbe',DBSIM_PATH,DBTREE_PATH,\
    treepalette,maxticksize,figxy,hratio,babundthreshold)
    vAbunds, _, vSpeciesColorDict, _, _, _, _ = \
    tree.speciesTreePlot(run_id,'virus',DBSIM_PATH,DBTREE_PATH,\
    treepalette,maxticksize,figxy,hratio,vabundthreshold)
    bAbunds = bAbunds[bAbunds.t <= max(vAbunds.t)]
    bAbunds = bAbunds[bAbunds['abundance']>0]
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
    microbe_stacked = bAbunds.pivot(index='t',columns='tree_bstrain_id',values='abundance')
    microbe_stacked.plot.area(ax = axes[0],stacked=True,legend=False, linewidth=0,color=bSpeciesColorDict,sort_columns=True)
    microbe_stacked.plot(stacked=True, ax=axes[0], legend=False, color='white',sort_columns=True,linewidth=.1)
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

if sys.argv[2] == 'meanintrinsicfitness':
    bAbunds, _, bSpeciesColorDict, _, _, _, _  = \
    tree.speciesTreePlot(run_id,'microbe',DBSIM_PATH,DBTREE_PATH,\
    treepalette,maxticksize,figxy,hratio,babundthreshold)
    vAbunds, _, vSpeciesColorDict, _, _, _, _ = \
    tree.speciesTreePlot(run_id,'virus',DBSIM_PATH,DBTREE_PATH,\
    treepalette,maxticksize,figxy,hratio,vabundthreshold)
    bAbunds = bAbunds[bAbunds.t <= max(vAbunds.t)]
    bAbunds = bAbunds[bAbunds['abundance']>0]
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
    microbe_stacked.plot(stacked=True, ax=axes[0], legend=False, color='white',sort_columns=True,linewidth=.1)
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

if sys.argv[2] == 'meanviraldeath':
    bAbunds, _, bSpeciesColorDict, _, _, _, _  = \
    tree.speciesTreePlot(run_id,'microbe',DBSIM_PATH,DBTREE_PATH,\
    treepalette,maxticksize,figxy,hratio,babundthreshold)
    vAbunds, _, vSpeciesColorDict, _, _, _, _ = \
    tree.speciesTreePlot(run_id,'virus',DBSIM_PATH,DBTREE_PATH,\
    treepalette,maxticksize,figxy,hratio,vabundthreshold)
    bAbunds = bAbunds[bAbunds.t <= max(vAbunds.t)]
    bAbunds = bAbunds[bAbunds['abundance']>0]
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
    singlematch = pd.read_sql_query("SELECT t, vstrain_id, bstrain_id\
        FROM bstrain_to_vstrain_matches WHERE match_length in (1)", conMatch)\
        .drop_duplicates()
    doublematch = pd.read_sql_query("SELECT t, vstrain_id, bstrain_id\
        FROM bstrain_to_vstrain_matches WHERE match_length in (2)", conMatch)\
        .drop_duplicates()
    triplematch = pd.read_sql_query("SELECT t, vstrain_id, bstrain_id\
        FROM bstrain_to_vstrain_matches WHERE match_length in (3)", conMatch)\
        .drop_duplicates()
    zeromatch.sort_values(by=['t', 'vstrain_id'])
    zeromatch = zeromatch.merge(vfreq,on=['t','vstrain_id'])\
                    .merge(bfreq,on=['t','bstrain_id'])
    zeromatch = zeromatch.groupby(['t','vstrain_id','vfreq'])\
                        .agg(bfreq=('bfreq','sum')).reset_index()
    zeromatch['meanfitness'] = zeromatch['vfreq']*zeromatch['bfreq']
    zeromatch = zeromatch.groupby(['t'])\
                        .agg(meanfitness=('meanfitness','sum')).reset_index()
    singlematch.sort_values(by=['t', 'vstrain_id'])
    singlematch = singlematch.merge(vfreq,on=['t','vstrain_id'])\
                    .merge(bfreq,on=['t','bstrain_id'])
    singlematch = singlematch.groupby(['t','vstrain_id','vfreq'])\
                        .agg(bfreq=('bfreq','sum')).reset_index()
    singlematch['meanfitness'] = singlematch['vfreq']*singlematch['bfreq']
    singlematch = singlematch.groupby(['t'])\
                        .agg(meanfitness=('meanfitness','sum')).reset_index()
    doublematch.sort_values(by=['t', 'vstrain_id'])
    doublematch = doublematch.merge(vfreq,on=['t','vstrain_id'])\
                    .merge(bfreq,on=['t','bstrain_id'])
    doublematch = doublematch.groupby(['t','vstrain_id','vfreq'])\
                        .agg(bfreq=('bfreq','sum')).reset_index()
    doublematch['meanfitness'] = doublematch['vfreq']*doublematch['bfreq']
    doublematch = doublematch.groupby(['t'])\
                        .agg(meanfitness=('meanfitness','sum')).reset_index()
    triplematch.sort_values(by=['t', 'vstrain_id'])
    triplematch = triplematch.merge(vfreq,on=['t','vstrain_id'])\
                    .merge(bfreq,on=['t','bstrain_id'])
    triplematch = triplematch.groupby(['t','vstrain_id','vfreq'])\
                        .agg(bfreq=('bfreq','sum')).reset_index()
    triplematch['meanfitness'] = triplematch['vfreq']*triplematch['bfreq']
    triplematch = triplematch.groupby(['t'])\
                        .agg(meanfitness=('meanfitness','sum')).reset_index()
    fig, ax = plt.subplots(2,sharex=True)
    axes = [ax[0], ax[0].twinx(), ax[1], ax[1].twinx()]
    microbe_stacked = bAbunds.pivot(index='t',columns='tree_bstrain_id',values='abundance')
    microbe_stacked.plot.area(ax = axes[0],stacked=True,legend=False, linewidth=0,color=bSpeciesColorDict,sort_columns=True)
    microbe_stacked.plot(stacked=True, ax=axes[0], legend=False, color='white',sort_columns=True,linewidth=.1)
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
    axes[3].plot(zeromatch['t'],zeromatch['meanfitness'],label='0 match',linewidth=1)
    axes[3].plot(singlematch['t'],singlematch['meanfitness'],label='1 match',linewidth=1)
    axes[3].plot(doublematch['t'],doublematch['meanfitness'],label='2 match',linewidth=1)
    axes[3].plot(triplematch['t'],triplematch['meanfitness'],label='3 match',linewidth=1)
    axes[3].yaxis.tick_left()
    axes[3].yaxis.set_label_position("left")
    axes[3].set_ylabel(ylabel ='Mean Fitness',labelpad=15,fontsize=7)
    axes[3].legend()
    fig.tight_layout()
    fig.savefig(os.path.join(PLOT_PATH,'mean-viral-death.png'),dpi=resolve)


if sys.argv[2] == 'MOI':
    bAbunds, _, bSpeciesColorDict, _, _, _, _  = \
    tree.speciesTreePlot(run_id,'microbe',DBSIM_PATH,DBTREE_PATH,\
    treepalette,maxticksize,figxy,hratio,babundthreshold)
    vAbunds, _, vSpeciesColorDict, _, _, _, _ = \
    tree.speciesTreePlot(run_id,'virus',DBSIM_PATH,DBTREE_PATH,\
    treepalette,maxticksize,figxy,hratio,vabundthreshold)
    bAbunds = bAbunds[bAbunds.t <= max(vAbunds.t)]
    bAbunds = bAbunds[bAbunds['abundance']>0]
    microbe_stacked = pd.read_sql_query("SELECT t,bstrain_id,abundance FROM babundance \
                        WHERE run_id = {}".format(run_id), conSim)
    microbe_stacked = microbe_stacked[microbe_stacked.t<=max(virus_total.t)]
    virus_stacked = pd.read_sql_query("SELECT t,vstrain_id,abundance FROM vabundance \
                        WHERE run_id = {}".format(run_id), conSim)
    bfreq = microbe_stacked.merge(microbe_total, on=['t'])
    bfreq = bfreq.rename(columns={"abundance": "babundance"})
    vfreq = virus_stacked.merge(virus_total, on=['t'])
    vfreq = vfreq.rename(columns={"abundance": "vabundance"})
    bfreq['bfreq'] = bfreq['babundance']/bfreq['Microbial Abundance']
    # vfreq['vfreq'] = vfreq['vabundance']/vfreq['Viral Abundance']
    zeromatch = pd.read_sql_query("SELECT t, vstrain_id, bstrain_id\
        FROM bstrain_to_vstrain_0matches", conMatch).drop_duplicates()
    zeromatch.sort_values(by=['t', 'bstrain_id'])
    zeromatch = zeromatch.merge(bfreq,on=['t','bstrain_id'])
    zeromatch = zeromatch.merge(vfreq,on=['t','vstrain_id'])
    zeromatch = zeromatch.groupby(['t','bstrain_id','bfreq','babundance'])\
                        .agg(vabundance=('vabundance','sum')).reset_index()
    zeromatch['MOI'] = np.array(zeromatch['vabundance'])/np.array(zeromatch['babundance'])
    zeromatch['expMOI'] = np.array(zeromatch['bfreq'])*np.array(zeromatch['MOI'])
    zeromatch = zeromatch.groupby(['t'])\
                        .agg(expMOI=('expMOI','sum')).reset_index()
    fig, ax = plt.subplots(2,sharex=True)
    axes = [ax[0], ax[0].twinx(), ax[1], ax[1].twinx()]
    microbe_stacked = bAbunds.pivot(index='t',columns='tree_bstrain_id',values='abundance')
    microbe_stacked.plot.area(ax = axes[0],stacked=True,legend=False, linewidth=0,color=bSpeciesColorDict,sort_columns=True)
    microbe_stacked.plot(stacked=True, ax=axes[0], legend=False, color='white',sort_columns=True,linewidth=.1)
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
    axes[3].plot(zeromatch['t'],zeromatch['expMOI'],linewidth=1)
    axes[3].yaxis.tick_left()
    axes[3].yaxis.set_label_position("left")
    axes[3].set_ylabel(ylabel ='Expected MOI',labelpad=15,fontsize=7)
    fig.tight_layout()
    fig.savefig(os.path.join(PLOT_PATH,'expected-MOI.png'),dpi=resolve)

if sys.argv[2] == 'MRCA':
    bAbunds, _, bSpeciesColorDict, _, _, _, _  = \
    tree.speciesTreePlot(run_id,'microbe',DBSIM_PATH,DBTREE_PATH,\
    treepalette,maxticksize,figxy,hratio,babundthreshold)
    vAbunds, _, vSpeciesColorDict, _, _, _, _ = \
    tree.speciesTreePlot(run_id,'virus',DBSIM_PATH,DBTREE_PATH,\
    treepalette,maxticksize,figxy,hratio,vabundthreshold)
    bAbunds = bAbunds[bAbunds.t <= max(vAbunds.t)]
    bAbunds = bAbunds[bAbunds['abundance']>0]
    microbe_stacked = pd.read_sql_query("SELECT t,bstrain_id,abundance FROM babundance \
                        WHERE run_id = {}".format(run_id), conSim)
    microbe_stacked = microbe_stacked[microbe_stacked.t<=max(virus_total.t)]
    virus_stacked = pd.read_sql_query("SELECT t,vstrain_id,abundance FROM vabundance \
                        WHERE run_id = {}".format(run_id), conSim)
    bfreq = microbe_stacked.merge(microbe_total, on=['t'])
    bfreq['bfreq'] = bfreq['abundance']/bfreq['Microbial Abundance']
    bfreq = bfreq.drop(columns=['abundance','Microbial Abundance'])
    vfreq = virus_stacked.merge(virus_total, on=['t'])
    vfreq['vfreq'] = vfreq['abundance']/vfreq['Viral Abundance']
    vfreq = vfreq.drop(columns=['abundance','Viral Abundance'])
    pairwiseMicrobe = pd.read_sql_query("SELECT t,bstrain_id_sample_1, bstrain_id_sample_2,\
                            lineage, MRCA_bstrain_id, t_creation \
                            FROM pairwise_MRCA_bstrains", conTree)
    bfreq = bfreq.rename(columns={'bstrain_id':'bstrain_id_sample_1'})
    pairwiseMicrobe = pairwiseMicrobe.merge(bfreq,on=['t','bstrain_id_sample_1'])
    pairwiseMicrobe = pairwiseMicrobe.rename(columns={'bfreq':'bfreq_sample_1'})
    bfreq = bfreq.rename(columns={'bstrain_id_sample_1':'bstrain_id_sample_2'})
    pairwiseMicrobe = pairwiseMicrobe.merge(bfreq,on=['t','bstrain_id_sample_2'])
    pairwiseMicrobe = pairwiseMicrobe.rename(columns={'bfreq':'bfreq_sample_2'})
    bfreq = bfreq.rename(columns={'bstrain_id_sample_2':'bstrain_id'})
    pairwiseMicrobe['prob'] = np.array(pairwiseMicrobe['bfreq_sample_1'])\
                                *np.array(pairwiseMicrobe['bfreq_sample_2'])
    norm = pairwiseMicrobe.groupby(['t','lineage'])\
                        .agg(norm=('prob','sum')).reset_index()
    pairwiseMicrobe = pairwiseMicrobe.merge(norm,on=['t','lineage'])
    pairwiseMicrobe['prob'] = np.array(pairwiseMicrobe['prob'])/np.array(pairwiseMicrobe['norm'])
    pairwiseMicrobe['t_creation_weighted'] = np.array(pairwiseMicrobe['prob'])\
                                                *np.array(pairwiseMicrobe['t_creation'])
    lineageMicrobe = pd.read_sql_query("SELECT t,bstrain_id,lineage \
                            FROM MRCA_bstrains", conTree)
    lineageMicrobe = lineageMicrobe.merge(bfreq,on=['t','bstrain_id'])
    lineageMicrobe = lineageMicrobe.groupby(['t','lineage'])\
                        .agg(lineage_bfreq=('bfreq','sum')).reset_index()
    pairwiseMicrobe = pairwiseMicrobe.merge(lineageMicrobe,on=['t','lineage'])
    pairwiseMicrobe['t_creation_weighted'] = np.array(pairwiseMicrobe['t_creation_weighted'])\
                                                *np.array(pairwiseMicrobe['lineage_bfreq'])
    expectedMicrobePairwiseMRCA = pairwiseMicrobe.groupby(['t'])\
                        .agg(exp_tCreation=('t_creation_weighted','sum')).reset_index()
    pairwiseVirus = pd.read_sql_query("SELECT t,vstrain_id_sample_1, vstrain_id_sample_2,\
                            lineage, MRCA_vstrain_id, t_creation \
                            FROM pairwise_MRCA_vstrains", conTree)
    vfreq = vfreq.rename(columns={'vstrain_id':'vstrain_id_sample_1'})
    pairwiseVirus = pairwiseVirus.merge(vfreq,on=['t','vstrain_id_sample_1'])
    pairwiseVirus = pairwiseVirus.rename(columns={'vfreq':'vfreq_sample_1'})
    vfreq = vfreq.rename(columns={'vstrain_id_sample_1':'vstrain_id_sample_2'})
    pairwiseVirus = pairwiseVirus.merge(vfreq,on=['t','vstrain_id_sample_2'])
    pairwiseVirus = pairwiseVirus.rename(columns={'vfreq':'vfreq_sample_2'})
    vfreq = vfreq.rename(columns={'vstrain_id_sample_2':'vstrain_id'})
    pairwiseVirus['prob'] = np.array(pairwiseVirus['vfreq_sample_1'])\
                                *np.array(pairwiseVirus['vfreq_sample_2'])
    norm = pairwiseVirus.groupby(['t','lineage'])\
                        .agg(norm=('prob','sum')).reset_index()
    pairwiseVirus = pairwiseVirus.merge(norm,on=['t','lineage'])
    pairwiseVirus['prob'] = np.array(pairwiseVirus['prob'])/np.array(pairwiseVirus['norm'])
    pairwiseVirus['t_creation_weighted'] = np.array(pairwiseVirus['prob'])\
                                            *np.array(pairwiseVirus['t_creation'])
    lineageVirus = pd.read_sql_query("SELECT t,vstrain_id,lineage \
                            FROM MRCA_vstrains", conTree)
    lineageVirus = lineageVirus.merge(vfreq,on=['t','vstrain_id'])
    lineageVirus = lineageVirus.groupby(['t','lineage'])\
                        .agg(lineage_vfreq=('vfreq','sum')).reset_index()
    pairwiseVirus = pairwiseVirus.merge(lineageVirus,on=['t','lineage'])
    pairwiseVirus['t_creation_weighted'] = np.array(pairwiseVirus['t_creation_weighted'])\
                                                *np.array(pairwiseVirus['lineage_vfreq'])
    expectedVirusPairwiseMRCA = pairwiseVirus.groupby(['t'])\
                    .agg(exp_tCreation=('t_creation_weighted','sum')).reset_index()
    fig, ax = plt.subplots(3,sharex=True)
    axes = [ax[0], ax[0].twinx(), ax[1], ax[1].twinx(), ax[2], ax[2].twinx()]
    microbe_stacked = bAbunds[bAbunds.t<=max(virus_total['t'])].pivot(index='t',columns='tree_bstrain_id',values='abundance')
    microbe_stacked.plot.area(ax = axes[0],stacked=True,legend=False, linewidth=0,color=bSpeciesColorDict,sort_columns=True)
    microbe_stacked.plot(stacked=True, ax=axes[0], legend=False, color='white',sort_columns=True,linewidth=.1)
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
    axes[3].plot(expectedVirusPairwiseMRCA['t'],\
                list(np.array(expectedVirusPairwiseMRCA['exp_tCreation'])),\
                    color = 'darkred',label='Virus',linewidth=1)
    axes[3].plot(expectedMicrobePairwiseMRCA['t'],\
                list(np.array(expectedMicrobePairwiseMRCA['exp_tCreation'])),\
                    color = 'darkblue',label='Host',linewidth=1)
    # axes[3].plot(expectedVirusPairwiseMRCA['t'],expectedVirusPairwiseMRCA['exp_tCreation'],\
    #                 color = 'darkred', label='Virus',linewidth=1)
    axes[3].yaxis.tick_left()
    axes[3].yaxis.set_label_position("left")
    axes[3].set_ylabel(ylabel ='Expected Time of MRCA',labelpad=15,fontsize=7)
    axes[3].legend()
    axes[3].set_xlabel(xlabel = 'Time t',fontsize=7)
    axes[4].fill_between(virus_total['t'],virus_total['Viral Abundance'], color='grey',alpha=0.6)
    lim = axes[2].get_ylim()
    axes[4].set_ylim(0,lim[1])
    axes[4].set_yticklabels([])
    axes[4].set_yticks([])
    axes[4].margins(x=0)
    axes[5].margins(x=0)
    axes[5].plot(expectedVirusPairwiseMRCA['t'],\
                list(np.array(expectedVirusPairwiseMRCA['t'])-np.array(expectedVirusPairwiseMRCA['exp_tCreation'])),\
                    color = 'darkred',label='Virus',linewidth=1)
    axes[5].plot(expectedMicrobePairwiseMRCA['t'],\
                list(np.array(expectedMicrobePairwiseMRCA['t'])-np.array(expectedMicrobePairwiseMRCA['exp_tCreation'])),\
                    color = 'darkblue',label='Host',linewidth=1)
    # axes[3].plot(expectedVirusPairwiseMRCA['t'],expectedVirusPairwiseMRCA['exp_tCreation'],\
    #                 color = 'darkred', label='Virus',linewidth=1)
    axes[5].yaxis.tick_left()
    axes[5].yaxis.set_label_position("left")
    axes[5].set_ylabel(ylabel ='Expected Time to MRCA ',labelpad=15,fontsize=7)
    axes[5].legend()
    axes[5].set_xlabel(xlabel = 'Time t',fontsize=7)
    fig.tight_layout()
    fig.savefig(os.path.join(PLOT_PATH,'expected-pairwise-mrca.png'),dpi=resolve)


if sys.argv[2] == 'viraltemporaladaptation':
    bAbunds, _, bSpeciesColorDict, _, _, _, _  = \
    tree.speciesTreePlot(run_id,'microbe',DBSIM_PATH,DBTREE_PATH,\
    treepalette,maxticksize,figxy,hratio,babundthreshold)
    vAbunds, _, vSpeciesColorDict, _, _, _, _ = \
    tree.speciesTreePlot(run_id,'virus',DBSIM_PATH,DBTREE_PATH,\
    treepalette,maxticksize,figxy,hratio,vabundthreshold)
    bAbunds = bAbunds[bAbunds.t <= max(vAbunds.t)]
    bAbunds = bAbunds[bAbunds['abundance']>0]
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
    singlematch = pd.read_sql_query("SELECT t, vstrain_id, bstrain_id\
        FROM bstrain_to_vstrain_matches WHERE match_length in (1)", conMatch)\
        .drop_duplicates()
    doublematch = pd.read_sql_query("SELECT t, vstrain_id, bstrain_id\
        FROM bstrain_to_vstrain_matches WHERE match_length in (2)", conMatch)\
        .drop_duplicates()
    triplematch = pd.read_sql_query("SELECT t, vstrain_id, bstrain_id\
        FROM bstrain_to_vstrain_matches WHERE match_length in (3)", conMatch)\
        .drop_duplicates()
    zeromatch.sort_values(by=['t', 'vstrain_id'])
    zeromatch = zeromatch.merge(vfreq,on=['t','vstrain_id'])\
                    .merge(bfreq,on=['t','bstrain_id'])
    zeromatch = zeromatch.groupby(['t','vstrain_id','vfreq'])\
                        .agg(bfreq=('bfreq','sum')).reset_index()
    zeromatch['meanfitness'] = zeromatch['vfreq']*zeromatch['bfreq']
    zeromatch = zeromatch.groupby(['t'])\
                        .agg(meanfitness=('meanfitness','sum')).reset_index()
    singlematch.sort_values(by=['t', 'vstrain_id'])
    singlematch = singlematch.merge(vfreq,on=['t','vstrain_id'])\
                    .merge(bfreq,on=['t','bstrain_id'])
    singlematch = singlematch.groupby(['t','vstrain_id','vfreq'])\
                        .agg(bfreq=('bfreq','sum')).reset_index()
    singlematch['meanfitness'] = singlematch['vfreq']*singlematch['bfreq']
    singlematch = singlematch.groupby(['t'])\
                        .agg(meanfitness=('meanfitness','sum')).reset_index()
    doublematch.sort_values(by=['t', 'vstrain_id'])
    doublematch = doublematch.merge(vfreq,on=['t','vstrain_id'])\
                    .merge(bfreq,on=['t','bstrain_id'])
    doublematch = doublematch.groupby(['t','vstrain_id','vfreq'])\
                        .agg(bfreq=('bfreq','sum')).reset_index()
    doublematch['meanfitness'] = doublematch['vfreq']*doublematch['bfreq']
    doublematch = doublematch.groupby(['t'])\
                        .agg(meanfitness=('meanfitness','sum')).reset_index()
    triplematch.sort_values(by=['t', 'vstrain_id'])
    triplematch = triplematch.merge(vfreq,on=['t','vstrain_id'])\
                    .merge(bfreq,on=['t','bstrain_id'])
    triplematch = triplematch.groupby(['t','vstrain_id','vfreq'])\
                        .agg(bfreq=('bfreq','sum')).reset_index()
    triplematch['meanfitness'] = triplematch['vfreq']*triplematch['bfreq']
    triplematch = triplematch.groupby(['t'])\
                        .agg(meanfitness=('meanfitness','sum')).reset_index()
    fig, ax = plt.subplots(2,sharex=True)
    axes = [ax[0], ax[0].twinx(), ax[1], ax[1].twinx()]
    microbe_stacked = bAbunds.pivot(index='t',columns='tree_bstrain_id',values='abundance')
    microbe_stacked.plot.area(ax = axes[0],stacked=True,legend=False, linewidth=0,color=bSpeciesColorDict,sort_columns=True)
    microbe_stacked.plot(stacked=True, ax=axes[0], legend=False, color='white',sort_columns=True,linewidth=.1)
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
    axes[3].plot(zeromatch['t'],zeromatch['meanfitness'],label='0 match',linewidth=1)
    axes[3].plot(singlematch['t'],singlematch['meanfitness'],label='1 match',linewidth=1)
    axes[3].plot(doublematch['t'],doublematch['meanfitness'],label='2 match',linewidth=1)
    axes[3].plot(triplematch['t'],triplematch['meanfitness'],label='3 match',linewidth=1)
    axes[3].yaxis.tick_left()
    axes[3].yaxis.set_label_position("left")
    axes[3].set_ylabel(ylabel ='Mean Fitness',labelpad=15,fontsize=7)
    axes[3].legend()
    fig.tight_layout()
    fig.savefig(os.path.join(PLOT_PATH,'viral-temporal-adaptation.png'),dpi=resolve)


print('Complete!')
