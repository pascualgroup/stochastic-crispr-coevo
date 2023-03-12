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
import tripartitegraph as tg
import treefunctions as tree
# import mpld3

run_id = sys.argv[1]
print('Compiling plots and networks for run ID: {0}'.format(sys.argv[1]))
times = [float(sys.argv[i]) for i in list(range(2,len(sys.argv)))]
print('...for times: {0}...'.format(', '.join(map(str,times))))
# threshold = float(sys.argv[-1])/100
print('...and with with {0}% abundances removed...'.format(float(sys.argv[-1])))
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
vabundthreshold = 0.05
hyperAnalyze = True
html = False
overlay = True
sSpacing = 2
vSpacing = 2
bSpacing = 2

SCRIPT_PATH = os.path.abspath(os.path.dirname(__file__)) # cluster
dir = 'crispr-sweep-7-2-2022/isolates/runID3297-c66-r47'
run = 'runID3297-c66-r47'

DBSIM_PATH = os.path.join(SCRIPT_PATH,'..','..','..','simulation','sweep_db_gathered.sqlite')
# DBSIM_PATH = os.path.join('/Volumes','Yadgah',dir,'{}.sqlite'.format(run))

conSim = sqlite3.connect(DBSIM_PATH)
curSim = conSim.cursor()
ID = curSim.execute('SELECT combo_id,replicate FROM runs WHERE run_id = {}'.format(run_id)).fetchall()
combo_id = ID[0][0]
replicate = ID[0][1]
RUN_DIR = os.path.join('runID{0}-c{1}-r{2}'.format(run_id,combo_id,replicate))

DBMATCH_PATH = os.path.join(SCRIPT_PATH,'..','..','isolates',RUN_DIR,'matches_output.sqlite') # cluster
DBCLADE_PATH = os.path.join(SCRIPT_PATH,'..','..','isolates',RUN_DIR,'clade-abundances_output.sqlite') # cluster
DBTREE_PATH = os.path.join(SCRIPT_PATH,'..','..','isolates',RUN_DIR,'trees_output.sqlite') # cluster
DBTRI_PATH = os.path.join(SCRIPT_PATH,'..','..','isolates',RUN_DIR,'tripartite-networks_output.sqlite') # cluster
PLOT_PATH = os.path.join(SCRIPT_PATH,'..','..','isolates',RUN_DIR,'tripartite-images')
# DBMATCH_PATH = os.path.join('/Volumes','Yadgah',dir,'matches_output.sqlite') # local
# DBCLADE_PATH = os.path.join('/Volumes','Yadgah',dir,'clade-abundances_output.sqlite') # local
# DBTREE_PATH = os.path.join('/Volumes','Yadgah',dir,'trees_output.sqlite') # local
# DBTRI_PATH = os.path.join('/Volumes','Yadgah',dir,'tripartite-networks_output.sqlite') # local
# PLOT_PATH = os.path.join('/Volumes','Yadgah','tripartite-images')
# PLOT_PATH = os.path.join('/Volumes','Yadgah','tripartite-check')
# PLOT_PATH = os.path.join('/Volumes','Yadgah')

if not os.path.isdir(PLOT_PATH):
    os.makedirs (PLOT_PATH)

conMatch = sqlite3.connect(DBMATCH_PATH)
curMatch = conMatch.cursor()
conTree = sqlite3.connect(DBTREE_PATH)
curTree = conTree.cursor()
conTri = sqlite3.connect(DBTRI_PATH)
curTri = conTri.cursor()

# print('SQLite Query: microbial abundance time series data')
# microbe_stacked = pd.read_sql_query("SELECT t,bstrain_id,abundance FROM babundance WHERE run_id = {}".format(run_id), conSim)
# microbe_stacked = microbe_stacked.pivot(index='t',columns='bstrain_id',values='abundance')
microbe_total = pd.read_sql_query("SELECT t,microbial_abundance FROM summary WHERE run_id = {}".format(run_id), conSim)\
.rename(columns={"microbial_abundance": "Microbial Abundance"})
# print('SQLite Query: viral abundance time series data')
# virus_stacked = pd.read_sql_query("SELECT t,vstrain_id,abundance FROM vabundance WHERE run_id = {}".format(run_id), conSim)
virus_total = pd.read_sql_query("SELECT t,viral_abundance FROM summary WHERE run_id = {}".format(run_id), conSim)\
.rename(columns={"viral_abundance": "Viral Abundance"})
# virus_stacked = virus_stacked.pivot(index='t',columns='vstrain_id',values='abundance')
# if len(microbe_stacked.index) > len(virus_stacked.index):
#     microbe_stacked.drop(index=list(set(microbe_stacked.index)-set(virus_stacked.index)),inplace=True)


virusColorDict, vcladeColorDict, vcladeDict, treevallstrains,figV, axesV = \
tree.virusCladeTreePlot(run_id,DBSIM_PATH,DBCLADE_PATH,DBTREE_PATH,\
colorpalClades,maxAbundTickSize,cladetreefigxy,cladetreehratio,\
lifeTimeThreshold,vabundthreshold)
microbeColorDict, mcladeColorDict, mcladeDict, figB, axesB = \
tree.microbeCladeTreePlot(run_id,DBSIM_PATH,DBCLADE_PATH,DBTREE_PATH,\
colorpalClades,maxAbundTickSize,cladetreefigxy,cladetreehratio)
figS, axesS = \
tree.susceptibleCladeTreePlot(run_id,DBSIM_PATH,DBCLADE_PATH,DBTREE_PATH,\
DBMATCH_PATH,colorpalClades,maxAbundTickSize,cladetreefigxy,cladetreehratio)

# if justTree == 1:
#     print('Compiling just trees...')
#     for imgType in treeImgTypes:
#         if overlay:
#             axesB.append(axesB[0].twinx())
#             axesB[2].plot(virus_total['t'],virus_total['Viral Abundance'],linewidth=0,color='grey')
#             axesB[2].fill_between(virus_total['t'],virus_total['Viral Abundance'], color='grey',alpha=0.6)
#             lim = axesB[2].get_ylim()
#             axesB[2].set_ylim(0,lim[1])
#             figB.tight_layout()
#             figB.savefig(os.path.join(PLOT_PATH,'microbe-clades-abundances-tree-time-slices_virus-overlay.{0}'\
#                 .format(imgType)),dpi=resolve)
#             axesV.append(axesV[0].twinx())
#             axesV[2].plot(microbe_total['t'],microbe_total['Microbial Abundance'],linewidth=0,color='grey')
#             axesV[2].fill_between(microbe_total['t'],microbe_total['Microbial Abundance'], color='grey',alpha=0.4)
#             lim = axesV[2].get_ylim()
#             axesV[2].set_ylim(0,lim[1])
#             figV.tight_layout()
#             figV.savefig(os.path.join(PLOT_PATH,'virus-clades-abundances-tree-time-slices_microbe-overlay.{0}'\
#                 .format(imgType)),dpi=resolve)
#         else:
#             figB.tight_layout()
#             figB.savefig(os.path.join(PLOT_PATH,'microbe-clades-abundances-tree-time-slices.{0}'\
#                 .format(imgType)),dpi=resolve)
#             figV.tight_layout()
#             figV.savefig(os.path.join(PLOT_PATH,'virus-clades-abundances-tree-time-slices.{0}'\
#                 .format(imgType)),dpi=resolve)
#         figS.tight_layout()
#         figS.savefig(os.path.join(PLOT_PATH,'susceptible-clades-abundances-tree-time-slices.{0}'\
#             .format(imgType)),dpi=resolve)


for time in times:
    print('Compiling tripartite and infection networks at time = {}'.format(time))
    tg.tripartiteGraph(time,run_id,sSpacing,vSpacing,bSpacing,html,\
                        colorpalSpacers,graphImgTypes,\
                        DBSIM_PATH,DBMATCH_PATH,DBTRI_PATH,PLOT_PATH)
    vTotal = virus_total[virus_total.t==time]['Viral Abundance'].values[0]
    bTotal = microbe_total[microbe_total.t==time]['Microbial Abundance'].values[0]
    tripartiteNetwork = pd.read_sql_query(
        "SELECT vmatch_id, spacer_id, bmatch_id \
        FROM single_match_tripartite_networks WHERE t = {0}".format(time),
        conTri)
    bstrainmatchIDs = pd.read_sql_query(
        "SELECT match_id, bstrain_id \
        FROM bmatches WHERE t = {0}".format(time),
        conTri).rename(columns={"match_id": "bmatch_id"})
    bstrainAbund = pd.read_sql_query("SELECT bstrain_id,abundance FROM babundance  \
        WHERE run_id = {0} AND t = {1}".format(run_id,time), conSim).rename(columns={"abundance": "bstrainAbund"})
    bmatchAbund = pd.read_sql_query("SELECT match_id,babundance \
            FROM bmatches_abundances  \
            WHERE t = {0} AND match_id in ({1})".format(time,', '
            .join(map(str,tripartiteNetwork['bmatch_id'].unique()))), conTri)\
            .rename(columns={"match_id": "bmatch_id"}).rename(columns={"babundance": "bmatchAbund"})
    tripartiteNetwork = tripartiteNetwork.merge(bstrainmatchIDs,on=['bmatch_id']).drop(columns='bmatch_id')
    # bstrainmatchIDs = bstrainmatchIDs.merge(bstrainAbund,on=['bstrain_id'])\
    #                 .merge(bmatchAbund,on=['bmatch_id']).sort_values(by=['bmatchAbund', 'bstrainAbund'])
    vAbunds = pd.read_sql_query("SELECT vstrain_id,abundance FROM vabundance \
    WHERE run_id = {0} AND t = {1}".format(run_id,time), conSim)
    bmatch0vstrain = bstrainmatchIDs[:]
    bstrainmatchIDs = bstrainmatchIDs.merge(bstrainAbund,on=['bstrain_id'])\
                .merge(bmatchAbund,on=['bmatch_id']).sort_values(by=['bstrain_id']) # bin values are sorted here (but actually one below...)
    bstrain0vstrain = pd.read_sql_query(
    "SELECT bstrain_id, vstrain_id \
    FROM bstrain_to_vstrain_0matches WHERE t = {0}".format(time), conMatch)
    bmatch0vstrain = bstrain0vstrain.merge(bmatch0vstrain,on=['bstrain_id'])\
                        .drop(columns='bstrain_id').drop_duplicates()\
                        .merge(vAbunds,on=['vstrain_id']).groupby(['bmatch_id'])\
                        .agg(MatchCollapseProb=('abundance','sum')).reset_index()
    bmatch0vstrain['MatchCollapseProb'] = 1/vTotal*bmatch0vstrain['MatchCollapseProb']
    bstrain0vstrain = bstrain0vstrain.merge(vAbunds,on=['vstrain_id'])\
    .groupby(['bstrain_id']).agg(CollapseProb=('abundance','sum')).reset_index()
    bstrain0vstrain['CollapseProb'] = 1/vTotal*bstrain0vstrain['CollapseProb']
    bstrain0vstrain = bstrainmatchIDs.merge(bstrain0vstrain,on=['bstrain_id'])\
                        .merge(bmatch0vstrain,on=['bmatch_id'])
    bstrain0vstrain['CollapseProb'] = list(np.array(bstrain0vstrain['CollapseProb'])\
                                    *1/bTotal*np.array(bstrain0vstrain['bstrainAbund']))
    bstrain0vstrain['MatchCollapseProb'] = list(np.array(bstrain0vstrain['MatchCollapseProb'])\
                                    *1/bTotal*np.array(bstrain0vstrain['bmatchAbund']))
    bstrain0vstrain = bstrain0vstrain.sort_values(by=['MatchCollapseProb', 'CollapseProb']) # bin values are sorte
    bstrainsOrdered = list(set(bstrainmatchIDs['bstrain_id']) - set(bstrain0vstrain['bstrain_id']))
    bstrainsOrdered = list(bstrainmatchIDs[bstrainmatchIDs.bstrain_id.isin(bstrainsOrdered)]\
                    .sort_values(by=['bmatchAbund', 'bstrainAbund'])['bstrain_id']) # bin values are sorted here
    bstrainsOrdered.extend(bstrain0vstrain['bstrain_id'])
    # bstrainmatchIDs = bstrainmatchIDs.sort_values(by=['MatchCollapseProb', 'CollapseProb']) # bin values are sorted here
    vmatchAbund = pd.read_sql_query("SELECT match_id,vabundance, bsusceptible \
            FROM vmatches_abundances  \
            WHERE t = {0} AND match_id in ({1}) \
            ORDER BY vabundance".format(time,', '
            .join(map(str,tripartiteNetwork['vmatch_id'].unique()))), conTri)\
            .rename(columns={"match_id": "vmatch_id"})\
            .rename(columns={"vabundance": "vmatchAbund",\
            "bsusceptible": "EscEventProb"})
    vmatchAbund['EscEventProb']= list(1/bTotal * np.array(vmatchAbund['EscEventProb'])\
    *1/vTotal* np.array(vmatchAbund['vmatchAbund']))
    vmatchAbund = vmatchAbund.sort_values(by=['EscEventProb']) # bin values are sorted here
    spacerIDs = list(tripartiteNetwork['spacer_id'].unique())
    spacerIDs.sort()
    spacerDict = {i:j for (i,j) in zip(spacerIDs,list(range(0,len(spacerIDs))))}
    print(spacerDict)
    tripartiteNetwork['spacer_id'] = [spacerDict[i] for i in tripartiteNetwork['spacer_id']]
    Mcmap = cm.get_cmap('{}'.format(colorpalSpacers))
    spacerCMAP = mc.LinearSegmentedColormap.from_list(
        'spacers', [Mcmap(i) for i in range(len(spacerIDs))], len(spacerIDs))
    tripartiteNetwork = tripartiteNetwork.pivot(index='vmatch_id',columns='bstrain_id',values='spacer_id')
    tripartiteNetwork = tripartiteNetwork.reindex(index = vmatchAbund['vmatch_id']) # reordered here
    tripartiteNetwork = tripartiteNetwork[bstrainsOrdered] # reordered here
    microbeBounds = [0]
    mtickLocation = []
    for i in bstrainsOrdered:
        j = bstrainmatchIDs[bstrainmatchIDs.bstrain_id == i]['bstrainAbund'].values[0]
        microbeBounds.append(microbeBounds[-1]+np.log(j)+1)
        mtickLocation.append(microbeBounds[-2]+(microbeBounds[-1]-microbeBounds[-2])/2)
    mtickLocation = list(1/max(microbeBounds)*np.array(mtickLocation))
    microbeBounds = list(1/max(microbeBounds)*np.array(microbeBounds))
    virusBounds = [0]
    vtickLocation = []
    for i in vmatchAbund['vmatchAbund']:
        # virusBounds.append(virusBounds[-1]+np.log(i)+1)
        print('for vmatch_id {}'.format(i))
        virusBounds.append(virusBounds[-1]+np.log(i)+1)
        print('virus bounds are {}'.format(virusBounds))
        vtickLocation.append(virusBounds[-2]+(virusBounds[-1]-virusBounds[-2])/2)
        print('virus tick locations are {}'.format(vtickLocation))
    vtickLocation = list(1/max(virusBounds)*np.array(vtickLocation))
    virusBounds = list(1/max(virusBounds)*np.array(virusBounds))
    print('final virus bounds are {}'.format(virusBounds))
    print('final virus tick locations are {}'.format(vtickLocation))
    fig, ax = plt.subplots(1)
    heat = ax.pcolormesh(microbeBounds, virusBounds, tripartiteNetwork,cmap=spacerCMAP,edgecolors='black',linewidth=0.01)
    ax.set_xticks(mtickLocation)
    ax.set_xticklabels(bstrainsOrdered, rotation=0,fontsize = 4)
    ax.set_yticks(vtickLocation)
    ax.set_yticklabels(vmatchAbund['vmatch_id'], rotation=0,fontsize = 3.5)
    cbticks = list(
    np.linspace(.5*(len(spacerIDs)-1)/len(spacerIDs),(len(spacerIDs)-0.5)
    *(len(spacerIDs)-1)/len(spacerIDs),len(spacerIDs)
    ))
    cbar = fig.colorbar(heat, ticks=cbticks)
    cbar.ax.set_yticklabels(spacerIDs)
    cbar.set_label("Spacer ID")
    ax.set_xlabel(xlabel ='Microbe Strain ID',labelpad=10)
    ax.set_ylabel(ylabel = 'Virus Match ID',labelpad=10)
    ax.set_title('runID{0}-c{1}-r{2}\nTripartite Match Network at t = {3} (contact-ordered)'.format(run_id,combo_id,replicate,time),pad=20,fontsize=10)
    fig.tight_layout()
    for imgType in imgTypes:
        fig.savefig(os.path.join(PLOT_PATH,'tripartite-match-network-contact-ordered-time{0}.{1}'.format(int(np.floor(time)),imgType)),dpi=resolve)
    plt.close(fig)
    vstrainmatchIDs = pd.read_sql_query(
        "SELECT match_id, vstrain_id \
        FROM vmatches WHERE t = {0}".format(time),
        conTri).rename(columns={"match_id": "vmatch_id"})
    nonMatchNetwork = pd.read_sql_query(
        "SELECT vstrain_id, bstrain_id \
        FROM bstrain_to_vstrain_0matches WHERE t = {0} \
        AND vstrain_id in ({1})".format(time,
        ', '.join(map(str,np.unique(vstrainmatchIDs['vstrain_id'])))), conMatch)\
        .merge(vstrainmatchIDs,on=['vstrain_id']).drop(columns=['vstrain_id']).drop_duplicates()
    bstrainAbund = pd.read_sql_query("SELECT bstrain_id,abundance FROM babundance  \
        WHERE run_id = {0} AND t = {1}".format(run_id,time), conSim)\
        .rename(columns={"abundance": "bstrainAbund"}).merge(nonMatchNetwork,on=['bstrain_id'])\
        .drop(columns=['vmatch_id']).drop_duplicates().sort_values(by=['bstrainAbund']) # bin values are sorted here
    bstrain0vstrain = pd.read_sql_query(
    "SELECT bstrain_id, vstrain_id \
    FROM bstrain_to_vstrain_0matches WHERE t = {0}".format(time), conMatch)
    bstrain0vstrain = bstrain0vstrain.merge(vAbunds,on=['vstrain_id'])\
    .groupby(['bstrain_id']).agg(CollapseProb=('abundance','sum')).reset_index()
    bstrain0vstrain['CollapseProb'] = 1/vTotal*bstrain0vstrain['CollapseProb']
    bstrain0vstrain = bstrainAbund.merge(bstrain0vstrain,on=['bstrain_id'])
    bstrain0vstrain['CollapseProb'] = list(np.array(bstrain0vstrain['CollapseProb'])\
                                    *1/bTotal*np.array(bstrain0vstrain['bstrainAbund']))
    bstrain0vstrain = bstrain0vstrain.sort_values(by=['CollapseProb']) # bin values are sorted here
    nonMatchNetwork['presence'] = nonMatchNetwork['vmatch_id'].size*[1]
    nonMatchNetwork = nonMatchNetwork.pivot(index='vmatch_id',columns='bstrain_id',values='presence')
    nonMatchNetwork= nonMatchNetwork.reindex(index = vmatchAbund['vmatch_id'])
    nonMatchNetwork = nonMatchNetwork[bstrain0vstrain['bstrain_id']]
    microbeBounds = [0]
    mtickLocation = []
    for i in bstrain0vstrain['bstrainAbund']:
        microbeBounds.append(microbeBounds[-1]+np.log(i)+1)
        mtickLocation.append(microbeBounds[-2]+(microbeBounds[-1]-microbeBounds[-2])/2)
    mtickLocation = list(1/max(microbeBounds)*np.array(mtickLocation))
    microbeBounds = list(1/max(microbeBounds)*np.array(microbeBounds))
    fig, ax = plt.subplots(1)
    ax.pcolormesh(microbeBounds, virusBounds, nonMatchNetwork, edgecolors='black',linewidth=0.01)
    ax.set_xticks(mtickLocation)
    ax.set_xticklabels(bstrain0vstrain['bstrain_id'], rotation=0,fontsize = 4)
    ax.set_yticks(vtickLocation)
    ax.set_yticklabels(vmatchAbund['vmatch_id'], rotation=0,fontsize = 3.5)
    ax.set_xlabel(xlabel ='Microbe Strain ID',labelpad=10)
    ax.set_ylabel(ylabel = 'Virus Match ID',labelpad=10)
    ax.set_title('runID{0}-c{1}-r{2}\nInfection Network at t = {3} (contact-ordered)'.format(run_id,combo_id,replicate,time),pad=20,fontsize=10)
    fig.tight_layout()
    for imgType in imgTypes:
        fig.savefig(os.path.join(PLOT_PATH,'infection-network-contact-ordered-time{0}.{1}'.format(int(np.floor(time)),imgType)),dpi=resolve)
    plt.close(fig)
    vstrains =  pd.read_sql_query("SELECT vstrain_id FROM vabundance  \
            WHERE run_id = {0} AND t = {1}".format(run_id,time), conSim)
    vstrains = vstrains[vstrains.vstrain_id != 1]
    treevstrains = {tree:actual for (tree,actual) in
    zip([vstrain_id[0]
            for vstrain_id in
                curTree.execute("SELECT tree_vstrain_id FROM tree_vstrain_order \
                WHERE vstrain_id in ({}) \
                ORDER BY vstrain_id".format(', '.join(map(str,vstrains['vstrain_id']))))],sorted(vstrains['vstrain_id'])
                )\
        if tree in treevallstrains
        }
    axesV[0].axvline(x=time, color='black', ls=':', lw=.7)
    axesV[1].axvline(x=time, color='black', ls=':', lw=.7)
    if hyperAnalyze:
        # texts = []
        for treestrain in treevstrains.keys():
            pID = [id[0] for id in curSim.execute("SELECT parent_vstrain_id FROM vstrains \
                WHERE run_id = {0} AND vstrain_id = {1}"\
                .format(run_id,treevstrains[treestrain]))][0]
            pSpacers = [id[0] for id in curSim.execute("SELECT spacer_id FROM vpspacers \
                        WHERE run_id = {0} \
                        AND vstrain_id = {1}".format(run_id,pID))]
            dSpacers = [id[0] for id in curSim.execute("SELECT spacer_id FROM vpspacers \
                        WHERE run_id = {0} \
                        AND vstrain_id = {1}".format(run_id,treevstrains[treestrain]))]
            frm = list(set(pSpacers) - set(dSpacers))
            to = list(set(dSpacers) - set(pSpacers))
            bID = [id[0] for id in curSim.execute("SELECT infected_bstrain_id \
                FROM vstrains \
                WHERE run_id = {0}\
                AND vstrain_id = {1}".format(run_id,treevstrains[treestrain]))][0]
            mID = vstrainmatchIDs[vstrainmatchIDs.vstrain_id\
            == treevstrains[treestrain]]['vmatch_id'].values[0]
            # axesV[1].text(time, treestrain,\
            # 'mID:{0}, vID:{1}, bID:{2}::[{3}]->[{4}]'\
            # .format(mID,treevstrains[treestrain],bID,', '.join(map(str,frm)),\
            # ', '.join(map(str,to))),\
            # fontsize=3)#, fontweight='bold')
            axesV[1].text(time, treestrain,\
            '{0},{1},{2}::[{3}]->[{4}]'\
            .format(mID,treevstrains[treestrain],bID,', '.join(map(str,frm)),\
            ', '.join(map(str,to))),\
            fontsize=3)#, fontweight='bold')
            # texts.append('mID:{0}, vID:{1}, bID:{2}::[{3}]->[{4}]'\
            # .format(mID,treevstrains[treestrain],bID,', '.join(map(str,frm)),\
            # ', '.join(map(str,to))))
        # points = axesV[1].plot([time]*len(treevstrains.keys()),treevstrains.keys(),alpha=0)
        # text = mpld3.plugins.PointHTMLTooltip(points[0], texts)
        # mpld3.plugins.connect(figV,text)
    else:
        for treestrain in treevstrains.keys():
            bID = [id[0] for id in curSim.execute("SELECT infected_bstrain_id FROM vstrains \
                WHERE run_id = {0}\
                AND vstrain_id = {1}".format(run_id,treevstrains[treestrain]))][0]
            mID = vstrainmatchIDs[vstrainmatchIDs.vstrain_id\
            == treevstrains[treestrain]]['vmatch_id'].values[0]
            axesV[1].text(time, treestrain,\
            '{1}::mID:{0}'.format(mID,treevstrains[treestrain]),\
            fontsize=3)#, fontweight='bold')
    bstrains =  pd.read_sql_query("SELECT bstrain_id FROM babundance  \
            WHERE run_id = {0} AND t = {1}".format(run_id,time), conSim)
    bstrains = bstrains[bstrains.bstrain_id != 1]
    treebstrains = {tree:actual for (tree,actual) in
    zip([bstrain_id[0]
            for bstrain_id in
                curTree.execute("SELECT tree_bstrain_id FROM tree_bstrain_order \
                WHERE bstrain_id in ({}) \
                ORDER BY bstrain_id".format(', '\
                .join(map(str,bstrains['bstrain_id']))))],\
                sorted(bstrains['bstrain_id'])
                )
        }
    axesB[0].axvline(x=time, color='black', ls=':', lw=.7)
    axesB[1].axvline(x=time, color='black', ls=':', lw=.7)
    axesS[0].axvline(x=time, color='black', ls=':', lw=.7)
    axesS[1].axvline(x=time, color='black', ls=':', lw=.7)
    if hyperAnalyze:
        for treestrain in treebstrains.keys():
            vID = [id[0] for id in curSim.execute("SELECT infecting_vstrain_id FROM bstrains \
                WHERE run_id = {0}\
                AND bstrain_id = {1}".format(run_id,treebstrains[treestrain]))][0]
            sIDs = [id[0] for id in curSim.execute("SELECT spacer_id FROM bspacers \
                WHERE run_id = {0}\
                AND bstrain_id = {1}\
                ORDER BY spacer_id".format(run_id,treebstrains[treestrain]))]
            axesB[1].text(time, treestrain, 'bID:{0}::[{1}], vID:{2}'\
            .format(treebstrains[treestrain],', '.join(map(str,sIDs)),vID),\
            fontsize=3)#, fontweight='bold')
            axesS[1].text(time, treestrain, 'bID:{0}::[{1}], vID:{2}'\
            .format(treebstrains[treestrain],', '.join(map(str,sIDs)),vID),\
            fontsize=3)#, fontweight='bold')
    else:
        for treestrain in treebstrains.keys():
            axesB[1].text(time, treestrain, '{0}'\
            .format(treebstrains[treestrain]), fontsize=3)#, fontweight='bold')
        for treestrain in treebstrains.keys():
            axesS[1].text(time, treestrain, '{0}'\
            .format(treebstrains[treestrain]), fontsize=3)#, fontweight='bold')

print('Compiling trees with time slices...')
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
        # html_str = mpld3.fig_to_html(figV)
        # Html_file= open(os.path.join(PLOT_PATH,'virus-clades-abundances-tree-time-slices_microbe-overlay.html'),"w")
        # Html_file.write(html_str)
        # Html_file.close()
    else:
        figB.tight_layout()
        figB.savefig(os.path.join(PLOT_PATH,'microbe-clades-abundances-tree-time-slices.{0}'\
            .format(imgType)),dpi=resolve)
        figV.tight_layout()
        figV.savefig(os.path.join(PLOT_PATH,'virus-clades-abundances-tree-time-slices.{0}'\
            .format(imgType)),dpi=resolve)
    figS.tight_layout()
    figS.savefig(os.path.join(PLOT_PATH,'susceptible-clades-abundances-tree-time-slices.{0}'\
        .format(imgType)),dpi=resolve)


print('Complete!')
