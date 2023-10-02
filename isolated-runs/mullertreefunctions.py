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
import seaborn as sns

def speciesMullerPlot(run_id,species,DBSIM_PATH,DBTREE_PATH,treepalette,maxticksize,figxy,hratio,abundthreshold,logScale):
    if species == 'virus':
        strains = 'vstrains'
        strain = 'vstrain'
        strain_id = 'vstrain_id'
        abundanceTitle = 'Viral Strain Abundances'
        strain_abundance = "viral_abundance"
        straintotal = "vtotal"
        abundance = 'vabundance'
        tree_strain_id = "tree_vstrain_id"
        tree_parent_strain_id = "tree_parent_vstrain_id"
        parent_strain_id = "parent_vstrain_id"
    else:
        strains = 'bstrains'
        strain = 'bstrain'
        strain_id = 'bstrain_id'
        abundanceTitle = 'Microbial Strain Abundances'
        strain_abundance = "microbial_abundance"
        straintotal = "btotal"
        abundance = 'babundance'
        tree_strain_id = "tree_bstrain_id"
        tree_parent_strain_id = "tree_parent_bstrain_id"
        parent_strain_id = "parent_bstrain_id"
    conSim = sqlite3.connect(DBSIM_PATH)
    curSim = conSim.cursor()
    ID = curSim.execute('SELECT combo_id,replicate FROM runs WHERE run_id = {}'.format(run_id)).fetchall()
    combo_id = ID[0][0]
    replicate = ID[0][1]
    RUN_DIR = os.path.join('runID{0}-c{1}-r{2}'.format(run_id,combo_id,replicate))
    species_stacked = pd.read_sql_query("SELECT t,{0},abundance FROM {1} WHERE run_id = {2}".format(strain_id,abundance,run_id), conSim)
    t = max(pd.read_sql_query("SELECT t FROM vabundance WHERE run_id = {}".format(run_id), conSim).t.values)
    species_stacked = species_stacked[species_stacked.t <= t]
    species_total = pd.read_sql_query("SELECT t,{0} FROM summary WHERE run_id = {1}".format(strain_abundance,run_id), conSim)\
    .rename(columns={strain_abundance: straintotal})
    species_total = species_total[species_total.t <= max(species_stacked.t)]
    maxIDs = species_stacked.set_index('t').groupby([strain_id]).agg(t = ('abundance','idxmax'),\
                                maxAbund = ('abundance','max')).reset_index()
    maxIDs = species_total.merge(maxIDs,on=['t'])
    maxIDs[straintotal] = abundthreshold*np.array(maxIDs[straintotal])
    keepStrains = list(maxIDs[maxIDs['maxAbund']>maxIDs[straintotal]][strain_id].values)
    conTree = sqlite3.connect(DBTREE_PATH)
    curTree = conTree.cursor()
    print('SQLite Query: tree data')
    parentTrack = keepStrains.copy()
    keepStrainsNew = keepStrains.copy()
    while len(parentTrack) != 0:
        parentTrack = list(pd.read_sql_query(
        "SELECT parent_{1} \
        FROM {0} WHERE {1} in ({2})".format(strains,strain_id,', '.join(map(str,parentTrack))), conSim)\
        [parent_strain_id].values)
        # print(parentTrack)
        keepStrainsNew.extend(parentTrack)
        keepStrainsNew = list(np.unique(keepStrainsNew))
        # print(keepStrainsNew)
        parentTrack = list(np.array(parentTrack)[np.array(parentTrack) != 0])
    keepStrainsNew = list(*np.array(keepStrainsNew)[keepStrainsNew!=0])
    keepStrainsNew.sort()
    species_stacked = species_stacked[[(i in keepStrainsNew) for i in species_stacked[strain_id]]].reset_index(drop=True)
    if logScale:
        species_stacked['abundance'] = np.log(species_stacked['abundance'])
    keepTreeStrainsDF = pd.read_sql_query(
    "SELECT tree_{0}, {0} \
    FROM tree_{1}_order WHERE {0} in ({2})".format(strain_id,strain,', '.join(map(str,keepStrainsNew))), conTree)
    keepTreeStrains = list(keepTreeStrainsDF['tree_{}'.format(strain_id)].values)
    keepTreeStrains.sort()
    strainTimes = pd.read_sql_query(
    "SELECT tree_{0}, t_creation, t_extinction, tree_parent_{0} \
    FROM tree_{1}_creation_extinction WHERE tree_{0} in ({2})".format(
    strain_id,strain,', '.join(map(str,keepTreeStrains))), conTree)
    strainTimesNoExt = pd.read_sql_query(
    "SELECT {0}, t_creation, parent_{0} \
    FROM {1} WHERE {0} in ({2})".format(
    strain_id,strains,', '.join(map(str,keepTreeStrainsDF['bstrain_id']))), conSim)\
    .merge(keepTreeStrainsDF,on=strain_id).drop(columns=[strain_id])\
    .rename(columns={'parent_{0}'.format(strain_id):strain_id})\
    .merge(keepTreeStrainsDF.rename(columns={tree_strain_id: 'tree_parent_{0}'.format(strain_id)}), on=strain_id)\
    .drop(columns=[strain_id])
    strainTimesNoExt = strainTimesNoExt[[(i not in list(strainTimes[tree_strain_id])) \
                                         for i in list(strainTimesNoExt[tree_strain_id])]]
    strainTimesNoExt['t_extinction'] = len(strainTimesNoExt)*[t]
    strainTimes = pd.concat([strainTimes, strainTimesNoExt]).reset_index(drop=True)
    treeAbundances = pd.read_sql_query(
    "SELECT t, tree_{0}, abundance \
    FROM tree_{1} WHERE tree_{0} in ({2})"\
    .format(strain_id,abundance,', '.join(map(str,keepTreeStrains))), conTree)
    newTreeOrder = {}
    newTreeID = 1
    for treeStrain in keepTreeStrains:
        newTreeOrder[treeStrain] = newTreeID
        newTreeID += 1
    treeAbundances.replace({tree_strain_id:newTreeOrder},inplace=True)
    strainTimes.replace({tree_strain_id:newTreeOrder,tree_parent_strain_id:newTreeOrder},inplace=True)
    numStrains = len(keepTreeStrains)
    Vcmap = cm.get_cmap(treepalette)
    # Vcmap = sns.color_palette("icefire",as_cmap=True)
    norm = Normalize(vmin=float(1), vmax=float(max(newTreeOrder.values())))
    speciesColorDict = {}
    if logScale:
        maxAbundanceSpecies = np.log(np.max(treeAbundances.abundance.values))
    else:
        maxAbundanceSpecies = np.max(treeAbundances.abundance.values)
    markerIncSpecies = maxticksize/maxAbundanceSpecies
    markerColorsSpecies = []
    hlinecSpecies = []
    vlinecSpecies = []
    hcolorsSpecies = []
    vcolorsSpecies = []
    print('Compiling stacked species abundances and tree plots')
    fig, ax = plt.subplots(2,sharex=True, figsize=figxy, gridspec_kw={'height_ratios': hratio})
    # fig.suptitle('Strain Abundances (run{0}-c{1}-r{2})'.format(run_id,combo_id,replicate))
    axes = [ax[0], ax[1]]
    for strainID in sorted(treeAbundances[tree_strain_id].unique()):
        # print(strainID)
        tCreate = strainTimes[strainTimes[tree_strain_id] == strainID].t_creation.values[0]
        tExtinct = strainTimes[strainTimes[tree_strain_id] == strainID].t_extinction.values[0]
        parent = strainTimes[strainTimes[tree_strain_id] == strainID][tree_parent_strain_id].values[0]
        hlinecSpecies.append([[tCreate, strainID],[tExtinct, strainID]])
        vlinecSpecies.append([[tCreate, parent],[tCreate, strainID]])
        speciesColorDict[strainID] = Vcmap(norm(strainID))
        hcolorsSpecies.append(speciesColorDict[strainID])
        vcolorsSpecies.append(speciesColorDict[strainID])
        markerColorsSpecies.append(speciesColorDict[strainID])
        axes[1].scatter(treeAbundances[treeAbundances[tree_strain_id]==strainID]['t'],
        treeAbundances[treeAbundances[tree_strain_id]==strainID][tree_strain_id],
        lw=1,
        s=(treeAbundances[treeAbundances[tree_strain_id]==strainID]['abundance'].values)*markerIncSpecies,
        color = speciesColorDict[strainID],marker='|')
    strainLineages = LineCollection(hlinecSpecies, linestyles='solid', colors=hcolorsSpecies,linewidths=(1))
    creationLines = LineCollection(vlinecSpecies, linestyles='solid', colors=vcolorsSpecies, linewidths=(1))
    axes[1].add_collection(strainLineages)
    axes[1].add_collection(creationLines)
    axes[1].set_xlabel(xlabel = 'Time t',fontsize=7)
    axes[1].set_ylim(0,np.max(strainTimes[tree_strain_id])+1)
    axes[1].set_xlim(0,t)
    axes[1].xaxis.set_minor_locator(ticker.MultipleLocator(25))
    axes[1].set_yticks([])

    species_stacked = species_stacked.merge(keepTreeStrainsDF,on=[strain_id])\
    .drop(columns=[strain_id]).replace({tree_strain_id:newTreeOrder})
    species_stacked.sort_values(by=['t',tree_strain_id])
    speciesDF = species_stacked.copy()
    species_total = species_stacked.groupby(['t']).agg(total=('abundance','sum')).reset_index()
    species_stacked = species_stacked.merge(species_total,on=['t'])
    species_stacked['freq'] = species_stacked['abundance']/species_stacked['total']
    species_stacked = species_stacked.drop(columns=['abundance','total'])
    species_total['total'] = np.log10(species_total['total'])
    species_stacked = species_stacked.merge(species_total,on=['t'])
    species_stacked['freq'] = species_stacked['freq']*species_stacked['total']
    # species_total.plot(ax=axes[0], legend=False, color='white',sort_columns=True,linewidth=.1)
    species_stacked = species_stacked.pivot(index='t',columns=tree_strain_id,values='freq')
    species_stacked.plot.area(ax = axes[0],stacked=True,legend=False, linewidth=0,color=speciesColorDict,sort_columns=True)
    species_stacked.plot(stacked=True, ax=axes[0], legend=False, color='white',sort_columns=True,linewidth=.1)
 
    axes[0].set_ylabel(ylabel =abundanceTitle,labelpad=15,fontsize=15)
    axes[0].set_xlabel(xlabel = 'Time t',fontsize=15)
    # axes[0].ticklabel_format(axis = 'y',style='sci',scilimits=(0,0))
    axes[0].xaxis.set_minor_locator(ticker.MultipleLocator(25))
    lim = axes[0].get_ylim()
    axes[0].set_ylim(0,lim[1])
    axes[0].set_xlim(0,t)
    keepTreeStrainsDF['new_tree_{}_id'.format(strain)] = keepTreeStrainsDF['tree_{}_id'.format(strain)].copy()
    keepTreeStrainsDF = keepTreeStrainsDF.replace({'new_tree_{}_id'.format(strain): newTreeOrder})
    return speciesDF, keepTreeStrainsDF,speciesColorDict, hlinecSpecies, vlinecSpecies, fig, axes






bAbunds,bKeepTreeStrainsDF, bSpeciesColorDict, _, _, figB, axesB = \
speciesMullerPlot(run_id,'microbe',DBSIM_PATH,DBTREE_PATH,\
treepalette,maxticksize,figxy,hratio,babundthreshold, False)