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
import seaborn as sns


def truncateTree(left,right,top,bottom,vlinecVirus,hlinecVirus):
    vline3d2 = np.array(vlinecVirus)
    first = True
    for i in range(len(np.array(vlinecVirus))):
        line = np.array(vlinecVirus)[i]
        if (line[0][1] < bottom) & (line[1][1] < bottom):
            continue
        elif (line[0][1] < bottom) & (line[1][1] > top):
            line[0][1] = bottom
            line[1][1] = top
            if first:
                vline3d2 = np.array([line.copy()])
                first = False
            else:
                vline3d2 = np.append(vline3d2, np.array([line]), axis=0)
        elif (line[0][1] < bottom) & (line[1][1] >= bottom) & (line[1][1] <= top):
            line[0][1] = bottom
            if first:
                vline3d2 = np.array([line.copy()])
                first = False
            else:
                vline3d2 = np.append(vline3d2, np.array([line]), axis=0)
        elif (line[0][1] > top) & (line[1][1] > top):
            continue
        elif (line[0][1] >= bottom) & (line[0][1] < top) & (line[1][1] > top):
            line[1][1] = top
            if first:
                vline3d2 = np.array([line.copy()])
                first = False
            else:
                vline3d2 = np.append(vline3d2, np.array([line]), axis=0)
        else:
            if first:
                vline3d2 = np.array([line.copy()])
                first = False
            else:
                vline3d2 = np.append(vline3d2, np.array([line]), axis=0)
    first = True
    vline3d = np.array(vline3d2)
    for i in range(len(np.array(vline3d2))):
        line = np.array(vline3d2)[i]
        # print(i)
        # print(first)
        if (line[0][0] < left) & (line[1][0] < left):
            continue
        elif (line[0][0] < left) & (line[1][0] > right):
            line[0][0] = left
            line[1][0] = right
            if first:
                vline3d = np.array([line.copy()])
                first = False
            else:
                vline3d = np.append(vline3d, np.array([line]), axis=0)
        elif (line[0][0] < left) & (line[1][0] >= left) & (line[1][0] <= right):
            line[0][0] = left
            if first:
                vline3d = np.array([line.copy()])
                first = False
            else:
                vline3d = np.append(vline3d, np.array([line]), axis=0)
        elif (line[0][0] > right) & (line[1][0] > right):
            continue
        elif (line[0][0] >= left) & (line[0][0] <= right) & (line[1][0] > right):
            line[1][0] = right
            if first:
                vline3d = np.array([line.copy()])
                first = False
            else:
                vline3d = np.append(vline3d, np.array([line]), axis=0)
        else:
            if first:
                vline3d = np.array([line.copy()])
                first = False
            else:
                vline3d = np.append(vline3d, np.array([line]), axis=0)
    hline3d2 = np.array(hlinecVirus)
    first = True
    for i in range(len(np.array(hlinecVirus))):
        line = np.array(hlinecVirus)[i]
        if (line[0][1] < bottom) & (line[1][1] < bottom):
            continue
        elif (line[0][1] < bottom) & (line[1][1] > top):
            line[0][1] = bottom
            line[1][1] = top
            if first:
                hline3d2 = np.array([line.copy()])
                first = False
            else:
                hline3d2 = np.append(hline3d2, np.array([line]), axis=0)
        elif (line[0][1] < bottom) & (line[1][1] >= bottom) & (line[1][1] <= top):
            line[0][1] = bottom
            if first:
                hline3d2 = np.array([line.copy()])
                first = False
            else:
                hline3d2 = np.append(hline3d2, np.array([line]), axis=0)
        elif (line[0][1] > top) & (line[1][1] > top):
            continue
        elif (line[0][1] >= bottom) & (line[0][1] < top) & (line[1][1] > top):
            line[1][1] = top
            if first:
                hline3d2 = np.array([line.copy()])
                first = False
            else:
                hline3d2 = np.append(hline3d2, np.array([line]), axis=0)
        else:
            if first:
                hline3d2 = np.array([line.copy()])
                first = False
            else:
                hline3d2 = np.append(hline3d2, np.array([line]), axis=0)
    first = True
    hline3d = np.array(hline3d2)
    for i in range(len(np.array(hline3d2))):
        line = np.array(hline3d2)[i]
        # print(i)
        # print(first)
        if (line[0][0] < left) & (line[1][0] < left):
            continue
        elif (line[0][0] < left) & (line[1][0] > right):
            line[0][0] = left
            line[1][0] = right
            if first:
                hline3d = np.array([line.copy()])
                first = False
            else:
                hline3d = np.append(hline3d, np.array([line]), axis=0)
        elif (line[0][0] < left) & (line[1][0] >= left) & (line[1][0] <= right):
            line[0][0] = left
            if first:
                hline3d = np.array([line.copy()])
                first = False
            else:
                hline3d = np.append(hline3d, np.array([line]), axis=0)
        elif (line[0][0] > right) & (line[1][0] > right):
            continue
        elif (line[0][0] >= left) & (line[0][0] <= right) & (line[1][0] > right):
            line[1][0] = right
            if first:
                hline3d = np.array([line.copy()])
                first = False
            else:
                hline3d = np.append(hline3d, np.array([line]), axis=0)
        else:
            if first:
                hline3d = np.array([line.copy()])
                first = False
            else:
                hline3d = np.append(hline3d, np.array([line]), axis=0)
    return vline3d, hline3d


# this function has diiversity/abundance and tree. Also inlet as separte image to be joined manually.
def speciesTreeDiv2Plot(run_id, species, DBSIM_PATH, DBTREE_PATH, treepalette, maxticksize, figxy, hratio, abundthreshold,t0,tf):
    if species == 'virus':
        strains = 'vstrains'
        strain = 'vstrain'
        strain_id = 'vstrain_id'
        abundanceTitle = 'Viral\nAbundance'
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
        abundanceTitle = 'Host\nAbundance'
        strain_abundance = "microbial_abundance"
        straintotal = "btotal"
        abundance = 'babundance'
        tree_strain_id = "tree_bstrain_id"
        tree_parent_strain_id = "tree_parent_bstrain_id"
        parent_strain_id = "parent_bstrain_id"
    conSim = sqlite3.connect(DBSIM_PATH)
    curSim = conSim.cursor()
    ID = curSim.execute(
        'SELECT combo_id,replicate FROM runs WHERE run_id = {}'.format(run_id)).fetchall()
    combo_id = ID[0][0]
    replicate = ID[0][1]
    RUN_DIR = os.path.join(
        'runID{0}-c{1}-r{2}'.format(run_id, combo_id, replicate))
    species_stacked = pd.read_sql_query(
        "SELECT t,{0},abundance FROM {1} WHERE run_id = {2}".format(strain_id, abundance, run_id), conSim)
    t = max(pd.read_sql_query(
        "SELECT t FROM vabundance WHERE run_id = {}".format(run_id), conSim).t.values)
    species_stacked = species_stacked[species_stacked.t <= t]
    species_total = pd.read_sql_query("SELECT t,{0} FROM summary WHERE run_id = {1}".format(strain_abundance, run_id), conSim)\
        .rename(columns={strain_abundance: straintotal})
    species_total = species_total[species_total.t <= max(species_stacked.t)]
    # speciesRichness = species_stacked.groupby(['t']).agg(
    #     richness=(strain_id, 'size')).reset_index()
    speciesShannon = species_stacked.merge(species_total, on=['t'])
    speciesShannon['shannon'] = speciesShannon['abundance'] / \
        speciesShannon[straintotal]
    speciesShannon = speciesShannon[speciesShannon.shannon > 0]
    speciesShannon['shannon'] = -1 * \
        np.array(speciesShannon['shannon']) * \
        np.log(np.array(speciesShannon['shannon']))
    speciesShannon = speciesShannon.groupby(['t']).agg(
        shannon=('shannon', 'sum')).reset_index()
    speciesShannon['shannon'] = np.exp(np.array(speciesShannon['shannon']))
    maxIDs = species_stacked.set_index('t').groupby([strain_id]).agg(t=('abundance', 'idxmax'),
                                                                     maxAbund=('abundance', 'max')).reset_index()
    maxIDs = species_total.merge(maxIDs, on=['t'])
    maxIDs[straintotal] = abundthreshold*np.array(maxIDs[straintotal])
    keepStrains = list(
        maxIDs[maxIDs['maxAbund'] > maxIDs[straintotal]][strain_id].values)
    conTree = sqlite3.connect(DBTREE_PATH)
    curTree = conTree.cursor()
    print('SQLite Query: tree data')
    parentTrack = keepStrains.copy()
    keepStrainsNew = keepStrains.copy()
    while len(parentTrack) != 0:
        parentTrack = list(pd.read_sql_query(
            "SELECT parent_{1} \
        FROM {0} WHERE {1} in ({2})".format(strains, strain_id, ', '.join(map(str, parentTrack))), conSim)
            [parent_strain_id].values)
        # print(parentTrack)
        keepStrainsNew.extend(parentTrack)
        keepStrainsNew = list(np.unique(keepStrainsNew))
        # print(keepStrainsNew)
        parentTrack = list(np.array(parentTrack)[np.array(parentTrack) != 0])
    keepStrainsNew = list(*np.array(keepStrainsNew)[keepStrainsNew != 0])
    keepStrainsNew.sort()
    species_stacked = species_stacked[[
        (i in keepStrainsNew) for i in species_stacked[strain_id]]].reset_index(drop=True)
    keepTreeStrainsDF = pd.read_sql_query(
        "SELECT tree_{0}, {0} \
    FROM tree_{1}_order WHERE {0} in ({2})".format(strain_id, strain, ', '.join(map(str, keepStrainsNew))), conTree)
    keepTreeStrains = list(
        keepTreeStrainsDF['tree_{}'.format(strain_id)].values)
    keepTreeStrains.sort()
    strainTimes = pd.read_sql_query(
        "SELECT tree_{0}, t_creation, t_extinction, tree_parent_{0} \
    FROM tree_{1}_creation_extinction WHERE tree_{0} in ({2})".format(
            strain_id, strain, ', '.join(map(str, keepTreeStrains))), conTree)
    treeAbundances = pd.read_sql_query(
        "SELECT t, tree_{0}, abundance \
    FROM tree_{1} WHERE tree_{0} in ({2})"
        .format(strain_id, abundance, ', '.join(map(str, keepTreeStrains))), conTree)
    newTreeOrder = {}
    newTreeID = 0
    for treeStrain in keepTreeStrains:
        newTreeOrder[treeStrain] = newTreeID
        newTreeID += 1
    treeAbundances.replace({tree_strain_id: newTreeOrder}, inplace=True)
    strainTimes.replace({tree_strain_id: newTreeOrder,
                        tree_parent_strain_id: newTreeOrder}, inplace=True)
    numStrains = len(keepTreeStrains)
    Vcmap = cm.get_cmap(treepalette)
    # Vcmap = sns.color_palette("icefire",as_cmap=True)
    norm = Normalize(vmin=float(1), vmax=float(max(newTreeOrder.values())))
    speciesColorDict = {}
    maxAbundanceSpecies = np.max(treeAbundances.abundance.values)
    markerIncSpecies = maxticksize/maxAbundanceSpecies
    markerColorsSpecies = []
    hlinecSpecies = []
    vlinecSpecies = []
    hcolorsSpecies = []
    vcolorsSpecies = []
    print('Compiling stacked species abundances and tree plots')
    fig, ax = plt.subplots(1, sharex=True, figsize=[figxy[0],figxy[0]/4]) # collapse plot
    fig2, ax2 = plt.subplots(2, sharex=True, figsize=figxy,
                           gridspec_kw={'height_ratios': hratio,'hspace':0}) # tree diversity plot
    axes = [ax, ax.twinx()] # collapse plot
    axes2 = [ax2[0], ax2[1]] # tree diversity plot
    for strainID in sorted(treeAbundances[tree_strain_id].unique()):
        tCreate = strainTimes[strainTimes[tree_strain_id]
                              == strainID].t_creation.values[0]
        tExtinct = strainTimes[strainTimes[tree_strain_id]
                               == strainID].t_extinction.values[0]
        parent = strainTimes[strainTimes[tree_strain_id]
                             == strainID][tree_parent_strain_id].values[0]
        hlinecSpecies.append([[tCreate, strainID], [tExtinct, strainID]])
        if parent != 0:
            vlinecSpecies.append([[tCreate, parent], [tCreate, strainID]])
        speciesColorDict[strainID] = Vcmap(norm(strainID))
        hcolorsSpecies.append(speciesColorDict[strainID])
        vcolorsSpecies.append(speciesColorDict[strainID])
        markerColorsSpecies.append(speciesColorDict[strainID])
        axes2[1].scatter(treeAbundances[treeAbundances[tree_strain_id] == strainID]['t'],
                        treeAbundances[treeAbundances[tree_strain_id]
                                       == strainID][tree_strain_id],
                        lw=1,
                        s=treeAbundances[treeAbundances[tree_strain_id] ==
                                         strainID]['abundance'].values*markerIncSpecies,
                        color=speciesColorDict[strainID], marker='|')
    strainLineages = LineCollection(
        hlinecSpecies, linestyles='solid', colors=hcolorsSpecies, linewidths=(.7))
    creationLines = LineCollection(
        vlinecSpecies, linestyles='solid', colors=vcolorsSpecies, linewidths=(.7))
    axes2[1].add_collection(strainLineages)
    axes2[1].add_collection(creationLines)
    axes2[1].set_ylim(0, np.max(strainTimes[tree_strain_id])+1)
    axes2[1].set_yticks([])
    species_stacked = species_stacked.merge(keepTreeStrainsDF, on=[strain_id])\
        .drop(columns=[strain_id]).replace({tree_strain_id: newTreeOrder})
    species_stacked.sort_values(by=['t', tree_strain_id])
    speciesDF = species_stacked.copy()
    species_stacked = species_stacked.pivot(
        index='t', columns=tree_strain_id, values='abundance')
    if species == 'virus':
        axes2[1].set_ylabel(ylabel='Viral Phylogeny', labelpad=15, fontsize=10)
        ## collapse plot
        microbe_total = pd.read_sql_query("SELECT t, microbial_abundance FROM summary WHERE run_id = {}".format(run_id), conSim)\
            .rename(columns={'microbial_abundance': 'btotal'})
        microbe_total = microbe_total[microbe_total.t <= max(
            species_stacked.index)]
        axes[0].plot(microbe_total['t'], microbe_total['btotal'],
                        linewidth=0, color='grey')
        axes[0].fill_between(
            microbe_total['t'], microbe_total['btotal'], color='grey', alpha=0.75)
        species_stacked.plot(ax=axes[1], stacked=False, legend=False,
                                linewidth=1.5, color=speciesColorDict, sort_columns=True)
        axes[1].set_ylabel(ylabel='Viral Strain\nAbundances',
                            labelpad=30, fontsize=10, rotation=270)
        axes[0].set_ylabel(ylabel='Host\nAbundance', labelpad=15, fontsize=10)
        ## tree diversity plot
        axes2[0].plot(microbe_total['t'], microbe_total['btotal'],
                     linewidth=0, color='grey')
        axes2[0].fill_between(
            microbe_total['t'], microbe_total['btotal'], color='grey', alpha=0.5)
        axes2[0].set_yticklabels([])
        axes2[0].set_yticks([])
        axes2.append(axes2[0].twinx())
        species_stacked.plot.area(ax=axes2[2], stacked=True, legend=False,
                                    linewidth=0, color=speciesColorDict, sort_columns=True, alpha=0.25)
        species_stacked.plot(stacked=True, ax=axes2[2], legend=False,
                                color='white', sort_columns=True, linewidth=.1)
        axes2.append(axes2[0].twinx())
        # axes2[3].plot(speciesRichness['t'], speciesRichness['richness'],
        #                 color='darkred', linewidth=1.5)
        # axes2[3].yaxis.tick_right()
        # axes2[3].yaxis.set_label_position("right")
        # axes2[3].set_ylabel(ylabel='Viral Richness',
        #                     labelpad=30, fontsize=10, rotation=270)
        axes2[3].plot(speciesShannon['t'], speciesShannon['shannon'],
                        color='darkblue', linewidth=1.5)
        axes2[3].yaxis.tick_left()
        axes2[3].yaxis.set_label_position("left")
        axes2[3].set_ylabel(
            ylabel='Viral\n Shannon Diversity', labelpad=15, fontsize=10)
    if species == 'microbe':
        axes2[1].set_ylabel(ylabel='Host Phylogeny',
                            labelpad=15, fontsize=10)
        ## collapse plot
        virus_total = pd.read_sql_query("SELECT t, viral_abundance FROM summary WHERE run_id = {}".format(run_id), conSim)\
            .rename(columns={'viral_abundance': 'vtotal'})
        virus_total = virus_total[virus_total.t <= max(species_stacked.index)]
        axes[1].plot(virus_total[(virus_total.t >= t0) & (virus_total.t <= tf)]['t'],
                     virus_total[(virus_total.t>=t0)&(virus_total.t<=tf)]['vtotal'],
                        linewidth=0, color='grey')
        axes[1].fill_between(
            virus_total[(virus_total.t >= t0) & (virus_total.t <= tf)]['t'], \
            virus_total[(virus_total.t >= t0) & (virus_total.t <= tf)]['vtotal'], color='grey', alpha=0.75)
        speciesDF[(speciesDF.t>=t0)&(speciesDF.t<=tf)]\
            .pivot(index='t', columns=tree_strain_id, values='abundance').plot(ax=axes[0], stacked=False, legend=False,
                             linewidth=1.5, color=speciesColorDict, sort_columns=True)
        axes[0].set_ylabel(
            ylabel='Host Strain\nAbundances', labelpad=15, fontsize=10)
        axes[1].set_ylabel(ylabel='Viral\nAbundance',
                        labelpad=30, fontsize=10, rotation=270)
        ## tree diversity plot
        species_stacked.plot.area(ax=axes2[0], stacked=True, legend=False,
                                  linewidth=0, color=speciesColorDict, sort_columns=True, alpha=0.25)
        species_stacked.plot(stacked=True, ax=axes2[0], legend=False,
                             color='white', sort_columns=True, linewidth=.1)
        axes2[0].set_yticklabels([])
        axes2[0].set_yticks([])
        axes2.append(axes2[0].twinx())
        axes2[2].plot(virus_total['t'], virus_total['vtotal'],
                      linewidth=0, color='grey')
        axes2[2].fill_between(
            virus_total['t'], virus_total['vtotal'], color='grey', alpha=0.5)
        axes2[2].set_yticklabels([])
        axes2[2].set_yticks([])
        axes2.append(axes2[0].twinx())
        axes2[3].plot(speciesShannon['t'], speciesShannon['shannon'],
                        color='darkblue', linewidth=1.5)
        axes2[3].yaxis.tick_left()
        axes2[3].yaxis.set_label_position("left")
        axes2[3].set_ylabel(
            ylabel='Host\n Shannon Diversity', labelpad=15, fontsize=10)
    ### settings for all plots
    # axes[0].set_xlim(0, t)
    axes[0].margins(x=0)
    axes[0].xaxis.set_minor_locator(ticker.MultipleLocator(25))
    lim = axes[0].get_ylim()
    axes[0].set_ylim(0, lim[1])
    axes[0].ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
    lim = axes[1].get_ylim()
    axes[1].set_ylim(0, lim[1])
    axes[1].margins(x=0)
    axes[1].xaxis.set_minor_locator(ticker.MultipleLocator(25))
    axes[1].ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
    axes[0].set_xlabel(xlabel='Time t', fontsize=10)
    lim = axes2[0].get_ylim()
    axes2[0].set_ylim(0, lim[1])
    axes2[0].set_xlim(0, t)
    axes2[0].margins(x=0)
    axes2[0].xaxis.set_minor_locator(ticker.MultipleLocator(25))
    axes2[1].set_xlabel(xlabel='Time t', fontsize=10)
    axes2[1].set_xlim(0, t)
    axes2[1].margins(x=0)
    axes2[1].xaxis.set_minor_locator(ticker.MultipleLocator(25))
    lim = axes2[2].get_ylim()
    axes2[2].set_ylim(0, lim[1])
    axes2[2].set_xlim(0, t)
    axes2[2].margins(x=0)
    axes2[2].xaxis.set_minor_locator(ticker.MultipleLocator(25))
    lim = axes2[3].get_ylim()
    axes2[3].set_ylim(0, lim[1])
    axes2[3].set_xlim(0, t)
    axes2[3].margins(x=0)
    axes2[3].xaxis.set_minor_locator(ticker.MultipleLocator(25))
    ###  
    keepTreeStrainsDF['new_tree_{}_id'.format(
        strain)] = keepTreeStrainsDF['tree_{}_id'.format(strain)].copy()
    keepTreeStrainsDF = keepTreeStrainsDF.replace(
        {'new_tree_{}_id'.format(strain): newTreeOrder})
    return speciesDF, keepTreeStrainsDF, speciesColorDict, hlinecSpecies, vlinecSpecies, fig, axes, fig2, axes2


# this function has diiversity/abundance, tree, and inlet but inlet is optimal
def speciesTreeDiv3Plot(run_id, species, DBSIM_PATH, DBTREE_PATH, treepalette, maxticksize, figxy, hratio, abundthreshold, t0, tf):
    if species == 'virus':
        strains = 'vstrains'
        strain = 'vstrain'
        strain_id = 'vstrain_id'
        abundanceTitle = 'Viral\nAbundance'
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
        abundanceTitle = 'Host\nAbundance'
        strain_abundance = "microbial_abundance"
        straintotal = "btotal"
        abundance = 'babundance'
        tree_strain_id = "tree_bstrain_id"
        tree_parent_strain_id = "tree_parent_bstrain_id"
        parent_strain_id = "parent_bstrain_id"
    conSim = sqlite3.connect(DBSIM_PATH)
    curSim = conSim.cursor()
    ID = curSim.execute(
        'SELECT combo_id,replicate FROM runs WHERE run_id = {}'.format(run_id)).fetchall()
    combo_id = ID[0][0]
    replicate = ID[0][1]
    RUN_DIR = os.path.join(
        'runID{0}-c{1}-r{2}'.format(run_id, combo_id, replicate))
    species_stacked = pd.read_sql_query(
        "SELECT t,{0},abundance FROM {1} WHERE run_id = {2}".format(strain_id, abundance, run_id), conSim)
    t = max(pd.read_sql_query(
        "SELECT t FROM vabundance WHERE run_id = {}".format(run_id), conSim).t.values)
    species_stacked = species_stacked[species_stacked.t <= t]
    species_total = pd.read_sql_query("SELECT t,{0} FROM summary WHERE run_id = {1}".format(strain_abundance, run_id), conSim)\
        .rename(columns={strain_abundance: straintotal})
    species_total = species_total[species_total.t <= max(species_stacked.t)]
    # speciesRichness = species_stacked.groupby(['t']).agg(
    #     richness=(strain_id, 'size')).reset_index()
    speciesShannon = species_stacked.merge(species_total, on=['t'])
    speciesShannon['shannon'] = speciesShannon['abundance'] / \
        speciesShannon[straintotal]
    speciesShannon = speciesShannon[speciesShannon.shannon > 0]
    speciesShannon['shannon'] = -1 * \
        np.array(speciesShannon['shannon']) * \
        np.log(np.array(speciesShannon['shannon']))
    speciesShannon = speciesShannon.groupby(['t']).agg(
        shannon=('shannon', 'sum')).reset_index()
    speciesShannon['shannon'] = np.exp(np.array(speciesShannon['shannon']))
    maxIDs = species_stacked.set_index('t').groupby([strain_id]).agg(t=('abundance', 'idxmax'),
                                                                     maxAbund=('abundance', 'max')).reset_index()
    maxIDs = species_total.merge(maxIDs, on=['t'])
    maxIDs[straintotal] = abundthreshold*np.array(maxIDs[straintotal])
    keepStrains = list(
        maxIDs[maxIDs['maxAbund'] > maxIDs[straintotal]][strain_id].values)
    conTree = sqlite3.connect(DBTREE_PATH)
    curTree = conTree.cursor()
    print('SQLite Query: tree data')
    parentTrack = keepStrains.copy()
    keepStrainsNew = keepStrains.copy()
    while len(parentTrack) != 0:
        parentTrack = list(pd.read_sql_query(
            "SELECT parent_{1} \
        FROM {0} WHERE {1} in ({2})".format(strains, strain_id, ', '.join(map(str, parentTrack))), conSim)
            [parent_strain_id].values)
        # print(parentTrack)
        keepStrainsNew.extend(parentTrack)
        keepStrainsNew = list(np.unique(keepStrainsNew))
        # print(keepStrainsNew)
        parentTrack = list(np.array(parentTrack)[np.array(parentTrack) != 0])
    keepStrainsNew = list(*np.array(keepStrainsNew)[keepStrainsNew != 0])
    keepStrainsNew.sort()
    species_stacked = species_stacked[[
        (i in keepStrainsNew) for i in species_stacked[strain_id]]].reset_index(drop=True)
    keepTreeStrainsDF = pd.read_sql_query(
        "SELECT tree_{0}, {0} \
    FROM tree_{1}_order WHERE {0} in ({2})".format(strain_id, strain, ', '.join(map(str, keepStrainsNew))), conTree)
    keepTreeStrains = list(
        keepTreeStrainsDF['tree_{}'.format(strain_id)].values)
    keepTreeStrains.sort()
    strainTimes = pd.read_sql_query(
        "SELECT tree_{0}, t_creation, t_extinction, tree_parent_{0} \
    FROM tree_{1}_creation_extinction WHERE tree_{0} in ({2})".format(
            strain_id, strain, ', '.join(map(str, keepTreeStrains))), conTree)
    treeAbundances = pd.read_sql_query(
        "SELECT t, tree_{0}, abundance \
    FROM tree_{1} WHERE tree_{0} in ({2})"
        .format(strain_id, abundance, ', '.join(map(str, keepTreeStrains))), conTree)
    newTreeOrder = {}
    newTreeID = 0
    for treeStrain in keepTreeStrains:
        newTreeOrder[treeStrain] = newTreeID
        newTreeID += 1
    treeAbundances.replace({tree_strain_id: newTreeOrder}, inplace=True)
    strainTimes.replace({tree_strain_id: newTreeOrder,
                        tree_parent_strain_id: newTreeOrder}, inplace=True)
    numStrains = len(keepTreeStrains)
    Vcmap = cm.get_cmap(treepalette)
    # Vcmap = sns.color_palette("icefire",as_cmap=True)
    norm = Normalize(vmin=float(1), vmax=float(max(newTreeOrder.values())))
    speciesColorDict = {}
    maxAbundanceSpecies = np.max(treeAbundances.abundance.values)
    markerIncSpecies = maxticksize/maxAbundanceSpecies
    markerColorsSpecies = []
    hlinecSpecies = []
    vlinecSpecies = []
    hcolorsSpecies = []
    vcolorsSpecies = []
    print('Compiling stacked species abundances and tree plots')
    fig2, ax2 = plt.subplots(3, figsize=figxy,
                             gridspec_kw={'height_ratios': hratio})  # tree diversity plot
    ax2[0].get_shared_x_axes().join(ax2[0], ax2[1])
    ax2[0].set_xticklabels([])
    axes2 = [ax2[0], ax2[1], ax2[2], ax2[2].twinx()]  # tree diversity plot
    for strainID in sorted(treeAbundances[tree_strain_id].unique()):
        tCreate = strainTimes[strainTimes[tree_strain_id]
                              == strainID].t_creation.values[0]
        tExtinct = strainTimes[strainTimes[tree_strain_id]
                               == strainID].t_extinction.values[0]
        parent = strainTimes[strainTimes[tree_strain_id]
                             == strainID][tree_parent_strain_id].values[0]
        hlinecSpecies.append([[tCreate, strainID], [tExtinct, strainID]])
        if parent != 0:
            vlinecSpecies.append([[tCreate, parent], [tCreate, strainID]])
        speciesColorDict[strainID] = Vcmap(norm(strainID))
        hcolorsSpecies.append(speciesColorDict[strainID])
        vcolorsSpecies.append(speciesColorDict[strainID])
        markerColorsSpecies.append(speciesColorDict[strainID])
        axes2[1].scatter(treeAbundances[treeAbundances[tree_strain_id] == strainID]['t'],
                         treeAbundances[treeAbundances[tree_strain_id]
                                        == strainID][tree_strain_id],
                         lw=1,
                         s=treeAbundances[treeAbundances[tree_strain_id] ==
                                          strainID]['abundance'].values*markerIncSpecies,
                         color=speciesColorDict[strainID], marker='|')
    strainLineages = LineCollection(
        hlinecSpecies, linestyles='solid', colors=hcolorsSpecies, linewidths=(.7))
    creationLines = LineCollection(
        vlinecSpecies, linestyles='solid', colors=vcolorsSpecies, linewidths=(.7))
    axes2[1].add_collection(strainLineages)
    axes2[1].add_collection(creationLines)
    axes2[1].set_ylim(0, np.max(strainTimes[tree_strain_id])+1)
    axes2[1].set_yticks([])
    species_stacked = species_stacked.merge(keepTreeStrainsDF, on=[strain_id])\
        .drop(columns=[strain_id]).replace({tree_strain_id: newTreeOrder})
    species_stacked.sort_values(by=['t', tree_strain_id])
    speciesDF = species_stacked.copy()
    species_stacked = species_stacked.pivot(
        index='t', columns=tree_strain_id, values='abundance')
    if species == 'virus':
        axes2[1].set_ylabel(ylabel='Viral Strains', labelpad=15, fontsize=10)
        # collapse plot
        microbe_total = pd.read_sql_query("SELECT t, microbial_abundance FROM summary WHERE run_id = {}".format(run_id), conSim)\
            .rename(columns={'microbial_abundance': 'btotal'})
        microbe_total = microbe_total[microbe_total.t <= max(
            species_stacked.index)]
        axes2[2].plot(microbe_total['t'], microbe_total['btotal'],
                     linewidth=0, color='grey')
        axes2[2].fill_between(
            microbe_total['t'], microbe_total['btotal'], color='grey', alpha=0.75)
        species_stacked.plot(ax=axes2[3], stacked=False, legend=False,
                             linewidth=1.5, color=speciesColorDict, sort_columns=True)
        axes2[3].set_ylabel(ylabel='Viral Strain\nAbundances',
                           labelpad=30, fontsize=10, rotation=270)
        axes2[2].set_ylabel(ylabel='Host\nAbundance', labelpad=15, fontsize=10)
        # tree diversity plot
        axes2[0].plot(microbe_total['t'], microbe_total['btotal'],
                      linewidth=0, color='grey')
        axes2[0].fill_between(
            microbe_total['t'], microbe_total['btotal'], color='grey', alpha=0.5)
        axes2[0].set_yticklabels([])
        axes2[0].set_yticks([])
        axes2.append(axes2[0].twinx())
        species_stacked.plot.area(ax=axes2[2], stacked=True, legend=False,
                                  linewidth=0, color=speciesColorDict, sort_columns=True, alpha=0.25)
        species_stacked.plot(stacked=True, ax=axes2[2], legend=False,
                             color='white', sort_columns=True, linewidth=.1)
        axes2.append(axes2[0].twinx())
        # axes2[3].plot(speciesRichness['t'], speciesRichness['richness'],
        #                 color='darkred', linewidth=1.5)
        # axes2[3].yaxis.tick_right()
        # axes2[3].yaxis.set_label_position("right")
        # axes2[3].set_ylabel(ylabel='Viral Richness',
        #                     labelpad=30, fontsize=10, rotation=270)
        axes2[3].plot(speciesShannon['t'], speciesShannon['shannon'],
                      color='darkblue', linewidth=1.5)
        axes2[3].yaxis.tick_left()
        axes2[3].yaxis.set_label_position("left")
        axes2[3].set_ylabel(
            ylabel='Viral\n Shannon Diversity', labelpad=15, fontsize=10)
    if species == 'microbe':
        axes2[1].set_ylabel(ylabel='Host Immune Strains',
                            labelpad=15, fontsize=10)
        # collapse plot
        virus_total = pd.read_sql_query("SELECT t, viral_abundance FROM summary WHERE run_id = {}".format(run_id), conSim)\
            .rename(columns={'viral_abundance': 'vtotal'})
        virus_total = virus_total[(virus_total.t >= t0) & (virus_total.t <= tf)]
        axes2[3].plot(virus_total['t'], virus_total['vtotal'],linewidth=0, color='grey')
        axes2[3].fill_between(virus_total['t'], virus_total['vtotal'], color='grey', alpha=0.75)
        speciesDF[(speciesDF.t>=t0)&(speciesDF.t<=tf)]\
            .pivot(index='t', columns=tree_strain_id, values='abundance').plot(ax=axes2[2], stacked=False, legend=False,
                             linewidth=1.5, color=speciesColorDict, sort_columns=True)
        axes2[2].set_ylabel(
            ylabel='Host Strain\nAbundances', labelpad=15, fontsize=10)
        axes2[3].set_ylabel(ylabel='Viral\nAbundance',
                           labelpad=30, fontsize=10, rotation=270)
        # tree diversity plot
        species_stacked.plot.area(ax=axes2[0], stacked=True, legend=False,
                                  linewidth=0, color=speciesColorDict, sort_columns=True, alpha=0.25)
        species_stacked.plot(stacked=True, ax=axes2[0], legend=False,
                             color='white', sort_columns=True, linewidth=.1)
        virus_total = pd.read_sql_query("SELECT t, viral_abundance FROM summary WHERE run_id = {}".format(run_id), conSim)\
            .rename(columns={'viral_abundance': 'vtotal'})
        virus_total = virus_total[virus_total.t <= max(species_stacked.index)]
        axes2[0].set_yticklabels([])
        axes2[0].set_yticks([])
        axes2.append(axes2[0].twinx())
        axes2[4].plot(virus_total['t'], virus_total['vtotal'],
                      linewidth=0, color='grey')
        axes2[4].fill_between(
            virus_total['t'], virus_total['vtotal'], color='grey', alpha=0.5)
        axes2[4].set_yticklabels([])
        axes2[4].set_yticks([])
        axes2.append(axes2[0].twinx())
        axes2[5].plot(speciesShannon['t'], speciesShannon['shannon'],
                      color='darkblue', linewidth=1.5)
        axes2[5].yaxis.tick_left()
        axes2[5].yaxis.set_label_position("left")
        axes2[5].set_ylabel(
            ylabel='Host\n Shannon Diversity', labelpad=15, fontsize=10)
    # settings for all plots
    axes2[2].set_xlim(0, t)
    axes2[2].margins(x=0)
    axes2[2].xaxis.set_minor_locator(ticker.MultipleLocator(25))
    lim = axes2[2].get_ylim()
    axes2[2].set_ylim(0, lim[1])
    # axes2[2].ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
    lim = axes2[3].get_ylim()
    axes2[3].set_ylim(0, lim[1])
    axes2[3].margins(x=0)
    axes2[3].xaxis.set_minor_locator(ticker.MultipleLocator(25))
    axes2[3].ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
    axes2[2].set_xlabel(xlabel='Time t', fontsize=10)
    lim = axes2[0].get_ylim()
    axes2[0].set_ylim(0, lim[1])
    axes2[0].set_xlim(0, t)
    axes2[0].margins(x=0)
    axes2[0].xaxis.set_minor_locator(ticker.MultipleLocator(25))
    axes2[1].set_xlabel(xlabel='Time t', fontsize=10)
    axes2[1].set_xlim(0, t)
    axes2[1].margins(x=0)
    axes2[1].xaxis.set_minor_locator(ticker.MultipleLocator(25))
    lim = axes2[2].get_ylim()
    axes2[2].set_ylim(0, lim[1])
    axes2[2].set_xlim(0, t)
    axes2[2].margins(x=0)
    axes2[2].xaxis.set_minor_locator(ticker.MultipleLocator(25))
    lim = axes2[3].get_ylim()
    axes2[3].set_ylim(0, lim[1])
    axes2[3].set_xlim(0, t)
    axes2[3].margins(x=0)
    axes2[3].xaxis.set_minor_locator(ticker.MultipleLocator(25))
    ###
    keepTreeStrainsDF['new_tree_{}_id'.format(
        strain)] = keepTreeStrainsDF['tree_{}_id'.format(strain)].copy()
    keepTreeStrainsDF = keepTreeStrainsDF.replace(
        {'new_tree_{}_id'.format(strain): newTreeOrder})
    return speciesDF, keepTreeStrainsDF, speciesColorDict, hlinecSpecies, vlinecSpecies, fig2, axes2



# this function outputs line abundances, tree, then diversity/stacked baundances
def speciesTreeDivPlot(run_id, species, DBSIM_PATH, DBTREE_PATH, treepalette, maxticksize, figxy, hratio, abundthreshold, stacked, overlay):
    if species == 'virus':
        strains = 'vstrains'
        strain = 'vstrain'
        strain_id = 'vstrain_id'
        abundanceTitle = 'Viral\nAbundance'
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
        abundanceTitle = 'Host\nAbundance'
        strain_abundance = "microbial_abundance"
        straintotal = "btotal"
        abundance = 'babundance'
        tree_strain_id = "tree_bstrain_id"
        tree_parent_strain_id = "tree_parent_bstrain_id"
        parent_strain_id = "parent_bstrain_id"
    conSim = sqlite3.connect(DBSIM_PATH)
    curSim = conSim.cursor()
    ID = curSim.execute(
        'SELECT combo_id,replicate FROM runs WHERE run_id = {}'.format(run_id)).fetchall()
    combo_id = ID[0][0]
    replicate = ID[0][1]
    RUN_DIR = os.path.join(
        'runID{0}-c{1}-r{2}'.format(run_id, combo_id, replicate))
    species_stacked = pd.read_sql_query(
        "SELECT t,{0},abundance FROM {1} WHERE run_id = {2}".format(strain_id, abundance, run_id), conSim)
    t = max(pd.read_sql_query(
        "SELECT t FROM vabundance WHERE run_id = {}".format(run_id), conSim).t.values)
    species_stacked = species_stacked[species_stacked.t <= t]
    species_total = pd.read_sql_query("SELECT t,{0} FROM summary WHERE run_id = {1}".format(strain_abundance, run_id), conSim)\
        .rename(columns={strain_abundance: straintotal})
    species_total = species_total[species_total.t <= max(species_stacked.t)]
    speciesRichness = species_stacked.groupby(['t']).agg(richness=(strain_id,'size')).reset_index()
    speciesShannon = species_stacked.merge(species_total,on=['t'])
    speciesShannon['shannon'] = speciesShannon['abundance']/speciesShannon[straintotal]
    speciesShannon = speciesShannon[speciesShannon.shannon>0]
    speciesShannon['shannon'] = -1*np.array(speciesShannon['shannon'])*np.log(np.array(speciesShannon['shannon']))
    speciesShannon = speciesShannon.groupby(['t']).agg(shannon=('shannon', 'sum')).reset_index()
    speciesShannon['shannon'] = np.exp(np.array(speciesShannon['shannon']))
    maxIDs = species_stacked.set_index('t').groupby([strain_id]).agg(t=('abundance', 'idxmax'),
                                                                     maxAbund=('abundance', 'max')).reset_index()
    maxIDs = species_total.merge(maxIDs, on=['t'])
    maxIDs[straintotal] = abundthreshold*np.array(maxIDs[straintotal])
    keepStrains = list(
        maxIDs[maxIDs['maxAbund'] > maxIDs[straintotal]][strain_id].values)
    conTree = sqlite3.connect(DBTREE_PATH)
    curTree = conTree.cursor()
    print('SQLite Query: tree data')
    parentTrack = keepStrains.copy()
    keepStrainsNew = keepStrains.copy()
    while len(parentTrack) != 0:
        parentTrack = list(pd.read_sql_query(
            "SELECT parent_{1} \
        FROM {0} WHERE {1} in ({2})".format(strains, strain_id, ', '.join(map(str, parentTrack))), conSim)
            [parent_strain_id].values)
        # print(parentTrack)
        keepStrainsNew.extend(parentTrack)
        keepStrainsNew = list(np.unique(keepStrainsNew))
        # print(keepStrainsNew)
        parentTrack = list(np.array(parentTrack)[np.array(parentTrack) != 0])
    keepStrainsNew = list(*np.array(keepStrainsNew)[keepStrainsNew != 0])
    keepStrainsNew.sort()
    species_stacked = species_stacked[[
        (i in keepStrainsNew) for i in species_stacked[strain_id]]].reset_index(drop=True)
    keepTreeStrainsDF = pd.read_sql_query(
        "SELECT tree_{0}, {0} \
    FROM tree_{1}_order WHERE {0} in ({2})".format(strain_id, strain, ', '.join(map(str, keepStrainsNew))), conTree)
    keepTreeStrains = list(
        keepTreeStrainsDF['tree_{}'.format(strain_id)].values)
    keepTreeStrains.sort()
    strainTimes = pd.read_sql_query(
        "SELECT tree_{0}, t_creation, t_extinction, tree_parent_{0} \
    FROM tree_{1}_creation_extinction WHERE tree_{0} in ({2})".format(
            strain_id, strain, ', '.join(map(str, keepTreeStrains))), conTree)
    treeAbundances = pd.read_sql_query(
        "SELECT t, tree_{0}, abundance \
    FROM tree_{1} WHERE tree_{0} in ({2})"
        .format(strain_id, abundance, ', '.join(map(str, keepTreeStrains))), conTree)
    newTreeOrder = {}
    newTreeID = 0
    for treeStrain in keepTreeStrains:
        newTreeOrder[treeStrain] = newTreeID
        newTreeID += 1
    treeAbundances.replace({tree_strain_id: newTreeOrder}, inplace=True)
    strainTimes.replace({tree_strain_id: newTreeOrder,
                        tree_parent_strain_id: newTreeOrder}, inplace=True)
    numStrains = len(keepTreeStrains)
    Vcmap = cm.get_cmap(treepalette)
    # Vcmap = sns.color_palette("icefire",as_cmap=True)
    norm = Normalize(vmin=float(1), vmax=float(max(newTreeOrder.values())))
    speciesColorDict = {}
    maxAbundanceSpecies = np.max(treeAbundances.abundance.values)
    markerIncSpecies = maxticksize/maxAbundanceSpecies
    markerColorsSpecies = []
    hlinecSpecies = []
    vlinecSpecies = []
    hcolorsSpecies = []
    vcolorsSpecies = []
    print('Compiling stacked species abundances and tree plots')
    fig, ax = plt.subplots(3, sharex=True, figsize=figxy,
                           gridspec_kw={'height_ratios': hratio})
    # fig.suptitle('Strain Abundances (run{0}-c{1}-r{2})'.format(run_id,combo_id,replicate))
    axes = [ax[0], ax[1], ax[2]]
    for strainID in sorted(treeAbundances[tree_strain_id].unique()):
        tCreate = strainTimes[strainTimes[tree_strain_id]
                              == strainID].t_creation.values[0]
        tExtinct = strainTimes[strainTimes[tree_strain_id]
                               == strainID].t_extinction.values[0]
        parent = strainTimes[strainTimes[tree_strain_id]
                             == strainID][tree_parent_strain_id].values[0]
        hlinecSpecies.append([[tCreate, strainID], [tExtinct, strainID]])
        if parent != 0:
            vlinecSpecies.append([[tCreate, parent], [tCreate, strainID]])
        speciesColorDict[strainID] = Vcmap(norm(strainID))
        hcolorsSpecies.append(speciesColorDict[strainID])
        vcolorsSpecies.append(speciesColorDict[strainID])
        markerColorsSpecies.append(speciesColorDict[strainID])
        axes[1].scatter(treeAbundances[treeAbundances[tree_strain_id] == strainID]['t'],
                        treeAbundances[treeAbundances[tree_strain_id]
                                       == strainID][tree_strain_id],
                        lw=.5,
                        s=treeAbundances[treeAbundances[tree_strain_id] ==
                                         strainID]['abundance'].values*markerIncSpecies,
                        color=speciesColorDict[strainID], marker='|')
    strainLineages = LineCollection(
        hlinecSpecies, linestyles='solid', colors=hcolorsSpecies, linewidths=(1))
    creationLines = LineCollection(
        vlinecSpecies, linestyles='solid', colors=vcolorsSpecies, linewidths=(1))
    axes[1].add_collection(strainLineages)
    axes[1].add_collection(creationLines)
    axes[1].set_xlabel(xlabel='Time t', fontsize=7)
    axes[1].set_ylim(0, np.max(strainTimes[tree_strain_id])+1)
    axes[1].set_xlim(0, t)
    axes[1].xaxis.set_minor_locator(ticker.MultipleLocator(25))
    axes[1].set_yticks([])
    species_stacked = species_stacked.merge(keepTreeStrainsDF, on=[strain_id])\
        .drop(columns=[strain_id]).replace({tree_strain_id: newTreeOrder})
    species_stacked.sort_values(by=['t', tree_strain_id])
    speciesDF = species_stacked.copy()
    species_stacked = species_stacked.pivot(
        index='t', columns=tree_strain_id, values='abundance')
    if (overlay) & (species == 'virus'):
        microbe_total = pd.read_sql_query("SELECT t, microbial_abundance FROM summary WHERE run_id = {}".format(run_id), conSim)\
            .rename(columns={'microbial_abundance': 'btotal'})
        microbe_total = microbe_total[microbe_total.t <= max(
            species_stacked.index)]
        if stacked:
            axes.append(axes[0].twinx())
            axes[0].plot(microbe_total['t'], microbe_total['btotal'],
                         linewidth=0, color='grey')
            axes[0].fill_between(
                microbe_total['t'], microbe_total['btotal'], color='grey', alpha=0.6)
            species_stacked.plot.area(ax=axes[3], stacked=True, legend=False,
                                      linewidth=0, color=speciesColorDict, sort_columns=True)
            species_stacked.plot(stacked=True, ax=axes[3], legend=False,
                                 color='white', sort_columns=True, linewidth=.1)
            axes[3].set_ylabel(ylabel=abundanceTitle,
                               labelpad=30, fontsize=10, rotation=270)
            species_stacked.plot.area(ax=axes[2], stacked=True, legend=False,
                                      linewidth=0, color=speciesColorDict, sort_columns=True, alpha=0.25)
            species_stacked.plot(stacked=True, ax=axes[2], legend=False,
                                 color='white', sort_columns=True, linewidth=.1)
            axes.append(axes[2].twinx())
            axes.append(axes[2].twinx())
            axes[2].set_yticklabels([])
            axes[2].set_yticks([])
            axes[2].margins(x=0)
            axes[4].plot(speciesRichness['t'],speciesRichness['richness'],color='darkred',linewidth=1.5)
            axes[4].yaxis.tick_right()
            axes[4].yaxis.set_label_position("right")
            axes[4].set_ylabel(ylabel ='Viral Richness',labelpad=30,fontsize=10,rotation=270)
            axes[5].plot(speciesShannon['t'], speciesShannon['shannon'],
                         color='darkblue', linewidth=1.5)
            axes[5].yaxis.tick_left()
            axes[5].yaxis.set_label_position("left")
            axes[5].set_ylabel(ylabel='Viral\n Shannon Diversity', labelpad=15,fontsize=10)
        else:
            axes.append(axes[0].twinx())
            axes[0].plot(microbe_total['t'], microbe_total['btotal'],
                         linewidth=0, color='grey')
            axes[0].fill_between(
                microbe_total['t'], microbe_total['btotal'], color='grey', alpha=0.5)
            species_stacked.plot(ax=axes[3], stacked=False, legend=False,
                                 linewidth=1.5, color=speciesColorDict, sort_columns=True)
            axes[3].set_ylabel(ylabel='Viral Strain\nAbundances',
                               labelpad=30, fontsize=10, rotation=270)
            species_stacked.plot.area(ax=axes[2], stacked=True, legend=False,
                                      linewidth=0, color=speciesColorDict, sort_columns=True, alpha=0.25)
            species_stacked.plot(stacked=True, ax=axes[2], legend=False,
                                 color='white', sort_columns=True, linewidth=.1)
            axes.append(axes[2].twinx())
            axes.append(axes[2].twinx())
            axes[2].set_yticklabels([])
            axes[2].set_yticks([])
            axes[2].margins(x=0)
            axes[4].plot(speciesRichness['t'],speciesRichness['richness'],color='darkred',linewidth=1.5)
            axes[4].yaxis.tick_right()
            axes[4].yaxis.set_label_position("right")
            axes[4].set_ylabel(ylabel ='Viral Richness',labelpad=30,fontsize=10,rotation=270)
            axes[5].plot(speciesShannon['t'], speciesShannon['shannon'],
                         color='darkblue', linewidth=1.5)
            axes[5].yaxis.tick_left()
            axes[5].yaxis.set_label_position("left")
            axes[5].set_ylabel(ylabel='Viral\n Shannon Diversity', labelpad=15,fontsize=10)
        axes[1].set_ylabel(ylabel='Viral Strains', labelpad=15, fontsize=10)
        axes[0].set_ylabel(ylabel='Host\nAbundance', labelpad=15, fontsize=10)
        lim = axes[3].get_ylim()
        axes[3].set_ylim(0, lim[1])
    if (overlay) & (species == 'microbe'):
        virus_total = pd.read_sql_query("SELECT t, viral_abundance FROM summary WHERE run_id = {}".format(run_id), conSim)\
            .rename(columns={'viral_abundance': 'vtotal'})
        virus_total = virus_total[virus_total.t <= max(species_stacked.index)]
        if stacked:
            axes.append(axes[0].twinx())
            axes[3].plot(virus_total['t'], virus_total['vtotal'],
                         linewidth=0, color='grey')
            axes[3].fill_between(
                virus_total['t'], virus_total['vtotal'], color='grey', alpha=0.75)
            species_stacked.plot.area(ax=axes[0], stacked=True, legend=False,
                                      linewidth=0, color=speciesColorDict, sort_columns=True)
            species_stacked.plot(stacked=True, ax=axes[0], legend=False,
                                 color='white', sort_columns=True, linewidth=.1)
            axes[0].set_ylabel(ylabel=abundanceTitle, labelpad=15, fontsize=10)
            species_stacked.plot.area(ax=axes[2], stacked=True, legend=False,
                                      linewidth=0, color=speciesColorDict, sort_columns=True, alpha=0.25)
            species_stacked.plot(stacked=True, ax=axes[2], legend=False,
                                 color='white', sort_columns=True, linewidth=.1)
            axes.append(axes[2].twinx())
            axes.append(axes[2].twinx())
            axes[2].set_yticklabels([])
            axes[2].set_yticks([])
            axes[2].margins(x=0)
            axes[4].plot(speciesRichness['t'],speciesRichness['richness'],color='darkred',linewidth=1.5)
            axes[4].yaxis.tick_right()
            axes[4].yaxis.set_label_position("right")
            axes[4].set_ylabel(ylabel ='Host Richness',labelpad=30,fontsize=10,rotation=270)
            axes[5].plot(speciesShannon['t'], speciesShannon['shannon'],
                         color='darkblue', linewidth=1.5)
            axes[5].yaxis.tick_left()
            axes[5].yaxis.set_label_position("left")
            axes[5].set_ylabel(ylabel='Host\n Shannon Diversity', labelpad=15,fontsize=10)
        else:
            axes.append(axes[0].twinx())
            axes[3].plot(virus_total['t'], virus_total['vtotal'],
                         linewidth=0, color='grey')
            axes[3].fill_between(
                virus_total['t'], virus_total['vtotal'], color='grey', alpha=0.5)
            species_stacked.plot(ax=axes[0], stacked=False, legend=False,
                                 linewidth=1.5, color=speciesColorDict, sort_columns=True)
            axes[0].set_ylabel(
                ylabel='Host Strain\nAbundances', labelpad=15, fontsize=10)
            species_stacked.plot.area(ax=axes[2], stacked=True, legend=False,
                                      linewidth=0, color=speciesColorDict, sort_columns=True, alpha=0.25)
            species_stacked.plot(stacked=True, ax=axes[2], legend=False,
                                 color='white', sort_columns=True, linewidth=.1)
            axes.append(axes[2].twinx())
            axes.append(axes[2].twinx())
            axes[2].set_yticklabels([])
            axes[2].set_yticks([])
            axes[2].margins(x=0)
            axes[4].plot(speciesRichness['t'],speciesRichness['richness'],color='darkred',linewidth=1.5)
            axes[4].yaxis.tick_right()
            axes[4].yaxis.set_label_position("right")
            axes[4].set_ylabel(ylabel ='Host Richness',labelpad=30,fontsize=10,rotation=270)
            axes[5].plot(speciesShannon['t'], speciesShannon['shannon'],
                         color='darkblue', linewidth=1.5)
            axes[5].yaxis.tick_left()
            axes[5].yaxis.set_label_position("left")
            axes[5].set_ylabel(ylabel='Host\n Shannon Diversity', labelpad=15,fontsize=10)
        axes[1].set_ylabel(ylabel='Host Immune Strains', labelpad=15, fontsize=10)
        axes[3].set_ylabel(ylabel='Viral\nAbundance',
                           labelpad=30, fontsize=10, rotation=270)
        lim = axes[5].get_ylim()
        axes[5].ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
    if (not overlay) & (stacked):
        species_stacked.plot.area(ax=axes[0], stacked=True, legend=False,
                                  linewidth=0, color=speciesColorDict, sort_columns=True)
        species_stacked.plot(
            stacked=True, ax=axes[0], legend=False, color='white', sort_columns=True, linewidth=.1)
        axes[0].set_ylabel(ylabel=abundanceTitle, labelpad=15, fontsize=10)
    if (not overlay) & (not stacked):
        species_stacked.plot(ax=axes[0], stacked=False, legend=False,
                             linewidth=1, color=speciesColorDict, sort_columns=True)
        if species == 'virus':
            straintitle = 'Viral Strain\nAbundances'
        else:
            straintitle = 'Host Strain\nAbundances'
        axes[0].set_ylabel(ylabel=straintitle, labelpad=15, fontsize=10)
    lim = axes[0].get_ylim()
    axes[0].set_ylim(0, lim[1])
    lim = axes[2].get_ylim()
    axes[2].set_ylim(0, lim[1])
    axes[0].set_xlim(0, t)
    axes[1].set_xlim(0, t)
    axes[2].set_xlim(0, t)
    axes[0].xaxis.set_minor_locator(ticker.MultipleLocator(25))
    axes[1].xaxis.set_minor_locator(ticker.MultipleLocator(25))
    axes[2].xaxis.set_minor_locator(ticker.MultipleLocator(25))
    axes[2].set_xlabel(xlabel='Time t', fontsize=10)
    axes[0].ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
    axes[1].ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
    lim = axes[3].get_ylim()
    axes[3].set_ylim(0, lim[1])
    lim = axes[4].get_ylim()
    axes[4].set_ylim(0, lim[1])
    lim = axes[5].get_ylim()
    axes[5].set_ylim(0, lim[1])
    axes[3].set_xlim(0, t)
    axes[4].set_xlim(0, t)
    axes[3].xaxis.set_minor_locator(ticker.MultipleLocator(25))
    axes[4].xaxis.set_minor_locator(ticker.MultipleLocator(25))
    axes[5].xaxis.set_minor_locator(ticker.MultipleLocator(25))
    axes[3].ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
    axes[4].ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
    keepTreeStrainsDF['new_tree_{}_id'.format(
        strain)] = keepTreeStrainsDF['tree_{}_id'.format(strain)].copy()
    keepTreeStrainsDF = keepTreeStrainsDF.replace(
        {'new_tree_{}_id'.format(strain): newTreeOrder})
    return speciesDF, keepTreeStrainsDF, speciesColorDict, hlinecSpecies, vlinecSpecies, fig, axes




def speciesTreePlot(run_id,species,DBSIM_PATH,DBTREE_PATH,treepalette,maxticksize,figxy,hratio,abundthreshold,stacked,overlay):
    if species == 'virus':
        strains = 'vstrains'
        strain = 'vstrain'
        strain_id = 'vstrain_id'
        abundanceTitle = 'Viral\nAbundance'
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
        abundanceTitle = 'Host\nAbundance'
        strain_abundance = "microbial_abundance"
        straintotal = "btotal"
        abundance = 'babundance'
        tree_strain_id = "tree_bstrain_id"
        tree_parent_strain_id = "tree_parent_bstrain_id"
        parent_strain_id = "parent_bstrain_id"
    conSim = sqlite3.connect(DBSIM_PATH)
    curSim = conSim.cursor()
    ID = curSim.execute(
        'SELECT combo_id,replicate FROM runs WHERE run_id = {}'.format(run_id)).fetchall()
    combo_id = ID[0][0]
    replicate = ID[0][1]
    RUN_DIR = os.path.join(
        'runID{0}-c{1}-r{2}'.format(run_id, combo_id, replicate))
    species_stacked = pd.read_sql_query(
        "SELECT t,{0},abundance FROM {1} WHERE run_id = {2}".format(strain_id, abundance, run_id), conSim)
    t = max(pd.read_sql_query(
        "SELECT t FROM vabundance WHERE run_id = {}".format(run_id), conSim).t.values)
    species_stacked = species_stacked[species_stacked.t <= t]
    species_total = pd.read_sql_query("SELECT t,{0} FROM summary WHERE run_id = {1}".format(strain_abundance, run_id), conSim)\
        .rename(columns={strain_abundance: straintotal})
    species_total = species_total[species_total.t <= max(species_stacked.t)]
    maxIDs = species_stacked.set_index('t').groupby([strain_id]).agg(t=('abundance', 'idxmax'),
                                                                     maxAbund=('abundance', 'max')).reset_index()
    maxIDs = species_total.merge(maxIDs, on=['t'])
    maxIDs[straintotal] = abundthreshold*np.array(maxIDs[straintotal])
    keepStrains = list(
        maxIDs[maxIDs['maxAbund'] > maxIDs[straintotal]][strain_id].values)
    conTree = sqlite3.connect(DBTREE_PATH)
    curTree = conTree.cursor()
    print('SQLite Query: tree data')
    parentTrack = keepStrains.copy()
    keepStrainsNew = keepStrains.copy()
    while len(parentTrack) != 0:
        parentTrack = list(pd.read_sql_query(
            "SELECT parent_{1} \
        FROM {0} WHERE {1} in ({2})".format(strains, strain_id, ', '.join(map(str, parentTrack))), conSim)
            [parent_strain_id].values)
        # print(parentTrack)
        keepStrainsNew.extend(parentTrack)
        keepStrainsNew = list(np.unique(keepStrainsNew))
        # print(keepStrainsNew)
        parentTrack = list(np.array(parentTrack)[np.array(parentTrack) != 0])
    keepStrainsNew = list(*np.array(keepStrainsNew)[keepStrainsNew != 0])
    keepStrainsNew.sort()
    species_stacked = species_stacked[[
        (i in keepStrainsNew) for i in species_stacked[strain_id]]].reset_index(drop=True)
    keepTreeStrainsDF = pd.read_sql_query(
        "SELECT tree_{0}, {0} \
    FROM tree_{1}_order WHERE {0} in ({2})".format(strain_id, strain, ', '.join(map(str, keepStrainsNew))), conTree)
    keepTreeStrains = list(
        keepTreeStrainsDF['tree_{}'.format(strain_id)].values)
    keepTreeStrains.sort()
    strainTimes = pd.read_sql_query(
        "SELECT tree_{0}, t_creation, t_extinction, tree_parent_{0} \
    FROM tree_{1}_creation_extinction WHERE tree_{0} in ({2})".format(
            strain_id, strain, ', '.join(map(str, keepTreeStrains))), conTree)
    treeAbundances = pd.read_sql_query(
        "SELECT t, tree_{0}, abundance \
    FROM tree_{1} WHERE tree_{0} in ({2})"
        .format(strain_id, abundance, ', '.join(map(str, keepTreeStrains))), conTree)
    newTreeOrder = {}
    newTreeID = 0
    for treeStrain in keepTreeStrains:
        newTreeOrder[treeStrain] = newTreeID
        newTreeID += 1
    treeAbundances.replace({tree_strain_id: newTreeOrder}, inplace=True)
    strainTimes.replace({tree_strain_id: newTreeOrder,
                        tree_parent_strain_id: newTreeOrder}, inplace=True)
    numStrains = len(keepTreeStrains)
    Vcmap = cm.get_cmap(treepalette)
    # Vcmap = sns.color_palette("icefire",as_cmap=True)
    norm = Normalize(vmin=float(1), vmax=float(max(newTreeOrder.values())))
    speciesColorDict = {}
    maxAbundanceSpecies = np.max(treeAbundances.abundance.values)
    markerIncSpecies = maxticksize/maxAbundanceSpecies
    markerColorsSpecies = []
    hlinecSpecies = []
    vlinecSpecies = []
    hcolorsSpecies = []
    vcolorsSpecies = []
    print('Compiling stacked species abundances and tree plots')
    fig, ax = plt.subplots(2, sharex=True, figsize=figxy,
                           gridspec_kw={'height_ratios': hratio})
    # fig.suptitle('Strain Abundances (run{0}-c{1}-r{2})'.format(run_id,combo_id,replicate))
    axes = [ax[0], ax[1]]
    for strainID in sorted(treeAbundances[tree_strain_id].unique()):
        tCreate = strainTimes[strainTimes[tree_strain_id]
                              == strainID].t_creation.values[0]
        tExtinct = strainTimes[strainTimes[tree_strain_id]
                               == strainID].t_extinction.values[0]
        parent = strainTimes[strainTimes[tree_strain_id]
                             == strainID][tree_parent_strain_id].values[0]
        hlinecSpecies.append([[tCreate, strainID], [tExtinct, strainID]])
        if parent != 0:
            vlinecSpecies.append([[tCreate, parent], [tCreate, strainID]])
        speciesColorDict[strainID] = Vcmap(norm(strainID))
        hcolorsSpecies.append(speciesColorDict[strainID])
        vcolorsSpecies.append(speciesColorDict[strainID])
        markerColorsSpecies.append(speciesColorDict[strainID])
        axes[1].scatter(treeAbundances[treeAbundances[tree_strain_id] == strainID]['t'],
                        treeAbundances[treeAbundances[tree_strain_id]
                                       == strainID][tree_strain_id],
                        lw=.5,
                        s=treeAbundances[treeAbundances[tree_strain_id] ==
                                         strainID]['abundance'].values*markerIncSpecies,
                        color=speciesColorDict[strainID], marker='|')
    strainLineages = LineCollection(
        hlinecSpecies, linestyles='solid', colors=hcolorsSpecies, linewidths=(1))
    creationLines = LineCollection(
        vlinecSpecies, linestyles='solid', colors=vcolorsSpecies, linewidths=(1))
    axes[1].add_collection(strainLineages)
    axes[1].add_collection(creationLines)
    axes[1].set_xlabel(xlabel='Time t', fontsize=7)
    axes[1].set_ylim(0, np.max(strainTimes[tree_strain_id])+1)
    axes[1].set_xlim(0, t)
    axes[1].xaxis.set_minor_locator(ticker.MultipleLocator(25))
    axes[1].set_yticks([])
    species_stacked = species_stacked.merge(keepTreeStrainsDF, on=[strain_id])\
        .drop(columns=[strain_id]).replace({tree_strain_id: newTreeOrder})
    species_stacked.sort_values(by=['t', tree_strain_id])
    speciesDF = species_stacked.copy()
    species_stacked = species_stacked.pivot(
        index='t', columns=tree_strain_id, values='abundance')
    if (overlay) & (species == 'virus'):
        microbe_total = pd.read_sql_query("SELECT t, microbial_abundance FROM summary WHERE run_id = {}".format(run_id), conSim)\
            .rename(columns={'microbial_abundance': 'btotal'})
        microbe_total = microbe_total[microbe_total.t <= max(species_stacked.index)]
        if stacked:
            axes.append(axes[0].twinx())
            axes[0].plot(microbe_total['t'],microbe_total['btotal'],linewidth=0,color='grey')
            axes[0].fill_between(microbe_total['t'],microbe_total['btotal'], color='grey',alpha=0.6)
            species_stacked.plot.area(ax=axes[2], stacked=True, legend=False,
                                      linewidth=0, color=speciesColorDict, sort_columns=True)
            species_stacked.plot(stacked=True, ax=axes[2], legend=False, \
                                    color='white', sort_columns=True, linewidth=.1)
            axes[2].set_ylabel(ylabel=abundanceTitle, labelpad=15, fontsize=10,rotation=270)
        else:
            axes.append(axes[0].twinx())
            axes[0].plot(microbe_total['t'],microbe_total['btotal'],linewidth=0,color='grey')
            axes[0].fill_between(microbe_total['t'],microbe_total['btotal'], color='grey',alpha=0.5)
            species_stacked.plot(ax=axes[2], stacked=False, legend=False,
                                 linewidth=1.5, color=speciesColorDict, sort_columns=True)
            axes[2].set_ylabel(ylabel='Viral Strain\nAbundances', labelpad=30, fontsize=10,rotation=270)
        axes[0].set_ylabel(ylabel='Host\nAbundance', labelpad=15, fontsize=10)
        lim = axes[2].get_ylim()
        axes[2].set_ylim(0, lim[1])
    if (overlay) & (species == 'microbe'):
        virus_total = pd.read_sql_query("SELECT t, viral_abundance FROM summary WHERE run_id = {}".format(run_id), conSim)\
            .rename(columns={'viral_abundance': 'vtotal'})
        virus_total = virus_total[virus_total.t <= max(species_stacked.index)]
        if stacked:
            axes.append(axes[0].twinx())
            axes[2].plot(virus_total['t'],virus_total['vtotal'],linewidth=0,color='grey')
            axes[2].fill_between(virus_total['t'],virus_total['vtotal'], color='grey',alpha=0.75)
            species_stacked.plot.area(ax=axes[0], stacked=True, legend=False,
                                      linewidth=0, color=speciesColorDict, sort_columns=True)
            species_stacked.plot(stacked=True, ax=axes[0], legend=False, \
                                    color='white', sort_columns=True, linewidth=.1)
            axes[0].set_ylabel(ylabel=abundanceTitle, labelpad=15, fontsize=10)
        else:
            axes.append(axes[0].twinx())
            axes[2].plot(virus_total['t'],virus_total['vtotal'],linewidth=0,color='grey')
            axes[2].fill_between(virus_total['t'],virus_total['vtotal'], color='grey',alpha=0.5)
            species_stacked.plot(ax=axes[0], stacked=False, legend=False,
                                  linewidth=1.5, color=speciesColorDict, sort_columns=True)
            axes[0].set_ylabel(ylabel='Host Strain\nAbundances', labelpad=15, fontsize=10)
        axes[2].set_ylabel(ylabel='Viral\nAbundance', labelpad=30, fontsize=10,rotation=270)
        lim = axes[2].get_ylim()
        axes[2].set_ylim(0, lim[1])
    if (not overlay) & (stacked):
        species_stacked.plot.area(ax=axes[0], stacked=True, legend=False,
                                linewidth=0, color=speciesColorDict, sort_columns=True)
        species_stacked.plot(
            stacked=True, ax=axes[0], legend=False, color='white', sort_columns=True, linewidth=.1)
        axes[0].set_ylabel(ylabel=abundanceTitle, labelpad=15, fontsize=10)
    if (not overlay) & (not stacked):
        species_stacked.plot(ax=axes[0], stacked=False, legend=False,
                                  linewidth=1, color=speciesColorDict, sort_columns=True)
        if species == 'virus':
            straintitle = 'Viral Strain\nAbundances'
        else:
            straintitle = 'Host Strain\nAbundances'
        axes[0].set_ylabel(ylabel=straintitle, labelpad=15, fontsize=10)
    lim = axes[0].get_ylim()
    axes[0].set_ylim(0, lim[1])
    axes[1].set_xlabel(xlabel='Time t', fontsize=10)
    axes[0].ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
    axes[0].xaxis.set_minor_locator(ticker.MultipleLocator(25))
    axes[1].ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
    axes[1].xaxis.set_minor_locator(ticker.MultipleLocator(25))
    axes[0].set_xlim(0, t)
    keepTreeStrainsDF['new_tree_{}_id'.format(
        strain)] = keepTreeStrainsDF['tree_{}_id'.format(strain)].copy()
    keepTreeStrainsDF = keepTreeStrainsDF.replace(
        {'new_tree_{}_id'.format(strain): newTreeOrder})
    return speciesDF, keepTreeStrainsDF, speciesColorDict, hlinecSpecies, vlinecSpecies, fig, axes


def speciesCladeTreePlot(run_id,species,colorSep,DBSIM_PATH,DBCLADE_PATH,DBTREE_PATH,treepalette,maxticksize,figxy,hratio,abundthreshold):
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
    conClade = sqlite3.connect(DBCLADE_PATH)
    curClade = conClade.cursor()
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
    species_stacked = species_stacked[[(i in keepStrainsNew) for i in species_stacked[strain_id]]].reset_index(drop=True)
    print('SQLite Query: clade data')
    speciesClades = pd.read_sql_query("SELECT DISTINCT clade_id, {0} \
    FROM {1}s WHERE {0} in ({2})".format(strain_id,abundance,', '.join(map(str,keepStrainsNew))), conClade)
    speciesCladeIDs = list(np.unique(speciesClades['clade_id'].copy().values))
    speciesCladeIDs.sort()
    keepTreeStrainsDF = pd.read_sql_query(
    "SELECT tree_{0}, {0} \
    FROM tree_{1}_order WHERE {0} in ({2})".format(strain_id,strain,', '.join(map(str,keepStrainsNew))), conTree)
    keepTreeStrains = list(keepTreeStrainsDF['tree_{}'.format(strain_id)].values)
    strainTimes = pd.read_sql_query(
    "SELECT tree_{0}, t_creation, t_extinction, tree_parent_{0} \
    FROM tree_{1}_creation_extinction WHERE tree_{0} in ({2})".format(
    strain_id,strain,', '.join(map(str,keepTreeStrains))), conTree)
    treeAbundances = pd.read_sql_query(
    "SELECT t, tree_{0}, abundance \
    FROM tree_{1} WHERE tree_{0} in ({2})"\
    .format(strain_id,abundance,', '.join(map(str,keepTreeStrains))), conTree)
    Vcmap = cm.get_cmap(treepalette)
    # Vcmap = sns.color_palette("icefire",as_cmap=True)
    maxCladeColorID = len(speciesCladeIDs) + 1
    speciesColorDict = {}
    cladeColorDict = {}
    colorInd = 0
    for clade_id in speciesCladeIDs:
        cladeColorDict[clade_id] = colorInd
        colorInd += 1
    print('Compiling stacked species clade abundances and tree plots')
    fig, ax = plt.subplots(2,sharex=True, figsize=figxy, gridspec_kw={'height_ratios': hratio})
    # fig.suptitle('Strain Abundances (run{0}-c{1}-r{2})'.format(run_id,combo_id,replicate))
    axes = [ax[0], ax[1]]
    hlinecSpecies = []
    vlinecSpecies = []
    hcolorsSpecies = []
    vcolorsSpecies = []
    maxAbundanceSpecies = np.max(treeAbundances.abundance.values)
    markerIncSpecies = maxticksize/maxAbundanceSpecies
    markerColorsSpecies = []
    treeStrainsAll = []
    newTreeOrder = {}
    for clade_id in speciesCladeIDs:
        strainsOfClade = curClade.execute("SELECT DISTINCT {0} FROM {1}s \
        WHERE clade_id = {2} AND {0} in ({3})".format(strain_id,abundance,clade_id,', '.join(map(str,keepStrainsNew))))
        strainsOfClade = [i[0] for i in strainsOfClade.fetchall()]
        treeStrainsOfClade = curTree.execute("SELECT DISTINCT tree_{0} FROM tree_{1}_order \
        WHERE {0} in ({2}) ORDER BY tree_{0}".format(strain_id,strain,', '.join(map(str,strainsOfClade))))
        treeStrainsOfClade = [i[0] for i in treeStrainsOfClade.fetchall()]
        treeStrainsAll.extend(treeStrainsOfClade)
    treeStrainsAll.sort()
    newTreeID = 1
    for treeStrain in treeStrainsAll:
        newTreeOrder[treeStrain] = newTreeID
        newTreeID += 1
    treeAbundances.replace({tree_strain_id:newTreeOrder},inplace=True)
    strainTimes.replace({tree_strain_id:newTreeOrder,tree_parent_strain_id:newTreeOrder},inplace=True)
    for clade_id in speciesCladeIDs:
        strainsOfClade = curClade.execute("SELECT DISTINCT {0} FROM {1}s \
        WHERE clade_id = {2} AND {0} in ({3})".format(strain_id,abundance,clade_id,', '.join(map(str,keepStrainsNew))))
        strainsOfClade = [i[0] for i in strainsOfClade.fetchall()]
        treeStrainsOfClade = curTree.execute("SELECT DISTINCT tree_{0} FROM tree_{1}_order \
        WHERE {0} in ({2}) ORDER BY tree_{0}".format(strain_id,strain,', '.join(map(str,strainsOfClade))))
        treeStrainsOfClade = [newTreeOrder[i[0]] for i in treeStrainsOfClade.fetchall()]
        numStrains = len(treeStrainsOfClade)
        colorID = cladeColorDict[clade_id]
        c1 = colorID/maxCladeColorID
        c2 = (colorID+1)/maxCladeColorID
        strainColorID = 0
        for strainID in treeStrainsOfClade[::-1]:
            tCreate = strainTimes[strainTimes[tree_strain_id] == strainID].t_creation.values[0]
            tExtinct = strainTimes[strainTimes[tree_strain_id] == strainID].t_extinction.values[0]
            parent = strainTimes[strainTimes[tree_strain_id] == strainID][tree_parent_strain_id].values[0]
            hlinecSpecies.append([[tCreate, strainID],[tExtinct, strainID]])
            vlinecSpecies.append([[tCreate, parent],[tCreate, strainID]])
            speciesColorDict[strainID] = Vcmap(c1 + (c2-c1)/(numStrains+colorSep)*strainColorID)
            print(c1 + (c2-c1)/(numStrains)*strainColorID)
            print(Vcmap(c1 + (c2-c1)/(numStrains)*strainColorID))
            hcolorsSpecies.append(speciesColorDict[strainID])
            vcolorsSpecies.append(speciesColorDict[strainID])
            markerColorsSpecies.append(speciesColorDict[strainID])
            strainColorID += 1
            axes[1].scatter(treeAbundances[treeAbundances[tree_strain_id]==strainID]['t'],
            treeAbundances[treeAbundances[tree_strain_id]==strainID][tree_strain_id],
            lw=.5,
            s=treeAbundances[treeAbundances[tree_strain_id]==strainID]['abundance'].values*markerIncSpecies,
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
    species_stacked = species_stacked.pivot(index='t',columns=tree_strain_id,values='abundance')
    species_stacked.plot.area(ax = axes[0],stacked=True,legend=False, linewidth=0,color=speciesColorDict,sort_columns=True)
    species_stacked.plot(stacked=True, ax=axes[0], legend=False, color='black',sort_columns=True,linewidth=.1)
    axes[0].set_ylabel(ylabel =abundanceTitle,labelpad=15,fontsize=15)
    axes[0].set_xlabel(xlabel = 'Time t',fontsize=15)
    axes[0].ticklabel_format(axis = 'y',style='sci',scilimits=(0,0))
    axes[0].xaxis.set_minor_locator(ticker.MultipleLocator(25))
    lim = axes[0].get_ylim()
    axes[0].set_ylim(0,lim[1])
    axes[0].set_xlim(0,t)
    keepTreeStrainsDF['new_tree_{}_id'.format(strain)] = keepTreeStrainsDF['tree_{}_id'.format(strain)].copy()
    keepTreeStrainsDF = keepTreeStrainsDF.replace({'new_tree_{}_id'.format(strain): newTreeOrder})
    return speciesDF, keepTreeStrainsDF, speciesColorDict, cladeColorDict, hlinecSpecies, vlinecSpecies, fig, axes


# def susceptibleCladeTreePlot(run_id,DBSIM_PATH,DBCLADE_PATH,DBTREE_PATH,DBMATCH_PATH,colorpalClades,maxTickSizeMicrobe,figxy,hratio):
#     conSim = sqlite3.connect(DBSIM_PATH)
#     curSim = conSim.cursor()
#     ID = curSim.execute('SELECT combo_id,replicate FROM runs WHERE run_id = {}'.format(run_id)).fetchall()
#     combo_id = ID[0][0]
#     replicate = ID[0][1]
#     RUN_DIR = os.path.join('runID{0}-c{1}-r{2}'.format(run_id,combo_id,replicate))
#     microbe_stacked = pd.read_sql_query("SELECT t,bstrain_id,abundance FROM babundance WHERE run_id = {}".format(run_id), conSim)
#     virus_stacked = pd.read_sql_query("SELECT t FROM vabundance WHERE run_id = {}".format(run_id), conSim)
#     microbe_stacked = microbe_stacked[microbe_stacked.t <= max(virus_stacked.t)]
#     microbe_stacked = microbe_stacked.pivot(index='t',columns='bstrain_id',values='abundance')
#     conClade = sqlite3.connect(DBCLADE_PATH)
#     curClade = conClade.cursor()
#     conTree = sqlite3.connect(DBTREE_PATH)
#     curTree = conTree.cursor()
#     conMatch = sqlite3.connect(DBMATCH_PATH)
#     curMatch = conMatch.cursor()
#     print('SQLite Query: clade data')
#     microbeClades = pd.read_sql_query("SELECT DISTINCT clade_id, bstrain_id \
#     FROM babundances", conClade)
#     microbeCladeIDs = pd.read_sql_query("SELECT DISTINCT clade_id \
#     FROM babundances", conClade)
#     print('SQLite Query: tree data')
#     bStrainTimes = pd.read_sql_query(
#     "SELECT tree_bstrain_id, t_creation, t_extinction, tree_parent_bstrain_id \
#     FROM tree_bstrain_creation_extinction", conTree)
#     bstrain0vstrain = pd.read_sql_query(
#     "SELECT t, bstrain_id, vstrain_id \
#     FROM bstrain_to_vstrain_0matches", conMatch)
#     vAbunds = pd.read_sql_query("SELECT t,vstrain_id,abundance FROM vabundance WHERE run_id = {}".format(run_id), conSim)
#     bstrain0vstrain = bstrain0vstrain.merge(vAbunds,on=['t','vstrain_id'])\
#     .groupby(['t','bstrain_id']).agg(abundance=('abundance','sum')).reset_index()
#     bTreeAbundances = pd.read_sql_query("SELECT tree_bstrain_id,bstrain_id \
#     FROM tree_bstrain_order \
#     ORDER BY bstrain_id",conTree)
#     bTreeAbundances = bstrain0vstrain.merge(bTreeAbundances,on=['bstrain_id']).drop(columns=['bstrain_id'])
#     vTreeAbundances = pd.read_sql_query(
#     "SELECT t \
#     FROM tree_vabundance", conTree)
#     if bTreeAbundances['t'][bTreeAbundances['t'].size-1] not in vTreeAbundances.t.values:
#         bTreeAbundances = bTreeAbundances[bTreeAbundances.t != bTreeAbundances['t'][bTreeAbundances['t'].size-1]]
#     Mcmap = cm.get_cmap('{}'.format(colorpalClades)) #usually cmap = 'turbo'
#     # Get normalize function (takes data in range [vmin, vmax] -> [0, 1])
#     Mnorm = Normalize(vmin=1, vmax=len(microbeCladeIDs))
#     microbeColorDict = {}
#     cladeColorDict = {}
#     cladeDict = {}
#     colorInd = 0
#     for clade_id in microbeCladeIDs.values:
#         # print(id[0])
#         cladeColorDict[clade_id[0]] = colorInd
#         cladeDict[clade_id[0]] = []
#         colorInd += 1
#     for strain in microbeClades.bstrain_id.values:
#         clade = microbeClades[microbeClades['bstrain_id']==strain]['clade_id'].values[0]
#         cladeDict[clade] = np.append(cladeDict[clade],strain)
#     columnOrder = []
#     for clade_id in microbeCladeIDs.values[::-1]:
#         columnOrder = np.append(columnOrder,cladeDict[clade_id[0]])
#     columnOrder = columnOrder.astype(int)
#     microbe_stacked = microbe_stacked[columnOrder]
#     for strain in microbe_stacked.columns.values:
#         clade = microbeClades[microbeClades['bstrain_id']==strain]['clade_id'].values[0]
#         microbeColorDict[strain] = Mcmap(Mnorm(np.arange(1, len(microbeCladeIDs)+1, 1)))[cladeColorDict[clade]]
#     print('Compiling stacked microbe clade abundances and tree plots')
#     fig, ax = plt.subplots(2,sharex=True, figsize=figxy, gridspec_kw={'height_ratios': hratio})
#     fig.suptitle('Strain Abundances (run{0}-c{1}-r{2})'.format(run_id,combo_id,replicate))
#     axes = [ax[0], ax[1]]
#     microbe_stacked.plot.area(ax = axes[0],stacked=True,legend=False, linewidth=0,color=microbeColorDict,sort_columns=True)
#     microbe_stacked.plot(stacked=True, ax=axes[0], legend=False, color='black',sort_columns=True,linewidth=.1)
#     axes[0].set_ylabel(ylabel ='Microbial Strain Abundances',labelpad=15,fontsize=7)
#     axes[0].set_xlabel(xlabel = 'Time t',fontsize=7)
#     axes[0].ticklabel_format(axis = 'y',style='sci',scilimits=(0,0))
#     axes[0].xaxis.set_minor_locator(ticker.MultipleLocator(25))
#     lim = axes[0].get_ylim()
#     axes[0].set_ylim(0,lim[1])
#     axes[0].set_xlim(0,np.max(bTreeAbundances.t.values))
#     axes[1].set_xlabel(xlabel = 'Time t',fontsize=7)
#     axes[1].set_ylim(0,np.max(bStrainTimes.tree_bstrain_id)+1)
#     axes[0].set_xlim(0,np.max(bTreeAbundances.t.values))
#     axes[1].xaxis.set_minor_locator(ticker.MultipleLocator(25))
#     axes[1].set_yticks([])
#     hlinecMicrobe = []
#     vlinecMicrobe = []
#     hcolorsMicrobe = []
#     vcolorsMicrobe = []
#     maxAbundanceMicrobe = np.max(bTreeAbundances.abundance.values)
#     markerIncMicrobe = maxTickSizeMicrobe/maxAbundanceMicrobe
#     markerColorsMicrobe = []
#     for clade_id in np.unique(microbeClades.clade_id):
#         strainsOfClade = curClade.execute("SELECT DISTINCT bstrain_id FROM babundances \
#         WHERE clade_id = {}".format(clade_id))
#         strainsOfClade = [i[0] for i in strainsOfClade.fetchall()]
#         treeStrainsOfClade = curTree.execute("SELECT DISTINCT tree_bstrain_id FROM tree_bstrain_order \
#         WHERE bstrain_id in ({}) ORDER BY tree_bstrain_id".format(', '.join(map(str,strainsOfClade))))
#         treeStrainsOfClade = [i[0] for i in treeStrainsOfClade.fetchall()]
#         axes[1].scatter(bTreeAbundances[bTreeAbundances['tree_bstrain_id'].isin(treeStrainsOfClade)]['t'],
#         bTreeAbundances[bTreeAbundances['tree_bstrain_id'].isin(treeStrainsOfClade)]['tree_bstrain_id'],
#         lw=.8,
#         s=bTreeAbundances[bTreeAbundances['tree_bstrain_id'].isin(treeStrainsOfClade)]['abundance'].values*markerIncMicrobe,
#         c = np.array([Mcmap(Mnorm(np.arange(1, len(microbeCladeIDs)+1, 1)))[cladeColorDict[clade_id]]]),marker='|')
#         for strain in treeStrainsOfClade:
#             tCreate = bStrainTimes[bStrainTimes.tree_bstrain_id == strain].t_creation.values[0]
#             tExtinct = bStrainTimes[bStrainTimes.tree_bstrain_id == strain].t_extinction.values[0]
#             parent = bStrainTimes[bStrainTimes.tree_bstrain_id == strain].tree_parent_bstrain_id.values[0]
#             hlinecMicrobe.append([[tCreate, strain],[tExtinct, strain]])
#             vlinecMicrobe.append([[tCreate, parent],[tCreate, strain]])
#             hcolorsMicrobe.append(Mcmap(Mnorm(np.arange(1, len(microbeCladeIDs)+1, 1)))[cladeColorDict[clade_id]])
#             vcolorsMicrobe.append(Mcmap(Mnorm(np.arange(1, len(microbeCladeIDs)+1, 1)))[cladeColorDict[clade_id]])
#             markerColorsMicrobe.append(Mcmap(Mnorm(np.arange(1, len(microbeCladeIDs)+1, 1)))[cladeColorDict[clade_id]])
#     strainLineages = LineCollection(hlinecMicrobe, linestyles='solid', colors=hcolorsMicrobe,linewidths=(0.35))
#     creationLines = LineCollection(vlinecMicrobe, linestyles='solid', colors=vcolorsMicrobe, linewidths=(0.35))
#     axes[1].add_collection(strainLineages)
#     axes[1].add_collection(creationLines)
#     return fig, axes
