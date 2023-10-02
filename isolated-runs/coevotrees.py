def coEvoTreePlot(run_id,DBSIM_PATH,DBTREE_PATH,treepaletteV,treepaletteB,maxticksize,figxy,hratio,abundthresholdV,abundthresholdB):
    conSim = sqlite3.connect(DBSIM_PATH)
    curSim = conSim.cursor()
    ID = curSim.execute(
        'SELECT combo_id,replicate FROM runs WHERE run_id = {}'.format(run_id)).fetchall()
    combo_id = ID[0][0]
    replicate = ID[0][1]
    RUN_DIR = os.path.join(
        'runID{0}-c{1}-r{2}'.format(run_id, combo_id, replicate))
    conTree = sqlite3.connect(DBTREE_PATH)
    curTree = conTree.cursor()
    fig, ax = plt.subplots(4, sharex=True, figsize=figxy,
                        gridspec_kw={'height_ratios': hratio})
    # fig.suptitle('Strain Abundances (run{0}-c{1}-r{2})'.format(run_id,combo_id,replicate))
    axes = [ax[0], ax[1], ax[2], ax[3]]
    for species in ['microbe','virus']:
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
            treepalette = treepaletteV
            idx0 = 2
            idx1 = 3
            abundthreshold = abundthresholdV
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
            treepalette = treepaletteB
            idx0 = 0
            idx1 = 1
            abundthreshold = abundthresholdB
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
            axes[idx1].scatter(treeAbundances[treeAbundances[tree_strain_id] == strainID]['t'],
                            treeAbundances[treeAbundances[tree_strain_id]
                                        == strainID][tree_strain_id],
                            lw=1,
                            s=treeAbundances[treeAbundances[tree_strain_id] ==
                                            strainID]['abundance'].values*markerIncSpecies,
                            color=speciesColorDict[strainID], marker='|')
        strainLineages = LineCollection(
            hlinecSpecies, linestyles='solid', colors=hcolorsSpecies, linewidths=(1))
        creationLines = LineCollection(
            vlinecSpecies, linestyles='solid', colors=vcolorsSpecies, linewidths=(1))
        axes[idx1].add_collection(strainLineages)
        axes[idx1].add_collection(creationLines)
        axes[idx1].set_xlabel(xlabel='Time t', fontsize=7)
        axes[idx1].set_ylim(0, np.max(strainTimes[tree_strain_id])+1)
        axes[idx1].set_xlim(0, t)
        axes[idx1].xaxis.set_minor_locator(ticker.MultipleLocator(25))
        axes[idx1].set_yticks([])
        species_stacked = species_stacked.merge(keepTreeStrainsDF, on=[strain_id])\
            .drop(columns=[strain_id]).replace({tree_strain_id: newTreeOrder})
        species_stacked.sort_values(by=['t', tree_strain_id])
        speciesDF = species_stacked.copy()
        species_stacked = species_stacked.pivot(
            index='t', columns=tree_strain_id, values='abundance')
        species_stacked.plot.area(ax=axes[idx0], stacked=True, legend=False,
                                linewidth=0, color=speciesColorDict, sort_columns=True)
        # species_stacked.plot(
        #     stacked=True, ax=axes[idx0], legend=False, color='white', sort_columns=True, linewidth=.1)
        axes[idx0].set_ylabel(ylabel=abundanceTitle, labelpad=20, fontsize=20)
        lim = axes[idx0].get_ylim()
        axes[idx0].set_ylim(0, lim[1])
        axes[idx1].set_xlabel(xlabel='Time t', fontsize=20)
        axes[idx0].ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
        axes[idx0].xaxis.set_minor_locator(ticker.MultipleLocator(25))
        axes[idx1].ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
        axes[idx1].xaxis.set_minor_locator(ticker.MultipleLocator(25))
        axes[idx0].set_xlim(0, t)
        axes[idx0].tick_params(axis='x', labelsize=20)
        axes[idx1].tick_params(axis='x', labelsize=20)
        axes[idx0].tick_params(axis='y', labelsize=20)
        axes[idx0].yaxis.get_offset_text().set_fontsize(20)
    return fig, axes


fig.savefig(os.path.join('/Users/armun/Desktop/abundTree.pdf'),dpi=resolve)

coEvoTreePlot2(3297,DBSIM_PATH,DBTREE_PATH,'viridis','turbo',100,(15,15*1.3),[1,1,3],vabundthreshold,babundthreshold)

fig, axes = coEvoTreePlot2(3297,DBSIM_PATH,DBTREE_PATH,'bone','turbo',100,(15,15),[1,1,3],vabundthreshold,babundthreshold)

def coEvoTreePlot2(run_id, DBSIM_PATH, DBTREE_PATH, treepaletteV, treepaletteB, maxticksize, figxy, hratio, abundthresholdV, abundthresholdB):
    conSim = sqlite3.connect(DBSIM_PATH)
    curSim = conSim.cursor()
    ID = curSim.execute(
        'SELECT combo_id,replicate FROM runs WHERE run_id = {}'.format(run_id)).fetchall()
    combo_id = ID[0][0]
    replicate = ID[0][1]
    RUN_DIR = os.path.join(
        'runID{0}-c{1}-r{2}'.format(run_id, combo_id, replicate))
    conTree = sqlite3.connect(DBTREE_PATH)
    curTree = conTree.cursor()
    fig, ax = plt.subplots(3, sharex=True, figsize=figxy,
                           gridspec_kw={'height_ratios': hratio})
    # fig.suptitle('Strain Abundances (run{0}-c{1}-r{2})'.format(run_id,combo_id,replicate))
    axes = [ax[0], ax[1], ax[2]]
    for species in ['microbe', 'virus']:
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
            treepalette = treepaletteV
            idx0 = 1
            idx1 = 2
            abundthreshold = abundthresholdV
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
            treepalette = treepaletteB
            idx0 = 0
            abundthreshold = abundthresholdB
        species_stacked = pd.read_sql_query(
            "SELECT t,{0},abundance FROM {1} WHERE run_id = {2}".format(strain_id, abundance, run_id), conSim)
        t = max(pd.read_sql_query(
            "SELECT t FROM vabundance WHERE run_id = {}".format(run_id), conSim).t.values)
        species_stacked = species_stacked[species_stacked.t <= t]
        species_total = pd.read_sql_query("SELECT t,{0} FROM summary WHERE run_id = {1}".format(strain_abundance, run_id), conSim)\
            .rename(columns={strain_abundance: straintotal})
        species_total = species_total[species_total.t <= max(
            species_stacked.t)]
        maxIDs = species_stacked.set_index('t').groupby([strain_id]).agg(t=('abundance', 'idxmax'),
                                                                         maxAbund=('abundance', 'max')).reset_index()
        maxIDs = species_total.merge(maxIDs, on=['t'])
        maxIDs[straintotal] = abundthreshold*np.array(maxIDs[straintotal])
        keepStrains = list(
            maxIDs[maxIDs['maxAbund'] > maxIDs[straintotal]][strain_id].values)
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
            parentTrack = list(np.array(parentTrack)[
                            np.array(parentTrack) != 0])
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
        newTreeOrder = {}
        speciesColorDict = {}
        markerColorsSpecies = []
        hlinecSpecies = []
        vlinecSpecies = []
        hcolorsSpecies = []
        vcolorsSpecies = []
        print('SQLite Query: tree data')
        strainTimes = pd.read_sql_query(
            "SELECT tree_{0}, t_creation, t_extinction, tree_parent_{0} \
        FROM tree_{1}_creation_extinction WHERE tree_{0} in ({2})".format(
                strain_id, strain, ', '.join(map(str, keepTreeStrains))), conTree)
        treeAbundances = pd.read_sql_query(
            "SELECT t, tree_{0}, abundance \
        FROM tree_{1} WHERE tree_{0} in ({2})"
            .format(strain_id, abundance, ', '.join(map(str, keepTreeStrains))), conTree)
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
        if species == 'virus':
            norm = Normalize(vmin=float(1), vmax=float(1.35*max(newTreeOrder.values())))
        else:
            norm = Normalize(vmin=float(1), vmax=float(max(newTreeOrder.values())))
        maxAbundanceSpecies = np.max(treeAbundances.abundance.values)
        markerIncSpecies = maxticksize/maxAbundanceSpecies
        print('Compiling stacked species abundances and tree plots')
        for strainID in sorted(treeAbundances[tree_strain_id].unique()):
            speciesColorDict[strainID] = Vcmap(norm(strainID))
            if species == 'virus':
                tCreate = strainTimes[strainTimes[tree_strain_id]
                                    == strainID].t_creation.values[0]
                tExtinct = strainTimes[strainTimes[tree_strain_id]
                                    == strainID].t_extinction.values[0]
                parent = strainTimes[strainTimes[tree_strain_id]
                                    == strainID][tree_parent_strain_id].values[0]
                hlinecSpecies.append([[tCreate, strainID], [tExtinct, strainID]])
                if parent != 0:
                    vlinecSpecies.append([[tCreate, parent], [tCreate, strainID]])
                hcolorsSpecies.append(speciesColorDict[strainID])
                vcolorsSpecies.append(speciesColorDict[strainID])
                markerColorsSpecies.append(speciesColorDict[strainID])
                axes[idx1].scatter(treeAbundances[treeAbundances[tree_strain_id] == strainID]['t'],
                                treeAbundances[treeAbundances[tree_strain_id]
                                                == strainID][tree_strain_id],
                                lw=1,
                                s=treeAbundances[treeAbundances[tree_strain_id] ==
                                                    strainID]['abundance'].values*markerIncSpecies,
                                color=speciesColorDict[strainID], marker='|')
        if species == 'virus':
            strainLineages = LineCollection(
                hlinecSpecies, linestyles='solid', colors=hcolorsSpecies, linewidths=(1))
            creationLines = LineCollection(
                vlinecSpecies, linestyles='solid', colors=vcolorsSpecies, linewidths=(1))
            axes[idx1].add_collection(strainLineages)
            axes[idx1].add_collection(creationLines)
            axes[idx1].set_xlabel(xlabel='Time t', fontsize=7)
            axes[idx1].set_ylim(0, np.max(strainTimes[tree_strain_id])+1)
            axes[idx1].set_xlim(0, t)
            axes[idx1].xaxis.set_minor_locator(ticker.MultipleLocator(25))
            axes[idx1].set_yticks([])
            axes[idx1].set_xlabel(xlabel='Time t', fontsize=20)
            axes[idx1].ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
            axes[idx1].xaxis.set_minor_locator(ticker.MultipleLocator(25))
            axes[idx1].tick_params(axis='x', labelsize=20)
        species_stacked = species_stacked.merge(keepTreeStrainsDF, on=[strain_id])\
            .drop(columns=[strain_id]).replace({tree_strain_id: newTreeOrder})
        species_stacked.sort_values(by=['t', tree_strain_id])
        speciesDF = species_stacked.copy()
        species_stacked = species_stacked.pivot(
            index='t', columns=tree_strain_id, values='abundance')
        species_stacked.plot.area(ax=axes[idx0], stacked=True, legend=False,
                                  linewidth=0, color=speciesColorDict, sort_columns=True)
        axes[idx0].set_ylabel(ylabel=abundanceTitle, labelpad=20, fontsize=20)
        lim = axes[idx0].get_ylim()
        axes[idx0].set_ylim(0, lim[1])
        axes[idx0].ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
        axes[idx0].xaxis.set_minor_locator(ticker.MultipleLocator(25))
        axes[idx0].set_xlim(0, t)
        axes[idx0].tick_params(axis='x', labelsize=20)
        axes[idx0].tick_params(axis='y', labelsize=20)
        axes[idx0].yaxis.get_offset_text().set_fontsize(20)
    return fig, axes
