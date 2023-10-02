def coTreePlot(run_id,DBSIM_PATH,DBMATCH_PATH,DBTREE_PATH,treepalette,maxticksize,figxy,hratio,abundthreshold,abundthreshold2,stacked,overlay):
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
    ####
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
    keepTreeStrainsDF['new_tree_{}_id'.format(
        strain)] = keepTreeStrainsDF['tree_{}_id'.format(strain)].copy()
    keepTreeStrainsDF = keepTreeStrainsDF.replace(
        {'new_tree_{}_id'.format(strain): newTreeOrder})
    keepTreeStrainsDF['color'] = keepTreeStrainsDF['new_tree_{}_id'.format(strain)].copy()
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
                        lw=1,
                        s=treeAbundances[treeAbundances[tree_strain_id] ==
                                            strainID]['abundance'].values*markerIncSpecies,
                        color=speciesColorDict[strainID], marker='|')
    keepTreeStrainsDF = keepTreeStrainsDF[keepTreeStrainsDF.vstrain_id != 0].copy()
    colors = []
    for strainID in keepTreeStrainsDF['new_tree_{}_id'.format(strain)]:
        colors.append(speciesColorDict[strainID])
    keepTreeStrainsDF['color'] = colors
    keepTreeStrainsDF = keepTreeStrainsDF.replace(
        {'color': speciesColorDict})
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
        .drop(columns=[strain_id, tree_strain_id,'color']).rename(columns={'new_tree_{}_id'.format(strain):tree_strain_id})
    species_stacked.sort_values(by=['t', tree_strain_id])
    speciesDF = species_stacked.copy()
    species_stacked = species_stacked.pivot(
        index='t', columns=tree_strain_id, values='abundance')
    if (overlay):
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
            axes[2].set_ylabel(ylabel=abundanceTitle, labelpad=20, fontsize=10,rotation=270)
        else:
            axes.append(axes[0].twinx())
            axes[0].plot(microbe_total['t'],microbe_total['btotal'],linewidth=0,color='grey')
            axes[0].fill_between(microbe_total['t'],microbe_total['btotal'], color='grey',alpha=0.5)
            species_stacked.plot(ax=axes[2], stacked=False, legend=False,
                                    linewidth=1.5, color=speciesColorDict, sort_columns=True)
            axes[2].set_ylabel(ylabel='Viral Strain\nAbundances', labelpad=30, fontsize=10,rotation=270)
        axes[0].set_ylabel(ylabel='Host\nAbundance', labelpad=20, fontsize=10)
        lim = axes[2].get_ylim()
        axes[2].set_ylim(0, lim[1])
    if (not overlay) & (stacked):
        species_stacked.plot.area(ax=axes[0], stacked=True, legend=False,
                                linewidth=0, color=speciesColorDict, sort_columns=True)
        species_stacked.plot(
            stacked=True, ax=axes[0], legend=False, color='white', sort_columns=True, linewidth=.1)
        axes[0].set_ylabel(ylabel=abundanceTitle, labelpad=20, fontsize=10)
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


    ######
    keepTreeStrainsDFvirus = keepTreeStrainsDF.copy()
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
    #####
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
    maxIDs[straintotal] = abundthreshold2*np.array(maxIDs[straintotal])
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
    keepTreeStrainsDF['new_tree_{}_id'.format(
        strain)] = keepTreeStrainsDF['tree_{}_id'.format(strain)].copy()
    keepTreeStrainsDF = keepTreeStrainsDF.replace(
        {'new_tree_{}_id'.format(strain): newTreeOrder})
    numStrains = len(keepTreeStrains)
    #
    speciesColorDict = {}
    maxAbundanceSpecies = np.max(treeAbundances.abundance.values)
    markerIncSpecies = maxticksize/maxAbundanceSpecies
    markerColorsSpecies = []
    hlinecSpecies = []
    vlinecSpecies = []
    hcolorsSpecies = []
    vcolorsSpecies = []
    ####
    conMatch = sqlite3.connect(DBMATCH_PATH)
    curMatch = conMatch.cursor()
    infectionNet = pd.read_sql_query(
        "SELECT t, bstrain_id,vstrain_id FROM bstrain_to_vstrain_0matches", conMatch)\
            .merge(keepTreeStrainsDFvirus[['vstrain_id','color']],on=['vstrain_id'])\
            .drop_duplicates().merge(\
                pd.read_sql_query(
                    "SELECT t, vstrain_id, abundance FROM vabundance", conSim),\
        on=['t', 'vstrain_id']).rename(columns={'abundance':'vabundance'}).merge(
                pd.read_sql_query(
                    "SELECT t, bstrain_id, abundance FROM babundance", conSim),
        on=['t', 'bstrain_id']).rename(columns={'abundance': 'babundance'})
    infectionNet = infectionNet.groupby(['bstrain_id'])\
                    .agg(max=('babundance', 'max')).reset_index()\
                    .merge(infectionNet, on=['bstrain_id'])
    rows = list(np.array(infectionNet['babundance']) ==
                    np.array(infectionNet['max']))
    infectionNet = infectionNet[rows]
    infectionNet = infectionNet.drop(columns=['t','max'])
    infectionNet = infectionNet.groupby(['bstrain_id'])\
                    .agg(max=('vabundance', 'max')).reset_index()\
                    .merge(infectionNet, on=['bstrain_id'])
    rows = list(np.array(infectionNet['vabundance']) ==
                    np.array(infectionNet['max']))
    infectionNet = infectionNet[rows]
    infectionNet = infectionNet.drop(columns=['max'])
    infectionNet = infectionNet.merge(keepTreeStrainsDF,on=['bstrain_id'])
    noMatchColor = mc.to_rgba('lightgrey')
    ####
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
        if strainID not in np.unique(infectionNet['new_tree_bstrain_id']):
            speciesColorDict[strainID] = noMatchColor
        else:
            speciesColorDict[strainID] = infectionNet[infectionNet.new_tree_bstrain_id == strainID]['color'].values[0]
        hcolorsSpecies.append(speciesColorDict[strainID])
        vcolorsSpecies.append(speciesColorDict[strainID])
        markerColorsSpecies.append(speciesColorDict[strainID])
        axes[1].scatter(treeAbundances[treeAbundances[tree_strain_id] == strainID]['t'],
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
    axes[1].add_collection(strainLineages)
    axes[1].add_collection(creationLines)
    axes[1].set_xlabel(xlabel='Time t', fontsize=7)
    axes[1].set_ylim(0, np.max(strainTimes[tree_strain_id])+1)
    axes[1].set_xlim(0, t)
    axes[1].xaxis.set_minor_locator(ticker.MultipleLocator(25))
    axes[1].set_yticks([])
    species_stacked = species_stacked.merge(keepTreeStrainsDF, on=[strain_id])\
        .drop(columns=[strain_id, tree_strain_id]).rename(columns={'new_tree_{}_id'.format(strain):tree_strain_id})
    species_stacked.sort_values(by=['t', tree_strain_id])
    speciesDF = species_stacked.copy()
    species_stacked = species_stacked.pivot(
        index='t', columns=tree_strain_id, values='abundance')
    if (overlay):
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
    return speciesDF, keepTreeStrainsDF, speciesColorDict, hlinecSpecies, vlinecSpecies, fig, axes
