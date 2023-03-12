bAbunds, bKeepTreeStrainsDF, bSpeciesColorDict, _, _, _, _  = \
tree.speciesTreePlot(run_id,'microbe',DBSIM_PATH,DBTREE_PATH,\
treepalette,maxticksize,figxy,hratio,babundthreshold)
p = curSim.execute('SELECT viral_mutation_rate,viral_burst_size,microbe_carrying_capacity, \
                    spacer_acquisition_prob, adsorption_rate, viral_decay_rate \
                    FROM param_combos WHERE combo_id = {}'.format(combo_id)).fetchall()
mu = p[0][0]
beta = p[0][1]
K = p[0][2]
q = p[0][3]
phi = p[0][4]
d = p[0][5]
tripartiteNets = pd.read_sql_query("SELECT t, vmatch_id, spacer_id, bmatch_id \
    FROM single_match_tripartite_networks", conTri)
vmatchstrains = pd.read_sql_query("SELECT t, match_id, vstrain_id \
    FROM vmatches", conTri).rename(columns={'match_id':'vmatch_id'})
bstrain0vstrain = pd.read_sql_query("SELECT t, vstrain_id, bstrain_id \
    FROM bstrain_to_vstrain_0matches", conMatch)
bstrainmatches = pd.read_sql_query("SELECT t, bstrain_id, match_id \
    FROM bmatches", conTri).rename(columns={"match_id":"bmatch_id"})
microbe_stacked = pd.read_sql_query("SELECT t, bstrain_id, abundance \
    FROM babundance", conSim).rename(columns={"abundance":"babundance"})
virus_stacked = pd.read_sql_query("SELECT t, vstrain_id, abundance \
    FROM vabundance", conSim).rename(columns={"abundance":"vabundance"})
bstrain0vstrain = bstrain0vstrain.merge(microbe_stacked,on=['t','bstrain_id'])\
                    .merge(virus_stacked,on=['t','vstrain_id'])
bstrain0vstrain = bstrain0vstrain.drop(columns=['vstrain_id']).drop_duplicates()\
                .groupby(['t']).agg(bsus=('babundance','sum')).reset_index()
bmatchAbunds = pd.read_sql_query("SELECT t, match_id, babundance \
    FROM bmatches_abundances", conTri).rename(columns={"match_id": "bmatch_id"})\
    .merge(microbe_total,on=['t'])
bmatchAbunds['bfreq'] = bmatchAbunds['babundance']/bmatchAbunds['Microbial Abundance']
vmatchAbunds = pd.read_sql_query("SELECT t, match_id, vabundance, bsusceptible, bimmune \
    FROM vmatches_abundances", conTri).rename(columns={"match_id": "vmatch_id"})\
    .merge(virus_total,on=['t'])
vmatchAbunds['vfreq'] = vmatchAbunds['vabundance']/vmatchAbunds['Viral Abundance']
vmatchAbunds = vmatchAbunds.merge(microbe_total,on=['t'])
vmatchAbunds['bsfreq'] = vmatchAbunds['bsusceptible']/vmatchAbunds['Microbial Abundance']
vmatchAbunds = vmatchAbunds.drop(columns=['Microbial Abundance'])
tripartiteNets = tripartiteNets.merge(bmatchAbunds,on=['t','bmatch_id'])
tripartiteTotal = tripartiteNets.drop(columns=['spacer_id','vmatch_id']).drop_duplicates()\
                .groupby(['t']).agg(btri=('babundance','sum')).reset_index()
tripartiteNets = tripartiteNets.merge(vmatchAbunds,on=['t','vmatch_id'])\
                .merge(tripartiteTotal,on=['t'])
matchLengths = tripartiteNets.drop(columns=['bmatch_id']).drop_duplicates()\
                .groupby(['t','vmatch_id']).agg(match_length=('spacer_id','size')).reset_index()
tripartiteNets = tripartiteNets.merge(matchLengths,on=['t','vmatch_id'])
spacerEscProb = pd.DataFrame().assign(\
            t=tripartiteNets['t'], vmatch_id=tripartiteNets['vmatch_id'],\
            vfreq=tripartiteNets['vfreq'],\
            vabundance=tripartiteNets['vabundance'],\
            match_length=tripartiteNets['match_length'],\
            spacer_id=tripartiteNets['spacer_id'],\
            bmatch_id=tripartiteNets['bmatch_id'],\
            bsus=tripartiteNets['bsusceptible'],\
            bim=tripartiteNets['bimmune'],\
            bsfreq=tripartiteNets['bsfreq'],\
            bfreq=tripartiteNets['bfreq'],\
            babundance=tripartiteNets['babundance'],\
            btri=tripartiteNets['btri'],\
            btotal=tripartiteNets['Microbial Abundance'])
spacerEscDiv = pd.DataFrame().assign(\
            t=tripartiteNets['t'], vmatch_id=tripartiteNets['vmatch_id'],\
            spacer_id=tripartiteNets['spacer_id'],\
            bmatch_id=tripartiteNets['bmatch_id'])
spacerEscDiv = spacerEscDiv.merge(bstrainmatches,on=['t','bmatch_id'])
times = []
strains = []
for t in sorted(spacerEscDiv.t.unique()):
    non0tri = {*spacerEscDiv[spacerEscDiv.t==t]['bstrain_id'].values} - {*bstrain0vstrain[bstrain0vstrain.t==t]['bstrain_id'].values}
    if len(non0tri) > 0:
        times.extend(len(non0tri)*[t])
        strains.extend(non0tri)
    else:
        continue
#
spacerEscDivNon0Tri = pd.DataFrame().assign(\
            t=times,\
            bstrain_id=strains)
spacerEscFreq = spacerEscDivNon0Tri.merge(microbe_stacked,on=['t','bstrain_id'])
spacerEscFreq = spacerEscFreq.groupby(['t']).agg(btotal=('babundance','sum')).reset_index()\
                .merge(spacerEscFreq,on=['t'])
spacerEscFreq['freq'] =  np.array(spacerEscFreq['babundance'])/np.array(spacerEscFreq['btotal'])
spacerEscFreq['freq'] = np.exp(-1*spacerEscFreq['freq']*np.log(spacerEscFreq['freq']))
spacerEscFreq = spacerEscFreq.groupby(['t']).agg(div=('freq','prod')).reset_index()
spacerEscDivNon0Tri = spacerEscDivNon0Tri.merge(bstrainmatches,on=['t','bstrain_id'])\
                        .drop(columns=['bstrain_id']).drop_duplicates()

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
axes[3].plot(spacerEscFreq['t'],spacerEscFreq['div'],color='darkblue',linewidth=1)
axes[3].yaxis.tick_left()
axes[3].yaxis.set_label_position("left")
axes[3].set_ylabel(ylabel ='Shannon Diversity\nof Single-Match Tripartite Network',\
                    labelpad=15,fontsize=7)
fig.tight_layout()
