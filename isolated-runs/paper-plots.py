PLOT_PATH = os.path.join('/Users/armun','Desktop') # local
fig, ax = plt.subplots(2,sharex=True)
axes = [ax[0], ax[0].twinx(), ax[1], ax[1].twinx()]
microbe_stacked = bAbunds[bAbunds.t<=max(virus_total['t'])]\
            .pivot(index='t',columns='tree_bstrain_id',values='abundance')
microbe_stacked.plot.area(ax = axes[0],stacked=True,legend=False, linewidth=0,\
                            color=bSpeciesColorDict,sort_columns=True)
microbe_stacked.plot(stacked=True, ax=axes[0], legend=False, color='black',\
                        sort_columns=True,linewidth=.1)
axes[0].set_ylabel(ylabel ='Host Abundance',labelpad=15,fontsize=10)
axes[0].set_xlabel(xlabel = 'Time t',fontsize=10)
axes[0].ticklabel_format(axis = 'y',style='sci',scilimits=(0,0))
axes[0].xaxis.set_minor_locator(ticker.MultipleLocator(25))
axes[1].plot(virus_total['t'],virus_total['Viral Abundance'],linewidth=0,\
                color='grey')
axes[1].fill_between(virus_total['t'],virus_total['Viral Abundance'],\
                        color='grey',alpha=0.6)
axes[1].set_ylabel(ylabel ='Viral Abundance',labelpad=15,fontsize=10)
axes[0].margins(x=0)
axes[1].margins(x=0)
lim = axes[1].get_ylim()
axes[1].set_ylim(0,lim[1])
axes[2].fill_between(virus_total['t'],virus_total['Viral Abundance'],\
                        color='grey',alpha=0.6)
axes[2].set_ylabel(ylabel ='Viral Abundance',labelpad=15,fontsize=10)
axes[2].margins(x=0)
axes[2].margins(x=0)
axes[2].xaxis.set_minor_locator(ticker.MultipleLocator(25))
axes[3].plot(sorted(spacerExpMax2['t'].unique()),np.array(spacerExpMax2['exp']),\
            color='darkblue',linewidth=1.5,label='2nd largest (scheme 1)')
axes[3].plot(sorted(spacerExpMax['t'].unique()),np.array(spacerExpMax['exp']),\
            color='darkgreen',linewidth=1.5,label='1st largest (scheme 1)')
axes[3].plot(sorted(spacerExp1['t'].unique()),np.array(spacerExp1['exp']),\
            color='darkorange',linewidth=1.5,label='fuel (scheme 1)')
axes[3].legend()
# axes[3].plot(sorted(spacerEscProb['t'].unique()),\
#         np.ones(len(spacerEscProb['t'].unique())),color='red',linewidth=1.0)
axes[3].set_ylabel(ylabel ='Expected Maximum of\nSingly-matched Host Abundances',\
            labelpad=15,fontsize=10)
axes[3].set_xlabel(xlabel = 'Time t',fontsize=10)
lim = axes[2].get_ylim()
axes[2].set_ylim(2,lim[1])
lim = axes[3].get_ylim()
axes[3].set_ylim(0,lim[1])
axes[3].yaxis.tick_left()
axes[3].yaxis.set_label_position("left")
axes[2].yaxis.tick_right()
axes[2].yaxis.set_label_position("right")
# lim = axes[2].get_ylim()
# axes[2].set_ylim(0,lim[1])
fig.tight_layout()
fig.savefig(os.path.join(PLOT_PATH,''),dpi=resolve)
###
###
###
fig, ax = plt.subplots(2,sharex=True)
axes = [ax[0], ax[0].twinx(), ax[1], ax[1].twinx()]
microbe_stacked = bAbunds[bAbunds.t<=max(virus_total['t'])]\
            .pivot(index='t',columns='tree_bstrain_id',values='abundance')
microbe_stacked.plot.area(ax = axes[0],stacked=True,legend=False, linewidth=0,\
                            color=bSpeciesColorDict,sort_columns=True)
microbe_stacked.plot(stacked=True, ax=axes[0], legend=False, color='black',\
                        sort_columns=True,linewidth=.1)
axes[0].set_ylabel(ylabel ='Host Abundance',labelpad=15,fontsize=10)
axes[0].set_xlabel(xlabel = 'Time t',fontsize=10)
axes[0].ticklabel_format(axis = 'y',style='sci',scilimits=(0,0))
axes[0].xaxis.set_minor_locator(ticker.MultipleLocator(25))
axes[1].plot(virus_total['t'],virus_total['Viral Abundance'],linewidth=0,\
                color='grey')
axes[1].fill_between(virus_total['t'],virus_total['Viral Abundance'],\
                        color='grey',alpha=0.6)
axes[1].set_ylabel(ylabel ='Viral Abundance',labelpad=15,fontsize=10)
axes[0].margins(x=0)
axes[1].margins(x=0)
lim = axes[1].get_ylim()
axes[1].set_ylim(0,lim[1])
axes[2].fill_between(virus_total['t'],virus_total['Viral Abundance'],\
                        color='grey',alpha=0.6)
axes[2].set_ylabel(ylabel ='Viral Abundance',labelpad=15,fontsize=10)
axes[2].margins(x=0)
axes[2].margins(x=0)
axes[2].xaxis.set_minor_locator(ticker.MultipleLocator(25))
axes[3].plot(sorted(spacerEscProb['t'].unique()),np.array(spacerExpMax['exp'])/K,\
            color='darkorange',linewidth=1.5,label='Max (scheme 1)')
axes[3].plot(sorted(spacerEscProb['t'].unique()),np.array(spacerExpStd['exp'])/K,\
            color='darkblue',linewidth=1.5,label='S.Dev. (scheme 1)')
axes[3].legend()
# axes[3].plot(sorted(spacerEscProb['t'].unique()),\
#         np.ones(len(spacerEscProb['t'].unique())),color='red',linewidth=1.0)
axes[3].set_ylabel(ylabel ='Expected Singly-matched Host Frequency (1/K)',\
            labelpad=15,fontsize=10)
axes[3].set_xlabel(xlabel = 'Time t',fontsize=10)
lim = axes[2].get_ylim()
axes[2].set_ylim(2,lim[1])
lim = axes[3].get_ylim()
axes[3].set_ylim(0,lim[1])
axes[3].yaxis.tick_left()
axes[3].yaxis.set_label_position("left")
axes[2].yaxis.tick_right()
axes[2].yaxis.set_label_position("right")
# lim = axes[2].get_ylim()
# axes[2].set_ylim(0,lim[1])
fig.tight_layout()
fig.savefig(os.path.join(PLOT_PATH,''),dpi=resolve)
###
###
###
fig, ax = plt.subplots(2,sharex=True)
axes = [ax[0], ax[0].twinx(), ax[1], ax[1].twinx(), ax[1].twinx()]
microbe_stacked = bAbunds[bAbunds.t<=max(virus_total['t'])]\
            .pivot(index='t',columns='tree_bstrain_id',values='abundance')
microbe_stacked.plot.area(ax = axes[0],stacked=True,legend=False, linewidth=0,\
                            color=bSpeciesColorDict,sort_columns=True)
microbe_stacked.plot(stacked=True, ax=axes[0], legend=False, color='black',\
                        sort_columns=True,linewidth=.1)
axes[0].set_ylabel(ylabel ='Host Abundance',labelpad=15,fontsize=10)
axes[0].set_xlabel(xlabel = 'Time t',fontsize=10)
axes[0].ticklabel_format(axis = 'y',style='sci',scilimits=(0,0))
axes[0].xaxis.set_minor_locator(ticker.MultipleLocator(25))
axes[1].plot(virus_total['t'],virus_total['Viral Abundance'],linewidth=0,\
                color='grey')
axes[1].fill_between(virus_total['t'],virus_total['Viral Abundance'],\
                        color='grey',alpha=0.6)
axes[1].set_ylabel(ylabel ='Viral Abundance',labelpad=15,fontsize=10)
axes[0].margins(x=0)
axes[1].margins(x=0)
lim = axes[1].get_ylim()
axes[1].set_ylim(0,lim[1])
axes[2].fill_between(virus_total['t'],virus_total['Viral Abundance'],\
                        color='grey',alpha=0.6)
axes[2].margins(x=0)
axes[2].margins(x=0)
axes[2].set_yticks([])
axes[2].xaxis.set_minor_locator(ticker.MultipleLocator(25))
l1, = axes[3].plot(sorted(spacerEscProb['t'].unique()),np.array(spacerExpMax['exp'])/K,\
            color='darkorange',linewidth=1.5,label='Max of Singly-matched Hosts (scheme 1)')
l2, = axes[3].plot(sorted(spacerEscProb['t'].unique()),np.array(spacerExpStd['exp'])/K,\
            color='darkgreen',linewidth=1.5,label='S.Dev of Singly-matched Hosts (scheme 1)')
l3, = axes[4].plot(sorted(spacerEscProb['t'].unique()),np.array(spacerExp7['exp']),\
            color='darkblue',linewidth=1.5,label='No. of Single Matches (scheme 1)')
axes[3].legend(handles=[l1,l2,l3])
# axes[3].plot(sorted(spacerEscProb['t'].unique()),\
#         np.ones(len(spacerEscProb['t'].unique())),color='red',linewidth=1.0)
axes[3].set_ylabel(ylabel ='Expected Singly-matched Host Frequency (1/K)',\
            labelpad=15,fontsize=10)
axes[4].set_ylabel(ylabel ='Expected Number of Single Matches\nper Viral Strain',\
            labelpad=15,fontsize=10)
axes[3].set_xlabel(xlabel = 'Time t',fontsize=10)
lim = axes[2].get_ylim()
axes[2].set_ylim(2,lim[1])
lim = axes[3].get_ylim()
axes[3].set_ylim(0,lim[1])
axes[3].yaxis.tick_left()
axes[3].yaxis.set_label_position("left")
axes[4].yaxis.tick_right()
axes[4].yaxis.set_label_position("right")
# lim = axes[2].get_ylim()
# axes[2].set_ylim(0,lim[1])
fig.tight_layout()
fig.savefig(os.path.join(PLOT_PATH,''),dpi=resolve)
###
###
fig, ax = plt.subplots(2,sharex=True)
axes = [ax[0], ax[0].twinx(), ax[1], ax[1].twinx()]
microbe_stacked = bAbunds[bAbunds.t<=max(virus_total['t'])]\
            .pivot(index='t',columns='tree_bstrain_id',values='abundance')
microbe_stacked.plot.area(ax = axes[0],stacked=True,legend=False, linewidth=0,\
                            color=bSpeciesColorDict,sort_columns=True)
microbe_stacked.plot(stacked=True, ax=axes[0], legend=False, color='black',\
                        sort_columns=True,linewidth=.1)
axes[0].set_ylabel(ylabel ='Host Abundance',labelpad=15,fontsize=10)
axes[0].set_xlabel(xlabel = 'Time t',fontsize=10)
axes[0].ticklabel_format(axis = 'y',style='sci',scilimits=(0,0))
axes[0].xaxis.set_minor_locator(ticker.MultipleLocator(25))
axes[1].plot(virus_total['t'],virus_total['Viral Abundance'],linewidth=0,\
                color='grey')
axes[1].fill_between(virus_total['t'],virus_total['Viral Abundance'],\
                        color='grey',alpha=0.6)
axes[1].set_ylabel(ylabel ='Viral Abundance',labelpad=15,fontsize=10)
axes[0].margins(x=0)
axes[1].margins(x=0)
lim = axes[1].get_ylim()
axes[1].set_ylim(0,lim[1])
axes[2].fill_between(virus_total['t'],virus_total['Viral Abundance'],\
                        color='grey',alpha=0.6)
axes[2].set_ylabel(ylabel ='Viral Abundance',labelpad=15,fontsize=10)
axes[2].margins(x=0)
axes[2].margins(x=0)
axes[2].xaxis.set_minor_locator(ticker.MultipleLocator(25))
axes[3].plot(sorted(spacerEscProb['t'].unique()),np.array(spacerExp1['exp'])/K,\
            color='darkorange',linewidth=1.5,label='Susceptibles (infection)')
axes[3].plot(sorted(spacerEscProb['t'].unique()),np.array(spacerExpMax['exp'])/K,\
            color='darkblue',linewidth=1.5,label='Max of Singly-matched Hosts (scheme 1)')
axes[3].legend()
# axes[3].plot(sorted(spacerEscProb['t'].unique()),\
#         np.ones(len(spacerEscProb['t'].unique())),color='red',linewidth=1.0)
axes[3].set_ylabel(ylabel ='Expected Host Frequency (1/K)',\
            labelpad=15,fontsize=10)
axes[3].set_xlabel(xlabel = 'Time t',fontsize=10)
lim = axes[2].get_ylim()
axes[2].set_ylim(2,lim[1])
lim = axes[3].get_ylim()
axes[3].set_ylim(0,lim[1])
axes[3].yaxis.tick_left()
axes[3].yaxis.set_label_position("left")
axes[2].yaxis.tick_right()
axes[2].yaxis.set_label_position("right")
# lim = axes[2].get_ylim()
# axes[2].set_ylim(0,lim[1])
fig.tight_layout()
fig.savefig(os.path.join(PLOT_PATH,''),dpi=resolve)
##
fig, ax = plt.subplots(2,sharex=True)
axes = [ax[0], ax[0].twinx(), ax[1], ax[1].twinx()]
microbe_stacked = bAbunds[bAbunds.t<=max(virus_total['t'])]\
            .pivot(index='t',columns='tree_bstrain_id',values='abundance')
microbe_stacked.plot.area(ax = axes[0],stacked=True,legend=False, linewidth=0,\
                            color=bSpeciesColorDict,sort_columns=True)
microbe_stacked.plot(stacked=True, ax=axes[0], legend=False, color='black',\
                        sort_columns=True,linewidth=.1)
axes[0].set_ylabel(ylabel ='Host Abundance',labelpad=15,fontsize=10)
axes[0].set_xlabel(xlabel = 'Time t',fontsize=10)
axes[0].ticklabel_format(axis = 'y',style='sci',scilimits=(0,0))
axes[0].xaxis.set_minor_locator(ticker.MultipleLocator(25))
axes[1].plot(virus_total['t'],virus_total['Viral Abundance'],linewidth=0,\
                color='grey')
axes[1].fill_between(virus_total['t'],virus_total['Viral Abundance'],\
                        color='grey',alpha=0.6)
axes[1].set_ylabel(ylabel ='Viral Abundance',labelpad=15,fontsize=10)
axes[0].margins(x=0)
axes[1].margins(x=0)
lim = axes[1].get_ylim()
axes[1].set_ylim(0,lim[1])
axes[2].fill_between(virus_total['t'],virus_total['Viral Abundance'],\
                        color='grey',alpha=0.6)
axes[2].set_ylabel(ylabel ='Viral Abundance',labelpad=15,fontsize=10)
axes[2].margins(x=0)
axes[2].margins(x=0)
axes[2].xaxis.set_minor_locator(ticker.MultipleLocator(25))
axes[3].plot(sorted(spacerEscProb['t'].unique()),np.array(spacerExp3['exp']),\
            color='darkorange',linewidth=1.5,label='R0 (infection)')
axes[3].plot(sorted(spacerEscProb['t'].unique()),np.array(spacerExpR1Max['exp']),\
            color='darkblue',linewidth=1.5,label='Max R1 (scheme 1)')
axes[3].legend()
axes[3].plot(sorted(spacerEscProb['t'].unique()),\
        np.ones(len(spacerEscProb['t'].unique())),color='red',linewidth=1.0)
axes[3].set_ylabel(ylabel ='Expected R',\
            labelpad=15,fontsize=10)
axes[3].set_xlabel(xlabel = 'Time t',fontsize=10)
lim = axes[2].get_ylim()
axes[2].set_ylim(2,lim[1])
lim = axes[3].get_ylim()
axes[3].set_ylim(0,lim[1])
axes[3].yaxis.tick_left()
axes[3].yaxis.set_label_position("left")
axes[2].yaxis.tick_right()
axes[2].yaxis.set_label_position("right")
# lim = axes[2].get_ylim()
# axes[2].set_ylim(0,lim[1])
fig.tight_layout()
fig.savefig(os.path.join(PLOT_PATH,''),dpi=resolve)











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
tripartiteNets = pd.read_sql_query("SELECT t, bstrain_id, vstrain_id, time_specific_match_id\
    FROM bstrain_to_vstrain_matches WHERE match_length = 1", conMatch)
spacerMatches = pd.read_sql_query("SELECT t, time_specific_match_id, spacer_id\
    FROM matches_spacers", conMatch)
tripartiteNets = tripartiteNets.merge(spacerMatches,on=['t','time_specific_match_id'])\
            .drop(columns=['time_specific_match_id'])
microbe_stacked = pd.read_sql_query("SELECT t,bstrain_id,abundance FROM babundance \
                    WHERE run_id = {}".format(run_id), conSim)
virus_stacked = pd.read_sql_query("SELECT t,vstrain_id,abundance FROM vabundance \
                    WHERE run_id = {}".format(run_id), conSim)
spacerEscProb = tripartiteNets.merge(microbe_stacked,on=['t','bstrain_id'])\
                    .rename(columns={'abundance':'babundance'})\
                    .merge(virus_stacked,on=['t','vstrain_id'])\
                    .rename(columns={'abundance':'vabundance'})\
                    .merge(virus_total,on=['t'])\
                    .merge(microbe_total,on=['t'])\
                    .rename(columns={'Viral Abundance':'vtotal',\
                                'Microbial Abundance':'btotal'})
bsusceptible = pd.read_sql_query("SELECT t, bstrain_id, vstrain_id\
                FROM bstrain_to_vstrain_0matches", conMatch)\
                .merge(microbe_stacked,on=['t','bstrain_id'])\
                .rename(columns={'abundance':'bsus'})\
                .groupby(['t','vstrain_id']).agg(bsus=('bsus','sum'))\
                .reset_index()
bimmune = pd.read_sql_query("SELECT t, bstrain_id, vstrain_id\
                FROM bstrain_to_vstrain_matches", conMatch)\
                .merge(microbe_stacked,on=['t','bstrain_id'])\
                .rename(columns={'abundance':'bim'})\
                .groupby(['t','vstrain_id']).agg(bim=('bim','sum'))\
                .reset_index()
spacerEscProb = spacerEscProb.merge(bsusceptible,on=['t','vstrain_id'])\
                    .merge(bimmune,on=['t','vstrain_id'])
bspacerAbunds = spacerEscProb[['t','vstrain_id','spacer_id','babundance']]\
                    .groupby(['t','vstrain_id','spacer_id'])\
                    .agg(babundance=('babundance','sum')).reset_index()
spacerEscProb = spacerEscProb.drop(columns=['bstrain_id','babundance'])\
                    .merge(bspacerAbunds,on=['t','vstrain_id','spacer_id'])\
                    .drop_duplicates()
matchLengths = spacerEscProb[['t','vstrain_id','spacer_id']]\
                .drop_duplicates().groupby(['t','vstrain_id'])\
                .agg(match_length=('spacer_id','size')).reset_index()
spacerEscProb = spacerEscProb.merge(matchLengths,on=['t','vstrain_id'])
max1 = spacerEscProb[['t','vstrain_id','spacer_id','babundance']]\
                        .groupby(['t','vstrain_id'])\
                        .agg(max=('babundance','max')).reset_index()\
                        .merge(spacerEscProb[['t','vstrain_id','spacer_id','babundance']],\
                        on=['t','vstrain_id'])
max1rows = list(np.array(max1['babundance']) == np.array(max1['max'])) # keep 1st rank
max1 = max1[max1rows] # keep 1st rank
#
max2 = spacerEscProb[['t','vstrain_id','spacer_id','babundance']]\
                        .merge(max1[['t','vstrain_id','max']].drop_duplicates(),\
                        on=['t','vstrain_id'])
max1rows = list(np.array(max2['babundance']) != np.array(max2['max'])) # remove 1st rank
max2 = max2[max1rows].drop(columns=['max']) # remove 1st rank
max2 = max2.groupby(['t','vstrain_id'])\
            .agg(max=('babundance','max')).reset_index()\
            .merge(max2,on=['t','vstrain_id'])
max2rows = list(np.array(max2['babundance']) == np.array(max2['max'])) # keep 2nd rank
max2 = max2[max2rows] # keep 2nd rank
#
max3 = spacerEscProb[['t','vstrain_id','spacer_id','babundance']]\
                        .merge(max1[['t','vstrain_id','max']],\
                        on=['t','vstrain_id'])
max1rows = list(np.array(max3['babundance']) != np.array(max3['max'])) # remove 1st rank
max3 = max3[max1rows].drop(columns=['max']) # remove 1st rank
max3 = max3.merge(max2[['t','vstrain_id','max']],on=['t','vstrain_id'])
max2rows = list(np.array(max3['babundance']) != np.array(max3['max'])) # remove 2nd rank
max3 = max3[max2rows].drop(columns=['max']) # remove 2nd rank
max3 = max3.groupby(['t','vstrain_id'])\
            .agg(max=('babundance','max')).reset_index()\
            .merge(max3,on=['t','vstrain_id'])
max3rows = list(np.array(max3['babundance']) == np.array(max3['max'])) # keep 3rd rank
max3 = max3[max3rows] # keep 3rd rank
# at least...
maxAtLeast2 = pd.concat([max1, max2]).sort_values(by=['t','vstrain_id'])
maxAtLeast3 = pd.concat([maxAtLeast2, max3]).sort_values(by=['t','vstrain_id'])
#
std = spacerEscProb[['t','vstrain_id','spacer_id','babundance']]\
                        .groupby(['t','vstrain_id'])\
                        .agg(std=('babundance','std')).reset_index()
## Probabilities and Normalizations
spacerEscProb['infectionProbs'] = \
    (np.array(spacerEscProb['vabundance'])/np.array(spacerEscProb['vtotal']))\
    *(np.array(spacerEscProb['bsus'])/np.array(spacerEscProb['btotal']))
spacerEscProb['infectionProbs2'] = \
    np.array(spacerEscProb['infectionProbs'])\
            *beta*mu*np.array(spacerEscProb['match_length'])
norm = spacerEscProb.groupby(['t']).agg(norm=('infectionProbs','sum'))\
            .reset_index()
spacerEscProb = spacerEscProb.merge(norm,on=['t'])
norm = spacerEscProb[['t','vstrain_id','infectionProbs']].drop_duplicates()\
        .groupby(['t']).agg(norm1=('infectionProbs','sum')).reset_index()
spacerEscProb = spacerEscProb.merge(norm,on=['t'])
norm = spacerEscProb[['t','vstrain_id','infectionProbs2']].drop_duplicates()\
        .groupby(['t']).agg(norm2=('infectionProbs2','sum')).reset_index()
spacerEscProb = spacerEscProb.merge(norm,on=['t'])
spacerEscProb['infectionProbsNormed'] = \
    np.array(spacerEscProb['infectionProbs'])\
    /np.array(spacerEscProb['norm'])
spacerEscProb['infectionProbsNormed2'] = \
    ((np.array(spacerEscProb['infectionProbs2'])\
    /np.array(spacerEscProb['norm2'])))/np.array(spacerEscProb['match_length'])
spacerEscProb['infectionProbsNormed1'] = \
    ((np.array(spacerEscProb['infectionProbs'])\
    /np.array(spacerEscProb['norm1'])))
spacerEscProb = spacerEscProb.drop(columns=['norm','norm1','norm2'])
#
max1 = spacerEscProb[['t','vstrain_id',\
            'infectionProbsNormed','infectionProbsNormed1']]\
                .drop_duplicates().merge(max1,on=['t','vstrain_id'])
max2 = spacerEscProb[['t','vstrain_id',\
            'infectionProbsNormed','infectionProbsNormed1']]\
                .drop_duplicates().merge(max2,on=['t','vstrain_id'])
max3 = spacerEscProb[['t','vstrain_id',\
            'infectionProbsNormed','infectionProbsNormed1']]\
                .drop_duplicates().merge(max3,on=['t','vstrain_id'])
maxAtLeast2 = spacerEscProb[['t','vstrain_id',\
            'infectionProbsNormed','infectionProbsNormed1']]\
                .drop_duplicates().merge(maxAtLeast2,on=['t','vstrain_id'])
maxAtLeast3 = spacerEscProb[['t','vstrain_id',\
            'infectionProbsNormed','infectionProbsNormed1']]\
                .drop_duplicates().merge(maxAtLeast3,on=['t','vstrain_id'])
###
norm =max1.groupby(['t']).agg(norm=('infectionProbsNormed1','sum')).reset_index()
max1 = max1.merge(norm,on=['t'])
max1['infectionProbsNormed1'] = ((np.array(max1['infectionProbsNormed1'])\
                                        /np.array(max1['norm'])))
max1['expMax'] = np.array(max1['max'])*\
                            np.array(max1['infectionProbsNormed1'])
norm =max2.groupby(['t']).agg(norm=('infectionProbsNormed','sum')).reset_index()
max2 = max2.merge(norm,on=['t'])
max2['infectionProbsNormed'] = ((np.array(max2['infectionProbsNormed'])\
                                        /np.array(max2['norm'])))
max2['expMax'] = np.array(max2['max'])*\
                            np.array(max2['infectionProbsNormed'])
norm =max3.groupby(['t']).agg(norm=('infectionProbsNormed','sum')).reset_index()
max3 = max3.merge(norm,on=['t'])
max3['infectionProbsNormed'] = ((np.array(max3['infectionProbsNormed'])\
                                        /np.array(max3['norm'])))
max3['expMax'] = np.array(max3['max'])*\
                            np.array(max3['infectionProbsNormed'])
norm = maxAtLeast2.groupby(['t']).agg(norm=('infectionProbsNormed','sum')).reset_index()
maxAtLeast2 = maxAtLeast2.merge(norm,on=['t'])
maxAtLeast2['infectionProbsNormed'] = \
            ((np.array(maxAtLeast2['infectionProbsNormed'])\
            /np.array(maxAtLeast2['norm'])))
maxAtLeast2['expMax'] = np.array(maxAtLeast2['max'])*\
                            np.array(maxAtLeast2['infectionProbsNormed'])
norm = maxAtLeast3.groupby(['t']).agg(norm=('infectionProbsNormed','sum')).reset_index()
maxAtLeast3 = maxAtLeast3.merge(norm,on=['t'])
maxAtLeast3['infectionProbsNormed'] = \
            ((np.array(maxAtLeast3['infectionProbsNormed'])\
            /np.array(maxAtLeast3['norm'])))
maxAtLeast3['expMax'] = np.array(maxAtLeast3['max'])*\
                            np.array(maxAtLeast3['infectionProbsNormed'])
#
expMax1 = max1[['t','vstrain_id','expMax']].drop_duplicates()\
                .groupby(['t']).agg(exp=('expMax','sum')).reset_index()
expMax2 = max2[['t','vstrain_id','expMax']].drop_duplicates()\
                .groupby(['t']).agg(exp=('expMax','sum')).reset_index()
expMax3 = max3[['t','vstrain_id','expMax']].drop_duplicates()\
                .groupby(['t']).agg(exp=('expMax','sum')).reset_index()
expMaxAtLeast2 = maxAtLeast2[['t','vstrain_id','expMax']].drop_duplicates()\
                .groupby(['t']).agg(exp=('expMax','sum')).reset_index()
expMaxAtLeast3 = maxAtLeast3[['t','vstrain_id','expMax']].drop_duplicates()\
                .groupby(['t']).agg(exp=('expMax','sum')).reset_index()


fig, ax = plt.subplots(2,sharex=True)
axes = [ax[0], ax[0].twinx(), ax[1], ax[1].twinx()]
microbe_stacked = bAbunds[bAbunds.t<=max(virus_total['t'])]\
            .pivot(index='t',columns='tree_bstrain_id',values='abundance')
microbe_stacked.plot.area(ax = axes[0],stacked=True,legend=False, linewidth=0,\
                            color=bSpeciesColorDict,sort_columns=True)
microbe_stacked.plot(stacked=True, ax=axes[0], legend=False, color='black',\
                        sort_columns=True,linewidth=.1)
axes[0].set_ylabel(ylabel ='Host Abundance',labelpad=15,fontsize=10)
axes[0].set_xlabel(xlabel = 'Time t',fontsize=10)
axes[0].ticklabel_format(axis = 'y',style='sci',scilimits=(0,0))
axes[0].xaxis.set_minor_locator(ticker.MultipleLocator(25))
axes[1].plot(virus_total['t'],virus_total['Viral Abundance'],linewidth=0,\
                color='grey')
axes[1].fill_between(virus_total['t'],virus_total['Viral Abundance'],\
                        color='grey',alpha=0.6)
axes[1].set_ylabel(ylabel ='Viral Abundance',labelpad=15,fontsize=10)
axes[0].margins(x=0)
axes[1].margins(x=0)
lim = axes[1].get_ylim()
axes[1].set_ylim(0,lim[1])
axes[2].fill_between(virus_total['t'],virus_total['Viral Abundance'],\
                        color='grey',alpha=0.6)
axes[2].set_ylabel(ylabel ='Viral Abundance',labelpad=15,fontsize=10)
axes[2].margins(x=0)
axes[2].margins(x=0)
axes[2].xaxis.set_minor_locator(ticker.MultipleLocator(25))
axes[3].plot(sorted(spacerExp1['t'].unique()),np.array(spacerExp1['exp'])/K,\
            color='darkorange',linewidth=1.5,label='Susceptibles')
axes[3].plot(sorted(expMax1['t'].unique()),np.array(expMax1['exp'])/K,\
            color='darkred',linewidth=1.5,label='Rank 1 Single Matches')
axes[3].plot(sorted(expMaxAtLeast2['t'].unique()),np.array(expMaxAtLeast2['exp'])/K,\
            color='darkgreen',linewidth=1.5,label='At Least Rank 2 Single Matches')
axes[3].plot(sorted(expMaxAtLeast3['t'].unique()),np.array(expMaxAtLeast3['exp'])/K,\
            color='darkblue',linewidth=1.5,label='At Least Rank 3 Single Matches')
axes[3].legend()
# axes[3].plot(sorted(spacerEscProb['t'].unique()),\
#         np.ones(len(spacerEscProb['t'].unique())),color='red',linewidth=1.0)
axes[3].set_ylabel(ylabel ='Expected Host Frequencies (1/K)',\
            labelpad=15,fontsize=10)
axes[3].set_xlabel(xlabel = 'Time t',fontsize=10)
lim = axes[2].get_ylim()
axes[2].set_ylim(2,lim[1])
lim = axes[3].get_ylim()
axes[3].set_ylim(0,lim[1])
axes[3].yaxis.tick_left()
axes[3].yaxis.set_label_position("left")
axes[2].yaxis.tick_right()
axes[2].yaxis.set_label_position("right")
# lim = axes[2].get_ylim()
# axes[2].set_ylim(0,lim[1])
fig.tight_layout()



microbe_stacked = pd.read_sql_query("SELECT t,bstrain_id,abundance FROM babundance \
                    WHERE run_id = {}".format(run_id), conSim)
microbe_stacked = microbe_stacked[microbe_stacked.abundance>0]
bsusceptible = pd.read_sql_query("SELECT t, bstrain_id, vstrain_id\
                FROM bstrain_to_vstrain_0matches", conMatch)\
                .drop(columns=['vstrain_id']).drop_duplicates()\
                .merge(microbe_stacked,on=['t','bstrain_id'])
microbe_stacked = microbe_stacked.groupby(['t']).agg(btotal=('abundance', 'sum'))\
    .reset_index().merge(microbe_stacked, on=['t'])
microbe_stacked['bfreq'] = np.array(microbe_stacked['abundance']) / \
    np.array(microbe_stacked['btotal'])
microbe_stacked['shannon'] = microbe_stacked['bfreq']*np.log(microbe_stacked['bfreq'])
bshannon = microbe_stacked.groupby(['t']).agg(shannon=('shannon', 'sum'))\
    .reset_index()
bshannon['shannon'] = np.exp(-1*np.array(bshannon['shannon']))
bshannon = bshannon[bshannon.t <= max(bshannonS.t)]
bsusceptible = bsusceptible.merge(microbe_stacked[['t','btotal']]\
                                    .drop_duplicates(), on=['t'])
bsusceptible['bfreq'] = np.array(bsusceptible['abundance']) / \
    np.array(bsusceptible['btotal'])
bsusceptible['shannon'] = bsusceptible['bfreq']*np.log(bsusceptible['bfreq'])
bshannonS = bsusceptible.groupby(['t']).agg(shannon=('shannon', 'sum'))\
    .reset_index()
bshannonS['shannon'] = np.exp(-1*np.array(bshannonS['shannon']))
bstotal = bsusceptible.groupby(['t']).agg(bfreq=('bfreq', 'sum'))\
            .reset_index()

fig, ax = plt.subplots(2,sharex=True)
axes = [ax[0], ax[0].twinx(), ax[1], ax[1].twinx()]
microbe_stacked = bAbunds[bAbunds.t<=max(virus_total['t'])]\
            .pivot(index='t',columns='tree_bstrain_id',values='abundance')
microbe_stacked.plot.area(ax = axes[0],stacked=True,legend=False, linewidth=0,\
                            color=bSpeciesColorDict,sort_columns=True)
microbe_stacked.plot(stacked=True, ax=axes[0], legend=False, color='black',\
                        sort_columns=True,linewidth=.1)
axes[0].set_ylabel(ylabel ='Host Abundance',labelpad=15,fontsize=10)
axes[0].set_xlabel(xlabel = 'Time t',fontsize=10)
axes[0].ticklabel_format(axis = 'y',style='sci',scilimits=(0,0))
axes[0].xaxis.set_minor_locator(ticker.MultipleLocator(25))
axes[1].plot(virus_total['t'],virus_total['Viral Abundance'],linewidth=0,\
                color='grey')
axes[1].fill_between(virus_total['t'],virus_total['Viral Abundance'],\
                        color='grey',alpha=0.6)
axes[1].set_ylabel(ylabel ='Viral Abundance',labelpad=15,fontsize=10)
axes[0].margins(x=0)
axes[1].margins(x=0)
lim = axes[1].get_ylim()
axes[1].set_ylim(0,lim[1])
axes[2].fill_between(virus_total['t'],virus_total['Viral Abundance'],\
                        color='grey',alpha=0.6)
axes[2].set_ylabel(ylabel ='Viral Abundance',labelpad=15,fontsize=10)
axes[2].margins(x=0)
axes[2].margins(x=0)
axes[2].xaxis.set_minor_locator(ticker.MultipleLocator(25))
# axes[3].plot(sorted(bshannon['t'].unique()),\
#             np.array(bshannonS['shannon'])/np.array(bshannon['shannon']),\
#             color='darkblue',linewidth=1.5,label='Susc. Shannon Div.')
axes[3].plot(sorted(bstotal['t'].unique()),
             np.array(bstotal['bfreq']),
             color='darkblue', linewidth=1.5)
# axes[3].legend()
axes[3].set_ylabel(ylabel ='Susceptible Frequency',\
            labelpad=15,fontsize=10)
axes[3].set_xlabel(xlabel = 'Time t',fontsize=10)
lim = axes[2].get_ylim()
axes[2].set_ylim(2,lim[1])
lim = axes[3].get_ylim()
axes[3].set_ylim(0,lim[1])
axes[3].yaxis.tick_left()
axes[3].yaxis.set_label_position("left")
axes[2].yaxis.tick_right()
axes[2].yaxis.set_label_position("right")
# lim = axes[2].get_ylim()
# axes[2].set_ylim(0,lim[1])
fig.tight_layout()


bAbunds, bKeepTreeStrainsDF, bSpeciesColorDict, _, _, _, _ = \
tree.speciesTreePlot(run_id, 'microbe',DBSIM_PATH,DBTREE_PATH,\
treepalette, maxticksize,figxy,hratio,babundthreshold)
bstrain0 = pd.read_sql_query("SELECT t, bstrain_id\
    FROM bstrain_to_vstrain_0matches", conMatch)\
        .drop_duplicates()
microbe_stacked = pd.read_sql_query("SELECT t,bstrain_id,abundance FROM babundance \
                    WHERE run_id = {}".format(run_id), conSim)
microbe_stacked = microbe_stacked[(\
    microbe_stacked.t >= 275) & (microbe_stacked.t <= 800)]
triRichness = []
for t in microbe_stacked.t.unique():
    numLeft = microbe_stacked[microbe_stacked.t==t]\
                .merge(tripartiteNets, on=['bstrain_id'])
    triRichness.append(len(numLeft['bstrain_id']))
    













bAbunds, bKeepTreeStrainsDF, bSpeciesColorDict, _, _, _, _, _  = \
tree.speciesCladeTreePlot(run_id,'microbe',colorjumpB,DBSIM_PATH,DBCLADE_PATH,DBTREE_PATH,\
treepalette,maxticksize,figxy,hratio,babundthreshold)
vAbunds, vKeepTreeStrainsDF, vSpeciesColorDict, _, hlinecVirus, vlinecVirus, _, _ = \
tree.speciesCladeTreePlot(run_id,'virus',colorjumpV,DBSIM_PATH,DBCLADE_PATH,DBTREE_PATH,\
treepalette,maxticksize,figxy,hratio,vabundthreshold)
bAbunds = bAbunds[bAbunds.t <= max(vAbunds.t)]
bAbunds = bAbunds[bAbunds['abundance']>0]
tripartiteNets = pd.read_sql_query("SELECT t, bstrain_id, vstrain_id, time_specific_match_id\
    FROM bstrain_to_vstrain_matches WHERE match_length = 1", conMatch)
spacerMatches = pd.read_sql_query("SELECT t, time_specific_match_id, spacer_id\
    FROM matches_spacers", conMatch)
tripartiteNets = tripartiteNets.merge(spacerMatches,on=['t','time_specific_match_id'])\
            .drop(columns=['time_specific_match_id'])
microbe_stacked = pd.read_sql_query("SELECT t,bstrain_id,abundance FROM babundance \
                    WHERE run_id = {}".format(run_id), conSim)
virus_stacked = pd.read_sql_query("SELECT t,vstrain_id,abundance FROM vabundance \
                    WHERE run_id = {}".format(run_id), conSim)
tripartiteNets = tripartiteNets.merge(microbe_stacked,on=['t','bstrain_id'])\
                    .rename(columns={'abundance':'babundance'})\
                    .merge(virus_stacked,on=['t','vstrain_id'])\
                    .rename(columns={'abundance':'vabundance'})\
                    .merge(virus_total,on=['t'])\
                    .merge(microbe_total,on=['t'])\
                    .rename(columns={'Viral Abundance':'vtotal',\
                                'Microbial Abundance':'btotal'})
spacerHeat = tripartiteNets[['t','bstrain_id','spacer_id','babundance']]\
                .drop_duplicates().groupby(['t','spacer_id'])\
                .agg(babundance=('babundance','sum')).reset_index()
spacerHeat = spacerHeat.sort_values(by=['t','babundance'])
newSpacerID = {}
idNew = 1
for spacerID in spacerHeat['spacer_id'].unique():
    if spacerID not in np.array([*newSpacerID.keys()]):
        newSpacerID[spacerID] = idNew
        idNew += 1



spacerHeat['spacer_id_new'] = spacerHeat['spacer_id'].copy()
spacerHeat = spacerHeat.replace({'spacer_id_new': newSpacerID})
spacerHeat = spacerHeat.sort_values(by=['t','spacer_id_new'])
extraTimes = list({*bAbunds['t']} - {*heat['t']})
df = pd.DataFrame({'t':extraTimes,\
        'spacer_id':len(extraTimes)*[spacerHeat['spacer_id'][0]],\
        'babundance':len(extraTimes)*[0],\
        'spacer_id_new':len(extraTimes)*[spacerHeat['spacer_id_new'][0]]})
spacerHeat = pd.concat([df,spacerHeat])
heat = spacerHeat.pivot_table(index='spacer_id_new', columns='t', values='babundance').fillna(0)
fig, ax = plt.subplots(2,sharex=True, gridspec_kw={'height_ratios': [1,6]})
axes = [ax[0], ax[0].twinx(), ax[1]]#, ax[2], ax[3]]
microbe_stacked = bAbunds.pivot(index='t',columns='tree_bstrain_id',values='abundance')
microbe_stacked.plot.area(ax = axes[0],stacked=True,legend=False, linewidth=0,color=bSpeciesColorDict,sort_columns=True)
axes[0].set_ylabel(ylabel ='Microbial Strain Abundances',labelpad=15,fontsize=7)
axes[0].set_xlabel(xlabel = 'Time t',fontsize=7)
axes[0].ticklabel_format(axis = 'y',style='sci',scilimits=(0,0))
axes[0].xaxis.set_minor_locator(ticker.MultipleLocator(25))
lim = axes[0].get_ylim()
axes[0].set_ylim(0,lim[1])
axes[0].set_xlim(0,1056)
axes[1].plot(virus_total['t'],virus_total['Viral Abundance'],linewidth=0,color='grey')
axes[1].fill_between(virus_total['t'],virus_total['Viral Abundance'], color='grey',alpha=0.6)
lim = axes[1].get_ylim()
axes[1].set_ylim(0,lim[1])
axes[1].set_xlim(0,1056)
spacersOrdered  = spacerHeat[['spacer_id','spacer_id_new']].drop_duplicates()\
                        .sort_values(by=['spacer_id_new'])
spacersOrdered = np.array(spacersOrdered['spacer_id'])
axes[2].pcolormesh(heat,cmap='turbo')
divider = make_axes_locatable(axes[2])
axes[2].set_xlim(0,1056)
axes[2].set_ylim(0,max(heat.index))
axes[2].margins(x=0)
axes[2].margins(y=0)
axes[2].set_yticks(np.arange(1-0.5,max(heat.index),1))
axes[2].set_yticklabels(spacersOrdered,fontsize=8)
[axes[2].axhline(y=i-0.5, linestyle='-',color='white',alpha=0.4,lw=.5) for i in heat.index]
axes[2].set_ylabel(ylabel ='Spacer ID',labelpad=15,fontsize=10)
axes[2].set_xlabel(xlabel ='Time t',labelpad=15,fontsize=10)
[axes[2].axvline(x=i, linestyle='--',color='white',alpha=1,lw=.5) for i in times]
[axes[1].axvline(x=i, linestyle='--',color='black',alpha=1,lw=.5) for i in times]
fig.tight_layout()
for t in times:
    strainHist =\
    tripartiteNets[tripartiteNets.t == t][['t','vstrain_id','spacer_id']]\
    .merge(vKeepTreeStrainsDF,on=['vstrain_id'])
    strainHist['spacer_id_new'] = strainHist['spacer_id'].copy()
    strainHist = strainHist.replace({'spacer_id_new': newSpacerID})
    strainHist = strainHist.sort_values(by=['new_tree_vstrain_id','spacer_id_new'])
    newStrainID = {}
    idNew = 1
    for strainID in strainHist['new_tree_vstrain_id'].unique():
        if strainID not in np.array([*newStrainID.keys()]):
            newStrainID[strainID] = idNew
            idNew += 1
    strainHist = strainHist.replace({'new_tree_vstrain_id': newStrainID})\
                    .sort_values(by=['new_tree_vstrain_id','spacer_id_new'])
    fig, ax = plt.subplots(1)
    ax.scatter(strainHist['new_tree_vstrain_id'],strainHist['spacer_id_new'],lw=1,
                color = 'blue',marker='o',s=50)
    spacersOrdered  = strainHist[['spacer_id','spacer_id_new']].drop_duplicates()\
                            .sort_values(by=['spacer_id_new'])
    ax.set_yticks(np.array(strainHist['spacer_id_new'].unique()))
    spacersOrdered = np.array(spacersOrdered['spacer_id'])
    ax.set_yticklabels(spacersOrdered,fontsize=8)
    ax.set_xticks(np.array(strainHist['new_tree_vstrain_id'].unique()))
    strainsOrdered = np.array(strainHist['vstrain_id'].unique())
    ax.set_xticklabels(strainsOrdered,fontsize=8)
    [ax.axhline(y=i, linestyle='-',color='grey',alpha=0.5,lw=.5) \
                for i in np.array(strainHist['spacer_id_new'].unique())]
    [ax.axvline(x=i, linestyle='-',color='grey',alpha=0.5,lw=.5) \
                for i in np.array(strainHist['new_tree_vstrain_id'].unique())]
    ax.set_xlim(0,max(strainHist['new_tree_vstrain_id'])+1)
    ax.set_ylim(0,max(strainHist['spacer_id_new'])+1)
    ax.margins(x=0)
    ax.margins(y=0)
    ax.set_ylabel(ylabel ='Spacer ID',labelpad=15,fontsize=15)
    ax.set_xlabel(xlabel = 'Viral Strain ID',labelpad=15,fontsize=15)
    ax.set_title('Viral Single Matches @t = {}'\
                    .format(t), pad = 25, fontsize = 15)
    fig.tight_layout()
    fig.savefig(os.path.join('/Users/armun/Desktop/spacers',\
        't{}_viral-single-matches.png'.format(t)),dpi=resolve)
    spacerHist =\
    tripartiteNets[tripartiteNets.t == t][['spacer_id','babundance']]\
        .drop_duplicates().groupby(['spacer_id'])\
        .agg(babundance=('babundance','sum')).reset_index()
    spacerHist['spacer_id_new'] = spacerHist['spacer_id'].copy()
    spacerHist = spacerHist.replace({'spacer_id_new': newSpacerID})
    spacerHist = spacerHist.sort_values(by=['spacer_id_new'])
    fig, ax = plt.subplots(1)
    ax.bar(spacerHist['spacer_id_new'],spacerHist['babundance']/K,color='darkred')
    ax.set_xticks(np.array(spacerHist['spacer_id_new'].unique()))
    ax.set_xticklabels(np.array(spacerHist['spacer_id'].unique()),fontsize=8)
    ax.set_ylabel(ylabel ='Host Frequency (1/K)',labelpad=15,fontsize=15)
    ax.set_xlabel(xlabel = 'Spacer ID',labelpad=15,fontsize=15)
    ax.set_title('Host Frequency Distribution\nof Spacers @t = {}'\
                    .format(t), pad = 25, fontsize = 15)
    fig.tight_layout()
    fig.savefig(os.path.join('/Users/armun/Desktop/spacers',\
        't{}_spacer-freq.png'.format(t)),dpi=resolve)
    spacerHist = spacerHist.sort_values(by=['babundance'])
    fig, ax = plt.subplots(1)
    ax.bar(np.arange(1,len(spacerHist['spacer_id_new'])+1,1)[::-1],\
            spacerHist['babundance']/K,color='darkgreen')
    ax.set_xticks(np.arange(1,len(spacerHist['spacer_id_new'])+1,1)[::-1])
    ax.set_xticklabels(np.array(spacerHist['spacer_id'].unique())[::-1],fontsize=8)
    ax.set_ylabel(ylabel ='Host Frequency (1/K)',labelpad=15,fontsize=15)
    ax.set_xlabel(xlabel = 'Spacer ID',labelpad=15,fontsize=15)
    ax.set_title('Ordered Host Frequency Distribution\nof Spacers @t = {}'\
                    .format(t), pad = 25, fontsize = 15)
    fig.tight_layout()
    fig.savefig(os.path.join('/Users/armun/Desktop/spacers',\
                't{}_ordered-spacer-freq.png'.format(t)),dpi=resolve)
