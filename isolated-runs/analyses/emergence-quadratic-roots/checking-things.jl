import matplotlib.cm as cm
import matplotlib.colors as mc
from scipy.special import lambertw
import pandas as pd
import numpy as np
import math
import matplotlib.pyplot as plt
import sys
import os
import seaborn as sns
import sqlite3
import matplotlib.ticker as ticker
import matplotlib.gridspec as grid_spec
# import plotly.graph_objects as go
# import plotly.express as px
# from mpl_toolkits.mplot3d import Axes3D
import itertools

run_id = 3297

resolve = 500
imgType = "png" #or "png"
vabundthreshold = .1

# dir = 'crispr-sweep-7-2-2022/isolates/runID3297-c66-r47'
run = 'runID3297-c66-r47'

# DBSIM_PATH = os.path.join('/Volumes','Yadgah',dir,'{}.sqlite'.format(run))
DBSIM_PATH = os.path.join('/Users/armun/Dropbox/Current/Projects\
/microbe-virus-crispr/stochastic-crispr/repository-2/isolated-runs/isolates\
/runID3297-c66-r47','{}.sqlite'.format(run))

conSim = sqlite3.connect(DBSIM_PATH)
curSim = conSim.cursor()
# Designating plot path from simulation data
ID = curSim.execute('SELECT combo_id,replicate FROM runs WHERE run_id = {}'.format(run_id)).fetchall()
combo_id = ID[0][0]
replicate = ID[0][1]
CC = curSim.execute('SELECT microbe_carrying_capacity FROM param_combos WHERE combo_id = {}'.format(combo_id)).fetchall()
CC = CC[0][0]
RUN_DIR = os.path.join('runID{0}-c{1}-r{2}'.format(run_id,combo_id,replicate))

# DBPROB_PATH = os.path.join('/Volumes','Yadgah',dir,'probability-of-emergence_output.sqlite') # local
# DBPROB_PATH = os.path.join('/Volumes','Yadgah',dir,'emergence-lambert-root_output.sqlite')
# DBMATCH_PATH = os.path.join('/Volumes','Yadgah',dir,'matches_output.sqlite') # local
# DBTREE_PATH = os.path.join('/Volumes','Yadgah',dir,'trees_output.sqlite') # local
# DBTRI_PATH = os.path.join('/Volumes','Yadgah',dir,'tripartite-networks_output.sqlite') # local
# PLOT_PATH = os.path.join('/Volumes','Yadgah') # local

DBMATCH_PATH = os.path.join('/Users/armun/Dropbox/Current/Projects\
/microbe-virus-crispr/stochastic-crispr/repository-2/isolated-runs/isolates\
/runID3297-c66-r47','matches_output.sqlite') # local
DBTREE_PATH = os.path.join('/Users/armun/Dropbox/Current/Projects\
/microbe-virus-crispr/stochastic-crispr/repository-2/isolated-runs/isolates\
/runID3297-c66-r47','trees_output.sqlite') # local
DBTRI_PATH = os.path.join('/Users/armun/Dropbox/Current/Projects\
/microbe-virus-crispr/stochastic-crispr/repository-2/isolated-runs/isolates\
/runID3297-c66-r47','tripartite-networks_output.sqlite') # local
DBCLADE_PATH = os.path.join('/Users/armun/Dropbox/Current/Projects\
/microbe-virus-crispr/stochastic-crispr/repository-2/isolated-runs/isolates\
/runID3297-c66-r47','clade-abundances_output.sqlite') # local
PLOT_PATH = os.path.join('/Users/armun/Desktop') # local

# conProb = sqlite3.connect(DBPROB_PATH)
# curProb = conProb.cursor()
conMatch = sqlite3.connect(DBMATCH_PATH)
curMatch = conMatch.cursor()
conTri = sqlite3.connect(DBTRI_PATH)
curTri = conTri.cursor()
conTree = sqlite3.connect(DBTREE_PATH)
curTree = conTree.cursor()

virus_stacked = pd.read_sql_query("SELECT t,vstrain_id,abundance FROM vabundance WHERE run_id = {}".format(run_id), conSim)
virus_total = pd.read_sql_query("SELECT t,viral_abundance FROM summary WHERE run_id = {}".format(run_id), conSim)\
.rename(columns={"viral_abundance": "vtotal"})
virus_total = virus_total[virus_total.t <= max(virus_stacked.t)]
# THIS IS FOR THE 3D DISTRIBUTIONS....
vMaxIDs = virus_stacked.set_index('t').groupby(['vstrain_id']).agg(t = ('abundance','idxmax'),\
                            vmaxAbund = ('abundance','max')).reset_index()
vMaxIDs = virus_total.merge(vMaxIDs,on=['t'])
vMaxIDs['vtotal'] = vabundthreshold*np.array(vMaxIDs['vtotal'])
keepStrains = list(vMaxIDs[vMaxIDs['vmaxAbund']>vMaxIDs['vtotal']]['vstrain_id'].values)
virus_stacked = virus_stacked[[(i in keepStrains) for i in virus_stacked.vstrain_id]].reset_index(drop=True)
vmatchAbund = virus_stacked[[(i in keepStrains) for i in virus_stacked.vstrain_id]].reset_index(drop=True)
virus_stacked = virus_stacked.pivot(index='t',columns='vstrain_id',values='abundance')

vmatchstrain = pd.read_sql_query("SELECT t, match_id, vstrain_id \
    FROM vmatches", conTri)\
    .rename(columns={"match_id":"vmatch_id"})
vmatchAbund = vmatchAbund.merge(vmatchstrain,on=['t','vstrain_id']).groupby(['t','vmatch_id'])\
.agg(abundance=('abundance','sum')).reset_index()
vMaxIDs = vmatchAbund.set_index('t').groupby(['vmatch_id']).agg(t = ('abundance','idxmax'),\
                            vmaxAbund = ('abundance','max')).reset_index()
vMaxIDs = virus_total.merge(vMaxIDs,on=['t'])
vMaxIDs['vtotal'] = vabundthreshold*np.array(vMaxIDs['vtotal'])
keepStrains = list(vMaxIDs[vMaxIDs['vmaxAbund']>vMaxIDs['vtotal']]['vmatch_id'].values)
vmatchAbund = vmatchAbund[[(i in keepStrains) for i in vmatchAbund.vmatch_id]].reset_index(drop=True)



######
vmatch_stacked = vmatchAbund.pivot(index='t',columns='vmatch_id',values='abundance')
pal = sns.color_palette("tab20c")
fig, ax = plt.subplots(1,sharex=True)
vmatch_stacked.plot.area(ax = ax,stacked=True,legend=False, linewidth=0,color=pal,sort_columns=True)


vmatchSabund = pd.read_sql_query("SELECT t,match_id,bsusceptible \
        FROM vmatches_abundances", conTri)\
        .rename(columns={"match_id":"vmatch_id"})\
        .rename(columns={"bsusceptible":"babundance"})
vmatchSabund = vmatchSabund.merge(vmatchAbund,on=['t','vmatch_id'])\
.drop(columns=['abundance']).drop_duplicates()
vmatchS_stacked = vmatchSabund.pivot(index='t',columns='vmatch_id',values='babundance')
pal = sns.color_palette("tab20c")
fig, ax = plt.subplots(1,sharex=True)
vmatchS_stacked.plot.area(ax = ax,stacked=True,legend=False, linewidth=0,color=pal,sort_columns=True)
######


tmin = 300
tmax = 600
vmatchSabund = pd.read_sql_query("SELECT t,match_id,bsusceptible \
        FROM vmatches_abundances", conTri)\
        .rename(columns={"match_id":"vmatch_id"})\
        .rename(columns={"bsusceptible":"babundance"})
vmatchSabund = vmatchSabund.merge(vmatchAbund,on=['t','vmatch_id'])\
.drop(columns=['abundance']).drop_duplicates()
vmatchS_stacked = vmatchSabund[(vmatchSabund.t >= tmin) & (vmatchSabund.t <= tmax)]\
                    .pivot(index='t',columns='vmatch_id',values='babundance')
vmatch_stacked = vmatchAbund[(vmatchAbund.t >= tmin) & (vmatchAbund.t <= tmax)]\
                    .pivot(index='t',columns='vmatch_id',values='abundance')
###
# vmatchID = 705
# vmatchS_stacked = vmatchSabund[\
#                     (vmatchSabund.t >= tmin) & (vmatchSabund.t <= tmax) & (vmatchSabund.vmatch_id == vmatchID)]\
#                     .pivot(index='t',columns='vmatch_id',values='babundance')
# vmatch_stacked = vmatchAbund[\
#                     (vmatchAbund.t >= tmin) & (vmatchAbund.t <= tmax) & (vmatchSabund.vmatch_id == vmatchID)]\
#                     .pivot(index='t',columns='vmatch_id',values='abundance')
###
pal = sns.color_palette("tab20c")
fig, ax = plt.subplots(1,sharex=True)
cmapb = cm.tab20b.colors                        # get a specific colormap
cmapc = cm.tab20c.colors                          # extract all colors
cmapNew = cmapc + cmapb              # repeat X times, here X = 2
cmapJoint = mc.LinearSegmentedColormap.from_list('cmap_new', np.array(cmapNew))
vmatch_stacked.plot.area(ax = ax,stacked=True,legend=True, linewidth=0,cmap=cmapJoint,sort_columns=True)
ax.legend(loc='center right',prop={'size': 2})
vmatch_stacked.plot(stacked=True, ax=axes[0], legend=False, color='white',sort_columns=True,linewidth=.2)
ax.set_ylabel(ylabel ='Viral MatchID Abundances',labelpad=10,fontsize=10)
ax.set_xlabel(xlabel = 'Time t',fontsize=15)
ax.ticklabel_format(axis = 'y',style='sci',scilimits=(0,0))
ax.xaxis.set_minor_locator(ticker.MultipleLocator(25))
ax.tick_params(axis='x', labelsize= 10)
ax.tick_params(axis='y', labelsize= 10)
fig.tight_layout()
fig.savefig(os.path.join(PLOT_PATH,'plot1'),dpi=resolve)
fig, ax = plt.subplots(1,sharex=True)
vmatchS_stacked.plot.area(ax = ax,stacked=True,legend=True, linewidth=0,cmap=cmapJoint,sort_columns=True)
ax.legend(loc='center right',prop={'size': 2})
vmatchS_stacked.plot(stacked=True, ax=ax, legend=False, color='white',sort_columns=True,linewidth=.2)
ax.set_ylabel(ylabel ='Susceptible-to-MatchID Abundances',labelpad=10,fontsize=10)
ax.set_xlabel(xlabel = 'Time t',fontsize=10)
ax.ticklabel_format(axis = 'y',style='sci',scilimits=(0,0))
ax.xaxis.set_minor_locator(ticker.MultipleLocator(25))
ax.tick_params(axis='x', labelsize= 10)
ax.tick_params(axis='y', labelsize= 10)
fig.tight_layout()
fig.savefig(os.path.join(PLOT_PATH,'plot2'),dpi=resolve)



vmatchID = 562
pal = sns.color_palette("tab20c")
microbe_total = pd.read_sql_query("SELECT t,microbial_abundance FROM summary WHERE run_id = {}".format(run_id), conSim)\
.rename(columns={"microbial_abundance": "btotal"})
# vmatchSabund = pd.read_sql_query("SELECT t,match_id,bsusceptible \
#         FROM vmatches_abundances", conTri)\
#         .rename(columns={"match_id":"vmatch_id"})\
#         .rename(columns={"bsusceptible":"babundance"})
vstrain0bstrain = pd.read_sql_query("SELECT t,vstrain_id,bstrain_id \
        FROM bstrain_to_vstrain_0matches", conMatch)
vabunds = pd.read_sql_query("SELECT t,vstrain_id,abundance \
        FROM vabundance", conSim)
babunds = pd.read_sql_query("SELECT t,bstrain_id,abundance \
        FROM babundance", conSim).merge(microbe_total,on=['t'])
babunds['bfreq'] = babunds['abundance']/babunds['btotal']
bstrainAbunds = vstrain0bstrain.merge(vabunds,on=['t','vstrain_id'])\
    .rename(columns={"abundance":"v_infect_total"})\
    .merge(babunds,on=['t','bstrain_id'])\
    .rename(columns={"abundance":"babundance"})\
    .groupby(['t','bstrain_id','babundance'])\
    .agg(v_infect_total = ('v_infect_total','sum')).reset_index()
vstrain0bstrain = vstrain0bstrain.merge(bstrainAbunds,on=['t','bstrain_id'])
vmatches = pd.read_sql_query("SELECT t,match_id,vstrain_id \
        FROM vmatches", conTri).rename(columns={"match_id":"vmatch_id"})
vstrain0bstrain = vstrain0bstrain.merge(vmatches,on=['t','vstrain_id'])\
                    .drop(columns=['vstrain_id']).drop_duplicates()
sampleSAbunds = vstrain0bstrain[vstrain0bstrain['vmatch_id'].isin([vmatchID])]\
                .sort_values(by=['t'])
vmatchAbunds = pd.read_sql_query("SELECT t,match_id,vabundance \
        FROM vmatches_abundances", conTri)\
        .rename(columns={"match_id":"vmatch_id"})
vmatchAbunds = vmatchAbunds[vmatchAbunds.vmatch_id==vmatchID]
# S0 = []
# for bstrainID in sorted(sampleSAbunds['bstrain_id'].unique()):
#     S0.append(sampleSAbunds[sampleSAbunds['bstrain_id'] == bstrainID]['babundance'].values[0])
# S0 = {'col1': sorted(sampleSAbunds['bstrain_id'].unique()), 'S0': S0}
# v0 = sampleSAbunds['v_infect_total'].values[0]
fxnSAbunds= []
for t in sorted(sampleSAbunds['t'].unique()):
    c = 0
    for bstrainID in sorted(sampleSAbunds[sampleSAbunds['t']==t]['bstrain_id']):
        v0 = sampleSAbunds[sampleSAbunds['bstrain_id']==bstrainID]['v_infect_total'].values[0]
        v = sampleSAbunds[(sampleSAbunds['t']==t) & (sampleSAbunds['bstrain_id']==bstrainID)]\
                ['v_infect_total'].values[0]
        # v0 = vmatchAbunds['vabundance'].values[0]
        # v = vmatchAbunds[vmatchAbunds['t']==t]['vabundance'].values[0]
        s0 = sampleSAbunds[sampleSAbunds['bstrain_id']==bstrainID]['babundance'].values[0]
        # print("v0 is {0}, v is {1}, s0 is s{2}".format(v0,v,s0))
        # c += (-v+v0)/50*babunds[(babunds.t==t) & (babunds.bstrain_id==bstrainID)]['bfreq'].values[0]\
                # +s0
        c += (-v+v0+s0*50)/(50)
        # c = (-1*.01*S0+(10**-7)-v*(10**-7)+S0*50*(10**-7))/(-1*.01+50*(10**-7))
        #l = (.00000000001*np.log(v) - .00000000001*np.log(v0))/(50*(10**-7))
    s = c #+ l
    fxnSAbunds.append(s)
# fxnSAbunds = []
# S0 = sampleSAbunds[sampleSAbunds['t']==381]['babundance'].values[0]
# v0 = samplevmatchAbund[samplevmatchAbund['t']==381]['abundance'].values[0]
# for t in sorted(sampleSAbunds['t'].unique()):
#     v = samplevmatchAbund[samplevmatchAbund['t']==t]['abundance'].values[0]
#     arg = (-1/.01)*np.exp(((10**-7)/.01)*(v-v0-50*S0))*S0*(v**(-1/.01))*(v0**(1/.01))*50*(10**-7)
#     s = -.01/(50*(10**-7))*lambertw(arg,k=0)
#     fxnSAbunds.append(s)
# fxnSAbunds = []
# S0 = sampleSAbunds[sampleSAbunds['t']==428]['babundance'].values[0]
# v0 = samplevmatchAbund[samplevmatchAbund['t']==428]['abundance'].values[0]
# # v0 = 2
# for t in sorted(sampleSAbunds['t'].unique()):
#     v = samplevmatchAbund[samplevmatchAbund['t']==t]['abundance'].values[0]
#     s = np.exp((1/.01)*(.01*np.log(S0)-np.log(v)+np.log(v0)))
#     fxnSAbunds.append(s)
btotal = sampleSAbunds.groupby(['t','vmatch_id'])\
.agg(btotal = ('babundance','sum')).reset_index()
fig, ax = plt.subplots(1,sharex=True)
ax.plot(sorted(sampleSAbunds['t'].unique()),btotal['btotal'],linewidth=1)
ax.plot(sorted(sampleSAbunds['t'].unique()),fxnSAbunds,linewidth=1)
fig.tight_layout()
fig.savefig(os.path.join(PLOT_PATH,'plot3'),dpi=resolve)


sampleSAbunds = sampleSAbunds.pivot(index='t',columns='vmatch_id',values='babundance')
fig, ax = plt.subplots(1,sharex=True)
ax.plot(sampleSAbunds,linewidth=1)

# vmatchIabund = pd.read_sql_query("SELECT t,match_id,bimmune \
#         FROM vmatches_abundances", conTri)\
#         .rename(columns={"match_id":"vmatch_id"})\
#         .rename(columns={"bimmune":"babundance"})
# vmatchIabund = vmatchIabund.merge(vmatchAbund,on=['t','vmatch_id'])\
# .drop(columns=['abundance']).drop_duplicates()
# vmatchI_stacked = vmatchIabund.pivot(index='t',columns='vmatch_id',values='babundance')
# pal = sns.color_palette("tab20c")
# fig, ax = plt.subplots(1,sharex=True)
# ax.plot(vmatchI_stacked,linewidth=1)
# vmatchI_stacked.plot.area(ax = ax,stacked=True,legend=False, linewidth=0,color=pal,sort_columns=True)

# fig, ax = plt.subplots(3,sharex=True)
# axes = [ax[0], ax[1], ax[2]]
# cmapb = cm.tab20b.colors                        # get a specific colormap
# cmapc = cm.tab20c.colors                          # extract all colors
# cmapNew = cmapc + cmapb              # repeat X times, here X = 2
# cmapJoint = mc.LinearSegmentedColormap.from_list('cmap_new', np.array(cmapNew))
# vmatch_stacked.plot.area(ax = axes[0],stacked=True,legend=True, linewidth=0,cmap=cmapJoint,sort_columns=True)
# vmatch_stacked.plot(stacked=True, ax=axes[0], legend=False, color='white',sort_columns=True,linewidth=.15)
# vmatchS_stacked.plot.area(ax = axes[1],stacked=True,legend=False, linewidth=0,cmap=cmapJoint,sort_columns=True)
# vmatchS_stacked.plot(stacked=True, ax=axes[1], legend=False, color='white',sort_columns=True,linewidth=.15)
# vmatchI_stacked.plot.area(ax = axes[2],stacked=True,legend=False, linewidth=0,cmap=cmapJoint,sort_columns=True)
# vmatchI_stacked.plot(stacked=True, ax=axes[2], legend=False, color='white',sort_columns=True,linewidth=.15)



# vfreq = pd.read_sql_query("SELECT t, vmatch_id, vfrequency \
#     FROM existing_vmatch_extinction ORDER BY t", conProb)\
#     .rename(columns={"vfrequency":"vfreq"})\
#     .merge(vmatchSabund,on=['t','vmatch_id'])
# vfreq['bexpected'] = vfreq['vfreq']*vfreq['babundance']
# vfreq = vfreq.groupby(['t'])\
#     .agg(b_exp=('bexpected','sum')).reset_index()
# fig, ax = plt.subplots(1,sharex=True)
# ax.plot(vfreq['t'],vfreq['b_exp'],linewidth=3,color='grey')












pEmergeExpected = pd.read_sql_query("SELECT t, plambert_emerge_weighted \
    FROM existing_vmatch_extinction ORDER BY t", conProb).groupby(['t'])\
    .agg(p_exp=('plambert_emerge_weighted','sum')).reset_index()
pal = sns.color_palette("tab20c")
fig, ax = plt.subplots(2,sharex=True)
axes = [ax[0], ax[0].twinx(), ax[1], ax[1].twinx(), ax[1].twinx()]
microbe_stacked.plot.area(ax = axes[0],stacked=True,legend=False, linewidth=0,color=pal,sort_columns=True)
axes[0].set_ylabel(ylabel ='Microbial Strain Abundances',labelpad=15,fontsize=15)
axes[0].set_xlabel(xlabel = 'Time t',fontsize=15)
axes[0].ticklabel_format(axis = 'y',style='sci',scilimits=(0,0))
axes[0].xaxis.set_minor_locator(ticker.MultipleLocator(25))
axes[0].tick_params(axis='x', labelsize= 10)
axes[0].tick_params(axis='y', labelsize= 10)
lim = axes[0].get_ylim()
axes[0].set_ylim(0,lim[1])
axes[1].plot(virus_total['t'],virus_total['vtotal'],linewidth=0,color='grey')
axes[1].fill_between(virus_total['t'],virus_total['vtotal'], color='grey',alpha=0.6)
lim = axes[1].get_ylim()
axes[1].set_ylim(0,lim[1])
axes[2].plot(pEmergeExpected['t'],pEmergeExpected['p_exp'],linewidth=1,color='darkblue')
axes[2].set_ylabel(ylabel ='Probability of Viral Strain Emergence',labelpad=15,fontsize=15)
axes[2].set_xlabel(xlabel = 'Time t',fontsize=15)
axes[2].tick_params(axis='x', labelsize= 10)
axes[2].tick_params(axis='y', labelsize= 10)
axes[3].plot(virus_total['t'],virus_total['vtotal'],linewidth=0,color='grey')
axes[3].fill_between(virus_total['t'],virus_total['vtotal'], color='grey',alpha=0.3)
axes[4].plot(vfreq['t'],vfreq['b_exp'],linewidth=1,color='red')
fig.tight_layout()









pal = sns.color_palette("tab20c")
fig, ax = plt.subplots(2,sharex=True)
axes = [ax[0], ax[1]]
vmatchS_stacked.plot.area(ax = axes[0],stacked=True,legend=False, linewidth=0,color=pal,sort_columns=True)
axes[0].set_ylabel(ylabel ='Susceptible Abundances',labelpad=15,fontsize=15)
axes[0].set_xlabel(xlabel = 'Time t',fontsize=15)
axes[0].ticklabel_format(axis = 'y',style='sci',scilimits=(0,0))
axes[0].xaxis.set_minor_locator(ticker.MultipleLocator(25))
axes[0].tick_params(axis='x', labelsize= 10)
axes[0].tick_params(axis='y', labelsize= 10)
lim = axes[0].get_ylim()
axes[0].set_ylim(0,lim[1])
vmatch_stacked.plot.area(ax = axes[1],stacked=True,legend=False, linewidth=0,color=pal,sort_columns=True)
lim = axes[1].get_ylim()
axes[1].set_ylim(0,lim[1])
axes[1].set_ylabel(ylabel ='Match Abundances',labelpad=15,fontsize=15)
axes[1].set_xlabel(xlabel = 'Time t',fontsize=15)
axes[1].tick_params(axis='x', labelsize= 10)
axes[1].tick_params(axis='y', labelsize= 10)
axes[1].ticklabel_format(axis = 'y',style='sci',scilimits=(0,0))
axes[1].xaxis.set_minor_locator(ticker.MultipleLocator(25))
fig.tight_layout()


pEmergence = pd.read_sql_query("SELECT t, vmatch_id, p_extinction_lambert \
    FROM existing_vmatch_extinction ORDER BY t, vmatch_id", conProb)
pEmergence['p_emergence'] = 1 - np.array(pEmergence['p_extinction_lambert'])
pEmergence = pEmergence.drop(columns=['p_extinction_lambert'])
pEmergence = pEmergence.merge(vmatchstrain,on=['t','vmatch_id'])
pEmergence = pEmergence.merge(virusAbund,on=['t','vstrain_id'])

pEmergenceStacked = pEmergence.drop(columns=['vmatch_id'])
pEmergenceStacked  = pEmergenceStacked[[(i in keepStrains) for i in pEmergenceStacked.vstrain_id]].reset_index(drop=True)
pEmergenceStacked = pEmergenceStacked.pivot_table(index='t', columns='vstrain_id', values='p_emergence').fillna(0)










pEmergenceSTD= pEmergence.drop(columns=['vmatch_id'])
pEmergeExpected = pd.read_sql_query("SELECT t, vmatch_id, p_extinction_lambert, vfrequency \
    FROM existing_vmatch_extinction ORDER BY t", conProb)\
    .merge(frequencies,on=['t','vmatch_id']).rename(columns={"vfrequency":"vfreq"})
pEmergeExpected['p_extinction_lambert'] = 1- np.array(pEmergeExpected['p_extinction_lambert'])
pEmergeExpected=pEmergeExpected.rename(columns={"p_extinction_lambert":"plambert_emerge_weighted"})
pEmergeExpected['plambert_emerge_weighted'] = \
    pEmergeExpected['plambert_emerge_weighted']*pEmergeExpected['vfreq']
pEmergeExpected = pEmergeExpected.groupby(['t'])\
    .agg(p_exp=('plambert_emerge_weighted','sum')).reset_index()
pEmergenceSTD = pEmergenceSTD.merge(pEmergeExpected,on=['t'])
pEmergenceSTD['p_std'] = (pEmergenceSTD['p_emergence']-pEmergenceSTD['p_exp'])**2
pStrains= pEmergenceSTD.groupby(['t']).agg(numStrains=('p_emergence','size')).reset_index()
pEmergenceSTD = pEmergenceSTD.groupby(['t']).agg(p_std=('p_std','sum')).reset_index()\
                .merge(pStrains,on=['t'])
pEmergenceSTD['p_std'] = pEmergenceSTD['p_std']/pEmergenceSTD['numStrains']
pEmergenceSTD['p_std'] = np.sqrt(pEmergenceSTD['p_std'])


# std of p_emergence
pal = sns.color_palette("tab20c")
fig, ax = plt.subplots(2,sharex=True)
axes = [ax[0], ax[1]]
virus_stacked.plot.area(ax = axes[0],stacked=True,legend=False, linewidth=0,color=pal,sort_columns=True)
axes[0].set_ylabel(ylabel ='Viral Strain Abundances',labelpad=15,fontsize=15)
axes[0].set_xlabel(xlabel = 'Time t',fontsize=15)
axes[0].ticklabel_format(axis = 'y',style='sci',scilimits=(0,0))
axes[0].xaxis.set_minor_locator(ticker.MultipleLocator(25))
axes[0].tick_params(axis='x', labelsize= 10)
axes[0].tick_params(axis='y', labelsize= 10)
lim = axes[0].get_ylim()
axes[0].set_ylim(0,lim[1])
axes[1].plot(pEmergenceSTD['t'],pEmergenceSTD['p_std'],linewidth=1,color='darkblue')
lim = axes[1].get_ylim()
axes[1].set_ylim(0,lim[1])
axes[1].set_ylabel(ylabel ='stand. dev. of p_emergence',labelpad=15,fontsize=15)
axes[1].set_xlabel(xlabel = 'Time t',fontsize=15)
axes[1].tick_params(axis='x', labelsize= 10)
axes[1].tick_params(axis='y', labelsize= 10)
axes[1].ticklabel_format(axis = 'y',style='sci',scilimits=(0,0))
axes[1].xaxis.set_minor_locator(ticker.MultipleLocator(25))
fig.tight_layout()





# stacked p_emergence fo strains
pal = sns.color_palette("tab20c")
fig, ax = plt.subplots(2,sharex=True)
axes = [ax[0], ax[1]]
virus_stacked.plot.area(ax = axes[0],stacked=True,legend=False, linewidth=0,color=pal,sort_columns=True)
axes[0].set_ylabel(ylabel ='Viral Strain Abundances',labelpad=15,fontsize=15)
axes[0].set_xlabel(xlabel = 'Time t',fontsize=15)
axes[0].ticklabel_format(axis = 'y',style='sci',scilimits=(0,0))
axes[0].xaxis.set_minor_locator(ticker.MultipleLocator(25))
axes[0].tick_params(axis='x', labelsize= 10)
axes[0].tick_params(axis='y', labelsize= 10)
lim = axes[0].get_ylim()
axes[0].set_ylim(0,lim[1])
pEmergenceStacked.plot.area(ax = axes[1],stacked=True,legend=False, linewidth=0,color=pal,sort_columns=True)
axes[1].set_ylabel(ylabel ='p_emergence',labelpad=15,fontsize=15)
axes[1].set_xlabel(xlabel = 'Time t',fontsize=15)
axes[1].ticklabel_format(axis = 'y',style='sci',scilimits=(0,0))
axes[1].xaxis.set_minor_locator(ticker.MultipleLocator(25))
axes[1].tick_params(axis='x', labelsize= 10)
axes[1].tick_params(axis='y', labelsize= 10)
lim = axes[1].get_ylim()
axes[1].set_ylim(0,lim[1])
fig.tight_layout()


# SPACER DIVERSITY PLOT
pal = sns.color_palette("tab20c")
fig, ax = plt.subplots(2,sharex=True)
axes = [ax[0], ax[0].twinx(), ax[1], ax[1].twinx()]
microbe_stacked.plot.area(ax = axes[0],stacked=True,legend=False, linewidth=0,color=pal,sort_columns=True)
axes[0].set_ylabel(ylabel ='Microbial Strain Abundances',labelpad=15,fontsize=15)
axes[0].set_xlabel(xlabel = 'Time t',fontsize=15)
axes[0].ticklabel_format(axis = 'y',style='sci',scilimits=(0,0))
axes[0].xaxis.set_minor_locator(ticker.MultipleLocator(25))
axes[0].tick_params(axis='x', labelsize= 10)
axes[0].tick_params(axis='y', labelsize= 10)
lim = axes[0].get_ylim()
axes[0].set_ylim(0,lim[1])
axes[1].plot(virus_total['t'],virus_total['vtotal'],linewidth=0,color='grey')
axes[1].fill_between(virus_total['t'],virus_total['vtotal'], color='grey',alpha=0.6)
lim = axes[1].get_ylim()
axes[1].set_ylim(0,lim[1])
axes[2].plot(expBFreq['t'],expBFreq['shanDiv'],linewidth=1,color='darkblue')
axes[2].set_ylabel(ylabel ='Spacer Shannon Diversity',labelpad=15,fontsize=15)
axes[2].set_xlabel(xlabel = 'Time t',fontsize=15)
axes[2].tick_params(axis='x', labelsize= 10)
axes[2].tick_params(axis='y', labelsize= 10)
axes[3].plot(virus_total['t'],virus_total['vtotal'],linewidth=0,color='grey')
axes[3].fill_between(virus_total['t'],virus_total['vtotal'], color='grey',alpha=0.3)
fig.tight_layout()


# TARGETED SPACER DIVERSITY PLOT
pal = sns.color_palette("tab20c")
fig, ax = plt.subplots(2,sharex=True)
axes = [ax[0], ax[0].twinx(), ax[1], ax[1].twinx()]
microbe_stacked.plot.area(ax = axes[0],stacked=True,legend=False, linewidth=0,color=pal,sort_columns=True)
axes[0].set_ylabel(ylabel ='Microbial Strain Abundances',labelpad=15,fontsize=15)
axes[0].set_xlabel(xlabel = 'Time t',fontsize=15)
axes[0].ticklabel_format(axis = 'y',style='sci',scilimits=(0,0))
axes[0].xaxis.set_minor_locator(ticker.MultipleLocator(25))
axes[0].tick_params(axis='x', labelsize= 10)
axes[0].tick_params(axis='y', labelsize= 10)
lim = axes[0].get_ylim()
axes[0].set_ylim(0,lim[1])
axes[1].plot(virus_total['t'],virus_total['vtotal'],linewidth=0,color='grey')
axes[1].fill_between(virus_total['t'],virus_total['vtotal'], color='grey',alpha=0.6)
lim = axes[1].get_ylim()
axes[1].set_ylim(0,lim[1])
axes[2].plot(expBFreq['t'],expBFreq['shanDivRatio'],linewidth=1,color='darkblue')
axes[2].set_ylabel(ylabel ='Proportion of Spacers Targeted',labelpad=15,fontsize=15)
axes[2].set_xlabel(xlabel = 'Time t',fontsize=15)
axes[2].tick_params(axis='x', labelsize= 10)
axes[2].tick_params(axis='y', labelsize= 10)
axes[3].plot(virus_total['t'],virus_total['vtotal'],linewidth=0,color='grey')
axes[3].fill_between(virus_total['t'],virus_total['vtotal'], color='grey',alpha=0.3)
fig.tight_layout()




pEmergeStacked = pEmergence.merge(tripartiteNets,on=['t','vmatch_id'])\
    .drop(columns=['spacer_id','bmatch_id','p_extinction']).drop_duplicates()
pEmergeStacked = pEmergeStacked[pEmergeStacked['p_emergence']>0]
# pEmergeStacked['p_emergence']= list(map(np.log,list(pEmergeStacked['p_emergence'])))
pEmergeStacked = pEmergeStacked.pivot_table(index='t', columns='vmatch_id', values='p_emergence')
pal = sns.color_palette("tab20c")
fig, ax = plt.subplots(1,sharex=True)
pEmergeStacked.plot(ax=ax,linewidth=1,color=pal)
pEmergeStacked.plot.area(ax = ax,stacked=True,legend=False, linewidth=0,color=pal,sort_columns=True)
ax.set_ylabel(ylabel ='P_emergence',labelpad=15,fontsize=15)
ax.set_xlabel(xlabel = 'Time t',fontsize=15)
ax.ticklabel_format(axis = 'y',style='sci',scilimits=(0,0))
ax.xaxis.set_minor_locator(ticker.MultipleLocator(25))
ax.tick_params(axis='x', labelsize= 10)
ax.tick_params(axis='y', labelsize= 10)
lim = ax.get_ylim()
ax.set_ylim(0,lim[1])




# NO P_EMERGENCE
spacerHeatPB_noV = tripartiteNets.drop(columns=['vmatch_id']).drop_duplicates()\
    .merge(bmatchstrain,on=['t','bmatch_id']).drop(columns=['bmatch_id'])\
    .merge(v0bmatches,on=['t','bstrain_id']).drop(columns=['bstrain_id'])\
    .drop_duplicates()
norm = spacerHeatPB_noV.groupby(['t']).agg(norm=('p_emergence','sum')).reset_index()
spacerHeatPB_noV = spacerHeatPB_noV.merge(norm,on=['t'])
spacerHeatPB_noV['p_emergence'] =  \
    [spacerHeatPB_noV['p_emergence'][i]/spacerHeatPB_noV['norm'][i] if spacerHeatPB_noV['norm'][i] != 0 else 0 \
    for i in range(0,len(spacerHeatPB_noV['p_emergence']))]
spacerHeatPB_noV = spacerHeatPB_noV.groupby(['t','spacer_id'])\
            .agg(p_spacer=('p_emergence','sum')).reset_index()
bfreq = tripartiteNets.drop(columns=['vmatch_id']).drop_duplicates()\
.merge(\
pd.read_sql_query("SELECT t,match_id,babundance \
    FROM bmatches_abundances",conTri).rename(columns={"match_id":"bmatch_id"}),\
on=['t','bmatch_id']).groupby(['t','spacer_id']).agg(bfreq=('babundance','sum'))\
.reset_index()
bfreq = bfreq.groupby(['t']).agg(btotal=('bfreq','sum')).reset_index()\
        .merge(bfreq,on=['t'])
bfreq['bfreq'] = bfreq['bfreq']/bfreq['btotal']
bfreq = bfreq.drop(columns=['btotal'])
spacerHeatPB_noV = spacerHeatPB_noV.merge(bfreq,on=['t','spacer_id'])
spacerHeatPB_noV['exp_freq'] = spacerHeatPB_noV['bfreq'].copy()
expBFreq = spacerHeatPB_noV.drop(columns=['p_spacer','bfreq']).groupby(['t'])\
            .agg(exp_freq=('exp_freq','sum')).reset_index()
shannon = []
for t in sorted(spacerHeatPB_noV['t'].unique()):
    div = list(spacerHeatPB_noV[(spacerHeatPB_noV['t']==t) & (spacerHeatPB_noV['exp_freq']!=0)]['exp_freq'])
    shannon.append(np.exp(-1*sum(np.array(div)*np.array(list(map(np.log,div))))))


expBFreq['shanDivBurn'] = shannon[:]
shannon = []
for t in sorted(spacerHeatPB_noV['t'].unique()):
    div = list(bfreq[(bfreq['t']==t)]['bfreq'])
    shannon.append(np.exp(-1*sum(np.array(div)*np.array(list(map(np.log,div))))))


expBFreq['shanDiv'] = shannon[:]

expBFreq['shanDivRatio'] = np.array(expBFreq['shanDivBurn'])/np.array(expBFreq['shanDiv'])


# spacerHeatPB_noV = spacerHeatPB_noV[spacerHeatPB_noV.exp_freq !=0 ]
# spacerHeatPB_noV['exp_freq'] = list(map(np.log,spacerHeatPB_noV.exp_freq))
bfreq = bfreq.replace({"spacer_id": newSpacerIDs})\
            .pivot_table(index='spacer_id', columns='t', values='bfreq').fillna(0)

spacerHeatPB_noV = spacerHeatPB_noV.replace({"spacer_id": newSpacerIDs})\
            .pivot_table(index='spacer_id', columns='t', values='p_spacer').fillna(0)



pal = sns.color_palette("tab20c")
fig, ax = plt.subplots(2,sharex=True, gridspec_kw={'height_ratios': [1,4]})
axes = [ax[0], ax[0].twinx(), ax[1]]#, ax[2], ax[3]]
microbe_stacked.plot.area(ax = axes[0],stacked=True,legend=False, linewidth=0,color=pal,sort_columns=True)
axes[0].set_ylabel(ylabel ='Microbial Strain Abundances',labelpad=15,fontsize=7)
axes[0].set_xlabel(xlabel = 'Time t',fontsize=7)
axes[0].ticklabel_format(axis = 'y',style='sci',scilimits=(0,0))
axes[0].xaxis.set_minor_locator(ticker.MultipleLocator(25))
lim = axes[0].get_ylim()
axes[0].set_ylim(0,lim[1])
axes[1].plot(virus_total[virus_total.t<=max(microbe_stacked.index)]['t'],virus_total[virus_total.t<=max(microbe_stacked.index)]['vtotal'],linewidth=0,color='grey')
axes[1].fill_between(virus_total[virus_total.t<=max(microbe_stacked.index)]['t'],virus_total[virus_total.t<=max(microbe_stacked.index)]['vtotal'], color='grey',alpha=0.6)
lim = axes[1].get_ylim()
axes[1].set_ylim(0,lim[1])
axes[2].imshow(bfreq,cmap='jet', aspect='auto',extent=[0, max(microbe_stacked.index), max(spacerHeatPB_noV.index), min(spacerHeatPB_noV.index)-1])
axes[2].set_yticks(np.arange(0.5,max(spacerHeatPB_noV.index)+0.5,1))
axes[2].set_yticklabels(spacersOrdered,fontsize=4)
axes[2].set_ylabel(ylabel ='Spacer IDs',labelpad=15,fontsize=10)
fig.tight_layout()



pal = sns.color_palette("tab20c")
fig, ax = plt.subplots(3,sharex=True, gridspec_kw={'height_ratios': [1,4,1]})
axes = [ax[0], ax[0].twinx(), ax[1], ax[2], ax[2].twinx()]#, ax[2], ax[3]]
microbe_stacked.plot.area(ax = axes[0],stacked=True,legend=False, linewidth=0,color=pal,sort_columns=True)
axes[0].set_ylabel(ylabel ='Microbial Strain Abundances',labelpad=15,fontsize=7)
axes[0].set_xlabel(xlabel = 'Time t',fontsize=7)
axes[0].ticklabel_format(axis = 'y',style='sci',scilimits=(0,0))
axes[0].xaxis.set_minor_locator(ticker.MultipleLocator(25))
lim = axes[0].get_ylim()
axes[0].set_ylim(0,lim[1])
axes[1].plot(virus_total[virus_total.t<=max(microbe_stacked.index)]['t'],virus_total[virus_total.t<=max(microbe_stacked.index)]['vtotal'],linewidth=0,color='grey')
axes[1].fill_between(virus_total[virus_total.t<=max(microbe_stacked.index)]['t'],virus_total[virus_total.t<=max(microbe_stacked.index)]['vtotal'], color='grey',alpha=0.6)
lim = axes[1].get_ylim()
axes[1].set_ylim(0,lim[1])
axes[2].imshow(bfreq,cmap='jet', aspect='auto',extent=[0, max(microbe_stacked.index), max(spacerHeatPB_noV.index), min(spacerHeatPB_noV.index)-1])
axes[2].set_yticks(np.arange(0.5,max(spacerHeatPB_noV.index)+0.5,1))
axes[2].set_yticklabels(spacersOrdered,fontsize=4)
axes[2].set_ylabel(ylabel ='Spacer IDs',labelpad=15,fontsize=10)
fig.tight_layout()

# axes[1].grid(which='major', color='w', linestyle='-', linewidth=.03)
# ax.set_yticks(np.arange(2, 10, 1), minor=True)
# axes[1].yaxis.grid(True)
# # axes[2].plot(expBFreq['t'],expBFreq['exp_freq'])
axes[3].plot(expBFreq['t'],expBFreq['shanDivRatio'],color='darkblue')
axes[3].set_ylabel(ylabel ='Proportion of Spacers Targeted',labelpad=15,fontsize=7)
axes[4].plot(virus_total[virus_total.t<=max(microbe_stacked.index)]['t'],virus_total[virus_total.t<=max(microbe_stacked.index)]['vtotal'],linewidth=0,color='grey')
axes[4].fill_between(virus_total[virus_total.t<=max(microbe_stacked.index)]['t'],virus_total[virus_total.t<=max(microbe_stacked.index)]['vtotal'], color='grey',alpha=0.5)






spacerHeatPB_noV = tripartiteNets.drop(columns=['vmatch_id']).drop_duplicates()\
    .merge(bmatchstrain,on=['t','bmatch_id']).drop(columns=['bmatch_id'])\
    .merge(v0bmatches,on=['t','bstrain_id']).drop(columns=['bstrain_id'])\
    .drop_duplicates().merge(pEmergence,on=['t','vstrain_id'])

pEmergeExpected = pd.read_sql_query("SELECT t, vmatch_id, p_emerge_weighted \
    FROM vmatch_extinction ORDER BY t", conProb)
tripartiteNets = pd.read_sql_query("SELECT t, vmatch_id, spacer_id, bmatch_id \
    FROM single_match_tripartite_networks", conTri)

pEmergeExpected = pEmergeExpected.merge(tripartiteNets,on=['t','vmatch_id']).
