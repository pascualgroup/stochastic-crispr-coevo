#!/usr/bin/env python3

import pandas as pd
import numpy as np
import scipy as sp
import math
import matplotlib.pyplot as plt
import sys
import os
import seaborn as sns
from scipy import stats
import sqlite3
# from scipy.interpolate import griddata
# from scipy.interpolate import interp1d
import matplotlib.ticker as ticker
import matplotlib.cm as cm
from matplotlib.colors import Normalize
import matplotlib.colors as mc
import colorsys

run_id = sys.argv[1]

SCRIPT_PATH = os.path.abspath(os.path.dirname(__file__)) # cluster

DBSIM_PATH = os.path.join(SCRIPT_PATH,'..','..','simulation','sweep_db_gathered.sqlite') # cluster
#home_dir = os.system("cd ~") # local
# DBSIM_PATH = os.path.join('/Volumes','Yadgah','sweep_db_gathered.sqlite') # local

# DBSIM_PATH = os.path.join('/home/armun/scratch/crispr-sweep-13-9-2021/data-analysis/match-abundances','..','..','simulation','sweep_db_gathered.sqlite')

# DBSIM_PATH = os.path.join('/Volumes','Yadgah','run_id1455_combo73_replicate15.sqlite')
# DBSIM_PATH = os.path.join('/Volumes','Yadgah','run_id1550_combo78_replicate10.sqlite')

# DBSIM_PATH = os.path.join('/Volumes','Yadgah','run_id1343_combo68_replicate3.sqlite')
# DBSIM_PATH = os.path.join('/Volumes','Yadgah','run_id1344_combo68_replicate4.sqlite')
# DBSIM_PATH = os.path.join('/Volumes','Yadgah','run_id1345_combo68_replicate5.sqlite')

DBMATCH_PATH = os.path.join(SCRIPT_PATH,'..','gathered-analyses','match-abundances','match-abundances.sqlite') # cluster
# DBMATCH_PATH = os.path.join('/Volumes','Yadgah','match-abundances.sqlite') # local
# DBMATCH_PATH = os.path.join('/Volumes','Yadgah','match-abundances_output.sqlite') # local. run_id fixed; for testing



DBCLADE_PATH = os.path.join(SCRIPT_PATH,'..','gathered-analyses','clade-abundances','clade-abundances.sqlite') # cluster
# DBCLADE_PATH = os.path.join('/Volumes','Yadgah','clade-abundances.sqlite') # local
# DBCLADE_PATH = os.path.join('/Volumes','Yadgah','clade-abundances_output.sqlite') # local. run_id fixed; for testing


# def adjust_lightness(color, amount=0.5):
#     try:
#         c = mc.cnames[color]
#     except:
#         c = color
#     c = colorsys.rgb_to_hls(*mc.to_rgb(c))
#     return colorsys.hls_to_rgb(c[0], max(0, min(1, amount * c[1])), c[2])




conSim = sqlite3.connect(DBSIM_PATH)
curSim = conSim.cursor()
print('SQLite Query: microbial abundance time series data')
microbeSim = pd.read_sql_query("SELECT t,bstrain_id,abundance FROM babundance WHERE run_id = {}".format(run_id), conSim)
microbe_stacked = microbeSim.pivot(index='t',columns='bstrain_id',values='abundance')
print('SQLite Query: viral abundance time series data')
virusSim = pd.read_sql_query("SELECT t,vstrain_id,abundance FROM vabundance WHERE run_id = {}".format(run_id), conSim)
virus_stacked = virusSim.pivot(index='t',columns='vstrain_id',values='abundance')
if microbeSim['t'][microbeSim['t'].size-1] not in virus_stacked.index:
    microbe_stacked.drop(microbeSim['t'][microbeSim['t'].size-1],inplace=True)


conClade = sqlite3.connect(DBCLADE_PATH)
curClade = conClade.cursor()
print('SQLite Query: clade data')
microbeClades = pd.read_sql_query("SELECT DISTINCT clade_id, bstrain_id \
FROM babundances WHERE run_id = {}".format(run_id), conClade)
microbeCladeIDs = pd.read_sql_query("SELECT DISTINCT clade_id \
FROM babundances WHERE run_id = {}".format(run_id), conClade)
microbeCladeAbundances = pd.read_sql_query("SELECT t, clade_id, abundance \
FROM clade_babundances WHERE run_id = {}".format(run_id), conClade)

virusClades = pd.read_sql_query("SELECT DISTINCT clade_id, vstrain_id \
FROM vabundances WHERE run_id = {}".format(run_id), conClade)
virusCladeIDs = pd.read_sql_query("SELECT DISTINCT clade_id \
FROM vabundances WHERE run_id = {}".format(run_id), conClade)
virusCladeAbundances = pd.read_sql_query("SELECT t, clade_id, abundance \
FROM clade_vabundances WHERE run_id = {}".format(run_id), conClade)



# conMatch = sqlite3.connect(DBMATCH_PATH)
# curMatch = conMatch.cursor()
# print('SQLite Query: match data')
# bmatchTypes = pd.read_sql_query("SELECT t, bstrain_match_length, babundance \
# FROM total_matched_bstrain_abundances WHERE run_id = {}".format(run_id), conMatch)
# v0abundances = pd.read_sql_query("SELECT t, bstrain_match_length, vabundance_0matches \
# FROM total_matched_bstrain_abundances WHERE run_id = {}".format(run_id), conMatch)
# v0abundances_stacked = v0abundances.pivot(index='t',columns='bstrain_match_length',
# values='vabundance_0matches')
# v0abundances_stacked.drop(0, axis=1, inplace=True)

# bmatchTypes = pd.read_sql_query("SELECT t, bstrain_match_length, babundance \
# FROM total_matched_bstrain_abundances", conMatch)
# v0abundances = pd.read_sql_query("SELECT t, bstrain_match_length, vabundance_0matches \
# FROM total_matched_bstrain_abundances", conMatch)
# v0abundances_stacked = v0abundances.pivot(index='t',columns='bstrain_match_length',
# values='vabundance_0matches')
# v0abundances_stacked.drop(0, axis=1, inplace=True)

# timeGrainSize = 8000
#
# typeDF = pd.DataFrame(data={'Match Length':[], 'Time':[], 'Microbial Abundance' :[]})
# newTime = np.linspace(0, microbeSim['t'].iloc[-1], num=timeGrainSize, endpoint=True)
# time = np.linspace(0,microbeSim['t'].iloc[-1],num=int(microbeSim['t'].iloc[-1])+1,endpoint=True)
#
#
# for matchType in np.unique(bmatchTypes["bstrain_match_length"]):
#     btypeAbundance = np.append(np.repeat(0,int(bmatchTypes[bmatchTypes['bstrain_match_length']==matchType]['t'].iloc[0])),
#     bmatchTypes[bmatchTypes['bstrain_match_length']==matchType]['babundance'])
#     btypeAbundance = np.append(btypeAbundance,np.repeat(0,int(microbeSim['t'].iloc[-1]-len(btypeAbundance)+1)))
#     bTypeInterp1 = interp1d(time, btypeAbundance)
#     bTypeInterp2 = interp1d(time, btypeAbundance, kind='cubic')
#     typeDF = typeDF.append(
#     pd.DataFrame(data={'Match Length':np.repeat(matchType,len(newTime)), 'Time': newTime, 'Microbial Abundance' : bTypeInterp2(newTime)}),
#     ignore_index=True)

conMatch = sqlite3.connect(DBMATCH_PATH)
curMatch = conMatch.cursor()

## for runs with run_id
# print('SQLite Query: match data')
# bcladematchTypesB = pd.read_sql_query(
# "SELECT t, bclade_id, match_type, matched_babundance \
# FROM bclades_match_type_abundances WHERE run_id = {}".format(run_id), conMatch)
# bcladematchTypesV = pd.read_sql_query(
# "SELECT t, bclade_id, match_type, matched_vabundance \
# FROM bclades_match_type_abundances WHERE run_id = {}".format(run_id), conMatch)
# bmatchTypes = pd.read_sql_query(
# "SELECT t, match_type, matched_babundance \
# FROM match_type_abundances WHERE run_id = {}".format(run_id), conMatch)
# bmatchTypes_stacked = bmatchTypes.pivot(index='t',columns='match_type',
# values='matched_babundance')
# vmatchTypes = pd.read_sql_query(
# "SELECT t, match_type, matched_vabundance \
# FROM match_type_abundances WHERE run_id = {}".format(run_id), conMatch)
# vmatchTypes_stacked = vmatchTypes.pivot(index='t',columns='match_type',
# values='matched_vabundance')


bcladematchTypesB = pd.read_sql_query(
"SELECT t, bclade_id, match_type, matched_babundance \
FROM bclades_match_type_abundances", conMatch)
bcladematchTypesV = pd.read_sql_query(
"SELECT t, bclade_id, match_type, matched_vabundance \
FROM bclades_match_type_abundances", conMatch)
bmatchTypes = pd.read_sql_query(
"SELECT t, match_type, matched_babundance \
FROM match_type_abundances", conMatch)
bmatchTypes_stacked = bmatchTypes.pivot(index='t',columns='match_type',
values='matched_babundance')
vmatchTypes = pd.read_sql_query(
"SELECT t, match_type, matched_vabundance \
FROM match_type_abundances", conMatch)
vmatchTypes_stacked = vmatchTypes.pivot(index='t',columns='match_type',
values='matched_vabundance')







# Designating plot path from simulation data
ID = curSim.execute('SELECT combo_id,replicate FROM runs WHERE run_id = {}'.format(run_id)).fetchall()
combo_id = ID[0][0]
replicate = ID[0][1]

#THIS PATH NEEDS TO BE LOCATED HERE BECAUSE OF COMBO_ID DEFINITION ABOVE!
PLOT_PATH = os.path.join(SCRIPT_PATH,'..','gathered-analyses','match-abundances','plots','c{}'.format(combo_id),'r{}'.format(replicate)) # cluster
# PLOT_PATH = os.path.abspath(os.path.dirname(__file__)) # local
# PLOT_PATH = os.getcwd() # local. run_id fixed; for testing


print('Compiling microbial clade abundance plots')
# Get a color map
Mcmap = cm.get_cmap('turbo')
# Get normalize function (takes data in range [vmin, vmax] -> [0, 1])
Mnorm = Normalize(vmin=1, vmax=len(microbeCladeIDs))
microbeColorDict = {}
cladeColorDict = {}
cladeDict = {}
colorInd = 0
for clade_id in microbeCladeIDs.values:
    # print(id[0])
    cladeColorDict[clade_id[0]] = colorInd
    cladeDict[clade_id[0]] = []
    colorInd += 1

for strain in microbe_stacked.columns.values:
    clade = microbeClades[microbeClades['bstrain_id']==strain]['clade_id'].values[0]
    cladeDict[clade] = np.append(cladeDict[clade],strain)

columnOrder = []
for clade_id in microbeCladeIDs.values:
    columnOrder = np.append(columnOrder,cladeDict[clade_id[0]])

columnOrder = columnOrder.astype(int)
microbe_stacked = microbe_stacked[columnOrder]

for strain in microbe_stacked.columns.values:
    clade = microbeClades[microbeClades['bstrain_id']==strain]['clade_id'].values[0]
    microbeColorDict[strain] = Mcmap(Mnorm(np.arange(1, len(microbeCladeIDs)+1, 1)))[cladeColorDict[clade]]


fig, ax = plt.subplots(1,sharex=True)
offset = 1
offsetMatches  = np.max(bcladematchTypesB['match_type']) + offset + 1
for clade_id in microbeCladeIDs.values:
    # if int(clade_id) <= 2:
    ax.scatter(bcladematchTypesB[bcladematchTypesB['bclade_id']==int(clade_id)]['t'],
    bcladematchTypesB[bcladematchTypesB['bclade_id']==int(clade_id)]['match_type'] + int((clade_id-1)*offsetMatches),
    c=bcladematchTypesB[bcladematchTypesB['bclade_id']==int(clade_id)]['matched_babundance'],
    lw=.4, s=100, cmap = 'gist_yarg',marker='|')
# new_list = range(-1, offsetMatches*len(microbeCladeIDs.values)+1)
# axes[1].yaxis.set_ticks(new_list)
ax.set_yticks(np.arange(-1, int(offsetMatches*len(microbeCladeIDs.values)), 1))
fig.canvas.draw()
matchNums = ['{}'.format(match) for match in range(0,np.max(bcladematchTypesB['match_type'])+1)]
axisSpaces = ['' for i in range(0,offset)]
matchAxis = ['']
for clade_id in microbeCladeIDs.values:
    if clade_id != microbeCladeIDs.values[-1]:
        matchAxis = np.append(matchAxis,matchNums)
        matchAxis = np.append(matchAxis,axisSpaces)
    else:
        matchAxis = np.append(matchAxis,matchNums)
        matchAxis = np.append(matchAxis,[''])

ax.set_yticklabels(matchAxis)
fig.canvas.draw()
FL = []
for clade_id in microbeCladeIDs.values:
    FL = np.append(FL,np.arange(0+int((clade_id-1)*offsetMatches),
    np.max(bcladematchTypesB['match_type'])+1+int((clade_id-1)*offsetMatches),1))

ax.yaxis.set_major_locator(plt.FixedLocator(FL))
fig.canvas.draw()
ax.xaxis.set_minor_locator(ticker.MultipleLocator(50))
for clade_id in microbeCladeIDs.values:
    ax.axhspan(-0.5+int((clade_id-1)*offsetMatches),
    np.max(bcladematchTypesB['match_type'])+.5+int((clade_id-1)*offsetMatches),
    facecolor=Mcmap(Mnorm(np.arange(1, len(microbeCladeIDs)+1, 1)))[cladeColorDict[int(clade_id)]], alpha=0.2)
fig.savefig(os.path.join('/Volumes/Yadgah','microbe-match-abundances.pdf'),dpi=1000)
plt.show()
# labels = [tick.get_text() for tick in ax.get_yticklabels()]




fig, ax = plt.subplots(2,sharex=True)
fig.suptitle('(run{0}-c{1}-r{2})'.format(run_id,combo_id,replicate))
axes = [ax[0], ax[1]]
pal = sns.color_palette("tab20b")
microbe_stacked.plot.area(ax = axes[0],stacked=True,legend=False, linewidth=0,color=microbeColorDict,sort_columns=True)
microbe_stacked.plot(stacked=True, ax=axes[0], legend=False, color='white',sort_columns=True,linewidth=.3)
axes[0].set_ylabel(ylabel ='Microbial Strain Abundances',labelpad=15,fontsize=7)
axes[0].set_xlabel(xlabel = 'Time t',fontsize=7)
axes[0].ticklabel_format(axis = 'y',style='sci',scilimits=(0,0))
axes[0].xaxis.set_minor_locator(ticker.MultipleLocator(50))
lim = axes[0].get_ylim()
axes[0].set_ylim(0,lim[1])
offset = 1
offsetMatches  = np.max(bcladematchTypesB['match_type']) + offset + 1
for clade_id in microbeCladeIDs.values:
    axes[1].scatter(bcladematchTypesB[bcladematchTypesB['bclade_id']==int(clade_id)]['t'],
    bcladematchTypesB[bcladematchTypesB['bclade_id']==int(clade_id)]['match_type'] + int((clade_id-1)*offsetMatches),
    c=bcladematchTypesB[bcladematchTypesB['bclade_id']==int(clade_id)]['matched_babundance'],
    lw=.75, s=10, cmap = 'gist_yarg',marker='|')
# new_list = range(-1, offsetMatches*len(microbeCladeIDs.values)+1)
# axes[1].yaxis.set_ticks(new_list)
axes[1].set_yticks(np.arange(-1, int(offsetMatches*len(microbeCladeIDs.values)), 1))
fig.canvas.draw()
matchNums = ['{}'.format(match) for match in range(0,np.max(bcladematchTypesB['match_type'])+1)]
axisSpaces = ['' for i in range(0,offset)]
matchAxis = ['']
for clade_id in microbeCladeIDs.values:
    if clade_id != microbeCladeIDs.values[-1]:
        matchAxis = np.append(matchAxis,matchNums)
        matchAxis = np.append(matchAxis,axisSpaces)
    else:
        matchAxis = np.append(matchAxis,matchNums)
        matchAxis = np.append(matchAxis,[''])

axes[1].set_yticklabels(matchAxis)
fig.canvas.draw()
FL = []
for clade_id in microbeCladeIDs.values:
    FL = np.append(FL,np.arange(0+int((clade_id-1)*offsetMatches),
    np.max(bcladematchTypesB['match_type'])+1+int((clade_id-1)*offsetMatches),1))

axes[1].yaxis.set_major_locator(plt.FixedLocator(FL))
fig.canvas.draw()
ax.xaxis.set_minor_locator(ticker.MultipleLocator(50))
for clade_id in microbeCladeIDs.values:
    axes[1].axhspan(-0.5+int((clade_id-1)*offsetMatches),
    np.max(bcladematchTypesB['match_type'])+.5+int((clade_id-1)*offsetMatches),
    facecolor=Mcmap(Mnorm(np.arange(1, len(microbeCladeIDs)+1, 1)))[cladeColorDict[int(clade_id)]], alpha=0.2)
size = fig.get_size_inches()
fig.set_size_inches(size[0]*2, size[1]*2)
fig.savefig(os.path.join('/Volumes/Yadgah','microbe-match-abundances-stacked.pdf'),dpi=1000)












# fig, ax = plt.subplots(1,sharex=True)
# ax.scatter(bcladematchTypesB[bcladematchTypesB['bclade_id']==clade_id]['t'],
# bcladematchTypesB[bcladematchTypesB['bclade_id']==clade_id]['match_type'] + offsetMatches,
# c=bcladematchTypesB[bcladematchTypesB['bclade_id']==clade_id]['matched_babundance'],
# lw=.25, s=200, cmap = 'RdPu',marker='|')
# plt.show()
#
# # new_list = range(-1, np.max(bmatchTypes["bstrain_match_length"])+2)
# axes[1].yaxis.set_ticks(new_list)
# # x_formatter = ticker.ScalarFormatter(useOffset=False)
# # axes[1].xaxis.set_major_formatter(x_formatter)
# axes[1].ticklabel_format(useOffset = False)
#
# fig.canvas.draw()
# fig.savefig(os.path.join('/Volumes/Yadgah','microbe-match-abundances.png'),dpi=1000)
# plt.show()










# fig, ax = plt.subplots(2,sharex=True)
# fig.suptitle('(run{0}-c{1}-r{2})'.format(run_id,combo_id,replicate))
# axes = [ax[0], ax[1]]
# pal = sns.color_palette("tab20b")
# microbe_stacked.plot.area(ax = axes[0],stacked=True,legend=False, linewidth=0,color=microbeColorDict,sort_columns=True)
# microbe_stacked.plot(stacked=True, ax=axes[0], legend=False, color='black',sort_columns=True,linewidth=.1)
# axes[0].set_ylabel(ylabel ='Microbial Strain Abundances',labelpad=15,fontsize=7)
# axes[0].set_xlabel(xlabel = 'Time t',fontsize=7)
# axes[0].ticklabel_format(axis = 'y',style='sci',scilimits=(0,0))
# axes[0].xaxis.set_minor_locator(ticker.MultipleLocator(50))
# lim = axes[0].get_ylim()
# axes[0].set_ylim(0,lim[1])
# # typeHeat = axes[1].scatter(typeDF['Time'], typeDF['Match Length'], c=typeDF['Microbial Abundance'],
# # lw=.2, s=200, cmap = 'RdPu',marker='|')
# typeHeat = axes[1].scatter(bmatchTypes['t'], bmatchTypes['bstrain_match_length'],
# c=bmatchTypes['babundance'],lw=.25, s=200, cmap = 'RdPu',marker='|')
# new_list = range(-1, np.max(bmatchTypes["bstrain_match_length"])+2)
# axes[1].yaxis.set_ticks(new_list)
# # x_formatter = ticker.ScalarFormatter(useOffset=False)
# # axes[1].xaxis.set_major_formatter(x_formatter)
# axes[1].ticklabel_format(useOffset = False)
# #cb = fig.colorbar(typeHeat, ax=axes[1], shrink = .1, location='bottom',pad=.9)
# # cb.ax.yaxis.set_offset_position('right')
# # offsety = cb.ax.yaxis.get_offset_text().get_text()
# # cb.ax.tick_params(labelsize=4)
# # cb.ax.set_xlabel(r'Microbial Abundance $/$ {}'.format(offsety),labelpad=.1,fontsize=4)
# # # cb.ax.yaxis.offsetText.set_visible(False)
# # cb.update_ticks()
# # v0abundances_stacked.plot(ax =axes[2], color=sns.color_palette("Accent"),linewidth=0.75)
# # axes[2].set_yscale('log')
# fig.canvas.draw()
# fig.savefig(os.path.join('/Volumes/Yadgah','microbe-match-abundances.png'),dpi=1000)
# # plt.show()
