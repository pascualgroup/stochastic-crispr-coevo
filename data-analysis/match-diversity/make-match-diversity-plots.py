#!/usr/bin/env python3

import pandas as pd
import numpy as np
import scipy as sp
import math
import matplotlib.pyplot as plt
import matplotlib.use
matplotlib.use('Agg')
import sys
import os
import seaborn as sns
from scipy import stats
import sqlite3
from scipy.interpolate import griddata
from scipy.interpolate import interp1d
import matplotlib.ticker as ticker

run_id = sys.argv[1]

SCRIPT_PATH = os.path.abspath(os.path.dirname(__file__)) # cluster

#DBSIM_PATH = os.path.join(SCRIPT_PATH,'..','..','simulation','sweep_db_gathered.sqlite') # cluster
#home_dir = os.system("cd ~") # local
#DBSIM_PATH = os.path.join('/Volumes','Yadgah','sweep_db_gathered.sqlite') # local

DBSIM_PATH = os.path.join('/Volumes','Yadgah','run_id1455_combo73_replicate15.sqlite')

#DBRICH_PATH = os.path.join(SCRIPT_PATH,'..','gathered-analyses','richness','richness.sqlite') # cluster
#DBRICH_PATH = os.path.join('/Volumes','Yadgah','richness.sqlite') # local

#DBSHAN_PATH = os.path.join(SCRIPT_PATH,'..','gathered-analyses','shannon','shannon.sqlite') # cluster
#DBSHAN_PATH = os.path.join('/Volumes','Yadgah','shannon.sqlite') # local

#DBPEAKS_PATH = os.path.join(SCRIPT_PATH,'..','gathered-analyses','walls','walls.sqlite') # cluster
#DBPEAKS_PATH = os.path.join('/Volumes','Yadgah','walls.sqlite') # local
#DBPEAKS_PATH = os.path.join('/Volumes','Yadgah','walls_output.sqlite') # local. run_id fixed; for testing

#DBMATCH_PATH = os.path.join(SCRIPT_PATH,'..','gathered-analyses','walls','match_diversity.sqlite') # cluster
#DBMATCH_PATH = os.path.join('/Volumes','Yadgah','match-diversity.sqlite') # local
DBMATCH_PATH = os.path.join('/Volumes','Yadgah','match-diversity_output.sqlite') # local. run_id fixed; for testing

conSim = sqlite3.connect(DBSIM_PATH)
curSim = conSim.cursor()
print('SQLite Query: microbial abundance time series data')
microbeSim = pd.read_sql_query("SELECT t,bstrain_id,abundance FROM babundance WHERE run_id = {}".format(run_id), conSim)
microbe_stacked = microbeSim.pivot(index='t',columns='bstrain_id',values='abundance')
print('SQLite Query: viral abundance time series data')
virusSim = pd.read_sql_query("SELECT t,vstrain_id,abundance FROM vabundance WHERE run_id = {}".format(run_id), conSim)
virus_stacked = virusSim.pivot(index='t',columns='vstrain_id',values='abundance')

conMatch = sqlite3.connect(DBMATCH_PATH)
curMatch = conMatch.cursor()
print('SQLite Query: virus shannon data')
bstrainMatches = pd.read_sql_query("SELECT t, bstrain_id, bstrain_frequency, bstrain_abundance \
FROM strain_matches WHERE run_id = {}".format(run_id), conMatch)
vstrainMatches = pd.read_sql_query("SELECT t, vstrain_id, vstrain_frequency, vstrain_abundance \
FROM strain_matches WHERE run_id = {}".format(run_id), conMatch)
bmatchTypes = pd.read_sql_query("SELECT t, match_type, matched_bfrequency, matched_babundance \
FROM match_abundances WHERE run_id = {}".format(run_id), conMatch)
vmatchTypes = pd.read_sql_query("SELECT t, match_type, matched_vfrequency, matched_vabundance \
FROM match_abundances WHERE run_id = {}".format(run_id), conMatch)
bspacers = pd.read_sql_query("SELECT t, spacer_id, bfrequency_with_spacer_id, babundance_with_spacer_id \
FROM spacer_match_frequencies WHERE run_id = {}".format(run_id), conMatch)
vspacers = pd.read_sql_query("SELECT t, spacer_id, vfrequency_with_spacer_id, vabundance_with_spacer_id \
FROM spacer_match_frequencies WHERE run_id = {}".format(run_id), conMatch)
pdi = pd.read_sql_query("SELECT t, PDI, maxPDI, PDI_degree1, PDI_degree2, PDI_degree3, PDI_degree4 \
FROM PDI WHERE run_id = {}".format(run_id), conMatch)


# Designating plot path from simulation data
ID = curSim.execute('SELECT combo_id,replicate FROM runs WHERE run_id = {}'.format(run_id)).fetchall()
combo_id = ID[0][0]
replicate = ID[0][1]

#THIS PATH NEEDS TO BE LOCATED HERE BECAUSE OF COMBO_ID DEFINITION ABOVE!
PLOT_PATH = os.path.join(SCRIPT_PATH,'..','gathered-analyses','match-diversity','plots','c{}'.format(combo_id),'r{}'.format(replicate)) # cluster
#PLOT_PATH = os.path.abspath(os.path.dirname(__file__)) # local
#PLOT_PATH = os.getcwd() # local. run_id fixed; for testing


print('Compiling microbial abundance and match diversity plots')




grid_kws = {"height_ratios": (.9, .9, .05), "hspace": (.1, .3)}

timeGrainSize = 8000

typeDF = pd.DataFrame(data={'Match Type':[], 'Time':[], 'Microbial Abundance' :[]})
newTime = np.linspace(0, microbeSim['t'].iloc[-1], num=timeGrainSize, endpoint=True)
time = np.linspace(0,microbeSim['t'].iloc[-1],num=int(microbeSim['t'].iloc[-1])+1,endpoint=True)


for matchType in np.unique(bmatchTypes["match_type"]):
    # time = np.append(np.linspace(0,bmatchTypes[bmatchTypes['match_type']==matchType]['t'].iloc[0],
    # num = int(bmatchTypes[bmatchTypes['match_type']==matchType]['t'].iloc[0]),endpoint=False),
    # bmatchTypes[bmatchTypes['match_type']==matchType]['t'])
    # time = np.append(time,np.linspace(bmatchTypes[bmatchTypes['match_type']==matchType]['t'].iloc[-1]+1,microbeSim['t'].iloc[-1],
    # num=int(microbeSim['t'].iloc[-1]-bmatchTypes[bmatchTypes['match_type']==matchType]['t'].iloc[-1])))
    btypeAbundance = np.append(np.repeat(0,int(bmatchTypes[bmatchTypes['match_type']==matchType]['t'].iloc[0])),
    bmatchTypes[bmatchTypes['match_type']==matchType]['matched_babundance'])
    btypeAbundance = np.append(btypeAbundance,np.repeat(0,int(microbeSim['t'].iloc[-1]-len(btypeAbundance)+1)))
    bTypeInterp1 = interp1d(time, btypeAbundance)
    bTypeInterp2 = interp1d(time, btypeAbundance, kind='cubic')
    typeDF = typeDF.append(
    pd.DataFrame(data={'Match Type':np.repeat(matchType,len(newTime)), 'Time': newTime, 'Microbial Abundance' : bTypeInterp2(newTime)}),
    ignore_index=True)


grid_kws = {"height_ratios": (.9, .9, .05), "hspace": (.1, .3)}

fig, ax = plt.subplots(2,sharex=True)
fig.suptitle('(run{0}-c{1}-r{2})'.format(run_id,combo_id,replicate))
axes = [ax[0], ax[1]]
pal = sns.color_palette("tab20b")
numColors = len(pal)
microbe_stacked.plot.area(ax = axes[0],stacked=True,legend=False, linewidth=0,color=pal)
axes[0].set_ylabel(ylabel ='Microbial Immune Abundances N_i',labelpad=15,fontsize=7)
axes[0].set_xlabel(xlabel = 'Time t',fontsize=7)
axes[0].ticklabel_format(style='sci',scilimits=(0,0))
lim = axes[0].get_ylim()
axes[0].set_ylim(0,lim[1])
typeHeat = axes[1].scatter(typeDF['Time'], typeDF['Match Type'], c=typeDF['Microbial Abundance'],
lw=.2, s=200, cmap = 'RdPu',marker='|')
new_list = range(0, np.max(bmatchTypes["match_type"])+2)
axes[1].yaxis.set_ticks(new_list)
x_formatter = matplotlib.ticker.ScalarFormatter(useOffset=False)
axes[1].xaxis.set_major_formatter(x_formatter)
axes[1].ticklabel_format(useOffset = False)
cb = fig.colorbar(typeHeat, ax=axes[1], shrink = 0.4, location='bottom')
cb.formatter.set_powerlimits((0, 0))
cb.ax.yaxis.set_offset_position('right')
cb.update_ticks()
fig.canvas.draw()
# typeDF = typeDF.pivot("Match Type", "Time", "Microbial Abundance")
# axes[1] = sns.heatmap(typeDF)
fig.savefig(os.path.join('/Volumes/Yadgah','microbe-match-diversity-abundances.png'),dpi=500)

fig, ax = plt.subplots(4,sharex=True, gridspec_kw=grid_kws)
fig.suptitle('(run{0}-c{1}-r{2})'.format(run_id,combo_id,replicate))
axes = [ax[0], ax[0].twinx()]
microbe_stacked.plot.area(ax = axes[0],stacked=True,legend=False, linewidth=0,color=pal)
axes[0].set_ylabel(ylabel ='Microbial Immune Abundances N_i',labelpad=15,fontsize=7)
axes[0].set_xlabel(xlabel = 'Time t',fontsize=7)
axes[0].ticklabel_format(style='sci',scilimits=(0,0))
lim = axes[0].get_ylim()
axes[0].set_ylim(0,lim[1])
bmatchTypes = sns.heatmap(typeDF, ax=ax)
axes[1].set_ylabel(ylabel ='Match Type',labelpad=15,fontsize=7)
axes[1].set_xlabel(xlabel = 'Time t',fontsize=7)
axes[1].ticklabel_format(style='sci',scilimits=(0,0))
axes[2].set_ylim(0,lim[1])
lim = axes[3].get_ylim()
fig.tight_layout()
fig.savefig(os.path.join(PLOT_PATH,'microbe-match-diversity-abundances.png'),dpi=500)




fig, ax = plt.subplots(3,sharex=True)
fig.suptitle('(run{0}-c{1}-r{2})'.format(run_id,combo_id,replicate))
axes = [ax[0], ax[1], ax[2]]
pal = sns.color_palette("tab20b")
microbe_stacked.plot.area(ax = axes[0],stacked=True,legend=False, linewidth=0,color=pal)
axes[0].set_ylabel(ylabel ='Microbial Immune Abundances N_i',labelpad=15,fontsize=7)
axes[0].set_xlabel(xlabel = 'Time t',fontsize=7)
axes[0].ticklabel_format(style='sci',scilimits=(0,0))
lim = axes[0].get_ylim()
axes[0].set_ylim(0,lim[1])
typeHeat = axes[1].scatter(typeDF['Time'], typeDF['Match Type'], c=typeDF['Microbial Abundance'],
lw=.2, s=200, cmap = 'RdPu',marker='|')
new_list = range(0, np.max(bmatchTypes["match_type"])+2)
axes[1].yaxis.set_ticks(new_list)
x_formatter = matplotlib.ticker.ScalarFormatter(useOffset=False)
axes[1].xaxis.set_major_formatter(x_formatter)
axes[1].ticklabel_format(useOffset = False)
cb = fig.colorbar(typeHeat, ax=axes[1], shrink = 0.2, location='bottom',pad=0.1)
cb.formatter.set_powerlimits((0, 0))
cb.ax.yaxis.set_offset_position('right')
offsety = cb.ax.yaxis.get_offset_text().get_text()
cb.ax.tick_params(labelsize=4)
cb.ax.set_xlabel(r'Microbial Abundance $/$ {}'.format(offsety),labelpad=.1,fontsize=4)
cb.ax.yaxis.offsetText.set_visible(False)
cb.update_ticks()
axes[2].plot(pdi['t'], pdi['PDI'],color=pal[0],linewidth=0.75)
fig.canvas.draw()
# typeDF = typeDF.pivot("Match Type", "Time", "Microbial Abundance")
# axes[1] = sns.heatmap(typeDF)
fig.savefig(os.path.join('/Volumes/Yadgah','microbe-match-diversity-abundances.png'),dpi=1000)
