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

run_id = sys.argv[1]

SCRIPT_PATH = os.path.abspath(os.path.dirname(__file__)) # cluster

DBSIM_PATH = os.path.join(SCRIPT_PATH,'..','..','simulation','sweep_db_gathered.sqlite') # cluster
#home_dir = os.system("cd ~") # local
#DBSIM_PATH = os.path.join('/Volumes','Yadgah','sweep_db_gathered.sqlite') # local

DBRICH_PATH = os.path.join(SCRIPT_PATH,'..','gathered-analyses','richness','richness.sqlite') # cluster
#DBRICH_PATH = os.path.join('/Volumes','Yadgah','richness.sqlite') # local

DBSHAN_PATH = os.path.join(SCRIPT_PATH,'..','gathered-analyses','shannon','shannon.sqlite') # cluster
#DBSHAN_PATH = os.path.join('/Volumes','Yadgah','shannon.sqlite') # local

DBPEAKS_PATH = os.path.join(SCRIPT_PATH,'..','gathered-analyses','walls','walls.sqlite') # cluster
#DBPEAKS_PATH = os.path.join('/Volumes','Yadgah','walls.sqlite') # local
#DBPEAKS_PATH = os.path.join('/Volumes','Yadgah','walls_output.sqlite') # local. run_id fixed; for testing

DBMATCH_PATH = os.path.join(SCRIPT_PATH,'..','gathered-analyses','walls','match_diversity.sqlite') # cluster
#DBMATCH_PATH = os.path.join('/Volumes','Yadgah','match_diversity.sqlite') # local
#DBMATCH_PATH = os.path.join('/Volumes','Yadgah','match_diversity_output.sqlite') # local. run_id fixed; for testing

conSim = sqlite3.connect(DBSIM_PATH)
curSim = conSim.cursor()
print('SQLite Query: microbial abundance time series data')
microbeSim = pd.read_sql_query("SELECT t,bstrain_id,abundance FROM babundance WHERE run_id = {}".format(run_id), conSim)
microbe_stacked = microbeSim.pivot(index='t',columns='bstrain_id',values='abundance')
print('SQLite Query: viral abundance time series data')
virusSim = pd.read_sql_query("SELECT t,vstrain_id,abundance FROM vabundance WHERE run_id = {}".format(run_id), conSim)
virus_stacked = virusSim.pivot(index='t',columns='vstrain_id',values='abundance')

conMatch = sqlite3.connect(DBMATCH_PATH)
curMatch = conRich.cursor()
print('SQLite Query: virus shannon data')
bstrainMatches = pd.read_sql_query("SELECT t, bstrain_id, bstrain_frequency, bstrain_abundance\\
FROM strain_matches WHERE run_id = {}".format(run_id), conMatch)
vstrainMatches = pd.read_sql_query("SELECT t, vstrain_id, vstrain_frequency, vstrain_abundance\\
FROM strain_matches WHERE run_id = {}".format(run_id), conMatch)
bmatchTypes = pd.read_sql_query("SELECT t, match_type, matched_bfrequency, matched_babundance\\
FROM match_abundances WHERE run_id = {}".format(run_id), conMatch)
vmatchTypes = pd.read_sql_query("SELECT t, match_type, matched_vfrequency, matched_vabundance\\
FROM match_abundances WHERE run_id = {}".format(run_id), conMatch)
bspacers = pd.read_sql_query("SELECT t, spacer_id, bfrequency_with_spacer_id, babundance_with_spacer_id\\
FROM spacer_match_frequencies WHERE run_id = {}".format(run_id), conMatch)
vspacers = pd.read_sql_query("SELECT t, spacer_id, vfrequency_with_spacer_id, vabundance_with_spacer_id\\
FROM spacer_match_frequencies WHERE run_id = {}".format(run_id), conMatch)



# Designating plot path from simulation data
ID = curSim.execute('SELECT combo_id,replicate FROM runs WHERE run_id = {}'.format(run_id)).fetchall()
combo_id = ID[0][0]
replicate = ID[0][1]

#THIS PATH NEEDS TO BE LOCATED HERE BECAUSE OF COMBO_ID DEFINITION ABOVE!
PLOT_PATH = os.path.join(SCRIPT_PATH,'..','gathered-analyses','match-diversity','plots','c{}'.format(combo_id),'r{}'.format(replicate)) # cluster
#PLOT_PATH = os.path.abspath(os.path.dirname(__file__)) # local
#PLOT_PATH = os.getcwd() # local. run_id fixed; for testing


print('Compiling microbial peak and shannon diversity plot')
fig, ax = plt.subplots(2,sharex=True)
#fig = plt.figure()
fig.suptitle('(run{0}-c{1}-r{2})'.format(run_id,combo_id,replicate))
axes = [ax[0], ax[0].twinx(), ax[1], ax[1].twinx()]
microbe_stacked.plot.area(ax = axes[0],stacked=True,legend=False, linewidth=0,color=pal)
axes[0].set_ylabel(ylabel ='Microbial Immune Abundances N_i',labelpad=15,fontsize=7)
axes[0].set_xlabel(xlabel = 'Time t',fontsize=7)
axes[0].ticklabel_format(style='sci',scilimits=(0,0))
lim = axes[0].get_ylim()
axes[0].set_ylim(0,lim[1])

microbeShannon.plot(x='t',ax = axes[2],legend=False,color=pal[0],linewidth=0.75)
axes[2].set_ylabel(ylabel ='Microbe Shannon Diversity',labelpad=15,fontsize=7)
axes[2].set_xlabel(xlabel = 'Time t',fontsize=7)
axes[2].ticklabel_format(style='sci',scilimits=(0,0))
axes[3].set_yticks([])
axes[3].fill_between(microbePeaks["t"],0,1, where=microbePeaks["peak_presence"]>0,interpolate=True, color="grey", alpha=0.5, linewidth=0,transform=axes[2].get_xaxis_transform())
lim = axes[2].get_ylim()
axes[2].set_ylim(0,lim[1])
lim = axes[3].get_ylim()
axes[3].set_ylim(0,lim[1])
fig.tight_layout()
fig.savefig(os.path.join(PLOT_PATH,'microbe-peaks-shannon.png'),dpi=500)


fig, ax = plt.subplots(2,sharex=True)
#fig = plt.figure()
fig.suptitle('(run{0}-c{1}-r{2})'.format(run_id,combo_id,replicate))
axes = [ax[0], ax[0].twinx(), ax[1], ax[1].twinx()]
virus_stacked.plot.area(ax = axes[0],stacked=True,legend=False, linewidth=0,color=pal)
axes[0].set_ylabel(ylabel ='Viral Strain Abundances V_i',labelpad=15,fontsize=7)
axes[0].set_xlabel(xlabel = 'Time t',fontsize=7)
axes[0].ticklabel_format(style='sci',scilimits=(0,0))
lim = axes[0].get_ylim()
axes[0].set_ylim(0,lim[1])

microbeShannon.plot(x='t',ax = axes[2],legend=False,color=pal[0],linewidth=0.75)
axes[2].set_ylabel(ylabel ='Microbe Shannon Diversity',labelpad=15,fontsize=7)
axes[2].set_xlabel(xlabel = 'Time t',fontsize=7)
axes[2].ticklabel_format(style='sci',scilimits=(0,0))
axes[3].set_yticks([])
axes[3].fill_between(microbePeaks["t"],0,1, where=microbePeaks["peak_presence"]>0,interpolate=True, color="grey", alpha=0.5, linewidth=0,transform=axes[2].get_xaxis_transform())
lim = axes[2].get_ylim()
axes[2].set_ylim(0,lim[1])
lim = axes[3].get_ylim()
axes[3].set_ylim(0,lim[1])
fig.tight_layout()
fig.savefig(os.path.join(PLOT_PATH,'microbe-peaks-shannon.png'),dpi=500)
