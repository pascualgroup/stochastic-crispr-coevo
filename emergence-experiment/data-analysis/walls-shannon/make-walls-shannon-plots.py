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
if len(sys.argv) == 2:
    threshold = 0
else:
    threshold = float(sys.argv[2])/100
resolve = 500

SCRIPT_PATH = os.path.abspath(os.path.dirname(__file__)) # cluster
DBSIM_PATH = os.path.join(SCRIPT_PATH,'..','..','simulation','sweep_db_gathered.sqlite') # cluster
#DBSIM_PATH = os.path.join('/Volumes','Yadgah','sweep_db_gathered.sqlite') # local
DBRICH_PATH = os.path.join(SCRIPT_PATH,'..','gathered-analyses','richness','richness.sqlite') # cluster
#DBRICH_PATH = os.path.join('/Volumes','Yadgah','richness.sqlite') # local
DBSHAN_PATH = os.path.join(SCRIPT_PATH,'..','gathered-analyses','shannon','shannon.sqlite') # cluster
#DBSHAN_PATH = os.path.join('/Volumes','Yadgah','shannon.sqlite') # local
DBPEAKS_PATH = os.path.join(SCRIPT_PATH,'..','gathered-analyses','walls-shannon','walls-shannon.sqlite') # cluster
#DBPEAKS_PATH = os.path.join('/Volumes','Yadgah','walls-shannon.sqlite') # local
#DBPEAKS_PATH = os.path.join('/Volumes','Yadgah','walls-shannon_output.sqlite') # local. run_id fixed; for testing

conSim = sqlite3.connect(DBSIM_PATH)
curSim = conSim.cursor()
ID = curSim.execute('SELECT combo_id,replicate FROM runs WHERE run_id = {}'.format(run_id)).fetchall()
combo_id = ID[0][0]
replicate = ID[0][1]
conShan = sqlite3.connect(DBSHAN_PATH)
curShan = conShan.cursor()
conPeaks = sqlite3.connect(DBPEAKS_PATH)
curPeaks = conPeaks.cursor()

PLOT_PATH = os.path.join(SCRIPT_PATH,'..','gathered-analyses','walls-shannon','plots','c{}'.format(combo_id),'r{}'.format(replicate)) # cluster
#PLOT_PATH = os.path.abspath(os.path.dirname(__file__)) # local
#PLOT_PATH = os.getcwd() # local. run_id fixed; for testing

print('SQLite Query: microbial abundance time series data')
microbe_stacked = pd.read_sql_query("SELECT t,bstrain_id,abundance \
    FROM babundance WHERE run_id = {}".format(run_id), conSim)
microbeThreshAbund = pd.read_sql_query("SELECT t, microbial_abundance \
    FROM summary WHERE run_id = {}".format(run_id), conSim)
microbeThreshAbund['microbial_abundance'] = \
    threshold*microbeThreshAbund['microbial_abundance']
microbeThreshAbund = \
    [microbeThreshAbund[microbeThreshAbund.t == t].microbial_abundance.values[0]
    for t in microbe_stacked.t]
microbe_stacked = microbe_stacked[microbe_stacked.abundance >= microbeThreshAbund]
microbe_stacked = microbe_stacked.pivot(index='t',columns='bstrain_id',values='abundance')

print('SQLite Query: viral abundance time series data')
virus_stacked = pd.read_sql_query("SELECT t,vstrain_id,abundance \
    FROM vabundance WHERE run_id = {}".format(run_id), conSim)
virusThreshAbund = pd.read_sql_query("SELECT t, viral_abundance \
    FROM summary WHERE run_id = {}".format(run_id), conSim)
virusThreshAbund['viral_abundance'] = \
    threshold*virusThreshAbund['viral_abundance']
virusThreshAbund = \
    [virusThreshAbund[virusThreshAbund.t == t].viral_abundance.values[0]
    for t in virus_stacked.t]
virus_stacked = virus_stacked[virus_stacked.abundance >= virusThreshAbund]
virus_stacked = virus_stacked.pivot(index='t',columns='vstrain_id',values='abundance')

if len(microbe_stacked.index) > len(virus_stacked.index):
    microbe_stacked.drop(index=list(set(microbe_stacked.index)-set(virus_stacked.index)),inplace=True)

# conRich = sqlite3.connect(DBRICH_PATH)
# curRich = conRich.cursor()
# print('SQLite Query: virus strain richness data')
# virusRichness = pd.read_sql_query("SELECT t, vrichness FROM richness WHERE run_id = {}".format(run_id), conRich)
# print('SQLite Query: microbe immune richness data')
# microbeRichness = pd.read_sql_query("SELECT t,brichness FROM richness WHERE run_id = {}".format(run_id), conRich)


print('SQLite Query: virus shannon data')
virusShannon = pd.read_sql_query("SELECT t, vshannon FROM shannon_diversity WHERE run_id = {}".format(run_id), conShan)
print('SQLite Query: microbe shannon data')
microbeShannon = pd.read_sql_query("SELECT t,bshannon FROM shannon_diversity WHERE run_id = {}".format(run_id), conShan)



thresholds = curPeaks.execute('SELECT upper_threshold,lower_threshold,shannon_threshold FROM threshold_values WHERE run_id = {}'.format(run_id)).fetchall()
upperThreshold = thresholds[0][0]
lowerThreshold = thresholds[0][1]
shannonThreshold = thresholds[0][2]

print('SQLite Query: microbial walls data')
microbePeaks = pd.read_sql_query("SELECT t,peak_presence FROM microbial_peak_series WHERE run_id = {}".format(run_id), conPeaks)
microbeWalls = pd.read_sql_query("SELECT t,wall_presence FROM microbial_wall_series WHERE run_id = {}".format(run_id), conPeaks)
numPeaks = pd.read_sql_query("SELECT num_peaks,num_walls FROM microbial_peakwall_count", conPeaks)
durations = pd.read_sql_query("SELECT begin_row,end_row,begin_t,end_t,peak_number,wall_number,duration FROM microbial_peakwall_durations", conPeaks)

pal = sns.color_palette("tab20b")

print('Compiling microbial peak plot')
fig, ax = plt.subplots()
fig.suptitle('Microbial Immune Abundances (run{0}-c{1}-r{2})'.format(run_id,combo_id,replicate))
axes = [ax, ax.twinx(), ax.twiny()]
microbe_stacked.plot.area(ax = axes[0],stacked=True,legend=False, linewidth=0,color=pal)
axes[0].set_ylabel(ylabel ='Abundances N_i',labelpad=15,fontsize=7)
axes[0].set_xlabel(xlabel = 'Time t',fontsize=7)
axes[0].ticklabel_format(style='sci',scilimits=(0,0))
axes[2].set_xticks([])
axes[2].axhline(upperThreshold, color='k', lw=.5, linestyle='dotted')
axes[2].axhline(lowerThreshold, color='k', lw=.5, linestyle='dotted')
axes[1].set_yticks([])
axes[1].fill_between(microbePeaks["t"],0,1, where=microbePeaks["peak_presence"]>0,color="grey", alpha=0.5, linewidth=0,transform=axes[0].get_xaxis_transform())
lim = axes[0].get_ylim()
axes[0].set_ylim(0,lim[1])
lim = axes[1].get_ylim()
axes[1].set_ylim(0,lim[1])
lim = axes[2].get_ylim()
axes[2].set_ylim(0,lim[1])
fig.tight_layout()
fig.savefig(os.path.join(PLOT_PATH,'microbe-peaks.png'),dpi=resolve)


print('Compiling microbial wall plot')
fig, ax = plt.subplots()
fig.suptitle('Microbial Immune Abundances (run{0}-c{1}-r{2})'.format(run_id,combo_id,replicate))
axes = [ax, ax.twinx(), ax.twiny()]
microbe_stacked.plot.area(ax = axes[0],stacked=True,legend=False, linewidth=0,color=pal)
axes[0].set_ylabel(ylabel ='Abundances N_i',labelpad=15,fontsize=7)
axes[0].set_xlabel(xlabel = 'Time t',fontsize=7)
axes[0].ticklabel_format(style='sci',scilimits=(0,0))
axes[2].set_xticks([])
axes[2].axhline(upperThreshold, color='k', lw=.5, linestyle='dotted')
axes[2].axhline(lowerThreshold, color='k', lw=.5, linestyle='dotted')
axes[1].set_yticks([])
axes[1].fill_between(microbeWalls["t"],0,1, where=microbeWalls["wall_presence"]>0,color="grey", alpha=0.5, linewidth=0,transform=axes[0].get_xaxis_transform())
lim = axes[0].get_ylim()
axes[0].set_ylim(0,lim[1])
lim = axes[1].get_ylim()
axes[1].set_ylim(0,lim[1])
lim = axes[2].get_ylim()
axes[2].set_ylim(0,lim[1])
fig.tight_layout()
fig.savefig(os.path.join(PLOT_PATH,'microbe-walls.png'),dpi=resolve)


print('Compiling microbial peak and shannon diversity plot')
fig, ax = plt.subplots(2,sharex=True)
fig.suptitle('(run{0}-c{1}-r{2})'.format(run_id,combo_id,replicate))
axes = [ax[0], ax[0].twinx(), ax[1], ax[1].twinx(), ax[0].twiny()]
microbe_stacked.plot.area(ax = axes[0],stacked=True,legend=False, linewidth=0,color=pal)
axes[0].set_ylabel(ylabel ='Microbial Immune Abundances N_i',labelpad=15,fontsize=7)
axes[0].set_xlabel(xlabel = 'Time t',fontsize=7)
axes[0].ticklabel_format(style='sci',scilimits=(0,0))
axes[4].set_xticks([])
axes[4].axhline(upperThreshold, color='k', lw=.5, linestyle='dotted')
axes[4].axhline(lowerThreshold, color='k', lw=.5, linestyle='dotted')
axes[4].axhline(0, color='k', lw=1)
axes[1].set_yticks([])
axes[1].fill_between(microbePeaks["t"],0,1, where=microbePeaks["peak_presence"]>0,interpolate=True, color="grey", alpha=0.5, linewidth=0,transform=axes[0].get_xaxis_transform())
lim = axes[0].get_ylim()
axes[0].set_ylim(0,lim[1])
lim = axes[1].get_ylim()
axes[1].set_ylim(0,lim[1])
lim = axes[4].get_ylim()
axes[4].set_ylim(0,lim[1])
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
fig.savefig(os.path.join(PLOT_PATH,'microbe-peaks-shannon.png'),dpi=resolve)



print('Compiling microbial walls and shannon diversity plot')
fig, ax = plt.subplots(2,sharex=True)
fig.suptitle('(run{0}-c{1}-r{2})'.format(run_id,combo_id,replicate))
axes = [ax[0], ax[0].twinx(), ax[1], ax[1].twinx(), ax[0].twiny(), ax[1].twiny()]
microbe_stacked.plot.area(ax = axes[0],stacked=True,legend=False, linewidth=0,color=pal)
axes[0].set_ylabel(ylabel ='Microbial Immune Abundances N_i',labelpad=15,fontsize=7)
axes[0].set_xlabel(xlabel = 'Time t',fontsize=7)
axes[0].ticklabel_format(style='sci',scilimits=(0,0))
axes[4].set_xticks([])
axes[4].axhline(upperThreshold, color='k', lw=.5, linestyle='dotted')
axes[4].axhline(lowerThreshold, color='k', lw=.5, linestyle='dotted')
axes[4].axhline(0, color='k', lw=1)
axes[1].set_yticks([])
axes[1].fill_between(microbeWalls["t"],0,1, where=microbeWalls["wall_presence"]>0,interpolate=True, color="grey", alpha=0.5, linewidth=0,transform=axes[0].get_xaxis_transform())
lim = axes[0].get_ylim()
axes[0].set_ylim(0,lim[1])
lim = axes[1].get_ylim()
axes[1].set_ylim(0,lim[1])
lim = axes[4].get_ylim()
axes[4].set_ylim(0,lim[1])
microbeShannon.plot(x='t',ax = axes[2],legend=False,color=pal[0],linewidth=0.75)
axes[2].set_ylabel(ylabel ='Microbe Shannon Diversity',labelpad=15,fontsize=7)
axes[2].set_xlabel(xlabel = 'Time t',fontsize=7)
axes[2].ticklabel_format(style='sci',scilimits=(0,0))
axes[3].set_yticks([])
axes[3].fill_between(microbeWalls["t"],0,1, where=microbeWalls["wall_presence"]>0,interpolate=True, color="grey", alpha=0.5, linewidth=0,transform=axes[2].get_xaxis_transform())
axes[5].axhline(shannonThreshold, color='k', lw=.5, linestyle='dotted')
axes[5].set_xticks([])
lim = axes[2].get_ylim()
axes[2].set_ylim(0,lim[1])
lim = axes[3].get_ylim()
axes[3].set_ylim(0,lim[1])
lim = axes[5].get_ylim()
axes[5].set_ylim(0,lim[1])
fig.tight_layout()
fig.savefig(os.path.join(PLOT_PATH,'microbe-walls-shannon.png'),dpi=resolve)



print('Compiling microbial wall plot projection on viral abundance and shannon diversity')
fig, ax = plt.subplots(2,sharex=True)
fig.suptitle('(run{0}-c{1}-r{2})'.format(run_id,combo_id,replicate))
axes = [ax[0], ax[0].twinx(), ax[1], ax[1].twinx()]
virus_stacked.plot.area(ax = axes[0],stacked=True,legend=False, linewidth=0,color=pal)
axes[0].set_ylabel(ylabel ='Viral Strain Abundances V_i',labelpad=15,fontsize=7)
axes[0].set_xlabel(xlabel = 'Time t',fontsize=7)
axes[0].ticklabel_format(style='sci',scilimits=(0,0))
axes[1].set_yticks([])
axes[1].fill_between(microbeWalls["t"],0,1, where=microbeWalls["wall_presence"]>0,interpolate=True, color="grey", alpha=0.5, linewidth=0,transform=axes[0].get_xaxis_transform())
lim = axes[0].get_ylim()
axes[0].set_ylim(0,lim[1])
lim = axes[1].get_ylim()
axes[1].set_ylim(0,lim[1])
virusShannon.plot(x='t',ax = axes[2],legend=False,color=pal[0],linewidth=0.75)
axes[2].set_ylabel(ylabel ='Virus Shannon Diversity',labelpad=15,fontsize=7)
axes[2].set_xlabel(xlabel = 'Time t',fontsize=7)
axes[2].ticklabel_format(style='sci',scilimits=(0,0))
axes[3].set_yticks([])
axes[3].fill_between(microbeWalls["t"],0,1, where=microbeWalls["wall_presence"]>0,interpolate=True, color="grey", alpha=0.5, linewidth=0,transform=axes[2].get_xaxis_transform())
lim = axes[2].get_ylim()
axes[2].set_ylim(0,lim[1])
lim = axes[3].get_ylim()
axes[3].set_ylim(0,lim[1])
fig.tight_layout()
fig.savefig(os.path.join(PLOT_PATH,'microbe-walls-on-virus-shannon.png'),dpi=resolve)


print('Compiling microbial peak plot projection on viral abundance and shannon diversity')
fig, ax = plt.subplots(2,sharex=True)
fig.suptitle('(run{0}-c{1}-r{2})'.format(run_id,combo_id,replicate))
axes = [ax[0], ax[0].twinx(), ax[1], ax[1].twinx()]
virus_stacked.plot.area(ax = axes[0],stacked=True,legend=False, linewidth=0,color=pal)
axes[0].set_ylabel(ylabel ='Viral Strain Abundances V_i',labelpad=15,fontsize=7)
axes[0].set_xlabel(xlabel = 'Time t',fontsize=7)
axes[0].ticklabel_format(style='sci',scilimits=(0,0))
axes[1].set_yticks([])
axes[1].fill_between(microbePeaks["t"],0,1, where=microbePeaks["peak_presence"]>0,interpolate=True, color="grey", alpha=0.5, linewidth=0,transform=axes[0].get_xaxis_transform())
lim = axes[0].get_ylim()
axes[0].set_ylim(0,lim[1])
lim = axes[1].get_ylim()
axes[1].set_ylim(0,lim[1])
virusShannon.plot(x='t',ax = axes[2],legend=False,color=pal[0],linewidth=0.75)
axes[2].set_ylabel(ylabel ='Virus Shannon Diversity',labelpad=15,fontsize=7)
axes[2].set_xlabel(xlabel = 'Time t',fontsize=7)
axes[2].ticklabel_format(style='sci',scilimits=(0,0))
axes[3].set_yticks([])
axes[3].fill_between(microbePeaks["t"],0,1, where=microbePeaks["peak_presence"]>0,interpolate=True, color="grey", alpha=0.5, linewidth=0,transform=axes[2].get_xaxis_transform())
lim = axes[2].get_ylim()
axes[2].set_ylim(0,lim[1])
lim = axes[3].get_ylim()
axes[3].set_ylim(0,lim[1])
fig.tight_layout()
fig.savefig(os.path.join(PLOT_PATH,'microbe-peaks-on-virus-shannon.png'),dpi=resolve)


print('Compiling microbial-virus stacked time series plots')
fig, ax = plt.subplots(2,sharex=True)
fig.suptitle('(run{0}-c{1}-r{2})'.format(run_id,combo_id,replicate))
axes = [ax[0], ax[0].twinx(), ax[1], ax[1].twinx()]
microbe_stacked.plot.area(ax = axes[0],stacked=True,legend=False, linewidth=0,color=pal)
axes[0].set_ylabel(ylabel ='Microbial Immune Abundances N_i',labelpad=15,fontsize=7)
axes[0].set_xlabel(xlabel = 'Time t',fontsize=7)
axes[0].ticklabel_format(style='sci',scilimits=(0,0))
axes[1].set_yticks([])
axes[1].fill_between(microbeWalls["t"],0,1, where=microbeWalls["wall_presence"]>0,interpolate=True, color="grey", alpha=0.5, linewidth=0,transform=axes[0].get_xaxis_transform())
lim = axes[0].get_ylim()
axes[0].set_ylim(0,lim[1])
lim = axes[1].get_ylim()
axes[1].set_ylim(0,lim[1])
virus_stacked.plot.area(ax = axes[2],stacked=True,legend=False, linewidth=0,color=pal)
axes[2].set_ylabel(ylabel ='Viral Strain Abundances V_i',labelpad=15,fontsize=7)
axes[2].set_xlabel(xlabel = 'Time t',fontsize=7)
axes[2].ticklabel_format(style='sci',scilimits=(0,0))
axes[3].set_yticks([])
axes[3].fill_between(microbeWalls["t"],0,1, where=microbeWalls["wall_presence"]>0,interpolate=True, color="grey", alpha=0.5, linewidth=0,transform=axes[2].get_xaxis_transform())
lim = axes[2].get_ylim()
axes[2].set_ylim(0,lim[1])
lim = axes[3].get_ylim()
axes[3].set_ylim(0,lim[1])
fig.tight_layout()
fig.savefig(os.path.join(PLOT_PATH,'walls-on-microbe-virus-stacked.png'),dpi=resolve)



print('Complete!')
