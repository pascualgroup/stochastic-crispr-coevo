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
#run_id = 1

SCRIPT_PATH = os.path.abspath(os.path.dirname(__file__))

DBSIM_PATH = os.path.join(SCRIPT_PATH,'..','..','simulation','sweep_db_gathered.sqlite') # cluster
#home_dir = os.system("cd ~") # local
#DBSIM_PATH = os.path.join('/Volumes','Yadgah','sweep_db_gathered.sqlite') # local
conSim = sqlite3.connect(DBSIM_PATH)
curSim = conSim.cursor()

ID = curSim.execute('SELECT combo_id,replicate FROM runs WHERE run_id = {}'.format(run_id)).fetchall()
combo_id = ID[0][0]
replicate = ID[0][1]

PLOT_PATH = os.path.join(SCRIPT_PATH,'..','gathered-analyses','richness','plots','c{}'.format(combo_id),'r{}'.format(replicate)) # cluster
#PLOT_PATH = os.path.abspath(os.path.dirname(__file__)) # local

print('SQLite Query: microbial abundance time series data')
microbeSim = pd.read_sql_query("SELECT t,bstrain_id,abundance FROM babundance WHERE run_id = {}".format(run_id), conSim)
microbe_stacked = microbeSim.pivot(index='t',columns='bstrain_id',values='abundance')

print('SQLite Query: viral abundance time series data')
virusSim = pd.read_sql_query("SELECT t,vstrain_id,abundance FROM vabundance WHERE run_id = {}".format(run_id), conSim)
virus_stacked = virusSim.pivot(index='t',columns='vstrain_id',values='abundance')


DBAnalysis_PATH = os.path.join(SCRIPT_PATH,'..','gathered-analyses','richness','richness.sqlite') # cluster
#DBAnalysis_PATH = os.path.join('/Volumes','Yadgah','richness.sqlite') # local
conAnalysis = sqlite3.connect(DBAnalysis_PATH)
curAnalysis = conAnalysis.cursor()

print('SQLite Query: virus strain richness data')
virusAnalysis = pd.read_sql_query("SELECT t, vrichness FROM richness WHERE run_id = {}".format(run_id), conAnalysis)

print('SQLite Query: microbe immune richness data')
microbeAnalysis = pd.read_sql_query("SELECT t,brichness FROM richness WHERE run_id = {}".format(run_id), conAnalysis)

print('Compiling viral strain richness plot')
pal = sns.color_palette("tab20b")
virusAnalysis.plot(x='t',xlabel = 'Time t',ylabel = 'Richness Sv', legend=False,color=pal[3])
plt.title('Viral Strain Richness (run{0}-c{1}-r{2})'.format(run_id,combo_id,replicate))
plt.tight_layout()
plt.savefig(os.path.join(PLOT_PATH,'virus-strain-richness.png'),dpi=500)
#plt.show()

print('Compiling microbial immune richness plot')
microbeAnalysis.plot(x='t',xlabel = 'Time t',ylabel = 'Richness Sn', legend=False, color=pal)
plt.title('Microbial Immune Richness (run{0}-c{1}-r{2})'.format(run_id,combo_id,replicate))
plt.tight_layout()
plt.savefig(os.path.join(PLOT_PATH,'microbe-strain-richness.png'),dpi=500)
#plt.show()

print('Compiling richness subplots')
fig, axs = plt.subplots(2,sharex=True)
axs[0].set(ylabel ='Microbial Immune Richness Sn')
axs[1].set(ylabel ='Viral Strain Richness Sv')
fig.suptitle('Richness (run{0}-c{1}-r{2})'.format(run_id,combo_id,replicate))
microbeAnalysis.plot(x='t',xlabel = 'Time t',ax = axs[0],legend=False,color=pal[0],linewidth=0.75)
axs[0].set_ylabel(ylabel ='Microbial Immune Richness Sn',labelpad=15,fontsize=7)
axs[0].set_xlabel(xlabel = 'Time t',fontsize=7)
axs[1].ticklabel_format(style='sci',scilimits=(0,0))
virusAnalysis.plot(x = 't',xlabel = 'Time t',ax = axs[1],legend=False,color=pal[3],linewidth=0.75)
axs[1].set_ylabel(ylabel ='Viral Strain Richness Sv',labelpad=15,fontsize=7)
axs[1].set_xlabel(xlabel = 'Time t',fontsize=7)
axs[1].ticklabel_format(style='sci',scilimits=(0,0))
plt.tight_layout()
plt.savefig(os.path.join(PLOT_PATH,'microbe-virus-richness-stacked.png'),dpi=500)

print('Compiling time series and richness subplots')
fig, axs = plt.subplots(2,sharex=True)
#axs[0,0].set(ylabel ='Microbial Immune Abundances N_i')
#axs[1,0].set(ylabel ='Microbial Immune Richness Sn')
#axs[0,1].set(ylabel ='Viral Strain Abundances V_i')
#axs[1,1].set(ylabel ='Viral Strain Richness Sv')
fig.suptitle('(run{0}-c{1}-r{2})'.format(run_id,combo_id,replicate))

microbe_stacked.plot.area(stacked=True, xlabel = 'Time t',ax = axs[0], legend=False, linewidth=0,color=pal)
axs[0].set_ylabel(ylabel ='Microbial Immune Abundances N_i',labelpad=15,fontsize=7)
axs[0].set_xlabel(xlabel = 'Time t',fontsize=7)
axs[0].ticklabel_format(style='sci',scilimits=(0,0))
microbeAnalysis.plot(x='t',xlabel = 'Time t',ax = axs[1],legend=False,color=pal[0],linewidth=0.75)
axs[1].set_ylabel(ylabel ='Microbial Immune Richness Sn',labelpad=15,fontsize=7)
axs[1].set_xlabel(xlabel = 'Time t',fontsize=7)
axs[1].ticklabel_format(style='sci',scilimits=(0,0))
plt.tight_layout()
plt.savefig(os.path.join(PLOT_PATH,'microbe-tSeries-richness-stacked.png'),dpi=500)

fig, axs = plt.subplots(2,sharex=True)
fig.suptitle('(run{0}-c{1}-r{2})'.format(run_id,combo_id,replicate))
virus_stacked.plot.area(stacked=True,xlabel = 'Time t',ax = axs[0], legend=False, linewidth=0,color=pal)
#axs[0,0].yaxis.set_label_position("right")
axs[0].set_ylabel(ylabel ='Viral Strain Abundances V_i',rotation=270,labelpad=15,fontsize=7)
axs[0].yaxis.set_label_position("right")
axs[0].yaxis.tick_right()
axs[0].tick_params(axis='x', labelsize=7)
axs[0].set_xlabel(xlabel = 'Time t',fontsize=7)
axs[0].ticklabel_format(style='sci',scilimits=(0,0))
virusAnalysis.plot(x = 't',xlabel = 'Time t',ax = axs[1],legend=False,color=pal[3],linewidth=0.75)
axs[1].yaxis.set_label_position("right")
axs[1].set_ylabel(ylabel ='Viral Strain Richness Sv',rotation=270,labelpad=15,fontsize=7)
axs[1].yaxis.tick_right()
axs[1].set_xlabel(xlabel = 'Time t',fontsize=7)
axs[1].ticklabel_format(style='sci',scilimits=(0,0))
plt.tight_layout()
plt.savefig(os.path.join(PLOT_PATH,'virus-tSeries-richness-stacked.png'),dpi=500)

print('Complete!')
