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

SCRIPT_PATH = os.path.abspath(os.path.dirname(__file__))

DBSIM_PATH = os.path.join(SCRIPT_PATH,'..','..','simulation','sweep_db_gathered.sqlite') # cluster
#home_dir = os.system("cd ~") # local
#DBSIM_PATH = os.path.join('/Volumes','Yadgah','sweep_db_gathered.sqlite') # local
conSim = sqlite3.connect(DBSIM_PATH)
curSim = conSim.cursor()

DBAnalysis_PATH = os.path.join(SCRIPT_PATH,'..','gathered-analyses','shannon','shannon.sqlite') # cluster
#DBAnalysis_PATH = os.path.join('/Volumes','Yadgah','shannon.sqlite') # local
#DBAnalysis_PATH = os.path.join('/Volumes','Yadgah','shannon_output.sqlite') # local
conAnalysis = sqlite3.connect(DBAnalysis_PATH)
curAnalysis = conAnalysis.cursor()

ID = curSim.execute('SELECT combo_id,replicate FROM runs WHERE run_id = {}'.format(run_id)).fetchall()
combo_id = ID[0][0]
replicate = ID[0][1]

PLOT_PATH = os.path.join(SCRIPT_PATH,'..','gathered-analyses','shannon','plots','c{}'.format(combo_id),'r{}'.format(replicate)) # cluster
#PLOT_PATH = os.path.abspath(os.path.dirname(__file__)) # local

#print('SQLite Query: virus shannon data') # local
#virusAnalysis = pd.read_sql_query("SELECT t, vshannon FROM shannon", conAnalysis) # local
#print('SQLite Query: microbe shannon data') # local
#microbeAnalysis = pd.read_sql_query("SELECT t,bshannon FROM shannon", conAnalysis) # local

print('SQLite Query: virus shannon data') # cluster
virusAnalysis = pd.read_sql_query("SELECT t, vshannon FROM shannon WHERE run_id = {}".format(run_id), conAnalysis) # cluster
print('SQLite Query: microbe shannon data') # cluster
microbeAnalysis = pd.read_sql_query("SELECT t,bshannon FROM shannon WHERE run_id = {}".format(run_id), conAnalysis) # cluster

print('SQLite Query: microbial abundance time series data')
microbeSim = pd.read_sql_query("SELECT t,bstrain_id,abundance FROM babundance WHERE run_id = {}".format(run_id), conSim)
microbe_stacked = microbeSim.pivot(index='t',columns='bstrain_id',values='abundance')

print('SQLite Query: viral abundance time series data')
virusSim = pd.read_sql_query("SELECT t,vstrain_id,abundance FROM vabundance WHERE run_id = {}".format(run_id), conSim)
virus_stacked = virusSim.pivot(index='t',columns='vstrain_id',values='abundance')

print('Compiling viral shannon entropy plot')
pal = sns.color_palette("tab20b")
virusAnalysis.plot(x='t',xlabel = 'Time t',ylabel = 'Shannon Entropy Hv', legend=False,color=pal[3])
plt.title('Viral Shannon Entropy (run{0}-c{1}-r{2})'.format(run_id,combo_id,replicate))
plt.tight_layout()
plt.savefig(os.path.join(PLOT_PATH,'virus-shannon.png'),dpi=500)
#plt.show()

print('Compiling microbial shannon plot')
microbeAnalysis.plot(x='t',xlabel = 'Time t',ylabel = 'Shannon Entropy Hn', legend=False, color=pal)
plt.title('Microbial Shannon Entropy (run{0}-c{1}-r{2})'.format(run_id,combo_id,replicate))
plt.tight_layout()
plt.savefig(os.path.join(PLOT_PATH,'microbe-shannon.png'),dpi=500)
#plt.show()

print('Compiling shannon subplots')
fig, axs = plt.subplots(2,sharex=True)
axs[0].set(ylabel ='Microbial Shannon Entropy Hn')
axs[1].set(ylabel ='Viral Shannon Entropy Hv')
fig.suptitle('Shannon Entropy (run{0}-c{1}-r{2})'.format(run_id,combo_id,replicate))
microbeAnalysis.plot(x='t',xlabel = 'Time t',ax = axs[0],legend=False,color=pal[0],linewidth=0.75)
axs[0].set_ylabel(ylabel ='Microbial Shannon Entropy Hn',labelpad=15,fontsize=7)
axs[0].set_xlabel(xlabel = 'Time t',fontsize=7)
axs[1].ticklabel_format(style='sci',scilimits=(0,0))
virusAnalysis.plot(x = 't',xlabel = 'Time t',ax = axs[1],legend=False,color=pal[3],linewidth=0.75)
axs[1].set_ylabel(ylabel ='Viral Shannon Entropy Hv',labelpad=15,fontsize=7)
axs[1].set_xlabel(xlabel = 'Time t',fontsize=7)
axs[1].ticklabel_format(style='sci',scilimits=(0,0))
plt.tight_layout()
plt.savefig(os.path.join(PLOT_PATH,'microbe-virus-shannon-stacked.png'),dpi=500)

print('Compiling time series and shannon subplots')
fig, axs = plt.subplots(2,sharex=True)
fig.suptitle('(run{0}-c{1}-r{2})'.format(run_id,combo_id,replicate))

microbe_stacked.plot.area(stacked=True, xlabel = 'Time t',ax = axs[0], legend=False, linewidth=0,color=pal)
axs[0].set_ylabel(ylabel ='Microbial Immune Abundances N_i',labelpad=15,fontsize=7)
axs[0].set_xlabel(xlabel = 'Time t',fontsize=7)
axs[0].ticklabel_format(style='sci',scilimits=(0,0))
microbeAnalysis.plot(x='t',xlabel = 'Time t',ax = axs[1],legend=False,color=pal[0],linewidth=0.75)
axs[1].set_ylabel(ylabel ='Microbial Shannon Entropy Hn',labelpad=15,fontsize=7)
axs[1].set_xlabel(xlabel = 'Time t',fontsize=7)
axs[1].ticklabel_format(style='sci',scilimits=(0,0))
plt.tight_layout()
plt.savefig(os.path.join(PLOT_PATH,'microbe-tSeries-shannon-stacked.png'),dpi=500)

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
axs[1].set_ylabel(ylabel ='Viral Shannon Entropy Hv',rotation=270,labelpad=15,fontsize=7)
axs[1].yaxis.tick_right()
axs[1].set_xlabel(xlabel = 'Time t',fontsize=7)
axs[1].ticklabel_format(style='sci',scilimits=(0,0))
plt.tight_layout()
plt.savefig(os.path.join(PLOT_PATH,'virus-tSeries-shannon-stacked.png'),dpi=500)

print('Complete!')
