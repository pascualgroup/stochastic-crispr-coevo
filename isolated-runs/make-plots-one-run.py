#!/usr/bin/env python3

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import seaborn as sns
import sqlite3

run = sys.argv[1]
DBSIM_PATH = os.path.join(SCRIPT_PATH,'isolates',run,'{}.sqlite'.format(run))
conSim = sqlite3.connect(DBSIM_PATH)
curSim = conSim.cursor()
ID = curSim.execute('SELECT combo_id,replicate FROM runs WHERE run_id = {}'.format(run_id)).fetchall()
combo_id = ID[0][0]
replicate = ID[0][1]

print('SQLite Query: virus abundance data')
virus = pd.read_sql_query("SELECT t,vstrain_id,abundance FROM vabundance", con)
virus_stacked = virus.pivot(index='t',columns='vstrain_id',values='abundance')

print('Compiling viral strain Abundances plot')
pal = sns.color_palette("tab20b")
fig, ax = plt.subplots(1)
virus_stacked.plot.area(ax = ax, stacked=True, legend=False, linewidth=0,color=pal)
plt.title('Viral Strain Abundances (run{0}-c{1}-r{2}):'.format(run_id,combo_id,replicate))
plt.xlabel('Time t')
plt.ylabel('Abundances V_i')
plt.tight_layout()
plt.savefig(os.path.join('virus-strain-abundances_plot.png'),dpi=500)
plt.close(fig)
#plt.show()

print('SQLite Query: microbe abundance data')
microbe = pd.read_sql_query("SELECT t,bstrain_id,abundance FROM babundance", con)
microbe_stacked = microbe.pivot(index='t',columns='bstrain_id',values='abundance')

print('Compiling microbial strain abundances plot')
fig, ax = plt.subplots(1)
microbe_stacked.plot.area(ax = ax, stacked=True, legend=False, linewidth=0,color=pal)
plt.title('Microbial Strain Abundances (run{0}-c{1}-r{2}):'.format(run_id,combo_id,replicate))
plt.xlabel('Time t')
plt.ylabel('Abundances N_i')
plt.tight_layout()
plt.savefig(os.path.join('microbe-strain-abundances_plot.png'),dpi=500)
plt.close(fig)
#plt.show()

print('Compiling microbial-virus stacked time series plots')
fig, ax = plt.subplots(2,sharex=True)
#fig = plt.figure()
fig.suptitle('Strain Abundances (run{0}-c{1}-r{2})'.format(run_id,combo_id,replicate))
axes = [ax[0], ax[1]]
microbe_stacked.plot.area(ax = axes[0],stacked=True,legend=False, linewidth=0,color=pal)
axes[0].set_ylabel(ylabel ='Microbial Immune Abundances N_i',labelpad=15,fontsize=7)
axes[0].set_xlabel(xlabel = 'Time t',fontsize=7)
axes[0].ticklabel_format(axis = 'y',style='sci',scilimits=(0,0))
axes[0].xaxis.set_minor_locator(ticker.MultipleLocator(25))
lim = axes[0].get_ylim()
axes[0].set_ylim(0,lim[1])
virus_stacked.plot.area(ax = axes[1],stacked=True,legend=False, linewidth=0,color=pal)
axes[1].set_ylabel(ylabel ='Viral Strain Abundances V_i',labelpad=15,fontsize=7)
axes[1].set_xlabel(xlabel = 'Time t',fontsize=7)
axes[1].ticklabel_format(axis = 'y',style='sci',scilimits=(0,0))
axes[1].xaxis.set_minor_locator(ticker.MultipleLocator(25))
lim = axes[1].get_ylim()
axes[1].set_ylim(0,lim[1])
fig.tight_layout()
fig.savefig(os.path.join(PLOT_PATH,'microbe-virus-stacked-abundances.png'),dpi=500)
plt.close(fig)
