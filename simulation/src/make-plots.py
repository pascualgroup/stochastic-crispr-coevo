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

#con = sqlite3.connect('/Volumes/Yadgah/sweep_db_gathered.sqlite')
SCRIPT_PATH = os.path.abspath(os.path.dirname(__file__))
DB_PATH = os.path.join(SCRIPT_PATH,'..','sweep_db_gathered.sqlite')
con = sqlite3.connect(DB_PATH)
cur = con.cursor()

ID = cur.execute('SELECT combo_id,replicate FROM runs WHERE run_id = {}'.format(run_id)).fetchall()
#ID = pd.read_sql_query("SELECT combo_id,replicate FROM runs WHERE run_id = {}".format(run_id), con)

combo_id = ID[0][0]
replicate = ID[0][1]

PLOT_PATH = os.path.join(SCRIPT_PATH,'..', 'plots','c{}'.format(combo_id),'r{}'.format(replicate))

#PLOT_PATH = os.path.abspath(os.path.dirname(__file__))

print('SQLite Query: virus abundance data')
virus = pd.read_sql_query("SELECT t,vstrain_id,abundance FROM vabundance WHERE run_id = {}".format(run_id), con)
virus_stacked = virus.pivot(index='t',columns='vstrain_id',values='abundance')

print('Compiling viral strain Abundances plot')
pal = sns.color_palette("tab20b")
virus_stacked.plot.area(stacked=True, legend=False, linewidth=0,color=pal)
plt.title('Viral Strain Abundances')
plt.xlabel('Time t')
plt.ylabel('Abundances V_i')
plt.tight_layout()
plt.savefig(os.path.join(PLOT_PATH,'virus-strain-abundances_plot.png'),dpi=500)
#plt.show()

print('SQLite Query: microbe abundance data')
microbe = pd.read_sql_query("SELECT t,bstrain_id,abundance FROM babundance WHERE run_id = {}".format(run_id), con)
microbe_stacked = microbe.pivot(index='t',columns='bstrain_id',values='abundance')

print('Compiling microbial strain abundances plot')
microbe_stacked.plot.area(stacked=True, legend=False, linewidth=0,color=pal)
plt.title('Microbial Strain Abundances')
plt.xlabel('Time t')
plt.ylabel('Abundances N_i')
plt.tight_layout()
plt.savefig(os.path.join(PLOT_PATH,'microbe-strain-abundances_plot.png'),dpi=500)
#plt.show()

print('Complete!')
