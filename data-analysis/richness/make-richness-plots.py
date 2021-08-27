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
DBSIM_PATH = os.path.join(SCRIPT_PATH,'..','..','simulation','sweep_db_gathered.sqlite')
conSim = sqlite3.connect(DBSIM_PATH)
curSim = conSim.cursor()

ID = curSim.execute('SELECT combo_id,replicate FROM runs WHERE run_id = {}'.format(run_id)).fetchall()
#ID = pd.read_sql_query("SELECT combo_id,replicate FROM runs WHERE run_id = {}".format(run_id), con)
combo_id = ID[0][0]
replicate = ID[0][1]

PLOT_PATH = os.path.join(SCRIPT_PATH,'plots','c{}'.format(combo_id),'r{}'.format(replicate))
#PLOT_PATH = os.path.abspath(os.path.dirname(__file__))

DBANALYSIS_PATH = os.path.join(SCRIPT_PATH,'..','data_analysis.sqlite')
conAn = sqlite3.connect(DBANALYSIS_PATH)
curAn = conAn.cursor()

print('SQLite Query: Virus Data')
virus = pd.read_sql_query("SELECT t, virus FROM richness WHERE run_id = {}".format(run_id), conAn)

print('Processing: Viral Strain Richness Plot')
pal = sns.color_palette("tab20b")
virus.plot(x='t',legend=False,color=pal)
plt.title('Viral Strain Richness (run {0}: c{1}-r{2})'.format(run_id,combo_id,replicate))
plt.xlabel('Time t')
plt.ylabel('Richness')
plt.tight_layout()
plt.savefig(os.path.join(PLOT_PATH,'virus-strain-richness.png'),dpi=500)
#plt.show()

print('SQLite Query: Microbe Immune Richness Data')
microbe = pd.read_sql_query("SELECT t,bstrain_id,abundance FROM babundance WHERE run_id = {}".format(run_id), conAn)

print('Processing: Microbial Immune Richness Plot')
microbe.plot(x='t',legend=False, color=pal)
plt.title('Microbial Immune Richness (run {0}: c{1}-r{2})'.format(run_id,combo_id,replicate))
plt.xlabel('Time t')
plt.ylabel('Richness')
plt.tight_layout()
plt.savefig(os.path.join(PLOT_PATH,'microbe-strain-richness.png'),dpi=500)
#plt.show()
