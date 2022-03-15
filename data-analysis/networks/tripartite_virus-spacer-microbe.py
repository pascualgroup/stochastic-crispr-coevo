#!/usr/bin/env python3

import pandas as pd
import numpy as np
import scipy as sp
import math
import matplotlib.pyplot as plt
from matplotlib import ticker
import sys
import os
import seaborn as sns
from scipy import stats
import sqlite3
from mpl_toolkits.axes_grid1 import make_axes_locatable

run_id = sys.argv[1]
time = sys.argv[2]

SCRIPT_PATH = os.path.abspath(os.path.dirname(__file__)) # cluster

DB_PATH = os.path.join(SCRIPT_PATH,'..','sweep_db_gathered.sqlite')
# DB_PATH = os.path.join('/Volumes','Yadgah','run_id1455_combo73_replicate15.sqlite')
# DB_PATH = os.path.join('/Volumes','Yadgah','crispr-sweep-7-2-2022/isolated-runs/
# run_id3297_combo66_replicate47/run_id3297_combo66_replicate47.sqlite')

DBMATCH_PATH = os.path.join(SCRIPT_PATH,'..','gathered-analyses','matches','matches.sqlite') # cluster
#DBMATCH_PATH = os.path.join('/Volumes','Yadgah','matches.sqlite') # local

DBCLADE_PATH = os.path.join(SCRIPT_PATH,'..','gathered-analyses','clade-abundances','clade-abundances.sqlite') # cluster
# DBCLADE_PATH = os.path.join('/Volumes','Yadgah','clade-abundances.sqlite') # local
# DBCLADE_PATH = os.path.join('/Volumes','Yadgah','clade-abundances_output.sqlite') # local. run_id fixed; for testing

DBTREE_PATH = os.path.join(SCRIPT_PATH,'..','gathered-analyses','trees','trees.sqlite') # cluster
# DBTREE_PATH = os.path.join('/Volumes','Yadgah','trees.sqlite') # local
# DBTREE_PATH = os.path.join('/Volumes','Yadgah','trees_output.sqlite') # local. run_id fixed; for testing

con = sqlite3.connect(DB_PATH)
cur = con.cursor()

conMatch = sqlite3.connect(DBMATCH_PATH)
curMatch = conMatch.cursor()

conClade = sqlite3.connect(DBCLADE_PATH)
curClade = conClade.cursor()

print('SQLite Query: clade data')
microbeClades = pd.read_sql_query("SELECT DISTINCT clade_id, bstrain_id \
FROM babundances", conClade)
microbeCladeIDs = pd.read_sql_query("SELECT DISTINCT clade_id \
FROM babundances", conClade)
microbeCladeAbundances = pd.read_sql_query("SELECT t, clade_id, abundance \
FROM clade_babundances", conClade)

virusClades = pd.read_sql_query("SELECT DISTINCT clade_id, vstrain_id \
FROM vabundances", conClade)
virusCladeIDs = pd.read_sql_query("SELECT DISTINCT clade_id \
FROM vabundances", conClade)
virusCladeAbundances = pd.read_sql_query("SELECT t, clade_id, abundance \
FROM clade_vabundances", conClade)


## for runs with run_id
# conTree = sqlite3.connect(DBTREE_PATH)
# curTree = conTree.cursor()
# print('SQLite Query: tree data')
# bStrainTimes = pd.read_sql_query(
# "SELECT tree_bstrain_id, t_creation, t_extinction \
# FROM tree_bstrain_creation_extinction WHERE run_id = {}".format(run_id), conTree)
# bTreeAbundances = pd.read_sql_query(
# "SELECT t, tree_bstrain_id, abundance \
# FROM tree_babundance WHERE run_id = {}".format(run_id), conTree)
# vStrainTimes = pd.read_sql_query(
# "SELECT tree_vstrain_id, t_creation, t_extinction \
# FROM tree_vstrain_creation_extinction WHERE run_id = {}".format(run_id), conTree)
# vTreeAbundances = pd.read_sql_query(
# "SELECT t, tree_vstrain_id, abundance \
# FROM tree_vabundance WHERE run_id = {}".format(run_id), conTree

conTree = sqlite3.connect(DBTREE_PATH)
curTree = conTree.cursor()
bStrainTimes = pd.read_sql_query(
"SELECT tree_bstrain_id, t_creation, t_extinction, tree_parent_bstrain_id \
FROM tree_bstrain_creation_extinction", conTree)
bTreeAbundances = pd.read_sql_query(
"SELECT t, tree_bstrain_id, abundance \
FROM tree_babundance", conTree)
vStrainTimes = pd.read_sql_query(
"SELECT tree_vstrain_id, t_creation, t_extinction, tree_parent_vstrain_id \
FROM tree_vstrain_creation_extinction", conTree)
vTreeAbundances = pd.read_sql_query(
"SELECT t, tree_vstrain_id, abundance \
FROM tree_vabundance", conTree)
if bTreeAbundances['t'][bTreeAbundances['t'].size-1] not in vTreeAbundances.t.values:
    bTreeAbundances = bTreeAbundances[bTreeAbundances.t != bTreeAbundances['t'][bTreeAbundances['t'].size-1]]

if vTreeAbundances['t'][vTreeAbundances['t'].size-1] not in bTreeAbundances.t.values:
    vTreeAbundances = vTreeAbundances[vTreeAbundances.t != vTreeAbundances['t'][vTreeAbundances['t'].size-1]]





# Designating plot path from simulation data
ID = curSim.execute('SELECT combo_id,replicate FROM runs WHERE run_id = {}'.format(run_id)).fetchall()
combo_id = ID[0][0]
replicate = ID[0][1]


# strainsOfClade = curClade.execute("SELECT DISTINCT bstrain_id FROM babundances \
# WHERE clade_id = {0} AND t = {1} AND run_id = {2} ORDER BY tree_bstrain_id".format(clade_id,time,run_id))




# singleMatchIDs = [i[0] for i in curMatch.execute("SELECT DISTINCT time_specific_match_id FROM matches_spacers \
# WHERE t = {0} AND match_type = 1 AND run_id = {1}".format(time,run_id)).fetchall()]
# bstrains0 = [i[0] for i in curMatch.execute("SELECT bstrain_id FROM bstrain_to_vstrain_0matches \
# WHERE t = {0} AND run_id = {1}".format(time,run_id))]
# bstrains1 = [i[0] for i in curMatch.execute("SELECT bstrains FROM bstrain_to_vstrain_matches \
# WHERE t = {0} AND time_specific_match_id in ({1}) AND run_id = {2}".format(time,', '.join(map(str,singleMatchIDs)),run_id)).fetchall()]
# vstrains0 = [i[0] for i in curMatch.execute("SELECT vstrain_id FROM bstrain_to_vstrain_0matches \
# WHERE t = {0} AND run_id = {1}".format(time,run_id))
# vstrains1 = [i[0] for i in curMatch.execute("SELECT vstrains FROM bstrain_to_vstrain_matches \
# WHERE t = {0} AND time_specific_match_id in ({1}) AND run_id = {2}".format(time,', '.join(map(str,singleMatchIDs)),run_id)).fetchall()]



singleMatchIDs = [i[0] for i in curMatch.execute("SELECT DISTINCT time_specific_match_id FROM matches_spacers \
WHERE t = {0} AND match_type = 1".format(time)).fetchall()]
bstrains0 = [i[0] for i in curMatch.execute("SELECT bstrain_id FROM bstrain_to_vstrain_0matches \
WHERE t = {0}".format(time))]
bstrains1 = [i[0] for i in curMatch.execute("SELECT bstrains FROM bstrain_to_vstrain_matches \
WHERE t = {0} AND time_specific_match_id in ({1})".format(time,', '.join(map(str,singleMatchIDs)))).fetchall()]
vstrains0 = [i[0] for i in curMatch.execute("SELECT vstrain_id FROM bstrain_to_vstrain_0matches \
WHERE t = {0}".format(time))
vstrains1 = [i[0] for i in curMatch.execute("SELECT vstrains FROM bstrain_to_vstrain_matches \
WHERE t = {0} AND time_specific_match_id in ({1})".format(time,', '.join(map(str,singleMatchIDs)))).fetchall()]






bStrainsOrdered = curTree.execute("SELECT bstrain_id FROM tree_bstrain_order \
WHERE bstrain_id in ({}) ORDER BY tree_bstrain_id".format(', '.join(map(str,treeBstrains))))

treeVstrains = curTree.execute("SELECT tree_vstrain_id FROM tree_vabundance \
WHERE t = {} ORDER BY tree_vstrain_id".format(time))

vStrainsOrdered = curTree.execute("SELECT vstrain_id FROM tree_vstrain_order \
WHERE bvtrain_id in ({}) ORDER BY tree_vstrain_id".format(', '.join(map(str,treeVstrains))))


triPartiteSingleMatch = {}
for i in bStrainsOrdered:



for clade_id in np.unique(microbeClades.clade_id):
    # strainsOfClade = curClade.execute("SELECT DISTINCT bstrain_id FROM babundances \
    # WHERE clade_id = {} AND run_id = {}".format(clade_id,run_id))
    strainsOfClade = curClade.execute("SELECT DISTINCT bstrain_id FROM babundances \
    WHERE clade_id = {0} AND t = {1} ORDER BY tree_bstrain_id".format(clade_id,time))
    strainsOfClade = [i[0] for i in strainsOfClade.fetchall()]
    for strain in treeStrainsOfClade:
        print(strain)
        tCreate = bStrainTimes[bStrainTimes.tree_bstrain_id == strain].t_creation.values[0]
        tExtinct = bStrainTimes[bStrainTimes.tree_bstrain_id == strain].t_extinction.values[0]
        parent = bStrainTimes[bStrainTimes.tree_bstrain_id == strain].tree_parent_bstrain_id.values[0]
        hlinecMicrobe.append([[tCreate, strain],[tExtinct, strain]])
        vlinecMicrobe.append([[tCreate, parent],[tCreate, strain]])
        hcolorsMicrobe.append(Mcmap(Mnorm(np.arange(1, len(microbeCladeIDs)+1, 1)))[cladeColorDict[clade_id]])
        vcolorsMicrobe.append(Mcmap(Mnorm(np.arange(1, len(microbeCladeIDs)+1, 1)))[cladeColorDict[clade_id]])
        markerColorsMicrobe.append(Mcmap(Mnorm(np.arange(1, len(microbeCladeIDs)+1, 1)))[cladeColorDict[clade_id]])
