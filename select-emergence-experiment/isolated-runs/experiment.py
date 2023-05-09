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


comboID = sys.argv[1]

SCRIPT_PATH = os.path.abspath(os.path.dirname(__file__)) # cluster
# SQLITE path to simulation data
DBSIM_PATH = os.path.join(SCRIPT_PATH,'../simulation/sweep_db_gathered.sqlite')
DBSIM_PATH = os.path.join(
    '/Volumes/Yadgah/crispr-sweep-7-2-2022/emergence-experiment/simulation/sweep_db_gathered.sqlite')  # local
# Call SQLITE simulation database and some related info
conSim = sqlite3.connect(DBSIM_PATH)
curSim = conSim.cursor()





