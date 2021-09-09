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

con = sqlite3.connect("output.sqlite")
cur = con.cursor()


print('SQLite Query: virus abundance data')
virus = pd.read_sql_query("SELECT t,vstrain_id,abundance FROM vabundance", con)
virus_stacked = virus.pivot(index='t',columns='vstrain_id',values='abundance')

print('Compiling viral strain Abundances plot')
pal = sns.color_palette("tab20b")
virus_stacked.plot.area(stacked=True, legend=False, linewidth=0,color=pal)
plt.title('Viral Strain Abundances:')
plt.xlabel('Time t')
plt.ylabel('Abundances V_i')
plt.tight_layout()
plt.savefig(os.path.join('virus-strain-abundances_plot.png'),dpi=500)
#plt.show()

print('SQLite Query: microbe abundance data')
microbe = pd.read_sql_query("SELECT t,bstrain_id,abundance FROM babundance", con)
microbe_stacked = microbe.pivot(index='t',columns='bstrain_id',values='abundance')

print('Compiling microbial strain abundances plot')
microbe_stacked.plot.area(stacked=True, legend=False, linewidth=0,color=pal)
plt.title('Microbial Strain Abundances:')
plt.xlabel('Time t')
plt.ylabel('Abundances N_i')
plt.tight_layout()
plt.savefig(os.path.join('microbe-strain-abundances_plot.png'),dpi=500)
#plt.show()
