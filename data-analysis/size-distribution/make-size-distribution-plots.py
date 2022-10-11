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
from scipy.signal import find_peaks, peak_widths
from scipy.interpolate import interp1d, CubicSpline
import scipy.integrate as integrate
import sqlite3
import matplotlib.ticker as ticker
import colorsys
from mpl_toolkits.axes_grid1 import make_axes_locatable
import math
from scipy.optimize import curve_fit

run_id = sys.argv[1]
resolve = 500
imgTypes = ["png"]

SCRIPT_PATH = os.path.abspath(os.path.dirname(__file__)) # cluster
dir = 'crispr-sweep-7-2-2022/isolates/runID3297-c66-r47' # local
run = 'runID3297-c66-r47' # local

SCRIPT_PATH = os.path.abspath(os.path.dirname(__file__)) # cluster
DBSIM_PATH = os.path.join(SCRIPT_PATH,'..','..','simulation','sweep_db_gathered.sqlite') # cluster
#DBSIM_PATH = os.path.join('/Volumes','Yadgah','sweep_db_gathered.sqlite') # local
# DBSIM_PATH = os.path.join('/Volumes','Yadgah',dir,'{}.sqlite'.format(run)) # local

# if file exists...
DBOUTPUT_PATH = os.path.join('size-distribution_output.sqlite')

conSim = sqlite3.connect(DBSIM_PATH)
curSim = conSim.cursor()
ID = curSim.execute('SELECT combo_id,replicate FROM runs WHERE run_id = {}'.format(run_id)).fetchall()
combo_id = ID[0][0]
replicate = ID[0][1]
CC = curSim.execute('SELECT microbe_carrying_capacity FROM param_combos WHERE combo_id = {}'.format(combo_id)).fetchall()
CC = CC[0][0]
RUN_DIR = os.path.join('runID{0}-c{1}-r{2}'.format(run_id,combo_id,replicate))

freqDisplaced =  pd.read_sql_query("SELECT t,microbial_abundance FROM summary \
WHERE run_id = {0}".format(run_id),conSim).rename(columns={"microbial_abundance": "bTotal"})

freqDisplaced['freq_displaced'] = (CC-np.array(freqDisplaced['bTotal']))/CC

freqDisplaced = freqDisplaced[freqDisplaced['t'] <= max(virus_stacked.index)]

# grid_kws = {"height_ratios": (1, 1, 1), "hspace": (.1, .1)}
pal = sns.color_palette("tab20b")
# fig, ax = plt.subplots(2,sharex=True)
# axes = [ax[0], ax[1]]
fig, ax = plt.subplots(2,sharex=True)
axes = [ax[0], ax[1]]

microbe_stacked.plot.area(ax = axes[0],stacked=True,legend=False, linewidth=0,color=pal,sort_columns=True)
axes[0].set_ylabel(ylabel ='Microbial Strain Abundances',labelpad=15,fontsize=7)
axes[0].set_xlabel(xlabel = 'Time t',fontsize=7)
axes[0].ticklabel_format(axis = 'y',style='sci',scilimits=(0,0))
axes[0].xaxis.set_minor_locator(ticker.MultipleLocator(25))
lim = axes[0].get_ylim()
axes[0].set_ylim(0,lim[1])
# axes[0].set_xlim(0,350.0)
# spacerHeat = pextWeighted.replace({"spacer_id": newSpacerIDs})
# # spacerHeat['p_spacer'] = list(map(np.log,spacerHeat.p_spacer))
# spacerHeat = spacerHeat.pivot_table(index='spacer_id', columns='t', values='p_spacer')
# axes[1].imshow(spacerHeat,cmap='viridis', aspect='auto')
# axes[1].set_xlim(0,350.0)
burnFreq.plot(x='t',ax = axes[1],legend=False,linewidth=0.75)
axes[1].plot(freqDisplaced['t'],freqDisplaced['freq_displaced'],linewidth=0.75)



peaks, _ = find_peaks(freqDisplaced['freq_displaced'], height=0)
peakWidths = peak_widths(freqDisplaced['freq_displaced'], peaks, rel_height=1.0)



plt.plot(freqDisplaced['freq_displaced'])
plt.plot(peaks, freqDisplaced['freq_displaced'][peaks],"x")
plt.hlines(*peakWidths[1:], color="C2")


interpBurnFreq = interp1d(freqDisplaced['t'],freqDisplaced['freq_displaced'])


epiSizes = []
for k in range(0,len(peakWidths[1:][1])):
    ti = round(peakWidths[1:][1][k] + min(freqDisplaced['t']),4)
    tj = round(peakWidths[1:][2][k] + min(freqDisplaced['t']),4)
    if tj == max(freqDisplaced['t']):
        epiSizes.extend(integrate.quad(interpBurnFreq, ti, tj-.001)[0])
        break
    size = integrate.quad(interpBurnFreq, ti, tj)[0]
    epiSizes.append(size)


# CDF

epiSizeData = pd.DataFrame()
counts = []
bin0 = 0.01
binf = 13
binSize = 0.01
bins = np.arange(bin0,binf,binSize)
for i in bins:
    sizes = list(filter(lambda x: x > i, epiSizes))
    counts.append(len(sizes))



epiSizeData['sizes'] = bins[:]
epiSizeData['count'] = counts[:]
fig, ax = plt.subplots(1)
ax.loglog(epiSizeData['sizes'], epiSizeData['count'],"o")


epiSizeData = pd.DataFrame()
epiSizeData['sizes'] = list(np.around(np.array(epiSizes), 1))
epiSizeData['count'] = epiSizeData.groupby('sizes')['sizes'].transform('count')
epiSizeData = epiSizeData.sort_values(by=['sizes'])
fig, ax = plt.subplots(1)
ax.loglog(epiSizeData['sizes'], epiSizeData['count'],"o")

ts = np.arange(min(burnFreq['t']), max(burnFreq['t']), 0.01)
fig, ax = plt.subplots(1)
ax.plot(ts, interpBurnFreq(ts))


conn = sqlite3.connect(DBOUTPUT_PATH)
epiSizeData.to_sql('displacement_sizes', conn, if_exists='replace', index=False)


conn = sqlite3.connect(DBOUTPUT_PATH)
df.to_sql('epidemic_sizes', conn, if_exists='replace', index=False)


# for k in range(0,len(peakWidths[1:][1])):
#     print(k)
#
# for (i,j) in zip(peakWidths[1:][1],peakWidths[1:][2]):
#     print(i)
#     print(j)
#
#
# for k in range(0,len(peakWidths[1:][1])):
#     ti = peakWidths[1:][1][k] + min(burnFreq['t'])
#     tj = peakWidths[1:][2][k] + min(burnFreq['t'])
#     print(ti,tj)



print('Complete!')
