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

DBSIM_PATH = os.path.join(SCRIPT_PATH,'..','..','..','sweep_db_gathered.sqlite') # cluster
# DBSIM_PATH = os.path.join('/Volumes','Yadgah',dir,'{}.sqlite'.format(run)) # local

conSim = sqlite3.connect(DBSIM_PATH)
curSim = conSim.cursor()
ID = curSim.execute('SELECT combo_id,replicate FROM runs WHERE run_id = {}'.format(run_id)).fetchall()
combo_id = ID[0][0]
replicate = ID[0][1]
CC = curSim.execute('SELECT microbe_carrying_capacity FROM param_combos WHERE combo_id = {}'.format(combo_id)).fetchall()
CC = CC[0][0]
RUN_DIR = os.path.join('runID{0}-c{1}-r{2}'.format(run_id,combo_id,replicate))

DBMATCH_PATH = os.path.join(SCRIPT_PATH,'..','..','isolates',,RUN_DIR,'matches_output.sqlite') # cluster
DBTRI_PATH = os.path.join(SCRIPT_PATH,'..','..','isolates',RUN_DIR,'tripartite-networks_output.sqlite') # cluster
DBPROB_PATH = os.path.join(SCRIPT_PATH,'..','..','isolates',RUN_DIR,'probability-of-emergence_output.sqlite') # cluster
# DBMATCH_PATH = os.path.join('/Volumes','Yadgah',dir,'matches_output.sqlite') # local
# DBTRI_PATH = os.path.join('/Volumes','Yadgah',dir,'tripartite-networks_output.sqlite') # local
# DBPROB_PATH = os.path.join('/Volumes','Yadgah',dir,'probability-of-emergence_output.sqlite') # local

# conMatch = sqlite3.connect(DBMATCH_PATH)
# curMatch = conMatch.cursor()
# conTri = sqlite3.connect(DBTRI_PATH)
# curTri = conTri.cursor()
conProb = sqlite3.connect(DBPROB_PATH)
curProb = conProb.cursor()

# print('SQLite Query: microbial abundance time series data')
# microbe_stacked = pd.read_sql_query("SELECT t,bstrain_id,abundance FROM babundance WHERE run_id = {}".format(run_id), conSim)
#
# print('SQLite Query: viral abundance time series data')
# virus_stacked = pd.read_sql_query("SELECT t,vstrain_id,abundance FROM vabundance WHERE run_id = {}".format(run_id), conSim)
#
# microbe_stacked = microbe_stacked[microbe_stacked.t>=max(virus_stacked.t)]
# microbe_stacked = microbe_stacked.pivot(index='t',columns='bstrain_id',values='abundance')
# virus_stacked = virus_stacked.pivot(index='t',columns='vstrain_id',values='abundance')

freqDisplaced =  pd.read_sql_query("SELECT t,microbial_abundance FROM summary \
WHERE run_id = {0}".format(run_id),conSim).rename(columns={"microbial_abundance": "bTotal"})
vTotal =  pd.read_sql_query("SELECT t,viral_abundance FROM summary \
WHERE run_id = {0}".format(run_id),conSim).rename(columns={"viral_abundance": "vTotal"})
# freqDisplaced =  pd.read_sql_query("SELECT t,viral_abundance FROM summary \
# WHERE run_id = {0}".format(run_id),conSim).rename(columns={"viral_abundance": "freq_displaced"})

freqDisplaced['freq_displaced'] = (CC-np.array(freqDisplaced['bTotal']))/CC
freqDisplaced = freqDisplaced[freqDisplaced['t'] <= max(vTotal.t)]

pal = sns.color_palette("tab20b")
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
        epiSizes.extend(integrate.quad(interpBurnFreq, ti, tj-.005)[0])
        break
    size = integrate.quad(interpBurnFreq, ti, tj)[0]
    epiSizes.append(size)


# CDF

epiSizeData = pd.DataFrame()
counts = []
bin0 = 0.01
binf = 13
binSize = 0.005
bins = np.arange(bin0,binf,binSize)
for i in bins:
    sizes = list(filter(lambda x: x > i, epiSizes))
    counts.append(len(sizes))




def power_law(x, a, b):
    return a*np.power(x, b)



# Fit the dummy power-law data
pars, cov = curve_fit(f=power_law, xdata=bins[:], ydata=counts[:], p0=[0, 0], bounds=(-np.inf, np.inf))
# Get the standard deviations of the parameters (square roots of the # diagonal of the covariance)
stdevs = np.sqrt(np.diag(cov))
# Calculate the residuals
res = y_dummy - power_law(x_dummy, *pars)

epiSizeData['sizes'] = bins[:]
epiSizeData['count'] = counts[:]
fig, ax = plt.subplots(1)
ax.loglog(epiSizeData['sizes'], epiSizeData['count'],"o")
ax.loglog(bins[:],power_law(bins[:],pars[0],pars[1]))
ax.loglog(bins[:],power_law(bins[:],1/2*pars[0],pars[1]-.5))
ax.loglog(bins[:],power_law(bins[:],pars[0],pars[1]+.5))


epiSizeData = pd.DataFrame()
epiSizeData['sizes'] = list(np.around(np.array(epiSizes), 1))
epiSizeData['count'] = epiSizeData.groupby('sizes')['sizes'].transform('count')
epiSizeData = epiSizeData.sort_values(by=['sizes'])
fig, ax = plt.subplots(1)
ax.loglog(epiSizeData['sizes'], epiSizeData['count'],"o")

ts = np.arange(min(burnFreq['t']), max(burnFreq['t']), 0.01)
fig, ax = plt.subplots(1)
ax.plot(ts, interpBurnFreq(ts))


for k in range(0,len(peakWidths[1:][1])):
    print(k)

for (i,j) in zip(peakWidths[1:][1],peakWidths[1:][2]):
    print(i)
    print(j)



for k in range(0,len(peakWidths[1:][1])):
    ti = peakWidths[1:][1][k] + min(burnFreq['t'])
    tj = peakWidths[1:][2][k] + min(burnFreq['t'])
    print(ti,tj)



print('Complete!')
