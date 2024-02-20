#!/usr/bin/env python3

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import math
import sqlite3
import matplotlib.cm as cm
import matplotlib.colors as mc
from matplotlib.colors import rgb2hex
import plotly.graph_objects as go

run = 'runID3297-c66-r47'

resolve = 500
imgTypes = ["pdf"]
figxy = (15,10) # setting for tree abundance figure
colorpalSpacers = 'tab20b' # Color palette for spacer in networks: use discrete color palette
sSpacing = 8
vSpacing = 8
bSpacing = 6

SCRIPT_PATH = os.path.abspath(os.path.dirname(__file__)) # cluster

# SQLITE paths, beware of manipulating these. Directories are structured accordingly
DBMATCH_PATH = os.path.join(SCRIPT_PATH,'..','generated-data','isolates',run,'matches_output.sqlite') 
DBTREE_PATH = os.path.join(SCRIPT_PATH,'..','generated-data','isolates',run,'trees_output.sqlite') 
DBTRI_PATH = os.path.join(SCRIPT_PATH,'..','generated-data','isolates',run,'tripartite-networks_output.sqlite')
DBSIM_PATH = os.path.join(SCRIPT_PATH,'..','generated-data','isolates',run,'{}.sqlite'.format(run))
PLOT_PATH = os.path.join(SCRIPT_PATH)
# local
dir = 'crispr-sweep-7-2-2022/isolates/runID3297-c66-r47'
DBPROB_PATH = os.path.join('/Volumes','Yadgah',dir,'emergence-lambert-root_output.sqlite') # local
DBMATCH_PATH = os.path.join('/Volumes','Yadgah',dir,'matches_output.sqlite') # local
DBTRI_PATH = os.path.join('/Volumes','Yadgah',dir,'tripartite-networks_output.sqlite') # local
DBSIM_PATH = os.path.join('/Volumes','Yadgah',dir,'{}.sqlite'.format(run))
PLOT_PATH = os.path.join('/Volumes','Yadgah')


conSim = sqlite3.connect(DBSIM_PATH)
curSim = conSim.cursor()
run_id = curSim.execute('SELECT DISTINCT run_id FROM summary').fetchall()
run_id = run_id[0][0]
ID = curSim.execute('SELECT combo_id,replicate FROM runs WHERE run_id = {}'.format(run_id)).fetchall()
combo_id = ID[0][0]
replicate = ID[0][1]
conMatch = sqlite3.connect(DBMATCH_PATH)
curMatch = conMatch.cursor()
conTri = sqlite3.connect(DBTRI_PATH)
curTri = conTri.cursor()

times = [275, 475, 750]
html = False

def tripartiteGraph(times, run_id, sSpacing, vSpacing, bSpacing, html, colorpalSpacers, imgTypes, DBSIM_PATH, DBMATCH_PATH, DBTRI_PATH, PLOT_PATH):
    conSim = sqlite3.connect(DBSIM_PATH)
    curSim = conSim.cursor()
    conMatch = sqlite3.connect(DBMATCH_PATH)
    curMatch = conMatch.cursor()
    ID = curSim.execute(
        'SELECT combo_id,replicate FROM runs WHERE run_id = {}'.format(run_id)).fetchall()
    combo_id = ID[0][0]
    replicate = ID[0][1]
    CC = curSim.execute(
        'SELECT microbe_carrying_capacity FROM param_combos WHERE combo_id = {}'.format(combo_id)).fetchall()
    CC = CC[0][0]
    RUN_DIR = os.path.join(
        'runID{0}-c{1}-r{2}'.format(run_id, combo_id, replicate))
    conTri = sqlite3.connect(DBTRI_PATH)
    curTri = conTri.cursor()
    abundNormV = 1
    abundNormB = 1
    for time in times:
        # Single Match Tripartite Network
        #
        tripartiteNetwork = pd.read_sql_query("SELECT bstrain_id, vstrain_id, time_specific_match_id \
                                FROM bstrain_to_vstrain_matches WHERE match_length = 1 \
                                    AND t = {0}".format(time), conMatch)
        spacerMatches = pd.read_sql_query("SELECT time_specific_match_id, spacer_id\
                                FROM matches_spacers WHERE t = {0}".format(time), conMatch)
        tripartiteNetwork = tripartiteNetwork.merge(spacerMatches, on=['time_specific_match_id'])\
            .drop(columns=['time_specific_match_id'])
        bstrainmatchIDs = pd.read_sql_query(
            "SELECT match_id, bstrain_id \
                FROM bmatches WHERE t = {0}".format(time),
            conTri).rename(columns={"match_id": "bmatch_id"})
        vstrainmatchIDs = pd.read_sql_query(
            "SELECT match_id, vstrain_id \
                FROM vmatches WHERE t = {0}".format(time),
            conTri).rename(columns={"match_id": "vmatch_id"})
        tripartiteNetwork = tripartiteNetwork.merge(
            bstrainmatchIDs, on=['bstrain_id'])\
            .merge(
            vstrainmatchIDs, on=['vstrain_id'])\
            .drop(columns=['bstrain_id', 'vstrain_id'])\
            .drop_duplicates()
        bstrainmatchIDs = bstrainmatchIDs.merge(
            tripartiteNetwork[['bmatch_id']].drop_duplicates(), on=['bmatch_id'])
        vstrainmatchIDs = vstrainmatchIDs.merge(
            tripartiteNetwork[['vmatch_id']].drop_duplicates(), on=['vmatch_id'])
        ###
        bstrainAbund = pd.read_sql_query("SELECT bstrain_id, abundance FROM babundance  \
            WHERE run_id = {0} AND t = {1} AND bstrain_id in ({2})".format(run_id, time, ', '
                                                                           .join(map(str, bstrainmatchIDs['bstrain_id']))), conSim)\
            .rename(columns={"abundance": "bstrainAbund"})
        bmatchAbund = pd.read_sql_query("SELECT match_id,babundance\
                FROM bmatches_abundances  \
                WHERE t = {0} AND match_id in ({1}) \
                ORDER BY babundance".format(time, ', '
                                            .join(map(str, tripartiteNetwork['bmatch_id'].unique()))), conTri)\
            .rename(columns={"match_id": "bmatch_id"})\
            .rename(columns={"babundance": "bmatchAbund"})
        bstrainmatchIDs = bstrainmatchIDs.merge(bstrainAbund, on=['bstrain_id'])\
            .merge(bmatchAbund, on=['bmatch_id']).sort_values(by=['bmatchAbund', 'bstrainAbund'])
        vstrainAbund = pd.read_sql_query("SELECT vstrain_id, abundance FROM vabundance  \
            WHERE run_id = {0} AND t = {1} AND vstrain_id in ({2})".format(run_id, time, ', '
                                                                           .join(map(str, vstrainmatchIDs['vstrain_id']))), conSim)\
            .rename(columns={"abundance": "vstrainAbund"})
        vmatchAbund = pd.read_sql_query("SELECT match_id,vabundance,bsusceptible \
                FROM vmatches_abundances  \
                WHERE t = {0} AND match_id in ({1}) \
                ORDER BY vabundance".format(time, ', '
                                            .join(map(str, tripartiteNetwork['vmatch_id'].unique()))), conTri)\
            .rename(columns={"match_id": "vmatch_id"})\
            .rename(columns={"vabundance": "vmatchAbund"})
        spacerAbund = tripartiteNetwork.drop(columns=['vmatch_id']).drop_duplicates()\
            .merge(bstrainmatchIDs.drop(columns=['bstrain_id', 'bstrainAbund']),
                   on=['bmatch_id']).groupby(['spacer_id']).\
            agg(bmatchAbund=('bmatchAbund', 'sum')).reset_index()
        if abundNormV < max(vmatchAbund['vmatchAbund']):
            abundNormV = max(vmatchAbund['vmatchAbund'])
        if abundNormB < max(spacerAbund['bmatchAbund']):
            abundNormB = max(spacerAbund['bmatchAbund'])
    for time in times:
        # Single Match Tripartite Network
        #
        tripartiteNetwork = pd.read_sql_query("SELECT bstrain_id, vstrain_id, time_specific_match_id \
                                FROM bstrain_to_vstrain_matches WHERE match_length = 1 \
                                    AND t = {0}".format(time), conMatch)
        spacerMatches = pd.read_sql_query("SELECT time_specific_match_id, spacer_id\
                                FROM matches_spacers WHERE t = {0}".format(time), conMatch)
        tripartiteNetwork = tripartiteNetwork.merge(spacerMatches, on=['time_specific_match_id'])\
            .drop(columns=['time_specific_match_id'])
        bstrainmatchIDs = pd.read_sql_query(
            "SELECT match_id, bstrain_id \
                FROM bmatches WHERE t = {0}".format(time),
            conTri).rename(columns={"match_id": "bmatch_id"})
        vstrainmatchIDs = pd.read_sql_query(
            "SELECT match_id, vstrain_id \
                FROM vmatches WHERE t = {0}".format(time),
            conTri).rename(columns={"match_id": "vmatch_id"})
        tripartiteNetwork = tripartiteNetwork.merge(
            bstrainmatchIDs, on=['bstrain_id'])\
            .merge(
            vstrainmatchIDs, on=['vstrain_id'])\
            .drop(columns=['bstrain_id', 'vstrain_id'])\
            .drop_duplicates()
        bstrainmatchIDs = bstrainmatchIDs.merge(
            tripartiteNetwork[['bmatch_id']].drop_duplicates(), on=['bmatch_id'])
        vstrainmatchIDs = vstrainmatchIDs.merge(
            tripartiteNetwork[['vmatch_id']].drop_duplicates(), on=['vmatch_id'])
        ###
        bstrainAbund = pd.read_sql_query("SELECT bstrain_id, abundance FROM babundance  \
            WHERE run_id = {0} AND t = {1} AND bstrain_id in ({2})".format(run_id, time, ', '
                                                                           .join(map(str, bstrainmatchIDs['bstrain_id']))), conSim)\
            .rename(columns={"abundance": "bstrainAbund"})
        bmatchAbund = pd.read_sql_query("SELECT match_id,babundance\
                FROM bmatches_abundances  \
                WHERE t = {0} AND match_id in ({1}) \
                ORDER BY babundance".format(time, ', '
                                            .join(map(str, tripartiteNetwork['bmatch_id'].unique()))), conTri)\
            .rename(columns={"match_id": "bmatch_id"})\
            .rename(columns={"babundance": "bmatchAbund"})
        bstrainmatchIDs = bstrainmatchIDs.merge(bstrainAbund, on=['bstrain_id'])\
            .merge(bmatchAbund, on=['bmatch_id']).sort_values(by=['bmatchAbund', 'bstrainAbund'])
        vstrainAbund = pd.read_sql_query("SELECT vstrain_id, abundance FROM vabundance  \
            WHERE run_id = {0} AND t = {1} AND vstrain_id in ({2})".format(run_id, time, ', '
                                                                           .join(map(str, vstrainmatchIDs['vstrain_id']))), conSim)\
            .rename(columns={"abundance": "vstrainAbund"})
        vmatchAbund = pd.read_sql_query("SELECT match_id,vabundance,bsusceptible \
                FROM vmatches_abundances  \
                WHERE t = {0} AND match_id in ({1}) \
                ORDER BY vabundance".format(time, ', '
                                            .join(map(str, tripartiteNetwork['vmatch_id'].unique()))), conTri)\
            .rename(columns={"match_id": "vmatch_id"})\
            .rename(columns={"vabundance": "vmatchAbund"})
        vstrainmatchIDs = vstrainmatchIDs.merge(vstrainAbund, on=['vstrain_id'])\
            .merge(vmatchAbund[['vmatch_id', 'vmatchAbund']], on=['vmatch_id'])\
            .sort_values(by=['vmatchAbund', 'vstrainAbund'])
        spacerIDs = list(tripartiteNetwork['spacer_id'].unique())
        spacerIDs.sort()
        spacerDict = {i: j for (i, j) in zip(
            spacerIDs, list(range(0, len(spacerIDs))))}
        print(spacerDict)
        bPhenos = pd.read_sql_query("SELECT bstrain_id, spacer_id \
            FROM bspacers WHERE run_id = {0} AND bstrain_id in ({1})"
                                    .format(run_id, ', '.join(map(str, bstrainmatchIDs.bstrain_id.values))), conSim)
        spacersOrdered = []  # order spacers by time acquired
        for bstrainID in sorted(bstrainmatchIDs['bstrain_id'].values):
            phenotype = bPhenos[bPhenos.bstrain_id ==
                                bstrainID].spacer_id.values
            phenotype = list(set(phenotype).intersection(set(spacerIDs)))
            print('phenotype is {0}'.format(phenotype))
            newSpacer = list(set([spacerDict[i]
                             for i in phenotype]) - set(spacersOrdered))
            print(newSpacer)
            if len(newSpacer) == 0:
                continue
            spacersOrdered.extend(newSpacer)
            print('spacersOrdered is {0}'.format(spacersOrdered))
        #
        tripartiteNetwork['spacer_id'] = [spacerDict[i]
                                          for i in tripartiteNetwork['spacer_id']]  # take note of this!!
        Mcmap = cm.get_cmap('{}'.format(colorpalSpacers))
        if len(spacerIDs) > 1:
            spacerCMAP = mc.LinearSegmentedColormap.from_list(
                'spacers', [Mcmap(i) for i in range(len(spacerIDs))], len(spacerIDs))
        else:
            spacerCMAP = Mcmap
        #
        # print('spacerIDs: {}'.format(spacerIDs))
        numV = len(tripartiteNetwork['vmatch_id'].unique())
        numS = len(tripartiteNetwork['spacer_id'].unique())
        numB = len(tripartiteNetwork['bmatch_id'].unique())
        ##
        scale = 20
        ##
        if (numV % 2) == 0:
            vOffset = -numV/2*vSpacing
        else:
            vOffset = -(numV-1)/2*vSpacing
        #
        if (numS % 2) == 0:
            sOffset = -numS/2*sSpacing
        else:
            sOffset = -(numS-1)/2*sSpacing
        #
        if (numB % 2) == 0:
            bOffset = -numB/2*bSpacing
        else:
            bOffset = -(numB-1)/2*bSpacing
        #
        vNodeDict = {}
        vNodeColor = []
        vNodeSize = []
        vNodeText = []
        vNodePositionsX = []
        vNodePositionsY = numV * [1]
        vstrainsOrdered = list(vstrainmatchIDs[['vmatch_id', 'vmatchAbund']]
                               .sort_values(by=['vmatchAbund'])['vmatch_id'].unique())  # bin values are sorted here
        for vID in vstrainsOrdered:
            vNodeColor.append('#888')
            vNodeSize.append(
                scale*vmatchAbund[vmatchAbund.vmatch_id == vID]['vmatchAbund'].values[0]/abundNormV)
            # vNodeSize.append(5*(np.log(*vmatchAbund[vmatchAbund.vmatch_id==vID]['vmatchAbund'])+1))
            vNodeText.append('{}'.format(vID))
            vNodeDict[vID] = [vOffset, vNodePositionsY[0]]
            # y = semicircle(vOffset,numS*sSpacing/2,0,.5)
            # vNodePositionsY.append(y)
            # vNodeDict[vmatchID] = [vOffset,y]
            vNodePositionsX.append(vOffset)
            vOffset += vSpacing
        sNodeDict = {}
        sNodeColor = {}
        sNodeSize = {}
        sNodeText = numS * [None]
        sNodePositionsX = {}
        sNodePositionsY = 0
        #
        spacerAbund = tripartiteNetwork.drop(columns=['vmatch_id']).drop_duplicates()\
            .merge(bstrainmatchIDs.drop(columns=['bstrain_id', 'bstrainAbund']),
                   on=['bmatch_id']).groupby(['spacer_id']).\
            agg(bmatchAbund=('bmatchAbund', 'sum')).reset_index()
        for sID in spacerAbund.sort_values(by=['spacer_id'])['spacer_id']:
            sNodeColor[sID] = rgb2hex(spacerCMAP(sID))
            print('Spacer node color: {0} for spacerID {1}'.format(
                rgb2hex(spacerCMAP(sID)), sID))
        spacersOrdered2 = [*spacersOrdered[1:], spacersOrdered[0]]
        for (sID, sID2) in zip(spacersOrdered, spacersOrdered2):
            sNodeSize[sID] = scale*spacerAbund[spacerAbund.spacer_id ==
                                               sID]['bmatchAbund'].values[0]/abundNormB
            # sNodeSize[sID] = 2
            # if sID == spacersOrdered[0]:
            #     sNodeSize[sID] = 50
            sNodeDict[sID] = [sOffset, sNodePositionsY]
            sNodePositionsX[sID] = sOffset
            sOffset += sSpacing
        #
        bNodeDict = {}
        bNodeColor = []
        bNodeSize = []
        bNodeText = []
        bNodePositionsX = []
        bNodePositionsY = numB * [-1]
        #
        bstrainsOrdered = list(bstrainmatchIDs[['bmatch_id', 'bmatchAbund']]
                               .sort_values(by=['bmatchAbund'])['bmatch_id'].unique())  # bin values are sorted here
        for bID in bstrainsOrdered:
            bNodeColor.append('#888')
            bNodeSize.append(
                scale*bmatchAbund[bmatchAbund.bmatch_id == bID]['bmatchAbund'].values[0]/abundNormB)
            bNodeText.append('{}'.format(bID))
            bNodeDict[bID] = [bOffset, bNodePositionsY[0]]
            bNodePositionsX.append(bOffset)
            bOffset += bSpacing
        node_trace = []
        node_trace.append(go.Scatter(
            x=vNodePositionsX + bNodePositionsX,
            y=vNodePositionsY + bNodePositionsY,
            mode='markers',
            hoverinfo='text',
            showlegend=False,
            marker=dict(
                showscale=False,
                color=vNodeColor+bNodeColor,
                size=vNodeSize+bNodeSize,
                line_width=.5)
        ))
        spacers = pd.read_sql_query(
            "SELECT spacer_id \
            FROM single_match_tripartite_networks WHERE t = {0}".format(time),
            conTri)
        spacers = sorted(spacers['spacer_id'].unique())
        for i in spacers[::-1]:
            node_trace.append(go.Scatter(
                x=[sNodePositionsX[spacerDict[i]]],
                y=[sNodePositionsY],
                mode='markers',
                hoverinfo='text',
                marker=dict(
                    showscale=False,
                    color=sNodeColor[spacerDict[i]],
                    size=sNodeSize[spacerDict[i]],
                    line_width=.5),
                showlegend=False,
                legendgroup='{}'.format(i),
            ))
        for i in spacers[::-1]:
            node_trace.append(go.Scatter(
                x=[None],
                y=[None],
                mode='markers',
                hoverinfo='text',
                marker=dict(
                    color=sNodeColor[spacerDict[i]],
                    size=10,
                    line_width=.5),
                legendgroup='{}'.format(i), showlegend=True, name='{}'.format(i),
            ))
        vtext_trace = go.Scatter(
            x=vNodePositionsX,
            y=vNodePositionsY,
            mode='text',
            text=vNodeText,
            textposition="top center",
            showlegend=False,
            textfont_size=4.5)
        btext_trace = go.Scatter(
            x=bNodePositionsX,
            y=bNodePositionsY,
            mode='text',
            text=bNodeText,
            textposition="bottom center",
            showlegend=False,
            textfont_size=4.5)
        edge_x = []
        edge_y = []
        edge_trace = []
        for sID in tripartiteNetwork['spacer_id'].unique():
            for edge in zip(tripartiteNetwork[tripartiteNetwork.spacer_id == sID]['vmatch_id'],
                            tripartiteNetwork[tripartiteNetwork.spacer_id == sID]['spacer_id']):
                x0 = vNodeDict[edge[0]][0]
                y0 = vNodeDict[edge[0]][1]
                x1 = sNodeDict[edge[1]][0]
                y1 = sNodeDict[edge[1]][1]
                edge_x.append(x0)
                edge_x.append(x1)
                edge_x.append(None)
                edge_y.append(y0)
                edge_y.append(y1)
                edge_y.append(None)
            for edge in zip(tripartiteNetwork[tripartiteNetwork.spacer_id == sID]['spacer_id'],
                            tripartiteNetwork[tripartiteNetwork.spacer_id == sID]['bmatch_id']):
                x0 = sNodeDict[edge[0]][0]
                y0 = sNodeDict[edge[0]][1]
                x1 = bNodeDict[edge[1]][0]
                y1 = bNodeDict[edge[1]][1]
                edge_x.append(x0)
                edge_x.append(x1)
                edge_x.append(None)
                edge_y.append(y0)
                edge_y.append(y1)
                edge_y.append(None)
            edge_trace.append(go.Scatter(
                x=edge_x, y=edge_y,
                line=dict(width=0.5, color=rgb2hex(spacerCMAP(edge[0]))),
                hoverinfo='none',
                showlegend=False,
                mode='lines'))
            edge_x = []
            edge_y = []
        fig = go.Figure(data=[*edge_trace, *node_trace, vtext_trace, btext_trace],
                        layout=go.Layout(
                        # paper_bgcolor='rgba(0,0,0,0)',
                        plot_bgcolor='rgba(0, 0, 0, 0)',
                        title='<br>Tripartite graph: runID{0}-c{1}-r{2} @ time = {3}'
                        .format(run_id, combo_id, replicate, time),
                        titlefont_size=16,
                        hovermode='closest',
                        margin=dict(b=20, l=5, r=5, t=40),
                        xaxis=dict(showgrid=False, zeroline=False,
                                   showticklabels=False),
                        yaxis=dict(showgrid=False, zeroline=False, showticklabels=False))
                        )
        for imgType in imgTypes:
            fig.write_image(os.path.join(PLOT_PATH, 'tripartite-graph-time{0}.{1}'
                                         .format(int(np.floor(time)), imgType)))
        if html:
            fig.show()
