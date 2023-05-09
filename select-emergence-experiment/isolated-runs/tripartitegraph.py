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

def tripartiteGraph(time,run_id,sSpacing,vSpacing,bSpacing,html,colorpalSpacers,imgTypes,DBSIM_PATH,DBMATCH_PATH,DBTRI_PATH,PLOT_PATH):
    conSim = sqlite3.connect(DBSIM_PATH)
    curSim = conSim.cursor()
    conMatch = sqlite3.connect(DBMATCH_PATH)
    curMatch = conMatch.cursor()
    ID = curSim.execute('SELECT combo_id,replicate FROM runs WHERE run_id = {}'.format(run_id)).fetchall()
    combo_id = ID[0][0]
    replicate = ID[0][1]
    CC = curSim.execute('SELECT microbe_carrying_capacity FROM param_combos WHERE combo_id = {}'.format(combo_id)).fetchall()
    CC = CC[0][0]
    RUN_DIR = os.path.join('runID{0}-c{1}-r{2}'.format(run_id,combo_id,replicate))
    conTri = sqlite3.connect(DBTRI_PATH)
    curTri = conTri.cursor()
    tripartiteNetwork = pd.read_sql_query(
        "SELECT vmatch_id, spacer_id, bmatch_id \
        FROM single_match_tripartite_networks WHERE t = {0}".format(time),
        conTri)
    bstrainmatchIDs = pd.read_sql_query(
        "SELECT match_id, bstrain_id \
        FROM bmatches WHERE t = {0}".format(time),
        conTri).rename(columns={"match_id": "bmatch_id"})
    bstrainAbund = pd.read_sql_query("SELECT bstrain_id,abundance FROM babundance  \
        WHERE run_id = {0} AND t = {1}".format(run_id,time), conSim).rename(columns={"abundance": "bstrainAbund"})
    bmatchAbund = pd.read_sql_query("SELECT match_id,babundance \
            FROM bmatches_abundances  \
            WHERE t = {0} AND match_id in ({1})".format(time,', '
            .join(map(str,tripartiteNetwork['bmatch_id'].unique()))), conTri)\
            .rename(columns={"match_id": "bmatch_id"}).rename(columns={"babundance": "bmatchAbund"})
    tripartiteNetwork = tripartiteNetwork.merge(bstrainmatchIDs,on=['bmatch_id']).drop(columns='bmatch_id')
    bstrainmatchIDs = bstrainmatchIDs.merge(bstrainAbund,on=['bstrain_id'])\
                    .merge(bmatchAbund,on=['bmatch_id']).sort_values(by=['bmatchAbund', 'bstrainAbund'])
    vmatchAbund = pd.read_sql_query("SELECT match_id,vabundance,bsusceptible \
            FROM vmatches_abundances  \
            WHERE t = {0} AND match_id in ({1}) \
            ORDER BY vabundance".format(time,', '
            .join(map(str,tripartiteNetwork['vmatch_id'].unique()))), conTri)\
            .rename(columns={"match_id": "vmatch_id"})\
            .rename(columns={"vabundance": "vmatchAbund",\
            "bsusceptible": "EscEventProb"})
    bTotal = curSim.execute("SELECT microbial_abundance FROM summary \
    WHERE run_id = {0} AND t = {1}".format(run_id,time)).fetchall()
    bTotal = bTotal[0][0]
    vTotal = curSim.execute("SELECT viral_abundance FROM summary \
    WHERE run_id = {0} AND t = {1}".format(run_id,time)).fetchall()
    vTotal = vTotal[0][0]
    bstrain0vstrain = pd.read_sql_query(
    "SELECT bstrain_id, vstrain_id \
    FROM bstrain_to_vstrain_0matches WHERE t = {0}".format(time), conMatch)
    vAbunds = pd.read_sql_query("SELECT vstrain_id,abundance FROM vabundance \
    WHERE run_id = {0} AND t = {1}".format(run_id,time), conSim)
    bstrain0vstrain = pd.read_sql_query(
    "SELECT bstrain_id, vstrain_id \
    FROM bstrain_to_vstrain_0matches WHERE t = {0}".format(time), conMatch)
    vAbunds = pd.read_sql_query("SELECT vstrain_id,abundance FROM vabundance \
    WHERE run_id = {0} AND t = {1}".format(run_id,time), conSim)
    bstrain0vstrain = bstrain0vstrain.merge(vAbunds,on=['vstrain_id'])\
    .groupby(['bstrain_id']).agg(CollapseProb=('abundance','sum')).reset_index()
    bstrain0vstrain['CollapseProb'] = 1/vTotal*bstrain0vstrain['CollapseProb']
    bstrain0vstrain = bstrainmatchIDs.merge(bstrain0vstrain,on=['bstrain_id'])
    vprobs = bstrain0vstrain['CollapseProb'][:]
    bstrain0vstrain['CollapseProb'] = np.array(bstrain0vstrain['CollapseProb'])\
                                    *1/bTotal*np.array(bstrain0vstrain['bstrainAbund'])
    bstrain0vstrain['MatchCollapseProb'] = vprobs*1/bTotal*np.array(bstrain0vstrain['bmatchAbund'])
    bstrain0vstrain = bstrain0vstrain.sort_values(by=['MatchCollapseProb', 'CollapseProb']) # bin values are sorted here
    vmatchAbund['EscEventProb']= list(1/bTotal * np.array(vmatchAbund['EscEventProb'])\
    *1/vTotal* np.array(vmatchAbund['vmatchAbund']))
    spacerIDs = list(tripartiteNetwork['spacer_id'].unique())
    spacerIDs.sort()
    spacerDict = {i:j for (i,j) in zip(spacerIDs,list(range(0,len(spacerIDs))))}
    print(spacerDict)
    bPhenos = pd.read_sql_query("SELECT bstrain_id, spacer_id \
        FROM bspacers WHERE run_id = {0} AND bstrain_id in ({1})"\
        .format(run_id,', '.join(map(str,bstrainmatchIDs.bstrain_id.values))), conSim)
    spacersOrdered = [] # order spacers by time acquired
    for bstrainID in sorted(bstrainmatchIDs['bstrain_id'].values):
        phenotype = bPhenos[bPhenos.bstrain_id == bstrainID].spacer_id.values
        phenotype = list(set(phenotype).intersection(set(spacerIDs)))
        print('phenotype is {0}'.format(phenotype))
        newSpacer = list(set([spacerDict[i] for i in phenotype]) - set(spacersOrdered))
        print(newSpacer)
        if len(newSpacer) == 0:
            continue
        spacersOrdered.extend(newSpacer)
        print('spacersOrdered is {0}'.format(spacersOrdered))
    print('fin')
    tripartiteNetwork['spacer_id'] = [spacerDict[i] for i in tripartiteNetwork['spacer_id']] # take note of this!!
    Mcmap = cm.get_cmap('{}'.format(colorpalSpacers))
    if len(spacerIDs) > 1:
        spacerCMAP = mc.LinearSegmentedColormap.from_list(
            'spacers', [Mcmap(i) for i in range(len(spacerIDs))], len(spacerIDs))
    else:
        spacerCMAP = Mcmap
    print('spacerIDs: {}'.format(spacerIDs))
    numV = len(tripartiteNetwork['vmatch_id'].unique())
    numS = len(tripartiteNetwork['spacer_id'].unique())
    numB = len(tripartiteNetwork['bstrain_id'].unique())
    if (numV % 2) == 0:
        vOffset = -numV/2*vSpacing
    else:
        vOffset = -(numV-1)/2*vSpacing
    if (numS % 2) == 0:
        sOffset = -numS/2*sSpacing
    else:
        sOffset = -(numS-1)/2*sSpacing
    if (numB % 2) == 0:
        bOffset = -numB/2*bSpacing
    else:
        bOffset = -(numB-1)/2*bSpacing
    vNodeDict = {}
    vNodeColor = []
    vNodeSize = []
    vNodeText = []
    vNodePositionsX = []
    vNodePositionsY = numV * [1]
    for vmatchID in vmatchAbund.sort_values(by=['EscEventProb'])['vmatch_id']:
        vNodeColor.append('#888')
        vNodeSize.append(np.log(*vmatchAbund[vmatchAbund.vmatch_id==vmatchID]['vmatchAbund'])+1)
        vNodeText.append('{}'.format(vmatchID))
        vNodeDict[vmatchID] = [vOffset,vNodePositionsY[0]]
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
    spacerAbund = tripartiteNetwork.drop(columns=['vmatch_id']).drop_duplicates()\
                    .merge(bstrainmatchIDs.drop(columns = ['bmatch_id','bmatchAbund']),\
                    on=['bstrain_id']).groupby(['spacer_id']).\
                    agg(bstrainAbund=('bstrainAbund','sum')).reset_index()
    spacerAbund = vmatchAbund.merge(tripartiteNetwork.drop(columns=['bstrain_id'])\
                    .drop_duplicates(),on='vmatch_id')\
                    .drop(columns=['EscEventProb'])\
                    .groupby(['spacer_id'])\
                    .agg(vmatchAbund=('vmatchAbund','sum')).reset_index()\
                    .merge(spacerAbund,on=['spacer_id'])
    print(tripartiteNetwork['spacer_id'].unique())
    print(spacerAbund)
    for sID in spacerAbund.sort_values(by=['spacer_id'])['spacer_id']:
        sNodeColor[sID] = rgb2hex(spacerCMAP(sID))
        print('Spacer node color: {0} for spacerID {1}'.format(rgb2hex(spacerCMAP(sID)),sID))
    for sID in spacersOrdered:
        sNodeSize[sID] = 20*spacerAbund[spacerAbund.spacer_id==sID]['bstrainAbund'].values[0]/CC
        # sNodeSize[sID] = np.log(*spacerAbund[spacerAbund.spacer_id==sID]['bstrainAbund'])+1
        print('frequency is {}'.format(spacerAbund[spacerAbund.spacer_id==sID]['bstrainAbund'].values[0]/CC))
        sNodeDict[sID] = [sOffset,sNodePositionsY]
        sNodePositionsX[sID] = sOffset
        sOffset += sSpacing
    bNodeDict = {}
    bNodeColor = []
    bNodeSize = []
    bNodeText = []
    bNodePositionsX = []
    bNodePositionsY = numB * [-1]
    # bNodePositionsY = []
    # def semicircle(x,r, h, k):
    #     # use numpy for array solving of the semicircle equation
    #     y = k + -1*np.sqrt(r**2 - (x - h)**2)
    #     return y
    bstrainsOrdered = list(set(bstrainmatchIDs['bstrain_id']) - set(bstrain0vstrain['bstrain_id']))
    bstrainsOrdered = list(bstrainmatchIDs[bstrainmatchIDs.bstrain_id.isin(bstrainsOrdered)]\
                    .sort_values(by=['bmatchAbund', 'bstrainAbund'])['bstrain_id']) # bin values are sorted here
    bstrainsOrdered.extend(bstrain0vstrain['bstrain_id'])
    for bstrainID in bstrainsOrdered:
        bNodeColor.append('#888')
        bNodeSize.append(bstrainmatchIDs[bstrainmatchIDs.bstrain_id==bstrainID]['bstrainAbund'].values[0]/CC)
        # bNodeSize.append(np.log(*bstrainmatchIDs[bstrainmatchIDs.bstrain_id==bstrainID]['bstrainAbund'])+1)
        print('frequency is {}'.format(bstrainmatchIDs[bstrainmatchIDs.bstrain_id==bstrainID]['bstrainAbund'].values[0]/CC))
        bNodeText.append('{}'.format(bstrainID))
        bNodeDict[bstrainID] = [bOffset,bNodePositionsY[0]]
        bNodePositionsX.append(bOffset)
        bOffset += bSpacing
    vNodeSize = list(20/max(vNodeSize)*np.array(vNodeSize)) # rescale node sizes
    # sNodeSize = list(20/max(sNodeSize)*np.array(sNodeSize)) # rescale node sizes
    # bNodeSize = list(20/max(bNodeSize)*np.array(bNodeSize)) # rescale node sizes
    # sNodeSize = list(20*np.array(sNodeSize)) # rescale node sizes
    bNodeSize = list(20*np.array(bNodeSize)) # rescale node sizes
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
        text = vNodeText,
        textposition="top center",
        showlegend=False,
        textfont_size=6)
    btext_trace = go.Scatter(
        x=bNodePositionsX,
        y=bNodePositionsY,
        mode='text',
        text = bNodeText,
        textposition="bottom center",
        showlegend=False,
        textfont_size=6)
    edge_x = []
    edge_y = []
    edge_trace = []
    for sID in tripartiteNetwork['spacer_id'].unique():
        for edge in zip(tripartiteNetwork[tripartiteNetwork.spacer_id == sID]['vmatch_id'],\
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
        for edge in zip(tripartiteNetwork[tripartiteNetwork.spacer_id == sID]['spacer_id'],\
                        tripartiteNetwork[tripartiteNetwork.spacer_id == sID]['bstrain_id']):
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
                    title='<br>Tripartite graph: runID{0}-c{1}-r{2} @ time = {3}'
                    .format(run_id,combo_id,replicate,time),
                    titlefont_size=16,
                    hovermode='closest',
                    margin=dict(b=20,l=5,r=5,t=40),
                    xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
                    yaxis=dict(showgrid=False, zeroline=False, showticklabels=False))
                    )
    for imgType in imgTypes:
        fig.write_image(os.path.join(PLOT_PATH,'tripartite-graph-time{0}.{1}'\
                            .format(int(np.floor(time)),imgType)))
    if html:
        fig.show()
