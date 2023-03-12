import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import sqlite3

def treetimes(run_id,curSim,conTri,axesV,axesB,times,hyperAnalyze,\
vAbunds,bAbunds,vKeepTreeStrainsDF,bKeepTreeStrainsDF):
    for time in times:
        bstrainmatchIDs = pd.read_sql_query(
            "SELECT match_id, bstrain_id \
            FROM bmatches WHERE t = {0}".format(time),
            conTri).rename(columns={"match_id": "bmatch_id"})
        vstrainmatchIDs = pd.read_sql_query(
            "SELECT match_id, vstrain_id \
            FROM vmatches WHERE t = {0}".format(time),
            conTri).rename(columns={"match_id": "vmatch_id"})
    # This applies info of interest on individual branches of viral tree
        axesV[0].axvline(x=time, color='black', ls=':', lw=.7)
        axesV[1].axvline(x=time, color='black', ls=':', lw=.7)
        if hyperAnalyze:
            for newtreestrainID in vAbunds[vAbunds['t']==time]['tree_vstrain_id']:
                strainID = vKeepTreeStrainsDF[vKeepTreeStrainsDF['new_tree_vstrain_id'] \
                                == newtreestrainID]['vstrain_id'].values[0]
                if strainID in vstrainmatchIDs['vstrain_id'].values:
                    matchID = vstrainmatchIDs[vstrainmatchIDs.vstrain_id\
                    == strainID]['vmatch_id'].values[0]
                else:
                    matchID = 0
                pID = [id[0] for id in curSim.execute("SELECT parent_vstrain_id FROM vstrains \
                    WHERE run_id = {0} AND vstrain_id = {1}"\
                    .format(run_id,strainID))][0]
                pSpacers = [id[0] for id in curSim.execute("SELECT spacer_id FROM vpspacers \
                            WHERE run_id = {0} \
                            AND vstrain_id = {1}".format(run_id,pID))]
                dSpacers = [id[0] for id in curSim.execute("SELECT spacer_id FROM vpspacers \
                            WHERE run_id = {0} \
                            AND vstrain_id = {1}".format(run_id,strainID))]
                frm = list(set(pSpacers) - set(dSpacers))
                to = list(set(dSpacers) - set(pSpacers))
                bID = [id[0] for id in curSim.execute("SELECT infected_bstrain_id \
                    FROM vstrains \
                    WHERE run_id = {0}\
                    AND vstrain_id = {1}".format(run_id,strainID))][0]
                # axesV[1].text(time, treestrain,\
                # 'mID:{0}, vID:{1}, bID:{2}::[{3}]->[{4}]'\
                # .format(mID,treevstrains[treestrain],bID,', '.join(map(str,frm)),\
                # ', '.join(map(str,to))),\
                # fontsize=3)#, fontweight='bold')
                axesV[1].text(time, newtreestrainID,\
                '{0},{1},bID:{2}::[{3}]->[{4}]'\
                .format(strainID,matchID,bID,', '.join(map(str,frm)),\
                ', '.join(map(str,to))),\
                fontsize=3)#, fontweight='bold')
                # texts.append('mID:{0}, vID:{1}, bID:{2}::[{3}]->[{4}]'\
                # .format(mID,treevstrains[treestrain],bID,', '.join(map(str,frm)),\
                # ', '.join(map(str,to))))
        else:
            for newtreestrainID in vAbunds[vAbunds['t']==time]['tree_vstrain_id']:
                strainID = vKeepTreeStrainsDF[vKeepTreeStrainsDF['new_tree_vstrain_id'] \
                                == newtreestrainID]['vstrain_id'].values[0]
                if strainID in vstrainmatchIDs['vstrain_id'].values:
                    matchID = vstrainmatchIDs[vstrainmatchIDs.vstrain_id\
                    == strainID]['vmatch_id'].values[0]
                else:
                    matchID = 0
                bID = [id[0] for id in curSim.execute("SELECT infected_bstrain_id FROM vstrains \
                    WHERE run_id = {0}\
                    AND vstrain_id = {1}".format(run_id,strainID))][0]
                mID = vstrainmatchIDs[vstrainmatchIDs.vstrain_id\
                == strainID]['vmatch_id'].values[0]
                axesV[1].text(time, newtreestrainID,\
                '{1},mID:{0}'.format(strainID,matchID),\
                fontsize=3)#, fontweight='bold')
        ###
        ##
        # This applies info of interest on individual branches of microbial tree
        axesB[0].axvline(x=time, color='black', ls=':', lw=.7)
        axesB[1].axvline(x=time, color='black', ls=':', lw=.7)
        # axesS[0].axvline(x=time, color='black', ls=':', lw=.7)
        # axesS[1].axvline(x=time, color='black', ls=':', lw=.7)
        if hyperAnalyze:
            for newtreestrainID in bAbunds[bAbunds['t']==time]['tree_bstrain_id']:
                if bAbunds[(bAbunds['t']==time) & \
                    (bAbunds['tree_bstrain_id'] == newtreestrainID)]\
                    .abundance.values[0] == 0:
                    continue
                strainID = bKeepTreeStrainsDF[bKeepTreeStrainsDF['new_tree_bstrain_id'] \
                                == newtreestrainID]['bstrain_id'].values[0]
                if strainID in bstrainmatchIDs['bstrain_id'].values:
                    matchID = bstrainmatchIDs[bstrainmatchIDs.bstrain_id\
                    == strainID]['bmatch_id'].values[0]
                else:
                    matchID = 0
                vID = [id[0] for id in curSim.execute("SELECT infecting_vstrain_id FROM bstrains \
                    WHERE run_id = {0}\
                    AND bstrain_id = {1}".format(run_id,strainID))][0]
                sIDs = [id[0] for id in curSim.execute("SELECT spacer_id FROM bspacers \
                    WHERE run_id = {0}\
                    AND bstrain_id = {1}\
                    ORDER BY spacer_id".format(run_id,strainID))]
                axesB[1].text(time, newtreestrainID, '{0},{1},[{2}]::vID:{3}'\
                .format(strainID,matchID,', '.join(map(str,sIDs)),vID),\
                fontsize=3)#, fontweight='bold')
        else:
            for newtreestrainID in bAbunds[bAbunds['t']==time]['tree_bstrain_id']:
                strainID = bKeepTreeStrainsDF[bKeepTreeStrainsDF['new_tree_bstrain_id'] \
                                == newtreestrainID]['bstrain_id'].values[0]
                if strainID in bstrainmatchIDs['bstrain_id'].values:
                    matchID = bstrainmatchIDs[bstrainmatchIDs.bstrain_id\
                    == strainID]['bmatch_id'].values[0]
                else:
                    matchID = 0
                axesB[1].text(time, newtreestrainID, '{0},mID:{1}'\
                .format(strainID,matchID), fontsize=3)#, fontweight='bold')
