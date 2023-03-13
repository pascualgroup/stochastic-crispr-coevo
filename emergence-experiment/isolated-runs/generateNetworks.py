import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import sqlite3
import matplotlib.cm as cm
import matplotlib.colors as mc

def generate(run_id,combo_id,replicate,\
conSim,conMatch,conTri,times,microbe_total,virus_total,\
colorpalSpacers,imgTypes,PLOT_PATH,resolve):
    for time in times:
        print('Compiling tripartite and infection networks at time = {}'.format(time))
        # tg.tripartiteGraph(time,run_id,sSpacing,vSpacing,bSpacing,html,\
        #                     colorpalSpacers,graphImgTypes,\
        #                     DBSIM_PATH,DBMATCH_PATH,DBTRI_PATH,PLOT_PATH)
        vTotal = virus_total[virus_total.t==time]['Viral Abundance'].values[0]
        vAbunds = pd.read_sql_query("SELECT vstrain_id,abundance FROM vabundance \
        WHERE run_id = {0} AND t = {1}".format(run_id,time), conSim)
        bTotal = microbe_total[microbe_total.t==time]['Microbial Abundance'].values[0]
        bAbunds = pd.read_sql_query("SELECT bstrain_id,abundance FROM babundance \
        WHERE run_id = {0} AND t = {1}".format(run_id,time), conSim)
        ###
        ## Single Match Tripartite Network
        #
        tripartiteNetwork = pd.read_sql_query(
            "SELECT vmatch_id, spacer_id, bmatch_id \
            FROM single_match_tripartite_networks WHERE t = {0}".format(time),
            conTri)
        ###
        ## The following code block orders the tripartite network microbe strains
        #
        bstrainmatchIDs = pd.read_sql_query(
            "SELECT match_id, bstrain_id \
            FROM bmatches WHERE t = {0}".format(time),
            conTri).rename(columns={"match_id": "bmatch_id"})
        bmatchAbund = pd.read_sql_query("SELECT match_id,babundance \
                FROM bmatches_abundances  \
                WHERE t = {0} AND match_id in ({1})".format(time,', '
                .join(map(str,tripartiteNetwork['bmatch_id'].unique()))), conTri)\
                .rename(columns={"match_id": "bmatch_id"}).rename(columns={"babundance": "bmatchAbund"})
        bstrainvstrain0 = pd.read_sql_query(
        "SELECT bstrain_id, vstrain_id \
        FROM bstrain_to_vstrain_0matches WHERE t = {0}".format(time), conMatch)
        # This adds up all of the viral abundances that can infect a microbe MATCH type
        bmatchcollapse = bstrainvstrain0.merge(bstrainmatchIDs,on=['bstrain_id'])\
                            .drop(columns='bstrain_id').drop_duplicates()\
                            .merge(vAbunds,on=['vstrain_id']).groupby(['bmatch_id'])\
                            .agg(prob=('abundance','sum')).reset_index()
        # And divides by total viral population size to compute a sampling probability
        bmatchcollapse['prob'] = 1/vTotal*bmatchcollapse['prob']
        # This adds column of associated match abundances
        bmatchcollapse = bmatchcollapse.merge(bmatchAbund,on=['bmatch_id'])
        bmatchcollapse['prob'] = bmatchcollapse['prob']*1/bTotal*np.array(bmatchcollapse['bmatchAbund'])
        bmatchcollapse = bmatchcollapse.sort_values(by=['prob']) # bin values are sorted
        # This adds the bmatchID that was not originally in the infection network
        bstrainsOrdered = list(set(bstrainmatchIDs['bmatch_id']) - set(bmatchcollapse['bmatch_id']))
        bstrainsOrdered = list(bmatchAbund[bmatchAbund.bmatch_id.isin(bstrainsOrdered)]\
                        .sort_values(by=['bmatchAbund'])['bmatch_id']) # bin values are sorted here
        bstrainsOrdered.extend(bmatchcollapse['bmatch_id'])
        ###
        ## The following code block orders the tripartite network viral strains
        #
        vmatchAbund = pd.read_sql_query("SELECT match_id,vabundance, bsusceptible \
                FROM vmatches_abundances  \
                WHERE t = {0} AND match_id in ({1}) \
                ORDER BY vabundance".format(time,', '
                .join(map(str,tripartiteNetwork['vmatch_id'].unique()))), conTri)\
                .rename(columns={"match_id": "vmatch_id"})\
                .rename(columns={"vabundance": "vmatchAbund",\
                "bsusceptible": "EscEventProb"})
        vmatchAbund['EscEventProb']= list(1/bTotal * np.array(vmatchAbund['EscEventProb'])\
        *1/vTotal* np.array(vmatchAbund['vmatchAbund']))
        vmatchAbund = vmatchAbund.sort_values(by=['EscEventProb']) # bin values are sorted here
        ###
        ## This assigns colors to (proto-)spacer matches
        #
        spacerIDs = list(tripartiteNetwork['spacer_id'].unique())
        spacerIDs.sort()
        spacerDict = {i:j for (i,j) in zip(spacerIDs,list(range(0,len(spacerIDs))))}
        print(spacerDict)
        tripartiteNetwork['spacer_id'] = [spacerDict[i] for i in tripartiteNetwork['spacer_id']]
        Mcmap = cm.get_cmap('{}'.format(colorpalSpacers))
        spacerCMAP = mc.LinearSegmentedColormap.from_list(
            'spacers', [Mcmap(i) for i in range(len(spacerIDs))], len(spacerIDs))
        ###
        ## This pivots tripartite network so it can be plotted as a matrix
        #
        tripartiteNetwork = tripartiteNetwork.pivot(index='vmatch_id',columns='bmatch_id',values='spacer_id')
        tripartiteNetwork = tripartiteNetwork.reindex(index = vmatchAbund['vmatch_id']) # reordered here
        tripartiteNetwork = tripartiteNetwork[bstrainsOrdered] # reordered here
        microbeBounds = [0]
        mtickLocation = []
        ###
        ## This creates bin sizes and tick marks for each match type
        #
        for i in bstrainsOrdered:
            j = bmatchAbund[bmatchAbund.bmatch_id == i]['bmatchAbund'].values[0]
            microbeBounds.append(microbeBounds[-1]+np.log(j)+1)
            mtickLocation.append(microbeBounds[-2]+(microbeBounds[-1]-microbeBounds[-2])/2)
        mtickLocation = list(1/max(microbeBounds)*np.array(mtickLocation))
        microbeBounds = list(1/max(microbeBounds)*np.array(microbeBounds))
        virusBounds = [0]
        vtickLocation = []
        for i in vmatchAbund['vmatchAbund']:
            # virusBounds.append(virusBounds[-1]+np.log(i)+1)
            # print('for vmatch_id {}'.format(i))
            virusBounds.append(virusBounds[-1]+np.log(i)+1)
            # print('virus bounds are {}'.format(virusBounds))
            vtickLocation.append(virusBounds[-2]+(virusBounds[-1]-virusBounds[-2])/2)
            # print('virus tick locations are {}'.format(vtickLocation))
        vtickLocation = list(1/max(virusBounds)*np.array(vtickLocation))
        virusBounds = list(1/max(virusBounds)*np.array(virusBounds))
        # print('final virus bounds are {}'.format(virusBounds))
        # print('final virus tick locations are {}'.format(vtickLocation))
        ###
        ## This plots the single match tripartite network as a matrix
        #
        fig, ax = plt.subplots(1)
        heat = ax.pcolormesh(microbeBounds, virusBounds, tripartiteNetwork,cmap=spacerCMAP,edgecolors='black',linewidth=0.01)
        ax.set_xticks(mtickLocation)
        ax.set_xticklabels(bstrainsOrdered, rotation=0,fontsize = 4)
        ax.set_yticks(vtickLocation)
        ax.set_yticklabels(vmatchAbund['vmatch_id'], rotation=0,fontsize = 3.5)
        cbticks = list(
        np.linspace(.5*(len(spacerIDs)-1)/len(spacerIDs),(len(spacerIDs)-0.5)
        *(len(spacerIDs)-1)/len(spacerIDs),len(spacerIDs)
        ))
        cbar = fig.colorbar(heat, ticks=cbticks)
        cbar.ax.set_yticklabels(spacerIDs)
        cbar.set_label("Spacer ID")
        ax.set_xlabel(xlabel ='Microbe Match ID',labelpad=10)
        ax.set_ylabel(ylabel = 'Virus Match ID',labelpad=10)
        ax.set_title('runID{0}-c{1}-r{2}\nTripartite Match Network at t = {3} (contact-ordered)'.format(run_id,combo_id,replicate,time),pad=20,fontsize=10)
        fig.tight_layout()
        for imgType in imgTypes:
            fig.savefig(os.path.join(PLOT_PATH,'tripartite-match-network-contact-ordered-time{0}.{1}'.format(int(np.floor(time)),imgType)),dpi=resolve)
        plt.close(fig)
        ###
        ## Infection Networks
        #
        # Infection Network with Match IDs
        vstrainmatchIDs = pd.read_sql_query(
            "SELECT match_id, vstrain_id \
            FROM vmatches WHERE t = {0}".format(time),
            conTri).rename(columns={"match_id": "vmatch_id"})
        # vmatch to bmatch infection network
        bmatchvmatch0 = pd.read_sql_query(
            "SELECT vstrain_id, bstrain_id \
            FROM bstrain_to_vstrain_0matches WHERE t = {0} \
            AND vstrain_id in ({1})".format(time,
            ', '.join(map(str,np.unique(vstrainmatchIDs['vstrain_id'])))), conMatch)\
            .merge(vstrainmatchIDs,on=['vstrain_id']).drop(columns=['vstrain_id']).drop_duplicates()\
            .merge(bstrainmatchIDs,on=['bstrain_id']).drop(columns=['bstrain_id']).drop_duplicates()
        bmatchvmatch0['presence'] = bmatchvmatch0['vmatch_id'].size*[1]
        bmatchvmatch0 = bmatchvmatch0.pivot(index='vmatch_id',columns='bmatch_id',values='presence')
        bmatchvmatch0 = bmatchvmatch0.reindex(index = vmatchAbund['vmatch_id'])
        bmatchvmatch0 = bmatchvmatch0[bmatchcollapse['bmatch_id'].values]

        microbeBounds = [0]
        mtickLocation = []
        for i in bmatchcollapse['bmatchAbund']:
            microbeBounds.append(microbeBounds[-1]+np.log(i)+1)
            mtickLocation.append(microbeBounds[-2]+(microbeBounds[-1]-microbeBounds[-2])/2)
        mtickLocation = list(1/max(microbeBounds)*np.array(mtickLocation))
        microbeBounds = list(1/max(microbeBounds)*np.array(microbeBounds))
        fig, ax = plt.subplots(1)
        ax.pcolormesh(microbeBounds, virusBounds, bmatchvmatch0, edgecolors='black',linewidth=0.01)
        ax.set_xticks(mtickLocation)
        ax.set_xticklabels(bmatchcollapse['bmatch_id'], rotation=0,fontsize = 4)
        ax.set_yticks(vtickLocation)
        ax.set_yticklabels(vmatchAbund['vmatch_id'], rotation=0,fontsize = 3.5)
        ax.set_xlabel(xlabel ='Microbe Match ID',labelpad=10)
        ax.set_ylabel(ylabel = 'Virus Match ID',labelpad=10)
        ax.set_title('runID{0}-c{1}-r{2}\nInfection Network of Match Types at t = {3} (contact-ordered)'.format(run_id,combo_id,replicate,time),pad=20,fontsize=10)
        fig.tight_layout()
        for imgType in imgTypes:
            fig.savefig(os.path.join(PLOT_PATH,'infection-network-match-types-contact-ordered-time{0}.{1}'.format(int(np.floor(time)),imgType)),dpi=resolve)
        plt.close(fig)
        ##
        # Infection Network with Strain IDs
        bstraincollapse = bstrainvstrain0.merge(vAbunds,on=['vstrain_id'])\
                            .groupby(['bstrain_id']).agg(prob=('abundance','sum'))\
                            .reset_index()
        bstraincollapse['prob'] = 1/vTotal*np.array(bstraincollapse['prob'])
        bstraincollapse = bstraincollapse.merge(bAbunds,on=['bstrain_id'])
        bstraincollapse['prob'] = 1/bTotal*np.array(bstraincollapse['prob'])*\
                                    np.array(bstraincollapse['abundance'])
        bstraincollapse = bstraincollapse.sort_values(by=['prob'])

        vstraincollapse = bstrainvstrain0.merge(bAbunds,on=['bstrain_id'])\
                            .groupby(['vstrain_id']).agg(prob=('abundance','sum'))\
                            .reset_index()
        vstraincollapse['prob'] = 1/bTotal*np.array(vstraincollapse['prob'])
        vstraincollapse = vstraincollapse.merge(vAbunds,on=['vstrain_id'])
        vstraincollapse['prob'] = 1/vTotal*np.array(vstraincollapse['prob'])*\
                                    np.array(vstraincollapse['abundance'])
        vstraincollapse = vstraincollapse.sort_values(by=['prob'])
        bstrainvstrain0['presence'] = bstrainvstrain0['vstrain_id'].size*[1]
        bstrainvstrain0 = bstrainvstrain0.pivot(index='vstrain_id',columns='bstrain_id',values='presence')
        bstrainvstrain0 = bstrainvstrain0.reindex(index = vstraincollapse['vstrain_id'])
        bstrainvstrain0 = bstrainvstrain0[bstraincollapse['bstrain_id'].values]

        microbeBounds = [0]
        mtickLocation = []
        for i in bstraincollapse['abundance']:
            microbeBounds.append(microbeBounds[-1]+np.log(i)+1)
            mtickLocation.append(microbeBounds[-2]+(microbeBounds[-1]-microbeBounds[-2])/2)
        mtickLocation = list(1/max(microbeBounds)*np.array(mtickLocation))
        microbeBounds = list(1/max(microbeBounds)*np.array(microbeBounds))
        virusBounds = [0]
        vtickLocation = []
        for i in vstraincollapse['abundance']:
            # virusBounds.append(virusBounds[-1]+np.log(i)+1)
            # print('for vmatch_id {}'.format(i))
            virusBounds.append(virusBounds[-1]+np.log(i)+1)
            # print('virus bounds are {}'.format(virusBounds))
            vtickLocation.append(virusBounds[-2]+(virusBounds[-1]-virusBounds[-2])/2)
            # print('virus tick locations are {}'.format(vtickLocation))
        vtickLocation = list(1/max(virusBounds)*np.array(vtickLocation))
        virusBounds = list(1/max(virusBounds)*np.array(virusBounds))
        fig, ax = plt.subplots(1)
        ax.pcolormesh(microbeBounds, virusBounds, bstrainvstrain0, edgecolors='black',linewidth=0.01)
        ax.set_xticks(mtickLocation)
        ax.set_xticklabels(bstraincollapse['bstrain_id'], rotation=0,fontsize = 4)
        ax.set_yticks(vtickLocation)
        ax.set_yticklabels(vstraincollapse['vstrain_id'], rotation=0,fontsize = 3.5)
        ax.set_xlabel(xlabel ='Microbe Strain ID',labelpad=10)
        ax.set_ylabel(ylabel = 'Virus Strain ID',labelpad=10)
        ax.set_title('runID{0}-c{1}-r{2}\nInfection Network of Strains at t = {3} (contact-ordered)'.format(run_id,combo_id,replicate,time),pad=20,fontsize=10)
        fig.tight_layout()
        for imgType in imgTypes:
            fig.savefig(os.path.join(PLOT_PATH,'infection-network-strains-contact-ordered-time{0}.{1}'.format(int(np.floor(time)),imgType)),dpi=resolve)
        plt.close(fig)

    return vstrainmatchIDs, bstrainmatchIDs
