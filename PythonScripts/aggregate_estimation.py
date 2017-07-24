from __future__ import division
import os
import sys
import pandas as pd
import numpy as np
import copy

PILOT = int(os.environ.get('PILOT'))
FINAL = int(os.environ.get('FINAL'))
SIMS = int(os.environ.get('SIMS'))
MODALITY = os.environ.get('MODALITY')
ADAPTIVE = os.environ.get('ADAPTIVE')
EXC = os.environ.get('EXC')
MODEL = os.environ.get('MODEL')
OUTDIR = os.environ.get('OUTDIR')
TABDIR = os.environ.get('TABDIR')

poweraim = 0.80

estimation_out = os.path.join(TABDIR,"estimation_"+MODALITY+'_'+ADAPTIVE+'_'+EXC+'_'+MODEL+'.csv')
prediction_out = os.path.join(TABDIR,"prediction_"+MODALITY+'_'+ADAPTIVE+'_'+EXC+'_'+MODEL+'.csv')
true_out = os.path.join(TABDIR,"true_"+MODALITY+'_'+ADAPTIVE+'_'+EXC+'_'+MODEL+'.csv')
conditional_out = os.path.join(TABDIR,"conditional_"+MODALITY+'_'+ADAPTIVE+'_'+EXC+'_'+MODEL+'.csv')

simw = [2,4,6,8]*4
cons = 47 if MODALITY=="HCP" else 16
simef = np.repeat(["050","080","100","120"],4)

estimation_all = []
prediction_all = []
true_all = []
conditional_all = []

for p in range(SIMS):
    # read file with estimation results
    estimation_file = os.path.join(OUTDIR,"estimation_"+MODALITY+"_"+str(p+1)+"%s.csv")%("" if MODALITY=="SIM" else "_0-48")
    if not os.path.isfile(estimation_file):
        continue
    estimation = pd.read_csv(estimation_file)
    estimation['sim']=p+1
    estimation = estimation.fillna(0)

    for c in range(cons):

        # read file with prediction of power
        if MODALITY=="SIM":
            post = "_w_"+str(simw[c])+"_e_"+str(simef[c])
        else:
            post = "_contrast_"+str(c)
        prediction_file = os.path.join(OUTDIR,"powpre_"+MODALITY+"_"+str(p+1)+str(post)+".csv")
        if not os.path.isfile(prediction_file):
            #print("not found")

            prediction = pd.DataFrame({
                'Unnamed:0':range(0,FINAL-PILOT),
                'BF':[0]*(FINAL-PILOT),
                'FDRc_predicted':[0]*(FINAL-PILOT),
                'RFT':[0]*(FINAL-PILOT),
                'UN':[0]*(FINAL-PILOT),
                'samplesize':range(PILOT,FINAL),
                'subjects':range(PILOT,FINAL)
                })

        else:
            prediction = pd.read_csv(prediction_file)
            if not 'BH' in prediction:
                    prediction['BH']=0
            prediction['subjects']=range(PILOT,FINAL)

        # replace FDR with predicted FDR

        prediction['BH'] = prediction['FDRc_predicted']

        # conditional: when is it > poweraim
        conditional = {}
        conditional["BFp"] = np.min(prediction.loc[prediction.BF>poweraim].subjects)
        conditional["RFp"] = np.min(prediction.loc[prediction.RFT>poweraim].subjects)
        conditional["UNp"] = np.min(prediction.loc[prediction.UN>poweraim].subjects)
        conditional["BHp"] = np.min(prediction.loc[prediction.BH>poweraim].subjects)

        # from wide to long
        prediction = pd.melt(prediction,id_vars=['subjects'],value_vars=['BF','UN','BH','RFT'])
        prediction.columns = ['subjects','mcp','power']
        prediction.mcp = [x if not x=="RFT" else "RF" for x in prediction.mcp]
        prediction['simulation'] = p+1
        prediction['condition'] = c+1

        true_file = os.path.join(OUTDIR,"powtru_"+MODALITY+"_"+str(p+1)+str(post)+".csv")
        if not os.path.isfile(true_file):
            true = pd.DataFrame({
                'Unnamed:0':range(0,FINAL-PILOT),
                'BF':[0]*(FINAL-PILOT),
                'BH':[0]*(FINAL-PILOT),
                'RFT':[0]*(FINAL-PILOT),
                'UN':[0]*(FINAL-PILOT),
                })
        else:
            true = pd.read_csv(true_file)
            if len(true)!=47:
                print(len(true))
                continue

        cols = true.columns
        cols = [x if not x == "RFT" else "RF" for x in cols]
        true.columns = cols
        true = true.iloc[:46]
        true['subjects']=range(PILOT,FINAL)

        # for simulated data: check FPR measures
        if MODALITY == 'SIM':
            true_FWE = copy.deepcopy(true)
            true_FPR = copy.deepcopy(true)
            true_FDR = copy.deepcopy(true)

            # read true file
            true["BH"] = true["BH_TP"]/(true["BH_TP"]+true["BH_FN"])
            true["BF"] = true["BF_TP"]/(true["BF_TP"]+true["BF_FN"])
            true["RF"] = true["RFT_TP"]/(true["RFT_TP"]+true["RFT_FN"])
            true["UN"] = true["UN_TP"]/(true["UN_TP"]+true["UN_FN"])
            true['subjects']=range(PILOT,FINAL)

            # true for FPR error rates
            true_FPR["BH"] = true_FPR["BH_FP"]/(true_FPR["BH_FP"]+true_FPR["BH_TN"])
            true_FPR["BF"] = true_FPR["BF_FP"]/(true_FPR["BF_FP"]+true_FPR["BF_TN"])
            true_FPR["RF"] = true_FPR["RFT_FP"]/(true_FPR["RFT_FP"]+true_FPR["RFT_TN"])
            true_FPR["UN"] = true_FPR["UN_FP"]/(true_FPR["UN_FP"]+true_FPR["UN_TN"])
            true_FPR['subjects']=range(PILOT,FINAL)

            # true for FDR error rates
            true_FDR["BH"] = true_FDR["BH_FP"]/(true_FDR["BH_TP"]+true_FDR["BH_FP"])
            true_FDR["BF"] = true_FDR["BF_FP"]/(true_FDR["BF_TP"]+true_FDR["BF_FP"])
            true_FDR["RF"] = true_FDR["RFT_FP"]/(true_FDR["RFT_TP"]+true_FDR["RFT_FP"])
            true_FDR["UN"] = true_FDR["UN_FP"]/(true_FDR["UN_TP"]+true_FDR["UN_FP"])
            true_FDR['subjects']=range(PILOT,FINAL)

            # true for FWE error rates
            true_FWE["BH"] = [1. if x>0 else 0. for x in true_FWE['BH_FP']]
            true_FWE["BF"] = [1. if x>0 else 0. for x in true_FWE['BF_FP']]
            true_FWE["RF"] = [1. if x>0 else 0. for x in true_FWE['RFT_FP']]
            true_FWE["UN"] = [1. if x>0 else 0. for x in true_FWE['UN_FP']]
            true_FWE['subjects']=range(PILOT,FINAL)

        # conditional: when is it > poweraim
        conditional["BFt"] = np.min(true.loc[true.BF>poweraim].subjects)
        conditional["RFt"] = np.min(true.loc[true.RF>poweraim].subjects)
        conditional["UNt"] = np.min(true.loc[true.UN>poweraim].subjects)
        conditional["BHt"] = np.min(true.loc[true.BH>poweraim].subjects)

        BFi = true.loc[true.subjects==conditional['BFp']].BF
        UNi = true.loc[true.subjects==conditional['UNp']].UN
        RFi = true.loc[true.subjects==conditional['RFp']].RF
        BHi = true.loc[true.subjects==conditional['BHp']].BH

        # conditional what is the power for predicted ss
        conditional["BFpow"] = 'nan' if len(BFi)==0 else min(BFi)
        conditional["RFpow"] = 'nan' if len(RFi)==0 else min(RFi)
        conditional["UNpow"] = 'nan' if len(UNi)==0 else min(UNi)
        conditional["BHpow"] = 'nan' if len(BHi)==0 else min(BHi)

        # wide to long
        true = pd.melt(true,id_vars=['subjects'],value_vars=['BF','UN','BH','RF'])
        true.columns = ['subjects','mcp','power']
        true['simulation'] = p+1
        true['condition'] = c+1

        # Conditional combine

        conditional_pd = pd.DataFrame(conditional,index=[0])
        # wide to long
        cond_melt = pd.melt(conditional_pd,id_vars=['subjects'],value_vars=conditional_pd.columns.tolist())
        cond_melt['mcp'] = [x[0:2] for x in cond_melt['variable']]
        cond_melt['val'] = [x[2:] for x in cond_melt['variable']]
        cond_melt_pd = cond_melt.pivot(index='mcp',columns='val',values='value')
        conditional_pd = cond_melt_pd.reset_index()
        conditional_pd.columns = ['mcp','predicted','power','true']
        conditional_pd['simulation']=p+1
        conditional_pd['condition']=c+1

        # for simulated data: check FPR measures
        if MODALITY == 'SIM':
            conditional_pd['FWE']=[min(true_FWE[x][true_FWE['subjects']==conditional[x+"p"]]) if conditional[x+"p"]>0 else 'nan' for x in ['BF','BH','RF','UN']]
            conditional_pd['FDR']=[min(true_FDR[x][true_FDR['subjects']==conditional[x+"p"]]) if conditional[x+"p"]>0 else 'nan' for x in ['BF','BH','RF','UN']]
            conditional_pd['FPR']=[min(true_FPR[x][true_FPR['subjects']==conditional[x+"p"]]) if conditional[x+"p"]>0 else 'nan' for x in ['BF','BH','RF','UN']]

        conditional_all.append(conditional_pd)

        prediction_all.append(prediction)

        # if c==14:
        #     np.load("stop")

        true_all.append(true)

    estimation_all.append(estimation)

# combine true and estimated
estimation_all = pd.concat(estimation_all,ignore_index=False)
estimation_all.to_csv(estimation_out)
prediction_all = pd.concat(prediction_all,ignore_index=False)
prediction_all.to_csv(prediction_out)
true_all = pd.concat(true_all,ignore_index=False)
true_all.to_csv(true_out)
conditional_all = pd.concat(conditional_all,ignore_index=False)
conditional_all.to_csv(conditional_out)
