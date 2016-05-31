from __future__ import division
import os
import sys
import pandas as pd
import numpy as np

PILOT = int(sys.argv[1])
FINAL = int(sys.argv[2])
SIMS = int(sys.argv[3])
MODALITY = sys.argv[4]
ADAPTIVE = sys.argv[5]
EXC = sys.argv[6]
MODEL = sys.argv[7]
OUTDIR = os.environ.get('OUTDIR')
TABDIR = os.environ.get('TABDIR')

estimation_out = os.path.join(TABDIR,"estimation_"+MODALITY+'_'+ADAPTIVE+'_'+EXC+'_'+MODEL+'.csv')
prediction_out = os.path.join(TABDIR,"prediction_"+MODALITY+'_'+ADAPTIVE+'_'+EXC+'_'+MODEL+'.csv')
true_out = os.path.join(TABDIR,"true_"+MODALITY+'_'+ADAPTIVE+'_'+EXC+'_'+MODEL+'.csv')
conditional_out = os.path.join(TABDIR,"conditional_"+MODALITY+'_'+ADAPTIVE+'_'+EXC+'_'+MODEL+'.csv')

simw = [2,4,6,8]*4
simef = np.repeat(['half','one','onehalf','two'],4)

estimation_all = []
prediction_all = []
true_all = []
conditional_all = []

for p in range(SIMS):
    print(p)

    # read file with estimation results
    estimation_file = os.path.join(OUTDIR,"estimation_"+MODALITY+"_"+str(p+1)+".csv")
    if not os.path.isfile(estimation_file):
        continue
    estimation = pd.read_csv(estimation_file)
    estimation['sim']=p+1
    estimation['condition']=[c+1 for c in range(len(estimation))]

    for c in range(16):

        # read file with prediction of power
        prediction_file = os.path.join(OUTDIR,"powpre_"+MODALITY+"_"+str(p+1)+"_w_"+str(simw[c])+"_e_"+str(simef[c])+".csv")
        if not os.path.isfile(prediction_file):
            continue
        prediction = pd.read_csv(prediction_file)
        if not 'BH' in prediction:
                prediction['BH']='nan'
        prediction['subjects']=range(PILOT,FINAL)

        # conditional: when is it > 0.8
        conditional = {}
        conditional["BFp"] = np.min(prediction.loc[prediction.BF>0.8].subjects)
        conditional["RFp"] = np.min(prediction.loc[prediction.RFT>0.8].subjects)
        conditional["UNp"] = np.min(prediction.loc[prediction.UN>0.8].subjects)
        conditional["BHp"] = np.min(prediction.loc[prediction.BH>0.8].subjects)

        # from wide to long
        prediction = pd.melt(prediction,id_vars=['subjects'],value_vars=['BF','UN','BH','RFT'])
        prediction.columns = ['subjects','mcp','power']
        prediction.mcp = [x if not x=="RFT" else "RF" for x in prediction.mcp]
        prediction['simulation'] = p+1
        prediction['condition'] = c+1

        true_file = os.path.join(OUTDIR,"powtru_"+MODALITY+"_"+str(p+1)+"_w_"+str(simw[c])+"_e_"+str(simef[c])+".csv")
        if not os.path.isfile(true_file):
            continue
        true = pd.read_csv(true_file)

        # read true file
        true["BH"] = true["BH_TP"]/(true["BH_TP"]+true["BH_FN"])
        true["BF"] = true["BF_TP"]/(true["BF_TP"]+true["BF_FN"])
        true["RF"] = true["RFT_TP"]/(true["RFT_TP"]+true["RFT_FN"])
        true["UN"] = true["UN_TP"]/(true["UN_TP"]+true["UN_FN"])
        true['subjects']=range(PILOT,FINAL)

        # conditional: when is it > 0.8
        conditional["BFt"] = np.min(true.loc[true.BF>0.8].subjects)
        conditional["RFt"] = np.min(true.loc[true.RF>0.8].subjects)
        conditional["UNt"] = np.min(true.loc[true.UN>0.8].subjects)
        conditional["BHt"] = np.min(true.loc[true.BH>0.8].subjects)

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

        prediction_all.append(prediction)
        true_all.append(true)
        conditional_all.append(conditional_pd)

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
