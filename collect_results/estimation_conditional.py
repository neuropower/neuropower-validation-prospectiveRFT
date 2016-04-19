from __future__ import division
import os
import sys
import pandas as pd
import numpy as np

pilot_sub = int(sys.argv[1])
final_sub = int(sys.argv[2])
sims = int(sys.argv[3])
conditions = int(sys.argv[4])
modality = sys.argv[5]
adaptive = sys.argv[6]
threshold = sys.argv[7]
basefolder = sys.argv[8]
outfolder = sys.argv[9]

folder = os.path.join(basefolder+'/interim',modality+'_'+adaptive,threshold)
prefile = os.path.join(outfolder,'powpred_'+modality+'_'+adaptive+'_'+threshold+'.csv')
obsfile = os.path.join(outfolder,'powtrue_'+modality+'_'+adaptive+'_'+threshold+'.csv')
outfile = os.path.join(outfolder,'conditional'+modality+'_'+adaptive+'_'+threshold+'.csv')
pre = pd.read_csv(prefile)
true = pd.read_csv(obsfile)

simw = [2,4,6,8]*4
simef = np.repeat(['half','one','onehalf','two'],4)

results = []
for p in range(sims):
    for q in range(conditions):
        for m in ['BF','UN','RF']:
            m1 = pre['mcp']==m
            m2 = pre['condition']==q+1
            m3 = pre['simulation']==p+1
            m4 = pre['power']>0.8
            index = np.min(np.where([all(tup) for tup in zip(m1,m2,m3,m4)]))
            sub_nec = pre.loc[index].subjects
            m1 = true['mcp']==m
            m2 = true['condition']==q+1
            m3 = true['simulation']==p+1
            m4 = true['subjects']==sub_nec
            index = np.where([all(tup) for tup in zip(m1,m2,m3,m4)])
            trueres = true.loc[index]
            result ={
            "mcp":trueres.mcp,
            "subjects":trueres.subjects,
            "sim":trueres.simulation,
            "condition":trueres.condition,
            "TPR":trueres.TPR,
            "FPR":trueres.FPR,
            "FDR":trueres.FDR,
            "FWER":trueres.FWER
            }
results.append(result)
results = pd.concat(results)
results.to_csv(outfile)
