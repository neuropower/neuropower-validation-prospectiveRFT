from __future__ import division
import os
import sys
import pandas as pd

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
outfile = os.path.join(outfolder,'powpred_'+modality+'_'+adaptive+'_'+threshold+'.csv')

results = []
for p in range(sims):
    for c in range(conditions):
        file = os.path.join(folder,"powpre_"+modality+"_"+str(p+1)+"_contrast_"+str(c)+".csv")
        if not os.path.isfile(file):
            continue
        res = pd.read_csv(file)
        res['sim']=p+1
        if not 'BH' in res:
                res['BH']='nan'
        res['subjects']=range(pilot_sub,final_sub)
        longres = pd.melt(res,id_vars=['subjects'],value_vars=['BF','UN','BH','RFT'])
        longres.columns = ['subjects','mcp','power']
        longres['simulation'] = p
        longres['condition'] = c
        results.append(longres)

results = pd.concat(results,ignore_index=False)
results.to_csv(outfile)
