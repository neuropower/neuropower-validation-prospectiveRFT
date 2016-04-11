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
outfile = os.path.join(outfolder,'powpred_'+modality+'_'+adaptive+'_'+threshold+'.csv')

simw = [2,4,6,8]*4
simef = np.repeat(['half','one','onehalf','two'],4)

results = []
for p in range(sims):
    for c in range(conditions):
        if modality == 'hcp':
            file = os.path.join(folder,"powpre_"+modality+"_"+str(p+1)+"_contrast_"+str(c)+".csv")
        else:
            file = os.path.join(folder,"powpre_"+modality+"_"+str(p+1)+"_w_"+str(simw[c])+"_e_"+str(simef[c])+".csv")
        if not os.path.isfile(file):
            print(file)
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

results = pd.concat(results)
results.to_csv(outfile)
