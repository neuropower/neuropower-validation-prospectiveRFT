from __future__ import division
import os
import sys
import pandas as pd

pilot_sub = int(sys.argv[1])
final_sub = int(sys.argv[2])
sims = int(sys.argv[3])
modality = sys.argv[4]
adaptive = sys.argv[5]
threshold = sys.argv[6]
basefolder = sys.argv[7]
outfolder = sys.argv[8]

folder = os.path.join(basefolder+'/interim',modality+'_'+adaptive,threshold)
outfile = os.path.join(outfolder,"estimation_"+modality+'_'+adaptive+'_'+threshold+'.csv')

results = []
for p in range(sims):
    file = os.path.join(folder,"estimation_"+modality+"_"+str(p+1)+".csv")
    if not os.path.isfile(file):
        continue
    res = pd.read_csv(file,header=None)
    res['sim']=p+1
    res['condition']=[c+1 for c in range(len(res))]
    results.append(res)

results = pd.concat(results,ignore_index=False)
if modality == "sim":
    results.columns=['eff','ac','pi1e','pi1t','effe','efft','effexp','sde','sdt','a','sim','condition']
else:
    results.columns=['condition','pi1t','pi1e','pi1c','effe','efft','effc','effexp','sde','sdt','sim','condition']
results.to_csv(outfile)
