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
outfile = os.path.join(outfolder,modality+'_'+adaptive+'_'+threshold+'.csv')

results = []
for p in range(sims):
    file = os.path.join(folder,"estimation_"+modality+"_"+str(p+1)+".csv")
    res = pd.read_csv(file,header=None)
    res['sim']=p+1
    results.append(res)

results = pd.concat(results,ignore_index=False)
results.columns=['condition','pi0t','pi0e','pi0c','effe','efft','effc','effest','sde','sdt','sim']
results.to_csv(outfile)
