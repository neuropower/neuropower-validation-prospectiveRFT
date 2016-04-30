from __future__ import division
import os
import sys
import pandas as pd
import numpy as np

pilot_sub = int(sys.argv[1])
final_sub = int(sys.argv[2])
simstart = int(sys.argv[3])
simstop = int(sys.argv[4])
conditions = int(sys.argv[5])
modality = sys.argv[6]
adaptive = sys.argv[7]
threshold = sys.argv[8]
basefolder = sys.argv[9]
outfolder = sys.argv[10]

# pilot_sub = 15
# final_sub = 60
# simstart = 3
# simstop = 4
# conditions = 47
# modality = 'hcp'
# adaptive = 'adaptive'
# threshold = 'u2'
# basefolder = "/Users/Joke/Documents/Onderzoek/Studie_4_propow/InterimPower_Results"
# outfolder = "/Users/Joke/Documents/Onderzoek/Studie_4_propow/InterimPower_Results/tables/"

folder = os.path.join(basefolder+'/interim',modality+'_'+adaptive,threshold)
prefile = os.path.join(outfolder,'powpred_'+modality+'_'+adaptive+'_'+threshold+'.csv')
obsfile = os.path.join(outfolder,'powtrue_'+modality+'_'+adaptive+'_'+threshold+'.csv')
outfile = os.path.join(outfolder,'conditional_simstart_'+str(simstart)+"_"+modality+'_'+adaptive+'_'+threshold+'.csv')
pre = pd.read_csv(prefile)
true = pd.read_csv(obsfile)

simw = [2,4,6,8]*4
simef = np.repeat(['half','one','onehalf','two'],4)
pow = 0.5

results = []
for p in np.arange(simstart,simstop+1):
    for q in range(conditions):
        print(q)
        for m in ['BF','UN','RF','BH']:
            m1 = pre['mcp']==m
            m2 = pre['condition']==q+1
            m3 = pre['simulation']==p
            m4 = pre['power']>pow
            if np.sum([all(tup) for tup in zip(m1,m2,m3,m4)])==0:
                continue
            index = np.min(np.where([all(tup) for tup in zip(m1,m2,m3,m4)]))
            sub_nec = pre.loc[index].subjects
            m1 = true['mcp']==(m if not m == "RF" else "RFT")
            m2 = true['condition']==q+1
            m3 = true['simulation']==p
            m4 = true['TPR']>pow
            if np.sum([all(tup) for tup in zip(m1,m2,m3,m4)])==0:
                continue
            index = np.min(np.where([all(tup) for tup in zip(m1,m2,m3,m4)]))
            sub_realnec = true.loc[index].subjects
            m4 = true['subjects']==sub_nec
            if np.sum([all(tup) for tup in zip(m1,m2,m3,m4)])==0:
                continue
            index = np.where([all(tup) for tup in zip(m1,m2,m3,m4)])
            trueres = true.loc[index]
            result ={
            "mcp":m,
            "subjects":int(trueres.subjects),
            "truesubneed":int(sub_realnec),
            "sim":int(trueres.simulation),
            "condition":int(trueres.condition),
            "TPR":float(trueres.TPR),
            }
            results.append(result)

results = pd.DataFrame(results)
results.to_csv(outfile)
