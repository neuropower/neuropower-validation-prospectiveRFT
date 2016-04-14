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

pre = pd.read_csv(prefile)
true = pd.read_csv(obsfile)

simw = [2,4,6,8]*4
simef = np.repeat(['half','one','onehalf','two'],4)

for p in range(sims):
for q in range(conditions):
for m in ['BF','UN','BH','RFT']:
m1 = pre['mcp']==m
m2 = pre['condition']==q+1
m3 = pre['simulation']==p+1
index = [all(tup) for tup in zip(m1,m2,m3)].index(True)

# nu eerst bepalen welke van de subjects >0.8
# daarna minimum bepalen

results = []
for p in range(sims):
    for c in range(conditions):

if modality == 'hcp':
    file1 = os.path.join(folder,"powpre_"+modality+"_"+str(p+1)+"_contrast_"+str(c)+".csv")
    file2 = os.path.join(folder,"powtru_"+modality+"_"+str(p+1)+"_contrast_"+str(c)+".csv")
else:
    file1 = os.path.join(folder,"powpre_"+modality+"_"+str(p+1)+"_w_"+str(simw[c])+"_e_"+str(simef[c])+".csv")
    file2 = os.path.join(folder,"powtru_"+modality+"_"+str(p+1)+"_w_"+str(simw[c])+"_e_"+str(simef[c])+".csv")
if not os.path.isfile(file1):
    print(file1)
    continue
pre = pd.read_csv(file1)
true = pd.read_csv(file2)
res['sim']=p+1
if not 'BH' in res:
        res['BH']='nan'
res['subjects']=range(pilot_sub,final_sub)
longres = pd.melt(res,id_vars=['subjects'],value_vars=['BF','UN','BH','RFT'])
longres.columns = ['subjects','mcp','power']
longres.mcp = [x if not x=="RFT" else "RF" for x in longres.mcp]
longres['simulation'] = p+1
longres['condition'] = c+1
results.append(longres)

results = pd.concat(results)
results.to_csv(outfile)


mcps <- c("UN","BH","BF","RF")
method_n <- c("Uncorrected","Benjamini-Hochberg","Bonferroni","Random Field Theory")

# compute minimal sample size


pre.ad.min <- pre.ad[pre.ad$power>0.8,]
minss <- ddply(pre.ad.min,
               ~simulation+condition+mcp,
               summarize,
               SS=min(subjects)
)


results <- data.frame()
for (l in 1:dim(minss)[1]){
  if(l%%100==0){print(l)}
  ind <- which(obs.ad$subjects==minss$SS[l] & obs.ad$condition==minss$condition[l] & obs.ad$simulation==minss$simulation[l] & obs.ad$mcp==minss$mcp[l])
  res <- unlist(obs.ad[ind,])
  results <- rbind(results,res)
}
