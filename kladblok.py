cons <- read.table("80_Unrelated_contrasts/contrasts.txt")$V1
paradigms <- read.table("80_Unrelated_contrasts/paradigms.txt")$V1
subID <- read.table("80_Unrelated_contrasts/subID.txt")$V1

pars <- c("SOCIAL","EMOTION","RELATIONAL","LANGUAGE","GAMBLING","WM","MOTOR")

commands <- c()
for (i in 1:80){
k <- 0
for (c in 1:47){
k <- k+1
file_in <- paste("80_Unrelated_contrasts/",subID[i],"/",pars[paradigms[c]],"/Level2/cope",cons[c],".feat/cope1.nii.gz",sep="")
file_out <- paste("80_Unrelated_cope/","subject_",i,"_contrast_",k,".nii.gz",sep="")
commands <- c(commands,paste("cp",file_in,file_out))
}
}

write.table(commands,"commands.sh",row.names=FALSE,col.names=FALSE,quote=FALSE)

=========================


import os
import sys
sys.path.append("/share/PI/russpold/software/anaconda/lib/python2.7/site-packages")
import numpy as np
import scipy
import math
import nibabel as nib
import shutil
import csv
import imp
import subprocess

RESDIR="/scratch/users/jdurnez/power_peak_HCP/"
TMPDIR="/local-scratch/jdurnez/6958122/"
HOMEDIR="/home/jdurnez/power_peak/"
DATADIR="/scratch/PI/russpold/data/HCP/"

os.chdir(TMPDIR)

for c in range(47):
  print(c)
  contrast = [line.rstrip('\n') for line in open(os.path.join(HOMEDIR,'HCP_contrasts.txt'))][c]
  paradigm = [line.rstrip('\n') for line in open(os.path.join(HOMEDIR,'HCP_paradigms.txt'))][c]
  subjects = [line.rstrip('\n') for line in open(os.path.join(HOMEDIR,'HCP_unrelated_subID.txt'))]
  disks = [line.rstrip('\n') for line in open(os.path.join(HOMEDIR,'HCP_unrelated_disk.txt'))]
  masks = []
  for sub in range(len(subjects)):
      masks.append(os.path.join(DATADIR,"Disk"+str(disks[sub])+"of5",subjects[sub],"MNINonLinear/Results/tfMRI_"+paradigm,"tfMRI_"+paradigm+"_hp200_s4_level2vol.feat","cope"+contrast+".feat","stats/mask.nii.gz"))
  combine_cmd = "fslmerge -t masks.nii.gz "+" ".join(masks)
  os.popen(combine_cmd).read()
  mean_cmd = "fslmaths masks.nii.gz -Tmean masks_mean.nii.gz"
  os.popen(mean_cmd).read()
  mask_cmd = "fslmaths masks_mean.nii.gz -thr 1 mask_contrast_"+str(c)+".nii.gz"
  os.popen(mask_cmd).read()


masks = []
for c in range(47):
    masks.append("mask_contrast_"+str(c)+".nii.gz")

combine_cmd = "fslmerge -t masks.nii.gz "+" ".join(masks)
os.popen(combine_cmd).read()
mean_cmd = "fslmaths masks.nii.gz -Tmean masks_mean.nii.gz"
os.popen(mean_cmd).read()
mask_cmd = "fslmaths masks_mean.nii.gz -thr 1 mask_final.nii.gz"
os.popen(mask_cmd).read()

=====================

histogram

from palettable.colorbrewer.qualitative import Paired_12
import matplotlib.pyplot as plt

est_eff
est_sd
bum['pi1']

twocol = Paired_12.mpl_colors

xn = np.arange(-10,10,0.01)
nul = neuropower.nulPDF(xn,exc=2,method="RFT")
alt = neuropower.altPDF(xn,mu=est_eff,sigma=est_sd,exc=2,method="RFT")
mix = neuropower.mixprobdens(xn,pi1=bum['pi1'],mu=est_eff,sigma=est_sd,exc=2,method="RFT")

plt.figure(figsize=(7,5))
plt.hist(peaks.peak,lw=0,facecolor=twocol[0],normed=True,bins=np.arange(2,10,0.3),label="observed distribution")
plt.xlim([0,10])
plt.ylim([0,3])
plt.plot(xn,nul,color=twocol[3],lw=3,label="null distribution")
plt.plot(xn,alt,color=twocol[5],lw=3, label="alternative distribution")
plt.plot(xn,mix,color=twocol[1],lw=3,label="total distribution")
plt.title("histogram")
plt.xlabel("peak height")
plt.ylabel("density")
plt.legend(loc="upper right",frameon=False)
plt.show()

colors=[twocol[0],twocol[2]]
xn = np.arange(0,1,0.01)
alt = [1-bum['pi1']]+scipy.stats.beta.pdf(xn, bum['a'], 1)
null = [1-bum['pi1']]*len(xn)
peks = [[peaks.pval[peaks.active==0]],[peaks.pval[peaks.active==1]]]
plt.figure(figsize=(7,5))
plt.hist(peks,lw=0,normed=True,stacked=True,fill=True,color=colors,bins=np.arange(0,1,0.05),label="observed distribution")
plt.ylim([0,3])
plt.plot(xn,null,color=twocol[1],lw=3,label="null distribution")
plt.plot(xn,alt,color=twocol[3],lw=3,label="null distribution")
plt.show()


===========================

from palettable.colorbrewer.qualitative import Paired_12
import matplotlib.pyplot as plt
twocol = Paired_12.mpl_colors



TRUEDIR = os.path.join(TMPDIR,"analysis_true_"+str(seed)+"/")
truefile = os.path.join(TRUEDIR,'out_zstat1.nii.gz')
tvals = nib.load(truefile).get_data()
tvals = [item for sublist in tvals for item in sublist]
tvals = [item for sublist in tvals for item in sublist]
tvals = [x for x in tvals if x != 0.0]


plt.figure(figsize=(7,5))
plt.hist(tvals,lw=0,facecolor=twocol[0],normed=True,bins=np.arange(-10,20,0.3),label="observed distribution")
plt.show()


==========================
from palettable.colorbrewer.qualitative import Paired_12
import matplotlib.pyplot as plt
twocol = Paired_12.mpl_colors

xn = np.arange(-10,10,0.01)
alt10 = neuropower.altPDF(xn,mu=est_eff,sigma=est_sd,exc=2,method="RFT")
est_eff_20 = est_eff/np.sqrt(9)*np.sqrt(19)
alt20 = neuropower.altPDF(xn,mu=est_eff_20,sigma=est_sd,exc=2,method="RFT")
est_eff_30 = est_eff/np.sqrt(9)*np.sqrxt(29)
alt30 = neuropower.altPDF(xn,mu=est_eff_30,sigma=est_sd,exc=2,method="RFT")
est_eff_50 = est_eff/np.sqrt(9)*np.sqrt(49)
alt50 = neuropower.altPDF(xn,mu=est_eff_50,sigma=est_sd,exc=2,method="RFT")

plt.figure(figsize=(7,5))
plt.hist(pks10.tolist(),lw=0,facecolor=twocol[0],normed=True,bins=np.arange(-10,20,0.3),label="observed distribution")
plt.plot(xn,alt10,color=twocol[1],lw=3, label="alternative distribution")
plt.hist(pks30.tolist(),lw=0,facecolor=twocol[4],normed=True,bins=np.arange(-10,20,0.3),label="observed distribution")
plt.plot(xn,alt30,color=twocol[5],lw=3, label="alternative distribution")
plt.hist(pks50.tolist(),lw=0,facecolor=twocol[8],normed=True,bins=np.arange(-10,20,0.3),label="observed distribution")
plt.plot(xn,alt50,color=twocol[9],lw=3, label="alternative distribution")
plt.xlim([2,11])
plt.show()

############################################
in R check which subjects have all contrasts
############################################

HOMEDIR <- "/home/jdurnez/power_peak/"
DATADIR <- "/scratch/PI/russpold/data/HCP/"

cons <- read.table(paste(HOMEDIR,'HCP_contrasts.txt',sep=""))$V1
pars <- read.table(paste(HOMEDIR,'HCP_paradigms.txt',sep=""))$V1
subs <- read.table(paste(HOMEDIR,'HCP_unrelated_subID.txt',sep=""))$V1
disks <- read.table(paste(HOMEDIR,'HCP_unrelated_disk.txt',sep=""))$V1

res <- array(NA,dim=c(47,206))

for (s in 1:206){
  for (c in 1:47){
    file <- paste(DATADIR,"Disk",disks[s],"of5/",subs[s],"/MNINonLinear/Results/tfMRI_",pars[c],"/tfMRI_",pars[c],"_hp200_s4_level2vol.feat/cope",cons[c],".feat/stats/cope1.nii.gz",sep="")
    res[c,s] <- file.exists(file)
  }
}

hasmissing <- which(apply(res,2,mean)!=1)

newsubs <- subs[-hasmissing]
newdisks <- disks[-hasmissing]

write.table(newsubs,paste(HOMEDIR,'HCP_unrelated_subID.txt',sep=""),col.names=FALSE,row.names=FALSE)
write.table(newdisks,paste(HOMEDIR,'HCP_unrelated_disk.txt',sep=""),col.names=FALSE,row.names=FALSE)

############################################
in R check which subjects have all conditions
############################################

RESDIR <- "/scratch/users/jdurnez/power_peak_SIM/"
a <- c()
for (s in c(1:500)){
  file <- paste(RESDIR,"estimation_sim_",s,".csv",sep="")
  if(dim(read.table(file,sep=",",header=FALSE))[1]!=16){a <- c(a,s)}
}

a <- c(a,410)
HOMEDIR="/home/jdurnez/power_peak/"
write.table(a,paste(HOMEDIR,"SIM-miss.txt",sep=""),row.names=FALSE,col.names=FALSE)
