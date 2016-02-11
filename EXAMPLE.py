############
# preamble #
############

import os
import sys
import numpy as np
import scipy
import math
import nibabel as nib
import shutil
import csv
import sys
import imp
import pandas as pd
from palettable.colorbrewer.qualitative import Paired_12
from palettable.colorbrewer.qualitative import Set1_9
import matplotlib.pyplot as plt
import pylab as plb
import matplotlib as mpl

os.chdir("../ProspectivePower-Functions/")
import BUM
import neuropower
import cluster
import peakdistribution
import simul_multisubject_fmri_dataset
import model

exc = 2
FIGDIR="/Users/Joke/Documents/Onderzoek/Studie_4_propow/ProspectivePower-Paper/Studie_4_v1.4/Figures/"

os.chdir("../ProspectivePower-Validation/")
maskfile = "EXAMPLE_mask.nii"

SPM = nib.load("EXAMPLE_zstat1.nii").get_data()
peaks = cluster.cluster(SPM,exc)

# compute P-values

pvalues = np.exp(-exc*(np.array(peaks.peak)-exc))
pvalues = [max(10**(-6),t) for t in pvalues]
peaks['pval'] = pvalues

# estimate model

bum = BUM.bumOptim(peaks['pval'].tolist(),starts=10)
modelfit = neuropower.modelfit(peaks.peak,bum['pi1'],exc=exc,starts=10,method="RFT")

# predict power

thresholds = neuropower.threshold(peaks.peak,peaks.pval,FWHM=8,mask=nib.load(maskfile),alpha=0.05,exc=exc)
effect_cohen = modelfit['mu']/np.sqrt(13)
power_predicted = []
newsubs = range(10,71)
for s in newsubs:
    projected_effect = effect_cohen*np.sqrt(s)
    powerpred =  {k:1-neuropower.altCDF(v,projected_effect,modelfit['sigma'],exc=exc,method="RFT") for k,v in thresholds.items() if v!='nan'}
    power_predicted.append(powerpred)

power_predicted_df = pd.DataFrame(power_predicted)

# figure modelfit

twocol = Paired_12.mpl_colors

xn = np.arange(-10,10,0.01)
nul = [1-bum['pi1']]*neuropower.nulPDF(xn,exc=exc,method="RFT")
alt = bum['pi1']*neuropower.altPDF(xn,mu=modelfit['mu'],sigma=modelfit['sigma'],exc=exc,method="RFT")
mix = neuropower.mixprobdens(xn,pi1=bum['pi1'],mu=modelfit['mu'],sigma=modelfit['sigma'],exc=2,method="RFT")

xn_p = np.arange(0,1,0.01)
alt_p = [1-bum['pi1']]*scipy.stats.beta.pdf(xn_p, bum['a'], 1)+1-bum['pi1']
null_p = [1-bum['pi1']]*len(xn_p)

mpl.rcParams['font.size']='11.0'

fig,axs=plt.subplots(1,2,figsize=(14,5))
fig.subplots_adjust(hspace=.5,wspace=0.3)
axs=axs.ravel()
axs[0].hist(peaks.pval,lw=0,normed=True,facecolor=twocol[0],bins=np.arange(0,1.1,0.1),label="observed distribution")
axs[0].set_ylim([0,3])
axs[0].plot(xn_p,null_p,color=twocol[3],lw=2,label="null distribution")
axs[0].plot(xn_p,alt_p,color=twocol[5],lw=2,label="alternative distribution")
axs[0].legend(loc="upper right",frameon=False)
axs[0].set_title("Distribution of 81 peak p-values \n $\pi_1$ = "+str(round(bum['pi1'],2)))
axs[0].set_xlabel("Peak p-values")
axs[0].set_ylabel("Density")

axs[1].hist(peaks.peak,lw=0,facecolor=twocol[0],normed=True,bins=np.arange(2,10,0.3),label="observed distribution")
axs[1].set_xlim([2,7])
axs[1].set_ylim([0,1])
axs[1].plot(xn,nul,color=twocol[3],lw=2,label="null distribution")
axs[1].plot(xn,alt,color=twocol[5],lw=2, label="alternative distribution")
axs[1].plot(xn,mix,color=twocol[1],lw=2,label="total distribution")
axs[1].set_title("Distribution of peak heights \n $\delta_1$ = "+str(round(modelfit['mu']/np.sqrt(13),2)))
axs[1].set_xlabel("Peak heights (z-values)")
axs[1].set_ylabel("Density")
axs[1].legend(loc="upper right",frameon=False)

os.chdir(FIGDIR)
plt.savefig("FIG_EX_pi0.pdf")
plt.close()

# Figure power curves

colset1 = Set1_9.mpl_colors

BFmin = np.min([i for i,elem in enumerate(power_predicted_df['BF']>0.8,1) if elem])+10
RFTmin = np.min([i for i,elem in enumerate(power_predicted_df['RFT']>0.8,1) if elem])+10
BHmin = np.min([i for i,elem in enumerate(power_predicted_df['BH']>0.8,1) if elem])+10
UNmin = np.min([i for i,elem in enumerate(power_predicted_df['UN']>0.8,1) if elem])+10

plt.figure(figsize=(7,5))
plt.plot([RFTmin,RFTmin],[0,power_predicted_df['RFT'][RFTmin-10]],color=colset1[2])
plt.plot([10,RFTmin],[power_predicted_df['RFT'][RFTmin-10],power_predicted_df['RFT'][RFTmin-10]],color=colset1[2])
plt.plot([BFmin,BFmin],[0,power_predicted_df['BF'][BFmin-10]],color=colset1[0])
plt.plot([10,BFmin],[power_predicted_df['BF'][BFmin-10],power_predicted_df['BF'][BFmin-10]],color=colset1[0])
plt.plot([BHmin,BHmin],[0,power_predicted_df['BH'][BHmin-10]],color=colset1[1])
plt.plot([10,BHmin],[power_predicted_df['BH'][BHmin-10],power_predicted_df['BH'][BHmin-10]],color=colset1[1])
plt.plot([UNmin,UNmin],[0,power_predicted_df['UN'][UNmin-10]],color=colset1[3])
plt.plot([10,UNmin],[power_predicted_df['UN'][UNmin-10],power_predicted_df['UN'][UNmin-10]],color=colset1[3])
plt.plot(newsubs,power_predicted_df['BF'],color=colset1[0],lw=2,label="Bonferroni")
plt.plot(newsubs,power_predicted_df['BH'],color=colset1[1],lw=2,label="Benjamini-Hochberg")
plt.plot(newsubs,power_predicted_df['RFT'],color=colset1[2],lw=2,label="Random Field Theory")
plt.plot(newsubs,power_predicted_df['UN'],color=colset1[3],lw=2,label="Uncorrected")
plt.text(RFTmin+1,0.02,str(RFTmin),color=colset1[2])
plt.text(BFmin-2.5,0.02,str(BFmin),color=colset1[0])
plt.text(UNmin-2.5,0.02,str(UNmin),color=colset1[3])
plt.text(BHmin-2.5,0.02,str(BHmin),color=colset1[1])
plt.ylim([0,1])
plt.title("Power with varying sample size")
plt.xlabel("Subjects")
plt.ylabel("Average power")
plt.legend(loc="center right",frameon=False)
plt.savefig("FIG_EX_power.pdf")
plt.close()
