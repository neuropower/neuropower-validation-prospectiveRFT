#### THIS MAKES NO SENSE

############
# preamble #
############

from __future__ import division
import os
import sys
sys.path.append("/share/PI/russpold/software/anaconda/lib/python2.7/site-packages")
pilot_sub = int(sys.argv[1])
final_sub = int(sys.argv[2])
seed = int(sys.argv[3])
HOMEDIR = sys.argv[4]
DATADIR = sys.argv[5]
TMPDIR = sys.argv[6]
os.chdir(HOMEDIR)

import numpy as np
import scipy
import math
import nibabel as nib
import shutil
import csv
import sys
import imp
import BUM
import neuropower
import cluster
import peakdistribution
import simul_multisubject_fmri_dataset
import model
import matplotlib.pyplot as plt
import pandas as pd

# results file

results = pd.DataFrame(columns=['eff','eff_est','sd','est_sd','pi1','est_pi1'])

# standard parameters
exc = 2.3
pi1 = 0.5
nopeaks = 500
eff = 4
sd = 1
sim = 1

# define distributions for these parameters
R = np.arange(exc,7+0.001,0.001)
pdf_null = neuropower.nulPDF(R,exc=exc)
pdf_alt = neuropower.altPDF(R,eff,sd,exc=exc)
cdf_null = neuropower.nulCDF(R,exc=exc)
cdf_alt = neuropower.altCDF(R,eff,sd,exc=exc)
pdf_mix = neuropower.mixprobdens(R,pi1,eff,sd,exc=exc)

# sample at random from distributions
samples = list()
truth = list()
for i in range(nopeaks):
    n = np.random.random()
    if i < (pi1*nopeaks):
        value = R[np.min(np.where(n<cdf_null))]
    else:
        value = R[np.min(np.where(n<cdf_alt))]
    samples.append(value)
    truth.append(0 if i < (pi1*nopeaks) else 1)

# make plot
#n,bins,patches = plt.hist(samples,20,normed=1,facecolor="black",alpha=0.5)
#plt.plot(R,pdf_mix,color="red",lw=3)
#ax = plt.axes()
#ax.set_xlim(exc,7)
#ax.set_ylim(0,1.05)
#plt.show()

# compute P-values

pvalues = np.exp(-exc*(np.array(samples)-exc))
pvalues = [max(10**(-6),t) for t in pvalues]

# estimate model

bum = BUM.bumOptim(pvalues,starts=10)
if bum['pi1'] == 0:
    est_eff = 'nan'
else:
    modelfit = neuropower.modelfit(samples,bum['pi1'],exc=exc,starts=20,method="RFT")
    est_eff = modelfit['mu']
    est_sd = modelfit['sigma']
    #tau = neuropower.TruncTau(est_eff,est_sd,exc)
    #est_exp_eff = est_eff + tau*est_sd

# write away estimation results with true values


estimation = {'eff':eff,
                'est_eff':est_eff,
                'sd':sd,
                'est_sd':est_sd,
                'pi1':pi1,
                'est_pi1':bum['pi1']
                }
results.loc[sim] = pd.Series(estimation)




estimation = [effectsize, str(wd_names[c]), bum['pi1'],true_pi1,est_eff,true_effectsize,est_exp_eff,est_sd,true_sd,bum['a']]
fd = open(resfile,"a")
wr = csv.writer(fd,quoting=csv.QUOTE_NONE)
wr.writerow(estimation)
fd.close()

if bum['pi1'] == 0:
    shutil.rmtree(TEMPDIR)
    continue

# predict power

thresholds = neuropower.threshold(peaks.peak,peaks.pval,FWHM=FWHM,mask=mask,alpha=0.05,exc=exc)
effect_cohen = modelfit['mu']/np.sqrt(pilot_sub)
power_predicted = []
for s in range(pilot_sub,final_sub):
    projected_effect = effect_cohen*np.sqrt(s)
    powerpred =  {k:1-neuropower.altCDF(v,projected_effect,modelfit['sigma'],exc=exc,method="RFT") for k,v in thresholds.items() if v!='nan'}
    power_predicted.append(powerpred)

shutil.rmtree(TEMPDIR)

#######################
# simulate final data #
#######################

power_true = []
for s in range(pilot_sub,final_sub):

    #simulate data
    os.mkdir(TEMPDIR)
    os.chdir(TEMPDIR)
    new_subj=s

    signal = nib.load(os.path.join(HOMEDIR,"SIM_activation"+str(width)+".nii")).get_data()
    signal = np.repeat(signal[:, :, :, np.newaxis], new_subj, axis=3)
    low_values_indices = signal < 0.1
    signal[low_values_indices] = 0
    high_values_indices = signal > 0
    signal[high_values_indices] = effectsize

    noise_final = simul_multisubject_fmri_dataset.surrogate_3d_dataset(n_subj=new_subj,mask=mask,sk=smooth_sigma,noise_level=1.0,width=5.0,out_text_file=None,out_image_file=None, seed=seed)
    noise_final = noise_final.transpose((1,2,3,0))
    noise = noise_final

    data = signal + noise
    img=nib.Nifti1Image(data,np.eye(4))
    img.to_filename(os.path.join("simulation.nii.gz"))

    # analyse data and extract peaks
    model.model(new_subj,TEMPDIR)
    fslcmd = 'flameo --copefile=simulation.nii.gz --covsplitfile=design.grp --designfile=design.mat --ld=stats --maskfile=%s --runmode=ols --tcontrastsfile=design.con' %(os.path.join(HOMEDIR,'SIM_mask.nii'))
    os.popen(fslcmd).read()

    SPM = nib.load("stats/zstat1.nii.gz").get_data()
    SPM = SPM[::-1,:,:]
    peaks = cluster.cluster(SPM,exc)
    pvalues = np.exp(-exc*(np.array(peaks.peak)-exc))
    pvalues = [max(10**(-6),t) for t in pvalues]
    peaks['pval'] = pvalues

    # compute true power for different procedures

    truth = []
    for i in range(len(peaks)):
        peak_act = activation[peaks.x[i],peaks.y[i],peaks.z[i]]
        truth.append(peak_act)

    truth = [0 if x == 0 else 1 for x in truth]
    peaks['active'] = truth
    thresholds = neuropower.threshold(peaks.peak,peaks.pval,FWHM=FWHM,mask=mask,alpha=0.05,exc=exc)
    res = {
        'UN_TP':'nan','BF_TP':'nan','RFT_TP':'nan','BH_TP':'nan',
        'UN_FP':'nan','BF_FP':'nan','RFT_FP':'nan','BH_FP':'nan',
        'UN_FN':'nan','BF_FN':'nan','RFT_FN':'nan','BH_FN':'nan',
        'UN_TN':'nan','BF_TN':'nan','RFT_TN':'nan','BH_TN':'nan'
        }
    for method in range(4):
        ind = ["UN","BF","RFT","BH"][method]
        if thresholds[ind] == 'nan':
            continue
        if np.sum(peaks.peak>exc) == 0:
            TP = FP = TN = FN = 0
        else:
            pos = peaks.peak>thresholds[ind]
            true = peaks.active==1
            TP = np.sum([a and b for a,b in zip(pos,true)])
            FP = np.sum([a and not b for a,b in zip(pos,true)])
            FN = np.sum([b and not a for a,b in zip(pos,true)])
            TN = np.sum([not a and not b for a,b in zip(pos,true)])
        res[str(ind+"_TP")] = TP
        res[str(ind+"_FP")] = FP
        res[str(ind+"_TN")] = TN
        res[str(ind+"_FN")] = FN
    power_true.append(res)
    shutil.rmtree(TEMPDIR)

# write away data

toCSV = power_predicted
keys = toCSV[0].keys()
with open(os.path.join(DATADIR,'powpre_sim_'+str(seed)+'_w_'+str(wd_names[c])+'_e_'+es_names[c]+'.csv'),'wb') as output_file:
    dict_writer = csv.DictWriter(output_file, keys)
    dict_writer.writeheader()
    dict_writer.writerows(toCSV)

toCSV = power_true
keys = toCSV[0].keys()
with open(os.path.join(DATADIR,'powtru_sim_'+str(seed)+'_w_'+str(wd_names[c])+'_e_'+es_names[c]+'.csv'),'wb') as output_file:
    dict_writer = csv.DictWriter(output_file, keys)
    dict_writer.writeheader()
    dict_writer.writerows(toCSV)
