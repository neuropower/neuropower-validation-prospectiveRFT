############
# preamble #
############

from __future__ import division
import os
import sys
import numpy as np
import scipy
import math
import nibabel as nib
import shutil
import csv
import imp
import subprocess
from neuropower import BUM, neuropowermodels, cluster
import model
import pandas
import uuid

exc = float(sys.argv[1])
pilot_sub = int(sys.argv[2])
final_sub = int(sys.argv[3])
seed = int(sys.argv[4])
FILEDIR = sys.argv[5]
DATADIR = sys.argv[6]
TMPDIR = sys.argv[7]
RESDIR = sys.argv[8]
ADAPTIVE = sys.argv[9]

TEMPDIR = os.path.join(TMPDIR,str(uuid.uuid4()))
os.mkdir(TEMPDIR)
os.chdir(TEMPDIR)

# parameters
all_sub = 180
true_sub = 100

resfile = os.path.join(RESDIR,"estimation_hcp_"+str(seed)+".csv")
if os.path.isfile(resfile):
    os.remove(resfile)

for c in range(47):
    # read mask, list with unique contrasts/paradigms for HCP, list of subject ID's
    maskfile = os.path.join(FILEDIR,'HCP_mask.nii.gz')
    mask = nib.load(maskfile).get_data()
    contrast = [line.rstrip('\n') for line in open(os.path.join(FILEDIR,'HCP_contrasts.txt'))][c]
    paradigm = [line.rstrip('\n') for line in open(os.path.join(FILEDIR,'HCP_paradigms.txt'))][c]
    subjects = [line.rstrip('\n') for line in open(os.path.join(FILEDIR,'HCP_unrelated_subID.txt'))]
    disks = [line.rstrip('\n') for line in open(os.path.join(FILEDIR,'HCP_unrelated_disk.txt'))]

    # select subjects for analysis
    all_subs = np.arange(all_sub)
    np.random.seed(seed=seed)
    suborder = np.random.choice(all_subs,len(all_subs),replace=False)
    if ADAPTIVE == "adaptive":
        pilot_subs = np.sort(suborder[0:pilot_sub])
        final_subs = np.sort(suborder[0:final_sub])
        true_subs = np.sort(suborder[final_sub:all_sub])
    elif ADAPTIVE == "predictive":
        pilot_subs = np.sort(suborder[0:pilot_sub])
        final_subs = np.sort(suborder[pilot_sub:(final_sub+pilot_sub)])
        true_subs = np.sort(suborder[(pilot_sub+final_sub):all_sub])

    ################################
    # select and analyze true data #
    ################################

    TRUEDIR = os.path.join(TEMPDIR,"analysis_true_"+str(seed)+"/")
    os.mkdir(TRUEDIR)
    os.chdir(TRUEDIR)

    true_cope = []
    for sub in true_subs:
        true_cope.append(os.path.join(DATADIR,"Disk"+str(disks[sub])+"of5",subjects[sub],"MNINonLinear/Results/tfMRI_"+paradigm,"tfMRI_"+paradigm+"_hp200_s4_level2vol.feat","cope"+contrast+".feat","stats/cope1.nii.gz"))

    cope_cmd = "fslmerge -t true_cope.nii.gz "+" ".join(true_cope)
    os.popen(cope_cmd).read()

    model.model(len(true_cope),TRUEDIR)
    fslcmd = 'flameo --copefile=true_cope.nii.gz --covsplitfile=design.grp --designfile=design.mat --ld=stats --maskfile=%s --runmode=ols --tcontrastsfile=design.con' %(os.path.join(FILEDIR,'HCP_mask.nii.gz'))
    os.popen(fslcmd).read()

    Res = np.sum(mask)/(5/2)**3
    zs = np.arange(exc,15,0.0001)
    pN_RFT = Res*np.exp(-zs**2/2)*zs**2
    cutoff_RFT = min(zs[pN_RFT<0.05])

    mask_cmd = 'fslmaths stats/zstat1.nii.gz -thr %s -bin thresh_zstat.nii.gz' %str(cutoff_RFT)
    os.popen(mask_cmd)

    #################################
    # select and analyze pilot data #
    #################################

    PILOTDIR = os.path.join(TEMPDIR,"analysis_pilot_"+str(seed)+"/")
    os.mkdir(PILOTDIR)
    os.chdir(PILOTDIR)

    pilot_cope = []
    for sub in pilot_subs:
        pilot_cope.append(os.path.join(DATADIR,"Disk"+str(disks[sub])+"of5",subjects[sub],"MNINonLinear/Results/tfMRI_"+paradigm,"tfMRI_"+paradigm+"_hp200_s4_level2vol.feat","cope"+contrast+".feat","stats/cope1.nii.gz"))

    cope_cmd = "fslmerge -t cope.nii.gz "+" ".join(pilot_cope)
    os.popen(cope_cmd).read()

    model.model(len(pilot_cope),PILOTDIR)
    fslcmd = 'flameo --copefile=cope.nii.gz --covsplitfile=design.grp --designfile=design.mat --ld=stats --maskfile=%s --runmode=ols --tcontrastsfile=design.con' %(os.path.join(FILEDIR,'HCP_mask.nii.gz'))
    os.popen(fslcmd).read()

    SPM = nib.load("stats/zstat1.nii.gz").get_data()

    peaks = cluster.cluster(SPM,exc)

    ###############################################################
    # estimate and compute model and estimate power on pilot data #
    ###############################################################

    # compute P-values

    pvalues = np.exp(-exc*(np.array(peaks.peak)-exc))
    pvalues = [max(10**(-6),t) for t in pvalues]
    peaks['pval'] = pvalues

    # estimate model

    bum = BUM.bumOptim(peaks['pval'].tolist(),starts=10)
    if bum['pi1'] == 0:
        est_eff = 'nan'
    else:
        modelfit = neuropowermodels.modelfit(peaks.peak,bum['pi1'],exc=exc,starts=10,method="RFT")
        est_eff = modelfit['mu']
        est_sd = modelfit['sigma']
        tau = neuropower.TruncTau(est_eff,est_sd,exc)
        est_exp_eff = est_eff + tau*est_sd

    # compute true parameters

    truefile = os.path.join(TRUEDIR,'thresh_zstat.nii.gz')
    activation = nib.load(truefile).get_data()

    truth = []
    for i in range(len(peaks)):
        peak_act = activation[peaks.x[i],peaks.y[i],peaks.z[i]]
        truth.append(peak_act)

    truth = [0 if x == 0 else 1 for x in truth]
    peaks['active'] = truth
    true_indices = [index for index,value in enumerate(truth) if value == 1]
    false_indices = [index for index,value in enumerate(truth) if value == 0]
    true_effectsize = np.mean(peaks.peak[true_indices])
    true_sd = np.std(peaks.peak[true_indices])
    true_pi1 = np.mean(truth)

    # correction pi0

    if np.sum(peaks['active'])==0:
        bumN = BUM.bumOptim(peaks['pval'][false_indices].tolist(),starts=100)
        pi0N = 1-bumN['pi1']
        PN0 = pi0N*np.sum(peaks['active']==0) #how many are null in the nonsignificant batch
        PN1 = (1-pi0N)*np.sum(peaks['active']==0)
        cor_pi1 = 1-PN0/len(peaks)
    elif np.sum(peaks['active'])==len(peaks):
        bumP = BUM.bumOptim(peaks['pval'][true_indices].tolist(),starts=100)
        pi0P = 1-bumP['pi1']
        PA0 = pi0P*np.sum(peaks['active']==1)
        PA1 = (1-pi0P)*np.sum(peaks['active']==1)
        cor_pi1 = 1-PA0/len(peaks)
    else:
        bumN = BUM.bumOptim(peaks['pval'][false_indices].tolist(),starts=100)
        pi0N = 1-bumN['pi1']
        bumP = BUM.bumOptim(peaks['pval'][true_indices].tolist(),starts=100)
        pi0P = 1-bumP['pi1']
        PN0 = pi0N*np.sum(peaks['active']==0)
        PA0 = pi0P*np.sum(peaks['active']==1)
        PN1 = (1-pi0N)*np.sum(peaks['active']==0)
        PA1 = (1-pi0P)*np.sum(peaks['active']==1)
        cor_pi1 = 1-(PN0+PA0)/len(peaks)

    # correction effect
    if np.sum(peaks['active'])==0:
        EN = np.mean(peaks['peak'][false_indices])
        EN1 = (EN - pi0N*(exc+1./exc))/(1-pi0N)
        cor_effectsize = (PN1*EN1)/(PN1)
    elif np.sum(peaks['active'])==len(peaks):
        EA = np.mean(peaks['peak'][true_indices])
        EA1 = (EA - pi0P*(exc+1/exc))/(1-pi0P)
        cor_effectsize = (PA1*EA1)/(PA1)
    else:
        EN = np.mean(peaks['peak'][false_indices])
        EA = np.mean(peaks['peak'][true_indices])
        EN1 = (EN - pi0N*(exc+1/exc))/(1-pi0N)
        EA1 = (EA - pi0P*(exc+1/exc))/(1-pi0P)
        cor_effectsize = (PN1*EN1 + PA1*EA1)/(PN1+PA1)

    # write away estimation results with true values

    if bum['pi1'] == 0:
        est_exp_eff = 'nan'
        est_sd = 'nan'

    estimation = [c, bum['pi1'],true_pi1,cor_pi1,est_eff,true_effectsize,cor_effectsize,est_exp_eff,est_sd,true_sd]
    fd = open(resfile,"a")
    wr = csv.writer(fd,quoting=csv.QUOTE_NONE)
    wr.writerow(estimation)
    fd.close()

    if bum['pi1'] == 0:
        shutil.rmtree(TRUEDIR)
        shutil.rmtree(PILOTDIR)
        continue

    # predict power

    thresholds = neuropower.threshold(peaks.peak,peaks.pval,FWHM=[2.5,2.5,2.5],voxsize=[1,1,1],nvox=np.product(SPM.size),alpha=0.05,exc=exc,method="RFT")
    effect_cohen = modelfit['mu']/np.sqrt(pilot_sub)
    power_predicted = []
    for s in range(pilot_sub,final_sub):
        projected_effect = effect_cohen*np.sqrt(s)
        powerpred =  {k:1-neuropower.altCDF(v,projected_effect,modelfit['sigma'],exc=exc,method="RFT") for k,v in thresholds.items() if v!='nan'}
        power_predicted.append(powerpred)

    shutil.rmtree(PILOTDIR)

    ######################
    # analyze final data #
    ######################

    power_true = []
    for s in range(pilot_sub,final_sub):
        #analyze data
        FINALDIR = os.path.join(TEMPDIR,"analysis_final_"+str(seed)+"/")
        os.mkdir(FINALDIR)
        os.chdir(FINALDIR)
        final_cope = []
        for sub in final_subs[range(s)]:
            final_cope.append(os.path.join(DATADIR,"Disk"+str(disks[sub])+"of5",subjects[sub],"MNINonLinear/Results/tfMRI_"+paradigm,"tfMRI_"+paradigm+"_hp200_s4_level2vol.feat","cope"+contrast+".feat","stats/cope1.nii.gz"))

        cope_cmd = "fslmerge -t cope.nii.gz "+" ".join(final_cope)
        os.popen(cope_cmd).read()

        model.model(len(final_cope),FINALDIR)
        fslcmd = 'flameo --copefile=cope.nii.gz --covsplitfile=design.grp --designfile=design.mat --ld=stats --maskfile=%s --runmode=ols --tcontrastsfile=design.con' %(os.path.join(FILEDIR,'HCP_mask.nii.gz'))
        os.popen(fslcmd).read()

        SPM = nib.load("stats/zstat1.nii.gz").get_data()
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
        thresholds = neuropower.threshold(peaks.peak,peaks.pval,FWHM=[2.5,2.5,2.5],voxsize=[1,1,1],nvox=np.product(SPM.size),alpha=0.05,exc=exc,method="RFT")
        TPR = {'UN':'nan','BF':'nan','RFT':'nan','BH':'nan'}
        for method in range(4):
            ind = ["UN","BF","RFT","BH"][method]
            if thresholds[ind] == 'nan':
                continue
            pos = peaks.peak>thresholds[ind]
            true = peaks.active==1
            TP = [a and b for a,b in zip(pos,true)]
            peaks['TP'] = TP
            if np.sum(true) == 0:
                continue
            TPR[ind] = float(np.sum(TP))/float(np.sum(true))
        power_true.append(TPR)
        shutil.rmtree(FINALDIR)

    toCSV = power_predicted
    keys = toCSV[0].keys()
    with open(os.path.join(RESDIR,'powpre_hcp_'+str(seed)+'_contrast_'+str(c)+'.csv'),'wb') as output_file:
        dict_writer = csv.DictWriter(output_file, keys)
        dict_writer.writeheader()
        dict_writer.writerows(toCSV)

    toCSV = power_true
    keys = toCSV[0].keys()
    with open(os.path.join(RESDIR,'powtru_hcp_'+str(seed)+'_contrast_'+str(c)+'.csv'),'wb') as output_file:
        dict_writer = csv.DictWriter(output_file, keys)
        dict_writer.writeheader()
        dict_writer.writerows(toCSV)

    shutil.rmtree(TRUEDIR)
shutil.rmtree(TEMPDIR)
