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
sys.path.append("/home1/03545/jdurnez/power/Functions/")
import model
import pandas
import uuid

EXC = float(sys.argv[1])
PILOT = int(sys.argv[2])
FINAL = int(sys.argv[3])
SEED = int(sys.argv[4])
ADAPTIVE = sys.argv[5]
MODEL = sys.argv[6]
startloop = int(sys.argv[7])
endloop = int(sys.argv[8])
TMPDIR = sys.argv[9]

FILEDIR = os.environ.get('FILEDIR')
HCPDIR = os.environ.get('HCPDIR')
OUTDIR = os.environ.get('OUTDIR')

TEMPDIR = TMPDIR

# parameters
all_sub = 180
true_sub = 100

for c in np.arange(startloop,endloop):

    resfile = os.path.join(OUTDIR,"estimation_HCP_"+str(SEED)+"_"+str(startloop)+"-"+str(endloop)+".csv")
    PredictionFile = os.path.join(OUTDIR,'powpre_HCP_'+str(SEED)+'_contrast_'+str(c)+'.csv')
    TrueFile = os.path.join(OUTDIR,'powtru_HCP_'+str(SEED)+'_contrast_'+str(c)+'.csv')
    # FileList = [resfile,PredictionFile,TrueFile]
    #
    # # Check if results files exist:
    #     # if sample is done => continue
    #     # if sample half done => remove files and restart
    #     # if not done => pass
    # exist = [os.path.isfile(x) for x in FileList]
    # if np.sum(exist) == 3:
    #     continue
    # elif np.sum(exist) == 0:
    #     pass
    # else:
    #     [os.remove(x) if os.path.isfile(x) else None for x in FileList]

    # read mask, list with unique contrasts/paradigms for HCP, list of subject ID's
    maskfile = os.path.join(FILEDIR,'HCP_mask.nii.gz')
    mask = nib.load(maskfile).get_data()
    contrast = [line.rstrip('\n') for line in open(os.path.join(FILEDIR,'HCP_contrasts.txt'))][c]
    paradigm = [line.rstrip('\n') for line in open(os.path.join(FILEDIR,'HCP_paradigms.txt'))][c]
    subjects = [line.rstrip('\n') for line in open(os.path.join(FILEDIR,'HCP_unrelated_subID.txt'))]
    disks = [line.rstrip('\n') for line in open(os.path.join(FILEDIR,'HCP_unrelated_disk.txt'))]

    # select subjects for analysis
    all_subs = np.arange(all_sub)
    np.random.seed(seed=SEED)
    suborder = np.random.choice(all_subs,len(all_subs),replace=False)
    if ADAPTIVE == "adaptive":
        pilot_subs = np.sort(suborder[0:PILOT])
        final_subs = np.sort(suborder[0:FINAL])
        true_subs = np.sort(suborder[FINAL:all_sub])
    elif ADAPTIVE == "predictive":
        pilot_subs = np.sort(suborder[0:PILOT])
        final_subs = np.sort(suborder[PILOT:(FINAL+PILOT)])
        true_subs = np.sort(suborder[(PILOT+FINAL):all_sub])

    ################################
    # select and analyze true data #
    ################################

    TRUEDIR = os.path.join(TEMPDIR,"analysis_true_"+str(SEED)+"/")
    os.popen("mkdir "+str(TRUEDIR))
    os.chdir(TRUEDIR)

    true_cope = []
    for sub in true_subs:
        true_cope.append(os.path.join(HCPDIR,"disk"+str(disks[sub]),subjects[sub],"MNINonLinear/Results/tfMRI_"+paradigm,"tfMRI_"+paradigm+"_hp200_s4_level2vol.feat","cope"+contrast+".feat","stats/cope1.nii.gz"))

    cope_cmd = "fslmerge -t true_cope.nii.gz "+" ".join(true_cope)
    os.popen(cope_cmd).read()

    model.model(len(true_cope),TRUEDIR)
    fslcmd = 'flameo --copefile=true_cope.nii.gz --covsplitfile=design.grp --designfile=design.mat --ld=stats --maskfile=%s --runmode=ols --tcontrastsfile=design.con' %(os.path.join(FILEDIR,'HCP_mask.nii.gz'))
    os.popen(fslcmd).read()

    Res = np.sum(mask)/(5/2)**3
    zs = np.arange(EXC,15,0.0001)
    pN_RFT = Res*np.exp(-zs**2/2)*zs**2
    cutoff_RFT = min(zs[pN_RFT<0.05])

    mask_cmd = 'fslmaths stats/zstat1.nii.gz -thr %s -bin thresh_zstat.nii.gz' %str(cutoff_RFT)
    os.popen(mask_cmd)

    #################################
    # select and analyze pilot data #
    #################################

    PILOTDIR = os.path.join(TEMPDIR,"analysis_pilot_"+str(SEED)+"/")
    os.popen("mkdir "+str(PILOTDIR))
    os.chdir(PILOTDIR)

    pilot_cope = []
    for sub in pilot_subs:
        pilot_cope.append(os.path.join(HCPDIR,"disk"+str(disks[sub]),subjects[sub],"MNINonLinear/Results/tfMRI_"+paradigm,"tfMRI_"+paradigm+"_hp200_s4_level2vol.feat","cope"+contrast+".feat","stats/cope1.nii.gz"))

    cope_cmd = "fslmerge -t cope.nii.gz "+" ".join(pilot_cope)
    os.popen(cope_cmd).read()

    model.model(len(pilot_cope),PILOTDIR)
    fslcmd = 'flameo --copefile=cope.nii.gz --covsplitfile=design.grp --designfile=design.mat --ld=stats --maskfile=%s --runmode=ols --tcontrastsfile=design.con' %(os.path.join(FILEDIR,'HCP_mask.nii.gz'))
    os.popen(fslcmd).read()

    SPM = nib.load("stats/zstat1.nii.gz").get_data()

    if MODEL == "RFT":
        peaks = cluster.PeakTable(SPM,EXC,mask)
    elif MODEL == "CS":
        peaks = cluster.PeakTable(SPM,-100,mask)

    ###############################################################
    # estimate and compute model and estimate power on pilot data #
    ###############################################################

    # compute P-values

    if MODEL == "RFT":
        pvalues = np.exp(-EXC*(np.array(peaks.peak)-EXC))
        pvalues = [max(10**(-6),t) for t in pvalues]
    elif MODEL == "CS":
        pvalues = 1-np.asarray(neuropowermodels.nulCDF(peaks.peak,method="CS"))
    peaks['pval'] = pvalues

    # estimate model

    bum = BUM.EstimatePi1(peaks['pval'].tolist(),starts=10)
    if bum['pi1'] == 0:
        est_eff = 'nan'
    else:
        if MODEL == "RFT":
            modelfit = neuropowermodels.modelfit(peaks.peak,bum['pi1'],exc=EXC,starts=20,method="RFT")
            est_eff = modelfit['mu']
            est_sd = modelfit['sigma']
            tau = neuropowermodels.TruncTau(est_eff,est_sd,EXC)
            est_exp_eff = est_eff + tau*est_sd
            mu = modelfit['mu']
        elif MODEL == "CS":
            modelfit = neuropowermodels.modelfit(peaks.peak,bum['pi1'],starts=5,method="CS")
            est_sd = 'nan'
            xn = np.arange(-10,30,0.01)
            alt = np.asarray(neuropowermodels.altPDF(xn,mu=modelfit['mu'],method="CS"))
            est_eff = xn[alt==np.max(alt)][0]
            est_exp_eff = 'nan'
            mu = modelfit['mu']

    # compute true parameters

    truefile = os.path.join(TRUEDIR,'thresh_zstat.nii.gz')
    activation = nib.load(truefile).get_data()

    truth = []
    for i in range(len(peaks)):
        peak_act = activation[int(peaks.x[i]),int(peaks.y[i]),int(peaks.z[i])]
        truth.append(peak_act)

    truth = [0 if x == 0 else 1 for x in truth]
    peaks['active'] = truth
    true_indices = [index for index,value in enumerate(truth) if value == 1]
    false_indices = [index for index,value in enumerate(truth) if value == 0]
    true_effectsize = np.mean(peaks.peak[true_indices])
    true_sd = np.std(peaks.peak[true_indices])
    true_pi1 = np.mean(truth)

    # correction pi0
    pi0P=0
    if np.sum(peaks['active'])==0:
        bumN = BUM.EstimatePi1(peaks['pval'][false_indices].tolist(),starts=100)
        pi0N = 1-bumN['pi1']
        PN0 = pi0N*np.sum(peaks['active']==0) #how many are null in the nonsignificant batch
        PN1 = (1-pi0N)*np.sum(peaks['active']==0)
        cor_pi1 = 1-PN0/len(peaks)
    elif np.sum(peaks['active'])==len(peaks):
        #bumP = BUM.EstimatePi1(peaks['pval'][true_indices].tolist(),starts=100)
        #pi0P = 1-bumP['pi1']
        PA0 = pi0P*np.sum(peaks['active']==1)
        PA1 = (1-pi0P)*np.sum(peaks['active']==1)
        cor_pi1 = 1-PA0/len(peaks)
    else:
        bumN = BUM.EstimatePi1(peaks['pval'][false_indices].tolist(),starts=100)
        pi0N = 1-bumN['pi1']
        #bumP = BUM.EstimatePi1(peaks['pval'][true_indices].tolist(),starts=100)
        #pi0P = 1-bumP['pi1']
        PN0 = pi0N*np.sum(peaks['active']==0)
        PA0 = pi0P*np.sum(peaks['active']==1)
        PN1 = (1-pi0N)*np.sum(peaks['active']==0)
        PA1 = (1-pi0P)*np.sum(peaks['active']==1)
        cor_pi1 = 1-(PN0+PA0)/len(peaks)

    # correction effect
    if np.sum(peaks['active'])==0:
        EN = np.mean(peaks['peak'][false_indices])
        EN1 = (EN - pi0N*(EXC+1./EXC))/(1-pi0N)
        cor_effectsize = (PN1*EN1)/(PN1)
    elif np.sum(peaks['active'])==len(peaks):
        EA = np.mean(peaks['peak'][true_indices])
        EA1 = (EA - pi0P*(EXC+1/EXC))/(1-pi0P)
        cor_effectsize = (PA1*EA1)/(PA1)
    else:
        EN = np.mean(peaks['peak'][false_indices])
        EA = np.mean(peaks['peak'][true_indices])
        EN1 = (EN - pi0N*(EXC+1/EXC))/(1-pi0N)
        EA1 = (EA - pi0P*(EXC+1/EXC))/(1-pi0P)
        cor_effectsize = (PN1*EN1 + PA1*EA1)/(PN1+PA1)

    # write away estimation results with true values

    if bum['pi1'] == 0:
        est_exp_eff = 'nan'
        est_sd = 'nan'

    estimation = [c, bum['pi1'],true_pi1,cor_pi1,est_eff,true_effectsize,cor_effectsize,est_exp_eff,est_sd,true_sd,bum['a']]
    fd = open(resfile,"a")
    wr = csv.writer(fd,quoting=csv.QUOTE_NONE)
    if c == 0:
        estnames = ['contrast','pi1e','pi1t','pi1c','ese','est','esc','esexp','sde','sdt','bumpar']
        wr.writerow(estnames)
    wr.writerow(estimation)
    fd.close()

    if bum['pi1'] == 0:
        shutil.rmtree(TRUEDIR)
        shutil.rmtree(PILOTDIR)
        continue

    # predict power

    thresholds = neuropowermodels.threshold(peaks.peak,peaks.pval,FWHM=[2.5,2.5,2.5],voxsize=[1,1,1],nvox=np.product(SPM.size),alpha=0.05,exc=EXC,method="RFT")
    if MODEL=="CS":
        thresholds["UN"] = 3.521390
    effect_cohen = modelfit['mu']/np.sqrt(PILOT)
    power_predicted = []
    for s in range(PILOT,FINAL):
        projected_effect = effect_cohen*np.sqrt(s)
        if MODEL == "RFT":
            powerpred =  {k:1-neuropowermodels.altCDF(v,projected_effect,modelfit['sigma'],exc=EXC,method="RFT") for (k,v) in thresholds.items() if v!='nan'}
        elif MODEL == "CS":
            xn = np.arange(-10,30,0.01)
            nul = np.asarray(neuropowermodels.nulPDF(xn,method="CS"))
            projected_effect = projected_effect-xn[nul==np.max(nul)][0]
            powerpred =  {k:[1-neuropowermodels.altCDF([v],projected_effect,method="CS")[0] for (k,v) in thresholds.items() if v!='nan']}
        power_predicted.append(powerpred)

    shutil.rmtree(PILOTDIR)

    ######################
    # analyze final data #
    ######################

    power_true = []
    for s in range(PILOT,FINAL):
        #analyze data
        FINALDIR = os.path.join(TEMPDIR,"analysis_final_"+str(SEED)+"/")
        os.popen("mkdir "+str(FINALDIR))
        os.chdir(FINALDIR)
        final_cope = []
        for sub in final_subs[range(s)]:
            final_cope.append(os.path.join(HCPDIR,"disk"+str(disks[sub]),subjects[sub],"MNINonLinear/Results/tfMRI_"+paradigm,"tfMRI_"+paradigm+"_hp200_s4_level2vol.feat","cope"+contrast+".feat","stats/cope1.nii.gz"))

        cope_cmd = "fslmerge -t cope.nii.gz "+" ".join(final_cope)
        os.popen(cope_cmd).read()

        model.model(len(final_cope),FINALDIR)
        fslcmd = 'flameo --copefile=cope.nii.gz --covsplitfile=design.grp --designfile=design.mat --ld=stats --maskfile=%s --runmode=ols --tcontrastsfile=design.con' %(os.path.join(FILEDIR,'HCP_mask.nii.gz'))
        os.popen(fslcmd).read()

        SPM = nib.load("stats/zstat1.nii.gz").get_data()
        if MODEL == "RFT":
            peaks = cluster.PeakTable(SPM,EXC,mask)
            pvalues = np.exp(-EXC*(np.array(peaks.peak)-EXC))
            pvalues = [max(10**(-6),t) for t in pvalues]
        elif MODEL == "CS":
            peaks = cluster.PeakTable(SPM,-100,mask)
            pvalues = 1-np.asarray(neuropowermodels.nulCDF(peaks.peak,method="CS"))
        peaks['pval'] = pvalues

        # compute true power for different procedures
        truth = []
        for i in range(len(peaks)):
            peak_act = activation[int(peaks.x[i]),int(peaks.y[i]),int(peaks.z[i])]
            truth.append(peak_act)
        truth = [0 if x == 0 else 1 for x in truth]
        peaks['active'] = truth
        thresholds = neuropowermodels.threshold(peaks.peak,peaks.pval,FWHM=[2.5,2.5,2.5],voxsize=[1,1,1],nvox=np.product(SPM.size),alpha=0.05,exc=EXC,method="RFT")
        if MODEL == "CS":
            thresholds["UN"] = 3.521390
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
    with open(PredictionFile,'wb') as output_file:
        dict_writer = csv.DictWriter(output_file, keys)
        dict_writer.writeheader()
        dict_writer.writerows(toCSV)

    toCSV = power_true
    keys = toCSV[0].keys()
    with open(TrueFile,'wb') as output_file:
        dict_writer = csv.DictWriter(output_file, keys)
        dict_writer.writeheader()
        dict_writer.writerows(toCSV)

    shutil.rmtree(TRUEDIR)
shutil.rmtree(TEMPDIR)
