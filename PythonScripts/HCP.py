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
sys.path.append(os.path.join(os.environ.get('HOMEDIR'),"Functions/"))
from neuropower_local import utils,poweranalysis,neuropowermodels,effectsize
from neuropower_local.utils import cluster,peakpvalues
import model
import pandas
import uuid
from scipy import integrate
import pandas as pd

EXC = float(sys.argv[1])
PILOT = int(sys.argv[2])
FINAL = int(sys.argv[3])
SEED = int(sys.argv[4])
ADAPTIVE = sys.argv[5]
MODEL = sys.argv[6]
startloop = int(sys.argv[7])
endloop = int(sys.argv[8])+1
TMPDIR = sys.argv[9]

FILEDIR = os.environ.get('FILEDIR')
HCPDIR = os.environ.get('HCPDIR')
OUTDIR = os.environ.get('OUTDIR')

TEMPDIR = TMPDIR

# parameters
all_sub = 180
true_sub = 100
print("startloop: "+str(startloop))

for c in np.arange(startloop,endloop):
    print("contrast tested: "+str(c))

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
        true_cope.append(os.path.join(HCPDIR,"Disk"+str(disks[sub])+"of5",subjects[sub],"MNINonLinear/Results/tfMRI_"+paradigm,"tfMRI_"+paradigm+"_hp200_s4_level2vol.feat","cope"+contrast+".feat","stats/cope1.nii.gz"))

    cope_cmd = "fslmerge -t true_cope.nii.gz "+" ".join(true_cope)
    os.popen(cope_cmd).read()

    model.model(len(true_cope),TRUEDIR)
    fslcmd = 'flameo --copefile=true_cope.nii.gz --covsplitfile=design.grp --designfile=design.mat --ld=stats --maskfile=%s --runmode=ols --tcontrastsfile=design.con --outputdof' %(os.path.join(FILEDIR,'HCP_mask.nii.gz'))
    os.popen(fslcmd).read()

    smoothest = os.popen('$FSLDIR/bin/smoothest -d $(cat stats/dof) -m stats/mask -r stats/res4d').read()
    VOLUMEcount = str(smoothest.split("\n")[1].split("VOLUME ")[1])
    RESELcount = str(smoothest.split("\n")[2].split("RESELS ")[1])
    Res = float(VOLUMEcount)/float(RESELcount)
    #Res = int(np.sum(mask)/(5/2)**3)
    cutoff_RFT = float(os.popen("ptoz 0.05 -g %d"%Res).read().split("\n")[0])

    # cmd = 'fdr -i stats/pstat1.nii.gz -m %s -q 0.05'%(os.path.join(FILEDIR,'HCP_mask.nii.gz'))
    # cutoff_FDR = float(os.popen(cmd).read().split("\n")[1])
    # import scipy.stats as st
    # st.norm.ppf(1-cutoff_FDR)

    # from locfdr import locfdr
    # zstat = nib.load("stats/zstat1.nii.gz").get_data().flatten()
    # zstat = zstat[zstat!=0]
    # res = locfdr.locfdr(zstat)

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
        pilot_cope.append(os.path.join(HCPDIR,"Disk"+str(disks[sub])+"of5",subjects[sub],"MNINonLinear/Results/tfMRI_"+paradigm,"tfMRI_"+paradigm+"_hp200_s4_level2vol.feat","cope"+contrast+".feat","stats/cope1.nii.gz"))

    cope_cmd = "fslmerge -t cope.nii.gz "+" ".join(pilot_cope)
    os.popen(cope_cmd).read()

    model.model(len(pilot_cope),PILOTDIR)
    fslcmd = 'flameo --copefile=cope.nii.gz --covsplitfile=design.grp --designfile=design.mat --ld=stats --maskfile=%s --runmode=ols --tcontrastsfile=design.con' %(os.path.join(FILEDIR,'HCP_mask.nii.gz'))
    os.popen(fslcmd).read()

    SPM = nib.load("stats/zstat1.nii.gz").get_data()

    ###############################################################
    # estimate and compute model and estimate power on pilot data #
    ###############################################################

    power = poweranalysis.power(spm=SPM,mask=mask,FWHM=[2.5,2.5,2.5],voxsize=1,alpha=0.05,samplesize=PILOT)
    power.estimate_model(exc=EXC)
    peaks = power.peaktable

    est_eff = power.mu
    est_sd = power.sigma
    pi1e = power.pi1
    tau = neuropowermodels.TruncTau(est_eff,est_sd,EXC)
    est_exp_eff = est_eff + tau*est_sd

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

    peaks_supra = peaks[peaks.peak>EXC]
    true_effectsize = np.mean(peaks_supra.peak[peaks_supra.active==1])
    true_sd = np.std(peaks_supra.peak[peaks_supra.active==1])
    true_pi1 = np.mean(peaks_supra.active)

    # correction pi0
    pi0P=0
    if np.sum(peaks_supra['active'])==0:
        pi1 = effectsize.Pi1()
        pi1.pvalues = peaks_supra['pvals']
        pi1.estimate(starts=10)
        pi0N = 1-pi1.pi1
        PN0 = pi0N*np.sum(peaks_supra['active']==0) #how many are null in the nonsignificant batch
        PN1 = (1-pi0N)*np.sum(peaks_supra['active']==0)
        cor_pi1 = 1-PN0/len(peaks_supra)
    elif np.sum(peaks_supra['active'])==len(peaks_supra):
        #bumP = BUM.EstimatePi1(peaks['pval'][true_indices].tolist(),starts=100)
        #pi0P = 1-bumP['pi1']
        PA0 = pi0P*np.sum(peaks_supra['active']==1)
        PA1 = (1-pi0P)*np.sum(peaks_supra['active']==1)
        cor_pi1 = 1-PA0/len(peaks_supra)
    else:
        pi1 = effectsize.Pi1()
        pi1.pvalues = peaks_supra['pvals'][peaks_supra['active']==0]
        pi1.estimate(starts=10)
        pi0N = 1-pi1.pi1
        #bumP = BUM.EstimatePi1(peaks['pval'][true_indices].tolist(),starts=100)
        #pi0P = 1-bumP['pi1']
        PN0 = pi0N*np.sum(peaks_supra['active']==0)
        PA0 = pi0P*np.sum(peaks_supra['active']==1)
        PN1 = (1-pi0N)*np.sum(peaks_supra['active']==0)
        PA1 = (1-pi0P)*np.sum(peaks_supra['active']==1)
        cor_pi1 = 1-(PN0+PA0)/len(peaks_supra)

    # correction effect
    pi1all = effectsize.Pi1()
    pi1all.pvalues = peaks_supra['pvals']
    pi1all.estimate(starts=10)
    pi0N = 1-pi1.pi1
    if np.sum(peaks_supra['active'])==0:
        EN = np.mean(peaks_supra['peak'][peaks_supra['active']==0])
        EN1 = (EN - pi0N*(EXC+1./EXC))/(1-pi0N)
        cor_effectsize = (PN1*EN1)/(PN1)
    elif np.sum(peaks_supra['active'])==len(peaks_supra):
        EA = np.mean(peaks_supra['peak'][peaks_supra['active']==1])
        EA1 = (EA - pi0P*(EXC+1/EXC))/(1-pi0P)
        cor_effectsize = (PA1*EA1)/(PA1)
    else:
        EN = np.mean(peaks_supra['peak'][peaks_supra['active']==0])
        EA = np.mean(peaks_supra['peak'][peaks_supra['active']==1])
        EN1 = (EN - pi0N*(EXC+1/EXC))/(1-pi0N)
        EA1 = (EA - pi0P*(EXC+1/EXC))/(1-pi0P)
        cor_effectsize = (PN1*EN1 + PA1*EA1)/(PN1+PA1)

    # write away estimation results with true values

    if pi1e == 0:
        est_exp_eff = 'nan'
        est_sd = 'nan'
        est_eff = 'nan'

    estimation = [c,pi1e,true_pi1,cor_pi1,est_eff,true_effectsize,cor_effectsize,est_exp_eff,est_sd,true_sd,'cohen']
    fd = open(resfile,"a")
    wr = csv.writer(fd,quoting=csv.QUOTE_NONE)
    if c == 0:
        estnames = ['contrast','pi1e','pi1t','pi1c','ese','est','esc','esexp','sde','sdt','cohen']
        wr.writerow(estnames)
    wr.writerow(estimation)
    fd.close()

    if pi1e == 0:
        shutil.rmtree(TRUEDIR)
        shutil.rmtree(PILOTDIR)
        continue

    # predict power
    power.compute_thresholds(EXC)
    power_predicted = power.powercurves(ssrange=range(PILOT,FINAL))

    os.chdir(FILEDIR)
    shutil.rmtree(PILOTDIR)

    ######################
    # analyze final data #
    ######################

    power_true = []
    for s in range(PILOT,(FINAL+1)):
        print("computing true power for "+str(s)+" subjects")
        #analyze data
        FINALDIR = os.path.join(TEMPDIR,"analysis_final_"+str(SEED)+"/")
        os.mkdir(FINALDIR)
        os.chdir(FINALDIR)
        final_cope = []
        for sub in final_subs[range(s)]:
            final_cope.append(os.path.join(HCPDIR,"Disk"+str(disks[sub])+"of5",subjects[sub],"MNINonLinear/Results/tfMRI_"+paradigm,"tfMRI_"+paradigm+"_hp200_s4_level2vol.feat","cope"+contrast+".feat","stats/cope1.nii.gz"))

        cope_cmd = "fslmerge -t cope.nii.gz "+" ".join(final_cope)
        os.popen(cope_cmd).read()

        model.model(len(final_cope),FINALDIR)
        fslcmd = 'flameo --copefile=cope.nii.gz --covsplitfile=design.grp --designfile=design.mat --ld=stats --maskfile=%s --runmode=ols --tcontrastsfile=design.con' %(os.path.join(FILEDIR,'HCP_mask.nii.gz'))
        os.popen(fslcmd).read()

        SPM = nib.load("stats/zstat1.nii.gz").get_data()

        power = poweranalysis.power(spm=SPM,mask=mask,FWHM=[2.5,2.5,2.5],voxsize=1,alpha=0.05,samplesize=PILOT)
        power.extract_peaks(exc=EXC)
        peaks = power.peaktable

        # compute true power for different procedures
        truth = []
        for i in range(len(peaks)):
            peak_act = activation[int(peaks.x[i]),int(peaks.y[i]),int(peaks.z[i])]
            truth.append(peak_act)
        truth = [0 if x == 0 else 1 for x in truth]
        peaks['active'] = truth

        # how many false negatives?
        if not np.sum(peaks['active'])==len(peaks):
            nonactid = np.where(peaks['active']==0)[0]
            pi1 = effectsize.Pi1()
            pi1.pvalues = peaks['pvals'][nonactid]
            pi1.estimate(starts=10)
            NoFalseNegatives = pi1.pi1*len(nonactid)

        power.compute_thresholds(exc=EXC)
        thresholds = power.thres.thresholds

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
            TPR[ind] = float(np.sum(TP))/float(np.sum(true)+NoFalseNegatives)

        power_true.append(TPR)
        os.chdir(TEMPDIR)
        shutil.rmtree(FINALDIR)

        # write away data
        predDF = pd.DataFrame(power_predicted)
        predDF.to_csv(PredictionFile)

        truDF = pd.DataFrame(power_true)
        truDF.to_csv(TrueFile)

        os.chdir(FILEDIR)

    shutil.rmtree(TRUEDIR)

shutil.rmtree(TEMPDIR)
