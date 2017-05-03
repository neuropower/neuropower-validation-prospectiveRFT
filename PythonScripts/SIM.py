############
# preamble #
############

from __future__ import division,print_function
import os
import sys
import numpy as np
import scipy
import math
import nibabel as nib
import shutil
import csv
import sys
sys.path.append(os.path.join(os.environ.get('HOMEDIR'),"Functions/"))
from neuropower_local import utils,poweranalysis,neuropowermodels,effectsize
from neuropower_local.utils import cluster,peakpvalues
import simul_multisubject_fmri_dataset
import model
import uuid
from scipy import integrate
from palettable.colorbrewer.qualitative import Paired_12,Set1_9
import matplotlib.pyplot as plt
import pandas as pd
import glob

EXC = float(os.environ.get("EXC"))
PILOT = int(os.environ.get("PILOT"))
FINAL = int(os.environ.get("FINAL"))
SEED = int(os.environ.get("SEED"))
ADAPTIVE = os.environ.get("ADAPTIVE")
MODEL = os.environ.get("MODEL")
startloop = 0
endloop = 16

SIMFILEDIR = os.environ.get('SIMFILEDIR')
RESDIR = os.environ.get('OUTDIR')
TMPDIR = os.environ.get('TMPDIR')
PEAKDIR = os.environ.get('PEAKDIR')

TEMPDIR = os.path.join(TMPDIR,str(uuid.uuid4()))

# conditions
effectsizes = np.repeat([0.5,1,1.5,2],4)
widths = [2,4,6,8]*4
es_names = np.repeat(["050","100","150","200"],4)
wd_names = [2,4,6,8]*4

resfile = os.path.join(RESDIR,"estimation_SIM_"+str(SEED)+".csv")

for c in np.arange(startloop,endloop):
#c=10
    os.popen("mkdir "+str(TEMPDIR))

    os.chdir(TEMPDIR)

    #####################
    # simulate all data #
    #####################

    #parameters
    smooth_FWHM = 3
    FWHM = [smooth_FWHM,smooth_FWHM,smooth_FWHM]
    smooth_sigma = smooth_FWHM/(2*math.sqrt(2*math.log(2)))
    mask = nib.load(os.path.join(SIMFILEDIR,"SIM_mask.nii"))
    width = widths[c]
    effectsize = effectsizes[c]

    #size of dataset
    if ADAPTIVE == "adaptive":
        total_sub = FINAL
    elif ADAPTIVE == "predictive":
        total_sub = FINAL+PILOT

    # generate noise + signal
    noise = simul_multisubject_fmri_dataset.surrogate_3d_dataset(
        n_subj=total_sub,
        mask=mask,
        sk=smooth_sigma,
        noise_level=1.0,
        width=5.0,
        out_text_file=None,
        out_image_file=None,
        seed=SEED)
    signal = nib.load(os.path.join(SIMFILEDIR,"SIM_activation"+str(width)+".nii")).get_data()
    signal = np.repeat(signal[:, :, :, np.newaxis], total_sub, axis=3)
    low_values_indices = signal < 0.1
    signal[low_values_indices] = 0
    high_values_indices = signal > 0
    signal[high_values_indices] = effectsize
    noise = noise.transpose((1,2,3,0))
    activation = signal[:,:,:,1]
    img=nib.Nifti1Image(activation,np.eye(4))
    img.to_filename("activation.nii.gz")
    data = signal + noise
    data_pilot = data[:,:,:,0:PILOT]
    img=nib.Nifti1Image(data_pilot,np.eye(4))
    img.to_filename("simulation.nii.gz")

    # perform second level OLS analysis on simulated data
    model.model(PILOT,TEMPDIR)
    fslcmd = 'flameo --copefile=simulation.nii.gz --covsplitfile=design.grp --designfile=design.mat --ld=stats --maskfile=%s --runmode=ols --tcontrastsfile=design.con' %(os.path.join(SIMFILEDIR,'SIM_mask.nii'))
    os.popen(fslcmd).read()

    # extract peaks
    SPM = nib.load("stats/zstat1.nii.gz").get_data()
    MASK = nib.load(os.path.join(SIMFILEDIR,'SIM_mask.nii')).get_data()
    SPM = SPM[::-1,:,:]

    ###############################################################
    # estimate and compute model and estimate power on pilot data #
    ###############################################################

    power = poweranalysis.power(spm=SPM,mask=MASK,exc=EXC,FWHM=smooth_FWHM,voxsize=1,alpha=0.05,samplesize=PILOT)
    power.estimate_model()

    ##################
    #nib.load('lala')
    # FigureCanvas
    # twocol = Paired_12.mpl_colors
    # xn = np.arange(-10,30,0.01)
    # nul = [1-power.pi1]*np.asarray(neuropowermodels.nulPDF(xn,exc=EXC))
    # alt = power.pi1*np.asarray(neuropowermodels.altPDF(xn,mu=power.mu,sigma=power.sigma,exc=EXC))
    # mix = neuropowermodels.mixPDF(xn,pi1=power.pi1,mu=float(power.mu),sigma=power.sigma,exc=EXC)
    # xn_p = np.arange(0,1,0.01)
    # alt_p = float(power.pi1)*scipy.stats.beta.pdf(xn_p, float(power.a), 1)+1-float(power.pi1)
    # null_p = [1-power.pi1]*len(xn_p)
    # fig,axs=plt.subplots(1,2,figsize=(14,5))
    # axs[0].hist(power.peaktable.pvals,lw=0,normed=True,facecolor=twocol[0],bins=np.arange(0,1.1,0.1),label="observed distribution")
    # axs[0].set_ylim([0,3])
    # axs[0].plot(xn_p,null_p,color=twocol[3],lw=2,label="null distribution")
    # axs[0].plot(xn_p,alt_p,color=twocol[5],lw=2,label="alternative distribution")
    # axs[0].legend(loc="upper right",frameon=False)
    # axs[0].set_title("Distribution of "+str(len(power.peaktable))+" peak p-values \n $\pi_1$ = "+str(round(float(power.pi1),2)))
    # axs[0].set_xlabel("Peak p-values")
    # axs[0].set_ylabel("Density")
    # axs[1].hist(power.peaktable.peak,lw=0,facecolor=twocol[0],normed=True,bins=np.arange(min(power.peaktable.peak),30,0.3),label="observed distribution")
    # axs[1].set_xlim([np.min(power.peaktable.peak),np.max(power.peaktable.peak)])
    # axs[1].set_ylim([0,1])
    # axs[1].plot(xn,nul,color=twocol[3],lw=2,label="null distribution")
    # axs[1].plot(xn,alt,color=twocol[5],lw=2, label="alternative distribution")
    # axs[1].plot(xn,mix,color=twocol[1],lw=2,label="total distribution")
    #
    # axs[1].set_xlabel("Peak heights (z-values)")
    # axs[1].set_ylabel("Density")
    # axs[1].legend(loc="upper right",frameon=False)
    ##################

    power.estimate_model()
    peaks = power.peaktable

    est_eff = power.mu
    est_sd = power.sigma
    pi1e = power.pi1
    tau = neuropowermodels.TruncTau(est_eff,est_sd,EXC)
    est_exp_eff = est_eff + tau*est_sd


    # compute true parameters
    truth = []
    for i in range(len(peaks)):
        peak_act = activation[int(peaks.x[i]),int(peaks.y[i]),int(peaks.z[i])]
        truth.append(peak_act)

    truth = [0 if x == 0 else 1 for x in truth]
    peaks['active'] = truth
    true_indices = [index for index,value in enumerate(truth) if value == 1]

    true_effectsize = np.mean(peaks.peak[true_indices])
    true_sd = np.std(peaks.peak[true_indices])
    true_pi1 = np.mean(truth)

    # subs = PILOT
    # actfile = os.path.join(PEAKDIR,'peaks_SIM_active_'+str(c)+'_'+str(subs)+'.csv')
    # nonactfile = os.path.join(PEAKDIR,'peaks_SIM_nonactive_'+str(c)+'_'+str(subs)+'.csv')
    # actpeaks = peaks.peak[peaks.active==0]
    # nonactpeaks = peaks.peak[peaks.active==1]
    # with open(actfile, 'a') as f:
    #     actpeaks.to_csv(f,header=False)
    # with open(nonactfile, 'a') as f:
    #     nonactpeaks.to_csv(f,header=False)
    #

    # write away estimation results with true values
    estimation = [effectsize, str(wd_names[c]), pi1e,true_pi1,est_eff,true_effectsize,est_exp_eff,est_sd,true_sd,'a','c']
    fd = open(resfile,"a")
    wr = csv.writer(fd,quoting=csv.QUOTE_NONE)
    if c == 0:
        estnames = ['es','activation','pi1e','pi1t','ese','est','esexp','sde','sdt','bumpar','cohen']
        wr.writerow(estnames)
    wr.writerow(estimation)
    fd.close()

    if pi1e == 0:
        shutil.rmtree(TEMPDIR)
        #continue

    # predict power

    power_predicted = power.powercurves(ssrange=range(PILOT,FINAL),exc_future=EXC)

    shutil.rmtree(TEMPDIR)

    #######################
    # simulate final data #
    #######################

    power_true = []
    for s in range(PILOT,FINAL):

        #simulate data
        os.mkdir(TEMPDIR)
        os.chdir(TEMPDIR)
        new_subj=s

        if ADAPTIVE == "adaptive":
            data_final = data[:,:,:,0:s]
        elif ADAPTIVE == "predictive":
            data_final = data[:,:,:,PILOT:PILOT+s]

        img=nib.Nifti1Image(data_final,np.eye(4))
        img.to_filename(os.path.join("simulation.nii.gz"))

        # analyse data and extract peaks
        model.model(new_subj,TEMPDIR)
        fslcmd = 'flameo --copefile=simulation.nii.gz --covsplitfile=design.grp --designfile=design.mat --ld=stats --maskfile=%s --runmode=ols --tcontrastsfile=design.con' %(os.path.join(SIMFILEDIR,'SIM_mask.nii'))
        os.popen(fslcmd).read()

        SPM = nib.load("stats/zstat1.nii.gz").get_data()
        SPM = SPM[::-1,:,:]
        # peaks = cluster(spm = SPM, exc = EXC, mask = MASK)
        # peaks = peakpvalues(peaks, exc = EXC)


        power = poweranalysis.power(spm=SPM,mask=MASK,exc=EXC,FWHM=smooth_FWHM,voxsize=1,alpha=0.05,samplesize=s)
        power.extract_peaks()
        power.compute_thresholds()
        thresholds = power.thres.thresholds
        peaks = power.peaktable

        # compute true power for different procedures
        truth = []
        for i in range(len(peaks)):
            peak_act = activation[int(peaks.x[i]),int(peaks.y[i]),int(peaks.z[i])]
            truth.append(peak_act)

        truth = [0 if x == 0 else 1 for x in truth]
        peaks['active'] = truth

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
            if MODEL == "RFT" and np.sum(peaks.peak>EXC) == 0:
                TP = FP = TN = FN = 0
            else:
                peaks['pos'] = peaks.peak>thresholds[ind]
                peaks['true'] = peaks.active==1
                if not (method == 1 or method == 3):
                    peaks = peaks[peaks.peak>EXC]
                TP = np.sum([a and b for a,b in zip(peaks.pos,peaks.true)])
                FP = np.sum([a and not b for a,b in zip(peaks.pos,peaks.true)])
                FN = np.sum([b and not a for a,b in zip(peaks.pos,peaks.true)])
                TN = np.sum([not a and not b for a,b in zip(peaks.pos,peaks.true)])
            res[str(ind+"_TP")] = TP
            res[str(ind+"_FP")] = FP
            res[str(ind+"_TN")] = TN
            res[str(ind+"_FN")] = FN
        power_true.append(res)
        shutil.rmtree(TEMPDIR)

        # subs = s
        # actfile = os.path.join(PEAKDIR,'peaks_SIM_active_'+str(c)+'_'+str(subs)+'.csv')
        # nonactfile = os.path.join(PEAKDIR,'peaks_SIM_nonactive_'+str(c)+'_'+str(subs)+'.csv')
        # actpeaks = peaks.peak[peaks.active==0]
        # nonactpeaks = peaks.peak[peaks.active==1]
        # with open(actfile, 'a') as f:
        #     actpeaks.to_csv(f,header=False)
        # with open(nonactfile, 'a') as f:
        #     nonactpeaks.to_csv(f,header=False)
        #
        # write away data
        predfile = os.path.join(RESDIR,'powpre_SIM_'+str(SEED)+'_w_'+str(wd_names[c])+'_e_'+es_names[c]+'.csv')
        predDF = pd.DataFrame(power_predicted)
        predDF.to_csv(predfile)

        trufile = os.path.join(RESDIR,'powtru_SIM_'+str(SEED)+'_w_'+str(wd_names[c])+'_e_'+es_names[c]+'.csv')
        truDF = pd.DataFrame(power_true)
        truDF.to_csv(trufile)
