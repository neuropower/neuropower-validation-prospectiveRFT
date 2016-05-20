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
import sys
import imp
from neuropower import BUM, neuropowermodels, cluster
import simul_multisubject_fmri_dataset
import model
import uuid
from palettable.colorbrewer.qualitative import Paired_12,Set1_9
import matplotlib.pyplot as plt

EXC = float(sys.argv[1])
PILOT = int(sys.argv[2])
FINAL = int(sys.argv[3])
SEED = int(sys.argv[4])
SIMFILEDIR = sys.argv[5]
RESDIR = sys.argv[6]
TMPDIR = sys.argv[7]
ADAPTIVE = sys.argv[8]
MODEL = sys.argv[9]

TEMPDIR = os.path.join(TMPDIR,str(uuid.uuid4()))

resfile = os.path.join(RESDIR,"estimation_sim_"+str(SEED)+".csv")
if os.path.isfile(resfile):
    os.remove(resfile)

for c in range(16):

    os.mkdir(TEMPDIR)
    os.chdir(TEMPDIR)

    #####################
    # simulate all data #
    #####################
    # conditions
    effectsizes = np.repeat([0.5,1,1.5,2],4)
    widths = [2,4,6,8]*4
    es_names = np.repeat(["half","one","onehalf","two"],4)
    wd_names = [2,4,6,8]*4

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
    if MODEL == "RFT":
        peaks = cluster.PeakTable(SPM,EXC,MASK)
    elif MODEL == "CS":
        peaks = cluster.PeakTable(SPM,-100,MASK)

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
            modelfit = neuropowermodels.modelfit(peaks.peak,bum['pi1'],EXC=EXC,starts=20,method="RFT")
            est_eff = modelfit['mu']
            est_sd = modelfit['sigma']
            tau = neuropowermodels.TruncTau(est_eff,est_sd,EXC)
            est_exp_eff = est_eff + tau*est_sd
        elif MODEL == "CS":
            modelfit = neuropowermodels.modelfit(peaks.peak,bum['pi1'],starts=5,method="CS")
            est_sd = 'nan'
            xn = np.arange(-10,30,0.01)
            alt = np.asarray(neuropowermodels.altPDF(xn,mu=modelfit['mu'],method="CS"))
            est_eff = xn[alt==np.max(alt)][0]
            est_exp_eff = 'nan'

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

    # FigureCanvas
    # twocol = Paired_12.mpl_colors
    # xn = np.arange(-10,30,0.01)
    # nul = [1-bum['pi1']]*np.asarray(neuropowermodels.nulPDF(xn,method="CS"))
    # alt = bum['pi1']*np.asarray(neuropowermodels.altPDF(xn,mu=est_eff,method="CS"))
    # mix = neuropowermodels.mixPDF(xn,pi1=bum['pi1'],mu=float(modelfit['mu']),method="CS")
    # xn_p = np.arange(0,1,0.01)
    # alt_p = float(bum['pi1'])*scipy.stats.beta.pdf(xn_p, float(bum['a']), 1)+1-float(bum['pi1'])
    # null_p = [1-bum['pi1']]*len(xn_p)
    # fig,axs=plt.subplots(1,2,figsize=(14,5))
    # axs[0].hist(peaks.pval,lw=0,normed=True,facecolor=twocol[0],bins=np.arange(0,1.1,0.1),label="observed distribution")
    # axs[0].set_ylim([0,3])
    # axs[0].plot(xn_p,null_p,color=twocol[3],lw=2,label="null distribution")
    # axs[0].plot(xn_p,alt_p,color=twocol[5],lw=2,label="alternative distribution")
    # axs[0].legend(loc="upper right",frameon=False)
    # axs[0].set_title("Distribution of "+str(len(peaks))+" peak p-values \n $\pi_1$ = "+str(round(float(bum['pi1']),2)))
    # axs[0].set_xlabel("Peak p-values")
    # axs[0].set_ylabel("Density")
    # axs[1].hist(peaks.peak,lw=0,facecolor=twocol[0],normed=True,bins=np.arange(min(peaks.peak),30,0.3),label="observed distribution")
    # axs[1].set_xlim([np.min(peaks.peak),np.max(peaks.peak)])
    # axs[1].set_ylim([0,1])
    # axs[1].plot(xn,nul,color=twocol[3],lw=2,label="null distribution")
    # axs[1].plot(xn,alt,color=twocol[5],lw=2, label="alternative distribution")
    # axs[1].plot(xn,mix,color=twocol[1],lw=2,label="total distribution")
    #
    # axs[1].set_xlabel("Peak heights (z-values)")
    # axs[1].set_ylabel("Density")
    # axs[1].legend(loc="upper right",frameon=False)

    # write away estimation results with true values
    estimation = [effectsize, str(wd_names[c]), bum['pi1'],true_pi1,est_eff,true_effectsize,est_exp_eff,est_sd,true_sd,bum['a']]
    fd = open(resfile,"a")
    wr = csv.writer(fd,quoting=csv.QUOTE_NONE)
    if c == 0:
        estnames = ['es','activation','pi1e','pi1t','ese','est','esexp','sde','sdt','bumpar']
        wr.writerow(estnames)
    wr.writerow(estimation)
    fd.close()

    if bum['pi1'] == 0:
        shutil.rmtree(TEMPDIR)
        continue

    # predict power
    thresholds = neuropowermodels.threshold(peaks.peak,peaks.pval,FWHM=[2.5,2.5,2.5],voxsize=[1,1,1],nvox=np.product(SPM.size),alpha=0.05,exc=EXC,method="RFT")
    effect_cohen = est_eff/np.sqrt(PILOT)
    power_predicted = []
    for s in range(PILOT,FINAL):
        projected_effect = effect_cohen*np.sqrt(s)
        if MODEL == "RFT":
            powerpred =  {k:1-neuropower.altCDF(v,projected_effect,modelfit['sigma'],exc=EXC,method="RFT") for k,v in thresholds.items() if v!='nan'}
        elif MODEL == "CS":
            powerpred =  {k:1-neuropowermodels.altCDF([v],projected_effect,method="CS")[0] for k,v in thresholds.items() if v!='nan'}
        power_predicted.append(powerpred)

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
        if MODEL == "RFT":
            peaks = cluster.PeakTable(SPM,EXC)
            pvalues = np.exp(-EXC*(np.array(peaks.peak)-EXC))
            pvalues = [max(10**(-6),t) for t in pvalues]
        elif MODEL == "CS":
            peaks = cluster.PeakTable(SPM,-100,MASK)
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
    with open(os.path.join(RESDIR,'powpre_sim_'+str(SEED)+'_w_'+str(wd_names[c])+'_e_'+es_names[c]+'.csv'),'wb') as output_file:
        dict_writer = csv.DictWriter(output_file, keys)
        dict_writer.writeheader()
        dict_writer.writerows(toCSV)

    toCSV = power_true
    keys = toCSV[0].keys()
    with open(os.path.join(RESDIR,'powtru_sim_'+str(SEED)+'_w_'+str(wd_names[c])+'_e_'+es_names[c]+'.csv'),'wb') as output_file:
        dict_writer = csv.DictWriter(output_file, keys)
        dict_writer.writeheader()
        dict_writer.writerows(toCSV)
