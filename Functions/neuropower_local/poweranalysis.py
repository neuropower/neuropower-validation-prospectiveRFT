from neuropowermodels import *
from effectsize import *
from utils import *
import numpy as np
import pandas as pd
from scipy import integrate

class power(object):
    '''
    Power calculator for fMRI imaging studies based on peaks.

    :param spm: parametric map
    :type spm: ndarray
    :param mask: mask of parametric map
    :type mask: ndarray
    '''

    def __init__(self,spm,mask=None,samplesize=None,FWHM=None,voxsize=None,alpha=None):

        self.spm = spm
        self.mask = mask
        self.samplesize = samplesize
        self.FWHM = FWHM
        self.voxsize = voxsize
        self.alpha = alpha
        self.cte = integrate.quad(lambda x:x*peakdens3D(x,1),-np.inf,np.inf)[0]

    def extract_peaks(self,exc):

        peaktable = cluster(spm = self.spm, exc = None, mask = self.mask)
        self.peaktable = peakpvalues(peaktable, exc = exc)

        return self

    def estimate_model(self,exc,starts=10):

        if not hasattr(self, 'peaktable'):
            self.extract_peaks(exc=exc)

        peakexc = self.peaktable[self.peaktable['peak']>exc]

        pi1 = Pi1(pvalues = peakexc['pvals'])
        pi1.estimate(starts=starts)
        self.pi1 = pi1.pi1
        self.lam = pi1.lam

        effect = Effect(peakexc['peak'])
        effect.pi1 = pi1.pi1
        effect.exc = exc
        effect.estimate()
        self.mu = effect.mu
        self.sigma = effect.sigma

        return self

    def compute_thresholds(self,exc):

        thres = thresholds(self.peaktable['peak'],self.peaktable['pvals'],FWHM=self.FWHM,voxsize=self.voxsize,alpha=self.alpha,nvox = np.sum(self.mask>0),exc=exc)
        thres.estimate()
        self.thres = thres

        return self

    def predict(self,newss,thresholds=None,exc=None,FDRpredict=False):

        cohen = (self.mu-self.cte)/np.sqrt(self.samplesize)
        mu_projected = cohen*np.sqrt(newss)+self.cte
        #mu_projected = self.mu/np.sqrt(self.samplesize)*np.sqrt(newss)

        if thresholds == None:
            if self.thres.thresholds:
                pass
            elif exc == None:
                raise ValueError('Excursion threshold needed for computing thesholds.  Set threshold or compute thresholds separately.')
            else:
                self.compute_thresholds(exc=EXC)
            if FDRpredict==True:
                self.thres.predict_FDR(mu_projected,self.sigma)
            thresholds = self.thres.thresholds
        if type(thresholds==dict):
            predicted = {k:1-altCDF(v,mu_projected,self.sigma,exc=exc) for (k,v) in thresholds.items() if not v == 'nan'}
        else:
            predicted = 1-altCDF(v,mu_projected,self.sigma,exc=exc)

        return predicted

    def powercurves(self,ssrange,thresholds=None,FDRpredict=True,exc=None):
        #FDRpredict needs a thresholds-object
        if thresholds == None:
            pass
        elif type(thresholds)==float:
            thresholds = {"power":thresholds}
        elif type(thresholds)==dict:
            thresholds = thresholds
        else:
            raise ValueError("unknown thresholds type")
        powercurves = pd.DataFrame({})

        for ss in ssrange:
            res = self.predict(thresholds=thresholds,newss=ss,FDRpredict=True,exc=exc)
            res['samplesize'] = ss
            powercurves = powercurves.append(res,ignore_index=True)

        return powercurves

class thresholds(object):
    '''
    Computes cutoffs.

    :param peaks: list of peaks
    :type peaks: list
    :param pvalues: list of pvalues (same order as peaks)
    :type pvalues: list
	:param FWHM: full width at half maximum
	:type FWHM: float or list of floats
    :param voxsize: size of the voxels
    :type voxsize: float or list of floats
    :param nvox: number of voxels to be evaluated
    :type nvox: integer
    :param alpha: alpha level
    :type alpha: float
    :param exc: the excursion threshold of the map, if set to None, no
    threshold is assumed
    :type exc: float
    '''

    def __init__(self,peaks,pvalues,FWHM,voxsize,nvox,alpha,exc):
        self.peaks = np.asarray(peaks)
        self.pvalues = np.asarray(pvalues)
        self.FWHM_vox = np.asarray(FWHM)/np.asarray(voxsize)
        self.nvox = nvox
        self.alpha = alpha
        self.exc = exc

    def estimate(self):
        'Computes cutoffs'

        # nuldistribution to get cutoffs based on pvalues
        peakrange = np.arange(self.exc,15,0.001)
    	pN = 1-nulCDF(np.array(peakrange),exc=self.exc)

    	# smoothness
    	resels = self.nvox/np.product(self.FWHM_vox)
    	pN_RFT = resels*np.exp(-peakrange**2/2)*peakrange**2

        # compute cutoffs
    	cutoff_UN = np.min(peakrange[pN<self.alpha])
    	cutoff_BF = np.min(peakrange[pN<self.alpha/len(self.peaks)])
    	cutoff_RFT = np.min(peakrange[pN_RFT<self.alpha])

    	#Benjamini-Hochberg
    	pvals_sortind = np.argsort(self.pvalues)
    	pvals_order = pvals_sortind.argsort()
    	FDRqval = pvals_order/float(len(self.pvalues))*0.05
    	reject = self.pvalues<FDRqval

    	if reject.any():
    		FDRc = np.max(self.pvalues[reject])
    	else:
    		FDRc = 0

    	cutoff_BH = 'nan' if FDRc==0 else min(peakrange[pN<FDRc])

        self.UN = cutoff_UN
        self.BF = cutoff_BF
        self.RFT = cutoff_RFT
        self.BH = cutoff_BH
        self.thresholds = {"UN":self.UN,
                           "BF":self.BF,
                           "RFT":self.RFT,
                           "BH":self.BH
                      }

        return self

    def predict_FDR(self,mu,sigma):
        'Predict FDR cutoff based on distributions'

        x = np.arange(self.exc,15,0.01)
        y1 = 1-altCDF(x,mu=mu,sigma=sigma,exc=None)
        y2 = 1-nulCDF(x,exc=self.exc)
        lfdr = y2/(y1+y2)

        if sum(lfdr<0.05)>0:
            cutoff_FDR = min(x[lfdr<0.05])
        else:
            cutoff_FDR = 'nan'

        self.FDRc_predicted = cutoff_FDR
        self.thresholds['FDRc_predicted']=cutoff_FDR

        return self
