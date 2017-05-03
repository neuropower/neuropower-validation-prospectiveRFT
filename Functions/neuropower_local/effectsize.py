from neuropowermodels import *
import numpy as np
from scipy.optimize import minimize

class Pi1(object):
    '''
    The proportion of active peaks in a statistical parametric map.
    '''

    def __init__(self,pi1=None,pvalues=None):

        self.pi1 = pi1
        self.pvalues = np.asarray(pvalues)

    def estimate(self, starts=1,seed=None):
        '''
        Searches the maximum likelihood estimator for the shape parameters of the BUM-model given a list of p-values.
        The BUM-model is introduced in pounds & Morris, 2003.

        :param starts: number of starting points
        :type starts: integer
        :param seed: seed for random start point
        :type seed: integer
        '''

        seed = np.random.uniform(0,1000,1) if not 'seed' in locals() else seed
    	a = np.random.uniform(0.05,0.95,(starts,))
    	l = np.random.uniform(0.05,0.95,(starts,))
    	best = []
    	par = []
    	x = [10**(-6) if y<= 10**(-6) else y for y in self.pvalues] #optimiser is stuck when p-values == 0

    	for i in range(0,starts):
    		pars = np.array((a[i],l[i]))
    		opt = minimize(self.fbumnLL,[pars[0],pars[1]],method='L-BFGS-B',args=(x,),jac=self.fpLL,bounds=((0.00001,1),(0.00001,1)))
    		best.append(opt.fun)
    		par.append(opt.x)
    	minind=best.index(np.nanmin(best))
    	bestpar=par[minind]

    	self.pi1 = 1-(bestpar[1] + (1-bestpar[1])*bestpar[0])
        self.a = bestpar[0]
        self.lam = bestpar[1]
        self.maxloglikelihood = best[minind]

        return self

    @staticmethod
    def fpLL(pars, x):
        '''
        Gradient function of the Beta Uniform Model.
        :param pars: [a,l]
        :type pars: list of 2
        :param x: p-value
        :type x: float
        :returns G: [dl,da]
        '''

    	a = pars[0]
    	l = pars[1]
    	dl = -sum((1-a*x**(a-1))/(a*(1-l)*x**(a-1)+l))
    	da = -sum((a*(1-l)*x**(a-1)*np.log(x)+(1-l)*x**(a-1))/(a*(1-l)*x**(a-1)+l))
        G = np.asarray([dl,da])
    	return G

    @staticmethod
    def fbumnLL(pars, x):
        '''
        Returns the negative sum of the loglikelihood of the Beta Uniform Model.
        :param pars: [a,l]
        :type pars: list of 2
        :param x: p-value
        :type x: float
        :returns LL: negative sum of the loglikelihood
        '''

    	a = pars[0]
    	l = pars[1]
    	L = l+(1-l)*a*x**(a-1)
    	negsumlog = -sum(np.log(L))

    	return negsumlog

class Effect(object):
    '''
    Effect size estimator.

    :param peaks: list of peaks
    :type peaks: list
    :param pi1: proportion of Ha
    :type pi1: float
	:param exc: the excursion threshold of the map, if set to None, no
	threshold is assumed
	:type exc: float
    '''

    def __init__(self,peaks,pi1=None,exc=None):
        self.peaks = np.asarray(peaks)
        self.pi1 = pi1
        self.exc = exc

    def estimate(self, starts=1,seed=None):
        '''
        Searches the maximum likelihood estimator for the effect size.

        :param starts: number of starting points
        :type starts: integer
        :param seed: seed for random start point
        :type seed: integer
        '''

    	seed = np.random.uniform(0,1000,1) if not 'seed' in locals() else seed
    	mus = np.random.uniform(self.exc+(1./self.exc),10,(starts,))
    	seed = np.random.uniform(0,1000,1) if not 'seed' in locals() else seed
    	sigmas = np.random.uniform(0.1,10,(starts,))

    	best = []
    	par = []

    	for i in range(0,starts):
    		opt = minimize(mixPDF_SLL,[mus[i],sigmas[i]],method='L-BFGS-B',args=(self.peaks,self.pi1,self.exc),bounds=((self.exc+(1./self.exc),50),(0.1,50)))
    		best.append(opt.fun)
    		par.append(opt.x)
    	minind=best.index(np.nanmin(best))

        self.maxloglikelihood = best[minind]
        self.mu = par[minind][0]
        self.sigma = par[minind][1]

    	return self
