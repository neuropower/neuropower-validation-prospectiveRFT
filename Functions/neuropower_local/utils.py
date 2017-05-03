import numpy as np
import pandas as pd
import scipy.stats as stats

def cluster(spm,exc,mask):
	""" Extract local maxima from a spm, return a csv file with variables:
	- x coordinate
	- y coordinate
	- z coordinate
	- peak height """

	# make a new array with an extra row/column/plane around the original array
	spm_newdim = tuple(map(lambda x: x+2,spm.shape))
	spm_ext = np.zeros((spm_newdim))
	msk_ext = np.zeros((spm_newdim))
	spm_ext.fill(-100)
	spm_ext[1:(spm.shape[0]+1),1:(spm.shape[1]+1),1:(spm.shape[2]+1)] = spm
	msk_ext[1:(spm.shape[0]+1),1:(spm.shape[1]+1),1:(spm.shape[2]+1)] = mask
	spm_ext = spm_ext * msk_ext
	shape = spm.shape
	spm = None
	# open peak csv
	labels = ['x','y','z','peak']
	peaks = pd.DataFrame(columns=labels)
	# check for each voxel whether it's a peak, if it is, add to table
	for m in range(1,shape[0]+1):
		for n in range(1,shape[1]+1):
			for o in range(1,shape[2]+1):
				surroundings = None
				res = None
				if spm_ext[m,n,o]>exc:
					surroundings=[spm_ext[m-1,n-1,o-1],
					spm_ext[m-1,n-1,o],
					spm_ext[m-1,n-1,o+1],
					spm_ext[m-1,n,o-1],
					spm_ext[m-1,n,o],
					spm_ext[m-1,n,o+1],
					spm_ext[m-1,n+1,o-1],
					spm_ext[m-1,n+1,o],
					spm_ext[m-1,n+1,o+1],
					spm_ext[m,n-1,o-1],
					spm_ext[m,n-1,o],
					spm_ext[m,n-1,o+1],
					spm_ext[m,n,o-1],
					spm_ext[m,n,o+1],
					spm_ext[m,n+1,o-1],
					spm_ext[m,n+1,o],
					spm_ext[m,n+1,o+1],
					spm_ext[m+1,n-1,o-1],
					spm_ext[m+1,n-1,o],
					spm_ext[m+1,n-1,o+1],
					spm_ext[m+1,n,o-1],
					spm_ext[m+1,n,o],
					spm_ext[m+1,n,o+1],
					spm_ext[m+1,n+1,o-1],
					spm_ext[m+1,n+1,o],
					spm_ext[m+1,n+1,o+1]]
					if spm_ext[m,n,o] > np.max(surroundings):
						res =pd.DataFrame(data=[[m-1,n-1,o-1,spm_ext[m,n,o]]],columns=labels)
						peaks=peaks.append(res)
	peaks = peaks.reset_index()

	return peaks

def peakpvalues(peaktable,exc):
	peaks = np.array(peaktable.peak)
	aboveID = np.where(peaks>exc)[0]
	pvals = np.ones(len(peaks))
	pvals[aboveID] = np.exp(-exc*(np.array(peaks[aboveID])-exc))
	peaktable['pvals'] = pvals
	return peaktable

def peakdens3D(x,k):
	# returns the PDF of a peak
	fd1 = 144*stats.norm.pdf(x)/(29*6**(0.5)-36)
	fd211 = k**2.*((1.-k**2.)**3. + 6.*(1.-k**2.)**2. + 12.*(1.-k**2.)+24.)*x**2. / (4.*(3.-k**2.)**2.)
	fd212 = (2.*(1.-k**2.)**3. + 3.*(1.-k**2.)**2.+6.*(1.-k**2.)) / (4.*(3.-k**2.))
	fd213 = 3./2.
	fd21 = (fd211 + fd212 + fd213)
	fd22 = np.exp(-k**2.*x**2./(2.*(3.-k**2.))) / (2.*(3.-k**2.))**(0.5)
	fd23 = stats.norm.cdf(2.*k*x / ((3.-k**2.)*(5.-3.*k**2.))**(0.5))
	fd2 = fd21*fd22*fd23
	fd31 = (k**2.*(2.-k**2.))/4.*x**2. - k**2.*(1.-k**2.)/2. - 1.
	fd32 = np.exp(-k**2.*x**2./(2.*(2.-k**2.))) / (2.*(2.-k**2.))**(0.5)
	fd33 = stats.norm.cdf(k*x / ((2.-k**2.)*(5.-3.*k**2.))**(0.5))
	fd3 = fd31 * fd32 * fd33
	fd41 = (7.-k**2.) + (1-k**2)*(3.*(1.-k**2.)**2. + 12.*(1.-k**2.) + 28.)/(2.*(3.-k**2.))
	fd42 = k*x / (4.*np.pi**(0.5)*(3.-k**2.)*(5.-3.*k**2)**0.5)
	fd43 = np.exp(-3.*k**2.*x**2/(2.*(5-3.*k**2.)))
	fd4 = fd41*fd42 * fd43
	fd51 = np.pi**0.5*k**3./4.*x*(x**2.-3.)
	f521low = np.array([-10.,-10.])
	f521up = np.array([0.,k*x/2.**(0.5)])
	f521mu = np.array([0.,0.])
	f521sigma = np.array([[3./2., -1.],[-1.,(3.-k**2.)/2.]])
	fd521,i = stats.mvn.mvnun(f521low,f521up,f521mu,f521sigma)
	f522low = np.array([-10.,-10.])
	f522up = np.array([0.,k*x/2.**(0.5)])
	f522mu = np.array([0.,0.])
	f522sigma = np.array([[3./2., -1./2.],[-1./2.,(2.-k**2.)/2.]])
	fd522,i = stats.mvn.mvnun(f522low,f522up,f522mu,f522sigma)
	fd5 = fd51*(fd521+fd522)
	out = fd1*(fd2+fd3+fd4+fd5)
	return out
