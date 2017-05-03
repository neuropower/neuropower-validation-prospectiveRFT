import scipy.stats
import numpy as np
import scipy.integrate as integrate

def altPDF(peaks, mu, sigma, exc=None):
	"""
	Probability density function of the peaks in a neuroimaging map
	under the alternative hypothesis.

	:param peaks: quantiles
	:type peaks: array_like
	:param mu: location parameter
	:type mu: float
	:param sigma: scale parameter
	:type sigma: float
	:param exc: the excursion threshold of the map, if set to None, no
	threshold is assumed
	:type exc: float

	:returns fa: probability density function evaluated at x
	"""

	peaks = np.asarray(peaks)
	ksi = (peaks - mu) / sigma
	num = 1 / sigma * scipy.stats.norm.pdf(ksi)
	if exc:
		alpha = (exc-mu)/sigma
		den = 1. - scipy.stats.norm.cdf(alpha)
	else:
		den = 1.
	fa = num/den
	return fa

def nulPDF(peaks,exc):
	"""
	Probability density function of the peaks in a neuroimaging map
	under the null hypothesis.

	:param peaks: quantiles
	:type peaks: array_like
	:param exc: the excursion threshold of the map
	:type exc: float

	:returns f0: probability density function evaluated at x
	"""

	peaks = np.asarray(peaks)
	f0 = exc*np.exp(-exc*(peaks-exc))
	return f0

def altCDF(peaks,mu,sigma,exc=None):
	"""
	Cumulative density function of the peaks in a neuroimaging map
	under the alternative hypothesis.

	:param peaks: quantiles
	:type peaks: array_like
	:param mu: location parameter
	:type mu: float
	:param sigma: scale parameter
	:type sigma: float
	:param exc: the excursion threshold of the map, if set to None, no
	threshold is assumed
	:type exc: float

	:returns Fa: cumulative density function evaluated at x
	"""

	peaks = np.asarray(peaks)
	ksi = (peaks-mu)/sigma
	if exc:
		alpha = (exc-mu)/sigma
		Fa = (scipy.stats.norm.cdf(ksi) - scipy.stats.norm.cdf(alpha))/(1-scipy.stats.norm.cdf(alpha))
	else:
		Fa = scipy.stats.norm.cdf(ksi)
	return Fa

def TruncTau(mu,sigma,exc):
	"""
	Help function to compute the expected peak height, given that the mean and expected
	value differ due to truncation.

	:param mu: location parameter
	:type mu: float
	:param sigma: scale parameter
	:type sigma: float
	:param exc: the excursion threshold of the map
	:type exc: float

	:returns tau: scaling factor
	"""
	num = scipy.stats.norm.cdf((exc-mu)/sigma)
	den = 1-scipy.stats.norm.pdf((exc-mu)/sigma)
	tau = num/den
	return tau

def nulCDF(peaks,exc):
	"""
	Cumulative density function of the peaks in a neuroimaging map
	under the alternative hypothesis.

	:param peaks: quantiles
	:type peaks: array_like
	:param exc: the excursion threshold of the map,
	:type exc: float

	:returns Fa: cumulative density function evaluated at x
	"""
	peaks = np.asarray(peaks)
	F0 = 1-np.exp(-exc*(peaks-exc))
	return F0

def mixPDF(peaks,pi1,mu,sigma,exc):
	"""
	Probability density function of the peaks in a neuroimaging map
	from a mix of the null and alternative hypothesis.

	:param peaks: quantiles
	:type peaks: array_like
	:param pi1: proportion of Ha in the total distribution
	:type pi1: float
	:param mu: location parameter
	:type mu: float
	:param sigma: scale parameter
	:type sigma: float
	:param exc: the excursion threshold of the map, if set to None, no
	threshold is assumed
	:type exc: float

	:returns Fa: cumulative density function evaluated at x
	"""
	peaks = np.array(peaks)
	f0=nulPDF(peaks,exc=exc)
	fa=altPDF(peaks,mu,sigma=sigma,exc=exc)
	f=[(1-pi1)*x + pi1*y for x, y in zip(f0, fa)]
	return f

def mixPDF_SLL(pars,peaks,pi1,exc=None):
	"""
	Returns the negative sum of the loglikelihood of mixPDF

	:param pars: [mu,sigma]
	:type pars: list of two objects
	:param peaks: quantiles
	:type peaks: array_like
	:param pi1: proportion of Ha in the total distribution
	:type pi1: float
	:param exc: the excursion threshold of the map
	:type exc: float

	:returns LL: negative sum of the loglikelihood of mixPDF
	"""
	mu = pars[0]
	sigma = pars[1]
	f = mixPDF(peaks,pi1=pi1,mu=mu,sigma=sigma,exc=exc)
	LL = -sum(np.log(f))
	return(LL)
