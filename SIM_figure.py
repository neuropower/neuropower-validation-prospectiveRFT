import matplotlib.pyplot as plt
import palettable.colorbrewer.qualitative as cb
import scipy.stats as stats

twocol = cb.Paired_12.mpl_colors
xn = np.arange(-10,10,0.01)
plt.figure(figsize=(7,5))
plt.hist(peaks.pval.tolist(),lw=0,facecolor=twocol[0],normed=True,bins=np.arange(0,1,0.1),label="observed distribution")
plt.hlines(1-bum["pi1"],0,1,color=twocol[1],lw=3,label="null part of distribution")
plt.plot(xn,stats.beta.pdf(xn,bum["a"],1)+1-bum["pi1"],color=twocol[3],lw=3,label="alternative part of distribution")
plt.xlim([0,1])
plt.ylim([0,8])
plt.title("histogram")
plt.xlabel("peak height")
plt.ylabel("density")
plt.legend(loc="upper right",frameon=False)
plt.show()

plt.figure(figsize=(7,5))
plt.hist(peaks.peak.tolist(),lw=0,facecolor=twocol[0],normed=True,bins=np.arange(exc,8,0.3),label="observed distribution")
plt.xlim([exc,8])
plt.ylim([0,1])
plt.plot(xn,neuropower.nulPDF(xn,exc)*(1-bum["pi1"]),color=twocol[3],lw=3,label="null distribution")
plt.plot(xn,neuropower.altPDF(xn,modelfit["mu"],modelfit["sigma"],exc)*(bum["pi1"]),color=twocol[5],lw=3, label="alternative distribution")
plt.plot(xn,neuropower.mixprobdens(xn,bum["pi1"],modelfit["mu"],modelfit["sigma"],exc=exc,method="RFT"))
plt.title("histogram")
plt.xlabel("peak height")
plt.ylabel("density")
plt.legend(loc="upper right",frameon=False)
plt.show()
