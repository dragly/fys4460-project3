# -*- coding: utf-8 -*-
"""
Created on Tue Apr  9 11:51:20 2013

@author: svenni
"""

from pylab import *
from scipy.ndimage.morphology import binary_closing
from scipy.ndimage import measurements
from percolation import clusterNumberDensity


L = 1024
Lx = L
Ly = L
nSamples = 100
nPVals = 12

pc = 0.59275
    
maxBinArea = Lx*Ly
nBins = 50 #min([100, int(log(maxBinArea))])


bins = array(sorted(set(logspace(0, log10(maxBinArea), nBins).astype(int64).tolist())))
pcbins, pcTotalBins = clusterNumberDensity(L,pc,nSamples=nSamples,bins=bins)

figure()
pVals = linspace(0.45, pc - 0.05, nPVals)
PVals = zeros(nPVals)
nBins = len(bins)

sxi = zeros(nPVals)

#figure(figsize=(16,5))
idString = "L" + str(L) + "-nsamples" + str(nSamples)
idString += "-frombelow"
for pIndex in range(len(pVals)):
    p = pVals[pIndex]
    Psum = 0
    Psamples = 0
    print "Calculating for p=" + str(p)
    bins, totalBins = clusterNumberDensity(L,p,nSamples=nSamples,bins=bins)
    totalBins /= pcTotalBins
    
#    plotbins = histogram(totalBins, bins=bins)[0]
    if p == pc:
        label = "p=p_c"
    else:
        label = "p=" + ("%.2f" % p)
    loglog(bins[:-1], totalBins, label=label)
    savetxt("results/1e/data-" + idString + "-p" + ("%.2f" % p) + ".dat", [bins[:-1], totalBins])
    finalElement = where(totalBins>0.5)[0][-1]
    sxi[pIndex] = (bins[finalElement] + bins[finalElement+1]) / 2 # pick the mean of the two elements around 0.5    
        
xlabel(r"$\alpha$")
ylabel(r"$n(s,p)/n(s,p_c)$")
#    if fromAbove == 0:
legend(prop={'size':14}, loc=3)
xlim(xmin=1)
title(r"Number densities near $p_c$")

savefig("results/1g/n-vs-s-" + idString + ".pdf")

figure()
xlabel(r"$|p-p_c|$")
ylabel(r"$s_\xi$")
plot(abs(pVals - pc), sxi)

savefig("results/1g/sxi-" + idString + ".pdf")

figure()
xlabel(r"$\log |p-p_c|$")
ylabel(r"$\log s_\xi$")
plot(log(abs(pVals - pc)), log(sxi))

savefig("results/1g/sxi-loglog-" + idString + ".pdf")

figure()
dlogsxi = diff(log(sxi)) / diff(log(abs(pVals - pc)))

xlabel(r"$\log |p - p_c|$")
ylabel(r"$d\log s_\xi / d\log|p-p_c|$")
plot(log(abs(pVals - pc))[:-1], dlogsxi)

savefig("results/1g/dsxi-loglog-" + idString + ".pdf")

sigma = -1 / ((log(sxi)[-1] - log(sxi)[0]) / (log(abs(pVals - pc))[:-1][-1] - log(abs(pVals - pc))[:-1][0]))

print "sigma = " + str(sigma)
savetxt("results/1g/sigma-" + idString + ".dat", array([sigma]))

show()
#show()
