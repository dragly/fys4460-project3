# -*- coding: utf-8 -*-
"""
Created on Tue Apr  9 11:51:20 2013

@author: svenni
"""

from pylab import *
from scipy.ndimage.morphology import binary_closing
from scipy.ndimage import measurements
from percolation import clusterNumberDensity


L = 256
Lx = L
Ly = L
nSamples = 100
nPVals = 4

pc = 0.59275
    
maxBinArea = Lx*Ly
nBins = 50 #min([100, int(log(maxBinArea))])

for fromAbove in range(2):
    figure()
    if fromAbove == 0:
        pVals = linspace(0.50, pc, nPVals)
    else:
        pVals = linspace(pc, 0.70, nPVals)
    PVals = zeros(nPVals)
    bins = array(sorted(set(logspace(0, log10(maxBinArea), nBins).astype(int64).tolist())))
    nBins = len(bins)
    
    #figure(figsize=(16,5))
    idString = "L" + str(L) + "-nsamples" + str(nSamples)
    if fromAbove == 0:
        idString += "-frombelow"
    else:
        idString += "-fromabove"
    for pIndex in range(len(pVals)):
        p = pVals[pIndex]
        Psum = 0
        Psamples = 0
        print "Calculating for p=" + str(p)
        bins, totalBins = clusterNumberDensity(L,p,nSamples=nSamples,bins=bins)
        
    #    plotbins = histogram(totalBins, bins=bins)[0]
        if p == pc:
            label = "p=p_c"
        else:
            label = "p=" + ("%.2f" % p)
        loglog(bins[:-1], totalBins, label=label)
        savetxt("results/1e/data-" + idString + "-p" + ("%.2f" % p) + ".dat", [bins[:-1], totalBins])
    xlabel(r"$s$")
    ylabel(r"$n(s,p)$")
#    if fromAbove == 0:
    legend(prop={'size':14}, loc=3)
    xlim(xmin=1)
    title(r"Number densities near $p_c$")
    
    savefig("results/1e/n-vs-s-" + idString + ".pdf")
show()
#show()
