# -*- coding: utf-8 -*-
"""
Created on Tue Apr  9 11:51:20 2013

@author: svenni
"""

from pylab import *
from scipy.ndimage.morphology import binary_closing
from scipy.ndimage import measurements



L = 256
Lx = L
Ly = L
nSamples = 10000
nPVals = 6
    
maxBinArea = Lx*Ly
nBins = 50 #min([100, int(log(maxBinArea))])

for fromAbove in range(2):
    figure()
    if fromAbove == 0:
        pVals = linspace(0.55, 0.59, nPVals)
    else:
        pVals = linspace(0.6, 0.65, nPVals)
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
        totalBins = zeros(nBins - 1)
    #    totalBins = zeros(maxBinArea)
        for sample in range(nSamples):
            r = rand(Lx,Ly)
            z = r<p
            lw, num = measurements.label(z)
            labelList = arange(num + 1)
            area = measurements.sum(z, lw, index=labelList) 
            maxArea = area.max()
            maxLabels = labelList[where(area == maxArea)]
            for label in maxLabels:
                sliced = measurements.find_objects(lw == label)
                if(len(sliced) > 0):
                    sliceX = sliced[0][1]
                    sliceY = sliced[0][0]
                    if sliceX.stop - sliceX.start >= Lx or sliceY.stop - sliceY.start >= Ly:
                        area[where(area == maxArea)] = 0 # remove the percolating cluster
            
            currentBins = histogram(area, bins=bins, weights=area)[0]
            currentBins /= maxBinArea
            totalBins += currentBins
    #        totalBins += totalAreas / float(maxBinArea)
        totalBins /= nSamples
    #    plotbins = histogram(totalBins, bins=bins)[0]
        loglog(bins[:-1], totalBins, label="p=" + ("%.2f" % p))
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
