# -*- coding: utf-8 -*-
"""
Created on Tue Apr  9 11:51:20 2013

@author: svenni
"""

from pylab import *
from scipy.ndimage.morphology import binary_closing
from scipy.ndimage import measurements

pc = 0.59275

kappa = arange(4,10)
Ls = 2**kappa
    
nSamples = 100000
figure()
idString = "nsamples" + str(nSamples)
for L in Ls:
    Lx = L
    Ly = L

    maxBinArea = Lx*Ly
    nBins = 50 #min([100, int(log(maxBinArea))])
    bins = array(sorted(set(logspace(0, log10(maxBinArea), nBins).astype(int64).tolist())))
    nBins = len(bins)
    
    #figure(figsize=(16,5))
    p = pc
    Psum = 0
    Psamples = 0
    print "Calculating for L=" + str(L)
    totalBins = zeros(nBins - 1)
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
    totalBins /= nSamples
    loglog(bins[:-1], totalBins, label="L=" + ("%d" % L))
    savetxt("data" + idString + "-L" + str(L) + ".dat", [bins[:-1], totalBins])

xlabel("s")
ylabel(r"$n(s,p)$")
legend()
xlim(xmin=1)
legend(prop={'size':14})

savefig("results/1f/n-vs-s-" + idString + ".pdf")
show()
#show()
