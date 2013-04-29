# -*- coding: utf-8 -*-
"""
Created on Wed Apr 24 16:03:19 2013

@author: svenni
"""

from pylab import *
from scipy.ndimage.morphology import binary_closing
from scipy.ndimage import measurements
from percolation import clusterNumberDensity

pc = 0.59275

L = 2048
    
nSamples = 10000
idString = "nsamples" + str(nSamples)
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

bins, totalBins = clusterNumberDensity(L,p,nSamples=nSamples,bins=bins)

figure(1)
loglog(bins[:-1], totalBins, label="L=" + ("%d" % L))
savetxt("results/1f/data" + idString + "-L" + str(L) + ".dat", [bins[:-1], totalBins])

figure(2)
clf()
xlabel(r"$\log s$")
ylabel(r"$d\log n(s,p)/d\log s$")
dlogn = diff(log(totalBins)) / diff(log(bins[:-1]))
plot(log(bins[1:-1]), dlogn)

savefig("results/1f/single-dlogn-vs-s-L" + str(L) + "-" + idString + ".pdf")
figure(1)
xlabel(r"$s$")
ylabel(r"$n(s,p)$")
legend()
xlim(xmin=1)
legend(prop={'size':14})

savefig("results/1f/single-n-vs-s-" + idString + ".pdf")

show()
#show()
