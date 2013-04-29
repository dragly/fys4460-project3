# -*- coding: utf-8 -*-
"""
Created on Wed Apr 24 14:48:01 2013

@author: svenni
"""

from pylab import *
from scipy.ndimage.morphology import binary_closing
from scipy.ndimage import measurements
from percolation import spanningClusterMass

pc = 0.59275

kappa = arange(4,12)
Ls = 2**kappa
    
nSamples = 10000
idString = "nsamples" + str(nSamples)
M = zeros(len(Ls))
i = 0
for L in Ls:
    
    #figure(figsize=(16,5))
    p = pc
    print "Calculating for L=" + str(L)

    M[i] = spanningClusterMass(L,p,nSamples=nSamples)
    i += 1

figure()
xlabel(r"$L$")
ylabel(r"$M$")
plot(Ls, M)
savefig("results/1h/M-vs-L-" + idString + ".pdf")

figure()
xlabel(r"$L$")
ylabel(r"$M$")
loglog(Ls, M)
savefig("results/1h/M-vs-L-loglog-" + idString + ".pdf")

figure()
dlogM = diff(log(M)) / diff(log(Ls))
xlabel(r"$\log L$")
ylabel(r"$d\log M / d\log L$")
plot(log(Ls)[:-1], dlogM)
savefig("results/1h/M-vs-L-derivative-" + idString + ".pdf")


show()
#show()
