# -*- coding: utf-8 -*-
"""
Created on Tue Apr  9 11:51:20 2013

@author: svenni
"""

from pylab import *
from scipy.ndimage.morphology import binary_closing
from scipy.ndimage import measurements
from percolation import clusterNumberDensity

pc = 0.59275

kappa = arange(5,12)
Ls = 2**kappa
    
nSamples = 1000
idString = "nsamples" + str(nSamples)

inclines = zeros(len(Ls))

i = 0
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
    
    savefig("results/1f/dlogn-vs-s-L" + str(L) + "-" + idString + ".pdf")
    
    lateBin = int(nBins * 0.75)
    inclines[i] = (log(totalBins[lateBin])- log(totalBins[0])) / (log(bins[lateBin]) - log(bins[0]))
    i += 1
figure(1)
xlabel(r"$s$")
ylabel(r"$n(s,p)$")
legend()
xlim(xmin=1)
legend(prop={'size':14})

savefig("results/1f/n-vs-s-" + idString + ".pdf")

figure(3)
plot(Ls[1:], -inclines[1:])
xlabel(r"$L$")
ylabel(r"$\tau$")

savefig("results/1f/inclines-" + idString + ".pdf")


plot(log(Ls[1:]), log(-inclines[1:]))

#from scipy.optimize import curve_fit
#def func(L, A, B, alpha):
#    if alpha > 0.05:
#        return 10000000
#    elif alpha < 1e-6:
#        return 10000000
#    return A - B * exp(-alpha * L)
#
#popt,pcov = curve_fit(func, Ls[1:], -inclines[1:], p0=[2,1e-6,1e-6], maxfev=1000)
#plot(Ls[1:], func(Ls[1:],*popt))

show()
#show()
