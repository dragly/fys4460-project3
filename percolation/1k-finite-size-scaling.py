# -*- coding: utf-8 -*-
"""
Created on Wed Apr 24 16:20:20 2013

@author: svenni
"""

from pylab import *
from scipy.ndimage.morphology import binary_closing
from scipy.ndimage import measurements
from percolation import spanningClusterMass
from percolation import PercolationSystem

pc = 0.59275

#kappa = arange(4,12)
#Ls = 2**kappa
#Ls = [800]
Ls = [25,50,100,200,400,800]
    
nSamples = 5000
M = zeros(len(Ls))
i = 0
nValues = 100
p = linspace(0.54,0.62,nValues)
Pi = zeros(nValues)

pPi08 = zeros(len(Ls))
pPi03 = zeros(len(Ls))
idString = "nsamples" + str(nSamples) + "-nvalues" + str(nValues)

for k in range(len(Ls)):
    if Ls[k] > 300:
        nSamples = 200
    L = Ls[k]
    print "Calculating for L=" + str(L)
    for i in range(len(p)):
        Pi[i] = 0.0
        for j in range(nSamples):
            system = PercolationSystem(L,p[i])
            Pi[i] += system.isPercolating
        Pi[i] /= nSamples
    plot(p,Pi, label="L=" + str(L))
    intersect08 = where(Pi < 0.8)[0][-1]
    intersect03 = where(Pi < 0.3)[0][-1]
    pPi08[k] = p[intersect08]
    pPi03[k] = p[intersect03]
    
legend()
xlabel(r"$p$")
ylabel(r"$\Pi$")
savefig("results/1i/Pi-vs-p-" + idString + ".pdf")

figure()
plot(Ls, pPi08, label=r"$\Pi = 0.8$")
plot(Ls, pPi03, label=r"$\Pi = 0.3$")
xlabel(r"$L$")
ylabel(r"$p_{\Pi=x}$")
legend()
savefig("results/1i/p-vs-L-" + idString + ".pdf")

figure()
pPiDiff = pPi08 - pPi03
plot(log(Ls), log(pPiDiff))
xlabel(r"$\log L$")
ylabel(r"$\log (p_{\Pi=0.8} - p_{\Pi=0.3})$")
savefig("results/1i/p-vs-L-loglog-" + idString + ".pdf")

figure()
pPiDiffDeriv = diff(log(pPiDiff)) / diff(log(Ls))
plot(log(Ls)[:-1], pPiDiffDeriv)
xlabel(r"$\log L$")
ylabel(r"$d \log (p_{\Pi=0.8} - p_{\Pi=0.3}) / d\log L$")
savefig("results/1i/p-vs-L-loglog-diff-" + idString + ".pdf")


show()


#show()