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
Ls = array([25,50,100,200,400,800])
    
nSamples = 500
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
        nSamples = 40
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
plot(log10(Ls), log10(pPiDiff),'-o')
xlabel(r"$\log10 L$")
ylabel(r"$\log10 (p_{\Pi=0.8} - p_{\Pi=0.3})$")
savefig("results/1i/p-vs-L-log10log10-" + idString + ".pdf")

figure()
pPiDiffDeriv = diff(log10(pPiDiff)) / diff(log10(Ls))
plot(log10(Ls)[:-1], pPiDiffDeriv,'-o')
xlabel(r"$\log10 L$")
ylabel(r"$d \log10 (p_{\Pi=0.8} - p_{\Pi=0.3}) / d\log10 L$")
savefig("results/1i/p-vs-L-log10log10-diff-" + idString + ".pdf")

nu = -1 / (log10(pPiDiff)[-1] - log10(pPiDiff)[0])/(log10(Ls)[-1] - log10(Ls)[0])
print nu
savetxt("results/1i/nu.txt", [nu])

exactNu = 1.3
figure()
plot(Ls**(-1/exactNu), pPi08, '-o', label=r"$\Pi = 0.8$")
plot(Ls**(-1/exactNu), pPi03, '-o', label=r"$\Pi = 0.3$")

savefig("results/1i/pPi-vs-Lpower-" + idString + ".pdf")

linreg08 = polyfit(Ls**(-1/exactNu), pPi08, deg=1)
pc08 = linreg08[1]
print "pc08=", pc08
linreg03 = polyfit(Ls**(-1/exactNu), pPi03, deg=1)
pc03 = linreg03[1]
print "pc08=", pc03

newX = linspace(0,0.1,10)

plot(newX, poly1d(linreg03)(newX), label="polyfit")
plot(newX, poly1d(linreg08)(newX), label="polyfit")
xlabel(r"$L^{-1/\nu}$")
ylabel(r"$p_{\Pi=x}$")
legend()


savefig("results/1i/pPi-vs-Lpower-polyfit-" + idString + ".pdf")

show()


#show()