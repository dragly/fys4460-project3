# -*- coding: utf-8 -*-
"""
Created on Wed May  1 14:48:22 2013

@author: svenni
"""

from pylab import *
from scipy.ndimage import measurements
from percolation import FIND_COND, MK_EQSYSTEM, sitetobond, coltomat
from matplotlib.colors import ListedColormap

morgenstemning = loadtxt("morgenstemning.txt")
colormap = ListedColormap(morgenstemning[::-1] / 255.)
colormap = "GnBu_r"

# First , find the backbone
# Generate spanning cluster (l - r spanning )
pc = 0.59275
Lvals = [25,50,75,100,150,200,250,300,400,600]
C = zeros(len(Lvals))
nSamplesMax = 500
idString = "nsamples" + str(nSamplesMax)
for i in range(len(Lvals)):
    L = Lvals[i]
    lx = L
    ly = L
    p = pc
    if L >= 300:
        nSamples = 25
    else:
        nSamples = nSamplesMax
    ncount = 0
    
    print "L=",L
    for j in range(nSamples):
        ncount = 0
        perc = []
    
        while (len(perc)==0):
            ncount = ncount + 1
            if (ncount > 1000):
                print "Couldn't make percolation cluster..."
                break
            
            z=rand(lx,ly)<p
        #        z = array([[1,1,1,1,1],[1,1,1,0,1],[0,1,1,0,1],[0,1,1,1,0],[1,0,1,0,1]])
            lw,num = measurements.label(z)
            perc_x = intersect1d(lw[0,:],lw[-1,:])
            perc = perc_x[where(perc_x > 0)]
#            print "Percolation attempt", ncount
        
        if len(perc) > 0:
            zz = asarray((lw == perc[0]))
            # zz now contains the spanning cluster
            # Transpose
            zzz = zz.T
            #    # Generate bond lattice from this
            g = sitetobond ( zzz )
            #    figure()
            #    imshow(g[:,0].reshape(lx,ly), interpolation='nearest')
            #    figure()
            #    imshow(g[:,1].reshape(lx,ly), interpolation='nearest')
            #    figure()
            #    imshow(zzz, interpolation='nearest')
            #    # Generate conductivity matrix
            P, c_eff = FIND_COND (g, lx, ly)
            C[i] += c_eff
    C[i] /= nSamples
    print "C=",C[i]
figure()        
plot(Lvals, C)
xlabel(r"$L$")
ylabel(r"$C$")
savefig("results/1p/L-vs-C-" + idString + ".pdf")

figure()
plot(log10(Lvals), log10(C))
xlabel(r"$\log L$")
ylabel(r"$\log C$")
savefig("results/1p/loglog-L-vs-C-" + idString + ".pdf")

etaR = (log10(C)[-1] - log10(C)[0]) / (log10(Lvals)[-1] - log10(Lvals)[0])
print -1/etaR

show()