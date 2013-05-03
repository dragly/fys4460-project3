# -*- coding: utf-8 -*-
"""
Created on Tue Apr 30 13:01:31 2013

@author: svenni
"""

# -*- coding: utf-8 -*-

from pylab import *
from scipy.ndimage import measurements
from percolation import FIND_COND, MK_EQSYSTEM, sitetobond, coltomat
from matplotlib.colors import ListedColormap

colormap = "GnBu_r"

# First , find the backbone
# Generate spanning cluster (l - r spanning )
pc = 0.59275
Lvals = [25,50,100,200,400,600,1000]
#Lvals = [25,50]
pVals = logspace(log10(0.63), log10(0.85), 6)
nSamplesMax = 600
idString = "nsamples" + str(nSamplesMax)
mu = zeros(len(Lvals))
for i in range(len(Lvals)):
    L = Lvals[i]
    lx = L
    ly = L
    C = zeros(len(pVals))
    for pIndex in range(len(pVals)):
        p = pVals[pIndex]
        if L >= 300:
            nSamples = nSamplesMax / 10
        else:
            nSamples = nSamplesMax
        ncount = 0
        
        print "L=",L, "p=",p
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
                Pvec, c_eff = FIND_COND (g, lx, ly)
                C[pIndex] += c_eff
                
        C[pIndex] /= nSamples
        print "C=", C[pIndex]
#        Msc[i] /= nSamples
#        Mbb[i] /= nSamples
#        Mde[i] /= nSamples
#        M[i] /= nSamples
    print "C=",C
    figure(1)        
    plot(pVals - pc, C, label="L " + str(L))
    xlabel(r"$|p-p_c|$")
    ylabel(r"$C$")
    l = legend()
    l.draggable(True)
    savefig("results/1p/p-vs-C-" + idString + ".pdf")
    
    figure(2)
    fitting = polyfit(log10(pVals - pc + 1e-6), log10(C + 1e-6), deg=1)
    mu[i] = fitting[0]
    plot(log10(pVals - pc), log10(C), 'o', label="L " + str(L) + " mu " + ("%.2f" % mu[i]))
#    plot(poly1d(fitting), label="L " + str(L) + " fit")
    print "mu=",mu[i]
    xlabel(r"$\log |p-p_c|$")
    ylabel(r"$\log C$")
    l = legend()
    l.draggable(True)
    savefig("results/1p/loglog-p-vs-C-" + idString + ".pdf")
figure()
plot(Lvals, mu, '-o')
xlabel(r"$L$")
ylabel(r"$\mu$")
savefig("results/1p/mu-vs-L-" + idString + ".pdf")
show()