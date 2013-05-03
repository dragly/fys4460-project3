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
Lvals = [25,50,75,100,150,200,250,300,400,600]
C = zeros(len(Lvals))
Msc = zeros(len(Lvals))
Mbb = zeros(len(Lvals))
Mde = zeros(len(Lvals))
M = zeros(len(Lvals))
nSamplesMax = 1000
idString = "nsamples" + str(nSamplesMax)
for i in range(len(Lvals)):
    L = Lvals[i]
    lx = L
    ly = L
    p = pc
    if L >= 300:
        nSamples = 50
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
            Pvec, c_eff = FIND_COND (g, lx, ly)
            C[i] += c_eff
            
            x = coltomat ( Pvec , lx , ly )
            P = x * zzz 
            g1 = g[:,0]
            g2 = g[: ,1]
            z1 = coltomat( g1 , lx , ly )
            z2 = coltomat( g2 , lx , ly )
            
            # Calculate flux from top to down (remember that flux is the negative of the pressure difference)
            f2 = zeros ( (lx , ly ))
            for iy in range(ly -1):
                f2[: , iy ] = ( P [: , iy ] - P [: , iy +1]) * z2 [: , iy ]
            
            # Calculate flux from left to right (remember that flux is the negative of the pressure difference)
            f1 = zeros ( (lx , ly ))
            for ix in range(lx-1):
                f1[ ix ,:] = ( P [ ix ,:] - P [ ix +1 ,:]) * z1 [ ix ,:]
            #    
            #    # Find the sum of absolute fluxes in and out of each site
            fn = zeros (( lx , ly ))
            fn = fn + abs ( f1 )
            fn = fn + abs ( f2 )
            # Add for each column, except the leftmost one, the up-down flux, but offset
            fn [: ,1: ly ] = fn [: ,1: ly ] + abs ( f2 [: ,0: ly -1])
            # For the left-most one, add the inverse pressure multiplied with the spanning cluster bool information
            fn [: ,0] = fn [: ,0] + abs (( P [: ,0] - 1.0)*( zzz [: ,0]))
            # For each row except the topmost one, add the left-right flux, but offset
            fn [1: lx ,:] = fn [1: lx ,:] + abs ( f1 [0: lx -1 ,:])
                        
            
            # if the sum of ingoing and outgoing flow of the system is
            # more than the maximum value of flow any place in the system
            # there is no singly connected bonds
            singlyDiff = fn.max() - (sum(fn[:,0]) + sum(fn[:,-1]))
#            print singlyDiff
            if abs(singlyDiff) > 1e-8:
                singlyConnected = False
                scLimit = inf
                colormap = matplotlib.colors.ListedColormap(["#084081", "#7BCCC4", "yellow"])
            else:
                singlyConnected = True
                scLimit = fn.max()
                colormap = matplotlib.colors.ListedColormap(["#084081", "#7BCCC4", "yellow", "red"])
#            print "Has singly connected bonds:", singlyConnected
            
            zsc = (fn > scLimit - 1e-6)
                
            zbb = (abs(fn) > 1e-6)
            zde = (zzz - zbb)
            
            Msc[i] += sum(zsc)
            Mbb[i] += sum(zbb)
            Mde[i] += sum(zde)
            M[i] += sum(zzz)
    C[i] /= nSamples
    Msc[i] /= nSamples
    Mbb[i] /= nSamples
    Mde[i] /= nSamples
    M[i] /= nSamples
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

figure()        
plot(Lvals, Msc)
xlabel(r"$L$")
ylabel(r"$M_{{SC}}$")
savefig("results/1p/L-vs-MSC-" + idString + ".pdf")

figure()
plot(log10(Lvals), log10(Msc))
xlabel(r"$\log L$")
ylabel(r"$\log M_{SC}$")
savefig("results/1p/loglog-L-vs-Msc-" + idString + ".pdf")

DSC = (log10(Msc)[-1] - log10(Msc)[0]) / (log10(Lvals)[-1] - log10(Lvals)[0])
print "DSC", DSC

figure()        
plot(Lvals, Mbb)
xlabel(r"$L$")
ylabel(r"$M_{{BB}}$")
savefig("results/1p/L-vs-MBB-" + idString + ".pdf")

figure()
plot(log10(Lvals), log10(Mbb))
xlabel(r"$\log L$")
ylabel(r"$\log M_{BB}$")
savefig("results/1p/loglog-L-vs-MBB-" + idString + ".pdf")

DBB = (log10(Mbb)[-1] - log10(Mbb)[0]) / (log10(Lvals)[-1] - log10(Lvals)[0])
print "DBB", DBB
figure()        
plot(Lvals, Mde)
xlabel(r"$L$")
ylabel(r"$M_{{DE}}$")
savefig("results/1p/L-vs-MDE-" + idString + ".pdf")

figure()
plot(log10(Lvals), log10(Mde))
xlabel(r"$\log L$")
ylabel(r"$\log M_{DE}$")
savefig("results/1p/loglog-L-vs-MDE-" + idString + ".pdf")

DDE = (log10(Mde)[-1] - log10(Mde)[0]) / (log10(Lvals)[-1] - log10(Lvals)[0])
print "DDE", DDE

figure()        
plot(Lvals, M)
xlabel(r"$L$")
ylabel(r"$M$")
savefig("results/1p/L-vs-M-" + idString + ".pdf")

figure()
plot(log10(Lvals), log10(M))
xlabel(r"$\log L$")
ylabel(r"$\log M$")
savefig("results/1p/loglog-L-vs-M-" + idString + ".pdf")

D = (log10(M)[-1] - log10(M)[0]) / (log10(Lvals)[-1] - log10(Lvals)[0])
print "D", D

etaR = -1/((log10(C)[-1] - log10(C)[0]) / (log10(Lvals)[-1] - log10(Lvals)[0]))
print "etaR", etaR

savetxt("results/1p/data.txt", [D, DSC, DBB, DDE, etaR])

show()