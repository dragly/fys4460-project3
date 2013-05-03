# -*- coding: utf-8 -*-
"""
Created on Fri Apr 26 11:19:57 2013

@author: svenni
"""

from pylab import *
from scipy.ndimage import measurements
from walk import walk

pc = 0.59275

# exwalk.py 
# Generate spanning cluster (l-r spanning)
L = [25,50,100,200,400,600]
#L = [25,50,75,100,125,150,200]
nSamples = 5000
p = pc
Msc = zeros(len(L))
idString = "nSamples" + str(nSamples)
for i in range(len(L)):
    print "L=" + str(L[i])
    lx = ly = L[i]
    if L[i] > 150:
        nSamples = 10
    for sample in range(nSamples):
        ncount = 0
        perc = []
        
        while (len(perc)==0):
            ncount = ncount + 1
            if (ncount >100):
                print "Couldn't make percolation cluster..."
                break
            
            z=rand(lx,ly)<p
            lw,num = measurements.label(z)
            perc_x = intersect1d(lw[0,:],lw[-1,:])
            perc = perc_x[where(perc_x > 0)]
#            print ncount
        
        if len(perc) > 0:
            zz = (lw == perc[0])
            # zz now contains the spanning cluster
    #            figure()
    #            imshow(zz, interpolation='nearest', origin='upper') # Display spanning cluster
    #            savefig("current.pdf")
            #show()
            #% Run walk on this cluster
            l,r = walk(zz)    
        #    figure()
        #    imshow(l, interpolation='nearest', origin='upper')
        #    figure()
        #    imshow(r, interpolation='nearest', origin='upper')
            zzz = (l*r > 0) # Find points where both l and r are non-zero
            Msc[i] += sum(zzz)
    Msc[i] /= nSamples
figure()
plot(L, Msc)
xlabel(r"$L$")
ylabel(r"$M_{SC}$")
savefig("results/1m/M-vs-L-" + idString + ".pdf")
figure()
plot(log10(L), log10(Msc))
xlabel(r"$\log L$")
ylabel(r"$\log M_{SC}$")
savefig("results/1m/loglog-M-vs-L-" + idString + ".pdf")
DSC = (log10(Msc)[-1] - log10(Msc)[0]) / (log10(L)[-1] - log10(L)[0])
print DSC
show()