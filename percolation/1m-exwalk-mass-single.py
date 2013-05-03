# -*- coding: utf-8 -*-
"""
Created on Fri Apr 26 11:19:57 2013

@author: svenni
"""

from pylab import *
from scipy.ndimage import measurements
from walk import walk

colormap = "GnBu_r"

pc = 0.59275

# exwalk.py 
# Generate spanning cluster (l-r spanning)
p = pc - 0.02
lx = ly = L = 256
ncount = 0
perc = []
idString = "L" + str(L)
while (len(perc)==0):
    ncount = ncount + 1
    if (ncount >100):
        print "Couldn't make percolation cluster..."
        break
    
    z=rand(lx,ly)<p
    lw,num = measurements.label(z)
    perc_x = intersect1d(lw[0,:],lw[-1,:])
    perc = perc_x[where(perc_x > 0)]
    print ncount

if len(perc) > 0:
    zz = (lw == perc[0])
    # zz now contains the spanning cluster
    figure()
    imshow(zz, interpolation='nearest', origin='upper', cmap=colormap) # Display spanning cluster
    savefig("results/1m/exwalk-" + idString + ".pdf", dpi=300)
    #show()
    #% Run walk on this cluster
    l,r = walk(zz)    
    figure()
    imshow(l, interpolation='nearest', origin='upper', cmap=colormap)
    savefig("results/1m/leftwalk-" + idString + ".pdf", dpi=300)
    figure()
    imshow(r, interpolation='nearest', origin='upper', cmap=colormap)
    savefig("results/1m/rightwalk-" + idString + ".pdf", dpi=300)
    zzz = (l*r > 0) # Find points where both l and r are non-zero
    figure()
    imshow(zzz * 2 + zz, interpolation='nearest', origin='upper', cmap=colormap)
    savefig("results/1m/bothwalk-" + idString + ".pdf", dpi=300)
    show()