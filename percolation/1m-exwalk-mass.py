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
L = [25,50,100,200,400,800]
L = [25,50,100]
p = linspace(0.58, 0.62, 20)
nSamples = 100

for i in range(len(L)):
    lx = ly = L[i]
    Msc = zeros(len(p))
    for j in range(len(p)):
        print "L=" + str(L[i]) + ", p=" + str(p[j])
        ncount = 0
        perc = []
        
        while (len(perc)==0):
            ncount = ncount + 1
            if (ncount >100):
                print "Couldn't make percolation cluster..."
                break
            
            z=rand(lx,ly)<p[j]
            lw,num = measurements.label(z)
            perc_x = intersect1d(lw[0,:],lw[-1,:])
            perc = perc_x[where(perc_x > 0)]
            print ncount
        
        if len(perc) > 0:
            labelList = arange(num + 1)
            area = measurements.sum(z, lw, index=labelList)
            areaImg = area[lw]
            maxArea = area.max()
            zz = (areaImg == area.max())
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
            zzz = l*r # Find points where both l and r are non-zero
            Msc[j] = sum(zzz > 0)
    plot(p - pc, Msc)
show()