# -*- coding: utf-8 -*-
"""
Created on Tue Apr  9 11:51:20 2013

@author: svenni
"""

from pylab import *
from scipy.ndimage.morphology import binary_closing
from scipy.ndimage import measurements



L = 512
Lx = L
Ly = L

nSamples = 500
nPVals = 50

pVals = linspace(0.59, 1, nPVals)
PVals = zeros(nPVals)

#figure(figsize=(16,5))
for pIndex in range(len(pVals)):
#    clf()
    p = pVals[pIndex]
    Psum = 0
    Psamples = 0
    print "Calculating for p=" + str(p)
    for sample in range(nSamples):
        r = rand(Lx,Ly)
        z = r<p
        
#        subplot(1,3,1)
#        imshow(z, origin='lower', interpolation='nearest')
#        colorbar()
#        title("Matrix")
        
        # Show image of labeled clusters (shuffled)
        lw, num = measurements.label(z)
#        subplot(1,3,2)
#        b = arange(lw.max() + 1) # create an array of values from 0 to lw.max() + 1
#        shuffle(b) # shuffle this array
#        shuffledLw = b[lw] # replace all values with values from b
#        imshow(shuffledLw, origin='lower', interpolation='nearest') # show image clusters as labeled by a shuffled lw
#        colorbar()
#        title("Labeled clusters")
        
        # Calculate areas
#        subplot(1,3,3)
        labelList = arange(lw.max() + 1)
        area = measurements.sum(z, lw, index=labelList)
#        areaImg = area[lw]
#        im3 = imshow(areaImg, origin='lower', interpolation='nearest')
#        colorbar()
#        title("Clusters by area")
        
        # Bounding boxes
        maxArea = area.max()
        maxLabels = labelList[where(area == maxArea)]
#        print "Found " + str(len(maxLabels)) + " clusters of size " + str(area.max())
        if area.max() <= 0:
            continue
        
        for label in maxLabels:
            sliced = measurements.find_objects(lw == label)
            if(len(sliced) > 0):
                sliceX = sliced[0][1]
                sliceY = sliced[0][0]
#                plotxlim=im3.axes.get_xlim()
#                plotylim=im3.axes.get_ylim()
#                plot([sliceX.start, sliceX.start, sliceX.stop, sliceX.stop, sliceX.start], \
#                                  [sliceY.start, sliceY.stop, sliceY.stop, sliceY.start, sliceY.start], \
#                                  color="red")
#                xlim(plotxlim)
#                ylim(plotylim)
                
                if sliceX.stop - sliceX.start >= Lx or sliceY.stop - sliceY.start >= Ly:
#                    print "Is spanning cluster"
                    P = maxArea / (Lx * Ly)
                else:
#                    print "Is not spanning cluster"
                    P = 0
#                print "p = " + str(p) + ", P = " + str(P)
                Psum += P
                Psamples += 1
    if Psamples > 0:
        PVals[pIndex] = Psum / Psamples
    else:
        PVals[pIndex] = 0

idString = "L" + str(L) + "-nsamples" + str(nSamples)        

figure()
plot(pVals, PVals)
title("L = " + str(L))
xlabel(r"$p$")
ylabel(r"$P(L,p)$")
grid()
savefig("results/1a/P-vs-p-" + idString + ".pdf")

figure()
pc = 0.59275
validpVals = pVals[where(pVals > pc)]
validPVals = PVals[where(pVals > pc)]
logpVals = log(validpVals - pc)
logPVals = log(validPVals)
plot(logpVals, logPVals)
title("loglog for L = " + str(L))
xlabel(r"$\log(p - p_c)$")
ylabel(r"$\log(P(L,p))$")
grid()
savefig("results/1a/P-vs-p-" + idString + "-loglog.pdf")

figure()
dPdp = (logPVals[1:] - logPVals[:-1]) / (logpVals[1:] - logpVals[:-1])
plot(logpVals[:-1], dPdp)
title("derivative of loglog for L = " + str(L))
xlabel(r"$\log(p - p_c)$")
ylabel(r"$d\log(P)/d\log(p-p_c)$")
grid()
savefig("results/1a/P-vs-p-" + idString + "-loglogderivative.pdf")

savetxt("results/1a/data-" + idString + ".dat", [pVals, PVals])

#show()
