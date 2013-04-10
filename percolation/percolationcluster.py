# -*- coding: utf-8 -*-

from pylab import *
from scipy.ndimage.morphology import binary_closing
from scipy.ndimage import measurements

#from numpy.random import seed
seed(101)

L = 100
Lx = L
Ly = L
r = rand(Lx,Ly)
p = 0.6
z = r<p

figure(figsize=(16,5))
subplot(1,3,1)
imshow(z, origin='lower', interpolation='nearest')
colorbar()
title("Matrix")

# Show image of labeled clusters (shuffled)
lw, num = measurements.label(z)
subplot(1,3,2)
b = arange(lw.max() + 1) # create an array of values from 0 to lw.max() + 1
shuffle(b) # shuffle this array
shuffledLw = b[lw] # replace all values with values from b
imshow(shuffledLw, origin='lower', interpolation='nearest') # show image clusters as labeled by a shuffled lw
colorbar()
title("Labeled clusters")

# Calculate areas
subplot(1,3,3)
labelList = arange(lw.max() + 1)
area = measurements.sum(z, lw, index=labelList)
areaImg = area[lw]
im3 = imshow(areaImg, origin='lower', interpolation='nearest')
colorbar()
title("Clusters by area")

# Bounding boxes
maxLabels = labelList[where(area == area.max())]
print "Found " + str(len(maxLabels)) + " clusters of size " + str(area.max())
for label in maxLabels:
    sliced = measurements.find_objects(lw == label)
    if(len(sliced) > 0):
        sliceX = sliced[0][1]
        sliceY = sliced[0][0]
        plotxlim=im3.axes.get_xlim()
        plotylim=im3.axes.get_ylim()
        plot([sliceX.start, sliceX.start, sliceX.stop, sliceX.stop, sliceX.start], \
                          [sliceY.start, sliceY.stop, sliceY.stop, sliceY.start, sliceY.start], \
                          color="red")
        xlim(plotxlim)
        ylim(plotylim)
        
        if sliceX.stop - sliceX.start >= Lx or sliceY.stop - sliceY.start >= Ly:
            print "Is spanning cluster"
            P = Lx * Ly - area.max()
        else:
            print "Is not spanning cluster"
            P = 0

show()