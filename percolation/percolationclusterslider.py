# -*- coding: utf-8 -*-
"""
Created on Sun Mar 24 13:00:15 2013

@author: svenni
"""

from pylab import *
from scipy.ndimage import measurements

L = 100
r = rand(L,L)
p = 0.4
z = r<p

figure(figsize=(18,6))
subplot(1,3,1)
im1 = imshow(z, origin='lower')
colorbar()
lw, num = measurements.label(z)

# Show image of labeled clusters (shuffled)
subplot(1,3,2)
im2 = imshow(z, origin='lower') # show image clusters as labeled by a shuffled lw
colorbar()

# Calculate areas
subplot(1,3,3)
im3 = imshow(z, origin='lower')
plotxlim=im3.axes.get_xlim()
plotylim=im3.axes.get_ylim()
ontopplot, = plot([0,0,50,50,0], [0,50,50,0,0], color="red")
xlim(plotxlim)
ylim(plotylim)
colorbar()

axp = axes([0.25, 0.02, 0.65, 0.03])
pSlider = Slider(axp, 'p', 0.01, 1.0, valinit=p)
def update(p):
    p = pSlider.val
    z = r<p
    im1.set_data(z)
    im1.set_clim(z.min(), z.max())
    lw, num = measurements.label(z)

    # labeled clusters
    b = arange(lw.max() + 1) # create an array of values from 0 to lw.max() + 1
    shuffle(b) # shuffle this array
    shuffledLw = b[lw] # replace all values with values from b
    im2.set_data(shuffledLw) # show image clusters as labeled by a shuffled lw    
    im2.set_clim(shuffledLw.min(), shuffledLw.max())
    
    # calculate area
    area = (measurements.sum(z, lw, index=range(lw.max() + 1))).astype(int)
    areaImg = area[lw]
    im3.set_data(areaImg)
    im3.set_clim(areaImg.min(), areaImg.max())
    sliced = measurements.find_objects(areaImg == areaImg.max())
    if(len(sliced) > 0):
        sliceX = sliced[0][1]
        sliceY = sliced[0][0]
        ontopplot.set_xdata([sliceX.start, sliceX.start, sliceX.stop, sliceX.stop, sliceX.start])
        ontopplot.set_ydata([sliceY.start, sliceY.stop, sliceY.stop, sliceY.start, sliceY.start])
    else:        
        ontopplot.set_xdata([0])
        ontopplot.set_ydata([0])
    
    
    draw()
pSlider.on_changed(update)

update(0)

show()