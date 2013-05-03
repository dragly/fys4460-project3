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

morgenstemning = loadtxt("morgenstemning.txt")
colormap = ListedColormap(morgenstemning[::-1] / 255.)
colormap = "GnBu_r"

# First , find the backbone
# Generate spanning cluster (l - r spanning )
L = 256
lx = L
ly = L
p = 0.585
ncount = 0
perc = []
#colormap = 'autumn'

idstring = "L" + str(L) + "-p" + ("%.2f" % p).replace(".", "_")

while (len(perc)==0):
    ncount = ncount + 1
    if (ncount > 100):
        print "Couldn't make percolation cluster..."
        break
    
    z=rand(lx,ly)<p
#        z = array([[1,1,1,1,1],[1,1,1,0,1],[0,1,1,0,1],[0,1,1,1,0],[1,0,1,0,1]])
    lw,num = measurements.label(z)
    perc_x = intersect1d(lw[0,:],lw[-1,:])
    perc = perc_x[where(perc_x > 0)]
    print "Percolation attempt", ncount

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
    p, c_eff = FIND_COND (g, lx, ly)
    #    # Transform this onto a nx x ny lattice
    x = coltomat ( p , lx , ly )
    P = x * zzz 
    g1 = g[:,0]
    g2 = g[: ,1]
    z1 = coltomat( g1 , lx , ly )
    z2 = coltomat( g2 , lx , ly )
    #    # Plotting
    figure()
    ax = subplot(111)
    ax.set_adjustable('box-forced')
    imshow(zzz, interpolation='nearest', cmap=colormap)
    title("Spanning cluster")
#    colorbar()
    grid(color="white")
    savefig("results/1o/spanning-cluster-" + idstring + ".pdf", dpi=300)
    
    #    subplot (2 ,2 ,1) , imagesc ( zzz )
    #    title ( " Spanning cluster ")
    #    axis equal
    figure()
    ax2 = subplot(111, sharex=ax, sharey=ax)
    ax2.set_adjustable('box-forced')
    imshow(P, interpolation='nearest', cmap=colormap, aspect=1)
    title("Pressure")
    colorbar()
    grid(color="white")
    savefig("results/1o/pressure-" + idstring + ".pdf", dpi=300)
    #    subplot (2 ,2 ,2) , imagesc ( P )
    #    title ( " Pressure " )
    #    axis equal
    
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
    figure()
    ax3 = subplot(111, sharex=ax, sharey=ax)
    ax3.set_adjustable('box-forced')
    imshow(fn, interpolation='nearest', cmap=colormap)
    title ( " Flux " )
    colorbar()
    grid(color="white")
    savefig("results/1o/flux-" + idstring + ".pdf", dpi=300)
    
    
    # if the sum of ingoing and outgoing flow of the system is
    # more than the maximum value of flow any place in the system
    # there is no singly connected bonds
    singlyDiff = fn.max() - (sum(fn[:,0]) + sum(fn[:,-1]))
    print singlyDiff
    if abs(singlyDiff) > 1e-8:
        singlyConnected = False
        scLimit = inf
        colormap = matplotlib.colors.ListedColormap(["#084081", "#7BCCC4", "yellow"])
    else:
        singlyConnected = True
        scLimit = fn.max()
        colormap = matplotlib.colors.ListedColormap(["#084081", "#7BCCC4", "yellow", "red"])
    print "Has singly connected bonds:", singlyConnected
    
    #print "fn"
    #print fn
    #    subplot (2 ,2 ,3) , imagesc ( fn )
    zsc = (fn > scLimit - 1e-6)
        
    zbb = (abs(fn) > 1e-6)
    zde = (zzz - zbb)
    #    zbb = ( zzz + 2* zfn )
    #    zbb = zbb / zbb.max()
    figure()
    ax4 = subplot(111, sharex=ax, sharey=ax)
    ax4.set_adjustable('box-forced')
    imshow(zsc * 3 + (zbb - zsc) * 2 + (zzz - zbb), interpolation='nearest', cmap=colormap)
    title ( "Singly connected, backbone and dead ends")
    savefig("results/1o/sc-bb-de-" + idstring + ".pdf", dpi=300)
#    colorbar()
    grid(color="white")
#    subplot(325, sharex=ax, sharey=ax)
#    imshow(zbb*2 + zzz, interpolation='nearest', cmap=colormap)
#    title ( "Backbone and dead ends")
#    colorbar()
#    grid(color="white")
    
#    figure()

#    imshow(c_eff, interpolation='nearest')
    print sum(c_eff)
    show()