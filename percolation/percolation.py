# -*- coding: utf-8 -*-
"""
Created on Tue Apr  9 11:49:58 2013

@author: svenni
"""

from pylab import *
from scipy.ndimage import measurements
from scipy.sparse import spdiags, dia_matrix, coo_matrix
from scipy.sparse.linalg import spsolve

class PercolationSystem():
    def __init__(self,L,p):
        self.r = rand(L,L)
        self.z = self.r<p
        self.lw, self.num = measurements.label(self.z)
        self.labelList = arange(self.num + 1)
        self.area = measurements.sum(self.z, self.lw, index=self.labelList)
        self.maxArea = self.area.max()
        self.maxLabels = self.labelList[where(self.area == self.maxArea)]
        lw = self.lw
        firstRow = lw[0,:]
        lastRow = lw[-1,:]
        firstColumn = lw[:,0]
        lastColumn = lw[:,-1]
        topBottomSpanningClusters = intersect1d(firstRow, lastRow)
        leftRightSpanningClusters = intersect1d(firstColumn, lastColumn)
        self.spanningClusters = union1d(topBottomSpanningClusters, leftRightSpanningClusters)
        self.spanningClusters = delete(self.spanningClusters, where(self.spanningClusters == 0))
        self.isPercolating = (len(self.spanningClusters) > 0)
        
        
        

def clusterNumberDensity(L,p,nSamples=100,bins=None):
    if bins == None:
        bins = array(sorted(set(logspace(0, log10(maxBinArea), nBins).astype(int64).tolist())))        
        
    totalBins = zeros(len(bins) - 1)
    maxBinArea = L*L
    for sample in range(nSamples):
        system = PercolationSystem(L,p)
        lw = system.lw
        area = system.area
        maxArea = system.maxArea
        maxLabels = system.maxLabels
        for label in maxLabels:
            sliced = measurements.find_objects(lw == label)
            if(len(sliced) > 0):
                sliceX = sliced[0][1]
                sliceY = sliced[0][0]
                if sliceX.stop - sliceX.start >= L or sliceY.stop - sliceY.start >= L:
                    area[where(area == maxArea)] = 0 # remove the percolating cluster
        
#            currentBins = histogram(area, bins=bins, weights=area)[0]
        binAreas = diff(bins)
        currentBins = histogram(area, bins=bins)[0].astype(float)
        currentBins /= binAreas # normalize to values of 1x1 s bins
        currentBins /= maxBinArea
        totalBins += currentBins
#        totalBins += totalAreas / float(maxBinArea)
    totalBins /= nSamples
    return bins, totalBins

def spanningClusterMass(L,p,nSamples=100):
    percolatingMass = 0
    for sample in range(nSamples):
        system = PercolationSystem(L,p)
        lw = system.lw
        for label in system.maxLabels:
            sliced = measurements.find_objects(lw == label)
            if(len(sliced) > 0):
                sliceX = sliced[0][1]
                sliceY = sliced[0][0]
                if sliceX.stop - sliceX.start >= L or sliceY.stop - sliceY.start >= L:
                    percolatingMass += system.maxArea
    percolatingMass /= nSamples
    return percolatingMass
    
    

#
# Written by Marin Soreng
# ( C ) 2004
#
# Calculates the effective flow conductance Ceff of the
# lattice A as well as the pressure P in every site .
def FIND_COND (A , X , Y ):
    P_in = 1.
    P_out = 0.
    # Calls MK_EQSYSTEM .
    B,C = MK_EQSYSTEM (A , X , Y )
    #print "B"
    #print B.todense()
    #print "C"
    #print C
    # Kirchhoff ’ s equations solve for P
    P = spsolve(B, C)
    # The pressure at the external sites is added
    # ( Boundary conditions )
    P = concatenate((P_in * ones (X), P,  P_out * ones (X)))
    # Calculate Ceff
    # second-last X elements of P multiplied with second-last elements of A
    # these are the second last column of the system
    # gives the conductivity of the system per row?
    Ceff = dot((P[-1-2*X:-1-X] - P_out).T, A[-1-2*X:-1-X, 1]) / ( P_in - P_out )
    #print "P" 
    #print P
    #print "Ceff"
    #print Ceff
    return P , Ceff

#
# Written by Marin S r e n g
# ( C ) 2004
#
# Sets up Kirchoff ’ s equations for the 2 D lattice A .
# A has X * Y rows and 2 columns . The rows indicate the site ,
# the first column the bond perpendicular to the flow direction
# and the second column the bond parallel to the flow direction .
#
# First we use that the flow into one site equals the pressure gradient
# between two sites. This is the same as saying that the potential difference
# is
#    U = RI
# and setting R = 1 so that U = I, for each site boundary.
#
# For each site, set up the flow into and out of the site as an equality.
# Multiply each flow with the probability for a flow to be there, i.e. 
# the 2xN matrix A in our case.
# Replace all flows with the pressure differences, I = P1 - P2.
# Factor out all the pressures, so you get P1(flowOutDown + flowOutRight + ...) = 0
# Notice that the right hand side now is zero or, for the left boundary, 1.
# Set this up as a matrix problem Dx = C, where x is the pressures, D is given
# by the flow elements (flowOutDown + flowOutRight + ...) and C is 0 or 1 when it
# is at a boundary.
#
# The equations are set up by noticing that all flow into each site
# must flow out of the site.
# 
#
# The return values are [B , C ] where B * x = C . This is solved
# for the site pressure by x = B \ C .

def MK_EQSYSTEM (A , X , Y ):
    # Total no of internal lattice sites
    sites = X *( Y - 2)
    #print "sites:", sites
    # Allocate space for the nonzero upper diagonals
    main_diag = zeros(sites)
    upper_diag1 = zeros(sites - 1)
    upper_diag2 = zeros(sites - X)
    # Calculates the nonzero upper diagonals
    main_diag = A[X:X*(Y-1), 0] + A[X:X*(Y-1), 1] + A[0:X*(Y-2), 1] + A[X-1:X*(Y-1)-1, 0]
    upper_diag1 = A [X:X*(Y-1)-1, 0]
    upper_diag2 = A [X:X*(Y-2), 1]
    main_diag[where(main_diag == 0)] = 1
    # Constructing B which is symmetric , lower = upper diagonals .
    B = dia_matrix ((sites , sites)) # B *u = t
    B = - spdiags ( upper_diag1 , -1 , sites , sites )
    B = B + - spdiags ( upper_diag2 ,-X , sites , sites )
    B = B + B.T + spdiags ( main_diag , 0 , sites , sites )
    # Constructing C
    C = zeros(sites)
    #    C = dia_matrix ( (sites , 1) )
    C[0:X] = A[0:X, 1]
    C[-1-X+1:-1] = 0*A [-1 -2*X + 1:-1-X, 1]
    return B , C

def sitetobond ( z ):
    #
    # Function to convert the site network z (L , L ) into a ( L *L ,2) bond
    # network
    # g [i,0] gives bond perpendicular to direction of flow
    # g [i,1] gives bond parallel to direction of flow
    # z [ nx , ny ] -> g [ nx * ny , 2]
    #
    nx = size (z ,1 - 1)
    ny = size (z ,2 - 1)
    N = nx * ny 
    # g = zeros (N ,2)
    gg_r = zeros ((nx , ny)) # First , find these
    gg_d = zeros ((nx , ny )) # First , find these
    gg_r [:, 0:ny - 1] = z [:, 0:ny - 1] * z [:, 1:ny]
    gg_r [: , ny  - 1] = z [: , ny  - 1]
    gg_d [0:nx - 1, :] = z [0:nx - 1, :] * z [1:nx, :]
    gg_d [nx - 1, :] = 0
    #print "gg_r"
    #print gg_r
    #print "gg_d"
    #print gg_d
    # Then , concatenate gg onto g
    g = zeros ((nx *ny ,2))
    g [:, 0] = gg_d.reshape(-1,order='F').T
    g [:, 1] = gg_r.reshape(-1,order='F').T
    return g
    
def coltomat (z, x, y):
    # Convert z ( x * y ) into a matrix of z (x , y )
    # Transform this onto a nx x ny lattice
    g = zeros ((x , y))
    #print "For"
    for iy in range(1,y):
        i = (iy - 1) * x + 1
        ii = i + x - 1
        #print iy, i, ii
        g[: , iy - 1] = z[ i - 1 : ii]
    return g