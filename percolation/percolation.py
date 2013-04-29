# -*- coding: utf-8 -*-
"""
Created on Tue Apr  9 11:49:58 2013

@author: svenni
"""

from pylab import *
from scipy.ndimage import measurements

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