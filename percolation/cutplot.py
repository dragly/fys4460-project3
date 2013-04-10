# -*- coding: utf-8 -*-
"""
Created on Tue Apr  9 16:55:30 2013

@author: svenni
"""

from pylab import *

k = 0.85337999999997971

beta = linspace(0,1,100)
f = beta*(0.9 - 0.59275)**(beta - 1)

plot(beta, f)
plot(beta, k*ones(100))
show()

beta = 0.441488