# -*- coding: utf-8 -*-
"""
Created on Tue Apr  9 18:24:43 2013

@author: svenni
"""

from pylab import *

z = rand(1e7)**(-3+1)
mybins = histogram(z, bins=logspace(0,log10(z.max()) / 2,100))
#gca().set_xscale("log")

cumulativeDistribution = cumsum(mybins[0]) / float(sum(mybins[0]))
figure()
zbins = mybins[1][:-1]
semilogx(zbins, cumulativeDistribution)
title("Cumulative distribution")
xlabel("z")
ylabel("P(Z > z)")

savefig("results/1c/cumulative-distribution.pdf")

figure()
dPdzbins = diff(cumulativeDistribution) / diff(zbins)
plot(log10(zbins[:-1]), log10(dPdzbins))
title("Distribution")
xlabel(r"$\log z$")
ylabel(r"$\log f_Z(z)$")

savefig("results/1c/distribution.pdf")

figure()
d2Pd2zbins = diff(log10(dPdzbins)) / (diff(log10(zbins[:-1])))
plot(log10(zbins[:-2]), d2Pd2zbins, label="Numerical derivative")
nMovAvg = 20
plot(log10(zbins[nMovAvg-1:-2]), movavg(d2Pd2zbins, n=nMovAvg), label="Trailing avg.")
title("Distribution derivative")
xlabel(r"$\log z$")
ylabel(r"$d \log f_Z(z) / d \log z$")
legend()
alpha = (log10(diff(cumulativeDistribution))[30] - log10(diff(cumulativeDistribution))[10]) / (log10((zbins))[29] - log10((zbins))[9])
print alpha

savefig("results/1c/distribution-derivative.pdf")

show()