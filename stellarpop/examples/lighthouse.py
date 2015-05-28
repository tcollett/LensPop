from estimator import Estimator

class Testimator(Estimator):

    def __init__(self,priors,data,model=None):
        self.data = data
        self.model = model

        Estimator.__init__(self,priors)


    def logLikelihood(self,pars):
        import numpy
        from math import pi
        y = pars['y']
        x = pars['x']
        resids = numpy.log(y/pi/((self.data-x)**2 + y**2))
        return resids.sum()



import pymc,numpy

""" FAKE DATA """
data = numpy.array([ 4.73,  0.45, -1.73,  1.09,  2.19,  0.12,
        1.31,  1.00,  1.32,  1.07,  0.86, -0.49, -2.59,  1.73,  2.11,
        1.61,  4.98,  1.71,  2.23,-57.20,  0.96,  1.25, -1.56,  2.45,
        1.19,  2.17,-10.66,  1.91, -4.16,  1.92,  0.10,  1.98, -2.51,
        5.55, -0.47,  1.91,  0.95, -0.78, -0.84,  1.72, -0.01,  1.48,
        2.70,  1.21,  4.41, -4.79,  1.33,  0.81,  0.20,  1.58,  1.29,
        16.19,  2.75, -2.38, -1.79,  6.50,-18.53,  0.72,  0.94,  3.64,
        1.94, -0.11,  1.57,  0.57])


"""
PRIORS -- The pymc object labels should be the same as the priors dictionary
    keys.
"""
xprior = pymc.Uniform('x',-2.,2.)
yprior = pymc.Uniform('y',0.,2.)

priors = {'x':{},'y':{}}
priors['x']['prior'] = xprior
priors['x']['gen'] = xprior.rand
priors['y']['prior'] = yprior
priors['y']['gen'] = yprior.rand


"""
Create the estimator object and set the number of observations, then iterate
  for 1000 iterations.
"""
est = Testimator(priors,data)
est.set_obs(1e2)
est.nested_sample(1000)

logZ = est.logZ


"""
Parameter estimation from the stored samples
"""
from math import exp
x = 0.
x2 = 0.
y = 0.
y2 = 0.
for sample in est.samples:
    w = exp(sample['logWt'] - logZ)
    x += w*sample['x']
    x2 += w*sample['x']**2
    y += w*sample['y']
    y2 += w*sample['y']**2
print "LogZ: %5.2f"%logZ
print "Information: %5.2f"%est.H
print "X: %4.2f +/- %f" % (x,(x2-x**2)**0.5)
print "Y: %4.2f +/- %f" % (y,(y2-y**2)**0.5)

#import pylab
#pylab.plot(est.logZtrace)
#pylab.show()
