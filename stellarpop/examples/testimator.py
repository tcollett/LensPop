from estimator import Estimator

class Testimator(Estimator):

    def __init__(self,priors,data,model):
        self.data = data
        self.model = model

        Estimator.__init__(self,priors)


    def logLikelihood(self,pars):
        chi2 = 0.

        for datum in self.data:
            model = self.model(datum,pars)
            chi2 += (model-datum['y'])**2/datum['err']**2
        return -0.5*chi2


"""
A linear fit to noisy data. The model is

    y = 1.2*x - 2.4

The model is passed to the likelihood function in this example to demonstrate
    one possible method for dealing with complicated models (although it is
    clearly unnecessary here). The upshot is that Testimator is a general class
    for *any* simple data/model chi-square.
"""

import pymc,numpy

""" FAKE DATA """
x = numpy.random.random(15)*10.
y = 1.2*x - 2.4
y += numpy.random.normal(0.,0.5,15)

data = []
for i in range(15):
    data.append({'x':x[i],'y':y[i],'err':0.5})


"""
PRIORS -- The pymc object labels should be the same as the priors dictionary
    keys.
"""
slope_prior = pymc.Uniform('slope',-1,2)
intercept_prior = pymc.Uniform('intercept',-5,5)

priors = {'slope':{},'intercept':{}}
priors['slope']['prior'] = slope_prior
priors['slope']['gen'] = slope_prior.rand
priors['intercept']['prior'] = intercept_prior
priors['intercept']['gen'] = intercept_prior.rand


""" MODEL """
def model(data,pars):
    a,b = pars['slope'],pars['intercept']
    return data['x']*a + b


"""
Create the estimator object and set the number of observations, then iterate
  for 1000 iterations.
"""
est = Testimator(priors,data,model)
est.set_obs(1e2)
est.nested_sample(2000)

logZ = est.logZ


"""
Parameter estimation from the stored samples
"""
from math import exp
slope = 0.
slope2 = 0.
intercept = 0.
intercept2 = 0.
for sample in est.samples:
    w = exp(sample['logWt'] - logZ)
    slope += w*sample['slope']
    slope2 += w*sample['slope']**2
    intercept += w*sample['intercept']
    intercept2 += w*sample['intercept']**2
print "LogZ: %5.2f"%logZ
print "Slope: %4.2f +/- %f" % (slope,(slope2-slope**2)**0.5)
print "Intercept: %4.2f +/- %f" % (intercept,(intercept2-intercept**2)**0.5)

import pylab
pylab.plot(est.logZtrace)
pylab.show()
