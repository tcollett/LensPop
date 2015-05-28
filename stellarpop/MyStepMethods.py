from pymc import AdaptiveMetropolis,Metropolis
from pymc.Node import ZeroProbability
import numpy as np
from math import exp

class MWAdaptiveMetropolis(AdaptiveMetropolis):
    """
    An overwrite of the step() function of PyMC.AdaptiveMetropolis.
    """
    def __init__(self, stochastic, cov=None, delay=1000, scales=None, interval=200, greedy=True, shrink_if_necessary=False, verbose=0,doLikelihood=False,Markovian=False,ignoreBounds=True):

        self.ordering = {}
        self.covorder = {}
        for i in range(len(stochastic)):
            self.ordering[stochastic[i]] = i
        if cov is not None:
            amcov = np.eye(cov.shape[0])
        else:
            amcov = None
        AdaptiveMetropolis.__init__(self, stochastic, amcov,delay,scales,interval,greedy,shrink_if_necessary,verbose)
        self.interval = interval
        self.doLikelihood = doLikelihood
        self.Markovian = Markovian
        self.ignoreBounds = ignoreBounds
        self.set_cov(cov)

    def set_cov(self,cov=None,scales=None):
        if cov is not None:
            self.C = cov
            return
        ord_sc = []
        for i in range(len(self.stochastics)):
            for s in self.stochastics:
                if self.ordering[s]==i:
                    break
            this_value = abs(np.ravel(s.value))
            if not this_value.any():
                this_value = [1.]
            for elem in this_value:
                ord_sc.append(elem)
        # print len(ord_sc), self.dim
        for i in xrange(len(ord_sc)):
            if ord_sc[i] == 0:
                ord_sc[i] = 1.
        self.C = np.eye(self.dim)*ord_sc/50.
        
    def step(self):
        """
        Draw samples from a set of normal distributions, rescaling the widths
            of the distributions after every iteration (such that the chain
            is highly non-Markovian but is still acceptable for evolving
            points). Retains the covariance matrix updating of
            pymc.AdaptiveMetropolis.

        Another (perhaps small) improvement is avoiding the initial calculation
            of loglike since this is (presumably!) 0, ie the current step is
            valid. This might not be true for all models.... The proposed
            loglikelihood is also determined before the sum of the stochastics
            to short-circuit their calculation in the event the loglikelihood
            is -inf.
        """
        # Probability and likelihood for stochastic's current value:
        logp = sum([stochastic.logp for stochastic in self.stochastics])
        loglike = 0
        if self.doLikelihood:
            loglike = self.loglike

        # Sample a candidate value              
        #arrayjump = np.dot(self.proposal_sd, np.random.normal(size=self.proposal_sd.shape[0]))
        arrayjump = np.random.multivariate_normal(np.zeros(self.C.shape[0]),self.C)
        for stochastic in self.stochastics:
            jump = arrayjump[self.ordering[stochastic]]
            if self.ignoreBounds is True:
                stochastic.value = stochastic.value + jump
                continue
            isbound = False
            if ['lower','upper']==stochastic.parent_names:
                isbound = True
                lb = stochastic.parents['lower']
                ub = stochastic.parents['upper']
            elif ['mu','tau','a','b']==stochastic.parent_names:
                isbound = True
                lb = stochastic.parents['a']
                ub = stochastic.parents['b']
            if isbound is True:
                if stochastic.value+jump<lb:
                    stochastic.value = lb
                elif stochastic.value+jump>ub:
                    stochastic.value = ub
                else:
                    stochastic.value = stochastic.value + jump
            else:
                stochastic.value = stochastic.value + jump

        # Metropolis acceptance/rejection test
        accept = False
        try:
            # Probability and likelihood for stochastic's proposed value:
            loglike_p = self.loglike
            logp_p = sum([stochastic.logp for stochastic in self.stochastics])

            if np.log(np.random.random()) < logp_p + loglike_p - logp - loglike:
                accept = True
                self.accepted += 1
            else:
                self.rejected += 1
        except ZeroProbability:
            self.rejected += 1
            logp_p = None
            loglike_p = None

        # Adjust the proposal width to drive the acceptance ratio to 50%
        if self.Markovian is True:
            if self.accepted>self.rejected:
                self.proposal_sd *= exp(1./self.accepted)
            else:
                self.proposal_sd /= exp(1./self.rejected)

        if self._current_iter == self.delay: 
            self.greedy = False
 
        if not accept:
            self.reject()

        if accept or not self.greedy:
            self.internal_tally()

        if self._current_iter>self.delay and self._current_iter%self.interval==0 and not self.Markovian:
            self.update_cov()
    
        self._current_iter += 1


class MWAMetropolis(AdaptiveMetropolis):
    """
    An overwrite of the step() function of PyMC.AdaptiveMetropolis.
    """
    def __init__(self, stochastic, cov=None):

        AdaptiveMetropolis.__init__(self, stochastic, None,1e99,None,1e99,True,False,0)
        self.cov = cov

    def step(self):
        """
        Draw samples from a set of normal distributions, rescaling the widths
            of the distributions after every iteration (such that the chain
            is highly non-Markovian but is still acceptable for evolving
            points). Retains the covariance matrix updating of
            pymc.AdaptiveMetropolis.

        Another (perhaps small) improvement is avoiding the initial calculation
            of loglike since this is (presumably!) 0, ie the current step is
            valid. This might not be true for all models.... The proposed
            loglikelihood is also determined before the sum of the stochastics
            to short-circuit their calculation in the event the loglikelihood
            is -inf.
        """
        # Probability and likelihood for stochastic's current value:
        logp = sum([stochastic.logp for stochastic in self.stochastics])
        loglike = self.loglike
            
        # Sample a candidate value       
#        arrayjump = np.dot(self.proposal_sd, np.random.normal(size=self.proposal_sd.shape[0]))
        arrayjump = np.random.multivariate_normal(np.zeros(self.cov.shape[0]),self.cov)
        for stochastic in self.stochastics:
            jump = arrayjump[self._slices[stochastic]]
            stochastic.value = stochastic.value + jump

        # Metropolis acceptance/rejection test
        accept = False
        try:
            # Probability and likelihood for stochastic's proposed value:
            loglike_p = self.loglike
            logp_p = sum([stochastic.logp for stochastic in self.stochastics])
#            print "current like",loglike
#            print "current priors",logp
#            print "prop like",loglike_p
#            print "prop priors",logp_p
            if np.log(np.random.random()) < logp_p + loglike_p - logp - loglike:
                accept = True
                self.accepted += 1
            else:
                self.rejected += 1
        except ZeroProbability:
            self.rejected += 1
            logp_p = None
            loglike_p = None

        # Adjust the proposal width to drive the acceptance ratio to 50%
        if not accept:
            self.reject()

        if accept:
            self.internal_tally()

        self._current_iter += 1

from numpy.random import normal as rnormal
from pymc.utils import symmetrize
class MWAMatrixMetropolis(Metropolis):
    """Metropolis sampler with proposals customised for symmetric positive definite matrices"""
    def __init__(self, stochastic, scale=1., proposal_sd=None, verbose=None, tally=True):
        Metropolis.__init__(self, stochastic, scale=scale, proposal_sd=proposal_sd, proposal_distribution="Normal", verbose=verbose, tally=tally)

    @staticmethod
    def competence(s):
        """
        The competence function for MatrixMetropolis
        """
        # MatrixMetropolis handles the Wishart family, which are valued as
        # _symmetric_ matrices.
        if any([isinstance(s,cls) for cls in [distributions.Wishart,distributions.InverseWishart,distributions.WishartCov]]):
            return 2
        else:
            return 0

    def propose(self):
        """
        Proposals for positive definite matrix using random walk deviations on the Cholesky
        factor of the current value.
        """

        # Locally store size of matrix
        dims = self.stochastic.value.shape

        # Add normal deviate to value and symmetrize
        dev =  rnormal(0, self.adaptive_scale_factor * self.proposal_sd, size=dims)
#        dev =  rnormal(0, self.proposal_sd, size=dims)
        symmetrize(dev)
        # Replace
        self.stochastic.value = dev + self.stochastic.value

class MVNMetropolis(Metropolis):
    def __init__(self, stochastic, cov=None, scale=1., proposal_sd=None, verbose=None, tally=True):
        Metropolis.__init__(self, stochastic, scale=scale, proposal_sd=proposal_sd, proposal_distribution="Normal", verbose=verbose, tally=tally)
        self.Cov = cov

    @staticmethod
    def competence(s):
        """
        The competence function for MatrixMetropolis
        """
        # MatrixMetropolis handles the Wishart family, which are valued as
        # _symmetric_ matrices.
        if any([isinstance(s,cls) for cls in [distributions.MvNormal,distributions.MvNormalCov]]):
            return 2
        else:
            return 0

    def propose(self):
        """
        Proposals for positive definite matrix using random walk deviations on the Cholesky
        factor of the current value.
        """

        size = self.stochastic.value.size
        dev = np.random.multivariate_normal(np.zeros(size),self.Cov)
        # Replace
        self.stochastic.value = dev + self.stochastic.value

class Binary(MWAMetropolis):
    """
    An overwrite of the step() function of PyMC.AdaptiveMetropolis.
    """
    def __init__(self, stochastic):

        MWAMetropolis.__init__(self, stochastic)

    def step(self):
        """
        Draw samples from a set of normal distributions, rescaling the widths
            of the distributions after every iteration (such that the chain
            is highly non-Markovian but is still acceptable for evolving
            points). Retains the covariance matrix updating of
            pymc.AdaptiveMetropolis.

        Another (perhaps small) improvement is avoiding the initial calculation
            of loglike since this is (presumably!) 0, ie the current step is
            valid. This might not be true for all models.... The proposed
            loglikelihood is also determined before the sum of the stochastics
            to short-circuit their calculation in the event the loglikelihood
            is -inf.
        """
        # Probability and likelihood for stochastic's current value:
        logp = sum([stochastic.logp for stochastic in self.stochastics])
        loglike = self.loglike

        # Sample a candidate value
        print stochastic
        low,hi = stochastic.parents['lower'],stochastic.parents['upper']
        jump = stochastic.value
        while jump==stochastic.value:
            jump = np.random.randint(low,hi+1)
        stochastic.value = jump

        # Metropolis acceptance/rejection test
        accept = False
        try:
            # Probability and likelihood for stochastic's proposed value:
            loglike_p = self.loglike
            logp_p = sum([stochastic.logp for stochastic in self.stochastics])
#            print "current like",loglike
#            print "current priors",logp
#            print "prop like",loglike_p
#            print "prop priors",logp_p
            if np.log(np.random.random()) < logp_p + loglike_p - logp - loglike:
                accept = True
                self.accepted += 1
            else:
                self.rejected += 1
        except ZeroProbability:
            self.rejected += 1
            logp_p = None
            loglike_p = None

        # Adjust the proposal width to drive the acceptance ratio to 50%
        if not accept:
            self.reject()

        if accept:
            self.internal_tally()

        self._current_iter += 1


def covfromtrace(sampler,pars):
    n = sampler.trace(pars[0].__name__,0)[:].size
    vals = numpy.empty((len(pars),n))
    for i in range(len(pars)):
        name = pars[i].__name__
        vals[i] = sampler.trace(name,0)[:]
    return numpy.cov(vals)

