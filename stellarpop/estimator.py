from nested_sampler import NestedSampler

class Estimator(NestedSampler):
    """
    A general estimator class to handle MCMC chains and nested sampling.
      Inheriting classes will need to re-define logLikelihood().
    """

    def __init__(self,priors):
        """
        priors is a dictionary with each key being the name of a
            prior/parameter and each value being a dictionary with two keys:
            'gen' and 'prior'. The 'prior' value is a PyMC stochastic object
            representing a prior, and the 'gen' value is a generating function
            from the prior (frequently the 'prior' object's rand() function).
        """
        self.priors = priors
        self.names = [name for name in self.priors.keys()]

        self._priors = {}
        for key in self.names:
            self._priors[key] = self.priors[key]['prior']

        self.constraint = None
        self._nevolve_steps = 20
        self._required_steps = 5
        self._proposal_cov = None


    def nested_sample(self,niter):
        """
        Wrapper to ensure that the observations have been set.
        """
        try:
            n = self.nobs
        except:
            print "No observations created; initializing 100 observations."
            self.obs(100)
        self._nested_sample(niter)


    def set_constraint(self):
        """
        Creates the constraint object needed for the nested sampling
          loops.
        """
        import pymc,numpy
        def do_likelihood(pars=self._priors):
            p = {}
            for key in self.names:
                p[key] = pars[key]
            if self.logLikelihood(p)>self.logLstar:
                return 0
            else:
                return -numpy.inf

        parents = {'pars':self._priors}
        constraint = pymc.Potential(logp=do_likelihood,name='constraint',
                        parents=parents,doc='Likelihood constraint',
                        verbose=0,cache_depth=2)

        return constraint


    def set_obs(self,nobs):
        """
        Set the number of active points to use, and create each point. Assumes
          no covariance between the priors.
        """
        from numpy import asarray
        NestedSampler.__init__(self)
        self.constraint = self.set_constraint()

        logL = []
        for i in range(nobs):
            pars = {}
            for par in self.names:
                pars[par] = self.priors[par]['gen']()
            self.obs.append(pars)
            logL.append(self.logLikelihood(pars))
        self.logL = asarray(logL)
        self.nobs = nobs


    def logLikelihood(self):
        """ Placeholder for subclass likelihood function """
        pass


    def mcmc(self,outname,niter,nburn=0,nthin=1):
        """
        A method for generating MCMC samples.

            outname - file to write the samples (and sampler state) to
            niter   - number of iterations (*including* burn iterations)
            nburn   - number of burn iterations
            nthin   - thinning interval of the chains
        """
        import pymc,numpy
        import MyStepMethods as MySM

        mcmcpars = {}
        covarpars = []
        for key in self.names:
            mcmcpars[key] = self.priors[key]['prior']
            if 'vals' not in self.priors[key].keys():
                covarpars.append(mcmcpars[key])

        @pymc.observed
        def loglikelihood(value=0.,pars=mcmcpars):
            p = {}
            for key in self.names:
                p[key] = pars[key]
            return self.logLikelihood(p)

        """
        Run the burn chain separately, creating a covariance matrix for the
            science run. We assume that npars*10 are the minimum number of
            iterations required for any meaningful covariance matrix, otherwise
            we use the length of the burn chain.
        """
        ""

        def covfromtrace(sampler,pars):
            n = sampler.trace(pars[0].__name__,0)[:].size
            vals = numpy.empty((len(pars),n))
            for i in range(len(pars)):
                name = pars[i].__name__
                v = sampler.trace(name,0)[:]
                if v.ndim==2 and v.shape[1]==1:
                    v = v[:,0]
                vals[i] = v.copy()
            return numpy.cov(vals)

        sampler = pymc.MCMC([loglikelihood,mcmcpars])
        if nburn<1000:
            nburn = 1000
        delay = nburn/2
        interval = nburn/4+1
        sampler.use_step_method(MySM.MWAdaptiveMetropolis,covarpars,delay=delay,interval=interval,greedy=False,doLikelihood=True)
        sampler.sample(nburn)

        cov = covfromtrace(sampler,covarpars)*4
        sampler = pymc.MCMC([loglikelihood,mcmcpars])
        sampler.use_step_method(MySM.MWAdaptiveMetropolis,covarpars,cov=cov,delay=delay,interval=interval,doLikelihood=True,Markovian=True)
        sampler.sample(nburn)

        cov = covfromtrace(sampler,covarpars)/4
        sampler = pymc.MCMC([loglikelihood,mcmcpars],db='pickle',dbname=outname)
        sampler.use_step_method(MySM.MWAdaptiveMetropolis,covarpars,cov=cov,delay=niter*10,interval=niter*10,doLikelihood=True,Markovian=True)
        sampler.sample(niter+nburn,nburn,nthin)

        sampler.db.commit()


    def evolve_obs(self,discard_index,weight):
        """
        Use MCMC to evolve observations from the nested sampler. This method
            also creates the samples from discarded points. The reason for this
            is that the details of the saved samples might depend on the
            problem at hand (and therefore cannot be blindly handled by the
            nested sampler).
        """
        import pymc
        from numpy.random import randint
        import MyStepMethods as MySM
        discard = {}
        trial = {}
        steps = {}

        """
        Copy the discarded object to the samples list
        """
        obs = self.obs[discard_index]
        for key in self.names:
            discard[key] = obs[key]
        discard['logL'] = self.logL[discard_index]
        discard['logWt'] = discard['logL'] + weight
        self.samples.append(discard)

        """
        Draw a random observation and evolve it using MCMC
        """
        index = randint(len(self.obs))
        # Make sure we don't use the discarded object
        if index==discard_index:
            index = (index+1)%len(self.obs)

        """
        Loop until the active observation has evolved to an `independent' point
            (this will be one iteration until the remaining volume is rather
            small) or through 8 iterations.
        """
        count = 0
        nsteps = self._nevolve_steps

        # Set the sampler values to the observation which will be evolved
        for key in self.names:
            self._priors[key].value = self.obs[index][key]

        sampler = pymc.MCMC([self.constraint,self._priors])

        """
        The MWAdaptiveMetropolis stepping algorithm uses the covariance matrix
            from previous runs for the proposal distribution. It also updates
            the width of the proposal distributions after *every iteration*
            following Skilling's method:

                if accepted: width *= exp(1./n_accepted)
                if rejected: width /= exp(1./n_rejected)

            There are other small speed improvements as well. Note: if the 
            model includes variables that are *not* stochastics, this will
            probably break.

        The delay is set to what *should* be a reasonable value -- after
            _nevolve_steps/3 **successful** iterations, the covariance matrix
            is updated.
        """
        sampler.use_step_method(MySM.MWAdaptiveMetropolis,self._priors.values(),cov=self._proposal_cov,delay=self._nevolve_steps/3)
        while count<8:
            sampler.sample(nsteps)

            """
            Ensure that each parameter has had at least two steps accepted;
                there must be better criteria for insuring that the evolved
                point is independent from the orignal, but this seems to work.
            """
            n = sampler.step_method_dict[self._priors.values()[0]][0].accepted
            if n>=self._required_steps:
                # The point has evolved!
                break

            """
            If the point is not sufficiently evolved, double the number of
                steps and try again. This is a poor approach, since it restarts
                the MCMC instead of continuing from where it left off.
            """
            nsteps = int(nsteps*2)
            count += 1

        """
        Grab the covariance matrix of the trace to inform the proposal
            distribution of the next iteration.
        """
        cov = sampler.step_method_dict[self._priors.values()[0]][0].C
        self.proposal_cov = cov

        """ Take the last sample from the chain to be the evolved point """
        obs = self.obs[discard_index]
        for key in self.names:
            obs[key] = sampler.trace(key)[:][-1]
        self.obs[discard_index] = obs
        self.logL[discard_index] = self.logLikelihood(obs)
