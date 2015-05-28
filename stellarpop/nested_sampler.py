
class NestedSampler:
    """
    This class implements the nested sampling algorithm as described by
        Skilling. The point of having a separate sampling class is to allow
        an implementation of the nested sampling that does not require details
        of the problem being estimated. This is a template class, and the
        inheriting class should define the following:

            a function to initialize the observations in the list obs, and set
                the number of samples, nobs (logL should also be initialized!)
            a function to determine the log-likelihood

        as well as re-defining the evolve_obs function.
    """
    def __init__(self):
        from math import log,exp

        self.nobs = None
        self.niter = 0
        self.logwidth = None
        self.logLstar = None

        self.obs = []
        self.logL = []
        self.samples = []

        self.H = 0.
        self.logZ = -1e300
        self.logZtrace = []
        self.Htrace = []

    def evolve_obs(self,index,logwidth):
        """ This should be re-defined by the inheriting class """
        pass


    def _nested_sample(self,niter):
        """
        An implementation of the nested sampling algorithm, with a data model
            that should allow straightforward extension to pausing the
            integration (to externally evaluate convergence, for example).
        """
        from math import log,exp

        """
        First check that there are valid observations and set logwidth if this
            is the first iteration.
        """
        if self.nobs is None or self.nobs==0:
            raise Exception,"The observations have not been initialized."

        if self.logwidth is None:
            logwidth = log(1.-exp(-1./self.nobs))
        else:
            logwidth = self.logwidth

        """ Short names... """
        H = self.H
        logZ = self.logZ

        def logadd(x,y):
            """ A helper function for log addition """
            from math import log,exp
            if x>y:
                return x+log(1.+exp(y-x))
            else:
                return y+log(1.+exp(x-y))


        """
        This is the nested sampling loop, as proposed by Skilling.
        """
        for iter in range(niter):
            index = self.logL.argmin()
            logL = self.logL[index]

            logWt = logL + logwidth
            logZ_new = logadd(logZ,logWt)

            H = exp(logWt-logZ_new)*logL + exp(logZ-logZ_new)*(H+logZ)
            H -= logZ_new
            logZ = logZ_new

            self.logZtrace.append(logZ)
            self.Htrace.append(H)

            self.logLstar = logL
            self.evolve_obs(index,logwidth)

            logwidth -= 1./self.nobs
            if iter>10.*self.nobs*H:
                print "Early exit condition reached after %d iterations:"%iter
                print "  (iter > 10*nobs*H)"
                break
            if iter>100 and self.logZtrace[-1]-self.logZtrace[-100]<0.01:
                print "Early exit condition reached after %d iterations:"%iter
                print "  (delta_logZ < 0.01 over 100 iterations)"
                break


        self.niter += iter
        self.logwidth = logwidth
        self.H = H
        self.logZ = logZ
