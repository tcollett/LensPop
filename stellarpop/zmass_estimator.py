from stellarpop.estimator import Estimator

class MassEstimator(Estimator):
    """
    An object used to determine estimates of stellar masses. This inherits
        from the base class NestedSampler, although this functionality is not
        necessary for simple MCMC chains.
    """

    def __init__(self,priors,data,model):
        """
        priors is a dictionary of dictionaries describing the prior for each
            parameter. Each parameter's dictionary contains two entries, `gen'
            and `prior'; the former is a function that generates samples from
            the prior and the latter is a PyMC #stochastic# object. Covariances
            are not yet supported. An example is something like this:

            priors = {'mass':{'gen':mass_genfunc,'prior':mass_prior}}

        data are the observational data, dictionaries containing (at least)
            four keys: filter, redshift, mag and error.
        model is an SPSmodel object
        """
        self.data = data
        self.model = model
        self.priors = priors
        self.names = priors.keys()


    def fastMCMC(self,niter,nburn,nthin=1):
        from Sampler import SimpleSample as sample
        from scipy import interpolate
        import pymc,numpy,time
        import ndinterp

        models = self.model.models
        data = self.data
        filters = data.keys()
        t = time.time()
        T1 = models[filters[0]]*0.
        T2 = 0.
        for f in filters:
            T1 += (models[f]-data[f]['mag'])/data[f]['sigma']**2
            T2 += 2.5/self.data[f]['sigma']**2
        M = T1/T2
        logp = 0.
        for f in filters:
            logp += -0.5*(-2.5*M+models[f]-data[f]['mag'])**2/data[f]['sigma']**2
        t = time.time()
        axes = {}
        i = 0
        ax = {}
        ind = numpy.unravel_index(logp.argmax(),logp.shape)
        best = []
        for key in self.model.axes_names:
            a = self.model.axes[key]['points']
            axes[i] = interpolate.splrep(a,numpy.arange(a.size),k=1,s=0)
            ax[key] = i
            best.append(a[ind[i]])
            i += 1

        print logp.max()
        logpmodel = ndinterp.ndInterp(axes,logp,order=1)
        massmodel = ndinterp.ndInterp(axes,M,order=1)

        pars = [self.priors[key] for key in self.names]

        doExp = []
        cube2par = []
        i = 0
        for key in self.names:
            if key.find('log')==0:
                pntkey = key.split('log')[1]
                #self.priors[key].value = numpy.log10(best[ax[pntkey]])
                doExp.append(True)
            else:
                pntkey = key
                doExp.append(False)
                #self.priors[key].value = best[ax[pntkey]]
            cube2par.append(ax[pntkey])
        doExp = numpy.array(doExp)==True
        par2cube = numpy.argsort(cube2par)

        logp -= logp.max()
        p = numpy.exp(logp)
        p /= p.sum()
        i = 0
        wmean = numpy.empty(p.ndim)
        axarr = []
        for key in self.model.axes_names:
            a = self.model.axes[key]['points']
            p0 = numpy.rollaxis(p,i,p.ndim)
            wmean[i] = (a*p0).sum()
            axarr.append(numpy.rollaxis(a+p0*0,p.ndim-1,i))
            i += 1
        cov = numpy.empty((p.ndim,p.ndim))
        #for i in range(p.ndim):
        #    for j in range(i,p.ndim):
        #        cov[i,j] = (p*(axarr[i]-wmean[i])*(axarr[j]-wmean[j])).sum()
        #        cov[j,i] = cov[i,j]
        for i in range(p.ndim):
            k = cube2par[i]
            for j in range(i,p.ndim):
                l = cube2par[j]
                cov[i,j] = (p*(axarr[k]-wmean[k])*(axarr[l]-wmean[l])).sum()
                cov[j,i] = cov[i,j]
        cov /= 1.-(p**2).sum()
        #for key in self.names:
        #    if key.find('log')==0:
        #        pntkey = key.split('log')[1]
        #        self.priors[key].value = numpy.log10(wmean[ax[pntkey]])
        #    else:
        #        self.priors[key].value = wmean[ax[key]]

        #self.priors['redshift'].value = 0.1
        pnt = numpy.empty((len(self.priors),1))
        @pymc.deterministic
        def mass_and_logp(value=0.,pars=pars):
            p = numpy.array(pars)
            p[doExp] = 10**p[doExp]
            p = numpy.atleast_2d(p[par2cube])
            mass = massmodel.eval(p)
            if mass==0.:
                return [0.,-1e200]
            logp = logpmodel.eval(p)
            return [mass,logp]

        @pymc.observed
        def loglikelihood(value=0.,lp=mass_and_logp):
            return lp[1]

        """
        logp -= logp.max()
        p = numpy.exp(logp)
        p /= p.sum()
        i = 0
        wmean = numpy.empty(p.ndim)
        for key in self.model.axes_names:
            a = self.model.axes[key]['points']
            p0 = numpy.rollaxis(p,i,p.ndim)
            wmean[i] = (a*p0).sum()
            i += 1

        

        cov = []
        for key in self.names:
            if key=='age':
                cov.append(0.5)
            elif key=='logage':
                cov.append(0.03)
            elif key=='tau':
                cov.append(0.1)
            elif key=='logtau':
                cov.append(0.03)
            elif key=='tau_V':
                cov.append(self.priors[key]['prior'].value/20.)
            elif key=='logtau_V':
                cov.append(0.1)
            elif key=='Z':
                cov.append(self.priors[key]['prior'].value/20.)
            elif key=='logZ':
                cov.append(0.03)
            elif key=='redshift':
                cov.append(0.1)
        cov = numpy.array(cov)
        """

        from SampleOpt import Sampler,AMAOpt
        S = AMAOpt(pars,[loglikelihood],[mass_and_logp],cov=cov)
        S.sample(nburn)
        logps,trace,dets = S.result()
        print logps.max()

        S = Sampler(pars,[loglikelihood],[mass_and_logp])
        S.setCov(cov)
        S.sample(nburn/2)

        logps,trace,dets = S.result()
        cov = numpy.cov(trace.T)

        S = Sampler(pars,[loglikelihood],[mass_and_logp])
        S.setCov(cov)
        S.sample(niter)

        logps,trace,dets = S.result()
        mass,logL = dets['mass_and_logp'][:,:,0].T
        o = {'logP':logps,'logL':logL,'logmass':mass}
        cnt = 0
        for key in self.names:
            o[key] = trace[:,cnt].copy()
            cnt += 1
        return o
        
        arg = logp.argmax()
        logp -= logp.max()
        p = numpy.exp(logp)
        p /= p.sum()
        print p.max()
        i = 0
        for key in self.model.axes_names:
            a = self.model.axes[key]['points']
            if key=='redshift':
                a = a[::5]
            p0 = numpy.rollaxis(p,i,p.ndim)
            print key,(a*p0).sum()
            i += 1

        print numpy.unravel_index(arg,logp.shape)
        logp -= max
        print (M*numpy.exp(logp)).sum()/numpy.exp(logp).sum()
        z = (M*0.+1)*self.model.axes['redshift']['points'][::5]
        print (z*numpy.exp(logp)).sum()/numpy.exp(logp).sum()
        f = open('check','wb')
        import cPickle
        cPickle.dump([M,logp],f,2)
        f.close()
        mod = ndinterp.ndInterp(self.models.axes,logp)




