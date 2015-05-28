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

        Estimator.__init__(self,priors)


    def logLikelihood(self,pars):
        """
        Given a set of pars, return the loglikelihood. This function uses the
            object's model and data members. pars is a dictionary of parameter
            values, eg.

            pars = {'mass':2.4,'age':2.39,'Z':0.008,...}

        If the parameters fall out of the range of the model (such that it
            returns 0), return -Inf. [The priors should exclude this from
            happening anyways.]
        """
        from math import log10
        chi2 = 0.

        if 'mass' in pars.keys():
            massoffset = -2.5*log10(pars['mass'])
        else:
            massoffset = -2.5*pars['logmass']

        for datum in self.data:
            filter,redshift = datum['filter'],datum['redshift']
            mag,error = datum['mag'],datum['error']
            pnt = {}
            for key in self.names:
                if key=='mass' or key=='logmass':
                    continue
                if key.find('log')==0:
                    pntkey = key.split('log')[1]
                    pnt[pntkey] = [10**pars[key]]
                elif key.find('disc')==0:
                    pntkey = key.split('disc')[1]
                    pkey = pars[key]
                    if type(pkey)!=type(1):
                        pkey = pkey[0]
                    pnt[pntkey] = [self.priors[key]['vals'][pkey]]
                else:
                    pnt[key] = [pars[key]]
            modelmag = self.model.eval(pnt,filter,redshift)[0]
            # Catch results out of the bounds of the model
            if modelmag==0:
                return -1e300

            modelmag += massoffset
            chi2 += (modelmag-mag)**2/error**2
        return -0.5*chi2


    def MCMC(self,niter,nburn,nthin=1,nupdates=4):
        from Sampler import SimpleSample as sample
        import pymc,numpy

        pars = []
        cov = []
        for key in self.names:
            pars.append(self.priors[key]['prior'])
            if key=='logmass':
                cov.append(0.03)
            elif key=='mass':
                cov.append(self.priors[key]['prior'].value/20.)
            elif key=='age':
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

        cov = numpy.array(cov)
        @pymc.observed
        def loglikelihood(value=0.,pars=pars):
            p = {}
            for key in self.names:
                p[key] = pars[len(p)]
            return self.logLikelihood(p)

        costs = [loglikelihood]
        logps,trace,dets = sample(pars,costs,[],nburn,cov=cov)
        for i in range(nupdates):
            logps,trace,dets = sample(pars,costs,[],nburn/nupdates,cov=cov)
            cov = numpy.cov(trace.T)
        logps,trace,dets = sample(pars,costs,[],niter,cov=cov)
        self.proposal_cov = cov
        self.trace= trace
        self.logp = logps
        cnt = 0
        o = {'logp':logps}
        for key in self.names:
            o[key] = trace[:,cnt].copy()
            cnt += 1
        return o


    def fastMCMC(self,niter,nburn,nthin=1):
        from Sampler import SimpleSample as sample
        import pymc,numpy
        from scipy import interpolate

        pmags = []
        for i in range(len(self.data)):
            f = self.data[i]['filter']
            z = self.data[i]['redshift']
            pmags.append(self.model.models[f][z].z)
            lpmodel = self.model.models[f][z]
        for i in range(len(self.data)-1):
            d1 = self.data[i]
            d2 = self.data[i+1]
            f1,m1,me1 = d1['filter'],d1['mag'],d1['error']
            f2,m2,me2 = d2['filter'],d2['mag'],d2['error']
            c = pmags[i]-pmags[i+1]
            dc = m1-m2
            if i==0:
                logp = -0.5*(c-dc)**2/(me1**2+me2**2)
            else:
                logp += -0.5*(c-dc)**2/(me1**2+me2**2)

        indx = numpy.unravel_index(logp.argmax(),logp.shape)

        m = 0.
        w = 0.
        for i in range(len(self.data)):
            d,e = self.data[i]['mag'],self.data[i]['error']
            M = (pmags[i][indx]-d)/2.5
            m += M/e**2
            w += 1./e**2
        m /= w

        import time
        t = time.time()
        s = [i for i in logp.shape]
        s.append(11)
        logp = numpy.zeros(s)
        s[:-1] = [1 for i in range(len(s)-1)]
        M = numpy.linspace(m-1.,m+1.,11)
        for i in range(len(self.data)):
            d1 = self.data[i]
            f,mag,me = d1['filter'],d1['mag'],d1['error']
            pmag = numpy.tile(numpy.expand_dims(pmags[i],-1),s)-2.5*M
            logp += -0.5*(pmag-mag)**2/me**2

        indx = numpy.unravel_index(logp.argmax(),logp.shape)
        lpmodel.z = logp
        lpmodel.axes[len(logp.shape)-1] = interpolate.splrep(M,numpy.arange(M.size),s=3)
        lpmodel.set_order(lpmodel.order)
        print time.time()-t

        axn = self.model.axes_names
        axes = self.model.axes

        pars = []
        cov = []
        for key in self.names:
            pars.append(self.priors[key]['prior'])
            if key=='logmass':
                pars[-1].value = m
                cov.append(0.03)
            elif key=='mass':
                pars[-1].value = 10**m
                cov.append(self.priors[key]['prior'].value/20.)
            elif key=='age':
                pars[-1].value = axes['age']['points'][indx[axn.index('age')]]
                cov.append(0.5)
            elif key=='logage':
                pars[-1].value = numpy.log10(axes['age']['points'][indx[axn.index('age')]])
                cov.append(0.03)
            elif key=='tau':
                pars[-1].value = axes['tau']['points'][indx[axn.index('tau')]]
                cov.append(0.1)
            elif key=='logtau':
                pars[-1].value = numpy.log10(axes['tau']['points'][indx[axn.index('tau')]])
                cov.append(0.03)
            elif key=='tau_V':
                pars[-1].value = axes['tau_V']['points'][indx[axn.index('tau_V')]]
                cov.append(self.priors[key]['prior'].value/20.)
            elif key=='logtau_V':
                pars[-1].value = numpy.log10(axes['tau_V']['points'][indx[axn.index('tau_V')]])
                cov.append(0.1)
            elif key=='Z':
                pars[-1].value = axes['Z']['points'][indx[axn.index('Z')]]
                cov.append(self.priors[key]['prior'].value/20.)
            elif key=='logZ':
                pars[-1].value = numpy.log10(axes['Z']['points'][indx[axn.index('Z')]])
                cov.append(0.2)


        axn.append('logmass')
        cov = numpy.array(cov)
        @pymc.observed
        def loglikelihood(value=0.,pars=pars):
            from math import log10
            points = numpy.zeros((1,len(pars)))
            i = 0
            for key in self.names:
                if key=='mass':
                    points[0,axn.index('logmass')] = log10(pars[i])
                elif key=='logmass':
                    points[0,axn.index('logmass')] = pars[i]
                elif key.find('log')==0:
                    key = key.split('log')[1]
                    points[0,axn.index(key)] = 10**pars[i]
                else:
                    points[0,axn.index(key)] = pars[i]
                i += 1
            lp = lpmodel.eval(points)
            if lp==0:
                return -1e300
            return lp            
            return self.logLikelihood(p)


        costs = [loglikelihood]
        logps,trace,dets = sample(pars,costs,[],nburn,cov=cov)
        cov = numpy.cov(trace[nburn/2:].T)
        print [cov[i,i]**0.5 for i in range(cov.shape[0])]
        logps,trace,dets = sample(pars,costs,[],niter,cov=cov)
        self.proposal_cov = cov
        self.trace= trace
        self.logp = logps
        cnt = 0
        o = {'logp':logps}
        for key in self.names:
            o[key] = trace[:,cnt].copy()
            cnt += 1
        return o



    def fasterMCMC(self,niter,nburn,nthin=1,grid=False):
        from Sampler import SimpleSample as sample
        import pymc,numpy,ndinterp

        axn = [a for a in self.model.axes_names]
        axes = self.model.axes

        pmags = []
        lpModel = None
        magModel = None
        Modaxes = None
        magUC = 0.
        order = None
        for i in range(len(self.data)):
            f = self.data[i]['filter']
            z = self.data[i]['redshift']
            pmags.append(self.model.models[f][z].z)
            if magModel is None:
                magModel = self.model.models[f][z].z*0.
                lpModel = magModel.copy()
                Modaxes = self.model.models[f][z].axes
                order = self.model.models[f][z].order
            m,me = self.data[i]['mag'],self.data[i]['error']**2
            d = (m-pmags[-1])
            lpModel += -0.5*d**2/me
            magModel += d/me
            magUC += -0.5/me

        for i in range(len(self.data)-1):
            d1 = self.data[i]
            d2 = self.data[i+1]
            f1,m1,me1 = d1['filter'],d1['mag'],d1['error']
            f2,m2,me2 = d2['filter'],d2['mag'],d2['error']
            c = pmags[i]-pmags[i+1]
            dc = m1-m2
            if i==0:
                logp = (c-dc)**2/(me1**2+me2**2)
            else:
                logp += (c-dc)**2/(me1**2+me2**2)

        indx = numpy.unravel_index(logp.argmin(),logp.shape)
        m = 0.
        w = 0.
        for i in range(len(self.data)):
            d,e = self.data[i]['mag'],self.data[i]['error']
            M = (pmags[i][indx]-d)/2.5
            m += M/e**2
            w += 1./e**2
        m /= w

        M = numpy.linspace(m-0.6,m+0.6,13)
        a = []
        for i in range(len(self.model.axes_names)):
            a.append(self.model.axes[self.model.axes_names[i]]['points'].copy())
        for key in self.names:
            p = self.priors[key]['prior']
            if key.find('mass')>=0:
                continue
            if key.find('log')==0:
                key = key[3:]
                a[axn.index(key)] = numpy.log10(a[axn.index(key)])
            i = axn.index(key)
            for j in range(len(a[i])):
                p.value = a[i][j]
                try:
                    a[i][j] = p.logp
                except:
                    a[i][j] = -1e300
        logp = lpModel+ndinterp.create_axes_array(a).sum(0)

        logp = numpy.expand_dims(logp,logp.ndim).repeat(M.size,logp.ndim)
        for i in range(M.size):
            M0 = M[i]*-2.5
            logp[:,:,:,:,i] += magModel*M0 +  M0**2*magUC
        logp -= logp.max()

        wt = numpy.exp(logp)
        wt /= wt.sum()
        a = []
        for i in range(len(self.model.axes_names)):
            a.append(self.model.axes[self.model.axes_names[i]]['points'])
        a.append(M)
        for key in self.names:
            if key.find('mass')>=0:
                if key.find('log')!=0:
                    a[-1] == 10**a[-1]
                axn.append('mass')
            elif key.find('log')==0:
                key = key[3:]
                a[axn.index(key)] = numpy.log10(a[axn.index(key)])
        vals = ndinterp.create_axes_array(a)
        if grid==True:
            m = (wt*vals[-1]).sum()
            st = ((wt*(vals[-1]-m)**2).sum())**0.5
            return m,st

        pars = []
        cov = []

        for key in self.names:
            pars.append(self.priors[key]['prior'])
            if key.find('log')==0:
                key = key[3:]
            i = axn.index(key)
            m = (wt*vals[i]).sum()
            st = ((wt*(vals[i]-m)**2).sum())**0.5
            pars[-1].value = m
            cov.append(st)

        for key in self.names:
            continue
            pars.append(self.priors[key]['prior'])
            if key=='logmass':
                pars[-1].value = m
                cov.append(0.03)
            elif key=='mass':
                pars[-1].value = 10**m
                cov.append(self.priors[key]['prior'].value/20.)
            else:
                if key.find('log')!=0:
                    i = axn.index(key)
                    mn = (vals[i]*wt).sum()
                    cov.append(((wt*(vals[i]-mn)**2).sum())**0.5)
                    pars[-1].value = mn
                else:
                    key = key.split('log')[1]
                    i = axn.index(key)
                    mn = (numpy.log10(vals[i])*wt).sum()
                    cov.append(((wt*(numpy.log10(vals[i])-mn)**2).sum())**0.5)
                    pars[-1].value = mn
        for key in self.names:
            continue
            pars.append(self.priors[key]['prior'])
            if key=='logmass':
                pars[-1].value = m
                cov.append(0.03)
            elif key=='mass':
                pars[-1].value = 10**m
                cov.append(self.priors[key]['prior'].value/20.)
            elif key=='age':
                pars[-1].value = axes['age']['points'][indx[axn.index('age')]]
                cov.append(0.5)
            elif key=='logage':
                pars[-1].value = numpy.log10(axes['age']['points'][indx[axn.index('age')]])
                cov.append(0.03)
            elif key=='tau':
                pars[-1].value = axes['tau']['points'][indx[axn.index('tau')]]
                cov.append(0.2)
            elif key=='logtau':
                pars[-1].value = numpy.log10(axes['tau']['points'][indx[axn.index('tau')]])
                cov.append(0.03)
            elif key=='tau_V':
                pars[-1].value = axes['tau_V']['points'][indx[axn.index('tau_V')]]
                cov.append(self.priors[key]['prior'].value/20.)
            elif key=='logtau_V':
                pars[-1].value = numpy.log10(axes['tau_V']['points'][indx[axn.index('tau_V')]])
                cov.append(0.1)
            elif key=='Z':
                pars[-1].value = axes['Z']['points'][indx[axn.index('Z')]]
                cov.append(self.priors[key]['prior'].value/20.)
            elif key=='logZ':
                pars[-1].value = numpy.log10(axes['Z']['points'][indx[axn.index('Z')]])
                cov.append(0.05)

        lpModel = ndinterp.ndInterp(Modaxes,lpModel,order)
        magModel = ndinterp.ndInterp(Modaxes,magModel,order)

        cov = numpy.array(cov)
        @pymc.observed
        def loglikelihood(value=0.,pars=pars):
            from math import log10
            points = numpy.zeros((1,len(pars)-1))
            i = 0
            for key in self.names:
                if key=='mass':
                    M = -2.5*log10(pars[i])
                elif key=='logmass':
                    M = -2.5*pars[i]
                elif key.find('log')==0:
                    key = key.split('log')[1]
                    points[0,axn.index(key)] = 10**pars[i]
                else:
                    points[0,axn.index(key)] = pars[i]
                i += 1
            lp = lpModel.eval(points)
            if lp==0:
                return -1e300  # Short circuit if out of range
            lp += magModel.eval(points)*M + M**2*magUC
            return lp


        costs = [loglikelihood]
        logps,trace,dets = sample(pars,costs,[],nburn/2,cov=cov,jump=[0.,0.])
#        cov = numpy.cov(trace.T)
        logps,trace,dets = sample(pars,costs,[],nburn/2,cov=cov,jump=[0.,0.])
        cov = numpy.cov(trace.T)
        logps,trace,dets = sample(pars,costs,[],niter,cov=cov,jump=[0.,0.])
        self.proposal_cov = cov
        self.trace= trace
        self.logp = logps
        cnt = 0
        o = {'logp':logps}
        for key in self.names:
            o[key] = trace[:,cnt].copy()
            cnt += 1
        return o

