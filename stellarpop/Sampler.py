import numpy

def SimpleSample(pars,costs,deterministics,niter,cov=None,jump=None):
    if jump is None:
        stretch,offset = 3.3,3.
    else:
        stretch,offset = jump
    nvars = len(pars)
    niter = int(niter)
    trace = numpy.empty((niter,nvars))
    logps = numpy.zeros(niter)
    dets = []
    if cov is None:
        widths = {'x':0.05,'y':0.05,'reff':0.1,'q':0.03,'pa':1.,'eta':0.03,'nu':0.03}
        cov = numpy.empty(nvars)
        for varIndx in xrange(nvars):
            name = pars[varIndx].__name__
            cov[varIndx] = widths[name.split('_')[0]]
    else:
        cov = numpy.asarray(cov)
    blank = numpy.zeros(nvars)
    for varIndx in xrange(nvars):
        trace[0,varIndx] = pars[varIndx].value
        logps[0] += pars[varIndx].logp
    for cost in costs:
        logps[0] += cost.logp
    dets.append([d.value for d in deterministics]) 

    for i in xrange(1,niter):
        z = 10**(numpy.random.random(nvars)*stretch-offset)
        if cov.ndim==1:
            W = numpy.random.randn(cov.size)*cov*z
        else:
            W = numpy.random.multivariate_normal(blank,cov)*z
        logp = 0.
        updates = trace[i-1].copy()+W
        bad = False
        for varIndx in xrange(nvars):
            pars[varIndx].value = updates[varIndx]
        for varIndx in xrange(nvars):
            try:
                logp += pars[varIndx].logp
            except:
                logp += -1e200
                bad = True
                break
        if bad==True:
            logps[i] = logps[i-1]
            trace[i] = trace[i-1].copy()
            dets.append(dets[-1])
            continue
        for cost in costs:
            logp += cost.logp
        if logp>logps[i-1]:
            logps[i] = logp
            trace[i] = updates
            dets.append([d.value for d in deterministics])
            continue
        if logp-logps[i-1]>numpy.log(numpy.random.random()):
            logps[i] = logp
            trace[i] = updates
            dets.append([d.value for d in deterministics])
        else:
            logps[i] = logps[i-1]
            trace[i] = trace[i-1].copy()
            dets.append(dets[-1])
    for varIndx in xrange(nvars):
        pars[varIndx].value = trace[-1][varIndx]

    return logps,trace,dets


def Optimizer(pars,costs,deterministics,niter,cov=None):
    from scipy.optimize import leastsq
    from numpy import exp

    nvars = len(pars)
    fake = numpy.ones(nvars*10)
    def optFunc(p):
        post = 0.
        for i in range(p.size):
            try:
                pars[i].value = p[i]
                post += pars[i].logp
            except:
                print 'blah'
                return fake*1e10
        for cost in costs:
            post += cost.logp
        return exp(post*-1)*fake

    inpar = numpy.empty(nvars)
    for i in range(nvars):
        inpar[i] = pars[i].value

    outpar,ier = leastsq(optFunc,inpar,epsfcn=1e-4)
    return outpar


def MCMCOpt(inpars,costs,deterministics,niter,cov=None,jump=None):
    if jump is None:
        stretch,offset = 3.3,3
    else:
        from math import log10
        lo,hi = jump
        lo,hi = log10(lo),log10(hi)
        stretch,offset = hi+lo,lo
    pars = []
    for par in inpars:
        try:
            tmp = par.logp
            pars.append(par)
        except:
            pass
    nvars = len(pars)
    trace = numpy.empty((niter,nvars))
    logps = numpy.zeros(niter)
    dets = []
    if cov is None:
        widths = {'x':0.05,'y':0.05,'re':0.1,'q':0.03,'pa':1.,'eta':0.03,'nu':0.03}
        cov = numpy.empty(nvars)
        for varIndx in xrange(nvars):
            name = pars[varIndx].__name__
            cov[varIndx] = widths[name.split('_')[0]]
    blank = numpy.zeros(nvars)
    for varIndx in xrange(nvars):
        trace[0,varIndx] = pars[varIndx].value
        logps[0] += pars[varIndx].logp
    for cost in costs:
        logps[0] += cost.logp
    dets.append([d.value for d in deterministics])

    for i in xrange(1,niter):
#        z = 10**(numpy.random.random(nvars)*stretch-offset)
        z = 10**(numpy.random.randn(nvars)*0.3)
        if cov.ndim==1:
            W = numpy.random.randn(cov.size)*cov*z
        else:
            W = numpy.random.multivariate_normal(blank,cov)*z
        logp = 0.
        updates = trace[i-1].copy()+W
        bad = False
        for varIndx in xrange(nvars):
            pars[varIndx].value = updates[varIndx]
        for varIndx in xrange(nvars):
            try:
                logp += pars[varIndx].logp
            except:
                logp = -1e300
                bad = True
                break
        if bad==True:
            logps[i] = logps[i-1]
            trace[i] = trace[i-1].copy()
            dets.append(dets[-1])
            continue

        for cost in costs:
            logp += cost.logp
        if logp>logps[i-1]:
            logps[i] = logp
            trace[i] = updates
            dets.append([d.value for d in deterministics])
            continue
        else:
            logps[i] = logps[i-1]
            trace[i] = trace[i-1].copy()
            dets.append(dets[-1])
    for varIndx in xrange(nvars):
        pars[varIndx].value = trace[-1][varIndx]

    return logps,trace,dets



    for i in xrange(1,niter):
        z = 10**(numpy.random.random(nvars)*3.3-3.)
        if cov.ndim==1:
            W = numpy.random.randn(cov.size)*cov*z
        else:
            W = numpy.random.multivariate_normal(blank,cov)*z
        logp = 0.
        updates = trace[i-1].copy()+W
        bad = False
        for varIndx in xrange(nvars):
            try:
                pars[varIndx].value = updates[varIndx]
                logp += pars[varIndx].logp
            except:
                logp += -1e200
                bad = True
                break
        if bad==True:
            logps[i] = logps[i-1]
            trace[i] = trace[i-1].copy()
            dets.append(dets[-1])
            continue
        for cost in costs:
            logp += cost.logp
        if logp>logps[i-1]:
            logps[i] = logp
            trace[i] = updates
            dets.append([d.value for d in deterministics])
            continue
        if logp-logps[i-1]>numpy.log(numpy.random.random()):
            logps[i] = logp
            trace[i] = updates
            dets.append([d.value for d in deterministics])
        else:
            logps[i] = logps[i-1]
            trace[i] = trace[i-1].copy()
            dets.append(dets[-1])
    for varIndx in xrange(nvars):
        pars[varIndx].value = trace[-1][varIndx]

    return logps,trace,dets

