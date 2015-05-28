import numpy

def SimpleSample(pars,costs,deterministics,niter,cov=None):
    nvars = len(pars)
    trace = numpy.empty((niter,nvars))
    logps = numpy.zeros(niter)
    dets = []
    blank = numpy.zeros(nvars)
    for varIndx in xrange(nvars):
        trace[0,varIndx] = pars[varIndx].value
        logps[0] += pars[varIndx].logp
    for cost in costs:
        logps[0] += cost.logp
    dets.append([d.value for d in deterministics]) 

    for i in xrange(1,niter):
        z = 10**(numpy.random.random(nvars)*0.8-0.5)
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

