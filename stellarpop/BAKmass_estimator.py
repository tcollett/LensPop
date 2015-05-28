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
