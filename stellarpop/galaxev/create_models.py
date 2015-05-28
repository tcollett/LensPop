
def create_csp_models(ised,tau=None,tau_V=None,mu=0.3,epsilon=0.):
    """
    Script to create B+C models for a given metallicity and a range of 
        exponential decays and dust models. Uses the following distributions
        by default:

        Tau (gyr) = U[0.04:9.1::35]
        Tau_V (reddening) = LU[0.01:2::21]

        mu (reddening) = 0.3
        epsilon (recycling) = 0 (ie no recycling)

    This code requires that the environment be setup to use B+C, including
        having the $bc03 environment variable set.
    """
    if (tau is None) or (tau_V is None):
        import numpy

    # Use default distributions if not provided
    if tau is None:
        tau = numpy.linspace(0.04,9.1,35)
    if tau_V is None:
        tau_V = numpy.logspace(-2.,numpy.log10(2.),21)

    for t in tau:
        for tV in tau_V:
            # Create the models
            _csp_model(ised,t,tV,mu,epsilon)


def _csp_model(ised,t,tV,mu,epsilon):
    """
    Worker function to create a BC03 CSP model given an ised, tau, tau_V, mu,
        and epsilon.
    """
    import os
    lookup = {'m22':0.0001,'m32':0.0004,'m42':0.004,'m52':0.008,
                'm62':0.02,'m72':0.05}

    imf = ised.split('_')[3]
    Zcode = ised.split('_')[2]
    Z = lookup[Zcode]

    tmpname = 'tmp_%s_%s.in'%(imf,Zcode)

    name = 'bc03_%s_Z=%6.4f_tV=%5.3f_mu=%3.1f_t=%5.3f_eps=%5.3f.ised' % (imf,Z,tV,mu,t,epsilon)

    file = open(tmpname,'w')

    # Metallicity/IMF
    file.write('%s\n'%ised)
    # Use dust with tau_V = tV and mu = mu
    file.write('Y\n')
    file.write('%f\n'%tV)
    file.write('%f\n'%mu)
    # Use exponential SFH with tau = t
    file.write('1\n')
    file.write('%f\n'%t)
    # Choose whether or not to recycle gas
    if epsilon!=0:
        file.write('Y\n')
        file.write('%f\n'%epsilon)
    else:
        file.write('N\n')
    # Cutoff star formation at 20 Gyr
    file.write('20\n')
    file.write('mySSP_%s\n'%Zcode)
    file.close()

    # Run bc03
    os.system('$bc03/csp_galaxev < %s'%tmpname)

    # Keep the output *.ised and *.4color files
    os.system('cp mySSP_%s.ised %s'%(Zcode,name))
    name = name.replace('ised','mass')
    os.system('cp mySSP_%s.4color %s'%(Zcode,name))

    # Clean up
    os.system('rm -f mySSP_%s*'%Zcode)


def create_seds(ised,outname,age=None,tau=None,tau_V=None,mu=0.3,epsilon=0.):
    """
    Script to create B+C SEDs over an array of age for a given metallicity and
        a range of exponential decays and dust models. If the BC03 CSP models
        have not yet been created, they will be. Uses the following
        distributions by default:

        Tau (gyr) = U[0.04:9.1::35]
        Tau_V (reddening) = LU[0.01:2::21]
        Age (gyr) = U[0.6:13.5::31]

        mu (reddening) = 0.3
        epsilon (recycling) = 0 (ie no recycling)

    This code requires that the environment be setup to use B+C, including
        having the $bc03 environment variable set.

    NOTE: The output with the default parameters uses >1GB of RAM (and a
        similar amount of disk space when the model is written to disk). The
        user is encouraged to set the age range to something reasonable for
        their system.
    """
    import os,glob,cPickle
    import numpy
    from scipy import interpolate
    from stellarpop import tools
    from math import log10

    # Use default distributions if not provided
    if tau is None:
        tau = numpy.linspace(0.04,9.1,35)
    if tau_V is None:
        tau_V = numpy.logspace(-2.,numpy.log10(2.),21)
    if age is None:
        age = numpy.linspace(0.6,13.5,31)

    lookup = {'m22':0.0001,'m32':0.0004,'m42':0.004,'m52':0.008,
                'm62':0.02,'m72':0.05}

    imf = ised.split('_')[3]
    Zcode = ised.split('_')[2]
    Z = lookup[Zcode]

    tmpname = 'tmp_%s.in'%Zcode


    # Create the mass normalization models
    files = glob.glob('bc03_%s_Z=%6.4f*mass'%(imf,Z))
    mass_models = {}
    for file in files:
        d = numpy.loadtxt(file)
        mass_models[file] = interpolate.splrep(d[:,0],d[:,6],k=3,s=0)


    output = {}
    for t in tau:
        for tV in tau_V:
            name = 'bc03_%s_Z=%6.4f_tV=%5.3f_mu=%3.1f_t=%5.3f_eps=%5.3f.mass' % (imf,Z,tV,mu,t,epsilon)
            if name not in files:
                _csp_model(ised,t,tV,mu,epsilon)

            input_ised = name.replace('mass','ised')

            # First use galaxevpl to create the SED text files
            file = open(tmpname,'w')
            for a in age:
                oname = name.replace('.mass','_age=%06.3f.sed'%a)
                # Select CSP model
                file.write('%s\n'%input_ised)
                # Use default wavelength ranges
                file.write('\n')
                # Set age
                file.write('%f\n'%a)
                # Set output name
                file.write('%s\n'%oname)
            file.close()

            os.system('$bc03/galaxevpl < %s'%tmpname)

            # Store the SEDs in output
            for a in age:
                oname = name.replace('.mass','_age=%06.3f.sed'%a)
                # Create the numpy SED
                sed = tools.makeUserSED(oname)
                wave = sed[0]
                # Renormalize the mass!
                logAge = log10(a)+9.
                mass = interpolate.splev(logAge,mass_models[name])
                sed = sed[1]/mass
                # Create the SPS object
                obj = {'sed':sed,'Z':Z,'age':a,'tau':t,'tau_V':tV,'mass':mass}
                output[oname] = obj
                os.system('rm %s'%oname)

    f = open(outname,'wb')
    cPickle.dump(output,f,2)
    cPickle.dump(wave,f,2)
    cPickle.dump(tau,f,2)
    cPickle.dump(tau_V,f,2)
    cPickle.dump(age,f,2)
    f.close()
