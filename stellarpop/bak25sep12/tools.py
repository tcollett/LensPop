"""
A module of tools to work with SED, spectral, and filter operations.
"""
import stellarpop as sp
filterpath = sp.__path__[0]+'/filters/'
SEDpath = sp.__path__[0]+'/templates/'

import glob
filterList = glob.glob(filterpath+"*res")
filterList.sort()
filterList = [filter.split('/')[-1].split('.res')[0] for filter in filterList]
SEDList = glob.glob(SEDpath+"*sed")
SEDList.sort()
SEDList = [SED.split('/')[-1].split('.sed')[0] for SED in SEDList]

from math import pi
BC03factor = 3.826e33/(4*pi*3.08568e19**2)

def filterfromfile(file):
    """
    Create a filter model from a file.
    """
    try:
        from numpy.lib.io import loadtxt
    except:
        from numpy import loadtxt
    from scipy.interpolate import splrep
    import scipy
    f = open(filterpath+file+".res")
    filter = loadtxt(f)
    f.close()
    return splrep(filter[:,0],filter[:,1],k=1,s=0)


def filterfrommodel(name):
    """
    Returns a pre-computed spline filter model for filter `name.'
    """
    import cPickle
    f = open(filterpath+name+".model")
    filter = cPickle.load(f)
    f.close()
    return filter


def getSED(name):
    """
    Returns a model of the SED, a tuple of (wave,data)
    """
    try:
        from numpy.lib.io import loadtxt
    except:
        from numpy import loadtxt
    f = open(SEDpath+name+".sed")
    SED = loadtxt(f)
    f.close()
    return SED[:,0],SED[:,1]


def makeUserSED(filename):
    """
    Returns a model of an SED supplied by the user
    """
    from numpy.lib.io import loadtxt
    SED = loadtxt(filename)
    return SED[:,0],SED[:,1]


def getUserSED(filename):
    """
    Loads a (binary) model of an SED from a pickle
    """
    import cPickle
    f = open(filename,'rb')
    sed = cPickle.load(f)
    f.close()
    return sed


def writeSEDtoFile(sed,filename):
    """
    Writes a binary pickle of the SED to file
    """
    import cPickle
    f = open(filename,'wb')
    cPickle.dump(sed,f,2)
    f.close()


def ABFilterMagnitude(filter,spectrum,redshift):
    """
    Determines the AB magnitude (up to a constant) given an input filter, SED,
        and redshift.
    """
    from scipy.interpolate import splev,splint,splrep
    from scipy.integrate import simps
    from math import log10
    sol = 299792452.

    wave = spectrum[0].copy()
    data = spectrum[1].copy()

    # Convert to f_nu
    data = data*wave**2/(sol*1e10)

    # Redshift the spectrum and determine the valid range of wavelengths
    wave *= (1.+redshift)
    wmin,wmax = filter[0][0],filter[0][-1]
    cond = (wave>=wmin)&(wave<=wmax)

    # Evaluate the filter at the wavelengths of the spectrum
    response = splev(wave[cond],filter)

    freq = sol*1e10/wave[cond]
    data = data[cond]*(1.+redshift)

    # Flip arrays
    freq = freq[::-1]
    data = data[::-1]
    response = response[::-1]

    # Integrate
    observed = splrep(freq,response*data/freq,s=0,k=1)
    flux = splint(freq[0],freq[-1],observed)

    bp = splrep(freq,response/freq,s=0,k=1)
    bandpass = splint(freq[0],freq[-1],bp)

    return -2.5*log10(flux/bandpass) - 48.6
ABFM = ABFilterMagnitude

def VegaFilterMagnitude(filter,spectrum,redshift):
    """
    Determines the Vega magnitude (up to a constant) given an input filter,
        SED, and redshift.
    """
    from scipy.interpolate import splev,splint,splrep
    from scipy.integrate import simps
    from math import log10

    wave = spectrum[0].copy()
    data = spectrum[1].copy()

    # Redshift the spectrum and determine the valid range of wavelengths
    wave *= (1.+redshift)
    data /= (1.+redshift)
    wmin,wmax = filter[0][0],filter[0][-1]
    cond = (wave>=wmin)&(wave<=wmax)

    # Evaluate the filter at the wavelengths of the spectrum
    response = splev(wave[cond],filter)

    # Determine the total observed flux (without the bandpass correction)
    observed = splrep(wave[cond],(response*data[cond]),s=0,k=1)
    flux = splint(wmin,wmax,observed)

    # Determine the magnitude of Vega through the filter
    vwave,vdata = getSED('Vega')
    cond = (vwave>=wmin)&(vwave<=wmax)
    response = splev(vwave[cond],filter)
    vega = splrep(vwave[cond],response*vdata[cond],s=0,k=1)
    vegacorr = splint(wmin,wmax,vega)

    return -2.5*log10(flux/vegacorr)#+2.5*log10(1.+redshift)


def FilterMagnitude(filter,spectrum,redshift,zp):
    from scipy.interpolate import splev,splint,splrep
    from math import log10

    wave = spectrum[0].copy()
    data = spectrum[1].copy()

    # Redshift the spectrum and determine the valid range of wavelengths
    wave *= (1.+redshift)
    data /= (1.+redshift)
    wmin,wmax = filter[0][0],filter[0][-1]
    cond = (wave>=wmin)&(wave<=wmax)

    # Evaluate the filter at the wavelengths of the spectrum
    response = splev(wave[cond],filter)

    # Determine the total observed flux (without the bandpass correction)
    observed = splrep(wave[cond],(response*data[cond]),s=0,k=1)
    flux = splint(wmin,wmax,observed)
    flux = (response*data[cond]).sum()
    return -2.5*log10(flux) + zp

