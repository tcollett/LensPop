import numpy,time
from scipy import interpolate,ndimage
"""
MINIMAL ERROR CHECKING!
"""
def cnts2mag(cnts,zp):
    from math import log10
    return -2.5*log10(cnts) + zp


class PixelizedModel:
    def __init__(self,image):
        self.image = image.copy()
        self.image /= image.sum()
        self.setCentroid()
        self.createModel()
        self.x = self.x0
        self.y = self.y0
        self.amp = 1.
        self.convolve = None

    def setCentroid(self):
        y,x = numpy.indices(self.image.shape).astype(numpy.float32)
        self.x0 = (x*self.image).sum()
        self.y0 = (y*self.image).sum()

    def createModel(self,order=1):
        if order==1:
            self.model = self.image.copy()
        else:
            self.model = ndimage.spline_filter(self.image,output=numpy.float64,order=order)
        self.order = order

    def pixeval(self,x,y):
        X = x-self.x+self.x0
        Y = y-self.y+self.y0
        psf = ndimage.map_coordinates(self.model,[Y,X],prefilter=False)
        psf /= psf.sum()
        return self.amp*psf


class GaussianModel:
    def __init__(self,parameters,isSDSS=False):
        self.x = 0.
        self.y = 0.
        self.amp = 1.
        self.pars = parameters
        self.isSDSS = isSDSS
        self.convolve = None

    def pixeval(self,x,y,window=10):
        X = x-self.x
        X2 = X**2
        Y = y-self.y
        Y2 = Y**2
        r2 = X2+Y2
        sig2x,sig2y,sigxy,corr,FWHM = self.pars
        c = (r2<(window*FWHM)**2)
        psf = x*0.
        psf[c] = numpy.exp(-0.5*(X2[c]/sig2x + Y2[c]/sig2y - 2.*corr**2*X[c]*Y[c]/sigxy))
        psf /= psf.sum()
        return self.amp*psf
