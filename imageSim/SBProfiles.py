import numpy,time
from scipy import interpolate
"""
MINIMAL ERROR CHECKING!
"""
def cnts2mag(cnts,zp):
    from math import log10
    return -2.5*log10(cnts) + zp


class Sersic:
    def __init__(self,x=None,y=None,q=None,pa=None,re=None,amp=None,n=None):
        self.x = x
        self.y = y
        self.q = q
        self.pa = pa
        self.re = re
        self.amp = amp
        self.n = n
        self.convolve = True

        #A flag to tell the code not to pixeval every step, if nothing changes
        self.NoFreeParams=False

    def setAmpFromMag(self,mag,zp):
        from math import exp,log10,pi
        from scipy.special import gamma
        cnts = 10**(-0.4*(mag-zp))
        n = self.n
        re = self.re
        k = 2.*n-1./3+4./(405.*n)+46/(25515.*n**2)
        self.amp = cnts/((re**2)*exp(k)*n*(k**(-2*n))*gamma(2*n)*2*pi)

    def eval(self,r):
        k = 2.*self.n-1./3+4./(405.*self.n)+46/(25515.*self.n**2)
        R = r/self.re
        return self.amp*numpy.exp(-k*(R**(1./self.n) - 1.))

    def pixeval(self,x,y,scale=1,csub=23):
        if self.NoFreeParams==True:
            try:
                return self.pixevalresult
            except:
                pass

        from scipy import interpolate
        from math import pi,cos as COS,sin as SIN
        shape = x.shape
        x = x.ravel()
        y = y.ravel()

        cos = COS(self.pa*pi/180.)
        sin = SIN(self.pa*pi/180.)
        xp = (x-self.x)*cos+(y-self.y)*sin
        yp = (y-self.y)*cos-(x-self.x)*sin
        r = (self.q*xp**2+yp**2/self.q)**0.5

        k = 2.*self.n-1./3+4./(405.*self.n)+46/(25515.*self.n**2)
        R = numpy.logspace(-5.,4.,451) # 50 pnts / decade
        s0 = numpy.exp(-k*(R**(1./self.n) - 1.))

        # Determine corrections for curvature
        rpow = R**(1./self.n - 1.)
        term1 = (k*rpow/self.n)**2
        term2 = k*(self.n-1.)*rpow/(R*self.n**2)
        wid = scale/self.re
        corr = (term1+term2)*wid**3/6.
        try:
            minR = R[abs(corr)<0.005].min()
        except:
            minR = 0

        # Evaluate model!
        model = interpolate.splrep(R,s0,k=3,s=0)
        model2 = interpolate.splrep(R,s0*R*self.re**2,k=3,s=0)
        R0 = r/self.re
        s = interpolate.splev(R0,model)*scale**2
        if self.n<=1. or minR==0:
            return self.amp*s.reshape(shape)
        coords = numpy.where(R0<minR)[0]
        c = (numpy.indices((csub,csub)).astype(numpy.float32)-csub/2)*scale/csub
        for i in coords:
            # The central pixels are tricky because we can't assume that we
            #   are integrating in delta-theta segments of an annulus; these
            #   pixels are treated separately by sub-sampling with ~500 pixels
            if R0[i]<3*scale/self.re:
                s[i] = 0.
                y0 = c[1]+y[i]
                x0 = c[0]+x[i]
                xp = (x0-self.x)*cos+(y0-self.y)*sin
                yp = (y0-self.y)*cos-(x0-self.x)*sin
                r0 = (self.q*xp**2+yp**2/self.q)**0.5/self.re
                s[i] = interpolate.splev(r0.ravel(),model).mean()*scale**2
                continue
            lo = R0[i]-0.5*scale/self.re
            hi = R0[i]+0.5*scale/self.re
            angle = (scale/self.re)/R0[i]
            s[i] = angle*interpolate.splint(lo,hi,model2)
            # The following code should no longer be needed
            """
            if lo<0:
                s[i] = ((interpolate.splint(0,abs(lo),model2)+interpolate.splint(0,hi,model2)))*pi*2
            else:
                s[i] = angle*interpolate.splint(lo,hi,model2)
            """

        self.pixevalresult=self.amp*s.reshape(shape)

        return self.pixevalresult


class deV(Sersic):

    def __init__(self,x=None,y=None,q=None,pa=None,re=None,amp=None):
        Sersic.__init__(self,x,y,q,pa,re,amp,4.)


class exp(Sersic):

    def __init__(self,x=None,y=None,q=None,pa=None,re=None,amp=None):
        Sersic.__init__(self,x,y,q,pa,re,amp,1.)


class Gauss:

    def __init__(self,x=None,y=None,q=None,pa=None,sigma=None,amp=None,r0=None):
        self.x = x
        self.y = y
        self.q = q
        self.pa = pa
        self.sigma = sigma
        self.amp = amp
        self.r0 = r0
        self.convolve = True

    def pixeval(self,x,y,factor=None,csub=None):
        from math import pi

        cos = numpy.cos(self.pa*pi/180.)
        sin = numpy.sin(self.pa*pi/180.)
        xp = (x-self.x)*cos+(y-self.y)*sin
        yp = (y-self.y)*cos-(x-self.x)*sin
        r2 = (self.q*xp**2+yp**2/self.q)
        if self.r0 is None:
            return self.amp*numpy.exp(-0.5*r2/self.sigma**2)
        return self.amp*numpy.exp(-0.5*(r2**0.5-self.r0)**2/self.sigma**2)


    def getMag(self,zp):
        from math import exp,pi
        if self.r0 is None:
            cnts = self.amp*(2*pi*self.sigma**2)
        else:
            from scipy.special import erf
            r0 = self.r0
            s = self.sigma
            r2pi = (2*pi)**0.5
            cnts = self.amp*pi*s*(r2pi*r0*(1.+erf(r0/(s*2**0.5)))+2*s*exp(-0.5*r0**2/s**2))
        return cnts2mag(cnts,zp)


    def eval(self,x,y):
        from math import pi
        try:
            cos = numpy.cos(self.theta)
            sin = numpy.sin(self.theta)
            xp = (x-self.x)*cos+(y-self.y)*sin
            yp = (y-self.y)*cos-(x-self.x)*sin
            r = (self.q*xp**2+yp**2/self.q)**0.5/self.sigma
            s = self.amp*numpy.exp(-0.5*r**self.n)/(2.*pi*self.sigma**2)**1.0
            return s
        except:
            return x*0.

