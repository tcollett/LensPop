from ndinterp import ndInterp
from math import pi

class MassProfile:
    """
    A generic mass model class, including a method to calculate deflection
        angles for power law models. This probably needs to be re-worked.
    """
    def align_coords(self,xin,yin,revert=False):
        from math import cos,sin,pi
        theta = self.theta-pi/2.
        ctheta = cos(theta)
        stheta = sin(theta)
        if revert:
            X = xin*ctheta-yin*stheta
            Y = yin*ctheta+xin*stheta
            x = X+self.x
            y = Y+self.y
            return x,y
        X = xin-self.x
        Y = yin-self.y
        x = X*ctheta+Y*stheta
        y = Y*ctheta-X*stheta
        return x,y


class PowerLaw(MassProfile):
    """
    A subclass for power-law mass models. The `power-law' aspect doesn't
        currently work, but it does work for SIE models.
    """
    def __init__(self,b=None,eta=1.,pa=None,q=None,x=None,y=None,load=False):
        self.b = b
        self.eta = eta
        self.pa = pa
        self.q = q
        self.x = x
        self.y = y
        self.nsteps = 1e3
        if load:
            import cPickle
            f = open('powerlaw.alphax','rb')
            self.xmodel = cPickle.load(f)
            f.close()
            f = open('powerlaw.alphay','rb')
            self.ymodel = cPickle.load(f)
            f.close()
        else:
            self.xmodel=None
            self.ymodel=None
        if eta!=1.:
            self.setGrid()


    def kappa(self,rho):
        return 0.5*rho**(self.eta-2.)

    def deflections(self,xin,yin):
        from numpy import ones,arctan as atan, arctanh as atanh
        from math import cos,sin,pi
        from numpy import arcsin as asin,arcsinh as asinh
        x,y = self.align_coords(xin,yin)
        q = self.q
        if q==1.:
            q = 1.-1e-7  # Avoid divide-by-zero errors
        eps = (1.-q**2)**0.5
        if self.eta==1.:
            # SIE models
            r = (x**2+y**2)**0.5
            r[r==0.] = 1e-9
            xout = self.b*asinh(eps*x/q/r)*q**0.5/eps
            yout = self.b*asin(eps*y/r)*q**0.5/eps
        else:
            from powerlaw import powerlawdeflections as pld
            b,eta = self.b,self.eta
            s = 1e-7
            if x.ndim>1:
                yout,xout = pld(-1*y.ravel(),x.ravel(),b,eta,s,q)
                xout,yout = xout.reshape(x.shape),-1*yout.reshape(y.shape)
            else:
                yout,xout = pld(-1*y,x,b,eta,s,q)
                yout = -1*yout
        theta = -(self.theta - pi/2.)
        ctheta = cos(theta)
        stheta = sin(theta)
        x = xout*ctheta+yout*stheta
        y = yout*ctheta-xout*stheta
        return x,y

    def caustic(self):
        if self.eta!=1:
            return None,None
        from numpy import ones,arctan as atan, arctanh as atanh
        from numpy import cos,sin,pi,linspace
        from numpy import arcsin as asin,arcsinh as asinh
        q = self.q
        if q==1.:
            q = 1.-1e-7  # Avoid divide-by-zero errors
        eps = (1.-q**2)**0.5
        theta = linspace(0,2*pi,5000)
        ctheta = cos(theta)
        stheta = sin(theta)
        delta = (ctheta**2+q**2*stheta**2)**0.5
        xout = ctheta*q**0.5/delta - asinh(eps*ctheta/q)*q**0.5/eps
        yout = stheta*q**0.5/delta - asin(eps*stheta)*q**0.5/eps
        xout,yout = xout*self.b,yout*self.b
        theta = -(self.theta - pi/2.)
        ctheta = cos(theta)
        stheta = sin(theta)
        x = xout*ctheta+yout*stheta
        y = yout*ctheta-xout*stheta
        return x+self.x,y+self.y



class ExtShear(MassProfile):

    def __init__(self,x=None,y=None,b=0.,theta=0.,pa=None):
        self.x = x
        self.y = y
        self.b = b
        if pa is None:
            self.theta = theta
        else:
            self.pa = pa

    def deflections(self,x,y):
        from math import sin,cos,pi
        x = x-self.x
        y = y-self.y
        theta = self.theta #- pi/2
        s = sin(2*theta)
        c = cos(2*theta)

        # From Kormann B1422 paper
        alpha_x = -self.b*(x*c+y*s)
        alpha_y = -self.b*(x*s-y*c)

        return alpha_x,alpha_y


class PointSource(MassProfile):
    def __init__(self,x=None,y=None,b=0.):
        self.x = x
        self.y = y

    def deflections(self,x,y):
        from math import pi 
        x = x-self.x
        y = y-self.y
        r2 = x**2+y**2
        b2 = self.b**2
        alpha_x = b2*x/r2
        alpha_y = b2*y/r2

        return alpha_x,alpha_y


class Sersic(MassProfile):
    def __init__(self,b=None,n=None,pa=None,q=None,x=None,y=None,re=None):
        self.b = b
        self.n = n
        self.pa = pa
        self.q = q
        self.x = x
        self.y = y
        self.re = re

    def deflections(self,xin,yin):
        from numpy import ones,arctan as atan, arctanh as atanh
        from math import cos,sin,pi,exp
        from numpy import arcsin as asin,arcsinh as asinh
        from scipy.special import gammainc,gamma
        x,y = self.align_coords(xin,yin)
        q = self.q
        if q==1.:
            q = 1.-1e-7  # Avoid divide-by-zero errors
        from sersic import sersicdeflections as sd
        b,n,re = self.b,self.n,self.re
        k = 2.*n-1./3+4./(405.*n)+46/(25515.*n**2)
        b = b/q**0.5
        re = re/q**0.5
        amp = b*b*k**(2*n)/(2*n*re**2*(gamma(2*n)*gammainc(2*n,k*(b/re)**(1/n))))
        if x.ndim>1:
            yout,xout = sd(-1*y.ravel(),x.ravel(),amp,re,n,q)
            xout,yout = xout.reshape(x.shape),-1*yout.reshape(y.shape)
        else:
            yout,xout = sd(-1*y,x,amp,re,n,q)
            yout = -1*yout
        theta = -(self.theta - pi/2.)
        ctheta = cos(theta)
        stheta = sin(theta)
        x = xout*ctheta+yout*stheta
        y = yout*ctheta-xout*stheta
        return x,y

    def getMass(self,sigCrit=1.):
        from math import pi
        from scipy.special import gammainc,gamma

        b,n,re = self.b,self.n,self.re
        k = 2.*n-1./3+4./(405.*n)+46/(25515.*n**2)
        amp = b*b*k**(2*n)/(2*n*re**2*(gamma(2*n)*gammainc(2*n,k*(b/re)**(1/n))))
        return 2*pi*sigCrit*re**2*amp*gamma(2*n)*n/k**(2*n)

    def getbFromMass(self,mass,sigCrit):
        from scipy import interpolate
        from scipy.special import gammainc,gamma
        import numpy

        n,re = self.n,self.re
        b = numpy.logspace(-3,1,401)*re
        k = 2.*n-1./3+4./(405.*n)+46/(25515.*n**2)
        amp = b*b*k**(2*n)/(2*n*re**2*(gamma(2*n)*gammainc(2*n,k*(b/re)**(1/n))))
        m = 2*pi*sigCrit*re**2*amp*gamma(2*n)*n/k**(2*n)
        model = interpolate.splrep(m,b)
        return interpolate.splev(mass,model)

    def setbFromMass(self,mass,sigCrit):
        self.b = self.getbFromMass(mass,sigCrit)



class SersicG(MassProfile):
    def __init__(self,b=None,n=None,pa=None,q=None,x=None,y=None,re=None):
        import numpy
        import pylens
        self.b = b
        self.n = n
        self.pa = pa
        self.q = q
        self.x = x
        self.y = y
        self.re = re
        path = pylens.__file__.split('pylens.py')[0]
        self.xmod,self.ymod = numpy.load('%s/serModelsHDR.dat'%path)

    def deflections(self,xin,yin):
        import numpy
        from math import cos,sin,pi
        from scipy.special import gammainc
        x,y = self.align_coords(xin,yin)
        q = self.q
        if q==1.:
            q = 1.-1e-7  # Avoid divide-by-zero errors
        b,n,re = self.b,self.n,self.re
        k = 2.*n-1./3+4./(405.*n)+46/(25515.*n**2)
        amp = (b/re)**2/gammainc(2*n,k*(b/re)**(1/n))
        r0 = re/q**0.5
        o = numpy.ones(x.size).astype(x.dtype)
        eval = numpy.array([abs(x).ravel()/r0,abs(y).ravel()/r0,n*o,q*o]).T
        xout = numpy.sign(x)*(amp*self.xmod.eval(eval)*r0).reshape(x.shape)
        yout = numpy.sign(y)*(amp*self.ymod.eval(eval)*r0).reshape(y.shape)
        theta = -(self.theta - pi/2.)
        ctheta = cos(theta)
        stheta = sin(theta)
        x = xout*ctheta+yout*stheta
        y = yout*ctheta-xout*stheta
        return x,y

    def getMass(self,sigCrit=1.):
        from math import pi
        from scipy.special import gammainc,gamma

        b,n,re = self.b,self.n,self.re
        k = 2.*n-1./3+4./(405.*n)+46/(25515.*n**2)
        amp = b*b*k**(2*n)/(2*n*re**2*(gamma(2*n)*gammainc(2*n,k*(b/re)**(1/n))))
        return 2*pi*sigCrit*re**2*amp*gamma(2*n)*n/k**(2*n)

    def getbFromMass(self,mass,sigCrit):
        from scipy import interpolate
        from scipy.special import gammainc,gamma
        import numpy

        n,re = self.n,self.re
        b = numpy.logspace(-3,1,401)*re
        k = 2.*n-1./3+4./(405.*n)+46/(25515.*n**2)
        amp = b*b*k**(2*n)/(2*n*re**2*(gamma(2*n)*gammainc(2*n,k*(b/re)**(1/n))))
        m = 2*pi*sigCrit*re**2*amp*gamma(2*n)*n/k**(2*n)
        model = interpolate.splrep(m,b)
        return interpolate.splev(mass,model)

    def setbFromMass(self,mass,sigCrit):
        self.b = self.getbFromMass(mass,sigCrit)



class sNFW(MassProfile):
    def __init__(self,b=None,rs=None,x=None,y=None):
        self.b = b
        self.rs = rs
        self.x = x
        self.y = y
        self.q = 1.
        self.pa = 0.
        self.theta = 0.

    def deflections(self,xin,yin):
        from numpy import arctanh,arctan,arctan2,log,sin,cos

        #x,y = self.align_coords(xin,yin)
        x = xin-self.x
        y = yin-self.y
        b,rs = self.b,self.rs
        X = b/rs
        if X<1.:
            amp = X**2/(8*arctanh(((1-X)/(1+X))**0.5)/(1-X**2)**0.5+4*log(X/2.))
        elif X==1:
            amp = 0.25/(1.+log(0.5))
        else:
            amp = X**2/(8*arctan(((X-1)/(1+X))**0.5)/(X**2-1)**0.5+4*log(X/2.))

        r2 = (x**2+y**2)/rs**2
        r = r2**0.5
        F = r*0.
        F[r<1.] = arctanh((1-r2[r<1.])**0.5)/(1-r2[r<1.])**0.5
        F[r==1.] = 1.
        F[r>1.] = arctan((r2[r>1.]-1.)**0.5)/(r2[r>1.]-1)**0.5

        dr = 4*amp*rs*(log(r/2)+F)/r
        A = arctan2(y,x)
        return dr*cos(A),dr*sin(A)

    def getbFromMass(self,mass,rvir,sigCrit):
        from numpy import arctanh,arctan,arctan2,log,sin,cos,pi,logspace
        from scipy import interpolate

        rs = self.rs
        c = rvir/self.rs

        x = logspace(-4,1,501)
        A = x*0.
        C = x<1.
        X = x[C].copy()
        A[C] = X**2/(8*arctanh(((1-X)/(1+X))**0.5)/(1-X**2)**0.5+4*log(X/2.))
        C = x==1.
        A[C] = 0.25/(1.+log(0.5))
        C = x>1.
        X = x[C].copy()
        A[C] = X**2/(8*arctan(((X-1)/(1+X))**0.5)/(X**2-1)**0.5+4*log(X/2.))

        m = 4*pi*A*rs**2*(log(1+c)-c/(1+c))*sigCrit
        model = interpolate.splrep(m,x*rs)
        return interpolate.splev(mass,model)

    def setbFromMass(self,mass,rvir,sigCrit):
        self.b = self.getbFromMass(mass,rvir,sigCrit)
