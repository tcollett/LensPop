from ndinterp import ndInterp
from math import pi

class _MassModel:
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


    def get_deflection_angles(self,x,y,**args):
        """
        Powerlaw deflection angles following Barkana
        """
        from scipy.integrate import quad
        from math import cos,sin
        import numpy
        shape = x.shape
        x = x.ravel()
        y = y.ravel()
        q = self.q
        x,y = self.align_coords(x,y)
        x2 = x**2
        y2 = y**2
        rho0 = (x2+y2/q**2)**0.5
        ell = (1.-q**2)**0.5

        rho = numpy.linspace(1e-11,1.,self.nsteps)
        rho = rho.repeat(x.size).reshape((self.nsteps,x.size))
        rho *= rho0

        rhoell = (rho*ell)**2
        d = ((rhoell+y2-x2)**2 + 4*x2*y2)**0.5
        w = ((d+x2+y2+rhoell)/(d+x2+y2-rhoell))**0.5

        kappa = self.kappa(rho,**args)
        integrand = rho*kappa*w/(x2+y2*w**4)
        xmap = 2.*x*q*integrand.sum(0)*rho0/self.nsteps
        ymap = 2.*y*q*(integrand*w**2).sum(0)*rho0/self.nsteps
        xmap[rho0==0] = 0.
        ymap[rho0==0] = 0.
        theta = self.theta
        ctheta = cos(theta)
        stheta = sin(theta)
        x = xmap*ctheta-ymap*stheta
        y = ymap*ctheta+xmap*stheta
        return x,y

        return self.align_coords(xmap,ymap,True)


    def setGrid(self):
        from scipy import interpolate
        x0 = numpy.logspace(-3,1,81)
        etas = numpy.linspace(0.,2.,21)
        qs = numpy.linspace(0.2,1.,17)
        grid1 = numpy.empty((x0.size,x0.size,etas.size,qs.size))
        grid2 = numpy.empty(grid1.shape)
        for i in range(qs.size):
            q = qs[i]
            q2 = q**2
            b = 1-q2
            for j in range(etas.size):
                eta = etas[j]
                g = 0.5*eta-1. # g = -1*gamma
                for k in range(x0.size):
                    x = x0[k]
                    for l in range(x0.size):
                        y = x0[l]
                        qb = ((2*x*y)/b)**g  # q_bar
                        qt = q*(x*y)**0.5/b  # q_tilde
                        sb = 0.5*(x/y - y/x) + s**2*b/(2*x*y)
                        nu1 = s**2*b/(2*x*y)
                        nu2 = nu1+ 0.5*b*(x/y + y/(x*q2))
                        nu = numpy.logspace(nu1,nu2,1001)
                        mu = nu-sb
                        t = (1+mu**2)**0.5
                        f1 = (t-mu)**0.5/t
                        f2 = (t+mu)**0.5/t
                        ng = nu**g
                        I1 = interpolate.splrep(nu,f1*ng)
                        I2 = interpolate.splrep(nu,f2*ng)
                        grid1[k,l,i,j] = qt*interpolate.splint(nu1,nu2,I1)
                        grid2[k,l,i,j] = qt*interpolate.splint(nu1,nu2,I2)
                pylab.imshow(grid1[:,:,i,j])
                pylab.show()

    def alphax(self,x,y):
        pass

    def alphay(self,x,y):
        pass


class _PowerLaw(_MassModel):
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

    def create_grid(self,x,y,q,eta):
        import numpy
        xgrid = numpy.empty((x.size,y.size,q.size,eta.size))
        ygrid = xgrid.copy()

        shape = (x.size,y.size)
        xnew = x.repeat(y.size)
        ynew = y.repeat(x.size).reshape((y.size,x.size)).T.flatten()

        for i in range(q.size):
            self.q = q[i]
            for j in range(eta.size):
                self.eta = eta[j]
                a,b = self.deflection_angles(xnew,ynew)
                xgrid[:,:,i,j] = a.reshape(shape)
                ygrid[:,:,i,j] = b.reshape(shape)
        ax = [x,y,q,eta]
        alphax = ndInterp()
        alphax.ndI_setup(ax,xgrid)
        alphay = ndInterp()
        alphay.ndI_setup(ax,ygrid)

        self.xmodel = alphax
        self.ymodel = alphay

    def deflections(self,xin,yin):
        if self.NoFreeParams==True:
            try:
                return self.deflx,self.defly
            except: pass


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
        self.deflx=x
        self.defly=y
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



class _ExtShear(_MassModel):

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


class _MassSheet(_MassModel):

    def __init__(self,x=None,y=None,b=0.):
        self.x = x
        self.y = y
        self.b = b

    def deflections(self,x,y):
        from math import sin,cos,pi
        x = x-self.x
        y = y-self.y

        #\vec(alpha)=-kappa*\vec(x)

        alpha_x = -self.b*x
        alpha_y = -self.b*y

        return alpha_x,alpha_y





class _NFW(_MassModel):
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

    def create_grid(self,x,y,q,eta):
        import numpy
        xgrid = numpy.empty((x.size,y.size,q.size,eta.size))
        ygrid = xgrid.copy()

        shape = (x.size,y.size)
        xnew = x.repeat(y.size)
        ynew = y.repeat(x.size).reshape((y.size,x.size)).T.flatten()

        for i in range(q.size):
            self.q = q[i]
            for j in range(eta.size):
                self.eta = eta[j]
                a,b = self.deflection_angles(xnew,ynew)
                xgrid[:,:,i,j] = a.reshape(shape)
                ygrid[:,:,i,j] = b.reshape(shape)
        ax = [x,y,q,eta]
        alphax = ndInterp()
        alphax.ndI_setup(ax,xgrid)
        alphay = ndInterp()
        alphay.ndI_setup(ax,ygrid)

        self.xmodel = alphax
        self.ymodel = alphay

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
