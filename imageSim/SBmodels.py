from __future__ import absolute_import, division, print_function
from math import pi

from .profiles import _Sersic, _Gauss
from .pointSource import PixelizedModel as PM, GaussianModel as GM

def cnts2mag(cnts,zp):
    from math import log10
    return -2.5*log10(cnts) + zp

_SersicPars = [['amp','n','pa','q','re','x','y'],
                ['logamp','n','pa','q','re','x','y'],
                ['amp','n','q','re','theta','x','y'],
                ['logamp','n','q','re','theta','x','y']]

class SBModel:
    def __init__(self,name,pars,convolve=0):
        self.keys = pars.keys()
        self.keys.sort()
        if self.keys not in self._SBkeys:
            import sys
            print('Not all parameters were defined!')
            sys.exit()
        self._baseProfile.__init__(self)
        self.vmap = {}
        self.pars = pars
        for key in self.keys:
            try:
                v = self.pars[key].value
                self.vmap[key] = self.pars[key]
            except:
                self.__setattr__(key,self.pars[key])
        self.setPars()
        self.name = name
        self.convolve = convolve


    def __setattr__(self,key,value):
        if key=='pa':
            self.__dict__['pa'] = value
            if value is not None:
                self.__dict__['theta'] = value*pi/180.
        elif key=='theta':
            if value is not None:
                self.__dict__['pa'] = value*180./pi
            self.__dict__['theta'] = value
        elif key=='logamp':
            if value is not None:
                self.__dict__['amp'] = 10**value
        else:
            self.__dict__[key] = value


    def setPars(self):
        for key in self.vmap:
            self.__setattr__(key,self.vmap[key].value)


class Sersic(SBModel, _Sersic):
    _baseProfile = _Sersic
    _SBkeys = [['amp','n','pa','q','re','x','y'],
                ['logamp','n','pa','q','re','x','y'],
                ['amp','n','q','re','theta','x','y'],
                ['logamp','n','q','re','theta','x','y']]

    def __init__(self,name,pars,convolve=0):
        SBModel.__init__(name,pars,convolve)

    def getMag(self,amp,zp):
        from scipy.special import gamma
        from math import exp,pi
        n = self.n
        re = self.re
        k = 2.*n-1./3+4./(405.*n)+46/(25515.*n**2)
        cnts = (re**2)*amp*exp(k)*n*(k**(-2*n))*gamma(2*n)*2*pi
        return cnts2mag(cnts,zp)

    def Mag(self,zp):
        return self.getMag(self.amp,zp)


class Gauss(_Gauss):

    def __init__(self,name,var=None,const=None,convolve=0):
        if const is None:
            const = {}
        if var is None:
            var = {}
        # Check for all keys to be set
        keys = var.keys()+const.keys()
        keys.sort()
        if 'r0' not in keys:
            keys.append('r0')
            keys.sort()
        if keys!=['amp','pa','q','r0','sigma','x','y']:
            import sys
            print("Not all parameters defined!")
            sys.exit()
        imageSim._Gauss.__init__(self)
        self.invar = var
        self.keys = keys
        self.values = {'r0':None}
        self.vmap = {}
        for key in var.keys():
            self.values[key] = None
            self.vmap[key] = var[key]
        for key in const.keys():
            self.__setattr__(key,const[key])
        self.setPars()
        self.name = name
        self.convolve = convolve

    def __setattr__(self,key,value):
        if key=='pa':
            self.__dict__['pa'] = value
            if value is not None:
                self.__dict__['theta'] = value*pi/180.
        elif key=='theta':
            if value is not None:
                self.__dict__['pa'] = value*180./pi
            self.__dict__['theta'] = value
        elif key=='logamp':
            if value is not None:
                self.__dict__['amp'] = 10**value
        elif key=='scale':
            self.__dict__['sigma'] = value
        else:
            self.__dict__[key] = value

    def setPars(self):
        for key in self.vmap:
            self.__setattr__(key,self.vmap[key].value)

    def getMag(self,amp,zp):
        from math import exp,pi
        if self.r0 is None:
            cnts = amp/(2*pi*self.sigma**2)
        else:
            from scipy.special import erf
            r0 = self.r0
            s = self.sigma
            r2pi = (2*pi)**0.5
            cnts = amp*pi*s*(r2pi*r0*(1.+erf(r0/(s*2**0.5)))+2*s*exp(-0.5*r0**2/s**2))
        return cnts2mag(cnts,zp)

    def Mag(self,zp):
        return self.getMag(self.amp,zp)



class PointSource(GM, PM):
    def __init__(self,name,model,var=None,const=None):
        if const is None:
            const = {}
        if var is None:
            var = {}
        keys = var.keys()+const.keys()
        keys.sort()
        if keys!=['amp','x','y']:
            print("Not all parameters defined!",keys)
            df
        self.keys = keys
        self.values = {}
        self.vmap = {}
        self.ispix = False
        for key in var.keys():
            self.values[key] = None
            self.vmap[var[key]] = key
        for key in const.keys():
            self.values[key] = const[key]
        if type(model)==type([]):
            GM.__init__(self,model)
        else:
            PM.__init__(self,model)
            self.ispix = True
        self.setValues()
        self.name = name
        self.convolve = None

    def __setattr__(self,key,value):
        if key=='logamp':
            if value is not None:
                self.__dict__['amp'] = 10**value
        else:
            self.__dict__[key] = value

    def pixeval(self,xc,yc,dummy1=None,dummy2=None,**kwargs):
        if self.ispix==True:
            return PM.pixeval(self,xc,yc)
        else:
            return GM.pixeval(self,xc,yc)

    def setValues(self):
        self.x = self.values['x']
        self.y = self.values['y']
        if 'amp' in self.keys:
            self.amp = self.values['amp']
        elif self.values['logamp'] is not None:
            self.amp = 10**self.values['logamp']

    def getMag(self,amp,zp):
        return cnts2mag(amp,zp)

    def Mag(self,zp):
        return self.getMag(self.amp,zp)

    def setPars(self,pars):
        for key in self.vmap:
            self.values[self.vmap[key]] = pars[key]
        self.setValues()
