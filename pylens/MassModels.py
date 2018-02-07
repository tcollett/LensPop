from __future__ import absolute_import, division, print_function

from .MassProfiles import PowerLaw as MP_PowerLaw, \
                          ExtShear as MP_ExtShear, \
                          PointSource as MP_PointSource, \
                          Sersic as MP_Sersic, \
                          SersicG as MP_SersicG, \
                          sNFW as MP_sNFW
from math import pi

class MassModel:
    def __init__(self,name,pars):
        self.keys = pars.keys()
        self.keys.sort()
        if self.keys not in self._MMkeys:
            import sys
            print('Not all (or too many) parameters were defined!')
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

    def __setattr__(self,key,value):
        if key=='pa':
            self.__dict__['pa'] = value
            if value is not None:
                self.__dict__['theta'] = value*pi/180.
        elif key=='theta':
            if value is not None:
                self.__dict__['pa'] = value*180./pi
            self.__dict__['theta'] = value
        else:
            self.__dict__[key] = value

    def setPars(self):
        for key in self.vmap:
            self.__setattr__(key,self.vmap[key].value)

        
class PowerLaw(MassModel,MP_PowerLaw):
    _baseProfile = MP_PowerLaw
    _MMkeys = [['b','eta','pa','q','x','y'],['b','eta','q','theta','x','y']]

    def __init__(self,name,pars):
        MassModel.__init__(self,name,pars)


class SIE(PowerLaw):
    def __init__(self,name,pars):
        pars['eta'] = 1.
        PowerLaw.__init__(self,name,pars)


class ExtShear(MassModel,MP_ExtShear):
    _baseProfile = MP_ExtShear
    _MMkeys = [['b','pa','x','y'],['b','theta','x','y']]

    def __init__(self,name,pars):
        MassModel.__init__(self,name,pars)


class PointSource(MassModel,MP_PointSource):
    _baseProfile = MP_PointSource
    _MMkeys = [['b','x','y'],['b','x','y']]

    def __init__(self,name,pars):
        MassModel.__init__(self,name,pars)


class Sersic(MassModel,MP_Sersic):
    _baseProfile = MP_Sersic
    _MMkeys = [['b','n','pa','q','reff','x','y'],
                ['b','n','pa','q','re','x','y'],
                ['b','n','q','reff','theta','x','y'],
                ['b','n','q','re','theta','x','y']]

    def __init__(self,name,pars):
        MassModel.__init__(self,name,pars)

    def __setattr__(self,key,value):
        if key=='pa':
            self.__dict__['pa'] = value
            if value is not None:
                self.__dict__['theta'] = value*pi/180.
        elif key=='theta':
            if value is not None:
                self.__dict__['pa'] = value*180./pi
            self.__dict__['theta'] = value
        elif key=='reff':
            if value is not None:
                self.__dict__['re'] = value
        else:
            self.__dict__[key] = value


class SersicG(MassModel,MP_SersicG):
    _baseProfile = MP_SersicG
    _MMkeys = [['b','n','pa','q','reff','x','y'],
                ['b','n','pa','q','re','x','y'],
                ['b','n','q','reff','theta','x','y'],
                ['b','n','q','re','theta','x','y']]

    def __init__(self,name,pars):
        MassModel.__init__(self,name,pars)

    def __setattr__(self,key,value):
        if key=='pa':
            self.__dict__['pa'] = value
            if value is not None:
                self.__dict__['theta'] = value*pi/180.
        elif key=='theta':
            if value is not None:
                self.__dict__['pa'] = value*180./pi
            self.__dict__['theta'] = value
        elif key=='reff':
            if value is not None:
                self.__dict__['re'] = value
        else:
            self.__dict__[key] = value


class sNFW(MassModel,MP_sNFW):
    _baseProfile = MP_sNFW
    _MMkeys = [['b','rs','x','y']]

    def __init__(self,name,pars):
        MassModel.__init__(self,name,pars)

    def __setattr__(self,key,value):
        self.__dict__[key] = value
