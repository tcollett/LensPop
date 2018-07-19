import models
from math import pi

_PowerLawPars = [['b','eta','pa','q','x','y'],['b','eta','q','theta','x','y']]
_ExtShearPars = [['b','pa','x','y'],['b','theta','x','y']]
_MassSheetPars = [['b','x','y'],['b','x','y']]
_NFWPars = [['b','x','y'],['b','x','y']]

class PowerLaw(models._PowerLaw):
    """
    A subclass for power-law mass models. The `power-law' aspect doesn't
        currently work, but it does work for SIE models.
    """
    def __init__(self,name,var=None,const=None):
        if const is None:
            const = {}
        if var is None:
            var = {}
        # Check for all keys to be set
        keys = var.keys()+const.keys()
        keys.sort()
        if keys not in _PowerLawPars:
            import sys
            print "Not all parameters defined!"
            sys.exit()
        models._PowerLaw.__init__(self)
        self.keys = keys
        self.values = {}
        self.vmap = {}
        for key in var.keys():
            self.values[key] = None
            self.vmap[key] = var[key]
        for key in const.keys():
            self.__setattr__(key,const[key])
        self.setPars()
        self.name = name
        self.NoFreeParams=False


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


class SIE(PowerLaw):
    def __init__(self,name,var=None,const=None):
        c = {}
        for key in const.keys():
            c[key] = const[key]
        c['eta'] = 1.
        PowerLaw.__init__(self,name,var,c)


class ExtShear(models._ExtShear):
    def __init__(self,name,var=None,const=None):
        if const is None:
            const = {}
        if var is None:
            var = {}
        # Check for all keys to be set
        keys = var.keys()+const.keys()
        keys.sort()
        if keys not in _ExtShearPars:
            import sys
            print "Not all parameters defined!",keys
            sys.exit()
        models._ExtShear.__init__(self)
        self.keys = keys
        self.values = {}
        self.vmap = {}
        for key in var.keys():
            self.values[key] = None
            self.vmap[key] = var[key]
        for key in const.keys():
            self.__setattr__(key,const[key])
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

class MassSheet(models._MassSheet):
    def __init__(self,name,var=None,const=None):
        if const is None:
            const = {}
        if var is None:
            var = {}
        # Check for all keys to be set
        keys = var.keys()+const.keys()
        keys.sort()
        if keys not in _MassSheetPars:
            import sys
            print "Not all parameters defined!",keys
            sys.exit()
        models._MassSheet.__init__(self)
        self.keys = keys
        self.values = {}
        self.vmap = {}
        for key in var.keys():
            self.values[key] = None
            self.vmap[key] = var[key]
        for key in const.keys():
            self.__setattr__(key,const[key])
        self.setPars()
        self.name = name

    def __setattr__(self,key,value):
        self.__dict__[key] = value

    def setPars(self):
        for key in self.vmap:
            self.__setattr__(key,self.vmap[key].value)




class NFW(models._NFW):
    """
    A subclass for power-law mass models. The `power-law' aspect doesn't
        currently work, but it does work for SIE models.
    """
    def __init__(self,name,var=None,const=None):
        if const is None:
            const = {}
        if var is None:
            var = {}
        # Check for all keys to be set
        keys = var.keys()+const.keys()
        keys.sort()
        if keys not in _PowerLawPars:
            import sys
            print "Not all parameters defined!"
            sys.exit()
        models._PowerLaw.__init__(self)
        self.keys = keys
        self.values = {}
        self.vmap = {}
        for key in var.keys():
            self.values[key] = None
            self.vmap[key] = var[key]
        for key in const.keys():
            self.__setattr__(key,const[key])
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


