"""
A module to compute cosmological distances, including:
    comoving_distance (Dc)
    angular_diameter_distance (Da)
    luminosity_distance (Dl)
    comoving_volume (volume)

"""
c = 299792458.
G = 4.3e-6
from math import pi
import warnings
warnings.warn("Default cosmology is Om=0.3,Ol=0.7,h=0.7,w=-1 and distance units are Mpc!",ImportWarning)

class Distance():
    def __init__(self,cosmo=[0.3,0.7,0.7]):
        self.OMEGA_M = cosmo[0]
        self.OMEGA_L = cosmo[1]
        self.h = cosmo[2]
        self.w = -1.
        self.wpars = None
        self.w_analytic = False
        self.Dc = self.comoving_distance
        self.Dt = self.comoving_transverse_distance
        self.Dm = self.comoving_transverse_distance
        self.Da = self.angular_diameter_distance
        self.Dl = self.luminosity_distance
        self.dm = self.distance_modulus
        self.volume = self.comoving_volume

    def set(self,cosmo):
        self.OMEGA_M = cosmo[0]
        self.OMEGA_L = cosmo[1]
        self.h = cosmo[2]

    def reset(self):
        self.OMEGA_M = 0.3
        self.OMEGA_L = 0.7
        self.h = 0.7
        self.w = -1.

    def age(self,z):
        from scipy import integrate
        f = lambda zp,m,l,k : (m/zp+k+l*zp**2)**-0.5
        om = self.OMEGA_M
        ol = self.OMEGA_L
        ok = 1.-om-ol
        return (9.778/self.h)*integrate.romberg(f,1e-300,1/(1.+z),(om,ol,ok))

    def comoving_distance(self,z1,z2=0.):
        from scipy import integrate
        if z2<z1:
            z1,z2 = z2,z1
        def fa(z):
            if self.w_analytic==True:
                return self.w(z,self.wpars)
            from math import exp
            wa = lambda z : (1.+self.w(z,self.wpars))/(1.+z)
            #return exp(3.*integrate.romberg(wa,0,z))
            return exp(3.*integrate.quad(wa,0,z)[0])
        if type(self.w)==type(self.comoving_distance) or type(self.w)==type(fa):
            f = lambda z,m,l,k : (m*(1.+z)**3+k*(1.+z)**2+l*fa(z))**-0.5
        elif self.w!=-1.:
            f = lambda z,m,l,k : (m*(1.+z)**3+k*(1.+z)**2+l*(1.+z)**(3.*(1.+self.w)))**-0.5
        else:
            f = lambda z,m,l,k : (m*(1.+z)**3+k*(1.+z)**2+l)**-0.5
        om = self.OMEGA_M
        ol = self.OMEGA_L
        ok = 1.-om-ol
#        return (c/self.h)*integrate.romberg(f,z1,z2,(om,ol,ok))/1e5
        return (c/self.h)*integrate.quad(f,z1,z2,(om,ol,ok))[0]/1e5

    def comoving_transverse_distance(self,z1,z2=0.):
        dc = 1e5*self.comoving_distance(z1,z2)/(c/self.h)
        ok = 1.-self.OMEGA_M-self.OMEGA_L
        if ok>0:
            from math import sinh,sqrt
            dtc = sinh(sqrt(ok)*dc)/sqrt(ok)
        elif ok<0:
            from math import sin,sqrt
            ok *= -1.
            dtc = sin(sqrt(ok)*dc)/sqrt(ok)
        else:
            dtc = dc
        return (c/self.h)*dtc/1e5

    def angular_diameter_distance(self,z1,z2=0.):
        if z2<z1:
            z1,z2 = z2,z1
        return self.comoving_transverse_distance(z1,z2)/(1.+z2)

    def luminosity_distance(self,z):
        return (1.+z)*self.comoving_transverse_distance(z)

    def comoving_volume(self,z1,z2=0.):
        from scipy import integrate
        if z2<z1:
            z1,z2 = z2,z1
        f = lambda z,m,l,k: (self.comoving_distance(0.,z)**2)/((m*(1.+z)**3+k*(1.+z)**2+l)**0.5)
        om = self.OMEGA_M
        ol = self.OMEGA_L
        ok = 1.-om-ol
        return 4*pi*(c/self.h)*integrate.romberg(f,z1,z2,(om,ol,ok))/1e5

    def rho_crit(self,z):
        H2 = (self.OMEGA_M*(1+z)**3 + self.OMEGA_L)*(self.h/10.)**2
        return 3*H2/(8.*pi*G)

    def distance_modulus(self,z):
        from math import log10
        if z>0:return 5*log10(self.luminosity_distance(z)*1e5)
        else: return 0
