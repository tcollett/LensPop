"""
A module to compute cosmological distances, including:
    comoving_distance (Dc)
    angular_diameter_distance (Da)
    luminosity_distance (Dl)
    comoving_volume (volume)

Assumes a default cosmology OMEGA_M = 0.27, OMEGA_L = 0.73, but these can be
    set in the functions or globally.

All distances have units of h**-1 Mpc.
"""
c = 299792458.

class Distance:
    def __init__(self,cosmo=[0.27,0.73,1.]):
        self.OMEGA_M = cosmo[0]
        self.OMEGA_L = cosmo[1]
        self.h = cosmo[2]
        self.Dc = self.comoving_distance
        self.Dt = self.comoving_transverse_distance
        self.Dm = self.comoving_transverse_distance
        self.Da = self.angular_diameter_distance
        self.Dl = self.luminosity_distance
        self.volume = self.comoving_volume

    def reset(self):
        self.OMEGA_M = 0.27
        self.OMEGA_L = 0.73
        self.h = 1.


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
        f = lambda z,m,l,k : (m*(1.+z)**3+k*(1.+z)**2+l)**-0.5
        om = self.OMEGA_M
        ol = self.OMEGA_L
        ok = 1.-om-ol
        return (c/self.h)*integrate.romberg(f,z1,z2,(om,ol,ok))/1e5

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
        f = lambda z,m,l: (self.comoving_distance(0.,z)**2)/((m*(1.+z)**3+l)**0.5)
        om = self.OMEGA_M
        ol = self.OMEGA_L
        return (c/self.h)*integrate.romberg(f,z1,z2,(om,ol))/1e5
