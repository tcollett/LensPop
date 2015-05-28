import numpy as N
import os

filters = ['u','g','r','i','z']
for filter in filters:
    #os.system("wget -nc http://www.sdss.org/dr7/instruments/imager/filters/%s.dat"%filter)
    data = N.loadtxt('%s.dat'%filter)
    w = data[:,0].copy()
    t = data[:,1].copy()
    a = data[:,4].copy()

    f = t*a
    o = open('%s_SDSS.res'%filter,'w')
    for i in range(f.size):
        o.write("%4.0f  %6.4f\n"%(w[i],f[i]))
    o.close()
    
