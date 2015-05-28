import glob,numpy,os

files = glob.glob('*UKIRT.res')
files.sort()
for f in files:
    wave,flux = numpy.loadtxt(f).T
    wave *= 10.
    A = wave.argsort()
    wave = wave[A]
    flux = flux[A]
    os.system('cp %s %s'%(f,f.replace('res','orig')))
    f = open(f,'w')
    for i in range(wave.size):
        f.write('%7.1f  %11.9f\n'%(wave[i],flux[i]))
    f.close()
