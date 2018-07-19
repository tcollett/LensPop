from __init__ import *
import cPickle
import pyfits
import sys,os
import pylab as plt
import glob

params = {
   'axes.labelsize': 16,
   'text.fontsize': 16,
   'legend.fontsize': 14,
   'xtick.labelsize': 14,
   'ytick.labelsize': 14,
   'text.usetex': False,
    'figure.figsize': [6, 4]
   }
plt.rcParams.update(params)

surveystoread=[]
surveystoread+=["Euclid"]
surveystoread+=["CFHT"]
surveystoread+=["CFHTa"]
surveystoread+=["DESc"]
surveystoread+=["DESb"]
surveystoread+=["DESa"]
surveystoread+=["LSSTc"]
surveystoread+=["LSSTb"]
surveystoread+=["LSSTa"]


weights={}
bl  ={} 
zs  ={}
zl  ={}
sigl={}
ql  ={}
rs  ={}
ms  ={}
mag ={}
for survey in surveystoread:
 for sourcepop in ["lsst"]:
   ss=str(survey)
   filename="../%s_%s_lists.pkl"%(survey,sourcepop)
   filename="%s_%s_lists.pkl"%(survey,sourcepop)
   try:
       f=open(filename,"rb")
       namelist=cPickle.load(f)
       f.close()
       weights[ss]=namelist[0]
       if survey=="CFHTa":
           for key in weights[ss].keys():
               weights[ss][key]=list(numpy.array(weights[ss][key])/100.)
       bl[ss]=namelist[1]
       zs[ss]=namelist[2]
       rs[ss]=namelist[3]
       ms[ss]=namelist[4]
       zl[ss]=namelist[5]
       sigl[ss]=namelist[6]
       ql[ss]=namelist[7]
       mag[ss]=namelist[8]
   except IOError:
       continue
for survey in surveystoread:
    print survey
    for key in ["resolved","rfpf"]:
        print key, numpy.array(weights[survey][key]).sum()


save=True#False
show=True

bson=numpy.array([2.66,1.24,1.27,2.39,1.41,1.27,1.00,1.3,1.0,1.19,1.22,1.36,1.76,1.19,1.29,1.56,1.04,0.85,1.10,1.23,1.16,0.93,1.03,1.4,0.74,1.21,1.14,1.74,2.03,1.23,2.55,1.05,1.51,4.36,0.94,0.93,3.11,1.79,0.96,1.40,1.3,0.81,1.95,1.66,1.55,1.07,1.06,1.38,0.52,2.16,1.40,1.44])
plt.hist(bson,bins=numpy.linspace(0,3,16),weights=bson*0+220./len(bson),fc="red",alpha=0.8,label="SL2S")
a,b=numpy.histogram(bl["CFHT"]["rfpf"],bins=numpy.linspace(0,3,31),weights=weights["CFHT"]["rfpf"])
a*=2#double since finer bins
plt.plot(b[:-1]+(b[1]-b[0])/2.,a,c="k",lw=3,ls="dashed",label="Prediction")
plt.xlabel(r"$\Theta_\mathrm{E}$ (arcsec)")
plt.ylabel(r"Lenses per $\Theta_\mathrm{E}$ bin")
plt.legend()
plt.tight_layout()
if save:plt.savefig("/home/ttemp/papers/LensPop/thetaE.pdf")
if save:plt.savefig("/home/ttemp/papers/LensPop/thetaE.png")
if show:plt.show()
plt.cla()




f=open('lsst.1sqdegree_catalog2.pkl','rb')
data=cPickle.load(f)
f.close()
glsstcat=data[:,3]
wl=numpy.ones(len(glsstcat))*numpy.sum(weights["LSSTc"]["resolved"])/float(len(glsstcat))
a,b=numpy.histogram(glsstcat,bins=numpy.linspace(20,30,21),weights=wl)
#a*=2
plt.hist(ms["LSSTc"]["resolved"],bins=numpy.linspace(20,30,21),weights=weights["LSSTc"]["resolved"],fc="grey",alpha=0.6,label="LSST-optimal")
plt.hist(ms["LSSTa"]["resolved"],bins=numpy.linspace(20,30,21),weights=weights["LSSTa"]["resolved"],fc="red",alpha=1,label="LSST-all")

plt.plot(b[:-1]+(b[1]-b[0])/2.,a,c="k",lw=3,ls="dashed",label="Catalogue")
plt.xlabel(r"Unlensed source $g$-band magnitude")
plt.ylabel(r"Lenses per bin")
plt.legend(loc=2)
plt.tight_layout()
if save:plt.savefig("/home/ttemp/papers/LensPop/LSSTg.pdf")
if show:plt.show()
plt.cla()

wc=numpy.array(weights["LSSTc"]["resolved"])
wa=numpy.array(weights["LSSTa"]["resolved"])
wb=numpy.array(weights["LSSTb"]["resolved"])
maskc=numpy.where((numpy.array(bl["LSSTc"]["resolved"]))>0)
maska=numpy.where((numpy.array(bl["LSSTa"]["resolved"]))>1)
maskb=numpy.where((numpy.array(bl["LSSTb"]["resolved"]))<1)
print numpy.sum(wc[maskc])- numpy.sum(wa[maska])- numpy.sum(wb[maskb])
plt.hist(bl["LSSTc"]["resolved"],bins=numpy.linspace(0,3,31),weights=weights["LSSTc"]["resolved"],fc="grey",alpha=1,label="LSST-optimal")
plt.hist(bl["LSSTa"]["resolved"],bins=numpy.linspace(0,3,31),weights=weights["LSSTa"]["resolved"],fc="red",alpha=0.6,label="LSST-all")
plt.hist(bl["LSSTb"]["resolved"],bins=numpy.linspace(0,3,31),weights=weights["LSSTb"]["resolved"],fc="blue",alpha=0.6,label="LSST-best")


a,b=numpy.histogram(bl["Euclid"]["resolved"],bins=numpy.linspace(0,3,31),weights=weights["Euclid"]["resolved"])
#a*=2
plt.plot(b[:-1]+(b[1]-b[0])/2.,a,c="k",lw=3,ls="dashed",label="Euclid")
plt.xlabel(r"$\Theta_\mathrm{E}$ (arcsec)")
plt.ylabel(r"Lenses per bin")
plt.legend()
plt.tight_layout()
if save:plt.savefig("/home/ttemp/papers/LensPop/LSSTb.pdf")
if show:plt.show()
plt.cla()


#a,b=numpy.histogram(mag["Euclid"]["resolved"],bins=numpy.linspace(0,100,31),weights=weights["Euclid"]["resolved"])
#a*=2#double for finer bins
#a=numpy.log10(a)
#plt.plot(b[:-1]+(b[1]-b[0])/2.,a,c="k",lw=3,ls="dashed")
#plt.xlabel(r"$\mu$")
#plt.ylabel(r"Lenses per mag bin")
#plt.tight_layout()
#plt.savefig("/home/ttemp/papers/LensPop/mags.pdf")
#plt.show()
#plt.cla()






plt.hist(bl["DESc"]["resolved"],bins=numpy.linspace(0,3,31),weights=weights["DESc"]["resolved"],fc="grey",alpha=1,label="DES-optimal")
plt.hist(bl["DESa"]["resolved"],bins=numpy.linspace(0,3,31),weights=weights["DESa"]["resolved"],fc="red",alpha=0.6,label="DES-all")
plt.hist(bl["DESb"]["resolved"],bins=numpy.linspace(0,3,31),weights=weights["DESb"]["resolved"],fc="blue",alpha=0.6,label="DES-best")
plt.xlabel(r"$\Theta_\mathrm{E}$ (arcsec)")
plt.ylabel(r"Lenses per bin")
plt.legend()
plt.tight_layout()
if save:plt.savefig("/home/ttemp/papers/LensPop/DESb.pdf")
if show:plt.show()
plt.cla()


plt.hist(zl["DESc"]["resolved"],bins=numpy.linspace(0,5.5,56),weights=weights["DESc"]["resolved"],fc="grey",alpha=1)
plt.hist(zs["DESc"]["resolved"],bins=numpy.linspace(0,5.5,56),weights=weights["DESc"]["resolved"],fc="red",alpha=0.6)
plt.xlabel(r"redshift")
plt.ylabel(r"Lenses perbin")
plt.tight_layout()
if save:plt.savefig("/home/ttemp/papers/LensPop/DESz.pdf")
if show:plt.show()
plt.cla()


plt.hist(zl["LSSTc"]["resolved"],bins=numpy.linspace(0,5.5,56),weights=weights["LSSTc"]["resolved"],fc="grey",alpha=1)
plt.hist(zs["LSSTc"]["resolved"],bins=numpy.linspace(0,5.5,56),weights=weights["LSSTc"]["resolved"],fc="red",alpha=0.6)
plt.xlabel(r"redshift")
plt.ylabel(r"Lenses perbin")
plt.tight_layout()
if save:plt.savefig("/home/ttemp/papers/LensPop/LSSTz.pdf")
if show:plt.show()
plt.cla()


plt.hist(zl["LSSTc"]["resolved"],bins=numpy.linspace(0,5.5,56),normed=True,fc="grey",alpha=1)
plt.hist(zs["LSSTc"]["resolved"],bins=numpy.linspace(0,5.5,56),normed=True,fc="red",alpha=0.6)
plt.xlabel(r"redshift")
plt.ylabel(r"Lenses perbin")
plt.tight_layout()
if save:plt.savefig("/home/ttemp/papers/LensPop/LSSTz.pdf")
#plt.show()
plt.cla()


plt.hist(zl["LSSTc"]["resolved"],bins=numpy.linspace(0,5.5,56),normed=True,fc="grey",alpha=1)
plt.hist(zs["LSSTc"]["resolved"],bins=numpy.linspace(0,5.5,56),normed=True,fc="red",alpha=0.6)
plt.xlabel(r"redshift")
plt.ylabel(r"Lenses per bin")
plt.tight_layout()
if save:plt.savefig("/home/ttemp/papers/LensPop/LSSTz.pdf")
#plt.show()
plt.cla()

a,bins=numpy.histogram(sigl["Euclid"]["resolved"],bins=numpy.linspace(100,400,81),normed=True)
b,bins1=numpy.histogram(sigl["CFHTa"]["resolved"],bins=numpy.linspace(100,400,81),normed=True)
c,bins2=numpy.histogram(sigl["LSSTc"]["resolved"],bins=numpy.linspace(100,400,81),normed=True)
#c*=2
d,bins3=numpy.histogram(sigl["DESc"]["resolved"],bins=numpy.linspace(100,400,81),normed=True)
from scipy.ndimage.filters import gaussian_filter as gf
idx=2
a[idx:]=gf(a,2)[idx:]
b=gf(b,4)
c=gf(c,4)
d=gf(d,4)

#d*=4
db=bins[1]-bins[0]
db1=bins1[1]-bins1[0]
db2=bins2[1]-bins2[0]
db3=bins3[1]-bins3[0]

plt.plot(bins[:-1]+db/2.,a,lw=3,c='k',label="Euclid")
plt.plot(bins2[:-1]+db2/2.,c,lw=3,c='b',label="LSST")
plt.plot(bins3[:-1]+db3/2.,d,lw=3,c='r',label="DES")
plt.plot(bins1[:-1]+db1/2.,b,lw=3,c='g',label="CFHT")
plt.legend()
plt.xlabel("$\sigma_{\mathrm{V}}$ (km s$^{-1}$)")
plt.ylabel("P($\sigma_{\mathrm{V}}$)")

#plt.ylabel("$\sigma_{\mathrm{V}}$ bin")
plt.tight_layout()
if save:plt.savefig("/home/ttemp/papers/LensPop/sigmadist.pdf")
if show:plt.show()
plt.cla()

azbins=numpy.linspace(0,5,81)
zbins=numpy.linspace(0,2,81)

a,bins=numpy.histogram(zl["Euclid"]["resolved"],bins=zbins,normed=True)
b,bins1=numpy.histogram(zl["CFHTa"]["resolved"],bins=zbins,normed=True)
c,bins2=numpy.histogram(zl["LSSTc"]["resolved"],bins=zbins,normed=True)
d,bins3=numpy.histogram(zl["DESc"]["resolved"],bins=zbins,normed=True)
aa,abins=numpy.histogram(zs["Euclid"]["resolved"],bins=azbins,normed=True)
ab,abins1=numpy.histogram(zs["CFHT"]["resolved"],bins=azbins,normed=True)
ac,abins2=numpy.histogram(zs["LSSTc"]["resolved"],bins=azbins,normed=True)
ad,abins3=numpy.histogram(zs["DESc"]["resolved"],bins=azbins,normed=True)


"""
from scipy.ndimage import gaussian_filter1d
import numpy as np
t = np.linspace(0, 1, len(bins[:-1]+db/2.))
t2 = np.linspace(0, 1, 100)
x2 = np.interp(t2, t, bins[:-1]+db/2.)
y2 = np.interp(t2, t, a)
sigma = 1
x3 = gaussian_filter1d(x2, sigma)
y3 = gaussian_filter1d(y2, sigma)
"""


from scipy.ndimage.filters import gaussian_filter as gf
idx=10
a[idx:]=gf(a,2)[idx:]
b[idx:]=gf(b,4)[idx:]
c[idx:]=gf(c,2)[idx:]
d[idx:]=gf(d,4)[idx:]
idx2=2
aa[idx2:]=gf(aa,2)[idx2:]
ab[idx2:]=gf(ab,4)[idx2:]
ac[idx2:]=gf(ac,2)[idx2:]
ad[idx2:]=gf(ad,4)[idx2:]

for item in [a,b,c,d]:
    item[:idx]=0
bins[:idx]=0


#x4 = np.interp(t, t2, x3)
#y4 = np.interp(t, t2, y3)



db=bins[1]-bins[0]
db1=bins1[1]-bins1[0]
db2=bins2[1]-bins2[0]
db3=bins3[1]-bins3[0]
adb=abins[1]-abins[0]
adb1=abins1[1]-abins1[0]
adb2=abins2[1]-abins2[0]
adb3=abins3[1]-abins3[0]
plt.plot(bins[:-1]+db/2.,a,lw=3,c='k',label="Euclid")
#plt.plot(bins[:-1]+db/2.,a,lw=3,c='k')
#plt.plot(x3,y3,lw=3,c='orange')

plt.plot(bins2[:-1]+db2/2.,c,lw=3,c='b',label="LSST")
plt.plot(bins3[:-1]+db3/2.,d,lw=3,c='r',label="DES")
plt.plot(bins1[:-1]+db1/2.,b,lw=3,c='g',label="CFHT")
plt.plot(abins[:-1]+adb/2.,aa,lw=3,c='k',ls="dotted")
plt.plot(abins2[:-1]+adb2/2.,ac,lw=3,c='b',ls="dotted")
plt.plot(abins3[:-1]+adb3/2.,ad,lw=3,c='r',ls="dotted")
plt.plot(abins1[:-1]+adb1/2.,ab,lw=3,c='g',ls="dotted")
plt.ylabel("P($z$)")
plt.xlabel("redshift")
#plt.ylabel("$\sigma_{\mathrm{V}}$ bin")
plt.legend()
plt.tight_layout()
if save:plt.savefig("/home/ttemp/papers/LensPop/zdist.pdf")
if show:plt.show()
plt.cla()
