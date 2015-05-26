from __init__ import *
import cPickle
import pyfits
import sys,os
import pylab as plt
import glob

params = {
   'axes.labelsize': 14,
   'text.fontsize': 14,
   'legend.fontsize': 10,
   'xtick.labelsize': 10,
   'ytick.labelsize': 10,
   'text.usetex': False,
    'figure.figsize': [6, 4]
   }
plt.rcParams.update(params)

sourcepops=["lsst"]

experiment="Euclid"
#experiment="CFHT"
#experiment="LSST"
#experiment="DES"

if len(sys.argv)>1:
    experiment=sys.argv[1]

surveystoread=[]
if experiment=="Euclid":
    surveystoread+=["Euclid"]
elif experiment=="CFHT":
    surveystoread+=["CFHT"]
elif experiment=="CFHTa":
    surveystoread+=["CFHTa"]
elif experiment=="DES":
    surveystoread+=["DESc"]
    surveystoread+=["DESb"]
    surveystoread+=["DESa"]
elif experiment=="LSST":
    surveystoread+=["LSSTc"]
    surveystoread+=["LSSTb"]
    surveystoread+=["LSSTa"]
else:
    surveystoread=[str(experiment)]
    experiment=experiment[:-1]

    
for survey in surveystoread:
  for sourcepop in sourcepops:
    if survey[-2]=="a":
        surveyname=survey[:-1]+"_full_coadd"
    elif survey[-2]=="b":
        surveyname=survey[:-1]+"_best_epoch"
    elif survey[-2]=="c":
        surveyname=survey[:-1]+"_optimal_coadd"
    else:
        surveyname=survey
    filename="%s_%s_lists.pkl"%(survey,sourcepop)
    lensparsfile="lenses_%s.txt"%survey
    f=open(lensparsfile,"w")
    print 
    #os.system("rm %s"%filename) #this line resets the read-in
    bl={}
    zs={}
    zl={}
    sigl={}
    ql={}
    rs={}
    ms={}
    mag={}
    weights={}
    for key in ["resolved","rfpf"]:
        bl[key]=[]
        zs[key]=[]
        rs[key]=[]
        ms[key]=[]
        zl[key]=[]
        sigl[key]=[]
        ql[key]=[]
        mag[key]=[]
        rs[key]=[]
        weights[key]=[]
       
    if experiment=="CFHT":
      frac=42000.*1./150.
      bands=["g_SDSS","r_SDSS","i_SDSS"]

    if experiment=="CFHTa":
      frac=42000.*1./150.
      bands=["g_SDSS","r_SDSS","i_SDSS"]

    elif experiment=="Euclid":
      frac=42000.*1./15000.
      bands=["VIS"]
                
    elif experiment=="DES":
      frac=42000.*1./5000.
      bands=["g_SDSS","r_SDSS","i_SDSS"]

    elif experiment=="LSST":
      frac=42000.*1./20000.
      bands=["g_SDSS","r_SDSS","i_SDSS"]


    filelist=glob.glob("LensStats/%s_%s_Lens_stats_*.pkl"%(experiment,sourcepop))

    chunki=0
    ilist=[]
    print survey
    for chunk in filelist:
        print chunki
        chunki+=1
        f2=open(chunk,"rb")
        fracsky,sspl=cPickle.load(f2)
        fract=frac*fracsky
        f2.close()
        I=0
        for i in sspl.keys():
            if i in ilist:  
                continue 
            else:
                try:
                    sspl[i]["seeing"][survey]
                except KeyError:
                    continue
                f.write("%.2f "%sspl[i]["zl"])
                f.write("%.2f "%sspl[i]["zs"][1])
                f.write("%.2f "%sspl[i]["b"][1])
                f.write("%.2f "%sspl[i]["sigl"])
                f.write("%.2f "%sspl[i]["ql"])
                f.write("%.2f "%sspl[i]["rl"]["g_SDSS"])
                for band in bands:
                    f.write("%.2f "%sspl[i]["ml"][band])
                f.write("%.2f "%sspl[i]["rl"]["g_SDSS"])
                f.write("%.2f "%sspl[i]["xs"][1])
                f.write("%.2f "%sspl[i]["ys"][1])
                f.write("%.2f "%sspl[i]["qs"][1])
                f.write("%.2f "%sspl[i]["ps"][1])
                f.write("%.2f "%sspl[i]["rs"][1])
                f.write("%.2f "%sspl[i]["mag"][1])
                for band in bands:
                    f.write("%.2f "%sspl[i]["seeing"][survey][band])
                    f.write("%.2f "%sspl[i]["SN"][survey][1][band][0])
                if survey!="Euclid":
                    f.write("%.2f "%sspl[i]["rfsn"][survey][1][0])
                f.write("\n")


                ilist.append(str(i))
                if sspl[i]["pf"][survey][1]==False:continue

                try:
                  if sspl[i]["resolved"][survey][1][sspl[i]["bestband"][survey][1]]:
                    bb=sspl[i]["bestband"][survey][1]
                    if sspl[i]["mag"][1]<3:continue
                    if sspl[i]["SN"][survey][1][bb][0]<20:continue
                    if sspl[i]["rs"][1]*sspl[i]["mag"][1]<sspl[i]["seeing"][survey][band]:continue
                    if sspl[i]["rs"][1]**2+sspl[i]["seeing"][survey][bb]**2\
                       < 2* sspl[i]["b"][1]**2 \
                       :continue


                    bl["resolved"].append(sspl[i]["b"][1])
                    weights["resolved"].append(1./frac)
                    zs["resolved"].append(sspl[i]["zs"][1])
                    rs["resolved"].append(sspl[i]["rs"][1])
                    zl["resolved"].append(sspl[i]["zl"])
                    sigl["resolved"].append(sspl[i]["sigl"])
                    ql["resolved"].append(sspl[i]["ql"])
                    mag["resolved"].append(sspl[i]["mag"][1])
                    ms["resolved"].append(sspl[i]["ms"][1]["g_SDSS"])

                    if sspl[i]["rfpf"][survey][1]:
                        if sspl[i]["rfsn"][survey][1][0]<20:continue


                        if experiment=="CFHT" or experiment=="CFHTa":
                            if sspl[i]["resolved"][survey][1]["RF"]==False:continue
                            if sspl[i]["zl"]>1:continue
                            if sspl[i]["zl"]<0.2:continue
                            if sspl[i]["ml"]["i_SDSS"]<17:continue
                            if sspl[i]["ml"]["i_SDSS"]>22:continue

                        bl["rfpf"].append(sspl[i]["b"][1])
                        weights["rfpf"].append(1./fract)
                        zs["rfpf"].append(sspl[i]["zs"][1])
                        rs["rfpf"].append(sspl[i]["rs"][1])
                        zl["rfpf"].append(sspl[i]["zl"])
                        sigl["rfpf"].append(sspl[i]["sigl"])
                        ql["rfpf"].append(sspl[i]["ql"])
                        mag["rfpf"].append(sspl[i]["mag"][1])
                        ms["rfpf"].append(sspl[i]["ms"][1]["g_SDSS"])


                except KeyError:
                  pass
    f.close()

    if survey[-2]=="a":
        surveyname=survey[:-1]+" (full coadd)"
    elif survey[-2]=="b":
        surveyname=survey[:-1]+" (best single epoch imaging)"
    elif survey[-2]=="c":
        surveyname=survey[:-1]+" (optimal coadd)"
    else:
        surveyname=survey

    print survey, "will find",
    print numpy.sum(numpy.array(weights["resolved"]).ravel()),
    print "lenses assuming poisson limited galaxy subtraction in all bands, or",
    print numpy.sum(numpy.array(weights["rfpf"]).ravel()), 
    print "lenses in the g-i difference images"

    f=open(filename,"wb")
    cPickle.dump([weights,bl,zs,rs,ms,zl,sigl,ql,mag],f,2)
    f.close()



bson=numpy.array([2.66,1.24,1.27,2.39,1.41,1.27,1.00,1.3,1.0,1.19,1.22,1.36,1.76,1.19,1.29,1.56,1.04,0.85,1.10,1.23,1.16,0.93,1.03,1.4,0.74,1.21,1.14,1.74,2.03,1.23,2.55,1.05,1.51,4.36,0.94,0.93,3.11,1.79,0.96,1.40,1.3,0.81,1.95,1.66,1.55,1.07,1.06,1.38,0.52,2.16,1.40,1.44])
plt.hist(bson,bins=numpy.linspace(0,3,16),weights=bson*0+220./len(bson),fc="grey",alpha=0.6)
a,b=numpy.histogram(bl["rfpf"],bins=numpy.linspace(0,3,31),weights=weights["rfpf"])
a*=2#double for finer bins
plt.plot(b[:-1]+(b[1]-b[0])/2.,a,c="k",lw=3,ls="dashed")
plt.xlabel(r"$\Theta_\mathrm{E}$ (arcsec)")
plt.ylabel(r"Lenses per $\Theta_\mathrm{E}$ bin")
plt.tight_layout()
plt.show()
