from __init__ import *
import cPickle
import pyfits
import sys
import pylab as plt
import time
sigfloor=200

L=LensSample(reset=False,sigfloor=sigfloor,cosmo=[0.3,0.7,0.7])

experiment="Euclid"
frac=0.1

a=20#SN threshold
b=3#Magnification threshold

c=1000
d=1000


#experiment="DES"
if len(sys.argv)>1:
    experiment=sys.argv[1]
    frac=float(sys.argv[2])
if len(sys.argv)>3:
    a=int(sys.argv[3])
    b=int(sys.argv[4])
    #c=int(sys.argv[5])
    #d=int(sys.argv[6])

firstod=1
nsources=1


surveys=[]

if experiment=="Euclid":
    surveys+=["Euclid"]
if experiment=="CFHT":
    surveys+=["CFHT"] #full coadd (Gaussianised)
if experiment=="CFHTa":
    surveys+=["CFHTa"] #dummy CFHT
if experiment=="DES":
    surveys+=["DESc"] #Optimal stacking of data
    surveys+=["DESb"] #Best Single epoch image
    surveys+=["DESa"] #full coadd (Gaussianised)
if experiment=="LSST":
    surveys+=["LSSTc"] #Optimal stacking of data
    surveys+=["LSSTb"] #Best Single epoch image
    surveys+=["LSSTa"] #full coadd (Gaussianised)
    #print "only doing LSSTc"


S={}
n={}
for survey in surveys:
    S[survey]=FastLensSim(survey,fractionofseeing=1)
    S[survey].bfac=float(2)
    S[survey].rfac=float(2)


t0=time.clock()

#for sourcepop in ["lsst","cosmos"]:
for sourcepop in ["lsst"]:
  chunk=0
  Si=0
  SSPL={}
  foundcount={}
  for survey in surveys:
      foundcount[survey]=0

  if sourcepop=="cosmos":
      nall=1100000
  elif sourcepop=="lsst":
      nall=12530000
  nall=int(nall*frac)

  for i in range(nall):
    if i%10000==0:
        print "about to load"
        L.LoadLensPop(i,sourcepop)
        print i,nall

    if i!=0:
        if i%10000==0 or i==100 or i==300 or i==1000 or i==3000:
            t1=time.clock()
            ti=(t1-t0)/float(i)
            tl=(nall-i)*ti
            tl/=60#mins
            hl=numpy.floor(tl/(60))
            ml=tl-(hl*60)
            print i,"%ih%im left"%(hl,ml)

    lenspars=L.lens[i]
    if lenspars["lens?"]==False:
        del L.lens[i]
        continue

    lenspars["rl"]["VIS"]=(lenspars["rl"]["r_SDSS"]+\
                           lenspars["rl"]["i_SDSS"]+lenspars["rl"]["z_SDSS"])/3
    for mi in [lenspars["ml"],lenspars["ms"][1]]:
        mi["VIS"]=(mi["r_SDSS"]+mi["i_SDSS"]+mi["z_SDSS"])/3

    


    #if lenspars["zl"]>1 or lenspars["zl"]<0.2 or lenspars["ml"]["i_SDSS"]<17 or lenspars["ml"]["i_SDSS"]>22:continue# this is a CFHT compare quick n dirty test

    lenspars["mag"]={}
    lenspars["msrc"]={}
    lenspars["mag"]={}
    lenspars["msrc"]={}
    lenspars["SN"]={}
    lenspars["bestband"]={}
    lenspars["pf"]={}
    lenspars["resolved"]={}
    lenspars["poptag"]={}
    lenspars["seeing"]={}
    lenspars["rfpf"]={}
    lenspars["rfsn"]={}

    lastsurvey="non"
    for survey in surveys:

        S[survey].setLensPars(lenspars["ml"],lenspars["rl"],lenspars["ql"],reset=True)
        for j in range(nsources):
            S[survey].setSourcePars(lenspars["b"][j+1],lenspars["ms"][j+1],\
                                    lenspars["xs"][j+1],lenspars["ys"][j+1],\
                                    lenspars["qs"][j+1],lenspars["ps"][j+1],\
                                    lenspars["rs"][j+1],sourcenumber=j+1    )

        if survey[:3]+str(i)!=lastsurvey:
            model=S[survey].makeLens(stochasticmode="MP")
            SOdraw=numpy.array(S[survey].SOdraw)
            if type(model)!=type(None):
                lastsurvey=survey[:3]+str(i)
            if S[survey].seeingtest=="Fail":
                lenspars["pf"][survey]={}
                lenspars["rfpf"][survey]={}
                for src in S[survey].sourcenumbers:
                    lenspars["pf"][survey][src]=False
                    lenspars["rfpf"][survey][src]=False
                continue#try next survey
        else: 
            S[survey].loadModel(model)
            S[survey].stochasticObserving(mode="MP",SOdraw=SOdraw)
            if S[survey].seeingtest=="Fail":
                lenspars["pf"][survey]={}
                for src in S[survey].sourcenumbers:
                    lenspars["pf"][survey][src]=False
                continue#try next survey
            S[survey].ObserveLens()

        mag,msrc,SN,bestband,pf=S[survey].SourceMetaData(SNcutA=a,magcut=b,SNcutB=[c,d])
        lenspars["SN"][survey]={}
        lenspars["bestband"][survey]={}
        lenspars["pf"][survey]={}
        lenspars["resolved"][survey]={}
        lenspars["poptag"][survey]=i
        lenspars["seeing"][survey]=S[survey].seeing
        rfpf={}
        rfsn={}
        for src in S[survey].sourcenumbers:
            rfpf[src]=False      
            rfsn[src]=[0]
            lenspars["mag"][src]=mag[src]
            lenspars["msrc"][src]=msrc[src]
            lenspars["SN"][survey][src]=SN[src]
            lenspars["bestband"][survey][src]=bestband[src]
            lenspars["pf"][survey][src]=pf[src]
            lenspars["resolved"][survey][src]=S[survey].resolved[src]
        if survey!="Euclid":
            if S[survey].seeingtest!="Fail":
                if survey not in ["CFHT","CFHTa"]:
                    S[survey].makeLens(noisy=True,stochasticmode="1P",SOdraw=SOdraw,MakeModel=False)
                    rfpf,rfsn=S[survey].RingFinderSN(SNcutA=a,magcut=b,SNcutB=[c,d],mode="donotcrossconvolve")
                else:
                    rfpf,rfsn=S[survey].RingFinderSN(SNcutA=a,magcut=b,SNcutB=[c,d],mode="crossconvolve")
        lenspars["rfpf"][survey]=rfpf
        lenspars["rfsn"][survey]=rfsn

        ###
        #This is where you can add your own lens finder
        #e.g.
        #found=Myfinder(S[survey].image,S[survey].sigma,\
        #                    S[survey].psf,S[survey].psfFFT)
        #NB/ image,sigma, psf, psfFFT are dictionaries 
        #    The keywords are the filters, e.g. "g_SDSS", "VIS" etc

        #then save any outputs you'll need to the lenspars dictionary:
        #lenspars["my_finder_result"]=found

        ###

        #If you want to save the images (it may well be a lot of data!):
        #folder="where_to_save_fits_images"
        #folder="%s/%i"%(folder,i)
        #for band in S[survey].bands:
            #img=S[survey].image[band]
            #sig=S[survey].sigma[band]
            #psf=S[survey].psf[band]
            #resid=S[survey].fakeResidual[0][band]#The lens subtracted

        #resid contains the lensed source, with the lens subtracted
        #assuming the subtraction is poisson noise limited (i.e. ideal)

            #pyfits.PrimaryHDU(img).writeto("%s/image_%s.fits"%(folder,band),\
                #                               clobber=True)
            #pyfits.PrimaryHDU(sig).writeto("%s/sigma_%s.fits"%(folder,band),\
                #                               clobber=True)
            #pyfits.PrimaryHDU(psf).writeto("%s/psf_%s.fits"%(folder,band),\
                #                               clobber=True)
            #pyfits.PrimaryHDU(resid).writeto("%s/galsub_%s.fits"%(folder,band),clobber=True)

        ###

        L.lens[i]=None #delete used data for memory saving
            
    accept=False
    for survey in surveys:
        if lenspars["pf"][survey][1]:
            accept=True

    if accept:
        #S[survey].display(band="VIS",bands=["VIS","VIS","VIS"])
        #if Si>100:exit()
        Si+=1
        SSPL[Si]=lenspars.copy() 
        if (Si+1)%1000==0:
            f=open("LensStats/%s_%s_Lens_stats_%i.pkl"%(experiment,sourcepop,chunk),"wb")
            cPickle.dump([frac,SSPL],f,2)
            f.close()
            SSPL={} # reset SSPL or memory fills up
            chunk+=1

    del L.lens[i]

  f=open("LensStats/%s_%s_Lens_stats_%i.pkl"%(experiment,sourcepop,chunk),"wb")
  cPickle.dump([frac,SSPL],f,2)
  f.close()
  print Si

bl=[]
for j in SSPL.keys():
    try: 
        if SSPL[j]["rfpf"][survey][1]:
            bl.append(SSPL[j]["b"][1])
    except KeyError:pass
