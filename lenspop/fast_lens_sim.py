from __future__ import absolute_import, division, print_function

import indexTricks as iT
import numpy
from six.moves import cPickle as pickle
import pylab as plt

from pylens import *
import imageSim
import distances as D
from .surveys import Survey

from StochasticObserving import SO
from SignaltoNoise import S2N


class FastLensSim(SO,S2N):
    def __init__(self,surveyname,fractionofseeing=1):
        #-----------------------------------------------------
        ### Read in survey
        self.surveyname=surveyname
        survey=Survey(surveyname)#This stores typical survey  in Surveys.Survey
        self.survey=survey
        self.pixelsize=survey.pixelsize
        self.side=survey.side
        self.readnoise=survey.readnoise
        self.nexposures=survey.nexposures
        self.f_sky=survey.f_sky

        self.bands=survey.bands
        self.strategy=survey.strategy
        self.strategyx=survey.strategyx

        self.exposuretimes={}
        self.zeropoints={}
        self.stochasticobservingdata={}
        self.gains={}
        self.seeing={}
        self.psfscale={}
        self.psf={}
        self.psfFFT={}

        self.ET={}
        self.SB={}


        for i in range(len(survey.bands)):
            self.exposuretimes[survey.bands[i]]=survey.exposuretimes[i]
            self.zeropoints[survey.bands[i]]=survey.zeropoints[i]
            self.gains[survey.bands[i]]=survey.gains[i]
            self.stochasticobservingdata[survey.bands[i]]=survey.stochasticobservingdata[i]
        self.zeroexposuretime=survey.zeroexposuretime
        #-----------------------------------------------------
        ###do some setup
        self.xl=(self.side-1.)/2.
        self.yl=(self.side-1.)/2.
        self.x,self.y = iT.coords((self.side,self.side))
        self.r2 = (self.x-self.xl)**2+(self.y-self.yl)**2

        self.pixelunits=False

        #-----------------------------------------------------
        self.Reset()

    def Reset(self):
        self.sourcenumbers=[]
        #Some objects that need pre-defining as dictionaries
        self.magnification={}
        self.image={}
        self.sigma={}
        self.residual={}
        self.zeromagcounts={}
        self.xs={}
        self.ys={}
        self.ms={}
        self.qs={}
        self.ps={}
        self.rs={}
        self.ns={}
        self.bl={}
        self.src={}
        self.galmodel={}
        self.sourcemodel={}
        self.model={}
        self.totallensedsrcmag={}
        self.fakeLens={}
        self.image={}
        self.sigma={}
        self.fakeResidual={}
        self.fakeResidual[0]={}
        self.SN={}
        self.SNRF={}
        self.convolvedsrc={}
        self.convolvedgal={}

#===========================================================================

    def trytoconvert(self,par,p):
        try:return par/p
        except NameError:print("warning one of the parameters is not defined")

#===========================================================================

    def setLensPars(self,m,r,q,n=4,pixelunits=False,reset=True,xb=0,xp=0,jiggle=0):
        if reset: self.Reset()
        self.rl={}
        if pixelunits==False:
            for band in r.keys():
                self.rl[band]=self.trytoconvert(r[band],self.pixelsize)
        self.ml=m
        self.ql=q

        self.deltaxl=(numpy.random.rand()-0.5)*2*jiggle
        self.deltayl=(numpy.random.rand()-0.5)*2*jiggle
        if jiggle!=0:
            self.deltap=0.+(numpy.random.rand()-0.5)*180
            n=(numpy.random.rand()+1)*4
        else:
            self.deltap=0.

        self.nl=n
        self.gal = imageSim.Sersic('gal', {'x':self.xl+self.deltaxl,
                                           'y':self.yl+self.deltayl,
                                           'q':self.ql,
                                           'pa':90+self.deltap,
                                           're':self.rl[band],
                                           'n':self.nl})

        self.xb=xb
        self.xp=xp


#===========================================================================

    def setSourcePars(self,b,m,x,y,q,p,r,n=1,pixelunits=False,sourcenumber=1):
        if pixelunits==False:
            x=self.trytoconvert(x,self.pixelsize)
            y=self.trytoconvert(y,self.pixelsize)
            r=self.trytoconvert(r,self.pixelsize)
            b=self.trytoconvert(b,self.pixelsize)
        self.xs[sourcenumber]=x+self.xl+self.deltaxl+0.000001
        self.ys[sourcenumber]=y+self.yl+self.deltayl+0.000001
        self.ms[sourcenumber]=m
        self.qs[sourcenumber]=q
        self.ps[sourcenumber]=p
        self.rs[sourcenumber]=r
        self.ns[sourcenumber]=n
        self.bl[sourcenumber]=b
        self.src[sourcenumber] = imageSim.Sersic('src%i'%sourcenumber,
            {'x':self.xs[sourcenumber],'y':self.ys[sourcenumber],\
             'q':self.qs[sourcenumber],'pa':self.ps[sourcenumber],\
             're':self.rs[sourcenumber],'n':self.ns[sourcenumber]})
        self.sourcenumbers.append(sourcenumber)
        self.sourcemodel[sourcenumber]={}
        self.totallensedsrcmag[sourcenumber]={}
        self.fakeResidual[sourcenumber]={}
        self.SN[sourcenumber]={}
        self.SNRF[sourcenumber]={}
        self.convolvedsrc[sourcenumber]={}

#===========================================================================

    def lensASource(self,sourcenumber,bands):
        src=self.src[sourcenumber]
        lens=massmodel.PowerLaw('lens',{},{'x':self.xl+self.deltaxl,'y':self.yl+self.deltayl,'q':self.ql,'pa':90+self.deltap,'b':self.bl[sourcenumber],'eta':1})
        es=massmodel.ExtShear('lens',{},{'x':self.xl+self.deltaxl,'y':self.yl+self.deltayl,'pa':self.xp,'b':self.xb})
        lenses=[lens,es]

        a=51
        ox,oy=iT.coords((a,a))
        ps=(self.rs[sourcenumber]*(10./a))
        ox=(ox-(a-1)/2.)*ps+(self.xs[sourcenumber])
        oy=(oy-(a-1)/2.)*ps+(self.ys[sourcenumber])

        unlensedsrcmodel=(src.pixeval(ox,oy,csub=5)*(ps**2)).sum()
        srcnorm=unlensedsrcmodel.sum()
        unlensedsrcmodel/=srcnorm

        srcmodel=pylens.lens_images(lenses,src,[self.x,self.y],getPix=True,csub=5)[0]
        srcmodel[srcmodel<0]=0
        srcmodel/=srcnorm

        self.magnification[sourcenumber]=(numpy.sum(numpy.ravel(srcmodel))/numpy.sum(numpy.ravel(unlensedsrcmodel)))
        sm={}
        for band in bands:
            unlensedtotalsrcflux=10**(-(self.ms[sourcenumber][band]-self.zeropoints[band])/2.5)
            sm[band]=srcmodel*unlensedtotalsrcflux

            if sm[band].max()>0:
                self.totallensedsrcmag[sourcenumber][band]=-2.5*numpy.log10(sm[band].sum())+self.zeropoints[band]
            else:
                self.totallensedsrcmag[sourcenumber][band]=99
        return sm

#===========================================================================

    def EvaluateGalaxy(self,light,mag,bands):
        model={}
        lightm=light.pixeval(self.x,self.y,csub=5)
        lightm[lightm<0]=0
        lightm/=lightm.sum()
        for band in bands:
            flux=10**(-(mag[band]-self.zeropoints[band])/2.5)
            model[band]=lightm*flux

        return model

#===========================================================================

    def MakeModel(self, bands):
        #did you know that self.gal is actually fixed for all bands currently?
        self.galmodel=self.EvaluateGalaxy(self.gal,self.ml,bands)
        for sourcenumber in self.sourcenumbers:
            self.sourcemodel[sourcenumber]=self.lensASource(sourcenumber,bands)

        for band in bands:
            self.model[band]=self.galmodel[band]*1
            for sourcenumber in self.sourcenumbers:
                self.model[band]+=self.sourcemodel[sourcenumber][band]

#===========================================================================

    def ObserveLens(self,noisy=True,bands=[]):
      if bands==[]:bands=self.bands
      for band in bands:
        if self.seeing[band]!=0:
          convolvedgal, self.psfFFT[band] = \
            imageSim.convolve(self.galmodel[band], self.psf[band], True)
          convolvedgal[convolvedgal<0]=0
          self.convolvedgal[band]=convolvedgal


          convolvedmodel=convolvedgal*1

          convolvedsrc={}

          for sourcenumber in self.sourcenumbers:
            convolvedsrc[sourcenumber] = \
              imageSim.convolve(self.sourcemodel[sourcenumber][band],
                                self.psfFFT[band], False)[0]
            convolvedsrc[sourcenumber][convolvedsrc[sourcenumber]<0]=0
            self.convolvedsrc[sourcenumber][band]=convolvedsrc[sourcenumber]
            convolvedmodel+=convolvedsrc[sourcenumber]

          self.zeromagcounts[band]=(10**(-(0-self.zeropoints[band])/2.5))

          exposurecorrection=((self.ET[band]*1./self.zeroexposuretime))*self.gains[band]
          convolvedmodel*=exposurecorrection

          #skybackground per second per square arcsecond
          background=(10**(-(self.SB[band]-self.zeropoints[band])/2.5))*(self.pixelsize**2)
          tot_bg=background*exposurecorrection


          sigma=((convolvedmodel+tot_bg)+self.nexposures*(self.readnoise**0.5)**2)**.5

          fakeLens=convolvedmodel*1.
          if noisy:fakeLens+=(numpy.random.randn(self.side,self.side)*(sigma))

          #convert back to ADU/second:
          fakeLens/=exposurecorrection
          sigma/=exposurecorrection

          self.image[band]=fakeLens*1
          self.fakeLens[band]=fakeLens*1
          self.sigma[band]=sigma*1
          self.fakeResidual[0][band]=fakeLens-convolvedgal
          for sourcenumber in self.sourcenumbers:
              self.SN[sourcenumber][band]=self.SNfunc(\
                  convolvedsrc[sourcenumber],sigma)
              self.fakeResidual[sourcenumber][band]=\
                  fakeLens-convolvedmodel+convolvedsrc[sourcenumber]

#===========================================================================
    def loadModel(self,ideallens):
        self.galmodel,self.sourcemodel,self.model,self.magnification,self.totallensedsrcmag=ideallens
        self.image=self.model

#===========================================================================

    def loadConvolvedModel(self,ideallens):
        self.galmodel,self.sourcemodel,self.model,self.magnification,self.totallensedsrcmag=ideallens
        self.image=self.model


#===========================================================================

    def makeLens(self,stochastic=True,save=False,noisy=True,stochasticmode="MP",SOdraw=[],bands=[],musthaveallbands=False,MakeModel=True):
        if stochastic==True:self.stochasticObserving(mode=stochasticmode,SOdraw=SOdraw,musthaveallbands=musthaveallbands)
        if self.seeingtest=="Fail":return None
        if bands==[]:bands=self.bands

        if MakeModel:
            self.MakeModel(bands)

        if self.strategy=="resolve":
            if stochastic==True:self.stochasticObserving(mode=stochasticmode,SOdraw=SOdraw) #have to rerun stochastic observing now we know the magnification

        self.ObserveLens(noisy=noisy)
        return [self.galmodel,self.sourcemodel,self.model,self.magnification,self.totallensedsrcmag]


#===========================================================================


    def makeColorLens(self,bands=["g_SDSS","r_SDSS","i_SDSS"],recolourize=True):
        if self.surveyname=="Euclid" and bands==["g_SDSS","r_SDSS","i_SDSS"]:
            bands=["VIS","VIS","VIS"]
        import colorImage
        goodbands=[]
        for band in bands:
            try:
                self.image[band]
                goodbands.append(band)
            except KeyError:
                pass
        bands=goodbands
        if len(bands)==1:
            bands=[bands[0],bands[0],bands[0]]
        if len(bands)==2:
            bands=[bands[0],"dummy",bands[1]]
            self.ml["dummy"]=(self.ml[bands[0]]+self.ml[bands[2]])/2
            self.image["dummy"]=(self.image[bands[0]]+self.image[bands[2]])/2
        if recolourize:
            self.color = colorImage.ColorImage()
            self.color.bMinusr=(self.ml[bands[0]]-self.ml[bands[2]])/4.
            self.color.bMinusg=(self.ml[bands[0]]-self.ml[bands[1]])/4.
            self.color.nonlin=4.
            self.colorimage = self.color.createModel(\
                 self.image[bands[0]],self.image[bands[1]],self.image[bands[2]])
        else:
            self.colorimage = self.color.colorize(\
                 self.image[bands[0]],self.image[bands[1]],self.image[bands[2]])

        return self.colorimage


#===========================================================================

    def display(self,band="g_SDSS",bands=["g_SDSS","r_SDSS","i_SDSS"]):
        if self.surveyname=="Euclid":bands=["VIS","VIS","VIS"]
        import pylab as plt
        plt.ion()
        plt.figure(1)
        plt.imshow(self.makeColorLens(bands=bands),interpolation="none")
        plt.figure(2)
        import colorImage
        self.color = colorImage.ColorImage()#sigma-clipped single band residual
        plt.imshow(self.color.createModel(self.fakeResidual[0][band],self.fakeResidual[0][band],self.fakeResidual[0][band])[:,:,0],interpolation="none")
        plt.figure(3)
        plt.imshow(self.fakeResidual[0][band],interpolation="none")
        try:
            self.fakeResidual[1]["RF"]
            plt.figure(4)
            plt.imshow(self.fakeResidual[1]["RF"],interpolation="none")
        except KeyError: pass

        plt.draw()
        raw_input()
        plt.ioff()

#===========================================================================

    def Rank(self,mode,band="g_SDSS",bands=["g_SDSS","r_SDSS","i_SDSS"]):
        import pylab as plt
        plt.ion()
        rank="d"
        while rank not in ["0","1","2","3","4","-1","-2","-3"]:
            if mode=="colour":
                plt.imshow(self.makeColorLens(bands=bands),interpolation="none")
                plt.draw()
            if mode=="rf":
                plt.imshow(self.fakeResidual[0]["RF"],interpolation="none")
                plt.draw()
            if mode=="best":
                plt.imshow(self.fakeResidual[0][band],interpolation="none")
                plt.draw()
            rank=raw_input()
            if rank=="":rank="0"
        plt.ioff()
        return rank
#===========================================================================
