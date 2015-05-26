import numpy

class S2N():
    def __init__(self):#This class can only be inherited from
        pass

    def imageRegions(image,sig,sigfloor=0.5):
        image[image/sig<significancefloor]=0
        masks, multiplicity = ndimage.measurements.label(image)
        labels=numpy.arange(1, multiplicity+1)


    def SNfunc(self,data,sig,significancefloor=0.5):
        D=data.ravel()
        S=sig.ravel()

        args=numpy.argsort(-D/S)
        D=numpy.take(D,args)
        S=numpy.take(S,args)
        Dsum=numpy.cumsum(D)
        Ssum=numpy.cumsum(S**2)**0.5
        SN=(Dsum/Ssum).max()

        #regional SN
        import scipy.ndimage as  ndimage
        data[data/sig<significancefloor]=0
        masks, multiplicity = ndimage.measurements.label(data)
        labels=numpy.arange(1, multiplicity+1)
        SNs=numpy.zeros(multiplicity+1)
        SNs[0]=SN
        for i in range(multiplicity):
            D=data[masks==i+1].ravel()
            S=sig[masks==i+1].ravel()
            args=numpy.argsort(-D/S)
            D=numpy.take(D,args)
            S=numpy.take(S,args)
            Dsum=numpy.cumsum(D)
            Ssum=numpy.cumsum(S**2)**0.5
            SNi=(Dsum/Ssum).max()
            SNs[i+1]=SNi
        SNs=-numpy.sort(-SNs)
        return SNs

    def SourceMetaData(self,SNcutA=15,magcut=3,SNcutB=[10,8]):
        self.mag={}
        self.msrc={}
        self.bestband={}
        self.passfail={}
        self.resolved={}

        for src in self.sourcenumbers:
            self.resolved[src]={}
            SNr={}
            self.mag[src]=self.magnification[src]
            self.msrc[src]={}
            for band in self.bands:
                if self.seeing[band]!=0:
                    self.msrc[src][band]=self.totallensedsrcmag[src][band]
                    if self.mag[src]*self.rs[src]>(self.seeing[band]/self.pixelsize):
                        self.resolved[src][band]=True
                        SNindex=0
                    else:
                        self.resolved[src][band]=False
                        SNindex=2
                    try:
                        SNr[band]=self.SN[src][band][SNindex]
                    except IndexError: SNr[band]=0
                    except KeyError: SNr[band]=0

                else:
                    self.SN[src][band]=[0,0,0]
                    self.msrc[src][band]=[99]
                    self.resolved[src][band]=False
                    SNr[band]=0


            self.bestband[src],dummy = max(SNr.iteritems(), key=lambda x:x[1])
            
            self.passfail[src]=False
            try:
                if self.SN[src][self.bestband[src]][2]>min(SNcutB) and \
                   self.SN[src][self.bestband[src]][1]>max(SNcutB):
                    self.passfail[src]=True
                    ltype=1
            except IndexError: 
                pass
                

            try:
                #print self.SN[src][self.bestband[src]][0],SNcutA,self.mag[src],magcut,self.resolved[src][self.bestband[src]]
                if self.SN[src][self.bestband[src]][0]>SNcutA and \
                   (self.mag[src]>magcut) and self.resolved[src][self.bestband[src]]:
                    self.passfail[src]=True
                    ltype=2
            except IndexError: 
                pass

            if  self.SeeingTest(src,self.bestband[src]) ==False:
                self.passfail[src]=False

        #debugger
        #for src in self.sourcenumbers:
        #    print src,self.passfail[src],self.SeeingTest(src,self.bestband[src]),self.SN[src][self.bestband[src]][0],self.resolved[src][self.bestband[src]],self.mag[src],ltype

        return self.mag,self.msrc,self.SN,self.bestband,self.passfail

#===========================================================================
    def RingFinderSN(self,bands=["g_SDSS","i_SDSS"],repair=True,mode="crossconvolve",SNcutA=15,magcut=3,SNcutB=[10,8],runringfinder=False,mustbeseen=False):
        self.rfpf={}
        for src in self.sourcenumbers:
            self.SNRF[src]=[0]
            self.rfpf[src]=False
        try: 
            if self.seeing[bands[0]]==0:
                return self.rfpf,self.SNRF
        except KeyError:
            return self.rfpf,self.SNRF
        try:
            if self.seeing[bands[1]]==0:
                return self.rfpf,self.SNRF
        except KeyError:     
            return self.rfpf,self.SNRF
        if mode=="crossconvolve":
            seeing=(self.seeing[bands[1]]**2+self.seeing[bands[0]]**2)**.5
            for band in bands:
                self.psfFFT[band]=None
                self.psfscale[band]=seeing/2.355
                self.psf[band]= numpy.exp(-0.5*self.r2/(self.psfscale[band]/self.pixelsize)**2)
                self.psf[band]/=numpy.sum(self.psf[band])
            self.ObserveLens(bands=bands)
        else: seeing=self.seeing[bands[0]]

        self.seeing["RF"]=seeing

        seen=False
        for src in self.sourcenumbers:
            if self.SeeingTest(src,"RF"):
                seen=True
        if mustbeseen:
            seen=True

        if seen==False:
            return [self.rfpf,self.SNRF]

        assert (self.psf[bands[0]]-self.psf[bands[1]]).sum()==0, "psf missmatch - can't run ringfinder"

        B=self.image[bands[0]]
        R=self.image[bands[1]]
        sB=self.sigma[bands[0]]
        sR=self.sigma[bands[1]]

        r=self.r2**0.5
        r*=self.pixelsize
        mask=((r<2.7) & (r>0.5))

        alpha=B[mask].sum()*1./R[mask].sum()

        self.D=B-alpha*R
        self.S=(sB**2+(alpha*sR)**2)**.5
        self.fakeResidual[0]["RF"]=self.D
        for src in self.sourcenumbers:
            self.SNRF[src]=self.SNfunc(self.convolvedsrc[src]["g_SDSS"]-alpha*self.convolvedsrc[src]["i_SDSS"],self.S)
            d=self.convolvedsrc[src]["g_SDSS"]-alpha*self.convolvedsrc[src]["i_SDSS"]
            d+=(numpy.random.randn(self.side,self.side)*(self.S))
            self.fakeResidual[src]["RF"]=d
            if self.mag[src]*self.rs[src]>(seeing/self.pixelsize):
                self.resolved[src]["RF"]=True
            else:
                self.resolved[src]["RF"]=False

            self.rfpf[src]=False
            try:
                if self.SNRF[src][2]>min(SNcutB) and \
                   self.SNRF[src][1]>max(SNcutB):
                    self.rfpf[src]=True
            except IndexError: pass
            try:
                if self.SNRF[src][0]>SNcutA \
                   and self.mag[src]>magcut \
                   and self.mag[src]*self.rs[src]>(seeing/self.pixelsize):
                    self.rfpf[src]=True
            except IndexError: pass
            if self.SeeingTest(src,"RF")==False:
                self.rfpf[src]=False       
                self.passfail[src]=False

            if runringfinder:
                import RingFinder
                RF=RingFinder.RingFinder(B,R,sB,sR,self.pixelsize,
                                         self.zeromagcounts["g_SDSS"],
                                         self.zeromagcounts["i_SDSS"])
                RFo=RF.ringfind()
                self.D=RF.D*1
                return RFo,self.rfpf,self.SNRF
        return self.rfpf,self.SNRF
