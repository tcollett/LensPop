import numpy,copy

class SO():
    def __init__(self):#This class can only be inherited from
        pass

    def drawPSFandSB(self,band):
        dat=self.stochasticobservingdata[band]
        k=numpy.random.randint(len(dat[:,0]))
        return dat[k,0],dat[k,1]       

    def CalculateETSB(self,sbs,band):
        et=self.exposuretimes[band]*(len(sbs)*(1./self.nexposures))
        sbf=10**(-(sbs)/2.5)
        sbf=sbf.mean()
        sb=-2.5*numpy.log10(sbf)
        return et,sb

    def PSFfloor(self,a=[],dummy=None):
        #function that encodes the seeing strategy
        mode=self.strategy
        x=self.strategyx
        if mode == "absolute":
            if x==0: 
                return 10
            else: return x
        if mode == "percentile":
            import scipy.stats
            return stats.scoreatpercentile(a,x)
        if mode == "best":
            a=numpy.sort(a)
            return a[int(x)]
        if mode == "resolve":
            #return (self.fos*self.bl[1])*self.pixelsize)
            try:
                if self.bfac*(self.bl[1])**2-self.rfac*(self.rs[1])**2<0:
                    floor1=0
                else:
                    floor1=(self.bfac*(self.bl[1])**2-self.rfac*(self.rs[1])**2)**0.5
            except FloatingPointError:
                floor1=0 # this should never get called

            try:
                floor2=self.rs[1]*self.magnification[1]
            except KeyError:
                floor2=999


            floor=numpy.min([floor1,floor2])

            return (floor)*self.pixelsize

        if mode == "resolveclever":
            print "warning: the resolveclever code isn't finsihed"
            (self.fos*self.bl[1]*self.pixelsize)
            numpy.sort(a)
            la=len(a)
            b=a[a<(self.fos*self.bl[1]*self.pixelsize)]
            lb=len(b)

            #definitely include seeings less than the source size:
            defoin=b[b>self.rs[1]]

            return a[int(x)-1]

            return (self.fos*self.bl[1]*self.pixelsize)


    def stochasticObserving(self,mode="MP",seeingstrategy="absolute",stratfloor=10,SOdraw=[],musthaveallbands=False):
        psfs={}
        sbs={}
        psfs2={}
        sbs2={}
        worstacceptedpsfband={}
        worstacceptedpsf=0
        for band in self.bands:
            worstacceptedpsfband[band]=0
            if SOdraw==[]:
                psfs[band]=numpy.zeros(self.nexposures)
                sbs[band]=numpy.zeros(self.nexposures)
                psfs2[band]=numpy.zeros(self.nexposures)
                sbs2[band]=numpy.zeros(self.nexposures)
                for i in range(self.nexposures):
                    a,b=self.drawPSFandSB(band)
                    psfs[band][i]=a
                    sbs[band][i]=b
                    psfs2[band][i]=a*1
                    sbs2[band][i]=b*1
                self.SOdraw=[psfs2,sbs2]
                psffloor=self.PSFfloor(psfs[band],dummy=band)
            else:
                psfs[band]=SOdraw[0][band]
                sbs[band]=SOdraw[1][band]
                psffloor=self.PSFfloor(psfs[band],dummy=band)

            for i in range(len(psfs[band])):
                if psfs[band][i]<psffloor and psfs[band][i]>worstacceptedpsf:
                    worstacceptedpsf=psfs[band][i]
                if psfs[band][i]<psffloor and psfs[band][i]>worstacceptedpsfband[band]:
                    worstacceptedpsfband[band]=psfs[band][i]

            sbs[band]=sbs[band][psfs[band]<psffloor]
            psfs[band]=psfs[band][psfs[band]<psffloor]
            if len(psfs[band][psfs[band]<psffloor])==0:
                sbs[band]=numpy.array([0])
                psfs[band]=numpy.array([0.01])

        #since images need to have same psf, we use the worst accepted:
        #for band in self.bands:###
        #    print band,worstacceptedpsfband[band],",",###
        #print worstacceptedpsfband###
        for band in self.bands:
            if mode=="1P":
                self.seeing[band]=worstacceptedpsf
            if mode=="MP":
                self.seeing[band]=worstacceptedpsfband[band]
            if self.seeing[band]==0:continue

            self.psfscale[band]=self.seeing[band]/2.355
            self.psf[band]= numpy.exp(-0.5*self.r2/(self.psfscale[band]/self.pixelsize)**2)
            self.psf[band]/=numpy.sum(self.psf[band])
            self.psfFFT[band]=None

        #print self.seeing###

        #now calculate new exposuretimes and skybrightnesses
        for band in self.bands:
            self.ET[band],self.SB[band]=self.CalculateETSB(sbs[band],band)


        self.seeingtest="Fail"
        for band in self.bands:
            if self.SeeingTest(1,band) ==True:
                self.seeingtest="Pass"
        if musthaveallbands:self.seeingtest="Pass"


    def SeeingTest(self,src,band):
        if self.bfac*(self.bl[src])**2<(self.rfac*(self.rs[src]))**2+(self.seeing[band]/self.pixelsize)**2:
            return False
        else:
            return True
