import scipy.ndimage as  ndimage
import pylab as plt
import moments
import indexTricks as iT
import numpy

class RingFinder():
    def __init__(self,B,R,sB,sR,pixelsize,zerofluxB,zerofluxR,xl=None,yl=None):
        self.pixelsize=pixelsize
        self.zerofluxB=zerofluxB
        self.zerofluxR=zerofluxR
        self.zeroflux=self.zerofluxB

        if xl==None:
            xl=B.shape[0]/2.
        if yl==None:
            yl=B.shape[1]/2.

        self.xl=xl
        self.yl=yl

        x,y=iT.coords(B.shape)
        x-=self.xl
        y-=self.yl
        r=(x**2+y**2)**0.5
        mask=((r<2.7) & (r>0.5))
        
        alpha=B[mask].sum()*1./R[mask].sum()


        self.D=B-alpha*R

        self.S=(sB**2+(alpha*sR)**2)**.5

        
        self.SN=self.D/self.S


    def ringfind(self,cmax=2.7,vb=False):
        self.RegionFind(vb=vb)#regionfind works
        self.MomentChecks(cmax=cmax,vb=vb)
        if len(self.goodlabels)==0:
            if vb:print "found nothing"
            return False
        elif len(self.goodlabels)==1:
            if self.q<0.7: 
                return True
            else:
                if vb:print "1 feature, not elongated"
                return False
        elif len(self.goodlabels)>1:
            return True
        else:
            if vb:print "why was I called?"
            return False

    def RegionFind(self,significancefloor=1.2,vb=False):

        self.significancefloor=significancefloor

        self.SN[self.SN<significancefloor]=0

        Amax=(7./self.pixelsize**2)

        regions, nlbl = ndimage.measurements.label(self.SN)

        lbls = numpy.arange(1, nlbl+1)

        for lbl in lbls:
            lblength= len(regions[regions==lbl])
            if lblength<10:
                self.SN[regions==lbl]=0
            elif lblength>Amax*1000:
                self.SN[regions==lbl]=0
            else:
                pass

        self.D[self.SN==0]=0

        self.regions, self.multiplicity = ndimage.measurements.label(self.D)
        self.labels=numpy.arange(1, self.multiplicity+1)  



    def MomentChecks(self,cmax=2.7,cmin=0.5,vb=False):
        import moments
        import indexTricks as iT
        x,y=iT.coords(self.D.shape)
        x-=self.xl
        y-=self.yl
        x*=self.pixelsize
        y*=self.pixelsize


        goodlabels=[]
        unalignedlabels=[]

        for label in self.labels:
            D2=self.D*1
            D2[self.regions!=label]=0

            #plt.imshow(D2,interpolation="None")
            #plt.show(block=True)



            com,(a,b,q), eigvecs, flux,mean= moments.main(D2)

            va=eigvecs[:,1]

            ca=(com[0]-self.xl,com[1]-self.yl)
            com_dist=((ca[0]**2+ca[1]**2)**0.5)*self.pixelsize
            ca=ca/((ca[0]**2+ca[1]**2)**0.5)
            

            dp=ca[0]*va[0]+ca[1]*va[1]
            theta=(numpy.arccos(dp)*180/3.14159)#-90

            #print theta,self.thetain
            
            #plt.imshow(D2,interpolation="none")
            #plt.show(block=True)

            if com_dist < cmin:#arcseconds 
                if vb:print 21
                continue 
            if com_dist > cmax:#arcseconds 
                if vb:print 22
                continue 
            if b*self.pixelsize<0.2:#arcseconds 
                if vb:print 23
                continue
            if self.FluxCheck(flux,mean)==False:
                if vb:print 24
                continue

            if numpy.abs(theta) >30 and (theta<150 or theta >210):#degrees
                unalignedlabels.append(label)
            else:
                goodlabels.append(label)

            self.q=q

        self.goodlabels=goodlabels
        self.unalignedlabels=unalignedlabels

    def FluxCheck(self, flux,mean):
        zeroflux=self.zeroflux

        if flux<0:return False
        
        mag=10**(2.5*(flux*1./self.zeroflux))
        mag=-2.5*numpy.log10(flux*1./self.zeroflux)

        mean/=self.pixelsize**2 #in arcsec^2

        meanSB=10**(2.5*(mean*1./self.zeroflux))
        meanSB=-2.5*numpy.log10(mean*1./self.zeroflux)

        if mag > 25.5:return False
        if meanSB>26.3: return False
        return True

    
if __name__ == '__main__':
    import pyfits
    imgg=pyfits.open("MockImages/Candidate_978_g_SDSS_image.fits")
    sigg=pyfits.open("MockImages/Candidate_978_g_SDSS_sigma.fits")
    psfg=pyfits.open("MockImages/Candidate_978_g_SDSS_psf.fits")
    imgi=pyfits.open("MockImages/Candidate_978_i_SDSS_image.fits")
    sigi=pyfits.open("MockImages/Candidate_978_i_SDSS_sigma.fits")
    psfi=pyfits.open("MockImages/Candidate_978_i_SDSS_psf.fits")
 
    RF=RingFinder(imgg,imgi,sigg,sigi,psfg,psfi,0.265,30,30)
    RF.thetain=theta
    RF.ringfind()


    continue
    exit()

    for xl in [30,40,50]:
        for yl in [40,50,60]:
            RF=RingFinder(galmodel*1,galmodel*1,sigmodel*1,sigmodel*1,1,30,30,xl,yl)
            print xl,yl
            RF.ringfind()














