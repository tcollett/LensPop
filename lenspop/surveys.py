from __future__ import absolute_import, division, print_function

from six.moves import cPickle as pickle
import numpy

class Survey():
    def  __init__(self,Name):
        self.zeroexposuretime=1
        self.strategy="resolve"
        self.strategyx=1
        if Name[:3]=="DES":
            self.pixelsize=0.263
            self.side=76
            self.bands=['g','r','i']
            self.zeropoints=[30,30,30]
            self.zeroexposuretime=90.
            self.skybrightnesses=[21.7,20.7,20.1]
            self.exposuretimes=[900,900,900]
            self.gains=[4.5,4.5,4.5]
            self.seeing=[.9,.9,.9]
            self.nexposures=10
            self.degrees_of_survey=5000
            self.readnoise=(10/4.5)
            twodg=pickle.load(open("2dpdfs/2dg_DES.pkl",'r'))
            twodr=pickle.load(open("2dpdfs/2dr_DES.pkl",'r'))
            twodi=pickle.load(open("2dpdfs/2di_DES.pkl",'r'))
            self.stochasticobservingdata=[twodg,twodr,twodi]
            if Name=="DESsv":
                self.degrees_of_survey=150
            if Name=="DESsv" or  Name =="DESa":
                self.strategy="absolute"
                self.strategyx=10
            if Name=="DESb":
                self.strategy="best"
                self.strategyx=1
            if Name=="DESdummy":
                self.strategy="absolute"
                self.strategyx=10
                dumg=numpy.array([[1.2,21.7],[1.2,21.7]])
                dumr=numpy.array([[0.95,20.7],[0.95,20.7]])
                dumi=numpy.array([[0.95,20.1],[0.95,20.1]])
                print("dummy seeing,strat")
                self.stochasticobservingdata=[dumg,dumr,dumi]
                self.strategy="absolute"
                self.strategyx=10



        elif Name[:4]=="LSST":
            self.pixelsize=0.18
            self.side=111
            self.bands=['g','r','i']
            self.zeropoints=[30,30,30]
            self.zeroexposuretime=25
            self.skybrightnesses=[21.7,20.7,20.1]
            self.exposuretimes=[3000,6000,6000]
            self.gains=[4.5,4.5,4.5]
            self.seeing=[.4,.4,.4]
            self.nexposures=100
            self.degrees_of_survey=18000
            self.readnoise=(10/4.5)
            twodg=pickle.load(open("2dpdfs/2dg_LSST.pkl",'r'))
            twodr=pickle.load(open("2dpdfs/2dr_LSST.pkl",'r'))
            twodi=pickle.load(open("2dpdfs/2di_LSST.pkl",'r'))
            self.stochasticobservingdata=[twodg,twodr,twodi]
            if Name[-1]=="a":
                self.strategy="absolute"
                self.strategyx=10
            if Name[-1]=="b":
                self.strategy="best"
                self.strategyx=1
            if Name[-1]=="c":
                self.strategy="resolve"
                self.strategyx=1



        elif Name=="CFHT" or Name=="CFHTa":
            self.pixelsize=0.187
            self.side=107
            self.bands=['g','r','i']
            self.zeropoints=[26.96,26.47,26.24]
            self.zeroexposuretime=1
            self.skybrightnesses=[21.9,20.6,19.2]
            self.exposuretimes=[3500,5500,5500]
            self.gains=[1.62,1.62,1.62]
            self.seeing=[.8,.8,0.8]
            self.nexposures=1 # this isn't actually true, but the 2d pdfs
                              # are for the CFHT coadds (approximately)
            self.degrees_of_survey=150
            self.readnoise=(5)
            twodg=pickle.load(open("2dpdfs/2dg_CFHT.pkl",'r'))
            twodr=pickle.load(open("2dpdfs/2dr_CFHT.pkl",'r'))
            twodi=pickle.load(open("2dpdfs/2di_CFHT.pkl",'r'))
            self.stochasticobservingdata=[twodg,twodr,twodi]
            self.strategy="absolute"
            self.strategyx=10

        elif Name=="HSC":
            self.pixelsize=0.17
            self.side=200
            self.bands=['g','r','i']
            self.zeropoints=[30,30]
            self.zeroexposuretime=90./(8.2/4)**2
            self.skybrightnesses=[21.9,19.2]
            self.exposuretimes=[600,600]
            self.gains=[4.5,4.5]
            self.seeing=[.8,.8]
            self.nexposures=10
            self.degrees_of_survey=1370
            self.readnoise=(10/4.5)
            twodg=numpy.array([[0.8,21.9],[0.8,21.9]])
            twodi=numpy.array([[0.8,19.2],[0.8,19.2]])
            self.stochasticobservingdata=[twodg,twodi]

        elif Name=="COSMOS":
            pass


        elif Name=="Euclid":
            self.pixelsize=0.1
            self.side=200
            self.bands=['VIS']
            self.zeropoints=[25.5]
            self.zeroexposuretime=1.
            self.skybrightnesses=[22.2]
            self.exposuretimes=[1610]
            self.gains=[1]
            self.seeing=[.2]
            self.nexposures=4
            self.degrees_of_survey=20000
            self.readnoise=(4.5)
            twodVIS=numpy.array([[0.17,22.2],[0.17,22.2]])
            self.stochasticobservingdata=[twodVIS]

        elif Name=="ideal":
            self.pixelsize=0.05
            self.side=400
            self.bands=['g','r','i']
            self.zeropoints=[24.5,24.5,24.5]
            self.zeroexposuretime=4.
            self.skybrightnesses=[220,220,220]
            self.exposuretimes=[22400000000,22400000000,22400000000]
            self.gains=[1,1,1]
            self.seeing=[.05,0.05,0.05]
            self.nexposures=1
            self.readnoise=(.005)
            twodr=numpy.array([[0.1,220],[0.1,220]])
            twodg=numpy.array([[0.1,220],[0.1,220]])
            twodi=numpy.array([[0.1,220],[0.1,220]])
            self.stochasticobservingdata=[twodg,twodr,twodi]
            self.degrees_of_survey=41253
        else:
            print("I don't know that survey")
            exit()



        #convert bandnames into the required formats
        for i in range(len(self.bands)):
            bandname=self.bands[i]
            if bandname=="g":
                self.bands[i]="g_SDSS"
            if bandname=="r":
                self.bands[i]="r_SDSS"
            if bandname=="i":
                self.bands[i]="i_SDSS"
            if bandname=="z":
                self.bands[i]="z_SDSS"
            if bandname=="F814":
                self.bands[i]="F814W_ACS"

        degrees_of_whole_sky=41253.
        self.f_sky=float(self.degrees_of_survey)/degrees_of_whole_sky
