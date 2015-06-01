The codes in this project will allow you to simulate observations of the galaxy-galaxy strong lensing population as seen by upcoming surveys.

===============================================================================
LICENSE

The code is open access, but please email thomas.collett@port.ac.uk to tell me that you are using it. Please cite Collett 2015, if you make use of these codes in your work.

The documentation is deliberately sparse - I'd much rather you email me asking to collaborate than use the code as a black box.

===============================================================================
LENS POPULATION

If you don't want to simulate lenses, but do want to know the properties of future samples, we make these available for each of the surveys. These are space seperated text files, with one lens per line. The columns are explained at the start of each file). 

The files are:
    Euclid_lenses.txt
    LSSTa_lenses.txt
    LSSTb_lenses.txt
    LSSTc_lenses.txt
    DESa_lenses.txt
    DESb_lenses.txt
    DESc_lenses.txt

a refers to lenses discoverable in the full co-add
b refers to lenses discoverable in the best single epoch imaging
c refers to lenses discoverable in the optimally stacked coadd

===============================================================================
INSTALLATION

First clone the repo with
   git clone "https://github.com/tcollett/LensPop.git" 

Make a few folders to put things in
   cd LensPop  
   mkdir idealisedlenses
   mkdir LensStats

The code is mostly python, so should work out of the box (if you have standard astrophysical libraries installed already - if you have any problems install anaconda python) except the deflection angles code which must be compiled using:

   cd pylens
   f2py -c -m powerlaw powerlaw.f
   
===============================================================================
HOW TO REPRODUCE COLLETT 2015 
(Timings are based on my ~2011 intel i5 desktop)
First generate idealized lenspopulation:
    python MakeLensPop.py (~7 hours, makes all the lenses on the sky)

Now observe the idealized lens population:
    (with very loose limits on SN and maglims)
    python ModelAll.py CFHT 0.1 (~8 hrs, 0.1 is fraction of sky, not survey fraction)
    python ModelAll.py Euclid 0.1 (~9 hrs)
    python ModelAll.py DES 0.1 (~16 hrs)
    python ModelAll.py LSST 0.1 (~22 hrs)

You can of course cut down shot noise by simulating more sky (upto 1), but it'll take longer. The 0.1 is hard coded into the next part (but can be changed trivially)

Now Make Results:
    python MakeResults.py CFHT
    python MakeResults.py DES
    python MakeResults.py LSST
    python MakeResults.py Euclid

Now make figures showing what you've found
    python MakeFigures34567.py


===============================================================================
HOW TO TEST YOUR OWN RINGFINDER

If your finder uses gaussianised coadds (or single epoch images), simply open ModelAll.py and add your finder at line ~180.

===============================================================================
HOW TO MAKE MOCK COADDS (OR BEST SINGLE EPOCH IMAGING)

Simply find somewhere to put all the data (it'll be a lot!), and uncomment  Lines 190 through 210.

===============================================================================
HOW TO MAKE ALL SINGLE EPOCH IMAGING.

Hack the code (it should be moderately easy) or email me: thomas.collett@port.ac.uk. This will be a huge amount of data [tens of terabytes for all the lenses in LSST, and that's before simulating non-lenses], so plan accordingly!
===============================================================================