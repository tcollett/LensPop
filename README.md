The codes in this project will allow you to simulate observations of the galaxy-galaxy strong lensing population as seen by upcoming surveys.

===============================================================================
LICENSE

The code is open access, but please email thomas.collett@port.ac.uk to tell me that you are using it. Please cite Collett 2015, if you make use of these codes in your work.

===============================================================================
LENS POPULATION

If you don't want to simulate lenses, but do want to know the properties of future samples, we make these available for each of the surveys. These are space seperated text files, with one lens per line:

Euclid_lenses.txt
LSST_lenses.txt
DES_lenses.txt

The columns in the tables are 








===============================================================================
INSTALLATION

The code is mostly python, so should work out of the box (if you have standard astrophysical libraries installed already) except the deflection angles code which must be compiled using:

   cd pylens
   f2py -c -m powerlaw powerlaw.f
   
===============================================================================
HOW TO REPRODUCE COLLETT 2015 

first generate idealized lenspopulation:
    python   MakeLensPop.py (~7 hours, makes all the lenses on the sky)

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



===============================================================================
HOW TO TEST YOUR OWN RINGFINDER

If your finder uses gaussianised coadds (or single epoch images), simply open ModelAll.py and add your finder at line ~180.

===============================================================================
HOW TO MAKE MOCK COADDS (OR BEST SINGLE EPOCH IMAGING)

Simply find somewhere to put all the data (it'll be a lot!), and uncomment  Lines 190 through 210.

===============================================================================
HOW TO MAKE SINGLE EPOCH IMAGING.

Hack the code (it should be pretty easy) or email me: thomas.collett@port.ac.uk
===============================================================================
