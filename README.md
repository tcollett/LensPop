# LensPop

Simulating observations of the galaxy-galaxy strong lensing population as expected to be seen in upcoming surveys.

----

## License, Credits etc

The code is open access, but please [email me]( mailto:thomas.collett@port.ac.uk) to tell him that you are using it. If you use any of the software or ideas in this repo in your own research, please cit "Collett (2015)."

The documentation is deliberately sparse - I'd much rather you email me asking to collaborate than use the code as a black box.

-- Tom Collett.

----

## Pre-made Lens Populations

If you don't want to simulate lenses, but do want to know the properties of future samples, we make these available for each of the surveys. These are space-separated text files, with one lens per line. The columns are explained at the start of each file).

The files are:
```
    Euclid_lenses.txt
    LSSTa_lenses.txt
    LSSTb_lenses.txt
    LSSTc_lenses.txt
    DESa_lenses.txt
    DESb_lenses.txt
    DESc_lenses.txt
```

* `a` refers to lenses discoverable in the full co-add.
* `b` refers to lenses discoverable in the best single epoch imaging.
* `c` refers to lenses discoverable in the optimally stacked coadd.

----

## Installation

### Users

The following _should_ work:
```
   pip install git+git://github.com/drphilmarshall/LensPop.git#egg=lenspop
```

### Developers

First clone the repo with
```
   git clone "https://github.com/tcollett/LensPop.git"
```
and then set up your path with
```
   python setup.py develop
```

### Notes

You'll need to make a few folders to put things in:
```
   mkdir idealisedlenses
   mkdir LensStats
```
The code is mostly python, so should work out of the box (if you have standard astrophysical libraries installed already - if you have any problems install anaconda python) except the deflection angles code which must be compiled using:
```
   cd pylens
   f2py -c -m powerlaw powerlaw.f
```

----

## Recipes

### How to Reproduce Collett (2015)

Timings are based on my ~2011 intel i5 desktop.
First generate idealized lenspopulation:
```
    python MakeLensPop.py
```
To make all the lenses on the sky, edit this script to use `fsky=1` - this will take about 7 hours.

Now observe the idealized lens population (with very loose limits on SN and maglims):
```
    python ModelAll.py CFHT 0.1 (~8 hrs, 0.1 is fraction of sky, not survey fraction)
    python ModelAll.py Euclid 0.1 (~9 hrs)
    python ModelAll.py DES 0.1 (~16 hrs)
    python ModelAll.py LSST 0.1 (~22 hrs)
```
You can of course cut down shot noise by simulating more sky (up to 1), but it'll take longer. The 0.1 is hard coded into the next part (but can be changed trivially).

Now Make Results:
```
    python MakeResults.py CFHT
    python MakeResults.py DES
    python MakeResults.py LSST
    python MakeResults.py Euclid
```
Finally, make the figures showing what you've found:
```
    python MakeFigures34567.py
```

----

### How to Test Your Own RingFinder

If your finder uses gaussianised coadds (or single epoch images), simply open `ModelAll.py` and add your finder at line ~180.

----

### How to Make Mock Coadds (or Best Single Epoch Imaging)

Simply find somewhere to put all the data (it'll be a lot!), and uncomment  Lines 190 through 210.

----

### How to Make All Single Epoch Imaging

Hack the code (it should be moderately easy) or email me: thomas.collett@port.ac.uk. This will be a huge amount of data [tens of terabytes for all the lenses in LSST, and that's before simulating non-lenses], so plan accordingly!
