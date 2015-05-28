#!/usr/bin/python -W ignore

USAGE = """
NAME
        mag2mag

PURPOSE
        Given a galaxy of type T at redshift z1 with magnitude m1 in 
        filter1, compute magnitude m2 in filter2 for the same galaxy at 
        redshift z2.
        
COMMENTS
        Tricks the hyper-z photometric redshift code into doing this.
        PGPLOT module requires /usr/local/bin/perl at KIPAC.

USAGE
        mag2mag [flags] [options]

FLAGS
        -u           Print this message [0]
        -q           quiet operation [0]
        -vega        Input magnitudes are in Vega system [def=AB]
        -convert     Convert from one magnitude system to the other [def=0]
        -H0          Hubble constant [70]
        -Om          Omega matter [0.3]
        -OL          Omega lambda [0.7]
        -plot        Illustrate results with a nice plot
        -noX         Do the plotting quietly (not to xwin)
        -test        Do an example run (the first one on the list below)
        -no-clean    Keep various output files for posterity
       
INPUTS
        -m1       f         input magnitude
        -f1       s         input filter
        -z1       f         input redshift
        -T        s         galaxy spectral type
        -f2       s         output filter
        -z2       f         output redshift
        stdin               Provide compulsory inputs when prompted [Not yet functional]

OPTIONAL INPUTS

OUTPUTS
        stdout       Useful information

EXAMPLES
        
        mag2mag -T CWW_Sbc_ext -m1 25 -f1 'Johnson_H'    -z1 0.6 \\
                 -plot -vega -convert -f2 'F814W_ACSWFC' -z2 1.4

        mag2mag -T CWW_E_ext -m1 -21.43 -f1 'Johnson_B' -z1 0.0 \\
                   -plot -vega -convert -f2 'F775W_ACSWFC' -z2 0.321

        mag2mag -T CWW_E_ext -m1 18.302625 -f1 'F775W_ACSWFC' -z1 0.321 \\
                            -plot -convert -f2 'Johnson_B'    -z2 0.0
"""

def parse_cmdline(args,template):
    options = {}
    parse = [i for i in args]

    del parse[0] # Delete program name
    while len(parse):
        flag = parse[0][1:]
        del parse[0]
        if flag.find('=')>=0:
            value = flag.split('=')[1]
        elif len(parse)>0 and parse[0][0]!='-':
            value = parse[0]
            del parse[0]
        else:
            value = True
        options[flag] = value
    flags = options.keys()

    output = []
    for item in template:
        tmp = item.split('=')
        option = tmp[0]
        if option not in flags:
            output.append(None)
            continue
        if len(tmp)>1:
            if tmp[1]=='f':
                output.append(float(options[option]))
            elif tmp[1]=='i':
                output.append(int(options[option]))
            elif tmp[1]=='s':
                output.append(options[option])
            else:
                output.append(True)
        else:
            output.append(True)

    return output

import sys

quiet,m1,f1,z1,SEDtype,f2,z2,vega,convert,H0,Om,OL,plot,noX,test,noclean,help = parse_cmdline(sys.argv,['q','m1=f','f1=s','z1=f','T=s','f2=s','z2=f','vega','convert','H0=f','Om=f','OL=f','plot','noX','test','no-clean','u'])

if help is not None:
    print USAGE
    sys.exit()

if H0 is None:
    H0 = 70.
if Om is None:
    Om = 0.3
if OL is None:
    OL = 0.7

if test is not None:
    SEDtype = 'Sbc_cww'
#    SEDtype = 'CWW_Sbc_ext'
    m1 = 25.
    f1 = 'H_Johnson'
    z1 = 0.6
    f2 = 'F814W_WFC'
    z2 = 1.4
    cmd = "%s -T %s -m1 %f -f1 %s -z1 %f -f2 %s -z2 %f -plot \
                -vega -convert" % (sys.argv[0],SEDtype,m1,f1,z1,f2,z2)
    import os
    os.system(cmd)
    sys.exit()

if vega is None:
    vega = False # No vega flag = use AB mags

if convert is None:
    convert = False

if m1 is None or f1 is None or z1 is None or SEDtype is None:
    print "Error: Incomplete input information. Use -u for usage."
    sys.exit()

from stellarpop import tools
f1 = tools.filterfromfile(f1)
if f2 is None:
    f2 = f1
else:
    f2 = tools.filterfromfile(f2)
if z2 is None:
    z2 = z1

SED = tools.getSED(SEDtype)

if vega:
    filtermag = tools.VegaFilterMagnitude(f1,SED,z1)
#    vegatmp = tools.getSED('Vega')
#    vegaAB = tools.ABFilterMagnitude(f1,vegatmp,0.)
else:
    filtermag = tools.ABFilterMagnitude(f1,SED,z1)
magoffset = m1-filtermag
#print magoffset

if vega:
    m2 = tools.VegaFilterMagnitude(f2,SED,z2)+magoffset
else:
    m2 = tools.ABFilterMagnitude(f2,SED,z2)+magoffset
if z1!=z2:
    from stellarpop import distances
    from math import log10
    dist = distances.Distance()
    dist.OMEGA_M = Om
    dist.OMEGA_L = OL
    dist.h = H0/100.
    if z2!=0.:
        dimming = dist.Dl(z1)/dist.Dl(z2)
    else:
        dimming = dist.Dl(z1)/1e-5
    m2 -= 5.*log10(dimming)

if convert:
    vegatmp = tools.getSED('Vega')
    vegatmp = (vegatmp[0],vegatmp[1])
    vegaAB = tools.ABFilterMagnitude(f2,vegatmp,0.)
    if vega:
        m2 += vegaAB
    else:
        m2 -= vegaAB
print m2

if plot is not None:
    import pylab
    from scipy.interpolate import splev
    pylab.plot(SED[0]*(1.+z1),SED[1]/SED[1].mean(),c='b')
    if z1!=z2:
        pylab.plot(SED[0]*(1.+z2),SED[1]/SED[1].mean(),c='r')
    wave = SED[0]*(1.+z1)
    cond = (wave>=f1[0][0])&(wave<=f1[0][-1])
    pylab.plot(wave[cond],splev(wave[cond],f1),c='b')
    if f1!=f2:
        wave = SED[0]*(1.+z2)
        cond = (wave>=f2[0][0])&(wave<=f2[0][-1])
        pylab.plot(wave[cond],splev(wave[cond],f2),c='r')
    pylab.show()
