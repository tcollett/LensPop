UVISfilters = """<option VALUE="F218W">F218W (ISM feature)
                <option VALUE="F200LP">F200LP (Clear)
                <option VALUE="F225W">F225W (UV)
                <option VALUE="F275W">F275W (UV)
                <option VALUE="F300X">F300X (Long Pass)
                <option VALUE="F336W">F336W (U)
                <option VALUE="F350LP">F350LP (Long Pass)
                <option VALUE="F390W">F390W (Washington C)
                <option VALUE="F438W">F438W (WFPC2 B)
                <option VALUE="F475W">F475W (SDSS g)
                <option VALUE="F475X">F475X (Long Pass)
                <option VALUE="F555W" SELECTED>F555W (WFPC2 V)
                <option VALUE="F600LP">F600LP (Long Pass)
                <option VALUE="F606W">F606W (WFPC2 Wide V)
                <option VALUE="F625W">F625W (SDSS r)
                <option VALUE="F775W">F775W (SDSS i)
                <option VALUE="F814W">F814W (WFPC2 Wide I)
                <option VALUE="F850LP">F850LP (SDSS z)"""

IRfilters = """ <option VALUE="F098M">F098M (Blue grism reference)
                <option VALUE="F105W">F105W (Fat Sloan z)
                <option VALUE="F110W" SELECTED>F110W (Wide J)
                <option VALUE="F125W">F125W (Broad J)
                <option VALUE="F140W">F140W (JH gap)
                <option VALUE="F126N">F126N (FeII)
                <option VALUE="F127M">F127M (Water/CH4 continuum)
                <option VALUE="F128N">F128N (Paschen &beta;)
                <option VALUE="F130N">F130N (Paschen &beta; continuum)
                <option VALUE="F132N">F132N (Paschen &beta; redshifted)
                <option VALUE="F139M">F139M (Water/CH4 line)
                <option VALUE="F153M">F153M (H20 & NH3)
                <option VALUE="F160W">F160W (Broad H)
                <option VALUE="F164N">F164N ([FeII])
                <option VALUE="F167N">F167N (FeII])"""

f = open('maketemplate.cl','w')
for line in UVISfilters.split('\n'):
    filter = line.split('"')[1].lower()
    f.write('calcband wfc3,uvis1,%s %s_thpt\n'%(filter,filter))
    f.write('tdump "%s_thpt" datafile="%s_WFC3UVIS.res" columns="wavelength,throughput"\n'%(filter,filter.upper()))


for line in IRfilters.split('\n'):
    filter = line.split('"')[1].lower()
    f.write('calcband wfc3,ir,%s %s_thpt\n'%(filter,filter))
    f.write('tdump "%s_thpt" datafile="%s_WFC3IR.res" columns="wavelength,throughput"\n'%(filter,filter.upper()))

f.close()
