import numpy

def coords(shape):
    return numpy.indices(shape).astype(numpy.float64)

def overSample(shape,factor):
    coords = numpy.indices((shape[0]*factor,shape[1]*factor)).astype(numpy.float64)/factor - 0.5*(factor-1.)/factor
    return coords[0],coords[1]

def resamp(a,factor,add=False):
    arr = a.copy()
    x = arr.shape[1]/factor
    y = arr.shape[0]/factor
    o = numpy.zeros((y,x))
    for i in range(factor):
        for j in range(factor):
            o += arr[i::factor,j::factor]
    if add==True:
        return o
    return o/factor**2

resample = resamp

def recube(a,factor):
    arr = a.copy()
    x = arr.shape[1]/factor
    y = arr.shape[0]/factor
    o = numpy.empty((x*y,factor,factor))
    for i in range(factor):
        for j in range(factor):
            o[:,i,j] = arr[i::factor,j::factor].ravel()
    return o

def mesh(x,y):
    grid=coords((x.size,y.size))
    grid[0]=grid[0]*(x[1]-x[0])+x[0]
    grid[1]=grid[1]*(y[1]-y[0])+y[0]
    return grid

