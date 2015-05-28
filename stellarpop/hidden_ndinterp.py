def create_axes_array(axes):
    """
    Given a list of N axes of length {a,b,c,...}, returns an N+1 dimension
        array of shape {N,a,b,c,...} describing the coordinates at each point
        in the grid.
    """
    import scipy
    ndim = len(axes)
    shape = [ndim]
    for i in axes:
        shape.append(i.size)

    coords = scipy.ones(shape)
    for i in range(ndim):
        coords[i] = scipy.rollaxis(scipy.rollaxis(coords[i],i,ndim)*axes[i],ndim-1,i)

    return coords


class ndInterp:
    """
    The ndInterp class is an interpolation model of an N-dimensional data cube.
        It is instantiated with a list of axes describing the dimensions of the
        cube and the cube itself. The model can be evaluated at discrete points
        within the cube -- points outside of the cube are evaluated as 0.
    """
    def __init__(self,axes,z,order=3):
        from scipy import ndimage
        import scipy
        self.axes = {}
        for key in axes.keys():
            self.axes[key] = axes[key]
        z = z.astype(scipy.float64)
        self.z = z.copy()
        if order==1:
            self.spline = z.copy()
        else:
            self.spline = ndimage.spline_filter(z,output=scipy.float64,order=order)
        self.order = order


    def evaluate_old(self,points):
        from scipy import interpolate,ndimage
        import numpy
        indices = []
        for i in range(len(self.axes.keys())):
            indices.append([])
        for i in range(len(points)):
            coords = points[i]
            for j in range(len(coords)):
                index = interpolate.splev(coords[j],self.axes[j])
                indices[j].append(index)
        return ndimage.map_coordinates(self.spline,indices,prefilter=False)

    def evaluate(self,points):
        from scipy import interpolate,ndimage
        import numpy
        points = numpy.array(points)
        if points.ndim==1:
            points = numpy.atleast_2d(points).T
        indices = numpy.empty((points.shape[1],points.shape[0]))
        for i in range(points.shape[-1]):
            indices[i] = interpolate.splev(points[:,i],self.axes[i])
        return ndimage.map_coordinates(self.spline,indices,prefilter=False)

    def eval(self,points):
        return self.evaluate(points)


    def set_order(self,order):
        from scipy import ndimage
        import scipy
        self.order = order
        if order==1:
            self.spline = self.z.copy()
            return
        self.spline = ndimage.spline_filter(self.z,output=scipy.float64,
                                                order=order)
