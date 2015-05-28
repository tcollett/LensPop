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

    TODO: It is possible to take N-dimensional slices from the cube instead of
        just points.
    """
    def __init__(self,axes,z,order=3):
        from scipy import ndimage
        import scipy
        self.axes = {}
        for key in axes.keys():
            self.axes[key] = axes[key]
        self.z = z.copy()
        self.spline = ndimage.spline_filter(z,output=scipy.float64,order=order)
        self.order = order


    def evaluate(self,points):
        from scipy import interpolate,ndimage
        indices = []
        for i in range(len(self.axes.keys())):
            indices.append([])
        for i in range(len(points)):
            coords = points[i]
            for j in range(len(coords)):
                index = interpolate.splev(coords[j],self.axes[j])
                indices[j].append(index)
        return ndimage.map_coordinates(self.spline,indices,prefilter=False)


    def eval(self,points):
        return self.evaluate(points)


    def set_order(self,order):
        from scipy import ndimage
        import scipy
        self.order = order
        self.spline = ndimage.spline_filter(self.z,output=scipy.float64,
                                                order=order)
