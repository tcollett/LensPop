
class MassModel:
    def __init__(self,files,axes):
        """
        All arguments (files,axes,filters,redshifts) are required.

        `files' is a list of filenames. Each filename points to a pickle'd
            set of SPS spectra.
        `axes' is a dictionary; the keys are each of the parameters of the
            SPS model. Each key points to a dictionary with two keys: `points'
            and `eval.' The `points' value is an array describing the point
            the axis is defined on, and `eval' is a spline model for the axis.
        """
        self.axes_names = axes.keys()
        self.naxes = len(axes.keys())

        self.axes = axes
        self.files = files

        self.mass_model,self.gas_model = self.create_model()

    def create_model(self):
        import scipy,cPickle
        from stellarpop import tools
        from stellarpop.ndinterp import ndInterp

        index = {}
        shape = []
        axes = {}
        axes_index = 0
        for key in self.axes_names:
            index[key] = {}
            shape.append(self.axes[key]['points'].size)
            axes[axes_index] = self.axes[key]['eval']
            axes_index += 1
            for i in range(self.axes[key]['points'].size):
                index[key][self.axes[key]['points'][i]] = i

        mass_model = scipy.empty(shape)*scipy.nan
        gas_model = mass_model.copy()

        for file in self.files:
            f = open(file,'rb')
            data = cPickle.load(f)
            wave = cPickle.load(f)
            f.close()
            for key in data.keys():
                obj = data[key]
                spec = obj['sed']
                ind = []
                for key in self.axes_names:
                    ind.append([index[key][obj[key]]])
                mass_model[ind] = obj['mass']
#                gas_model[ind] = obj['gas']

        mass_model = ndInterp(axes,mass_model)
#        gas_model = ndInterp(axes,gas_model)
        return mass_model,gas_model

    def stars(self,points):
        pnts = []
        npoints = len(points[self.axes_names[0]])
        for i in range(npoints):
            p = []
            for key in self.axes_names:
                p.append(points[key][i])
            pnts.append(p)

        return self.mass_model.eval(pnts)

    def gas(self,points):
        pnts = []
        npoints = len(points[self.axes_names[0]])
        for i in range(npoints):
            p = []
            for key in self.axes_names:
                p.append(points[key][i])
            pnts.append(p)

        return self.gas_model.eval(pnts)

