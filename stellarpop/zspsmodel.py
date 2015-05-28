import numpy

class zSPSModel:
    def __init__(self,files,axes,filters,oldModels=None):
        """
        All arguments (files,axes,filters,redshifts) are required.

        `files' is a list of filenames. Each filename points to a pickle'd
            set of SPS spectra.
        `axes' is a dictionary; the keys are each of the parameters of the
            SPS model. Each key points to a dictionary with two keys: `points'
            and `eval.' The `points' value is an array describing the point
            the axis is defined on, and `eval' is a spline model for the axis.
        `filters' is a string (or list of strings) containing the name(s) of
            the filter(s) to be evaulated.
        `redshifts' is the redshift (or list of redshifts) at which the SEDs
            should be evaluated.
        """
        from stellarpop import tools
        self.filters = {}
        if type(filters)==type([]):
            for filter in filters:
                self.filters[filter] = tools.filterfromfile(filter)
        else:
            self.filters[filter] = tools.filterfromfile(filter)
        self.filter_names = self.filters.keys()

        self.axes_names = axes.keys()
        self.naxes = len(axes.keys())

        self.axes = axes
        self.files = files

        self.models = {}
        self.models = self.create_models(oldModels)


    def luminosity_correction(self):
        """
        The galaxev templates are in L_sun/Angstrom, with L_sun = 3.826e26W
          so we change to the appropriate flux at each relevent redshift. We
          also multiply by 1e11, so that the flux is appropriate for a galaxy
          of mass 1e11 M_sun.
        
        Our SED code gives correct magnitudes if the SED is in units of
          ergs/s/cm/cm/Angstrom, and our distance class gives distances in
          Mpc/h, so we will need to use the conversion Mpc = 3.08568e24cm
          and assume a value for h.

        The conversion factor is conv = flux_density/luminosity_density and
          the templates should therefore be multiplied by conv to get the
          correct flux_density. The factor of (1+z) accounts for the fact that
          these are _densities_.

                    1e11   3.826e33
            conv =  ----- ----------
                    (1+z) 4 pi Dl**2
        """
        from stellarpop import distances
        from math import pi
        dist = distances.Distance()
        dist.h = 0.7
        dist.OMEGA_M = 0.3
        dist.OMEGA_L = 0.7
        cm_per_Mpc = 3.08568e24
        corr = []
        for z in self.redshifts:
            if z==0:
                dl = 1e-5*cm_per_Mpc
            else:
                dl = dist.Dl(z)*cm_per_Mpc
            conv = 1/(4*pi*dl**2)
            conv *= 3.826e33
            corr.append(conv)

        return corr


    def create_models(self,old=None):
        import numpy,cPickle
        from stellarpop import tools
        from ndinterp import ndInterp

        index = {}
        shape = []
        axes = {}
        inferRedshift = True
        axes_index = 0
        # Only axes with more than one element are interesting
        for key in self.axes_names:
            if self.axes[key]['points'].size<=1:
                continue
            index[key] = {}
            shape.append(self.axes[key]['points'].size)
            axes[axes_index] = self.axes[key]['eval']
            axes_index += 1
            for i in range(self.axes[key]['points'].size):
                index[key][self.axes[key]['points'][i]] = i
        outKeys = index.keys()
        self.redshifts = self.axes['redshift']['points']
        if len(self.redshifts)==1:
            inferRedshift = False
        self.corrections = self.luminosity_correction()

        if old is not None:
            models = old
            for f in self.filter_names:
                model = models[f].copy()
                if numpy.isnan(model).any():
                    models[f] = None
                else:
                    models[f] = ndInterp(axes,model,order=1)
            return models

        zindex = self.axes_names.index('redshift')

        models = {}
        model = numpy.empty(shape)*numpy.nan
        for f in self.filter_names:
            models[f] = model.copy()

        for file in self.files:
            f = open(file,'rb')
            data = cPickle.load(f)
            wave = cPickle.load(f)
            f.close()
            for key in data.keys():
                obj = data[key]
                jj = key
                spec = obj['sed']
                ind = []
                for key in self.axes_names:
                    if key=='redshift' or key not in outKeys:
                        continue
                    try:
                        ind.append([index[key][obj[key]]])
                    except:
                        print key,index[key]
                        print obj
                        df
                if inferRedshift:
                    ind.insert(zindex,None)
                for f in self.filter_names:
                    for i in range(self.redshifts.size):
                        z = self.redshifts[i]
                        # correction is the units correction factor
                        correction = self.corrections[i]
                        sed = [wave,spec*correction]
                        mag = tools.ABFilterMagnitude(self.filters[f],sed,z)
                        if numpy.isnan(mag)==True:
                            df
                        if inferRedshift:
                            ind[zindex] = i
                        models[f][ind] = mag
        self.interpAxes = axes
        for f in self.filter_names:
            model = models[f].copy()
            if numpy.isnan(model).any():
                models[f] = None
            else:
                models[f] = ndInterp(axes,model,order=1)
        return models

    def eval(self,points,filter):
        pnts = []
        npoints = len(points[self.axes_names[0]])
        for i in range(npoints):
            p = []
            for key in self.axes_names:
                p.append(points[key][i])
            pnts.append(p)

        return self.models[filter].eval(pnts)
