
class SPSModel:
    def __init__(self,spex,axes,filters,redshift=None):
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

        self.axes = axes
        self.redshift = redshift
        self.spex = spex

        self.models = {}
        self.models = self.create_models()
        del self.spex


    def luminosity_correction(self):
        """
        The FSPS templates are in L_sun/Angstrom, with L_sun = 3.83e26W
          so we change to the appropriate flux at each relevent redshift.

        Our SED code gives correct magnitudes if the SED is in units of
          ergs/s/cm/cm/Angstrom, and our distance class gives distances in
          Mpc/h, so we will need to use the conversion Mpc = 3.08568e24cm
          and assume a value for h.

        The conversion factor is conv = flux_density/luminosity_density and
          the templates should therefore be multiplied by conv to get the
          correct flux_density. The factor of (1+z) accounts for the fact that
          these are _densities_.

                        3.83e33
            conv =  ----------------
                    (1+z) 4 pi Dl**2
        """
        import distances
        from math import pi
        z = self.redshift
        dist = distances.Distance()
        cm_per_Mpc = 3.08568e24
        if z==0:
            dl = 1e-5*cm_per_Mpc
        else:
            dl = dist.Dl(z)*cm_per_Mpc
        conv = 1./(4*pi*dl**2)
        conv *= 3.83e33
        return conv


    def create_models(self):
        import numpy,cPickle
        from stellarpop import tools
        from ndinterp import ndInterp

        wave = self.axes[2]
        self.axes_names = self.axes[1]
        nmodels = self.spex.size/wave.size
        spex = self.spex.reshape(nmodels,wave.size)
        spex = spex*self.luminosity_correction()

        axes = self.axes[0]

        out = {}
        for F in self.filter_names:
            filt = self.filters[F]
            mags = numpy.zeros(nmodels)
            for i in range(nmodels):
                mags[i] = tools.ABFM(filt,[wave,spex[i]],self.redshift)
            mags = mags.reshape(self.spex.shape[:-1])
            out[F] = ndInterp(axes,mags)
        return out


    def eval(self,points,filter,redshift):
        return self.models[filter].eval(points)
        pnts = []
        npoints = len(points[self.axes_names[0]])
        for i in range(npoints):
            p = []
            for key in self.axes_names:
                p.append(points[key][i])
            pnts.append(p)

        return self.models[filter][redshift].eval(pnts)
