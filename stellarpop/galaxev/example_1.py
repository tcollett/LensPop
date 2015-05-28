######################################################################
#
# EXAMPLE 1: An example to create a SPS model from the BC03 Chabrier models
#
######################################################################


#
# Part 1: Create SEDs (this is a one-time step that will take hours)
#
import create_models as galaxev

iseds = ['bc2003_lr_m%d2_chab_ssp.ised'%i for i in range(2,8)]

# First create the CSP models
for ised in iseds:
    galaxev.create_csp_models(ised)

lookup = {'m22':0.0001,'m32':0.0004,'m42':0.004,'m52':0.008,'m62':0.02,
            'm72':0.05}

# Then create the SED models
for ised in iseds:
    Z = lookup[ised.split('_')[2]]

    outname = 'chabrier_Z=%6.4f.dat'%Z
    galaxev.create_seds(ised,outname)



#
# Part 2: Create the SPS model for a particular galaxy (takes ~3-15minutes).
#
# It is probably a good idea to split this out into it's own script.
#

# A redshift (or list of redshifts) and filter set must be specified. Also
#   specifiy an output name for pickle'ing.
redshift = 0.6459
filters = ['u_SDSS','g_SDSS','r_SDSS','i_SDSS','z_SDSS']
outname = 'mySPSmodel.dat'

import glob,cPickle,numpy
from stellarpop import spsmodel
from scipy import interpolate

# Grab the SED models
files = glob.glob('chabrier_Z*dat')
axes = {'Z':[],'age':None,'tau':None,'tau_V':None}

for file in files:
    z = float(file.split('=')[1].split('_')[0])
    axes['Z'].append(z)
axes['Z'] = numpy.sort(numpy.unique(numpy.asarray(axes['Z'])))

# We could open any of the files; they'll all have the same tau, tau_V, age
#   arrays.
f = open(file)
tmp = cPickle.load(f) # Trash read of output
del tmp
tmp = cPickle.load(f) # Trash read of wave
del tmp
t = cPickle.load(f)
tV = cPickle.load(f)
a = cPickle.load(f)
f.close()

axes['age'] = a
axes['tau'] = t
axes['tau_V'] = tV

# The interpolation order of the input data. Use linear (order=1) for the
#   crazy 'Z' sampling chosen by BC03.
order = {'Z':1,'age':3,'tau':3,'tau_V':3}
axes_models = {}
interpmodels = {}

for key in axes.keys():
    arr = axes[key]
    axes_models[key] = {}
    axes_models[key]['points'] = numpy.sort(arr)
    axes_models[key]['eval'] = interpolate.splrep(arr,numpy.arange(arr.size),k=order[key],s=0)

model = spsmodel.SPSModel(files,axes_models,redshift,filters)

f = open(outname,'wb')
cPickle.dump(model,f,2)
f.close()

