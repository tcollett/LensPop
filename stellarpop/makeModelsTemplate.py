import glob,cPickle,numpy
from stellarpop import zspsmodel
from scipy import interpolate
import sys,time

# Grab the SED models
files = glob.glob('/local/mauger/STELLARPOP/bc03models/chabrier_Z*dat')
axes = {'Z':[],'age':None,'tau':None,'tau_V':None}

for file in files:
    z = float(file.split('=')[1].split('.dat')[0])
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
order = {'Z':1,'age':1,'tau':1,'tau_V':1}
axes_models = {}
interpmodels = {}

for key in axes.keys():
    arr = numpy.sort(axes[key])
    axes_models[key] = {}
    axes_models[key]['points'] = arr
    axes_models[key]['eval'] = interpolate.splrep(arr,numpy.arange(arr.size),k=order[key],s=0)

arr = numpy.arange(0.01,5.001,0.05)
axes_models['redshift'] = {'points':arr.copy()}
axes_models['redshift']['eval'] = interpolate.splrep(arr,numpy.arange(arr.size),k=1,s=0)

outname = "B1608_chabBC03_z5.model"
filters = ['g_gunn','r_gunn','i_gunn','F606W_ACS','F814W_ACS']

print "Creating model: ",outname
t = time.time()
f = open(outname,'wb')
model = zspsmodel.zSPSModel(files,axes_models,filters)
cPickle.dump(model,f,2)
f.close()
print "Time elapsed: ",(time.time()-t)

