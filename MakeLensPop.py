import distances
import lenspop

# Simulate all lenses on the sky - takes about 7 hours:
# fsky = 1
# Fast test - under a minute:
fsky = 0.001

D = distances.Distance()

Lpop = lenspop.LensPopulation(reset=True, sigfloor=100, zlmax=2, D=D)

Ndeflectors = Lpop.Ndeflectors(2, zmin=0, fsky=fsky)

L = lenspop.LensSample(reset=False, sigfloor=100, cosmo=[0.3,0.7,0.7],
               sourcepop="lsst")

L.Generate_Lens_Pop(int(Ndeflectors), firstod=1, nsources=1, prunenonlenses=True)
