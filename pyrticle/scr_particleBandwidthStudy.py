import time
import cPickle as pkl
import particleSolver
reload(particleSolver)
from particleSolver import *

def cutoffStep(x, delta, axis = 0, gamma = .95):
	"""
	find the step size for which we retain a fraction gamma of the energy in x
	along axis axis, which has samples spaced delta apart
	"""
	nfft = 2**13
	X = abs(rfft(x,nfft, axis = axis))**2
	e = cumsum(X, axis = axis)/expand_dims(sum(X, axis = axis), axis)
	
	I = argmin(abs(e-.9), axis = axis)
	i = double(np.max(I))
	return delta/(i/nfft)

	
if __name__== "__main__":
	pdb.set_trace()
	runOn = True
	if runOn:
		#fix the solver parameters at something conservative
		alpha = .5
		thres = .0001
		maxSteps = 10000
		
		#we'll vary the temporal and spatial resolution for some 
		#diameters and injection velocities of interest
		thermalTimestep = linspace(1e-5, 1e-4, 5)
		nNodes = arange(5,100, 5)
		d = linspace(10e-6, 150e-6, 5)
		v0 = -1*linspace(1, 20,5)

		#create the plasmaField
		try:
			PF.density(array([0,0,0]))
		except NameError:
			PF = plasmaField()
		
		#looking for the smallest temporal step to resolve the drag force, transient temperature
		#and smallest volumetric step size we can use for the temperature profile
		forceCutoffStep = zeros((d.shape[0], v0.shape[0], thermalTimestep.shape[0], nNodes.shape[0]))
		thermalCutoffStep = zeros_like(forceCutoffStep)
		thermalSpatialStep = zeros_like(forceCutoffStep)
		
		#for every diamter and velocity, sweep over resolutions and find the cutoff step sizes
		rg = lambda x_: range(x_.shape[0])
		cnt = 0
		for i in rg(d):
			for j in rg(v0):
				for k in rg(thermalTimestep):
					for m in rg(nNodes):
						print ('progress: '+str(double(cnt)/forceCutoffStep.size)+ 'd: '+ str(d[i])
								+ 'v0: '+ str(v0[j])
								+ 'dt: '+ str(thermalTimestep[k])
								+ 'N: '+ str(nNodes[m]))
								
						#particle with this iteration's diameter and spatio-temporal resolution
						p = particle(diameter = d[i], T =300,thermalTimestep = thermalTimestep[k],nNodes = nNodes[m], nStep = maxSteps)
						
						#create the solver
						PS = particleSolver(p, maxSteps = maxSteps)
						
						#inject at the given velocity
						PS.p.inject(PF,x0 = array([0.001,.005,0.]), v0 = array([0.,v0[j],0.]))
						
						#solve it
						PS.solve()
						
						#get the spacetime domain drag force and temperature
						f = array(PS.p.forceHistory)
						T = PS.getSoln()
						
						forceCutoffStep[i,j,k,m] = cutoffStep(f, thermalTimestep[k])
						thermalCutoffStep[i,j,k,m] = cutoffStep(T, thermalTimestep[k], axis = 1)
						thermalSpatialStep[i,j,k,m] =  cutoffStep(T, d[i]**3/nNodes[m], axis = 1)
						cnt+=1
	theFile = open('particleBandwidth.pkl','wb')
	D = dict(zip(['d', 'v0', 'thermalTimestep', 'nNodes', 'forceTemporalStep', 'thermalTemporalStep', 'thermalSpatialStep'],
						[d, v0, thermalTimestep, nNodes, forceCutoffStep, thermalCutoffStep, thermalSpatialStep]))
	pkl.dump(D,theFile, -1)
	theFile.close()

