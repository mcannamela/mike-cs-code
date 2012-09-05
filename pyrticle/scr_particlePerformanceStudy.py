import time
import cPickle as pkl
import particleSolver
reload(particleSolver)
from particleSolver import *

runOn = True
if runOn:

	# alpha = exp(linspace(log(.1), log(1), 15))
	# thres = exp(linspace(log(.0005), log(.025), 15))
	# thermalTimestep = exp(linspace(log(1e-6), log(1e-3), 30))

	alpha = linspace(.25, 1, 15)
	thres = linspace(.001, .01, 10)
	thermalTimestep = linspace(1e-6, 1e-4, 50)


	nNodes = 20
	maxSteps = 10000
	PS = particleSolver(particle(diameter = 50e-6, T =300,thermalTimestep = 5e-5,nNodes = nNodes, nStep = maxSteps), 
										maxSteps = maxSteps, verbose = False)
	PS.p.inject(plasmaField(),x0 = array([0.001,.005,0.]), v0 = array([0.,-10.,0.]))

	nTrials = 1
	elapsed = zeros((alpha.shape[0], thres.shape[0],thermalTimestep.shape[0]))
	T = zeros((alpha.shape[0], thres.shape[0],thermalTimestep.shape[0], nNodes))
	for i in range(alpha.shape[0]):
		PS.alpha = alpha[i]
		for j in range(thres.shape[0]):
			PS.stepErrorThres = thres[j]
			for k in range(thermalTimestep.shape[0]):
				PS.p.thermalTimestep = thermalTimestep[k]
			
				start = time.time()
				PS.solve()
				elapsed[i,j,k] = (time.time()-start)/nTrials
				T[i,j,k,:]=squeeze(PS.solvedSteps[-1])
				PS.reset()
				print str([i,j,k, elapsed[i,j,k]])
else:
	theFile = open('particlePerformance.pkl','rb')
	d = cPickle.load(theFile)
	alpha = d['alpha']
	thres = d['thres']
	thermalTimestep = d['thermalTimestep']
	elapsed = d['elapsed']
	T = d['T']
	theFile.close()
	
f = mlm.fig('time trials')
ml.clf(f)
mlm.cont4(alpha, thres, thermalTimestep, elapsed, c = linspace(np.min(elapsed), np.max(elapsed), 20).tolist(), f=f)
mlm.axe(['alpha', 'thres', 'thermalTimestep,s'], mlm.rgs(alpha,thres, thermalTimestep), f)

# D = zeros((T.shape[0]*T.shape[1], T.shape[0]*T.shape[1]))
# for i in range(T.shape[0]-1):
	# for j in range(T.shape[1]-1):
		# for k in range(i+1, T.shape[0]):
			# for m in range(j+1, T.shape[1]):
				# D[i*T.shape[1]+j, k*T.shape[1]+m] = max(abs(T[i,j,:]-T[k,m,:]))


# f2 = mlm.fig('distance matrix')
# ml.clf(f2)
# mlm.surf3(arange(D.shape[0]),arange(D.shape[1]),D, f2)
# mlm.axe(['ij', 'km', 'distance'], mlm.rgs(arange(D.shape[0]),arange(D.shape[1]), D), f2)
# mlm.out()

# figure(figsize = (20,22))
# hist(D[D!=0].flatten(),100) 

f3 = mlm.fig('mean temperature')
ml.clf(f3)
mlm.cont4(alpha, thres, thermalTimestep, mean(T,axis = 3),c = linspace(np.min(mean(T,axis = 3)), np.max(mean(T,axis = 3)), 20).tolist(), f= f3)
mlm.axe(['alpha', 'thres', 'thermalTimestep'], mlm.rgs(alpha,thres, thermalTimestep), f3)
mlm.out()
        
theFile = open('particlePerformance.pkl','wb')
d = dict(zip(['alpha', 'thres', 'thermalTimestep', 'T', 'elapsed'],[alpha, thres, thermalTimestep, T, elapsed]))
pkl.dump(d,theFile, -1)
theFile.close()