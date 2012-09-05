import time
import cPickle as pkl
import particleSolver
reload(particleSolver)
from particleSolver import *
import mlabMacros as mlm

close('all')
if __name__== "__main__":
	runOn = True
	if runOn:
		#fix the solver parameters at something conservative
		alpha = .5
		thres = .0001
		maxSteps = 100000
		
		#we'll vary the temporal and spatial resolution for moderate diameter 
		#and a high injection velocity
		thermalTimestep = (linspace((1e-6),(1e-4), 20))
		kinematicTimestep = outer( thermalTimestep,1./arange(101.,0.,-4 ))
		thermalTimestep = tile(thermalTimestep[:,newaxis], (1,kinematicTimestep.shape[1]))
		skip = kinematicTimestep<1e-7
		
		
		#kinematicTimestep = thermalTimestep.copy()
		
		nNodes = 50
		d = 50e-6
		v0 = -18.

		#create the plasmaField, with cubic interpolation
		try:
			PF.density(array([0,0,0]))
		except NameError:
			PF = plasmaField(fastOn = True)
		
		cnt = 0
		
		# ############# highest rez case ###################
		
		#make the particle and shoot it into the plasma field
		PS = particleSolver(particle(diameter = d, 
														T = 300, 
														thermalTimestep = thermalTimestep[0,0], 
														kinematicTimestep = kinematicTimestep[0,0],
														nNodes = nNodes, 
														nStep = maxSteps), 
										 maxSteps = maxSteps)
		PS.p.inject(PF,x0 = array([0.001,.005,0.]), v0 = array([0.,v0,0.]))
		PS.p.splitScales = False
		start = time.time()
		PS.solve()
		elapsed = time.time()-start

		PS.plotSoln()
		
		xStar = PS.getTrajectory()[:-1]
		TStar = PS.getSoln()
		hStar = sum(PS.p.mass*PS.p.enthalpy())
		
		# ################################
		
		#track error in trajectory, temperature, and enthalpy
		m,n = (thermalTimestep.shape[0], kinematicTimestep.shape[1])
		ex =zeros((m,n))
		ey =zeros((m,n))
		et =zeros((m,n))
		eh = zeros((m,n))
		elapsed = zeros((m,n))
		ratio = zeros((m,n))
		X = [xStar]
		T = [TStar]
		
		matchLength = lambda  slave,master: (scipy.interpolate.interp1d(linspace(0,master.shape[0]-1, slave.shape[0]), transpose(slave))(arange(master.shape[0])))
		rg = lambda x_: range(x_.shape[0])
		
		
		cnt = 0
		N = thermalTimestep.size
		for i in rg(thermalTimestep):
			for j in range(kinematicTimestep.shape[1]):
				if skip[i,j]:
					cnt+=1
					continue
				print str(cnt)+'/'+str(N)+ '    '+ str(round(sum(elapsed)))+'s elapsed so far'
				PS = particleSolver(particle(diameter = d, 
														T = 300, 
														thermalTimestep = thermalTimestep[i,j], 
														kinematicTimestep = kinematicTimestep[i,j],
														nNodes = nNodes, 
														nStep = maxSteps), 
										 maxSteps = maxSteps)
				PS.p.inject(PF,x0 = array([0.001,.005,0.]), v0 = array([0.,v0,0.]))
				if i == j:
					PS.p.splitScales = False
				else:
					PS.p.splitScales = True
					
				# PS.p.thermalTimestep = thermalTimestep[i]
				# PS.p.kinematicTimestep = thermalTimestep[i]/round(thermalTimestep[i]/kinematicTimestep[j])
			
				start = time.time()
				PS.solve()
				elapsed[i,j] = time.time()-start
				#elapsed[j,i] = elapsed[i,j]
				#get the trajectory in space and the spatiotemporal temperature profile
				x = PS.getTrajectory()[:-1]
				t = PS.getSoln()
				
				#error in the trajectory
				ex[i,j] = sum(abs(xStar[0]-matchLength(x[0],xStar[0])))/xStar[0].shape[0]
				#ex[j,i] = ex[i,j]
				ey[i,j] = sum(abs(xStar[1]-matchLength(x[1],xStar[1])))/xStar[1].shape[0]
				#ey[j,i] = ey[i,j]
				#error in the final temperature profile
				et[i,j] =sum(abs(TStar[-1] - matchLength(t[-1],TStar[-1])  ))/TStar[-1].shape[0]
				#et[j,i] = et[i,j]
				#error in total enthalpy
				eh[i,j] = ( hStar- sum(PS.p.mass*PS.p.enthalpy()) )/hStar
				#eh[j,i] = eh[i,j]
				
				ratio[i,j] = double(PS.p.timestepRatio)
				#ratio[j,i] = ratio[i,j]
				
				X+=[x]
				T+=[t]
				cnt+=1
		
		elapsed[elapsed==0] = np.max(elapsed)
		theFile = open('particleConvergence.pkl','wb')
		D = dict(zip(['d', 'v0', 'thermalTimestep', 'kinematicTimestep', 'X','T','ex', 'ey', 'et', 'eh','elapsed','ratio'],
							[d, v0, thermalTimestep, kinematicTimestep, X,T, ex, ey, et, eh, elapsed,ratio ]))
		pkl.dump(D,theFile, -1)
		theFile.close()
		print 'saved!'
	
	else:
		theFile = open('particleConvergenceSplit.pkl','rb')
		D = pkl.load(theFile)
		theFile.close()
		
	cf = lambda x_: contourf(D['kinematicTimestep'], D['thermalTimestep'], x_, linspace(np.min(x_), np.max(x_),30))
	
	fig = lambda: figure(figsize = (20,22))
	
	mf = lambda x_, zLab_, fName_: mlm.mesh3(D['thermalTimestep'], 
																			D['kinematicTimestep'], 	
																			x_, f = mlm.fig(fName_), 
																			axeLab = ['thermal timestep, s', 'kinematic timestep, s', zLab_])
	
	mf(D['ex'], 'error, m', 'error in axial location')
	mf(D['ey'], 'error, m', 'error in vertical location')
	mf(D['et'], 'error, K', 'error in final mean temperature')
	mf(D['eh'], 'pct error', 'error in final total enthalpy')
	mf(log2(D['elapsed']/np.min(D['elapsed'])), 'log2(runtime/min(runtime))', 'simulation runtime, fastest run = '+str(np.min(D['elapsed'])))
	mf(D['ratio'], 'runtime, s', 'timestep ratio')
	
	# fig()
	# cf(D['ex'])
	# colorbar()
	# title('error in x(t)')
	# xlabel('kinematic timestep, s')
	# ylabel('thermal timestep, s')
	
	# fig()
	# cf(D['ey'])
	# colorbar()
	# title('error in y(t)')
	# xlabel('kinematic timestep, s')
	# ylabel('thermal timestep, s')
	
	# fig()
	# cf(D['et'])
	# colorbar()
	# title('error in final temperature profile')
	# xlabel('kinematic timestep, s')
	# ylabel('thermal timestep, s')
	
	# fig()
	# cf(D['eh'])
	# colorbar()
	# title('error in final enthalpy')
	# xlabel('kinematic timestep, s')
	# ylabel('thermal timestep, s')
	
	# fig()
	# cf(D['elapsed'])
	# colorbar()
	# title('elapsed time, s')
	# xlabel('kinematic timestep, s')
	# ylabel('thermal timestep, s')