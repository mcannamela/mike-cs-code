from particleSolver import particleSolver

PS = particleSolver(particle(diameter = array([1e-6,50e-6]), 
														T =300,
														thermalTimestep = 4e-5,
														kinematicTimestep = 2e-5,
														nNodes = 50, 
														nStep = maxSteps))
#create the plasmaField
try:
	PF.density(array([0,0,0]))
except NameError:
	PF = plasmaField(fastOn = True)
		
PS.p.inject(PF,x0 = array([0.001,.005,0.]), v0 = array([0.,-10.,0.]))

nTrials = 20
for i in range(nTrials):
	PS.reset()
	PS.solve()
