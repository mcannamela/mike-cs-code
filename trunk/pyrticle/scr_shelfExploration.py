import time
start = time.time()
import particleSolver as ps
#reload(ps)
import particle as part
import gasMixture as gm
#reload(part)
import mlabMacros as mlm
import pdb
import cPickle as pkl
print 'import time: ' +str(time.time()-start)

fsz = 20
xlab = lambda s_: xlabel(s_, fontsize = fsz)
ylab = lambda s_: ylabel(s_, fontsize = fsz)
tit = lambda s_: title(s_, fontsize = fsz)
sf = lambda figName_:  savefig(figName_+'.png', format = 'png', transparent = True)
fig = lambda n_:  figure(n_, figsize = (20,22))

def axisFontSize(ax = gca(), fsz = 20):
	ax = gca()
	for tick in ax.xaxis.get_major_ticks():
		tick.label1.set_fontsize(fsz)
	for tick in ax.yaxis.get_major_ticks():
		tick.label1.set_fontsize(fsz)

start = time.time()

#mesh initial velocity and diameter
v0 = 7
nD = 5

d = linspace(10e-6,100e-6, nD)

nNodes = 50

#create the plasmaField
sim =  "cur-700_flo-50_gas-argon68helium32_make-SG100"
try:
	PF.density(array([0,0,0]))
except NameError:
	PF = part.plasmaField(simName = sim, fastOn = True, temperatureScaling = 1.20)
# PF = part.plasmaField(simName = sim, 
	# fastOn = True, temperatureScaling = 1.20)
	
P = []
g = gm.gasMixture('zirconiaGauss')
injectionPosition = array([.0001, .0035, 0])

for i in range(d.shape[0]):
	p = part.particle(diameter = d[i], 
					T =300,
					thermalTimestep = 4e-5,
					kinematicTimestep = 2e-5,
					nNodes = nNodes,
					material = g)
	p.inject(PF, injectionPosition,array([0, -v0, 0]))
	
	P+=[p]


print 'init time: ' +str(time.time()-start)

PS = []
for p in P:
	PS+= [ps.particleSolver(p)]
	PS[-1].solve()
	print 'solved!'

close('all')

#plot the trajectories
fig(1)
PF.showField(velocity = False)
i = 0 
for ps in PS:#injection velocity
	x = ps.getTrajectory()
	plot(x[0], x[1],'.' ,
			color = [(nD-i)/nD, (nD-i)/nD, (nD-i)/nD], 
			linewidth = 2,alpha = .5, 
			label = 'd = '+str(ps.p.diameter()*1e6)+' microns')
	i+=1
axis((0,.1, -.025,.025))
axisFontSize()
	
xlab('z, m')
ylab('r, m')
tit('Temperature Field with Trajectories\ninjection velocity = '+str(v0)+' m/s, diameters in microns ')
legend()
sf('shelfExplorationTrajectoryFig')


F = []
T = []
R = []
for ps in PS:
	di = ps.p.diameter()*1e6
	
	
	t = ps.getSoln()
	R+= [ps.getSolnCoords()[1][-1]]
	T+= [t[-1]]
	
	#uncomment to see histories
	#F+= [mlm.fig('Temperature History, '+str(di)+ 'microns')]
	#ml.clf(F[-1])
	#mlm.mesh3(ps.getSolnCoords()[0],ps.getSolnCoords()[1], t, f = F[-1])

fig(2)	
N = double(d.shape[0])
for i in range(N):
	#plot(R[i]/np.max(R[i]), T[i],linewidth = 2,color = [(N-i)/N, 0, i/N], label = 'd = ' +str(d[i]*1e6)+ ' microns')
	
	plot(R[i]/np.max(R[i]), (T[i]-2950.)/75,linewidth = 2,color = [(N-i)/N, 0, i/N], label = 'd = ' +str(d[i]*1e6)+ ' microns')
	
# subplot(121)	
# axisFontSize()
legend(loc = 'upper left')
# xlab('dimensionless radius, r/R')
# ylab('Temperature, K')
# tit('Final Temperature Profiles for\n Different Diameter Particles')

# subplot(122)
axisFontSize()
legend(loc = 'lower left')
xlab('dimensionless radius, r/R')
ylab('molten fraction\n (T-T$_{solidus}$)/(T$_{solidus}$-T$_{liquidus}$)')
tit('Molten Fraction Profiles for\n Different Diameter Particles, v0 = '+str(v0)+ ' m/s')	
sf('shelfMoltenFractionProfiles')

