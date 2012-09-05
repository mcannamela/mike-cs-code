from numpy import *
v0 = linspace(3,20, 10)
d = array([.1e-6, 45e-6])
nNodes = 50

import time
start = time.time()
import particleSolver as ps
reload(ps)
import particle as part
reload(part)
import mlabMacros as mlm
import pdb
print 'import time: ' +str(time.time()-start)

xlab = lambda s_: xlabel(s_, fontsize = 24)
ylab = lambda s_: ylabel(s_, fontsize = 24)
tit = lambda s_: title(s_, fontsize = 20)
sf = lambda figName_:  savefig(figName_+'.png', format = 'png', transparent = True)
fig = lambda:  figure(figsize = (20,22))

start = time.time()
p = part.particle(diameter = d, 
														T =300,
														thermalTimestep = 4e-5,
														kinematicTimestep = 2e-5,
														nNodes = nNodes)

#create the plasmaField
sim =  "cur-700_flo-50_gas-argon68helium32_make-SG100"
try:
	PF.density(array([0,0,0]))
except NameError:
	PF = part.plasmaField(simName = sim, fastOn = True)
# PF = part.plasmaField(simName = sim, 
	# fastOn = True)
p.inject(PF, array([.0001, .0035, 0]),array([0, -v0[0], 0]))
PS = ps.particleSolver(p)
T = zeros( (v0.shape[0],nNodes))
r = zeros_like(T)
h = zeros(v0.shape[0])
t = zeros_like(h)
v = zeros_like(h)
X = []

print 'init time: ' +str(time.time()-start)
for i in range(v0.shape[0]):
	PS.p = p
	PS.reset(v0 = array([0, -v0[i], 0]))
	PS.solve()
	print 'solved!'
	try:
		T[i] = PS.getSoln()[-1]
		r[i] = PS.getSolnCoords()[1][-1]
		h[i] = PS.getEnthalpy()
		t[i] = sum(PS.p.thermalTimestepHistory)
		X+=[PS.getTrajectory()]
		v[i] = sqrt(sum(PS.p.velocity()**2))
	except:
		pdb.set_trace()	

# f = mlm.fig('Temperature Profiles for Different Injection Velocities')
# mlm.mesh3(v0[:,newaxis]*ones_like(r), r, T, f = f, axeLab = ['injection velocity, m/s', 'radius, um', 'temperature, K'])
close('all')
#fig()
PF.showField(velocity = False)
N = double(v0.shape[0])
for i in range(v0.shape[0]):
	plot(X[i][0], X[i][1],'.-', color = [(N-i)/N, (N-i)/N, (N-i)/N], linewidth = 2,alpha = .5, label = str(v0[i]))
fontsize = 20
ax = gca()
for tick in ax.xaxis.get_major_ticks():
	tick.label1.set_fontsize(fontsize)
for tick in ax.yaxis.get_major_ticks():
	tick.label1.set_fontsize(fontsize)
xlab('z, m')
ylab('r, m')
tit('Temperature Field with Trajectories\nwhite is 3m/s, black is 20m/s')
sf('trajectoryFig')

fig()
plot(v0, T[:,-1],'ko-', linewidth = 2)
xlabel('injection velocity, m/s')
ylabel('temperature, K')
title('surface temperature vs injection velocity')

fig()
plot(v0, h,'ko-', linewidth = 2)
xlabel('injection velocity, m/s')
ylabel('enthalpy, J/kg')
title('net enthalpy aquired vs injection velocity')

fig()
plot(v0,t,'ko-', linewidth = 2)
xlabel('injection velocity, m/s')
ylabel('travel time, s')
title('travel time vs injection velocity')

fig()
plot(v0, v, 'ko-', linewidth = 2)
xlabel('injection velocity, m/s')
ylabel('final speed, m/s')
title('final speed vs injection velocity')