

from numpy import *
seterr(invalid='raise')
seterr(divide='raise')
import time
import particleSolver as ps
import particle as part
import mlabMacros as mlm
import pdb
from plasmaField import *
from pylab import *
from scipy import integrate

import cPickle as pkl
import warnings
warnings.simplefilter('error')


xlab = lambda s_: xlabel(s_, fontsize = 24)
ylab = lambda s_: ylabel(s_, fontsize = 24)
tit = lambda s_: title(s_, fontsize = 20)
sf = lambda figName_:  savefig(figName_+'.png', format = 'png', transparent = True)
fig = lambda:  figure(figsize = (20,22))

def plotOverSpace(X, F, skip):
    N = double(len(X[0]))
    for i in range(len(X[0])):
        if mod(i,skip)==0:        
            plot(F[i], '.-', color = [(N-i)/N, (N-i)/N, (N-i)/N], 
                 linewidth = 2,alpha = .5, label = str(i))
    


start = time.time()

#create the plasmaField
#with open('/media/raidArray/CODE/jetModel/jet_sg100_ArPilot_40_500_DGMesh.pkl') as f:
with open('/media/raidArray/CODE/jetModel/jet_arPilot_40_500_highDensity.pkl') as f:
    PF = pkl.load(f)

N = 10
#v0 = linspace(4.1,14.3, N)
v0 = linspace(5.6,23.32, N)#Qc=6
t0 = .00135+linspace(0,3e-4, N)
theta = linspace(-15.,15.,N)*pi/180.
d = array([.1e-6, 35e-6])
nNodes = 20

p = part.timeDependentJetParticle(      diameter=d,
                        x0=float32([0.,0.,0.]),
                        v0=float32([0.,0.,0.]),
                        T=float32(300.), 
                        nNodes = nNodes,
                        thermalTimestep = float32(1.e-6),
                        kinematicTimestep = float32(1.e-6),
                        nStep = uint16(2000),
                        solidus_liquidus = float32([2912.,  3062.]),
                        vaporizationTemperature = float32(4600.),
                        heatOfVaporization = 5227000., 
                        vaporPressureFun = lambda T_: 101325.*exp(-4.5176e4/T_+10.088),
                        materialName = 'zirconiaGauss',
                        material = None,
                        plasmaField = PF,
                        t0=0.0)

p.adaptiveTimestepping=True
p.splitScales = 0

p.thermalStep   = float32(200)
p.inject(PF, array([0,  0.000,.004 ]),array([0, 0, -v0[0]]))
PS = ps.explicitParticleSolver(p, xMax = .11)
PS.alpha = .4
PS.stepErrorThres = 2e-5
PS.maxIter = 200
T = zeros( (v0.shape[0],nNodes))
r = zeros_like(T)
h = zeros(v0.shape[0])
timeHistory = []
v = zeros_like(h)
X = []
VMF = []
Tp = []
Ts = []
Tc = []
TsMinus = []
SS = []
R = []
xEnd = .1
print 'init time: ' +str(time.time()-start)

for i in range(v0.shape[0]):
    PS.p = p
#    PS.reset(x0 = array([0, .004, 0]), v0 = array([v0[i]*sin(10*pi/180), -v0[i]*cos(10*pi/180), 0  ]), t0 = .0011)
    PS.reset(x0 = array([0, .004, 0]), v0 = array([mean(v0)*sin(10*pi/180), -mean(v0)*cos(10*pi/180), 0  ]), t0 = t0[i])
#    PS.reset(x0 = array([0, .004*cos(theta[i]), .004*sin(theta[i])]), v0 = array([0, -mean(v0)*cos(theta[i]), -mean(v0)*sin(theta[i])  ]), t0 = .0011)
#    PS.reset(x0 = array([0, .004, 0]), v0 = array([0, -mean(v0)*cos(theta[i]), -mean(v0)*sin(theta[i])  ]), t0 = .0011)
    PS.solve()
    print 'solved!'

    r[i] = PS.getSolnCoords()[1][-1]
    R+=[PS.getSolnCoords()[1][:,-1]]
    h[i] = PS.getEnthalpy()
    timeHistory += [array(PS.p.t[:(PS.p.stepCount-1)]).copy()]
    X+=[PS.getTrajectory()]
    VMF+= [array(PS.p.vmfHistory)]
    Tp+= [array(PS.p.plasmaTemperatureHistory)]
    v[i] = sqrt(sum(PS.p.velocity()**2))
    SS+=[PS.getSoln()]
    S = PS.getSoln()
    
    Ts+=[S[:,-1]]
    Tc+=[S[:,0]]
    
    
    j = flatnonzero(X[-1][0]<xEnd)[-1]
    x = X[-1][0][j]
    try:
        xp = X[-1][0][j+1]
    except IndexError:
        T[i] = S[j]
        continue
        

    T[i] = S[j]+(xEnd-x)*(S[j+1]-S[j])/(xp-x)
    
   

#f = mlm.fig('Temperature Profiles for Different Injection Velocities')
#mlm.mesh3(v0[:,newaxis]*ones_like(r), r, T, f = f, 
           #axeLab = ['injection velocity, m/s', 'radius, um', 'temperature, K'])
#T08 = []
#for i, Ti in enumerate(T):
    
close('all')
PF.showField(velocity = False, planeNormal = 'z', t = t0[1], withColor = True)

N = double(v0.shape[0])
cu = cm.copper
for i in range(v0.shape[0]):
    plot(X[i][0], X[i][1],'.-', color = cu(double(i/N)),#'k',#[(N-i)/N, (N-i)/N, (N-i)/N], 
         linewidth = 2,alpha = .6, label = str(v0[i]))

xlab('x, m')
ylab('y, m')
tit('Temperature field with y trajectories')
xlim((0,.1))

PF.showField(velocity = False, planeNormal = 'y', t = t0[1], withColor = True)
for i in range(v0.shape[0]):
    plot(X[i][0], X[i][2],'.-', color = cu(double(i/N)),#'k',#[(N-i)/N, (N-i)/N, (N-i)/N], 
         linewidth = 2,alpha = .6, label = str(v0[i]))
xlab('x, m')
ylab('z, m')
tit('Temperature field with z trajectories')
xlim((0,.1))

#figure()
#for i in range(v0.shape[0]):
#    plot(X[i][0], VMF[i],'.-', color = [(N-i)/N, (N-i)/N, (N-i)/N], 
#         linewidth = 2,alpha = .5, label = str(v0[i]))
#xlab('x, m')
#ylab('vapor mass flux, kg/s')
#tit('vapor mass flux')

figure()
for i in range(v0.shape[0]):
    plot(X[i][0], Tp[i],':', color = cm.Blues(double(i/N)),#[(N-i)/N, (N-i)/N, (N-i)/N], 
         linewidth = 2,alpha = .8, label = str(v0[i]))
    plot(X[i][0], Ts[i],'--', color = cm.Reds(double(i/N)),#,[(N-i)/N, (N-i)/N, (N-i)/N], 
         linewidth = 2,alpha = .8, label = str(v0[i]))
    plot(X[i][0], Tc[i],'-', color = cm.Greens(double(i/N)),#[1-(N-i)/N, (N-i)/N, (N-i)/N], 
         linewidth = 2,alpha = .8, label = str(v0[i]))
xlab('x, m')
ylab('Temperature, K')
tit('Plasma and particle surface/center temperature histories')

figure()
for i in range(v0.shape[0]):
    plot(r[i], T[i],'.-', color = cu(double(i/N)),#[(N-i)/N, (N-i)/N, (N-i)/N], 
         linewidth = 2,alpha = 1, label = str(v0[i]))
xlab(r'radius, $\mu$m')
ylab('Temperature, K')
tit('Particle temperature profiles at x = 10 cm')



#figure()
#plot(v0, T[:,-1],'ko-', linewidth = 2)
#xlabel('injection velocity, m/s')
#ylabel('temperature, K')
#title('surface temperature vs injection velocity')
#
#figure()
#plot(v0, h,'ko-', linewidth = 2)
#xlabel('injection velocity, m/s')
#ylabel('enthalpy, J/kg')
#title('net enthalpy aquired vs injection velocity')
#
#figure()
#plot(v0,[t[-1] for t in timeHistory],'ko-', linewidth = 2)
#xlabel('injection velocity, m/s')
#ylabel('travel time, s')
#title('travel time vs injection velocity')
#
#figure()
#plot(v0, v, 'ko-', linewidth = 2)
#xlabel('injection velocity, m/s')
#ylabel('final speed, m/s')
#title('final speed vs injection velocity')

Ucl = mean(PF.U[18:120,0,:, PF.U.shape[3]/2, PF.U.shape[4]/2], axis = 0)
x = PF.X[0][:]*1000
gasTransitTime = integrate.cumtrapz(1./Ucl, x*.001)

figure()
for i in range(v0.shape[0]):
    plot(X[i][0]*1000, (timeHistory[i]-timeHistory[i][0])*1000,'--', color = cm.Reds(double(i+1)/(N+1)),         linewidth = 2,alpha = 1, label = str(v0[i]))
plot(x[:-1], gasTransitTime*1000, 'k-')
xlabel('axial distance, mm')
ylabel('transit time, ms')
title('transit time for particles (red dashed) \n and plasma gas (black solid)')
   
sly = slice(17,41)
Tin = PF.T[120:238,0, sly, sly]
TbarIn = mean(mean(Tin, axis = -1), axis = -1)
TmaxIn = amax(amax(Tin, axis = -1), axis = -1)
TbarP = array([mean((r[i]**3)*TT) for i,TT in enumerate(T)])
a = array([amin(TmaxIn), amax(TmaxIn)])
TbarP-= amin(TbarP)
TbarP*= (a[1]-a[0])/amax(TbarP)
TbarP+= a[0]
figure()
plot(PF.t[120:238]*1e6, TmaxIn, 'k-' )
for i,tt in enumerate(t0):
    plot(array([tt,tt])*1e6, a, 'r--')
    plot(tt*1e6, TbarP[i], 'ro')

xlabel(r'time, $\mu s$')
ylabel('max jet inlet temperature, K')
title('inlet max temperature shown with \n particle injection times marked')
