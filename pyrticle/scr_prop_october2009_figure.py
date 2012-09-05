import time
start = time.time()
import particleSolver as ps
reload(ps)
import particle as part
import gasMixture as gm
reload(part)
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
nV = 15
nD = 5

v0 = linspace(3,18, nV)
d = linspace(25e-6,75e-6, nD)
G = meshgrid(v0,d)

plotIdx = (atleast_1d(int32(nV/2)), int32(linspace(0,nD-1,3)))

nNodes = 20

#create the plasmaField
sim =  "cur-700_flo-50_gas-argon68helium32_make-SG100"
try:
    PF.density(array([0,0,0]))
except NameError:
    caselabel = 'I=500_V=60_Ar=40_H2=10'
    PF = part.lavaField(simName = 'lavaJet '+caselabel+'.pkl', gasPickle = 'lavaJet_gas '+caselabel+'.pkl')
# PF = part.plasmaField(simName = sim, 
    # fastOn = True, temperatureScaling = 1.15)
    
P = []
g = gm.gasMixture('zirconiaGauss')
injectionPosition = array([.001, .006, 0])
for i in range(G[0].flatten().shape[0]):
    p = part.particle(diameter = G[1].flatten()[i], 
                    T =300,
                    thermalTimestep = 4e-5,
                    kinematicTimestep = 2e-5,
                    nNodes = nNodes,
                    material = g)
    p.inject(PF, injectionPosition,array([0, -G[0].flatten()[i], 0]))
    
    P+=[p]
P = array(P).reshape(G[0].shape)

nTraj = 50
T = zeros( G[0].shape+ (nNodes,) )
r = zeros_like(T)
X = zeros( (2,)+G[0].shape+ (nTraj,) )
h = zeros_like(G[0])
y = zeros_like(h)
mi = zeros_like(h)

print 'init time: ' +str(time.time()-start)

# PS = ps.threadedSolve(P,ps.particleSolver)
# print 'done solving'
for i in range(G[0].flatten().shape[0]):
    idx = np.unravel_index(i,G[0].shape)
    
    PS = ps.particleSolver(P[idx])
    
    PS.solve()
    print 'solved! '+str(i) +' of ' +str(G[0].ravel().shape[0])

    T[idx] = PS.getSoln()[-1]
    r[idx] = PS.getSolnCoords()[1][-1]
    h[idx] = PS.getEnthalpy()
    y[idx] = PS.p.verticalPosition()
    mi[idx] = PS.p.totalMass/PS.p.initialMass
    x = PS.getTrajectory()[:2]
    s = hstack(  (  zeros(1)  , cumsum(sqrt(sum(diff(x,axis = 1)**2, axis = 0)))  )  )
    s/=s[-1]
    X[0][idx]= scipy.interp( linspace(0,1,nTraj),s, x[0])
    X[1][idx]= scipy.interp( linspace(0,1,nTraj),s, x[1])

        
close('all')




#plot the trajectories
#fig(1)
PF.showField(velocity = False)
M = plotIdx[0].shape[0]
N = plotIdx[1].shape[0]
sty = ['.', '*','o']
for j in range(M):#injection velocity
    for i in range(N):#diameter
        J = plotIdx[0][j]
        I = plotIdx[1][i]
        
        plot(X[0][I,J], X[1][I,J],sty[i]+'-', 
            color = [0,0,0], 
            linewidth = 2,alpha = .5, 
            label = '(v, D) = ('+str(round(G[0][I,J]))+', '+ str(round(G[1][I,J]*1e6))+')')
axis((0,.1, -.025,.025))
axisFontSize()
    
xlab('z, m')
ylab('r, m')
tit('Temperature Field with Trajectories\ninjection velocity in m/s, diameters in microns ')
legend()
sf('proposal_trajectoryFig')


#plot vertical travel vs injection velocity
fig(2)
fig(3)
mu = v0[plotIdx[0][0]]
sig = 2
for i in range(N):
    I = plotIdx[1][i]
    m = (y[I,-1]-y[I,0])/(v0[-1]-v0[0])
    mu_mapped = y[I,plotIdx[0][0]]
    sig_mapped = abs(m*sig)
    supp = linspace(mu_mapped-3*sig_mapped, mu_mapped+3*sig_mapped, 100)
    
    figure(2)
    plot(v0, y[I,:], linewidth = 2, label = 'D = '+str(round(G[1][I,0]*1e6)))
    
    figure(3)
    plot(supp, scipy.stats.norm.pdf(supp, loc = mu_mapped, scale = sig_mapped),linewidth = 2, label ='D = '+str(round(G[1][I,0]*1e6)) )
    
figure(2)
axisFontSize()    
xlab('injection velocity, m/s')
ylab('final vertical location, m')
tit('Spread of Trajectories for Different Size Particles\ndiameters in microns ')
legend()
sf('proposal_spreadFig1')

figure(3)
axisFontSize()    
xlab('final vertical location, m')
ylab('probability density')
tit('Posterior Distributions of Final Locations \nfor Particles of Increasing Diameter')
legend(loc = 'upper left')
sf('proposal_spreadFig2')

fig(4)
contourf(G[1]*1e6, G[0], mi, 30) 
axisFontSize()    
xlab('d, microns')
ylab('v0, m/s')
tit('ratio of final to initial mass')

mlm.surf3(G[1]*1e6,G[0], T[...,-1], f = mlm.fig('final surface temperature'))

D = dict(zip(['v0','d','G','T'],[v0,d,G,T]))
f = open('T_v0'+str(v0[0])+'-'+str(v0[-1])+'_d'+str(d[0]*1e6)+'-'+str(d[-1]*1e6)+'.pkl','wb')
pkl.dump(D,f,-1)
f.close()


