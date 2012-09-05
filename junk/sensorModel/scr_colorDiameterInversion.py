from diameterInversionTools import *
from scipy.stats import norm as normRand 
from scipy.optimize import brent 
import time
import pdb 
from plotMacros import *
import cPickle as pkl
from scipy.integrate import cumtrapz

#################### script settings  ########################
logarithmicOn = False
computeOn = False
N = 5000
stdThresLog = .25
stdThres = 6
noiseFactor = 1e-9

#file where measured diameter and blur number are tabulated
#theFile = 'D_B_lorenz_color.txt'
#theFile = 'D_B_lorenz_color_snap.txt'
theFile = 'D_B_lorenz_color_snap_extended.txt'
obsFile = 'fakeObsfile.txt'
pPkl = 'pPkl'
pklname = 'obsList_D.pkl'
##############################################################

def zPrior(z, sigz):
    return gauss(z,0,sigz)

def p2map(d,p):
#    f = interpolate.interp1d(d, -p, kind = 'cubic', bounds_error=False,fill_value = 0 )
#    try:
#        xmin, fval, iter, funcalls = brent(f, brack = (d[0], d[-1]), full_output = True)
#    except:
#        pdb.set_trace()
#    return xmin
    return d[argmax(p)]


def generateObservations(pd, mu, sig, n):
    """
    use the interpolant dictionary  of a 
    forwardDiameterMapper object to 
    pick n measurement-observation pairs for inversion
    by sampling from gauss-distributed measurements
    mu[0] - diameter mean, pixels
    sig[0] - diameter standard dev, pixels
    mu[1] - z mean, mm
    sig[1] = z std dev, mm
    """
    m = float64([normRand(mu[i], sig[i]).rvs(n) for i in range(2)])
    obs = [pd['Dmeas'][i](m) for i in range(3)]
    obs += [pd['B'][i](m) for i in range(3)]
    return (m, float64(obs))
    
fm = forwardDiameterMapper(theFile, cutoff = -1)


if logarithmicOn:
    Dref = 7.
    Bref = .15
    fm.model[0] = log(fm.model[0]/Dref)
    fm.meas[:3] = log(fm.meas[:3]/Dref)
    fm.meas[3:] = log(fm.meas[3:]/Bref)
    fm.refreshPolDict()
    fm.headerCallback()
    sigs = 3*r_[ones(3)*.05,ones(3)*.05 ]
else:
    sigs = r_[ones(3)*3.,ones(3)*.03 ]
    
idm = inverseDiameterMapper(fm,sigs )

##########################
####### fake data ########

if logarithmicOn:
    mu = [log(12./Dref), 0.]#[mean diameter, mean z location
    sig = [.1, 4.]#[diameter std dev, z std dev]
    pPkl+= '_muD_%d_sigD_%d_N_%d_logOn'%(mu[0], sig[0],N)
else:
    mu = [20., 0.]#[mean diameter, mean z location
    sig = [4., 4.]#[diameter std dev, z std dev]
    pPkl+= '_muD_%d_sigD_%d_N_%d'%(mu[0], sig[0],N)


models, observations = generateObservations(fm.polDict, mu, sig, N)
obsClean = observations.copy()
for i in range(observations.shape[0]):
    observations[i]+= normRand(0., sigs[i]).rvs(observations.shape[1])

##########################

p = zeros((fm.arrDict['D'].shape[0],fm.arrDict['z'].shape[0]))


if computeOn:
    start = time.time()
    cnt = 0
    P = []
    Pz = []
    for i in range(observations.shape[1]):
        if mod(cnt,100)==0 and cnt>0:
            print 'completed %d in %.2e s'%(cnt, time.time()-start )
#            figure()
#            contourf(fm.arrDict['z'], fm.arrDict['D'], p, 30, cmap = cm.Accent)
            
        p = idm.modelProbability(observations.T[i],pFuns = [None, lambda z_: zPrior(z_,sig[1])  ])
        #p = idm.modelProbability(observations.T[i],pFuns = [None, None ])
        Pz+= [trapz(p, fm.arrDict['D'], axis = 0)]
        P+= [trapz(p, fm.arrDict['z'], axis = 1)]
        cnt+=1
    with open(pPkl, 'wb') as f:
        pkl.dump((P,Pz,p),f,-1)
else:
    with open(pPkl, 'rb') as f:
        (P,Pz,p) = pkl.load(f)
    
P = float64(P)
Pz = float64(Pz)


pMap = array([p2map(fm.arrDict['D'],P[i]) for i in range(P.shape[0])])
pMean = trapz(P*fm.arrDict['D'], fm.arrDict['D'], axis = 1) 
pStd = trapz(P*(fm.arrDict['D']-pMean[...,newaxis])**2, fm.arrDict['D'], axis = 1)**.5
pSkew = trapz(P*(fm.arrDict['D']-pMean[...,newaxis])**3, fm.arrDict['D'], axis = 1)/pStd**3
if logarithmicOn:
    idx = pStd<stdThresLog
else:
    idx = pStd<stdThres
print 'discarding %d of %d due to wide spread'%(pStd.shape[0]-flatnonzero(idx).shape[0], N)

pzMean = trapz(Pz*fm.arrDict['z'], fm.arrDict['z'], axis = 1) 
hz, bz = histogram(pzMean, 50, normed = True)
z = linspace(mu[1]-3*sig[1], mu[1]+3*sig[1], 100)


h,b = histogram(pMap[idx], 50, normed = True)
hm,bm = histogram(pMean[idx], 50, normed = True)
d = linspace(mu[0]-3*sig[0], mu[0]+3*sig[0], 100)




figure()
hist(observations[0], 50, normed = True, facecolor = 'r', alpha = .35, label = 'noisy raw estimate')
hist(obsClean[0], 50, normed = True, facecolor = 'y', alpha = .35, label = 'clean raw estimate')
#plot(b[:-1]+diff(b)/2, h, 'b', drawstyle = 'steps',alpha = .35, label = 'map estimate')
plot(bm[:-1]+diff(bm)/2, hm, 'g', drawstyle = 'steps', label = 'mean estimate')
plot(d,gauss(d, mu[0], sig[0]), 'k', label = 'target')
plot(fm.arrDict['D'], sum(P[idx], axis = 0)/P[idx].shape[0], 'y', label = 'pdf disjunction estimate')
legend()
xlabel('diameter, pixels')
if logarithmicOn:
    title('reconstructed diameter distributions, logarithmic land')
else:
    title('reconstructed diameter distributions')


#if logarithmicOn:
#    figure()
#    efun = lambda x_:Dref*exp(x_)
#    hist(efun(observations[0]), 50, normed = True, facecolor = 'r', alpha = .35, label = 'noisy raw estimate')
#    hist(efun(obsClean[0]), 50, normed = True, facecolor = 'y', alpha = .35, label = 'clean raw estimate')
#    hist(efun(pMean), 50, normed = True, facecolor = 'g', alpha = .35, label = 'mean estimate')
#    F = r_[0,cumtrapz(sum(P[idx], axis = 0)/P[idx].shape[0], fm.arrDict['D'])]
#    pp = r_[0, diff(F)/diff(efun(fm.arrDict['D']))]
#    plot(efun(fm.arrDict['D']), pp/trapz(pp,efun(fm.arrDict['D'])), 'y', label = 'pdf disjunction estimate')
#    legend()
#    xlabel('diameter, pixels')
#    title('reconstructed diameter distributions')

figure()
plot(bz[:-1]+diff(bz)/2, hz, 'g', drawstyle = 'steps', label = 'mean estimate')
plot(z,gauss(z, mu[1], sig[1]), 'k', label = 'target')
plot(fm.arrDict['z'], sum(Pz, axis = 0)/Pz.shape[0], 'y', label = 'pdf disjunction estimate')
xlabel('z, mm')
title('reconstructed position distribution')

figure()
xlabel('standard error of diameter distributions')
hist(pStd, 50, normed = True)

figure()
xlabel('skewness of diameter distributions')
hist(pSkew, 50, normed = True)


figure()
plot(observations[0], observations[1], 'b.', alpha = .15)
xlabel(r'D$_r$')
xlabel(r'D$_g$')
title('scatter of measure diameters')

figure()
plot(observations[3], observations[4], 'b.', alpha = .15)
xlabel(r'B$_r$')
xlabel(r'B$_g$')
title('scatter of blur numbers')

show()

