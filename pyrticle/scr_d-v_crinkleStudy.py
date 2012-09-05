from numpy import *
import numpy as np
from plotMacros import *
import mlabMacros as mlm
import cPickle as pkl
import particleModelEnsemble as pmEns
from plasmaField import lavaField
import pdb
import os
from scipy.ndimage import gaussian_filter1d as smooth
from scipy import interp
from ridgeFinder import blobFinder2d, basicRidgeNetwork, scaleSpace2d, neigh3d_26connected
from helpers import matchImageSize

RF, RN, SS = blobFinder2d, basicRidgeNetwork, scaleSpace2d  


def detectSeam(x, s = (.25, .25)):
    mx0 = zeros_like(x)
    mx1 = zeros_like(x)
    mn0 = zeros_like(x)
    mn1 = zeros_like(x)
    y = smooth(x,s[1], axis = 1)
    y = smooth(y,s[0], axis = 0)
    
    for i in range(1,x.shape[1]-1):
        gtm = y[:,i]>y[:,i-1] 
        gtp = y[:,i]>y[:,i+1]
        
        ltm = y[:,i]<y[:,i-1] 
        ltp =  y[:,i]<y[:,i+1]
        
        mx0[:,i] = logical_and(gtm, gtp)
        mn0[:,i] = logical_and(ltm, ltp)
        
    for i in range(1,x.shape[0]-1):
        gtm = y[i,:]>y[i-1,:] 
        gtp = y[i,:]>y[i+1,:]
        
        ltm = y[i,:]<y[i-1,:] 
        ltp = y[i,:]<y[i+1,:]
        
        mx1[i,:] = logical_and(gtm, gtp)
        mn1[i,:] = logical_and(ltm, ltp)
        
    mn = logical_or(mn0, mn1)
    mx = logical_or(mx0, mx1)
    return (mn, mx)
    



##########################    
##select the case to study#####
##########################
close('all')
#name = 'lavaJet_9MB_I=600_Ar=38.5_H=5.5_V=64.5_.pkl'
name = 'sharpTemp_600_38.pkl'
with open(name, 'rb') as f: 
            D=pkl.load(f)

##########################    
##find the half width#####
##########################
hw = argmin(abs(D['axialVelocity']-.5*D['axialVelocity'][:,0][:,newaxis]), axis = 1)
hwT = argmin(abs(D['temperature']-.5*D['temperature'][:,0][:,newaxis]), axis = 1)
D['axialCoord']*=1e-2
D['radialCoord']*=1e-2

##########################    
##solve the ensemble, if need be#####
##########################
#PE = pmEns.basicParticleEnsemble(pmEns.basicGriddedParticleParameterGenerator, name+'Ens')
#
#PE.generate(dMinMax = [3e-6,125e-6],
#            vMinMax = [1, 35],
#            x0 =array([0.00635,.00792,0.]),
#            n = (100,99))
#            
#PE.ensembleSolve(plasmaField = 
#    lavaField(simName = name, simDict = D, 
#              fastOn = False, gasPickle = None))
# LF = lavaField(simName = name, simDict = D, fastOn = False, gasPickle = None)
# LF.showField()
##########################    
##bounce ensemble to disk#####
##########################
#with open('crinkle-SharpT.pkl', 'wb') as f:
##with open('crinkle-6.pkl', 'wb') as f:
## with open('crinkleWide.pkl', 'wb') as f:
#   pkl.dump(PE, f, -1)


##########################    
##load a saved ensemble#####
##########################
with open('crinkle-SharpT.pkl', 'rb') as f:
#with open('crinkle-6.pkl', 'rb') as f:
#with open('crinkleWide.pkl', 'rb') as f:
    PE = pkl.load(f)

##########################    
##summary for ensemble, plot if desired#####
##########################
S = PE.summary(plotOn = True)


##########################    
##now we want to locate the ridges#####
##########################
#nscales = 75
#maxFeatureSize = .3
#im = matchImageSize(zeros((128, 128)),S['surfaceTemperature'].T)
#
#ss = SS(im, 
#            nscales = nscales, maxFeatureSize = maxFeatureSize )
#rf = RF(ss)
#mn = rf('dark')
#mx = rf('bright')
#
#print 'ridges found!'
#fig()
#subplot(1,3,1)
#imshow(sum(mn, axis = 0))
#title('mins')
#subplot(1,3,2)
#imshow(sum(mx, axis = 0), alpha = .25)
#title('maxs')
#subplot(1,3,3)
#imshow(im)
#
#pdb.set_trace()
#rnb = RN(rf, kind = 'bright', neighFun = neigh3d_26connected, fname = None)
#rnd = RN(rf, kind = 'dark', neighFun = neigh3d_26connected, fname = None)
mn, mx = detectSeam(S['surfaceEnthalpy'], s = (2., 2.))
G = float64(mgrid[:mn.shape[0],:mn.shape[1]])
g = float64(mgrid[:mn.shape[0],:mn.shape[1]])
g[0] = float64(mn.shape[0] - g[0])
g[1] = float64(mn.shape[1] - g[1])
#R = ((mn.shape[0] -g[0])**2+ ( mn.shape[1]-g[1])**2)**.5
R = ((g[0])**2+ ( g[1])**2)**.5

mnThres = [88,96]
mxThres1 = [90, 107]
mxThres2 = [75, 83]
cutoff = 65

sidx = lambda thres, cut: logical_and(logical_and(R>thres[0],R<thres[1]),
                                          G[1]<cut)

SI = [flatnonzero(mx[sidx(mxThres1, cutoff)]),
    flatnonzero(mn[sidx(mnThres, 90)]),
      flatnonzero(mx[sidx(mxThres2, cutoff)]),
    ]
#pdb.set_trace()
fig()
subplot(1,3,1)
imshow(mn)
imshow(sidx(mnThres, 90), alpha = .25)
title('mins')
subplot(1,3,2)
imshow(mx)
imshow(sidx(mxThres1, cutoff), alpha = .25)
imshow(sidx(mxThres2, cutoff), alpha = .25)
title('maxs')
subplot(1,3,3)
imshow(S['surfaceEnthalpy'])
#
##show()
#
#mnmx = transpose(array(detectSeam(S['surfaceTemperature'].T)), (0,2,1))
#
#
##
#
##
d = PE.parameterVectors[...,1]
v = abs(PE.parameterVectors[..., -2])
z = D['axialCoord']
r = D['radialCoord']
##
###########################    
###extract quantities of interest from the ensemble#####
###########################
thermalStep = array(zeros(len(SI)), dtype = 'object')
surfTemp =array(zeros(len(SI)), dtype = 'object')
traj = array(zeros(len(SI)), dtype = 'object')
force = array(zeros(len(SI)), dtype = 'object')
bcTemp = array(zeros(len(SI)), dtype = 'object')
for i in range(len(SI)):
    funs = [lambda m: diff(m.PS.getSoln()[:,-1]),
            lambda m: m.PS.getSoln()[:-1,-1],
            lambda m: m.PS.getTrajectory()[:-1],
            lambda m: array(m.PS.p.forceHistory),
            lambda m: array(m.PS.bcHistory)]
    
    thermalStep[i],surfTemp[i], traj[i], force[i] , bcTemp[i]  = PE.arrayGet(
                                            SI[i], funs, flatten = True)
                                            

cutPlane = .02
cutFn = lambda m: interp(cutPlane, m.PS.getTrajectory()[:-1][0],m.PS.getTrajectory()[:-1][1])
yy = PE.arrayGet(getFunList = [cutFn])[0]

Y = [[] for T in traj]                                            
hy = [[] for T in traj]                                            
by = [[] for T in traj]                                            
for i in range(traj.shape[0]):
    for t in traj[i]:
        Y[i]+= [interp(cutPlane, t[0], t[1])]
    Y[i] = array(Y[i])
    hy[i], by[i] = histogram(Y[i], 15)
    by[i] = by[i][:-1]+diff(by[i])/2
        
#
###########################    
###plot the thermalStep size vs surface temperature
### and the force-plasma temeprature scatter#####
###########################

c = ['r', 'b', 'c']
mark = ['s', '.', 'x']
labs = ['bask', 'whisk', 'sweet']
##for i in range(len(thermalStep)):
##    for j in range(thermalStep[i].shape[0]):
##        figure(3)
##        plot(surfTemp[i][j], thermalStep[i][j], mark[i], color =c[i], alpha = .45)
##        figure(4)
##        plot(force[i][j], bcTemp[i][j], mark[i], color =c[i], alpha = .45)
##
##figure(3)
##xlab('surface temp, K')
##ylab('next thermal step, K')
##tit('Step size vs temperature \n r - basking, b - whisked, cy - sweet')
##axisFontSize()
##
##figure(4)
##xlab('force, N')
##ylab('plasma temperature, K')
##tit('scatter of heat and momentum forcing \n r - basking, b - whisked, cy - sweet')
##axisFontSize()
##
##########################    
##plot the trajectories on the background of the T-V field#####
##########################
fig()

subplot(1,2,1)
##derivatives of T and V fields##
dt = abs(diff(D['temperature'][15:], axis = 1).T/diff(r[:2]*1e3))
dv = abs(diff(D['axialVelocity'][15:], axis = 1).T/diff(r[:2]*1e3))
##################################

##contour of the derivative fields##
contourf(z[15:], r_[-r[-2:0:-1], r[:-1]],r_[dt[-1:0:-1], dt] , 30, 
                             cmap = cm.copper, linewidths = None)
contour(z[15:], r_[-r[-2:0:-1], r[:-1]],r_[dv[-1:0:-1], dv] ,15, 
                            cmap = cm.cool, linewidths = 2, alpha = .5)
##         ##                   
#plot(z, r[hw], 'y', linewidth = 3, label = 'half width V')
#plot(z, r[hwT],'g', linewidth = 3, label = 'half width T')
#legend()

#for i in range(len(traj)):
#    for j in range(traj[i].shape[0]):
#        plot(traj[i][j][0], traj[i][j][1], '.-', color = c[i], alpha = .2, ms = 8)

xlim([.0055,.03])
ylim([-.002,.006])
xlabel('z, m')
ylabel('r, m')
#tit('trajectories over temperature and velocity gradient fields \n'+ c[0]+' - basking, '+c[1]+' - whisked, '+c[2]+' - sweet\n'+r'filled contours - $\frac{dt}{dr}$, contours - $\frac{dv}{dr}$' )
title('temperature and velocity gradient fields \n'+r'filled contours - $\frac{dt}{dr}$, contours - $\frac{dv}{dr}$' )
axisFontSize()

subplot(1,2,2)
contourf(d[:,0]*1e6,v[0],S['surfaceEnthalpy'].T, 30)
CS = contour(d[:,0]*1e6,v[0],yy.T*1e3, 20, cmap = cm.cool)
clabel(CS, CS.levels, inline=False,fmt = '%.1f' , fontsize=12)

xlabel(r'diameter, $\mu$m')
ylabel('injection velocity, m/s')
title('particle surface enthalpy and \n vertical location at axial cutplane z = %.1f cm \n'%(cutPlane*100)+'filled contours - surface enthalpy \n contours - vertical location' )
subplots_adjust(top = .85, left = .06, right = .94)
colorbar()
#for i in range(len(Y)):
#    barh(by[i], hy[i], height = ones_like(by[i])*diff(by[i][:2]), 
#         align = 'center', alpha = .3, facecolor = c[i])

#ylim([-.002,.006])

fig()
plot(z, r[hw], 'y', linewidth = 3, label = 'half width V')
plot(z, r[hwT],'g', linewidth = 3, label = 'half width T')
title('Thermal and Momentum Jet Spreading' )
legend()
axisFontSize()
