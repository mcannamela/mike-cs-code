from diameterInversionTools import *
import time
import cPickle as pkl

def getMapIdx(p, N):
    return unravel_index(argmax(p.ravel()), N)

def getMeanIdx(x,y,p):
    return r_[sum(x[..., newaxis]*p), sum(y*p)]

def writeTable(fname, Dobs, Bobs, X):
    with open(fname, 'w') as f:
        f.write('Dobs \n Bobs \n Dfac \n')
        [f.write('%.3e\t'%d) for d in tuple(Dobs)]
        f.write('\n')
        [f.write('%.3e\t'%d) for d in tuple(Bobs)]
    
        for i in range(X.shape[0]):
            f.write('\n')
            for j in range(X.shape[1]):
                f.write('%.3e\t'%X[i,j])
    
    
#get input for blur model from measurements
start = array([3065, 3099, 3135, 3108, 3085])
stop = array([3242, 3213, 3194, 3233, 3279])
z = linspace(-2, 2, 5)

###gaussian##
##assume we measured 4*stdDev, convert to 6*stdDev
#W = 6.*(stop - start)/4.


##lorentz##
#assume we measured 6*stdDev
W = (stop - start)



#plumeWidth/6 in mm
sigz = 4.

#file where measured diameter and blur number are tabulated
#theFile = 'D_B_gauss_75.txt'
theFile = 'D_B_lorenz_75_GiantGrid.txt'
#theFile = 'D_B_lorenz_75_bigGrid.txt'
pklname = 'lorenzBigLUT.pkl'

#init blur objects and forward model
blur = blurModel(z,W, sigz)
fm = forwardMap(theFile, sigs = r_[1.,  .02 ])

#mesh the observations space
Bobs = linspace(0, .55, 149)
Dobs = linspace(5., 75,100)
#Bobs = linspace(0, .3, 29)
#Dobs = linspace(2., 50,30)

#allocate for diameter and intensity correctors
Dfac = zeros((Dobs.shape[0], Bobs.shape[0]))
Ifac = zeros_like(Dfac)

DfacMean = zeros((Dobs.shape[0], Bobs.shape[0]))
IfacMean = zeros_like(Dfac)

#allocate for the max indices and mean D, w values


#dimensions of the input tables
N = (fm.D.shape[0], fm.w.shape[0])


#begin computation of the tables
start = time.time()
cnt = 0
for i in range(Dobs.shape[0]):
    if mod(i,10)==0:
        print 'completed %d of %d in %.2e s'%(cnt,Dobs.shape[0]*Bobs.shape[0], time.time()-start )
    
    for j in range(Bobs.shape[0]):
        p = fm.modelProbability(r_[Dobs[i], Bobs[j]], pFuns = [None, blur.pdf])
        
        idx = getMapIdx(p, N)
        mu = getMeanIdx(fm.D, fm.w, p)
        
        Dfac[i,j] = fm.D[idx[0]]#/Dobs[i]
        Ifac[i,j] = fm.I[idx]/fm.Imeas[idx]
        
        DfacMean[i,j] = mu[0]#/fm.DmeasFun(mu)
        IfacMean[i,j] = fm.IFun(mu)/fm.ImeasFun(mu)
        
        cnt+=1

elapsed  = time.time()-start
print 'time elapsed: %.2e s'%elapsed

keys = ['Bobs', 'Dobs', 'Dfac', 'Ifac', 'DfacMean', 'IfacMean']
vals = [Bobs, Dobs, Dfac, Ifac, DfacMean, IfacMean]
pkldict = dict(zip(keys, vals))
with open(pklname, 'wb') as f:
    pkl.dump(pkldict, f, -1)

close('all')

fig()
plot(fm.Dmeas.ravel(), fm.B.ravel(),'.', alpha = .2)
ylab(r'B$_{obs}$')
xlab(r'D$_{obs}$, pix')

fig()
subplot(1,2,1)
contourf(Bobs, Dobs, Dfac, 30, cmap = cm.Accent); colorbar()
xlab(r'B$_{obs}$')
ylab(r'D$_{obs}$, pix')
tit(r'diameter factor D$^*$/D$_{obs}$')
subplot(1,2,2)
contourf(Bobs, Dobs, DfacMean, 30, cmap = cm.Accent); colorbar()
xlab(r'B$_{obs}$')
ylab(r'D$_{obs}$')
tit('diameter factor from mean \n'+r'D$^*$/D$_{obs}$')

fig()
subplot(1,2,1)
contourf(Bobs, Dobs, Ifac, 30, cmap = cm.Accent); colorbar()
xlab(r'B$_{obs}$')
ylab(r'D$_{obs}$')
tit('Intensity factor from max \n'+ r'I$^*$/I$_{obs}$')

subplot(1,2,2)
contourf(Bobs, Dobs, IfacMean, 30, cmap = cm.Accent); colorbar()
xlab(r'B$_{obs}$')
ylab(r'D$_{obs}$')
tit('Intensity factor from mean \n'+r'I$^*$/I$_{obs}$')



print 'writing tables to disk...'
writeTable('diameterCorrector.lut', Dobs, Bobs, Dfac/Dobs[...,newaxis])
writeTable('intensityCorrector.lut', Dobs, Bobs, Ifac)
writeTable('diameterCorrectorMean.lut', Dobs, Bobs, DfacMean/fm.DmeasFun(DfacMean))
writeTable('intensityCorrectorMean.lut', Dobs, Bobs, IfacMean)



print 'done!'
        
#show()

mlm.surf3(fm.D,fm.w,  fm.B, f = mlm.fig('blur #'))
mlm.surf3(fm.D, fm.w, fm.Dmeas, f = mlm.fig('D_meas'))
mlm.surf3(fm.D, fm.w, fm.I/fm.Imeas, f = mlm.fig('I_meas'))
mlm.surf3(fm.D, fm.w, fm.I, f = mlm.fig('I'))

mlm.surf3(Dobs, Bobs, Dfac, f = mlm.fig('Dfac'))
mlm.surf3(Dobs, Bobs, DfacMean, f = mlm.fig('DfacMean'))
mlm.surf3(Dobs, Bobs, Ifac, f = mlm.fig('Ifac'))
mlm.surf3(Dobs, Bobs, IfacMean, f = mlm.fig('IfacMean'))



mlm.ml.show()