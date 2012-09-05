from scr_dumbHDFOutPlot import *
import os
import mlabMacros as mlm
from scipy import ndimage
import cPickle as pkl

def hdfoutAccumulate(dirname, nArrMax = None):
    X = []
    theTime = []
    cnt = 0
    for fname in os.walk(dirname).next()[2]:
        if cnt%100==0:
            print "completed %d of %d"%(cnt, len(os.walk(dirname).next()[2]))
        cnt+=1
        if cnt%4!=0:
            continue

        if fname[-3:]!='dat':
            continue
            
        

        (x,y,arrs, t) = HDFOutRead(os.path.join(dirname,fname), nArrMax)
        
        for i in range(arrs.shape[0]):
            arrs[i] = cookList[i](arrs[i])
        X+= [arrs.copy()]
        theTime+=[t]

        
    if all(array(theTime)==0):
        theTime = arange(len(theTime))
    else:
        theTime = array(theTime)
    print "all done"
    
    idx = argsort(theTime)
    return (x,y,theTime[idx],array(X)[idx,...])


if __name__=='__main__':
    dirname = 'E:\\labCODE\\LAVA\\550_67_42_6_1e-7'
    if not os.path.isdir(dirname):
        dirname = 'D:\\CODE\\LAVA\\startupRun'
    
    compute = True
    
    dumpMean = False
    if not os.path.isfile(os.path.join(dirname, os.path.split(dirname)[1]+'.pkl')):
        compute =True
        
    if compute:
        z,r,t,X = hdfoutAccumulate(dirname)
        with open(os.path.join(dirname, os.path.split(dirname)[1]+'.pkl'),'wb') as f:
            pkl.dump((z,r,t,X), f, -1)
    else:
        with open(os.path.join(dirname, os.path.split(dirname)[1]+'.pkl'),'rb') as f:
            z,r,t,X = pkl.load(f)
     
    if dumpMean:
        n = (300,200)
        nCut = (60,40)
        Y = mean(X, axis = 0)
        Yreg = zeros((Y.shape[0],)+n)
        for i in range(Yreg.shape[0]):
            zReg, rReg, Yreg[i] = regularize2d(z,r,Y[i],n,nCut)
            if i<3:
                mlm.surf3(x = zReg, y = rReg, z = Yreg[i],axeLab = ['axial, cm','radial, cm',labs[i]], f = mlm.fig(labs[i]))
        dumplist = [Yreg[i] for i in range(5)]+[Yreg[5]+Yreg[6]]
        with open('e:\\labcode\\pyrticle\\lavaJet.pkl','wb') as f:
            pkl.dump(dict(zip(['radialCoord', 'axialCoord']+labs[:5]+['airMols'],                      [cookDist(rReg),cookDist(zReg)]+dumplist )),f, -1)
         
    Tcl = ndimage.gaussian_filter(squeeze(X[:,0,:,0,...]), 1)
    Vcl = ndimage.gaussian_filter(squeeze(X[:,1,:,0,...]),1)
#    
#    ONratio =squeeze(X[:, 6,...]/X[:,5,...])
#    ONratio[logical_or(isinf(ONratio), isnan(ONratio))] = .25
#    ONratio[ONratio>.5] = .5
#    ONratio[ONratio<.05] = .05
# 
#    AirRatio = sum(X[:, 5:,...], axis = 1)/sum(X[:,3:,...], axis = 1)
#    AirRatio[logical_or(isinf(AirRatio), isnan(AirRatio))] = 0
# 
#    hist(log(ONratio.ravel()), 200)
#    xlabel('ox/ni ratio')
#    figure()
#    hist(AirRatio.ravel(), 200)
#    xlabel('air ratio')
#         
#    mlm.mesh3(x =x, y = y, z = ONratio[0],axeLab = ['radial, cm','axial, cm','oxygen/nitrogen ratio'], f = mlm.fig('air composition'))
#    mlm.mesh3(x =x, y = y, z = AirRatio[0],axeLab = ['radial, cm','axial, cm','air mole fraction'], f = mlm.fig('air'))
    mlm.mesh3(x =t, y = z, z = Tcl,axeLab = ['time, s','axial, cm',labs[0]], f = mlm.fig(labs[0]))
    mlm.mesh3(x =t, y = z, z = Vcl,axeLab = ['time, s','axial, cm',labs[1]], f = mlm.fig(labs[1]))
#    mlm.ml.show()
    
    