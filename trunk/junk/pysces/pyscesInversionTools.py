from inversePysces import *
import shutil

defaultFileDic = {
    'pklName':'inverseSave.pkl',
    'dataFile':'timeseries.csv',
    'pscFile':'model.psc', 
    'speciesDictFile':'species.dic',
    'measuredSpeciesFile':'measured.txt',
    'referenceParameterFile':'referenceParameters.dic'
    }

class pyscesInversionTrace(object):
    def __init__(self, kRef, sRef, logPMax,tK, tS, tD, tLogP):
        self.K = tK.copy()
        self.S = tS.copy()
        self.logP = tLogP.copy()
        self.D = tD.copy()
        self.kRef = copy(kRef)
        self.sRef = copy(sRef)
        self.logPMax = copy(logPMax)
        
        self.svdFigNum = 100
        
    def __add__(self, other):
        self.convertRef(other.kRef, other.sRef)
        self.K = r_['1,2', self.K, other.K]
        self.S = r_['1,2', self.S, other.S]
        self.logP = r_['0,1', self.logP, other.logP]
        self.D = r_['0', self.D, other.D]
        
    def copy(self):    
        cop = pyscesInversionTrace(self.kRef, self.sRef, self.logPMax, 
                                   self.K, self.S, self.D, self.logP)
        return cop
        
    def convertRef(self, kRef, sRef):
        assert all(kRef>0) and all(sRef>0), "negative reference values detected"
        K = self.unLogK()
        S = self.unLogS()
        self.kRef = kRef
        self.sRef = sRef
        
        try:
            self.K = self.logK(K)
        except:
            pdb.set_trace()
        self.S = self.logS(S)
        
        
    def SVD(self):
        A = r_['0,2', self.K, self.S]        
        U, s, Vt = svd(A, full_matrices = False)
        
        self.modeshapes = U.T
        self.singularValues = s
        self.modalAmplitudes = Vt
        
        return (self.modeshapes.copy(), 
                self.singularValues.copy(), 
                self.modalAmplitudes.copy())
    
    def unLogK(self):
        return (2**self.K)*self.kRef[...,newaxis]
    def unLogS(self):
        return (2**self.S)*self.sRef[...,newaxis]
    
    def logK(self, K):
        return log2(K/self.kRef[...,newaxis])
        
    def logS(self, S):
        return log2(S/self.sRef[...,newaxis])
        
    def plotSVD(self):
        if 'singularValues' not in self.__dict__.keys():
            self.SVD()
        figure(self.svdFigNum)
        plot(self.singularValues, '.-')
        xlabel('index')
        title('Singular Values of Parameter Trace Matrix')
        
class pyscesInversionCase(object):
    def __init__(self, modelFolder, fileDic = defaultFileDic, nP = 2):
        
        self.nP = nP
        self.modelFolder = modelFolder
        self.fileDic = fileDic.copy()
        for k in self.fileDic.keys():
            self.fileDic[k] = os.path.join(modelFolder,self.fileDic[k])
          
        #\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
        #             initialize data
        #\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
        print "data initializing..."
        try:
            self.D = [ret_erk_experiment(self.fileDic['dataFile'],
                                     self.fileDic['speciesDictFile'])
                                     for i in range(self.nP)]    
        except IOError:
            print ('either datafile or species file not found, '
                    'perhaps you intend to load a case? continuing without data...')
        
        print "done..."
        #\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

    def __call__(self, N = 10, burn = 0, thin = 1,
                 logRateConstantStdDev = 1.5, 
                 logInitialConcentrationsStdDev = 1.5,
                 verbose = False):
                     
        self.logRateConstantStdDev = logRateConstantStdDev
        self.logInitialConcentrationsStdDev = logInitialConcentrationsStdDev
        
        
        #-------------------------------------------------------
        #             initialize the pysces models
        #-------------------------------------------------------
        print "initializing forward models..."
        self.FM = [pyscesForwardModel(self.fileDic['pscFile'], 
                     self.fileDic['speciesDictFile'], 
                     self.fileDic['measuredSpeciesFile'],
                     self.fileDic['referenceParameterFile'],
                     self.D[i].erkTime(),
                     self.D[i].constraintFunctionList()
                     ) for i in range(self.nP)]
        print "done."
        #-------------------------------------------------------
        
        print "will now initialize inverse model list..."        
        #initialize inverse model
        self.IMList = []
        
        for i in range(self.nP):
            print "initializing %d of %d"%(i+1,self.nP)
            self.IMList += [inversePyscesModel( self.FM[i], 
                                                self.D[i], 
                                                logRateConstantStdDev = self.logRateConstantStdDev, 
                                                logInitialConcentrationsStdDev = self.logInitialConcentrationsStdDev, 
                                                verbose = verbose)]
        print "list initialized, setting up the pool..."
        self.IMPool = inverseModelPool(self.IMList, self.fileDic['pklName'])
        print "done. now for the tough part..."
        
        #sample from the posterior    
        start = time.time()   
        print "here we go!"
        print 'the time is '+time.ctime()
        self.IMPool(N, burn = burn, thin = thin)
        elapsed = time.time()-start
        print 'sampled %d in %.2f s: %.2e samples/s'%(N*self.nP,elapsed,self.nP*N/elapsed  )
        
        tK = self.IMPool.tK
        tS = self.IMPool.tS
        tD = self.IMPool.tD
        tLogP = self.IMPool.tLogP
        logNormConst = self.IMPool.invModList[0].D_obs.log_normConst()
        
        gidx = flatnonzero(tLogP>-300)
        if gidx.size==0:
            gidx = flatnonzero(tLogP>median(tLogP))
        
        print "discarding %d of %d due to extreme bad fit"%(len(tLogP)-len(gidx), len(tLogP))
        
        tK = tK[:,gidx]
        tS = tS[:,gidx]
        tD = tD[gidx]
        tLogP = tLogP[gidx]
        
        self.T = pyscesInversionTrace(self.FM[0].k0, self.FM[0].s0, logNormConst,
                                          tK, tS, tD, tLogP )
                                          
        return self.T.copy()
    
    def trace(self):
        return self.T.copy()
    
    def updateReferenceParameters(self, kRef = None, sRef = None):
        if kRef == None or sRef == None:
            bidx = argmax(self.T.logP)
            kRef = self.T.unLogK()[:,bidx]
            sRef = self.T.unLogS()[:,bidx]
        self.T.convertRef(kRef, sRef)
        if any(isnan(kRef)) or any(isnan(self.T.kRef)):
            pdb.set_trace()
        self.writeRefParFile(os.path.join(os.curdir,'refPar.temp'), kRef, sRef)
        [FM.initReferenceValues(os.path.join(os.curdir,'refPar.temp')) for FM in self.FM]
    
    def bounceModel(self, modelPath):
                
        if not os.path.isdir(modelPath):
            os.mkdir(modelPath)
        for k in self.fileDic.keys():
            try:
                shutil.copy(self.fileDic[k], modelPath)
            except:
                pass
            
        refFile = os.path.join(modelPath, 
                               os.path.split(self.fileDic['referenceParameterFile'])[1])
        self.writeRefParFile(refFile, self.T.kRef, self.T.sRef)
        self.save(os.path.join(modelPath, os.path.split(self.pklName())[1]))
    
    def writeRefParFile(self, refFilePath, kRef, sRef):
        assert not any(isnan(kRef)) and not any(isnan(sRef)), "nans in the refPars, writing failed"
        with open(refFilePath, 'wb') as f:
            [f.write('k%d %.3e\n'%(i+1, k))  for i,k in enumerate(kRef)]
            [f.write('s%d %.3e\n'%(i+1, s))  for i,s in enumerate(sRef)]
        

    def save(self, pklName = None):
        if pklName ==None:
            pklName = self.pklName()
        with open(pklName, 'wb') as f:
            pkl.dump(self, f, -1)
    
    def load(self, pklName = None):
        if pklName ==None:
            pklName = self.pklName()
        with open(pklName, 'rb') as f:
            obj = pkl.load(f)
            for k in obj.__dict__.keys():
                self.__dict__[k] = obj.__dict__[k]
                
    def pklName(self, baseName = None):
        if baseName == None:
            baseName = self.fileDic['pklName']
        return baseName+'_case.pkl'
    
    def plotFuns(self):
        return ('plotMarginals plotModel plotIntersampleDistance plotCorr ' 
                'plotGoodness plotTraces plotFitCorr'.split())
    
    def plotListOfFuns(self, plotFunList):
        for pf in plotFunList:
            self.__getattribute__(pf)()            

    def plotMarginals(self):
        
    #////////////////////////////////////////////////////////
    #               posterior histograms
    #////////////////////////////////////////////////////////

        subN = (3,3)
        xK = linspace(-3*self.logRateConstantStdDev, 3*self.logRateConstantStdDev, 200)
        for i in range(self.T.K.shape[0]):
            idx = mod(i, prod(subN))        
            if idx==0:
                figure()
            subplot(subN[0],subN[1],idx+1)
            
            hist(self.T.K[i],50, label = 'posterior', normed = True)
            plot(xK, stats.norm.pdf(xK,loc = 0, scale = self.logRateConstantStdDev ), label = 'prior')
            title('k%d, ref %.1f'%(i+1, self.T.kRef[i]), fontsize = 16)
            xlabel(r'log(k/k$_{ref}$)', fontsize = 14)
            

        xS = linspace(-3*self.logInitialConcentrationsStdDev, 3*self.logInitialConcentrationsStdDev, 200)
        for i in range(self.T.S.shape[0]):
            idx = mod(i, prod(subN))        
            if idx==0:
                figure()
            subplot(subN[0],subN[1],idx+1)
            
            hist(self.T.S[i],50, label = 'posterior', normed = True)
            plot(xS, stats.norm.pdf(xS,loc = 0, scale = self.logInitialConcentrationsStdDev ), label = 'prior')
            title(r'$s$'+'%d, ref %.1e'%(i+1, self.T.sRef[i]), fontsize = 16)
            xlabel(r'log($s/s_{ref}$)', fontsize = 14)
     #////////////////////////////////////////////////////////   

    def plotModel(self):
        #-------------------------------------------------------
        #             random models
        #------------------------------------------------------
        didx = randint(0,self.T.D.shape[0], 9)
        print "randomly selected trace points are:"    
        print didx
        figure()
        for i, idx in enumerate(didx):
            subplot(3,3,i+1) 
            plot(self.D[0].erkTime(), self.D[0].val().T, 'k')
            plot(self.D[0].erkTime(), squeeze(self.T.D[idx]).T, 'y')
            title('log(P(D|m)) = %.2f'%self.T.logP[idx])
        #------------------------------------------------------

        #-------------------------------------------------------
        #             best model
        #------------------------------------------------------
        bidx = argmax(self.T.logP)
        figure()
        plot(self.D[0].erkTime(), self.D[0].val().T, 'k')
        plot(self.D[0].erkTime(), squeeze(self.T.D[bidx]).T, 'y')
        title('log(P(D|m)) = %.2f'%self.T.logP[bidx])
        self.IMPool.invModList[0].bounceReferenceParameters(self.T.K[:,bidx], self.T.S[:,bidx], 
                                self.fileDic['bestRefParameterFile'])
            #------------------------------------------------------    
    
    def plotIntersampleDistance(self):
        #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
        #:::::::::::::   intersample distances  :::::::::::::::::::::::::
        #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
        dk = diff(self.T.K[:], axis = 1)
        dkMag = sum(dk**2, axis = 0)**.5/dk.shape[0]**.5
        figure()
        subplot(1,2,1)
        plot(dkMag)
        xlabel('transition number')
        ylabel('normalized distance')
        title('change in position per transition for k vector\n normalized by unit change in every component, N**.5')
        subplot(1,2,2)
        try:
            hist(log10(dkMag[dkMag>0]), 75)
        except:
            print "hist of intersample distance failed!"
        xlabel('log10(distance)')
        title('per transtion log distance (unnormalized)\n%d of %d zero values discarded'%(sum(dkMag==0), len(dkMag)))
        #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::    
    
    def plotCorr(self):
        ########################################################    
        ################### correlation coeffs##################
        ########################################################
        figure()
        pcolor(arange(1.,self.T.K.shape[0]+1)-.5, 
               arange(1.,self.T.K.shape[0]+1)-.5,
               corrcoef(self.T.K[:])-eye(self.T.K.shape[0]))
        colorbar()
        xlabel('constant index')
        ylabel('constant index')
        title('correlation coefficient between rate constants')
        ########################################################
    
    def plotGoodness(self):
        #.......................................................
        #................goodness of fit histograms.............
        #.......................................................
        figure()
        hist(self.T.logP, 100)
        xlabel('log(P(D|m))')
        title('goodness of fit, max possible log(P(D|m)) = %.2e'%self.T.logPMax)
        print "best fitting model has goodness of fit %.2e of a possible %.2e"%(amax(self.T.logP),self.T.logPMax) 
        #.......................................................    
    def plotTraces(self):
        figure()
        plot(self.T.K.T)
        title('traces of rate constants')
    
    def plotFitCorr(self):
        figure()
        subplot(1,3,1)
        plot(mean(abs(self.T.K), axis = 0), self.T.logP, '.', alpha = .2)
        xlabel('mean abs rate constant')
        ylabel('goodness of fit')
        subplot(1,3,2)
        plot(median(abs(self.T.K), axis = 0), self.T.logP, '.', alpha = .2)
        xlabel('median abs rate constant')
        ylabel('goodness of fit')
        subplot(1,3,3)
        plot(amax(abs(self.T.K), axis = 0), self.T.logP, '.', alpha = .2)
        xlabel('max abs rate constant')
        ylabel('goodness of fit')
        
if __name__=="__main__":
    
    modelPath = '/media/ubuntuData/CODE/pysces/ret_erk_simple'
    if not os.path.isdir(modelPath): 
        print 'model path not found...'
        modelPath = 'D:\\CODE\\pysces\\ret_erk_simple'
        print modelPath +' is the new path'
        
    fileDic = {
    'pklName':'pyscesInversionToolsTest',
    'dataFile':'ret_erk_peak_timeseries_2.csv',
    'pscFile':'ret_erk_s.psc', 
    'speciesDictFile':'ret_erk_species.dic',
    'measuredSpeciesFile':'ret_erk_MeasuredSpecies.txt',
    'referenceParameterFile':'best_referenceParameters_pilot_2.dic',
    'bestRefParameterFile':'best_referenceParameters.dic'
    }

    
    PIC = pyscesInversionCase(modelPath, fileDic = fileDic, nP = 1)
    PIC(N = 100, burn = 0, thin = 1,
                 logRateConstantStdDev = .1, 
                 logInitialConcentrationsStdDev = .1,
                 verbose = True)
                 
    PIC.save()
    PIC.load()
                 
    plotFuns = ('plotMarginals plotModel plotIntersampleDistance plotCorr ' 
                'plotGoodness plotTraces plotFitCorr'.split())
    for pf in plotFuns:
        PIC.__getattribute__(pf)()