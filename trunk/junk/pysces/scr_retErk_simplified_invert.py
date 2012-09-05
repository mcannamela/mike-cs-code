from inversePysces import *

if __name__=="__main__": 
    ########################################################    
    ###################MCMC parameters######################
    ########################################################    
    N = 2000
    burn = 500
    thin = 10
    thinMore = 10
    
    logRateConstantStdDev = 1.5
    logInitialConcentrationsStdDev = 1.5
    ########################################################
    
    force = False
    nP = 4
    histOn = False
    corrOn = False
    modelsOn = False
    traceOn = True
    distOn = False
    fitOn = True
    fitCorrOn = True
    
    
    modelPath = '/media/ubuntuData/CODE/pysces/ret_erk_simple'
    if not os.path.isdir(modelPath): 
        print 'model path not found...'
        modelPath = 'D:\\CODE\\pysces\\ret_erk_simple'
        print modelPath +' is the new path'
    fileDic = {
        'pklName':'simple_bigRun_series_2_0',
        'dataFile':'ret_erk_peak_timeseries_2.csv',
        'pscFile':'ret_erk_s.psc', 
        'speciesDictFile':'ret_erk_species.dic',
        'measuredSpeciesFile':'ret_erk_MeasuredSpecies.txt',
        'referenceParameterFile':'best_referenceParameters_pilot_2.dic',
        'bestRefParameterFile':'best_referenceParameters.dic'
        }
    for k in fileDic.keys():
        fileDic[k] = os.path.join(modelPath, fileDic[k])
    
    
    
    #\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    #             initialize data
    #\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    print "data initializing..."
    D = [ret_erk_experiment(fileDic['dataFile'], fileDic['speciesDictFile'])for i in range(nP)]    
    measuredSpecies = ['R_star', 'ERK_half_star']
    print "done..."
    #\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    
    #-------------------------------------------------------
    #             initialize the pysces model
    #-------------------------------------------------------


    
    print "initializing forward models..."
    FM = [pyscesForwardModel(fileDic['pscFile'], 
                 fileDic['speciesDictFile'], 
                 fileDic['measuredSpeciesFile'],
                 fileDic['referenceParameterFile'],
                 D[i].erkTime(),
                 D[i].constraintFunctionList()
                 ) for i in range(nP)]
    print "done. prepare to be patient."
    #-------------------------------------------------------
    




    #////////////////////////////////////////////////////////
    #               get initial parameters
    #////////////////////////////////////////////////////////
        #determine number of k's    
    k0 = FM[0].k0
    s0 = FM[0].s0
    #////////////////////////////////////////////////////////
    
    
    ########################################################    
    ########### read from disk or run simulation############
    ########################################################
    if not force:
        try:
            IMPool = savedInverseModelPool(fileDic['pklName'])
            IMPool.trace()   
        except ValueError:
            print "load failed! forcing computation."
            force = True
        
    if force:            
        print "computation is forced! will now initialize inverse model list..."        
        #initialize inverse model
        IMList = []
        
        for i in range(nP):
            print "initializing %d of %d"%(i+1,nP)
            IMList += [inversePyscesModel(FM[i], 
                                    D[i], 
                                    logRateConstantStdDev = logRateConstantStdDev, 
                                    logInitialConcentrationsStdDev = logInitialConcentrationsStdDev, 
                                    verbose = False)]
        print "list initialized, setting up the pool..."
        IMPool = inverseModelPool(IMList, fileDic['pklName'])
        print "done. now for the tough part..."
        
        #sample from the posterior    
        start = time.time()    
        print "here we go!"
        IMPool(N, burn = burn, thin = thin)
        elapsed = time.time()-start
        print 'sampled %d in %.2f s: %.2e samples/s'%(N*nP,elapsed,nP*N/elapsed  )
        

        
    #break out results for easy plotting
    thinSlice = slice(0, None, thinMore)
    tK = IMPool.tK.copy()[:,thinSlice]
    tS = IMPool.tS.copy()[:,thinSlice]
    tD = IMPool.tD.copy()[...,thinSlice]
    tLogP = IMPool.tLogP.copy()[thinSlice]
    logNormConst = IMPool.invModList[0].D_obs.log_normConst()
    
    gidx = flatnonzero(tLogP>-300)
    if gidx.size==0:
        gidx = flatnonzero(tLogP>median(tLogP))
    
    print "discarding %d of %d due to extreme bad fit"%(len(tLogP)-len(gidx), len(tLogP))
    
    tK = tK[:,gidx]
    tS = tS[:,gidx]
    tD = tD[gidx]
    tLogP = tLogP[gidx]
        
    
    ########################################################
    
    
    
    
    #////////////////////////////////////////////////////////
    #               posterior histograms
    #////////////////////////////////////////////////////////
    if histOn:
        subN = (3,3)
        nFig = ceil(k0.shape[0]/prod(subN)) 
        x = linspace(-3*logRateConstantStdDev, 3*logRateConstantStdDev, 200)
        for i in range(tK.shape[0]):
            idx = mod(i, prod(subN))        
            if idx==0:
                figure()
            subplot(subN[0],subN[1],idx+1)
            
            hist(tK[i],50, label = 'posterior', normed = True)
            plot(x, stats.norm.pdf(x,loc = 0, scale = logRateConstantStdDev ), label = 'prior')
            title('k%d, ref %.1f'%(i+1, k0[i]), fontsize = 16)
            xlabel(r'log(k/k$_{ref}$)', fontsize = 14)
            
        subN = (3,3)    
        nFig = ceil(s0.shape[0]/prod(subN)) 
        x = linspace(-3*logInitialConcentrationsStdDev, 3*logInitialConcentrationsStdDev, 200)
        for i in range(tS.shape[0]):
            idx = mod(i, prod(subN))        
            if idx==0:
                figure()
            subplot(subN[0],subN[1],idx+1)
            
            hist(tS[i],50, label = 'posterior', normed = True)
            plot(x, stats.norm.pdf(x,loc = 0, scale = logInitialConcentrationsStdDev ), label = 'prior')
            title(r'$s$'+'%d, ref %.1e'%(i+1, s0[i]), fontsize = 16)
            xlabel(r'log($s/s_{ref}$)', fontsize = 14)
     #////////////////////////////////////////////////////////   
#        
#    
#    
    #-------------------------------------------------------
    #             random models
    #------------------------------------------------------
    if modelsOn:
        didx = randint(0,tD.shape[0], 4)
        print "randomly selected trace points are:"    
        print didx
        figure()
        for i, idx in enumerate(didx):
            subplot(2,2,i+1) 
            plot(D[0].erkTime(), D[0].val().T, 'k')
            plot(D[0].erkTime(), squeeze(tD[idx]).T, 'y')
            title('log(P(D|m)) = %.2f'%tLogP[idx])
        #------------------------------------------------------
    #
    #
        #-------------------------------------------------------
        #             best model
        #------------------------------------------------------
        bidx = argmax(tLogP)
        figure()
        plot(D[0].erkTime(), D[0].val().T, 'k')
        plot(D[0].erkTime(), squeeze(tD[bidx]).T, 'y')
        title('log(P(D|m)) = %.2f'%tLogP[bidx])
        IMPool.invModList[0].bounceReferenceParameters(tK[:,bidx], tS[:,bidx], 
                                fileDic['bestRefParameterFile'])
        #------------------------------------------------------    
    
    
    #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    #:::::::::::::   intersample distances  :::::::::::::::::::::::::
    #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    if distOn:
        dk = diff(tK[:], axis = 1)
        dkMag = sum(dk**2, axis = 0)**.5/dk.shape[0]**.5
        figure()
        subplot(1,2,1)
        plot(dkMag)
        xlabel('transition number')
        ylabel('normalized distance')
        title('change in position per transition for k vector\n normalized by unit change in every component, N**.5')
        subplot(1,2,2)
        hist(log10(dkMag[dkMag>0]), 75)
        xlabel('log10(distance)')
        title('per transtion log distance (unnormalized)\n%d of %d zero values discarded'%(sum(dkMag==0), len(dkMag)))
    #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::    
    
    
    ########################################################    
    ################### correlation coeffs##################
    ########################################################
    if corrOn:
        figure()
        pcolor(arange(1.,tK.shape[0]+1)-.5, 
               arange(1.,tK.shape[0]+1)-.5,
               corrcoef(tK[:])-eye(tK.shape[0]))
        colorbar()
        xlabel('constant index')
        ylabel('constant index')
        title('correlation coefficient between rate constants')
    ########################################################
    
    #.......................................................
    #................goodness of fit histograms.............
    #.......................................................
    if fitOn:
        figure()
        hist(tLogP, 100)
        xlabel('log(P(D|m))')
        title('goodness of fit, max possible log(P(D|m)) = %.2e'%logNormConst)
    print "best fitting model has goodness of fit %.2e of a possible %.2e"%(amax(tLogP),logNormConst) 
    #.......................................................    
    
    if traceOn:
        figure()
        plot(tK.T)
        title('traces of rate constants')
    
    if fitCorrOn:
        figure()
        subplot(1,3,1)
        plot(mean(abs(tK), axis = 0), tLogP, '.', alpha = .2)
        xlabel('mean abs rate constant')
        ylabel('goodness of fit')
        subplot(1,3,2)
        plot(median(abs(tK), axis = 0), tLogP, '.', alpha = .2)
        xlabel('median abs rate constant')
        ylabel('goodness of fit')
        subplot(1,3,3)
        plot(amax(abs(tK), axis = 0), tLogP, '.', alpha = .2)
        xlabel('max abs rate constant')
        ylabel('goodness of fit')
        
#    figure()
#    U, s, Vh = svd(tK)
#    stem(arange(len(s)), s)
#    title('singular values of rate constant trace')
#    