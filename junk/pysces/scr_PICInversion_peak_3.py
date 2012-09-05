from pyscesInversionTools import *


if __name__=='__main__':
        
    modelPath = '/media/ubuntuData/CODE/pysces/ret_erk_peak_3_bounce'
    if not os.path.isdir(modelPath): 
        print 'model path not found...'
        modelPath = 'D:\\CODE\\pysces\\ret_erk_simple'
        print modelPath +' is the new path'
        
    fileDic = {
    'pklName':'series_3',
    'dataFile':'ret_erk_peak_timeseries_3.csv',
    'pscFile':'ret_erk_s.psc', 
    'speciesDictFile':'ret_erk_species.dic',
    'measuredSpeciesFile':'ret_erk_MeasuredSpecies.txt',
    'referenceParameterFile':'best_referenceParameters_3_pilot.dic',
    'bestRefParameterFile':'best_referenceParameters_3.dic'
    }
    
    compute = False    
    bounce = False
    bounceDir = os.path.join(os.path.split(modelPath)[0],'ret_erk_peak_3_fullBounce' )
    N = 2500
    burn = 500
    thin = 10
    loadPickle = None
    
    logRateConstantStdDev = 1.0
    logInitialConcentrationsStdDev = 1.0
    
    
    PIC = pyscesInversionCase(modelPath, fileDic = fileDic, nP = 4)
    
    if compute:
        PIC(N,burn,thin,
                logRateConstantStdDev = logRateConstantStdDev, 
                logInitialConcentrationsStdDev = logInitialConcentrationsStdDev,
                verbose = False)
        PIC.save()
                
    else:
        PIC.load(loadPickle)
    
    if bounce:
        PIC.updateReferenceParameters()        
        PIC.bounceModel(bounceDir)
    
#    plotFuns =(' plotModel plotMarginals plotIntersampleDistance ' 
#                'plotGoodness plotTraces plotFitCorr'.split())
    plotFuns =['plotModel']
    PIC.plotListOfFuns(plotFuns)
    