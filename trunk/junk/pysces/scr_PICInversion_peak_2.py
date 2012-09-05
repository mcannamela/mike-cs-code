from pyscesInversionTools import *


if __name__=='__main__':
        
    modelPath = '/media/ubuntuData/CODE/pysces/ret_erk_peak_2'
    if not os.path.isdir(modelPath): 
        print 'model path not found...'
        modelPath = 'D:\\CODE\\pysces\\ret_erk_simple'
        print modelPath +' is the new path'
        
    fileDic = {
    'pklName':'series_2',
    'dataFile':'ret_erk_peak_timeseries_2.csv',
    'pscFile':'ret_erk_s.psc', 
    'speciesDictFile':'ret_erk_species.dic',
    'measuredSpeciesFile':'ret_erk_MeasuredSpecies.txt',
    'referenceParameterFile':'best_referenceParameters_pilot_2.dic',
    'bestRefParameterFile':'best_referenceParameters_2.dic'
    }
    
    compute = False    
    bounce = False
    bounceDir = os.path.join(os.path.split(modelPath)[0],'ret_erk_peak_2_bounce' )
    N = 5000
    burn = 500
    thin = 10
    loadPickle = None
    
    logRateConstantStdDev = .75
    logInitialConcentrationsStdDev = .75
    
    
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
    
#    plotFuns =(' plotModel plotMarginals ' 
#                'plotGoodness plotTraces plotFitCorr'.split())
    plotFuns =['plotModel']                
    PIC.plotListOfFuns(plotFuns)
    