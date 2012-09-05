from inversePysces import *

def getSpecies(FM):
    return FM.pscMod.data_sim.getSpecies().T[1:].copy()
def initSpecies(FM):
    return float64([FM.pscMod.__dict__[k+'_init'] for k in FM.pscMod.species])
    
def cook(S):
    """
    after solving, pick out the measured species and do any processing to
    make the result comparable to any measured data
    """
    s = atleast_2d(S)
    for i in range(s.shape[0]):
        s[i] /= amax(s[i])
        s[i] /= s[i][0]
    return squeeze(s)

if __name__== "__main__":   
    
    modelPath = '/media/ubuntuData/CODE/pysces/ret_erk_scratch_model'
    if not os.path.isdir(modelPath):
        modelPath = 'D:\\CODE\\pysces\\ret_erk_scratch_model'
    
    fileDic = {
        'pklName':'scratch_pERK.0.pkl',
        'dataFile':'ret_erk_timeseries.csv',
        'pscFile':'R-ERK_s.psc', 
        'speciesDictFile':'R-ERK_s_species.dic',
        'measuredSpeciesFile':'R-ERK_s_MeasuredSpecies.txt',
        'referenceParameterFile':'R-ERK_s_referenceParameters.dic'
        }
    for k in fileDic.keys():
        fileDic[k] = os.path.join(modelPath, fileDic[k])
    
        
    #\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    #             initialize data
    #\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    D = ret_erk_experiment(fileDic['dataFile'], fileDic['speciesDictFile'])
    
    measuredSpecies = ['RetStar', 'ppERKHalf']
    #\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    
    #-------------------------------------------------------
    #             initialize the pysces model
    #-------------------------------------------------------
    FM = pyscesForwardModel(fileDic['pscFile'], 
                 fileDic['speciesDictFile'], 
                 fileDic['measuredSpeciesFile'],
                 fileDic['referenceParameterFile'],
                 D.erkTime(),
                 D.constraintFunctionList()
                 )
                 
    pm = constrainedPysMod(os.path.split(fileDic['pscFile'])[1], os.path.split(fileDic['pscFile'])[0])
    #pm = pysces.model(os.path.split(fileDic['pscFile'])[1], os.path.split(fileDic['pscFile'])[0])
    
    for k in pm.parameters:
        pm.__dict__[k]*=1
        
    pm.sim_end = 6000
    pm.sim_points = 200
    pm.Simulate(D.constraintFunctionList()[0])
    #pm.Simulate()
    
    S = pm.data_sim.getSpecies().T[:,3:]
    figure()
    subplot(1,2,1)
    plot(S[0], S[[1,-1]].T)
    title("raw")
    subplot(1,2,2)
    plot(S[0], cook(S[[1,-1]]).T)
    title("cooked")
#
#    show()
    #-------------------------------------------------------
    
    #////////////////////////////////////////////////////////
    #               get initial parameters
    #////////////////////////////////////////////////////////
        #determine number of k's    
    k0 = FM.k0
    s0 = FM.s0
    #////////////////////////////////////////////////////////
    
    D_mod = FM(0*k0,0*s0)
    Sfm = FM.S.copy()
    t = FM.t()

    for i in range(1):
        figure()
        subplot(1,2,1)
        plot(t, Sfm[i,-1])
        subplot(1,2,2)
        plot(t, cook(Sfm[i,-1]))
    
    
#    for nD in range(4):
#        figure()
#        plot(t, S[nD,0], 'y-')
#        plot(D.t, D.RET[nD][1], 'k.-')
#        plot(FM.F[nD].T, FM.F[nD].S,'b')
#        plot(FM.F[nD].t, FM.F[nD].s,'g')
