from particle import *
import cPickle as pkl
import particleSolver as ps
import particleModelEnsemble  as pmEns

thedir = os.path.join(os.curdir, 'lavaJetBatch_pwvel=4_pwtemp=10_')


PE = pmEns.basicParticleEnsemble(pmEns.basicGriddedParticleParameterGenerator)
PE.generate(dMinMax = [3e-6,150e-6],
                                            vMinMax = [1, 35],
                                            x0 =array([0.00635,.00792,0.]),
                                            n = (75,74))
cnt = 0
for fname in os.walk(thedir).next()[2]:
	if fname.find('gas')==-1:
		print 'solving jet %d'%cnt
		cnt+=1
		n = fname.split('_')
		jet = lavaField(simName = os.path.join(thedir, fname), fastOn = False, gasPickle = None)
		PE.ensembleSolve(plasmaField = jet)
		with open(os.path.join(os.curdir, 'jet-particleBatch_'+n[2]+'_'+n[3]+'_.pkl'), 'wb') as f:
			pkl.dump(PE, f, -1)
