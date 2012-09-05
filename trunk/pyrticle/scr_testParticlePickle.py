from particle import *
import cPickle as pkl
import particleSolver as ps
import particleModelEnsemble  as pmEns
def objectBlacklist(obj):
    """
    find out which of an ojects attributes cannot be pickled
    """
    d = obj.__dict__
    black = []
    for k in d.keys():
        try: 
            pkl.dumps(d[k])
        except:
            black+=[k]
    return black
    
PF = plasmaField(fastOn = True)
p = particle()
ep = emittingParticle()
p.inject(PF)
PS = ps.particleSolver(p)
PS.solve()
PE = pmEns.basicParticleEnsemble(pmEns.basicGriddedParticleParameterGenerator, 'pickleTestEns')
PE.generate(n = (2,2))

PE.ensembleSolve()
print 'done!'

L = [ p,ep, PE,  PS, gm.gasMixture()]
n = [ 'particle','emitting particle', 'particle ensemble', 'particle solver', 'gas mixture']
print "begin pickling..."
for i in range(len(L)):
    try:
        pkl.dumps(L[i],-1)
        print n[i]+'...success!'
    except:
        print n[i]+'...failed.'
        
print "begin unpickling..."
for i in range(len(L)):
    try:
        pkl.loads(pkl.dumps(L[i],-1))
        print n[i]+'...success!'
    except:
        print n[i]+'...failed.'

print '\n testing reconstituted particles...'
pe = pkl.loads(pkl.dumps(PE))
pp = pe.modelArray.flatten()[0].PS.p
print 'specific heat after reconstitution:'
print pp.specificHeat()
pp.T -= 500.
print 'values should now be changed:'
print pp.specificHeat()

print 'plasma temp:'
print pp.plasmaTemperature
pp.X[pp.stepCount-1] = array([.001,.001,.001])
print 'values should now be changed:'
pp.refreshCache()
print pp.plasmaTemperature

print 'plasma Temperature'

#pfBlack = objectBlacklist(PF)
#psBlack = objectBlacklist(PS)
#gasBlack = objectBlacklist(gm.gasMixture())
#peBlack = objectBlacklist(PE)
#
#print 'blacklists for pf, ps, gas, pe:'
#print pfBlack
#print psBlack
#print gasBlack
#print peBlack

    
