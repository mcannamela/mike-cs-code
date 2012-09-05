from modelValidation import *
import gasMixture as gm
from pylab import *
from scipy.integrate import trapz


g = gm.defaultMaterial




timestep =1e-5
nSteps = 2000
nNodes = 50
diameter = [.1e-6, 50e-6] 

#x0 = array([0.00, .003537, 0.0])#T=3400
#x0 = array([0.00, .003437, 0.0])#T=3900
x0 = array([0.00, .003377, 0.0])#T=4200
#x0 = array([0.00, .0020, 0.0])#T=8500
#x0 = array([0.00, .0000, 0.0])#T=8500

try:
    PF.density(array([0,0,0]))
except NameError:
    PF = part.pf.plasmaField(fastOn = False)
    
print 'Tinf is %d'%PF.temperature(x0)

Tinf = PF.temperature(x0)
T0 = 1000

H0 = g.Enthalpy(T0)
Hinf = g.Enthalpy(Tinf)

    
p = rigidParticle(diameter = diameter, 
                            T =T0,
                            thermalTimestep = timestep,
                            kinematicTimestep = timestep,
                            nNodes = nNodes,
                            nStep = nSteps,
                            material = g)

m0 = double(p.totalMass)
print 'initial mass is %.3e'%p.totalMass
                            
p.adaptiveTimestepping=False
p.splitScales = 0

p.inject(PF, x0,array([0, 0, 0]))
        
PS = explicitRigidParticleSolver(particle = p, verbose = False)
PS.alpha = .7
PS.stepErrorThres = 1e-4
PS.maxIter = 50
PS.solve()

mf = double(p.totalMass)
print 'final mass is %.3e'%p.totalMass

qConv, qRad, qVap = [array(x) for x in (PS.convHistory, PS.radHistory, PS.vapHistory)]
Qconv, Qrad, Qvap = [sum(x)*timestep for x in [qConv, qRad, qVap]]

S = PS.getSolnH()
plot(p.material.Temperature(S[:,0]),'r.-')
plot(p.material.Temperature(S[:,-1]),'k.-')

figure()
plot(qConv, 'r.-', label = 'conv')
plot(qRad, 'b.-', label = 'rad')
plot(qVap, 'g.-', label = 'vap')
plot(qConv-qVap-qRad, 'k.-', label = 'net')

Tbar = sum((p.mesh.volume()*p.T)/sum(p.mesh.volume()))
Hbar = p.material.Temperature(Tbar)

print 'total change in particle enthalpy should be %.2e'%((Hinf-H0)*m0)
print 'but only T= %d was acheived, and thus the change is %.2e'%(Tbar, (p.material.Enthalpy(Tbar)-H0)*m0)
print 'sum of heat in and out is %.2e'%(Qconv-Qrad-Qvap)


print 'heat to vaporize the missing mass is %.2e'%(p.heatOfVaporization*(m0-mf))
print 'and the total heat lost to vaporization is %.2e'%Qvap
