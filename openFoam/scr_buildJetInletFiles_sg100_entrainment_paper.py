from paraFoamDataReader import *
from gasses import Ar
from scipy.integrate import trapz
r    = float64([  .001,    .2,    .4,  .6,   .8, 1.0])*.01
Tin  = float64([12700, 11800, 5000, 800, 300, 300])
Hin = Ar.TFuns['enthalpy'](Tin)
Uin0 = float64([1550,  1400,  450,  100, 20,  0])

R = r_[-flipud(r[1:]), r]
theta = linspace(0,2*pi, 101)[:-1]

p = r_['0,2', (0*r*cos(theta)[...,newaxis]).flatten(), 
               (r*cos(theta)[...,newaxis]).flatten(), 
                (r*sin(theta)[...,newaxis]).flatten()]

blank = ones_like(0*r*cos(theta)[...,newaxis])

t = float64([0,1e-3])
T = (blank*Tin).flatten()
H = (blank*Hin).flatten()
U = array([(blank*Uin0).flatten(), (0*blank).flatten(), (0*blank).flatten()])

#DTO = dummyTorchOutlet(p, t, r_['1', 300+.001*T[..., newaxis],T[..., newaxis]], 
#                             r_['1', 305000+.001*H[..., newaxis],H[..., newaxis]],
#                             r_['2', .001*U[..., newaxis],U[..., newaxis]])
#
#DTO.setup(300, 156096, .0805, N = 300)
#
#outFile = '/media/raidArray/SG100_entrainment_paper_jet/inlet'
#    
#TLW = torchLUTWriter()
#TLW(DTO, outFile)

rho0 = 1.2
T0 = 300.
sly = slice(0,3)
Afac =trapz(r[sly],r[sly])
Tbar = trapz(r[sly]*Tin[sly], r[sly])/Afac
rTbar = trapz(r[sly]/Tin[sly], r[sly])/Afac
rhoUBar = trapz(rho0*T0*r[sly]*Uin0[sly]/Tin[sly],r[sly])/Afac
uBar = trapz(r[sly]*Uin0[sly],r[sly])/Afac
rhoBar = rho0*T0*rTbar
uBar_rho = trapz(rho0*T0*r[sly]*Uin0[sly]/Tin[sly],r[sly])/rhoBar


print 'Tbar is %d'%Tbar
print '1/Tbar is %.2e'%(1.0/Tbar)
print '(1/T)bar is %.2e'%(rTbar)


print 'incoming density and temperature are %f and %d'%(rho0, T0)
print 'mean outgoing density is %.2e'%(rhoBar)

print 'mean rho*U: %.2e'%rhoUBar
print 'mean U: %.2e'%uBar
print 'mean rho * mean U: %.2e'%(uBar*rhoBar)
