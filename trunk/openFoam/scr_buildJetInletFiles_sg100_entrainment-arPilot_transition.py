from paraFoamDataReader import *
from gasses import Ar
r    = float64([  .001,    .2,    .4,  .6,   .8, 1.0])*.01
Tin  = float64([12700, 11800, 5000, 800, 300, 300])
Hin = Ar.TFuns['enthalpy'](Tin)
Uin0 = float64([1550,  1400,  450,  100, 20,  0])

T_target = 11858.0
U_target = 513.0


R = r_[-flipud(r[1:]), r]
theta = linspace(0,2*pi, 101)[:-1]

p = r_['0,2', (0*r*cos(theta)[...,newaxis]).flatten(), 
               (r*cos(theta)[...,newaxis]).flatten(), 
                (r*sin(theta)[...,newaxis]).flatten()]

blank = ones_like(0*r*cos(theta)[...,newaxis])

t = float64([0,.5e-3])
T = (blank*Tin).flatten()
H = (blank*Hin).flatten()
U = array([(blank*Uin0).flatten(), (0*blank).flatten(), (0*blank).flatten()])

DTO = dummyTorchOutlet(p, t, r_['1', T[..., newaxis],(T_target/Tin[0])*T[..., newaxis]], 
                             r_['1', 305000+.001*H[..., newaxis],H[..., newaxis]],
                             r_['2', U[..., newaxis],(U_target/Uin0[0])*U[..., newaxis]])

DTO.setup(300, 156096, .0805, N = 300)

outFile = '/home/wichtelwesen/liveOpenFoamRuns/jet_sg100_ArPilot_40_500_realizableKE_transition/inlet'
    
TLW = torchLUTWriter()
TLW(DTO, outFile)