from plasmaField import *

#runFolder = "/media/bigBoy/openFoamRuns/jet_sg100_ar_40_500_sigmaTe_final"
runFolder = "/media/fatMan/taggedOpenFoamRuns/jet_sg100_ArPilot_40_500_highDensity"
    
with open('/media/raidArray/CODE/gasModel/argonPropertiesDict.pkl', 'rb') as f:
    DAr = pkl.load(f)
    DAr['enthalpy']-=DAr['enthalpy'][0]
Ar = gas(40.0, DAr, 10000)


with open('/media/raidArray/CODE/gasModel/airPropertiesDict.pkl', 'rb') as f:
    DA = pkl.load(f)
    DA['enthalpy']-=DA['enthalpy'][0]
Air = gas(29.0, DA, 10000)

G = binaryGasMixture(Ar, Air)


tMin = .00300
tMax = .00600
xMin = array([0.0, -.02, -.02])
xMax = array([.1, .02, .02])

P = openFoamPlasma(runFolder, tMin, tMax, xMin,xMax, G)
#with open('jet_ar_40_500_sigmaTe.pkl','wb') as f:
with open('jet_arPilot_40_500_highDensity.pkl','wb') as f:
    pkl.dump(P, f, -1)