from numpy  import *
from pylab import *
import mlabMacros as mlm
import cPickle as pkl
import particleModelEnsemble as pmEns
from plasmaField import lavaField
import pdb
import os
import scipy

def showField(plasmaField, n = 10, temperature = True, velocity = True):
    """
    make contour plots of the field
    """
    z = plasmaField.D['axialCoord']
    r = hstack((-flipud(plasmaField.D['radialCoord']),plasmaField.D['radialCoord']))
    T = vstack((flipud(transpose(plasmaField.D['temperature'])), transpose(plasmaField.D['temperature'])))
    try:
        V = vstack((flipud(transpose(plasmaField.D['velocity'])), transpose(plasmaField.D['velocity'])))
    except:
        V = vstack((flipud(transpose(plasmaField.D['axialVelocity'])), transpose(plasmaField.D['axialVelocity'])))
    figure(figsize = (20, 8))
    if temperature and velocity:
        figure()        
        subplots_adjust(hspace = .5)
        subplot(2,1,1)
        cs = contourf(z*1e3,r*1e3, T,n, cmap = cm.Greys); colorbar(format = '%d')
        xlabel('x, mm')
        ylabel('y, mm')
        title('temperature of the jet, in Kelvin')
        subplot(2,1,2)
        cs = contourf(z*1e3,r*1e3, V,n, cmap = cm.Greys); colorbar(format = '%d')
        xlabel('x, mm')
        ylabel('y, mm')
        title('speed of the jet (to the right), in m/s')
        
        
    elif temperature:
        cs = contour(z,r, T,n, colors = 'k'); clabel(cs, inline = True, fmt = '%d')
        
    elif velocity:
        cs = contour(z,r, V,n, colors = 'k'); clabel(cs, inline = True, fmt = '%d')
        

name = 'sharpTemp_600_38.pkl'
LF = lavaField(simName = os.path.join(os.pardir,'pyrticle',name),  fastOn = False, gasPickle = None)

showField(LF)

with open(os.path.join(os.pardir,'pyrticle','twoPeaksMappingSharpT.pkl'), 'rb') as f:
    DD = pkl.load(f)

D = DD['D']
V = DD['V']
d = DD['d']
v = DD['v']

z = DD['z']
r = DD['r']
S = DD['S']



#diameter varies with row, injection velocity with column
Tmap = scipy.ndimage.gaussian_filter(S['surfaceTemperature'],1.,mode = 'nearest')
Hmap = scipy.ndimage.gaussian_filter(S['enthalpy']/(4*((D/2)**3)/(3*pi)),1.,mode = 'nearest')


pts = ['A','B', 
       'C', 'D',
       'E','F',
       'G','H']

dAC = [20,100]
vEG = [5, 20]

figure()
plot(v, Tmap[argmin(abs(d*1e6-dAC[0]))], '-k', label = 'small motes, d = %d'%dAC[0]+r' $\mu$m' )
plot(v, Tmap[argmin(abs(d*1e6-dAC[1]))], '-.k', label = 'big motes, d = %d'%dAC[1]+r' $\mu$m' )
xlabel('how fast we shot the dust mote into the jet, m/s')
ylabel('temperature of the dust when it hit the target, in Kelvin')
legend(loc = 'upper left')

figure()
plot(d, Tmap[:,argmin(abs(v-vEG[0]))], '-k', label = ' slow motes, v = %d m/s'%vEG[0])
plot(d, Tmap[:,argmin(abs(v-vEG[1]))], '-.k', label = 'fast motes, v = %d m/s'%vEG[1])
xlabel(r'size of dust mote, $\mu$m')
ylabel('temperature of the dust when it hit the target, in Kelvin')
legend()

figure()
cs = contour(d*1e6,v,Tmap.T, 30, colors = 'k')
xlabel(r'size of dust mote, $\mu$m')
ylabel('how fast we shot the dust mote into the jet, m/s')
title('temperature of the dust when it hit the target, in Kelvin' )
#clabel(cs,  fmt = '%d', manual = True)

#show()
