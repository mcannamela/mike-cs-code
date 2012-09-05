from arduinoLog import listeningLogger
from numpy import *
from pylab import *
import serial

L = listeningLogger(commPortNumber = '/dev/ttyACM0', silent = True)

t, X, dt = L(30)


figure()
[plot(3.3*convolve(x, boxcar(4)/4.0)/4096.0,'.-') for x in X]
print "sampling period in us is %d"%dt
print "that's %.2e Hz"%(1e6/dt)

I = '10V'
Q = 'test'

with open('sg100_voltage_200kHz_10'+Q+'_'+I+'.pkl', 'wb') as f:
    pkl.dump(X[1:], f, -1)

#show()

with open('sg100_voltage_200kHz_10'+Q+'_'+I+'.pkl', 'rb') as f:
    Y = pkl.load(f)
    
figure()
[plot(3.3*convolve(x, boxcar(4)/4.0)/4096.0,'.-') for x in Y]

show()