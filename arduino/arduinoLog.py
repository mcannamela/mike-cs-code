import serial
import time
import sys
import pdb
import scipy
from pylab import *
from numpy import *

class listeningLogger(object):
    def __init__(self, commPortNumber = 10, 
                         logfileName = 'maple.log',
                         silent = True, 
                         plotInterval = 1.0):
        self.commPortNumber = commPortNumber
        self.logfileName = logfileName
        self.silent = silent
        
        self.S = serial.Serial(commPortNumber, timeout = 3)
        self.f = open(logfileName, 'w')
        self.X = []
        self.t = []
        self.dt = 0
        #self.fig = figure()
        self.plotInterval = plotInterval
        self.lastPlot = 0
        
    def __call__(self, logDuration = inf):
        
        
        try:
            self.t0 = time.time()
            while (time.time()-self.t0)<logDuration:
                t, L = self.rawQuery()
                
                dt, X = self.cook(L)
                
                if X.shape[0] == 0:
                    continue
                
                self.dt = dt
                self.X += [X]
                self.t += [t]
                
                self.log(self.cookedLine(t, X))
                if (time.time()-self.lastPlot)>self.plotInterval:
                    self.plotCurrent()

            self.close()      
        except KeyboardInterrupt:
            self.close()

            
        return (array(self.t), self.X, self.dt)
    def plotCurrent(self):
        pass
    def rawQuery(self):
        if self.S.inWaiting()>0:
            L = self.S.readline()
            t = time.time()
            return (t, L)
        else:
            return (time.time(), '')
        
    def cook(self, L):
        sL = L.replace('\n','').split()
        try:
            dt = uint32(sL[0])
        except IndexError:
            dt = 0
            X = array([])
            return (dt, X)
        try:
            X = int16(sL[1:])
        except ValueError:
            pdb.set_trace()
        return (dt, X)
        
    def cookedLine(self, t, X):
        L = ''
        for x in X:
            L+= '    %d'%x
            
        return "%.2e"%(t-self.t0) + L+'\n'
        
    def log(self, cookedLine):
        self.f.write(cookedLine.replace('\n', ''))
        self.f.write('\n')
        if not self.silent:
            print cookedLine
    
    def close(self):
        self.S.close()
        self.f.close()  
    
class queryLogger(object):
    def __init__(self, commPortNumber = 6, 
                 logfileName = 'arduino.log.', 
                 sampleInterval = 1, 
                 silent = True):
        
        self.commPortNumber = commPortNumber
        self.logfileName = logfileName
        self.sampleInterval = sampleInterval
        self.silent = silent
        
        self.S = serial.Serial(commPortNumber, timeout = 3)
        self.f = open(logfileName, 'w')
        self.X = []
        self.t = []
     
    def __call__(self, logDuration = Inf):
        try:
            t0 = time.time()
            while (time.time()-t0)<logDuration:
                t, L = self.rawQuery()
                
                X = self.cook(L)
                
                if X.shape[0] == 0:
                    continue
                
                self.X += [X]
                self.t += [t]
                
                self.log(self.cookedLine(t, X))
                time.sleep(self.sampleInterval)
                
            self.close()
            
                    
                    
        except KeyboardInterrupt:
            self.close()

            
        return (array(self.t), array(self.X))
            
        
    def rawQuery(self):
        self.S.write('t')
        L = self.S.readline()       
        t = time.time()
        return (t, L)
        
    def cook(self, L):
        X = int32(L[:-1].split())
        return X
        
    def cookedLine(self, t, X):
        L = ''
        for x in X:
            L+= '    '+ str(x)
            
        return time.asctime(time.localtime(t)).replace(' ', '_') + L+'\n'
        
    def log(self, cookedLine):
        self.f.write(cookedLine.replace('\n', ''))
        self.f.write('\n')
        if not self.silent:
            print cookedLine
    
    def close(self):
        self.S.close()
        self.f.close()

class floatLogger(queryLogger):
    def cook(self, L):
        try:
            X = float64(L[:-1].split())
        except:
            pdb.set_trace()
        return X
      
class htm2500Logger(queryLogger):
    def __init__(self, commPortNumber = 6, 
                 logfileName = 'arduino.log.', 
                 sampleInterval = 1,
                 silent = True):
        queryLogger.__init__(self, 
                             commPortNumber, 
                             logfileName, 
                             sampleInterval,
                             silent)
        
        self.logHeader = 'timestamp    humidity_1    temp_1    humidity_2   temp_2    humidity_3   temp_3'
        self.f.write(self.logHeader+'\n')

        #from the htm2500 datasheet
        R = flipud(float64([169149, 125546,94143, 71172, 54308, 41505, 32014, 25011, 19691, 15618,
                     12474, 10000, 8080, 6569, 5372, 4424, 3661, 3039, 2536, 2128]))
        T = flipud(linspace(-30, 65, 20))
        
        VH = float64([925,1080 , 1235, 1390, 1540, 1685, 1825, 1960, 2090, 2220, 2350,  
                      2480,  2605, 2730, 2860, 2990, 3125, 3260, 3405, 3555, 3705 ])
        RH = linspace(0, 100, 21)
        
        self.thermistorBallastResistance = 10000. #ohms of the resistor in the voltage divider
        self.thermistorSupplyVoltage = 5000
         
        self.humidityLookup = lambda x: interp(x, VH, RH)
        self.temperatureLookup = lambda x: interp(x, R, T) 
        
    def cook(self, L):
        x = queryLogger.cook(self, L)
        if x.shape[0] == 0:
            return x
        tSly = slice(1, None, 2)
        hSly = slice(0, None, 2)
        
        VT = float64(x[tSly])#millivolts across thermistor
        VH = float64(x[hSly])# millivolts out on humidity pin
        
        
        T = self.cookTemperature(VT)
        H = self.cookHumidity(VH, T)
        
        X = zeros(x.shape[0])
        X[tSly] = T
        X[hSly] = H
        
        return X
        
     
    def thermistorResistance(self, VT):
        return self.thermistorBallastResistance*VT/(self.thermistorSupplyVoltage-VT)
        
    def cookTemperature(self, VT):
        R = self.thermistorResistance(VT)
        return array([round(T, 1) for T in self.temperatureLookup(R)])
        
    def cookHumidity(self, VH, T): 
#        a = float64([-1.9206e-3, 1.437e-5, 3.421e-3, -12.4])
#        Tfactor = (1+(T-23)*2.4e-3)**-1
#        x = flipud(float64([VH**i for i in range(4)]))
        
#        return sum(x*a)*Tfactor

        return array([round(H, 1) for H in self.humidityLookup(VH)])
        



if __name__=='__main__':
#    logger = queryLogger()
    logger = floatLogger(0, silent = False, logfileName = 'brew_temps_feb_20.log')

        
    t, x = logger()
    logger.close()
#    HL = htm2500Logger(9,silent = False)
#    t,x = HL()