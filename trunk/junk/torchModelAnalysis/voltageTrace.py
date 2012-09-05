from numpy import *
from pylab import *
from scipy import fft
from scipy.signal import hamming, boxcar, hanning, convolve, gaussian
from scipy.interpolate import interp1d
from scipy.integrate import trapz
import pdb

def readTrace(fname):
    with open(fname, 'r') as f:
        F = f.readlines()[1:]
    try:
        X = float64([L.split(',') for L in F]).T
    except:
        X = float64([L.split(';') for L in F]).T
    
    t = X[2]
    V = X[0]
    
    return t.copy(), V.copy()
    
def zCross(x):
    return flatnonzero(abs(diff(sign(x)))>0)
    
class voltageTrace(object):
    def __init__(self, t, V):
        self.tRaw = t.copy()
        self.VRaw = V.copy()
        self.DCV = trapz(self.VRaw, self.tRaw)/(self.tRaw[-1]-self.tRaw[0])
        self.process(len(t))
    
    def process(self, N):
        
        zc = zCross(self.VRaw-self.DCV)
        sly = slice(zc[0], zc[-1])
        self.DCV = trapz(self.VRaw[sly], self.tRaw[sly])/(self.tRaw[zc[-1]]-self.tRaw[zc[0]])
        
        t0 = self.tRaw[zc[0]]
        t1 = self.tRaw[zc[-1]]
        self.t = linspace(0, t1-t0, zc[-1]-zc[0])
        
        self.V = interp1d(self.tRaw, self.VRaw-self.DCV)(t0+self.t)
        
        self.fNy = .5/diff(self.t[:2])
        
    def plot(self, N = None):
        if N!=None:
            self.process(N)
        
        figure()
        plot(self.t*1e3, self.V+self.DCV)
        xlabel('time, ms')
        ylabel('voltage, V')
        
    def fft(self, N = None, winFun = None):
        if N!=None:
            self.process(N)
        
        if winFun == None:
            winFun = boxcar
        
        W = winFun(len(self.t))
        W/=sum(W)
        
        M = uint32(2**ceil(log2(double(len(self.t)))))
        
        self.Vhat = fftshift(fft(self.V*W, M))
        self.f = linspace(-1,1, M)*self.fNy
        
    def plotPowerSpectrum(self, nSmooth = 0, N = None, winFun = None):
        
        
        self.fft(N, winFun)
        
        y = self.Vhat.copy()
        if nSmooth >1:
            y = convolve(y, boxcar(nSmooth)/nSmooth, mode = 'same')
        
        figure()
        db = 10*log10(abs(y)**2)
        
        plot(self.f, db)
        xlim(0,self.f[-1] )
        xlabel('frequency, Hz')
        ylabel('power, dB')
        

if __name__=="__main__":
#    fname = '/media/ubuntuData/Dropbox/voltageTrace_blowout.csv'
#    fname = '/media/ubuntuData/Dropbox/sg100_40_700_voltageTrace0.0.csv'
#    fname = '/media/ubuntuData/Dropbox/sg100_voltageTrace_40_700_conductiveLayer0.0.csv'
#    fname = '/home/wichtelwesen/Dropbox/sg100_ArH20_40_700_voltageTrace0.0.csv'
#    fname = '/home/wichtelwesen/Dropbox/sg100_ar_40_700_voltageTrace0.0.csv'
    fname = r'D:\My Dropbox\sg100_ar_40_700_voltageTrace0.0.csv'
    fname = r'D:\My Dropbox\sg100_ArH20_40_700_voltageTrace0.0.csv'
    t,V = readTrace(fname)
    
    #V = mean(V)+convolve(V-mean(V),gaussian(12,1.5)/sum(gaussian(12,1.5)), mode = 'same' )
    vt = voltageTrace(t,V)
    
    vt.plot()
    vt.plotPowerSpectrum()
#    plot(t,V)
