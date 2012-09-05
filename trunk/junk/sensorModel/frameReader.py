from numpy import *
from pylab import *
from scipy import signal
import pdb

class abstractFrameReader(object):
    def __init__(self):
        pass
        
    def byteRead(self, fname):
        with open(fname, 'rb') as f:
            self.byteArray = float64(array(map(ord,f.read())))
    def parse(self):
        pass
        
    def read(self, fname):
        pass
    

class basicFrameReader(abstractFrameReader):
    def __init__(self):
        self.bwLenIdx = range(16,20)
        self.colorLenIdx =  range(8183,8187)
        self.frameLen = 20432
    def getLen(self, idx):
        return uint32(sum(flipud(2**(arange(array(idx).shape[0])*8))*self.byteArray[idx]))
    def read(self, fname):
        self.byteRead(fname)
        self.nFrame = self.byteArray.shape[0]/self.frameLen
        self.bwLen = self.getLen(self.bwLenIdx)
        self.colorLen = self.getLen(self.colorLenIdx)
        bwStart = self.bwLenIdx[-1]+1
        colorStart = self.colorLenIdx[-1]+1
        self.bw = uint32(zeros((self.nFrame, self.bwLen)))
        self.color = uint32(zeros((self.nFrame, 3, self.colorLen)))
        
        
        for i in range(self.nFrame):
            self.bw[i] = self.byteArray[bwStart+i*self.frameLen: (bwStart+self.bwLen+i*self.frameLen)]
            self.color[i] = self.byteArray[colorStart+i*self.frameLen:(colorStart+3*self.colorLen+i*self.frameLen)].reshape(3,self.colorLen)


if __name__=="__main__":
    fname = "E:\\CODE\\PointSpread\\PointSpread\\PatternSimpleLens_0-00"
    FR = basicFrameReader()
    FR.read(fname)
    plot(sum(FR.bw, axis = 0))
    figure()
    plot(sum(FR.color, axis=0).T)