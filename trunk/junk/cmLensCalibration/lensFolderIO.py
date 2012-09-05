from numpy import *
from pylab import *
import struct
import os
import pdb

class colorFrameReader(object):
    """
    reads just the color from a flux sentinel frame
    
    _TODO__TODO__TODO__TODO__TODO__TODO__TODO__TODO__TODO_
    _TODO__TODO__TODO__TODO__TODO__TODO__TODO__TODO__TODO_
   
         will need to change this class when the new sentinel
         is up and running. this should really be sub-classed from 
         an abstract base class, but is probably suitable for 
         subclassing when the new sentinel is in full use
    _TODO__TODO__TODO__TODO__TODO__TODO__TODO__TODO__TODO_
    
    """
    def __init__(self, offset = 48664):
        """
        to initialze, just tell us how many bytes to skip at the beginning
        offset - should be the byte where [offset:offset+4]
        """
        self. offset = offset
        #we're using singles so stride is 4        
        self.stride = 4
        #slice containing an int32 with the length of the array
        self.lenSlice = slice(self.offset,self.offset+4)
        
    def __call__(self, fname):
        """
        open the file fname and read the frame
        """
        self.fname = fname
        with open(self.fname,'rb') as f:
            self.byteArray = f.read()
        L = self.getLen()
        self.frame = zeros(3, L)
        idx =self.offset+4  
        
        x = zeros(3*L)
        for i in range(L*3):
            x[i] = self.hexStr2float(self.byteArray[slice(idx, idx+self.stride)])
            idx+= self.stride
        self.frame = x.reshape(3,L)
        return self.frame.copy()
            
    def getLen(self):
        """
        return the length of the array as read from the binary
        """
        return self.intArr2int32(self.hexStr2intArr(self.byteArray[self.lenSlice]))
        
    def hexStr2float(self, harr):
        """
        convert a string of 4 bytes to a float32
        """
        return struct.unpack('!f', harr)[0]
        
    def hexStr2intArr(self, harr):
        """
        convert a string of 4 bytes to an array of int32s between 0 and 255
        """
        return int32([ord(h) for h in harr])
        
    def intArr2int32(self, arr):
        """
        convert from an array of 4 int8s to one int32
        """
        pw = flipud(arange(0,4))
        return sum(arr*256**pw)
 
class scanFolderReader(object):
    """
    this class loops over the files in a 'scan folder' which is one scan through
    the z-coordinate for one location on the array and reads in the images there
    """
    def __init__(self, frameReader = None):
        if frameReader ==None:
            self.FR = colorFrameReader()

    def __call__(self, scanFolder):
        self.scanFolder = scanFolder
        
        #read the nominal array position from the folder name        
        self.parseScanFolderName()
        
        #initialize lists
        self.scanFiles = []
        z = []
        F = []
        for fname in  os.walk(self.scanFolder).next()[2]:
            #check for -00 in case there are stray files in the directory             
            if fname[-3:]!='-00':
                continue
            self.scanFiles+=[fname]
            z += [self.parseScanName(fname)]
            F += [self.FR(self.scanFileName(fname))]
        
        #the nominal z pos
        self.z = float64(z)
        self.frames = zeros((len(F), 3, F[0].shape[1]))
        IDX = argsort(z)
        self.z = self.z[IDX]
        for i,idx  in enumerate(IDX):
            self.frames[i] = F[idx]
        
        return (self.nominalArrayPosition, self.z.copy(), self.frames.copy())
        
    def scanFileName(self, fname):
        return os.path.join(self.scanFolder, fname)
        
    def parseScanFolderName(self):
        self.nominalArrayPosition = int32(os.path.split(self.scanFolder)[1])
        
    def parseScanName(self, fname):
                
        nominalArrayPosition = int32(fname[:2])
        sledPosition = float64(fname[3:8])
        assert nominalArrayPosition == self.nominalArrayPosition, "nominal array position mismatch between file and folder, please verify that the files are in the correct folder"
        return sledPosition
        
        
class lensFolderReader(object):
    """
    this is the big boy. it reads folders (one per x location) of folders 
    (one per z location) and builds up an image for each (x,z)  pair
    """
    def __init__(self, lensFolder, scanReader = None):
        self.lensFolder = lensFolder
        if scanReader == None:
            self.SR = scanFolderReader()
            
        scanFolders = os.walk(self.lensFolder).next()[1]
		
		self.scanFolders = []
		for f in scanFolders:
			if os.path.split(f)[1][0]=='.':
				continue
			else:
				self.scanFolders+=[f]
				
        self.cnt = 0

    def __iter__(self):
        return self
            
    def next(self):
        if not self.checkIteration():
            self.cnt = 0
            raise StopIteration
        
        scanFolder = self.scanFolders[self.cnt]
        self.cnt+=1

        #nominal x position (scalar), Z positions (1d array of len nZ), 
        #and stepImages (nZ x 3 x cameraImageLength array ) for the current scan folder
        x, Z, stepImages = self.SR(os.path.join(self.lensFolder, scanFolder))

        return (x, Z, stepImages)
            
    def checkIteration(self):
         return self.cnt<len(self.scanFolders)
        
            
   
            
if __name__ == '__main__':
    fname = 'blah.hex'
    scanFolder = os.path.join(os.path.curdir, "Initial Scans", "00")
    lensFolder = os.path.join(os.path.curdir, "Initial Scans")

#    print "reading single frame"
#    FR = colorFrameReader()
#    y = FR(fname)
#    plot(y.T)
#    
    print "reading z-scans at one array location"
    SR = scanFolderReader()
    x,z,Y = SR(scanFolder)
    figure()
    for yy in Y:
        plot(yy[0])
    
#    print "reading all scans in succession"
#    LR = lensFolderReader(lensFolder)
#    YY = [L for L in LR]
    
#    figure()
#    for yy in YY[1][-1]:
#        plot(yy[0])
    
    

    