import cPickle as pkl
from numpy import *
import os

class bigPickle(object):
    """
    breaks a big array into chunks and pickles them individually, to reduce
    memory spikes when loading. 
    """
    

    def cutPoints(self,shape, n):
        return int32(linspace(0, prod(array(shape)), n))
        
    def dump(self, arr,fname,  n = 10):
        """
        decompose array arr in to n pieces and pickle them.
        a dummy file named fname will contain the relative path 
        to a directory containing the pickles, numbered sequentially
        arr - numpy array to be pickled 
        fname - string, name of file to dump to
        n - number of pieces to pickle
        """
        self.cnt = 0;
        idx = self.cutPoints(arr.shape, n)
        with open(fname, 'w') as f:
            f.write(self.dirName(fname)+'\n')
            f.write(str(n)+'\n')
            f.write(str(arr.shape).strip('(').strip(')').replace(',', ' ' ))
        
        if not os.path.exists(self.dirName(fname)):
            os.mkdir(self.dirName(fname))
            
        for i in range(n-1):
            p = os.path.join(self.dirName(fname), self.pklNames())
            print 'dumping '+p+'...'
            with open(p, 'wb') as g:
                pkl.dump(arr.ravel()[idx[i]:idx[i+1]], g, -1)
                
                    
        
        
    def load(self,fname):
        """
        open the dummy file created by dump, unpickle all 
        sub pickles and reconstruct into an array
        """
        self.cnt = 0;
        with open(fname, 'r') as f:
            theDir = f.readline().strip('\n')
            n = int32(f.readline())
            N = int32(f.readline().split())
        
        idx = self.cutPoints(N, n)
        
        for i in range(n-1):
            p = os.path.join(theDir, self.pklNames())
            print 'loading '+p+'...'
            if i==0:
                with open(p, 'rb') as g:
                    a = pkl.load(g)
                    typ = a.dtype
                    A = zeros(tuple(N), dtype = typ)
                    A.ravel()[idx[i]:idx[i+1]] = a.copy()
            else:
                with open(p, 'rb') as g:
                    A.ravel()[idx[i]:idx[i+1]] = pkl.load(g).copy()
        
        return A
                    
    def pklNames(self):
        self.cnt+=1
        return str(self.cnt-1)+'.pkl'
        
    def dirName(self, fname):
        return fname+'.bigPickle'
    

if __name__ == '__main__':
    a = arange(1000)
    bigPkl().dump(a, 'bigPick')
    b = bigPkl().load('bigPick')
    any((b-a)!=0)
    
