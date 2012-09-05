from numpy import *
import pdb
from helpers import virtualWarning
import time
#define abstract container classes for models and collections
#of models

class model(object):
    """
    base class for providing interfaces to models
    """
    def __init__(self):
        """
        construct the model
        """
        virtualWarning()

    def parameterVector(self):
        """
        return an array of parameters that fully specify the model
        """
        virtualWarning()
    def parameterKeys(self):
        """
        list the names of the parameters in the order that
        parameterVector keeps them
        """
        virtualWarning()
        
    def forwardSolve(self):
        """
        run the model, computing the forward solution
        """
        virtualWarning()
        return zeros(1)
    def summarize(self):   
            virtualWarning()


def solve(M, **kwargs):
    M.forwardSolve(**kwargs)
    print 'solved!'
    
class modelEnsemble(object):
    """
    base class for collections of models
    """
    def __init__(self):
        """
        set up the object
        """
        virtualWarning()
    
    def ensembleSolve(self,**kwargs):
        """
        solve all the models in the ensemble
        default behavior is looping over the model array, calling
        their forwardSolve() method
        """
        print 'solving ensemble...'
        cnt = 0
        bigStart = time.time()
        for idx, M in ndenumerate(self.modelArray):
            if cnt%500==0:
                print 'solved %d of %d in %d min'%(cnt, self.modelArray.size, (time.time()-bigStart)/60)
            start = time.time()
            M.forwardSolve(**kwargs)
            self.solutionTime[idx] = time.time()-start
            cnt+=1
    
    def summary(self):
        """
        print/plot basics about the ensemble
        """
        virtualWarning()
    
    def modelConstructor(self, parameterVector):
        """
        return a model given a parameter vector
        """
        virtualWarning()    
    def parameterGenerator(self, shape = (1,)):
        """
        generate an array of parameter vectors upon request
        """
        virtualWarning()
        
    def generate(self, **kwargs):
        """
        make an array of model parameter vectors,
        then call construct to build an array of models
        """
        self.parameterVectors = self.parameterGenerator(**kwargs)
        pv = self.parameterVectors
        
        self.modelArray = empty(pv.shape[:-1],dtype(object))
        
        for i in arange(self.modelArray.size):
            idx = unravel_index(i,self.modelArray.shape)
            try:
                self.modelArray[idx] = self.modelConstructor(pv[idx])
            except:
                pdb.set_trace()
                
        self.solutionTime = zeros_like(self.modelArray)
    def arrayGet(self, selectionTuple = None, getFunList = [lambda m_: m_], flatten = False):
        """
        apply getFun() to all models in self.modelArray[selectionTuple]
        """
        
        #use the whole array if no selection tuple is passed
        if selectionTuple ==None:
            M = self.modelArray                
        else:
            if flatten:
                M = atleast_1d(self.modelArray.ravel()[selectionTuple])
            else:
                M = atleast_1d(self.modelArray[selectionTuple])

        
        
        x = [array(zeros(M.shape), dtype = type(fn(M.ravel()[0]))) for fn in getFunList]
        for idx, m, in ndenumerate(M):
            cnt = 0
            for fn in getFunList:
                x[cnt][idx] = fn(m)
                cnt+=1
        
        return tuple(x)
        
        
        
        
        
        