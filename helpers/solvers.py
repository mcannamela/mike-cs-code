from numpy import *
from helpers import *

#print "importing mlabMacros"
#
#import mlabMacros as mlm
import pdb
import numpy as np
import scipy
import sys
scipy.pkgload('weave')
from scipy import weave
print "import complete"

def coeffs2mat(c):
    """
    make a full matrix from an array of the 3 diagonals
    """
    return diag(c[0,0,1:],-1)+diag(c[0,1])+diag(c[0,2,:-1],1)


class oneDMesh(object):
    """
    a base class for a one d mesh of nodes and interfaces
    """
    def __init__(s, x0 = float32(.5), xN = float32(5.5), nNodes = 5):
        """
        initialize the mesh, computing node/interface spacing and geometric quantities
        """
        if type(x0) == type(float32(0)):
            s.cast = s.toSingle
        elif type(x0) == type(float64(0)):
            s.cast = s.toDouble
        else:
            s.cast = s.noCast
        
        #basic linear mesh
        s.nNodes = int(nNodes)
        
        #neighboring nodes of all inner interfaces
        s.interfaceNeighborNodes = vstack((arange(s.nNodes-1), arange(1,s.nNodes)))
        
        #create the node/interface location arrays
        s.xInt, s.xNode = s.initNodes(x0,xN)
        s.dxInt = diff(s.xInt)
        s.dxNode = diff(s.xNode)
        
        s.xInt = s.cast(s.xInt)
        s.xNode = s.cast(s.xNode)
        
        #area and volume
        s.interfaceArea = s.area() 
        s.nodeVolume = s.volume()  
        s.totalVolume = sum(s.nodeVolume)
    

    def __getstate__(self):
        """
        need to exclude all lambda objects
        """
        black = ['cast']
        S = []
        K = []
        for k in self.__dict__.keys():
            if k in black:
                continue
            K+= [k]
            S+= [self.__dict__[k]]
            
        return dict(zip(K,S))
        
    def __setstate__(self, D):
        """
        gonna have to call the methods that reconstruct the lambdas, which cannot be pickled
        """
        self.__dict__ = D.copy()
        
        if type(self.xNode) == type(float32(0)):
            self.cast = self.toSingle
        elif type(self.xNode) == type(float64(0)):
            self.cast = self.toDouble
        else:
            self.cast = self.noCast
            
        
    def initNodes(s, x0,xN):
        """
        make the node and interface spacings
        """
        virtualWarning()
        print '\n return a tuple of arrays, first of interface locations then of node locations'
        return (linspace(x0,xN, s.nNodes), x0+diff(linspace(x0, xN,s.nNodes))/2)
        
    def area(s, idx = None):
        """
        area of the interfaces
        virtual
        """
        virtualWarning()
        print '\n return an array of interface areas'
    
    def volume(s, idx = None):
        """
        volume of the nodes
        """
        virtualWarning()
        print '\n return an array of node volumes'
        
    def interfaceProperty(s, nodeProps, boundaryInterfaceProps):
        """
        given properties at the nodes, compute harmonic mean properties at the interface
        """
        return hstack((boundaryInterfaceProps[0], 
                                    harmonicMean(nodeProps[s.interfaceNeighborNodes], s.interfaceNeighborWeights()),
                                    boundaryInterfaceProps[1]))#nNodes
        
    def interfaceNeighborWeights(s):
        """
        weight each neighboring node by the distance of the OTHER node to the interface, divided by the internodal distance
        """
        xNeigh = s.xNode[s.interfaceNeighborNodes] #nNodes-1
        dxNeigh = (xNeigh[1]-xNeigh[0]) #nNodes-1
        dxInt = diff(s.xInt)#nNodes
        return vstack((dxInt[1:]/(2*dxNeigh), dxInt[:-1]/(2*dxNeigh)))#nNodes-1
        
    def weightedAdjacencyMatrix(s, interfaceWeights, BCWeights = None):
        """
        use geometry and the supplied interface weights to get weights between each node and its neighbors, as 
        well as between a node and itself. boundary nodes return a weight to their respective boundary interfaces
        """
        #print "this has changed"
        if BCWeights == None:
            #must return the tuple (backwardWeights, forwardWeights)
            internodalWeights = s.internodalWeight(hstack((s.xInt[0], s.xNode,s.xInt[-1])), interfaceWeights)
            
            #must return the array selfWeights
            selfWeights = s.selfWeight(internodalWeights)
            
            return vstack((-internodalWeights[0], selfWeights,-internodalWeights[1]))
        
        bcR = []
        for b in BCWeights:
            if b==0.:
                bcR+=[Inf]
            else:
                bcR+=[1/b]
            
        
        bR = (1/s.xInt[:-1] - 1/s.xNode)*(4*pi*interfaceWeights)**-1
        fR = (1/s.xNode     - 1/s.xInt[1:])*(4*pi*interfaceWeights)**-1
        forwardWeight = (fR+hstack((bR[1:],bcR[1])))**-1
        backwardWeight = (bR+ hstack((bcR[0],fR[:-1])))**-1
#        print "bR: "+str(bR)
#        print "fR: "+str(fR)
#        print "fwdWeight: "+str(forwardWeight)
#        print "backWeight: "+str(backwardWeight)
        
        return vstack((-backwardWeight, forwardWeight+backwardWeight, -forwardWeight))
    
    def internodalWeight(s):
        """
        weights between neighbor nodes
        virtual
        """
        virtualWarning()
        print '\n return a tuple of 1d arrays (backwardWeights, forwardWeights)'
        return (zeros(s.nNodes), zeros(s.nNodes))
    def selfWeight(s, internodalWeights):
        """
        weight between the node and itself
        """
        virtualWarning()
        print '\n return an array of one weight per node'
        return ones_like(interNodalWeights)
    
    def toSingle(self,c):
        return float32(c)
        
    def toDouble(self,c):
        return float64(c)
        
    def noCast(self,c):
        return c

class marchingSolver(object):
    """
    abstract base class
    solve a 1-d line of a finite difference equation with nonlinear coefficients, and march this solution one way in time or 
    space
    """
    def __init__(self, IC = None, nNodes = 5, nVars = 1, verbose = False, sema = None):
        """
        initialize data members that every solver must have. subclasses should write their own init methods to take care
        of initial conditions
        """
        self.verbose = verbose
        if IC == None:
            IC = zeros(nVars,nNodes)
        
        #number of nodes, number of variables
        self.nNodes = nNodes
        self.nVars = nVars
        
        #in case we want to keep track of time
        self.currentStep = 0
        
        #right hand side, the known side of the equation
        self.rhs = zeros((nVars,nNodes))
        
        #tridiagonal matrix of coefficients multiplying the unknown vector
        self.coeffs = zeros((nVars,3, nNodes))
        
        #the current and previous step solutions
        self.soln = IC.copy()
        self.prevIterSoln = IC.copy()
        
        if self.soln.ndim ==1:
            self.soln = expand_dims(self.soln,0)
            self.prevIterSoln = expand_dims(self.prevIterSoln,0)
            
        #steps we've already solved
        self.solvedSteps = (self.soln.copy(),)
        self.bcHistory = [] #provide a place to record the BC's
        
        #iteration parameters and counters
        self.itercount = [0.]#how many iterations at each step
        self.E = [0] #keep track of step error
        self.stepErrorThres = .002#threshold for terminating iteration
        self.maxIter = 100#timeout on a step after this many iterations
        self.alpha = 1#underrelaxation parameter
        
        self.sema = sema #semaphore, if used
        
        
    def nSteps(self):
        return len(self.itercount)
        
    def reset(self):
        """
        return to initial conditions and prep to start the sim over again
        """
        self.currentStep = 0
        self.solvedSteps = (self.solvedSteps[0].copy(),)
        self.soln = self.solvedSteps[0].copy()
        self.prevIterSoln = self.solvedSteps[0].copy()
        self.itercount = [0.]#how many iterations at each step
        self.E = [0] #keep track of step error
        self.bcHistory = []
        
    def solve(self):
        """
        solve steps then advance, marching through the domain until the stopping criterion is met.
        """
        #lock the semaphore
        if self.sema!=None:
            self.sema.acquire()
            
        self.updateBCs()
        while self.stoppingCriterion()!=True:
            self.solveStep()
            if self.verbose:
                print "the current solution is: "+str(self.soln)
            self.advance()
        
        self.E = array(self.E)
        if self.verbose:
            print 'mean iterations: ' +str(self.itercount/self.nSteps())
        
        #unlock the semaphore
        if self.sema!=None:
            self.sema.release()
        
        return array(self.solvedSteps)
        
    def solveStep(self):
        """
        iteratively solve the current step until it is self consistent
        """
        
        i = 0
        while self.stepError()>self.stepErrorThres and i<self.maxIter or i==0:
            self.updateCoeffs()
            self.prevIterSoln = self.soln.copy()
            for j in arange(self.nVars):
                deltaT = (self.tdmaSolve(self.coeffs[j], self.rhs[j])-self.soln[j])
                self.soln[j] += self.alpha*deltaT
                #print str(mean(deltaT))+'   ' +str(np.max(deltaT))
                
            i+=1
            if i>=self.maxIter:
                print 'iteration limit reached'
            
        self.itercount+=[i]
                
        
    def advance(self, timestep = 1):
        """
        boot the current solution to the solved steps list, obtain new boundary conditions and update all coeffs and rhs accordingly
        """
        self.solvedSteps+= (self.soln.copy(),)
        self.updateBCs()
        self.currentStep+=1
        if self.verbose:
            print "advance at "+str(self.currentStep)
    
    def tdmaSolve(self, coeffs,rhs):
        """
        use the tridiagonal matrix algorithm to solve the equations  given by coefficients in coeffs, and right hand side in rhs
        """
               
        a = coeffs[0].copy()
        b = coeffs[1].copy()
        c = coeffs[2].copy()
        
        d = rhs.copy()
        n = a.shape[0]
        x = zeros_like(a)
        
#        inlineOn = True
        try:
#            if inlineOn:
            code = """
                        c[0]/=b[0];
                        d[0]/=b[0];
                        for(int i=1; i<(int(n)-1);i++){
                            c[i]/=b[i]-c[i-1]*a[i];
                            d[i] = (d[i]-d[i-1]*a[i])/(b[i] -c[i-1]*a[i]);
                            }
                        
                        x[int(n)-1] = (d[int(n)-1]-d[int(n)-2]*a[int(n)-1])/(b[int(n)-1] -c[int(n)-2]*a[int(n)-1]);
                      
                        for(int i=(int(n)-2); i>=0;i--){
                            x[i] = d[i]-c[i]*x[i+1];
                            }
                            
                       return_val = 0;
                        """
            
            
            if sys.platform.find('linux')!=-1:
                retval = weave.inline(code, ['a', 'b', 'c', 'd', 'x', 'n'])
            else:
                try:
                    retval = weave.inline(code, ['a', 'b', 'c', 'd', 'x', 'n'], compiler = 'msvc')
                except:
                    retval = weave.inline(code, ['a', 'b', 'c', 'd', 'x', 'n'], compiler = 'mingw32')
        except:
            print "weave inline failed, we'll be doing it the slow way...pure python"
            #else:
                   
            c[0]/=b[0]
            d[0]/=b[0]
            
            for i in range(1, a.shape[0]-1):
                c[i]/=b[i]-c[i-1]*a[i]
                d[i]= (d[i]-d[i-1]*a[i])/(b[i] -c[i-1]*a[i])
                
            x[-1] = (d[-1]-d[-2]*a[-1])/(b[-1] -c[-2]*a[-1])
            for i in range(a.shape[0]-2,-1,-1):
                x[i] = d[i]-c[i]*x[i+1]
       
        return x
        
    def stepError(self):
        """
        compute an error between the last solution and current solution for this step
        """
        e = np.max(sum(abs(self.soln-self.prevIterSoln), axis = 0)/sum(abs(self.prevIterSoln), axis = 0))
        if e!=self.E[-1]:
            self.E+=[e]
        
        return e
    
    def stoppingCriterion(self):
        """
        return true if the simulation is finished
        """
        virtualWarning()
        print '\n return true if the simulation should be stopped'
        return True
        
    def updateCoeffs(self):
        """
        use the current and solved steps to generate new coefficients and rhs for obtaining the next step
        """
        virtualWarning()
        print '\n must assign new values to self.coeffs'
        
    def updateBCs(self):
        """
        obtain new boundary conditions for solving the next step
        """
        virtualWarning()
        print '\n must set appropriate values for the first and last equations'
    
class basicOneDMesh(oneDMesh):
    def initNodes(s, x0, xN):
        xInt = linspace(x0,xN, s.nNodes+1)
        xNode = xInt[:-1]+diff(xInt/2)
        return (xInt, xNode)
    def area(s):
        return ones_like(s.xInt)
    def volume(s):
        return ones_like(s.xNode)
    def internodalWeight(s, xNode, interfaceWeight):
        w = diff(xNode)*interfaceWeight
        return (w[:-1], w[1:])
    def selfWeight(s, internodalWeights):
        return internodalWeights[0]+internodalWeights[1]
    
class radial1DMesh(oneDMesh):
        """
        a one dimensional mesh of spherical shells
        """
        def initNodes(s, x0,xN):
            """
            equal volume shells with nodes in the middle of each shell
            equal radius shells with nodes in the middle of each shell
            x0 - inner diameter
            xN - outer diameter
            """
            #interfaces, equal volume nodes
            #nextInterface = lambda r_: pow(  ((xN/2)**3-(x0/2)**3)/(s.nNodes)+r_**3, 1./3.)
            nextInterface = lambda r_: r_+(xN/2-x0/2)/s.nNodes
            
            xInt = zeros(s.nNodes+1)
            xInt[0] = x0/2
            for i in range(1,s.nNodes+1):
                xInt[i] = nextInterface(xInt[i-1])
            
            #nodes halfway between interfaces
            xNode = zeros(s.nNodes)
            xNode = xInt[:-1]+diff(xInt)/2
            
            return (xInt, xNode)
            
        def area(s, idx = slice(None, None)):
            return 4*pi*s.xInt[idx]**2
            
        def volume(s,idx = None):
            if idx==None:
                v = abs((4./3)*pi*(s.xInt[:-1]**3-s.xInt[1:]**3))
            elif idx<1:
                v = abs((4./3)*pi*(s.xInt[idx-1]**3-s.xInt[idx]**3))
            else:
                v = abs((4./3)*pi*(s.xInt[idx]**3-s.xInt[idx+1]**3))
            return v
            
        def internodalWeight(s,xNode,  interfaceWeights):
            """
            weighting between nodes that results from geometry alone
            """
            w = ( abs(diff(1/xNode))**-1 )*(4*pi)*interfaceWeights
            
            return (w[:-1], w[1:])
            
        def selfWeight(s,internodalWeights):
            """
            self weight is the sum of the internodal weights
            """
            return internodalWeights[0]+internodalWeights[1]
            
class basicSolver(marchingSolver):
    def __init__(self, IC = 5, mesh = basicOneDMesh() , maxIter = 10):
        self.maxIter = 10
        self.mesh = mesh
        marchingSolver.__init__(self,IC*ones(mesh.nNodes),  nNodes = mesh.nNodes)
                
    def stoppingCriterion(self):
        return self.currentStep>50
        
    def updateCoeffs(self):
        self.coeffs = self.mesh.weightedAdjacencyMatrix(self.mesh.interfaceProperty(.1*ones(self.mesh.nNodes), array([0,0])))
        self.coeffs[0,0] = 0
        self.coeffs[1]+=1
        self.rhs = self.solvedSteps[-1][0].copy()
        self.rhs[-1]-=self.coeffs[2,-1]*self.BC
        
        #coeffs and rhs must be 3-d
        self.coeffs = expand_dims(self.coeffs,0)
        self.rhs = expand_dims(self.rhs,0)
        
    def updateBCs(self):
        self.BC = 10
        
    def plotSoln(self):
        ss = array(self.solvedSteps)[:,0,:]
        x = self.mesh.xNode
        y = arange(ss.shape[0])
#        f = mlm.fig('theAnswer')
#        mlm.surf3(y,x ,ss,f)
#        mlm.axe([ 'step','x',  'solution'], mlm.rgs(y,x,ss),f)
#        mlm.out()
#        return f
        return None

    
if __name__== "__main__":
    print "init solver..."
    BS = basicSolver(arange(5))
    print '...init complete'
    print 'solving...'
    BS.solve()
    print 'complete'
#    f = BS.plotSoln()
#    mlm.ml.show()
    
    
    
    
