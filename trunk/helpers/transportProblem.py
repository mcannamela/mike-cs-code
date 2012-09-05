##here we implement the streamlined simplex method to solve the transportation 
##problem. though the transportProblem class is general, the purpose behind it 
##is to compute the earth mover's distance between two sets of vector samples
from numpy import *
import networkx
import pdb
from scipy.cluster import vq
import scipy.stats
from pylab import *
import time


def cloud2features(cloud, nFeatures):
    """
    use vector quantization to transform a cloud of points into a vector of 
    feature weights
    
    cloud - nPts x nDim
    nFeatures - number of features
    
    returns:
        centroids - nFeatures x nDim, centroids of the feature vectors
        weights - nFeatures x 1, weight for each feature
    """
    
    #compute std dev to transform the features back into real units
    s = std(float64(cloud), axis = 0)
    
    #we must whiten the data since they may be very different in magnitude
    w = vq.whiten(cloud)
    
    assert not any(isnan(w)), 'nans detected in whitened features, check for zero variance features!'
    
    #compute the whitened centroids
    m, d = vq.kmeans(w, nFeatures)
    
    #un-whiten the centroids
    centroids = m*s
    
    if centroids.shape[0] != nFeatures:
        
        weight, codes = centroids2weights(cloud, centroids) 
        idx = argsort(weight)
        for i in idx[(centroids.shape[0]-nFeatures):]:
            newCentroids, ww = cloud2features(cloud[flatnonzero(codes==i)], 2)
            if newCentroids.shape[0]!=2:
                pdb.set_trace()
            centroids[i] = newCentroids[0].copy()
            centroids = r_['0,2', centroids, newCentroids[1].copy()]
        
    assert centroids.shape[0] !=0, 'centroids should not have zero length' 
    assert centroids.shape[0] == nFeatures, 'we should have the right number of centroids!'
    
    weights, codes = centroids2weights(cloud, centroids)
    
    return (centroids, weights)
    
def centroids2weights(cloud, centroids):
    #quantize the cloud to the centroids
    codes, d  = vq.vq(cloud, centroids)
    
    nFeatures = centroids.shape[0]
    
    #weights for each feature will be the number of samples falling 
    #in the voronoi cell of that feature
    weights = float64(zeros(nFeatures))
    for i in range(nFeatures):
            weights[i] = float64(sum(codes==i))
    return weights, codes
    
class cloud(object):
    """
    represent a point cloud and endow it with a notion of distance
    """
    def __init__(self, X, nFeatures = 10 ):
        """
        X - nPts x nDim
        nFeatures - scalar, number of centroids to find in vector quantization
        """
        
        self.X = X.copy()
        # print self.X.shape
        self.nFeatures = nFeatures
    
    def __getitem__(self, y):
        return self.X[y]
        
    def __add__(self, other):
        X = r_['0,2', self.X, other.X]
        return cloud(X, max([self.nFeatures, other.nFeatures]))
        
    def featureSequence(self, scales):
        self.whiten(scales)
        self.computeFeatures()
    
    def whiten(self, scales):
        """
        normalize each column of the cloud by scales to get a whitened cloud
        
        scales - nDim, each element is the scale for that column of X
        """
        self.s = scales.copy()
        self.Y = self.X/self.s
    
    def computeScales(self):
        """
        return the scale of each column of the cloud
        """
        return std(self.X, axis = 0)
        
    def computeFeatures(self):
        #compute the whitened centroids
        self.fWhite, d = vq.kmeans(self.Y, self.nFeatures)
        
        try:
            assert self.fWhite.shape[0]==self.nFeatures, 'too few features found'
        except:
            pdb.set_trace()
    
        #un-whiten the centroids
        self.fColor = self.fWhite*self.s
        self.w, self.codes = centroids2weights(self.Y, self.fWhite)
    
    def pointDistance(self, x,y):
        """
        the distances between two points of one of the clouds, x and y
        
        x - m x nDim
        y - n x nDim
        
        for convenience this function will return an m x n array of distances in the case that 
        x and y are themselves clouds. however, m and n should be kept small to avoid swamping memory. 
        
        """
        
        X = atleast_2d(x)
        Y = atleast_2d(y)
        
        return sum((X[:,newaxis,...]-Y[newaxis, :, ...])**2, axis = -1)**.5
    
    def distance(self, other):
        C = self.pointDistance(self.fWhite, other.fWhite)
        return sum(transportProblem()(C, float64(self.w)/sum(self.w), float64(other.w)/sum(other.w)))
    

class transportProblem(object):
    """
    the transportation problem is a special case of linear programming where 
    we seek a basis for transforming one vector into another using the least 
    amount of work. 
    
    conceptually, this is the problem of shipping units of supply from N 
    factories to M consumers. the cost to ship one unit from supplier i to 
    consumer j is C[i,j]. C is called the cost matrix. We must find flows 
    F_i_j such that all supplies are exhausted and all demand is met while 
    minimizing the total shipping cost, the sum over i and j of (F_i_j*C_i_j). 
    """
    
    def __init__(self,  maxiter = 1000, eps = 1e-10, verbose = False):
        """
        maxiter - maximum number of times to pivot
        eps - floating point numbers x with abs(x)<eps will be considered zero
        verbose - whether to print little messages along the way
        """
        self.eps = eps
        self.verbose = verbose
        self.maxiter = maxiter
        
    def __call__(self,cost, supply, demand):
        """
        main loop to compute the solution
        """
        #set the parameters of the problem
        self.reset(cost, supply, demand)
        
#        print self.c
        
        #check for special cases of one supplier or one consumer
        if self.s.shape[0]==1:
            self.f = self.d.copy()
            return sum(self.c*self.d)
        if self.d.shape[0]==1:
            self.f = self.s.copy()
            return sum(self.c*self.s)
        
        assert abs(sum(supply) - sum(demand))<self.eps, 'supply is not equal to demand! this should be a balanced transport problem'
        #construct a basic feasible solution
        self.bfs()
        
        if flatnonzero(self.f>self.eps).shape[0] ==1:
            return self.objective()
            
        
        #construct a graph of the non-zero cells ('basic' cells)
        self.buildBasicGraph()
        
        #compute the reduced cost from the graph of basic cells
        self.reducedCost()
        
        #print starting total cost
        if self.verbose:
            print 'objective after basic feasible solution is %.5e'%self.objective()
        
        #check if we need to optimize the bfs
        if self.isOptimal():
            return self.objective()
        
        #the main optimization loop. each time around the loop we make one pivot 
        #of the simplex
        cnt =0
        while not self.isOptimal():
            self.pivot()
            self.reducedCost()
            if cnt==self.maxiter:
                print 'max iterations reached, breaking main loop...'
                break
            cnt+=1
        
        #message if we were able to attain the optimum
        if self.verbose:
            print  'terminated after %d pivots'%cnt
            print 'objective after optimization is %.2e'%self.objective()
            
        return self.objective()
        
    def reset(self, cost, supply, demand):
        #cost matrix, supply and demand arrays
        ridx = flatnonzero(supply>self.eps)
        cidx = flatnonzero(demand>self.eps)
        

        self.c = cost[ridx][:,cidx]

        self.s = atleast_1d(supply[ridx])
        self.d = atleast_1d(demand[cidx])
        self.f = zeros((self.s.shape[0], self.d.shape[0]))
        self.idxArr = arange(self.f.size).reshape(self.f.shape)
        
        #alias the dims since we use them often
        nS = self.s.shape[0]
        nD = self.d.shape[0]
        
        #error checking
        # assert self.c.shape == (nS, nD), 'cost matrix must be nSupply x nDemand'        
        assert (sum(supply)-sum(demand))<=self.eps, 'supply and demand are not equal, this is a balanced transportation problem'
        
        #preallocate for duals
        self.u = zeros_like(self.s)
        self.v = zeros_like(self.d)
        
    def objective(self):
        """
        the objective function is the amount of stuff moved (f) times
        the distance through which it is moved (c)
        """
        return sum(self.c*self.f)
        
        
    def duals(self):
        """
        the dual values are defined for basic cells by the equations 
        u_i + v_j = c_i_j
        
        they are used in the computation of the reduced cost, 
        rc_i_j = c_i_j - (u_i + v_j) 
        
        we choose arbitrarily u_0 = 0 since we'll have an extra equation
        and work from there
        """
        #aliases       
        u = self.u
        v = self.v
        G = self.basicGraph
        c = self.c
        
        #there will always be nRows + nColumns -1 cells in the basis so
        #the dual equations are underconstrained. hence we can pick one value
        #arbitrarily
        u[0] = 0
        
        #keep arrays marking whether the duals have been computed
        done = (zeros_like(u), zeros_like(v))
        done[0][0] = 1
        
        #basic nodes
        B  = G.nodes()
        
        #nodes for which one dual is known, meaning the other can 
        #be computed
        available = set() #this will act as a queue during our loop
        for b in B:
            if b[0]==0:
                available.add(b)
        
        if self.verbose:
            print 'finding duals...'
        
        #iterate until we have computed all the duals u_i and v_j
        while not (all(done[0]>0) and all(done[1])>0):
            rem = []#these nodes will be removed from the queue on each loop
            toAdd = []#these will be added to the queue on each lop
            
            #loop over nodes in the queue
            for nde in (available):
                
                ##################################################
                ########## three cases are possible 
                ########## 1) we know both u and v for this node
                ########## 2) we know u for this node and can solve for v
                ########## 3) we know v for this node and can solve for u
                
                #both u and v are known for this node
                if done[0][nde[0]] and done[1][nde[1]]:
                    continue
                
                if done[0][nde[0]]:#u is known for this row
                    v[nde[1]] = c[nde]-u[nde[0]]
                    done[1][nde[1]] = 1
                            
                if done[1][nde[1]]:#v is known for this row
                    u[nde[0]] = c[nde]-v[nde[1]]
                    done[0][nde[0]] = 1 
                ###################################################
                
                #we know this node now so we remove it
                rem+=[nde] 
                
                #since we know this node, it's neighbors become available
                neigh = G.neighbors(nde)
                for n in neigh:
                        toAdd+=[n]
            
            #update the queue
            [available.add(nde) for nde in toAdd]
            [available.remove(nde) for nde in rem]
            
        
        
    def reducedCost(self):
        """
        the reduced cost is the difference between the cost and 
        the sum of the duals for that cell
        
        it represents the change in the objective function that can be obtained
        by bringing that each cell into the basis
        """
        self.duals()
        self.RC = self.c-(self.u[...,newaxis]+self.v)
        
    def isOptimal(self):
        """
        if all reduced costs are positive, then the solution must be optimal
        """
        return all(self.RC>=-self.eps)
    
    def isConnected(self, n,m):
        """
        test for connectivity between two matrix entries (2-tuples) n
        and m
        """
        return n!=m and (n[0]==m[0] or n[1]==m[1])
        
    def augmentGraph(self,  G, n):
        """
        add a matrix entry to the the graph G as a new node
        and hook up appropriate edges
        """
        G.add_node(n)
        for m in G.nodes():
            if self.isConnected(m,n):
                G.add_edge(m,n)
                
        
    def pivot(self):
        """
        perform one pivot of the simplex method
        
        we choose the empty cell ('non-basic') with the most negative 
        reduced cost. this is called the entering cell, since we are 
        bringing it into the solution. 
        
        for network theoretic reasons, there is guaranteed to be exactly
        one loop through the entering cell that we can use to compensate 
        for the added flow to that cell. 
        
        by inspecting the loop, we can determine how much flow is allowed 
        to place in the entering cell, and also which cell in the loop will 
        leave the basis. the cell that leaves the basis (i.e. it's flow becomes 
        zero) is called the exiting cell. 
        """
        entering = tuple(unravel_index(argmin(self.RC), self.f.shape))
        G = self.basicGraph
        self.augmentGraph(G, entering)
        cb = networkx.cycle_basis(G)
        
        
        L = []
        for c in cb:
            good = array([False, False])
            
            #must have an even number of cells in the loop
            if len(c)%2!=0:
                continue
                
            #the entering variable must appear in the loop
            if all([n!=entering for n in c]):
                continue
            
            #we must have balancing members of the loop in the row
            #and column of the entering variable
            for n in c:
                for i in range(2):
                    if n!=entering and n[i]==entering[i]:
                        good[i] = True
            if all(good):
                L = c
                break

            assert len(L)>=4, 'loop not found or loop too short'

        
        
        idx = flatnonzero(array([n == entering for n in L]))
       # print idx
        try:
            reIdx = [i%len(L) for i in range(idx[0], idx[0]+len(L))]
        except IndexError:
            print idx
            pdb.set_trace()
        L = [L[i] for i in reIdx]
        #print L
        even = [L[i] for i in range(0, len(L), 2)]
        odd = [L[i] for i in range(1, len(L), 2)]
        
        pivot = odd[argmin([self.f[n] for n in odd])]
        delta = self.f[pivot]
        
        for n in even:
            self.f[n]+=delta
        for n in odd:
            self.f[n]-=delta          
        
        
        zeroAndOdd = [odd[m] for m in 
                flatnonzero(array([self.f[n]<=self.eps for n in odd]))]
        if len(zeroAndOdd)>1:
            G.remove_node(zeroAndOdd[scipy.stats.randint(0, len(zeroAndOdd)).rvs(1)])
        else:
            G.remove_node(pivot)
            
        
        
    def bfs(self):
        """
        find a basic feasible solution (BFS) to the problem
        """
        C = self.c.copy()
        s = self.s.copy()
        d = self.d.copy()
        rtabu = ones_like(s)
        ctabu = ones_like(d)
        idxArr = arange(self.f.size).reshape(self.f.shape)
        
        while any(r_[rtabu, ctabu]==1):
#            print r_[rtabu, ctabu]
#            print s
#            print d
            
            
            #ignore rows and columns that are exhausted
            rt = flatnonzero(rtabu)
            ct = flatnonzero(ctabu)
            
            try:
                rp, cp = self.penalties(C[rt][:,ct])
            except:
                pdb.set_trace()
            
            #find the cell of greatest penalty and least cost
            if amax(rp)>amax(cp):#max row penalty dominates
                ridx = argmax(rp)
                cidx = argmin(C[rt][:,ct][ridx])
            else:#max column penalty dominates
                cidx = argmax(cp)
                ridx = argmin(C[rt][:,ct][:,cidx])
            
#            print (ridx, cidx)
            
            fidx = idxArr[rt][:,ct][ridx, cidx]
            ri = rt[ridx]
            ci = ct[cidx]
            
            
            
            #put the maximum possible flow into that cell and 
            #mark the according row or column as exhausted
            #pdb.set_trace()
            
            #demand exceeds supply
#            print 'making assignment'
            
            if s[ri]<d[ci]:
                self.f.ravel()[fidx] = s[ri] #set the flow to the supply
                d[ci]-=s[ri]#decrement the demand by the supply
                rtabu[ri] = 0  #this supplier is now exhausted
                
                if d[ci]<=self.eps:
                    ctabu[ci] = 0
#                print 'supplier exhausted'
                
            elif s[ri]>d[ci]:#supply exceeds demand
                self.f.ravel()[fidx] = d[ci]
                s[ri]-=d[ci]
                ctabu[ci]= 0
#                print 'consumer exhausted'
                
                if s[ri]<=self.eps:
                    rtabu[ri] = 0
                    
            elif s[ri]==d[ci]:
                self.f.ravel()[fidx] = d[ci]
                ctabu[ci]= 0
                rtabu[ri] = 0
                d[ci]-=s[ri]
                s[ri]-=d[ci]
                
#                print 'supplier and consumer exhausted'
                
        
    
    def buildBasicGraph(self):
        """
        build a graph of the basic variables.
        if the graph is not connected, find a non-basic 
        variable s.t. adding it to the basis makes the graph connected
        """
        self.basic = nonzero(self.f>self.eps)
        self.basicGraph = networkx.Graph()
        [self.basicGraph.add_node((self.basic[0][i], self.basic[1][i])) 
            for i in range(self.basic[0].shape[0])]
            
        for n in self.basicGraph.nodes():
            for m in self.basicGraph.nodes():
                if n==m:
                    continue
                elif self.basicGraph.has_edge(n,m):
                    continue
                elif n[0] == m[0] or n[1] == m[1]:
                    self.basicGraph.add_edge(n,m)
        
        
        N = self.basic[0].shape[0]
        
        #alias the dims since we use them often
        nS = self.s.shape[0]
        nD = self.d.shape[0]
        
        
        if N==(nS+nD-1):
            done = True
            if self.verbose:
                print 'network has the proper number of nodes'
        else:
            assert N<(nS+nD-1), 'too many basic cells, something went wrong'
            done = False
            print 'network has too few nodes, will now add zero cells until connected'
        
        
        
        G = self.basicGraph
        nonbasic = nonzero(self.f<=self.eps)
        nonbasic = set([(nonbasic[0][i], nonbasic[1][i]) for i in range(nonbasic[0].shape[0])])
        basic = set(G.nodes())
        while len(G.nodes())<(nS+nD-1) and not networkx.is_connected(G):
            print 'adding...'
            gg = networkx.connected_component_subgraphs(G)
            

            cnt = 1
            for n in gg[0].nodes():
                try:
                    for m in gg[1].nodes():
                        nd = [(n[0], m[1]), (n[1], m[0])]
                        for i in range(2):
                            if nd[i] in nonbasic and nd[i] not in basic:
                                self.augmentGraph(G, nd[i])
                                nonbasic.remove(nd[i])
                                basic.add(nd[i])
                                cnt-=1
                                break
                        if cnt==0:
                            break
                except:
                    pdb.set_trace()
                if cnt==0:
                        break

    
    def supplyConstraints(self):
        return sum(self.f, axis = 1)<=self.s
        
    def demandConstraints(self):
        return sum(self.f, axis = 0)<=self.d
    
    
    def penalties(self, C):
        c = C.copy()
        big = amax(c)
        smidx = argmin(c, axis = 0)
        sm = c[smidx, arange(c.shape[1])].copy()
        c[smidx,arange(c.shape[1]) ]=big
        nsm = amin(c, axis = 0)
        columnPenalty = nsm-sm
        
        
        c = C.copy()
        smidx = argmin(c, axis = 1)
        sm = c[arange(c.shape[0]),smidx ].copy()
        c[arange(c.shape[0]),smidx ]=big
        nsm = amin(c, axis = 1)
        rowPenalty = nsm-sm
        
        
#        return dict(zip(['rowPenalty', 'columnPenalty'], [rowPenalty, columnPenalty]))
        return (rowPenalty, columnPenalty)
        

class earthMoverDistance(object):
    def __init__(self, features, verbose = False, dFun = None):
        """
        features - (f1, f2), two vectors of features that we can use to 
            form the cost matrix
            
        verbose - control message printing
        
        dFun - function to compute distance between features[i] and features[j]
        """
        self.features = features
        self.verbose = verbose
        self.TP = transportProblem()
        if dFun==None:
            self.dFun = self.defaultDistance
        else:
            self.dFun = dFun
        
        self.C = self.dFun(self.features[0],
                           self.features[1])
                           
    def __call__(self, weights):
        """
        determine the earthmover's distance between two point clouds
        
        weights - (w1, w1), two vectors of weights for the features used by this emd
        cloud1 - nPts1 x nDim
        
        cloud2 - nPts2 x nDim
        
        dFun - distance between two points, defaults to 2-norm
        """
                
        self.emd = self.TP(self.C, weights[0], weights[1])
        return self.emd
    
    def plot(self):
        
        figure()
        imshow(self.C, cmap = cm.gray)
        xlabel('demand node')
        ylabel('supply node')
        
 
        
    def defaultDistance(self, p1, p2):
        """
        p1 - nPts1 x nDim
        p2 - nPts2 x nDim
        
        """
        return euclidNorm(p1,p2)

def euclidNorm(p1, p2):
    """
    p1 - nPts1 x nDim
    p2 - nPts2 x nDim
    
    """
    if p1.ndim ==1:
        return abs(p1[:,newaxis,...]-p2)
    else:
        return sum((p1[:,newaxis,...]-p2)**2, axis = -1)**.5


if __name__=='__main__':
    tpTestOn = False
    emdTestOn = True
    
    if tpTestOn:
        s = float64(arange(1,4))
        d = float64(arange(1,5))
        
        print 'supply: '+ str(s)
        print 'demand: '+ str(d)
        c = s[...,newaxis]*d
        s[0] += sum(d)-sum(s)
        TP = transportProblem(verbose = True)
        TP(c,s,d)
        print TP.f
    if emdTestOn:
        N = 10000.
        R = (scipy.stats.norm(loc = 0.).rvs( N),
             scipy.stats.norm(loc = 1.).rvs( N),
             scipy.stats.norm(loc = 2.).rvs( N))
        
        #F[0] is centroids, F[1] is weights
        F = [cloud2features(r, 10) for r in R]
        
        print 'finding EMD...'
        D = zeros((len(R),len(R)))
        start = time.time()
        for i in range(len(R)):
            for j in range(i, len(R)):
                D[i,j] = earthMoverDistance((F[i][0], F[j][0]))((F[i][1]/N, F[j][1]/N))
                D[j,i] = D[j,i]
        
        print '...done. %d s'%(time.time()-start)
        
        print 'means are 0, 1, 2, emd between [(0,0), (0,1), (1,2), (0,2)] is:'
        print (D[0,0], D[0,1], D[1,2], D[0,2])