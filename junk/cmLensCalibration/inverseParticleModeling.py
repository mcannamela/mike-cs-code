from forwardParticleModeling import *
import os
try:
    from scipy.spatial import cKDTree
    canHasKDTree = True
except:
    canHasKDTree = False
    
canHasKDTree = False

#if os.name == 'posix':
#    import mlabMacros as mlm

class forwardMapper(object):
    """
    for a given set of scans (x and z positions), calling forwardMapper on 
    a grid of particle states (D and z) will compute the blurry images and 
    observed rgb diameters for every (D,z,x) tuple in the grid
    
    the return value is a forwardMapPostProcessor object that holds all the 
    computed values and has facilities for writing the result to disk and 
    transforming the observed diameters to other parameters, if necessary, e.g. 
    difference of observed diamters between colors
    """
    def __init__(self, lensFolder, verbose = False):
        """
        initialize the forward mapper for a given folder of scan positions
        """
        self.lensFolder = lensFolder        
        self.verbose = verbose
        self.LFR = lensFolderReader(lensFolder)
        self.KE = kernelExtractor()
        
        self.sPrint("    will now read the lens folders... ")
            
        self.initStepImages()
                
        self.sPrint("    ...extracting kernels...")
            
        self.initKernels()
                
        self.sPrint("done")
        
        #make a blurrifier for each x location
        self.initBlurrifiers()
        
    def __call__(self, PG):
        """
        for each particle in PG, for each x-position, 
        obtain the appropriate kernel for that particle's z-position and 
        use it to produce a blurry image. then measure the observed diamters 
        from the blurry image
        """
        
        nX = len(self.X)

        #dimension of the arrays will be the (numberOfXPositions, numberOfDiameters, 
        #                                      numberOfZPositions)
        N = (nX,)+PG.zeroArr().shape        
        
        #IM[i,j,k] holds the blurry image at x[i], PG.D()[j], PG.z()[k]        
        #each entry of IM is itself a 3 x nImage[i,j,k] array
        IM = zeros(N, dtype = 'object')
        
        #like IM but for the kernels. each entry of K is a 3 x nKern[i,j,k] array
        K = zeros(N, dtype = 'object')
        
        #like IM but for the measured rgb diameters. D_obs[i,j,k] gives a length
        # 3 array of the observed diamters there
        D_obs = float32(zeros(N+(3,)))
        
        self.sPrint("blurrifying")
        start = time.time()
        cnt = 0
        for k in range(nX):
            for i,j,z,im in PG:
                IM[k,i,j], D_obs[k,i,j], K[k,i,j] = self.B[k](im, z)
                cnt += 1
        elapsed = time.time()-start
        self.sPrint("done. time was %d s for %d images"%(elapsed, cnt))
        self.sPrint("that's %.2e s/image!"%(elapsed/double(cnt)))
        
        #return value spawns a post processor for the forward map
        return forwardMapPostProcessor(IM, D_obs, K, PG, array([mean(x) for x in self.X]) )
    
        
    def initStepImages(self):
        
        #list of nominal x values
        self.Xraw = []
        
        #list (for each x) of arrays of z positions for the step images
        self.Zraw = []
        
        #list (for each x) of lists (for each z) of arrays holding the raw step images
        self.stepIM = []
        
        #for every scan folder, append to the master lists. lensFolderReader objects
        #have a next method, so we can use it like an iterator
        for T in self.LFR:
            self.Xraw+=[T[0]]
            self.Zraw+=[T[1]]
            self.stepIM+=[float32(T[2])]
            
            
    def initKernels(self):
        #list (for each x) of lists (for each z) of arrays holding the kernels
        self.K = []
        
        #list (for each x) of lists (for each z) of the pixel location of that kernel's
        #peak. mean(X[i]) should give the x position of the step for a particular scan
        #folder (scan of z-positions)
        self.X = []
        
        #for each x position (scan folder)
        for i,IM in enumerate(self.stepIM):
            #list of (stepPixel, fwhm, kernelArray) tuples for each z
            T = [self.KE(im) for im in IM]

            #append list of kernels to the master list
            self.K += [[float32(t[-1]) for t in T]]
            
            #append list of specific x-positions to the master list
            self.X += [int32([t[0] for t in T])]
            
        #list (for each x) of arrays (for each z) of kernel variances (scalar)
        self.V = []
        for K in self.K:
            self.V+=[self.kernelVariance(K)]
        
        #list (for each x) of arrays (for each z) of corrected z values based on
        #minimum kernel variance
        self.Z = []
        self.Zb = []
        for i,v in enumerate(self.V):
            #smooth the kernel variance curves by nsm before finding the min            
            #the minimum variance must be roughly in the middle of the array, 
            #otherwise, THIS WILL NOT WORK!
            nsm = 50.

            #smoothed red variance                
            smv = convolve(v[0].copy(), ones(nsm)/nsm, mode = 'same')
            
            #smoothed blue variance
            smv_b = convolve(v[2].copy(), ones(nsm)/nsm, mode = 'same')
            
            #correct for boundary effects
            smv[:ceil(nsm*.5)] = amax(v[0])
            smv[-ceil(nsm*.5):] = amax(v[0])
            
            smv_b[:ceil(nsm*.5)] = amax(v[2])
            smv_b[-ceil(nsm*.5):] = amax(v[2])
            
            #find the point of minimum (red) variance and call that 0
            z0Idx = argmin(smv)
            zCorrected = (self.Zraw[i] - self.Zraw[i][z0Idx])
            
            #find point of minumum blue variance
            zBIdx = argmin(smv_b)
            bOffset = zCorrected[zBIdx]
            
            #append to the list of z positions
            self.Z += [zCorrected]
            
            #append to list of blue offsets
            self.Zb+=[bOffset]

    
    def initBlurrifiers(self):
        #list (for each x) of blurrifiers
        self.B = []
        for i in range(len(self.X)):
            self.B+= [blurrifier(self.K[i], self.Z[i])]
    
    def zLim(self, offset = .5):
        """
        later we will be interested in the z limits of the dataset. 
        this function returns those limits, which are
        offset mm less than the red focal plane, and offset mm greater than 
        the blue focal plane
        
        the red focal plane is, by definition, at z = 0, hence the lower limit
        will always be -offset mm

        """
        zb = [-offset,max(self.Zb)+offset]
        return zb
            
    def kernelVariance(self,K):
        """
        compute the variance of each kernel (per z location)
        
        K - a list of kernels, one per z location
        
        returns:
            v - array (one per z), the width of each kernel, taken to be its variance
        """
        v = zeros((3, len(K)))
        for j, k in enumerate(K):
            for i in range(K[0].shape[0]):
                y = k[i]/trapz(k[i]) 
                x = linspace(-k[i].shape[0]/2., k[i].shape[0]/2., k[i].shape[0])
                mu = trapz(x*y)
                v[i,j] = 6*sqrt(trapz(((x-mu)**2)*y))
        return v
    def sPrint(self, msg):
        if self.verbose:
            print msg
        
class forwardMapPostProcessor(object):
    """
    this class takes a forward map as computed for a particular grid of particle 
    states and set of scans (as produced by a forwardMapper object) and writes it 
    to disk. 
    it also performs any transformation of the observed diameters to e.g. difference
    of observed diameters between colors
    """
    def __init__(self,IM, D_obs, K, PG, X ):
        """
        initialize with:
            IM - IM[i,j,k] holds the blurry image at X[i], PG.D()[j], PG.z()[k]        
                  each entry of IM is itself a 3 x nImage[i,j,k] array 
            D_obs - D_obs[i,j,k] is a length 3 array holding the rgb observed diamters
            K - K[i,j,k] is a 3 x nKern[i,j] array holding the blurring kernels
            PG - particleGrid object holding the particle states D and z
            X - array of x positions corresponding to the first dimension of IM, D_obs, K
        """
        self.IM = IM
        self.D_obs = D_obs
        self.K = K
        
        #2-d arrays are the grids of the particle states
        self.D = PG.D()
        self.Z = PG.Z()

        self.X = X
        
        self.transform()
        
    def __call__(self):
        """
        return the x position of the particles, the transformed image metrics, 
        and the true diameter in a tuple
        
        D_trans is nX x nD x nZ x 3
        
        D is nD x nZ
        """
        return (array(self.X), self.D_trans, self.D)
        
    def observedDiameters(self, flat = False):
        """
        return an nX x nD x nZ x 3 array of observed diameters if flat is False,
        otherwise return a 3 x (nX*nD*nZ) array  
        """
        if not flat:
            return self.D_obs.copy()
        else:
            dobs = zeros((3,prod(self.D_obs.shape[:-1])))
            for i in range(3):
                dobs[i] = self.D_obs[...,i].flatten().copy()
            return dobs
    def transformedDiameters(self, flat = False):
        """
        return an nX x nD x nZ x 3 array of observed diameters if flat is False,
        otherwise return a 3 x (nX*nD*nZ) array  
        """
        if not flat:
            return self.D_trans.copy()
        else:
            dtrans = zeros((3,prod(self.D_trans.shape[:-1])))
            for i in range(3):
                dtrans[i] = self.D_trans[...,i].flatten().copy()
            return dtrans

        
    def writeForwardMap(self, path):
        with open(path, 'w') as f:
             for (i,j,k), im in ndenumerate(self.IM):
                T = (i,j,k,int(mean(self.X[i])), self.Z[j,k], self.D[j,k])
                T+= tuple([self.D_obs[i,j,k,m] for m in range(3)])
                f.write('%d %d %d %d %.2e %.2e %.2e %.2e %.2e\n'%T)

                
    def transform(self):
        """
        this particular transform takes the nX x nD x nZ array D_obs which 
        holds r,g,b observed diamters in its [i,j,k]'th entry, and returns
        an array of the same size but having redObservedDiameter, green-redDifference, green-blueDifference
        as its metrics instead of r,g,b diamters
        """
        R, G, B = tuple([self.D_obs[...,i] for i in range(3)])
        self.D_trans = r_['3,4', R[...,newaxis], 
                          (G-R)[...,newaxis], 
                            (B-G)[...,newaxis]]
        
class inverseMapper(object):
    def __init__(self,x, D_trans, D):
        """
        x - vector of x positions indexing the first dimension of D_trans, D
        D_trans - nX x nD x nZ x 3 array of blurry image metrics, 
        D - nD x nZ array of true diameters corresponding to the middle two dimensions of 
            D_trans
        """
        self.X = x
        self.D = D.flatten()
        
        #pre-allocate for the self.D_trans array, which will be nX x 3 x nParticleStates
        #thus, self.D_trans[i,j] will give a vector of the j'th image metric at x[i]
        self.D_trans = zeros((D_trans.shape[0], 3, prod(D.shape)))
        for i in range(D_trans.shape[0]):
            for j in range(3):
                self.D_trans[i,j] = D_trans[i,..., j].flatten()
        
    def __call__(self, x, D_obs):
        D, e = self.invert(x,D_obs)
        return (D,e)
        
    def invert(self, x,D_obs):
        pass
    def transform(self, D_obs):
        D_trans = r_['0,1', atleast_1d(D_obs[0]), atleast_1d(D_obs[1]-D_obs[0]), atleast_1d(D_obs[2]-D_obs[1])]
        return D_trans
    
        
class brutalInverseMapper(inverseMapper):
    """
    invert the mapping (xPosition, diameter, zPosition)->(xObs, D_obs_r, Dobs_g, D_obs_b)
    by a binned but nonetheless brutal nearest neighbor search
    """
    def __init__(self, x, D_trans, D, nBins = 100):
        """
        x - vector of x positions indexing the first dimension of D_trans, D
        D_trans - nX x nD x nZ x 3 array of blurry image metrics, 
        D - nD x nZ array of true diameters corresponding to the middle two dimensions of 
            D_trans
        nBins - number of bins to use for teh red observed diameter
        """
        
        #if the difference between color diameters is greater than this, the 
        #point will be skipped
        self.D_trans_thr = 7
              
        self.nBins = nBins        
        
        self.X = array(x)
        
        #pre-allocate for the self.D_trans array, which will be nX x 3 x nParticleStates
        #thus, self.D_trans[i,j] will give a vector of the j'th image metric at x[i]
        self.D_trans = zeros((D_trans.shape[0], 3, prod(D.shape)))
        
        #list (for each x) of indirect sorts of observed red diameter
        self.ridx = []
        #for each x position        
        for i in range(D_trans.shape[0]):
            #sort D_trans by the red observed diamter
            self.ridx += [argsort(D_trans[i,...,0].flatten())]

            #for each image metric, sort indirectly using idx            
            for j in range(3):
                self.D_trans[i,j] = D_trans[i,...,j].flatten()[self.ridx[i]]


        self.D = D.flatten()
        self.bin_D_red()

    def invert(self, x, Dobs):
        Dt = self.transform(Dobs)        
        if abs(Dt[1]>self.D_trans_thr) or abs(Dt[2]>self.D_trans_thr):
            return (0, 30)
        
        #under linux, we can has KDTree, otherwise we have to use a slower method
        if canHasKDTree:
            xIdx, bidx = self.locateBin(x, Dobs[0])
            d,idx = self.binTrees[xIdx][bidx].query(Dt)
            return (self.D[idx], d)
        else:
            idx, Dtrans , D = self.sly(x, Dobs[0])
            d = self.distance(Dtrans, Dobs)
            i = argmin(d)
            return (D[i], d[i])
        
        
    def bin_D_red(self):
        """
        compute bins for the red observed diameter
        """
        self.H = []
        self.B = []
        self.binTrees = []
        self.binIdx = zeros(self.X.shape[0], dtype = 'object')
        
        #for each xPosition, perform binning on the set of observed red diameters        
        for i in range(self.X.shape[0]):
            #b contains equally spaced bins for the range of the i'th set of red diameters
            #h contains the histograms, but we don't really care about them
            h_, b_ = histogram(self.D_trans[i][0], self.nBins)
            
            h,b = self.mergeBins(h_,b_)
            
            self.H+=[h.copy()]
            self.B+=[b.copy()]
            
            #find indices where the bins would fit
            idx = searchsorted(self.D_trans[i,0],b)
            
            #list (for each x) of lists (for each bin) of indices belonging to 
            #that bin
            self.binIdx[i] = []
            
            
            binTrees = []
            for j in range(len(idx)-1):
                self.binIdx[i]+= [arange(idx[j],
                            min([idx[j+1],self.D_trans[i,0].shape[0]-1] ))]
                if canHasKDTree:
                    binTrees+= [cKDTree(self.D_trans[i][:,self.binIdx[i][-1]].T.copy())]
                    
            self.binTrees+=[[bt for bt in binTrees]]
    def mergeBins(self, H,B):
        h = H.copy()
        b = B.copy()
        while True:
            zidx = flatnonzero(h<1)
            if len(zidx)==0:
                break
            h = r_['0,1', h[:zidx[0]], h[(zidx[0]+1):]]
            b = r_['0,1', b[:zidx[0]], b[(zidx[0]+1):]]
        return h, b
            
                
        
    def locateBin(self, x, Dred):
        #find nearest xPosition        
        xIdx = argmin(abs(self.X-x))
        
        #find which red diameter bin (for the given x) we are in
        bidx = max([searchsorted(self.B[xIdx], Dred)-1,0])
        bidx = min([bidx,len(self.binIdx[xIdx])-1 ])        
        return (xIdx, bidx)
        
    def sly(self, x, Dred):
        xIdx, bidx = self.locateBin(x, Dred)
        
        #indices into the big array of red diameters
        idx = self.binIdx[xIdx][bidx]

        #return (indicesIntoD_trans, D_transInBin, DInBin)
        return (idx, self.D_trans[xIdx][:,idx], self.D[self.ridx[xIdx][idx]])
        

    def distance(self,Dtrans, Dobs):
        """
        compute the distance between a transformed 
        """
        Dt = array([Dobs[1]-Dobs[0], Dobs[2]-Dobs[1]])
        return sum((Dt[...,newaxis]-Dtrans[1:])**2, axis = 0)**.5
        

class cKDTreeInverseMapper(brutalInverseMapper):
        def __init__(self,x, D_trans, D):
            inverseMapper.__init__(self, x, D_trans, D)
            self.initKDTree()
        
        def initKDTree(self):
            self.trees = []
            for i in range(len(self.X)):
                self.trees+=[cKDTree(self.D_trans[i].T.copy())]
                
        def invert(self, x, D_obs):
            xIdx = argmin(abs(self.X-x))
            D_trans = self.transform(D_obs)
            e, idx = self.trees[xIdx].query(D_trans)
            return (self.D[idx], e)
            
        
    
class inverseMapPostProcessor(object):
    def __init__(self, X, observedDiameterGrid, recoveredDiameter, residual):
        """
        X - nX array of scan locations
        observedDiameterGrid - 3 x nR_obs x nG_obs x nB_obs grid of observed r,g,b 
                                diameters
        recoveredDiameter - nR_obs x nG_obs x nB_obs array of recovered diameter values
        residual - nR_obs x nG_obs x nB_obs array of whatever quality metric you used when 
                    when recovering the diameters. lower should be better
        """
        self.DOG = observedDiameterGrid
                
        self.D = recoveredDiameter
        if self.D.ndim == 3:
            self.D = self.D[newaxis,...]
        self.E = residual
        if self.E.ndim == 3:
            self.E = self.E[newaxis,...]
        self.X = X
                
    def reset(self,X, recoveredDiameter, residual ):
        self.D = recoveredDiameter
        if self.D.ndim == 3:
            self.D = self.D[newaxis,...]
        self.E = residual
        if self.E.ndim == 3:
            self.E = self.E[newaxis,...]
        self.X = X
  
    def writeHeader(self, fname): 
        b = ([amin(d) for d in self.DOG], [amax(d) for d in self.DOG])
        n = self.DOG[0].shape
        with open(fname+'.HEAD', 'wb') as f:
            f.write('dObsMin dObsMax dObsNPoints\n')
            pre = 'R: G: B:'.split()
            for i in range(3):
                f.write(pre[i]+' '+str(b[0][i])+' '+str(b[1][i])+' '+str(n[i]))
                if i!=2:
                    f.write('\n')
                
             
    def write(self, fname, withHeader = False):
        """
        bounce the inverse map to disc
        """
        if withHeader:
            self.writeHeader(fname)
            
        R, G, B = tuple([self.DOG[i].ravel() for i in range(self.DOG.shape[0])])
        
        for i, x in enumerate(self.X):
            with open(fname+'.%d.LUT'%x, 'wb') as f:
                for j, d in enumerate(self.D[i].ravel()):
                    f.write('%.1f %.1f %.1f %.1f %.1f'%(R[j], G[j], B[j],d, self.E[i].ravel()[j]))
                    if j!=(self.D[i].size-1):
                        f.write('\n')

        
        
if __name__ == "__main__":
#    lensFolder = os.path.join(os.path.curdir, "Initial Scans")    
    lensFolder = os.path.join(os.path.curdir, "oneScan")    
   # lensFolder = "C:\\Documents and Settings\\Matt Tuttle\\Desktop\\CMS\\Users\\Matt\\Plasma\\LensCalibration\\trash"    
    print "initializing the forward map."    
    FM = forwardMapper(lensFolder)
    zb = FM.zLim()
    PG = particleGrid(D_tup = (.1, 40, 200), Z_tup = (zb[0], zb[1], 200))

    print "begin..."
    #forwardMapPostProcessor spawned by FM
    PP = FM(PG)
    
    print "done. bouncing thrward map to disk"
    PP.writeForwardMap(os.path.join(os.path.curdir, "forwardMap"))
    print "initializing the inverse map"
    BIM = brutalInverseMapper(PP()[0], PP()[1], PP()[2])
#    BIM = kdTreeInverseMapper(PP()[0], PP()[1], PP()[2])
#    BIM = cKDTreeInverseMapper(PP()[0], PP()[1], PP()[2])
    
    #use this slice to generate the observed diamter grid
    DSly = slice(.1, 30, 100j)
    
    #grid of observed rgb diameters
    DOG = mgrid[[DSly for i in range(3)]]
    
    #preallocate for the (inverted) true diameters and error metric
    D = zeros(BIM.X.shape+DOG[0].shape)
    E = zeros(BIM.X.shape+DOG[0].shape)
    
    print "forward mapping complete. inverse map initialized. begin inversion!"
    start = time.time()
    cnt = 0
    for i,x in enumerate(BIM.X):
        for idx, blah in ndenumerate(DOG[0]):
            D[i][idx], E[i][idx] = BIM(x, DOG[(slice(None, None),)+idx])
            cnt+=1
            if mod(cnt, (prod(DOG[0].shape)*BIM.X.shape[0])/100)==0:
                print 'completed %d of %d in %d s'%(cnt, prod(DOG[0].shape)*BIM.X.shape[0], time.time()-start) 
    elapsed = time.time()-start
    print "completed %d inversions in %.2e s"%(cnt, elapsed)
    print "that's %.2e s/inversion"%(float64(elapsed)/cnt)
    
    IMPP = inverseMapPostProcessor(BIM.X, DOG, D, E)
    IMPP.write('inverse')
    figure()
    hist(E.flatten()[D.flatten()>0], 100)
    title("errors")
    
#    figure()
#    hist(D.flatten()[D.flatten()>0], 100)
#    title("recovered diameters")
    
    figure()
    dobs = r_['0,1', DOG[0].flatten(), DOG[0].flatten(), DOG[0].flatten()]
    d = D.flatten()
    plot(dobs[d[0]>0], d[0][d[0]>0], '.', alpha = .2)
    xlabel(r'$D^{obs}_{R}$, pixels')
    ylabel(r'$D^{inverted}$')
    
#    if os.name == 'posix':
#        mlf = mlm.fig('R-G-B Diameters')
#        dobs = PP.observedDiameters(flat = True)
#        mlm.scat(dobs, f=mlf, axeLab = 'R G B'.split())
#        
#        mlf = mlm.fig('transformed Observations')
#        dtrans = PP.transformedDiameters(flat = True)
#        mlm.scat(dtrans, f=mlf, axeLab = 'R G-R G-B'.split())
#          
#    
#    