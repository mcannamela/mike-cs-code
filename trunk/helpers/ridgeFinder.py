from scaleSpace import *
from matplotlib.widgets import Slider, Button

def graph2sly(g):
        """
        from a graph g, retrieve its nodes and put them in slice format
        """
        return tuple([array(g.nodes())[:,i] for i in range(3)])


def boundingBox(x, s, b):
    """
    for the set of indices x and std devs s, find the bounding box 3 
    of x plus 3 standard devs in any direction with in the confines of
    the bounds b
    
    x - nDim x nPts, each column is a vector in the set
    s - nPts, each is a standard deviation from the points in x
    b - nDim x 2, each row is (min, max) for that dimension
    """
    mn = amin(x, axis = 1)
    mx = amax(x, axis = 1)
    
    S = 3*amax(s)
    mn-= S
    mx+= S
    
    mn = amax(r_['0,2', mn, b[:,0]], axis = 0)
    mx = amin(r_['0,2', mx, b[:,1]], axis = 0)
    
    bb = zip( mn.tolist(), mx.tolist())
    

    return  bb
    
def mergeBlobs(blobNet):
    """
    merge the connected blobs in the graph into one blob object
    """
    nd = []
    for b in blobNet.nodes():
        nd+=b.bpNet.nodes()
    
    bpNet = networkx.Graph()
    bpNet.add_nodes_from(nd)
    
    b = blob(blobNet.nodes()[0].im, bpNet, blobNet.nodes()[0].allScales)
    b.computeFootprint()
    
    return b

def pickleName(name, addendum):
    return name[:-3]+addendum+'.pkl'

def summaryName(name, addendum):
    return name[:-3]+addendum+'.txt'
    
    
class blob(object):
    """
    an image with a network of connected indices and scale for each index defines
    the blob
    """
    def __init__(self, im, blobPointNetwork,scales, imCopy = False):
        """
        set the main data members 
        
        im - image to which indices in blobPointNetwork.nodes() refer
        blobPointNetwork - a networkx network of connected blob points
        """
        if imCopy:
            self.im = im.copy()
        else:
            self.im = im
            
        self.bpNet = networkx.Graph(blobPointNetwork)
        self.sly = graph2sly(self.bpNet)
        self.t = self.sly[0]
        self.x = array(self.sly[1:])
        self.allScales = scales.copy()
        self.scales = scales[self.t]
        
        self.meanSig = mean(self.scales**.5)
        
        self.maxArrSize = uint64(10.**9/4.)
        
        self.b = zeros((self.im.ndim, 2))
        self.b[:,1] = array(self.im.shape) - 1
        

    def computeFootprint(self,  maskTemp = None, gaussTemp = None):
        
        #get the bounding box for the cloud of blob points
        bb = boundingBox(self.x,self.scales**.5,self.b)

        #this is the part of the image domain 
        gSlices = tuple([slice(b[0], b[1]+1) for b in bb])

            
        g = mgrid[gSlices]
        x = float32(r_['0,2', g[0].ravel(), g[1].ravel()])
        
        if maskTemp == None:
            maskTemp = zeros_like(x[0])>0
        if gaussTemp == None:
            gaussTemp = float32(zeros_like(x[0]))
        
        nPts = self.scales.shape[0]
        maxArrSize = nPts*x.shape[0]*x.shape[1]
        
        if maxArrSize>self.maxArrSize:
            n = range(ceil(float64(maxArrSize)/self.maxArrSize))
            print 'blob has too many points! breaking it up into %d chunks of size %d'%(len(n), self.maxChunkSize)
            for i in n:
                sly = slice(i*self.maxArrSize, min([(i+1)*self.maxArrSize,nPts]) )
                mt, gt = multiGauss32(x, self.x[:,sly], self.scales[sly]**.5)

                maskTemp[mt] += 1

                gaussTemp[mt]+= gt
            
            gaussTemp /= sum(gaussTemp)
            
            idx = flatnonzero(maskTemp)
            self.density = gaussTemp[idx]
            
            
        else:
            
            (idx, self.density) = multiGauss32(x, self.x, self.scales**.5)
            
            
        self.footprintSly = (uint32(x[0][idx]), uint32(x[1][idx]))
        self.X = float32(self.footprintSly)
        
        
        self.centroid()

        self.momentOfInertia()

    
    def centroid(self):
        self.xc = sum(self.density*self.X, axis = 1)
        self.xPrime = self.X-self.xc[...,newaxis]
        return self.xc
        
    def momentOfInertia(self):
        try:
            xp = self.xPrime
        except:
            self.centroid()
            xp = self.xPrime

        diag = sum(self.density*xp**2, axis = 1)
        cross = sum(self.density*xp[0]*xp[1])
        
        self.I = r_[diag.copy(), cross.copy()]

    def isConnected(self, other):
        s1 = set(self.footprintSly[0]*self.im.shape[1]+self.footprintSly[1])
        s2 = set(other.footprintSly[0]*other.im.shape[1]+other.footprintSly[1])
        return len(s1.intersection(s2))>0
    
    
    def writeArray(self):
        meanBrightness = sum(self.density*self.im[self.footprintSly])
        arr = r_[self.xc, self.I, self.x[0].shape[0], 
                 self.X[0].shape[0], meanBrightness, 
                 self.meanSig]
        return arr
        
    def writeHeader(self):
        s = ('centroid_0 centroid_1 '+
        'moment_00 moment_11 moment_10 '+
        'nBlobPoints nFootprintPoints meanBrightness '+
        'meanBlobPointRadius')
        return s
     
    
class ridgeFinder2d(object):
    """
    find ridges in a scaleSpace2d object 
    """
    def __init__(self, ss):
        self.ss = ss
        
        self.D2 = zeros((2,)+ss.ss.shape)
        
        print "updating derivatives..."
        self.ss.updateDerivatives()
        print "finding principal curvatures..."
        self.pqDerivatives()
        
        print "quantizing directions..."
        self.quantizePrincipalCurvatureAngle()
        print "ridge finder initialized."
        
    
    def __call__(self, kind = 'bright'):
        """
        return masks that are True for ridge points
        """
        self.ridgeCriteria2d()
        if kind=='bright':
            return self.brightMask
        else:
            return self.darkMask
            
        
    def pqDerivatives(self):
        """
        compute the derivatives in the p-q coordinate frame 
        (aligned with principal curvatures)
        """
        ht = r_['0,2',self.ss.dget([2,0]).flatten(), self.ss.dget([1,1]).flatten()]
        hb = r_['0,2',self.ss.dget([1,1]).flatten(), self.ss.dget([0,2]).flatten()]
        H = r_['1', ht[:,newaxis,...], hb[:,newaxis,...]]
        
        p = zeros((ht.shape[-1],2))
        q = zeros((ht.shape[-1],2))
        k = zeros((ht.shape[-1],2))
        
        start = time.time()
        e, v = symeig2x2(H)
        
        k = e.T
        p = v[:,0].T
        q = v[:,1].T
        
        
        print 'principal curvatures complete in %d s'%((time.time()-start))
        P = zeros((2,)+ self.ss.ss.shape)
        Q = zeros((2,)+ self.ss.ss.shape)
        K = zeros((2,)+ self.ss.ss.shape)
        
        reList = [p,q,k]
        targ = [P,Q,K]
        sh = self.ss.ss.shape
        for i in range(len(reList)):
            for j in range(2):
                targ[i][j] = reList[i][:,j].reshape(sh)
       
        self.P = P.copy()
        self.Q = Q.copy()
        
        #D2[0] is greatest curvature
        self.D2[0] = K[0].copy()#in the beta direction
        self.D2[1] = K[1].copy()#in the beta - 90 direction
        

        
    def quantizePrincipalCurvatureAngle(self):
        self.pMask = abs(self.D2[0])>abs(self.D2[1])
        qMask = logical_not(self.pMask)
        ang = zeros_like(self.ss.ss)
        self.d = uint8(zeros_like(self.ss.ss))
        
        ang[self.pMask] = arctan2(self.P[1][self.pMask], self.P[0][self.pMask])
        ang[qMask] = arctan2(self.P[1][qMask], self.P[0][qMask])-pi/2
        a = array([i*pi/4 for i in range(4)]+[i*pi/4-pi for i in range(4)]+ [i*pi/4+pi for i in range(4)])
        id = r_[arange(4),arange(4),arange(4)]
        a0 = 0.
        
        self.d = uint8(id[argmin(ang[...,newaxis]-a, axis = -1)])

        self.stencils = [[[1,0],[1,1],[0,1],[1,-1]],
                        [[-1,0],[-1,-1],[0,-1],[-1,1]]]
    
    def ridgeCriteria2d(self):
        """
        compute a mask that is true for each point in scale space that satisfies 
        the ridge criteria
        find bright ridges by default.
        if dark ==True, find dark ridges
        """
        self.localStrengthMaxes()
        self.darkMask = self.sMask.copy()
        self.brightMask = self.sMask.copy()
        
        idces = transpose(nonzero(self.sMask))
        for j in range(idces.shape[0]):
            #idx = unravel_index(m, self.sMask.shape)
            idx = idces[j]
            
            #exclude points on the boundary of the image
            if any(idx==0) or any([idx[i]==(self.sMask.shape[i]-1) for i in range(len(idx))]):
                self.darkMask[tuple(idx)] = False
                self.brightMask[tuple(idx)] = False
                continue
            
            #get the directions of principal curvature 
            d = self.d[tuple(idx)]
            
            nidx = [r_[0,array(self.stencils)[i][d]]+idx for i in range(2)]
            

            D = diff(self.ss.ss[tuple(r_['0,2',nidx[0], idx, nidx[1]].T)])

            if sign(D[0])==sign(D[1]):
                self.darkMask[tuple(idx)] = False
                self.brightMask[tuple(idx)] = False
                continue
            else:
                if self.D2[1-self.pMask[tuple(idx)]][tuple(idx)]<0.:
                    self.brightMask[tuple(idx)] = True
                    self.darkMask[tuple(idx)] = False
                elif self.D2[1-self.pMask[tuple(idx)]][tuple(idx)]>0.:
                    self.brightMask[tuple(idx)] = False
                    self.darkMask[tuple(idx)] = True
                else:
                    self.brightMask[tuple(idx)] = False
                    self.darkMask[tuple(idx)] = False
            
 
        
    def ridgeStrength(self):
        dxy = self.ss.dget([1,1])
        dyy = self.ss.dget([0,2])
        dxx = self.ss.dget([2,0]) 
        
        self.s = ((dxx-dyy)**2 +4*dxy**2) #+ (256.-self.ss.ss)/15.
        del dxy, dyy, dxx
        gc.collect()
        
    def localStrengthMaxes(self):
        """
        compute a mask that is true for each pixel that is 
        a local max in ridgeStrength in the scale direction
        """
        self.ridgeStrength()
        self.sMask = zeros_like(self.s)==1
        
        prev = (self.s[0]>self.s[1])
        next = (self.s[-1]>self.s[-2])
        self.sMask[0] = prev.copy()
        self.sMask[1] = next.copy()
        for i in range(1,self.sMask.shape[0]-1):
            prev = self.s[i-1]
            next = self.s[i+1]
            self.sMask[i] = logical_and(self.s[i]>prev,self.s[i]>next)
            

class blobFinder2d(ridgeFinder2d):
    def ridgeStrength(self):
        k0 = abs(self.D2[0])
        k1 = abs(self.D2[1])
        self.s = (k0**2+k1**2)**.5

class bareBonesRidgeFinder2d(ridgeFinder2d):
    """
    retain the essentials of the ridge finder and trim away everything else 
    to reduce memory footprint
    """
    def __init__(self, rf2d):
        """
        copy the appropriate data members from rf2d
        
        rf2d - a ridgeFinder2d object that has already been run
        """
        self.ss = bareBonesScaleSpace2d(rf2d.ss)
        self.P = float32(rf2d.P.copy())
        self.D2 = float32(rf2d.D2.copy())
        self.darkMask = rf2d.darkMask.copy()
        self.brightMask = rf2d.brightMask.copy()
        self.s = float32(rf2d.s.copy())
        
class bareBonesBlobFinder2d(blobFinder2d):
    def __init__(self, bf2d):
        self.ss = bareBonesScaleSpace2d(bf2d.ss)
        self.P = float32(bf2d.P.copy())
        self.D2 = float32(bf2d.D2.copy())
        self.darkMask = bf2d.darkMask.copy()
        self.brightMask = bf2d.brightMask.copy()
        self.s = float32(bf2d.s.copy())
        
        
    
class basicBlobNetwork(object):
    """
    aggregate the results of blob detection and select the most important blobs
    """
    def __init__(self, blobFinder, kind = 'dark', neighFun = neigh3d_124connected, fname = None):
        """
        takes a finished blobFinder object and a function that returns
        the neighborhood of a pixel and uses them to form a network 
        of blob points and to distinguish individual blobs
        """
        #assume the file name (for the image) ends with .xxx, e.g. .tif
        if fname==None:
            self.fname = 'blobNet.tif'
        else:
            self.fname = str(fname)
        
        #having a blob finder is essential to forming a blob network...
        self.bf = bareBonesBlobFinder2d(blobFinder)
        
        self.graph2sly = graph2sly
        
        #make a graph connecting individual blob points
        if kind=='dark':
            self.pointGraph = mask2graph(self.bf.darkMask)(neighFun)
        elif kind =='bright':
            self.pointGraph = mask2graph(self.bf.brightMask)(neighFun)
        
        #identify connected components in the point graph as individual blobs
        self.buildComponentGraphs()
    
    def buildComponentGraphs(self):
        """
        individual blobs are identified here as connected components of 
        the blob point network
        """
        #now find groups of blob points
        self.componentGraphs = array(networkx.connected_component_subgraphs(self.pointGraph), dtype = object)
        self.blobInit()
        
        #compute mean scale of each group and discard large blobs
        #self.trim()
        self.slys = [self.graph2sly(g) for g in self.componentGraphs]
        
        self.refreshBlobFeatures()
        
        
    def blobInit(self):
        #create blob objects
        self.B = array([blob( self.bf.ss.im, g, self.bf.ss.scales) for g in self.componentGraphs])
        [b.computeFootprint( ) for b in self.B]
        self.nBlobs = self.B.shape[0]

    def refreshBlobFeatures(self):
        self.featureKeys = ['meanStdDev', 'meanBrightness', 
                            'minBrightness','medianBrightness',
                            'nBlobPoints','nFootprintPoints',
                            'meanCurvature', 'areaFraction']
                            
        self.nFeatures = len(self.featureKeys)
        
        self.blobFeatures = zeros((self.nBlobs, len(self.featureKeys)))
        bf = self.blobFeatures
        
        s = array([mean(self.bf.s[sly]) for sly in self.slys])
        idx = argsort(s)
        
        self.B = self.B[idx]
        self.componentGraphs = self.componentGraphs[idx]
        self.slys = [self.graph2sly(g) for g in self.componentGraphs]
        
        bf[:,0] = array([b.meanSig for b in self.B])
        bf[:,1] = array([sum(b.density*b.im[b.footprintSly]) for b in self.B])
        bf[:,2] = array([median(b.im[b.footprintSly]) for b in self.B])
        bf[:,3] = array([amin(b.im[b.footprintSly]) for b in self.B])
        bf[:,4] = array([len(g.nodes()) for g in self.componentGraphs])
        bf[:,5] = array([b.footprintSly[0].shape[0] for b in self.B])
        bf[:,6] = array([mean(self.bf.s[sly]) for sly in self.slys])
        bf[:,7] = array([float64(b.footprintSly[0].shape[0])/b.im.size for b in self.B])
        
        
        valList = [bf[:,i] for i in range(self.nFeatures)]
        
        self.featureDict = dict(zip(self.featureKeys, valList)); del valList
    def pickleName(self):
        return pickleName(self.fname,'blobNet')
        
    def pickle(self):
        with open(self.pickleName, 'wb') as f:
            f.dump(self, f, -1)
    
    

class interactiveBlobNetwork(basicBlobNetwork):
    def __call__(self, sliderKeys = None, sliderSense= None):
        """
        set up the figures, attach callbacks and begin the interaction
        
        sliderKeys - list of strings into featureKeys specifying which features 
                    to enable as sliders
        sliderSense - corresponding list of +/- 1 telling whether blob features
                        must be greater than (+1) or less than (-1) the 
                        threshold
        """
        self.sliderKeys = sliderKeys
        self.sliderSense = sliderSense
        self.nActiveFeat = len(self.sliderKeys)
        
        self.activeFeatures = array([self.featureDict[k] for k in self.sliderKeys])
        
        self.isVisible = zeros(self.nBlobs)==0
        self.mask = self.isVisible.copy()
        
        self.idxIm = arange(self.bf.ss.im.size).reshape(self.bf.ss.im.shape)
        self.footMask = zeros_like(self.bf.ss.im)
        
        self.fig = figure()
        subplot(121)
        imshow(self.bf.ss.im)
        self.ax = subplot(122)
        imshow(self.bf.ss.im, cmap = cm.gray)
        subplots_adjust(top = .95, bottom=0.25)
        
        self.ellipsePlot(idx = arange(self.mask.shape[0]))
        self.thresholdInit()
        self.sliderInit()
       # print 'merging visible blobs now...'
#        self.merge(None)
        self.nMergedBlobs = 0
        self.mergedFlag = False
        self.update(self.sliders[0].val)
        
        show()
    
        
        
    def merge(self, event):
        print 'merging...'
        print 'the number of visible blobs is %d'%sum(self.mask!=0)
        self.visibleBlobNetwork = networkx.Graph()
        G = self.visibleBlobNetwork 
        
        for b in self.B[self.mask]:
            G.add_node(b)
            
        for b in self.B[self.mask]:
            for g in self.B[self.mask]:
                if G.has_edge(b,g):
                   continue
                elif b.isConnected(g):
                    G.add_edge(b,g)
                    
        self.connectedBlobs = networkx.connected_component_subgraphs(G)
        
        self.mergedBlobs = [mergeBlobs(g) for g in self.connectedBlobs]
        
        self.nMergedBlobs = len(self.mergedBlobs)
        print '...merge complete, will now update.'
        
        self.update(None)
        self.mergedFlag = True
        print 'writing summary to text file...'
        self.write()
        print 'done.'
        

        
        
    def sliderInit(self):
        """
        generate the sliders for the thresholds
        """
        axcolor = 'lightgoldenrodyellow'
        sliderAxes = [axes([0.25, 0.05+.05*i, 0.65, 0.015], axisbg=axcolor)
                        for i in range(self.nActiveFeat)]
      
        slf = self.sliderKeys
        smn = self.sliderMins
        smx = self.sliderMaxs
        svl = (smx+smn)/2
        
        self.sliders = [Slider(sliderAxes[i], slf[i], smn[i], smx[i], valinit=svl[i])
                            for i in range(self.nActiveFeat)]
        [s.on_changed(self.update) for s in self.sliders]
        
        buttonax = axes([0.5, 0.3, 0.1, 0.04])
        button = Button(buttonax, 'Merge', color=axcolor, hovercolor='0.975')
        button.on_clicked(self.merge)
        
        
    def thresholdInit(self):
        """
        set thresholds as loosely as possible to begin with
        """
        nF = self.nActiveFeat
        self.sliderMins = array([amin(self.featureDict[k]) for k in self.sliderKeys])
        self.sliderMaxs = array([amax(self.featureDict[k]) for k in self.sliderKeys])
        self.thresholds = zeros(nF)
        for i in range(nF):
            assert self.sliderSense[i]!=0, 'slider sense must be a positive or negative number!'
            if self.sliderSense[i]>0:#this slider is a lower bound
                self.thresholds[i] = self.sliderMins[i]
            elif self.sliderSense[i]<0:#this slider is an upper
                self.thresholds[i] = self.sliderMaxs[i]
        


    def update(self, val):
        """
        callback for the sliders' on_change() event
        gets the slider values, puts them into the thresholds
        makes the threshold mask, then updates the visibility of 
        each footprintPlot 
        """
        print 'updating...'
        self.thresholds = array([s.val for s in self.sliders])
        self.getVisibility()
        self.thresholdMask()
        print '...mask complete, will now set visibility'
        self.setVisibility()
        self.updateTitle()
        
    def updateTitle(self):
        
        self.footMask*=0
        for b in self.B[self.mask]:
            self.footMask[b.footprintSly]+=1

        areaFraction = float64(flatnonzero(self.footMask).shape[0])/self.bf.ss.im.size
                
        self.fig.sca(self.ax)
        title('Showing %d of %d connected components\n comprising %d components after last merge \n having area fraction %.2f'
                %(sum(self.mask!=0), self.B.shape[0], self.nMergedBlobs, areaFraction ))
        draw()
        
    
        
    def getVisibility(self):
        """
        loop over footprintPlots and get their visibility state
        """
        
        self.isVisible = array([L.get_visible() for L in self.ax.collections])
            
        
    def setVisibility(self):
        """
        set the visibility of the footprintPlots according to the current
        thresholdMask
        """
        iv = self.isVisible
        L = self.ax.collections
        [L[i].set_visible(not iv[i]) for i in self.changed]
        
    def thresholdMask(self):
        """
        compare the blob features to the thresholds and make a mask
        that is true for those blobs that meet all threshold criteria
        """
        self.mask = ones_like(self.isVisible)
        print 'thresholds are '+str(self.thresholds)
        for i in range(self.nActiveFeat):
            assert self.sliderSense[i]!=0, 'slider sense must be a positive or negative number!'
            if self.sliderSense[i]>0:#this slider is a lower bound
                self.mask = logical_and(self.mask, 
                                    self.activeFeatures[i]>self.thresholds[i])
            elif self.sliderSense[i]<0:#this slider is an upper bound
                self.mask = logical_and(self.mask, 
                                    self.activeFeatures[i]<self.thresholds[i])
                
        
        self.mask = self.mask!=0
        self.changed = flatnonzero(self.isVisible!=self.mask)
        print '%d have changed!'%self.changed.shape[0]
    
    def summaryName(self):
        return summaryName(self.fname, 'blobNet')
    def write(self):
        if not self.mergedFlag:
            self.merge(None)
            
        with open(self.summaryName(), 'wb') as f:
            f.write(self.mergedBlobs[0].writeHeader()+'\n')
            for b in self.mergedBlobs:
                f.write(str(['%.2e'%s 
                    for s in b.writeArray()])[1:-1].replace(',','').replace('\'','')+'\n')
                
        
        
    def ellipsePlot(self, ax = None, idx = None):
        if idx == None:
            idx = flatnonzero(self.mask)        
        if ax == None:
            ax = self.ax
        
        cnt = 0
        for gg in self.componentGraphs[idx]:
            x,y = (array(gg.nodes())[:,1], array(gg.nodes())[:,2])
            t = array(gg.nodes())[:,0]
            
            
            c = (cnt-1)/float64(self.componentGraphs[idx].shape[0])
            ec = EllipseCollection(
                        2*self.bf.ss.scales[t]**.5,
                        2*self.bf.ss.scales[t]**.5,
                        0*x,
                        units='x',
                        offsets=r_['0,2', y,x].T,
                        transOffset=ax.transData,
                        facecolor = plt.cm.RdYlGn(c), 
                        edgecolor = plt.cm.RdYlGn(c), 
                        alpha = .35 )

            ax.add_collection(ec)
            cnt+=1


    def footprintPlot(self,ax = None, idx = None):
        if idx == None:
            idx = flatnonzero(self.mask)
        if ax == None:
            figure()
            ax = subplot(1,2,1)
            imshow(self.B[0].im)
            
        self.ellipsePlot(ax=ax, idx=idx)
        
        subplot(1,2,2)
        imshow(self.B[0].im)
        
        cnt = 0
        for b in self.B[idx]:
            c = (cnt-1)/float64(self.B[idx].shape[0])
            plot(b.X[1], b.X[0], '.',ms = 2, color = plt.cm.RdYlGn(c), alpha = .55)
            cnt+=1
        
        xlim(0, self.B[0].im.shape[1])
        ylim(0, self.B[0].im.shape[0])
        show()
        
            
            
        
        
            
    
class basicRidgeNetwork(basicBlobNetwork):
    def __init__(self, ridgeFinder, kind = 'dark', neighFun = neigh3d_124connected, fname = None):
        """
        takes a finished blobFinder object and a function that returns
        the neighborhood of a pixel and uses them to form a network 
        of blob points and to distinguish individual blobs
        """
        #assume the file name (for the image) ends with .xxx, e.g. .tif
        if fname==None:
            self.fname = 'blobNet.tif'
        else:
            self.fname = str(fname)
        
        #having a blob finder is essential to forming a blob network...
        self.bf = bareBonesRidgeFinder2d(ridgeFinder)
        
        self.graph2sly = graph2sly
        
        #make a graph connecting individual blob points
        if kind=='dark':
            self.pointGraph = mask2graph(self.bf.darkMask)(neighFun)
        elif kind =='bright':
            self.pointGraph = mask2graph(self.bf.brightMask)(neighFun)
        
        #identify connected components in the point graph as individual blobs
        self.buildComponentGraphs()
        
class crudeBlobNetwork(basicBlobNetwork):
    def __init__(self, blobFinder, 
                     neighFun = neigh3d_124connected, 
                     sizeCutoff= 5.5, 
                     nKeep = 500, 
                     pctKeep = 1., 
                     saliencyPctKeep = 1.):
        
        self.darknessWeight = 1e-3
        
        #only keep blob groups with a mean scale less than cutoff
        self.cutoff = sizeCutoff
        
        #only display those groups with all ranks higher than these
        self.nKeep = nKeep
        self.pctKeep = pctKeep
        self.saliencyPctKeep = saliencyPctKeep
    
        basicBlobNetwork.__init__(self, blobFinder, neighFun)
        
        #compute the saliency of each blob
        self.blobSaliency()
        
        #sort by saliency
        self.sort()
        
        #compute additional stats for each group
        self.refreshGraphStats()
        
        #compute the rank of each group in a number, fraction, and saliency fraction sense
        self.rank()
    
    def buildComponentGraphs(self):
        #now find groups of blob points
        self.componentGraphs = array(networkx.connected_component_subgraphs(self.pointGraph), dtype = object)
        self.blobInit()
        
        #compute mean scale of each group and discard large blobs
        self.trim()
        self.slys = [self.graph2sly(g) for g in self.componentGraphs]
        
        self.refreshGraphStats()
    
    def refreshBlobFeatures(self):
        self.meanSig = array([b.meanSig for b in self.B])
        self.meanBrightness = array([sum(b.density*b.im[b.footprintSly]) for b in self.B])
        
    def refreshGraphStats(self):
        """
        compute any other stats for blob groups here
        """
        self.nPoints = array([len(g.nodes()) for g in self.componentGraphs])
        self.meanVal = array([mean(self.bf.ss.ss[sly]) for sly in self.slys])
        self.minVal = array([amin(self.bf.ss.ss[sly]) for sly in self.slys])
    def blobSaliency(self):
        """
        compute a measure of the importance of each blob group, the saliency
        """
        self.refreshBlobFeatures()
        self.s = array([mean(self.bf.s[sly]) for sly in self.slys])

    
    def sort(self):
        """
        sort the groups by saliency 
        """
        idx = argsort(self.s)
        self.componentGraphs = self.componentGraphs[idx]
        self.s = self.s[idx]
        
        self.slys = [self.graph2sly(g) for g in self.componentGraphs]
    
    def rank(self):
        """
        compute a rank for each group in terms of number, number percentile, and 
        saliency percentile
        """
        self.nRank = flipud(arange(self.componentGraphs.shape[0]))
        self.pctRank = float64(arange(self.componentGraphs.shape[0]))/self.componentGraphs.shape[0]
        self.saliencyPctRank = cumsum(self.s)/sum(self.s)
    def trim(self):
        self.refreshBlobFeatures()
        self.componentGraphs = self.componentGraphs[self.meanSig<self.cutoff]
        self.B = self.B[self.meanSig<self.cutoff]
        self.slys = [self.graph2sly(g) for g in self.componentGraphs]
        self.refreshBlobFeatures()
        self.refreshGraphStats()
        
    def retrim(self, cutoff):
        self.cutoff = cutoff
        self.buildComponentGraphs()
        
        self.blobSaliency()
        self.sort()
        self.refreshGraphStats()
        self.rank()    
        
    def snip(self):
        """
        apply the cutoffs to the rankings and choose the most stringent requirement
        """
        N = self.s.shape[0]
        snips = [N-self.nKeep, 
                 argmin(abs((1-self.pctKeep)-self.pctRank)),
                 argmin(abs((1-self.saliencyPctKeep)-self.saliencyPctRank))]
        print 'snip parameters are ' +str(array([self.nKeep, self.pctKeep, self.saliencyPctKeep]))
        print 'snipping at the most stringent point'+str(N-array(snips))
        self.snidx = arange(max(snips), N)
        
    def resnip(self, keepTuple = (None,None, None)):
        if keepTuple[0]!=None:
            self.nKeep = keepTuple[0]
        if keepTuple[1]!=None:
            self.pctKeep = keepTuple[1]
        if keepTuple[2]!=None:
            self.saliencyPctKeep = keepTuple[2]
        
        self.snip()

    def plot(self):
        self.snip()
        idx = self.snidx
        
        im = self.bf.ss.ss[0]
        
        
        
        #############################################
        #######  blob overlay figure ################
        #############################################
        ellipseFig = figure()  
        subplot(1,2,1)
        imshow(im, cmap = cm.gray)     
        M = zeros_like(self.bf.brightMask[0,...])
        cnt = 1.
        
        subplot(1,2,2)
        ax = ellipseFig.gca()
        imshow(im, cmap = cm.gray)

        for gg in self.componentGraphs[idx]:
            x,y = (array(gg.nodes())[:,1], array(gg.nodes())[:,2])
            t = array(gg.nodes())[:,0]
            M[x,y] = cnt

            
            if mean(self.bf.ss.scales[t]**.5<self.cutoff):# or blobTestOn and not realTestOn:
                c = (cnt-1)/float64(idx.shape[0])
                ec = EllipseCollection(
                            2*self.bf.ss.scales[t]**.5,
                            2*self.bf.ss.scales[t]**.5,
                            0*x,
                            units='x',
                            offsets=r_['0,2', y,x].T,
                            transOffset=ax.transData,
                            facecolor = plt.cm.RdYlGn(c), 
                            edgecolor = plt.cm.RdYlGn(c), 
                            alpha = .18
                            )
    
                ax.add_collection(ec)
             
            cnt+=1.
        
        title('image with most salient \nconnected components overlaid\ncircle size - mean scale of connected components')
        ########################################################
        ########################################################
      
        #///////////////////////////////////////////////////////
        #////////////  histograms and cdfs /////////////////////
        histFig = figure()
        
        
        L = self.nPoints
        N = self.nPoints.shape[0]
        
        subplot(2,2,1)
        hist(log10(float64(L)), 100)
        xlabel('nr nodes')
        title('log10(connected component sizes)')
        
        subplot(2,2,2)
        plot(self.nRank,1-float64(cumsum((L)))/sum(L))
        ylabel('fraction total blob points')
        xlabel('component rank')
        ylim(0, 1.1)
        xlim(0, N)
    
        subplot(2,2,3)
        hist(log10(float64(self.s)), 100)
        xlabel('log10(total node strength)')
        title('connected component saliency')
        
        subplot(2,2,4)
        plot(self.nRank,1-self.saliencyPctRank)
        xlabel('component rank')
        ylabel('saliency percentile')
        ylim(0, 1.1)
        xlim(0, N)
        #///////////////////////////////////////////////////////
        #///////////////////////////////////////////////////////
        
        #\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
        #\\\\\\\\\\\\\\\\\\  scatter plots \\\\\\\\\\\\\\\\\\\\\
        scatFig = figure()
        subplot(1,3,1)
        plot(self.nPoints, self.s, 'k.')
        xlabel('nr nodes')
        ylabel('saliency')
        
        subplot(1,3,2)
        plot(self.meanSig, self.s, 'k.')
        xlabel('meanScale')
        ylabel('saliency')
        
        subplot(1,3,3)
        plot(self.nPoints, self.meanSig, 'k.')
        xlabel('meanScale')
        ylabel('nr nodes')
        
        #\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
        #\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
        return ( ellipseFig, histFig, scatFig)
        
        
    def overlayPlot(self, idx):
        im = self.bf.ss.ss[0]
        #############################################
        #######  blob overlay figure ################
        #############################################
        ellipseFig = figure()  
        subplot(1,2,1)
        imshow(im, cmap = cm.gray)     
        M = zeros_like(self.bf.brightMask[0,...])
        cnt = 1.
        
        subplot(1,2,2)
        ax = ellipseFig.gca()
        imshow(im, cmap = cm.gray)

        for gg in self.componentGraphs[idx]:
            x,y = (array(gg.nodes())[:,1], array(gg.nodes())[:,2])
            t = array(gg.nodes())[:,0]
            M[x,y] = cnt

            
            if mean(self.bf.ss.scales[t]**.5<self.cutoff):# or blobTestOn and not realTestOn:
                c = (cnt-1)/float64(idx.shape[0])
                ec = EllipseCollection(
                            2*self.bf.ss.scales[t]**.5,
                            2*self.bf.ss.scales[t]**.5,
                            0*x,
                            units='x',
                            offsets=r_['0,2', y,x].T,
                            transOffset=ax.transData,
                            facecolor = plt.cm.RdYlGn(c), 
                            edgecolor = plt.cm.RdYlGn(c), 
                            alpha = .18 )
    
                ax.add_collection(ec)
             
            cnt+=1.
        
        title('image with most salient \nconnected components overlaid\ncircle size - mean scale of connected components')
        ########################################################
        ########################################################
    def footprintPlot(self, idx):
        im = self.bf.ss.ss[0]
        
        #############################################
        #######  blob overlay figure ################
        #############################################
        ellipseFig = figure()  
        subplot(1,2,1)
        imshow(im, cmap = cm.gray)     
        M = zeros_like(self.bf.brightMask[0,...])
        cnt = 1.
        
        subplot(1,2,2)
        ax = ellipseFig.gca()
        imshow(im, cmap = cm.gray)

        for b in self.B[idx]:
            x,y = (int32(b.X[0]), int32(b.X[1]))
            M[x,y] = cnt

            plot(y,x, '.', color = plt.cm.RdYlGn(float64(cnt)/idx.shape[0]))
            cnt+=1
            
        title('image with most salient \nconnected components overlaid\ncircle size - mean scale of connected components')
        ########################################################
        ########################################################
        

    
        
if __name__ == '__main__':
    diagnosticOn = False
    realTestOn = False
    blobTestOn = True
    plotDerivativesOn = False
    plotCrossSectionOn = False
    plotMasksOn = False
    plotSummaryOn = True
    plotDirectionsOn = False
    computeOn = True
    buildGraphOn = True
    
    ##::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    ##::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    ##::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    if diagnosticOn:
        x = linspace(-1,1, 49)
        y = linspace(-1,1, 51)
        im = ones_like(y)[...,newaxis]*gauss(x, 0,.1)
        #im = ones_like(y)*gauss(x, 0,.1)[...,newaxis]
        im+=normRand(0,.9).rvs(im.size).reshape(im.shape)
        print "matplotlib backend is "+mpl.get_backend()
    ##::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    ##::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    ##::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    ##::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    if realTestOn:
        computeOn = True
        fname = 'cracksRF.pkl'
        if computeOn:
            sz = 512
            im_ = float64(imread('cracks.tif')[...,0])
            im = match(zeros((im_.shape[0]*sz/im_.shape[1],sz )), im_)
        else:
            fname = 'cracksRF.pkl'
            print "compute flag not thrown, loading pickle..."
            start = time.time()
            with open(fname, 'rb') as f:
                rf = pkl.load(f)
                ss = rf.ss
                im = rf.ss.ss[0]
            print "done. %d"%(time.time()-start)
            
    if blobTestOn and not realTestOn:
        print "blob test is on and real test is not!"
        computeOn = False
        fname = 'blobsRF.pkl'
        if computeOn:
            im = float64(imread('blobs.tif'))
            
        else:
            print "compute flag not thrown, loading pickle..."
            start = time.time()
            fname = 'blobsRF.pkl'
            with open(fname, 'rb') as f:
                rf = pkl.load(f)
                ss = rf.ss
                im = rf.ss.ss[0]
            print "done. %d"%(time.time()-start)
        
    ##::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    ##::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    ##::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    if computeOn:
        print "initializing scale space..."
        start = time.time()
        if blobTestOn and not realTestOn:
            ss = scaleSpace2d(im, nscales = 100, maxFeatureSize = .65)
        else:
            ss = scaleSpace2d(im, nscales = 100, maxFeatureSize = .35)
        print "done. %d"%(time.time()-start)
        
#        print "computing derivatives..."
#        start = time.time()
#        ss.updateDerivatives()
#        print "done. %d"%(time.time()-start)
        
        print "initializing ridge finder..."
        start = time.time()
        if blobTestOn and not realTestOn:
            print "belay that, it's a blob finder!"
            rf = blobFinder2d(ss)
        else:
            rf = ridgeFinder2d(ss)
        print "done. %d"%(time.time()-start)
        
        print "finding ridges..."
        start = time.time()
        r = rf()
        print "done. %d"%(time.time()-start)
        
        print 'bouncing to pickle...'
        start = time.time()
        with open(fname, 'wb') as f:
            pkl.dump(rf, f, -1)
        print "done. %d"%(time.time()-start)
    
    ##::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    ##::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    ##::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    if buildGraphOn:
        #//////////// build the graph /////////////////
        print "building graph..."
        start = time.time()
        #g = mask2graph(any(rf.darkMask, axis = 0))()
        g = mask2graph(rf.darkMask)()
        print "done. %d"%(time.time()-start)
        
        
        #//////////// find connected components/////////////////
        print "finding connected subgraphs..."
        start = time.time()
        G = array(networkx.connected_component_subgraphs(g), dtype = object)
        print "done. %d"%(time.time()-start)
        
        
        #//////////// number of nodes and ridge saliency/////////////////
        L = array([len(gg.nodes()) for gg in G])
        V = array([sum(rf.s[tuple([array(gg.nodes())[:,i] for i in range(3)])])**.5 for gg in G])
        #V = array([sum((256.-rf.ss.ss[tuple([array(gg.nodes())[:,i] for i in range(3)])])/25.+ rf.s[tuple([array(gg.nodes())[:,i] for i in range(3)])])**.5 for gg in G])
        meanVal = array([mean(rf.ss.ss[tuple([array(gg.nodes())[:,i] for i in range(3)])]) for gg in G])
        meanScale = array([mean(rf.ss.scales[array(gg.nodes())[:,0]]**.5) for gg in G])
        idx = flatnonzero(L>1)
        
        V = V[meanVal<median(meanVal)]
        
        
        snip = V.shape[0]-1000
        sidx = argsort(V)[snip:]
        idx = sidx
        
        
        M = zeros_like(rf.brightMask[0,...])
        cnt = 1.
        
        oFig = figure()
#        ax = subplot(1,2,1)
        ax = oFig.gca()
        imshow(im, cmap = cm.gray)
#        subplot(1,2,2)
#        imshow(any(rf.darkMask, axis = 0), cmap = cm.Accent)
#        for gg in G:
#            #x,y = (array(gg.nodes())[:,0], array(gg.nodes())[:,1])
#            x,y = (array(gg.nodes())[:,1], array(gg.nodes())[:,2])
#            plot(y,x, 'k.', ms = 4)
        
        
        
        for gg in G[idx]:
            #x,y = (array(gg.nodes())[:,0], array(gg.nodes())[:,1])
            x,y = (array(gg.nodes())[:,1], array(gg.nodes())[:,2])
            t = array(gg.nodes())[:,0]
            M[x,y] = cnt
            
            
            if blobTestOn and not realTestOn:
                cutoff = 5.5
            else:
                cutoff = 3.5
            if mean(rf.ss.scales[t]**.5<cutoff):# or blobTestOn and not realTestOn:
                c = (cnt-1)/float64(idx.shape[0])
                ec = EllipseCollection(
                            rf.ss.scales[t]**.5,
                            rf.ss.scales[t]**.5,
                            0*x,
                            units='x',
                            offsets=r_['0,2', y,x].T,
                            transOffset=ax.transData,
                            facecolor = plt.cm.RdYlGn(c), 
                            edgecolor = plt.cm.RdYlGn(c), 
                            alpha = .25 )
    
                ax.add_collection(ec)
            
#            subplot(1,2,2)
#            plot(y,x, 'y.', ms = 4)
            
            
            cnt+=1.
        
#        
#        subplot(1,2,1)
        title('image with most salient \nconnected components overlaid\ncircle size - mean scale of connected components')
        #legend()
#        subplot(1,2,2)
#        title('all ridge pixels with retained connected components overlaid\n y - large, k - all')

        
        
        

        figure()
        subplot(1,2,1)
        hist(log10(float64(L)), 100)
        xlabel('nr nodes')
        title('log10(connected component sizes)')
        
        subplot(1,2,2)
        plot(1- float64(cumsum(sort(L)))/sum(L))
        ylabel('1 - fraction total edge length')
        xlabel('component rank')
        ylim(0, 1.1)
        xlim(0, L.shape[0])
    
        figure()
        subplot(1,2,1)
        hist(log10(float64(V)), 100)
        xlabel('log10(total node strength)')
        title('connected component saliency')
        
        subplot(1,2,2)
        plot(1- cumsum(sort(V))/sum(V))
        xlabel('component rank')
        ylabel('1 - fraction total saliency')
        ylim(0, 1.1)
        xlim(0, V.shape[0])
    ##::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    
    

    ##::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    ##::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    ##::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    if plotDerivativesOn:
        figure()
        subplot(2,2,1)
        imshow(im, cmap = cm.copper)
        title('the image')
        titles = [r'im$_{yy}$',r'im$_{xy}$',r'im$_{xx}$']
        for i in range(3):
            subplot(2,2,i+2)
            imshow(ss.D[1][i,10], cmap = cm.copper); colorbar()
            title(titles[i])
        
        figure()
        titles = [r'im$_{y}$',r'im$_{x}$']
        for i in range(2):
            subplot(1,2,i+1)
            imshow(ss.D[0][i,0], cmap = cm.copper); colorbar()
            title(titles[i])
        
        sly = (slice(None, None), x.shape[0]/2, slice(None, None))
        #sly = (slice(None, None), slice(None, None),y.shape[0]/2)
    ##::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    
    
        
    ##::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    ##::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    ##::::::::::::::::::::::::::::::::::::::::::::::::::::::::::    
    if plotCrossSectionOn:
        figure()
        subplot(2,2,1)
        imshow( ss.ss[sly], cmap = cm.gray)
        title('scale space representation')
        
        subplot(2,2,2)
        imshow( rf.s[sly], cmap = cm.gray)
        title('ridge strength')
        
        subplot(2,2,3)
        imshow( rf.sMask[sly], cmap = cm.gray)
        title('strength mask')
        
        subplot(2,2,4)
        imshow( rf.d[sly], cmap = cm.gray)
        title('direction')
    ##::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    
    
    
    ##::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    ##::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    ##::::::::::::::::::::::::::::::::::::::::::::::::::::::::::    
    if plotMasksOn:
        figure()
        subplot(2,2,1)
        imshow( sum(rf.sMask, axis = 0), cmap = cm.gray)
        title('strength mask')
        
        subplot(2,2,2)
        imshow( any(rf.brightMask,axis = 0), cmap = cm.gray)
        title('bright mask')
        
        subplot(2,2,3)
        imshow( any(rf.darkMask,axis = 0), cmap = cm.gray)
        title('dark mask')
        
        subplot(2,2,4)
        quiver(rf.P[1,0], rf.P[0,0], color = 'k')
        quiver(rf.Q[1,0], rf.Q[0,0], color = 'y')
        title('P in black, Q in yellow')
    ##::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
        
        
    ##::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    ##::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    ##::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    if plotDirectionsOn:
        figure()
        for i in range(4):
            subplot(2,2,i+1)
            imshow(rf.d[i*rf.d.shape[0]/4], cmap = cm.gray); colorbar()
            title('i = %d'%(i*rf.d.shape[0]/4))
    
    
    
    
    ##::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    ##::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    ##::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    if plotSummaryOn:
        figure()
        subplot(1,2,1)
        imshow(im, cmap = cm.gray)
        title('the image')
        subplot(1,2,2)
        #imshow(sum(rf.ss.scales[:,newaxis, newaxis]*rf.darkMask, axis = 0 ))
        #imshow(log(1+sum(rf.darkMask, axis = 0 )), cmap = cm.Accent)
        imshow(any(rf.darkMask, axis = 0 ), cmap = cm.gray)
        title('dark ridge points')
    ##::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    
        
    
    
    
    show()
