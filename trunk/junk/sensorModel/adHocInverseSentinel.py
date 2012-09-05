import numpy as np
from numpy import *
import scipy
scipy.pkgload('weave')
from scipy import weave
import cPickle as pkl
import time
import pdb
from pylab import *
import particle as part
import sensor as sens
import particleModelEnsemble as pmEns
from plotMacros import *
import mlabMacros as mlm
from interpolators import zerothOrderInterpolator

def alignArrays(master, slave, tol):
    """
    match each element of master to the nearest element of slave within a distance tol of 
    the master element. elements of master that are not matched are discarded
    
    return the matched 
    """
    
    
    mIdx = []
    sIdx = []
    J = 0
    
    for i in range(master.shape[0]):
        theMin = -1
        for j in range(J,slave.shape[0]):
            #if we are left of the ROI, continue, if we are right of the ROI, break
            if slave[j]<(master[i]-tol):
                continue
            elif slave[j]>(master[i]+tol):
                J = j
                break
            #if we are in the ROI, direct search for the minimum value
            else:
                #first time, initialize theMin
                d = abs(master[i]-slave[j])
                dSlavePlus = abs(master[i]-slave[min([j+1, slave.shape[0]-1])])
                dMasterPlus = abs(master[min([i+1, master.shape[0]-1])]-slave[j])
                
                if d<=dSlavePlus and d<=dMasterPlus:
                        J = min([j+1, slave.shape[0]-1])
                        mIdx+=[i]
                        sIdx+=[j]
#                        print mIdx
#                        print sIdx
#                        print master
#                        print slave
                        break
    
    return (mIdx,sIdx)       

class adHocInverseSentinel(object):
    """
    use clever guesses to obtain from a sensor image the parameters of the particles that made it
    """
    def __init__(self, fluxSentinel = sens.basicFluxSentinel(), LUTGenerate = False, highPassOn = True, lowPassOn = True):
        
        self.FS = fluxSentinel
        self.sensorImage = zeros_like(self.FS.sensorOutput['sensorImage'])
        
        self.peakThreshold =  3
        self.peakAlignmentTol = 10
        self.peakDropTol = .75
        self.highPassConstant = .99
        self.setIntensityConstant()
        
        if LUTGenerate:
            self.generateLUTs()
        else:
            self.generateBlackbodyCurve()
            self.setLUT()
        
        self.highPassOn = highPassOn
        self.lowPassOn = lowPassOn
        self.setFilterKernels()
        
    def generateLUTs(self, **kwargs):
        pass
    
    def computeTemperature(self):
        pass
    def lookupTemperature(self, color):
        pass
    def computeColor(self, arr):
        pass
        
    def setIntensityConstant(self):
        """
        set the constant factor in I = C(D^2T^4/V)
        """
        boltzmannConstant = 5.67e-8
        emissivity = 1
        self.intensityConstant = (self.FS.gains[0]*
                (self.FS.pixelSize/self.FS.lens.magnification)*
                (.5*self.FS.lens.diameter/self.FS.lens.objectPlane)*
                .5*pi*emissivity*boltzmannConstant)
                
    def setLUT(self, LUTfiles = ['tLUT.pkl', 'rLUT.pkl']):
        """
        load default temperature and residual lookup tables if LUTfiles is None
        if filenames are passed in a list, load the LUT's from these
        """
        ft = open(LUTfiles[0], 'rb')
        fr = open(LUTfiles[1], 'rb')
        self.tLUT = pkl.load(ft)
        self.rLUT = pkl.load(fr)
        
    def setFilterKernels(self):
        """
        build any FIR kernels that will be used to prefilter the image
        """
        self.lowPassKernel = ones(3.)
        proto = ones_like(self.lowPassKernel)
        pad = lambda master, slave: r_[zeros( (slave.shape[0]-1)/2), master,zeros((slave.shape[0]-1)/2)]
        for i in range(5):
            self.lowPassKernel = scipy.ndimage.convolve1d(pad(self.lowPassKernel,proto), proto)
        self.lowPassKernel/=sum(self.lowPassKernel)
    
    def generateBlackbodyCurve(self,p = part.emittingParticle(), tMinMax = [1500, 5000], n = 350):
        """
        use the passed emittingParticle to generate a curve of its color vs temperature
        """
        
        T = linspace(1500, 5000, n)
        rgb = zeros((T.shape[0],3))
        
        for i in range(T.shape[0]):
            rgb[i] = squeeze(sum(p.intensityProfile(1, 1e-9*self.FS.wavelength, T = T[i])[newaxis,:]*self.FS.spectralSensitivity[:,newaxis,...], axis = 2))[1:]
            rgb[i]/= sqrt(sum(rgb[i]**2))
        
        self.blackbodyCurve = r_['1,2,0',T, rgb].T
        plotOn = False
        if plotOn:
            c = ['r','g','b']
            for i in range(1,4):
                plot(self.blackbodyCurve[0], self.blackbodyCurve[i], linewidth = 2, color = c[i-1])
       
    def invert(self, pEns = None, img = None):
        """
        given a particle ensemble pEns, use the fluxSentinel to generate an image, then invert that image
        or given and image img, set this as the sensorImage and invert it
        """
        assert not (pEns==None and img ==None), "both args to invert() can't be none, give an image or a particle ensemble!"
        
        if pEns !=None:
            self.FS.sense(pEns)
            self.sensorImage = self.FS.sensorOutput['sensorImage'].copy()
        else:
            self.sensorImage = img
            
         ##################### FILTERING ######################### 
        if self.lowPassOn:
            self.sensorImage = self.lowPass(self.sensorImage)
       
        if self.highPassOn:
            for i in range(self.sensorImage.shape[0]):
                self.sensorImage[i] = self.highPass(self.sensorImage[i])
        #######################################################
        
        ############# PEAK SEARCH AND PROCESSING ##################
        self.peakFind()
        
        #allocate the array for the answer
        self.parameterVectors = zeros((self.peakIdx[0].shape[0],9))
        ############################################################
        
        ################ DIAMETER ESTIMATION######################
        D = []
        for sly in self.slices.tolist():
            D+=[self.inferDiameter(self.sensorImage[0][sly])]
        self.setDiameters(D)
        ###########################################################
        
        #################### POSITION ESTIMATION #####################
        self.setPositions((self.peakIdx[0] - self.FS.sensorShape[1]/2)*self.FS.pixelSize/self.FS.magnification)
        ##############################################################
        
        ####################TEMPERATURE ESTIMATION####################
        self.computeTemperature()
        #################################################################
        
        ################### INTENSITY AND VELOCITY ESTIMATION ##############
        self.intensity = array([sum(self.sensorImage[0][sly]) for sly in self.slices])
        self.maxIntensity =  array([np.max(self.sensorImage[0][sly]) for sly in self.slices])
        self.computeVelocity()
        #####################################################################
       
       #post process the answer into a particle ensemble
        self.generateEnsemble()
    
    def particleWiseInvert(self, pEns):
        """
        process and invert each particle's image individually, producing a 1 to 1 matching ensemble
        """
        self.FS.sense(pEns)
        imgs = self.FS.sensorOutput['particleImages']
        
        lowPassPad = lambda a: r_['-1', zeros(a.shape[:-1]+(1+(self.lowPassKernel.shape[0]-1)/2,)), a, zeros(a.shape[:-1]+(1+(self.lowPassKernel.shape[0]-1)/2,))]
        
        self.particleImages = empty(imgs.shape, dtype = 'object')
        self.parameterVectors  = zeros(imgs.shape+(pEns.parameterVectors.shape[-1],) )
        D = zeros(imgs.shape)
        T = zeros_like(D)
        self.temperatureResiduals = zeros_like(D)
        y = zeros_like(D)
        self.intensity = zeros_like(D)
        self.maxIntensity = zeros_like(D)
        for idx, x in ndenumerate(imgs):
            im = self.particleImages[idx]
            im =  self.lowPass(lowPassPad(x))
            for i in range(im.shape[0]):
                im[i] = self.highPass(im[i])
                
            D[idx] = self.inferDiameter(im[0])
            
            T[idx], self.temperatureResiduals[idx] = self.lookupTemperature(self.computeColor(im[1:]))
            y[idx] = pEns.modelArray[idx].PS.p.verticalPosition()
            self.intensity[idx] = sum(im[0])
            self.maxIntensity[idx] = np.max(im[0])
            self.particleImages[idx] = im.copy()
        self.setDiameters(D)
        self.setTemperatures(T)
        self.setPositions(y)
        self.computeVelocity()
        self.generateEnsemble()
        
    
    def lowPass(self, arr):
        return scipy.ndimage.convolve1d(arr, self.lowPassKernel, mode = 'constant')
    
    def highPass(self, arr):
        n = arr.shape[0]
        y = zeros(n)
        c = self.highPassConstant
        code = """
                    y[0] = 0;
                    for(int i =1; i<n;i++){
                        y[i] = c*( arr[i]- arr[i-1] + y[i-1]) ;
                        if (y[i]<0)
                            y[i] = 0;
                        };
                    return_val = 0;
                """
        retval = weave.inline(code, ['arr', 'c', 'y', 'n'], compiler = 'mingw32')    
        return y
                    
    def peakFind(self):
        """
        find peaks and valleys in the image above a certain threshold, put their locations into peakIdx, valleyIdx.
        check for alignment between channels and throw out any peaks that aren't in all channels
        """
        nChan = self.sensorImage.shape[0]
        self.peakIdx = empty(nChan, dtype = 'object')
        self.valleyIdx = empty(nChan, dtype = 'object')
        
        d = diff(int32(sign(diff(self.sensorImage, axis = 1))))
        for i in range(0,nChan):
            self.peakIdx[i] = sort(1+flatnonzero(logical_or(d[i]==-2, d[i]==-1)))
            self.valleyIdx[i] = sort(1+flatnonzero(logical_or(d[i]==2, d[i]==1)))
            
            self.peakIdx[i] = self.peakIdx[i][flatnonzero(diff(r_[self.peakIdx[i], 
                                atleast_1d(self.sensorImage.shape[1]+2)])>1)]
            self.peakIdx[i] = self.peakIdx[i][flatnonzero(self.sensorImage[i][
                                                self.peakIdx[i]]>self.peakThreshold)]
        
        self.parsePeaks()
        self.alignPeaks()
        
    def parsePeaks(self):
        """
        parse the peaks according to their neighboring valleys
        """
      #aliases 
        p = r_[self.peakIdx[0], atleast_1d(self.sensorImage.shape[1]+2)]
        v = self.valleyIdx[0]
        s = self.sensorImage[0]
        
        P = [0]
        Vleft = []
        Vright = []
        J = 0
        
        for i in range(p.shape[0]-1):
            for j in range(J,v.shape[0]-1):
                #if the current valley is: less than the current peak, greater than the last peak
                #and the next valley is greater than the current peak, less than the next peak
                cond = v[j]<p[i] and v[j+1]>p[i] and v[j]>P[-1] and v[j+1]<p[i+1]
                #print "condition is  "+ str(cond)+" at (%d,%d)"%(i,j)
                if cond:
                    if s[v[j]]<=.25*s[p[i]] and s[v[j+1]]<=.25*s[p[i]]:
                        P+=[p[i]]
                        Vleft+= [v[j]]
                        Vright+=[v[j+1]]
                        J = j+1
                        break
        self.peakIdx[0] = array(P[1:])
        
        self.slices = empty(len(Vleft), dtype = 'object')
        for i in range(len(Vleft)):
            self.slices[i] = int32(arange(Vleft[i], Vright[i]+1))
    
    def alignPeaks(self):
        """
        verify that the peak appears on all channels
        """
        #check bw and b
        bwIdx, bIdx = alignArrays(self.peakIdx[0], self.peakIdx[3], self.peakAlignmentTol)
        self.peakIdx[0] = self.peakIdx[0][bwIdx]
        self.slices = self.slices[bwIdx]
        #check bw and g
        bwIdx, gIdx = alignArrays(self.peakIdx[0], self.peakIdx[2], self.peakAlignmentTol)
        self.peakIdx[0] = self.peakIdx[0][bwIdx]
        self.slices = self.slices[bwIdx]
        #check bw and r
        bwIdx, rIdx = alignArrays(self.peakIdx[0], self.peakIdx[1], self.peakAlignmentTol)
        self.peakIdx[0] = self.peakIdx[0][bwIdx]
        self.slices = self.slices[bwIdx]
        #re-check b and g now that bw has been cut down to final size
        bwIdx, bIdx = alignArrays(self.peakIdx[0], self.peakIdx[3], self.peakAlignmentTol)
        bwIdx, gIdx = alignArrays(self.peakIdx[0], self.peakIdx[2], self.peakAlignmentTol)
        
        self.peakIdx[1]= self.peakIdx[1][rIdx]
        self.peakIdx[2]= self.peakIdx[2][gIdx]
        self.peakIdx[3]= self.peakIdx[3][bIdx]
        
        
        
        self.slices = self.slices[bwIdx]
        
  
    def inferDiameter(self, arr):
        #peak index and midpoint
        midx = argmax(arr)
        p = arr[midx]
        mid = midx
        
        #left and right slices
        sL = float64(arange(0,midx+1))
        sR = float64(flipud(arange(midx,arr.shape[0])))
        
        #left and right images
        left = arr[int32(sL)]
        right =flipud(arr[flipud(int32(sR))])
        
        #alias for readability
        interp = scipy.interpolate.interp1d
        M = self.FS.lens.magnification
        pix = self.FS.pixelSize
        
        #left and right diameter estimates
        dL = 2*pix*(mid - interp(left, sL, kind = 'linear')(.5*p))/M
        dR = 2*pix*(interp(right, sR, kind = 'linear')(.5*p)-mid)/M
        
        return min([dL,dR])
        
    
    
    def computeVelocity(self):
        """
        from temperature, diameter, and intensity estimates, estimate the velocity
        """
        self.setVelocities(self.intensityConstant*self.diameter()**2
                            *self.temperature()**4/self.intensity)
        
        
    def generateEnsemble(self):
        """
        create a particle ensemble from the parameter array that has been built up
        """
        self.outputEnsemble = pmEns.emittingParticleEnsemble(pmEns.dummyEmittingParticleGenerator)
        self.outputEnsemble.generate(pVec = self.parameterVectors)
        
    def setDiameters(self, D):
        """
        use the list or array of diameters D to set the proper column of parameterVectors
        """
        self.parameterVectors[...,0] = .01*array(D)
        self.parameterVectors[...,1] = array(D)
        
    def setTemperatures(self,T):
        """
        use the list or array of temperatures T to set the proper column of parameterVectors
        """
        self.parameterVectors[...,2] = array(T)
        
    def setPositions(self, y,z = None):
        """
        use the list or array of positions y and z to set the proper columns of parameterVectors
        """
        self.parameterVectors[...,3] = self.FS.axialPos
        self.parameterVectors[...,4] = array(y)
        if z!=None:
            self.parameterVectors[...,5] = array(z)
            
    def setVelocities(self, u,v = None,w = None):
        """
        use the list or array of velocities u to set the proper columns of parameterVectors
        """
        self.parameterVectors[...,6] = array(u)
        
        if v != None:
            self.parameterVectors[...,7] = array(v)
        if w != None:
            self.parameterVectors[...,8] = array(w)
            
    def diameter(self):
        return self.parameterVectors[...,1]
    def temperature(self):
        return self.parameterVectors[...,2]
    def pixelTemperature(self):
        rgb = uint8(255*self.sensorImage[1:]/np.max(self.sensorImage[1:]))
        return self.tLUT[tuple(rgb.tolist())]
        
    def verticalPosition(self):
        return self.parameterVectors[...,4]
    def axialVelocity(self):
        return self.parameterVectors[...,6]
    def show(self):
        figure(figsize = (20,20))
        c = ['k','r', 'g', 'b']
        #for i in range(self.sensorImage.shape[0]):
        for i in range(1):
            plot(self.sensorImage[i], color = c[i], linewidth = 1.25)
        for i in range(self.slices.shape[0]):
            plot(self.slices[i],i*ones_like(self.slices[i]), 's-',color = 'y', linewidth = 1.5)
        ylim(0,256)
        axisFontSize()
    
class adHocInverseSentinel_3dLUT(adHocInverseSentinel):
    def generateLUTs(self, fnames = ['tLUT','rLUT']):
        """
        use the spectral radiance of the passed emitting particle class and the 
        spectral sensitivity of the attached fluxSentinel to generate 
        lookup tables for temperature and residual
        """
        print 'generating LUT, please be patient'
                
        try:
            print "blackbody curve has %d entries"%(self.blackbodyCurve[0].shape[0])
        except:
            "generating blackbody curve with default values"
            self.generateBlackbodyCurve()
           
        print 'blackbody curve generation complete, real work begins now...'
        self.tLUT = zeros((256,256,256))
        self.rLUT = zeros_like(self.tLUT)
        T = self.blackbodyCurve[0]
        rgb = self.blackbodyCurve[1:].T
        n = T.shape[0]
        
        
        try:
            tLUT = self.tLUT
            rLUT = self.rLUT
            code = """
                    unsigned long idx = 0;
                    double mres = 0;
                    double r = 0;
                    const int N = 256;
                    double norm = 0;
                    std::cout<< "we are in c now\\n";
                    for(int i = 0; i<N; i++){
                        if (i%20==0)
                            std::cout << "completed " << i*N*N << " of " << N*N*N << " colors \\n";
                        for(int j = 0; j<N; j++){
                            for(int k = 0; k<N; k++){
                                idx = N*(N*i+j)+k;
                                mres = 0;
                                norm = sqrt(pow(i,2)+pow(j,2)+pow(k,2));
                                for(int m = 0; m<n; m++){
                                    r = i*rgb[m*3]/norm+j*rgb[m*3+1]/norm+k*rgb[m*3+2]/norm;
                                    if(r>mres){
                                        mres = r;
                                        tLUT[idx] = T[m];
                                        rLUT[idx] = r;
                                    }
                                }
                            }
                        }
                    }
                    return_val = 0;
                    """
            print 'beginning inline generation of table at '
            print time.ctime()
            start = time.time()
            retval = weave.inline(code, ['n', 'tLUT', 'rLUT','T','rgb'], compiler = 'mingw32', headers = ['<iostream>'], verbose = 0, force = 0)
        
        except:
            print 'something went wrong with inline, doing it the slow way'
            d = zeros_like(T)
            idx = [0,0,0]
            r = zeros(3)
            
            start = time.time()
            for i in arange(256.):
                r[0] = i
                idx[0] = int32(i)
                elapsed = time.time()-start
                print 'completed %d of %d colors: %2.f percent'%(i*256**2,256**3,float64(i)/256.)
                print 'elapsed time is %d seconds or %d minutes'%(elapsed, elapsed/60)
                for j in arange(256.):
                    r[1] = j
                    idx[1] = int32(j)
                    # if mod(j,20)==0:
                        # print str(j)
                    for k in arange(256.):
                        r[2] = k
                        r/=sqrt(sum(r**2))
                        idx[2] = int32(k)
                        d = sqrt(sum((r*rgb)**2, axis = 1))
                        midx = argmin(d)
                        self.tLUT[tuple(idx)] = T[midx]
                        self.rLUT[tuple(idx)] = d[midx]

        elapsed = time.time()-start
        print 'generation ended at '
        print time.ctime()
        print 'elapsed time is %d s'%(elapsed)
        
        print 'generation complete! saving...'
        ft = open(fnames[0]+'.pkl','wb')
        fr = open(fnames[1]+'.pkl','wb')
        pkl.dump(self.tLUT,ft,-1)
        pkl.dump(self.rLUT,fr,-1)
        ft.close()
        fr.close()
        
    def lookupTemperature(self, color):
        return (self.tLUT[tuple(uint8(color))], self.rLUT[tuple(uint8(color))])
        
    def computeColor(self, arr):
        return squeeze(mean(arr,axis = 1))
        
    def computeTemperature(self):
        """
        loop over the peaks in peakIdx and use the color information to look up a temperature and a residual for each
        """
        T = zeros(self.parameterVectors.shape[0])
        self.temperatureResiduals = zeros_like(T)
        self.rgb = zeros((T.shape[0],3))
        rgb = float64(zeros(3))
        normConst = 255./np.max(self.sensorImage[1:])
        for j in range(len(self.slices)):
            sly = self.slices[j]
            p = self.peakIdx[0][j]
            img = self.sensorImage
            rgb = normConst*computeColor(r_['0,2', img[1][sly-(p-self.peakIdx[1][j])], img[2][sly-(p-self.peakIdx[2][j])]  , img[3][sly-(p-self.peakIdx[3][j])]  ])
           
            if any(rgb<1):
                rgb/=sum(rgb)/100
            T[j], self.temperatureResiduals[j] = self.lookupTemperature(color)
            self.rgb[j] = rgb
        self.setTemperatures(T)
        
class adHocInverseSentinel_2dLUT(adHocInverseSentinel):
    """
    use a different LUT method for this sentinel based on the ratios of green to red, blue to red
    """
    def setLUT(self, LUTfiles = ['tLUTratio.pkl','rLUTratio.pkl', 'bbCurve.pkl']):
        adHocInverseSentinel.setLUT(self, LUTfiles[:-1])
        f = open(LUTfiles[-1],'rb')
        self.blackbodyCurve = pkl.load(f)
        f.close()
        
        self.gb = self.blackbodyCurve[2:]/self.blackbodyCurve[1]
        self.tLookup = zerothOrderInterpolator( self.gb[:,0],self.gb[:,-1], self.tLUT, precision = 'single')
        self.rLookup = zerothOrderInterpolator(self.gb[:,0],self.gb[:,-1], self.rLUT, precision = 'single')
        
    def generateLUTs(self, fnames = ['tLUTratio','rLUTratio', 'bbCurve'], N = 2000):
        """
        use the spectral radiance of the passed emitting particle class and the 
        spectral sensitivity of the attached fluxSentinel to generate 
        lookup tables for temperature and residual
        """
        print 'generating rational LUT, please be patient'
        
        try:
            print "blackbody curve has %d entries"%(self.blackbodyCurve[0].shape[0])
        except:
            print "generating blackbody curve with default values"
            self.generateBlackbodyCurve(n = 2500)
           
        print 'blackbody curve generation complete, real work begins now...'
        self.tLUT = zeros((N,N))
        self.rLUT = zeros_like(self.tLUT)
        T = self.blackbodyCurve[0]
        gb = self.blackbodyCurve[2:]/self.blackbodyCurve[1]
        self.gb = gb
        g = linspace(gb[0,0], gb[0,-1], N)
        b = linspace(gb[1,0], gb[1,-1], N)
        G = gb[0]
        B = gb[1]
        n = T.shape[0]
        
        
        tLUT = self.tLUT
        rLUT = self.rLUT
        code = """
                unsigned long idx = 0;
                double minD = 0;
                double dist = 0;
                std::cout<< "we are in c now\\n";
                for(int i = 0; i<N; i++){
                    if (i%20==0)
                        std::cout << "completed " << i*N << " of " << N*N << " colors \\n";
                    for(int j = 0; j<N; j++){
                        minD = 10;
                        idx = (N*i+j); 
                        for(int m = 0; m<n; m++){
                            dist = pow(g[i]-G[m],2)+pow(b[j]-B[m],2);
                            if(  dist < minD  ){
                                minD = dist;
                                tLUT[idx] = T[m];
                                rLUT[idx] = dist;
                            }
                        }
                    }
                }
                return_val = 0;
                """
        print 'beginning inline generation of table at '
        print time.ctime()
        start = time.time()
        retval = weave.inline(code, ['n', 'N','tLUT', 'rLUT','T','g','b','G','B'], compiler = 'mingw32', headers = ['<iostream>'], verbose = 1, force = 0)
        self.rLUT = sqrt(self.rLUT)

        elapsed = time.time()-start
        print 'generation ended at '
        print time.ctime()
        print 'elapsed time is %d s'%(elapsed)
        
        
        self.tLookup = zerothOrderInterpolator( gb[:,0],gb[:,-1], self.tLUT, precision = 'single')
        self.rLookup = zerothOrderInterpolator(gb[:,0],gb[:,-1], self.rLUT, precision = 'single')
        print 'generation complete! saving...'
        ft = open(fnames[0]+'.pkl','wb')
        fr = open(fnames[1]+'.pkl','wb')
        fb = open(fnames[2]+'.pkl','wb')
        pkl.dump(self.tLUT,ft,-1)
        pkl.dump(self.rLUT,fr,-1)
        pkl.dump(self.blackbodyCurve,fb,-1)
        ft.close()
        fr.close()
        fb.close()
   
    def lookupTemperature(self, color):
        return (self.tLookup(color[...,newaxis]), self.rLookup(color[...,newaxis]))
    
    def computeColor(self, arr):
        rgb = mean(arr, axis = 1)
        rgb/= rgb[0]
            
        if rgb[1]>self.gb[0,-1]:
            rgb[1]=self.gb[0,-1]
        elif rgb[1]<self.gb[0,0]:
            rgb[1]=self.gb[0,0]
        
        if rgb[2]>self.gb[1,-1]:
            rgb[2]=self.gb[1,-1]
        elif rgb[2]<self.gb[1,0]:
            rgb[2]=self.gb[1,0]
        
        return rgb[1:]
        
    def computeTemperature(self):
        """
        do the lookup with the new kind of LUT
        """
        T = zeros(self.parameterVectors.shape[0])
        self.temperatureResiduals = zeros_like(T)
        self.rgb = zeros((T.shape[0],2))
        rgb = float64(zeros(3))
        
        for j in range(len(self.slices)):
            sly = self.slices[j]
            p = self.peakIdx[0][j]
            img = self.sensorImage
            gb = self.computeColor(r_['0,2',   img[1][sly-(p-self.peakIdx[1][j])], 
                                                                 img[2][sly-(p-self.peakIdx[2][j])], 
                                                                 img[3][sly-(p-self.peakIdx[3][j])]  ])
            
            
            T[j], self.temperatureResiduals[j] = self.lookupTemperature(gb)
            self.rgb[j] = gb
            
        self.setTemperatures(T)
        
    
if __name__== "__main__":
    alignTestOn = False
    lutGenerateTestOn = False
    invertTestOn = True
    newLutGenerateTestOn = False
    timeTrialOn = False
    
    if alignTestOn:
        a = sort(unique(np.random.randint(0,200,10)))
        b = sort(unique(a+np.random.randint(-5,5,a.shape[0]) ))
        
        
        plot(a,'.', markerfacecolor = None, label= 'master')
        plot(b, '.', markerfacecolor = None,label = 'slave')
        
        aIdx, bIdx = alignArrays(a,b,3)
        print 'aIdx:'+ str(aIdx)
        print 'bIdx:'+ str(bIdx)
        print 'a:'+ str(a)
        print 'b:'+ str(b)
    
    if lutGenerateTestOn:
        IS = adHocInverseSentinel_2dLUT(sens.basicFluxSentinel(), LUTGenerate = True)
        idx = flatnonzero(logical_and(IS.tLUT.ravel()!=0, IS.tLUT.ravel()!=5000))
        close('all')
        figure()
        hist(log(IS.rLUT.ravel()[idx]), 100)
        figure()
        hist(IS.tLUT.ravel()[idx], 100)
        rgb = int32(IS.blackbodyCurve[1:]*200)
        T = IS.blackbodyCurve[0]
        Tlook = IS.tLUT[rgb[0], rgb[1],rgb[2]]
        figure()
        plot(T, linewidth = 2)
        plot(Tlook, 'x')
    
    if newLutGenerateTestOn:
        close('all')
        print "running new lookup table test"
        print "inverse sensor initializing"
        IS = adHocInverseSentinel_2dLUT(LUTGenerate = True)
        r = IS.blackbodyCurve[1]
        g = IS.blackbodyCurve[2]/r
        b = IS.blackbodyCurve[3]/r
        T = IS.blackbodyCurve[0]
        
#        print "plotting in 3d"
#        mlm.surf3(x = g, y =b, z = IS.tLUT, f = mlm.fig('temperature lookup'), axeLab = ['g', 'b', 'T'])
#        mlm.surf3(x = g, y =b, z = IS.rLUT, f = mlm.fig('distance lookup'), axeLab = ['g', 'b', 'd'])
#        mlm.plot3(r_['0,2', g, b, T], f = mlm.fig('blackbody curve'), axeLab = ['g', 'b', 'T'])
        figure()
        hist(IS.tLUT.ravel(), 100)
        xlabel('Temperature, K')
        figure()
        hist(IS.rLUT.ravel(), 100)
        xlabel('Distance')
        
        print "looking up blackbody temperatures"
        Tlook = IS.tLookup(r_['0,2', g,b])
        print "max error: %.2f"%np.max(abs(Tlook-T))
        print "mean error: %.2f"%np.mean(abs(Tlook-T))
        figure()
        plot(T, linewidth = 2)
        plot(Tlook, 'x')
        show()
    
    if timeTrialOn:
        import cProfile as prof
        import pstats
        
        PE = pmEns.emittingParticleEnsemble(pmEns.marginalRandomEmittingParticleGenerator)
        PE.generate(N = 1e3)
        IS = adHocInverseSentinel_2dLUT()
        prof.run('IS.invert(PE)', 'prof')
        S = pstats.Stats('prof').strip_dirs()
        S.sort_stats('time').print_stats()
        S.sort_stats('cumulative').print_stats()
        
        
    if invertTestOn:
        print "running model inversion test"
        print "genrating inverse sensor model: IS"
        IS = adHocInverseSentinel_2dLUT()
        print "generating input ensemble: PE"
    #    PE = pmEns.emittingParticleEnsemble(pmEns.basicGriddedEmittingParticleGenerator)
        PE = pmEns.emittingParticleEnsemble(pmEns.marginalRandomEmittingParticleGenerator)
#        PE.generate(dMinMax = [20e-6, 60e-6],
#                    tMinMax = [2700, 4000], 
#                    zMinMax = [0,1e-3], 
#                    n = (3,3,1,1), nPixels = IS.FS.sensorShape[1])
        PE.generate(N = 1e2)
        print "beginning inversion"
        start = time.time()
        IS.invert(PE)
        elapsed = time.time()-start
        print "inversion completed in %f s"%elapsed
        print "mean, std T in: %d , %d"%(mean(PE.parameterVectors[...,2]), std(PE.parameterVectors[...,2]))
        print "mean, std T out: %d , %d"%(mean(IS.temperature()), std(IS.temperature()))
        print "mean, std  D in: %d , %d"%(mean(PE.parameterVectors[...,1]*1e6), std(PE.parameterVectors[...,1]*1e6))
        print "mean, std D out: %d , %d"%(mean(IS.diameter()*1e6),std(IS.diameter()*1e6))
        
        IS.show()
        figure()
        plot(IS.FS.sensorOutput['sensorImage'].T)
        ylim(0,256)
        
        figure(figsize = (20,20))
        subplot(1,2,1)
        hist(PE.parameterVectors[...,2], 100, label = 'input')
        hist(IS.temperature(), 100, label = 'output')
        xlabel('temperature, K')
        legend()
        subplot(1,2,2)
        hist(PE.parameterVectors[...,1]*1e6, 100, label = 'input')
        hist(IS.diameter()*1e6, 100, label = 'output')
        xlabel('diameter, K')
        legend()
        
        show()
        
    
    
    