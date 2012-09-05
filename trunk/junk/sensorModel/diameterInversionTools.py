from numpy import *
import numpy as np
from plotMacros import *
import mlabMacros as mlm
from scipy import integrate
from scipy import interpolate
from interpolators import cubicInterpolator as pol
from scipy.ndimage import median_filter as smoothMed
from scipy.ndimage import gaussian_filter as smoothGauss
from scipy.stats import norm as normRand 
import pdb
import time
import os
import cPickle as pkl
smooth = lambda im: smoothMed(im, size = 8)
smoothg = lambda im: smoothGauss(im, sigma = 3.)


def gauss(x,mu,sig):
    return exp(-(x-mu)**2/(2*sig**2))/sqrt(2*pi*sig**2)

def polate(x,y):
    return interpolate.interp1d(x,y, fill_value = 0., bounds_error = False)
    
def generateObservations(pd, mu, sig, n):
    """
    use the interpolant dictionary  of a 
    forwardDiameterMapper object to 
    pick n measurement-observation pairs for inversion
    by sampling from gauss-distributed measurements
    mu[0] - diameter mean, pixels
    sig[0] - diameter standard dev, pixels
    mu[1] - z mean, mm
    sig[1] = z std dev, mm
    """
    m = float64([normRand(mu[i], sig[i]).rvs(n) for i in range(2)])
    obs = [pd['Dmeas'][i](m) for i in range(3)]
    obs += [pd['B'][i](m) for i in range(3)]
    return (m, float64(obs))


#this class was hastily written to test the idea
# and is now depreciated since we changed our whole scheme to accomodate color
class blurModel(object):
    """
    """
    def __init__(self, z,W,sigz,zeroMean = True, cutoff = Inf):
        #out of plane distance in mm
        self.z = z 
        
        #standard deviation of out of plane distance for a particle in mm
        self.sigz = sigz
        
        #width of the blur kernel in pixels
        self.W = W
        
        #(change in kernel width)/(change in z distance)
        self.m = mean(abs(diff(self.W)/diff(self.z))) 
        
        #standard dev of the kernel width
        self.sigw = self.m*self.sigz
        
        #can enforce zero blur at zero out of plane distance
        if zeroMean:
            self.muw = 0.
        else:
            self.muw = min(self.W)
        
        #mesh over kernel width
        self.w = linspace(self.muw, 4*self.sigw, 500)
        
        #pdf of kernel width
        self.p = 2*gauss(self.w, self.muw, self.sigw)
        self.pFun = polate(self.w, self.p)
        
        #cdf of kernel width
        self.f = r_[0, integrate.cumtrapz(self.p, self.w)]
        self.fFun = polate(self.w, self.f)
        
        #statistics of the kernel width
        self.meanKernelWidth = integrate.trapz(self.w*self.p, self.w)
        self.medianKernelWidth = self.w[argmin(abs(self.f-.5))]
        self.q1KernelWidth = self.w[argmin(abs(self.f-.25))]#1st quartile
        self.q3KernelWidth = self.w[argmin(abs(self.f-.75))]#2nd quartile
         
    
    def pdf(self, w):
        return self.pFun(w)
        
    def cdf(self, w):
        return self.fFun(w)
    
    def quartiles(self):
        return r_[self.q1KernelWidth, self.medianKernelWidth,self.q3KernelWidth ]
        
    def whisker(self):
        f = fig()
        q = self.quartiles()
        plot(atleast_1d([self.meanKernelWidth]),atleast_1d([1]),'+',ms=10, figure = f)
        plot(q, ones_like(q), '.-', ms = 8, figure = f)
        return f
#this class is not used but still hanging around
class observationFileReader(object):
    def __init__(self, fname):
        self.fname = fname
        self.read()
        self.cnt = 0
    def setFile(self,fname):
        self.__init__(fname)
        
    def next(self):
        if self.cnt>=self.X.shape[0]:
            return None
        else:
            return self.X[self.cnt]
            
    def read(self):
        with open(self.fname, 'r') as f:
            X = []
            for L in f.readlines():
                X+=[float64(L.replace('\t', ' ').split())]
        self.X = array(X)
        
#base class for implementation of the forward mapping, model->observation
class forwardMappingFileReader(object):
    """
    read a file that whose lines are (modelParameter1 \t modelparameter2 ...\t observation1 \t observation2...)
    and set up interpolators for the functions observations(models)
    this class is meant to be subclassed for specific applications
    """
    
    def __init__(self, fname, cutoff = -1):
        """
        set file name, read in the data, call the format specific headerCallback()
        """
        self.fname = fname
        self.read()
        self.headerCallback()
    def headerCallback(self):
        """
        make dictionaries for ease of use that reference case-specific parameter names
        """
        pass
        
    def setFile(self, fname, cutoff = -1):
        self.__init__(fname, cutoff)
        
    def read(self):
        #pull in the tabulated functions Dmeas(D,w), B(D,w)
        with open(self.fname) as f:
            lines = f.readlines()
            X = []
            
            for L in lines:
                try:
                    X+=[float64(L.replace('\t',' ').split())]
                    
                except ValueError:
                    print 'passed'
                    pass
            
            self.X = array(X).T
        #D (pixels), z position (mm),  Dmeas r-g-b (microns), blur numbers r-g-b
        
        #not used any more but if you want to chop the table you can hard code it here
        cutoff = -1
        
        #first and second columns are presumed to be the values of the model parameters
        self.model = [sort(unique(self.X[0])),sort(unique(self.X[1]))]
        
        #since we will tranpose everything this will be the final shape of the array
        n = (self.model[1].shape[0], self.model[0].shape[0])
        
        #define a macro to reshape the linear array from the file into a 2d array
        f = lambda x: smoothg((reshape(x.copy(),n).T))[:cutoff]
        
        #concatenate 2d arrays of observations into one big 3d array
        self.meas = float64([f(self.X[i]) for i in range(2,self.X.shape[0])])
        
        #apply the cutoff to the model vector as well
        self.model[0] = self.model[0][:cutoff]
        
        #this will create the dictionary of interpolators
        self.refreshPolDict()
    
    def refreshPolDict(self):
        """
        make an interpolator for every observation
        """
        self.mn = r_[self.model[0][0], self.model[1][0]]
        self.mx = r_[self.model[0][-1], self.model[1][-1]]
        
        #interpolators 
        self.pol = [pol(self.mn, self.mx, self.meas[i]) 
                    for i in range(self.meas.shape[0])]

#this class includes a method specific to our data format                    
class forwardDiameterMapper(forwardMappingFileReader):
    def headerCallback(self):
        """
        define dictionaries for the interpolators and the arrays underlying them
        not strictly necessary but will make the code that follows easier to track
        """
        self.headers = ['D','z', 'Dmeas', 'B']
        self.arrDict = dict(zip(self.headers, [self.model[0], self.model[1],
                            self.meas[:3],self.meas[3:6] ]))
        
        self.polDict = dict(zip(self.headers[2:], [[self.pol[i] for i in range(3)],
                            [self.pol[i] for i in range(3,6)]]))
                            
#this class performs the probabilistic inversion according to tarantola
#it solves for the map observations -> model by computing the probability of 
#every model for a given observation
class inverseDiameterMapper(object):
    def __init__(self,forwardMapper, sigs = r_[2.,2., 2.,.02,.02, .02 ]):
        """
        forwardMapper - object that computes the forward map
        sigs - standard devs representing uncertainty in measurement of D and B
        """

        self.sigs = sigs
        self.fm = forwardMapper
        self.computeFlag = True
        
        
    def __call__(self, D,z):
        """
        does nothing more than prepare the observation vectors corresponding to D,z
        in the form that modelLikelihood() expects
        """
        #we can get the observations for any model vectors using interpolation, but if we 
        #want to use the default mesh of the model space, we don't need to interpolate anything and we 
        #can make it much faster
        if self.computeFlag:
            g = meshgrid(z,D)
            
            #these are vectors in the "model space" of tarantola
            #cache for later use
            self.modelVec = r_['0,3', g[1].copy(), g[0].copy()]
            
            return r_['0,3', float64([self.fm.polDict['Dmeas'][i](self.modelVec) for i in range(3)]),
                            float64([self.fm.polDict['B'][i](self.modelVec) for i in range(3)])]
        else:
            return self.fm.meas
    
 
        
        
    
    def modelLikelihood(self, obsVec, D = None, z = None):
        """
        model likelihood factor for all true diameters D (1d) 
        and plume locations z (1d),
        for given observed values of diameters and blur numbers
        """
        #by default just use the grids we have stored in the forward mapper
        if z==None:
            z = self.fm.arrDict['z'].copy()
        if D ==None:
            D = self.D.copy()
        
        assert (D.ndim==1 and z.ndim==1)
        
        self.dataVec = self(D,z)
        
        #we are using gaussian measurement uncertainties with no covariance, so in the end it's just a product
        #of individual measurement uncertainties
        self.L = ones((D.shape[0], z.shape[0]))
        for i in range(self.dataVec.shape[0]):
            self.L*= gauss(self.dataVec[i], obsVec[i], self.sigs[i])
            

        return self.L.copy()
        
    def modelProbability(self, obsVec, 
                             D = None, z = None, 
                             pFuns = None):
        """
        compute the probability of every model D,z  (each 1d)
        for a given Dobs and Bobs (each 0d)
        and functions pFuns that compute the prior 
        probabilities for each model
        """
        
        #in this case we can be lazy since we have the observation values cached
        if z==None and D==None:
            self.computeFlag = False
        else:
            self.computeFlag = True
            
        if z==None:
            z = self.fm.arrDict['z'].copy()
        if D ==None:
            D = self.fm.arrDict['D'].copy()
        
        #use uniform priors if we are not passed a prior function
        if pFuns==None:
           pFuns = []
           for i in range(obsVec.shape[0]):
               pFuns+= lambda x: ones_like(x)
        else:
            for i in range(len(pFuns)):
                if pFuns[i]==None:
                    pFuns[i] = lambda x: ones_like(x)
                    
        assert (D.ndim==1 and z.ndim==1)
        
        #preallocate for the probability density
        self.p = ones((D.shape[0], z.shape[0]))
        
        #apply the prior probabilities
        self.p*= pFuns[0](D)[...,newaxis]
        self.p*= pFuns[1](z)
        
        #compute and apply the model likelihood
        self.modelLikelihood( obsVec, D , z )
        self.p*= self.L
        
        #normalize the pdf
        self.p/=integrate.trapz(integrate.trapz(self.p,z, axis = 1), D )
        
        #catch nans and infs
        if any(isnan(self.p)) or any(isinf(self.p)):
            self.p = zeros_like(self.p)
        return self.p.copy()

#used for proving out the concept on the black and white version, this class is depreciated
class forwardMap(object):
    def __init__(self,fname, sigs = r_[2.,  .02 ]):
        
        #standard devs representing uncertainty in measurement
        #of D and B
        self.sigs = sigs
        
        #pull in the tabulated functions Dmeas(D,w), B(D,w)
        with open(fname) as f:
            lines = f.readlines()[1:]
            X = []
            
            for L in lines:
                try:
                    X+=[float64(L.split('\t'))]
                except ValueError:
                    pass
            
            self.X = array(X).T
        #D (pixels), kernel width (6sigma, pixels), blur number, Dmeas (microns)
        #self.X = self.X.T
        
        cutoff = -1
        if fname == 'D_B_lorenz_75_GiantGrid.txt':
            cutoff = 122
        
        self.D = sort(unique(self.X[0]))
        self.w = sort(unique(self.X[1]))
        
        self.B = smoothg(smooth(reshape(self.X[2].copy(),
                    (  self.w.shape[0], self.D.shape[0])).T))[:cutoff]
        self.Dmeas = smoothg(smooth(reshape(self.X[3].copy(),
                    (  self.w.shape[0], self.D.shape[0])).T/3.3))[:cutoff]
        #self.Imeas = reshape(self.X[5].copy(),(  self.w.shape[0], self.D.shape[0])).T
        #self.Imeas[:,0]=1.00001*self.Imeas[:,1]
        self.Imeas = (smooth(reshape(self.X[5],
                    (  self.w.shape[0], self.D.shape[0])).T))[:cutoff]
        #self.I = reshape(self.X[4].copy(),(  self.w.shape[0], self.D.shape[0])).T
        self.I = smoothg(smooth(reshape(self.X[4],
                    (  self.w.shape[0], self.D.shape[0])).T))[:cutoff]
       
        
        self.Brough = ((reshape(self.X[2],
                    (  self.w.shape[0], self.D.shape[0])).T))[:cutoff]
        self.Dmeasrough = (reshape(self.X[3],
                    (  self.w.shape[0], self.D.shape[0])).T/3.3)[:cutoff]
        
        
        self.D = self.D[:cutoff]
        
        self.mn = r_[self.D[0], self.w[0]]
        self.mx = r_[self.D[-1], self.w[-1]]
        
        #interpolators 
        self.DmeasFun = pol(self.mn, self.mx, self.Dmeas)
        self.BFun = pol(self.mn, self.mx, self.B)
        self.ImeasFun = pol(self.mn, self.mx, self.Imeas)
        self.IFun = pol(self.mn, self.mx, self.I)
        
        
        
    def __call__(self, D,w):
        """
        use the functional mapping from particle modeling
        """
        g = meshgrid(w,D)
        
        #these are vectors in the "model space" of tarantola
        #cache for later use
        self.modelVec = r_['0,3', g[1].copy(), g[0].copy()]
        
        return r_['0,3', self.DmeasFun(self.modelVec),  self.BFun(self.modelVec)]
    def intensityCorrector(D,w):
        g = meshgrid(w,D)
        
        #these are vectors in the "model space" of tarantola
        #cache for later use
        m = r_['0,3', g[1].copy(), g[0].copy()]

        J = self.I(m)/self.Imeas(m)
        J[isinf(J)] = 0.
        
        return J
        
        
    
    def modelLikelihood(self, obsVec, D = None, w = None):
        """
        model likelihood factor for all true diameters D (1d) 
        and kernel widths w (1d),
        for given observed values of Dobs (0d) and Bobs (0d)
        """
        if w==None:
            w = self.w.copy()
        if D ==None:
            D = self.D.copy()
        
        assert (D.ndim==1 and w.ndim==1)
        
        
        self.dataVec = self(D,w)
        
        #cache it for later
        self.L = ones((D.shape[0], w.shape[0]))
        for i in range(self.dataVec.shape[0]):
            self.L*= gauss(self.dataVec[i], obsVec[i], self.sigs[i])
        
        
        return self.L.copy()
        
    def modelProbability(self, obsVec, 
                             D = None, w = None, 
                             pFuns = None):
        """
        compute the probability of every model D,w  (each 1d)
        for a given Dobs and Bobs (each 0d)
        and functions pwFun(w), pdFun(D) that compute the prior 
        probabilities for each model
        """
        if w==None:
            w = self.w.copy()
        if D ==None:
            D = self.D.copy()
        if pFuns==None:
           pFuns = []
           for i in range(obsVec.shape[0]):
               pFuns+= lambda x: ones_like(x)
        else:
            for i in range(len(pFuns)):
                if pFuns[i]==None:
                    pFuns[i] = lambda x: ones_like(x)
                    
        assert (D.ndim==1 and w.ndim==1)
        
        self.p = ones((D.shape[0], w.shape[0]))
        
        self.p*= pFuns[0](D)[...,newaxis]
        self.p*= pFuns[1](w)
        
        self.modelLikelihood( obsVec, D , w )
        
        self.p*= self.L
        self.p/=integrate.trapz(integrate.trapz(self.p,w, axis = 1), D )
        if any(isnan(self.p) or isinf(self.p)):
            pdb.set_trace()
        return self.p.copy()

#this class streamlines the management of the classes that do that actual work, creates synthetic target distributions
#for testing, and includes the basic post-processing. it's used in scripts and such, and is the highest level object
#of those in this file

class diameterInversionCase(object):
    """
    organize all the parameters needed for the inference of diameter based on image measurements and a 
    particular forward map
    """
    def __init__(self,inputTableName = 'D_B_lorenz_color_snap_extended.txt' , #the table containing the forward map
                 pklPrefix = 'ppkl', #will save the results of the inversion as a pickle under this name
                 measurementStdDevs = [.05, .05], #error in the image measurements, in terms of log-variables by default
                 noiseMultiplier = 1e-12, #for the case of a synthetic target distribution, the amount of noise to add to the samples
                 targetMeans = [12., 0.], #mean in pixels, mm for true diameter and z location of the target distribution
                 targetStdDevs = [.2,4. ], #std dev in log(D/Dref), mm for log diameter and z location of the target distribution
                 posteriorStdThreshold = .25, #we will discard in post processing all samples with posterior variance that is too large
                 N = 500, #number of samples to run
                 Dref = 12.,#reference diameter for log transform
                 Bref = .15,#reference blur number for log transform
                 DPrior = None, #pass a prior probability for diameter
                 logarithmic = True, #whether or not to use the log transform. keeping it true is probably best.
                 forceCompute = True #if this flag is false and there exists a matching pickle on disk, just load the answer without computing it
                 ):
         
        #set all the parameters
        self.inputTableName =inputTableName
        self.pklPrefix = pklPrefix
        self.measurementStdDevs = measurementStdDevs
        self.noiseMultiplier = noiseMultiplier
        self.targetMeans = targetMeans
        self.targetStdDevs = targetStdDevs
        self.posteriorStdThreshold = posteriorStdThreshold
        self.N = N
        self.Dref = Dref
        self.Bref = Bref
        self.DPrior = DPrior
        self.logarithmic = logarithmic
        self.forceCompute = forceCompute
        self.refs = r_[ones(3)*self.Dref, ones(3)*self.Bref]
        
        if self.logarithmic:
            self.targetMeans[0] = log(self.targetMeans[0]/self.Dref)
            
    

    def __call__(self,samples = None, clean = False, nD = None):
        """
        samples - can pass our own observation vector samples from e.g. an experiment
        clean - run with no noise in the samples if clean is True
        nD - number of meshpoints in the diameter direction, used to control the resolution in diameter 
                independent of the table that created it
        """
        
        #get prepped for the inversion
        self.buildForwardMapper(nD = nD)
        self.generateSamples(samples)
        self.buildInverseMapper()
        
        #preallocate for the marginal distributions of D and z
        self.pD = zeros((self.N,self.fm.arrDict['D'].shape[0]))
        self.pz = zeros((self.N,self.fm.arrDict['z'].shape[0]))
        
        #this will be a list of selected posterior distributions
        P = []
        
        #control how often we save a whole posterior distribution and print progress
        if self.N>10000:
            pMod = 1000
        elif self.N>1000:
            pMod = 500
        else:
            pMod = 100
        
        #switch whether using clean or noisy samples
        if clean:
            obs = self.cleanObsSamples.T
        else:
            obs = self.obsSamples.T
        
        #the main loop!
        if self.forceCompute or not os.path.exists(self.pklFileName()):
            start = time.time()
            for i in range(self.N):
                #compute the joint posterior distribution of D and z. this is the answer!
                p = self.im.modelProbability(obs[i], pFuns = [self.DPrior, self.zPrior])
                
                #compute marginals by integrating out the other variable
                self.pD[i] = trapz(p, self.fm.arrDict['z'], axis = 1)
                self.pz[i] = trapz(p, self.fm.arrDict['D'], axis = 0)
                
                #print progress every so often
                if mod(i,pMod)==0:
                    print 'completed %d in %.2e s'%(i, time.time()-start )
                    #retain some joint distributions for later viewing
                    #looking at a few of them can be very helpful for diagnosing
                    #problems and adjusting measurement variances, and understanding
                    #what's really going on
                    P+=[p.copy()]
            self.P = float64(P)
            
            #save the results
            self.pickle()
        else:
            #if this case has already been run, we can just load the results
            #useful if you just want to mess with the post processing
            self.P, self.pD, self.pz = self.unpickle()
        
        #post-process the results
        self.post()
        
    def pklFileName(self):
        """
        give the pickle file a name based on the case parameters
        """
        if self.logarithmic:
            return self.pklPrefix+'_muD_%d_sigD_%d_N_%d_nois_%.2e_sigMeas_%.2e_logOn'%(
                                    self.targetMeans[0], self.targetStdDevs[0],self.N, self.noiseMultiplier, 
                                    self.measurementStdDevs[0])
        else:
            return self.pklPrefix+'_muD_%d_sigD_%d_N_%d_nois_%.2e_sigMeas_%.2e_logOff'%(
                                    self.targetMeans[0], self.targetStdDevs[0],self.N, self.noiseMultiplier, 
                                    self.measurementStdDevs[0])
    def pickle(self):
        """
        save the posterior distributions
        """
        with open(self.pklFileName(), 'wb') as f:
            pkl.dump((self.P, self.pD, self.pz),f,-1)
            
    def unpickle(self):
        """
        load the posterior distributions
        """
        with open(self.pklFileName(), 'rb') as f:
            return pkl.load(f)
            
    def zPrior(self, z):
        """
        we use the target distribution as a prior for z,
        but this function could easily be rewritten to use a different prior
        """
        return gauss(z, self.targetMeans[1], self.targetStdDevs[1])
        
    def buildForwardMapper(self, nD = None):
        """
        instantiate a forward mapping object and perform necessary conversions for 
        the log transform case
        """
        self.fm = forwardDiameterMapper(self.inputTableName, cutoff = -1)
        if not self.logarithmic:
            return 
        else:
            def flog(x,r):
                return log(x/r)
            def rlog(x,r):
                return r*exp(x)
            
            fm = self.fm
            if nD ==None:
                nD = fm.arrDict['D'].shape[0]
            
            d0 = self.Dref
            b0 = self.Bref
            
    
            D = rlog(linspace(flog(fm.arrDict['D'][0],d0),flog(fm.arrDict['D'][-1],d0), nD),d0)
            newz, newd = meshgrid(fm.arrDict['z'],D)
            polArg = r_['0,3', newd, newz]
            
            newMeas = float64([fm.pol[i](polArg) for i in range(len(fm.pol))])
#            
#            figure()
#            subplot(1,2,1)
#            contourf(fm.arrDict['z'].copy(),fm.arrDict['D'].copy(), fm.meas[0].copy(),30)
#            subplot(1,2,2)
#            contourf(fm.arrDict['z'].copy(),D.copy(), newMeas[0].copy(),30)
#            
#            print newMeas.shape
#            print fm.meas.shape
            

            fm.model[0] = flog(D,d0)
            for i in range(newMeas.shape[0]):
                newMeas[i] = flog(newMeas[i],self.refs[i])
            fm.meas = newMeas.copy()
            fm.refreshPolDict()
            fm.headerCallback()
                

    def generateSamples(self, samples = None):
        """
        if we are not passed samples, generate some from a synthetic target distribution and add the specified amount
        of noise
        if we are passed samples, just set the necessary data members and move on
        """
        if samples == None:
            self.synthFlag = True
            self.modSamples, self.cleanObsSamples = generateObservations(self.fm.polDict, 
                                    self.targetMeans, self.targetStdDevs, self.N)
            
            devs = r_[ones(3)*self.measurementStdDevs[0], ones(3)*self.measurementStdDevs[1]]
            self.obsSamples = self.cleanObsSamples+float64([normRand(0,self.noiseMultiplier*devs[i]).rvs(self.N)
                                            for i in range(self.cleanObsSamples.shape[0])])
        else:
            assert samples.shape[0] ==6
            self.synthFlag = False
            self.cleanObsSamples = samples.copy()
            self.obsSamples = samples.copy()
            self.N = self.obsSamples.shape[1]
            
       
    def buildInverseMapper(self):
        """
        instantiate the inverse mapping object
        """
        self.im = inverseDiameterMapper(self.fm, r_[ones(3)*self.measurementStdDevs[0],
                                                    ones(3)*self.measurementStdDevs[1]])
    
    def targetDiameterDist(self):
        """
        function factory for the target distribution in diameter, for use outside this object
        """
        return lambda d: gauss(d, self.targetMeans[0], self.targetStdDevs[0])
        
        
    def post(self):
        """
        post process the data, 
        creating an index of samples with successful inversion as measured by the posterior variance
        """
        
        #preallocate for the posterior means and std devs
        self.postModelMean = zeros((2,self.N))
        self.postModelStd = zeros((2,self.N))
        
        #alias for readability
        D = self.fm.arrDict['D']
        z = self.fm.arrDict['z']
        
        #compute statistics of the posterior distributions
        for i in range(self.N):
            self.postModelMean[0,i] = trapz(D*self.pD[i],D)
            self.postModelMean[1,i] = trapz(z*self.pz[i],z)

            self.postModelStd[0,i] = trapz((D-self.postModelMean[0,i])**2*self.pD[i],D)**.5
            self.postModelStd[1,i] = trapz((z-self.postModelMean[1,i])**2*self.pz[i],z)**.5
        
        #wrap it up into a nice neat dictionary
        self.postDict = dict(zip(['muD', 'muz', 'sigD', 'sigz'],
                                 [self.postModelMean[0],self.postModelMean[1],
                                 self.postModelStd[0],self.postModelStd[1]]))
        
        #find samples with good posterior variance. posterior variance might be zero if the posterior distribution itself 
        #was identically zero everywhere. this can happen if the observation vector is particularly outrageous or the 
        #measurement standard devs are too small (e.g. smaller than the grid size!)
        self.keeperIdx = flatnonzero(logical_and(self.postDict['sigD']<self.posteriorStdThreshold,self.postDict['sigD']>0.))
        
    def pdfDisjunction(self):
        """
        estimate the target distribution as the disjunction (i.e. mean) of the individual 
        distributions
        """
        return sum(self.pD[self.keeperIdx], axis = 0)/self.pD[self.keeperIdx].shape[0]
        
    def plot(self, keepersOnly = True):
        """
        make a bunch of plots, using only the samples from self.keeperIdx
        if keepersOnly is true 
        """
        
        #use keepers or everything
        if keepersOnly:
            idx = self.keeperIdx
        else:
            idx = arange(self.N)
        
        #the compliment of the keepers is the set of particles with bad variance
        nidx = sort(array(list(set(range(self.N))-set(idx))))
        
        #histograms of the mean estimates of D an z
        hm, bm =  histogram(self.postDict['muD'][idx], 50, normed = True)
        hz, bz =  histogram(self.postDict['muz'][idx], 50, normed = True)
        
        #pdf disjunction estimates of the posterior distributions
        pdfD = sum(self.pD[idx], axis = 0)/self.pD[idx].shape[0]
        pdfZ = sum(self.pz[idx], axis = 0)/self.pz[idx].shape[0]
        
        #alias the target means and std devs
        mu = self.targetMeans
        sig = self.targetStdDevs
        
        #make arrays for the domains of the target distributions
        z = linspace(mu[1]-3*sig[1], mu[1]+3*sig[1], 100)
        d = linspace(mu[0]-3*sig[0], mu[0]+3*sig[0], 100)

        #########################################################################################
        ############   big daddy figure, showing posterior and target distributions #############
        #########################################################################################
        postFig = figure()
        
        #dont plot the target distribution unless this is a synthetic run
        if self.synthFlag:
            hist(self.obsSamples[0], 50, label = 'noisy raw estimate', normed = True, facecolor = 'r', alpha = .35)
        hist(self.cleanObsSamples[0], 50, label = 'clean raw estimate',normed = True, facecolor = 'y', alpha = .35 )
        plot(bm[:-1]+diff(bm)/2, hm, 'g', label = 'mean estimate', drawstyle = 'steps')
        if self.synthFlag:
            plot(d,gauss(d, mu[0], sig[0]), 'k', label = 'target')
        plot(self.fm.arrDict['D'], pdfD, 'y',label = 'pdf disjunction estimate' )
        legend()
        ylabel('probability density')
        
        #adjust plot labels 
        if self.logarithmic:
            dstr = r'log(D/D$_{ref})$'
            title('posterior diameter distributions, logarithmic space')
        else:
            dstr = 'diameter, pixels'
            title('posterior diameter distributions')
        xlabel(dstr)
        #########################################################################################
        
        ##################################################
        #### plot posterior distribution of z for kix ####
        postFigZ = figure()
        plot(bz[:-1]+diff(bz)/2, hz, 'g', drawstyle = 'steps', label = 'mean estimate')
        plot(z,gauss(z, mu[1], sig[1]), 'k', label = 'target')
        plot(self.fm.arrDict['z'], pdfZ, 'y', label = 'pdf disjunction estimate')
        xlabel('z, mm')
        title('reconstructed position distribution')
        ##################################################
        
        #################################################################
        ############  plots to help diagnose the inversion  #############
        #################################################################
        statsFig = figure()
        
        #std devs of the posterior diameters
        subplot(2,3,1)
        hist(self.postDict['sigD'], 50)
        xlabel('diameter posterior std dev, '+dstr)
        
        #std devs of the posterior positions
        subplot(2,3,4)
        hist(log10(self.postDict['sigz']), 50)
        xlabel('position posterior std dev, '+r'log$_{10}$(z/1mm)')
        
        #axis limits for the next plots
        dlim = [[np.min(self.obsSamples[i]),np.max(self.obsSamples[i])] for i in range(2)]
        blim = [[np.min(self.obsSamples[i]),np.max(self.obsSamples[i])] for i in range(3,5)]
        
        ##scatterplots of some of the observations##
        subplot(2,3,2)
        #green vs red diameters, all samples
        plot(self.cleanObsSamples[0], self.cleanObsSamples[1], 'b.',label = 'clean', alpha = .1)
        if self.synthFlag:
            plot(self.obsSamples[0], self.obsSamples[1], 'r+',label = 'noisy', alpha = .1)
        xlim(dlim[0])
        ylim(dlim[1])
        xlabel(r'D$_r$')
        ylabel(r'D$_g$')
        title('measured diameters, all')
        legend(loc = 'upper left')
        
        subplot(2,3,3)
        #green vs red  diameters with low posterior variance
        plot(self.cleanObsSamples[0][idx], self.cleanObsSamples[1][idx], 'g.',label = 'clean', alpha = .1)
        xlim(dlim[0])
        ylim(dlim[1])
        xlabel(r'D$_r$')
        ylabel(r'D$_g$')
        title('measured diameters,\n low posterior variance')
        
        subplot(2,3,5)
        #green vs red blur numbers, all samples
        plot(self.cleanObsSamples[3], self.cleanObsSamples[4], 'b.',label = 'clean', alpha = .1)
        if self.synthFlag:
            plot(self.obsSamples[3], self.obsSamples[4], 'r+',label = 'noisy', alpha = .1)
        xlim(blim[0])
        ylim(blim[1])
        xlabel(r'D$_r$')
        ylabel(r'D$_g$')
        title('blur numbers, all')
        legend(loc = 'upper left')
        
        subplot(2,3,6)
        #green vs red blur numbers, low posterior variance
        plot(self.cleanObsSamples[3][idx], self.cleanObsSamples[4][idx], 'g.',label = 'clean', alpha = .1)
        xlim(blim[0])
        ylim(blim[1])
        xlabel(r'D$_r$')
        ylabel(r'D$_g$')
        title('blur numbers,\n low posterior variance')
        #################################################################
        
        
        return (postFig, postFigZ, statsFig)
        
        
    def pdfPlot(self):
        """
        plot the joint posterior distributions. it will make many plots so it's  broken out from plot()
        """
        for i in range(self.P.shape[0]):
            figure()
            contourf(self.fm.arrDict['z'], self.fm.arrDict['D'], self.P[i], 30, cmap = cm.Accent); colorbar()
            
        
        
if __name__== '__main__':
    
    #get input for blur model from measurements
    start = array([3065, 3099, 3135, 3108, 3085])
    stop = array([3242, 3213, 3194, 3233, 3279])
    z = linspace(-2, 2, 5)
    
    #assume we measured 4*stdDev, convert to 6*stdDev
    W = 6.*(stop - start)/4.
    
    #plumeWidth/6 in mm
    sigz = 4.
    
    #file where measured diameter and blur number are tabulated
    #theFile = 'D_B_gauss_75.txt'
    theFile = 'D_B_lorenz_75_bigGrid.txt'
    
    
    #init blur objects and forward model
    blur = blurModel(z,W, sigz)
    fm = forwardMap(theFile, sigs = r_[2.,  .02 ])
    
    
#    #test observations, gauss kernel
#    Dobs = [40.,30., 60., 82., 70., 20]
#    Bobs = [.14, .18, .16, .168, .14, .16 ]
    
    #test observations, lorentz kernel
    Dobs = [50., 59., 48.47, 60.4, 80., 18]
    Bobs = [.2, .286, .288, .294, .1, .2 ]
    
    

    labs = ['good',              'dense',  'questionable', 
            'indistinguishable', 'impossible', 'slow-pitch softball' ]
    ms = ['bo', 'r+', 'rx', 'k+', 'kx', 'go']
    P = []
    for i in range(len(Dobs)):
        P +=[ fm.modelProbability(r_[Dobs[i], Bobs[i]], pFuns = [None, blur.pdf])]
    
    fig()
    plot(fm.Dmeas.ravel(), fm.B.ravel(),'y.', alpha = .75)
    for i in range(len(Dobs)):
        plot(Dobs[i], Bobs[i], ms[i],ms = 12, label = labs[i])
    legend(loc = 'lower right')
    axisFontSize()
    xlab(r'D$_{meas}$, pix')
    ylab(r'blur number, $\frac{A_{tail}}{A_{central}}$')
    tit('Topography of the Observation Vector Space')
    
    
    ########################################################
    ########### plot probabilities #########################
    ########################################################
    pfig = fig()
    subplot(2,3,1)
    pfig.subplots_adjust(left  = 0.125/2, right = 0.99, bottom = 0.075, 
                         top = 0.95, wspace = 0.4 , hspace = 0.4)    

    for i in range(len(P)):
        subplot(2,3,i+1)
        contourf( fm.w,fm.D/Dobs[i], P[i], 30, cmap = cm.Accent)
        xlim(fm.w[0], fm.w[-1]) 
        ylim((fm.D/Dobs[i])[0], (fm.D/Dobs[i])[-1]) 
        ylab(r'diameter corrector,$\frac{D_{true}}{D_{meas}}$')
        xlab( 'blur kernel width, pix')
        tit(labs[i])
        axisFontSize()
    ########################################################


########################################################
########### plot quantities of interest#################
########################################################
    qfig = fig()
    subplot(2,2,1)
    qfig.subplots_adjust(left  = 0.125/2, right = 0.99, bottom = 0.05, 
                         top = 0.925, wspace = 0.15 , hspace = 0.35)
    plotList = [fm.Dmeas, fm.B, fm.I/fm.Imeas, fm.D[...,newaxis]/fm.Dmeas]
    labList = [r'D$_{meas}$, pix', r'blur number, $\frac{A_{tail}}{A_{central}}$',
                r'intesity corrector,$\frac{I_{true}}{I_{meas}}$',
                r'diameter corrector,$\frac{D_{true}}{D_{meas}}$']
    for i in range(len(plotList)):
        subplot(2,2,i+1)
        contourf( fm.w,fm.D, plotList[i], 50, cmap = cm.copper); colorbar()
        contour(fm.w,fm.D, plotList[i], [1.], colors = 'r',linewidths = 2)
        xlim(fm.w[0], fm.w[-1]) 
        ylim(fm.D[0], fm.D[-1]) 
        xlab('blur kernel width, pix')
        ylab('true diameter, pix')
        tit(labList[i])
#########################################################
    
        #mlm.surf3(fm.D,fm.w,  P[i],axeLab = ['true D, pix', 'blur kernel width, pix', 'p(D,w)'], f = mlm.fig(labs[i]))
    
#    d = fm(fm.D, fm.w)
    
#    
#    whisker = blur.whisker()
#    
#    fig()
#    plot(blur.w, blur.p)
#    fig()
#    plot(blur.w, blur.f)
#    
#    mlm.surf3(fm.D,fm.w,  fm.B, f = mlm.fig('blur #'))
#    mlm.surf3(fm.D, fm.w, fm.Dmeas, f = mlm.fig('D_meas'))
#    mlm.surf3(fm.D, fm.w, fm.Imeas, f = mlm.fig('I_meas'))
#    mlm.surf3(fm.D, fm.w, fm.I, f = mlm.fig('I'))
    
    #mlm.scat(r_[fm.X[0], fm.X[1], fm.X[5], fm.X[5]], f = mlm.fig('scat I'))
#    
#    mlm.surf3(fm.D,fm.w,  fm.Brough - fm.B, f = mlm.fig('blur # err'))
#    mlm.surf3(fm.D, fm.w, fm.Dmeasrough-fm.Dmeas, f = mlm.fig('D_meas err'))

#    mlm.surf3(fm.D,fm.w,  d[1] - fm.B, f = mlm.fig('blur # interp err'))
#    mlm.surf3(fm.D, fm.w, d[0]-fm.Dmeas, f = mlm.fig('D_meas interp err'))
#  
    
    #mlm.scat(r_['0,2', fm.X[0],fm.X[1],fm.X[2], fm.X[2] ], f = mlm.fig('blur scat'))
    
#    mlm.ml.show()
    show()
    
    
    