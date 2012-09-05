from numpy import *
from scipy import interpolate
import gasMixture as gm
import pdb
import os
import cPickle as pkl
from pylab import *

class lavaPropertyGenerator(object):
    """
    generate a text file containing a FORTRAN data statement suitable for use as properties
    with LAVA
    """
    def __init__(self):
        """
        define temperatures at which to evaluate props, plus some strings we'll need to write the 
        files properly
        """
        self.propList = ['conductivity', 'enthalpy', 'radiationLoss', 'viscosity']
        
        self.propKeys = dict(zip(  self.propList, ['CND','HK', 'SSUBR', 'VIS']   ))
        self.conversionFactors = dict(zip(self.propList, [1e5, 2.39e-7, 1e7, 10]  ))        
        self.fileExtensions = dict(zip(self.propList, ['.K', '.H', '.R', '.V']))        
        
        self.lavaTemp = linspace(0, 20000, 201)
        self.dataInterval = array([1,50])
        self.si2cgs = True #if you supply properties in SI, set this to true since LAVA works in CGS
        
    
    def __call__(self, propArr, prop, propIdx, gasName, convert = True, molarMass = None, path = None):
        """
        propArr - 2xN array, propArr[0] specifies temperatures, propArr[1] the value of the property
        prop - key into propList, names the property we are working with 
        gasName - string naming the gas will become the filename
        convert - boolean specifying whether to convert the units of propArr
        molarMass - if we are specifying enthalpy in SI units, LAVA expects kcal/mole so we need to know molar mass
        path - where to put the new file 
        """
        assert not (molarMass==None and prop=='enthalpy' and convert==True), "Error, conversion specified for enthalpy but molarMass not given. Either pass the array in kcal/mole and specify convert =False or pass the molarMass"
        if molarMass!=None and prop=='enthalpy':
            self.conversionFactors['enthalpy']*=molarMass
        C = dict(zip([True, False],[self.conversionFactors[prop], 1]))
        p = C[convert]*interpolate.interp1d(propArr[0], propArr[1], kind = 'cubic')(self.lavaTemp)
        
        if path!=None:
            fname = os.path.join(path, gasName+self.fileExtensions[prop])
        else:
            fname = gasName+self.fileExtensions[prop]
            
        with open(fname, 'w') as f:
            idx = arange(5)
            cnt = 0
            for i in range(4):
                f.write(self.dataStatement(prop, propIdx, i))
                f.write(self.dataLine(p[5*cnt+idx], isFirst=True )     )
                cnt+=1
                for j in range(8):
                    f.write(self.dataLine(p[5*cnt+idx]))
                    cnt+=1
                
                if i==3:
                    f.write(self.dataLine(p[5*cnt+idx])     )
                    f.write(self.dataLine(atleast_1d(p[-1]), isLast=True )     )
                else:
                    f.write(self.dataLine(p[5*cnt+idx], isLast=True )     )
  

    def dataStatement(self, prop = 'conductivity', propIdx = None, interval = 0):
        """
        generate the string that makes 
        """
        if interval == 3:
            ival = 50*interval+self.dataInterval+array([0,1])
        else:
            ival = 50*interval+self.dataInterval
        
        assert propIdx!=None, "propIdx not specified, it must match component index in LAVA!" 
        return '      DATA ('+self.propKeys[prop]+'(N,%d),N=%d,%d)\n'%(propIdx, ival[0],ival[1])

    def dataLine(self, arr, isFirst = False, isLast = False):
        """
        format the numbers in array arr into a string suitable for the file
        """
        assert arr.shape[0]<=5, "argument to arr must be pre-sliced to have at most length 5!"
        f = lambda n:('%.4E'%n).replace('E','D')+', '
        if isFirst:
            L = '     & /'
        else:
            L = '     &  '
            
        for x in arr.tolist():
            L+= f(x)
            
        if isLast:
            L=L[:-2]
            L+='/\n'
        else:
            L+='\n'
        return L
        
def cookTemperature(x):
    return x
def cookVel(x):
    return x/100.
def cookDist(x):
    return x/100.
def cookDensity(x):
    return x*1000
def cookConductivity(x):
    return x*1e-5
def cookViscosity(x):
    return .1*x
def cookSpeciesDensity(x):
    return x*1000
def cookSpeciesMoleFraction(x):
    return x
    
class lavaPropertyReader(object):
    def __init__(self, infile,prop,  cookFun):
        self.T = linspace(0,20000,201)
        self.prop = prop
        X = []
        with open(infile,'r') as f:
            for L in f.readlines():
                L= L.strip().replace('D','e').replace(',',' ').replace('/','')
                if L[0]!='&':
                    continue
                else:
                    X+= L[1:].split()
                pdb.set_trace()
        self.X = cookFun(float64(X))
        
    def __call__(self, outFile):
        with open(outFile, 'wb') as f:
            pkl.dump(dict(zip( ['pressure', 'temperature', self.prop],[1.,self.T, self.X])), f,-1)
        return (self.T, self.X)
    
if __name__=='__main__':
    propGen = False
    propRead = True
    
    if propRead:
        indir = 'E:\\labCODE\\LAVA\\XPTY'
        outdir = 'E:\\labCODE\\gasModel'
        infiles = ['Air.V' ]
        outfiles = ['airViscosity.pkl']
        props = ['viscosity']
        cookfuns = [cookViscosity]
        for i in range(len(infiles)):
            LPR = lavaPropertyReader(os.path.join(indir,infiles[i]), props[i], cookfuns[i])
            T,X = LPR(os.path.join(outdir,outfiles[i]))
            figure()
            plot(T,X)
            ylabel(props[i])
            xlabel('T')
        
    if propGen:
        #fetch up the properties of argon, hydrogen, air
        g = gm.gasMixture()
        h = gm.gasMixture(['hydrogen'])
        with open(os.path.join(os.pardir, 'gasModel', "airConductivity_p=1_nT=201.pkl"),'rb') as f:
            D = pkl.load(f)
            
        #construct the generator
        LPG = lavaPropertyGenerator()
        
        #generate the fortran files
        LPG(r_['0,2',g.mixProps[0], g.mixProps[1]], 'conductivity', 1,'mikesArgon')
        LPG(r_['0,2',h.mixProps[0], h.mixProps[1]], 'conductivity', 2,'mikesHydrogen')
        LPG(r_['0,2',D['temperature'], D['conductivity']], 'conductivity', 3,'mikesAir')
    
