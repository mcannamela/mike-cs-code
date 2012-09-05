from numpy import *
import scipy
scipy.pkgload()
# Importing Mayavi's mlab as 'ml'
#try:
#    from mayavi import mlab as ml
#except ImportError:
#    from enthought.mayavi import mlab as ml
    
from pylab import *

import pdb
import os
import helpers
from helpers import indicialOneD_LUT
    
class gasMixture(object):
    """
    model the transport properties of a blend of gasses
    """

    def __init__(self,gasNames=['argon'], volumeFractions=[1], molFractions = None):
        """
        initialize a default gas
        """
        
        self.pIdx = {'Temperature':0, 'Conductivity':1, 'Density':2,'Enthalpy':3,'SpecificHeat':4,'Viscosity':5}
        self.gasMolarMass = {'argon':39.95, 'hydrogen':1.0079,'helium':4.,'nitrogen':2*14.01,
                                                'oxygen':2*16., 'air':28.966, 'zirconia':131.,
												'zirconiaSmooth':131.,'zirconiaGauss':131., 
                                                'zirconiaConstant':131., 'simpleGas':131.}
        self.gasVanDerWaalsRadius = {'argon':188e-12, 'hydrogen':1.91e-10,'helium':2.11e-10,'nitrogen':2.25e-10,
                                                'oxygen':2.06e-10, 'air':2.21e-10, 'zirconia':159.e-12,
												'zirconiaSmooth':159.e-12,'zirconiaGauss':159.e-12, 
                                                'zirconiaConstant':159.e-12, 'simpleGas':159.e-12}
        self.gasNames = gasNames
        self.firstCall = True
        if volumeFractions!=None:
            self.blend(gasNames, volumeFractions)
        else:
            self.blend(gasNames, None, molFractions)
        
    def __getstate__(self):
        """
        need to exclude all lambda objects
        """
        black =['TemperaturePol', 'Viscosity', 'Enthalpy', 'SpecificHeat', 'Density', 'Conductivity']
        S = []
        K = []
        for k in self.__dict__.keys():
            if k in black:
                continue
            K+=[k]
            S+=[self.__dict__[k]]
            
        return dict(zip(K,S))
        
    def __setstate__(self, D):
        """
        gonna have to call the methods that reconstruct the lambdas, which cannot be pickled
        """
        self.__dict__=D.copy()
        self.setCallbacks()
    
    def meanMolarMass(self):
        m = 0.
        cnt = 0
        for g in self.gasNames:
            m+=self.molFractions[cnt]*self.gasMolarMass[g]
            cnt+=1
        return m
        
    def meanMolecularRadius(self):
        v = 0.
        cnt = 0
        for g in self.gasNames:
            v+=4.*pi*self.molFractions[cnt]*self.gasVanDerWaalsRadius[g]**3/3
            cnt+=1
        return (3.*v/(4*pi))**(1./3)
        
    def blend(self, gasNames, volumeFractions, molFractions = None):
        """
        blend the specified gasses together in the specified proportions
        
        gasNames - list of gas names
        volumeFractions - arrays of volume fractions of the various gasses
        """
        #ensure the single gas properties are read in
        if all(array(gasNames)==array(self.gasNames)) and not self.firstCall:
            readOn = False
        else:
            readOn = True
            if type(gasNames)==type(str()):
                self.gasNames = [gasNames]
            else:
                    self.gasNames = gasNames
        
        if readOn:
            self.read()

        assert not(volumeFractions ==None and molFractions ==None)
        if volumeFractions !=None:
            self.volumeFractions = float64(volumeFractions)
            self.volumeFractions/=sum(self.volumeFractions)
            self.molFractions = self.vol2mol(self.volumeFractions)
        else:
            self.molFractions = molFractions
            self.molFractions/=sum(self.molFractions)
            self.volumeFractions= self.mol2vol(self.molFractions)
        
        #weighted average is fine for most properties
        self.mixProps = average(self.gasProps, 2, array(self.molFractions*self.density))
        self.mixProps[1,:] = 0;
        self.mixProps[-1,:] = 0;
        
        #deal with viscosity and conductivity
        for i in range(self.gasProps.shape[1]):
            #viscosity
            self.mix('Viscosity',i)
            #conductivity
            self.mix('Conductivity', i)
            
        self.mixProps = float32(self.mixProps)
        self.setCallbacks()
        self.firstCall = False
     
    def vol2mol(self, volFrac):
        """
        convert volumeFractions to mole fractions
        """
        self.setDensity()
        
        #molar masses
        self.molarMass = array([self.gasMolarMass[g] for g in self.gasNames] )
        
        #mole fractions from volume fraction, density, and molar mass
        self.molFractions = self.density*volFrac/self.molarMass; #mole fraction
        self.molFractions/=sum(self.molFractions)
        return self.molFractions
        
    def mol2vol(self, volFrac):
        """
        convert mole fractions to volume fractions
        """
        self.setDensity()
        
        #molar masses
        self.molarMass = array([self.gasMolarMass[g] for g in self.gasNames] )
        
        #mole fractions from volume fraction, density, and molar mass
        try:
            self.volumeFractions= self.molFractions*self.molarMass/self.density 
        except:
            pdb.set_trace()
        return self.volumeFractions
        
    def setDensity(self):
         #get the density at room temp
        self.density = zeros(len(self.gasNames))
        I = self.pIdx
        for i in range(len(self.gasNames)):
            f = scipy.interpolate.interp1d(self.gasProps[I['Temperature'],:,i],self.gasProps[I['Density'],:,i] )
            self.density[i] = f(300)
            
    def setCallbacks(self, mixProps = None):
        #callbacks
        I = self.pIdx
        T = self.mixProps[I['Temperature'],:]
        if mixProps ==None:
            mixProps = self.mixProps
            
       
        tempInterpolator = scipy.interpolate.interp1d
        #theInterpolator = indicialOneD_LUT
        theInterpolator = scipy.interpolate.interp1d

        
#        self.TemperaturePol = tempInterpolator(r_[-1e10, mixProps[I['Enthalpy'],:],1e37 ],r_[0,T,T[-1]] )
        self.TemperaturePol = tempInterpolator(mixProps[I['Enthalpy'],:],T )
        self.Conductivity = theInterpolator(T,mixProps[I['Conductivity'],:])
        self.Density = theInterpolator(T,mixProps[I['Density'],:])
        self.Enthalpy = theInterpolator(T,mixProps[I['Enthalpy'],:])
        self.SpecificHeat =  theInterpolator(T,mixProps[I['SpecificHeat'],:])
        self.Viscosity =  theInterpolator(T,mixProps[I['Viscosity'],:])
    
    def Temperature(self, *args):
        return float32(self.TemperaturePol(*args))
            
    def getCallbacks(self):
        """
        return a dictionary of the property fetching callbacks
        """
        return dict(zip(['Temperature', 'Conductivity', 'Density', 'Enthalpy', 'SpecificHeat', 'Viscosity'],
                                  [self.Temperature, self.Conductivity, self.Density, self.Enthalpy, self.SpecificHeat, self.Viscosity]))
    
    def phi(self, mu):
        """
        Wilke 1949, "Gas Mixture Viscosities"
        eqn. 14
        
        compute matrix phi needed find the conductivities and viscosities of the mixture
        """
        n = len(self.gasNames)
        Phi = zeros(n)
        I = self.pIdx
        
        #make column vectors of component properties and molar masses
        M_i = mat(self.molarMass).T
        mu_i = mat(mu).T
        
        #need matrices of ratios between component props and molar masses
        M_ji = (1/M_i)*M_i.T  #numerator in the columns
        M_ij = M_i*(1/M_i).T  #denom in the columns
        mu_ij = mu_i*(1/mu_i).T # numerator in the columns
        
        #compute the numerator 
        Phi = pow(            1 + pow( array(mu_ij)  ,.5)*pow( array(M_ji)  ,.25)        ,2)
        
        #divide by the denominator
        Phi/= (4/pow(2,.5))*pow(1+ array(M_ij),.5)
        
        #enforce zeros on the diagonal
        Phi[range(Phi.shape[0]), range(Phi.shape[0])] = 0;
        
        return mat(Phi)
        
    def mix(self, mixProp, i):
        """
        use Wilke's formula to compute the properties of the mixture
        """
        mu = squeeze(self.gasProps[self.pIdx[mixProp], i,:])
        Phi = self.phi(mu)
        denom = 1+squeeze(array(Phi*mat(self.molFractions).T))/self.molFractions
        self.mixProps[self.pIdx[mixProp], i] = sum(mu/denom)
        
    def read(self):
        """
        read the component gas properties from text files 
        """
        gasDirs = ["E:\\labCODE\\gasModel\\",  "D:\\CODE\\gasModel\\",  
                   "/media/SW_Preload/CODE/gasModel",
                   "/media/raidArray/CODE/gasModel",
                   "M:\\CODE\\gasModel\\"]

        for gasDir in gasDirs:
            try:
                g = []
              #  print gasDir
           #     print os.path.exists(gasDir)
                for gas in self.gasNames:
                    f = open(os.path.join(gasDir,gas+'.txt'))
                    g += [array([array(double(ss.replace('\n','').replace(',','').split())) for ss in f.readlines()])]
                    f.close()
                break
            except IOError:
                continue
        
        G = zeros([g[0].shape[0], g[0].shape[1], len(g)])
        
        for i in range(len(g)):
            G[:,:,i] = g[i]
            
        self.gasProps = array(G)
                
    def write(self):
        """
        write the mixture properties to their own file, suitable for passing to the jet model
        """
        fname = self.gasName()+'.gas'
        f = open(fname,'w')
        f.write(self.jetInpFileStr())
        f.close()
        return fname
    
    def gasName(self):
        """
        build a name for this gas mixture by concatenating the component gas names and their volume fractions
        """
        fname = ""
        for i in range(len(self.gasNames)):
            fname += self.gasNames[i]+str(self.volumeFractions[i]).replace('.','')
        return fname
    
    def jetInpFileStr(self):
        """
        format the property arrays for the jet model input file
        """
        return str(self.mixProps).strip(']').replace('\n','').replace('[','').replace(']','\n')
    
    def propPlot(self):
        """
        plot the mixture properties
        """
        I = dict(zip(self.pIdx.values(), self.pIdx.keys()))
        for i in range(self.mixProps.shape[0]-1):
            figure()
            plot(self.mixProps[self.pIdx['Temperature'],:], self.mixProps[i+1,:])
            xlabel('Temperature, K')
            ylabel(I[i+1])
        show()

defaultMaterial = gasMixture('zirconiaGauss')        

if __name__== "__main__":
    close('all')
    gm = gasMixture()
    gm.blend(['argon', 'helium'], [.68, .32])
    gm.write()
    f = open('argon068helium032.gas')
    g = array([array(double(ss.replace('\n','').replace(',','').split())) for ss in f.readlines()])
    f.close()
    labs = ['temperature', 'conductivity', 'density', 'enthalpy', 'specific heat', 'viscosity']
    for i in range(g.shape[0]):
        figure(i+1)
        #plot(gm.mixProps[i,:], g[i,:])
        plot(gm.mixProps[0], gm.mixProps[i])
        xlabel('temperature')
        ylabel(labs[i])
    #show()
        
        