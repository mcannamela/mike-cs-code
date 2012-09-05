from numpy import *
from pylab import *
from numpy import round as around
import pdb
import xlrd
import xlwt

carbonFactor_gas = 11.7*(1./100)*(1/0.2931)*(1./2000) # lb/therm * therm/kbtu * kbtu/kWh * ton/lb
carbonFactor_oil = 20.4*(1./140)*(1/0.2931)*(1./2000) # lb/gal * gal/kbtu * kbtu/kWh * ton/lb
carbonFactor_elec = .9/2000# lb/kWh * ton/lb

heatingIntensity_office =    70*0.2931/0.09290304    #kbtu/ft2 * kWh/kbtu / m2/ft2
heatingIntensity_apartment=  50*0.2931/0.09290304     #kbtu/ft2* kWh/kbtu / m2/ft2

##well-insulated house
#heatingIntensity_house=      26*0.2931/0.09290304      #kbtu/ft2* kWh/kbtu / m2/ft2

#typical house
heatingIntensity_house=      40*0.2931/0.09290304      #kbtu/ft2* kWh/kbtu / m2/ft2

##leaky house
#heatingIntensity_house=      75*0.2931/0.09290304      #kbtu/ft2* kWh/kbtu / m2/ft2


elecIntensity_office = 20/.093 #kWh/ft2 /m2/ft2
elecIntensity_apartment = 8/.093 #kWh/ft2 /m2/ft2

##good house
#elecIntensity_house = 2.6/.093 #kWh/ft2 /m2/ft2

##typical house
elecIntensity_house = 3.8/.093 #kWh/ft2 /m2/ft2

##efficient house
#elecIntensity_house = 4.4/.093 #kWh/ft2 /m2/ft2


minDayLength = 9	
maxDayLength = 15.25
equinoxSunrise = 	6.75



        
    

def sunriseSunset(minDay = minDayLength, 
                 maxDay = maxDayLength, 
                 equinoxSunrise = equinoxSunrise):
        eqDay = 81.0
        t = linspace(0,1,365)-eqDay/365
        T = sin(t*2*pi)
        A = .25*(maxDay-minDay)
        
        sr = equinoxSunrise-A*T
        ss = equinoxSunrise+12+A*T
        return array([sr, ss])

class diurnalTemperatureModel(object):
    def __init__(self, SR):
        """
        SR - 2 x 365 array of sunrise SR[0], and sunset SR[1], hours starting 
            with Jan 1
        """
        
        self.sr = SR[0]
        self.ss = SR[1]
        
    def __call__(self, T_hi, T_lo,  d = None):
        """
        from Cesaraccio et al, "An improved method for determining degree-day values 
        from daily temperatures", Int. J. Biometeorology
        """
        
        if d==None:
            d = arange(len(T_hi))
        
        #important hours of the day
        Hn = around(self.sr[int32(d)])#sunrise
        H0 = around(self.ss[int32(d)])#sunset
        Hx = H0-4#max temperature
        Hp = Hn+24#next sunrise
        
        #min, max, and next day min temperatures
        Tn = T_lo
        Tx = T_hi
        Tp = r_[T_lo[1:], T_lo[0]]
        
        #empirical model shape parameter
        c = .39
        
        #sunset temperature
        T0 = Tx-c*(Tx-Tp)
        
        #shorthand for other model parameters
        a = Tx-Tn
        R = Tx-T0
        b = (Tp-T0)/(Hp-H0)**.5
        
        def hourlyTemperature(D):
            i = int32(D)
            
            t = arange(Hn[i], Hp[i])
            T = zeros_like(t)
            
            I = t<=Hx[i]
            J = (t>Hx[i])*(t<=H0[i])
            K = t>H0[i]
            
            T[I] = Tn[i]+a[i]*sin(.5*pi*(t[I]-Hn[i])/(Hx[i]-Hn[i]))
            T[J] = T0[i]+R[i]*sin(pi*.5*(1+.25*(t[J]-Hx[i])))
            T[K] = T0[i] + b[i]*(t[K]-H0[i])**.5
            
            return T.copy()
        
        TT = array([hourlyTemperature(D) for D in d])
        
        T = self.flattenT(TT)
        
        return r_[ T[-Hn[0]:], T[:-Hn[0]] ]
        
    def flattenT(self, TT):
        T_hourly = array([])
        for tt in TT:
            T_hourly = r_['0,1', T_hourly, tt]
        return T_hourly
        

def readData(fname):
    wb = xlrd.open_workbook(fname)
    sh = wb.sheet_by_index(0)
    
    Tlo = zeros(365)
    Thi = zeros(365)
    Tlo_next = zeros(365)
    E = zeros(sh.nrows-1)
    
    for i in range(1, sh.nrows):
        r = sh.row_values(i)
        
        if i<=365:
            Tlo[i-1] = r[1]        
            Thi[i-1] = r[2]
            
        
        E[i-1] = r[0]
        
    DTM = diurnalTemperatureModel(sunriseSunset())
    
    
    T = DTM(Thi, Tlo)
    
    return E, T

class energyCarbonBudget(object):
    def __init__(self, fractionSeries, annualUseDensity):
        """
        fractionSeries - a time series that will define our usage profile, 
                        ultimately it will be non-dimensionalized and will 
                        determine what portion of annual use is parceled out 
                        for each chunk of time
        annualUseDensity - the annual energy usage, in kWh/yr/m^2
        """
        self.fractionSeries = fractionSeries.copy()
        self.annualUseDensity = annualUseDensity
        
    
    def fraction(self):
        """
        
        """
        return self.fractionSeries/max([sum(self.fractionSeries),0])
        
    def energy(self):
        return self.fraction()*self.annualUseDensity
    
    def carbon(self, carbonFactor = carbonFactor_elec):
        return self.energy()*carbonFactor

class heatingCarbonBudget(energyCarbonBudget):
    def __init__(self, fractionSeries, annualUseDensity, heatingThreshold = 65):
        """
        heatingThreshold - temperature in deg F below which we will turn on the heat
        """        
        energyCarbonBudget.__init__(self,fractionSeries, annualUseDensity)
        self.heatingThreshold = heatingThreshold   
        
    def fraction(self):
        T = self.fractionSeries
        Th = self.heatingThreshold        
        isCold = float64((Th-T)>0)
        return (Th-T)*isCold/sum((Th-T)*isCold)
        
    def carbon(self, fuel = 'gas',
                   carbonFactor_oil = carbonFactor_oil,
                       carbonFactor_gas = carbonFactor_gas ):
        cf = {'oil':carbonFactor_oil, 'gas':carbonFactor_oil}[fuel]
        return cf*self.energy()
        
if __name__=="__main__":
    t = linspace(0,1, 100)
    Th = 63
    T = Th+5*sin(4*pi*t)
    Q = 100
    
    E = 12+10*sin(4*pi*t)+random(size = len(t))
    J = 100
    
    SR =  sunriseSunset()
    
    DTM = diurnalTemperatureModel(SR)
    
    T_lo = T-5*random(size = len(t))
    T_hi = T+5*random(size = len(t))
    T_hourly = DTM(T_hi, T_lo)

        
    th = arange(24*len(t))
    
    figure()
    plot(th, T_hourly)
    title('hourly temperature variation')
    
    
    
    ECB = energyCarbonBudget(E, J)
    HCB = heatingCarbonBudget(T,Q, Th)
    
#    figure()
#    plot(t, ECB.energy())
#    title('Electricity use per square meter')
#    
#    figure()
#    plot(t, ECB.carbon())
#    title('Carbon use per square meter due to electricity')
#    
#    figure()
#    plot(t, HCB.energy())
#    title('Heating per square meter')
#    
#    figure()
#    plot(t, HCB.carbon(fuel = 'gas'))
#    title('Carbon use per square meter due to gas heating')
#    
#    figure()
#    plot(t, HCB.carbon(fuel = 'oil'))
#    title('Carbon use per square meter due to oil heating')
    
#    figure()
#    plot(SR[0],'y')
#    plot(SR[1],'k')
#    title('sunrise sunset')
    
    show()
        
    