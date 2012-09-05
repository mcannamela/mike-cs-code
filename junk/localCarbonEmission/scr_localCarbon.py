from localCarbonEmission import *

dataFile = 'data.xls'

EL, T = readData(dataFile)



heatingThresholdTemperature = 65 #deg F

ECB_office = energyCarbonBudget(EL, elecIntensity_office)
ECB_apt = energyCarbonBudget(EL, elecIntensity_apartment)
ECB_house = energyCarbonBudget(EL, elecIntensity_house)

HCB_office = heatingCarbonBudget(T,
                                 heatingIntensity_office, 
                                 heatingThresholdTemperature)
HCB_apartment = heatingCarbonBudget(T,
                                 heatingIntensity_apartment, 
                                 heatingThresholdTemperature)
                                 
HCB_house = heatingCarbonBudget(T,
                                 heatingIntensity_house, 
                                 heatingThresholdTemperature)

L = [ECB_office, ECB_apt      , ECB_house,
     HCB_office, HCB_apartment, HCB_house]

E = array([x.energy() for x in L]).T
C = array([x.carbon() for x in L]).T

headers = ('officeElec_energyDensity aptElec_energyDensity houseElec_energyDensity'.split()+
          'officeHeating_energyDensity aptHeating_energyDensity houseHeating_energyDensity'.split()+ 
          'officeElec_carbon aptElec_carbon houseElec_carbon'.split()+
          'officeHeating_carbon aptHeating_carbon houseHeating_carbon'.split())

units = ('kWh/m^2 kWh/m^2 kWh/m^2 '.split()+
        'kWh/m^2 kWh/m^2 kWh/m^2 '.split()+
        'tonsC/kWh/yr tonsC/kWh/yr tonsC/kWh/yr'.split()+
        'tonsC/kWh/yr tonsC/kWh/yr tonsC/kWh/yr'.split())

fname = 'localCarbonEmissionCurves'
with open(fname+'.csv', 'w') as f:
    [f.write(h+', ') for h in headers]
    f.write('\n')
    [f.write(u+', ') for u in units]
    f.write('\n')
    for i,e in enumerate(E):
        [f.write(str(x)+', ') for x in e]
        [f.write(str(x)+', ') for x in C[i]]
        f.write('\n')
        
with open(fname+'.totals','w') as f:
    f.write('Total annual energy intensities, in kWh/m^2\n')
    f.write(':::::::::::::::::::::::Heating:::::::::::::::::::\n')
    f.write('office %d, apartment %d, house %d\n'%(heatingIntensity_office,
                                                    heatingIntensity_apartment,
                                                    heatingIntensity_house))
    f.write(':::::::::::::::::::::::Electricity:::::::::::::::\n')
    f.write('office %d, apartment %d, house %d\n'%(elecIntensity_office,
                                                    elecIntensity_apartment,
                                                    elecIntensity_house))
        
figure()
series = 'officeElec aptElec houseElec heatOffice heatApt heatHouse'.split()
m = '-r -k -y :r :k :y'.split()

midnightWeekly = float64(mod(arange(len(E.T[0])), 7*24)==0)
weeks = 7*24*arange(sum(midnightWeekly))

for i,e in enumerate(E.T):
    plot(e, m[i], label = series[i])
    
stem(weeks, ones_like(weeks))
ylim(0, amax(E))
xlabel('hour')
ylabel(r'energy density, $\frac{kWh}{m^2}$')
title('Hourly Energy Use \nblue lines mark midnight at the start of each week')
legend()

figure()
for i,e in enumerate(C.T):
    plot(e, m[i], label = series[i])
stem(weeks, ones_like(weeks))
ylim(0, amax(C))
xlabel('hour')
ylabel(r'carbon density, $\frac{tons C}{m^2}$')
title('Hourly Carbon Emission\nblue lines mark midnight at the start of each week')
legend()
