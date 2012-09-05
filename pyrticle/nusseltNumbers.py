from numpy import *
from pylab import *

def nu_nonEvap0(Re, Pr = .7):
    return 2+(.4*Re**.5-.06*Re**.6666)*Pr**.3333
    
def nu_nonEvap1(Re_, Pr = .7):
    Re = array(Re_)
    f = ones_like(Re)
    f[Re>1] = Re**.077
    return 1+f*(1+Re*Pr)**.33333
    
def nuEvap0(Re, Pr = .7):
    return 2+.57*Re**.5*Pr**.3333
    
def evapFactor0(B_f):
    return 1./(1.+B_f)**.7
    
def generalSpaldingMassNumber(T_gas,T_surf, cp_vapor, L_v, q_conv, q_vap, _q_rad):
    return (cp_vapor*(T_gas-T_surf)/L_v)*(1-(q_conv-q_vap-q_rad)/q_conv)
    
def drag0(Re):
    return (24/Re)*(1+.1935*Re**.6305)
def drag1(Re):
    return (24/Re)*(1+.325*Re**.474)
def drag2(Re):
    return (24/Re)*(1+(Re**(.666))/6)
    
if __name__=="__main__":
    Re = linspace(1, 30, 50)
    figure()
    plot(Re, nu_nonEvap0(Re), label = 'nonEvap0')
    plot(Re, nu_nonEvap1(Re), label = 'nonEvap1')
    plot(Re, nuEvap0(Re), label = 'evap0')
    legend()
    xlabel('Re')
    ylabel('Nu')
    
    z = linspace(1e-10, 3, 1000)
    
    figure()
    plot(z, z/(exp(z)-1))
    xlabel('z')
    ylabel('Nu exp evaporation factor')
    
    Bf = linspace(0, 3, 1000)
    figure()
    plot(Bf, evapFactor0(Bf))
    xlabel('spalding heat transfer number')
    ylabel('nu evap factor')
    
    figure()
    plot(Re, drag0(Re), label = 'drag0')
    plot(Re, drag1(Re), label = 'drag1')
    plot(Re, drag2(Re), label = 'drag2')
    xlabel('Re')
    ylabel('C_d')
    legend()