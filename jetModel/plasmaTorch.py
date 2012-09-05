# ########################################################################################
#this module defines jetSim objects, which are meant to create, run and manage simulations of 
#turbulent plasma jets, which are carried out in c++ 
# ########################################################################################
import gasMixture
reload(gasMixture)
from gasMixture import *
from helpers import *
from numpy import *
import subprocess
import scipy
import pdb

class plasmaTorch(object):
	"""
	base class for torch models. they initialize with basic parameters 
	
	current    - DC current in the torch
	flowrate  - volumetric flowrate in slm
	gas            - an object that returns properties of the gas
	make          - designate SG-100, 9-MB, etc
	
	internally, they simulate a particular torch and compute the exit conditions for the
	given torch parameters
	
	derived objects must define the members rx, Vx, Tx, the temperature and velocity 
	at the exit plane locations rx
	"""
	def __init__(self,current, flowrate, gas, make = 'SG100'):
		
		self.current = current
		self.flowrate = flowrate
		self.gas = gas
		self.make = make
		self.inputConditionsString = "cur-" + str(current) +"_flo-"+str(flowrate)+ "_gas-"+gas.gasName()+ "_make-"+make
		#gonna need the exit diameter
		self.d = {
			'SG100' : .008,
			'9MB' : .008
		}[self.make]
		
		#make the default mesh
		self.rx = linspace(0,self.d/2,200)
		self.drx = diff(self.rx[:2])
		self.areax = 2*pi*(self.rx+self.drx/2)*self.drx 	#exit node area, m^2
		
		self.Vx = 1000+0*self.rx				#exit velocity, m/s		
		self.Tx = 10000+0*self.rx    		#exit temperature, K
	
		self.inputMassFlowrate = self.flowrate*self.gas.Density(300)/6e4  #[L/min]*[1m^3/1000L]*[1min/60s]*[kg/m^3] = kg/s
		
	def stateVector(self):
		return (zip(['current', 'flowrate']+self.gas.gasNames+['make'],
			[self.current, self.flowrate]+ self.gas.volumeFractions.tolist() +[self.make]))
			
	def jetInpFileStr(self, rx):
		"""
		return a tuple of (location, temperature, velocity), one location per line, in string form 
		describing the conditions at the exit
		"""
		
		r = array(self.rx.tolist()+[self.rx[-1]+diff(self.rx[:2]), 1.1*max(rx)])
		T = array(self.Tx.tolist()+[300, 300])
		V = array(self.Vx.tolist()+[0, 0])
		Tx = scipy.interpolate.interp1d(r,T)(rx)
		Vx = scipy.interpolate.interp1d(r,V)(rx)
		return (str(rx).replace('[','').replace(']','').split(),str(Tx).replace('[','').replace(']','').split(),str(Vx).replace('[','').replace(']','').split())
		
	

class ramshawAndChangTorch(plasmaTorch):
	"""
	a crude jet model based on the paper by ramshaw and chang, 
	"Numerical Simulations of Argon Plasma Jets Flowing into Cold Air"
	Plasma Chemistry and Plasma Processing, vol. 13, no. 2, 1993
	"""
	def __init__(self, current  = 900, flowrate = 35.4, gas = gasMixture() , make = 'SG100' , X0 = array([12658, 1096]), power = 12100 ):
		"""
		set defaults and some hardwired parameters. compute torch power.
		pass None for X0 and power to compute these from fixed model parameters
		"""
		
		#call the base class constructor
		plasmaTorch.__init__(self, current, flowrate, gas, make)
		self.X0 = X0
		
		# hardwired default parameters
		self.Twall = 700 #wall temp, from ramshaw and chang
		
		#obtained for pfender's case, ar68he32, T0 = 15000, v0 = 1000
		self.m = 3.23				#temperature exponent
		self.n = 5.47				#velocity exponent
		
		#estimate torch power from current if it is not given
		if power ==None:
			P = {
				'SG100' : array([0, 16500, 19250, 22080, 28000]),
				'9MB' : array([0, 16500, 19250, 22080, 28000])
			}[self.make]
			I = array([0, 600, 700, 800, 1000])
			
			#power and flowrate are used to fit the exit profile
			self.inputPower = scipy.interpolate.interp1d(I,P)(self.current)
		else:
			self.inputPower=power
		
		
		
	def solve(self):
		"""
		optimize the free parameters and use the optimized values to compute the temperature and 
		velocity profiles at the exit
		"""
		#pdb.set_trace()
		if self.X0 == None:
			fn = lambda x0: self.centerlineError(x0)
			xOpt, fxOpt, gr, fgr = scipy.optimize.brute(fn, ((5000,19999),(100,2000)), Ns = 100, full_output = True)
			self.X0 = array(xOpt)
			#xOpt = scipy.optimize.fmin_l_bfgs_b(fn,array([10000,1000]), approx_grad=1,bounds = [(5000,19999),(100,2000)],  factr=1e3, epsilon = 1e-8 )
			#self.X0 = xOpt[0]
		else:
			fn = lambda N: self.parameterError(N)
			xOpt = scipy.optimize.fmin_l_bfgs_b(fn,array([5,5]), approx_grad=1,bounds = [(1,inf),(1,inf)]  )
			self.m = xOpt[0][0]
			self.n = xOpt[0][1]
			
		self.Tx = (self.X0[0] - self.Twall)*(1-pow(self.rx/(self.d/2),self.m))
		self.Vx = self.X0[1]*(1-pow(self.rx/(self.d/2),self.n))
	
	def kineticPower(self,Tx,Vx):
		return sum(self.gas.Density(Tx)*self.areax*Vx**3/2)
	def enthalpicPower(self,Tx,Vx):
		return sum(self.gas.Density(Tx)*self.areax*Vx*self.gas.Enthalpy(Tx))
	def massFlowrate(self, Tx, Vx):
		return sum(self.gas.Density(Tx)*self.areax*Vx)
	def exitProfiles(self, T0,V0):
		s = tuple([self.rx.shape[0]] + ones(array(T0).ndim).tolist())
		rx = self.rx
		rx.reshape(s)
		#curve fits of ramshaw and chang
		Tx = (T0 - self.Twall)*(1-pow(rx/(self.d/2),self.m))
		Vx = V0*(1-pow(rx/(self.d/2),self.n))
		return (Tx,Vx)
	
	def centerlineError(self, X0):
		"""
		for fixed parameters m and n and trial values X0 of the centerline temperature and velocity, 
		compute the error in flowrate and energy
		"""
		#curve fits of ramshaw and chang
		self.Tx = (X0[0] - self.Twall)*(1-pow(self.rx/(self.d/2),self.m))
		self.Vx = X0[1]*(1-pow(self.rx/(self.d/2),self.n))
		
		#compute the error in flowrate and energy terms
		try:
			flowrateError = 		abs(self.inputMassFlowrate-sum(self.gas.Density(self.Tx)*self.areax*self.Vx))/self.inputMassFlowrate
			kineticPower = 		sum(self.gas.Density(self.Tx)*self.areax*self.Vx**3/2)
			enthalpicPower = 	sum(self.gas.Density(self.Tx)*self.areax*self.Vx*self.gas.Enthalpy(self.Tx))
			powerError = abs(self.inputPower - kineticPower-enthalpicPower)/self.inputPower
		except:
			pdb.set_trace()
		#flowrate is more reliable so weight it more heavily with alpha>.5
		alpha = .7
		e = (alpha*flowrateError + (1-alpha)*(powerError))
		
		return e
	
	def parameterError(self, N):
		"""
		for fixed centerline temperature and velocity, given trial values N of parameters m and n, 
		compute the error in flowrate and power
		"""
		
		#curve fits of ramshaw and chang
		self.Tx = (self.X0[0] - self.Twall)*(1-pow(self.rx/(self.d/2),N[0]))
		self.Vx = self.X0[1]*(1-pow(self.rx/(self.d/2),N[1]))
		
		#compute the error in flowrate and energy terms
		flowrateError = 		abs(self.inputMassFlowrate-sum(self.gas.Density(self.Tx)*self.areax*self.Vx))/self.inputMassFlowrate
		kineticPower = 		sum(self.gas.Density(self.Tx)*self.areax*self.Vx**3/2)
		enthalpicPower = 	sum(self.gas.Density(self.Tx)*self.areax*self.Vx*self.gas.Enthalpy(self.Tx))
		powerError = abs(self.inputPower - kineticPower-enthalpicPower)/self.inputPower
		
		
		#flowrate is more reliable so weight it more heavily with alpha>.5
		alpha = .7
		e = alpha*flowrateError + (1-alpha)*(powerError)
		
		return e
		
	def plotProfile(self):
		figure()
		
		subplot(1,2,1)
		plot(self.rx*1000, self.Tx)
		ylabel('exit temperature, K')
		xlabel('radius, mm')
		
		subplot(1,2,2)
		plot(self.rx*1000, self.Vx)
		ylabel('exit velocity, m/s')
		xlabel('radius, mm')
#-----------------------------------------------------------------

if __name__== "__main__":
	g = gasMixture()
	g.blend(['argon', 'helium'],[68, 32])

	finkeTorch = ramshawAndChangTorch()
	finkeTorch.solve()

	pfenderTorch = ramshawAndChangTorch(800, 69, g,'SG100' , array([15000, 1000]))
	pfenderTorch.solve()
	
	
		