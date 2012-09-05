from numpy import *
import scipy

class gaussCompositePdf_1d(object):
	"""
	1-d gaussian mixture
	"""
	
	def __init__(s, pars = dict(zip(['E','S','P'],[linspace(-5,5,3), linspace(.5,1.5,3), array([.6, .3, .1])])) ):
		"""
		initialize given a dictionary of parameters:
		
		E - expected value
		S - std dev
		P - prior probability (scale)
		
		each value of the dictionary is an array of parameters; the length of those arrays must match and is 
		equal to the number of components in the mixture
		"""
		
		#set parameters, update other data members
		s.set_pars(pars)
		
				
	def set_pars(s, p):
		"""
		set new values of the parameters from the dictionary p
		"""
		assert type(p)==type(dict()), 'ERROR: pars must be set with a dictionary with keys E,S,P'
		
		s.pars = p
		s.N = s.E().shape[0]
		s.set_support()
		s.set_randomVars()
	
	def set_support(s):
		"""
		update the support of the distribution based on the current parameters
		"""
		mn = argmin(s.E())
		mx = argmax(s.E())
		s.supp = (s.E()[mn]-3*(s.S()[mn]), s.E()[mx]+3*(s.S()[mx]))
	
	def set_randomVars(s):
		"""
		build up a list of normal random variables using the current parameters
		"""
		s.rv = [scipy.stats.norm(loc = s.E()[i], scale = (s.S()[i])) for i in range(s.N)]
	
	def E(s, val = None):
		"""
		gets/sets the array of means
		"""
		if val ==None:
			return s.pars['E']
		else:
			assert val.shape ==s.pars['E'].shape, 'ERROR: cannot change the number of parameters with this function, use set_pars() instead'
			s.pars['E'] = val
		s.set_pars(s.pars)
		
	def S(s, val = None):
		"""
		gets the array of variances
		"""
		if val ==None:
			return s.pars['S']
		else:
			assert val.shape ==s.pars['S'].shape, 'ERROR: cannot change the number of parameters with this function, use set_pars() instead'
			s.pars['S'] = val
		s.set_pars(s.pars)
	def P(s, val = None):
		"""
		gets the array of priors
		"""
		if val ==None:
			return s.pars['P']
		else:
			assert val.shape ==s.pars['P'].shape, 'ERROR: cannot change the number of parameters with this function, use set_pars() instead'
			s.pars['P'] = val
		s.set_pars(s.pars)
	def parsVec(s, parsVec = None):
		"""
		cat the parameters into one vector and return it if parsVec is None, 
		else set the parameters according to parsVec
		"""
		if parsVec==None:
			return hstack([s.E(), s.S(), s.P()])
		else:
			assert mod(parsVec.shape[0],3)==0, 'ERROR: parameter vector must be divisible by 3!'
			n = parsVec.shape[0]/3
			p = reshape(parsVec,(3,n))
			s.set_pars(dict(zip(['E','S','P'],[p[0], p[1], p[2]])))
	def pdf(s, x = None):
		"""
		evaluate the pdf at the points x
		"""
		if x == None:
			x = linspace(s.supp[0], s.supp[1], 100)
		
		pdf = zeros_like(x)
		i = 0
		for rv in s.rv:
			pdf+= rv.pdf(x)*s.P()[i]
			i+=1
		return (pdf,x)
	def componentPdf(s,x = None):
		"""
		evaluate the component pdfs at x
		"""
		if x == None:
			x = linspace(s.supp[0], s.supp[1], 100)
		
		p = zeros((x.shape[0], s.N))
		i = 0
		for rv in s.rv:
			p[:,i] = s.P()[i]*rv.pdf(x)
			i+=1
		return (p,x)
		
	def sample(s, n):
		"""
		approximate sampling from the composite distribution by sampling the consituent distributions independently
		"""
		x = array([])
		for i in range(s.N):
			x=hstack((x,scipy.stats.norm(loc = s.E()[i], scale = s.S()[i]).rvs(ceil(n*s.P()[i]))   ))
		
		return x[:n]
	def rvs(s, n):
		"""
		alias for sample to be compatible with other distribution objects
		"""
		s.sample(n)
		
if __name__== "__main__":
	close("all")
	gm = gaussCompositePdf_1d()
	
	print gm.supp
	print gm.N
	print gm.E()
	print gm.S()
	print gm.P()
	
	x = None
	p,x = gm.pdf(x)
	ps,x = gm.componentPdf(x)
	xSam = gm.sample(10000)
	
	
	f = figure(figsize = (15,20))
	subplot(211)
	plot(x,p, 'k', linewidth = 2)
	plot(x,ps, 'x')
	xlim(gm.supp)
	xlabel('x')
	ylabel('p(x)')
	
	title('a 1-d gaussian mixture pdf\n means: '+ str(gm.E()) +
			'\n stdDev: '+str((gm.S()))+
			'\n priors: '+str(gm.P()))
	subplot(212)
	hist(xSam,100)
	xlim(gm.supp)
	xlabel('x')
	ylabel('count')
	title('histogram of samples from the mixture')
	f.subplots_adjust(top = .8)
	
	show()
	
	