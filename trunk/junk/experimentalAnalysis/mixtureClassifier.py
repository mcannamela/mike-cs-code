from compositePdf import *
from numpy import *
from numpy.random import rand
import scipy
import pdb

class mixtureClassifier_1d(object):
	"""
	break down 1d data into classes using a mixture model
	"""
	def __init__(s,x, N):
		"""
		set the data, initialize the compositePdf
		"""
		#must know number of components
		s.N = N
		
		#data should be in an array
		s.x = x.flatten()
		
		#support of the data
		s.supp = (min(s.x), max(s.x))
		
		#take random initial values for means, std devs, and prior distributions, encapsulating in a composite pdf
		s.pdf = gaussCompositePdf_1d(dict(zip(['E','S','P'],[s.supp[0]+s.supp[1]*rand(N),ones(N)*diff(s.supp)/6,ones(N)/N])))
		
		#use simple clustering to update initial values
		s.kMeans()
		
		
	def kMeans(s):
		"""
		perform some iterations of the kMeans algorithm to get estimates for the mean and std devs
		"""
		mu = s.pdf.E()
		
		while True:
			mu_old = mu
			idxTup = s.nearestMeanClassify(mu)
			mu = array([mean(s.x[idxTup[i]])   for i in range(s.N)])
			if sum(abs(mu-mu_old))<.01*diff(s.supp)[0]:
				break
		
		s.kMeansClusters = idxTup
		dev  = array([std(s.x[idxTup[i]])   for i in range(s.N)])
		pri = array([double(sum(idxTup[i]))/s.x.shape[0] for i in range(s.N)])
		#pdb.set_trace()
		s.pdf.parsVec(hstack((mu, dev, pri)))
		
	def nearestMeanClassify(s,mu):
		"""
		classify each point in the data according to the nearest mean
		"""
		D = abs(expand_dims(s.x,1)-mu)
		idx = argmin(D, axis = 1)
		return tuple([idx==i for i in range(s.N)])
		
	def minusLogLike(s,pars):
		"""
		compute the negative log-likelihood of the data given the parameters in pars
		"""
		s.oldPars = s.pdf.parsVec()
		s.pdf.parsVec(pars)
		p = s.pdf.pdf(s.x)[0]
		#print sum(s.pdf.parsVec()[-3:])
		return -sum(log(p))+1000*s.x.shape[0]*abs((1-sum(pars[-3:])))
	
	def optimize(s):
		"""
		find the MAP estimate of the parameters through estimation
		"""
		#pdb.set_trace()
		parsOpt = scipy.optimize.fmin(s.minusLogLike, s.pdf.parsVec(), maxfun = 1e4) 
		s.pdf.parsVec(parsOpt)
	
if __name__== "__main__":
	close("all")
	
	#construct a target distribution
	gm = gaussCompositePdf_1d()
	gm.E(array([-3,0,3]))
	gm.S(array([.5,1.5,1]))
	
	#draw samples from the target distribution and take their histogram
	xSam = gm.sample(1000)
	h,b = histogram(xSam, 100, normed=True)
	
	#create a classifier
	mc = mixtureClassifier_1d(xSam,gm.N)
	
	#print the true parameters and the guesses obtained from clustering
	print 'true means: '+str(gm.E()) + '\n kMeans: '+ str(mc.pdf.E())+'\n'
	print 'true devs: '+str(gm.S()) + '\n kMeansDevs: '+ str(mc.pdf.S())+'\n'
	print 'true priors: '+str(gm.P()) + '\n kMeansPriors: '+ str(mc.pdf.P())+'\n'
	
	#perform MAP estimate of the parameters
	mc.optimize()
	
	#print the MAP estimates
	print 'optimized means: ' +str(mc.pdf.E())
	print 'optimized devs: ' +str(mc.pdf.S())
	print 'optimized priors: ' +str(mc.pdf.P())
	
	#plot the results
	f = figure(figsize = (15,20))
	bar(b[:-1],h, width = diff(b[:2])[0], alpha = .5, label = 'samples')
	plot(gm.pdf()[1], gm.pdf()[0], 'k', linewidth = 2, label = 'truth')
	plot(mc.pdf.pdf()[1], mc.pdf.pdf()[0],  'r', linewidth = 2, label = 'guess')
	xlabel('x')
	ylabel('p(x)')
	legend()
	title('True posterior distribution and its MAP estimate from samples: ')
	
	
	
	
	
		