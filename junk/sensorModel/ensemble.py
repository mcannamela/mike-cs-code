import scipy
import helpers as hp
reload(hp)
import particle as part
reload(part)
import  mixtureClassifier as mixClass


normDist = scipy.stats.norm
lognormDist = scipy.stats.lognorm

meanD = 20e-6
stdD = 10e-6

class emittingParticleSet(list)
	"""
	dress up the list class with cached vectors of particle states, T,V,D
	"""
	def __init__(self, P):
		assert type(P[0]) ==type(part.particle()), 'particle set must be initialized with a list of particles!'
		list.__init__(P)
		self.T = array([p.T[-1] for p in P])
		self.V = transpose(array([p.velocity() for p in P]))
		self.D = transpose(array([p.mesh.xInt[[0,-1]] for p in P]))
		
	def hist(self, normOn = False)
		"""
		open up a figure, break out subplots, make 100 bin histograms of T,axialVelocity, and diameter
		"""
		#make sure to specify whether histogram is normed using the normOn argument
		pass 
		
	def temperature(self):
		return self.T.copy()
	def velocity(self):
		return self.V.copy()
	def lateralVelocity(self)
		return self.V[2].copy()
	def axialVelocity(self)
		return self.V[0].copy()
	def verticalVelocity(self)
		return self.V[1].copy()
	def diameter(self)
		return self.D[1].copy()
	def innerDiameter(self)
		return(self.D[0]).copy()
	
	
		

		
class inverseParticleDistributionSet(particleDistributionSet)
		"""
		infer distributions from a measured or computed sample of particles
		"""
		def __init__(particleList):
			"""
			use the particle list and a the inferDistribution method to back out distributions for the particle states
			
			particleList - a list of particle objects
			"""
			self.P = emittingParticleSet(particleList)
			
		def inferDistribution(self):
			"""
			estimate the pdf of the particle states using the particleList
			"""
			#for now, assume independently distributed states and use a mixtureClassifier_1d object to get the densities
			#initialize a mixClass.mixtureClassifier_1d object for each component of the particle state
			# T, V, D, d
			#note you have to pick the number of component distributions N, start with 3
			
			#call the optimize method of the classifer objects
			
			#call particleDistributionSet.__init__ method passing the 
			#pdf datamembers of the classifier objects as arguments 
			pass
			
class particleEnsemble(object):
	"""
	a collection of particles and their statistics
	"""
	def __init__(self, particleList = None, particleDistributions = None):
		"""
		set up the datamembers for the ensemble
		
		particleList - list of emittingParticle objects
		particleDistributions - a particleDistributionSet object
		
		one of the args should be None
		"""
		
		assert (particleList!=None or particleDistributions!=None)#use an assert statement to check that not both args ar None
		
		#if particleList is None, assert that particleDistributions is a particleDistributionSet object
			#define a datamember distributions = particleDistributions
		
		#if particleDistributions is None, use particleList to construct an inverseParticleDistributionSet
			#define a datamember P = particleList
			#define a datamemeber distributions = inverseParticleDistributionSet(particleList)
		pass
	def draw(self, N):
		"""
		use the draw method of the distributions datamember to get N random particles
		
		N - an integer number of particles to draw
		"""
		#set self.P = emittingParticleSet(self.distributions.draw(N) )
		pass
		
	def plot(self):
		"""
		show  histograms of the ensemble overlaid on plots of their respective pdfs
		"""
		#should be as simple as calling the hist method of the particle list then the plot method of 
		#the distributions datamember
	
	def distributionDistance(self, other):
		"""
		compute the distance between the distributions of other and self
		other - a particleEnsemble object
		"""
		#use the L1 norm for functions. 
		
		#evaluate the distributions of self and other on the same grid
		
		#subtract the evaluated distributions and take the absolute value
		
		#sum the absolute values of the differences of the distributions to get the distance
	
	def ensembleDistance(self, other)
		"""
		compute the ensemble distance between self and other
		other - a particleEnsemble object
		"""
		S = sloppySetComparator(self.P, other.P, self.P[0].norm)
		return S()
	def particleList(self):
		"""
		get the list of particles
		"""
		return self.P
		
	def setImages(self, I):
		"""
		set the array of images
		"""
		self.I = array(I)
	def frameImage(self)
		"""
		spit the sum of the individual particle images
		"""
		return sum(I,2)
		
class sloppySetComparator(object):
	"""
	class for computing various distances between two sets with non-strictly equivalent (sloppy) elements
	"""
	def __init__(self, S1, S2, elementNorm = None)
		"""
		initialize with two lists, S1 and S2, whose elements are all of the same type
		
		S1 - first set
		S2 - second set
		elementNorm - function that returns the norm of an element of the sets
		
		"""
		self.S = [S1, S2]
		self.N = [len(self.S1), len(self.S2)]
		self.I = [range(self.N[0]), range(self.N[1])]
		if elementNorm!=None:
			self.elementNorm = elementNorm
		else:
			self.elementNorm = lambda x: sqrt(sum((x)**2))
		self.Iprime = []
		self.setDistanceMatrix()
		
		self.reductionLoop()
	def __call__(self):
		return self.setDistance.copy()
	
	def reductionLoop(self)
		"""
		recursively reduce the distance matrix D until it has vanished
		then record the set difference and the error between the identified elements of the set
		"""
		
		for i in range(min(self.N)):
			self.reduce()
			
		self.setDifference = []
		if self.N[0]>self.N[1]:
			idx = 0
		elif self.N[1]>self.N[0]:
			idx = 1
		else:
			idx = None
		
		
			
		self.loss = 0
		if idx!=None:
			self.setDifference = set(self.I[idx]).difference(self.Iprime[idx])
			for s in self.setDifference:
				self.loss+=elementNorm(s)
		
		self.error = 0
		for I in self.Iprime:
			self.error+=self.distanceMatrix[I[0],I[1]]
			
		self.setDistance = self.loss+self.error
			
		
	def setDistanceMatrix(self)
		"""
		compute the distance matrix between the elements of the two sets
		"""
		self.D = zeros(self.N[0], self.N[1])
		for i in range(self.N[0]):
			for j in range(self.N[1])
				self.D[i,j] = self.elementNorm(self.S[0][i]-self.S[1][j])
		self.distanceMatrix = self.D.copy()
		
	def reduce(self):
		"""
		reduce the distance matrix by deleting the row and column corresponding to the minimum element
		"""
		i,j = unravel_index(np.argmin(self.D), self.D.shape)
		self.Iprime +=[(self.I[0][i], self.I[1][j])]
		
		s1 = range(i)+range(i+1,self.D.shape[0])
		s2 = range(j)+range(j+1,self.D.shape[1])
		
		self.D = self.D[s1,s2]
		self.I[0] = self.I[0][s1]
		self.I[1] = self.I[1][s2]
		
	
		
		
		
	
	