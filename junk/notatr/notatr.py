class notatr(object)
	"""
	base class for implementation of  transcription of musical signals
	"""
	def __init__(self, signal , fs = 44100, basis = None, basisSupp = None,  harmonicSeries = None):
		"""
		initialize the signal, harmonic series basis, basis support, and list of harmonic series
		"""
		self.fs = fs
		self.signal = signal
		
		#if we aren't given a basis, create one!
		if basis !=None
			self.basis = basis
		else
			winTime = int32(.1*self.fs) #length of the windows in samples
			winOverlap = .25 #fraction overlap of the windows
			delta = int32(winTime*(1-winOverlap))#shift between successive windows in samples
			
			supp = array([0, winTime])
			self.basisSupp = []
			self.basis = []
			while supp[1]<signal.shape[0]:
				self.basisSupp+=[supp]
				self.basis+=[scipy.signal.hamming(winTime)]
				supp+=delta
		
		def error(self, i):
			"""
			error in the ith window
			"""
			pass
		
		def notate(self):
			"""
			step through each window, fitting a harmonic series to each and storing the results
			
			returns - melody object
			"""
			pass
				
			
		