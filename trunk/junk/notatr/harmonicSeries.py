class harmonicSeries(object):
	"""
	fundamental frequencies, number of partials, partial detunings, and amplitudes define a harmonic series
	"""
	def __init__(self, fundamentalFrequencies = atleast_1d(0), 
							nPartials = atleast_1d(1), 
							detunes = None, 
							amplitudes = None, 
							maxNotes = 1,
							maxPartials = 10):
		"""
		set the datamembers:
		fundamentalFrequencies - length K vector of fundamental frequencies, one per note, omega_k
		nPartials - length K vector giving the number of partials in each note m_k
		detunes - maxNotes x maxPartials matrix of detuning frequencies for note k, partial m
		amplitudes - maxNotes x maxPartials matrix of amplitudes for each note k, partial m
		"""
		#set the fundamental frequencies and number of partials for each freq
		self.fundamentalFrequencies = fundamentalFrequencies
		self.nPartials = nPartials 
		
		assert self.maxNotes>self.fundamentalFrequencies.shape[0], 'nr. of fundamental frequencies exceeds maxNotes!'
		
		#set the maximum number of notes and partials allowed for this series
		self.maxNotes = maxNotes
		self.maxPartials = maxPartials
		
		#if we are not passed detunes and amplitudes matrices, allocate and fill them with zeros
		sz = (self.maxNotes, self.maxPartials)
		if detunes!=None
			assert detunes.shape == sz, 'detunes matrix dimension mismatch! \nshould be maxNotes x maxPartials'
			self.detunes = detunes
		else
			self.detunes = zeros((self.maxNotes, self.maxPartials))
		
		if amplitudes!=None
		assert amplitudes.shape == sz, 'amplitudes matrix dimension mismatch! \nshould be maxNotes x maxPartials'
			self.amplitudes = amplitudes
		else
			self.amplitudes = zeros((self.maxNotes, self.maxPartials))
		
		#preallocate for the mask of the frequency matrix
		self.partialMask = zeros((self.maxNotes(), self.maxPartials()))==1
		
	def M(self, M = None):
		"""
		max number of partials
		"""
		if M ==None:
			return np.max(self.nPartials)
			
	def K(self):
		"""
		number of notes
		"""
		return self.fundamentalFrequencies.shape[0]
	
	def mkSlice(self):
		"""
		slice out the valid K notes and at most M partials from any full maxNotes x maxPartials matrix
		"""
		return (slice(self.K()), slice(self.M) )
		
	def mask(self):
		"""
		return the boolean matrix that is true where the kth note, mth partial should have nonzero amplitudes
		"""
		#initialize the mask to true
		self.partialMask[:] = True
		
		#set mask false where there's no mode
		self.partialMask[self.nPartials[:,newaxis]<(arange(self.M())+1)] = False
		self.partialMask[self.K():] = False
		
		return self.partialMask
		
	def detunesMask(self):
		"""
		to make sure the first detune is always zero
		"""
		#first column of detunes must always be zero
		return tile(arange(self.M())[newaxis,:], (self.K(), 1))>0
		
		
	def freqs(self):
		"""
		return kth note, mth partial's frequency in a matrix
		"""
		#compute the actual frequencies of the series
		return self.fundamentalFrequencies[:,newaxis]*(arange(self.M())+1)+detunes[self.mkSlice()]*self.detunesMask
	
	def spectrum(self):
		"""
		return the complex spectrum of the series
		"""
		virtualWarn()
		return None
		
	def spectralPower(self):
		"""
		return the spectral power of the series
		"""
		virtualWarn()
		
		return None
		
	def freqInHz(self, freqs):
		"""
		given 
		freqs - array of frequencies in the harmonicSeries object's units
		return freqs in units of Hz
		"""
		virtualWarn()
		return None
	def timeSeries(self):
		"""
		return a time series representing the signal
		"""
class integralHarmonicSeries(harmonicSeries):
	"""
	harmonic series projected on a discrete basis of sinusoids
	"""
	def __init__(self, samplingFrequency = 44100, nfft = 2**14, windowFn = scipy.signal.hamming, winTime = .1,winOverlap = .25, **kwargs):
		"""
		samplingFrequency - in Hz
		nfft - number of points to use in the fft, determines resolution in frequency domain
		windowFn() - the fft of the window function, shifted by (m*fundamentalFreqs+detunes) will be used to 
						represent the sinusoids
		winTime - length of the window in s
		"""
		
		harmonicSeries.__init__(self, **kwargs)
		self.fs = samplingFrequency
		self.nfft = nfft
		self.winTime = winTime
		self.windowFn = windowFn
		self.setWindowFFTs()
		
	def setWindowFFTs(self):
		"""
		pre-compute the fft's of the shifted windows
		"""
		#nr of samples in the window
		nWin = self.winTime*self.fs
		#create the window
		win = self.windowFn(nWin)
		
		#fft the window
		winfft = fftshift(fft(win))
		
		#window energy
		e = cumsum(real(conj(winfft)*winfft))
		e/=e[-1]
		
		#select the central 98% by energy and construct to have odd length
		idx = np.argmin(abs(e-.01))
		width = round(double(winfft.shape[0]-2*idx)/2)
		center = np.argmax(abs(winfft))
		w = winfft
		w[0:center-width] = 0
		w[center+width+1:] = 0
		#the selected region of the window's fft
		w = winfft(range(center-width, center+width+1))
		
		#preallocate the shifted windows array
		self.shiftedWindows = zeros((self.nfft/2,self.nfft/2))
		
		#construct the initial slice
		iSly = int32(arange(self.nfft/2-1, self.nfft))
		
		#compute and cache the shifted windows
		for i in range(self.nfft):
			#from iSly, construct shiftedSly and winSly 
			shiftedSly = iSly-i
					
			#set the ith shifted window
			self.shiftedWindowFFTs[i] = w[winSly]
		
	def shifts(self):
		"""
		return indices into the rows of self.shiftedWindowFFTs for each 
		frequency in the sum
		"""
		#assumes freqs and detunes are in indicial units of the fft
		return int32(self.freqs)
	
	def spectrum(self):
		"""
		return the complex spectrum of the series
		"""
		idx = self.shifts()
		valid = self.partialMask[self.mkSlice()].nonzero()
		return sum(self.amplitudes[self.mkSlice()][valid][:,newaxis]*self.shiftedWindowFFTs[idx[valid]])
	
	def spectralPower(self):
		"""
		return the spectral power of this series
		"""
		S = self.spectrum()
		
		return S*conj(S)
		
		
