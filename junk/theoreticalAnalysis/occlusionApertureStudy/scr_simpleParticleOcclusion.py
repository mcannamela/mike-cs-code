#compute the fraction of occluding particles for a given sensor distance and vertical aperture
from numpy import *
import scipy
import pdb
sf = lambda figName_:  savefig(figName_+'.png', format = 'png', transparent = True)

class occlusionCalculator(object):
	"""
	"""
	def __init__(self, sensorDistance = 3., aperture = 1. , gridDim = 500):
		"""
		sensorDistance - plume widths to the sensor
		aperture - diameter of the aperture in plume widths
		"""
		self.plumeWidth = 6 #width of plume in standard devs
		
		#center of plume to sensor in standard devs
		self.sensorDistance = self.plumeWidth*sensorDistance
		
		#height of the lens in standard devs, divided by 2 since we use symmetry
		self.aperture =  self.plumeWidth*aperture/2
		
		#grid for computation of the pdf
		self.g = meshgrid(linspace(-3.5,3.5, gridDim), linspace(0,3.5,round(gridDim/2)))
		self.r = sqrt(self.g[0]**2+self.g[1]**2)
		
		self.setPdf(scipy.stats.norm.pdf(self.r))
		
	
	def setPdf(self, halfPlaneDensity):
		"""
		given an unnormalized density in the x-y upper half plane, normalize it and set the data member
		"""
		
		self.p = halfPlaneDensity
		self.p/= 4*sum(self.p)*self.aperture
	
	
	def __call__(self, xp = arange(-3,4)):
		"""
		return the probability of an occlusion for a given particle location xp
		"""
		#apply the weighted mask to the pdf and sum, multiplying by 4 to account for symmetry
		return 4*squeeze(sum(sum(self.mask(xp)*self.p[...,newaxis], axis = 0), axis = 0))
	
	def mask(self, xp):
		"""
		true where a particle could occlude a second particle located at (xp, 0)
		"""
		#need 3 dims for this
		Xp = xp[newaxis,newaxis,:]
		
		#slope of the line  x =  f(y) defining the zone of occlusion
		m = (self.sensorDistance-Xp)/self.aperture
		
		#if the x coordinate is greater than f(y), then the grid point is within the zone of occlusion
		bMask = self.g[0][...,newaxis]>(Xp+m*self.g[1][...,newaxis])
		
		#weight for thickness of the cone of occlusion
		theArg = (  (self.g[0][...,newaxis]-Xp)/m  )**2-self.g[1][...,newaxis]**2
		theArg[theArg<0]=0
		coneThickness = sqrt(theArg)
		
		
		return bMask*coneThickness
	
def plumeParticleCount(aperture = 2, 
			massFlux = .5e-3, 
			meanVelocity = 150, 
			meanVolume = 4*pi*(50e-6)**3/3, 
			density = 6100,
			plumeWidth = 30e-3):
	"""
	number of particles we can expect in the plume at any one time
	aperture - lens diameter in fraction of plume widths
	massFlux - flow of particles in kg/s
	meanVelocity  - average particle velocity m/s
	meanVolume - average volume of a particle, m^3
	density - particle density in kg/m^3
	plumeWidth - width of the plume in m
	"""
	particleFlux = massFlux/(density*meanVolume)
	linearParticleDensity = particleFlux/meanVelocity
	samplingVolumeLength = plumeWidth*aperture
	return linearParticleDensity*samplingVolumeLength-1

class weightedOcclusionCalculator(occlusionCalculator):
	"""
	compute the expectation of the blockage
	"""
	def __init__(self, sensorDistance = 3., aperture = 1. , gridDim = 500):
		occlusionCalculator.__init__(self, sensorDistance, aperture , gridDim)
	def mask(self, xp):
		"""
		weight the mask by 1/conicSectionalArea so we can get the fraction of blockage just by multiplying the result 
		by the average particle area
		"""
		#need 3 dims for this
		Xp = xp[newaxis,newaxis,:]
		
		#slope of the line  x =  f(y) defining the zone of occlusion
		m = (self.sensorDistance-Xp)/self.aperture
		
		#if the x coordinate is greater than f(y), then the grid point is within the zone of occlusion
		bMask = self.g[0][...,newaxis]>=(Xp+m*self.g[1][...,newaxis])
		
		#weight for thickness of the cone of occlusion
		R2 = (  (self.g[0][...,newaxis]-Xp)/m  )**2
		theArg = R2-self.g[1][...,newaxis]**2
		theArg[theArg<0]=0
		coneThickness= sqrt(theArg)
		
		#fractional blockage for a 25 micron particle in a 30 mm plume
		fractionalBlockage = (5e-3)**2/R2
		fractionalBlockage[fractionalBlockage==Inf] = 0
		fractionalBlockage[fractionalBlockage>1] = 1
		return bMask*coneThickness*fractionalBlockage
	
	
class uniformOcclusionCalculator(occlusionCalculator):
	def __init__(self, sensorDistance = 3., aperture = 1.):
		occlusionCalculator.__init__(self, sensorDistance, aperture)
		p = double(self.r<=(self.plumeWidth/2))
		self.setPdf(self, p)

	
		
		
if __name__=='__main__':
	plotsOn =True
	
	xp = linspace(-3,3,100)
	N = 20
	aperture = linspace(.1, 3, N) #aperture in fraction plume width
	plumeWidth = 30e-6 #m
	
	
	close('all')
	figure(0,figsize = (20,30))
	p = []
	r =[]
	for i in arange(N):
		print 'beginning aperture '+str(i)
		oc = weightedOcclusionCalculator(aperture = aperture[i])
		pOcclusion = oc(xp)
		relativeOcclusion = pOcclusion/pOcclusion[round(xp.shape[0]/2)]
		
		apStr = str(double(round(100*aperture[i]))/100)
		oStr = 'Blockage'
		pStr = 'fractional blockage'
		p+=[pOcclusion]
		r+=[relativeOcclusion]
		
		
		if plotsOn:
			figure(num = 0,figsize = (20,30))
			subplot(2,1,1)
			plot(xp, pOcclusion, linewidth = 1.5, color = (double(i)/N,0,double(N-i)/N), label = apStr)
			subplot(2,1,2)
			plot(xp, relativeOcclusion, linewidth = 1.5, color = (double(i)/N,0,double(N-i)/N))
			
			
			figure(num = 1, figsize = (20,30))
			subplot(2,2,1)
			subplots_adjust( hspace = .3, wspace = .2, left = .075, right = .95)
			title(oStr+'  in a Normal Density Plume,\n aperture = '+ apStr + ' plume widths')
			plot(xp, pOcclusion)
			xlabel('lateral particle location (z-score)')
			ylabel(pStr)
			
			
			subplot(2,2,2)
			plot(xp, relativeOcclusion)
			xlabel('lateral particle location (z-score)')
			ylabel('relative ' +pStr)
			title(oStr+' Relative to Center of Plume')
			
			
			subplot(2,2,3)
			title('Plume Particle Density')
			contourf(oc.g[0], oc.g[1], oc.p)
			colorbar()
			axis([-3,3,0,3])
			
			
			xlabel('plume lateral z-score')
			ylabel('plume vertical z-score')
			
			subplot(2,2,4)
			title('Zones of Occlusion')
			m = double(sum(oc.mask(xp)>0,axis = 2))
			contourf(oc.g[0], oc.g[1], m, range(0,np.max(m), 10))
			axis([-3,3,0,3])
			colorbar()
			xlabel('plume lateral z-score')
			ylabel('plume vertical z-score')
			sf('simpleOcclusion_aperture='+apStr)
			close(1)
	
	p = array(p)
	r = array(r)
	if plotsOn:
		figure(0)
		subplots_adjust( hspace = .3, wspace = .2, left = .075, right = .95)
		subplot(2,1,1)
		xlabel('plume lateral z-score')
		ylabel(pStr)
		legend()
		title('comparison of '+oStr+' across a range of apertures\naperture expressed in fraction plume width')
		subplot(2,1,2)
		xlabel('plume lateral z-score')
		ylabel('relative ' +pStr)
		title(oStr+' Relative to Center of Plume')
		sf('simpleOcclusionApertureComparison')
		
		figure(figsize = (20,30))
		subplot(2,1,1)
		subplots_adjust( hspace = .3, wspace = .2, left = .075, right = .95)
		plot(aperture, p[:,0], linewidth = 2)
		xlabel('aperture, fraction plume width')
		ylabel(pStr)
		title(oStr +' at Rear of Plume')
		
		subplot(2,1,2)
		plot(aperture, r[:,0], linewidth = 2)
		xlabel('aperture, fraction plume width')
		ylabel('relative ' + pStr)
		title(oStr+' Relative to Center of Plume')
		
		sf('backsideOcclusionVsAperture')
	
	meanVolumeDiameter = 50e-6
	Nstar = plumeParticleCount(aperture, meanVolume = 4*pi*(meanVolumeDiameter)**3/3)
	
	figure(figsize = (20,30))
	subplot(2,1,1)
	subplots_adjust( hspace = .3, wspace = .2, left = .075, right = .95)
	plot(aperture, Nstar)
	xlabel('aperture, fraction plume width')
	ylabel('nr. expected particles')
	title('Particles In Sensing Volume\n (meanParticleVolume)$^{1/3}$ = '+str(meanVolumeDiameter*1e6)+' microns')
	
	subplot(2,1,2)
	plot(aperture, squeeze(p[:,0])*Nstar, label = 'back')
	plot(aperture, squeeze(p[:,0])*Nstar/squeeze(r[:,0]), label = 'center')
	xlabel('aperture, fraction plume width')
	ylabel('expected '+pStr)
	title('Expected '+oStr+'Over a Range of Apertures')
	legend()
	sf('expectedOcclusions')
	