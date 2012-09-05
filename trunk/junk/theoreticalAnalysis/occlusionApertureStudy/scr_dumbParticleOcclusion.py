# assuming the plume is a uniformly dense circle of radius 1, estimate the fraction of the plume that can occlude a particle located at the vertical center of the plume for various 
#lateral locations given distance to the sensor and the vertical aperture of the sensor

if __name__=='__main__':
	R = 1
	sensorDistance = 5 #center of plume to sensor aperture
	verticalAperture = 3 #height of the lens

	linearParticleDensity = 1e6 #particles per meter
	sensingThickness = 10e-6 #width of a pixel


	x = linspace(-R,.999*R,100) #particle location

	occludedAreaFraction = (R-x)**2*arctan2(verticalAperture/2,sensorDistance-x,)/(pi*R**2)
	relativeOcclusion = occludedAreaFraction.copy()/occludedAreaFraction[round(x.shape[0]/2)]

	close('all')
	figure()
	subplot(2,1,1)
	plot(x, occludedAreaFraction)
	xlabel('lateral particle location')
	ylabel('fraction area occluding')
	title('approximate occlusion in a uniformly dense plume')
	
	subplot(2,1,2)
	plot(x, relativeOcclusion)
	xlabel('lateral particle location')
	ylabel('occlusion relative to center of plume')