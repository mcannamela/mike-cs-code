from numpy import *
def carrierGasParameters(flowrate = array(5)):
	"""
	for flowrate in slm, return the mean and std dev of particle ejection velocities
	"""
	#obtain from: Finke, 'the influence of injector geometry and carrier gas flow rate on spray pattern'.
	#proceedings of the united thermal spray conference, 1997.
	mu=flowrate*1.76  + 3.9 #Finke, figure 4
	s = flowrate*.416 + .46 #finke, figure 5
	return (mu, s)