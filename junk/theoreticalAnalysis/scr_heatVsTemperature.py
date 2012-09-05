from numpy import *
import scipy
scipy.pkgload()
import pdb
import mlabMacros as mlm
from plotMacros import *


class piecewiseLinear(object):
	"""
	a piecewise linear curve
	"""
	def __init__(s, xi = arange(3), mi = array([2,.5,1]), y0 = 0):
		"""
		initialize given 
		xi - changepoints for the slopes
		mi[i] - slope between xi[i] and x[i+1]
		y0 - value of y at the first x
		"""
		#set number of changepoints, changepoints in the domain, slopes
		s.N = xi.shape[0]
		s.xi = xi
		s.mi = mi
		
		#compute the changepoints in the range
		d = diff(s.xi)
		s.yi = cumsum(array([y0]+[d[i]*s.mi[i] for i in range(s.N-1)] ))
		
		s.xbins = hstack((-Inf, s.xi,Inf))
		s.ybins = hstack((-sign(s.mi[0])*Inf, s.yi,sign(s.mi[-1])*Inf))
	
	def __call__(s,x):
		"""
		evaluate the curve at x
		"""
		IDX  = [(x>=s.xbins[i])*(x<s.xbins[i+1])!=0 for i in range(s.N+1)]
		IDX[1] = (IDX[1]+IDX[0])!=0; IDX = IDX[1:]
		y = zeros_like(x)
		for i in range(s.N):
			y[IDX[i]] =s.yi[i]+s.mi[i]*(x[IDX[i]]-s.xi[i]) 
		
		return y
		
	def gaussMap(s,mu,sig, n = 500):
		"""
		map the normal distributions over the domain x described by mu, sig to distributions over y
		"""
		normDist = lambda mu_,sig_,x_: scipy.stats.norm(loc = mu_, scale = sig_).pdf(x_)
		suppx = (mu-3*sig, mu+3*sig)
		x = linspace(suppx[0],suppx[1],n)
		px = normDist(mu, sig, x)
		y = s.__call__(x) 
		IDX  = [(y>=s.ybins[i])*(y<s.ybins[i+1])!=0 for i in range(s.N+1)]
		IDX[1] = (IDX[1]+IDX[0])!=0; IDX = IDX[1:]
		py = zeros_like(y)
		
		for j in range(s.N):
			py[IDX[j]]= normDist(s.yi[j]+(mu-s.xi[j])*s.mi[j], sig*s.mi[j],y[IDX[j]])
		
		return (x,px,y,py)
	
	def histMap(s,inpPdf):
		"""
		map any input distribution over x to an output distribution over y using histograms 
		"""
		N = 100000
		nHist = 1000
		x = inpPdf.rvs(N)
		hi, bi = histogram(x, nHist, normed = False)
		y = s.__call__(x)
		ho, bo = histogram(y, nHist, normed = False)
		return (bi[:-1]+diff(bi)/2, hi, bo[:-1]+diff(bo)/2, ho)
	
if __name__ =="__main__":
	close("all")
	
	#------------uncomment toy problem to test----------------
	# PL = piecewiseLinear()
	# mu = .5
	# sig = .5
	# x = linspace(-1, 4,400)
	# y = PL(x)
	# plot(x,y)
	
	# x,px, y, py = PL.gaussMap(mu,sig)
	# bi, hi, bo, ho = PL.histMap(scipy.stats.norm(loc = mu, scale = sig))
	
	# figure(figsize = (20,30))
	# plot(x,y, label = 'y(x)')
	# plot(x, mean(y)*px/max(px), label = 'p(x) (scaled)')
	# plot(mean(x)*py/max(py), y, label = 'p(y) (scaled)')
	# plot(bi, mean(y)*hi/max(hi), label = 'input histogram')
	# plot(mean(x)*ho/max(ho), bo, label = 'output histogram')
	# legend()
	#------------------------------------------------------------------
	
	#temperature dependence of specific heat
	T = double(array([0,300,900,1500,2100,2700,  2925, 2951, 3026, 3050, 3326, 5000]))
	Cp = double(array([444, 444, 607, 676, 944, 1227, 1227, 9417, 9417, 1381, 1381, 1400]))
	Q0 = 0
	
	#create the mapping between heat and temperature
	TQ = piecewiseLinear(T,Cp,Q0)
	QT = piecewiseLinear(TQ.yi, 1./TQ.mi, T[0])
	
	
	
	
	#make a mesh of the parameters of the heat distribution
	mu = linspace(.5e6,4e6,40)
	sig = .6e6
	n = 500
	
	idx = 32
	
	#preallocate for binsIn, histIn, binsOut, histOut
	q = zeros((mu.shape[0], n))
	pq = zeros((mu.shape[0], n))
	t = zeros((mu.shape[0], n))
	pt = zeros((mu.shape[0], n))
	
	
	for i in range(mu.shape[0]):
		q[i,:],pq[i,:],t[i,:],pt[i,:] = QT.gaussMap(mu[i], sig, n)
	

	#check the mapping 
	fig()
	plot(QT.xi, QT.yi, 'k', linewidth = 2, label = 'mapping')
	xlab('heat input, J/kg')
	ylab('Temperature, K')
	tit('Mapping from Heat Input to Temperature')
	plot(q[idx], pq[idx]*.5*np.max(QT.yi)/np.max(pq[idx]), 'r', linewidth = 1.5, label = 'heat pdf')
	plot(pt[idx]*.5*np.max(QT.xi)/np.max(pt[idx]), t[idx], 'g', linewidth = 1.5, label = 'temperature pdf')
	legend()
	axisFontSize()
	sf('heat-temp_distributionMap')
	
	with open('heat-temp.csv', 'w') as f:
		f.write('heat, pHeat, temp, pTemp\n')
		for i in range(t[idx].shape[0]):
			f.write(str(q[idx][i])+ ', ')
			f.write(str(pq[idx][i])+ ', ')
			f.write(str(t[idx][i])+ ', ')
			f.write(str(pt[idx][i])+ ', ')
			f.write('\n')
	
		
	
	f = mlm.fig('Transformation of Unimodal Heat Distribution to Bimodal Temperature Distribution')
	
	ml.mesh(mlm.sc(t), tile(expand_dims(mlm.sc(mu),1), (1, t.shape[1])), log(mlm.sc(pt)), extent = mlm.ex, figure = f)
	mlm.axe(['Temperature', 'mean heat input', 'p(T)'], mlm.rgs(t, mu, pt), f)
	mlm.out()
	
	
	
	