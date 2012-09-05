#study the distribution and behavior of the blackbody residuals
import helpers as hp
reload(hp)
from helpers import *

from numpy import *
import pdb
from psMCMC import mcFluxSentinel as FS
from scipy.ndimage.filters import uniform_filter1d as  mvavg
smooth = lambda x,n: mvavg(x,n, axis = 0)


#\\\\\\\\\\\\\\\\\\\\\get the LUT's  and blackbody curve from the sensor object
try:
	print fs.tLUT.shape
except:
	fs = FS()
	fs.loadTemperatureLUT('tLUT.pkl')
	# bbRbg is normalized RGB corresponding to temperature tbb
	bbRbg, tbb = fs.blackBodyRgb()
	
	#choose a magnitude for the rgb vectors
	rgbMag = 100
	bbRgb = (rgbMag*bbRbg)
	
	#cutoff temperature
	tCutoff = 500
	bbRgb = bbRgb[tbb>tCutoff]
	tbb = tbb[tbb>tCutoff]
	
	#plot the rgb curves
	c = ['r', 'g', 'b']
	figure(figsize = (100, 120))
	[plot(tbb,bbRgb[:,i], c[i]+'-', linewidth = 2) for i in range(3)]
	xlabel('temperature, K')
	ylabel('rgb')
	title('Variation of color with temperature')
	
	#derivative of the temperature
	figure(figsize = (100, 120))
	plot(tbb[:-1], smooth(double(diff(tbb))/pnorm(diff(bbRgb,axis = 0),1),3) , linewidth = 2)
	xlabel('temperature, K')
	ylabel('dT/ds, K per pixel')
	title('derivative of temperature wrt rgb arclength\n rgb magnitude = '+str(rgbMag))
	
	#eliminate duplicate rgb tuples
	uIdx = uniqueRowsIdx(int32(bbRgb))
	bbRgb = bbRgb[uIdx]
	tbb = tbb[uIdx]
	
	#rgb grid
	G = int32(mgrid[:255,:255,:255])
	#\\\\\\\\\\\\\\\\\\\\\

#size of the neighborhood
rgbRadius = 5

#size of rgb vectors / size of noise vectors
snr = rgbMag/sqrt(3.*rgbRadius**2)




#helpers to find neighbors and their locations
displacementFn = lambda rgb_: transpose(array([excise(G[i],rgbRadius,int32(rgb_))[1].ravel()-rgb_[i] for i in range(3)]))
dFn = lambda rgb_: sum(displacementFn(rgb_)**2,1)

#bound the number of neighbors with a cube
nCols = (2*rgbRadius)**3

#thin the temperatures for a cleaner plot
thinIdx = arange(0,tbb.shape[0],1)

#initialize arrays of interest: temperature, residual, temperature error, rgb displacement, and rgb direction
T = transpose(tile(tbb[thinIdx], (nCols,1)))
B = [transpose(tile(bbRgb[thinIdx,i], (nCols,1))) for i in range(3)]
tRes = zeros((T.shape[0],nCols ))-1
tErr = tRes.copy()
displacement = zeros((T.shape[0], nCols,3))-1
direction = displacement.copy()
curveDirection = displacement.copy()
cdNorm = tRes.copy()
D = G[0].copy()
for i in range(T.shape[0]):
	#retrieve the neighboring residuals and slicing object 
	sly, ttRes = excise(fs.tResLUT, rgbRadius, int32(bbRgb[i]))
	n = ttRes.ravel().shape[0]	
	#look up the residual and error in temperature
	tRes[i,:n] = ttRes.ravel()
	tErr[i,:n] = fs.tLUT[sly].ravel()-tbb[i]
	
	#compute displacement and direction of the neighbors
	displacement[i,:n, :]= expand_dims(displacementFn(bbRgb[i]),0)
	dirNorm = expand_dims(pnorm(displacement[i,:n,:],1),1); dirNorm[dirNorm==0]=Inf
	direction[i,:n,:] = double(displacement[i,:n,:])/dirNorm
	try:
		cDirNorm = pnorm((bbRgb[i+1]-bbRgb[i]))
		curveDirection[i,:n,:] = expand_dims((bbRgb[i+1]-bbRgb[i])/cDirNorm,0)
		cdNorm[i,:n] = cDirNorm
	except:
		cDirNorm = pnorm((bbRgb[i]-bbRgb[i-1]))
		curveDirection[i,:n,:] = expand_dims((bbRgb[i]-bbRgb[i-1])/cDirNorm,0)
	print i


#some macros to aid plotting, range, scale, quiver, scatter, figure
rg = lambda x_: double(array([min(x_), max(x_)]))
sc = lambda x_: double(x_)/diff(rg(x_))
ex = (-1,1,  -1,1,  -1,1)
quiv = lambda x_,y_,z_,u_,v_,w_,f_: ml.quiver3d(sc(x_),sc(y_),sc(z_),sc(u_),sc(v_),sc(w_), 
																				extent = ex, scale_factor = .05, figure = f_, mode ='cone')
scat = lambda x_,y_,z_,s_,f_: ml.points3d(sc(x_),sc(y_),sc(z_),sc(s_), extent = ex, scale_factor = .01, scale_mode = 'none', figure = f_)
fig = lambda nm_: ml.figure(name = nm_, size = (1024, 720))

#pick out only valid neighbors from all bb points
idx = find( (tRes>0).ravel())


t = T.ravel()[idx]
res = arccos(tRes.ravel()[idx])
err = tErr.ravel()[idx]
dsp = pnorm(array([displacement[:,:,i].ravel()[idx] for i in range(3)]),0)
dire = [direction[:,:,i].ravel()[idx] for i in range(3)]
cdire = [curveDirection[:,:,i].ravel()[idx] for i in range(3)]
dirScalar = 0*t.copy()
for i in range(3):
	dirScalar+= dire[i]*cdire[i] 


tsly = permutation(res.shape[0])[:ceil(res.shape[0]/75.)]
qsly = permutation(res.shape[0])[:ceil(res.shape[0]/100.)]

figure()
hist(res[tsly], 100)
xlabel('residual, radians')
ylabel('count')
title('histogram of residuals')


#scatter: residual - error - temperature - rgbDisplacement
labs = ['Residual', 'Error', 'Temperature', 'Displacement']
rgs = array( [rg(res), rg(err), rg(t)]).ravel()
f = [fig(labs[0]+labs[1]+labs[2]+labs[3]+' Scatterplot: snr = '+str(snr))]; fnr = 0
scat( (res[tsly]),		err[tsly],		t[tsly],		(dsp[tsly]),		f[fnr])
ml.axes(xlabel = labs[0], ylabel = labs[1], zlabel = labs[2],
		ranges =rgs, extent = ex,  figure = f[fnr])
ml.outline(extent = ex, figure = f[fnr])


labs = ['Residual', 'Error', 'Displacement', 'Temperature']
rgs = array( [rg(res), rg(err), rg(dsp)]).ravel()
f+=[fig(labs[0]+labs[1]+labs[2]+labs[3]+' Scatterplot: snr = '+str(snr))]; fnr+=1
scat( (res[tsly]),		err[tsly],		(dsp[tsly]),  t[tsly], 	f[fnr])
ml.axes(xlabel = labs[0], ylabel = labs[1], zlabel = labs[2],
		ranges =rgs, extent = ex,  figure = f[fnr])
ml.outline(extent = ex, figure = f[fnr])


labs = ['Residual', 'Error', 'Off-Tangent Angle', 'Temperature']
rgs = array( [rg(res), rg(err),  rg(180*arccos(dirScalar[tsly])/pi)]).ravel()
#scatter: residual - error - temperature - rgbDisplacement
f+=[fig(labs[0]+labs[1]+labs[2]+labs[3]+' Scatterplot: snr = '+str(snr))]; fnr+=1
scat( (res[tsly]),		err[tsly],		arccos(dirScalar[tsly]) ,t[tsly],		f[fnr])
ml.axes(xlabel = labs[0], ylabel = labs[1], zlabel = labs[2],
		ranges =rgs, extent = ex,  figure = f[fnr])	
ml.outline(extent = ex, figure = f[fnr])


#quiver: the blackbody curve
# f+=[fig('Blackbody Curve Tangent vectors: snr = '+str(snr))]; fnr+=1
# quiv(B[0].ravel()[idx], B[1].ravel()[idx], B[2].ravel()[idx], 
	# cdire[0]*cdNorm.ravel()[idx],cdire[1]*cdNorm.ravel()[idx],cdire[2]*cdNorm.ravel()[idx],f[fnr])
# ml.axes(xlabel = 'red', ylabel = 'green', zlabel = 'blue',
		# ranges =array( [rg(B[0].ravel()), rg(B[1].ravel()), rg(B[2].ravel())]).ravel(), extent = ex,  figure = f[fnr])	
# ml.outline(extent = ex, figure = f[fnr])
