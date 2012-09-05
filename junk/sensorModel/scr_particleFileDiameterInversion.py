from diameterInversionTools import *

#read in a tuttle particle file and break off the diameter estimates and blur numbers into an array
def partFile2Obs(fname):
    """
    grab the last six columns from a particle file, transpose and return them as an array
    """
    with open(fname, 'r') as f:
        X = []
        for L in f.readlines():
            try:
                X+=[float64(L.split('\t')[-6:])] #build up the array X into nParticles x 6 array
            except:
                continue
        return float64(X).T #return a 6 x nParticles array

def partFileAppend(oldfname,newfname, idx, d):
    """
    make a new particle file that includes the estimated diameter and z-location
    """
    cnt = 0
    goodCnt = 0
    with open(oldfname, 'r') as f:
        with open(newfname, 'w') as g:
            for L in f.readlines():
                if cnt>(idx.shape[0]-1):
                    break
                if idx[cnt]:
                    g.write(L.strip('\n')+'\t%.3e\t%.3e\n'%(d.Dref*exp(d.postDict['muD'][goodCnt]),d.postDict['muz'][goodCnt]))
                    goodCnt+=1
                else:
                    g.write(L.strip('\n')+'\tnan\tnan\n')
                cnt+=1
    

   
#pick the particle file you want here, make sure it's in the working directory
particleFile = 'Particles - 502.6.2'

#read in all the particle measurements
samples = partFile2Obs(particleFile)

#retain only those samples with small but nonzero blur numbers, and with estimated diameters less than 200
idx = all(logical_and(logical_and(samples[3:]<.3, samples[3:]>0.0), 6.6*samples[:3]<200), axis = 0)

#obs now holds the selected samples
obs = samples[:,idx]

#use N to control how many samples to run 
N = 11000


#refer to diameter inversion tools if you want to see how this object works
dic = diameterInversionCase(
                 inputTableName = 'D_B_lorenz_color_snap_extended.txt' , 
                 pklPrefix = 'particleRun_502.6.2',
                 measurementStdDevs = [.06, .08],#this is error we expect in log(measuredDiameter/Dref) and log(blurNr/Bref). your guess is as good as mine.
                 posteriorStdThreshold = .25,#only samples whose posterior standard deviations are less than this value will be plotted
                 Dref = 5., #reference values when transforming to log variables
                 Bref = .15,#reference values when transforming to log variables
                 N=N,
                 logarithmic = True,
                 forceCompute = False)

#retrieve the reference values for the log transform                 
refs = r_[ones(3)*dic.Dref, ones(3)*dic.Bref]

#transform the samples into logarithmic variables
obs = float64([log(obs[i]/refs[i]) for i in range(refs.shape[0])])

#calling the dic object performs the inversion
dic(samples = obs[:,:N], nD = 120, clean = True)

#make the summarizing figures
dic.plot()

#note we also discard samples with very low (think 0.0) likelihood, which can happen if measurementStdDevs is small
print 'discarded %d due to low likelihood, %d due to high variance, of %d'%(
            sum((dic.postDict['sigD']==0.)),
            sum((dic.postDict['sigD']>dic.posteriorStdThreshold)),
            N)


#plotting the distributions in micron (non-log) variables is more familiar
figure()
hist(samples[1,idx]*6.66, 50,label = 'raw estimate, green',facecolor = 'y',normed = True,  alpha = .35)
hist(exp(dic.postDict['muD'])*dic.Dref*6.66, 50,label = 'posterior mean estimate',facecolor = 'k',normed = True,  alpha = .35)
xlabel('diameter, $\mu$m')

##uncomment to plot the joint pdfs
#dic.pdfPlot()

##note that dic.postDict holds the posterior statistics you are interested in

#this will write you a new particle file with the estimates tacked on the end
partFileAppend(particleFile, particleFile+'_plus', idx[:N], dic)



