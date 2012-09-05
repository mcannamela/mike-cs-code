from __future__ import with_statement
from stepwatch import *

def afunfac(activityThreshold = 100, lowThreshold = 2019, moderateThreshold = 6000,
                  highThreshold = Inf, nRank = 5):
                      
    def meanStepsPerDay(X, isOn):
        """
        mean number of steps taken per day
        """
        if X.size == 0:
            return 0
        else:
            return mean(sum(X, axis = 1))
        
    def stdStepsPerDay(X, isOn):
        """
        standard error of number of steps per day
        """
        if X.size == 0:
            return 0
        else:
            
            return std(sum(X, axis = 1), ddof = 1)
            
    def meanHrWornPerDay(X, isOn):
        """
        mean minutes worn per day
        """
        if X.size == 0:
            return 0
        else:
            return mean(sum(isOn, axis = 1)/60.0)
    def stdHrWornPerDay(X, isOn):
        """    
        standard error of minutes worn per day
        """
        if X.size == 0:
            return 0
        else:
            return std(sum(isOn, axis = 1)/60.0, ddof = 1)
    
    def meanHrLowPerDay(X, isOn, lb = activityThreshold, ub = lowThreshold):
        """
        mean minutes worn per day
        """
        if X.size == 0:
            return 0
        else:
            h = 0
            for i,x in enumerate(X):
                h+= sum((x[isOn[i]]>=lb)*(x[isOn[i]]<ub))/60.0
            return h/X.shape[0]
    def meanHrModPerDay(X, isOn, lb = lowThreshold, ub = moderateThreshold):
        """
        mean minutes worn per day
        """
        if X.size == 0:
            return 0
        else:
            
            h = 0
            for i,x in enumerate(X):
                h+= sum((x[isOn[i]]>=lb)*(x[isOn[i]]<ub))/60.0
            return h/X.shape[0]
    def meanHrVigorousPerDay(X, isOn, lb = moderateThreshold, ub = highThreshold):
        """
        mean minutes worn per day
        """
        if X.size == 0:
            return 0
        else:
            
            h = 0
            for i,x in enumerate(X):
                h+= sum((x[isOn[i]]>=lb)*(x[isOn[i]]<ub))/60.0
            return h/X.shape[0]
            
    def rankedPeakSteps(X, isOn):
        if X.size == 0:
            return 0
        else:
            return mean(array([mean(sort(x[isOn[i]])[-int(nRank):]) for i,x in enumerate(X)]))
    
    activityStatsFunList = [ meanStepsPerDay, stdStepsPerDay,
                            meanHrWornPerDay, stdHrWornPerDay, 
                            meanHrLowPerDay, meanHrModPerDay,
                            meanHrVigorousPerDay, rankedPeakSteps
                           ]
                           
    activityStatsHeaders = [ 'Mean Activity Counts per Day', 'Std. Err in AC per Day', 
                            'Mean Hours Worn per Day', 'Std. Err in Hours Worn', 
                            'Mean Low Hours per Day', 'Mean Moderate Hours per Day', 
                            'Mean Vigorous Hours per Day',  '%d avg. peak activity'%nRank
                           ]     
    return activityStatsFunList, activityStatsHeaders

activityStatsFunList, activityStatsHeaders = afunfac()


class wheyReader(stepwatchReader):
    """
    read in data from the whey study
    """
    def __init__(self, dataStartRow = 11, columnHeaders ='date time activity steps'.split(), returnSteps=False ):
        """
        returnSteps - normally we are interested in activity counts. throw this flag true to switch to step counts
        """
        self.columnHeaders = columnHeaders
        self.dataStartRow = dataStartRow
        self.returnSteps = returnSteps
        
    def __call__(self, fname):
        """
        pass a path to the file you want to read
        first try it as an excel, then as tab delimited 
        if it fails
        
        returns:
            tt - nMin array of times
            uD - nDays list of date tuples
            activities - nDays x nTime array of activity counts, by date and time
            steps - nDays x nTime array of steps, by date and time
            tStr - list of times in string form
        """
        
        #pull the data in one long array spanning multiple days
        t, X, D = self.tabDelimitedRead(fname, ',')
        
        #find unique days in the data
        uD, n = self.uniqueDates(D)
        
        #need to make the long list into 2d with nDays x nMin
        tt = arange(int32(amin(t)), 1+int32(amax(t)))

        #had better have every minute of the day!
        assert len(tt)==60*24, "bad times, check the file: "+fname
        
        #time in string form
        tStr = [self.min2timeStr(T) for T in tt]
        
        #preallocate for steps and activity counts        
        steps = zeros((len(uD), len(tt)))
        activities = zeros((len(uD), len(tt)))
        
        for i,x in enumerate(X[0]):
            activities[n[i], argmin(abs(t[i]-tt))] = x
        for i,x in enumerate(X[1]):
            steps[n[i], argmin(abs(t[i]-tt))] = x
        
        if self.returnSteps:
            return (tt, uD, steps, tStr)
        else:
            return (tt, uD, activities, tStr)
    
    def excelRead(self, fname):
        pass
    
    def tabDelimitedRead(self, fname, splitchar = '\t'):
        with open(fname, 'rb') as f:
            R = []
            rawLines = f.readlines()
            if len(rawLines)==1:
                LINES = rawLines[0].split('\r')
            else:
                LINES = rawLines
            
            for L in LINES:
                R+=[L.replace('\r','').replace('\n','').split(splitchar)]
              
        D = []
        X = []
        t = []
        for i,r in enumerate(R[self.dataStartRow:]):
            
            try:
                assert len(r) == 4, 'bad row at %d, skipping!'%i
            except AssertionError:
                print "bad row at %d, length must be 4, skipping!"%i
                continue
            
            
            D+= [self.dateStr2tup(r[0])]
            
            t+= [self.timeStr2min(r[1])]
            try:
                X+= [float64(r[2:])]
            except ValueError:
                X+= [zeros(2)]
        

        return (float64(t), float64(X).T, list(D))
        
    def uniqueDates(self, D):
        """
        D - list of dates in tuple form
        
        returns:
            uD - list of unique members of D
            n  - array of integer labels for every member of D 
        """
        uD = [D[0]]
        n = zeros(len(D))
        cnt = 0
        d_ = D[0]
        for i,d in enumerate(D):
            if d==d_:
                n[i] += cnt
            else:
                d_= d
                uD+=[d]
                cnt+=1
                n[i] += cnt
                
        return (uD, n)
        
    def timeStr2min(self, tstr):
        """
        convert a timestamp string to a time in minutes
        """
        h,m,s = tuple(int32(tstr.split(':')))
            
        return float64(60*h+m)
        
class activityFolderProcessor(stepwatchFolderProcessor):
    def __init__(self, SWR  = wheyReader(),
                       SWA  = stepwatchAnalyzer(),
                       SWW  = individualStepwatchWriter(),
                       verbose = False,
                       statsFunList = activityStatsFunList,
                       statsHeaders = activityStatsHeaders
                       ):
        self.SWR = SWR
        self.SWA = SWA
        self.SWW = SWW
        self.verbose = verbose
        self.statsFunList = statsFunList
        self.statsHeaders = statsHeaders
    def subjectID(self,fname):
        return os.path.split(fname)[1][:-4]


# this class depreciated in favor of activityBlockParser
class activityBlockAnalyzer(object):
    def __init__(self, featureFunList = None):
        if featureFunList == None:
            self.featureFunList = [self.blockDuration, self.blockIntensity]
        else:
            self.featureFunList = featureFunList
    def __call__(self, XX, WAS_WORN, WORN_DIC_LISTS):
        """
        WAS_WORN - nSubjects array of: wasWorn - 1d length nDays array, true for days of valid signal
        WORN_DIC_LISTS - nSubjects array of: wornDictList - length nDays list of dictionaries holding masks
                            for the signals' ON and OFF time, see stepwatchValidator.wasWorn()
        XX - nSubjects array of: nDays x nMinutes arrays of data       
        
        returns:
            F - nSubjects array of: nFeatures x nBlocks array of features for 
                every block of activity performed by a subject on all valid days. 
        
        """
        
        nSubj = len(XX)
        F = zeros(nSubj, dtype = 'object')
        cnt = 0
        for i,X in enumerate(XX):
            ww = WAS_WORN[i]
            wdl = WORN_DIC_LISTS[i]
            f = []
            if not any(ww):
                f+= [zeros((1,len(self.featureFunList)))]
            else:
                for j in flatnonzero(ww):
                    for b in wdl[j]['onBlocks']:
                        x = X[j][b]
                        f+=[float64([fn(x) for fn in self.featureFunList])]
                
            F[cnt]= float64(f).T
            cnt+=1
        return F
    
    def blockDuration(self, x):
        return len(x)
        
    def blockIntensity(self, x):
        return mean(x)

class activityBlockParser(object):
    """
    class should parse out blocks of activity for valid data
    """
    def __init__(self, threshold, lowTime, featureFunList = None, daysWornMax = None):
        """
        threshold - activity counts above threshold will start the block
        lowTime - activity counts below the threshold for lowTime consecutive minutes will stop the block
        """
        self.thres = threshold
        self.lowTime = lowTime
        if featureFunList == None:
            self.featureFunList = [self.blockDuration, self.blockIntensity]
        else:
            self.featureFunList = featureFunList
        
        if daysWornMax!=None:
            daysWornMax = int(daysWornMax)
        self.sly = slice(0, daysWornMax)
        
    def __call__(self, XX, WAS_WORN, WORN_DIC_LISTS): 
        """
        WAS_WORN - nSubjects array of: wasWorn - 1d length nDays array, true for days of valid signal
        WORN_DIC_LISTS - nSubjects array of: wornDictList - length nDays list of dictionaries holding masks
                            for the signals' ON and OFF time, see stepwatchValidator.wasWorn()
        XX - nSubjects array of: nDays x nMinutes arrays of data       
        
        returns:
            F - nSubjects array of: nFeatures x nBlocks array of features for 
                every block of activity performed by a subject on all valid days. 
        
        """
        nSubj = len(XX)
        F = zeros(nSubj, dtype = 'object')
        cnt = 0
        for i,X in enumerate(XX):
            Y = X[self.sly]
            ww = WAS_WORN[i][self.sly]
            wdl = WORN_DIC_LISTS[i][self.sly]
            f = []
            if not any(ww):
                f+= [zeros((1,len(self.featureFunList)))]
            else:
                for j in flatnonzero(ww):
                    onBlocks = self.parse(Y[j])
                    for b in onBlocks:
                        x = Y[j][b]
                        f+=[float64([fn(x) for fn in self.featureFunList])]
                
            F[cnt]= float64(f).T
            cnt+=1
        return F
        
    def parse(self, X):
        low = X < self.thres 
        rawBreaks = self.getBlocks(low)

        breaks = []
        for b in rawBreaks:
            if len(X[b])>=self.lowTime or amin(X[b])<1:
                breaks+=[b]
        
        high = ones_like(X)==1        
        for b in breaks:
            high[b] = False
            
        activeBlocks = self.getBlocks(high)
        return activeBlocks

    def blockDuration(self, x):
        return len(x)
        
    def blockIntensity(self, x):
        return mean(x)
        
    def getBlocks(self, x):
        """
        find the blocks of x where x is true, and return them in a list of slices
        """
        x = r_['0,1', atleast_1d(0),x,atleast_1d(0)]        
        d = diff(x)

        start = flatnonzero(d == 1)
        stop = flatnonzero(d == -1)
        
        stop[stop>len(x)] = len(x)
        
        blocks = [slice(start[i], stop[i]) for i in range(len(start))]
        
        return blocks


def activityBlockHistogram(F, fb, N):
    """
    F - nSubjects array of: nFeatures x nBlocks array of features for 
                every block of activity performed by a subject on all valid days.
                
    fb - 2 x nFeatures array of lower bounds fb[0] and upper bounds fb[1] 
        for making the histograms of the features in F
    
    N - nFeatures array of number of bins for each feature
    
    returns: 
        H - nSubjects array of: multidimensional histograms of dimension according to N
        B - length nFeatures list of arrays of bin centers for each dimension of the histogram
    """
    
    B = [linspace(fb[0][i], fb[1][i], n) for i,n in enumerate(N)]
    H = zeros(len(F), dtype = 'object')
    for i,f in enumerate(F):
        H[i], blah = histogramdd(f.T, B)
        H[i]/=sum(H[i])
        
    return (H,[b[:-1] + .5*diff(b) for b in B])

class blockFeatureWriter(object):
    def __init__(self):
        self.wb = xlwt.Workbook()
        self.SH = [self.wb.add_sheet('features0')]
        self.wrapStyle = xlwt.easyxf('alignment: wrap true, horizontal center; font: bold true; borders: bottom medium')
        self.centerStyle = xlwt.easyxf('alignment: horizontal center')
        self.dataStartRow = 1
    def addSheet(self):
        self.SH+=[self.wb.add_sheet('features%d'%len(self.SH))]
        return self.SH[-1]
        
    def __call__(self, fname, F, featureHeaders, SID, H = None, B = None):
        """
        F - nSubjects array of: nFeatures x nBlocks array of features for 
                every block of activity performed by a subject on all valid days. 
        featureHeaders - nFeatures list of feature names
        SID - nSubjects array of subject ID strings
        
        H - nSubjects array of: multidimensional histograms of dimension according to N
        B - length nFeatures list of arrays of bin edges for each dimension of the histogram
        """
        sh = self.SH[0]
        self.writeHeader(featureHeaders, sh)
        self.featureHeaders = featureHeaders
        self.SID = SID
        rw = self.dataStartRow
        
        for i,f in enumerate(F):
            for k in range(f.shape[1]):      
                try:
                    sh.write(rw, 0, SID[i], self.centerStyle)
                    [sh.write(rw, j+1, double(f[j][k]), self.centerStyle) for j in range(f.shape[0])]
                except ValueError:
                    print "too many rows for one excel sheet "
                    print "subject id is "+SID[i]+", the %dth subject of %d"%(i, len(SID))
                rw+= 1
                if rw>=65536:
                    sh = self.addSheet()
                    self.writeHeader(featureHeaders, sh)
                    rw = self.dataStartRow
                    print "row overflow, adding features sheet"
                    
        if H != None:
            self.writeHistogram(H,B)
        
        self.wb.save(fname)
            
        
    def writeHeader(self, H, sh):
        sh.write(0,0,'subject ID',self.wrapStyle)
        for i,h in enumerate(H):
            sh.write(0,i+1, h, self.wrapStyle)
    
    def writeHistogram(self, H, B ):
        sh = self.wb.add_sheet('histogram')
        sh.write(0,0,'subject ID',self.wrapStyle)
        [sh.write(i,1,h,self.wrapStyle) for i,h in enumerate(self.featureHeaders)]
        
        G = mgrid[[slice(0, len(b)) for b in B]]
        
        for i,b in enumerate(B):
            g = G[i].flatten()
            [sh.write(i,j+2, b[g[j]], self.wrapStyle) for j in range(len(g))]
            
        for i, h in enumerate(H):
            x = h.flatten()            
            sh.write(i+len(B), 0, self.SID[i], self.centerStyle)
            [sh.write(i+len(B), j+2, double(x[j]),self.centerStyle)
                for j in range(len(x))]
            
      
        
    
 
if __name__=="__main__":
    folderTestOn = True
    singleTestOn = False
    
    
    if folderTestOn:
        testFolder = os.path.join(os.getcwd(), 'activityTestFolder')
        outFolder = os.path.join(testFolder, 'out')
        
        SWR = wheyReader(returnSteps=False)
        SWFP = activityFolderProcessor(SWR = SWR, verbose = True)
        STATS, WAS_WORN, WORN_DIC_LISTS, statsHeaders, X, SID = SWFP(testFolder, individualOutput = False )
        
        ABA = activityBlockAnalyzer()
        ABP = activityBlockParser(100,3)
        BFW = blockFeatureWriter()
#        F = ABA(X, WAS_WORN, WORN_DIC_LISTS )
        F = ABP(X, WAS_WORN, WORN_DIC_LISTS)
        H,B = activityBlockHistogram(F, array([[0,0],[1200,500]]),3*ones(2) ) 
        BFW(os.path.join(outFolder, 'featureTest.xls'), F,
            'duration intensity'.split(), SID, H, B)
        
        summaryFile = SWFP.summaryName(outFolder)
        SWSS = stepwatchSummarySummarizer()
        #SWSS(summaryFile, minimumDaysWorn = 1)
        
     #   figure()
    #    for i, f in enumerate(F):
    #        plot(f[0], f[1], '.', alpha = .5)
            
     #   xlabel('duration, min')
     #   ylabel(r'intensity, $\frac{counts}{min}$')
        
    if singleTestOn:    
        fname = '5.csv'
        outName = fname[:-3]+'out.xls'
        
        SWR = wheyReader()
        
        print 'reading file'
        t, D, activities, tStr = SWR(fname)
        
        SWV = stepwatchValidator(stepThreshold = 100, 
                 zeroStepRunLength = 60,
                 subThresholdRunLength = 3,
                 hoursWornThreshold = 10)
                 
        SWA = stepwatchAnalyzer(SWV)
        
        print 'analyzing'
        wasWorn, wornDicList, S = SWA(activities)
        print 'done. will now write results.'
        
        SWW = individualStepwatchWriter()
        SWW(outName, wasWorn, t, D, activities, S = S, wornDicList = wornDicList )
       # print 'complete. begin plotting.'
        
       # isOn = array([w['isOn'] for w in wornDicList])    
       # for i in range(activities.shape[0]):
       #     figure()
        #    subplot(2,1,1)
        #    plot(t, activities[i],'.-')
            #ylim(99.5, 120)
        #    title('wasWorn = '+str(wasWorn[i]))
        #    subplot(2,1,2)
        #    plot(t, isOn[i])
        #    ylim(-.1, 1.1)
            
   # show()
