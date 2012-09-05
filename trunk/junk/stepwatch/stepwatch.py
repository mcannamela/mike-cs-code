from __future__ import with_statement
import xlrd 
import xlwt
import os
from numpy import *
#from pylab import *
import pdb


#threshold for steps/min for activity to be considered high
stepCountThreshold = 0

#threshold for steps/min for activity to be considered high
hiThres = 100
#minimum length for a block of activity above hiThres to be considered sustained
minHiBlockLength = 10
#maximum allowed break in a period of high activity
maxHiBreakLength = 2

#number of ranked peak steps to average for peak steps per day
nRank = 5

def funFac(stepCountThreshold, hiThres = hiThres, nRank = nRank):
    def stddev(X):
        s = std(X*(X>stepCountThreshold), ddof = 1)
        return s
    
    def meanStepsPerDay(X, isOn):
        """
        mean number of steps taken per day
        """
        if X.size == 0:
            return 0
        else:
            return mean(sum(X*(X>stepCountThreshold), axis = 1))
        
    def stdStepsPerDay(X, isOn):
        """
        standard error of number of steps per day
        """
        if X.size == 0:
            return 0
        else:
            s = std(sum(X*(X>stepCountThreshold), axis = 1), ddof = 1)
            
            return s
       
    def meanMinWornPerDay(X, isOn):
        """
        mean minutes worn per day
        """
        if X.size == 0:
            return 0
        else:
            return mean(sum(isOn, axis = 1))
    def stdMinWornPerDay(X, isOn):
        """    
        standard error of minutes worn per day
        """
        if X.size == 0:
            return 0
        else:
            return std(sum(isOn, axis = 1), ddof = 1)
    
    
    def meanMinHiPerDay(X, isOn):
        """
        mean minutes of high activity per day, where high activity is 
        minutes of steps/min above a threshold
        """
        if X.size == 0:
            return 0
        else:    
            return mean(sum(X*(X>stepCountThreshold)>=hiThres, axis = 1))    
    def stdMinHiPerDay(X, isOn):
        if X.size == 0:
            return 0
        else:    
            return std(sum(X*(X>stepCountThreshold)>=hiThres, axis = 1), ddof = 1)
    
    
    
    def sustainedMinHiPerDay(X, maxBreak):
        """
        minutes of sustained perious of high activity, defined as 
        having steps/min above a threshold for a minimum number of consecutive minutes
        within such a block of high activity, there may be breaks of up to maxBreak minutes
        where the steps/min falls below the threshold
        """
        if X.size == 0:
            return zeros(X.shape[0])
        else:
            BV = blockValidator(minHiBlockLength,maxBreak)
            minutesHi = zeros(X.shape[0])
            for i,x in enumerate(X):
                (isBig, blocks, isSmall, breaks) = BV(x, hiThres)
                minutesHi[i] = sum(isBig)
                
            return minutesHi
        
    def meanSustainedMinHi(X, isOn):
        if X.size == 0:
            return 0
        else:      
            return mean(sustainedMinHiPerDay(X, 0))
    def stdSustainedMinHi(X, isOn):
        if X.size == 0:
            return 0
        else:
            return std(sustainedMinHiPerDay(X, 0), ddof = 1)
        
    def meanSustainedMinHi_withBreaks(X, isOn):
        if X.size == 0:
            return 0
        else:
            return mean(sustainedMinHiPerDay(X, maxHiBreakLength))
    def stdSustainedMinHi_withBreaks(X, isOn):
        if X.size == 0:
            return 0
        else:
            return std(sustainedMinHiPerDay(X, maxHiBreakLength), ddof = 1)
            
    def rankedPeakSteps(X, isOn):
        if X.size == 0:
            return 0
        else:
            return mean(array([mean(sort(x[isOn[i]])[-int(nRank):]) for i,x in enumerate(X)]))
                
                 

    defaultStatsFunList = [ meanStepsPerDay, stdStepsPerDay,
                        meanMinWornPerDay, stdMinWornPerDay, 
                       meanMinHiPerDay, stdMinHiPerDay,
                       meanSustainedMinHi, stdSustainedMinHi,
                       meanSustainedMinHi_withBreaks, stdSustainedMinHi_withBreaks, rankedPeakSteps]
                       
    defaultStatsHeaders = [ 'Mean Steps per Day', 'Std. Err in Steps per Day', 
                        'Mean Minutes Worn per Day', 'Std. Err in Minutes Worn', 
                       'Mean Minutes Above %d per Day'%hiThres, 'Std. Error in Minutes Above %d'%hiThres,
                       'Mean Minutes Sustained Above %d for %d min, no breaks'%(hiThres, minHiBlockLength),
                       'Std. Error in Minutes Sustained Above %d for %d min, no breaks'%(hiThres, minHiBlockLength),
                       'Mean Minutes Sustained Above %d for %d min, with breaks of at most %d min'%(hiThres, minHiBlockLength, maxHiBreakLength),
                       'Std. Error in Minutes Sustained Above %d for %d min, with breaks of at most %d min'%(hiThres, minHiBlockLength, maxHiBreakLength),
                        '%d avg. peak steps'%nRank
                       ]
    return defaultStatsFunList, defaultStatsHeaders
    
defaultStatsFunList, defaultStatsHeaders = funFac(stepCountThreshold)

    
class stepwatchReader(object):
    """
    object in charge of reading the pedometer data, either from an excel
    or from a tab delimited text file
    """
    def __init__(self, dateRow = 29, dataStartRow = 31):
        """
        dateRow - row in the sheet where we can read the date
        dataStartRow - row in the sheet where the data begins
        
        """
        self.dateRow = dateRow
        self.dataStartRow = 31
        
    def __call__(self, fname):
        """
        pass a path to the file you want to read
        first try it as an excel, then as tab delimited 
        if it fails
        """
        try:
            t, X, D = self.excelRead(fname)
        except xlrd.biffh.XLRDError:
            t, X, D = self.tabDelimitedRead(fname)
        
        tStr = [self.min2timeStr(T) for T in t]
        
        return (t, D, X.T,  tStr)
        
    def excelRead(self, fname):
        """
        use xlrd to get the data out of an excel sheet
        
        returns:
            T - nMin array of times
            X - nDays x nMin array of steps, by date and time
        
        """
        wb = xlrd.open_workbook(fname)
        sh = wb.sheet_by_index(0)
        
        #dates are listed in this row
        dateRow = sh.row_values(self.dateRow)[1:]
        
        #convert dates to list of tuples
        D = [xlrd.xldate_as_tuple(d, wb.datemode)[:3] for d in dateRow]
        
        #loop of data rows, accumulating steps and times
        X = []
        T = []
        for i in range(self.dataStartRow, sh.nrows):
            r = sh.row_values(i)
            h,m,s = xlrd.xldate_as_tuple(r[0], wb.datemode)[-3:]
            if h==0 and m==0:
                t = 24*60
            else:
                t = 60*h+m
            T+=[t]
            X+=[float64(r[1:])]

        return (float64(T), float64(X), D)
        
    def tabDelimitedRead(self, fname):
        """
        as excelRead but for a tab delimited format
        """
        with open(fname, 'rb') as f:
            R = []
            for L in f.readlines():
                R+=[L.replace('\r\n','').split('\t')]
                
        D = [self.dateStr2tup(d) for d in R[self.dateRow][1:]]

        X = []
        t = []
        for r in R[self.dataStartRow:]:
            t+= [self.timeStr2min(r[0])]
            X+=[float64(r[1:])]
        
        
        return (float64(t), float64(X), D)
    def dateStr2tup(self, d):
        """
        convert a date in string form to a tuple
        """
        strTup = d.split(r'/')
        if len(strTup)==1:
            strTup = d.split(r'-')
            assert len(strTup)==3, 'date string parsing failed, looked for "/" and "-" but string was: '+d
            
            sep = '-'
        else:
            sep = '/'
            
        if sep == '/':
            t = (int32(strTup[2]), int32(strTup[0]), int32(strTup[1]))
        elif sep == '-':
            t = (int32(strTup[0]), int32(strTup[1]), int32(strTup[2]))
        return t
        
    def timeStr2min(self, tstr):
        """
        convert a timestamp string to a time in minutes
        """
        T = tstr.split()
        try:
            h,m,s = tuple(int32(T[0].split(':')))
        except ValueError:
            pdb.set_trace()
        if T[1] == 'AM' and h==12 and m>0:
            h = 0
        elif T[1] == 'AM' and h==12 and m==0:
            h = 24
        elif T[1] =='PM' and h!=12:
            h+=12
            
        return float64(60*h+m)

    def t2hms(self,t):
        """
        convert a time in minutes to hrs, min, sec tuple
        """
        h = floor(floor(t/60))
        tm = t-60*h
        m = floor(tm)
        s = (tm-m)*60
        return (h,m,s)
        
    def min2timeStr(self, t):
        """
        convert a time in min to a timestamp string
        """
        h,m,s = self.t2hms(t)
        return '%.2d:%.2d:%.2d'%(h,m,s)
        

class blockValidator(object):
    """
    this class analyzes a timeseries to find "valid" blocks
    of signal, where the signal has been above a threshold for 
    a certain length of time. 
    
    criteria for a valid block:
        -signal must be high (above a threshold) for a minimum number of 
            consecutive samples
        -small breaks are allowed, according to some break parameter
        
    """
    def __init__(self, minBlockLength = 10 , maxBreakLength = 2):
        """
        sets the parameters that determine a valid block
        
        minBlockLength - length for blocks of super-threshold values
        maxBreakLength - allow this many sub-threshold values in a block of 
                        superthreshold signal
        """
        self.minBlockLength = 10
        self.maxBreakLength = 2
    
    def __call__(self, x, threshold = 0):
        
        #find where x is below the threshold
        isSmall = x<threshold        
        
        #these blocks are breaks including small ones
        rawBreaks = self.getBlocks(isSmall)
        
        #retain only breaks long enough to break a streak
        breaks = self.validateBlockLength(rawBreaks,isSmall, self.maxBreakLength)
        
        #start out with all points assumed above the threshold
        isBig = (zeros_like(x)==0)
        
        #set isBig to false wherever there is a break long enough to break a streak
        self.setBlocks(breaks, isBig, val = False)
        
        #get all blocks where isBig is above the threshold, ignoring small breaks
        rawBlocks = self.getBlocks(isBig)
        
        #retain only blocks longer than minBlockLength
        blocks = self.validateBlockLength(rawBlocks, isBig, self.minBlockLength)
        
        #re-build isBig using only the valid blocks        
        isBig = (zeros_like(x)>0)
        self.setBlocks(blocks, isBig, val = True)
        
        return (isBig, blocks, logical_not(isBig), breaks)
    
    def setBlocks(self, B, x, val = True):
        """
        sets values of an array x in-place, according to 
        the list of blocks B
        
        B - list of slices into the array x
        x - 1d boolean array
        """
        for b in B:
            x[b] = val  
        
    def validateBlockLength(self, B, x, validBlocklength):
        """
        filters a list of blocks (slices into x) by the length 
        of the block
        
        B - list of slices into the array x
        x - 1d boolean array
        validBlockLength - blocks must be >= this length to be valid
        
        returns:
            validBlocks - list of blocks that pass the validity criteria
        """
        validBlocks = []        
        for b in B:
            if len(x[b]) < validBlocklength:
                continue
            else:
                validBlocks+=[b]
                
        return validBlocks
           
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
        
    def blockedArray(self, blocks, n):
        """
        use a list of blocks to make a boolean array of length n
        
        blocks - list of slices where the array will be true
        n - length of the array
        
        returns:
            x - 1d array of length n that is true for all slices in blocks
        """
        x = zeros(n)
        for b in blocks:
            x[b]=1 
        return x    
        
class stepwatchValidator(object):
    """
    this class analyzes stepwatch data for a series of days, 
    and determines whether each day is a valid datapoint, i.e. whether the 
    monitor was worn that day
    """
    def __init__(self, stepThreshold = 10, 
                 zeroStepRunLength = 90,
                 subThresholdRunLength = 3,
                 hoursWornThreshold = 10
                 ):
        #must have > stepThreshold steps in one minute to count as nonzero steps
        self.stepThreshold = stepThreshold
        
        #must have zeroStepRunLength zero-step minutes in a row to count as notWorn
        self.zeroStepRunLength = zeroStepRunLength 
        
        #if subThresholdRunLength or more consecutive minutes are observed with 
        # <stepThreshold steps, this can break a run of zero-step minutes
        self.subThresholdRunLength = subThresholdRunLength 
        
        #must have worn the device for at least hoursWornThreshold to count the day
        self.hoursWornThreshold = hoursWornThreshold
        
        #number of hours, minutes in a day
        self.nHours = 24.0
        self.nMin = 24.0*60
        
        #convert hoursWornThreshold to minutes
        self.minWornThreshold = self.hoursWornThreshold*60
        
    def __call__(self, X):
        """
        determines whether the device was worn for each day of data
        
        X - nDays x nMin, arrays holding the stepwatch signal for each day
        
        returns:
            wasWorn - length nDays boolean array, true if the device was worn that day
            wornDictList - list of dictionaries with masks and block slices for
                            each day, as returned by wasWorn()
        """
        wornDicList = [self.wasWorn(x) for x in X]
        wasWorn = array([w['wasWorn'] for w in wornDicList])
        
        return (wasWorn, wornDicList)
        
    
    def wasWorn(self, x):
        """
        apply the rules to one day's worth of data and label each minute as worn (true)
        or notWorn (false)
        x - length nMin array holding the number of steps for each minute of the day
        
        returns:
            dictionary with the keys:
                isOn - 1d array of len(x), true where the signal meets criteria for being ON
                onBlocks - list of slices into x where the signal is ON
                isOff - 1d array of len(x), true where the signal is not ON
                offBlocks - list of slices into x where the signal is not ON
                wasWorn - boolean, true if enough ON time was accumulated
        """
        assert len(x) == self.nMin, 'length of array is different from number of minutes in a day'
        
        
        #true if x nonzero
        nz = (x>0)
        
        #get all nonzero blocks
        nzBlocks = self.getBlocks(nz)
        
        #filter out short or sub-threshold blocks
        nzB = self.validateNZBlocks(nzBlocks, x)
        
        #initialize estimate of whether the signal is on/off
        isOn = self.blockedArray(nzB,self.nMin )
        isOff = logical_not(isOn)
        
        #now get zero-blocks from our first estimate
        zBlocks = self.getBlocks(isOff)
        
        #filter out short zero blocks
        zB = self.validateZBlocks(zBlocks)
        
        #now we know for sure when the device is off, and therefore
        #also when it is on
        isOff = self.blockedArray(zB,self.nMin)
        isOn = logical_not(isOff)
        
        #retrieve valid blocks of on-signal and off-signal
        onBlocks = self.getBlocks(isOn)
        offBlocks = self.getBlocks(isOff)
        
        #can now check whether the device was worn that day
        wasWorn = sum(isOn)>=self.minWornThreshold
        
        return dict(zip('isOn onBlocks isOff offBlocks wasWorn'.split(), 
                        [isOn, onBlocks, isOff, offBlocks, wasWorn]))
    def validateNZBlocks(self, B, x):
        """
        make checks on a list of nonzero blocks of x, accumulating those
        blocks which are valid.
        valid blocks must be either super-threshold in value (stepThreshold)
            or longer than some maximum low activity run length (subThresholdRunLength)
        
        B - list of slices into x
        x - 1d array representing the signal we are thresholding
        
        returns:
            validBlocks - list of slices meeting validty criteria
        """
        validBlocks = []        
        for b in B:
            if all(x[b]<self.stepThreshold) and len(x[b])<self.subThresholdRunLength:
                continue
            else:
                validBlocks+=[b]
        return validBlocks
    
    def validateZBlocks(self, B):
        """
        validate zero blocks, which must be of a certain length (zeroStepRunLength)
        
        B - list of slices representing regions of low signal
        
        returns:
            validBlocks - list of slices meeting validty criteria
        """
        validBlocks = []   
        for b in B:
            if (b.stop-b.start)>=self.zeroStepRunLength:
                validBlocks+=[b]
                
        return validBlocks
        
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
        
    def blockedArray(self, blocks, n):
        x = zeros(n)
        for b in blocks:
            x[b]=1 
        return x
        
        
class stepwatchAnalyzer(object):
    """
    equipped with a validator, this class calls user supplied stats functions
    on the valid parts of the data
    """
    def __init__(self, SWV = stepwatchValidator(), daysWornMax = None):
        """
        data will be masked according to the output of the stepwatchValidator
        """
        self.SWV = SWV
        if daysWornMax!=None:
            daysWornMax = int(daysWornMax)
        self.sly = slice(0, daysWornMax)
    
    def __call__(self, X, statsFunList = defaultStatsFunList):
        """
        determine the validity of the data in X
        compute the statistics specified by the statsFunList
        
        X - nDays x nMin, arrays holding the stepwatch signal for each day 
        statsFunList - list of statistical functions to apply to X, must have the 
                       call signature (X_valid, isOn), where X_valid is the signal 
                       from valid days in X and isOn a mask specifying valid
                       points in X.
                       
        returns:
            wasWorn - 1d length nDays array, true for days of valid signal
            wornDictList - length nDays list of dictionaries holding masks
                            for the signals' ON and OFF time, see stepwatchValidator.wasWorn()
            
        """
        if statsFunList ==None:
            statsFunList = defaultStatsFunList
        self.wasWorn, self.wornDicList = self.SWV(X)
        
        self.isOn = array([w['isOn'] for w in self.wornDicList])[self.wasWorn]
        
        S = self.computeStats(X[self.wasWorn], statsFunList)
        
        return (array(self.wasWorn), [dict(wd) for wd in self.wornDicList], S)
        
    def computeStats(self, X, statsFunList):
        """
        compute statistics of X for all functions in statsFunList
        """
        S = zeros(len(statsFunList), dtype = 'object')
        for j, fn in enumerate(statsFunList):
            x = X[self.sly]
            on = self.isOn[self.sly]
            S[j] = fn(x, on)
        return S
        
class individualStepwatchWriter(object):
    """
    write detailed output of the analysis of a single
    stepwatch file
    """
    def __init__(self, fmt = 'excel'):
        """
        can specify the output format
        for now excel is your only choice
        """
        self.fmt = 'excel'
        
        
    def __call__(self, fname, 
                 wasWorn, 
                 t, D, X,
                 S = None,
                 statsHeaders = defaultStatsHeaders, 
                 wornDicList = None):
        """
        fname - path to write the output
        wasWorn - length nDays array, true for valid days of data
        t, D, X - output from stepwatchReader, they are the time in minutes (nMin)
                    list of data tuples (nDays), and signal (nDays x nMin)
        S - the statistics of the data as computed by stepwatchAnalyzer
        statsHeaders - the names of the statistics that will appear in the output
        wornDictList - list of dictionaries (nDays) that hold masks specifying when
                        the signal is on/off (see stepwatchValidator.wasWorn())
        """
                     
        self.wasWorn = wasWorn      
        
        self.wb = xlwt.Workbook()
        
        self.writeValidDataSheet(t,D,X)
            
        
        if wornDicList != None:
            self.writeDiagnosticSheet(t, D,X, wornDicList)
            self.writeBlockSheet(wornDicList)
            
        if S!=None:
            self.writeStatSheet(self.subjectID(fname), S, statsHeaders)
            
        self.wb.save(fname)
            
    def writeValidDataSheet(self, t,D,X):
        """
        write a sheet with just data from valid days
        """
        sh = self.wb.add_sheet('validData')
        
        headerStyle = xlwt.easyxf('font: bold true; alignment: horizontal center; borders: bottom medium')
        centerStyle = xlwt.easyxf('alignment: horizontal center')
        L = self.formDateLine(D)
        [sh.row(0).write(i, d, headerStyle) for i,d in enumerate(L)]
        
        for i in range(len(t)):
            [sh.row(i+1).write(j, x, centerStyle) for j,x in enumerate(self.formDataLine(t[i], X[:,i]))]
            
    def writeDiagnosticSheet(self, t, D, X, wornDicList):
        """
        write a sheet with all the data, plus diagnostics
        like total minutes of on time, whether the day was flagged as valid
        and a side-by-side of threshold flags by the minute
        """
        allWorn = (self.wasWorn==self.wasWorn)
        isOn = array([w['isOn'] for w in wornDicList])
        sh = self.wb.add_sheet('diagnostics')
        
        dateRow = 0
        totalMinutesRow = 2
        dayValidRow = 1
        columnHeaderRow = 3
        dataStartRow = 4
        
        boldStyle = xlwt.easyxf('font: bold true')
        centerStyle = xlwt.easyxf('alignment: horizontal center')
        mergedCenterStyle = xlwt.easyxf('alignment: horizontal center; borders: right medium, left medium')
        dataHeaderStyle = xlwt.easyxf('alignment: horizontal center; font: italic true; borders: bottom medium, left medium')
        flagHeaderStyle = xlwt.easyxf('alignment: horizontal center; font: italic true; borders: bottom medium, right medium')
        
        #write the dates in merged cells 2 wide
        L = self.formDateLine(D, wasWorn = allWorn)
        sh.row(dateRow).write(0, L[0],boldStyle)
        sh.write(columnHeaderRow,0, 'time', 
                 xlwt.easyxf('font: bold true; alignment: horizontal center; borders: bottom medium'))
        for i in range(len(L)-1):
            sh.write_merge(dateRow,dateRow,2*i+1,2*i+2, L[i+1],mergedCenterStyle)
            sh.row(columnHeaderRow).write(2*i+1, 'data', dataHeaderStyle)
            sh.row(columnHeaderRow).write(2*i+2, 'wornFlag', flagHeaderStyle)
            
        dataCols = range(1, 2*X.shape[0], 2)
        flagCols = range(2, 2*(X.shape[0])+1, 2)

        #write total minutes worn that day
        sh.row(totalMinutesRow).write(0,'total minutes worn:', boldStyle)
        [sh.write_merge(totalMinutesRow,totalMinutesRow,c,c+1,sum(isOn[i]),mergedCenterStyle) for i,c in enumerate(dataCols)]
        
        #write whether the day was considered valid          
        sh.row(dayValidRow).write(0,'day valid:', boldStyle)
        [sh.write_merge(dayValidRow,dayValidRow,c,c+1, int(self.wasWorn[i]), 
                        mergedCenterStyle) 
                        for i,c in enumerate(dataCols)]
        
        dataStyle = xlwt.easyxf('alignment: horizontal center; borders: left medium')
        flagStyle = xlwt.easyxf('alignment: horizontal center; borders: right medium')
        for i in range(len(t)):
            L = self.formDataLine(t[i], X[:,i], wasWorn = allWorn)
            sh.row(dataStartRow+i).write(0,L[0], )
            [sh.row(dataStartRow+i).write(c, L[j+1],
                         dataStyle,
                         ) for j,c in enumerate(dataCols)]
            [sh.row(dataStartRow+i).write(c, int(isOn[j, i]), 
                        flagStyle,
                        ) for j,c in enumerate(flagCols)]
    
    def writeStatSheet(self, subjectID, S, headers):
        """
        write a sheet with the summary statistics from the valid days of data
        """
        sh = self.wb.add_sheet('summaryStatistics')
        wrapStyle = xlwt.easyxf('alignment: wrap true, horizontal center; font: bold true; borders: bottom medium')
        centerStyle = xlwt.easyxf('alignment: horizontal center')
        sh.row(0).write(0, 'subjectID', wrapStyle )
        sh.row(1).write(0, subjectID, centerStyle)
        sh.row(0).write(1, 'days worn', wrapStyle)
        [sh.row(0).write(i+2,h, wrapStyle) for i,h in enumerate(headers)]
        [sh.row(1).write(i+2,s, centerStyle) for i,s in enumerate(S)]
        sh.row(1).write(1, sum(self.wasWorn), centerStyle)
        for i in range(len(headers)+2):
             sh.col(i).width = int(sh.col(i).width*1.5)
        
    
    def writeBlockSheet(self, wornDicList):
        pass

    def formDataLine(self, t, x, wasWorn = None):
        """
        macro to make a time in minutes t and an array
        of data x into a list palatable for xlwt, while also masking for 
        valid days
        """
        if wasWorn==None:
            wasWorn = self.wasWorn
        x_ = []
        for i,xx in enumerate(x):
            if wasWorn[i]:
                x_+=[xx]
        L = [self.tStr(t)]+[int(xx) for xx in x_]

        return L
        
    def formDateLine(self, D, wasWorn = None):
        if wasWorn==None:
            wasWorn = self.wasWorn
        D_ = []        
        for i,d in enumerate(D):
            if wasWorn[i]:
                D_+=[d]
        L = ['date:']+[self.dateStr(d) for d in D_]
        return L
        
    
    def tStr(self, t):
        h,m,s = self.t2hms(t)
        return '%.2d:%.2d:%.2d'%(h,m,s)
    
    def t2hms(self,t):
        h = floor(floor(t/60))
        tm = t-60*h
        m = floor(tm)
        s = (tm-m)*60
        return (h,m,s)
  
    def dateStr(self, d):
        return str(d[0])+'/'+str(d[1])+'/'+str(d[2])
    def subjectID(self, fname):
        return os.path.split(fname)[1][:7]
        
class stepwatchStatsWriter(object):
    """
    class to write statistics for all subjects in one sheet
    """
    def __init__(self, fname, headers):
        """
        initialize the worksheet
        
        fname - path to write the output
        headers - list of strings to put in the column headers
        """
        self.wb = xlwt.Workbook()
        self.sh = self.wb.add_sheet('summaryStatistics')
         
        wrapStyle = xlwt.easyxf('alignment: wrap true, horizontal center; font: bold true; borders: bottom medium')
         
        self.sh.row(0).write(0, 'subjectID',wrapStyle)
        self.sh.row(0).write(1, 'days worn',wrapStyle)
        [self.sh.row(0).write(i+2,h,wrapStyle) for i,h in enumerate(headers)]
        for i in range(len(headers)+2):
            self.sh.col(i).width = int(self.sh.col(i).width*1.5)
        self.cnt = 1
        self.fname = fname
         
    def __call__(self, subjectID, daysWorn, S):
        """
        write a line for the given subject
        
        subjectID - int identifying the subject
        daysWorn - # of valid days for this subject
        S - list or array of summary statistics
        """
        centerStyle = xlwt.easyxf('alignment: horizontal center')
        self.sh.row(self.cnt).write(0, subjectID, centerStyle)
        self.sh.row(self.cnt).write(1, int(daysWorn), centerStyle)
        [self.sh.row(self.cnt).write(i+2, s, centerStyle) for i,s in enumerate(S)]
        self.cnt+=1
        
    def finalize(self):
        self.wb.save(self.fname)
        
        
class stepwatchFolderProcessor(object):
    """
    case management class, controls looping over the data for many subjects
    and writing individual and population output
    """
    def __init__(self, SWR  = stepwatchReader(),
                       SWA  = stepwatchAnalyzer(),
                       SWW  = individualStepwatchWriter(),
                       verbose = False, 
                       statsFunList = defaultStatsFunList,
                       statsHeaders = defaultStatsHeaders
                       ):
        self.SWR = SWR
        self.SWA = SWA
        self.SWW = SWW
        self.verbose = verbose
        self.statsFunList = statsFunList
        self.statsHeaders = statsHeaders
        
        
    def __call__(self, folder, 
                       outFolder = None, 
                       individualOutput = False):
        """
        folder - path to the data files
        outFolder - can re-direct the output if desired
        statsHeaders - 
        
        returns:
            STATS - nSubjects x nStats array, the stats for each individual
            WAS_WORN - nSubjects array of: wasWorn - 1d length nDays array, true for days of valid signal
            WORN_DIC_LISTS - nSubjects array of: wornDictList - length nDays list of dictionaries holding masks
                            for the signals' ON and OFF time, see stepwatchValidator.wasWorn()
            statsHeaders - names of all the stats computed and stored in STATS 
            XX - nSubjects array of: nDays x nMinutes arrays of data
            SID - nSubjects array of: strings representing the subject ID
        """
        
        if outFolder ==None:
            outFolder = os.path.join(folder, 'out')
        
        if not os.path.isdir(outFolder):
            os.mkdir(outFolder)
        
        self.SWSW = stepwatchStatsWriter(self.summaryName(outFolder), 
                                         self.statsHeaders)

        fnames = os.walk(folder).next()[2]
        self.sPrint( "walking folder, there are %d files here"%len(fnames))
        cnt = 0
        STATS = zeros((len(fnames),len(self.statsHeaders)))
        WAS_WORN = zeros(len(fnames), dtype = 'object')
        WORN_DIC_LISTS = zeros(len(fnames), dtype = 'object')
        XX = zeros(len(fnames), dtype = 'object')
        SID = zeros(len(fnames), dtype = 'object')
        try:
            for fn in fnames:
                if mod(cnt, max([1,len(fnames)/100]))==0:
                    self.sPrint( 'completed %d '%(100*cnt/float64(len(fnames)))+r'%'+': %d of %d'%(cnt, len(fnames)))
                    
                if fn[-3:]!='xls' and fn[-3:]!='csv':
                    continue
                    
                fname = os.path.join(folder, fn)
                outName = os.path.join(outFolder, fn[:-3]+'out.xls')
                
                t, D, X, tStr = self.SWR(fname)
                
                
                wasWorn, wornDicList, S = self.SWA(X, self.statsFunList)
                
                if individualOutput:
                    self.SWW(outName, wasWorn, t, D, X, 
                             S = S, 
                             wornDicList = wornDicList,
                             statsHeaders = self.statsHeaders)
                
                self.SWSW(self.subjectID(fname), sum(wasWorn), S)         
                STATS[cnt] = S.copy()
                WAS_WORN[cnt] = wasWorn.copy()
                WORN_DIC_LISTS[cnt] = [w.copy() for w in wornDicList]
                XX[cnt] = X.copy()
                SID[cnt] = self.subjectID(fname)
                cnt+=1
            self.SWSW.finalize()
            
        except:
            print "fatal error while processing "+fname
            raise
        
        return (STATS[:cnt], WAS_WORN[:cnt], WORN_DIC_LISTS[:cnt], self.statsHeaders, XX[:cnt], SID[:cnt] )
    

        
    def subjectID(self, fname):
        return os.path.split(fname)[1][:7]
        
    def summaryName(self, outFolder):
        return os.path.join(outFolder, 'statsSummary.xls')
        
    def sPrint(self, msg):
        if self.verbose:
            print msg
        
class stepwatchSummarySummarizer(object):
    """
    even higher level than the summary statistics, 
    gives statistics of the statistics, as well
    as separating subjects into good and bad based on 
    whether they produced valid data
    """
    
    def __init__(self):
        self.wb = xlwt.Workbook()
        self.gsh = self.wb.add_sheet('good subjects')
        self.bsh = self.wb.add_sheet('bad subjects')
        
        self.wrapStyle = xlwt.easyxf('alignment: wrap true, horizontal center; font: bold true; borders: bottom medium')
        self.centerStyle = xlwt.easyxf('alignment: horizontal center')
        
    def __call__(self, summaryFile, minimumDaysWorn = 3):
        
        #obtain the summary statistics for all subjects
        sID, dw, x = self.load(summaryFile)
        
        #split good and bad indices by number of days worn
        gidx = dw>=minimumDaysWorn
        gList = [int(f) for f in flatnonzero(gidx)]
        
        bidx = logical_not(gidx)
        bList = [int(f) for f in flatnonzero(bidx)]
        
        #entries with enough days
        gx = x[gidx]
        gsID = array(sID)[gidx]
        gdw = dw[gidx]
        
        dicKeys = ['subjectID', 'daysWorn', 'stats', 'statsHeaders']
        
        goodDic = dict(zip(dicKeys, [gsID, gdw, gx, self.rawHeaders[2:]]))
        
        #entries with too few days
        bx = x[bidx]
        bsID = array(sID)[bidx]
        bdw = dw[bidx]
        
        badDic = dict(zip(dicKeys, [bsID, bdw, bx, self.rawHeaders[2:]]))
        
        #write to the good and bad summary sheets
        self.writeSheet(gsID, gdw, gx)
        self.writeSheet(bsID, bdw, bx, isGoodSheet = False)
        
        self.wb.save(self.outFile(summaryFile))
        
        return (goodDic, badDic)
        
        
    def load(self, summaryFile):
        wb = xlrd.open_workbook(summaryFile)
        sh = wb.sheet_by_index(0)
        self.rawHeaders = sh.row_values(0)
        Xraw = [sh.row_values(i) for i in range(1,sh.nrows)]
        
        assert self.rawHeaders[0] == "subjectID", "first column must be 'subjectID'! double check to make sure you are using a valid summary file."
        assert self.rawHeaders[1] == "days worn", "second column must be 'days worn'! double check to make sure you are using a valid summary file."
        
        X = []
        subjectID = []
        daysWorn = []
        for x in Xraw:
            subjectID +=[x[0]]
            daysWorn +=[int32(x[1])]  
            X += [float64(x[2:])]

            
        return (subjectID, array(daysWorn), array(X))
        
    def writeSheet(self, sID, dw, x, isGoodSheet = True):
        if isGoodSheet:
            sh = self.gsh
        else:
            sh = self.bsh
        
        #set up row placement for various things
        headerRow = 1
        dataRow = 2
        
        sIDHeaderRow = 4
        sIDStartRow = 5
        
        #simple declaration of how many subjects make up this sheet
        sh.write(0,0, 'Number of subjects in this sheet: ', self.wrapStyle)
        sh.write(0,1, len(sID),self.centerStyle)

        #write the column headers for the statistics
        [sh.write(headerRow,i,h,self.wrapStyle) for i,h in enumerate(self.headers())]

        
        #need to know what statistical functions to apply to the columns
        fnList = self.functions()[1]
        nF = len(fnList)
        
        #summarize number of days worn
        for i in range(3):
            fn = fnList[mod(i,nF)]
            sh.write(dataRow, i, fn(dw), self.centerStyle)
        
        #summarize other statistics
        for i in range(3, len(self.headers())):
            fn = fnList[mod(i,nF)]
            idx = int(floor(i/nF))-1
            sh.write(dataRow, i, fn(x[:,idx]),self.centerStyle )

                
            
        #list the subjects making up this sheet
        sh.write(sIDHeaderRow, 0, 'Subjects in this sheet',self.wrapStyle)
        for i,s in enumerate(sID):
            sh.write(sIDStartRow+i, 0, s,self.centerStyle )
        
        #widen the columns a bit
        for i in range(len(self.headers())):
             sh.col(i).width = int(sh.col(i).width*1.5)
            
        
    def headers(self):
        RH = self.rawHeaders[1:]
        
        H = []
        for rh in RH:
            for k in self.functions()[0]:
                H+=[k+'( '+rh+' )']
        
        return H
            
    def functions(self):
        h = 'mean median stdDev'.split()
        f = [mean, median, std]        
        return (h,f)
    
    def outFile(self, fname):
        return fname[:-3]+'summary.xls'
        
if __name__=="__main__":

    folderTestOn = True
    singleTestOn = False
    
    
    if folderTestOn:
        testFolder = os.path.join(os.getcwd(), 'testFolder')
        outFolder = os.path.join(testFolder, 'out')
        
        SWFP = stepwatchFolderProcessor( verbose = True)
        (STATS, WAS_WORN, WORN_DIC_LISTS, statsHeaders, XX, SID) = SWFP(testFolder, individualOutput = False)
        
        summaryFile = SWFP.summaryName(outFolder)
        SWSS = stepwatchSummarySummarizer()
        SWSS(summaryFile, minimumDaysWorn = 1)
        
    if singleTestOn:    
        #fname = 'MB00004_1 [x2].xls'
        fname = 'MB00005_1 [x2].xls'
        #fname = 'MB00004_1_tabbed.xls'
        #fname = 'MB00099_1 [x2].xls'
        
        outName = fname[:-3]+'out.xls'
        
        SWR = stepwatchReader()
        
        print 'reading file'
        t, D, X, tStr = SWR(fname)
        
        SWA = stepwatchAnalyzer()
        
        print 'analyzing'
        wasWorn, wornDicList, S = SWA(X)
        print 'done. will now write results.'
        
        SWW = individualStepwatchWriter()
        SWW(outName, wasWorn, t, D, X, S = S, wornDicList = wornDicList )
        #print 'complete. begin plotting.'
        
        #isOn = array([w['isOn'] for w in wornDicList])    
        #for i in range(X.shape[0]):
         #   figure()
          #  subplot(2,1,1)
           # plot(t, X[i],'.-')
           # #ylim(99.5, 120)
           # title('wasWorn = '+str(wasWorn[i]))
           # subplot(2,1,2)
           # plot(t, isOn[i])
        

    
    
    
    
