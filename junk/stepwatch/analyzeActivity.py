from __future__ import with_statement
from activityMonitor import *
import time
import getopt
import sys

###############################################################################
###########default arguments and command line flags defined here###############
###############################################################################
shortOpts = '-f         -o                 -r                      -s               -z                      -w                   -d                   -D            -P             -L               -M                   -H               -u                    -n                     -N          -t             -m          -i                 -S             -v         -h'.split()
longOpts = '--folder   --outputFolder --subThresholdRunLength --activityThreshold  --zeroActivityRunLength --hoursWornThreshold --daysWornThreshold  --daysWornMax --avgPeakSteps --lowThreshold  --moderateThreshold  --highThreshold  --durationUpperBound  --intensityUpperBound  --binCount  --histogram    --summary  --individualOutput  --returnSteps --verbose  --help'.split()
defaultDic = dict(zip(shortOpts,
        [os.path.curdir,     '',              '3',                 '100',                 '60',                  '10',               '3',                '100',        '5',         '2019',         '6000',               'Inf',                100,                  3000,            3,           '',           '',            '',            '',          '' ,       '']))
        
shortOptsStr = 'f:o:r:s:z:w:d:D:L:M:H:u:n:N:tmiSvh'
###############################################################################                        
                        
long2short = dict(zip(longOpts, shortOpts))                 


def argHelp():
    print "Welcome to the Pythonic Activity Analyzer, PyAA."
    print "You can set the following options, read each entry as '-shortFlag, --longFlag, defaultValue':"
    print "-f, --folder, . : this folder must contain only the xls files (which may be tab delimited or true xls format) for analysis"
    print "-o, --outputFolder, /inputFolder/out/: folder where the output will be written"    
    
    print "-r, --subThresholdRunLength, 3: if subThresholdRunLength or more consecutive minutes are observed with < activityThreshold activity counts, this can break a run of zero-activity minutes"    
    print "-s, --activityThreshold, 100: must have > activityThreshold activity counts in one minute to count as nonzero activity counts for the purpose of determining run length"    
    print "-z, --zeroActivityRunLength, 60: must have zeroActivityRunLength zero-activity minutes in a row to count as notWorn"
    print "-w, --hoursWornThreshold, 10: must have worn the device for at least hoursWornThreshold to count the day"
    print "-d, --daysWornThreshold, 3: must have worn the device for at least this many days to use the subject in pan-subject statistics"
    print "-D, --daysWornMax, 100: will accext no more than this many days"
    print "-P, --avgPeakSteps, 5: average of the top P steps in each day"
    
    print "-L, --lowThreshold, 2019: activity counts between activityThreshold and lowThreshold are classified as low activity"
    print "-M, --moderateThreshold, 6000: activity counts between lowThreshold and moderateThreshold are classified as moderate activity"
    print "-H, --highThreshold, Inf: activity counts between moderateThreshold and highThreshold are classified as vigorous activity"
    
    print "-u, --durationUpperBound, 1200: upper bound on duration for purposes of specifying histogram bins"
    print "-n, --intensityUpperBound, 500: upper bound on intensity for purposes of specifying histogram bins"
    print "-N, --binCount, 3: number of bins to use in each dimension for the histogram"
    
    print "-t, --histogram: takes no argument. If present will write a histogram of the features for each subject"
    print "-m, --summary: takes no argument. If present will write a summary of the summary statistics"    
    print "-i, --individualOutput: takes no argument. If present will write individual output files"    
    print "-S, --returnSteps: takes no argument. If present will use step counts instead of activity counts"    
    print "-v, --verbose: takes no argument. If present, causes messages to be printed during operation"
    
    print "-h, --help: displays this message. But you already knew that because you're reading this."

if __name__ == '__main__':
#    try:    
#        opts, args = getopt.getopt(sys.argv[1:], 'f:o:r:s:z:w:ivh',
#                               'folder  outputFolder subThresholdRunLength activityThreshold zeroActivityRunLength hoursWornThreshold individualOutput verbose  help'.split())
#    except:
#        print "argument fail. try 'python analyzeActivity.py -h' for help."
#        pdb.set_trace()
#        raise ValueError, "check your flags and arguments."
        
    opts, args = getopt.getopt(sys.argv[1:], shortOptsStr,
                               [L.replace('--', '') for L in longOpts])
    
    #convert any long flags to short ones
    optDic = dict(opts)
    for k in optDic.keys():
        if k in longOpts:
            optDic[long2short[k]] = optDic[k]
    
    
    #turn on or off verbosity
    verbose = '-v' in optDic.keys()
        
    def sPrint(msg):
        if verbose:
            print msg
            
    #turn on or off summary of the summary 
    summarySummaryOn = '-m' in optDic.keys()
    
    #turn on or off individual output
    individualOutputOn = '-i' in optDic.keys()
    
    #turn on or off individual output
    returnSteps = '-i' in optDic.keys()
    
    #turn on or off default output folder
    nonDefaultOutFolder = '-o' in optDic.keys()
    
    #turn on or off histogram output
    histogramOn = '-t' in optDic.keys()
                
    #show help or run the script
    if '-h' in optDic.keys():
        argHelp()
    else:    
        for k in defaultDic.keys():
            if k not in optDic.keys() and k not in ['-i -o -v -h'.split()]:
                optDic[k] = defaultDic[k]
                
        folder = optDic['-f']
        
        #setup output folder
        hf, tf = os.path.split(folder)
        
        if nonDefaultOutFolder:
            outFolder = optDic['-o']
        else:
            outFolder = os.path.join(folder, 'out')            
            
        if not os.path.isdir(outFolder):
            os.mkdir(outFolder)
    
    #initialize the file reader
    SWR  = wheyReader(returnSteps = returnSteps)
    
##########################################################    
################initialize the validator##################
##########################################################
    
    #must have > activityThreshold activity counts in one minute to count as nonzero activity counts
    #for the purpose of determining run length
    activityThreshold =  int(optDic[long2short['--activityThreshold']])
    
    #must have zeroActivityRunLength zero-activity minutes in a row to count as notWorn
    zeroActivityRunLength =  int(optDic[long2short['--zeroActivityRunLength']]) 
    
    #if subThresholdRunLength or more consecutive minutes are observed with 
    # <activityThreshold activity counts, this can break a run of zero-activity minutes
    subThresholdRunLength =  int(optDic[long2short['--subThresholdRunLength']]) 
    
    #must have worn the device for at least hoursWornThreshold to count the day
    hoursWornThreshold =  double(optDic[long2short['--hoursWornThreshold']])

    #need daysWornThreshold days of wearing the monitor to count the subject
    daysWornThreshold =  double(optDic[long2short['--daysWornThreshold']])  
    daysWornMax =  int(optDic[long2short['--daysWornMax']])  
    
    #boundaries for counting activity as low, moderate, or vigorous. these are upper bounds on 
    #the catagories, with activityThreshold forming the lower bound of the low catagory
    lowThreshold = double(optDic[long2short['--lowThreshold']])  
    moderateThreshold = double(optDic[long2short['--moderateThreshold']])  
    highThreshold = double(optDic[long2short['--highThreshold']]) 
    
    nRank = int(optDic[long2short['--avgPeakSteps']])
    
    activityStatsFunList, activityStatsFunHeaders = afunfac(activityThreshold,
                                                            lowThreshold,
                                                            moderateThreshold,
                                                            highThreshold,
                                                            nRank)
#    lowFun = lambda x,y: meanHrLowPerDay(x, y, activityThreshold, lowThreshold)
#    modFun = lambda x,y: meanHrModPerDay(x, y, lowThreshold, moderateThreshold)
#    highFun = lambda x,y: meanHrVigorousPerDay(x, y, moderateThreshold, highThreshold)
#    
#    activityStatsFunList[-3:] = [lowFun, modFun, highFun]

    #parameters for the histogram
    durationUpperBound = double(optDic[long2short['--durationUpperBound']])
    intensityUpperBound = double(optDic[long2short['--intensityUpperBound']])
    binCount = double(optDic[long2short['--binCount']])
    
    
    
        
    SWV = stepwatchValidator(stepThreshold = activityThreshold, 
                             zeroStepRunLength = zeroActivityRunLength,
                             subThresholdRunLength = subThresholdRunLength,
                             hoursWornThreshold = hoursWornThreshold)
##########################################################                             

    #initialize the analyzer and individual output writer
    SWA  = stepwatchAnalyzer(SWV = SWV, daysWornMax = daysWornMax)
    SWW  = individualStepwatchWriter()
                       
    SWFP = activityFolderProcessor(SWR = SWR,
                                    SWA = SWA,
                                    SWW = SWW,
                                    verbose = verbose,
                                    statsFunList = activityStatsFunList,
                                    statsHeaders = activityStatsHeaders)
    
    ########## ANALYZE! #################
    sPrint("beginning analysis...")
    start = time.time()    
    (STATS, WAS_WORN, WORN_DIC_LISTS, statsHeaders, X, SID) = SWFP(folder,
                                                           outFolder = outFolder, 
                                                           individualOutput = individualOutputOn)
    sPrint("done. elapsed time is %d minutes"%((time.time()-start)/60.0))
    ######################################
    
#    ABA = activityBlockAnalyzer()
    ABP = activityBlockParser(activityThreshold,subThresholdRunLength, daysWornMax = daysWornMax)
    BFW = blockFeatureWriter()
#    F = ABA(X, WAS_WORN, WORN_DIC_LISTS )
    F = ABP(X, WAS_WORN, WORN_DIC_LISTS )
    
    
    if histogramOn:    
        H,B = activityBlockHistogram(F, array([[0,0],
                                        [durationUpperBound,intensityUpperBound]]),
                                        (binCount+1)*ones(2) ) 
    else:
        H,B = (None, None)
        
    BFW(os.path.join(outFolder, 'features.xls'), F,
            'duration intensity'.split(), SID, H, B)    
    
    if summarySummaryOn:
        summaryFile = SWFP.summaryName(outFolder)
        SWSS = stepwatchSummarySummarizer()
        goodDic, badDic = SWSS(summaryFile, minimumDaysWorn = daysWornThreshold)
