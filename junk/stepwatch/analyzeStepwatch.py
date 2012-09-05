from __future__ import with_statement
from stepwatch import *
from activityMonitor import activityBlockParser,activityBlockHistogram,blockFeatureWriter
import time
import getopt
import sys


###############################################################################
###########default arguments and command line flags defined here###############
###############################################################################
shortOpts = '-f         -o                 -r                      -s           -S                   -z                  -T                 -w                       -d               -D             -P             -u                    -n                     -N          -t             -m             -i             -v         -h'.split()
longOpts = '--folder   --outputFolder --subThresholdRunLength --stepThreshold  --stepCountThreshold --zeroStepRunLength --sustainedHiThres --hoursWornThreshold --daysWornThreshold  --daysWornMax  --avgPeakSteps --durationUpperBound  --intensityUpperBound  --binCount  --histogram    --summary --individualOutput --verbose  --help'.split()
defaultDic = dict(zip(shortOpts,
        [os.path.curdir,     '',              '3',                 '10',               '0',             '90',             '100',                '10',               '3',                   '100',       '5',           100,                  100,                  3,           '',            '',           '',           '' ,      '']))
        
shortOptsStr = 'f:o:r:s:S:z:T:w:d:D:P:u:n:N:tmivh'
###############################################################################                        
                        
long2short = dict(zip(longOpts, shortOpts))                 


def argHelp():
    print "Welcome to the Pythonic Stepwatch Analyzer, PySA."
    print "You can set the following options, read each entry as '-shortFlag, --longFlag, defaultValue':"
    print "-f, --folder, . : this folder must contain only the xls files (which may be tab delimited or true xls format) for analysis"
    print "-o, --outputFolder, /inputFolder/out/: folder where the output will be written"    
    
    print "-r, --subThresholdRunLength, 3: if subThresholdRunLength or more consecutive minutes are observed with < stepThreshold steps, this can break a run of zero-step minutes"    
    print "-s, --stepThreshold, 10: must have > stepThreshold steps in one minute to count as nonzero steps for the purpose of determining run length"    
    print "-S, --stepCountThreshold, 0: must have > stepCountThreshold steps in one minute to count as nonzero steps for the purpose of computing statistics"    
    print "-z, --zeroStepRunLength, 90: must have zeroStepRunLength zero-step minutes in a row to count as notWorn"
    print "-T, --sustainedHiThres, 100: sustained high activity must have at least this many steps"
    print "-w, --hoursWornThreshold, 10: must have worn the device for at least hoursWornThreshold to count the day"
    print "-d, --daysWornThreshold, 3: must have worn the device for at least this many days to use the subject in pan-subject statistics"
    print "-D, --daysWornMax, 100: will accept no more than this many days"
    print "-P, --avgPeakSteps, 5: average of the top P steps in each day"
    
    print "-m, --summary: takes no argument. If present will a summary of the summary statistics"    
    print "-i, --individualOutput: takes no argument. If present will write individual output files"    
    print "-v, --verbose: takes no argument. If present, causes messages to be printed during operation"
    
    print "-h, --help: displays this message. But you already knew that because you're reading this."

if __name__ == '__main__':
#    try:    
#        opts, args = getopt.getopt(sys.argv[1:], 'f:o:r:s:z:w:ivh',
#                               'folder  outputFolder subThresholdRunLength stepThreshold zeroStepRunLength hoursWornThreshold individualOutput verbose  help'.split())
#    except:
#        print "argument fail. try 'python analyzeStepwatch.py -h' for help."
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
    SWR  = stepwatchReader()
    
##########################################################    
################initialize the validator##################
##########################################################
    
    #must have > stepThreshold steps in one minute to count as nonzero steps
    #for the purpose of determining run length
    stepThreshold =  int(optDic[long2short['--stepThreshold']])
    sustainedHiThres = int(optDic[long2short['--sustainedHiThres']])
    
    stepCountThreshold =  int(optDic[long2short['--stepCountThreshold']])
    nRank = int(optDic[long2short['--avgPeakSteps']])
    statsFunList, statsHeaders = funFac(stepCountThreshold, sustainedHiThres, nRank)
        
    #must have zeroStepRunLength zero-step minutes in a row to count as notWorn
    zeroStepRunLength =  int(optDic[long2short['--zeroStepRunLength']]) 
    
    #if subThresholdRunLength or more consecutive minutes are observed with 
    # <stepThreshold steps, this can break a run of zero-step minutes
    subThresholdRunLength =  int(optDic[long2short['--subThresholdRunLength']]) 
    
    #must have worn the device for at least hoursWornThreshold to count the day
    hoursWornThreshold =  double(optDic[long2short['--hoursWornThreshold']])

    
    daysWornThreshold =  double(optDic[long2short['--daysWornThreshold']])
    daysWornMax =  int(optDic[long2short['--daysWornMax']]) 
        
    SWV = stepwatchValidator(stepThreshold = stepThreshold, 
                             zeroStepRunLength = zeroStepRunLength,
                             subThresholdRunLength = subThresholdRunLength,
                             hoursWornThreshold = hoursWornThreshold)
                             
    #parameters for the histogram
    durationUpperBound = double(optDic[long2short['--durationUpperBound']])
    intensityUpperBound = double(optDic[long2short['--intensityUpperBound']])
    binCount = double(optDic[long2short['--binCount']])
##########################################################                             

    #initialize the analyzer and individual output writer
    SWA  = stepwatchAnalyzer(SWV = SWV, daysWornMax = daysWornMax)
    SWW  = individualStepwatchWriter()
                       
    SWFP = stepwatchFolderProcessor(SWR = SWR,
                                    SWA = SWA,
                                    SWW = SWW,
                                    verbose = verbose,
                                    statsFunList = statsFunList,
                                    statsHeaders = statsHeaders,)
    
    ########## ANALYZE! #################
    sPrint("beginning analysis...")
    start = time.time()    
    (STATS, WAS_WORN, WORN_DIC_LISTS, statsHeaders, XX, SID) = SWFP(folder,
                                                           outFolder = outFolder, 
                                                           individualOutput = individualOutputOn)
    sPrint("done. elapsed time is %d minutes"%((time.time()-start)/60.0))
    ######################################
    
    ABP = activityBlockParser(stepThreshold,subThresholdRunLength, daysWornMax = daysWornMax)
    BFW = blockFeatureWriter()
    F = ABP(XX, WAS_WORN, WORN_DIC_LISTS )
    
    
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
