from arduinoLog import htm2500Logger 
import sys
import getopt



defaultDict = dict(zip(['-c','-f', '-d', '-t', '-v'],
    ['9', 'htmLog.log', 'Inf', '1','']))

def argFail():
    print 'command not understood'
    print ('usage: python HTMLog.py  \
-c COMM_PORT_NUMBER  \
-f LOG_FILE_NAME  \
-d LOG_DURATION  \
-t LOG_INTERVAL')

    print 'can also use flag -v for verbose logging'
    print 'when passing Inf as LOG_DURATION (the default), \
press ctrl-c to stop logging'

    print 'default call is: '
    
    callstr = 'HTMLog'
    for k in defaultDict.keys():
        callstr+= ' '+k+' '+defaultDict[k]
        
    print callstr

if __name__=='__main__':
    
 
    try:
        opts, args = getopt.getopt(sys.argv[1:], 'c:f:d:t:v')
        
            
        D = dict(opts)
        
        if '-v' in D.keys():
            silent = False
        else:
            silent = True
            
            
        for k in D.keys():
            if k!='-v':
                defaultDict[k] = D[k]
                
        d = defaultDict
            
            
            
        logger = htm2500Logger(commPortNumber = int(d['-c']), 
                     logfileName = d['-f'], 
                     sampleInterval = float(d['-t']),
                     silent = silent)
                         
        logger(float(d['-d']))
    except:
        argFail()
        
    
    