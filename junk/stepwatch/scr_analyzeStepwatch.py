from stepwatch import stepwatchFolderProcessor
import os
import time

if __name__=="__main__":
    runFolder = os.path.join(os.getcwd(), 'runFolder')
    SWFP = stepwatchFolderProcessor()
    print "beginning analysis..."
    start = time.time()    
    (STATS, WAS_WORN, WORN_DIC_LISTS, statsHeaders) = SWFP(runFolder,individualOutput = True)
    print "done. elapsed time is %d minutes"%((time.time()-start)/60.0)