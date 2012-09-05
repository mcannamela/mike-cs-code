from inverseParticleModeling import forwardMapper, particleGrid, brutalInverseMapper, inverseMapPostProcessor
import os, getopt, sys, time
from numpy import *
import pdb


###############################################################################
###########default arguments and command line flags defined here###############
###############################################################################
shortOpts = '-f -d -D -l -z -m -o -O -n -v -h'.split()
longOpts = '--folder   --D_min --D_max  --nD --z_offset --nZ --D_obs_min --D_obs_max --nD_obs --verbose  --help'.split()
defaultDic = dict(zip(shortOpts,
        [os.path.curdir, '.1',   '40',  '200', '.5',    '200',    '1',      '30',     '150',      '' ,      '']))
###############################################################################                        
                        
long2short = dict(zip(longOpts, shortOpts))                 


def argHelp():
    print "Welcome to the Pythonic LookUp Table Imputation System, pyLUTIS."
    print "You can set the following options, format is '-shortFlag, --longFlag, defaultValue':"
    print "-f, --folder, . : lookup table will be generated from this folder of folders of scans"
    print "-d, --D_min, .1: minimum true diameter to use in the forward map, in pixels"
    print "-D, --D_max, 40: maximum true diameter to use in the forward map, in pixels"
    print "-l, --nD, 200: number of meshpoints for true diameter in forward map"
    print "-z, --z_offset, .5: z, the optical axis coordinate, will be meshed from z_offset nearer than the red focal plane and further than the blue focal plane"
    print "-m, --nZ, 200: number of meshpoints for z in the forward map"
    print "-o, --D_obs_min, 1: minimum observed diameter to use in the inverse map, in pixels"
    print "-O, --D_obs_max, 30: maximum observed diameter to use inthe inverse map, in pixels"
    print "-n, --nD_obs, 150: "
    print "-v, --verbose: does not take an argument. If present, causes messages to be printed during operation"
    
    print "-h, --help: displays this message. But you already knew that because you're reading this."

if __name__ == '__main__':
    try:    
        opts, args = getopt.getopt(sys.argv[1:], 'f:d:D:l:z:m:o:O:n:vh',
                               'folder  D_min D_max nD z_offset nZ D_obs_min D_obs_max nD_obs verbose  help'.split())
    except:
        print "argument fail. try 'python scr_LUT_inversify.py -h' for help, which you clearly need."
        raise ValueError
    
    optDic = dict(opts)
    for k in optDic.keys():
        if k in longOpts:
            optDic[long2short[k]] = optDic[k]
    
    if '-v' in optDic.keys():
        verbose = True
    else:
        verbose = False
        
    def sPrint(msg):
        if verbose:
            print msg
               
    if '-h' in optDic.keys():
        argHelp()
    else:    
        for k in defaultDic.keys():
            if k not in optDic.keys() and k not in ['-v -h'.split()]:
                optDic[k] = defaultDic[k]
        
        lensFolder = optDic['-f']

        ########### DEBUG DEBUG DEBUG #######################
       # lensFolder = os.path.join(os.path.curdir, "oneScan") 
        #####################################################
        
        #setup output folder
        hf, tf = os.path.split(lensFolder)
        outFolder = os.path.join(hf, tf+'_out')
        if not os.path.isdir(outFolder):
            os.mkdir(outFolder)
    
        D_min = double(optDic[long2short['--D_min']])
        D_max = double(optDic[long2short['--D_max']])
        nD = int(optDic[long2short['--nD']])    
        
        z_offset = double(optDic[long2short['--z_offset']])  
        nZ = int(optDic[long2short['--nZ']])  
    
        D_obs_min = double(optDic[long2short['--D_obs_min']])
        D_obs_max = double(optDic[long2short['--D_obs_max']])
        nD_obs = int(optDic[long2short['--nD_obs']])
            
    
        sPrint("initializing the forward map...")
        
        FM = forwardMapper(lensFolder, verbose)
        zb = FM.zLim(offset = z_offset)
        PG = particleGrid(D_tup = (D_min, D_max, nD), Z_tup = (zb[0], zb[1], nZ))
    
        sPrint( "forward map initialized, begin the mapping...")
        #forwardMapPostProcessor spawned by FM
        PP = FM(PG)
        
        sPrint( "done. bouncing the forward map to disk")
        PP.writeForwardMap(os.path.join(os.path.curdir, "forwardMap"))
        sPrint( "initializing the inverse map...")
        BIM = brutalInverseMapper(PP()[0], PP()[1], PP()[2])
        X = BIM.X.copy()
        
        #use this slice to generate the observed diamter grid
        DSly = slice(D_obs_min, D_obs_max, nD_obs*1j)
        
        #grid of observed rgb diameters
        DOG = mgrid[[DSly for i in range(3)]]
        
        #preallocate for the (inverted) true diameters and error metric
        D = zeros(DOG[0].shape)
        E = zeros(DOG[0].shape)
        
        sPrint("done. forward mapping complete. inverse map initialized. begin inversion!")
        start = time.time()
        cnt = 0
        for i,x in enumerate(X):
            for idx, blah in ndenumerate(DOG[0]):
                D[idx], E[idx] = BIM(x, DOG[(slice(None, None),)+idx])
                cnt+=1
                if mod(cnt, (prod(DOG[0].shape)*X.shape[0])/100)==0:
                    sPrint('completed %d of %d in %d s'%(cnt, prod(DOG[0].shape)*X.shape[0], time.time()-start)) 
            if i==0:
                IMPP = inverseMapPostProcessor(X[i:i+1], DOG, D, E)
                IMPP.write(os.path.join(outFolder,'inverse'), withHeader = True)
            else:
                IMPP.reset(X[i:i+1], D, E)
                IMPP.write(os.path.join(outFolder,'inverse'), withHeader = False)
                
        elapsed = time.time()-start
        sPrint("completed %d inversions in %.2e s"%(cnt, elapsed))
        sPrint("that's %.2e s/inversion (including bounce time!)"%(float64(elapsed)/cnt))
        #raw_input('strike return to continue')
    