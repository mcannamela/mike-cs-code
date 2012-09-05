from scaleSpace import *
plotSummaryOn = True

computeOn = True
buildGraphOn = True
plotStrengthHistograms = False

nscales = 75
maxFeatureSize = .30
cutoff = 5.5
##::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
##::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
##::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

nKeep = 200
#imFile = 'blobsCropped.tif'
imFile = 'blobs.tif'
neighFun = neigh3d_124connected

fname = 'blobsRF.pkl'
if computeOn:
    im = float64(imread(imFile))
    
else:
    print "compute flag not thrown, loading pickle..."
    start = time.time()
    fname = 'blobsRF.pkl'
    with open(fname, 'rb') as f:
        rf = pkl.load(f)
        ss = rf.ss
        im = rf.ss.ss[0]
    print "done. %d"%(time.time()-start)
    
##::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
##::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
##::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
if computeOn:
    print "initializing scale space..."
    start = time.time()
    ss = scaleSpace2d(im, nscales = nscales, maxFeatureSize = maxFeatureSize)
    print "done. %d"%(time.time()-start)
    
    print "initializing blob finder..."
    start = time.time()
    rf = blobFinder2d(ss)
    print "done. %d"%(time.time()-start)
    
    print "finding blobs..."
    start = time.time()
    r = rf()
    print "done. %d"%(time.time()-start)
    
    print 'bouncing to pickle...'
    start = time.time()
    with open(fname, 'wb') as f:
        pkl.dump(rf, f, -1)
    print "done. %d"%(time.time()-start)

##::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
##::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
##::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
if buildGraphOn:
    #//////////// build the graph /////////////////
    print "building graph..."
    start = time.time()
    #g = mask2graph(any(rf.darkMask, axis = 0))()
    g = mask2graph(rf.darkMask)(neighFun)
    print "done. %d"%(time.time()-start)
    
    
    #//////////// find connected components/////////////////
    print "finding connected subgraphs..."
    start = time.time()
    G = array(networkx.connected_component_subgraphs(g), dtype = object)
    print "done. %d"%(time.time()-start)
    
    
    #//////////// number of nodes and ridge saliency/////////////////
    meanScale = array([mean(rf.ss.scales[array(gg.nodes())[:,0]]**.5) for gg in G])
    szIdx = flatnonzero(meanScale<cutoff)

    G = G[szIdx]
    meanScale = array([mean(rf.ss.scales[array(gg.nodes())[:,0]]**.5) for gg in G])
    L = array([len(gg.nodes()) for gg in G])
    #V = array([sum(rf.s[tuple([array(gg.nodes())[:,i] for i in range(3)])])**.5 for gg in G])
    V = array([mean(rf.s[tuple([array(gg.nodes())[:,i] for i in range(3)])])**.5 for gg in G])
    #V = array([sum((256.-rf.ss.ss[tuple([array(gg.nodes())[:,i] for i in range(3)])])/25.+ rf.s[tuple([array(gg.nodes())[:,i] for i in range(3)])])**.5 for gg in G])
    meanVal = array([mean(rf.ss.ss[tuple([array(gg.nodes())[:,i] for i in range(3)])]) for gg in G])
    minVal = array([amin(rf.ss.ss[tuple([array(gg.nodes())[:,i] for i in range(3)])]) for gg in G])
    
    

    szIdx = arange(V.shape[0])    
#    V = V[meanVal< median(meanVal[meanVal<median(meanVal)])]
#    V = V[meanVal< median(meanVal)]
#    V = V[minVal< median(minVal)]
    
    
    snip = V.shape[0]-nKeep
#    sidx = argsort(V+2*(amax(meanVal)-meanVal)/amax(meanVal) )[snip:]
    
    idx = argsort(V)[snip:]
    
    
    M = zeros_like(rf.brightMask[0,...])
    
    
    oFig = figure()
#        ax = subplot(1,2,1)
    ax = oFig.gca()
    imshow(im, cmap = cm.gray)
#        subplot(1,2,2)
#        imshow(any(rf.darkMask, axis = 0), cmap = cm.Accent)
#        for gg in G:
#            #x,y = (array(gg.nodes())[:,0], array(gg.nodes())[:,1])
#            x,y = (array(gg.nodes())[:,1], array(gg.nodes())[:,2])
#            plot(y,x, 'k.', ms = 4)
    
    
    cnt = 1.
    for gg in G[idx]:
        #x,y = (array(gg.nodes())[:,0], array(gg.nodes())[:,1])
        x,y = (array(gg.nodes())[:,1], array(gg.nodes())[:,2])
        t = array(gg.nodes())[:,0]
        M[x,y] = cnt
        
        
        
        if mean(rf.ss.scales[t]**.5)<cutoff:
            c = (cnt-1)/float64(idx.shape[0])
#            c = meanVal[idx[cnt-1]]/amax(meanVal[idx])
            ec = EllipseCollection(
                        2*rf.ss.scales[t]**.5,
                        2*rf.ss.scales[t]**.5,
                        0*x,
                        units='x',
                        offsets=r_['0,2', y,x].T,
                        transOffset=ax.transData,
                        facecolor = plt.cm.RdYlGn(c), 
                        edgecolor = plt.cm.RdYlGn(c), 
                        alpha = .35 )

            ax.add_collection(ec)
        
#            subplot(1,2,2)
#            plot(y,x, 'y.', ms = 4)
        
        
        cnt+=1.
    
#        
#        subplot(1,2,1)
    title('image with most salient \nconnected components overlaid\ncircle size - mean scale of connected components')
    #legend()
#        subplot(1,2,2)
#        title('all ridge pixels with retained connected components overlaid\n y - large, k - all')

    
    
    
if plotStrengthHistograms:
    figure()
    subplot(1,2,1)
    hist(log10(float64(L)), 100)
    xlabel('nr nodes')
    title('log10(connected component sizes)')
    
    subplot(1,2,2)
    plot(1- float64(cumsum(sort(L)))/sum(L))
    ylabel('1 - fraction total edge length')
    xlabel('component rank')
    ylim(0, 1.1)
    xlim(0, L.shape[0])

    figure()
    subplot(1,2,1)
    hist(log10(float64(V)), 100)
    xlabel('log10(total node strength)')
    title('connected component saliency')
    
    subplot(1,2,2)
    plot(1- cumsum(sort(V))/sum(V))
    xlabel('component rank')
    ylabel('1 - fraction total saliency')
    ylim(0, 1.1)
    xlim(0, V.shape[0])
##::::::::::::::::::::::::::::::::::::::::::::::::::::::::::



##::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
##::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
##::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
if plotSummaryOn:
    figure()
    subplot(1,2,1)
    imshow(im, cmap = cm.gray)
    title('the image')
    subplot(1,2,2)
    #imshow(sum(rf.ss.scales[:,newaxis, newaxis]*rf.darkMask, axis = 0 ))
    #imshow(log(1+sum(rf.darkMask, axis = 0 )), cmap = cm.Accent)
    imshow(any(rf.darkMask, axis = 0 ), cmap = cm.gray)
    title('dark ridge points')
##::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    



show()