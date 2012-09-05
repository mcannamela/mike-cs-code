import os 
import scipy
from numpy import *
import numpy as np
import cPickle as pkl
from plotMacros import *
import mlabMacros as mlm
import trellesHelpers as TH
import qhullCommands as qh
pl = mlm.ml.pipeline
ml = mlm.ml
import pdb

N = 400
var = 'Th'
x,el,bx,bel = TH.fetchMesh(elements = True, boundaries = True)
t,y = TH.fetchVar(var, N = N)
U,s,Vh,yBar = TH.fetchSVD(var,N = N, compute = False)


bidx = int32([])
bels = int32([])
belFace = int32([])
blab = array([])

for i in range(len(bx)):
    bidx = r_[bidx,bx[i]]
    blab = r_[blab, ones_like(bx[i])*i]
    bels = r_[bels,bel[i][0]]
    belFace =  r_[belFace,bel[i][1]]


def bel2tri(el, bel, bface):
    #pdb.set_trace()
    f = array([  el[  bel[i]   ][  TH.faces[  bface[i]  ]  ]     for i in range(bel.shape[0])])
    tri = r_['0', f[:,:3], f[:,[0,2,3]]]
    return tri
    
tri = bel2tri(el, bels,belFace)
sly = slice(0,None,10)

#
xend = x[:-1,bx[-2]]*1e3
mend = Vh[0,bx[-2]]
otri,idx, yend = qh.delaunay(xend)

ft = mlm.fig('Outlet Mode 0')
ml.triangular_mesh(mlm.sc(xend[0]), mlm.sc(xend[1]), mlm.sc(mend), otri,extent = mlm.ex, figure = ft)
ml.axes(figure = ft, ranges = [np.min(xend[0]),np.max(xend[0]),
                                      np.min(xend[1]),np.max(xend[1]),
                                     -1,1], extent = mlm.ex, nb_labels = 5, 
                                     xlabel = 'x, mm', ylabel = 'y, mm', zlabel = '0th spatial mode, scaled')
ml.outline(figure = ft, extent = mlm.ex)

gt = mlm.fig('Outlet Solution at t = 0')
ml.triangular_mesh(mlm.sc(xend[0]), mlm.sc(xend[1]), mlm.sc(y[0,bx[-2]]), otri,extent = mlm.ex, figure =gt)
ml.axes(figure = gt, ranges = [np.min(xend[0]),np.max(xend[0]),
                                      np.min(xend[1]),np.max(xend[1]),
                                     np.min(y[0,bx[-2]]),np.max(y[0,bx[-2]])], extent = mlm.ex, nb_labels = 5, 
                                     xlabel = 'x, mm', ylabel = 'y, mm', zlabel = 'temperature, K')
ml.outline(figure = gt, extent = mlm.ex)

#

#f = mlm.fig('mode 0')
#trisrc = pl.triangular_mesh_source(x[0], x[1], x[2], tri[:tri.shape[0]/2])
#edge = pl.extract_edges(trisrc)
#msh = pl.surface(edge, figure = f)
#
#src = pl.scalar_scatter(x[0,sly], x[1,sly], x[2,sly],Vh[0,sly], figure = f)
#splat = pl.gaussian_splatter(src)
#scp = pl.scalar_cut_plane(splat)
#
#theFig = f
#execfile("mayascr_ROM_proposal.py")
##
#g = mlm.fig('Full Solution')
#
#gtrisrc = pl.triangular_mesh_source(x[0], x[1], x[2], tri[:tri.shape[0]/2], figure = g)
#gedge = pl.extract_edges(gtrisrc)
#gmsh = pl.surface(gedge, figure = g)
#
#gsrc = pl.scalar_scatter(x[0,sly], x[1,sly], x[2,sly],y[0,sly], figure = g)
#gsplat = pl.gaussian_splatter(gsrc)
#gscp = pl.scalar_cut_plane(gsplat)
##
#theFig = g
#execfile("mayascr_ROM_proposal.py")
ml.show()
