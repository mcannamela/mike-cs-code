#define macros for plot formatting
from pylab import *
fsz = 26
xlab = lambda s_: xlabel(s_, fontsize = fsz)
ylab = lambda s_: ylabel(s_, fontsize = fsz)
tit = lambda s_: title(s_, fontsize = fsz)
sf = lambda figName_:  savefig(figName_+'.png', format = 'png', transparent = False)
def fig(n = None):
    if n==None:
        return figure(figsize = (14,10))
    else:
        return figure(n, figsize = (14,10))


def axisFontSize(ax = gca(), fsz = 22):
	for tick in ax.xaxis.get_major_ticks():
		tick.label1.set_fontsize(fsz)
	for tick in ax.yaxis.get_major_ticks():
		tick.label1.set_fontsize(fsz)

close('all')
