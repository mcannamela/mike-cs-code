# -*- coding: utf-8 -*-
"""
Created on Fri Feb 25 01:16:51 2011

@author: wichtelwesen
"""

import pysces
import cProfile
import pstats 
import time
pcMod = pysces.model('pysces_test_moiety1') 
pcMod.mode_integrator = 'LSODA'
pcMod.sim_end = 5.0
pcMod.sim_points = 40

N = 10.

def nRun():
    for i in range(N):
        pcMod.Simulate()
        

#cProfile.run('nRun()', 'prof')
#p = pstats.Stats('prof')
#p.strip_dirs().sort_stats(-1).print_stats()

start = time.time()
nRun()
elapsed  = time.time()-start

print "elapsed time is %d s, \n time per call: %.2e"%(elapsed, elapsed/N)