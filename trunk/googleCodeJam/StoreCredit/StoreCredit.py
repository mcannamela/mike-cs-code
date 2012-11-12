# -*- coding: utf-8 -*-
"""
Created on Fri Nov 09 22:09:57 2012

@author: Michael
"""
from numpy import *
from CodeJamIO import *

class StoreCreditInputReader(CodeJamInputReader):
    def initSpecial(self):
        self.linesPerCase = 3
    def readCase(self):
        caseLines = self.fLines[self.lineCounter:self.lineCounter+self.linesPerCase]
        self.lineCounter+=self.linesPerCase
        C = int(caseLines[0])
        nItems = int(caseLines[1])
        L = int32(caseLines[2].split())
        assert len(L)==int(caseLines[1]), "L should have %d elements but got %d"%(nItems, len(L))
        return (C,L)
            
        

class StoreCreditProblem(object):
    def __init__(self, C, L):
        """
        C - integer value of store credit
        L - list or array of integer item values
        """
        self.C = C
        self.L = array(L)
        self.sortIdx = argsort(L)
        self.N = len(L)
        
        self.bBound = C-1
        self.aBound = self.L[0]
        self.bBoundIdx = self.N-1
        self.aBoundIdx = 0
        
    def solve(self):
        
        L = self.L[self.sortIdx]
        for i in arange(self.N-1, 0, -1):
            if L[i]>self.C-1:
                continue
            else:
                bIdx = i
                break
        
        cnt = 0
        aIdx = 0        
        while cnt<self.N**2:
            for i in arange(aIdx, bIdx):
                if (L[i]+L[bIdx])<self.C:
                    aIdx+=1
                if (L[i]+L[bIdx])==self.C:
                    aIdx = i
                    idx = [self.sortIdx[aIdx]+1,self.sortIdx[bIdx]+1]
                    idx.sort()
                    return idx 
                elif (L[i]+L[bIdx]>self.C):
                    bIdx-= 1
                    break
            cnt+=1
        return (-1,-1)

if __name__=="__main__":
#    reader = StoreCreditInputReader("StoreCreditSampleInput.txt")
    reader = StoreCreditInputReader("A-large-practice.in")
    cases = reader.getCases()
    print "cases:"
    for c in cases:
        print c
    
    soln = [StoreCreditProblem(c[0], c[1]).solve() for c in cases]
    
    writer =  CodeJamOutputWriter("StoreCreditOutput.txt")
    print "\noutput"
    for i in range(len(cases)):
        writer.writeCase(i,soln[i])
        print soln[i]
    
    
    
