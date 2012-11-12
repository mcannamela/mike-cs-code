# -*- coding: utf-8 -*-
"""
Created on Mon Nov 12 12:36:56 2012

@author: Michael
"""
from CodeJamIO import *
import re


class AlienLanguageInputReader(CodeJamInputReader):
    def initSpecial(self):
        self.linesPerCase = 1
        
        x = [int(i) for i in self.fLines[0].split()]
        self.L = x[0]
        self.D = x[1]
        self.N = x[2]

        self.lineCounter = 1+self.D
        
        self.words = []
        for i in range(1,1+self.D):
            word = self.fLines[i].replace('\n','')
            self.words+= [word]
            
    
    def getWordDict(self):
        return (self.L, dict(zip(self.words, [True for i in range(self.D)])))
        
    def readCase(self):
        caseLines = self.fLines[self.lineCounter:self.lineCounter+self.linesPerCase]
        self.lineCounter+=self.linesPerCase
        W = caseLines[0].replace('\n','')
        
#        assert len(W)==self.L, "W %d elements but got %d"%(self.L, len(W))
        return W
        
class AlienLanguageProblem(object):
    def __init__(self, L, wordDict):
        self.L = L
        self.D = wordDict
    
    def solve(self, W):
        self.reg = re.compile(fixupRegexp(W))
        cnt = 0
        for d in self.D.keys():
            
            if self.reg.match(d)!=None:
                cnt+=1
        return cnt
    
 
    
def fixupRegexp(W):
    w = ""
    ex = ''
    for c in W:
        if c=='(':
            ex='|'
            w+=c
            continue
        elif c==')':
            ex=''
        w+=c+ex
    return w.replace('|)',')')
    
if __name__=="__main__":
    reader = AlienLanguageInputReader("A-large-practice.in")
    cases = reader.getCases()
    print "cases:"
    for c in cases:
        print c
        
    P = AlienLanguageProblem(reader.getWordDict()[0],reader.getWordDict()[1] )
    
    soln = [P.solve(c) for c in cases]
    
    writer =  CodeJamOutputWriter("sample.out")
    print "\noutput"
    for i in range(len(cases)):
        writer.writeCase(i,[soln[i]])
        print soln[i]    
        
    if False:
        s = "(ab)c(de)"
        S = fixupRegexp(s)
        print S
        reg = re.compile(S)
        shouldMatch = ["acd", "bcd", "ace", "bce"]
        shouldNotMatch = "ccd add bde".split()
        
        matched = [reg.match(s)!=None for s in shouldMatch]
            
        notMatched = [reg.match(s)==None for s in shouldNotMatch]
            
            