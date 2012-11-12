# -*- coding: utf-8 -*-
"""
Created on Fri Nov 09 22:19:10 2012

@author: Michael
"""

class CodeJamOutputWriter(object):
    def __init__(self, fname):
        self.fname = fname
        open(self.fname,'w')
        
    def writeCase(self, caseNumber, caseOutput):
        with open(self.fname, 'a') as f:
            f.write(self.caseLine(caseNumber, caseOutput))
    
    def caseNumberStringPrefix(self, caseNumber):
        return "Case #%d: "%(caseNumber+1)
    def outputString(self,caseOutput):
        s = ""
        for c in caseOutput:
            s+=str(c)+" "
        return s[:-1]+"\n"
        
    def caseLine(self, caseNumber, caseOutput):
        return self.caseNumberStringPrefix(caseNumber)+self.outputString(caseOutput)
        
class CodeJamInputReader(object):
    def __init__(self,fname):
        self.fname = fname
        self.lineCounter = 1
        self.initLines()
        self.initSpecial()
        self.initCases()
    
    def getCases(self):
        return self.cases
        
    def initLines(self):
        with open(self.fname, 'r') as f:
            self.fLines = f.readlines()
    
    def initSpecial(self):
        pass
    
    def initCases(self):
        self.cases = []
        for i in range(self.nCases()):
            self.cases+=[self.readCase()]
    
    def nCases(self):
        return self.N
    
    def readCase(self):
        return None