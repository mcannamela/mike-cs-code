import sys,pdb
from PyQt4 import QtCore, QtGui
from math import ceil
from timeWise_masterfulScheduling_ui import Ui_MainWindow
import os

try:
    import xlwt
except ImportError:
    xlwt = None
    
try:
    _fromUtf8 = QtCore.QString.fromUtf8
except AttributeError:
    _fromUtf8 = lambda s: s
qstr = QtCore.QString
    
conn = QtCore.QObject.connect
disconn = QtCore.QObject.disconnect
sig = QtCore.SIGNAL

frozenPolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Fixed,
                                              QtGui.QSizePolicy.Expanding)
unfrozenPolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Expanding,
                                              QtGui.QSizePolicy.Expanding)

teamBreakString = "  -----"
qTItem = QtGui.QTableWidgetItem


def find(L):
    """
    macro for finding the indices of the true elements of a list 
    """
    t = []
    for i,x in enumerate(L):
        if x:
            t.append(i)
    return t
def adjustTableColumns(tab, w):
    """
    pass a QTable tab and a width w to set the column widths
    """
    [tab.setColumnWidth(i,int(w)) for i in range(tab.columnCount())]

def it2float(it):
    return float(str(it.text()))
    
def it2int(it):
    return int(ceil(it2float(it)))
    
def dummyDollarsPerSTUFun():
    return 50000.0

class studentBody(QtGui.QWidget):
    def setup(self, totalEnrollment, 
                         selfContained, 
                         attendanceRate, 
                         dropoutsPerYear,
                         effectiveEnrollment,
                         nSections,
                         classSize,
                         updateButton):
        """
        totalEnrollment - list of spin boxes containing raw enrollment values 
                            per grade; last element is the total across grades
        
        selfContained - as totalEnrollment for self-contained student pop
        
        attendanceRate - list of spin  boxes for attendance rate (in %) by grade;
                            last element is the number-weighted mean across grades
        dropoutsPerYear - as totalEnrollment for the number of dropouts
        
        effectiveEnrollment - list of lcds displaying the effective enrollment by grade; 
                                last element is the total effective enrollment
            
        nSections - as effectiveEnrollment for the number of simultaneous sections
                    required based on the effective enrollment and class sizes
        
        classSize - as totalEnrollment for the class size; last element is the 
                    weighted average class size 
        
        updateButton - an invisible button used to signal an update to the student 
                        body
        """
        self.TE = totalEnrollment
        self.SC = selfContained
        self.AR = attendanceRate
        self.DPY = dropoutsPerYear
        self.EE = effectiveEnrollment
        self.nS = nSections
        self.CS = classSize
        self.idx = range(len(self.TE))
        self.updateButton = updateButton
        
    
    def getVal(self, wList):
        """
        macro for getting the values from a list of widgets
        """
        return [w.value() for w in wList]
    
    ###### getters ########
    def totalEnrollment(self):
        return self.getVal(self.TE)
    def selfContained(self):
        return self.getVal(self.SC)
    def attendanceRate(self):
        return self.getVal(self.AR)
    def dropoutsPerYear(self):
        return self.getVal(self.DPY)
    def classSize(self):
        return self.getVal(self.CS) 
    ########################
    
    def nSections(self):
        """
        compute the number of simultaneous sections we'll need 
        based on effective enrollment and class size
        """
        return [int(ceil(self.effectiveEnrollment()[i]/
                    max([1, float(self.classSize()[i])]))) for i in self.idx]
    def fullNSections(self):
        """
        compute the number of simultaneous sections we'll need
        based on the raw enrollment and class size
        """
        return [int(ceil(self.totalEnrollment()[i]/
                    max([1, float(self.classSize()[i])]))) for i in self.idx]

    def effectiveEnrollment(self):
        """
        compute the effective enrollment from the other enrollment values
        """
        return [round(float(self.totalEnrollment()[i]
                            -self.selfContained()[i]
                            -self.dropoutsPerYear()[i]/2.0)*self.attendanceRate()[i]/100.0)
                            for i in self.idx]
    def updateEffectiveEnrollment(self):
        """
        update the effective enrollment whenever some enrollment value changes
        """
        [w.display(int(self.effectiveEnrollment()[i])) for i,w in enumerate(self.EE)]
        self.updateTotals()
        self.updateButton.emit(sig('pressed()'))
        
        
    def updateNSections(self):
        """
        update the number of sections we need whenever it changes
        """
        [w.display(int(self.nSections()[i])) for i,w in enumerate(self.nS)]
        self.updateTotals()
        self.updateButton.emit(sig('pressed()'))
        
    def updateTotals(self):
        """
        update the total (across grade) enrollment values whenever a value changes
        """
        T = [0 for i in range(7)]
        A = [self.totalEnrollment, self.selfContained , self.attendanceRate , 
             self.dropoutsPerYear , self.effectiveEnrollment , self.nSections ,
             self.classSize ]
        arIdx = 2
        csIdx = 6
        for i in self.idx[:-1]:
            for j in range(len(T)):
                if j == arIdx:
                    T[j]+= self.totalEnrollment()[i]*A[j]()[i]
                elif j== csIdx:
                    T[j]+= self.nSections()[i]*A[j]()[i]
                else:
                    T[j]+=A[j]()[i]
        
        T[arIdx]/=float(T[0])
        T[csIdx]/=float(T[-2])
        
        X = [self.TE, self.SC , self.AR ,  self.DPY , 
                 self.EE , self.nS , self.CS ]
        for i,x in enumerate(X):
            try:
                x[-1].setValue(T[i])
            except AttributeError:
                x[-1].display(T[i])

    def write(self):
        """
        write the state of the object to an excel workbook
        """
        #if we don't have xlwt we will not write the output
        if xlwt==None:
            return 
            
        
            
        cHeaderStyle = xlwt.easyxf('alignment: wrap true, horizontal center; font: bold true; borders: bottom medium') 
        rHeaderStyle = xlwt.easyxf('alignment: wrap true, horizontal center; font: bold true; borders: right medium')     
        centerStyle = xlwt.easyxf('alignment: wrap true, horizontal center')        
        
        
        wb = self.parent().wb         
        shName = 'studentBody'
        try:
            shidx = wb._Workbook__worksheet_idx_from_name[shName.lower()]
            sh = wb.get_sheet(shidx)
            print 'sheet exists and was indexed'
        except:
            print 'sheet not found, creating it!'
            sh = wb.add_sheet(shName)
            sh.col(0).width = int(sh.col(0).width*1.5)

        sh._cell_overwrite_ok = True        
        
        
        
        rH = r'total enrollment; self-contained; dropouts/yr; attendance rate, %; nr. sections; class size'.split(';')
        cH = r'grade 9 10 11 12 total/avg.'.split()
        [sh.row(0).write(i, h, cHeaderStyle) for i,h in enumerate(cH)]
        [sh.write(i+1,0,h, rHeaderStyle) for i,h in enumerate(rH)]
        
        
        fnList = [self.__getattribute__(k) for k in 
                    'totalEnrollment selfContained dropoutsPerYear attendanceRate nSections classSize'.split()]
        for i,f in enumerate(fnList):
            for j,x in enumerate(f()):
                sh.write(i+1, j+1, x, centerStyle)
                
        p = os.path.join(os.getcwd(), self.parent().saveName()+'.xls')
        print 'writing excel to the path '+ p
        wb.save(p)
                
        
            
            
        
class subjectsParser(QtGui.QWidget):
    """
    used to get subjects out of the text boxes in the scheduling splitters, 
    so that they can be tablulated 
    """
    def setup(self, #subjectSplitters, 
                 subjectsTable, 
                 subjectsTableCredits,
                 periodLengthFun,
                 periodLengthCreditsFun):
        """
        subjectSplitters - list of splitters for each grade, each one holding 
                            nPeriods text boxes with the subject titles in them
        
        subjectsTable - the QTable object used to display the real hours spent
                        per day for each subject and grade

        subjectsTableCredits - as subjectsTable, but displays the time spent per
                                subject in credits

        periodLengthFun - function to tell us the length of each period in minutes

        periodLengthCreditsFun - function to tell us the length of each period in credits
        """
        
        
        #the splitters hold the text box widgets where the subjects 
        #will be typed
#        self.SP = subjectSplitters
#        print self.SP
        
        #store the text boxes in a list of lists, one list per grade
        self.TB = [[] for i in range(4)]
#        self.TB = [self.getWidgets(sp) for sp in self.SP]
#        print len(self.TB)
#        print [len(tb) for tb in self.TB]
        
        #the subjects table has a subject in each row and the time 
        #allocated to that subject (per grade) in each column
        self.ST = subjectsTable
        self.STC = subjectsTableCredits
        
        #the periodLengths calback returns the current length of 
        #each block (same across grades) so that we can figure out 
        #how much airtime each period is getting
        self.periodLengths = periodLengthFun
        self.periodLengthsCredits = periodLengthCreditsFun
        
        #need a list to hold any timers we create
        self.timerID = []
        
        #hook up the text box signals to our parsing function
#        for W in self.TB:
#            for w in W:
#                conn(w, sig("textChanged()"), self.textChanged)
    def addPeriod(self, TB):
        for i,tb in enumerate(TB):
            self.TB[i]+=[tb]
            conn(tb,sig("textChanged()"), self.textChanged)
    def removePeriod(self):
        for i in range(len(self.TB)):
            self.TB[i] = self.TB[i][:-1]
    def getWidgets(self, sp):
        """
        macro to return the widgets contained by a splitter
        """
        return [sp.widget(i) for i in range(sp.count())]
    def block(self, b):
        """
        we'll need to block signals from the text boxes sometimes
        """
        for W in self.TB:
            for w in W:
                w.blockSignals(b)
        
    def textChanged(self):
        """
        whenever the text changes in one of the boxes, we're just going
        to re-parse them all...
        to make this approach efficient we block any further signals that might
        call this function and set a timer to allow time for the user to finish 
        typing ~1 word
        """
        self.block(True)
        self.timerID += [self.startTimer(500)]
        
    
    def timerEvent(self, E):
        """
        here we perform the actual update from the textChanged() callback
        """
        #we need to kill all timers to avoid bigtime errors
        [self.killTimer(tid) for tid in self.timerID] 
        self.timerID = []
        
        #parse the text boxes for subjects and the time 
        #allotted to each
        S,P,PC = self.getSubjectsAndLengths()
                
        #find unique subjects and total the time allotted to each
        D,DC, Dtot, DCtot = self.reduceSubjects(S,P,PC)
                
        #guard against an empty table
        if D!=None:
            #if we have some subjects, we'll update the tables
            self.setTable(D,DC, Dtot, DCtot)
        else:
            self.blankTable()
            
        #unblock the text box signals so we can catch future updates
        self.block(False)
        
             
    def getSubjectsAndLengths(self):
        """
        parse through every text box, and make two lists of lists
        S - S[i] holds a list of subjects for grade i
        P - P[i] holds a list of hours allotted to the subjects in S[i]
        PC - as P but measured in credits instead of hours
        
        thus S[i][j] gives the jth subject of the ith grade
        """
        S = []
        P = []
        PC = []
        #for each grade's text boxes
        for j in range(len(self.TB)):
            S += [[]]
            P += [[]]
            PC += [[]]
            #for each period in each grade
            for i,tb in enumerate(self.TB[j]):
                
                #s is a list of the subjects in the box
                s = self.parseBox(tb)
                L = self.periodLengths()[i]
                LC = self.periodLengthsCredits()[i]
                
                #could have multiple subjects in one box, so we have to divide
                #the time of the block amongst those subjects for proper
                #accounting
                S[j]+= list(s) #must call list() to copy s! 
                if s==[]:
                    P[j] += []
                    PC[j]+= []
                else:
                    P[j]+= [L/float(len(s))/60.0 for i in range(len(s))]
                    PC[j]+= [LC/float(len(s)) for i in range(len(s))]
        
        return (S,P,PC)
                
    def reduceSubjects(self, S,P,PC):
        """
        compute a unique list of subjects and tally the hours/credits spent on each
        """
        #accumulate a master list of subjects
        allSubjects = []
        for SS in S:
            for s in SS:
                allSubjects+=[s]
        
        if allSubjects == []:
            return (None, None, None, None)
        #use set to get uniques keys
        uniqueSubjects = list(set(allSubjects))
        
        #make one dictionary for each grade, initialize the time to 0
        D = [dict(zip(uniqueSubjects, [0 for u in uniqueSubjects])) 
                for j in range(len(S))]
        DC = [dict(zip(uniqueSubjects, [0 for u in uniqueSubjects])) 
                for j in range(len(S))]
        
        
        #accumulate the time in each grade's dictionary
        #for each grade's subjects
        for j,ss in enumerate(S):
            #for each block of the day
            for i, s in enumerate(ss):
                D[j][s]+= P[j][i] 
                DC[j][s]+= PC[j][i] 
        
        #sum up the times over all grades
        Dtot = dict(zip(uniqueSubjects, [0 for u in uniqueSubjects]))
        DCtot = dict(zip(uniqueSubjects, [0 for u in uniqueSubjects]))
        for k in Dtot.keys():
            for d in D:
                Dtot[k]+=d[k]
            for d in DC:
                DCtot[k]+=d[k]
                
        return (D,DC, Dtot, DCtot)
                
    def setTable(self, D,DC,  Dtot, DCtot):
        """
        use the dictionaries to set the tables
        D - D[subjectName][i] gives the real hours per day spent on subject subjectName for 
            grade i
        DC - as D but uses credits instead of real hours
        
        Dtot - Dtot[subjectName] gives the total hours per day across all grades 
               for subject subjectName
        DCtot - as Dtot but uses credits instead of real hours
        """
        
        #list of all subjects
        K = Dtot.keys()
        
        #make sure we list them in order
        K.sort()
        
        #!!!!!!!!!!!!!!!!! DEBUG !!!!!!!!!!!!!!!!!
        if K==None:
            #this should never execute
            pdb.set_trace()
        #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            
        
        #update the table headers and the totals
        for i, k in enumerate(K):
            self.updateTableHeader(i, k)
            self.updateTable(i,len(D), ["%.1f"%Dtot[k],"%.1f"%DCtot[k]])
        
        #update the per grade table entries
        for j,d in enumerate(D):
            for i, k in  enumerate(K):
                self.updateTable(i, j, ["%.1f"%d[k],"%.1f"%DC[j][k]])
                
        #ensure that the tables have no excess entries      
        self.trimSubjectsTable(Dtot)
        
        #total hours down the columns
        self.updateColumnTotals()
        
    
    def updateColumnTotals(self):
        """
        sum up the hours/credits for each table column and 
        write it to the appropriate cell
        """
        
        #T[0][j] - total hours per day for each grade j
        #T[1][j] - total credits per day for each grade j
        T = [[],[]]
        
        #gather up both tables for easy looping
        tabs = [self.ST, self.STC]
        
        #for each table k
        for k, tt in enumerate(tabs):
            #for each column j
            for j in range(tt.columnCount()):
                #initialize the jth total at 0
                T[k].append(0)
                #loop over rows and accumulate the sum
                for i in range(tt.rowCount()):
                    try:
                        #add the hours/credits from the ith row 
                        T[k][j]+=tt.item(i,j).text().toDouble()[0]
                    except AttributeError:
                        #this adds 0 whenever the cell is empty
                        T[k][j]+=0.0
        
        #expand the table one row and label it 'Total'
        self.updateTableHeader(self.ST.rowCount(), 'Total')
        
        #fill in the total values
        [self.updateTable(self.ST.rowCount()-1, j,["%.1f"%T[0][j],"%.1f"%T[1][j]]) 
            for j,t in enumerate(T[0])]
        

    def blankTable(self):
        """
        write us an empty table
        """
        self.ST.setRowCount(0)
        self.STC.setRowCount(0)
                
    def parseBox(self, tb):
        """
        macro to obtain a list of python strings from a text box
        """
        tb.blockSignals(True)        
        s = str(tb.toPlainText())
        tb.blockSignals(False)   
        return s.split()
        
    def updateTable(self, i,j,S):
        """
        update both the hours (ST) and the credits (STC) table
        with the strings in the list S
        S[0] - update for ST
        S[1] - update for STC
        """
        T = [self.ST, self.STC]
        for k,s in enumerate(S):           
            t = qTItem(s)
            if (i+1)>T[k].rowCount():
                T[k].setRowCount(i+1)
            T[k].setItem(i,j,t)
        
    def updateTableHeader(self, i,s):
        """
        update the ith row label with the string s
        """
        T = [self.ST, self.STC]
        for k,tt in enumerate(T):
            t = qTItem(s)
            
            if (i+1)>T[k].rowCount():
                T[k].setRowCount(i+1)
            T[k].setVerticalHeaderItem(i,t)
        
    def trimSubjectsTable(self, D):
        """
        make sure our table has no extra rows
        """
        T = [self.ST, self.STC]
        n = len(D.keys())
        for k,s in enumerate(T):
            if (n)<T[k].rowCount():
                T[k].setRowCount(n)
                
    def getItem(self, T, i,j):
        try:
            return str(T.item(i,j).text())
        except AttributeError:
            return str(0)
    def getVerticalHeader(self, T, i):
        try:
            return str(T.verticalHeaderItem(i).text())
        except AttributeError:
            return str('')
    def getHorizontalHeader(self, T, i):
        try:
            return str(T.horizontalHeaderItem(i).text())
        except AttributeError:
            return str('')
            
    def write(self):
        """
        write the state of the object to an excel workbook
        """
        self.block(True)
        #if we don't have xlwt we will not write the output
        if xlwt==None:
            return 
        
        cHeaderStyle = xlwt.easyxf('alignment: wrap true, horizontal center; font: bold true; borders: bottom medium') 
        rHeaderStyle = xlwt.easyxf('alignment: wrap true, horizontal center; font: bold true; borders: right medium')     
        rTotalHeaderStyle = xlwt.easyxf('alignment: wrap true, horizontal center; font: bold true; borders: right medium, top medium')     
        centerStyle = xlwt.easyxf('alignment: wrap true, horizontal center') 
        
        cTotalStyle  =xlwt.easyxf('alignment: wrap true, horizontal center; borders: top medium')
        rTotalStyle  =xlwt.easyxf('alignment: wrap true, horizontal center; borders: left medium')
        
        wb = self.parent().wb
        
        shName = 'courseRequirements'
        try:
            shidx = wb._Workbook__worksheet_idx_from_name[shName.lower()]
            sh = wb.get_sheet(shidx)
            print 'sheet exists and was indexed'
        except:
            print 'sheet not found, creating it!'
            sh = wb.add_sheet(shName)
            sh.col(0).width = int(sh.col(0).width*2)
            
        sh._cell_overwrite_ok = True
            

        nR = self.ST.rowCount()
        nC = self.ST.columnCount()
        
#        pdb.set_trace()
        sh.write_merge(0,0,1,nC, 'Hours Per Day',cHeaderStyle)
        sh.write_merge(0,0,nC+2,2*nC+1, 'Credits',cHeaderStyle)
        
        rH = [self.getVerticalHeader(self.ST,i) for i in range(nR)]
        cH = [self.getHorizontalHeader(self.ST,i) for i in range(nC)]
        
        
        headerRow = 1
        
        [sh.row(headerRow).write(i+1, h, cHeaderStyle) for i,h in enumerate(cH)]
        [sh.row(headerRow).write(i+2+nC, h, cHeaderStyle) for i,h in enumerate(cH)]
        
        styles = [rHeaderStyle for i in range(nR)]
        styles[-1] = rTotalHeaderStyle
        [sh.write(i+headerRow+1,0,h, styles[i]) for i,h in enumerate(rH)]
        
        for i in range(nR):
            for j in range(nC):
                if i==nR-1:
                    style = cTotalStyle
                elif j==nC-1:
                    style = rTotalStyle
                else:
                    style = centerStyle
                    
                sh.write(i+headerRow+1, j+1, 
                         self.getItem(self.ST, i,j), style)
                sh.write(i+headerRow+1, j+2+nC, 
                         self.getItem(self.STC, i,j), style)
                         
                
        p = os.path.join(os.getcwd(), self.parent().saveName()+'.xls')
        print 'writing excel to the path '+ p
        wb.save(p)            
        self.block(False)

class FTETabulator(QtGui.QWidget):
    """
    this object tabulates full time equivalent (FTE) from the 
    table of subjects an the hours spent on them, the number of simultaneous sections,
    and the staffing pattern
    """
    def setup(self, fteTable, 
                 subjectsTable, 
                 studentBodyObj, 
                 lcdTotalFTE,
                 periodLengthFun,
                 prepsPerDayFun,
                 blocklist):
        """
        fteTable - the Qtable we'll use for displaying the FTE we compute 
        subjectsTable - table of the hours spent on each subject, by grade
        studentBodyObj - object representing the student body will provide us with
                          enrollment info
        lcdTotalFTE - lcd we'll use to display the upper bound on FTE 
                        (100% attendance, no dropouts or self-contained)
        periodLengthFun - function returns a list of period lengths in minutes
        prepsPerDayFun - function tells us how many periods per day are prep periods
        blocklist - a list of objects that need to be blocked whenever we update this 
                    object
        """
                     
        self.FT = fteTable
        self.ST = subjectsTable
        self.SB = studentBodyObj
        self.lcd = lcdTotalFTE
        
        self.periods = periodLengthFun
        self.ppd = prepsPerDayFun
        self.blocklist = blocklist
        
        self.button_fteChanged = QtGui.QPushButton(self)
        self.button_fteChanged.setVisible(False)
        
        #fte changes whenever the number of sections changes or the allocation 
        #of subjects changes
        [conn(n,sig('valueChanged(int)') ,self.fteChanged) for n in self.SB.nS]
        conn(self.ST,sig('cellChanged(int,int)') ,self.fteChanged)
        
        #begin by updating the fte        
        self.fteChanged()
    
    def fteChanged(self):
        """
        callback that updates the FTE table 
        blocks signals, updates the table and the lcd, unblocks
        """
        [b.blockSignals(True) for b in self.blocklist]
        self.matchTable()
        self.lcd.display(self.fullEnrollmentFTE())
        [b.blockSignals(False) for b in self.blocklist]
    

    def updateTable(self, i,j,s):
        """
        update one table entry i,j with the string s
        """
        t = qTItem(s)
        if (i+1)>self.FT.rowCount():
            self.FT.setRowCount(i+1)
        self.FT.setItem(i,j,t)
        
    def updateTableHeader(self, i,s):
        """
        update the ith row label with the string s
        """
        t = qTItem(s)
        if (i+1)>self.FT.rowCount():
            self.FT.setRowCount(i+1)
        self.FT.setVerticalHeaderItem(i,t)
        
    def getHours(self, i,j):
        """
        return the hours spent on subject i, grade j
        """
        try:
            return self.ST.item(i,j).text().toDouble()[0]
        except AttributeError:
            return 0
            
    def getHeaders(self):
        """
        copy the row labels from the subjects table
        """
        H = []
        for i in range(self.ST.rowCount()):
            try:
                H.append(str(self.ST.verticalHeaderItem(i).text()))
            except:
                H.append('')
        return H
        
    def matchTable(self):
        """
        match the fte to the hours reported in the subjects table
        """
        self.FT.blockSignals(True)
        #resize the FTE table        
        self.FT.setRowCount(self.ST.rowCount())
        
        #match table headers
        H = self.getHeaders()
        for i,h in enumerate(H):
            self.updateTableHeader(i,h)
        
        #make indices for the rows and columns
        R = range(self.ST.rowCount())
        C = range(self.ST.columnCount())
        
        #get some things we need to compute FTE:
            #a list of period lengths in min, 
            #number of preps per day,
            #and the mean period length
        P = self.periods()
        ppd = self.ppd()
        meanPeriod = max([1e-1, round(sum(P)/float(len(P))/60.0,1)])
        
        #we're gonna need the number of sections per grade and the total number 
        #of sections
        N = self.SB.nSections()
        N.append(sum(N))
        N.append(sum(N))
        
        
        #loop over table entries and compute the FTE for each
        for i in R:
            for j in C:
                
                #number of times students must take this course
                nCred = ceil(self.getHours(i,j)/meanPeriod)
                
                #fte is the number of sections times the number of times the course 
                #is taken divided by the number of sections one teacher can teach per day
                fte = N[j]*nCred/max([1,float(len(P)-ppd)])
                
                #in the case of totals, we will compute them directly from the 
                #current table values
                if j==C[-1]:
                    fte = sum([self.FT.item(i,jj).text().toDouble()[0] 
                        for jj in C[:-1]])
                if i==R[-1]:
                    fte = sum([self.FT.item(ii,j).text().toDouble()[0] 
                        for ii in R[:-1]])
                        
                #make the update to the table
                self.updateTable(i,j,"%.1f"%fte)
        self.FT.blockSignals(False)
        self.button_fteChanged.emit(sig('pressed()'))

    def fullEnrollmentFTE(self):
        """
        compute the FTE for an ideal student body with 100% attendance, no 
        self-containeds and no dropouts
        """
        try:
            fteEff = self.FT.item(self.FT.rowCount()-1,self.FT.columnCount()-1).text().toDouble()[0]
        except AttributeError:
            fteEff = 0

        nS = float(self.SB.nSections()[-1])
        nSFull= self.SB.fullNSections()[-1]
        fteFull = fteEff*(nSFull/nS)
        
        return fteFull
        
    def getItem(self, T, i,j):
        try:
            return str(T.item(i,j).text())
        except AttributeError:
            return str(0)
    def getVerticalHeader(self, T, i):
        try:
            return str(T.verticalHeaderItem(i).text())
        except AttributeError:
            return str('')
    def getHorizontalHeader(self, T, i):
        try:
            return str(T.horizontalHeaderItem(i).text())
        except AttributeError:
            return str('')
            
    def write(self):
        """
        write the state of the object to an excel workbook
        """
        
        #if we don't have xlwt we will not write the output
        if xlwt==None:
            return 
        
        cHeaderStyle = xlwt.easyxf('alignment: wrap true, horizontal center; font: bold true; borders: bottom medium') 
        rHeaderStyle = xlwt.easyxf('alignment: wrap true, horizontal center; font: bold true; borders: right medium')     
        rTotalHeaderStyle = xlwt.easyxf('alignment: wrap true, horizontal center; font: bold true; borders: right medium, top medium')     
        centerStyle = xlwt.easyxf('alignment: wrap true, horizontal center') 
        
        cTotalStyle  =xlwt.easyxf('alignment: wrap true, horizontal center; borders: top medium')
        rTotalStyle  =xlwt.easyxf('alignment: wrap true, horizontal center; borders: left medium')
        
        wb = self.parent().wb
        
        shName = 'FTE'
        try:
            shidx = wb._Workbook__worksheet_idx_from_name[shName.lower()]
            sh = wb.get_sheet(shidx)
            print 'sheet exists and was indexed'
        except:
            print 'sheet not found, creating it!'
            sh = wb.add_sheet(shName)
            sh.col(0).width = int(sh.col(0).width*2)
            
        sh._cell_overwrite_ok = True
            

        nR = self.FT.rowCount()
        nC = self.FT.columnCount()
        
        rH = [self.getVerticalHeader(self.FT,i) for i in range(nR)]
        cH = [self.getHorizontalHeader(self.FT,i) for i in range(nC)]
        
        
        headerRow = 0
        
        [sh.row(headerRow).write(i+1, h, cHeaderStyle) for i,h in enumerate(cH)]
        
        
        styles = [rHeaderStyle for i in range(nR)]
        styles[-1] = rTotalHeaderStyle
        [sh.write(i+headerRow+1,0,h, styles[i]) for i,h in enumerate(rH)]
        
        for i in range(nR):
            for j in range(nC):
                if i==nR-1:
                    style = cTotalStyle
                elif j==nC-1:
                    style = rTotalStyle
                else:
                    style = centerStyle
                    
                sh.write(i+headerRow+1, j+1, 
                         self.getItem(self.FT, i,j), style)
                
 
        p = os.path.join(os.getcwd(), self.parent().saveName()+'.xls')
        print 'writing excel to the path '+ p
        wb.save(p)            
        

class budgetTable(QtGui.QWidget):
    """
    implements functionality for the broken-out budget tables in the
    budget tab
    """
    
    def setup(self, tableWidget, 
                         addButton, 
                         deleteButton, 
                         spinBox_deleteLine, 
                         lcd_totalDollars,
                         lcd_totalSTU,
                         dollarsPerSTUFun,
                         buttonBudgetChanged):
        """
        tableWidget - the table of interest
        addButton - button that will add a row
        deleteButton - button that will delete a row
        spinBox_deleteLine - spinBox that tells us which row to delete
        lcd_totalDollars - lcd that displays the total dollars spent for this budget table
        lcd_totalSTU - displays the total cost of this budget table in STU
        dollarsPerSTUFun - returns the conversion rate between dollars and STU
        """
        self.timerID = []
        self.T = tableWidget
        self.AB = addButton
        self.DB = deleteButton
        self.sbDelete = spinBox_deleteLine
        self.lcdDollars = lcd_totalDollars
        self.lcdSTU = lcd_totalSTU
        self.dollarsPerSTU = dollarsPerSTUFun
        self.buttonBudgetChanged = buttonBudgetChanged
        
        self.lcdDollars.setDigitCount(7)
        adjustTableColumns(self.T, 110)
        
        #column dict
        self.cD = dict(zip(['number', 'unitCostDollars', 
                            'unitCostSTU',  'costDollars','costSTU'],
                           [1, 2, 3, 4, 5 ]))

        self.tableChanged()
        self.sbDelete.setValue(self.T.rowCount()-1)

        
        conn(self.AB, sig('pressed()'), self.addLine)
        conn(self.DB, sig('pressed()'), self.deleteLine)
        conn(self.T, sig('cellChanged(int,int)'), self.tableChanged)
        
    def tableChanged(self):
        self.T.blockSignals(True)
        self.buttonBudgetChanged.blockSignals(True)
        self.timerID += [self.startTimer(250)]
        
    def timerEvent(self, tid):
        #we need to kill all timers to avoid bigtime errors
        [self.killTimer(tid) for tid in self.timerID] 
        self.timerID = []
        self.cacheTable()
        self.crunchTable()
        self.updateTable()
        self.buttonBudgetChanged.blockSignals(False)
        self.T.blockSignals(False)
        
        
    def tableIter(self):
        return (range(self.T.rowCount()), range(self.T.columnCount()))
        
    def getEntry(self,i,j):
        """
        return the ith, jth entry of the table as a string
        """
        try:
            
            return str(self.T.item(i,j).text())
        except AttributeError:
            return '0'
            
    def cacheTable(self):
        """
        read the table and put it into lists for easy manipulation
        """
        R,C = self.tableIter()
        T = []
        for i in R:
            T+=[[]]
            for j in C:
                if j==0:
                    T[i]+=[self.getEntry(i,j)]
                else:
                    T[i]+=[float(self.getEntry(i,j))]
        
        self.listTable = list(T)
        return T
        
    def crunchTable(self):
        """
        use the number of staff and dollar unit cost columns to compute values
        for the other columns and totals
        """
        
        T = self.listTable
        DpSTU = float(self.dollarsPerSTU())
#        print DpSTU
        R,C = self.tableIter()
        cD = self.cD
        
        n = [T[i][cD['number']] for i in R]
        uD = [T[i][cD['unitCostDollars']] for i in R]
        
        for i in R:
            T[i][cD['unitCostSTU']] = uD[i]/max([1.0, DpSTU])
            T[i][cD['costDollars']] = n[i]*T[i][cD['unitCostDollars']]            
            T[i][cD['costSTU']] = n[i]*T[i][cD['unitCostSTU']]
            
            
#        print T
            
        totalSTU = sum([T[i][cD['costSTU']] for i in R])
        self.updateTotals(totalSTU)
        
    def updateTable(self):
        """
        update the table using the (presumably changed) cached values
        """
#        print self.listTable
        R,C = self.tableIter()
        floatCols = [3,5]
        for i in R:
            for j in C[1:]:
                x = {True:'%.2f'%self.listTable[i][j],
                     False:'%d'%self.listTable[i][j]}[j in floatCols]             
                self.setEntry(i,j,x )
                
    def updateTotals(self, totalSTU):
        totalDollars = totalSTU*self.dollarsPerSTU()
        self.lcdDollars.display(totalDollars)
        self.lcdSTU.display(totalSTU)
    
    def setEntry(self, i, j, s):
        """
        update ith row, jth column of the table with the string s
        """
        t = qTItem(s)
        self.T.setItem(i,j,t)
    
    def addLine(self):
        """
        append a row to the table
        """
        self.T.setRowCount(self.T.rowCount()+1)
        t = qTItem(str(self.T.rowCount()-1))
        self.T.setVerticalHeaderItem(self.T.rowCount()-1,t)
        self.setSpinBoxMax()
        
    
    def deleteLine(self):
        """
        delete a row from the table
        """
        self.T.removeRow(self.sbDelete.value())
        self.setSpinBoxMax()
    
    def setSpinBoxMax(self):
        """
        adjust the maximum value of the deleteLine spinbox based
        on how many lines we have
        """
        self.sbDelete.setMaximum(self.T.rowCount()-1)
    def getItem(self, T, i,j):
        try:
            return str(T.item(i,j).text())
        except AttributeError:
            return str(0)
    def getVerticalHeader(self, T, i):
        try:
            return str(T.verticalHeaderItem(i).text())
        except AttributeError:
            return str('')
    def getHorizontalHeader(self, T, i):
        try:
            return str(T.horizontalHeaderItem(i).text())
        except AttributeError:
            return str('')    
    def write(self, sh, startRow):
        """
        bounce the table to a spreadsheet
        sh - the spreadsheet we'll write to
        startRow - row in the sheet to begin writing this table
        
        returns:
            endRow - last row in the sheet that we wrote to
        """
        #if we don't have xlwt we will not write the output
        if xlwt==None:
            return 
        
        cHeaderStyle = xlwt.easyxf('alignment: wrap true, horizontal center; font: bold true; borders: bottom medium') 
        rHeaderStyle = xlwt.easyxf('alignment: wrap true, horizontal center; font: bold true; borders: right medium')     
        rTotalHeaderStyle = xlwt.easyxf('alignment: wrap true, horizontal center; font: bold true; borders: right medium, top medium')     
        centerStyle = xlwt.easyxf('alignment: wrap true, horizontal center') 
        
        cTotalStyle  =xlwt.easyxf('alignment: wrap true, horizontal center; borders: top medium')
        rTotalStyle  =xlwt.easyxf('alignment: wrap true, horizontal center; borders: left medium')
        
        sh._cell_overwrite_ok = True
            

        nR = self.T.rowCount()
        nC = self.T.columnCount()
        
        rH = [self.getVerticalHeader(self.T,i) for i in range(nR)]
        cH = [self.getHorizontalHeader(self.T,i) for i in range(nC)]
        
        [sh.row(startRow).write(i+1, h, cHeaderStyle) for i,h in enumerate(cH)]
              
        styles = [rHeaderStyle for i in range(nR)]
        styles[-1] = rTotalHeaderStyle
        [sh.write(i+startRow+1,0,h, styles[i]) for i,h in enumerate(rH)]
        
        for i in range(nR):
            for j in range(nC):
                if i==nR-1:
                    style = cTotalStyle
                elif j==nC-1:
                    style = rTotalStyle
                else:
                    style = centerStyle
                    
                sh.write(i+startRow+1, j+1, 
                         self.getItem(self.T, i,j), style)
                         
        return startRow+1+nR
                 
class budgetTab(QtGui.QWidget):
    """
    class to manage the budget tab, and communicate between the various staff budget tables
    """
    def setup(self, budgetTableList,
                 sbMeanSalary,
                 sbPctOverhead,
                 tableStartingBudget,
                 buttonBudgetChanged,
                 lcd_teachingSTU,
                 lcd_totalSTU,
                 lcdListNonPersonnel,
                 lcdListStartingBudget,
                 lcdListRegularTeaching,
                 lcdListCertifiedNonTeaching,
                 lcdListClassified,
                 lcdListFinalBudget):
        """
             budgetTableList - list of budget tables whose signals we must catch
             sbMeanSalary - spinBox holding the mean secondary teacher salary (1 STU)
             sbPctOverhead - spinBox giving the percent of school funded budget
                             that is non-personnel cost
             tableStartingBudget - table holding the starting values for the school 
                                 funded and centrally funded budgets
             buttonBudgetChanged - will press this (invisible) button whenever the budget updates
             lcd_teachingSTU - display the STU allocated to teaching on the scheduling panel
             lcd_totalSTU - display total STU in budget on the scheduling tab
             
             ---- all following lists have the format: [lcd_costInDollars, lcd_costInSTU -------
             lcdListNonPersonnel - for the nonPersonnel costs
             lcdListStartingBudget - for the starting staffing budget available
             lcdListRegularTeaching - for the regular teaching staff
             lcdListCertifiedNonTeaching - for the certified non teaching staff
             lcdListClassified - for the classified staff
             lcdListFinalBudget - budget remaining after staffing costs
        """
        self.timerID = []
        self.BTList = budgetTableList
        self.sbMeanSalary = sbMeanSalary
        self.sbPctOverhead = sbPctOverhead
        self.tableStartingBudget = tableStartingBudget
        self.buttonBudgetChanged = buttonBudgetChanged
        self.lcd_teachingSTU = lcd_teachingSTU
        self.lcd_totalSTU = lcd_totalSTU
        self.lcdListNonPersonnel = lcdListNonPersonnel
        self.lcdListStartingBudget = lcdListStartingBudget
        self.lcdListRegularTeaching = lcdListRegularTeaching
        self.lcdListCertifiedNonTeaching = lcdListCertifiedNonTeaching
        self.lcdListClassified = lcdListClassified
        self.lcdListFinalBudget = lcdListFinalBudget
        
        self.lcdDict = dict(zip('nonPersonnel starting regular certified classified final'.split(),
                                [lcdListNonPersonnel, lcdListStartingBudget,
                                    lcdListRegularTeaching,lcdListCertifiedNonTeaching,
                                    lcdListClassified, lcdListFinalBudget]))
                                    
        [L[0].setDigitCount(8) for L in self.lcdDict.values()]
                                    
        self.lcdStaffDict = dict(zip('regular certified classified'.split(),
                                [lcdListRegularTeaching,lcdListCertifiedNonTeaching,
                                    lcdListClassified]))
                                    
        for b in self.BTList:
            b.dollarsPerSTU = self.dollarsPerSTU
                                    
        [conn(T.T, sig('cellChanged(int,int)'), self.budgetChanged) for T in self.BTList]
        
        conn(self.sbMeanSalary, sig('valueChanged(double)'), self.tableChanged)
        conn(self.tableStartingBudget, sig('cellChanged(int,int)'), self.budgetChanged)
        conn(self.sbPctOverhead, sig('valueChanged(double)'), self.budgetChanged)
        
        
#        self.budgetChanged()
                                    
    def budgetChanged(self):
#        print "budgetChanged"
        self.tableStartingBudget.blockSignals(True)
        self.timerID += [self.startTimer(400)]
#        print self.timerID
        
    def tableChanged(self):
#        print 'table changed'
        self.tableStartingBudget.blockSignals(True)
        self.timerID += [self.startTimer(250)]
#        print self.timerID
        
        
    def timerEvent(self, tid):
#        print 'budget changed timer fired!'
        #we need to kill all timers to avoid bigtime errors
        [self.killTimer(TID) for TID in self.timerID] 
        self.timerID = []
        
        self.updateStartingBudgetSTU()
        self.updatePersonnelCosts()
        self.updateStartingBudget()
        self.updateFinalBudget()
        
        self.lcd_teachingSTU.display(self.lcdDict['regular'][1].value())
                                        #+self.lcdDict['final'][1].value())
        self.lcd_totalSTU.display(self.lcdDict['starting'][1].value())
        
        
        for T in self.BTList:
            T.T.blockSignals(True)
            T.cacheTable()
            T.crunchTable()
            T.updateTable()
            T.T.blockSignals(False)
        
        
        
        self.buttonBudgetChanged.emit(sig('pressed()'))
        self.tableStartingBudget.blockSignals(False)
    
    def updateStartingBudgetSTU(self):
#        print 'updating starting stu'
        s = self.schoolFundedBudget()/self.dollarsPerSTU()
        c = self.centralFundedBudget()/self.dollarsPerSTU()
        self.tableStartingBudget.setItem(0,1,qTItem('%0.f'%c))
        self.tableStartingBudget.setItem(1,1,qTItem('%1.f'%s))
        
        
        
    def updatePersonnelCosts(self):
        self.lcdDict['nonPersonnel'][0].display(self.overhead())
        self.lcdDict['nonPersonnel'][1].display(self.overhead()/self.dollarsPerSTU())
        
    def updateStartingBudget(self):
        s = self.startingBudget()
        self.lcdDict['starting'][0].display(s)
        self.lcdDict['starting'][1].display(s/self.dollarsPerSTU())
    
    def updateFinalBudget(self):
        f = self.finalBudget()
        self.lcdDict['final'][0].display(f)
        self.lcdDict['final'][1].display(f/self.dollarsPerSTU())

    def schoolFundedBudget(self):
        try:
            return self.tableStartingBudget.item(1,0).text().toDouble()[0]
        except AttributeError:
            return 0
            
        
    def centralFundedBudget(self):
        try:
            return self.tableStartingBudget.item(0,0).text().toDouble()[0]
        except AttributeError:
            return 0
            
    def startingBudget(self):
        return self.schoolFundedBudget()+self.centralFundedBudget()-self.overhead()
    
    def finalBudget(self):
        return self.startingBudget() - self.totalStaffingCosts()
    
    def totalStaffingCosts(self):
        return sum([self.lcdStaffDict[k][0].value() for k in self.lcdStaffDict.keys()])

    def dollarsPerSTU(self):
        """
        exchange rate between dollars and Secondary Teacher Units
        """
        dps = float(max([1,self.sbMeanSalary.value()]))
#        print 'dps is %0.f'%dps
        return dps
        
    def overhead(self):
        """
        compute non-personell costs in dollars
        """
        return .01*self.sbPctOverhead.value()*self.schoolFundedBudget()
    def getItem(self, T, i,j):
        try:
            return str(T.item(i,j).text())
        except AttributeError:
            return str(0)
    def getVerticalHeader(self, T, i):
        try:
            return str(T.verticalHeaderItem(i).text())
        except AttributeError:
            return str('')
    def getHorizontalHeader(self, T, i):
        try:
            return str(T.horizontalHeaderItem(i).text())
        except AttributeError:
            return str('')
            
    def write(self):
        """
        write the state of the object to an excel workbook
        """
        #if we don't have xlwt we will not write the output
        if xlwt==None:
            return 
        
        cHeaderStyle = xlwt.easyxf('alignment: wrap true, horizontal center; font: bold true; borders: bottom medium') 
        rHeaderStyle = xlwt.easyxf('alignment: wrap true, horizontal center; font: bold true; borders: right medium')     
        rTotalHeaderStyle = xlwt.easyxf('alignment: wrap true, horizontal center; font: bold true; borders: right medium, top medium')     
        centerStyle = xlwt.easyxf('alignment: wrap true, horizontal center') 
        
        
        wb = self.parent().wb
        
        shName = 'budget'
        try:
            shidx = wb._Workbook__worksheet_idx_from_name[shName.lower()]
            sh = wb.get_sheet(shidx)
            print 'sheet exists and was indexed'
        except:
            print 'sheet not found, creating it!'
            sh = wb.add_sheet(shName)
            sh.col(0).width = int(sh.col(0).width*1.5)
            
        sh._cell_overwrite_ok = True
            
        sh.write_merge(0,0,0,1,'Summary',cHeaderStyle)
        
        sh.write(1,0,'mean secondary teacher salary, $',rHeaderStyle)
        sh.write(1,1, self.sbMeanSalary.value() ,centerStyle)
        
        sh.write(2,0,'non-personnel cost, % of school funded budget',rHeaderStyle)
        sh.write(2,1, self.sbPctOverhead.value() ,centerStyle)
        
        
        
        nR = self.tableStartingBudget.rowCount()
        nC = self.tableStartingBudget.columnCount()
        
        rH = [self.getVerticalHeader(self.tableStartingBudget,i) for i in range(nR)]
        cH = [self.getHorizontalHeader(self.tableStartingBudget,i) for i in range(nC)]
        
        
        ########### write the starting budget table #########
        headerRow = 5
        sh.write_merge(headerRow-1,headerRow-1, 0, 3,'Starting Budget', cHeaderStyle)
        [sh.row(headerRow).write(i+1, h, cHeaderStyle) for i,h in enumerate(cH)]
        styles = [rHeaderStyle for i in range(nR)]
        styles[-1] = rTotalHeaderStyle
        [sh.write(i+headerRow+1,0,h, styles[i]) for i,h in enumerate(rH)]
        
        for i in range(nR):
            for j in range(nC):
                style = centerStyle
                sh.write(i+headerRow+1, j+1, 
                self.getItem(self.tableStartingBudget, i,j), style)
        #######################################################  
        
        rw = headerRow+nC+1
        sh.write(rw, 0, 'non-personnel cost, ($, STU)', rHeaderStyle)
        [sh.write(rw, i+1, self.lcdListNonPersonnel[i].value(), centerStyle) for i in range(2)]
        
        sh.write(rw+1, 0, 'staffing budget, ($, STU)', rHeaderStyle)
        [sh.write(rw+1, i+1, self.lcdListStartingBudget[i].value(), centerStyle) for i in range(2)]
        
        sh.write_merge(rw+3,rw+3, 0, 6, 'Staffing Breakdown', cHeaderStyle)
        
        staffTitles = 'certified non-teaching, regular teaching, classified'.split(',')
        
        startRow = rw+5
        for i,st in enumerate(staffTitles):
            sh.write_merge(startRow,startRow, 0,3, st,cHeaderStyle )
            endRow = self.BTList[i].write(sh, startRow+1)
            startRow = endRow+2
            
        sh.write_merge(startRow,startRow, 0, 6, 'Budget Remaining After Staffing', cHeaderStyle)
        sh.write(startRow+1,0, 'final budget, ($, stu)', rHeaderStyle)
        [sh.write(startRow+1, i+1, self.lcdListFinalBudget[i].value(), centerStyle) for i in range(2)]
            
        p = os.path.join(os.getcwd(), self.parent().saveName()+'.xls')
        print 'writing excel to the path '+ p
        wb.save(p)   

class teacherPool(QtGui.QTableWidget):
    """
    the teacherPool class represents all the teachers we'll need to 
    cover the schedule designed in the scheduling tab. from the teacher pool, 
    we assign teachers to superteams, which are teachers who share a prep block
    """
    def setup(self, FTETab, STArray):
        """
        FTETab - an FTETabulator object that will keep us informed on the 
            number and kind of teacher, by grade
        STArray - a superteamArray object that will hold the superteams we 
                create
        """
        self.fteTable = FTETab.FT
        self.STArray = STArray
        
        
        #dummy button to be pushed whenever the teams change
        self.button_teamsChanged = QtGui.QPushButton(self)
        self.button_teamsChanged.setVisible(False)
        
        #whenever the composition of teachers changes, we reset the pool
        conn(FTETab.button_fteChanged, sig("pressed()"), self.setPool )
        #whenever a cell is doubleClicked, we draft a teacher into the active superteam
        conn(self, sig("cellDoubleClicked(int,int)"), self.draft)
        #whenever a superteam is added, we need to connect some signals from it 
        conn(self.STArray.button_teamAdded, sig('pressed()'), self.connectTeams)
        self.connectTeams()
        
        self.setGeometry(QtCore.QRect(10, 60, 321, 291))
        self.setRowCount(1)
        self.setColumnCount(4)
        adjustTableColumns(self, 70)
        self.setObjectName(_fromUtf8("table_teacherPool"))
        
        #columns are the grade
        [self.setHorizontalHeaderItem(i,
            qTItem(str(i+9))) for i in range(4)]
        
        #need to keep us from entering the table if we click it
        self.setSelectionMode(0)
        self.setEditTriggers(QtGui.QAbstractItemView.NoEditTriggers)
        self.setPool()
        
        
    def connectTeams(self):
        [conn(t.L, sig("itemDoubleClicked(QListWidgetItem*)"), self.undraft) for t in self.STArray.teams()]
        
    def setSubjects(self):
        """
        set row labels to the subject the teacher teaches
        """
        [self.setVerticalHeaderItem(i,qTItem(s)) 
            for i,s in  enumerate(self.getSubjects())]
            
    def setPool(self):
        """
        set the number of teachers by subject and grade according to the fte
        table
        """
        self.setRowCount(self.fteTable.rowCount()-1)
        self.setSubjects()
        for r in self.rows():
            for c,g in self.grades():
                self.setItem(r,c,
                             qTItem(str(it2int(self.fteTable.item(r,c).clone()))))
                
        self.STArray.clearTeams() 
        
    def draft(self, r,c):
        """
        decrement the teacher count in the cell that was doubleClicked
        add that teacher to the active team 
        """
        memberLabel =  str(self.grade(c))+"_"+self.subj(r)
        self.STArray.draftMember(memberLabel)
        self.setItem(r,c,
            qTItem(str(it2int(self.item(r,c))-1)))
        
        self.button_teamsChanged.emit(sig('pressed()'))
        
        
    def undraft(self, it):
        """
        if a teacher can be drafted into a superteam, we can undraft them back
        into the pool
        
         it - the item from the superteam's list that was doubleclicked
        """
        if str(it.text())!=teamBreakString:
            g, s = str(it.text()).split("_")
            c = int(g)-9
            r = find([j==s for j in self.getSubjects()])[0]
            self.setItem(r,c,
                qTItem(str(it2int(self.item(r,c))+1)))
        self.button_teamsChanged.emit(sig('pressed()'))
    
    def getSubjects(self):
        """
        return a list of subjects from the fteTable
        """
        s = []
        for i in self.rows():
            s+= [str(self.fteTable.verticalHeaderItem(i).text())]
        return s
        
    def subj(self, r):
        return self.getSubjects()[r]
    def grade(self, c):
        return self.grades()[c][1]
        
    def rows(self):
        return range(self.fteTable.rowCount()-1)
    def grades(self):
        return zip(range(4), range(9,13))
    
class superteam(QtGui.QGroupBox):
    """
    represents a group of teachers who will prep together
    it consists of a list of teachers in the superteam, and 
    divisions within the superteam that are parsed into teams
    """
    def setup(self, teamNr):
        """
        teamNr - an int identifying the superteam
        """
        self.setGeometry(QtCore.QRect(130, 130, 231, 161))
        self.setObjectName(qstr("superteam"+str(teamNr)))
        self.setTitle("superteam "+str(teamNr))
        self.setSizePolicy(QtGui.QSizePolicy.Preferred, QtGui.QSizePolicy.Preferred)
        self.L = QtGui.QListWidget(self)
        self.L.setGeometry(QtCore.QRect(10, 20, 101, 201))
        self.L.setObjectName(qstr("superteamList"+str(teamNr)))
        self.L.setDragEnabled(True)
        self.L.setDragDropMode(QtGui.QAbstractItemView.InternalMove)
        self.L.setSizePolicy(QtGui.QSizePolicy.Preferred, QtGui.QSizePolicy.Preferred)
        
        #clicking this button will make it the active superteam
        self.activateButton = QtGui.QPushButton(self)
        self.activateButton.setGeometry(QtCore.QRect(130, 20, 60, 23))
        self.activateButton.setObjectName(qstr("superteamActivateButton"+str(teamNr)))
        self.activateButton.setText("activate")
        self.activateButton.setCheckable(True)
        
        #clicking this button will insert a team break into the superteam
        self.splitButton = QtGui.QPushButton(self)
        self.splitButton.setGeometry(QtCore.QRect(130, 50, 60, 23))
        self.splitButton.setObjectName(qstr("superteamSplitButton"+str(teamNr)))
        self.splitButton.setText("split")
        
        #this button will alert other objects whenever the team changes
        self.button_teamChanged = QtGui.QPushButton(self)
        self.button_teamChanged.setVisible(False)
        
        #double clicking an item on the list will remove it
        conn(self.L, sig("itemDoubleClicked(QListWidgetItem*)"), self.delMember)
        #clicking an item on the list makes this the active superteam
        conn(self.L, sig("itemClicked(QListWidgetItem*)"), self.activate)
        conn(self.L, sig("clicked(QModelIndex&)"), self.activate)
        
        conn(self.splitButton, sig('pressed()'), self.insertTeamBreak)
        
        self.timerID = []
    
    def insertTeamBreak(self):
        """
        the insertion of a special string into the list will be used to designate
        teams
        """
        idx = self.L.selectedIndexes()
        if len(idx)>0:
            self.L.insertItem(idx[0].row(), QtGui.QListWidgetItem(teamBreakString))
            self.button_teamChanged.emit(sig('pressed()'))
        
    def addMember(self, memberLabel):
        """
        memberLabel - a str that identifies a teacher's grade and subject
        """
        it = QtGui.QListWidgetItem()
        it.setText(qstr(memberLabel))
        self.L.addItem(it)
        
    def delMember(self, it):
        """
        it - an item from the superteam list to be deleted
        """
        self.activate()
        if str(it.text())!=teamBreakString:
            grade, subj = str(it.text()).split("_")
            print grade+", " +subj + " deleted"
        self.tmp_it = it
        self.timerID += [self.startTimer(50)]
        
    def timerEvent(self, E):
        #we need to kill all timers to avoid bigtime errors
        [self.killTimer(tid) for tid in self.timerID] 
        self.timerID = []
        
        self.L.takeItem(self.L.row(self.tmp_it))
        
    def clear(self):
        """
        remove all items from the list
        """
        while self.L.count()>0:
            self.L.takeItem(self.L.count()-1)
    def isActive(self):
        return self.activateButton.isChecked()
    def activate(self):
        self.activateButton.click()
    def parseTeams(self):
        """
        divide the superteam into teams based on the locations of certain strings 
        in the list
        """
        S = [str(self.L.item(i).text()) for i in range(self.L.count())]
        S+=[teamBreakString]
        T = []
        t = []
        for s in S:
            if s==teamBreakString:
                if t==[]:
                    continue
                else:
                    T+=[t]
                    t=[]
            else:
                t+=[s]
        return T
        
class superteamArray(QtGui.QWidget):
    def setup(self, nr, periodsPerDayFun):
        
        self.nPeriods = periodsPerDayFun
        self.button_teamAdded = QtGui.QPushButton(self)
        self.button_teamAdded.setVisible(False)
        self.button_teamDeleted = QtGui.QPushButton(self)
        self.button_teamDeleted.setVisible(False)
        self.button_teamsChanged = QtGui.QPushButton(self)
        self.button_teamsChanged.setVisible(False)
        
     
        #print "setting up superteam array layout"
        self.setGeometry(QtCore.QRect(10, 50, 680, 680))
        self.setObjectName(_fromUtf8("superteamArray"+str(nr)))
        self.gridSuperteams = QtGui.QGridLayout(self)
        self.gridSuperteams.setMargin(0)
        self.gridSuperteams.setObjectName(_fromUtf8("gridSuperteams"))
        
       #print "adding teams, %d"%self.nPeriods()
        self.teamCount = 0
        self.active = 0
        self.buttonGroup = QtGui.QButtonGroup(self)
        [self.addTeam(i) for i in range(self.nPeriods())]
        self.team(0).activateButton.setChecked(True)
        
    def draftMember(self, memberLabel):
        self.activeTeam().addMember(memberLabel)
            
    def addTeam(self, idx):
        #print "adding superteam"
        ST = superteam()
        
        ST.setup(idx)
        self.buttonGroup.addButton(ST.activateButton)
        
        conn(ST.button_teamChanged, sig('pressed()'), self.teamsChanged)
            
        x = self.idx2sub(idx)
        self.gridSuperteams.addWidget(ST, x[0],x[1])
        
        self.teamCount+=1
        self.button_teamAdded.emit(sig('pressed()'))
        self.teamsChanged()

    def removeTeam(self):
        #print "removing team"
        I = self.gridSuperteams.takeAt(self.teamCount-1)
        w = I.widget()
        self.gridSuperteams.removeItem(I)
        w.deleteLater()
        
        self.teamCount-=1
        self.button_teamDeleted.emit(sig('pressed()'))
        self.teamsChanged()

    def teamsChanged(self):
        self.button_teamsChanged.emit(sig('pressed()'))        
    def clearTeams(self):
        [t.clear() for t in self.teams()]
        self.teamsChanged()
        
    def updateNTeams(self):
        if self.teamCount<self.nPeriods():
            self.addTeam(self.teamCount)
            self.updateNTeams()
        elif self.teamCount>self.nPeriods():
            self.removeTeam()
            self.updateNTeams()
   
                
    def team(self, idx):
        return self.gridSuperteams.itemAt(idx).widget()
        
    def teams(self):
        return [self.team(i) for i in range(self.teamCount)]
        
    def activeTeam(self):
        for i, t in enumerate(self.teams()):
            if t.isActive():
                return t

    def idx2sub(self, idx):
        return (idx/3, idx-3*(idx/3))
    
class team(QtGui.QGroupBox):
    def setup(self, superteamNr, teamNr, roster):
        self.superteamNr = superteamNr
        self.teamNr = teamNr
        
        self.setGeometry(QtCore.QRect(110, 90, 191, 271))
        self.setObjectName(qstr("team"+self.idStr()))
        self.setTitle("superteam %d, team %d"%(self.superteamNr, self.teamNr))
        self.setSizePolicy(QtGui.QSizePolicy.Preferred, QtGui.QSizePolicy.Preferred)
        
        self.L = QtGui.QListWidget(self)
        self.L.setGeometry(QtCore.QRect(10, 70, 171, 192))
        self.L.setObjectName(_fromUtf8("list_team"))
        self.L.setObjectName(qstr("teamList"+self.idStr()))
        self.L.setDragEnabled(True)
        self.L.setDragDropMode(QtGui.QAbstractItemView.InternalMove)
        self.L.setSizePolicy(QtGui.QSizePolicy.Preferred, QtGui.QSizePolicy.Preferred)
        for s in roster:
            self.L.addItem(QtGui.QListWidgetItem(s))
        
        self.lineEdit_teamName = QtGui.QLineEdit(self)
        self.lineEdit_teamName.setGeometry(QtCore.QRect(40, 40, 113, 22))
        self.lineEdit_teamName.setObjectName(_fromUtf8("lineEdit_teamName"))
        self.lineEdit_teamName.setText("team"+self.idStr())
        
        self.label_teamName = QtGui.QLabel(self)
        self.label_teamName.setGeometry(QtCore.QRect(60, 20, 71, 16))
        self.label_teamName.setObjectName(_fromUtf8("label_teamName"))
        
    def idStr(self):
        return "_"+str(self.superteamNr)+"_"+str(self.teamNr)
    
    def clear(self):
       while self.L.count()>0:
           self.L.takeItem(self.L.count()-1)
    def teamName(self):
        return str(self.lineEdit_teamName.text())
    def nMembers(self):
        return self.L.count()

class teamArray(QtGui.QFrame):
    def setup(self, STArray, tPool):
        self.STArray = STArray
        self.tPool = tPool
        #self.frame_teamArray = QtGui.QFrame(self.teamTab)
        self.setGeometry(QtCore.QRect(10, 10, 800, 800))
        self.setFrameShape(QtGui.QFrame.StyledPanel)
        self.setFrameShadow(QtGui.QFrame.Raised)
        self.setObjectName(_fromUtf8("frame_teamArray"))
        
        self.button_teamsChanged = QtGui.QPushButton(self)
        self.button_teamsChanged.setVisible(False)
        
        self.button_namesChanged = QtGui.QPushButton(self)
        self.button_namesChanged.setVisible(False)
        
        self.teamCounts = []
        self.timerID = []
        
        self.gridTeams = QtGui.QGridLayout(self)
        self.gridTeams.setMargin(0)
        self.gridTeams.setObjectName(_fromUtf8("gridTeams"))
        
#        conn(self.STArray.button_teamAdded, sig('pressed()'), self.refresh)
#        conn(self.STArray.button_teamDeleted, sig('pressed()'), self.refresh)
        conn(self.tPool.button_teamsChanged, sig('pressed()'), self.refresh)
        conn(self.STArray.button_teamsChanged, sig('pressed()'), self.refresh)
        
        self.refresh()
    
    def refresh(self):
        if self.timerID ==[]:
            self.timerID += [self.startTimer(100)]
        
    def timerEvent(self, E):  
        #we need to kill all timers to avoid bigtime errors
        [self.killTimer(tid) for tid in self.timerID] 
        self.timerID = []
        
        self.removeTeams()
        for i,ST in enumerate(self.STArray.teams()):
            T = ST.parseTeams()
            self.addTeams(i,T)
            
        self.button_teamsChanged.emit(sig('pressed()'))
        
    def addTeams(self, superteamNr, teamList):
        T = [team() for t in teamList]
        [T[i].setup(superteamNr, i, t) for i,t in enumerate(teamList)]
        
        [self.gridTeams.addWidget(t, t.superteamNr,t.teamNr) for t in T]
        [conn(t.lineEdit_teamName, sig('editingFinished()'), self.teamNamesChanged) for t in T]
        self.button_teamsChanged.emit(sig('pressed()'))
        
    def removeTeams(self):
        while(self.gridTeams.count()>0):
            I = self.gridTeams.takeAt(self.gridTeams.count()-1)
            w = I.widget()
            self.gridTeams.removeItem(I)
            w.deleteLater()
    def teamNamesChanged(self):
        self.button_namesChanged.emit(sig('pressed()'))
    def teams(self):
        """
        return a tuple of lists:
            (teamNames, superteamNumbers, teamNumbers)
        """
        nm = []
        st = []
        nr = []
        nSec = []
        for i in range(self.gridTeams.count()):
            w = self.gridTeams.itemAt(i).widget()
            nm+=[w.teamName()]
            st+=[w.superteamNr]
            nr+=[w.teamNr]
            nSec+=[w.nMembers()]
        return (nm, st, nr, nSec)
        
        
class studentPool(QtGui.QTableWidget):
    def setup(self, SB, SAT, nPeriodsFun):
        """
        SB - studentBody object will tell us how many sections of each
            grade we have
        SAT - studentAllocationTable object tracks which teams we have assigned
                students to
        """
        self.SB = SB
        self.SAT = SAT
        self.nPeriods = nPeriodsFun
        
        self.setGeometry(QtCore.QRect(300, 70, 1381, 151))
        
        self.setRowCount(4)
        [self.setVerticalHeaderItem(i,qTItem(str(i+9))) for i in range(4)]
        
        self.button_refreshed = QtGui.QPushButton(self)
        self.button_refreshed.setVisible(False)
        
        self.dialog_allocate = QtGui.QInputDialog(self)
        self.dialog_error = QtGui.QErrorMessage(self)
        
        conn(self, sig("cellDoubleClicked(int,int)"), self.draft)
        conn(self.SB.updateButton, sig('pressed()'), self.refresh)
        
        
        self.setSelectionMode(0)
        self.setEditTriggers(QtGui.QAbstractItemView.NoEditTriggers)
        
        self.refresh()
        
    def refresh(self):
        N  =  self.SB.nSections()
        self.setColumnCount(self.nPeriods())
        for i in range(self.nPeriods()):
            for j in range(4):
                self.setItem(j,i,qTItem(str(N[j])))
        [self.setColumnWidth(i,170) for i in range(self.columnCount())]
        self.button_refreshed.emit(sig('pressed()'))
    
    def SATCleared(self):
        self.button_refreshed.blockSignals(True)
        self.refresh()
        self.button_refreshed.blockSignals(False)
        
    def undraft(self,per, grade, n):
        c = per-1
        r = grade-9
        N = it2int(self.item(r,c))
        N+=n
        self.setItem(r,c, qTItem(str(N)))
        
    def draft(self, r,c):
        n = self.dialog_allocate.getInt(self, 'allocate student sections', 
                                        'sections to allocate?', 0,
                                        0, it2int(self.item(r,c)) )
        per = c+1
        grade = r+9
        if not self.SAT.isPrep(per,grade):
            valid, nMax, nAlloc = self.SAT.allocationValid(per, grade, n[0])
            if valid:
                m = n[0]
            elif(nMax-nAlloc)>0:
                m = nMax-nAlloc
            else:
                self.dialog_error.showMessage("attempt to overbook a team; team can only handle as many sections in one period as it as members")
                return
            N = it2int(self.item(r,c))-m
            self.setItem(r,c, qTItem(str(N)))
            self.SAT.allocateSection(per, grade,  m)
                
        
        
class studentAllocationTable(QtGui.QTableWidget):
    def setup(self, SP, TA, nPeriodsFun):
        """
        SP - studentPool object will keep track of sections for us
        TA - teamArray object tells us what teams are available
        """
        self.SP = SP
        self.TA = TA
        self.nPeriods = nPeriodsFun
        self.STPreps = []
        
        self.draftingRow = 0
        
        self.setGeometry(QtCore.QRect(150, 230, 1521, 661))
        
        self.button_cleared = QtGui.QPushButton(self)
        self.button_cleared.setVisible(False)
        
        
        conn(self, sig('cellDoubleClicked(int,int)'), self.clearCell)
#        conn(self,sig("cellClicked(int,int)"), self.setDraftingCell)
        conn(self.TA.button_teamsChanged, sig('pressed()'), self.refresh)
        conn(self.TA.button_namesChanged, sig('pressed()'), self.refreshNames)
        conn(self.SP.button_refreshed, sig('pressed()'), self.SPRefreshed)
        
        self.setEditTriggers(QtGui.QAbstractItemView.NoEditTriggers)
        
        self.refresh()
        
    def clear(self):
        for i in range(self.rowCount()):
            for j in range(self.columnCount()):
                self.setItem(i,j,qTItem(''))
        
        self.button_cleared.emit(sig('pressed()'))
    
    def SPRefreshed(self):
        self.button_cleared.blockSignals(True)
        self.clear()
        self.button_cleared.blockSignals(False)
           
    def refresh(self):
        self.setColumnHeaders()
        nm, st, nr, nSec= self.TA.teams()
        self.setRowHeaders(nm)
        self.clear()
        
        [self.setItem(i,0,qTItem(str(n))) for i,n in enumerate(nSec)]
        
        self.STPreps = []
        for i,s in enumerate(st):
            self.STPreps+=[int(s)]
            self.setItem(i, s+1, qTItem('PREP'))
        
            
    def refreshNames(self):
        nm, st, nr, nSec= self.TA.teams()
        print 'refreshing team names in student allocation'
        print nm
        self.setRowHeaders(nm)
    def isPrep(self, per, grade):
        
        try:
            s = str(self.item(self.currentRow(), per).text())
        except AttributeError:
            s = ''
        if s=='PREP':
            return True
        else:
            return False
        
    def allocateSection(self,per, grade, n):
        r = self.currentRow()
        N = self.parseCell(r,per)
        N[grade-9]+=n
        self.setItem(r,per, qTItem(self.cellStr(N)))

    def allocationValid(self, per, grade, n):
        r = self.currentRow()
        nMax = it2int(self.item(r,0))
        N = self.parseCell(r,per)
        nAlloc = 0
        for m in N:
            nAlloc+=m
        if n>(nMax-nAlloc):
            valid = False
        else:
            valid = True
        return (valid, nMax, nAlloc)
        
    def clearCell(self, r,c):
        try: 
            s = str(self.item(r,c).text())
        except AttributeError:
            return
        if s!='PREP':
            N = self.parseCell(r,c)
            [self.SP.undraft(c, i+9,n) for i,n in enumerate(N)]
            self.setItem(r,c,qTItem(''))
        
    def parseCell(self, r,c):
        try:
            S = str(self.item(r,c).text())
        except AttributeError:
            S = ''
        if S=='':
            return [0,0,0,0]
        else:
            L = S.split()
            return [int(s.split('x')[0]) for s in L]
    def cellStr(self, N):
        L = [str(n)+'x'+str(i+9) for i,n in enumerate(N)]
        S = ''
        for s in L:
            S+=s+' '
        return S
            
    def setColumnHeaders(self):
        headers = 'sec/period'.split()
        for i in range(self.nPeriods()):
            headers+=[str(i+1)]
        
        self.setColumnCount(len(headers))
#        print 'column headers are '
#        print headers
        [self.setColumnWidth(i,170) for i in range(1,self.columnCount())]
        [self.setHorizontalHeaderItem(i,qTItem(h)) for i,h in enumerate(headers)]
    
    def setRowHeaders(self, nm):
        self.setRowCount(len(nm))
        [self.setVerticalHeaderItem(i,qTItem(n)) for i,n in enumerate(nm)]
    

        
class startMasterfulScheduling(QtGui.QMainWindow):
    """
    this is the main application window
    it contains all elements of the ui, and (currently) the functionality of
    the main scheduling apparatus. 
    it creates and hooks in the studentBody, subjectsParser, and FTETabulator 
    objects to complete its functionality
    """
    def __init__(self, parent=None):
        ############# initialize some miscellany ###########        
        self.spinCount = 0
        self.timerID = []
        self.wb = xlwt.Workbook()
        ####################################################
        
        #usual Qt setup 
        QtGui.QWidget.__init__(self, parent)
        self.ui = Ui_MainWindow()
        self.ui.setupUi(self)
        
        
        #dummy buttons for updates of totalEnrollment 
        #and number of sections
        self.TEValueChangedButton = QtGui.QPushButton()
        self.nSValueChangedButton = QtGui.QPushButton()
        self.buttonBudgetChanged = QtGui.QPushButton()
        self.buttonNPeriodsChanged = QtGui.QPushButton()
        self.buttonNPeriodsChanged.setVisible(False)
        
        
        #################################################################
        ############# setup scroll area ################################
        self.sa = QtGui.QScrollArea(self)
        self.sa.setGeometry(QtCore.QRect(10, 10, 800, 600))
        self.sa.setWidgetResizable(False)
        self.sa.setObjectName("scrollArea")
        
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.sa.sizePolicy().hasHeightForWidth())
       
        
        self.ui.allWidget.setParent(self.sa)
        self.sa.setWidget(self.ui.allWidget)
       
        self.setCentralWidget(self.sa)
        self.ui.centralwidget.deleteLater()
        #################################################################
        
        #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
        #::::::::::::::fix table column widths ::::::::::::::::::::::::::
        adjustTableColumns(self.ui.table_courseRequirements, 50)
        adjustTableColumns(self.ui.table_courseRequirements_credits, 50)
        adjustTableColumns(self.ui.table_FTE, 77)
        #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
        
        
        #################### initialize major ui components ############
        self.initStudentBody()
        self.initSubjectsParser()
        self.initFTETabulator()
        self.initBudget()
        
        #################################################################
        
        #>>>>>>>>>>>>>>>>>>>>  connect callbacks <<<<<<<<<<<<<<<<<<<<<<<<
        [self.addPeriod() for i in range(5)]
        

        
        self.connectInit()
        self.connectButtons()
        self.connectPrepPossible()
        self.connectFreezePeriods()
        self.connectSpinners()

        #>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        
        
        #call any necessary updates for the first time
        self.initSuperteamTab()
        self.initTeamTab()
        self.initStudentAllocationTab()
        self.updateDayLength()
        
    def initStudentBody(self):
        """
        initialize the student body section of the ui. 
        it deals with all computations involving the number of students and how
        they are counted
        """
        
        #total enrollment        
        idx = range(1,6)
        self.TE = [self.ui.__dict__['spinBox_totalEnrollment_%d'%i] for i in idx]

        #self contained        
        self.SC = [self.ui.__dict__['spinBox_selfContained_%d'%i] for i in idx]

        #attendance rate
        self.AR = [self.ui.__dict__['spinBox_attendanceRate_%d'%i] for i in idx]

        #dropouts per year
        self.DPY = [self.ui.__dict__['spinBox_dropoutsPerYear_%d'%i] for i in idx]
        
        self.EE = [self.ui.__dict__['lcd_effectiveEnrollment_%d'%i] for i in idx]
        
        self.nS = [self.ui.__dict__['lcd_nSections_%d'%i] for i in idx]
        
        self.CS = [self.ui.__dict__['spinBox_classSize_%d'%i] for i in idx]
        
        self.studentBody = studentBody(self)
        self.studentBody.setup(self.TE, 
                                       self.SC, 
                                       self.AR, 
                                       self.DPY, 
                                       self.EE, 
                                       self.nS, 
                                       self.CS,
                                       self.nSValueChangedButton)
        S = [self.TE, self.SC, self.AR, self.DPY, self.CS]
        
        
        #connect the update callbacks to all the spinboxes in the object
        for i in range(5):
            for s in S:
                conn(s[i], sig("valueChanged(int)"), self.studentBody.updateEffectiveEnrollment)
                conn(s[i], sig("valueChanged(int)"), self.studentBody.updateNSections)
        self.studentBody.updateEffectiveEnrollment()
        self.studentBody.updateNSections()
        
        conn(self.ui.actionSave_to_excel, sig("triggered(bool)"), self.studentBody.write)
        conn(self.ui.button_save, sig("pressed()"), self.studentBody.write)
    
    def initSubjectsParser(self):
        self.subParser = subjectsParser(self)
        self.subParser.setup(#self.sp()[-4:],
                 self.ui.table_courseRequirements, 
                 self.ui.table_courseRequirements_credits,
                 self.periodLengths,
                 self.periodLengthsCredits)
        SP = self.sp()[-4:]
        
        self.subParser.addPeriod([sp.widget(0) for sp in SP])
        conn(self.ui.startTime, sig("timeChanged(QTime)"), 
             self.subParser.textChanged)
        conn(self.ui.stopTime, sig("timeChanged(QTime)"), 
             self.subParser.textChanged)
             
        conn(self.ui.button_save, sig("pressed()"), self.subParser.write)
        
            
            
                 
    def initFTETabulator(self):
        self.fteTab= FTETabulator(self)
        self.fteTab.setup(self.ui.table_FTE, 
                 self.ui.table_courseRequirements, 
                 self.studentBody,
                 self.ui.lcd_totalFTE,
                 self.periodLengths,
                 self.prepsPerDay,
                 [self.nSValueChangedButton])
        
        conn(self.TEValueChangedButton, sig('pressed()'),self.fteTab.fteChanged)
        conn(self.nSValueChangedButton, sig('pressed()'),self.fteTab.fteChanged)
        conn(self.ui.button_save, sig("pressed()"), self.fteTab.write)
        
    def initBudget(self):
        K = 'certifiedNonTeaching regularTeaching classifiedStaff'.split()
        BT = [self.ui.__dict__['tableWidget_'+k] for k in K]
        AB = [self.ui.__dict__['button_'+k+'AddLine'] for k in K]
        DB = [self.ui.__dict__['button_'+k+'DeleteLine'] for k in K]
        spD = [self.ui.__dict__['spinBox_'+k+'DeleteLine'] for k in K]
        lcdD = [self.ui.__dict__['lcdNumber_'+k+'Dollars'] for k in K]
        lcdSTU = [self.ui.__dict__['lcdNumber_'+k+'STU'] for k in K]
        
#        [adjustTableColumns(tab, 120) for tab in BT]
        

        self.budgetBoxes = []
        self.budgetBoxKeys = K
        for i, bt in enumerate(BT):
            self.budgetBoxes += [budgetTable()]
            self.budgetBoxes[i].setup(bt, AB[i], DB[i],spD[i],
                                            lcdD[i],lcdSTU[i],
                                            dummyDollarsPerSTUFun,
                                            self.buttonBudgetChanged)
                                            
        self.budgetTab = budgetTab(self)
        
        self.budgetTab.setup(self.budgetBoxes,
                             self.ui.doubleSpinBox_meanTeacherSalary,
                             self.ui.doubleSpinBox_nonPersonnelCostPct,
                             self.ui.tableWidget_startingBudget,
                             self.buttonBudgetChanged,
                             self.ui.lcdNumber_teachingSTU,
                             self.ui.lcdNumber_totalSTU,
                             [self.ui.lcdNumber_startingBudgetNonPersonellDollars,self.ui.lcdNumber_startingBudgetNonPersonellSTU],
                             [self.ui.lcdNumber_startingBudgetStaffingDollars, self.ui.lcdNumber_startingBudgetStaffingSTU],
                             [lcdD[1], lcdSTU[1]],
                             [lcdD[0], lcdSTU[0]],
                             [lcdD[2], lcdSTU[2]],
                             [self.ui.lcdNumber_finalBudgetDollars, self.ui.lcdNumber_finalBudgetSTU])
        
                             
        conn(self.ui.button_save, sig("pressed()"), self.budgetTab.write)
                             
    def initSuperteamTab(self):
        
        self.STArr = superteamArray(self.ui.frame_superteams)
        self.STArr.setup(0, self.nPeriods)
        conn(self.buttonNPeriodsChanged, sig("pressed()"), self.STArr.updateNTeams)
        
        self.table_teacherPool = teacherPool(self.ui.frame_teacherPool)
        self.table_teacherPool.setup(self.fteTab, self.STArr)
    def initTeamTab(self):
        self.TArray = teamArray(self.ui.teamTab)
        self.TArray.setup(self.STArr, self.table_teacherPool)
    
    def initStudentAllocationTab(self):
        self.SPool = studentPool(self.ui.studentAllocationTab)
        self.SATable = studentAllocationTable(self.ui.studentAllocationTab)
        
        self.SPool.setup(self.studentBody, self.SATable, self.nPeriods)
        self.SATable.setup(self.SPool, self.TArray, self.nPeriods)
        conn(self.SATable.button_cleared, sig('pressed()'), self.SPool.SATCleared)

    def connectInit(self):
        conn(self.ui.spinBox_prepsPerDay, sig("valueChanged(int)"), 
             self.updateTeachingEfficiency)
        conn(self.ui.startTime, sig("timeChanged(QTime)"), 
             self.runtimeChanged)
        conn(self.ui.stopTime, sig("timeChanged(QTime)"), 
             self.runtimeChanged)
        conn(self.ui.spinBox_lunchtime, sig("valueChanged(int)"), 
             self.runtimeChanged)
        conn(self.ui.spinBox_passingPeriod, sig("valueChanged(int)"), 
             self.runtimeChanged)
        
        conn(self.TEValueChangedButton, sig('pressed()'),self.fteTab.fteChanged)
        
        self.runtimeChanged()

        self.connectSplitters()
        self.connectPrepPossible()
        self.connectFreezePeriods()
        
        idx = range(1,6)
        nS = [self.ui.__dict__['lcd_nSections_%d'%i] for i in idx]
        [conn(self.nSValueChangedButton, sig("pressed()"), self.fteTab.fteChanged) for w in nS]
        
    
    
    
    def connectSplitters(self):
        for i,sp in enumerate(self.sp()):
            conn(sp, sig("splitterMoved(int,int)"),self.__getattribute__('sp%d'%i))
            conn(sp, sig("splitterMoved(int,int)"),self.updateSpinners)    
            conn(sp, sig("splitterMoved(int,int)"),self.subParser.textChanged)
   
    def connectSpinners(self):
        sp = self.ui.splitter_periodLength
        for i in range(sp.count()):
            conn(sp.widget(i), sig("valueChanged(int)"), self.spinnerChanged)
            conn(sp.widget(i), sig("valueChanged(int)"), self.subParser.textChanged)
    def blockSpinners(self, b = True):
        sp = self.ui.splitter_periodLength
        sp.blockSignals(b)
        [sp.widget(i).blockSignals(b) for i in range(sp.count())]
            
    def connectButtons(self):
        conn(self.ui.button_addPeriod, sig("clicked()"), self.addPeriod)
        conn(self.ui.button_removePeriod, sig("clicked()"), self.removePeriod)
        conn(self.ui.button_evenPeriods, sig("clicked()"), self.evenPeriods)
        
    def connectPrepPossible(self):
        sp = self.ui.splitter_prepPossible
        for i in range(sp.count()):
            conn(sp.widget(i), sig("stateChanged(int)"), self.updateTeachingEfficiency)
    def connectFreezePeriods(self):
#        sp = self.ui.splitter_freezePeriods
#        for i in range(sp.count()):
#            conn(sp.widget(i), sig("stateChanged(int)"), self.__getattribute__('fr%d'%i))
        pass
        
   
    
    def addPeriod(self):
        self.updateMaxPrepPeriodsPerDay(self.nPeriods()+1)        
        newWidgets = [f(self.sp()[i]) for i,f in enumerate(self.newPeriodFuns())]
        font = QtGui.QFont()
        font.setPointSize(14)
        font.setWeight(75)
        font.setBold(True)
        [w.setFont(font) for w in newWidgets]
        
        self.connectSpinners()
        self.connectPrepPossible()
        
        self.evenPeriods()
        
        TB = newWidgets[-4:]
        self.subParser.addPeriod(TB)
#        self.initSubjectsParser()
        self.cacheTotalSize()
        self.nPeriodsChanged()
        
    def removePeriod(self):
        if self.nPeriods()>1:
            self.updateMaxPrepPeriodsPerDay(self.nPeriods()-1)
            cnt = self.sp()[2].count()-1
            [sp.widget(cnt).deleteLater() for sp in self.sp()]
        self.timerID += [self.startTimer(50)]
        
    def nPeriodsChanged(self):
        self.buttonNPeriodsChanged.emit(sig('pressed()'))
        
    def timerEvent(self, E):
        [self.killTimer(tid) for tid in self.timerID] 
        self.timerID = []
        self.evenPeriods()
#        self.initSubjectsParser()
        self.subParser.removePeriod()
        self.subParser.textChanged()
        self.cacheTotalSize()
        self.nPeriodsChanged()
            
    def periodDestroyed(self):
        """
        need to perform some cleanup once we actually destroy the widget...
        """
#        print "period destroyed!"
        self.evenPeriods()
#        print "periods evened"
#        self.initSubjectsParser()
        self.cacheTotalSize()
        
    
    def evenPeriods(self):
        active = [not x for x in self.isFrozen()]
        if not any(active):
            return
        idx = find(active)
        total = sum([self.sizes()[i] for i in idx])
        
        even = round(float(total)/sum(active))
        newSizes = self.sizes()
        for i in idx:
            newSizes[i] = int(even)
        newSizes[idx[-1]] = int(total-int(even)*(sum(active)-1))
        self.updateSplitters(newSizes)
        self.updateSpinners()
        
    def freezePeriod(self, idx):
        W = [sp.widget(idx) for sp in self.sp()]
        [self.freezeWidget(w) for w in W]
    def unfreezePeriod(self, idx):
        W = [sp.widget(idx) for sp in self.sp()]
        [self.unfreezeWidget(w) for w in W]
            

    def splitterMoved(self, sp):
        self.updateSplitters(sp.sizes())
        self.updateSpinners()
        self.subParser.textChanged()
        self.cacheSizes()

    def spinnerChanged(self, val = None):
#        print "spinnerChanged: %d"%self.spinCount
        self.spinCount+=1
        self.blockSpinners(True)

        sp = self.ui.splitter_periodLength
        
        spCheck = self.ui.splitter_freezePeriods
        sb = self.spinnerVals()
        P = self.periodLengths()
        chidx = find([sb[i]!=p for i,p in enumerate(P)])
        
        if chidx==[]:
            return
        idx = chidx[0]
        d = sb[idx]-P[idx]
        spCheck.widget(idx).setChecked(True)
        
        active = [not x for x in self.isFrozen()]
        aidx = find(active)

        if aidx==[]:
            sp.widget(idx).setValue(P[idx])
            spCheck.widget(idx).setChecked(False)
            return
            
        sb[aidx[0]]-= d
        
        newSizes = [self.totalSize*float(s)/sum(sb) for s in sb]
        
        self.updateSplitters(newSizes)
        self.updateSpinners()
        self.blockSpinners(False)
        
    
        
    def updateSplitters(self, newSizes = None):
        if newSizes ==None:
            newSizes = self.ui.splitter_periodLength.sizes()
        for sp in self.sp():
            sp.blockSignals(True)
            sp.setSizes(newSizes)
            sp.blockSignals(False)
        self.cacheSizes()
        self.updateTeachingEfficiency()
        
    def updateSpinners(self):
        P = self.periodLengths()
        self.blockSpinners(True)
        sp = self.ui.splitter_periodLength
        [sp.widget(i).setValue(p)  for i,p in enumerate(P)]
        self.blockSpinners(False)

    def updateDayLength(self, time = None):
        self.ui.lcd_duration.display(self.dayLength())
        self.updateSpinners()
    def updateMaxPrepPeriodsPerDay(self,mx = None):
        if mx==None:
            mx = self.nPeriods()
        sb = self.ui.spinBox_prepsPerDay
        if sb.value()>mx:
            sb.setValue(mx)
        self.ui.spinBox_prepsPerDay.setMaximum(mx)
        
    def updateTeachingEfficiency(self):
        self.ui.lcd_teachingEfficiency.display(self.teachingEfficiency())
        self.TEValueChangedButton.emit(sig('pressed()'))
    
        
        
    def runtimeChanged(self):
        self.updateDayLength()
        self.updateClasstime()
        self.updateSpinners()
        self.updateTeachingEfficiency()
    def updateClasstime(self):
        self.ui.lcd_classtime.display(self.classTime())
        
    def cacheSizes(self):
        self.lastSize = self.sizes()
    def cacheTotalSize(self):
        self.totalSize = sum(self.sizes())

    def sp0(self, pos = None, idx = None):
        self.splitterMoved(self.sp()[0])
    def sp1(self, pos = None, idx = None):
        self.splitterMoved(self.sp()[1])
    def sp2(self, pos = None, idx = None):
        self.splitterMoved(self.sp()[2])
    def sp3(self, pos = None, idx = None):
        self.splitterMoved(self.sp()[3])
    def sp4(self, pos = None, idx = None):
        self.splitterMoved(self.sp()[4])
    def sp5(self, pos = None, idx = None):
        self.splitterMoved(self.sp()[5])
    def sp6(self, pos = None, idx = None):
        self.splitterMoved(self.sp()[6])
        
    def fr0(self, state = None):
        if state:
            self.freezePeriod(0)
        else:
            self.unfreezePeriod(0)
    def fr1(self, state = None):
        if state:
            self.freezePeriod(1)
        else:
            self.unfreezePeriod(1)
    def fr2(self, state = None):
        if state:
            self.freezePeriod(2)
        else:
            self.unfreezePeriod(2)
    def fr3(self, state = None):
        if state:
            self.freezePeriod(3)
        else:
            self.unfreezePeriod(3)
    def fr4(self, state = None):
        if state:
            self.freezePeriod(4)
        else:
            self.unfreezePeriod(4)
    def fr5(self, state = None):
        if state:
            self.freezePeriod(5)
        else:
            self.unfreezePeriod(5)
    def fr6(self, state = None):
        if state:
            self.freezePeriod(6)
        else:
            self.unfreezePeriod(6)
   
    def sp(self):
        return [self.ui.splitter_prepPossible, self.ui.splitter_freezePeriods,
                self.ui.splitter_periodLength,
                self.ui.splitter_subjects_9,self.ui.splitter_subjects_10,
                self.ui.splitter_subjects_11,self.ui.splitter_subjects_12]
                
    def newPeriodFuns(self):
        return [self.newPrepPossible,QtGui.QCheckBox,
                self.newPeriodLength,
                QtGui.QPlainTextEdit,QtGui.QPlainTextEdit,
                QtGui.QPlainTextEdit,QtGui.QPlainTextEdit]
    def newPrepPossible(self,parent):
        cb = QtGui.QCheckBox(parent)
        cb.setChecked(True)
        return cb
    def newPeriodLength(self,parent):
        sb = QtGui.QSpinBox(parent)
        sb.setMaximum(999)
        return sb
    
    def freezeWidget(self, w):
        w.setSizePolicy(frozenPolicy)
    def unfreezeWidget(self, w):
        w.setSizePolicy(unfrozenPolicy)
        
    def nPeriods(self):
        return self.sp()[2].count()
    def sizes(self):
        return self.sp()[2].sizes()
    def periodLengths(self):
        return self.sizes2SpinVals()
    def periodLengthsCredits(self):
        return [round(2*p/(60*self.hoursPerCredit()))/2.0 for p in self.periodLengths()]
    def spinnerVals(self):
        sp = self.ui.splitter_periodLength
        return [sp.widget(i).value() for i in range(sp.count())]
    def sizes2SpinVals(self, sizes = None):
        if sizes == None:
            sizes = self.sizes()
        return [int(round(self.classTimeMin()*sz/max([1,float(sum(sizes))]))) for sz in sizes]

    def isFrozen(self):
        sp = self.ui.splitter_freezePeriods        
        return [sp.widget(i).isChecked() for i in range(self.nPeriods())]
        
    def prepPossible(self):
        sp = self.ui.splitter_prepPossible
        return [sp.widget(i).isChecked() for i in range(sp.count())]
    
    def meanPrepLengthMin(self):
        pp = self.prepPossible()
        P = self.periodLengths()
        active = find(pp)
        if active==[]:
            return 0
        return int(round(float(sum([P[i] for i in active]))/sum(pp)))
    def prepsPerDay(self):
        return self.ui.spinBox_prepsPerDay.value()
        
    def teachingEfficiency(self):
        t = (self.classTimeMin()-self.meanPrepLengthMin()*self.prepsPerDay())/60.0
        e = 100*t/self.classTime()
        return int(e)
    def hoursPerCredit(self):
        m = min(self.periodLengths())
        return m/60.0
        
        
    def dayLength(self):
        start = self.ui.startTime.time()
        stop = self.ui.stopTime.time()
        return stop.hour()+stop.minute()/60.0-start.hour()-start.minute()/60.0 
        
    def dayLengthMin(self):
        return self.dayLength()*60
        
    def functionalTime(self):
        return (self.ui.spinBox_lunchtime.value()
                    +self.ui.spinBox_passingPeriod.value()*self.nPeriods())/60.0
                    
    def classTime(self):
        return self.dayLength()-self.functionalTime()
        
    def classTimeMin(self):
        return self.classTime()*60
        
    def teachingTime(self):
        t = (self.classTimeMin()-self.meanPrepLengthMin()*self.ui.spinBox_prepsPerDay.value())/60.0
        return t
        
    def saveName(self):
        return str(self.ui.plainTextEdit_saveName.toPlainText()).strip().split()[0]
    
  
    
    

if __name__ == "__main__":
    app = QtGui.QApplication(sys.argv)
    myapp = startMasterfulScheduling()
    myapp.show()
    sys.exit(app.exec_())