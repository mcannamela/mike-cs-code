/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package blendedlearningprogram;
import java.lang.Math;

/**
 *
 * @author mcannamela
 */
public class ProgramSize {
    private BlendedLearningModelInterface blendedLearningModel;
    private int nrStudents;
    private int nrPeriods;
    

    public int getNrPeriods() {
        return nrPeriods;
    }
    public int getNrStudents() {
        return nrStudents;
    }
 
    
    public void setNrStudents(int nrStudents) {
        this.nrStudents = nrStudents;
    }
    public void setNrPeriods(int nrPeriods) {
        this.nrPeriods = nrPeriods;
    }
    
    
    public int getNrTeachers(){
        return (int) Math.ceil((double)nrStudents/
                       (double)blendedLearningModel.getStudentToTeacherRatio());
    }
    public int getNrSimultaneousSections(){
        return (int) Math.ceil((double)getNrTeachers()/
                       (double)getNrPeriods());
    }
    
}
