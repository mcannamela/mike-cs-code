/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package blendedlearningprogram;

/**
 *
 * @author mcannamela
 */
public class ProgramSize {
    private BlendedLearningModelInterface blendedLearningModel;
    private int nrStudents;
    private int nrPeriods;

    public ProgramSize() {
        this.blendedLearningModel = new StandardBlendedLearningModel();
        this.nrStudents = 100;
        this.nrPeriods = 8;
    }

    public ProgramSize(BlendedLearningModelInterface blendedLearningModel, int nrStudents, int nrPeriods) {
        this.blendedLearningModel = blendedLearningModel;
        this.nrStudents = nrStudents;
        this.nrPeriods = nrPeriods;
    }

    public int getNrPeriods() {
        return nrPeriods;
    }
    public int getNrStudents() {
        return nrStudents;
    }

    public BlendedLearningModelInterface getBlendedLearningModel() {
        return blendedLearningModel;
    }
    
 
    
    public void setNrStudents(int nrStudents) {
        this.nrStudents = nrStudents;
    }
    public void setNrPeriods(int nrPeriods) {
        this.nrPeriods = nrPeriods;
    }

    public void setBlendedLearningModel(BlendedLearningModelInterface blendedLearningModel) {
        this.blendedLearningModel = blendedLearningModel;
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
