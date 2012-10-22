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
    private BaseBlendedLearningModel blendedLearningModel;
    private int nrStudents;
    private int nrPeriods;

    public ProgramSize() {
        this.blendedLearningModel = new StandardBlendedLearningModel();
        this.nrStudents = 100;
        this.nrPeriods = 8;
    }

    public ProgramSize(BaseBlendedLearningModel blendedLearningModel, int nrStudents, int nrPeriods) {
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

    public BaseBlendedLearningModel getBlendedLearningModel() {
        return blendedLearningModel;
    }
    
 
    
    public void setNrStudents(int nrStudents) {
        this.nrStudents = nrStudents;
    }
    public void setNrPeriods(int nrPeriods) {
        this.nrPeriods = nrPeriods;
    }

    public void setBlendedLearningModel(BaseBlendedLearningModel blendedLearningModel) {
        this.blendedLearningModel = blendedLearningModel;
    }
    
    public int getNrSections(){
        return (int) Math.ceil((double)getNrStudents()/
                       (double)blendedLearningModel.getStudentToTeacherRatio());
    }
    public int getNrTeachers(){
        return (int) Math.ceil((double)getNrStudents()/
                       (double)(blendedLearningModel.getStudentToTeacherRatio()*
                                getNrPeriods()));
    }
    public int getNrSimultaneousSections(){
        return (int) Math.ceil((double)getNrSections()/
                       (double)getNrPeriods());
    }
    
}
