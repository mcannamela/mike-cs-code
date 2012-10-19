/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package blendedlearningprogram;

/**
 *
 * @author mcannamela
 */
public class ExtremeBlendedLearningModel implements BlendedLearningModelInterface {
    @Override
    public int getStudentToTeacherRatio() {
        return 1000;
    }

    @Override
    public String getDescription() {
        return "This is a model to serve an extreme number of students, for instance by webcasting a lecture.";
    }
}
