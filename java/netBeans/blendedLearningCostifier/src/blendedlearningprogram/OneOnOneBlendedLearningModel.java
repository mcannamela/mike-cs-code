/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package blendedlearningprogram;

/**
 *
 * @author mcannamela
 */
public class OneOnOneBlendedLearningModel implements BlendedLearningModelInterface {
    @Override
    public int getStudentToTeacherRatio() {
        return 1;
    }

    @Override
    public String getDescription() {
        return "This is a very inefficent way to blend learning!";
    }
}
