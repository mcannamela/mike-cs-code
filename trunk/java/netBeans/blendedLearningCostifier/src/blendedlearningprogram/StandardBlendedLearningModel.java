/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package blendedlearningprogram;

/**
 *
 * @author wichtelwesen
 */
public class StandardBlendedLearningModel implements BlendedLearningModelInterface {

    public StandardBlendedLearningModel() {
    }
    
    @Override
    public int getStudentToTeacherRatio() {
        return 25;
    }
    
}
