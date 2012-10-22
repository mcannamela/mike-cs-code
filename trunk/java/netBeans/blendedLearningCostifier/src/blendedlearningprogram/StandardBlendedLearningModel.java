/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package blendedlearningprogram;

/**
 *
 * @author wichtelwesen
 */
public class StandardBlendedLearningModel extends BaseBlendedLearningModel {

    public StandardBlendedLearningModel() {
        description = "This is a standard schmandard model with a standard student/teacher ratio";
        type = BlendedLearningModelEnum.STANDARD.toString();
        studentToTeacherRatio = 25;
    }
    
    
    
   
    
}
