/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package blendedlearningprogram;

/**
 *
 * @author mcannamela
 */
public class ExtremeBlendedLearningModel extends BaseBlendedLearningModel {
    
    public ExtremeBlendedLearningModel(){
        description = "This is a model to serve an extreme number of students, for instance by webcasting a lecture.";
        type = BlendedLearningModelEnum.EXTREME.toString();
        studentToTeacherRatio = 1000;
    }
            
    
}
