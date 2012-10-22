/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package blendedlearningprogram;

/**
 *
 * @author mcannamela
 */
public class OneOnOneBlendedLearningModel extends BaseBlendedLearningModel {

    public OneOnOneBlendedLearningModel() {
        description = "This is a very inefficent way to blend learning!";
        type = BlendedLearningModelEnum.ONE_ON_ONE.toString();
        studentToTeacherRatio = 1;
    }
    

   
}
