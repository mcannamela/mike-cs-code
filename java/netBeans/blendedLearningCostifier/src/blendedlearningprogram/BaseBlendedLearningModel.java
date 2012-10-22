/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package blendedlearningprogram;

/**
 *
 * @author wichtelwesen
 */
public class BaseBlendedLearningModel {
    protected  String description = "";
    protected String type;
    protected int studentToTeacherRatio;
   
    public int getStudentToTeacherRatio(){
        return studentToTeacherRatio;
    }
    public String getType(){
        return type;
    }
    
    public String getDescription() {
        return description;
    }

    public void setDescription(String description) {
        this.description = description;
    }

    

    
    
    
}
