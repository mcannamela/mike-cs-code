/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package blendedlearningprogram;

import java.util.ArrayList;

/**
 *
 * @author mcannamela
 */
public enum ScalingLawTypeEnum {
    CONSTANT ,
        LINEAR_IN_STUDENTS,
        LINEAR_IN_TEACHERS;
        
    public ArrayList<String> getScalingLawTypes(){
        ArrayList<String> types = new ArrayList();
        for(ScalingLawTypeEnum type: ScalingLawTypeEnum.values()){
            types.add(type.toString());
        }
        return types;
    }
}
