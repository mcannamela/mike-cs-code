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
public enum BlendedLearningModelEnum {
    STANDARD ,
        ONE_ON_ONE,
        EXTREME;
        
    public ArrayList<String> getScalingLawTypes(){
        ArrayList<String> types = new ArrayList();
        for(ScalingLawTypeEnum type: ScalingLawTypeEnum.values()){
            types.add(type.toString());
        }
        return types;
    }
}
