/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package blendedlearningprogram;

import org.apache.commons.configuration.Configuration;

/**
 *
 * @author mcannamela
 */
public class BlendedLearningModelFactory {
    
    private static final String TYPE_KEY = "type";
    
    public BlendedLearningModelInterface makeBlendedLearningModel(Configuration config){
        String type;
        type = config.getString(TYPE_KEY, BlendedLearningModelEnum.STANDARD.toString());
        return makeBlendedLearningModel(type);
    }
    public BlendedLearningModelInterface makeBlendedLearningModel(String type){
        BlendedLearningModelInterface model;
        if (type.equals(BlendedLearningModelEnum.ONE_ON_ONE.toString())){
            model =  new OneOnOneBlendedLearningModel();
        }
        else if (type.equals(BlendedLearningModelEnum.STANDARD.toString())){
            model =  new StandardBlendedLearningModel();
        }
        else if (type.equals(BlendedLearningModelEnum.EXTREME.toString())){
            model =  new ExtremeBlendedLearningModel();
        }
        else{
            throw new UnsupportedOperationException("blended learning model is not of a known type.");
        }
        return model;
    }
    
}
