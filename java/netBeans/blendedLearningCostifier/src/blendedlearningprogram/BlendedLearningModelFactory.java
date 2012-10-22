/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package blendedlearningprogram;

import costoptiontree.CostOptionNodeFactory;
import org.apache.commons.configuration.Configuration;

/**
 *
 * @author mcannamela
 */
public class BlendedLearningModelFactory {
    
    private static final String TYPE_KEY = "modelType";
    private static final String DESCRIPTION_KEY = "modelDescription";
    
    public BaseBlendedLearningModel makeBlendedLearningModel(Configuration config){
        String type, description;
        String[] descriptionTokens;
        type = config.getString(TYPE_KEY, BlendedLearningModelEnum.STANDARD.toString());
        descriptionTokens = config.getStringArray(DESCRIPTION_KEY);
        description = CostOptionNodeFactory.reconstituteCommaConfigString(descriptionTokens);
        return makeBlendedLearningModel(type);
    }
    public BaseBlendedLearningModel makeBlendedLearningModel(String type){
        BaseBlendedLearningModel model;
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
    
    public BaseBlendedLearningModel makeBlendedLearningModel(String type, String description){
        BaseBlendedLearningModel model;
        model = makeBlendedLearningModel(type);
        model.setDescription(description);
        return model;
    }
    
    
}
