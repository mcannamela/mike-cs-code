/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package costoptiontree;

import blendedlearningprogram.*;
import org.apache.commons.configuration.*;
/**
 *
 * @author wichtelwesen
 */
public class CostOptionFactory {
    private ProgramSize programSize;
    

    
    private static final String labelKey = "label"; 
    private static final String descriptionKey = "description";
    private static final String minCostKey = "minCost";
    private static final String maxCostKey = "maxCost";
    private static final String selectedCostKey = "selectedCost";
    private static final String scalingLawTypeKey = "scalingLawType";

    public CostOptionFactory() {
        this.programSize = new ProgramSize();
    }
    
    public CostOptionFactory(ProgramSize programSize) {
        this.programSize = programSize;
    }
    
    
    
    
    public PropertiesConfiguration defaultConfig(){
        PropertiesConfiguration config = new PropertiesConfiguration();
        config.addProperty(labelKey,  "optionLabel"); 
        config.addProperty(descriptionKey, "this is a description of the option, "
                                            +"which is very long and will probably be "
                                            +"displayed on multiple lines");
        config.addProperty(minCostKey, 10);
        config.addProperty(maxCostKey, 100);
        config.addProperty(selectedCostKey,60);
        config.addProperty(scalingLawTypeKey, ScalingLawTypeEnum.LINEAR_IN_STUDENTS.toString());
                
        return config;
    }
    public CostOption makeCostOption(){
        return makeCostOption(defaultConfig());
    }
    /**
     *
     * @param config the configuration object, presumably from a file, that 
     * describes this cost option.
     * @return returns the configured CostOption object
     */
    public CostOption makeCostOption(AbstractConfiguration config){
        String label = config.getString(labelKey);
        config.setListDelimiter('0');
        
        String[] descriptions = config.getStringArray(descriptionKey);
        
        String description = CostOptionNodeFactory.reconstituteCommaConfigString(descriptions);
        
        
        config.setListDelimiter(',');
        int minCost = config.getInt(minCostKey);
        int maxCost = config.getInt(maxCostKey);
        int selectedCost = config.getInt(selectedCostKey, (minCost+maxCost)/2 );
        
        
        String scalingLawType = config.getString(scalingLawTypeKey);
              
        AbstractScalingLaw scalingLaw;
        if (scalingLawType.equals(ScalingLawTypeEnum.CONSTANT.toString())){
            scalingLaw = new ConstantScalingLaw(programSize);
        }       
        else if (scalingLawType.equals(ScalingLawTypeEnum.LINEAR_IN_STUDENTS.toString())){
            scalingLaw = new LinearInStudentsScalingLaw(programSize);
        }
         
        else if (scalingLawType.equals(ScalingLawTypeEnum.LINEAR_IN_TEACHERS.toString())){
            scalingLaw = new LinearInTeachersScalingLaw(programSize);
        }
        else {
            throw new UnsupportedOperationException("scaling law is not of a known type.");
        }
        return new CostOption(label,  description, minCost, maxCost,
                         selectedCost,  scalingLaw);
    }
    
}
