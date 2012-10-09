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
    
    private enum ScalingLawTypes{
        CONSTANT ,
        LINEAR_IN_STUDENTS,
        LINEAR_IN_TEACHERS;
    }

    public CostOptionFactory(ProgramSize programSize) {
        this.programSize = programSize;
    }
    
    
    
    /**
     *
     * @param config the configuration object, presumably from a file, that 
     * describes this cost option.
     * @return returns the configured CostOption object
     */
    public CostOption makeCostOption(Configuration config){
        String label = config.getString("label");
        String description = config.getString("description");
        int minCost = config.getInt("minCost");
        int maxCost = config.getInt("maxCost");
        int selectedCost = config.getInt("selectedCost", (minCost+maxCost)/2 );
        
        ScalingLawTypes scalingLawTypes;
        String scalingLawType = config.getString("scalingLawType");
              
        AbstractScalingLaw scalingLaw;
        if (scalingLawType.equals(ScalingLawTypes.CONSTANT.toString())){
            scalingLaw = new ConstantScalingLaw(programSize);
        }       
        else if (scalingLawType.equals(ScalingLawTypes.LINEAR_IN_STUDENTS.toString())){
            scalingLaw = new LinearInStudentsScalingLaw(programSize);
        }
         
        else if (scalingLawType.equals(ScalingLawTypes.LINEAR_IN_TEACHERS.toString())){
            scalingLaw = new LinearInTeachersScalingLaw(programSize);
        }
        else {
            throw new UnsupportedOperationException("scaling law is not of a known type.");
        }
        return new CostOption(label,  description, minCost, maxCost,
                         selectedCost,  scalingLaw);
    }
    
}
