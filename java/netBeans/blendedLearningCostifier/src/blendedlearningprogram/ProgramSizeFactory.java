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
public class ProgramSizeFactory {
    private BlendedLearningModelFactory modelFactory= new BlendedLearningModelFactory();
    public static final String NR_STUDENTS_KEY = "nrStudents";
    public static final String NR_PERIODS_KEY = "nrPeriods";
    
    public ProgramSize makeProgramSize(Configuration config){
        BaseBlendedLearningModel model;
        int nrStudents, nrPeriods; 
                
        nrPeriods = config.getInt(NR_PERIODS_KEY, 8);
        nrStudents= config.getInt(NR_STUDENTS_KEY, 100);
        model = modelFactory.makeBlendedLearningModel(config);
        
        return new ProgramSize(model, nrStudents, nrPeriods);
        
    }
    
}
