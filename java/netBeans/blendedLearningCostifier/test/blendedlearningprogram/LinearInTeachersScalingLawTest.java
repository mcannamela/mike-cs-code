/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package blendedlearningprogram;

import org.junit.After;
import static org.junit.Assert.*;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;

/**
 *
 * @author Michael
 */
public class LinearInTeachersScalingLawTest {
    private static ProgramSize programSize;
    
    private static BlendedLearningModelInterface blendedLearningModel;
    private static int nrStudents;
    private static int nrPeriods;
    
    public LinearInTeachersScalingLawTest() {
    }
    
    @BeforeClass
    public static void setUpClass() {
        blendedLearningModel = new StandardBlendedLearningModel();
        nrStudents = 300;
        nrPeriods = 6;
        programSize = new ProgramSize(blendedLearningModel, nrStudents, nrPeriods);
    }

    

    /**
     * Test of scale method, of class LinearInTeachersScalingLaw.
     */
    @Test
    public void testScale() {
        System.out.println("scale");
        int cost = 10;
        LinearInTeachersScalingLaw instance = new LinearInTeachersScalingLaw(programSize);
        int expResult = 20;
        int result = instance.scale(cost);
        assertEquals(expResult, result);
        
    }
}
