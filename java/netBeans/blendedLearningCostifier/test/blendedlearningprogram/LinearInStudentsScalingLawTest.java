/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package blendedlearningprogram;

import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;
import static org.junit.Assert.*;

/**
 *
 * @author Michael
 */
public class LinearInStudentsScalingLawTest {
    
    private static ProgramSize programSize;
    
    private static BaseBlendedLearningModel blendedLearningModel;
    private static int nrStudents;
    private static int nrPeriods;
    
    public LinearInStudentsScalingLawTest() {
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
        LinearInStudentsScalingLaw instance = new LinearInStudentsScalingLaw(programSize);
        int expResult = 3000;
        int result = instance.scale(cost);
        assertEquals(expResult, result);
        
    }
}
