/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package blendedlearningprogram;

import org.junit.AfterClass;
import org.junit.BeforeClass;
import org.junit.Test;
import static org.junit.Assert.*;

/**
 *
 * @author wichtelwesen
 */
public class StandardBlendedLearningModelTest {
    
    public StandardBlendedLearningModelTest() {
    }
    

    /**
     * Test of getStudentToTeacherRatio method, of class StandardBlendedLearningModel.
     */
    @Test
    public void testGetStudentToTeacherRatio() {
        System.out.println("getStudentToTeacherRatio");
        StandardBlendedLearningModel instance = new StandardBlendedLearningModel();
        int expResult = 25;
        int result = instance.getStudentToTeacherRatio();
        assertEquals(expResult, result);
        
    }
}
